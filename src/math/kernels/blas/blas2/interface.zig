//! The level-2 interface: `interface/{gemv,zgemv,ger,zger,symv,zsymv,zhemv}.c`'s
//! checks, `beta` scaling and pointer set-up, as `blas1/interface.zig` takes
//! its arguments. Scratch is on the stack up to `stack_bytes` (`STACK_ALLOC`),
//! else allocated. Above OpenBLAS's thresholds the routines thread on `io`,
//! `symv` and `hemv` through partial vectors folded in chunk order.

const std = @import("std");
const sysinfo = @import("builtin");
const blas = @import("../__root__.zig");
const types = @import("../types.zig");
const half_precision = @import("../half_precision.zig");
const util = @import("../util.zig");
const sym = @import("symmetric.zig");
const portable = @import("../portable.zig");
const bf = @import("../ext/bf16/__root__.zig");
const accel = @import("../apple_accel.zig");
const threading = @import("../threading.zig");
const Ops = @import("ops.zig").Ops;
const Op = portable.Op;

const Allocator = std.mem.Allocator;
const Io = std.Io;
const BlasLong = types.BlasLong;
const Trans = types.Trans;
const Uplo = types.Uplo;

/// Bad arguments, by the reference's `xerbla` argument names; shared by every
/// level-2 interface file.
pub const ArgError = error{
    InvalidUplo,
    InvalidTrans,
    InvalidDiag,
    InvalidM,
    InvalidN,
    InvalidK,
    InvalidLda,
    InvalidIncx,
    InvalidIncy,
};
/// What a level-2 routine returns.
pub const Error = ArgError || Allocator.Error || Io.Cancelable;
/// What `gemvHalf`, which allocates nothing, returns.
pub const HalfError = ArgError || Io.Cancelable;

const isComplex = types.isComplex;
const Float = types.Float;

const raw = util.raw;

const first = util.first;

const isOne = util.isOne;

const isZero = util.isZero;

// ---- scratch buffer ----------------------------------------------------------

/// `MAX_STACK_ALLOC`, in every level-2 driver.
pub const stack_bytes = util.stack_bytes;
/// The reference's stack buffers are `aligned(0x20)`; no ported kernel needs
/// more than 16 (`buffers.md`).
pub const stack_align = 32;
/// `blas_memory_alloc`'s buffers are pages; the generic symv and hemv lay
/// their copies out on page boundaries from it.
pub const heap_align = 4096;

/// Scratch for `F` elements at the level-2 alignments.
fn Scratch(comptime F: type) type {
    return util.Scratch(F, stack_align, heap_align);
}

// ---- gemv --------------------------------------------------------------------

/// `?GEMV` (`interface/gemv.c`, `zgemv.c`): `y := alpha*op(A)*x + beta*y`,
/// `A` being `m` by `n`. `y` is scaled by `beta` first (a zero `beta`
/// overwrites NaN), and nothing more happens for a zero `alpha`.
pub fn gemv(allocator: Allocator, io: Io, comptime T: type, trans: Trans, m: isize, n: isize, alpha: T, a: [*]T, lda: isize, x: [*]T, incx: isize, beta: T, y: [*]T, incy: isize) Error!void {
    if (!isComplex(T) and trans.conjX()) return error.InvalidTrans;
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, m)) return error.InvalidLda;
    if (incx == 0) return error.InvalidIncx;
    if (incy == 0) return error.InvalidIncy;
    if (m == 0 or n == 0) return;
    const tr = trans.transposed();
    const lenx = if (tr) m else n;
    const leny = if (tr) n else m;
    if (!util.fitsMatrix(m, n, lda, @sizeOf(T))) return error.InvalidLda;
    if (!util.fits(lenx, incx, @sizeOf(T))) return error.InvalidIncx;
    if (!util.fits(leny, incy, @sizeOf(T))) return error.InvalidIncy;
    if (try accel.gemv(io, T, trans, m, n, alpha, a, lda, x, incx, beta, y, incy)) return;
    if (!isOne(T, beta)) Ops(T).init().scal(leny, beta, y, @intCast(@abs(incy)));
    if (isZero(T, alpha)) return;
    switch (trans) {
        inline else => |t| {
            // `r`..`d` have no real kernel; `gemv` returned above for those
            if (comptime !isComplex(T) and t.conjX()) unreachable;
            try gemvRun(allocator, io, T, t, m, n, alpha, a, lda, first(T, x, lenx, incx), incx, first(T, y, leny, incy), incy);
        },
    }
}

/// Chunks of a threaded level-2 routine start on multiples of this many rows
/// or columns, so a kernel's row or column tail falls where one call puts
/// it, and chunks of `y` or `A` stay off each other's cache lines.
const quantum = 64;

/// `gemv`'s row chunks without a transpose start on multiples of this:
/// x86_64's kernels take rows in blocks of up to 4096 (`NBMAX`), the last
/// block's remainder first (`gemv/driver.zig`), so only whole blocks split.
const gemv_rows = 4096;

fn gemvRun(allocator: Allocator, io: Io, comptime T: type, comptime t: Trans, m: isize, n: isize, alpha: T, a: [*]T, lda: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) (Allocator.Error || Io.Cancelable)!void {
    const F = Float(T);
    const complex = comptime isComplex(T);
    if (comptime !complex and t.conjX()) return;
    const op = comptime portable.gemvOp(if (complex) t else if (t.transposed()) .t else .n);
    const cs = if (complex) 2 else 1;
    // `gemv_thread.c` splits `y`: rows without a transpose, columns with one
    const tr = comptime t.transposed();
    const mu: usize = @intCast(m);
    const nu: usize = @intCast(n);
    const p = threading.plan(.gemv, T, threading.product(&.{ mu, nu }), if (tr) nu else mu, .even, if (tr) quantum else gemv_rows);
    if (p.workers > 1) {
        const G = struct {
            m: isize,
            n: isize,
            alpha: T,
            a: [*]T,
            lda: isize,
            x: [*]T,
            incx: isize,
            y: [*]T,
            incy: isize,

            fn dims(g: *const @This(), c: threading.Chunk) [2]usize {
                const k = c.end - c.start;
                return if (tr) .{ @intCast(g.m), k } else .{ k, @intCast(g.n) };
            }

            /// Each chunk's buffer on a cache line of its own.
            fn len(g: *const @This(), c: threading.Chunk) usize {
                const d = g.dims(c);
                return std.mem.alignForward(usize, sym.gemvLen(F, cs, d[0], d[1]), 64);
            }

            fn one(g: *const @This(), c: threading.Chunk, part: []F) void {
                const d = g.dims(c);
                const s: isize = @intCast(c.start);
                const pa: [*]T = @ptrCast(util.at(T, g.a, if (tr) s * g.lda else s));
                const py: [*]T = @ptrCast(util.at(T, g.y, s * g.incy));
                Ops(T).init().gemvK(op, @intCast(d[0]), @intCast(d[1]), g.alpha, pa, g.lda, g.x, g.incx, py, g.incy, part.ptr);
            }
        };
        const g: G = .{ .m = m, .n = n, .alpha = alpha, .a = a, .lda = lda, .x = x, .incx = incx, .y = y, .incy = incy };
        return threading.forEachPart(allocator, io, &p, F, heap_align, &g, G.len, G.one);
    }
    var stack: [stack_bytes]u8 align(stack_align) = undefined;
    const sc = try Scratch(F).init(allocator, &stack, sym.gemvLen(F, cs, mu, nu));
    defer sc.deinit(allocator);
    // the reference clears its stack buffer for every `trans` but `n` on
    // x86_64, for `cgemv_t.S` (<https://github.com/OpenMathLib/OpenBLAS/issues/746>,
    // `interface/zgemv.c:250`, v0.3.34-317), which no core here takes; no
    // kernel here lets what the buffer held reach `y` (`buffers.md`)
    Ops(T).init().gemvK(op, m, n, alpha, a, lda, x, incx, y, incy, sc.ptr);
}

// ---- half-precision gemv -----------------------------------------------------

/// Which half-precision `gemv` (`half_precision.zig`).
pub const HalfPrecision = half_precision.HalfPrecision;

/// `SBGEMV`, `BGEMV`, `SHGEMV` (`interface/sbgemv.c`): `y := alpha*op(A)*x +
/// beta*y`, `A` being `m` by `n`; `trans` is `n`/`r` or `t`/`c`. Unlike
/// `gemv`, the kernel applies `beta` itself, and a zero `alpha` only scales
/// `y` (when `beta` is not 1) without reading `A` or `x`. The kernels take no
/// scratch buffer, so there is no allocator; the reference's threaded
/// drivers are not here.
pub fn gemvHalf(io: Io, comptime kind: HalfPrecision, trans: Trans, m: isize, n: isize, alpha: kind.Out(), a: [*]kind.In(), lda: isize, x: [*]kind.In(), incx: isize, beta: kind.Out(), y: [*]kind.Out(), incy: isize) HalfError!void {
    _ = io;
    const In = kind.In();
    const Out = kind.Out();
    if (trans.conjX()) return error.InvalidTrans;
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, m)) return error.InvalidLda;
    if (incx == 0) return error.InvalidIncx;
    if (incy == 0) return error.InvalidIncy;
    if (m == 0 or n == 0) return;
    const tr = trans.transposed();
    const lenx = if (tr) m else n;
    const leny = if (tr) n else m;
    if (!util.fitsMatrix(m, n, lda, @sizeOf(In))) return error.InvalidLda;
    if (!util.fits(lenx, incx, @sizeOf(In))) return error.InvalidIncx;
    if (!util.fits(leny, incy, @sizeOf(Out))) return error.InvalidIncy;
    // `SBF16TOS_K` of the scalars, for bfloat16
    const alpha_f = halfToF32(Out, alpha);
    const beta_f = halfToF32(Out, beta);
    if (alpha_f == 0) {
        if (beta_f != 1) {
            const inc: BlasLong = @intCast(@abs(incy));
            if (comptime kind == .bf16) {
                _ = blas.kernel(.bscal_k)(leny, 0, 0, beta, y, inc, y, 0, y, 0);
            } else {
                _ = blas.kernel(.sscal_k)(leny, 0, 0, beta, y, inc, y, 0, y, 0);
            }
        }
        return;
    }
    const px = first(In, x, lenx, incx);
    const py = first(Out, y, leny, incy);
    _ = if (tr)
        blas.kernel(comptime portable.halfGemvSlot(kind.gemm(), true))(m, n, alpha, a, lda, px, incx, beta, py, incy)
    else
        blas.kernel(comptime portable.halfGemvSlot(kind.gemm(), false))(m, n, alpha, a, lda, px, incx, beta, py, incy);
}

/// A `gemvHalf` scalar as `f32`.
fn halfToF32(comptime T: type, v: T) f32 {
    return if (T == f32) v else bf.scalarToF32(v);
}

// ---- ger ---------------------------------------------------------------------

/// `?GER`, `?GERU` (`interface/ger.c`, `zger.c`): `A := alpha*x*y^T + A`.
pub fn ger(allocator: Allocator, io: Io, comptime T: type, m: isize, n: isize, alpha: T, x: [*]T, incx: isize, y: [*]T, incy: isize, a: [*]T, lda: isize) Error!void {
    return gerImpl(allocator, io, T, false, m, n, alpha, x, incx, y, incy, a, lda);
}

/// `?GERC`: `A := alpha*x*y^H + A`. Complex only.
pub fn gerc(allocator: Allocator, io: Io, comptime T: type, m: isize, n: isize, alpha: T, x: [*]T, incx: isize, y: [*]T, incy: isize, a: [*]T, lda: isize) Error!void {
    comptime util.requireComplex(T, "gerc");
    return gerImpl(allocator, io, T, true, m, n, alpha, x, incx, y, incy, a, lda);
}

fn gerImpl(allocator: Allocator, io: Io, comptime T: type, comptime conj: bool, m: isize, n: isize, alpha: T, x: [*]T, incx: isize, y: [*]T, incy: isize, a: [*]T, lda: isize) Error!void {
    const F = Float(T);
    const complex = comptime isComplex(T);
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    if (incx == 0) return error.InvalidIncx;
    if (incy == 0) return error.InvalidIncy;
    if (lda < @max(1, m)) return error.InvalidLda;
    if (m == 0 or n == 0) return;
    if (!util.fitsMatrix(m, n, lda, @sizeOf(T))) return error.InvalidLda;
    if (!util.fits(m, incx, @sizeOf(T))) return error.InvalidIncx;
    if (!util.fits(n, incy, @sizeOf(T))) return error.InvalidIncy;
    if (isZero(T, alpha)) return;
    if (try accel.ger(io, T, conj, m, n, alpha, x, incx, y, incy, a, lda)) return;
    const op: Op = comptime if (complex) (if (conj) .gerc_k else .geru_k) else .ger_k;
    const px = first(T, x, m, incx);
    const py = first(T, y, n, incy);
    var stack: [stack_bytes]u8 align(stack_align) = undefined;
    // a strided `x`'s copy; no kernel reads the buffer for a unit `x`
    // (`buffers.md`), where the reference may pass NULL
    const sc: Scratch(F) = if (incx == 1) .{ .ptr = @ptrCast(&stack), .len = 0 } else try Scratch(F).init(allocator, &stack, (if (complex) 2 else 1) * @as(usize, @intCast(m)));
    defer sc.deinit(allocator);
    // the reference passes over a column whose `y` is zero, so NaN in `x`
    // stays out of it; the kernels get the runs between zeros, a strided `x`
    // copied once for all of them
    var g: Ger(T, op) = .{ .m = m, .n = n, .alpha = alpha, .x = px, .incx = incx, .y = py, .incy = incy, .a = a, .lda = lda, .buf = sc.ptr };
    // `ger_thread.c` splits the columns
    const p = threading.plan(.ger, T, threading.product(&.{ @intCast(m), @intCast(n) }), @intCast(n), .even, quantum);
    if (p.workers > 1) {
        if (incx != 1) g.copyX();
        return threading.forEach(io, &p, @as(*const Ger(T, op), &g), Ger(T, op).chunk);
    }
    g.columns(0, n);
}

/// `ger`'s runs of columns, whole or in chunks.
fn Ger(comptime T: type, comptime op: Op) type {
    return struct {
        const Self = @This();

        m: isize,
        n: isize,
        alpha: T,
        x: [*]T,
        incx: isize,
        y: [*]T,
        incy: isize,
        a: [*]T,
        lda: isize,
        buf: [*]Float(T),

        fn nonzero(g: *const Self, j: isize) bool {
            return !isZero(T, util.at(T, g.y, j * g.incy).*);
        }

        /// `x` copied to the buffer, taken at unit stride from then on.
        fn copyX(g: *Self) void {
            _ = portable.fetch(T, .copy_k, blas.vtable)(g.m, raw(T, g.x), g.incx, g.buf, 1);
            g.x = @ptrCast(g.buf);
            g.incx = 1;
        }

        /// Where the part of the run from `r0` that starts at or past `c`
        /// begins: a multiple of `quantum` past `r0`, so a run's parts make
        /// its kernel call's column blocks.
        fn cut(c: isize, r0: isize) isize {
            if (c <= r0) return r0;
            return r0 + @divFloor(c - r0 + quantum - 1, quantum) * quantum;
        }

        fn chunk(g: *const Self, c: threading.Chunk) void {
            var h = g.*;
            h.columns(@intCast(c.start), @intCast(c.end));
        }

        /// The parts of every run that `cut` puts in columns `c0..c1`; `copyX`
        /// first when they are chunks of a call.
        fn columns(g: *Self, c0: isize, c1: isize) void {
            var j = c0;
            // the run `c0` is in starts before it
            if (j < g.n and g.nonzero(j)) {
                while (j > 0 and g.nonzero(j - 1)) j -= 1;
            }
            while (j < c1) {
                while (j < g.n and !g.nonzero(j)) j += 1;
                const r0 = j;
                while (j < g.n and g.nonzero(j)) j += 1;
                if (j == r0) break;
                const lo = @max(cut(c0, r0), r0);
                const hi = @min(cut(c1, r0), j);
                if (lo >= hi) continue;
                if (g.incx != 1) g.copyX();
                const ys: [*]T = @ptrCast(util.at(T, g.y, lo * g.incy));
                const as: [*]T = @ptrCast(util.at(T, g.a, lo * g.lda));
                const k = portable.fetch(T, op, blas.vtable);
                _ = @call(.auto, k, .{ g.m, hi - lo, 0 } ++ util.scalarArgs(T, g.alpha) ++ .{ raw(T, g.x), g.incx, raw(T, ys), g.incy, raw(T, as), g.lda, g.buf });
            }
        }
    };
}

// ---- symv, hemv --------------------------------------------------------------

/// `?SYMV` (`interface/symv.c`, `zsymv.c`): `y := alpha*S*x + beta*y` for the
/// symmetric `S` whose `uplo` triangle (`u` or `l`) is stored in `a`.
pub fn symv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, alpha: T, a: [*]T, lda: isize, x: [*]T, incx: isize, beta: T, y: [*]T, incy: isize) Error!void {
    if (uplo == .v or uplo == .m) return error.InvalidUplo;
    return mv(allocator, io, T, false, uplo, n, alpha, a, lda, x, incx, beta, y, incy);
}

/// `?HEMV` (`interface/zhemv.c`): the same for the Hermitian `S`; `v` and
/// `m` use the upper and lower triangles for `conj(S)`. Complex only.
pub fn hemv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, alpha: T, a: [*]T, lda: isize, x: [*]T, incx: isize, beta: T, y: [*]T, incy: isize) Error!void {
    comptime util.requireComplex(T, "hemv");
    return mv(allocator, io, T, true, uplo, n, alpha, a, lda, x, incx, beta, y, incy);
}

fn mv(allocator: Allocator, io: Io, comptime T: type, comptime herm: bool, uplo: Uplo, n: isize, alpha: T, a: [*]T, lda: isize, x: [*]T, incx: isize, beta: T, y: [*]T, incy: isize) Error!void {
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, n)) return error.InvalidLda;
    if (incx == 0) return error.InvalidIncx;
    if (incy == 0) return error.InvalidIncy;
    if (n == 0) return;
    if (!util.fitsMatrix(n, n, lda, @sizeOf(T))) return error.InvalidLda;
    if (!util.fits(n, incx, @sizeOf(T))) return error.InvalidIncx;
    if (!util.fits(n, incy, @sizeOf(T))) return error.InvalidIncy;
    if (accel.symv(T, herm, uplo, n, alpha, a, lda, x, incx, beta, y, incy)) return;
    if (!isOne(T, beta)) Ops(T).init().scal(n, beta, y, @intCast(@abs(incy)));
    if (isZero(T, alpha)) return;
    const px = first(T, x, n, incx);
    const py = first(T, y, n, incy);
    switch (uplo) {
        inline else => |u| try mvRun(allocator, io, T, herm, u, n, alpha, a, lda, px, incx, py, incy),
    }
}

fn mvRun(allocator: Allocator, io: Io, comptime T: type, comptime herm: bool, comptime uplo: Uplo, n: isize, alpha: T, a: [*]T, lda: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) (Allocator.Error || Io.Cancelable)!void {
    const F = Float(T);
    const complex = comptime isComplex(T);
    // symv has no `V` and `M` kernels; `symv` returned above for those
    const rev = comptime uplo == .v or uplo == .m;
    if (comptime !herm and rev) unreachable;
    const op: Op = comptime switch (uplo) {
        .u => if (herm) .hemv_U else .symv_U,
        .l => if (herm) .hemv_L else .symv_L,
        .v => .hemv_V,
        .m => .hemv_M,
    };
    const lower = comptime uplo == .l or uplo == .m;
    // the portable kernel on a core without its own gemv runs plain loops,
    // which need no buffer (`symmetric.Auto`), as does any type without slots
    const direct = if (comptime portable.slot(T, op)) |s| blk: {
        const kind: sym.Kind = comptime if (!herm) .symmetric else if (rev) .hermitian_rev else .hermitian;
        const auto = sym.Auto(portable.letter(T).?, T, lower, kind);
        break :blk @field(blas.table().*, s) == (if (complex) &auto.cplx else &auto.real) and !auto.native();
    } else true;
    const nu: usize = @intCast(n);
    // `symv_thread.c` splits the columns on multiples of 4, balanced by area
    const p = threading.plan(if (herm) .hemv else .symv, T, nu, nu, if (lower) .falling else .rising, 4);
    if (p.chunks.count > 1) return mvChunks(allocator, io, T, op, lower, direct, &p, n, alpha, a, lda, x, incx, y, incy);
    var none: [1]F = undefined;
    const sc: Scratch(F) = if (direct) .{ .ptr = &none, .len = 0 } else blk: {
        const bytes = mvBytes(F, complex, nu, incx == 1, incy == 1);
        break :blk try Scratch(F).heapOnly(allocator, (bytes + @sizeOf(F) - 1) / @sizeOf(F));
    };
    defer sc.deinit(allocator);
    // x86_64's SSE `symv` and `hemv` read past their copies of `x` and `y`,
    // and NaN there reaches `y`
    if (sysinfo.target.cpu.arch == .x86_64) @memset(sc.slice(), 0);
    Ops(T).init().mv(op, n, n, alpha, a, lda, x, incx, y, incy, sc.ptr);
}

/// Bytes of buffer the core's `symv` or `hemv` kernel of `m` needs, whichever
/// it is (`buffers.md`).
fn mvBytes(comptime F: type, comptime complex: bool, m: usize, unit_x: bool, unit_y: bool) usize {
    const cs = if (complex) 2 else 1;
    return @max(sym.composeBytes(F, cs, m, unit_x, unit_y), sym.sseBytes(F, cs, m));
}

/// `symv_thread.c`'s split over `p`'s column chunks: each chunk's columns
/// through the kernel into a partial vector of its own (rows `0..end` of the
/// upper triangle, `start..n` of the lower), with a page-aligned buffer of its
/// own, then the partials added to `y` in chunk order (`threading.Partials`).
fn mvChunks(allocator: Allocator, io: Io, comptime T: type, comptime op: Op, comptime lower: bool, direct: bool, p: *const threading.Plan, n: isize, alpha: T, a: [*]T, lda: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) (Allocator.Error || Io.Cancelable)!void {
    const F = Float(T);
    const nu: usize = @intCast(n);
    const G = struct {
        n: usize,
        alpha: T,
        a: [*]T,
        lda: isize,
        x: [*]T,
        incx: isize,
        direct: bool,
        parts: *const threading.Partials(T),

        fn rows(g: *const @This(), c: threading.Chunk) threading.Rows {
            return if (lower) .{ .lo = c.start, .hi = g.n } else .{ .lo = 0, .hi = c.end };
        }

        /// Kernel rows of the chunk: its partial's length.
        fn m(g: *const @This(), c: threading.Chunk) usize {
            return if (lower) g.n - c.start else c.end;
        }

        /// Bytes of the chunk's buffer.
        fn len(g: *const @This(), c: threading.Chunk) usize {
            const bytes = if (g.direct) 0 else mvBytes(F, comptime isComplex(T), g.m(c), g.incx == 1, true);
            return std.mem.alignForward(usize, bytes, heap_align);
        }

        fn one(g: *const @This(), c: threading.Chunk, b: []u8) void {
            const part = g.parts.take(c.index);
            if (sysinfo.target.cpu.arch == .x86_64) @memset(b, 0);
            const c0: isize = @intCast(c.start);
            const off: isize = if (lower) c0 else 0;
            Ops(T).init().mv(op, @intCast(g.m(c)), @intCast(c.end - c.start), g.alpha, @ptrCast(util.at(T, g.a, off * (g.lda + 1))), g.lda, @ptrCast(util.at(T, g.x, off * g.incx)), g.incx, part.ptr, 1, @ptrCast(@alignCast(b.ptr)));
        }
    };
    var g: G = .{ .n = nu, .alpha = alpha, .a = a, .lda = lda, .x = x, .incx = incx, .direct = direct, .parts = undefined };
    const parts: threading.Partials(T) = try .init(allocator, p, @as(*const G, &g), G.rows);
    defer parts.deinit(allocator);
    g.parts = &parts;
    try threading.forEachPart(allocator, io, p, u8, heap_align, @as(*const G, &g), G.len, G.one);
    try parts.fold(io, y, incy);
}
