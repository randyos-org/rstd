//! `?TRMV`, `?TRSV`, `?TBMV`, `?TBSV`, `?TPMV`, `?TPSV`: the checks of
//! `interface/{trmv,trsv,tbmv,tbsv,tpmv,tpsv}.c` over
//! `driver/level2/{,z}{trmv,trsv,tbmv,tbsv,tpmv,tpsv}_{U,L}.c`, loops over the
//! copy, axpy, dot and gemv kernels. `_U` also does the transposed lower
//! matrix, as upstream; the products thread (`Split`), the solves are serial.

const std = @import("std");
const sysinfo = @import("builtin");
const blas = @import("../__root__.zig");
const types = @import("../types.zig");
const util = @import("../util.zig");
const Vtable = @import("../Vtable.zig");
const portable = @import("../portable.zig");
const accel = @import("../apple_accel.zig");
const threading = @import("../threading.zig");
const interface = @import("interface.zig");
const drivers = @import("triangular_driver.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const BlasLong = types.BlasLong;
const Trans = types.Trans;
const Uplo = types.Uplo;
const Diag = types.Diag;

pub const Error = interface.Error;

/// Work space kept on the stack; bigger problems use the allocator.
pub const stack_bytes = interface.stack_bytes;

/// `dtb_entries` of a hand-built vtable that has none.
const default_dtb = 64;

// ---- element helpers -------------------------------------------------------

const isComplex = types.isComplex;
const Float = types.Float;

const raw = util.raw;

const at = util.at;

const off = util.off;

const first = util.first;

const Num = @import("ops.zig").Num;
const Ops = @import("ops.zig").Ops;

/// Real types have no `r` and `c`: they are `n` and `t`.
pub fn canon(comptime T: type, op: Trans) Trans {
    if (isComplex(T)) return op;
    return if (op.transposed()) .t else .n;
}

// ---- drivers ---------------------------------------------------------------

const File = drivers.File;
const Drive = drivers.Drive;

// ---- interface -------------------------------------------------------------

pub const Kind = enum { trmv, trsv, tbmv, tbsv, tpmv, tpsv };

/// Elements of `F` a copy of `x` takes, rounded to a 64-byte multiple so
/// the `gemv` buffer after it stays aligned (`f80` on x86 is 12 bytes).
fn copyElems(comptime T: type, n: usize) usize {
    const per = 64 / std.math.gcd(64, @sizeOf(Float(T)));
    const comp: usize = if (isComplex(T)) 2 else 1;
    return std.mem.alignForward(usize, n * comp, per);
}

fn run(
    comptime T: type,
    comptime kind: Kind,
    comptime file: File,
    comptime op: Trans,
    comptime unit: bool,
    o: Ops(T),
    n: BlasLong,
    k: BlasLong,
    a: [*]T,
    lda: BlasLong,
    x: [*]T,
    incx: BlasLong,
    buf: [*]Float(T),
) void {
    const D = Drive(T, file, op, unit);
    var b = x;
    var g = buf;
    if (incx != 1) {
        b = @ptrCast(@alignCast(buf));
        g = buf + copyElems(T, @intCast(n));
        o.copy(n, x, incx, b, 1);
    }
    switch (kind) {
        .trmv => D.trmv(o, n, a, lda, b, g),
        .trsv => D.trsv(o, n, a, lda, b, g),
        .tbmv => D.tbmv(o, n, k, a, lda, b),
        .tbsv => D.tbsv(o, n, k, a, lda, b),
        .tpmv => D.tpmv(o, n, a, b),
        .tpsv => D.tpsv(o, n, a, b),
    }
    if (incx != 1) o.copy(n, b, 1, x, incx);
}

/// The public routines: Accelerate where it runs them (`apple_accel.zig`),
/// after `call`'s checks, else `call` on the selected kernels.
fn public(allocator: Allocator, io: Io, comptime T: type, comptime kind: Kind, uplo: Uplo, trans: Trans, diag: Diag, n: isize, k: isize, a: [*]T, lda: isize, x: [*]T, incx: isize) Error!void {
    if (comptime accel.covers(T)) {
        if (!try checks(T, kind, uplo, trans, n, k, lda, incx)) return;
        if (accel.tri(T, @field(accel.TriKind, @tagName(kind)), uplo, trans, diag, n, k, a, lda, x, incx)) return;
    }
    return call(allocator, io, T, kind, blas.vtable(), uplo, trans, diag, n, k, a, lda, x, incx);
}

/// The reference's checks, and nothing past the address space; false for
/// nothing to do.
fn checks(comptime T: type, comptime kind: Kind, uplo: Uplo, trans: Trans, n: isize, k: isize, lda: isize, incx: isize) interface.ArgError!bool {
    const banded = kind == .tbmv or kind == .tbsv;
    const packed_ = kind == .tpmv or kind == .tpsv;
    if (uplo.isRev()) return error.InvalidUplo;
    if (trans.conjX()) return error.InvalidTrans;
    if (n < 0) return error.InvalidN;
    if (banded) {
        if (k < 0) return error.InvalidK;
        if (lda < k + 1) return error.InvalidLda;
    } else if (!packed_ and lda < @max(1, n)) return error.InvalidLda;
    if (incx == 0) return error.InvalidIncx;
    if (n == 0) return false;
    // nothing past the address space (`util.fits`)
    if (packed_) {
        if (!util.fitsMatrix(n, n, n, @sizeOf(T))) return error.InvalidN;
    } else if (!util.fitsMatrix(if (banded) lda else n, n, lda, @sizeOf(T))) return error.InvalidLda;
    if (!util.fits(n, incx, @sizeOf(T))) return error.InvalidIncx;
    return true;
}

/// The checks, work space and driver choice common to the six routines;
/// `vt` is the kernels' vtable.
pub fn call(
    allocator: Allocator,
    io: Io,
    comptime T: type,
    comptime kind: Kind,
    vt: *const Vtable,
    uplo: Uplo,
    trans: Trans,
    diag: Diag,
    n: isize,
    k: isize,
    a: [*]T,
    lda: isize,
    x: [*]T,
    incx: isize,
) Error!void {
    const F = Float(T);
    const banded = kind == .tbmv or kind == .tbsv;
    const packed_ = kind == .tpmv or kind == .tpsv;
    if (!try checks(T, kind, uplo, trans, n, k, lda, incx)) return;

    const kn: BlasLong = n;
    const o: Ops(T) = .{ .vt = vt, .dtb = if (vt.dtb_entries > 0) vt.dtb_entries else default_dtb };
    const px = first(T, x, kn, incx);
    if (comptime kind == .trmv or kind == .tbmv or kind == .tpmv) {
        const nu: usize = @intCast(n);
        // `{trmv,tbmv,tpmv}_thread.c` split the columns on multiples of 8,
        // balanced by length but for the band
        const shape: threading.Shape = if (kind == .tbmv) .even else if (uplo == .l) .falling else .rising;
        const p = threading.plan(if (kind == .trmv) .trmv else .tpmv, T, nu, nu, shape, 8);
        if (p.chunks.count > 1) return split(allocator, io, T, kind, o, &p, uplo, trans, diag, kn, k, a, lda, px, incx);
    }

    // x's copy, then `gemv`'s buffer (the full routines, past one block); the
    // 40-element pad is ours (`interface/ztrmv.c:245` pads its small case by 40)
    const un: usize = @intCast(n);
    var elems: usize = if (incx != 1) copyElems(T, un) else 0;
    if (!banded and !packed_ and kn > o.dtb) {
        const comp: usize = if (isComplex(T)) 2 else 1;
        elems += comp * (2 * un + @as(usize, @intCast(o.dtb))) + 40;
    }
    // the copy of `x` is rounded to 64 bytes, so the gemv buffer after it
    // keeps the stack's alignment
    var stack: [stack_bytes]u8 align(64) = undefined;
    const sc = try util.Scratch(F, 64, 64).init(allocator, &stack, elems);
    defer sc.deinit(allocator);
    const buf = sc.slice();
    // `zgemv_t.S` and its siblings return NaN if the buffer holds NaN or
    // infinity; cleared as `interface.gemv` does
    if (comptime isComplex(T) and sysinfo.target.cpu.arch == .x86_64) @memset(buf, 0);

    // `_U` for upper untransposed or lower transposed
    const fl: File = if ((uplo == .u) != trans.transposed()) .u else .l;
    switch (fl) {
        inline else => |f| switch (trans) {
            inline .n, .t, .r, .c => |t| switch (diag) {
                inline else => |d| run(T, kind, f, canon(T, t), d == .u, o, kn, k, a, lda, px, incx, buf.ptr),
            },
            .o, .u, .s, .d => unreachable,
        },
    }
}

/// The products over more than one chunk of columns: each chunk's part of
/// `op(A) * x` into a partial vector of its own, the partials then folded into
/// `x` in chunk order (`threading.Partials`). With a transpose a chunk's part
/// is its own rows, so nothing overlaps.
fn split(allocator: Allocator, io: Io, comptime T: type, comptime kind: Kind, o: Ops(T), p: *const threading.Plan, uplo: Uplo, trans: Trans, diag: Diag, n: BlasLong, k: BlasLong, a: [*]T, lda: BlasLong, x: [*]T, incx: BlasLong) (Allocator.Error || Io.Cancelable)!void {
    const F = Float(T);
    const un: usize = @intCast(n);
    const comp: usize = if (isComplex(T)) 2 else 1;
    // a strided `x`'s copy, then each chunk's `gemv` buffer, as `run` lays
    // them out
    const xe: usize = if (incx != 1) copyElems(T, un) else 0;
    const ge: usize = if (kind == .trmv) copyElems(F, comp * (2 * un + @as(usize, @intCast(o.dtb))) + 40) else 0;
    const buf = try allocator.alignedAlloc(F, .@"64", xe + p.chunks.count * ge);
    defer allocator.free(buf);
    // as `call` clears it
    if (comptime isComplex(T) and sysinfo.target.cpu.arch == .x86_64) @memset(buf, 0);
    const xs: [*]T = if (incx != 1) @ptrCast(buf.ptr) else x;
    if (incx != 1) o.copy(n, x, incx, xs, 1);
    switch (uplo == .l) {
        inline else => |lower| switch (trans) {
            inline .n, .t, .r, .c => |t| switch (diag) {
                inline else => |d| {
                    const S = Split(T, kind, lower, canon(T, t), d == .u);
                    var s: S = .{ .o = o, .n = n, .k = k, .a = a, .lda = lda, .x = xs, .g = buf.ptr + xe, .ge = ge, .parts = undefined };
                    const parts: threading.Partials(T) = try .init(allocator, p, @as(*const S, &s), S.rows);
                    defer parts.deinit(allocator);
                    s.parts = &parts;
                    try threading.forEach(io, p, @as(*const S, &s), S.one);
                    // `x` is the output: zeroed, then the partials added
                    for (0..un) |i| at(T, x, @as(BlasLong, @intCast(i)) * incx).* = std.mem.zeroes(T);
                    try parts.fold(io, x, incx);
                },
            },
            .o, .u, .s, .d => unreachable,
        },
    }
}

/// A chunk of columns of `op(A) * x`, out of place: the `trmv_kernel`s of
/// `driver/level2/trmv_thread.c`, `tbmv_thread.c` and `tpmv_thread.c`.
fn Split(comptime T: type, comptime kind: Kind, comptime lower: bool, comptime op: Trans, comptime unit: bool) type {
    const F = Float(T);
    const N = Num(T);
    const tr = comptime op.transposed();
    const cj = comptime op.conjA();
    return struct {
        const Self = @This();

        o: Ops(T),
        n: BlasLong,
        k: BlasLong,
        a: [*]T,
        lda: BlasLong,
        /// Unit stride.
        x: [*]T,
        /// The `gemv` buffers, `ge` elements each.
        g: [*]F,
        ge: usize,
        parts: *const threading.Partials(T),

        /// The rows of `op(A) * x` columns `c.start..c.end` reach.
        fn rows(s: *const Self, c: threading.Chunk) threading.Rows {
            if (tr) return .{ .lo = c.start, .hi = c.end };
            const nu: usize = @intCast(s.n);
            const k: usize = if (kind == .tbmv) @intCast(s.k) else nu;
            return if (lower) .{ .lo = c.start, .hi = @min(nu, c.end + k) } else .{ .lo = c.start -| k, .hi = c.end };
        }

        fn one(s: *const Self, c: threading.Chunk) void {
            const y = s.parts.take(c.index).ptr;
            const yb: BlasLong = @intCast(s.parts.rows[c.index].lo);
            const c0: BlasLong = @intCast(c.start);
            const c1: BlasLong = @intCast(c.end);
            if (kind == .trmv) s.full(y, yb, c0, c1, s.g + c.index * s.ge) else s.columns(y, yb, c0, c1);
        }

        /// `trmv_thread.c`: blocks of `dtb` columns from `c0`, each a `gemv`
        /// off the block and `axpy` or `dot` calls inside it.
        fn full(s: *const Self, y: [*]T, yb: BlasLong, c0: BlasLong, c1: BlasLong, g: [*]F) void {
            const o = s.o;
            const a = s.a;
            const lda = s.lda;
            const x = s.x;
            const m = s.n;
            var is = c0;
            while (is < c1) : (is += o.dtb) {
                const min_i = @min(c1 - is, o.dtb);
                if (!lower and is > 0) {
                    const pa = off(T, a, is * lda);
                    if (tr) o.gemv(op, is, min_i, 1, pa, lda, x, off(T, y, is - yb), g) else o.gemv(op, is, min_i, 1, pa, lda, off(T, x, is), off(T, y, -yb), g);
                }
                var i = is;
                while (i < is + min_i) : (i += 1) {
                    const col = off(T, a, i * lda);
                    const yi = at(T, y, i - yb);
                    const xi = at(T, x, i).*;
                    if (!lower and i > is) {
                        if (tr) yi.* = N.add(yi.*, o.dot(cj, i - is, off(T, col, is), off(T, x, is))) else o.axpy(cj, i - is, xi, off(T, col, is), off(T, y, is - yb));
                    }
                    yi.* = N.add(yi.*, if (unit) xi else N.mul(at(T, col, i).*, xi, cj));
                    if (lower and is + min_i > i + 1) {
                        const len = is + min_i - i - 1;
                        if (tr) yi.* = N.add(yi.*, o.dot(cj, len, off(T, col, i + 1), off(T, x, i + 1))) else o.axpy(cj, len, xi, off(T, col, i + 1), off(T, y, i + 1 - yb));
                    }
                }
                if (lower and m > is + min_i) {
                    const pa = off(T, a, is + min_i + is * lda);
                    if (tr) o.gemv(op, m - is - min_i, min_i, 1, pa, lda, off(T, x, is + min_i), off(T, y, is - yb), g) else o.gemv(op, m - is - min_i, min_i, 1, pa, lda, off(T, x, is), off(T, y, is + min_i - yb), g);
                }
            }
        }

        /// `tbmv_thread.c`, `tpmv_thread.c`: one column at a time, the
        /// off-diagonal part before the diagonal above it, after it below.
        fn columns(s: *const Self, y: [*]T, yb: BlasLong, c0: BlasLong, c1: BlasLong) void {
            const n = s.n;
            const x = s.x;
            var i = c0;
            while (i < c1) : (i += 1) {
                const room = if (lower) n - 1 - i else i;
                const len = if (kind == .tbmv) @min(room, s.k) else room;
                const col: BlasLong = if (kind == .tbmv) i * s.lda else if (lower) i * n - @divExact(i * (i - 1), 2) else @divExact(i * (i + 1), 2);
                const dp = off(T, s.a, col + (if (lower) 0 else if (kind == .tbmv) s.k else i));
                const yi = at(T, y, i - yb);
                const xi = at(T, x, i).*;
                if (!lower and len > 0) s.part(y, yb, yi, xi, len, off(T, dp, -len), i - len);
                yi.* = N.add(yi.*, if (unit) xi else N.mul(dp[0], xi, cj));
                if (lower and len > 0) s.part(y, yb, yi, xi, len, off(T, dp, 1), i + 1);
            }
        }

        /// Column `i`'s `len` off-diagonal elements at `ap`, rows `r0..`.
        fn part(s: *const Self, y: [*]T, yb: BlasLong, yi: *T, xi: T, len: BlasLong, ap: [*]T, r0: BlasLong) void {
            if (tr) yi.* = N.add(yi.*, s.o.dot(cj, len, ap, off(T, s.x, r0))) else s.o.axpy(cj, len, xi, ap, off(T, y, r0 - yb));
        }
    };
}

/// `?TRMV` (`interface/trmv.c`, `ztrmv.c`): `x := op(A) * x`, `A` upper or
/// lower triangular of `n`, column-major with leading dimension `lda`. The
/// other triangle (and the diagonal for `unit`) is not read. A negative
/// `incx` walks `x` from its end.
pub fn trmv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, diag: Diag, n: isize, a: [*]T, lda: isize, x: [*]T, incx: isize) Error!void {
    return public(allocator, io, T, .trmv, uplo, trans, diag, n, 0, a, lda, x, incx);
}

/// `?TRSV` (`interface/trsv.c`, `ztrsv.c`): solves `op(A) * y = x`, `y`
/// replacing `x`. No singularity check, as in the reference.
pub fn trsv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, diag: Diag, n: isize, a: [*]T, lda: isize, x: [*]T, incx: isize) Error!void {
    return public(allocator, io, T, .trsv, uplo, trans, diag, n, 0, a, lda, x, incx);
}

/// `trsv` over `vt`'s kernels at unit stride, never routed to Accelerate:
/// LAPACK's `getrs` with one right-hand side (`lapack/getrs/getrs_single.c`),
/// which calls the level-3 drivers alike. The caller checks the arguments.
pub fn trsvTable(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, uplo: Uplo, trans: Trans, diag: Diag, n: usize, a: [*]const T, lda: usize, x: [*]T) (Allocator.Error || Io.Cancelable)!void {
    call(allocator, io, T, .trsv, vt, uplo, trans, diag, @intCast(n), 0, @constCast(a), @intCast(lda), x, 1) catch |e| switch (e) {
        error.OutOfMemory, error.Canceled => |f| return f,
        else => unreachable,
    };
}

/// `?TBMV` (`interface/tbmv.c`, `ztbmv.c`): `x := op(A) * x`, `A` triangular
/// band with `k` off-diagonals in LAPACK band storage: `A(i,j)` is
/// `a[k + i - j + j*lda]` (upper) or `a[i - j + j*lda]` (lower).
pub fn tbmv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, diag: Diag, n: isize, k: isize, a: [*]T, lda: isize, x: [*]T, incx: isize) Error!void {
    return public(allocator, io, T, .tbmv, uplo, trans, diag, n, k, a, lda, x, incx);
}

/// `?TBSV` (`interface/tbsv.c`, `ztbsv.c`): solves `op(A) * y = x` for band `A`.
pub fn tbsv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, diag: Diag, n: isize, k: isize, a: [*]T, lda: isize, x: [*]T, incx: isize) Error!void {
    return public(allocator, io, T, .tbsv, uplo, trans, diag, n, k, a, lda, x, incx);
}

/// `?TPMV` (`interface/tpmv.c`, `ztpmv.c`): `x := op(A) * x`, `A` packed by
/// columns: `A(i,j)` is `ap[j*(j+1)/2 + i]` (upper, `i <= j`) or
/// `ap[j*(2n-j+1)/2 + i - j]` (lower, `i >= j`).
pub fn tpmv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, diag: Diag, n: isize, ap: [*]T, x: [*]T, incx: isize) Error!void {
    return public(allocator, io, T, .tpmv, uplo, trans, diag, n, 0, ap, 1, x, incx);
}

/// `?TPSV` (`interface/tpsv.c`, `ztpsv.c`): solves `op(A) * y = x` for packed `A`.
pub fn tpsv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, diag: Diag, n: isize, ap: [*]T, x: [*]T, incx: isize) Error!void {
    return public(allocator, io, T, .tpsv, uplo, trans, diag, n, 0, ap, 1, x, incx);
}
