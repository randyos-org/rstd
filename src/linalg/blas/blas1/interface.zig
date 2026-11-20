//! The level-1 interface: `interface/*.c`'s checks and pointer set-up, for
//! any float `T` or its `Complex` (`portable.kernel`). A negative stride
//! starts at the last element; kernels get `n >= 1`; a vector past the
//! address space is its stride's `ArgError`. Above OpenBLAS's thresholds
//! most routines thread on `io`, reductions folded in chunk order.

const std = @import("std");
const blas = @import("../__root__.zig");
const threading = @import("../threading.zig");
const types = @import("../types.zig");
const util = @import("../util.zig");
const portable = @import("../portable.zig");
const bf = @import("../ext/bf16/__root__.zig");
const accel = @import("../apple_accel.zig");
const rot_root = @import("rot/__root__.zig");
const rot_reals = @import("rot_reals.zig");
const extremes = @import("interface_extreme.zig");
const half = @import("interface_bf16.zig");

const Io = std.Io;
const BlasLong = types.BlasLong;
const BF16 = types.BF16;
const XDouble = types.XDouble;

const isComplex = types.isComplex;
const Float = types.Float;

const raw = util.raw;

/// Bad arguments, by the reference's names.
pub const ArgError = error{ InvalidIncx, InvalidIncy };
/// What a level-1 routine returns.
pub const Error = ArgError || Io.Cancelable;

const first = util.first;
const zero = util.zero;

/// Refuses an `x` that could not exist in memory (`util.fits`), so no index
/// overflows.
pub fn spanX(comptime T: type, n: isize, incx: isize) ArgError!void {
    if (!util.fits(n, incx, @sizeOf(T))) return error.InvalidIncx;
}

/// `spanX` for `x` of `X`, then `y` of `Y`.
pub fn spanXY(comptime X: type, comptime Y: type, n: isize, incx: isize, incy: isize) ArgError!void {
    try spanX(X, n, incx);
    if (!util.fits(n, incy, @sizeOf(Y))) return error.InvalidIncy;
}

// ---- threading ---------------------------------------------------------------

/// Chunks start on multiples of this many elements: every chunk but the last
/// is whole kernel blocks, the last ends where the vector does, and each
/// starts at the vector's alignment, so every element takes the path, and
/// gets the bits, one call would give it.
pub const quantum = 1024;

/// `family`'s plan for `n` elements, if it has more than one chunk.
pub fn chunked(comptime family: threading.Family, comptime T: type, n: isize) ?threading.Plan {
    const len: usize = @intCast(n);
    if (threading.serial(family, T, len)) return null;
    const p = threading.plan(family, T, len, len, .even, quantum);
    return if (p.chunks.count > 1) p else null;
}

/// `family`'s plan for `n` elements, if it takes more than one worker; an
/// element-wise routine otherwise makes one kernel call, with the same bits.
fn spread(comptime family: threading.Family, comptime T: type, n: isize) ?threading.Plan {
    const p = chunked(family, T, n) orelse return null;
    return if (p.workers > 1) p else null;
}

/// `call(ctx, len, at)` for each chunk of `p`, `at` its first element.
fn each(io: Io, p: *const threading.Plan, ctx: anytype, comptime call: fn (@TypeOf(ctx), isize, isize) void) Io.Cancelable!void {
    const S = struct {
        fn one(c: @TypeOf(ctx), ch: threading.Chunk) void {
            call(c, @intCast(ch.end - ch.start), @intCast(ch.start));
        }
    };
    return threading.forEach(io, p, ctx, S.one);
}

/// `part(ctx, len, at)` for each chunk of `p`, folded in chunk order.
pub fn fold(io: Io, p: *const threading.Plan, comptime R: type, ctx: anytype, comptime part: fn (@TypeOf(ctx), isize, isize) R, comptime combine: fn (R, R) R) Io.Cancelable!R {
    const S = struct {
        fn one(c: @TypeOf(ctx), ch: threading.Chunk) R {
            return part(c, @intCast(ch.end - ch.start), @intCast(ch.start));
        }
    };
    return threading.reduce(io, p, R, ctx, S.one, combine);
}

/// Element `at` of a vector whose element 0 is at `p` (as `first` gives
/// it), at stride `inc`.
pub fn shift(comptime T: type, p: [*]T, at: isize, inc: isize) [*]T {
    return @ptrCast(util.at(T, p, at * inc));
}

/// Two vectors from their element 0, and the routine's scalars.
pub fn Pair(comptime T: type, comptime Extra: type) type {
    return struct { x: [*]T, incx: isize, y: [*]T, incy: isize, e: Extra };
}

/// Whether the sums thread for `T`: `f32` and `f64` and their `Complex`.
/// Wider, they give the bits of OpenBLAS's x87 kernels (`readme.md`), which
/// no fold keeps; the other widths run portable kernels only.
fn sums(comptime T: type) bool {
    const F = Float(T);
    return (F == f32 or F == f64) and !types.isExtended(F);
}

/// What the partial sums of `T` fold in: `f64` for `f32`, so the result
/// rounds to `f32` once.
fn Wide(comptime T: type) type {
    const W = if (Float(T) == f32) f64 else Float(T);
    return if (isComplex(T)) std.math.Complex(W) else W;
}

/// The fold of the partial sums; on the x87 each sum is rounded to its
/// float (`util.rounded`), so the bits do not hang on what stays in a register.
fn plus(comptime W: type) fn (W, W) W {
    return struct {
        fn f(a: W, b: W) W {
            if (comptime !isComplex(W)) return util.rounded(W, a + b);
            const F = @TypeOf(a.re);
            return .init(util.rounded(F, a.re + b.re), util.rounded(F, a.im + b.im));
        }
    }.f;
}

/// `W` (`Wide(T)`) rounded to `T`.
fn rounded(comptime T: type, w: Wide(T)) T {
    if (comptime isComplex(T)) return .init(@floatCast(w.re), @floatCast(w.im));
    return @floatCast(w);
}

/// A chunk's `dsdot`: `f32` products summed in `f64`.
fn dsdotPart(q: Pair(f32, void), len: isize, at: isize) f64 {
    return blas.kernel(.dsdot_k)(len, shift(f32, q.x, at, q.incx), q.incx, shift(f32, q.y, at, q.incy), q.incy);
}

// ---- copy, swap, scal --------------------------------------------------------

/// `?COPY` (`interface/copy.c`): `y := x`.
pub fn copy(io: Io, comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) Error!void {
    if (n <= 0) return;
    try spanXY(T, T, n, incx, incy);
    if (try accel.copy(io, T, n, x, incx, y, incy)) return;
    const c: Pair(T, void) = .{ .x = first(T, x, n, incx), .incx = incx, .y = first(T, y, n, incy), .incy = incy, .e = {} };
    const S = struct {
        fn run(q: @TypeOf(c), len: isize, at: isize) void {
            _ = portable.fetch(T, .copy_k, blas.vtable)(len, raw(T, shift(T, q.x, at, q.incx)), q.incx, raw(T, shift(T, q.y, at, q.incy)), q.incy);
        }
    };
    // a zero stride runs serially, as OpenBLAS's do: the elements depend on each other
    if (incx != 0 and incy != 0) if (spread(.copy, T, n)) |p| return each(io, &p, c, S.run);
    S.run(c, n, 0);
}

/// `?SWAP` (`interface/swap.c`, `zswap.c`): `x <-> y`.
pub fn swap(io: Io, comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) Error!void {
    if (n <= 0) return;
    try spanXY(T, T, n, incx, incy);
    if (try accel.swap(io, T, n, x, incx, y, incy)) return;
    const c: Pair(T, void) = .{ .x = first(T, x, n, incx), .incx = incx, .y = first(T, y, n, incy), .incy = incy, .e = {} };
    const S = struct {
        fn run(q: @TypeOf(c), len: isize, at: isize) void {
            const px = raw(T, shift(T, q.x, at, q.incx));
            const py = raw(T, shift(T, q.y, at, q.incy));
            if (comptime isComplex(T)) {
                _ = portable.fetch(T, .swap_k, blas.vtable)(len, 0, 0, 0, 0, px, q.incx, py, q.incy, py, 0);
            } else {
                _ = portable.fetch(T, .swap_k, blas.vtable)(len, 0, 0, 0, px, q.incx, py, q.incy, py, 0);
            }
        }
    };
    if (incx != 0 and incy != 0) if (spread(.swap, T, n)) |p| return each(io, &p, c, S.run);
    S.run(c, n, 0);
}

/// `?SCAL` (`interface/scal.c`, `zscal.c`): `x := alpha*x`. Nothing for
/// `incx <= 0`, and nothing to do for `alpha == 1`. With the flag the
/// interface passes (1), a zero `alpha` leaves NaN and infinity as NaN, as
/// `alpha*x` would, and finite elements 0.
pub fn scal(io: Io, comptime T: type, n: isize, alpha: T, x: [*]T, incx: isize) Error!void {
    if (n <= 0 or incx <= 0) return;
    try spanX(T, n, incx);
    if (comptime isComplex(T)) {
        if (alpha.re == 1 and alpha.im == 0) return;
    } else if (alpha == 1) return;
    if (try accel.scal(io, T, n, alpha, x, incx)) return;
    const c: Pair(T, T) = .{ .x = x, .incx = incx, .y = x, .incy = incx, .e = alpha };
    const S = struct {
        fn run(q: @TypeOf(c), len: isize, at: isize) void {
            const px = raw(T, shift(T, q.x, at, q.incx));
            _ = @call(.auto, portable.fetch(T, .scal_k, blas.vtable), .{ len, 0, 0 } ++ util.scalarArgs(T, q.e) ++ .{ px, q.incx, px, 0, px, 1 });
        }
    };
    if (spread(.scal, T, n)) |p| return each(io, &p, c, S.run);
    S.run(c, n, 0);
}

/// `CSSCAL`/`ZDSCAL`: complex `x` scaled by a real `alpha`.
pub fn scalReal(io: Io, comptime T: type, n: isize, alpha: Float(T), x: [*]T, incx: isize) Error!void {
    comptime util.requireComplex(T, "scalReal");
    return scal(io, T, n, .init(alpha, 0), x, incx);
}

// ---- dot ---------------------------------------------------------------------

/// `?DOT`/`?DOTU` (`interface/dot.c`, `zdot.c`): `x . y`, 0 for `n <= 0`.
pub fn dot(io: Io, comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) Error!T {
    return dotImpl(io, T, false, n, x, incx, y, incy);
}

/// `?DOTC`: `conj(x) . y`. Complex only.
pub fn dotc(io: Io, comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) Error!T {
    comptime util.requireComplex(T, "dotc");
    return dotImpl(io, T, true, n, x, incx, y, incy);
}

fn dotImpl(io: Io, comptime T: type, comptime conj: bool, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) Error!T {
    if (n <= 0) return zero(T);
    try spanXY(T, T, n, incx, incy);
    if (try accel.dot(io, T, conj, n, x, incx, y, incy)) |r| return r;
    // `kernel/x86_64/ddot.c`, `zdot.c`: a zero stride runs serially
    if (comptime sums(T)) if (incx != 0 and incy != 0) if (chunked(.dot, T, n)) |p| {
        const c: Pair(T, void) = .{ .x = first(T, x, n, incx), .incx = incx, .y = first(T, y, n, incy), .incy = incy, .e = {} };
        const S = struct {
            fn part(q: @TypeOf(c), len: isize, at: isize) Wide(T) {
                // `sdot` sums in `f64` and rounds once (`readme.md`)
                if (comptime T == f32) return dsdotPart(q, len, at);
                const a = raw(T, shift(T, q.x, at, q.incx));
                const b = raw(T, shift(T, q.y, at, q.incy));
                if (comptime isComplex(T)) {
                    const r = types.complexOf(portable.fetch(T, if (conj) .dotc_k else .dotu_k, blas.vtable)(len, a, q.incx, b, q.incy));
                    return .init(r.re, r.im);
                }
                return portable.fetch(T, .dot_k, blas.vtable)(len, a, q.incx, b, q.incy);
            }
        };
        return rounded(T, try fold(io, &p, Wide(T), c, S.part, plus(Wide(T))));
    };
    const px = raw(T, first(T, x, n, incx));
    const py = raw(T, first(T, y, n, incy));
    if (comptime isComplex(T)) {
        return types.complexOf(portable.fetch(T, if (conj) .dotc_k else .dotu_k, blas.vtable)(n, px, incx, py, incy));
    }
    return portable.fetch(T, .dot_k, blas.vtable)(n, px, incx, py, incy);
}

/// `DSDOT` (`interface/dsdot.c`): the dot of two `f32` vectors, summed and
/// returned in `f64`.
pub fn dsdot(io: Io, n: isize, x: [*]f32, incx: isize, y: [*]f32, incy: isize) Error!f64 {
    if (n <= 0) return 0;
    try spanXY(f32, f32, n, incx, incy);
    if (try accel.dsdot(io, n, x, incx, y, incy)) |r| return r;
    const c: Pair(f32, void) = .{ .x = first(f32, x, n, incx), .incx = incx, .y = first(f32, y, n, incy), .incy = incy, .e = {} };
    if (incx != 0 and incy != 0) if (chunked(.dot, f32, n)) |p| return fold(io, &p, f64, c, dsdotPart, plus(f64));
    return dsdotPart(c, n, 0);
}

/// `SDSDOT` (`interface/sdsdot.c`): `alpha + x . y`, summed in `f64` and
/// returned as `f32`; `alpha` alone for `n <= 0`.
pub fn sdsdot(io: Io, n: isize, alpha: f32, x: [*]f32, incx: isize, y: [*]f32, incy: isize) Error!f32 {
    if (n <= 0) return alpha;
    return @floatCast(try dsdot(io, n, x, incx, y, incy) + alpha);
}

// ---- norms and sums ----------------------------------------------------------

/// `?NRM2`/`SCNRM2`/`DZNRM2` (`interface/nrm2.c`): the Euclidean norm. For
/// `incx == 0` it is `sqrt(n) * |x[0]|`, worked out here as the reference
/// does; a negative stride is the same vector backwards. Serial: chunks'
/// norms would fold with a rounding each, where `s` and `c` round once and
/// `d` and `z` keep Blue's three sums (`threading.md`).
pub fn nrm2(io: Io, comptime T: type, n: isize, x: [*]T, incx: isize) Error!Float(T) {
    _ = io;
    const F = Float(T);
    if (n <= 0) return 0;
    // the count and its root in the accumulator type: `n` may be past `F`'s range
    const A = types.Acc(F);
    const fn_: A = @floatFromInt(n);
    if (comptime !isComplex(T)) {
        if (n == 1) return @abs(x[0]);
        if (incx == 0) return @floatCast(@sqrt(fn_) * @abs(@as(A, x[0])));
    } else if (incx == 0) {
        const fr: A = @abs(x[0].re);
        const fi: A = @abs(x[0].im);
        const lo = @min(fr, fi);
        const hi = @max(fr, fi);
        if (hi == 0) return 0;
        if (hi == lo) return @floatCast(@sqrt(fn_) * @sqrt(@as(A, 2)) * hi);
        return @floatCast(@sqrt(fn_) * hi * @sqrt(1 + (lo / hi) * (lo / hi)));
    }
    try spanX(T, n, incx);
    if (accel.nrm2(T, n, x, incx)) |r| return r;
    return portable.fetch(T, .nrm2_k, blas.vtable)(n, raw(T, first(T, x, n, incx)), incx);
}

/// `?ASUM`/`SCASUM`/`DZASUM` (`interface/asum.c`): the sum of `|x|`, of
/// `|re| + |im|` for complex. 0 for `n <= 0`; the kernel gives 0 for
/// `incx <= 0`.
pub fn asum(io: Io, comptime T: type, n: isize, x: [*]T, incx: isize) Error!Float(T) {
    if (n <= 0) return 0;
    try spanX(T, n, incx);
    if (try accel.asum(io, T, n, x, incx)) |r| return r;
    return summed(io, T, .asum_k, n, x, incx);
}

/// `asum` or `sum` (`op`) over chunks of `n`, folded in chunk order.
fn summed(io: Io, comptime T: type, comptime op: portable.Op, n: isize, x: [*]T, incx: isize) Io.Cancelable!Float(T) {
    const F = Float(T);
    // `kernel/x86_64/dasum.c` and the rest: a stride `<= 0` runs serially
    if (comptime sums(T)) if (incx > 0) if (chunked(.asum, T, n)) |p| {
        const c: Pair(T, void) = .{ .x = x, .incx = incx, .y = x, .incy = incx, .e = {} };
        const S = struct {
            fn part(q: @TypeOf(c), len: isize, at: isize) Wide(F) {
                return portable.fetch(T, op, blas.vtable)(len, raw(T, shift(T, q.x, at, q.incx)), q.incx);
            }
        };
        return rounded(F, try fold(io, &p, Wide(F), c, S.part, plus(Wide(F))));
    };
    return portable.fetch(T, op, blas.vtable)(n, raw(T, x), incx);
}

/// `?SUM` (`interface/sum.c`, an OpenBLAS extension): the sum of `x`, of
/// `re + im` for complex.
pub fn sum(io: Io, comptime T: type, n: isize, x: [*]T, incx: isize) Error!Float(T) {
    if (n <= 0) return 0;
    try spanX(T, n, incx);
    return summed(io, T, .sum_k, n, x, incx);
}

// ---- extremes ----------------------------------------------------------------

pub const Extreme = extremes.Extreme;
pub const index = extremes.index;
pub const extreme = extremes.extreme;
pub const Best = extremes.Best;

// ---- rotations ---------------------------------------------------------------

/// `?ROT`, `CSROT`/`ZDROT` (`interface/rot.c`, `zrot.c`): the plane rotation
/// `x := c*x + s*y`, `y := c*y - s*x`, with real `c` and `s`. Long complex
/// vectors at unit strides take the real kernel on x86_64 (`rot_reals.zig`).
pub fn rot(io: Io, comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize, c: Float(T), s: Float(T)) Error!void {
    if (n <= 0) return;
    try spanXY(T, T, n, incx, incy);
    if (try accel.rot(io, T, n, x, incx, y, incy, c, s)) return;
    const Scalars = struct { c: Float(T), s: Float(T), reals: bool };
    // the kernel is chosen on the whole `n`, before any split: the two differ in bits
    const v: Pair(T, Scalars) = .{ .x = first(T, x, n, incx), .incx = incx, .y = first(T, y, n, incy), .incy = incy, .e = .{ .c = c, .s = s, .reals = rot_reals.routed(T, n, incx, incy) } };
    const S = struct {
        fn run(q: @TypeOf(v), len: isize, at: isize) void {
            const px = shift(T, q.x, at, q.incx);
            const py = shift(T, q.y, at, q.incy);
            if (comptime rot_reals.serves(T)) if (q.e.reals) return rot_reals.run(T, len, px, py, q.e.c, q.e.s);
            _ = portable.fetch(T, .rot_k, blas.vtable)(len, raw(T, px), q.incx, raw(T, py), q.incy, q.e.c, q.e.s);
        }
    };
    if (incx != 0 and incy != 0) if (spread(.rot, T, n)) |p| return each(io, &p, v, S.run);
    S.run(v, n, 0);
}

/// `?ROTM` (`interface/rotm.c`): the modified rotation `param` describes
/// (`param[0]` the flag). The reference passes everything to the kernel,
/// which handles negative strides itself; nothing for `n <= 0`.
pub fn rotm(io: Io, comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize, param: *const [5]T) Error!void {
    _ = io;
    comptime util.requireReal(T, "rotm");
    if (n <= 0) return;
    try spanXY(T, T, n, incx, incy);
    if (accel.rotm(T, n, x, incx, y, incy, param)) return;
    _ = portable.fetch(T, .rotm_k, blas.vtable)(n, x, incx, y, incy, @constCast(param));
}

// ---- axpy, axpby -------------------------------------------------------------

/// `?AXPY` (`interface/axpy.c`, `zaxpy.c`): `y := alpha*x + y`. Nothing to
/// do for `alpha == 0`; both strides zero is `y[0] += n*alpha*x[0]`.
pub fn axpy(io: Io, comptime T: type, n: isize, alpha: T, x: [*]T, incx: isize, y: [*]T, incy: isize) Error!void {
    return axpyImpl(io, T, false, n, alpha, x, incx, y, incy);
}

/// `?AXPYC`: `y := alpha*conj(x) + y`. Complex only.
pub fn axpyc(io: Io, comptime T: type, n: isize, alpha: T, x: [*]T, incx: isize, y: [*]T, incy: isize) Error!void {
    comptime util.requireComplex(T, "axpyc");
    return axpyImpl(io, T, true, n, alpha, x, incx, y, incy);
}

fn axpyImpl(io: Io, comptime T: type, comptime conj: bool, n: isize, alpha: T, x: [*]T, incx: isize, y: [*]T, incy: isize) Error!void {
    if (n <= 0) return;
    if (comptime isComplex(T)) {
        if (alpha.re == 0 and alpha.im == 0) return;
        if (incx == 0 and incy == 0) {
            // the reference does this unconjugated for `axpyc` too
            // the count in the accumulator type: `n` may be past `F`'s range
            const A = types.Acc(Float(T));
            const cn: A = @floatFromInt(n);
            y[0].re += @floatCast(cn * (alpha.re * x[0].re - alpha.im * x[0].im));
            y[0].im += @floatCast(cn * (alpha.im * x[0].re + alpha.re * x[0].im));
            return;
        }
        try spanXY(T, T, n, incx, incy);
        // CBLAS has no `axpyc`
        if (!conj and try accel.axpy(io, T, n, alpha, x, incx, y, incy)) return;
    } else {
        if (alpha == 0) return;
        if (incx == 0 and incy == 0) {
            y[0] += @floatCast(@as(types.Acc(T), @floatFromInt(n)) * alpha * x[0]);
            return;
        }
        try spanXY(T, T, n, incx, incy);
        if (try accel.axpy(io, T, n, alpha, x, incx, y, incy)) return;
    }
    const c: Pair(T, T) = .{ .x = first(T, x, n, incx), .incx = incx, .y = first(T, y, n, incy), .incy = incy, .e = alpha };
    const S = struct {
        fn run(q: @TypeOf(c), len: isize, at: isize) void {
            const px = shift(T, q.x, at, q.incx);
            const py = shift(T, q.y, at, q.incy);
            const op: portable.Op = comptime if (conj and isComplex(T)) .axpyc_k else .axpy_k;
            _ = @call(.auto, portable.fetch(T, op, blas.vtable), .{ len, 0, 0 } ++ util.scalarArgs(T, q.e) ++ .{ raw(T, px), q.incx, raw(T, py), q.incy, raw(T, py), 0 });
        }
    };
    if (incx != 0 and incy != 0) if (spread(.axpy, T, n)) |p| return each(io, &p, c, S.run);
    S.run(c, n, 0);
}

/// `?AXPBY` (`interface/axpby.c`, `zaxpby.c`, an OpenBLAS extension):
/// `y := alpha*x + beta*y`.
pub fn axpby(io: Io, comptime T: type, n: isize, alpha: T, x: [*]T, incx: isize, beta: T, y: [*]T, incy: isize) Error!void {
    _ = io;
    if (n <= 0) return;
    try spanXY(T, T, n, incx, incy);
    if (accel.axpby(T, n, alpha, x, incx, beta, y, incy)) return;
    const px = raw(T, first(T, x, n, incx));
    const py = raw(T, first(T, y, n, incy));
    _ = @call(.auto, portable.fetch(T, .axpby_k, blas.vtable), .{n} ++ util.scalarArgs(T, alpha) ++ .{ px, incx } ++ util.scalarArgs(T, beta) ++ .{ py, incy });
}

// ---- bfloat16 ----------------------------------------------------------------

pub const sbdot = half.sbdot;
pub const bscal = half.bscal;
pub const sbtobf16 = half.sbtobf16;
pub const bf16to = half.bf16to;

// ---- rotation set-up ---------------------------------------------------------

pub const rotg = rot_root.rotg;
pub const zrotg = rot_root.zrotg;
pub const rotmg = rot_root.rotmg;

// ---- complex magnitude -------------------------------------------------------

/// `SCABS1`/`DCABS1` (`kernel/generic/cabs.c`; `cabs.S` on x86, x86_64,
/// POWER and SPARC): `|re| + |im|`, the key the complex extremes compare,
/// not the modulus. Any complex type; infinite where the sum overflows.
pub fn cabs1(comptime T: type, z: T) Float(T) {
    comptime util.requireComplex(T, "cabs1");
    return util.abs(z.re) + util.abs(z.im);
}
