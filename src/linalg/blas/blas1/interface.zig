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

/// Refuses an `x` that could not exist in memory (`util.fits`), so no index
/// overflows.
fn spanX(comptime T: type, n: isize, incx: isize) ArgError!void {
    if (!util.fits(n, incx, @sizeOf(T))) return error.InvalidIncx;
}

/// `spanX` for `x` of `X`, then `y` of `Y`.
fn spanXY(comptime X: type, comptime Y: type, n: isize, incx: isize, incy: isize) ArgError!void {
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
fn chunked(comptime family: threading.Family, comptime T: type, n: isize) ?threading.Plan {
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
fn fold(io: Io, p: *const threading.Plan, comptime R: type, ctx: anytype, comptime part: fn (@TypeOf(ctx), isize, isize) R, comptime combine: fn (R, R) R) Io.Cancelable!R {
    const S = struct {
        fn one(c: @TypeOf(ctx), ch: threading.Chunk) R {
            return part(c, @intCast(ch.end - ch.start), @intCast(ch.start));
        }
    };
    return threading.reduce(io, p, R, ctx, S.one, combine);
}

/// Element `at` of a vector whose element 0 is at `p` (as `first` gives
/// it), at stride `inc`.
fn shift(comptime T: type, p: [*]T, at: isize, inc: isize) [*]T {
    return @ptrCast(util.at(T, p, at * inc));
}

/// Two vectors from their element 0, and the routine's scalars.
fn Pair(comptime T: type, comptime Extra: type) type {
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

/// Which extreme: by magnitude (`|re| + |im|` for complex) or by value.
pub const Extreme = enum {
    amax,
    amin,
    max,
    min,

    /// The slot of the kernel giving the extreme's value.
    pub fn op(e: Extreme) portable.Op {
        return switch (e) {
            .amax => .amax_k,
            .amin => .amin_k,
            .max => .max_k,
            .min => .min_k,
        };
    }

    /// The slot of the kernel giving its index.
    pub fn indexOp(e: Extreme) portable.Op {
        return switch (e) {
            .amax => .iamax_k,
            .amin => .iamin_k,
            .max => .imax_k,
            .min => .imin_k,
        };
    }
};

/// `I?AMAX`, `I?AMIN`, `I?MAX`, `I?MIN` (`interface/imax.c`): the 1-based
/// index of the first extreme element; 0 for `n <= 0` (and, from the
/// kernel, for `incx <= 0`). `max` and `min` are real only.
pub fn index(io: Io, comptime T: type, comptime which: Extreme, n: isize, x: [*]T, incx: isize) Error!isize {
    comptime realOnly(T, which);
    if (n <= 0) return 0;
    try spanX(T, n, incx);
    const B = Best(T, which);
    if (incx > 0) if (chunked(.index, T, n)) |p| {
        const c: Pair(T, void) = .{ .x = x, .incx = incx, .y = x, .incy = incx, .e = {} };
        const S = struct {
            fn part(q: @TypeOf(c), len: isize, at: isize) B.Pick {
                const ret = portable.fetch(T, comptime which.indexOp(), blas.vtable)(len, raw(T, shift(T, q.x, at, q.incx)), q.incx);
                const i = at + @as(isize, @intCast(ret)) - 1;
                return B.of(shift(T, q.x, i, q.incx)[0], i);
            }
        };
        return (try fold(io, &p, B.Pick, c, S.part, B.pick)).at + 1;
    };
    const ret: BlasLong = portable.fetch(T, comptime which.indexOp(), blas.vtable)(n, raw(T, x), incx);
    return @intCast(@min(ret, n));
}

/// `?AMAX`, `?AMIN`, `?MAX`, `?MIN` (`interface/max.c`, OpenBLAS extensions):
/// the extreme value itself, `|x|` for `amax`/`amin`. 0 for `n <= 0`; for
/// `incx == 0`, `x[0]` (by magnitude for `amax`/`amin`).
pub fn extreme(io: Io, comptime T: type, comptime which: Extreme, n: isize, x: [*]T, incx: isize) Error!Float(T) {
    comptime realOnly(T, which);
    if (n <= 0) return 0;
    try spanX(T, n, incx);
    if (incx == 0) {
        if (comptime isComplex(T)) return util.abs(x[0].re) + util.abs(x[0].im);
        return if (which == .amax or which == .amin) @abs(x[0]) else x[0];
    }
    const c: Pair(T, void) = .{ .x = x, .incx = incx, .y = x, .incy = incx, .e = {} };
    const S = struct {
        fn part(q: @TypeOf(c), len: isize, at: isize) Float(T) {
            return portable.fetch(T, comptime which.op(), blas.vtable)(len, raw(T, shift(T, q.x, at, q.incx)), q.incx);
        }
    };
    if (incx > 0) if (chunked(.index, T, n)) |p| return fold(io, &p, Float(T), c, S.part, Best(T, which).value);
    return S.part(c, n, 0);
}

/// The extremes' fold, exact: the first NaN, else the first extreme, as one
/// kernel call over all of `x` finds them (`extremum.zig`). Where complex
/// keys overflow, a chunk's kernel has picked by `|re|/4 + |im|/4`, and so
/// does the fold between two such chunks.
fn Best(comptime T: type, comptime which: Extreme) type {
    const F = Float(T);
    return struct {
        /// An element's index from 0, its key and its quarter key.
        const Pick = struct { at: isize, key: F, q: F };

        fn better(k: F, best: F) bool {
            return switch (which) {
                .amax, .max => k > best,
                .amin, .min => k < best,
            };
        }

        fn value(a: F, b: F) F {
            if (a != a) return a;
            if (b != b) return b;
            return if (better(b, a)) b else a;
        }

        fn of(v: T, at: isize) Pick {
            if (comptime isComplex(T)) return .{ .at = at, .key = util.abs(v.re) + util.abs(v.im), .q = util.abs(v.re) / 4 + util.abs(v.im) / 4 };
            const k = if (which == .amax or which == .amin) @abs(v) else v;
            return .{ .at = at, .key = k, .q = k };
        }

        fn pick(a: Pick, b: Pick) Pick {
            if (a.key != a.key) return a;
            if (b.key != b.key) return b;
            if (std.math.isInf(a.key) and std.math.isInf(b.key)) return if (better(b.q, a.q)) b else a;
            return if (better(b.key, a.key)) b else a;
        }
    };
}

fn realOnly(comptime T: type, comptime which: Extreme) void {
    if (isComplex(T) and (which == .max or which == .min)) @compileError(@tagName(which) ++ ": " ++ @typeName(T) ++ " has no ordering; use ." ++ (if (which == .max) "amax" else "amin") ++ ", which compares |re| + |im|");
}

// ---- rotations ---------------------------------------------------------------

/// `?ROT`, `CSROT`/`ZDROT` (`interface/rot.c`, `zrot.c`): the plane rotation
/// `x := c*x + s*y`, `y := c*y - s*x`, with real `c` and `s`.
pub fn rot(io: Io, comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize, c: Float(T), s: Float(T)) Error!void {
    if (n <= 0) return;
    try spanXY(T, T, n, incx, incy);
    if (try accel.rot(io, T, n, x, incx, y, incy, c, s)) return;
    const v: Pair(T, [2]Float(T)) = .{ .x = first(T, x, n, incx), .incx = incx, .y = first(T, y, n, incy), .incy = incy, .e = .{ c, s } };
    const S = struct {
        fn run(q: @TypeOf(v), len: isize, at: isize) void {
            _ = portable.fetch(T, .rot_k, blas.vtable)(len, raw(T, shift(T, q.x, at, q.incx)), q.incx, raw(T, shift(T, q.y, at, q.incy)), q.incy, q.e[0], q.e[1]);
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

/// `SBDOT` (`interface/bf16dot.c`): the dot of two bfloat16 vectors, 0 for
/// `n <= 0`, returned in `f32`.
pub fn sbdot(io: Io, n: isize, x: [*]BF16, incx: isize, y: [*]BF16, incy: isize) Error!f32 {
    _ = io;
    if (n <= 0) return 0;
    try spanXY(BF16, BF16, n, incx, incy);
    return blas.kernel(.sbdot_k)(n, first(BF16, x, n, incx), incx, first(BF16, y, n, incy), incy);
}

/// `BSCAL` (`interface/scal.c` built `BGEMM`): `x := alpha*x` over bfloat16.
/// Nothing for `incx <= 0`, and nothing to do for `alpha == 1`; like `scal`
/// it passes the kernel flag 1.
pub fn bscal(io: Io, n: isize, alpha: BF16, x: [*]BF16, incx: isize) Error!void {
    _ = io;
    if (incx <= 0 or n <= 0) return;
    try spanX(BF16, n, incx);
    if (bf.scalarToF32(alpha) == 1) return;
    _ = blas.kernel(.bscal_k)(n, 0, 0, alpha, x, incx, x, 0, x, 1);
}

/// `SBSTOBF16` and `SBDTOBF16` (`interface/tobf16.c`): `out := bf16(in)`,
/// rounding to nearest even; `f32` or `f64` in. The errors name `in` as
/// `x` and `out` as `y`.
pub fn sbtobf16(io: Io, comptime F: type, n: isize, in: [*]F, inc_in: isize, out: [*]BF16, inc_out: isize) Error!void {
    _ = io;
    if (n <= 0) return;
    const kern = if (comptime F == f32)
        blas.kernel(.sbstobf16_k)
    else if (comptime F == f64)
        blas.kernel(.sbdtobf16_k)
    else
        @compileError("sbtobf16: " ++ @typeName(F) ++ " is not f32 or f64");
    try spanXY(F, BF16, n, inc_in, inc_out);
    kern(n, first(F, in, n, inc_in), inc_in, first(BF16, out, n, inc_out), inc_out);
}

/// `SBF16TOS` and `DBF16TOD` (`interface/bf16to.c`): `out := in` widened to
/// `f32` or `f64`; errors as `sbtobf16`'s.
pub fn bf16to(io: Io, comptime F: type, n: isize, in: [*]BF16, inc_in: isize, out: [*]F, inc_out: isize) Error!void {
    _ = io;
    if (n <= 0) return;
    const kern = if (comptime F == f32)
        blas.kernel(.sbf16tos_k)
    else if (comptime F == f64)
        blas.kernel(.dbf16tod_k)
    else
        @compileError("bf16to: " ++ @typeName(F) ++ " is not f32 or f64");
    try spanXY(BF16, F, n, inc_in, inc_out);
    kern(n, first(BF16, in, n, inc_in), inc_in, first(F, out, n, inc_out), inc_out);
}

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

// ---- tests -------------------------------------------------------------------

const check = @import("../_test/check.zig");
const reference = @import("../_test/reference.zig");
const blas2_interface = @import("../blas2/interface.zig");
const Complex32 = types.Complex32;
const Complex64 = types.Complex64;

const reals = .{ f32, f64 } ++ check.x_real;
const all = check.all_types;

const zero = util.zero;

fn same(comptime T: type, want: []const T, got: []const T) !void {
    try check.expectCloseSlices(T, want, got, check.tolerance(Float(T), 4));
}

fn sameInt(want: isize, got: isize) !void {
    if (want == got) return;
    check.say("got {d}, want {d}: ", .{ got, want });
    return error.TestExpectedEqual;
}

test "interface: n <= 0 changes nothing and returns zero" {
    inline for (all) |T| {
        var x: [4]T = undefined;
        var y: [4]T = undefined;
        check.fill(T, &x, 1);
        check.fill(T, &y, 2);
        const x0 = x;
        const y0 = y;
        for ([_]isize{ 0, -1 }) |n| {
            try copy(check.io, T, n, &x, 1, &y, 1);
            try swap(check.io, T, n, &x, 1, &y, 1);
            try scal(check.io, T, n, zero(T), &x, 1);
            // the reference has no extended axpby
            if (comptime !types.isExtended(Float(T))) try axpby(check.io, T, n, zero(T), &x, 1, zero(T), &y, 1);
            try rot(check.io, T, n, &x, 1, &y, 1, 0, 1);
            try same(T, &x0, &x);
            try same(T, &y0, &y);
            try check.expectClose(T, zero(T), try dot(check.io, T, n, &x, 1, &y, 1), 0);
            try check.expectClose(Float(T), 0, try nrm2(check.io, T, n, &x, 1), 0);
            try check.expectClose(Float(T), 0, try asum(check.io, T, n, &x, 1), 0);
            try check.expectClose(Float(T), 0, try sum(check.io, T, n, &x, 1), 0);
            try check.expectClose(Float(T), 0, try extreme(check.io, T, .amax, n, &x, 1), 0);
            try sameInt(0, try index(check.io, T, .amax, n, &x, 1));
        }
    }
    try check.expectClose(f32, 2.5, try sdsdot(check.io, 0, 2.5, undefined, 1, undefined, 1), 0);
}

test "interface: a vector past the address space is its stride's error, x's first" {
    const testing = std.testing;
    // two elements of 2 bytes or more at this stride overflow
    const big = std.math.maxInt(isize) / 2;
    const bad = [_]isize{ big, -big, std.math.minInt(isize) };
    inline for (all) |T| {
        var x: [2]T = undefined;
        var y: [2]T = undefined;
        check.fill(T, &x, 1);
        check.fill(T, &y, 2);
        const two: T = if (comptime isComplex(T)) .init(2, 0) else 2;
        const extended = comptime types.isExtended(Float(T));
        for (bad) |b| {
            const cases = [_]struct { ix: isize, iy: isize, err: ArgError }{
                .{ .ix = b, .iy = 1, .err = error.InvalidIncx },
                .{ .ix = 1, .iy = b, .err = error.InvalidIncy },
                .{ .ix = b, .iy = b, .err = error.InvalidIncx },
            };
            for (cases) |c| {
                try testing.expectError(c.err, copy(check.io, T, 2, &x, c.ix, &y, c.iy));
                try testing.expectError(c.err, swap(check.io, T, 2, &x, c.ix, &y, c.iy));
                try testing.expectError(c.err, dot(check.io, T, 2, &x, c.ix, &y, c.iy));
                try testing.expectError(c.err, rot(check.io, T, 2, &x, c.ix, &y, c.iy, 0.6, 0.8));
                try testing.expectError(c.err, axpy(check.io, T, 2, two, &x, c.ix, &y, c.iy));
                if (comptime !extended) try testing.expectError(c.err, axpby(check.io, T, 2, two, &x, c.ix, two, &y, c.iy));
                if (comptime isComplex(T)) {
                    try testing.expectError(c.err, dotc(check.io, T, 2, &x, c.ix, &y, c.iy));
                    try testing.expectError(c.err, axpyc(check.io, T, 2, two, &x, c.ix, &y, c.iy));
                } else {
                    try testing.expectError(c.err, rotm(check.io, T, 2, &x, c.ix, &y, c.iy, &.{ -1, 2, 3, -1, 0.5 }));
                }
            }
            // `scal` and `bscal` skip a negative stride before checking it
            if (b > 0) try testing.expectError(error.InvalidIncx, scal(check.io, T, 2, two, &x, b));
            try testing.expectError(error.InvalidIncx, nrm2(check.io, T, 2, &x, b));
            try testing.expectError(error.InvalidIncx, asum(check.io, T, 2, &x, b));
            try testing.expectError(error.InvalidIncx, sum(check.io, T, 2, &x, b));
            try testing.expectError(error.InvalidIncx, index(check.io, T, .amax, 2, &x, b));
            try testing.expectError(error.InvalidIncx, extreme(check.io, T, .amax, 2, &x, b));
            // the early returns come first
            try copy(check.io, T, 0, &x, b, &y, b);
            try testing.expectEqual(zero(T), try dot(check.io, T, 0, &x, b, &y, b));
            if (comptime !isComplex(T)) try rotm(check.io, T, 0, &x, b, &y, b, &.{ -1, 2, 3, -1, 0.5 });
        }
        // nothing was written
        var x0: [2]T = undefined;
        var y0: [2]T = undefined;
        check.fill(T, &x0, 1);
        check.fill(T, &y0, 2);
        try same(T, &x0, &x);
        try same(T, &y0, &y);
    }
    var f = [_]f32{ 1, 2 };
    var d = [_]f64{ 1, 2 };
    var h = [_]BF16{ b16(1), b16(2) };
    for (bad) |b| {
        try testing.expectError(error.InvalidIncx, dsdot(check.io, 2, &f, b, &f, 1));
        try testing.expectError(error.InvalidIncy, dsdot(check.io, 2, &f, 1, &f, b));
        try testing.expectError(error.InvalidIncx, sdsdot(check.io, 2, 0.5, &f, b, &f, b));
        try testing.expectError(error.InvalidIncx, sbdot(check.io, 2, &h, b, &h, 1));
        try testing.expectError(error.InvalidIncy, sbdot(check.io, 2, &h, 1, &h, b));
        try testing.expectError(error.InvalidIncx, sbtobf16(check.io, f32, 2, &f, b, &h, 1));
        try testing.expectError(error.InvalidIncy, sbtobf16(check.io, f64, 2, &d, 1, &h, b));
        try testing.expectError(error.InvalidIncx, bf16to(check.io, f32, 2, &h, b, &f, 1));
        try testing.expectError(error.InvalidIncy, bf16to(check.io, f64, 2, &h, 1, &d, b));
        if (b > 0) try testing.expectError(error.InvalidIncx, bscal(check.io, 2, b16(2), &h, b));
    }
}

test "interface: a negative stride walks the vector from its end" {
    inline for (all) |T| {
        const n = 7;
        var x: [n * 2]T = undefined;
        var y: [n * 3]T = undefined;
        check.fill(T, &x, 3);
        check.fill(T, &y, 4);
        // `xr` is x as stride -2 sees it, stored forwards at stride 1
        var xr: [n]T = undefined;
        for (0..n) |i| xr[i] = x[(n - 1 - i) * 2];
        var y1 = y;
        var y2 = y;
        try copy(check.io, T, n, &x, -2, &y1, 3);
        try copy(check.io, T, n, &xr, 1, &y2, 3);
        try same(T, &y2, &y1);
        try check.expectClose(T, try dot(check.io, T, n, &xr, 1, &y, 3), try dot(check.io, T, n, &x, -2, &y, 3), check.tolerance(Float(T), 4 * n));
        try check.expectClose(Float(T), try nrm2(check.io, T, n, &xr, 1), try nrm2(check.io, T, n, &x, -2), check.tolerance(Float(T), 4 * n));
    }
}

test "interface: scal does nothing for incx <= 0 or alpha 1; a zero alpha keeps NaN" {
    inline for (all) |T| {
        var x: [4]T = undefined;
        check.fill(T, &x, 5);
        const x0 = x;
        const two: T = if (comptime isComplex(T)) .init(2, 0) else 2;
        const one: T = if (comptime isComplex(T)) .init(1, 0) else 1;
        try scal(check.io, T, 4, two, &x, 0);
        try scal(check.io, T, 4, two, &x, -1);
        try scal(check.io, T, 4, one, &x, 1);
        try same(T, &x0, &x);
        const nan = std.math.nan(Float(T));
        x[1] = if (comptime isComplex(T)) .init(nan, nan) else nan;
        try scal(check.io, T, 4, zero(T), &x, 1);
        const at1 = if (comptime isComplex(T)) x[1].re else x[1];
        if (!std.math.isNan(at1)) return error.TestExpectedNan;
        x[1] = zero(T);
        try same(T, &@as([4]T, @splat(zero(T))), &x);
    }
    var z = [_]Complex64{ .init(1, 2), .init(-3, 4) };
    try scalReal(check.io, Complex64, 2, 2, &z, 1);
    try same(Complex64, &.{ .init(2, 4), .init(-6, 8) }, &z);
}

test "interface: dotc conjugates x, and dsdot/sdsdot sum in f64" {
    var x = [_]Complex32{ .init(1, 2), .init(3, -1) };
    var y = [_]Complex32{ .init(2, 1), .init(-1, 1) };
    // conj(1+2i)(2+i) + conj(3-i)(-1+i) = (4-3i) + (-4+2i)
    try check.expectClose(Complex32, .init(0, -1), try dotc(check.io, Complex32, 2, &x, 1, &y, 1), 1e-6);
    // (1+2i)(2+i) + (3-i)(-1+i) = 5i + (-2+4i)
    try check.expectClose(Complex32, .init(-2, 9), try dot(check.io, Complex32, 2, &x, 1, &y, 1), 1e-6);
    // 1e8 + 1 - 1e8 is 1 in f64 but 0 in f32
    var a = [_]f32{ 1e4, 1, -1e4 };
    var b = [_]f32{ 1e4, 1, 1e4 };
    try check.expectClose(f64, 1, try dsdot(check.io, 3, &a, 1, &b, 1), 0);
    try check.expectClose(f32, 1.5, try sdsdot(check.io, 3, 0.5, &a, 1, &b, 1), 0);
}

test "interface: nrm2 of a zero stride is sqrt(n) * |x[0]|" {
    var x = [_]f64{-3};
    try check.expectClose(f64, 3 * @sqrt(4.0), try nrm2(check.io, f64, 4, &x, 0), 1e-12);
    var z = [_]Complex64{.init(3, -4)};
    try check.expectClose(f64, 5 * @sqrt(4.0), try nrm2(check.io, Complex64, 4, &z, 0), 1e-12);
    var w = [_]Complex64{.init(2, -2)};
    try check.expectClose(f64, @sqrt(9.0) * @sqrt(2.0) * 2, try nrm2(check.io, Complex64, 9, &w, 0), 1e-12);
}

test "interface: extremes, by magnitude and by value" {
    var x = [_]f64{ 3, -7, 7, 1, -2 };
    try sameInt(2, try index(check.io, f64, .amax, 5, &x, 1));
    try sameInt(4, try index(check.io, f64, .amin, 5, &x, 1));
    try sameInt(3, try index(check.io, f64, .max, 5, &x, 1));
    try sameInt(2, try index(check.io, f64, .min, 5, &x, 1));
    try sameInt(0, try index(check.io, f64, .amax, 5, &x, 0));
    try check.expectClose(f64, 7, try extreme(check.io, f64, .amax, 5, &x, 1), 0);
    try check.expectClose(f64, 1, try extreme(check.io, f64, .amin, 5, &x, 1), 0);
    try check.expectClose(f64, 7, try extreme(check.io, f64, .max, 5, &x, 1), 0);
    try check.expectClose(f64, -7, try extreme(check.io, f64, .min, 5, &x, 1), 0);
    try check.expectClose(f64, -7, try extreme(check.io, f64, .max, 5, x[1..], 0), 0);
    try check.expectClose(f64, 7, try extreme(check.io, f64, .amax, 5, x[1..], 0), 0);
    var z = [_]Complex32{ .init(1, 1), .init(-3, 2), .init(0, 4) };
    try sameInt(2, try index(check.io, Complex32, .amax, 3, &z, 1));
    try check.expectClose(f32, 5, try extreme(check.io, Complex32, .amax, 3, &z, 1), 0);
}

test "interface: rot, rotm and axpby match a plain loop" {
    inline for (reals) |T| {
        const n = 9;
        var x: [n * 2]T = undefined;
        var y: [n]T = undefined;
        check.fill(T, &x, 6);
        check.fill(T, &y, 7);
        var wx = x;
        var wy = y;
        for (0..n) |i| {
            const a = wx[(n - 1 - i) * 2];
            const b = wy[i];
            wx[(n - 1 - i) * 2] = 0.6 * a + 0.8 * b;
            wy[i] = 0.6 * b - 0.8 * a;
        }
        try rot(check.io, T, n, &x, -2, &y, 1, 0.6, 0.8);
        try same(T, &wx, &x);
        try same(T, &wy, &y);

        const param = [5]T{ -1, 2, 3, -1, 0.5 };
        var mx = x;
        var my = y;
        for (0..n) |i| {
            const a = mx[i * 2];
            const b = my[i];
            mx[i * 2] = 2 * a + -1 * b;
            my[i] = 3 * a + 0.5 * b;
        }
        try rotm(check.io, T, n, &x, 2, &y, 1, &param);
        try same(T, &mx, &x);
        try same(T, &my, &y);

        if (comptime types.isExtended(T)) continue;
        var ay = y;
        for (0..n) |i| ay[i] = 1.5 * x[i * 2] - 0.5 * y[i];
        try axpby(check.io, T, n, 1.5, &x, 2, -0.5, &y, 1);
        try same(T, &ay, &y);
    }
}

test "interface: axpy returns for n <= 0 or zero alpha, and collapses zero strides" {
    var x = [_]f64{ 1, 2 };
    var y = [_]f64{ 5, 6 };
    try axpy(check.io, f64, 0, 1, &x, 1, &y, 1);
    try axpy(check.io, f64, -3, 1, &x, 1, &y, 1);
    try axpy(check.io, f64, 2, 0, &x, 1, &y, 1);
    try same(f64, &.{ 5, 6 }, &y);
    var cx = [_]Complex64{.init(1, 1)};
    var cy = [_]Complex64{.init(5, 6)};
    try axpy(check.io, Complex64, 1, .init(0, 0), &cx, 1, &cy, 1);
    try same(Complex64, &.{.init(5, 6)}, &cy);

    var x1 = [_]f64{3};
    var y1 = [_]f64{1};
    try axpy(check.io, f64, 4, 2, &x1, 0, &y1, 0);
    try same(f64, &.{25}, &y1);
    var cz = [_]Complex64{.init(0, 0)};
    try axpy(check.io, Complex64, 2, .init(2, 3), &cx, 0, &cz, 0);
    // 2 * ((2+3i)(1+i)) = 2 * (-1+5i)
    try same(Complex64, &.{.init(-2, 10)}, &cz);
}

test "interface: axpy with negative strides, and axpyc" {
    var x = [_]f64{ 1, 2, 3 };
    var y = [_]f64{ 0, 0, 0 };
    try axpy(check.io, f64, 3, 1, &x, -1, &y, 1);
    try same(f64, &.{ 3, 2, 1 }, &y);
    // both negative: the same pairing as both positive
    var y2 = [_]f64{ 0, 0, 0 };
    try axpy(check.io, f64, 3, 1, &x, -1, &y2, -1);
    try same(f64, &x, &y2);
    var cx = [_]Complex32{ .init(1, 0), .init(2, 0), .init(3, 0) };
    var cy = [_]Complex32{ .init(0, 0), .init(0, 0), .init(0, 0) };
    try axpy(check.io, Complex32, 3, .init(1, 0), &cx, 1, &cy, -1);
    try same(Complex32, &.{ .init(3, 0), .init(2, 0), .init(1, 0) }, &cy);

    var zx = [_]Complex64{.init(1, 1)};
    var zy = [_]Complex64{.init(0, 0)};
    try axpyc(check.io, Complex64, 1, .init(2, 3), &zx, 1, &zy, 1);
    // (2+3i) * conj(1+i) = 5 + i
    try same(Complex64, &.{.init(5, 1)}, &zy);

    var qx = [_]XDouble{ 1, 2, 3 };
    var qy = [_]XDouble{ 1, 1, 1 };
    try axpy(check.io, XDouble, 3, 2, &qx, 1, &qy, 1);
    try same(XDouble, &.{ 3, 5, 7 }, &qy);
}

/// Level 1 over `T`, which has no vtable slots, against an `f128` reference:
/// `x` at stride 2 and `y` at stride -1.
fn anyFloat(comptime T: type) !void {
    const W = reference.Wide(T);
    const n = 37;
    const ulps = 2 * (n + 2);
    var x: [2 * n]T = undefined;
    var y: [n]T = undefined;
    check.fill(T, &x, 11);
    check.fill(T, &y, 12);
    var wx: [n]W = undefined;
    var wy: [n]W = undefined;
    for (0..n) |i| {
        wx[i] = reference.wide(T, x[2 * i]);
        wy[i] = reference.wide(T, y[n - 1 - i]);
    }
    const alpha: T = if (comptime isComplex(T)) .init(1.25, -0.5) else 1.25;
    const beta: T = if (comptime isComplex(T)) .init(-0.75, 0.25) else -0.75;
    const wa = reference.wide(T, alpha);
    const wb = reference.wide(T, beta);

    // reductions
    var d = zero(W);
    var dc = zero(W);
    var dscale: f128 = 0;
    var ssq: f128 = 0;
    var s1: f128 = 0;
    var best: usize = 0;
    for (0..n) |i| {
        d = reference.addW(d, reference.mulW(wx[i], wy[i]));
        dc = reference.addW(dc, reference.mulW(reference.conjW(wx[i]), wy[i]));
        dscale += reference.abs1(wx[i]) * reference.abs1(wy[i]);
        ssq += if (comptime isComplex(T)) wx[i].re * wx[i].re + wx[i].im * wx[i].im else wx[i] * wx[i];
        s1 += reference.abs1(wx[i]);
        if (reference.abs1(wx[i]) > reference.abs1(wx[best])) best = i;
    }
    try reference.expectNear(T, d, try dot(check.io, T, n, &x, 2, &y, -1), ulps, dscale);
    if (comptime isComplex(T)) try reference.expectNear(T, dc, try dotc(check.io, T, n, &x, 2, &y, -1), ulps, dscale);
    const norm = reference.sqrtW(ssq);
    try reference.expectNear(Float(T), norm, try nrm2(check.io, T, n, &x, 2), ulps + reference.rootUlps(Float(T)), norm);
    try reference.expectNear(Float(T), s1, try asum(check.io, T, n, &x, 2), ulps, s1);
    try reference.expectNear(Float(T), reference.abs1(wx[best]), try extreme(check.io, T, .amax, n, &x, 2), 4, s1);
    // the first of near-equal keys may differ; its key must be the largest
    const at = try index(check.io, T, .amax, n, &x, 2);
    try reference.expectNear(Float(T), reference.abs1(wx[best]), @floatCast(reference.abs1(wx[@intCast(at - 1)])), 4, reference.abs1(wx[best]));

    // updates, element by element
    var y2 = y;
    try axpy(check.io, T, n, alpha, &x, 2, &y2, -1);
    for (0..n) |i| {
        const want = reference.addW(reference.mulW(wa, wx[i]), wy[i]);
        try reference.expectNear(T, want, y2[n - 1 - i], 8, reference.abs1(wa) * reference.abs1(wx[i]) + reference.abs1(wy[i]));
    }
    y2 = y;
    try axpby(check.io, T, n, alpha, &x, 2, beta, &y2, -1);
    for (0..n) |i| {
        const want = reference.addW(reference.mulW(wa, wx[i]), reference.mulW(wb, wy[i]));
        try reference.expectNear(T, want, y2[n - 1 - i], 8, reference.abs1(wa) * reference.abs1(wx[i]) + reference.abs1(wb) * reference.abs1(wy[i]));
    }
    var x2 = x;
    try scal(check.io, T, n, alpha, &x2, 2);
    for (0..n) |i| try reference.expectNear(T, reference.mulW(wa, wx[i]), x2[2 * i], 8, reference.abs1(wa) * reference.abs1(wx[i]));
    const c: Float(T) = 0.6;
    const sn: Float(T) = 0.8;
    x2 = x;
    y2 = y;
    try rot(check.io, T, n, &x2, 2, &y2, -1, c, sn);
    const wc = reference.wide(Float(T), c);
    const ws = reference.wide(Float(T), sn);
    for (0..n) |i| {
        const mag = reference.abs1(wx[i]) + reference.abs1(wy[i]);
        const rx = if (comptime isComplex(T)) W.init(wc * wx[i].re + ws * wy[i].re, wc * wx[i].im + ws * wy[i].im) else wc * wx[i] + ws * wy[i];
        const ry = if (comptime isComplex(T)) W.init(wc * wy[i].re - ws * wx[i].re, wc * wy[i].im - ws * wx[i].im) else wc * wy[i] - ws * wx[i];
        try reference.expectNear(T, rx, x2[2 * i], 8, mag);
        try reference.expectNear(T, ry, y2[n - 1 - i], 8, mag);
    }
    x2 = x;
    y2 = y;
    try swap(check.io, T, n, &x2, 2, &y2, -1);
    try copy(check.io, T, n, &x2, 2, &y2, -1);
    for (0..n) |i| {
        try reference.expectNear(T, wy[i], x2[2 * i], 0, 0);
        try reference.expectNear(T, wy[i], y2[n - 1 - i], 0, 0);
    }
    if (comptime isComplex(T)) return;

    var ps = [5]T{ -1, 2, 3, -1, 0.5 };
    x2 = x;
    y2 = y;
    try rotm(check.io, T, n, &x2, 2, &y2, -1, &ps);
    for (0..n) |i| {
        const mag = 2 * @abs(wx[i]) + @abs(wy[i]);
        try reference.expectNear(T, 2 * wx[i] - wy[i], x2[2 * i], 8, mag);
        try reference.expectNear(T, 3 * wx[i] + 0.5 * wy[i], y2[n - 1 - i], 8, mag);
    }
    try reference.expectNear(T, sumOf(T, &wx), try sum(check.io, T, n, &x, 2), ulps, s1);
    var ga: T = 3;
    var gb: T = 4;
    var gc: T = undefined;
    var gs: T = undefined;
    rotg(T, &ga, &gb, &gc, &gs);
    try reference.expectNear(T, 5, ga, 4, 5);
    try reference.expectNear(T, 0.6, gc, 4, 1);
    try reference.expectNear(T, 0.8, gs, 4, 1);
}

fn sumOf(comptime T: type, w: []const reference.Wide(T)) f128 {
    var acc: f128 = 0;
    for (w) |v| acc += v;
    return acc;
}

test "interface: level 1 at any float type runs the portable kernels built for it" {
    inline for (check.any_floats) |F| {
        inline for (.{ F, std.math.Complex(F) }) |T| anyFloat(T) catch |e| {
            check.say("{s}\n", .{@typeName(T)});
            return e;
        };
    }
}

test "interface: f16 reductions accumulate wide, and counts past its range stay finite" {
    const H = f16;
    const C = std.math.Complex(H);
    const al = std.testing.allocator;
    // past 2048 an f16 sum of ones stalls; past 65504 it overflows
    const n = 5000;
    const ones = try al.alloc(H, n);
    defer al.free(ones);
    @memset(ones, 1);
    try check.expectClose(H, 5000, try asum(check.io, H, n, ones.ptr, 1), 0);
    try check.expectClose(H, 5000, try sum(check.io, H, n, ones.ptr, 1), 0);
    try check.expectClose(H, 5000, try dot(check.io, H, n, ones.ptr, 1, ones.ptr, 1), 0);
    const big = try al.alloc(H, 4000);
    defer al.free(big);
    for (big, 0..) |*v, i| v.* = if (i < 2000) 40 else -40;
    try check.expectClose(H, 0, try sum(check.io, H, 4000, big.ptr, 1), 0);
    const cones = try al.alloc(C, 3000);
    defer al.free(cones);
    @memset(cones, .init(1, 0));
    try check.expectClose(C, .init(3000, 0), try dot(check.io, C, 3000, cones.ptr, 1, cones.ptr, 1), 0);
    try check.expectClose(C, .init(3000, 0), try dotc(check.io, C, 3000, cones.ptr, 1, cones.ptr, 1), 0);
    // gemv's transposed sums, real and complex
    var y = [_]H{0};
    try blas2_interface.gemv(al, check.io, H, .t, n, 1, 1, ones.ptr, n, ones.ptr, 1, 0, &y, 1);
    try check.expectClose(H, 5000, y[0], 0);
    var cy = [_]C{.init(0, 0)};
    try blas2_interface.gemv(al, check.io, C, .c, 3000, 1, .init(1, 0), cones.ptr, 3000, cones.ptr, 1, .init(0, 0), &cy, 1);
    try check.expectClose(C, .init(3000, 0), cy[0], 0);
    // a zero stride: sqrt(n) |x[0]| with n past f16's range
    var x = [_]H{0.01};
    const want: f64 = @sqrt(70000.0) * @as(f64, @as(H, 0.01));
    try check.expectClose(f64, want, try nrm2(check.io, H, 70000, &x, 0), 4e-3);
    var cx = [_]C{.init(0.01, 0)};
    try check.expectClose(f64, want, try nrm2(check.io, C, 70000, &cx, 0), 4e-3);
    var yz = [_]H{0};
    try axpy(check.io, H, 70000, 1, &x, 0, &yz, 0);
    try check.expectClose(f64, 700, yz[0], 1);
}

fn b16(v: f32) BF16 {
    return bf.toBf16(v);
}

test "interface: sbdot, the conversions and bscal do nothing for n <= 0" {
    var x = [_]BF16{ b16(1), b16(2) };
    var f = [_]f32{ 7, 8 };
    var d = [_]f64{ 7, 8 };
    for ([_]isize{ 0, -3 }) |n| {
        try std.testing.expectEqual(@as(f32, 0), try sbdot(check.io, n, &x, 1, &x, 1));
        try sbtobf16(check.io, f32, n, &f, 1, &x, 1);
        try sbtobf16(check.io, f64, n, &d, 1, &x, 1);
        try bf16to(check.io, f32, n, &x, 1, &f, 1);
        try bf16to(check.io, f64, n, &x, 1, &d, 1);
        try bscal(check.io, n, b16(0), &x, 1);
    }
    try std.testing.expectEqual(@as(f32, 7), f[0]);
    try std.testing.expectEqual(@as(f64, 8), d[1]);
    try std.testing.expectEqual(@as(f32, 2), bf.fromBf16(x[1]));
}

test "interface: a negative stride reads and writes bfloat16 vectors from their ends" {
    const n = 5;
    var f: [n * 2]f32 = undefined;
    check.fill(f32, &f, 21);
    var h: [n * 3]BF16 = @splat(b16(0));
    // stride -3 puts element 0 last: the last element of the buffer region
    try sbtobf16(check.io, f32, n, &f, 2, &h, -3);
    for (0..n) |i| {
        try std.testing.expectEqual(@backingInt(b16(f[i * 2])), @backingInt(h[(n - 1 - i) * 3]));
    }
    var back: [n]f64 = undefined;
    try bf16to(check.io, f64, n, &h, -3, &back, 1);
    for (0..n) |i| try std.testing.expectEqual(@as(f64, bf.fromBf16(b16(f[i * 2]))), back[i]);
    // a reversed read of the same storage dots against itself forwards
    var fwd: [n]BF16 = undefined;
    for (0..n) |i| fwd[i] = h[(n - 1 - i) * 3];
    try std.testing.expectEqual(try sbdot(check.io, n, &fwd, 1, &fwd, 1), try sbdot(check.io, n, &h, -3, &h, -3));
    var want: f64 = 0;
    for (0..n) |i| want += @as(f64, bf.fromBf16(fwd[i])) * bf.fromBf16(fwd[i]);
    try check.expectClose(f32, @floatCast(want), try sbdot(check.io, n, &fwd, 1, &fwd, 1), check.tolerance(f32, 4 * n));
}

test "interface: bscal skips alpha 1 and a non-positive stride, and scales the rest" {
    var x = [_]BF16{ b16(1.5), b16(-2), b16(4) };
    try bscal(check.io, 3, b16(1), &x, 1);
    try bscal(check.io, 3, b16(2), &x, 0);
    try bscal(check.io, 3, b16(2), &x, -1);
    try std.testing.expectEqual(@as(f32, -2), bf.fromBf16(x[1]));
    try bscal(check.io, 3, b16(2), &x, 1);
    try std.testing.expectEqual(@as(f32, 3), bf.fromBf16(x[0]));
    try std.testing.expectEqual(@as(f32, -4), bf.fromBf16(x[1]));
    try std.testing.expectEqual(@as(f32, 8), bf.fromBf16(x[2]));
    // flag 1: a zero alpha zeroes the finite elements and keeps NaN
    x[2] = b16(std.math.nan(f32));
    try bscal(check.io, 3, b16(0), &x, 1);
    try std.testing.expectEqual(@as(f32, 0), bf.fromBf16(x[0]));
    try std.testing.expect(std.math.isNan(bf.fromBf16(x[2])));
}

test "interface: cabs1 is |re| + |im|, overflowing to infinity" {
    inline for (.{ std.math.Complex(f16), Complex32, Complex64 } ++ check.x_complex) |C| {
        const F = Float(C);
        const inf = std.math.inf(F);
        const max = std.math.floatMax(F);
        const tiny = std.math.floatTrueMin(F);
        const cases = [_][3]F{
            .{ 3, -4, 7 },
            .{ -0.5, 0.25, 0.75 },
            .{ -0.0, -0.0, 0 },
            .{ tiny, -tiny, 2 * tiny },
            .{ max, max, inf },
            .{ -inf, 1, inf },
        };
        for (cases) |k| {
            const got = cabs1(C, .init(k[0], k[1]));
            try std.testing.expectEqual(k[2], got);
            try std.testing.expect(!std.math.signbit(got));
        }
        try std.testing.expect(std.math.isNan(cabs1(C, .init(1, std.math.nan(F)))));
        try std.testing.expect(std.math.isNan(cabs1(C, .init(inf, -std.math.nan(F)))));
    }
}
