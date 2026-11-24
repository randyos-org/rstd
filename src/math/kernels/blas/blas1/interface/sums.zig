//! The level-1 interface's norms and sums: `nrm2`, `asum`, `sum` and their `.fast` forms.

const std = @import("std");
const blas = @import("../../__root__.zig");
const threading = @import("../../threading.zig");
const types = @import("../../types.zig");
const util = @import("../../util.zig");
const portable = @import("../../portable.zig");
const accel = @import("../../apple_accel.zig");
const iroot = @import("../interface.zig");
const idot = @import("dot.zig");

const Io = std.Io;
const isComplex = types.isComplex;
const Float = types.Float;
const raw = util.raw;
const Error = iroot.Error;
const Pair = iroot.Pair;
const Wide = iroot.Wide;
const chunked = iroot.chunked;
const first = iroot.first;
const fold = iroot.fold;
const plus = iroot.plus;
const rounded = iroot.rounded;
const shift = iroot.shift;
const spanX = iroot.spanX;
const sums = iroot.sums;
const Lanes = idot.Lanes;
const single = idot.single;
const wide = idot.wide;

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
/// `incx <= 0`. Single precision sums in `f64` (`Lanes`).
pub fn asum(io: Io, comptime T: type, n: isize, x: [*]T, incx: isize) Error!Float(T) {
    return asumBy(io, T, .plain, n, x, incx);
}

/// `asum` for `Summation.fast`: single precision on the core's own kernel.
pub fn asumFast(io: Io, comptime T: type, n: isize, x: [*]T, incx: isize) Error!Float(T) {
    return asumBy(io, T, .fast, n, x, incx);
}

fn asumBy(io: Io, comptime T: type, comptime lanes: Lanes, n: isize, x: [*]T, incx: isize) Error!Float(T) {
    if (n <= 0) return 0;
    try spanX(T, n, incx);
    if (comptime wide(T, lanes)) return summed(io, T, .asum_wide_k, n, x, incx);
    if (comptime !single(T)) if (try accel.asum(io, T, n, x, incx)) |r| return r;
    return summed(io, T, .asum_k, n, x, incx);
}

/// `asum` or `sum` (`op`) over chunks of `n`, folded in chunk order. The
/// serial test first, at compile time's cost: a serial call never holds a
/// `threading.Plan`, whose 536 bytes cost it a frame (issue #677).
inline fn summed(io: Io, comptime T: type, comptime op: portable.Op, n: isize, x: [*]T, incx: isize) Io.Cancelable!Float(T) {
    // `kernel/x86_64/dasum.c` and the rest: a stride `<= 0` runs serially
    if (comptime sums(T)) if (incx > 0) if (!threading.serial(.asum, T, @intCast(n))) return summedChunks(io, T, op, n, x, incx);
    return portable.fetch(T, op, blas.vtable)(n, raw(T, x), incx);
}

noinline fn summedChunks(io: Io, comptime T: type, comptime op: portable.Op, n: isize, x: [*]T, incx: isize) Io.Cancelable!Float(T) {
    const F = Float(T);
    if (chunked(.asum, T, n)) |p| {
        const c: Pair(T, void) = .{ .x = x, .incx = incx, .y = x, .incy = incx, .e = {} };
        const S = struct {
            fn part(q: @TypeOf(c), len: isize, at: isize) Wide(F) {
                return portable.fetch(T, op, blas.vtable)(len, raw(T, shift(T, q.x, at, q.incx)), q.incx);
            }
        };
        return rounded(F, try fold(io, &p, Wide(F), c, S.part, plus(Wide(F))));
    }
    return portable.fetch(T, op, blas.vtable)(n, raw(T, x), incx);
}

/// `?SUM` (`interface/sum.c`, an OpenBLAS extension): the sum of `x`, of
/// `re + im` for complex. Single precision sums in `f64` (`Lanes`).
pub fn sum(io: Io, comptime T: type, n: isize, x: [*]T, incx: isize) Error!Float(T) {
    if (n <= 0) return 0;
    try spanX(T, n, incx);
    return summed(io, T, if (comptime wide(T, .plain)) .sum_wide_k else .sum_k, n, x, incx);
}

/// `sum` for `Summation.fast`: single precision on the core's own kernel.
pub fn sumFast(io: Io, comptime T: type, n: isize, x: [*]T, incx: isize) Error!Float(T) {
    if (n <= 0) return 0;
    try spanX(T, n, incx);
    return summed(io, T, .sum_k, n, x, incx);
}
