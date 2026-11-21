//! The moments of a slice: mean, variance and standard deviation, through
//! the kernel layer where it has the sums (`kernels.blas.sumWith`,
//! `kernels.stats.deviationSums`).

const std = @import("std");
const blas = @import("../math/kernels/blas/__root__.zig");
const blas_types = @import("../math/kernels/blas/types.zig");
const kernels_stats = @import("../math/kernels/stats/__root__.zig");
const util = @import("util.zig");

const Io = std.Io;
const Complex = std.math.Complex;
const Wide = util.Wide;
const Elem = util.Elem;
const isComplex = util.isComplex;
const Real = util.Real;
const to = util.to;
const checkPair = util.checkPair;
const checkFloat = util.checkFloat;

/// The variance's divisor: `n - 1` (MATLAB's default) or `n`. A sample of
/// one has variance 0 either way, as MATLAB's `var`.
pub const Normalize = enum { sample, population };

pub const MomentError = error{ EmptyInput, Canceled };

/// How the sums are taken (`blas.Summation`): as the BLAS's kernels do,
/// compensated (Neumaier's), or `.fast`, `f32` lanes for single precision;
/// the loops here, for integers and mixed types, sum it as `.plain`.
pub const Summation = blas.Summation;

/// The `*With` functions' options; each default is the plain function's.
pub const MomentOptions = struct {
    summation: Summation = .plain,
};

/// The mean of `xs` (a slice or array of integers, floats or `Complex`), as
/// `R`: a float for real elements, a `Complex` for complex ones. Floats of
/// `R`'s own type sum through the BLAS `sum`.
pub fn mean(io: Io, comptime R: type, xs: anytype) MomentError!R {
    return meanWith(io, R, xs, .{});
}

/// `mean`, summed as `options` says; complex elements of `R`'s own type as
/// two strided BLAS sums over their interleaved parts.
pub fn meanWith(io: Io, comptime R: type, xs: anytype, options: MomentOptions) MomentError!R {
    const E = Elem(@TypeOf(xs));
    comptime checkPair(R, E);
    if (xs.len == 0) return error.EmptyInput;
    const s = options.summation;
    const n: Real(R) = @floatFromInt(xs.len);
    if (comptime isComplex(E)) {
        if (comptime E == R and kernelFloat(Real(R))) {
            const F = Real(R);
            const sl: []const E = xs;
            const fl = @as([*]const F, @ptrCast(sl.ptr))[0 .. 2 * sl.len];
            const re = try blas.sumWith(io, F, .{ .data = fl, .n = sl.len, .stride = 2 }, s);
            const im = try blas.sumWith(io, F, .{ .data = fl[1..], .n = sl.len, .stride = 2 }, s);
            return .{ .re = re / n, .im = im / n };
        }
        const W = Wide(Real(R));
        var re: Acc(W) = .{};
        var im: Acc(W) = .{};
        for (xs) |x| {
            re.add(to(W, x.re), s);
            im.add(to(W, x.im), s);
        }
        return .{ .re = @as(Real(R), @floatCast(re.value())) / n, .im = @as(Real(R), @floatCast(im.value())) / n };
    }
    return try total(io, R, xs, s) / n;
}

/// The corrected two-pass variance (Chan, Golub and LeVeque): the mean `m`,
/// then `(Σd² − |Σd|²/n) / (n − 1)` for `d = x − m` (`|x − m|²` for complex
/// `xs`). `R` is a float.
pub fn variance(io: Io, comptime R: type, xs: anytype, normalize: Normalize) MomentError!R {
    return varianceWith(io, R, xs, normalize, .{});
}

/// `variance`, both passes summed as `options` says: floats of `R` (or
/// `Complex(R)`) through `kernels.stats.deviationSums`, other elements in
/// order here, in `f64` or wider.
pub fn varianceWith(io: Io, comptime R: type, xs: anytype, normalize: Normalize, options: MomentOptions) MomentError!R {
    const E = Elem(@TypeOf(xs));
    comptime checkFloat(R);
    const W = Wide(R);
    const M = if (comptime isComplex(E)) Complex(R) else R;
    const m = try meanWith(io, M, xs, options);
    var sum_sq: W = undefined;
    var square_of_sum: W = undefined;
    if (comptime E == M and kernelFloat(R)) {
        const sl: []const E = xs;
        const d = try kernels_stats.deviationSums(io, E, options.summation, m, .of(sl));
        sum_sq = d.sum_sq;
        square_of_sum = if (comptime isComplex(E)) @as(W, d.sum.re) * d.sum.re + @as(W, d.sum.im) * d.sum.im else @as(W, d.sum) * d.sum;
    } else {
        const s = options.summation;
        var sq: Acc(W) = .{};
        var re: Acc(W) = .{};
        var im: Acc(W) = .{};
        for (xs) |x| {
            if (comptime isComplex(E)) {
                const dr = to(W, x.re) - @as(W, m.re);
                const di = to(W, x.im) - @as(W, m.im);
                re.add(dr, s);
                im.add(di, s);
                sq.add(dr * dr + di * di, s);
            } else {
                const d = to(W, x) - @as(W, m);
                re.add(d, s);
                sq.add(d * d, s);
            }
        }
        sum_sq = sq.value();
        square_of_sum = re.value() * re.value() + im.value() * im.value();
    }
    const n = xs.len;
    // a sample of one: 0, or NaN from a NaN or an infinity
    const div = if (normalize == .sample and n > 1) n - 1 else n;
    return @floatCast((sum_sq - square_of_sum / @as(W, @floatFromInt(n))) / @as(W, @floatFromInt(div)));
}

pub fn stdDev(io: Io, comptime R: type, xs: anytype, normalize: Normalize) MomentError!R {
    return @sqrt(try variance(io, R, xs, normalize));
}

/// `stdDev`, summed as `options` says.
pub fn stdDevWith(io: Io, comptime R: type, xs: anytype, normalize: Normalize, options: MomentOptions) MomentError!R {
    return @sqrt(try varianceWith(io, R, xs, normalize, options));
}

/// The float types the kernel layer's sums take.
fn kernelFloat(comptime F: type) bool {
    return F == f32 or F == f64 or F == blas_types.XDouble;
}

/// The sum of real `xs` as `R`: the BLAS's where the elements are `R` and
/// it has a `sum` for them, else in order here.
fn total(io: Io, comptime R: type, xs: anytype, s: Summation) error{Canceled}!R {
    const E = Elem(@TypeOf(xs));
    if (comptime E == R and kernelFloat(R)) {
        const sl: []const R = xs;
        return blas.sumWith(io, R, .of(sl), s);
    }
    const W = Wide(R);
    var acc: Acc(W) = .{};
    for (xs) |x| acc.add(to(W, x), s);
    return @floatCast(acc.value());
}

/// A running sum, plain or with Neumaier's compensation.
fn Acc(comptime W: type) type {
    return struct {
        sum: W = 0,
        c: W = 0,

        fn add(a: *@This(), x: W, s: Summation) void {
            const t = a.sum + x;
            if (s == .compensated) {
                if (@abs(a.sum) >= @abs(x)) a.c += (a.sum - t) + x else a.c += (x - t) + a.sum;
            }
            a.sum = t;
        }

        fn value(a: @This()) W {
            return a.sum + a.c;
        }
    };
}
