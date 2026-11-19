//! The moments of a slice: mean, variance and standard deviation, through
//! the BLAS where it has the sum.

const std = @import("std");
const blas = @import("../linalg/blas/__root__.zig");
const blas_types = @import("../linalg/blas/types.zig");
const util = @import("util.zig");

const Io = std.Io;
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

/// How the sums are taken: as the BLAS's kernels do, or compensated
/// (Neumaier's, `blas.Summation`).
pub const Summation = blas.Summation;

/// The mean of `xs` (a slice or array of integers, floats or `Complex`), as
/// `R`: a float for real elements, a `Complex` for complex ones. Floats of
/// `R`'s own type sum through the BLAS `sum`.
pub fn mean(io: Io, comptime R: type, xs: anytype) MomentError!R {
    return meanWith(io, R, xs, .plain);
}

/// `mean`, summed as `s` says.
pub fn meanWith(io: Io, comptime R: type, xs: anytype, s: Summation) MomentError!R {
    const E = Elem(@TypeOf(xs));
    comptime checkPair(R, E);
    if (xs.len == 0) return error.EmptyInput;
    const n: Real(R) = @floatFromInt(xs.len);
    if (comptime isComplex(E)) {
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

/// Two passes: the mean, then the squared deviations (`|x - mean|²` for
/// complex `xs`). `R` is a float.
pub fn variance(io: Io, comptime R: type, xs: anytype, normalize: Normalize) MomentError!R {
    return varianceWith(io, R, xs, normalize, .plain);
}

/// `variance`, both passes summed as `s` says; the deviations' squares
/// summed in order, in `f64` or wider.
pub fn varianceWith(io: Io, comptime R: type, xs: anytype, normalize: Normalize, s: Summation) MomentError!R {
    const E = Elem(@TypeOf(xs));
    comptime checkFloat(R);
    const W = Wide(R);
    var acc: Acc(W) = .{};
    if (comptime isComplex(E)) {
        const m = try meanWith(io, std.math.Complex(R), xs, s);
        for (xs) |x| {
            const dr = to(W, x.re) - @as(W, m.re);
            const di = to(W, x.im) - @as(W, m.im);
            acc.add(dr * dr + di * di, s);
        }
    } else {
        const m = try meanWith(io, R, xs, s);
        for (xs) |x| {
            const d = to(W, x) - @as(W, m);
            acc.add(d * d, s);
        }
    }
    const n = xs.len;
    // a sample of one: 0, or NaN from a NaN or an infinity
    const d = if (normalize == .sample and n > 1) n - 1 else n;
    return @floatCast(acc.value() / @as(W, @floatFromInt(d)));
}

pub fn stdDev(io: Io, comptime R: type, xs: anytype, normalize: Normalize) MomentError!R {
    return @sqrt(try variance(io, R, xs, normalize));
}

/// `stdDev`, summed as `s` says.
pub fn stdDevWith(io: Io, comptime R: type, xs: anytype, normalize: Normalize, s: Summation) MomentError!R {
    return @sqrt(try varianceWith(io, R, xs, normalize, s));
}

/// The sum of real `xs` as `R`: the BLAS's where the elements are `R` and
/// it has a `sum` for them, else in order here.
fn total(io: Io, comptime R: type, xs: anytype, s: Summation) error{Canceled}!R {
    const E = Elem(@TypeOf(xs));
    if (comptime E == R and (R == f32 or R == f64 or R == blas_types.XDouble)) {
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
