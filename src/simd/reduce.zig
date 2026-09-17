const std = @import("std");
const zigstd = @import("../zigstd.zig");
const ReduceOp = std.builtin.ReduceOp;

fn Elem(comptime ArrType: type) type {
    return switch (@typeInfo(ArrType)) {
        .array => |info| info.child,
        .vector => |info| info.child,
        .pointer => |info| info.child,
        else => @compileError("expects a fixed-size array, a vector, or a slice, found " ++ @typeName(ArrType)),
    };
}

/// Combines two scalars, or two equal-shape vectors, the same way `@reduce(op, ...)`
/// would combine two lanes -- this is the piece `@reduce` doesn't expose on its own,
/// needed to fold separate SIMD chunks (and the scalar tail) back together.
fn combine(comptime op: ReduceOp, a: anytype, b: @TypeOf(a)) @TypeOf(a) {
    return switch (op) {
        .Max => @max(a, b),
        .Min => @min(a, b),
        .Add => a + b,
        .Mul => a * b,
        .And => a & b,
        .Or => a | b,
        .Xor => a ^ b,
    };
}

pub fn reduce(comptime op: ReduceOp, arr: anytype) Elem(@TypeOf(arr)) {
    return switch (@typeInfo(@TypeOf(arr))) {
        .array => |info| @reduce(op, @as(@Vector(info.len, info.child), arr)),
        .vector => @reduce(op, arr),
        .pointer => |info| sliceReduce(op, info.child, arr),
        else => unreachable,
    };
}

fn scalarReduce(comptime op: ReduceOp, comptime T: type, arr: []const T) T {
    var x = arr[0];
    for (arr[1..]) |v| x = combine(op, x, v);
    return x;
}

fn sliceReduce(comptime op: ReduceOp, comptime T: type, arr: []const T) T {
    const chunk_len = comptime std.simd.suggestVectorLength(T) orelse 1;

    if (chunk_len < 2 or arr.len < chunk_len) return scalarReduce(op, T, arr);

    var x: @Vector(chunk_len, T) = arr[0..chunk_len].*;
    var i: usize = chunk_len;
    while (i + chunk_len <= arr.len) : (i += chunk_len) {
        const next: @Vector(chunk_len, T) = arr[i..][0..chunk_len].*;
        x = combine(op, x, next);
    }

    const vec_result = @reduce(op, x);
    if (i == arr.len) return vec_result;
    return combine(op, vec_result, scalarReduce(op, T, arr[i..]));
}

pub fn arrayMax(arr: anytype) Elem(@TypeOf(arr)) {
    return reduce(.Max, arr);
}

pub fn arrayMin(arr: anytype) Elem(@TypeOf(arr)) {
    return reduce(.Min, arr);
}

pub fn arraySum(arr: anytype) Elem(@TypeOf(arr)) {
    return reduce(.Add, arr);
}

pub fn arrayProduct(arr: anytype) Elem(@TypeOf(arr)) {
    return reduce(.Mul, arr);
}

pub fn arrayAnd(arr: anytype) Elem(@TypeOf(arr)) {
    return reduce(.And, arr);
}

pub fn arrayOr(arr: anytype) Elem(@TypeOf(arr)) {
    return reduce(.Or, arr);
}

pub fn arrayXor(arr: anytype) Elem(@TypeOf(arr)) {
    return reduce(.Xor, arr);
}

/// A tuple type of `n` copies of `ArgT` - the args-tuple `mapFn` gets called
/// with, built by hand (rather than `std.meta.ArgsTuple(@TypeOf(mapFn))`)
/// since `mapFn` is expected to have `anytype` parameters, which
/// `ArgsTuple` explicitly rejects.
fn ArgsOf(comptime n: usize, comptime ArgT: type) type {
    var types: [n]type = undefined;
    for (&types) |*t| t.* = ArgT;
    return @Tuple(&types);
}

fn mapAt(comptime T: type, comptime mapFn: anytype, slices: anytype, i: usize) T {
    var args: ArgsOf(slices.len, T) = undefined;
    inline for (slices, 0..) |s, idx| args[idx] = s[i];
    return @call(.auto, mapFn, args);
}

fn mapScalarReduce(comptime op: ReduceOp, comptime T: type, comptime mapFn: anytype, slices: anytype, start: usize, end: usize) T {
    var x: T = mapAt(T, mapFn, slices, start);
    var i = start + 1;
    while (i < end) : (i += 1) x = combine(op, x, mapAt(T, mapFn, slices, i));
    return x;
}

fn mapVecAt(comptime T: type, comptime N: usize, comptime mapFn: anytype, slices: anytype, offset: usize) @Vector(N, T) {
    var args: ArgsOf(slices.len, @Vector(N, T)) = undefined;
    inline for (slices, 0..) |s, idx| args[idx] = @as(@Vector(N, T), s[offset..][0..N].*);
    return @call(.auto, mapFn, args);
}

/// Generalizes `sliceReduce` to slices that need a per-index transform before
/// folding, e.g. `sum(|re[i]| + |im[i]|)` or `sum(xr[i]*yr[i] - xi[i]*yi[i])`
/// - the shape of `blas1.wasum`/`wdot`'s hand-rolled SIMD loops. `slices` is
/// a tuple of equal-length `[]const T` slices (e.g. `.{ xr, xi }`); `mapFn`
/// combines one element from each into a single `T`.
///
/// `mapFn` must be written with `anytype` parameters, e.g.:
/// ```
/// fn cabs1(re: anytype, im: anytype) @TypeOf(re) {
///     return @abs(re) + @abs(im);
/// }
/// ```
/// so the exact same body type-checks both against a scalar `T` (the tail)
/// and a `@Vector(N, T)` chunk (the SIMD body) - every operator a map like
/// this would use (`+`, `-`, `*`, `@abs`, ...) is already defined
/// element-wise for both, so there is no separate vector/scalar copy to keep
/// in sync.
pub fn mapReduce(
    comptime op: ReduceOp,
    comptime T: type,
    comptime mapFn: anytype,
    slices: anytype,
) T {
    const n = slices[0].len;
    inline for (slices) |s| std.debug.assert(s.len == n);

    const chunk_len = comptime std.simd.suggestVectorLength(T) orelse 1;
    if (chunk_len < 2 or n < chunk_len) return mapScalarReduce(op, T, mapFn, slices, 0, n);

    var x: @Vector(chunk_len, T) = mapVecAt(T, chunk_len, mapFn, slices, 0);
    var i: usize = chunk_len;
    while (i + chunk_len <= n) : (i += chunk_len) {
        x = combine(op, x, mapVecAt(T, chunk_len, mapFn, slices, i));
    }

    const vec_result = @reduce(op, x);
    if (i == n) return vec_result;
    return combine(op, vec_result, mapScalarReduce(op, T, mapFn, slices, i, n));
}

// ---

fn cabs1(re: anytype, im: anytype) @TypeOf(re) {
    return @abs(re) + @abs(im);
}

fn unconjRe(xr: anytype, xi: @TypeOf(xr), yr: @TypeOf(xr), yi: @TypeOf(xr)) @TypeOf(xr) {
    return xr * yr - xi * yi;
}

test "mapReduce(.Add, cabs1) matches a plain scalar loop, scalar-only path" {
    const xr = [_]f64{ 1.0, -2.0, 3.0 };
    const xi = [_]f64{ -0.5, 1.5, 0.0 };
    const got = mapReduce(.Add, f64, cabs1, .{ &xr, &xi });

    var want: f64 = 0.0;
    for (xr, xi) |re, im| want += cabs1(re, im);
    try zigstd.testExpectEqual(want, got);
}

test "mapReduce(.Add, cabs1) matches a plain scalar loop, crossing a SIMD-width boundary" {
    const chunk_len = comptime std.simd.suggestVectorLength(f64) orelse 1;
    const n = chunk_len * 2 + 3;

    var xr: [n]f64 = undefined;
    var xi: [n]f64 = undefined;
    for (0..n) |k| {
        xr[k] = @floatFromInt(k + 1);
        xi[k] = -@as(f64, @floatFromInt(k));
    }

    const got = mapReduce(.Add, f64, cabs1, .{ &xr, &xi });

    var want: f64 = 0.0;
    for (xr, xi) |re, im| want += cabs1(re, im);
    try zigstd.testExpectEqual(want, got);
}

test "mapReduce(.Add, unconjRe) reduces four parallel slices" {
    const chunk_len = comptime std.simd.suggestVectorLength(f64) orelse 1;
    const n = chunk_len + 1;

    var xr: [n]f64 = undefined;
    var xi: [n]f64 = undefined;
    var yr: [n]f64 = undefined;
    var yi: [n]f64 = undefined;
    for (0..n) |k| {
        xr[k] = @floatFromInt(k + 1);
        xi[k] = @floatFromInt(2 * k);
        yr[k] = 1.0;
        yi[k] = -1.0;
    }

    const got = mapReduce(.Add, f64, unconjRe, .{ &xr, &xi, &yr, &yi });

    var want: f64 = 0.0;
    for (xr, xi, yr, yi) |a, b, c, d| want += unconjRe(a, b, c, d);
    try zigstd.testExpectEqual(want, got);
}
