const std = @import("std");
const ReduceOp = std.builtin.ReduceOp;

fn Elem(comptime ArrType: type) type {
    return switch (@typeInfo(ArrType)) {
        .array => |info| info.child,
        .pointer => |info| info.child,
        else => @compileError("expects a fixed-size array or a slice, found " ++ @typeName(ArrType)),
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
