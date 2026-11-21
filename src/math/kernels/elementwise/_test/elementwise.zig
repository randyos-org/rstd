//! The elementwise kernels against Zig's scalar builtins, bit for bit, over
//! `@Vector` blocks and strided views, in place, and at any worker count.

const std = @import("std");
const types = @import("../../blas/types.zig");
const strided = @import("../../blas/strided.zig");
const threading = @import("../../blas/threading.zig");
const check = @import("../../blas/_test/check.zig");
const test_io = @import("../../blas/_test/io.zig");
const scalar = @import("../../../scalar.zig");
const ew = @import("../__root__.zig");

const testing = std.testing;
const Io = std.Io;
const ConstStrided = strided.ConstStrided;
const Strided = strided.Strided;
const Complex32 = types.Complex32;
const Complex64 = types.Complex64;

const Unary = enum { sqrt, exp, exp2, log, log2, log10, sin, cos, tan, abs, floor, ceil, trunc, round };
const Binary = enum { add, sub, mul, div };

fn want1(comptime op: Unary, v: anytype) @TypeOf(v) {
    return switch (op) {
        .sqrt => @sqrt(v),
        .exp => @exp(v),
        .exp2 => @exp2(v),
        .log => @log(v),
        .log2 => @log2(v),
        .log10 => @log10(v),
        .sin => @sin(v),
        .cos => @cos(v),
        .tan => @tan(v),
        .abs => @abs(v),
        .floor => scalar.floor(v),
        .ceil => scalar.ceil(v),
        .trunc => @trunc(v),
        .round => scalar.round(v),
    };
}

fn want2(comptime T: type, comptime op: Binary, a: T, b: T) T {
    if (comptime types.isComplex(T)) return switch (op) {
        .add => a.add(b),
        .sub => a.sub(b),
        .mul => a.mul(b),
        .div => a.div(b),
    };
    return switch (op) {
        .add => a + b,
        .sub => a - b,
        .mul => a * b,
        .div => a / b,
    };
}

fn bits(comptime T: type, v: T) u128 {
    if (comptime types.isComplex(T)) return (@as(u128, bits(@TypeOf(v.re), v.re)) << 64) | bits(@TypeOf(v.im), v.im);
    return @as(std.meta.Int(.unsigned, @bitSizeOf(T)), @bitCast(v));
}

/// Values in the functions' domains: positive, from 1/16 to about 2.
fn fillPositive(comptime T: type, buf: []T, seed: u64) void {
    check.fill(T, buf, seed);
    for (buf) |*v| v.* = @abs(v.*) + 0.0625;
}

/// The span a view of `n` elements at `stride` needs.
fn span(n: usize, stride: isize) usize {
    return if (n == 0) 0 else (n - 1) * @abs(stride) + 1;
}

const sizes = [_]usize{ 0, 1, 3, 15, 16, 17, 64, 1000, 3001 };
const strides = [_][2]isize{ .{ 1, 1 }, .{ 2, 1 }, .{ 1, -3 }, .{ -2, -1 } };

test "one input: each function's builtin, bit for bit, at any stride" {
    const al = testing.allocator;
    inline for (.{ f32, f64 }) |T| {
        const bx = try al.alloc(T, 3 * 3001);
        defer al.free(bx);
        const by = try al.alloc(T, 3 * 3001);
        defer al.free(by);
        fillPositive(T, bx, 1);
        inline for (comptime std.enums.values(Unary)) |op| {
            for (sizes) |n| for (strides) |s| {
                const x = try ConstStrided(T).init(bx[0..span(n, s[0])], n, s[0]);
                const y = try Strided(T).init(by[0..span(n, s[1])], n, s[1]);
                try @field(ew, @tagName(op))(check.io, T, x, y);
                for (0..n) |i| try testing.expectEqual(bits(T, want1(op, x.get(i))), bits(T, y.get(i)));
            };
        }
    }
}

test "two inputs: exact arithmetic, real and complex, at any stride" {
    const al = testing.allocator;
    inline for (.{ f32, f64, Complex32, Complex64 }) |T| {
        const bx = try al.alloc(T, 3 * 3001);
        defer al.free(bx);
        const by = try al.alloc(T, 3 * 3001);
        defer al.free(by);
        const bz = try al.alloc(T, 3 * 3001);
        defer al.free(bz);
        check.fill(T, bx, 2);
        check.fill(T, by, 3);
        inline for (comptime std.enums.values(Binary)) |op| {
            for (sizes) |n| for (strides) |s| {
                const x = try ConstStrided(T).init(bx[0..span(n, s[0])], n, s[0]);
                const y = try ConstStrided(T).init(by[0..span(n, s[1])], n, s[1]);
                const z = try Strided(T).init(bz[0..span(n, s[1])], n, s[1]);
                try @field(ew, @tagName(op))(check.io, T, x, y, z);
                for (0..n) |i| try testing.expectEqual(bits(T, want2(T, op, x.get(i), y.get(i))), bits(T, z.get(i)));
            };
        }
    }
}

test "abs of complex: the modulus, without overflow" {
    const big = std.math.floatMax(f64) / 2;
    const x = [_]Complex64{ .init(3, 4), .init(big, big), .init(-0.0, 0), .init(std.math.inf(f64), std.math.nan(f64)) };
    var y: [4]f64 = undefined;
    try ew.abs(check.io, Complex64, .of(&x), .of(&y));
    try testing.expectEqual(5, y[0]);
    try testing.expectEqual(std.math.hypot(big, big), y[1]);
    try testing.expect(std.math.isFinite(y[1]));
    try testing.expectEqual(0, y[2]);
    try testing.expect(std.math.isInf(y[3]));
}

test "in place, and a length mismatch" {
    var x = [_]f64{ 0.25, 4, 9, 16, 25, 36, 49, 64, 81, 100 };
    try ew.sqrt(check.io, f64, .of(&x), .of(&x));
    try testing.expectEqualSlices(f64, &.{ 0.5, 2, 3, 4, 5, 6, 7, 8, 9, 10 }, &x);
    var a = [_]f32{ 1, 2, 3 };
    try ew.mul(check.io, f32, .of(&a), .of(&a), .of(&a));
    try testing.expectEqualSlices(f32, &.{ 1, 4, 9 }, &a);
    var y: [2]f32 = undefined;
    try testing.expectError(error.LengthMismatch, ew.exp(check.io, f32, .of(&a), .of(&y)));
    try testing.expectError(error.LengthMismatch, ew.add(check.io, f32, .of(&a), .of(a[0..2]), .of(&a)));
    try testing.expectError(error.LengthMismatch, ew.add(check.io, f32, .of(&a), .of(&a), .of(&y)));
}

test "round: ties away from zero; floor, ceil and trunc of negatives" {
    const x = [_]f64{ 0.5, 1.5, 2.5, -0.5, -2.5, -1.25 };
    var y: [6]f64 = undefined;
    try ew.round(check.io, f64, .of(&x), .of(&y));
    try testing.expectEqualSlices(f64, &.{ 1, 2, 3, -1, -3, -1 }, &y);
    try ew.floor(check.io, f64, .of(&x), .of(&y));
    try testing.expectEqualSlices(f64, &.{ 0, 1, 2, -1, -3, -2 }, &y);
    try ew.ceil(check.io, f64, .of(&x), .of(&y));
    try testing.expectEqualSlices(f64, &.{ 1, 2, 3, -0.0, -2, -1 }, &y);
    try ew.trunc(check.io, f64, .of(&x), .of(&y));
    try testing.expectEqualSlices(f64, &.{ 0, 1, 2, -0.0, -2, -1 }, &y);
}

test "elementwise: one worker's bits at any worker count" {
    if (comptime !test_io.hosted) return error.SkipZigTest;
    const al = testing.allocator;
    var runs: test_io.WorkerRuns = .init(al);
    defer runs.deinit();
    const n = (1 << 15) + 37;
    const x = try al.alloc(f64, n);
    defer al.free(x);
    fillPositive(f64, x, 9);
    var outs: [test_io.worker_runs.len][]f64 = undefined;
    for (&outs) |*o| o.* = try al.alloc(f64, 2 * n);
    defer for (outs) |o| al.free(o);
    for (outs, 0..) |o, c| {
        const io = runs.io(c);
        // the same chunks at every count: only the workers differ
        threading.test_serial_max = 4096;
        try ew.log(io, f64, .of(x), .of(o[0..n]));
        try ew.div(io, f64, .of(x), .of(o[0..n]), .of(o[n..]));
    }
    try test_io.expectSameBytes(outs);
}
