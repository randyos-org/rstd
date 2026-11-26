//! `blas1/interface.zig`'s tests: the checks, strides and quick returns.

const std = @import("std");
const types = @import("../../types.zig");
const util = @import("../../util.zig");
const bf = @import("../../ext/bf16/__root__.zig");
const interface = @import("../../blas1/interface.zig");

const BFloat16 = types.BFloat16;
const XDouble = types.XDouble;
const isComplex = types.isComplex;
const Float = types.Float;
const ArgError = interface.ArgError;
const copy = interface.copy;
const swap = interface.swap;
const scal = interface.scal;
const scalReal = interface.scalReal;
const dot = interface.dot;
const dotc = interface.dotc;
const dsdot = interface.dsdot;
const sdsdot = interface.sdsdot;
const nrm2 = interface.nrm2;
const asum = interface.asum;
const sum = interface.sum;
const index = interface.index;
const extreme = interface.extreme;
const rot = interface.rot;
const rotm = interface.rotm;
const axpy = interface.axpy;
const axpyc = interface.axpyc;
const axpby = interface.axpby;
const sbdot = interface.sbdot;
const bscal = interface.bscal;
const sbtobf16 = interface.sbtobf16;
const bf16to = interface.bf16to;
const cabs1 = interface.cabs1;

const check = @import("../check.zig");
const reference = @import("../reference.zig");
const blas2_interface = @import("../../blas2/interface.zig");
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
            try check.expectClose(Float(T), 0, try extreme(check.io, T, .largest_magnitude, n, &x, 1), 0);
            try sameInt(0, try index(check.io, T, .largest_magnitude, n, &x, 1));
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
            try testing.expectError(error.InvalidIncx, index(check.io, T, .largest_magnitude, 2, &x, b));
            try testing.expectError(error.InvalidIncx, extreme(check.io, T, .largest_magnitude, 2, &x, b));
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
    var h = [_]BFloat16{ b16(1), b16(2) };
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
    try sameInt(2, try index(check.io, f64, .largest_magnitude, 5, &x, 1));
    try sameInt(4, try index(check.io, f64, .smallest_magnitude, 5, &x, 1));
    try sameInt(3, try index(check.io, f64, .largest, 5, &x, 1));
    try sameInt(2, try index(check.io, f64, .smallest, 5, &x, 1));
    try sameInt(0, try index(check.io, f64, .largest_magnitude, 5, &x, 0));
    try check.expectClose(f64, 7, try extreme(check.io, f64, .largest_magnitude, 5, &x, 1), 0);
    try check.expectClose(f64, 1, try extreme(check.io, f64, .smallest_magnitude, 5, &x, 1), 0);
    try check.expectClose(f64, 7, try extreme(check.io, f64, .largest, 5, &x, 1), 0);
    try check.expectClose(f64, -7, try extreme(check.io, f64, .smallest, 5, &x, 1), 0);
    try check.expectClose(f64, -7, try extreme(check.io, f64, .largest, 5, x[1..], 0), 0);
    try check.expectClose(f64, 7, try extreme(check.io, f64, .largest_magnitude, 5, x[1..], 0), 0);
    var z = [_]Complex32{ .init(1, 1), .init(-3, 2), .init(0, 4) };
    try sameInt(2, try index(check.io, Complex32, .largest_magnitude, 3, &z, 1));
    try check.expectClose(f32, 5, try extreme(check.io, Complex32, .largest_magnitude, 3, &z, 1), 0);
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

fn b16(v: f32) BFloat16 {
    return bf.toBf16(v);
}

test "interface: sbdot, the conversions and bscal do nothing for n <= 0" {
    var x = [_]BFloat16{ b16(1), b16(2) };
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
    var h: [n * 3]BFloat16 = @splat(b16(0));
    // stride -3 puts element 0 last: the last element of the buffer region
    try sbtobf16(check.io, f32, n, &f, 2, &h, -3);
    for (0..n) |i| {
        try std.testing.expectEqual(@backingInt(b16(f[i * 2])), @backingInt(h[(n - 1 - i) * 3]));
    }
    var back: [n]f64 = undefined;
    try bf16to(check.io, f64, n, &h, -3, &back, 1);
    for (0..n) |i| try std.testing.expectEqual(@as(f64, bf.fromBf16(b16(f[i * 2]))), back[i]);
    // a reversed read of the same storage dots against itself forwards
    var fwd: [n]BFloat16 = undefined;
    for (0..n) |i| fwd[i] = h[(n - 1 - i) * 3];
    try std.testing.expectEqual(try sbdot(check.io, n, &fwd, 1, &fwd, 1), try sbdot(check.io, n, &h, -3, &h, -3));
    var want: f64 = 0;
    for (0..n) |i| want += @as(f64, bf.fromBf16(fwd[i])) * bf.fromBf16(fwd[i]);
    try check.expectClose(f32, @floatCast(want), try sbdot(check.io, n, &fwd, 1, &fwd, 1), check.tolerance(f32, 4 * n));
}

test "interface: bscal skips alpha 1 and a non-positive stride, and scales the rest" {
    var x = [_]BFloat16{ b16(1.5), b16(-2), b16(4) };
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
