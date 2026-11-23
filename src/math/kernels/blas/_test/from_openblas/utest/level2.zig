//! `test_trmv.c` and `test_gemv.c`: `trmv` at unit stride, and `gemv` through NaN and infinity.

const std = @import("std");
const ob = @import("../__root__.zig");
const l2 = @import("../../../blas2/core.zig");
const ut = @import("../utest.zig");

const C32 = std.math.Complex(f32);
const C64 = std.math.Complex(f64);
const Want = ut.Want;
const is = ut.is;
const nears = ut.nears;
const single_eps = ut.single_eps;

// ---- test_trmv.c --------------------------------------------------------------------

pub fn strmvLower() !void {
    // row-major lower, not transposed: column-major upper, transposed
    var lower = [9]f32{ 2.0, 0.0, 0.0, 3.0, 4.0, 0.0, 5.0, 6.0, 7.0 };
    var v = [3]f32{ 11.0, 13.0, 17.0 };
    try l2.trmv(std.testing.allocator, ob.io, f32, .u, .t, .n, 3, &lower, 3, &v, 1);
    try nears(f32, &.{ 22.0, 85.0, 252.0 }, &v, single_eps, "strmv");
}

test "utest: strmv, a small lower triangle at unit stride" {
    try strmvLower();
}

// ---- test_gemv.c --------------------------------------------------------------------

pub fn gemvNanInf() !void {
    const al = std.testing.allocator;
    const n = 17;
    inline for (.{ f32, f64 }) |F| for ([_]isize{ 1, 2 }) |inc| {
        var a: [n * n]F = @splat(0);
        var x: [n]F = @splat(0);
        var y: [2 * n - 1]F = @splat(0);
        const ny: usize = if (inc == 1) n else 2 * n - 1;
        for (0..n) |i| y[i * @as(usize, @intCast(inc))] = if (i % 2 == 0) std.math.nan(F) else std.math.inf(F);
        try l2.gemv(al, ob.io, F, .n, n, n, 0, &a, n, &x, 1, 0, &y, inc);
        for (y[0..ny]) |v| try ob.expect(v == 0, "sgemv 0_nan_inf");
    };
}

test "utest: gemv, zero alpha and beta overwrite NaN and infinity" {
    try gemvNanInf();
}

pub fn cgemvNanInf() !void {
    const al = std.testing.allocator;
    const n = 17;
    inline for (.{ C32, C64 }) |C| {
        const F = ob.Float(C);
        const Fn = std.math.nan(F);
        const Fi = std.math.inf(F);
        // beta; and y's even and odd elements' parts, before and after
        const Case = struct { name: []const u8, beta: C, even: C, odd: C, want_even: [2]Want, want_odd: [2]Want };
        const cases = [_]Case{
            .{ .name = "0_nan_inf", .beta = .init(0, 0), .even = .init(Fn, Fn), .odd = .init(Fi, Fi), .want_even = .{ .zero, .zero }, .want_odd = .{ .zero, .zero } },
            .{ .name = "0_2_nan_1_inf_1", .beta = .init(0, 2), .even = .init(Fn, 1), .odd = .init(Fi, 1), .want_even = .{ .nan, .nan }, .want_odd = .{ .nan, .inf } },
            .{ .name = "2_0_nan_1_inf_1", .beta = .init(2, 0), .even = .init(Fn, 1), .odd = .init(Fi, 1), .want_even = .{ .nan, .nan }, .want_odd = .{ .inf, .nan } },
        };
        for (cases) |c| for ([_]isize{ 1, 2 }) |inc| {
            var a: [n * n]C = @splat(.init(0, 0));
            var x: [n]C = @splat(.init(0, 0));
            var y: [2 * n]C = @splat(.init(0, 0));
            const s: usize = @intCast(inc);
            for (0..n) |i| y[i * s] = if (i % 2 == 0) c.even else c.odd;
            try l2.gemv(al, ob.io, C, .n, n, n, .init(0, 0), &a, n, &x, 1, c.beta, &y, inc);
            for (y[0 .. n * s], 0..) |v, k| {
                if (k % s != 0) {
                    try ob.expect(v.re == 0 and v.im == 0, c.name);
                    continue;
                }
                const w = if ((k / s) % 2 == 0) c.want_even else c.want_odd;
                try ob.expect(is(w[0], v.re) and is(w[1], v.im), c.name);
            }
        };
    }
}

test "utest: cgemv, zgemv, NaN and infinity through beta 0, 2i and 2" {
    try cgemvNanInf();
}
