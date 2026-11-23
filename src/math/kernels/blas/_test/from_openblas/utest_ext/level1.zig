//! The level-1 extensions: `axpby`, `axpyc`, complex `rot`, `scal` and `rotg`, and `rotmg`.

const std = @import("std");
const ob = @import("../__root__.zig");
const blas = @import("../../../__root__.zig");
const l1 = @import("../../../blas1/interface.zig");
const util = @import("util.zig");

const testing = std.testing;
const Float = util.Float;
const isComplex = util.isComplex;
const C32 = util.C32;
const C64 = util.C64;
const single_eps = util.single_eps;
const double_eps = util.double_eps;
const tol = util.tol;
const Case = util.Case;
const val = util.val;
const fill = util.fill;
const span = util.span;
const diffNorm = util.diffNorm;
const expectNorm = util.expectNorm;
const conjugate = util.conjugate;

pub fn axpby() !void {
    var prng = std.Random.DefaultPrng.init(0xA8);
    const al = testing.allocator;
    inline for (.{ f32, f64, C32, C64 }) |T| {
        const cases: []const Case = if (comptime isComplex(T)) &axpby_complex else &axpby_real;
        for (cases) |c| {
            const x = try al.alloc(T, span(c.n, c.incx));
            defer al.free(x);
            const y = try al.alloc(T, span(c.n, c.incy));
            defer al.free(y);
            fill(T, x, prng.random());
            fill(T, y, prng.random());
            const yv = try al.dupe(T, y);
            defer al.free(yv);
            // scal takes only positive strides
            try l1.scal(ob.io, T, c.n, val(T, c.beta), yv.ptr, @intCast(@abs(c.incy)));
            try l1.axpy(ob.io, T, c.n, val(T, c.alpha), x.ptr, c.incx, yv.ptr, c.incy);
            try l1.axpby(ob.io, T, c.n, val(T, c.alpha), x.ptr, c.incx, val(T, c.beta), y.ptr, c.incy);
            try expectNorm(c.name, T, diffNorm(T, c.n, y, yv, c.incy), tol(T, single_eps, double_eps));
        }
    }
}

test "utest_ext: axpby against scal and axpy" {
    try axpby();
}

pub fn axpyc() !void {
    var prng = std.Random.DefaultPrng.init(0xA9);
    const al = testing.allocator;
    inline for (.{ C32, C64 }) |T| for (axpyc_cases) |c| {
        const x = try al.alloc(T, span(c.n, c.incx));
        defer al.free(x);
        const y = try al.alloc(T, span(c.n, c.incy));
        defer al.free(y);
        fill(T, x, prng.random());
        fill(T, y, prng.random());
        const xv = try al.dupe(T, x);
        defer al.free(xv);
        const yv = try al.dupe(T, y);
        defer al.free(yv);
        conjugate(T, xv);
        try l1.axpy(ob.io, T, c.n, val(T, c.alpha), xv.ptr, c.incx, yv.ptr, c.incy);
        try l1.axpyc(ob.io, T, c.n, val(T, c.alpha), x.ptr, c.incx, y.ptr, c.incy);
        try expectNorm(c.name, T, diffNorm(T, c.n, y, yv, c.incy), tol(T, single_eps, double_eps));
    };
}

test "utest_ext: axpyc against conjugation and axpy" {
    try axpyc();
}

pub fn complexRot() !void {
    var prng = std.Random.DefaultPrng.init(0xAA);
    const al = testing.allocator;
    inline for (.{ C32, C64 }) |T| for (rot_cases) |c| {
        const x = try al.alloc(T, span(c.n, c.inc_x));
        defer al.free(x);
        const y = try al.alloc(T, span(c.n, c.inc_y));
        defer al.free(y);
        fill(T, x, prng.random());
        fill(T, y, prng.random());
        const xv = try al.dupe(T, x);
        defer al.free(xv);
        const yv = try al.dupe(T, y);
        defer al.free(yv);
        // c*x + s*y, then -conj(s)*x + c*y
        try l1.axpby(ob.io, T, c.n, val(T, c.s), y.ptr, c.inc_y, val(T, c.c), xv.ptr, c.inc_x);
        try l1.axpby(ob.io, T, c.n, val(T, .{ -c.s[0], c.s[1] }), x.ptr, c.inc_x, val(T, c.c), yv.ptr, c.inc_y);
        try l1.rot(ob.io, T, c.n, x.ptr, c.inc_x, y.ptr, c.inc_y, @floatCast(c.c[0]), @floatCast(c.s[0]));
        const norm = (diffNorm(T, c.n, x, xv, c.inc_x) + diffNorm(T, c.n, y, yv, c.inc_y)) / 2;
        try expectNorm(c.name, T, norm, tol(T, single_eps, double_eps));
    };
}

test "utest_ext: complex rot against axpby" {
    try complexRot();
}

pub fn complexScal() !void {
    var prng = std.Random.DefaultPrng.init(0xAB);
    const al = testing.allocator;
    inline for (.{ C32, C64 }) |T| for (scal_cases) |c| {
        const x = try al.alloc(T, span(c.N, c.inc));
        defer al.free(x);
        fill(T, x, prng.random());
        const xv = try al.dupe(T, x);
        defer al.free(xv);
        const a = val(T, c.alpha);
        for (0..@intCast(c.N)) |i| xv[i * @as(usize, @intCast(c.inc))] = a.mul(xv[i * @as(usize, @intCast(c.inc))]);
        try l1.scal(ob.io, T, c.N, a, x.ptr, c.inc);
        try expectNorm(c.name, T, diffNorm(T, c.N, x, xv, c.inc), tol(T, single_eps, double_eps));
    };
}

test "utest_ext: complex scal against the plain product" {
    try complexScal();
}

pub fn complexRotg() !void {
    const Rotg = struct { name: []const u8, a: [2]f64, b: [2]f64, c: f64, s: [2]f64, r: [2]f64 };
    // the double cases' expected values; the single ones' are the same to 5 places
    const cases = [_]Rotg{
        .{ .name = "zero_a", .a = .{ 0, 0 }, .b = .{ 1, 1 }, .c = 0.0, .s = .{ 0.70710678118655, -0.70710678118655 }, .r = .{ 1.41421356237310, 0.0 } },
        .{ .name = "zero_b", .a = .{ 1, 1 }, .b = .{ 0, 0 }, .c = 1.0, .s = .{ 0.0, 0.0 }, .r = .{ 1.0, 1.0 } },
        .{ .name = "zero_real", .a = .{ 0, 1 }, .b = .{ 0, 1 }, .c = 0.70710678118654, .s = .{ 0.70710678118654, 0.0 }, .r = .{ 0.0, 1.41421356237309 } },
        .{ .name = "positive_real_positive_img", .a = .{ 3, 4 }, .b = .{ 4, 6 }, .c = 0.56980288229818, .s = .{ 0.82051615050939, -0.04558423058385 }, .r = .{ 5.26497863243527, 7.01997150991369 } },
        .{ .name = "negative_real_positive_img", .a = .{ -3, 4 }, .b = .{ -4, 6 }, .c = 0.56980288229818, .s = .{ 0.82051615050939, 0.04558423058385 }, .r = .{ -5.26497863243527, 7.01997150991369 } },
        .{ .name = "positive_real_negative_img", .a = .{ 3, -4 }, .b = .{ 4, -6 }, .c = 0.56980288229818, .s = .{ 0.82051615050939, 0.04558423058385 }, .r = .{ 5.26497863243527, -7.01997150991369 } },
        .{ .name = "negative_real_negative_img", .a = .{ -3, -4 }, .b = .{ -4, -6 }, .c = 0.56980288229818, .s = .{ 0.82051615050939, -0.04558423058385 }, .r = .{ -5.26497863243527, -7.01997150991369 } },
    };
    inline for (.{ C32, C64 }) |T| for (cases) |c| {
        var a = val(T, c.a);
        var cc: Float(T) = undefined;
        var s: T = undefined;
        blas.zrotg(T, &a, val(T, c.b), &cc, &s);
        const t = tol(T, single_eps, double_eps);
        const got = [_]f64{ cc, s.re, s.im, a.re, a.im };
        const want = [_]f64{ c.c, c.s[0], c.s[1], c.r[0], c.r[1] };
        for (got, want) |g, w| if (!(@abs(g - w) <= t)) {
            ob.say("{s} {s}: got {d}, want {d}\n", .{ @typeName(T), c.name, g, w });
            return error.TestUnexpectedResult;
        };
    };
}

test "utest_ext: complex rotg" {
    try complexRotg();
}

pub fn rotmg() !void {
    const Rotmg = struct { name: []const u8, in: [4]f64, out: [4]f64, param: [5]f64 };
    const cases = [_]Rotmg{
        .{ .name = "y1_zero", .in = .{ 2, 2, 8, 0 }, .out = .{ 2, 2, 8, 0 }, .param = .{ -2, 0, 0, 0, 0 } },
        .{ .name = "d1_negative", .in = .{ -1, 2, 8, 8 }, .out = .{ 0, 0, 0, 8 }, .param = .{ -1, 0, 0, 0, 0 } },
        .{ .name = "d1_positive_d2_positive_x1_zero", .in = .{ 2, 2, 0, 8 }, .out = .{ 2, 2, 8, 8 }, .param = .{ 1, 0, 0, 0, 0 } },
        .{ .name = "scaled_y_greater_than_scaled_x", .in = .{ 1, -2, 8, 8 }, .out = .{ 0, 0, 0, 8 }, .param = .{ -1, 0, 0, 0, 0 } },
    };
    inline for (.{ f32, f64 }) |F| for (cases) |c| {
        var d1: F = @floatCast(c.in[0]);
        var d2: F = @floatCast(c.in[1]);
        var x1: F = @floatCast(c.in[2]);
        const y1: F = @floatCast(c.in[3]);
        var param: [5]F = @splat(0);
        blas.rotmg(F, &d1, &d2, &x1, y1, &param);
        const t = tol(F, single_eps, double_eps);
        const got = [_]f64{ d1, d2, x1, y1, param[0], param[1], param[2], param[3], param[4] };
        const want = c.out ++ c.param;
        for (got, want) |g, w| if (!(@abs(g - w) <= t)) {
            ob.say("{s} {s}: got {d}, want {d}\n", .{ @typeName(F), c.name, g, w });
            return error.TestUnexpectedResult;
        };
    };
}

test "utest_ext: rotmg" {
    try rotmg();
}

/// `test_daxpby.c (test_saxpby.c)`.
const axpby_real = [_]Case{
    .{ .name = "inc_x_1_inc_y_1_N_100", .kind = .run, .n = 100, .incx = 1, .incy = 1, .alpha = .{ 1.0, 0 }, .beta = .{ 1.0, 0 } },
    .{ .name = "inc_x_2_inc_y_1_N_100", .kind = .run, .n = 100, .incx = 2, .incy = 1, .alpha = .{ 2.0, 0 }, .beta = .{ 1.0, 0 } },
    .{ .name = "inc_x_1_inc_y_2_N_100", .kind = .run, .n = 100, .incx = 1, .incy = 2, .alpha = .{ 1.0, 0 }, .beta = .{ 2.0, 0 } },
    .{ .name = "inc_x_2_inc_y_2_N_100", .kind = .run, .n = 100, .incx = 2, .incy = 2, .alpha = .{ 3.0, 0 }, .beta = .{ 4.0, 0 } },
    .{ .name = "inc_x_neg_1_inc_y_2_N_100", .kind = .run, .n = 100, .incx = -1, .incy = 2, .alpha = .{ 5.0, 0 }, .beta = .{ 4.0, 0 } },
    .{ .name = "inc_x_2_inc_y_neg_1_N_100", .kind = .run, .n = 100, .incx = 2, .incy = -1, .alpha = .{ 1.0, 0 }, .beta = .{ 6.0, 0 } },
    .{ .name = "inc_x_neg_2_inc_y_neg_1_N_100", .kind = .run, .n = 100, .incx = -2, .incy = -1, .alpha = .{ 7.0, 0 }, .beta = .{ 3.5, 0 } },
    .{ .name = "inc_x_1_inc_y_1_N_100_alpha_zero", .kind = .run, .n = 100, .incx = 1, .incy = 1, .alpha = .{ 0.0, 0 }, .beta = .{ 1.0, 0 } },
    .{ .name = "inc_x_1_inc_y_2_N_100_alpha_zero", .kind = .run, .n = 100, .incx = 1, .incy = 2, .alpha = .{ 0.0, 0 }, .beta = .{ 1.0, 0 } },
    .{ .name = "inc_x_1_inc_y_1_N_100_beta_zero", .kind = .run, .n = 100, .incx = 1, .incy = 1, .alpha = .{ 1.0, 0 }, .beta = .{ 0.0, 0 } },
    .{ .name = "inc_x_2_inc_y_1_N_100_beta_zero", .kind = .run, .n = 100, .incx = 2, .incy = 1, .alpha = .{ 1.0, 0 }, .beta = .{ 0.0, 0 } },
    .{ .name = "inc_x_1_inc_y_2_N_100_beta_zero", .kind = .run, .n = 100, .incx = 1, .incy = 2, .alpha = .{ 1.0, 0 }, .beta = .{ 0.0, 0 } },
    .{ .name = "inc_x_2_inc_y_2_N_100_beta_zero", .kind = .run, .n = 100, .incx = 2, .incy = 2, .alpha = .{ 1.0, 0 }, .beta = .{ 0.0, 0 } },
    .{ .name = "inc_x_1_inc_y_1_N_100_alpha_beta_zero", .kind = .run, .n = 100, .incx = 1, .incy = 1, .alpha = .{ 0.0, 0 }, .beta = .{ 0.0, 0 } },
    .{ .name = "inc_x_1_inc_y_2_N_100_alpha_beta_zero", .kind = .run, .n = 100, .incx = 1, .incy = 2, .alpha = .{ 0.0, 0 }, .beta = .{ 0.0, 0 } },
    .{ .name = "check_n_zero", .kind = .run, .n = 0, .incx = 1, .incy = 1, .alpha = .{ 1.0, 0 }, .beta = .{ 1.0, 0 } },
};

/// `test_zaxpby.c (test_caxpby.c)`.
const axpby_complex = [_]Case{
    .{ .name = "inc_x_1_inc_y_1_N_100", .kind = .run, .n = 100, .incx = 1, .incy = 1, .alpha = .{ 1.0, 1.0 }, .beta = .{ 1.0, 1.0 } },
    .{ .name = "inc_x_2_inc_y_1_N_100", .kind = .run, .n = 100, .incx = 2, .incy = 1, .alpha = .{ 2.0, 1.0 }, .beta = .{ 1.0, 1.0 } },
    .{ .name = "inc_x_1_inc_y_2_N_100", .kind = .run, .n = 100, .incx = 1, .incy = 2, .alpha = .{ 1.0, 1.0 }, .beta = .{ 2.0, 1.0 } },
    .{ .name = "inc_x_2_inc_y_2_N_100", .kind = .run, .n = 100, .incx = 2, .incy = 2, .alpha = .{ 3.0, 1.0 }, .beta = .{ 4.0, 3.0 } },
    .{ .name = "inc_x_neg_1_inc_y_2_N_100", .kind = .run, .n = 100, .incx = -1, .incy = 2, .alpha = .{ 5.0, 2.2 }, .beta = .{ 4.0, 5.0 } },
    .{ .name = "inc_x_2_inc_y_neg_1_N_100", .kind = .run, .n = 100, .incx = 2, .incy = -1, .alpha = .{ 1.0, 1.0 }, .beta = .{ 6.0, 3.0 } },
    .{ .name = "inc_x_neg_2_inc_y_neg_1_N_100", .kind = .run, .n = 100, .incx = -2, .incy = -1, .alpha = .{ 7.0, 2.0 }, .beta = .{ 3.5, 1.3 } },
    .{ .name = "inc_x_1_inc_y_1_N_100_alpha_zero", .kind = .run, .n = 100, .incx = 1, .incy = 1, .alpha = .{ 0.0, 0.0 }, .beta = .{ 1.0, 1.0 } },
    .{ .name = "inc_x_1_inc_y_1_N_100_beta_zero", .kind = .run, .n = 100, .incx = 1, .incy = 1, .alpha = .{ 1.0, 1.0 }, .beta = .{ 0.0, 0.0 } },
    .{ .name = "inc_x_1_inc_y_1_N_100_alpha_beta_zero", .kind = .run, .n = 100, .incx = 1, .incy = 1, .alpha = .{ 0.0, 0.0 }, .beta = .{ 0.0, 0.0 } },
    .{ .name = "inc_x_1_inc_y_2_N_100_alpha_beta_zero", .kind = .run, .n = 100, .incx = 1, .incy = 2, .alpha = .{ 0.0, 0.0 }, .beta = .{ 0.0, 0.0 } },
    .{ .name = "check_n_zero", .kind = .run, .n = 0, .incx = 1, .incy = 1, .alpha = .{ 1.0, 1.0 }, .beta = .{ 1.0, 1.0 } },
};

/// `test_zaxpyc.c (test_caxpyc.c)`.
const axpyc_cases = [_]Case{
    .{ .name = "conj_strides_one", .kind = .run, .n = 100, .incx = 1, .incy = 1, .alpha = .{ 5.0, 2.2 } },
    .{ .name = "conj_incx_one", .kind = .run, .n = 100, .incx = 1, .incy = 2, .alpha = .{ 5.0, 2.2 } },
    .{ .name = "conj_incy_one", .kind = .run, .n = 100, .incx = 2, .incy = 1, .alpha = .{ 5.0, 2.2 } },
    .{ .name = "conj_strides_two", .kind = .run, .n = 100, .incx = 2, .incy = 2, .alpha = .{ 5.0, 2.2 } },
};

/// `test_zrot.c (test_crot.c); `inc_x_0_inc_y_0` is disabled there`.
const rot_cases = [_]Case{
    .{ .name = "inc_x_1_inc_y_1", .kind = .run, .n = 100, .inc_x = 1, .inc_y = 1, .c = .{ 1.0, 0.0 }, .s = .{ 1.0, 0.0 } },
    .{ .name = "inc_x_neg_1_inc_y_neg_1", .kind = .run, .n = 100, .inc_x = -1, .inc_y = -1, .c = .{ 1.0, 0.0 }, .s = .{ 1.0, 0.0 } },
    .{ .name = "inc_x_2_inc_y_1", .kind = .run, .n = 100, .inc_x = 2, .inc_y = 1, .c = .{ 3.0, 0.0 }, .s = .{ 2.0, 0.0 } },
    .{ .name = "inc_x_neg_2_inc_y_1", .kind = .run, .n = 100, .inc_x = -2, .inc_y = 1, .c = .{ 1.0, 0.0 }, .s = .{ 1.0, 0.0 } },
    .{ .name = "inc_x_1_inc_y_2", .kind = .run, .n = 100, .inc_x = 1, .inc_y = 2, .c = .{ 1.0, 0.0 }, .s = .{ 1.0, 0.0 } },
    .{ .name = "inc_x_1_inc_y_neg_2", .kind = .run, .n = 100, .inc_x = 1, .inc_y = -2, .c = .{ 2.0, 0.0 }, .s = .{ 1.0, 0.0 } },
    .{ .name = "inc_x_2_inc_y_2", .kind = .run, .n = 100, .inc_x = 2, .inc_y = 2, .c = .{ 1.0, 0.0 }, .s = .{ 2.0, 0.0 } },
    .{ .name = "inc_x_neg_2_inc_y_neg_2", .kind = .run, .n = 100, .inc_x = -2, .inc_y = -2, .c = .{ 1.0, 0.0 }, .s = .{ 1.0, 0.0 } },
    .{ .name = "inc_x_2_inc_y_2_c_zero", .kind = .run, .n = 100, .inc_x = 2, .inc_y = 2, .c = .{ 0.0, 0.0 }, .s = .{ 1.0, 0.0 } },
    .{ .name = "inc_x_2_inc_y_2_s_zero", .kind = .run, .n = 100, .inc_x = 2, .inc_y = 2, .c = .{ 1.0, 0.0 }, .s = .{ 0.0, 0.0 } },
    .{ .name = "check_n_zero", .kind = .run, .n = 0, .inc_x = 1, .inc_y = 1, .c = .{ 1.0, 0.0 }, .s = .{ 1.0, 0.0 } },
};

/// `test_zscal.c (test_cscal.c)`.
const scal_cases = [_]Case{
    .{ .name = "alpha_r_zero_alpha_i_not_zero", .kind = .run, .N = 100, .inc = 1, .alpha = .{ 0.0, 1.0 } },
    .{ .name = "alpha_r_zero_alpha_i_zero_inc_2", .kind = .run, .N = 100, .inc = 2, .alpha = .{ 0.0, 0.0 } },
};
