//! OpenBLAS's unit tests: `utest/test_amax.c`, `test_amin.c`, `test_ismin.c`,
//! `test_min.c`, `test_axpy.c`, `test_axpby.c`, `test_dnrm2.c`,
//! `test_dotu.c`, `test_dsdot.c`, `test_rot.c`, `test_rotmg.c`,
//! `test_swap.c`, `test_trmv.c`, `test_zscal.c`, `test_gemv.c`, and the
//! pending `test_potrs.c` and `test_kernel_regress.c`, OpenBLAS v0.3.34-317.
//! Inputs and expected values are the sources' own.
//!
//! `test_trmv.c`'s row-major call is its column-major equivalent. An
//! expected infinity must be met exactly, where `ASSERT_DBL_NEAR_TOL` would
//! let NaN through.

const std = @import("std");
const ob = @import("__root__.zig");
const l1 = @import("../../blas1/interface.zig");

const C32 = std.math.Complex(f32);
const C64 = std.math.Complex(f64);
const nan = std.math.nan(f64);
const inf = std.math.inf(f64);

/// `SINGLE_EPS` and `DOUBLE_EPS`.
const single_eps = 1e-4;
const double_eps = 1e-13;

/// `ASSERT_DBL_NEAR_TOL`.
fn near(want: f64, got: f64, tol: f64, what: []const u8) !void {
    const ok = if (std.math.isFinite(want)) @abs(want - got) <= tol else want == got;
    if (ok) return;
    ob.say("{s}: got {e}, want {e}\n", .{ what, got, want });
    return error.TestUnexpectedResult;
}

fn nears(comptime F: type, want: []const f64, got: []const F, tol: f64, what: []const u8) !void {
    for (want, got) |w, g| try near(w, g, tol, what);
}

/// `F`'s elements, `n` complex ones as `2n` floats.
fn flat(comptime T: type, p: []T) []ob.Float(T) {
    const F = ob.Float(T);
    return @as([*]F, @ptrCast(p.ptr))[0 .. p.len * (@sizeOf(T) / @sizeOf(F))];
}

fn arr(comptime F: type, comptime n: usize, v: [n]f64) [n]F {
    var out: [n]F = undefined;
    for (&out, v) |*o, x| o.* = @floatCast(x);
    return out;
}

/// `n` complex elements from `2n` floats.
fn carr(comptime C: type, comptime n: usize, v: [2 * n]f64) [n]C {
    var out: [n]C = undefined;
    for (&out, 0..) |*o, i| o.* = .init(@floatCast(v[2 * i]), @floatCast(v[2 * i + 1]));
    return out;
}

fn eps(comptime F: type) f64 {
    return if (F == f32) single_eps else double_eps;
}

// ---- test_amax.c, test_amin.c, test_min.c, test_ismin.c -------------------------

/// Every case here but those that wait, for `check.eachCore`.
pub fn all() !void {
    try amaxAmin();
    try minMax();
    try isminIsmax();
    try axpyWithZeroStrides();
    try axpbyStrides();
    try swapWithZeroStrides();
    try rotAndRotm();
    try drotmg();
    try dnrm2Cases();
    try zdotuDsdot();
    try strmvLower();
    try scalNanInf();
    try gemvNanInf();
    try cgemvNanInf();
}

fn amaxAmin() !void {
    inline for (.{ f32, f64 }) |F| {
        var x = arr(F, 9, .{ -1.1, 2.2, -3.3, 4.4, -5.5, 6.6, -7.7, 8.8, -9.9 });
        try near(3.3, try l1.extreme(F, .amax, 3, &x, 1), eps(F), "amax");
        try near(1.1, try l1.extreme(F, .amin, 3, &x, 1), eps(F), "amin");
    }
    inline for (.{ C32, C64 }) |C| {
        var x = carr(C, 9, .{ -1.1, 2.2, -3.3, 4.4, -5.5, 6.6, -7.7, 8.8, -9.9, 10.10, -1.1, 2.2, -3.3, 4.4, -5.5, 6.6, -7.7, 8.8 });
        try near(20.0, try l1.extreme(C, .amax, 9, &x, 1), eps(ob.Float(C)), "camax");
        try near(3.3, try l1.extreme(C, .amin, 9, &x, 1), eps(ob.Float(C)), "camin");
    }
}

const l2 = @import("../../blas2/interface.zig");
const tri = @import("../../blas2/triangular.zig");

test "utest: amax, amin" {
    try amaxAmin();
}

fn minMax() !void {
    const Case = struct { []const u8, [3]f64, f64 };
    const mins = [_]Case{ .{ "smin_negative", .{ -1.1, -2.2, -3.3 }, -3.3 }, .{ "dmin_positive", .{ 1.1, 0.0, 3.3 }, 0.0 }, .{ "smin_zero", .{ 1.1, 2.2, 0.0 }, 0.0 } };
    const maxs = [_]Case{ .{ "smax_negative", .{ -1.1, -2.2, -3.3 }, -1.1 }, .{ "dmax_positive", .{ 1.1, 0.0, 3.3 }, 3.3 }, .{ "smax_zero", .{ -1.1, -2.2, 0.0 }, 0.0 } };
    inline for (.{ f32, f64 }) |F| {
        for (mins) |c| {
            var x = arr(F, 3, c[1]);
            try near(c[2], try l1.extreme(F, .min, 3, &x, 1), eps(F), c[0]);
        }
        for (maxs) |c| {
            var x = arr(F, 3, c[1]);
            try near(c[2], try l1.extreme(F, .max, 3, &x, 1), eps(F), c[0]);
        }
    }
}

test "utest: min, max" {
    try minMax();
}

fn isminIsmax() !void {
    const n = 50;
    const inc = 2;
    var x: [n * inc]f32 = undefined;
    for (&x, 0..) |*v, i| v.* = @floatFromInt(i + 1000);
    x[8 * inc] = 0;
    try ob.expect(try l1.index(f32, .min, n, &x, inc) == 9, "ismin positive_step_2");
    for (&x, 0..) |*v, i| v.* = -@as(f32, @floatFromInt(i)) - 1000;
    x[8 * inc] = -123456.0;
    try ob.expect(try l1.index(f32, .min, n, &x, inc) == 9, "ismin negative_step_2");
    for (&x, 0..) |*v, i| v.* = @floatFromInt(i + 1000);
    x[8 * inc] = 123456.0;
    try ob.expect(try l1.index(f32, .max, n, &x, inc) == 9, "ismax positive_step_2");
    for (&x, 0..) |*v, i| v.* = -@as(f32, @floatFromInt(i)) - 1000;
    x[8 * inc] = 0;
    try ob.expect(try l1.index(f32, .max, n, &x, inc) == 9, "ismax negative_step_2");
}

test "utest: ismin, ismax at stride 2" {
    try isminIsmax();
}

// ---- test_axpy.c, test_axpby.c, test_swap.c, test_rot.c ---------------------------

const xs = [8]f64{ 1.0, 3.0, 5.0, 7.0, 1.0, 3.0, 5.0, 7.0 };
const ys = [8]f64{ 2.0, 4.0, 6.0, 8.0, 2.0, 4.0, 6.0, 8.0 };

fn axpyWithZeroStrides() !void {
    // the reference's tolerance is DOUBLE_EPS throughout
    inline for (.{ f32, f64 }) |F| {
        var x = arr(F, 8, xs);
        var y = arr(F, 8, ys);
        try l1.axpy(F, 8, 0.25, &x, 0, &y, 0);
        try nears(F, &xs, &x, double_eps, "axpy_inc_0 x");
        try nears(F, &.{ 4.0, 4.0, 6.0, 8.0, 2.0, 4.0, 6.0, 8.0 }, &y, double_eps, "axpy_inc_0 y");
    }
    inline for (.{ C32, C64 }) |C| {
        var x = carr(C, 4, xs);
        var y = carr(C, 4, ys);
        try l1.axpy(C, 4, .init(0.25, 0.5), &x, 0, &y, 0);
        try nears(ob.Float(C), &xs, flat(C, &x), double_eps, "zaxpy_inc_0 x");
        try nears(ob.Float(C), &.{ -3.0, 9.0, 6.0, 8.0, 2.0, 4.0, 6.0, 8.0 }, flat(C, &y), double_eps, "zaxpy_inc_0 y");
        x = carr(C, 4, xs);
        y = carr(C, 4, ys);
        try l1.axpy(C, 4, .init(0.25, 0.5), &x, 0, &y, 1);
        try nears(ob.Float(C), &xs, flat(C, &x), double_eps, "zaxpy_incx_0 x");
        try nears(ob.Float(C), &.{ 0.75, 5.25, 4.75, 9.25, 0.75, 5.25, 4.75, 9.25 }, flat(C, &y), double_eps, "zaxpy_incx_0 y");
    }
}

test "utest: axpy with zero strides" {
    try axpyWithZeroStrides();
}

const bx = [9]f64{ 1.0, 3.0, 5.0, 7.0, 1.0, 3.0, 5.0, 7.0, 9.0 };
const by = [9]f64{ 2.0, 4.0, 6.0, 8.0, 2.0, 4.0, 6.0, 8.0, 10.0 };

fn axpbyStrides() !void {
    inline for (.{ f32, f64 }) |F| {
        const e = eps(F);
        var x = arr(F, 9, bx);
        var y = arr(F, 9, by);
        try l1.axpby(F, 9, 1.0, &x, 0, 2.0, &y, 0);
        try nears(F, &bx, &x, e, "axpby_inc_0 x");
        try nears(F, &.{ 1535.0, 4.0, 6.0, 8.0, 2.0, 4.0, 6.0, 8.0, 10.0 }, &y, e, "axpby_inc_0 y");
        x = arr(F, 9, bx);
        y = arr(F, 9, by);
        try l1.axpby(F, 9, 0.25, &x, 1, 0.75, &y, 1);
        try nears(F, &.{ 1.75, 3.75, 5.75, 7.75, 1.75, 3.75, 5.75, 7.75, 9.75 }, &y, e, "axpby_inc_1 y");
        var x2 = arr(F, 18, bx ++ bx);
        var y2 = arr(F, 18, by ++ by);
        try l1.axpby(F, 9, 0.25, &x2, 2, 0.75, &y2, 2);
        try nears(F, &(bx ++ bx), &x2, e, "axpby_inc_2 x");
        try nears(F, &.{ 1.75, 4.00, 5.75, 8.00, 1.75, 4.00, 5.75, 8.00, 9.75, 2.00, 3.75, 6.00, 7.75, 2.00, 3.75, 6.00, 7.75, 10.00 }, &y2, e, "axpby_inc_2 y");
    }
    inline for (.{ C32, C64 }) |C| {
        const F = ob.Float(C);
        const e = eps(F);
        var x = carr(C, 9, bx ++ bx);
        var y = carr(C, 9, by ++ by);
        try l1.axpby(C, 9, .init(1.0, 2.0), &x, 0, .init(2.0, 1.0), &y, 0);
        try nears(F, &(bx ++ bx), flat(C, &x), e, "zaxpby_inc_0 x");
        try nears(F, &.{ 9355.0, -8865.0, 6.0, 8.0, 2.0, 4.0, 6.0, 8.0, 10.0, 2.0, 4.0, 6.0, 8.0, 2.0, 4.0, 6.0, 8.0, 10.0 }, flat(C, &y), e, "zaxpby_inc_0 y");
        x = carr(C, 9, bx ++ bx);
        y = carr(C, 9, by ++ by);
        try l1.axpby(C, 9, .init(0.25, 0.25), &x, 1, .init(0.75, 0.75), &y, 1);
        try nears(F, &.{ -2.0, 5.5, -2.0, 13.5, -2.0, 5.5, -2.0, 13.5, 8.0, 11.5, -2.0, 9.5, 6.0, 9.5, -2.0, 9.5, -2.0, 17.5 }, flat(C, &y), e, "zaxpby_inc_1 y");
        var x2 = carr(C, 18, bx ++ bx ++ bx ++ bx);
        var y2 = carr(C, 18, by ++ by ++ by ++ by);
        try l1.axpby(C, 9, .init(0.25, 0.25), &x2, 2, .init(0.75, 0.75), &y2, 2);
        try nears(F, &(bx ++ bx ++ bx ++ bx), flat(C, &x2), e, "zaxpby_inc_2 x");
        try nears(F, &.{
            -2.0, 5.5,  6.0,  8.0,  -2.0, 5.5,  6.0,  8.0,  8.0,
            11.5, 4.0,  6.0,  6.0,  9.5,  4.0,  6.0,  -2.0, 17.5,
            2.0,  4.0,  -2.0, 13.5, 2.0,  4.0,  -2.0, 13.5, 10.0,
            2.0,  -2.0, 9.5,  8.0,  2.0,  -2.0, 9.5,  8.0,  10.0,
        }, flat(C, &y2), e, "zaxpby_inc_2 y");
    }
}

test "utest: axpby at strides 0, 1 and 2" {
    try axpbyStrides();
}

fn swapWithZeroStrides() !void {
    inline for (.{ f32, f64 }) |F| {
        var x = arr(F, 4, xs[0..4].*);
        var y = arr(F, 4, ys[0..4].*);
        try l1.swap(F, 4, &x, 0, &y, 0);
        try nears(F, xs[0..4], &x, eps(F), "swap_inc_0 x");
        try nears(F, ys[0..4], &y, eps(F), "swap_inc_0 y");
    }
    inline for (.{ C32, C64 }) |C| {
        var x = carr(C, 4, xs);
        var y = carr(C, 4, ys);
        try l1.swap(C, 4, &x, 0, &y, 0);
        try nears(ob.Float(C), &xs, flat(C, &x), eps(ob.Float(C)), "zswap_inc_0 x");
        try nears(ob.Float(C), &ys, flat(C, &y), eps(ob.Float(C)), "zswap_inc_0 y");
    }
}

test "utest: swap with zero strides" {
    try swapWithZeroStrides();
}

fn rotAndRotm() !void {
    inline for (.{ f32, f64 }) |F| {
        const e = eps(F);
        var x = arr(F, 4, xs[0..4].*);
        var y = arr(F, 4, ys[0..4].*);
        try l1.rot(F, 4, &x, 0, &y, 0, 0.25, 0.5);
        try nears(F, &.{ -0.21484375000000, 3.0, 5.0, 7.0 }, &x, e, "rot_inc_0 x");
        try nears(F, &.{ 0.03906250000000, 4.0, 6.0, 8.0 }, &y, e, "rot_inc_0 y");
        x = arr(F, 4, xs[0..4].*);
        y = arr(F, 4, ys[0..4].*);
        try l1.rot(F, 4, &x, 1, &y, 1, 1.0, 1.0);
        try nears(F, &.{ 3.0, 7.0, 11.0, 15.0 }, &x, e, "rot_inc_1 x");
        try nears(F, &.{ 1.0, 1.0, 1.0, 1.0 }, &y, e, "rot_inc_1 y");
        const twelve = [12]f64{ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0 };
        var mx = arr(F, 12, twelve);
        var my = arr(F, 12, twelve);
        const param = arr(F, 5, .{ 1.0, 2.0, 3.0, 4.0, 5.0 });
        try l1.rotm(F, 12, &mx, 1, &my, 1, &param);
        try nears(F, &.{ 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0, 36.0 }, &mx, e, "rotm_inc_1 x");
        try nears(F, &.{ 4.0, 8.0, 12.0, 16.0, 20.0, 24.0, 28.0, 32.0, 36.0, 40.0, 44.0, 48.0 }, &my, e, "rotm_inc_1 y");
    }
    inline for (.{ C32, C64 }) |C| {
        var x = carr(C, 4, xs);
        var y = carr(C, 4, ys);
        try l1.rot(C, 4, &x, 0, &y, 0, 0.25, 0.5);
        try nears(ob.Float(C), &.{ -0.21484375000000, -0.45703125000000, 5.0, 7.0, 1.0, 3.0, 5.0, 7.0 }, flat(C, &x), eps(ob.Float(C)), "zdrot_inc_0 x");
        try nears(ob.Float(C), &.{ 0.03906250000000, 0.17187500000000, 6.0, 8.0, 2.0, 4.0, 6.0, 8.0 }, flat(C, &y), eps(ob.Float(C)), "zdrot_inc_0 y");
    }
}

test "utest: rot and rotm" {
    try rotAndRotm();
}

// ---- test_rotmg.c, test_dnrm2.c, test_dotu.c, test_dsdot.c -------------------------

fn drotmg() !void {
    const Case = struct { name: []const u8, in: [4]f64, param: [5]f64, out: [4]f64, want: [5]f64 };
    const cases = [_]Case{
        // libGoto's bug, fixed by the 2014 rewrite; netlib's values
        .{ .name = "rotmg", .in = .{ 0.21149573940783739, 0.046892057172954082, -0.42272687517106533, 0.42211309121921659 }, .param = @splat(0), .out = .{ 0.17320483687975, 0.03840233915037, -0.51618034832329, 0.42211309121922 }, .want = .{ 0.0, 0.0, 0.99854803659786, -0.22139439665872, 0.0 } },
        // issue 1452, gonum's values
        .{ .name = "rotmg_issue1452", .in = .{ 5.9e-8, 5.960464e-8, 1.0, 150.0 }, .param = .{ 0, 0, 0, 1.0 / 4096.0, 0 }, .out = .{ 0.99995592822897, 0.98981219860583, 0.03662270484346, 150.000000000000 }, .want = .{ -1.0, 0.00000161109346, -0.00024414062500, 0.00024414062500, 0.00000162760417 } },
        .{ .name = "rotmg_D1eqD2_X1eqX2", .in = .{ 2, 2, 8, 8 }, .param = @splat(0), .out = .{ 1.0, 1.0, 16.0, 8.0 }, .want = .{ 1.0, 1.0, 0.0, 0.0, 1.0 } },
        .{ .name = "drotmg_D1_big_D2_big_flag_zero", .in = .{ 1600000000.0, 800000000.0, 8.0, 7.0 }, .param = .{ 1.0, 4096.0, -4096.0, 1.0, 4096.0 }, .out = .{ 68.96627824858757, 34.483139124293785, 45312.0, 7.0 }, .want = .{ -1.0, 4096.0, -3584.0, 1792.0, 4096.0 } },
    };
    for (cases) |c| {
        var d1 = c.in[0];
        var d2 = c.in[1];
        var x1 = c.in[2];
        const y1 = c.in[3];
        var param = c.param;
        l1.rotmg(f64, &d1, &d2, &x1, y1, &param);
        try nears(f64, &c.out, &.{ d1, d2, x1, y1 }, double_eps, c.name);
        try nears(f64, &c.want, &param, double_eps, c.name);
    }
}

test "utest: drotmg" {
    try drotmg();
}

fn dnrm2Cases() !void {
    var x: [29]f64 = @splat(0);
    x[10] = -inf;
    try near(inf, try l1.nrm2(f64, 28, &x, 1), double_eps, "dnrm2_inf");
    @memset(&x, 7.457008414e-310);
    try near(0, try l1.nrm2(f64, 28, &x, 1), double_eps, "dnrm2_tiny");
    var y: [5]f64 = @splat(10.0);
    try near(@sqrt(500.0), try l1.nrm2(f64, 5, &y, -1), double_eps, "dnrm2_neg_incx");
}

test "utest: dnrm2 of infinity, of denormals, and at a negative stride" {
    try dnrm2Cases();
}

fn zdotuDsdot() !void {
    var x1 = [_]C64{.init(1.0, 1.0)};
    var y1 = [_]C64{.init(1.0, 2.0)};
    const d1 = try l1.dot(C64, 1, &x1, 1, &y1, 1);
    try nears(f64, &.{ -1.0, 3.0 }, &.{ d1.re, d1.im }, double_eps, "zdotu_n_1");
    // one double past the start: (2, 3) and (6, 7)
    var x2 = [_]f64{ 1.0, 2.0, 3.0, 4.0 };
    var y2 = [_]f64{ 5.0, 6.0, 7.0, 8.0 };
    const px: [*]C64 = @ptrCast(&x2[1]);
    const py: [*]C64 = @ptrCast(&y2[1]);
    const d2 = try l1.dot(C64, 1, px, 1, py, 1);
    try nears(f64, &.{ -9.0, 32.0 }, &.{ d2.re, d2.im }, double_eps, "zdotu_offset_1");
    var sx = [_]f32{0.172555164};
    var sy = [_]f32{-0.0138700781};
    try near(-0.00239335360107, try l1.dsdot(1, &sx, 1, &sy, 1), double_eps, "dsdot_n_1");
}

test "utest: zdotu, dsdot" {
    try zdotuDsdot();
}

// ---- test_trmv.c --------------------------------------------------------------------

fn strmvLower() !void {
    // row-major lower, not transposed: column-major upper, transposed
    var lower = [9]f32{ 2.0, 0.0, 0.0, 3.0, 4.0, 0.0, 5.0, 6.0, 7.0 };
    var v = [3]f32{ 11.0, 13.0, 17.0 };
    try tri.trmv(f32, std.testing.allocator, .upper, .t, .non_unit, 3, &lower, 3, &v, 1);
    try nears(f32, &.{ 22.0, 85.0, 252.0 }, &v, single_eps, "strmv");
}

test "utest: strmv, a small lower triangle at unit stride" {
    try strmvLower();
}

// ---- test_zscal.c -------------------------------------------------------------------

const Want = enum { nan, inf, zero };

fn is(w: Want, v: anytype) bool {
    return switch (w) {
        .nan => std.math.isNan(v),
        .inf => std.math.isInf(v),
        .zero => v == 0,
    };
}

fn scalNanInf() !void {
    // alpha, the value x holds; every result NaN; at strides 1 and 2
    const real = [_][2]f64{ .{ 0, nan }, .{ nan, 0 }, .{ 0, inf }, .{ inf, 0 }, .{ nan, inf }, .{ inf, nan } };
    inline for (.{ f32, f64 }) |F| for (real) |c| for ([_]isize{ 1, 2 }) |inc| {
        var x: [20]F = @splat(@floatCast(c[1]));
        try l1.scal(F, 9, @floatCast(c[0]), &x, inc);
        try ob.expect(std.math.isNan(x[0]) and std.math.isNan(x[8]), "sscal");
    };
    // alpha, x's elements, and what raw floats 0, 1, 16, 17 become
    const Case = struct { name: []const u8, alpha: [2]f64, x: [2]f64, x2: [2]f64, want: [2]Want };
    const cases = [_]Case{
        .{ .name = "i_nan", .alpha = .{ 0, 1 }, .x = .{ nan, 0 }, .x2 = .{ nan, 0 }, .want = .{ .nan, .nan } },
        .{ .name = "nan_i", .alpha = .{ nan, 0 }, .x = .{ 0, 1 }, .x2 = .{ 0, 1 }, .want = .{ .nan, .nan } },
        .{ .name = "i_inf", .alpha = .{ 0, 1 }, .x = .{ inf, 0 }, .x2 = .{ inf, 0 }, .want = .{ .nan, .inf } },
        .{ .name = "inf_i", .alpha = .{ inf, 0 }, .x = .{ 0, 1 }, .x2 = .{ 0, 1 }, .want = .{ .nan, .inf } },
        .{ .name = "i_0inf", .alpha = .{ 0, 1 }, .x = .{ 0, inf }, .x2 = .{ 0, inf }, .want = .{ .inf, .nan } },
        // the stride-2 case holds (0, NaN)
        .{ .name = "i00_NAN", .alpha = .{ 0, 0 }, .x = .{ nan, 0 }, .x2 = .{ 0, nan }, .want = .{ .nan, .nan } },
    };
    inline for (.{ C32, C64 }) |C| for (cases) |c| for ([_]isize{ 1, 2 }) |inc| {
        const v = if (inc == 1) c.x else c.x2;
        var x: [18]C = @splat(.init(@floatCast(v[0]), @floatCast(v[1])));
        try l1.scal(C, 9, .init(@floatCast(c.alpha[0]), @floatCast(c.alpha[1])), &x, inc);
        const f = flat(C, &x);
        try ob.expect(is(c.want[0], f[0]) and is(c.want[1], f[1]) and is(c.want[0], f[16]) and is(c.want[1], f[17]), c.name);
    };
}

test "utest: scal of NaN and infinity" {
    try scalNanInf();
}

// ---- test_gemv.c --------------------------------------------------------------------

fn gemvNanInf() !void {
    const al = std.testing.allocator;
    const n = 17;
    inline for (.{ f32, f64 }) |F| for ([_]isize{ 1, 2 }) |inc| {
        var a: [n * n]F = @splat(0);
        var x: [n]F = @splat(0);
        var y: [2 * n - 1]F = @splat(0);
        const ny: usize = if (inc == 1) n else 2 * n - 1;
        for (0..n) |i| y[i * @as(usize, @intCast(inc))] = if (i % 2 == 0) std.math.nan(F) else std.math.inf(F);
        try l2.gemv(al, F, .n, n, n, 0, &a, n, &x, 1, 0, &y, inc);
        for (y[0..ny]) |v| try ob.expect(v == 0, "sgemv 0_nan_inf");
    };
}

test "utest: gemv, zero alpha and beta overwrite NaN and infinity" {
    try gemvNanInf();
}

fn cgemvNanInf() !void {
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
            try l2.gemv(al, C, .n, n, n, .init(0, 0), &a, n, &x, 1, c.beta, &y, inc);
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

// ---- pending ------------------------------------------------------------------------

test "utest: potrf, potrs (waits for LAPACK)" {
    return error.SkipZigTest;
}

test "utest: kernel_regress skx_avx (waits for LAPACK's dgesvd)" {
    return error.SkipZigTest;
}
