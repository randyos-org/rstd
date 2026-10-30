//! OpenBLAS's unit tests: `utest/test_amax.c`, `test_amin.c`, `test_ismin.c`,
//! `test_min.c`, `test_axpy.c`, `test_axpby.c`, `test_dnrm2.c`,
//! `test_dotu.c`, `test_dsdot.c`, `test_rot.c`, `test_rotmg.c`,
//! `test_swap.c`, `test_trmv.c`, `test_zscal.c`, `test_gemv.c`,
//! `test_potrs.c` and `test_kernel_regress.c` (through
//! `rstd.linalg.lapack`), OpenBLAS v0.3.34-317. Inputs and expected values
//! are the sources' own.
//!
//! `test_trmv.c`'s row-major call is its column-major equivalent. An
//! expected infinity must be met exactly, where `ASSERT_DBL_NEAR_TOL` would
//! let NaN through.

const std = @import("std");
const ob = @import("__root__.zig");
const blas = @import("../../__root__.zig");
const lapack = @import("../../../lapack/__root__.zig");

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

/// `test_potrs.c`, which runs with the level-3 suite: LAPACK reaches the
/// level-3 drivers, too big for the level-1 suites' rig images.
pub fn lapackAll() !void {
    try potrfBug695();
    try potrfSmoke();
}

fn amaxAmin() !void {
    inline for (.{ f32, f64 }) |F| {
        var x = arr(F, 9, .{ -1.1, 2.2, -3.3, 4.4, -5.5, 6.6, -7.7, 8.8, -9.9 });
        try near(3.3, try blas.extreme(ob.io, F, .amax, 3, &x, 1), eps(F), "amax");
        try near(1.1, try blas.extreme(ob.io, F, .amin, 3, &x, 1), eps(F), "amin");
    }
    inline for (.{ C32, C64 }) |C| {
        var x = carr(C, 9, .{ -1.1, 2.2, -3.3, 4.4, -5.5, 6.6, -7.7, 8.8, -9.9, 10.10, -1.1, 2.2, -3.3, 4.4, -5.5, 6.6, -7.7, 8.8 });
        try near(20.0, try blas.extreme(ob.io, C, .amax, 9, &x, 1), eps(ob.Float(C)), "camax");
        try near(3.3, try blas.extreme(ob.io, C, .amin, 9, &x, 1), eps(ob.Float(C)), "camin");
    }
}

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
            try near(c[2], try blas.extreme(ob.io, F, .min, 3, &x, 1), eps(F), c[0]);
        }
        for (maxs) |c| {
            var x = arr(F, 3, c[1]);
            try near(c[2], try blas.extreme(ob.io, F, .max, 3, &x, 1), eps(F), c[0]);
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
    try ob.expect(try blas.index(ob.io, f32, .min, n, &x, inc) == 9, "ismin positive_step_2");
    for (&x, 0..) |*v, i| v.* = -@as(f32, @floatFromInt(i)) - 1000;
    x[8 * inc] = -123456.0;
    try ob.expect(try blas.index(ob.io, f32, .min, n, &x, inc) == 9, "ismin negative_step_2");
    for (&x, 0..) |*v, i| v.* = @floatFromInt(i + 1000);
    x[8 * inc] = 123456.0;
    try ob.expect(try blas.index(ob.io, f32, .max, n, &x, inc) == 9, "ismax positive_step_2");
    for (&x, 0..) |*v, i| v.* = -@as(f32, @floatFromInt(i)) - 1000;
    x[8 * inc] = 0;
    try ob.expect(try blas.index(ob.io, f32, .max, n, &x, inc) == 9, "ismax negative_step_2");
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
        try blas.axpy(ob.io, F, 8, 0.25, &x, 0, &y, 0);
        try nears(F, &xs, &x, double_eps, "axpy_inc_0 x");
        try nears(F, &.{ 4.0, 4.0, 6.0, 8.0, 2.0, 4.0, 6.0, 8.0 }, &y, double_eps, "axpy_inc_0 y");
    }
    inline for (.{ C32, C64 }) |C| {
        var x = carr(C, 4, xs);
        var y = carr(C, 4, ys);
        try blas.axpy(ob.io, C, 4, .init(0.25, 0.5), &x, 0, &y, 0);
        try nears(ob.Float(C), &xs, flat(C, &x), double_eps, "zaxpy_inc_0 x");
        try nears(ob.Float(C), &.{ -3.0, 9.0, 6.0, 8.0, 2.0, 4.0, 6.0, 8.0 }, flat(C, &y), double_eps, "zaxpy_inc_0 y");
        x = carr(C, 4, xs);
        y = carr(C, 4, ys);
        try blas.axpy(ob.io, C, 4, .init(0.25, 0.5), &x, 0, &y, 1);
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
        try blas.axpby(ob.io, F, 9, 1.0, &x, 0, 2.0, &y, 0);
        try nears(F, &bx, &x, e, "axpby_inc_0 x");
        try nears(F, &.{ 1535.0, 4.0, 6.0, 8.0, 2.0, 4.0, 6.0, 8.0, 10.0 }, &y, e, "axpby_inc_0 y");
        x = arr(F, 9, bx);
        y = arr(F, 9, by);
        try blas.axpby(ob.io, F, 9, 0.25, &x, 1, 0.75, &y, 1);
        try nears(F, &.{ 1.75, 3.75, 5.75, 7.75, 1.75, 3.75, 5.75, 7.75, 9.75 }, &y, e, "axpby_inc_1 y");
        var x2 = arr(F, 18, bx ++ bx);
        var y2 = arr(F, 18, by ++ by);
        try blas.axpby(ob.io, F, 9, 0.25, &x2, 2, 0.75, &y2, 2);
        try nears(F, &(bx ++ bx), &x2, e, "axpby_inc_2 x");
        try nears(F, &.{ 1.75, 4.00, 5.75, 8.00, 1.75, 4.00, 5.75, 8.00, 9.75, 2.00, 3.75, 6.00, 7.75, 2.00, 3.75, 6.00, 7.75, 10.00 }, &y2, e, "axpby_inc_2 y");
    }
    inline for (.{ C32, C64 }) |C| {
        const F = ob.Float(C);
        const e = eps(F);
        var x = carr(C, 9, bx ++ bx);
        var y = carr(C, 9, by ++ by);
        try blas.axpby(ob.io, C, 9, .init(1.0, 2.0), &x, 0, .init(2.0, 1.0), &y, 0);
        try nears(F, &(bx ++ bx), flat(C, &x), e, "zaxpby_inc_0 x");
        try nears(F, &.{ 9355.0, -8865.0, 6.0, 8.0, 2.0, 4.0, 6.0, 8.0, 10.0, 2.0, 4.0, 6.0, 8.0, 2.0, 4.0, 6.0, 8.0, 10.0 }, flat(C, &y), e, "zaxpby_inc_0 y");
        x = carr(C, 9, bx ++ bx);
        y = carr(C, 9, by ++ by);
        try blas.axpby(ob.io, C, 9, .init(0.25, 0.25), &x, 1, .init(0.75, 0.75), &y, 1);
        try nears(F, &.{ -2.0, 5.5, -2.0, 13.5, -2.0, 5.5, -2.0, 13.5, 8.0, 11.5, -2.0, 9.5, 6.0, 9.5, -2.0, 9.5, -2.0, 17.5 }, flat(C, &y), e, "zaxpby_inc_1 y");
        var x2 = carr(C, 18, bx ++ bx ++ bx ++ bx);
        var y2 = carr(C, 18, by ++ by ++ by ++ by);
        try blas.axpby(ob.io, C, 9, .init(0.25, 0.25), &x2, 2, .init(0.75, 0.75), &y2, 2);
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
        try blas.swap(ob.io, F, 4, &x, 0, &y, 0);
        try nears(F, xs[0..4], &x, eps(F), "swap_inc_0 x");
        try nears(F, ys[0..4], &y, eps(F), "swap_inc_0 y");
    }
    inline for (.{ C32, C64 }) |C| {
        var x = carr(C, 4, xs);
        var y = carr(C, 4, ys);
        try blas.swap(ob.io, C, 4, &x, 0, &y, 0);
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
        try blas.rot(ob.io, F, 4, &x, 0, &y, 0, 0.25, 0.5);
        try nears(F, &.{ -0.21484375000000, 3.0, 5.0, 7.0 }, &x, e, "rot_inc_0 x");
        try nears(F, &.{ 0.03906250000000, 4.0, 6.0, 8.0 }, &y, e, "rot_inc_0 y");
        x = arr(F, 4, xs[0..4].*);
        y = arr(F, 4, ys[0..4].*);
        try blas.rot(ob.io, F, 4, &x, 1, &y, 1, 1.0, 1.0);
        try nears(F, &.{ 3.0, 7.0, 11.0, 15.0 }, &x, e, "rot_inc_1 x");
        try nears(F, &.{ 1.0, 1.0, 1.0, 1.0 }, &y, e, "rot_inc_1 y");
        const twelve = [12]f64{ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0 };
        var mx = arr(F, 12, twelve);
        var my = arr(F, 12, twelve);
        const param = arr(F, 5, .{ 1.0, 2.0, 3.0, 4.0, 5.0 });
        try blas.rotm(ob.io, F, 12, &mx, 1, &my, 1, &param);
        try nears(F, &.{ 3.0, 6.0, 9.0, 12.0, 15.0, 18.0, 21.0, 24.0, 27.0, 30.0, 33.0, 36.0 }, &mx, e, "rotm_inc_1 x");
        try nears(F, &.{ 4.0, 8.0, 12.0, 16.0, 20.0, 24.0, 28.0, 32.0, 36.0, 40.0, 44.0, 48.0 }, &my, e, "rotm_inc_1 y");
    }
    inline for (.{ C32, C64 }) |C| {
        var x = carr(C, 4, xs);
        var y = carr(C, 4, ys);
        try blas.rot(ob.io, C, 4, &x, 0, &y, 0, 0.25, 0.5);
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
        blas.rotmg(f64, &d1, &d2, &x1, y1, &param);
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
    try near(inf, try blas.nrm2(ob.io, f64, 28, &x, 1), double_eps, "dnrm2_inf");
    @memset(&x, 7.457008414e-310);
    try near(0, try blas.nrm2(ob.io, f64, 28, &x, 1), double_eps, "dnrm2_tiny");
    var y: [5]f64 = @splat(10.0);
    try near(@sqrt(500.0), try blas.nrm2(ob.io, f64, 5, &y, -1), double_eps, "dnrm2_neg_incx");
}

test "utest: dnrm2 of infinity, of denormals, and at a negative stride" {
    try dnrm2Cases();
}

fn zdotuDsdot() !void {
    var x1 = [_]C64{.init(1.0, 1.0)};
    var y1 = [_]C64{.init(1.0, 2.0)};
    const d1 = try blas.dot(ob.io, C64, 1, &x1, 1, &y1, 1);
    try nears(f64, &.{ -1.0, 3.0 }, &.{ d1.re, d1.im }, double_eps, "zdotu_n_1");
    // one double past the start: (2, 3) and (6, 7)
    var x2 = [_]f64{ 1.0, 2.0, 3.0, 4.0 };
    var y2 = [_]f64{ 5.0, 6.0, 7.0, 8.0 };
    const px: [*]C64 = @ptrCast(&x2[1]);
    const py: [*]C64 = @ptrCast(&y2[1]);
    const d2 = try blas.dot(ob.io, C64, 1, px, 1, py, 1);
    try nears(f64, &.{ -9.0, 32.0 }, &.{ d2.re, d2.im }, double_eps, "zdotu_offset_1");
    var sx = [_]f32{0.172555164};
    var sy = [_]f32{-0.0138700781};
    try near(-0.00239335360107, try blas.dsdot(ob.io, 1, &sx, 1, &sy, 1), double_eps, "dsdot_n_1");
}

test "utest: zdotu, dsdot" {
    try zdotuDsdot();
}

// ---- test_trmv.c --------------------------------------------------------------------

fn strmvLower() !void {
    // row-major lower, not transposed: column-major upper, transposed
    var lower = [9]f32{ 2.0, 0.0, 0.0, 3.0, 4.0, 0.0, 5.0, 6.0, 7.0 };
    var v = [3]f32{ 11.0, 13.0, 17.0 };
    try blas.trmv(std.testing.allocator, ob.io, f32, .u, .t, .n, 3, &lower, 3, &v, 1);
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
        try blas.scal(ob.io, F, 9, @floatCast(c[0]), &x, inc);
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
        try blas.scal(ob.io, C, 9, .init(@floatCast(c.alpha[0]), @floatCast(c.alpha[1])), &x, inc);
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
        try blas.gemv(al, ob.io, F, .n, n, n, 0, &a, n, &x, 1, 0, &y, inc);
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
            try blas.gemv(al, ob.io, C, .n, n, n, .init(0, 0), &a, n, &x, 1, c.beta, &y, inc);
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

// ---- test_potrs.c ------------------------------------------------------------------

/// `potrf.bug_695`: `cpotrf` of `A1`, `zpotrs` on `A2`, then `cpotrf` of
/// `A3` (`A1`'s data again) must leave no NaN in `A3[91]`, the one element
/// upstream checks.
fn potrfBug695() !void {
    const al = std.testing.allocator;
    var a1 = carr(C32, 100, potrf_a1);
    try lapack.potrf(al, ob.io, C32, .u, 10, &a1, 10, null);
    var a2 = carr(C64, 100, potrs_a2);
    var b = carr(C64, 20, potrs_b);
    try lapack.potrs(al, ob.io, C64, .l, 10, 2, &a2, 10, &b, 10);
    var a3 = carr(C32, 100, potrf_a1);
    try lapack.potrf(al, ob.io, C32, .u, 10, &a3, 10, null);
    try ob.expect(!std.math.isNan(a3[91].re) and !std.math.isNan(a3[91].im), "cpotrf bug 695: NaN in A3[91]");
}

test "utest: potrf bug 695, cpotrf and zpotrs" {
    try potrfBug695();
}

/// `potrf.smoketest_trivial`: the 2 x 2 factor of each precision, both
/// triangles, multiplied back (`L Lᴴ` or `Uᴴ U`) within 1e-5 or 1e-12.
fn potrfSmoke() !void {
    const al = std.testing.allocator;
    inline for (.{ f32, f64, C32, C64 }) |T| for ([_]lapack.Uplo{ .l, .u }) |uplo| {
        const complex = comptime ob.isComplex(T);
        const a1: [4]T = if (complex) carr(T, 4, .{ 2, 0, 0.3, 0.1, 0.3, -0.1, 3, 0 }) else arr(T, 4, .{ 2, 0.3, 0.3, 3 });
        var a = a1;
        try lapack.potrf(al, ob.io, T, uplo, 2, &a, 2, null);
        // the other triangle
        a[if (uplo == .l) 2 else 1] = if (complex) .init(0, 0) else 0;
        var b: [4]T = undefined;
        const one: T = if (complex) .init(1, 0) else 1;
        const zero: T = if (complex) .init(0, 0) else 0;
        const t1: blas.Trans = if (uplo == .l) .n else .c;
        const t2: blas.Trans = if (uplo == .l) .c else .n;
        try blas.gemm(al, ob.io, T, t1, t2, 2, 2, 2, one, &a, 2, &a, 2, zero, &b, 2);
        const tol: f64 = if (ob.Float(T) == f32) 1e-5 else 1e-12;
        for (a1, b) |w, g| {
            const err: f64 = if (complex) std.math.hypot(@as(f64, w.re - g.re), @as(f64, w.im - g.im)) else @abs(@as(f64, w - g));
            try ob.expect(err <= tol, "potrf smoketest: difference");
        }
    };
}

test "utest: potrf smoketest, 2 x 2 in each precision" {
    try potrfSmoke();
}

// ---- test_kernel_regress.c ----------------------------------------------------------

/// glibc's `rand` after `srand(seed)` (its `TYPE_3` additive generator,
/// `stdlib/random_r.c`), which `test_kernel_regress.c` fills its matrix
/// from on OpenBLAS's Linux CI; `RAND_MAX` is `2^31 - 1`.
const GlibcRand = struct {
    r: [34]u32,
    i: usize = 0,

    fn init(seed: u32) GlibcRand {
        var g: GlibcRand = .{ .r = undefined };
        var w: [344]u32 = undefined;
        w[0] = if (seed == 0) 1 else seed;
        for (1..31) |i| {
            // `16807 * r mod (2^31 - 1)`, by Schrage's method as glibc does
            const r: i64 = @as(i32, @bitCast(w[i - 1]));
            var word = 16807 * @rem(r, 127773) - 2836 * @divTrunc(r, 127773);
            if (word < 0) word += 2147483647;
            w[i] = @intCast(word);
        }
        for (31..34) |i| w[i] = w[i - 31];
        for (34..344) |i| w[i] = w[i - 31] +% w[i - 3];
        @memcpy(&g.r, w[344 - 34 ..]);
        return g;
    }

    fn next(g: *GlibcRand) u32 {
        // the last 34 values, a ring: r[k] = r[k - 31] + r[k - 3]
        const k = g.i % 34;
        const v = g.r[(k + 34 - 31) % 34] +% g.r[(k + 34 - 3) % 34];
        g.r[k] = v;
        g.i += 1;
        return v >> 1;
    }
};

/// `kernel_regress.skx_avx`: `dgesvd` with jobs `A` of a random 100 x 100
/// in [0, 10] (`srand(0)`), then `‖U Σ Vᵀ - A‖₂` (Frobenius) within 1e-10.
/// The row-major buffer read column-major is `Aᵀ`, whose decomposition the
/// same check covers: the residual's norm is the transpose's.
fn kernelRegress() !void {
    const al = std.testing.allocator;
    const n = 100;
    const m = try al.alloc(f64, n * n);
    defer al.free(m);
    const tmp = try al.alloc(f64, n * n);
    defer al.free(tmp);
    const u = try al.alloc(f64, n * n);
    defer al.free(u);
    const vt = try al.alloc(f64, n * n);
    defer al.free(vt);
    const x = try al.alloc(f64, n * n);
    defer al.free(x);
    var s: [n]f64 = undefined;
    var g: GlibcRand = .init(0);
    // glibc's first two values after `srand(0)`, which is `srand(1)`
    var g2 = g;
    try ob.expect(g2.next() == 1804289383 and g2.next() == 846930886, "kernel_regress: glibc rand");
    for (m, tmp) |*a, *t| {
        a.* = (@as(f64, @floatFromInt(g.next())) + 0.0) / 2147483647.0 * 10;
        t.* = a.*;
    }
    try lapack.gesvd(al, ob.io, f64, .a, .a, n, n, m.ptr, n, &s, u.ptr, n, vt.ptr, n);
    for (0..n) |j| for (0..n) |i| {
        u[i + j * n] *= s[j];
    };
    try blas.gemm(al, ob.io, f64, .n, .n, n, n, n, 1, u.ptr, n, vt.ptr, n, 0, x.ptr, n);
    for (x, tmp) |*a, t| a.* -= t;
    const norm = try blas.nrm2(ob.io, f64, n * n, x.ptr, 1);
    try near(0, norm, 1e-10, "kernel_regress skx_avx: |U S Vt - A|");
}

test "utest: kernel_regress skx_avx, dgesvd" {
    try kernelRegress();
}

// ---- test_potrs.c's data -------------------------------------------------------------
// Transcribed mechanically from `utest/test_potrs.c` (OpenBLAS v0.3.34-317), as
// real and imaginary parts in turn; re-derive rather than edit. Its `A3` is `A1`.

/// `A1`.
const potrf_a1 = [_]f64{
    5.8525753,   0.0,         -0.79540455,  -0.7066077,   0.98274714,  -1.3824869,  2.619998,     -1.8532984,
    -1.8306153,  1.2336911,   0.32275113,   -0.015575029, 2.1968813,   -1.0640624,  0.27894387,   -0.97911835,
    3.0476584,   -0.18548489, 0.3842994,    -0.7050991,   -0.79540455, 0.7066077,   8.313246,     0.0,
    -1.8076122,  0.8882447,   0.47806996,   -0.48494184,  0.5096429,   0.5395974,   -0.7285097,   0.10360408,
    -1.1760061,  2.7146957,   -0.4271084,   -0.042899966, -1.7228563,  -2.8335886,  1.8942566,    -0.6389735,
    0.98274714,  1.3824869,   -1.8076122,   -0.8882447,   9.367975,    0.0,         -0.1838578,   -0.6468568,
    -1.8338387,  -0.7064959,  0.041852742,  0.6556877,    2.5673025,   -1.9732997,  -1.1148382,   0.15693812,
    2.4704504,   1.0389464,   1.0858271,    1.298006,     2.619998,    1.8532984,   0.47806996,   0.48494184,
    -0.1838578,  0.6468568,   3.1117508,    0.0,          -1.956626,   -0.22825956, 0.07081801,   0.31801307,
    0.3698375,   0.5400855,   0.80686307,   -1.5315914,   1.5649154,   1.6229297,   -0.112077385, -1.2014246,
    -1.8306153,  -1.2336911,  0.5096429,    -0.5395974,   -1.8338387,  0.7064959,   -1.956626,    0.22825956,
    3.6439795,   0.0,         -0.2594722,   -0.48786148,  -0.47636223, 0.27821827,  -0.61608654,  2.01858,
    -2.7767487,  -1.7693765,  0.048102796,  0.9741874,    0.32275113,  0.015575029, -0.7285097,   -0.10360408,
    0.041852742, -0.6556877,  0.07081801,   -0.31801307,  -0.2594722,  0.48786148,  3.624376,     0.0,
    -1.6697118,  -0.4017511,  -1.4397877,   0.7550918,    -0.31456697, 1.0403451,   -0.31978557,  -0.13701046,
    2.1968813,   1.0640624,   -1.1760061,   -2.7146957,   2.5673025,   1.9732997,   0.3698375,    -0.5400855,
    -0.47636223, -0.27821827, -1.6697118,   0.4017511,    6.8273163,   0.0,         -0.10051322,  -0.24303961,
    1.4415971,   -0.29750675, 1.221786,     0.85654986,   0.27894387,  0.97911835,  -0.4271084,   0.042899966,
    -1.1148382,  -0.15693812, 0.80686307,   1.5315914,    -0.61608654, -2.01858,    -1.4397877,   -0.7550918,
    -0.10051322, 0.24303961,  3.4057708,    0.0,          -0.5856801,  1.0203559,   0.7103452,    -0.8422135,
    3.0476584,   0.18548489,  -1.7228563,   2.8335886,    2.4704504,   -1.0389464,  1.5649154,    -1.6229297,
    -2.7767487,  1.7693765,   -0.31456697,  -1.0403451,   1.4415971,   0.29750675,  -0.5856801,   -1.0203559,
    7.005772,    0.0,         -0.9617417,   1.2486815,    0.3842994,   0.7050991,   1.8942566,    0.6389735,
    1.0858271,   -1.298006,   -0.112077385, 1.2014246,    0.048102796, -0.9741874,  -0.31978557,  0.13701046,
    1.221786,    -0.85654986, 0.7103452,    0.8422135,    -0.9617417,  -1.2486815,  3.4629636,    0.0,
};

/// `A2`.
const potrs_a2 = [_]f64{
    3.0607147216796875,   0.0,                   -0.5905849933624268,  -0.29020825028419495,  0.321084201335907,     0.45168760418891907,   0.8387917876243591,   -0.644718587398529,
    -0.3642411530017853,  0.051274992525577545,  0.8071482181549072,   0.33944568037986755,   0.013674172572791576,  0.21422699093818665,   0.35476258397102356,  0.42408594489097595,
    -0.5991537570953369,  -0.23082709312438965,  -0.0600702166557312,  -0.2113417387008667,   -0.7954045534133911,   0.7066076993942261,    2.807175397872925,    0.0,
    -0.1691000759601593,  0.313548743724823,     -0.30911174416542053, 0.7447023987770081,    -0.22347848117351532,  0.03316075727343559,   -0.4088296890258789,  -1.0214389562606812,
    -0.2344931811094284,  0.08056317269802094,   0.793269693851471,    -0.17507623136043549,  0.03163455054163933,   0.20559945702552795,   0.13581633567810059,  -0.2110036462545395,
    0.9827471375465393,   1.3824869394302368,    -1.8076121807098389,  -0.8882446885108948,   2.3277781009674072,    0.0,                   0.830405056476593,    -0.19296252727508545,
    0.1394239068031311,   -0.5260677933692932,   1.239942193031311,    -0.09915469586849213,  0.06731037050485611,   -0.059320636093616486, 0.11507681757211685,  -0.1984301060438156,
    -0.6843825578689575,  0.4647614359855652,    1.213119387626648,    -0.7757048010826111,   2.619997978210449,     1.8532984256744385,    0.4780699610710144,   0.48494184017181396,
    -0.18385779857635498, 0.6468567848205566,    2.0811400413513184,   0.0,                   -0.035075582563877106, 0.09732913225889206,   0.27337002754211426,  -0.9032229781150818,
    -0.8374675512313843,  0.0479498989880085,    0.6916252374649048,   0.45711082220077515,   0.1883818507194519,    0.06482727080583572,   -0.32384994626045227, 0.05857187137007713,
    -1.8306152820587158,  -1.2336910963058472,   0.5096428990364075,   -0.5395973920822144,   -1.833838701248169,    0.7064958810806274,    -1.956626057624817,   0.22825956344604492,
    1.706615924835205,    0.0,                   -0.2895336151123047,  0.17579378187656403,   -0.923172116279602,    -0.4530014097690582,   0.5040621757507324,   -0.37026339769363403,
    -0.2824432849884033,  -1.0374568700790405,   0.1399831622838974,   0.4977008104324341,    0.32275113463401794,   0.015575028955936432,  -0.7285097241401672,  -0.10360407829284668,
    0.041852742433547974, -0.655687689781189,    0.07081800699234009,  -0.318013072013855,    -0.25947219133377075,  0.4878614842891693,    1.5735365152359009,   0.0,
    -0.2647853195667267,  -0.26654252409935,     -0.6190430521965027,  -0.24699924886226654,  -0.6288471221923828,   0.48154571652412415,   0.02446540631353855,  -0.2611822783946991,
    2.1968812942504883,   1.0640623569488525,    -1.1760060787200928,  -2.714695692062378,    2.5673024654388428,    1.9732997417449951,    0.3698374927043915,   -0.54008549451828,
    -0.4763622283935547,  -0.27821826934814453,  -1.6697118282318115,  0.4017511010169983,    1.2674795389175415,    0.0,                   0.3079095482826233,   -0.07258892804384232,
    -0.5929520130157471,  -0.038360968232154846, 0.04388086497783661,  -0.025549031794071198, 0.27894386649131775,   0.9791183471679688,    -0.42710840702056885, 0.0428999662399292,
    -1.1148382425308228,  -0.1569381207227707,   0.8068630695343018,   1.5315914154052734,    -0.6160865426063538,   -2.0185799598693848,   -1.439787745475769,   -0.7550917863845825,
    -0.10051321983337402, 0.24303960800170898,   0.9066106081008911,   0.0,                   0.05315789580345154,   -0.06136537343263626,  -0.21304509043693542, 0.6494344472885132,
    3.0476584434509277,   0.1854848861694336,    -1.7228562831878662,  2.8335886001586914,    2.4704504013061523,    -1.0389463901519775,   1.564915418624878,    -1.6229296922683716,
    -2.7767486572265625,  1.769376516342163,     -0.314566969871521,   -1.0403450727462769,   1.4415971040725708,    0.29750674962997437,   -0.5856801271438599,  -1.0203559398651123,
    0.5668219923973083,   0.0,                   0.033351436257362366, -0.07832501083612442,  0.3842993974685669,    0.7050991058349609,    1.894256591796875,    0.6389734745025635,
    1.085827112197876,    -1.2980060577392578,   -0.11207738518714905, 1.2014245986938477,    0.04810279607772827,   -0.9741873741149902,   -0.31978556513786316, 0.13701045513153076,
    1.2217860221862793,   -0.856549859046936,    0.7103452086448669,   0.84221351146698,      -0.9617416858673096,   -1.2486815452575684,   0.0756804421544075,   0.0,
};

/// `B`.
const potrs_b = [_]f64{
    -0.21782716937787788, -0.9222220085490986, -0.7620356655676837, 0.15533508334193666,  -0.905011814118756,    0.2847570854574069,  -0.3451346708401685, 1.076948486041297,
    0.25336108035924787,  0.975317836492159,   0.11192755545114,    -0.1603741874112385,  -0.20604111555491242,  0.10570814584017311, -1.0568488936791578, -0.06025820467086475,
    -0.6650468984506477,  -0.5000967284800251, -1.0509472322215125, 0.5022165705328413,   -0.727775859267237,    0.50638268521728,    0.39947219167701153, -0.4576746001199889,
    -0.7122162951294634,  -0.630289556702497,  0.9870834574024372,  -0.2825689605519449,  0.0628393808469436,    -0.1253397353973715, 0.8439562576196216,  1.0850814110398734,
    0.562377322638969,    -0.2578030745663871, 0.12696236014017806, -0.09853584666755086, -0.023682508769195098, 0.18093440285319276, -0.7264975746431271, 0.31670415674097235,
};
