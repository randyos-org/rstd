//! Tests: products, reductions, transposes, shapes.

const sysinfo = @import("builtin");
const std = @import("std");
const Allocator = std.mem.Allocator;
const Io = std.Io;
const blas = @import("../../../math/kernels/blas/__root__.zig");
const ext = @import("../../../math/kernels/blas/ext/interface.zig");
const l3 = @import("../../../math/kernels/blas/blas3/core.zig");
const l2 = @import("../../../math/kernels/blas/blas2/core.zig");
const l1 = @import("../../../math/kernels/blas/blas1/interface.zig");
const lapack = @import("../../../math/kernels/lapack/__root__.zig");
const vector = @import("../../vector.zig");
const matrix = @import("../../matrix.zig");
const ops = @import("../../../math/kernels/lapack/ops.zig");

const View = matrix.View;
const ConstView = matrix.ConstView;
const Mat = matrix.Mat;
const isComplex = vector.isComplex;
const Real = vector.Real;
const eMul = matrix.eMul;
const eMulAdd = matrix.eMulAdd;

const zigstd = @import("../../../zigstd.zig");
const bitwise = @import("../../../math/kernels/blas/_test/bitwise.zig");
const testing = std.testing;
const sameBitsOrX87 = bitwise.sameBitsOrX87;
const C64 = std.math.Complex(f64);
const test_io = Io.failing;

const factor = @import("../factor.zig");
const chol = factor.chol;
const inverse = @import("../inverse.zig");
const det = inverse.det;
const inv = inverse.inv;
const invInto = inverse.invInto;
const solve = inverse.solve;
const solveInto = inverse.solveInto;
const products = @import("../products.zig");
const ctranspose = products.ctranspose;
const mul = products.mul;
const mulInto = products.mulInto;
const transpose = products.transpose;
const transposeInto = products.transposeInto;
const reductions = @import("../reductions.zig");
const axpy = reductions.axpy;
const dot = reductions.dot;
const max = reductions.max;
const min = reductions.min;
const norm = reductions.norm;
const scale = reductions.scale;
const sum = reductions.sum;
const util = @import("../util.zig");
const copyOf = util.copyOf;

pub fn randomElem(comptime T: type, rnd: std.Random) T {
    if (comptime isComplex(T)) return .init(rnd.float(f64) * 2 - 1, rnd.float(f64) * 2 - 1);
    return std.math.ldexp(rnd.float(T) * 2 - 1, rnd.intRangeAtMost(i32, -4, 4));
}

pub fn fill(comptime T: type, rnd: std.Random, s: []T) void {
    for (s) |*x| x.* = randomElem(T, rnd);
}

/// `mulInto` on a `Mat` against `linalg.matrix.mul` on `StaticMat`s of the
/// same numbers, `R` by `K` times `K` by `C`.
pub fn expectStaticBits(comptime T: type, comptime r: usize, comptime k: usize, comptime c: usize, rnd: std.Random) !void {
    const al = testing.allocator;
    for (0..20) |_| {
        var a: matrix.StaticMat(r, k, T) = undefined;
        var b: matrix.StaticMat(k, c, T) = undefined;
        fill(T, rnd, a.view().data);
        fill(T, rnd, b.view().data);
        const want = matrix.mul(a, b);
        const got = try mul(al, test_io, &a, &b);
        defer got.deinit(al);
        var w = want;
        try testing.expect(sameBitsOrX87(T, w.view().data, got.data));
        // a column: `M v`
        const v: [k]T = b.data[0];
        var col: [r]T = undefined;
        try mulInto(al, test_io, &col, &a, &v);
        try testing.expect(sameBitsOrX87(T, &col, &@as([r]T, w.data[0])));
    }
}

test "up to the inline limit, a Mat multiplies to a StaticMat's bits" {
    var prng: std.Random.DefaultPrng = .init(0x7e57);
    const shapes = .{ .{ 1, 1, 1 }, .{ 2, 3, 4 }, .{ 4, 4, 4 }, .{ 3, 8, 5 }, .{ 8, 8, 8 }, .{ matrix.inline_limit, matrix.inline_limit, matrix.inline_limit } };
    inline for (.{ f32, f64, C64 }) |T| inline for (shapes) |s| try expectStaticBits(T, s[0], s[1], s[2], prng.random());
}

/// `Σ_p a(i, p) b(p, j)` in `f64` (or its complex), against `got` within
/// `k` roundings of `T`.
pub fn expectNearProduct(comptime T: type, a: ConstView(T), b: ConstView(T), got: ConstView(T)) !void {
    const eps: f64 = std.math.floatEps(Real(T));
    for (0..got.cols) |j| for (0..got.rows) |i| {
        var re: f64 = 0;
        var im: f64 = 0;
        var mag: f64 = 0;
        for (0..a.cols) |p| {
            const x = a.at(i, p);
            const y = b.at(p, j);
            if (comptime isComplex(T)) {
                re += x.re * y.re - x.im * y.im;
                im += x.re * y.im + x.im * y.re;
                mag += (@abs(x.re) + @abs(x.im)) * (@abs(y.re) + @abs(y.im));
            } else {
                re += @as(f64, x) * y;
                mag += @abs(@as(f64, x) * y);
            }
        }
        const g = got.at(i, j);
        const tol = 2 * @as(f64, @floatFromInt(a.cols + 1)) * eps * mag;
        if (comptime isComplex(T)) {
            try testing.expect(@abs(g.re - re) <= tol and @abs(g.im - im) <= tol);
        } else try testing.expect(@abs(@as(f64, g) - re) <= tol);
    };
}

pub fn expectBlasProduct(comptime T: type, m: usize, n: usize, k: usize, rnd: std.Random) !void {
    const al = testing.allocator;
    const a = try Mat(T).init(al, m, k);
    defer a.deinit(al);
    const b = try Mat(T).init(al, k, n);
    defer b.deinit(al);
    fill(T, rnd, a.data);
    fill(T, rnd, b.data);
    // undefined storage may hold NaNs: `beta = 0` must not read them
    const out = try Mat(T).init(al, m, n);
    defer out.deinit(al);
    const nan = std.math.nan(Real(T));
    @memset(out.data, if (comptime isComplex(T)) T.init(nan, nan) else nan);
    try mulInto(al, test_io, out, a, b);
    try expectNearProduct(T, a.constView(), b.constView(), out.constView());
}

test "past the inline limit, mul is the BLAS's gemm or gemv" {
    var prng: std.Random.DefaultPrng = .init(0xb1a5);
    inline for (.{ f32, f64, C64 }) |T| {
        try expectBlasProduct(T, matrix.inline_limit + 1, 7, 5, prng.random());
        try expectBlasProduct(T, 3, 4, 20, prng.random());
        try expectBlasProduct(T, 40, 33, 17, prng.random());
        try expectBlasProduct(T, 30, 1, 12, prng.random());
    }
}

pub fn expectVectorOps(comptime T: type, n: usize, rnd: std.Random) !void {
    const al = testing.allocator;
    // `x` is a row of a 3-row matrix, stride 3; `y` a column
    const xs = try al.alloc(T, 3 * n);
    defer al.free(xs);
    fill(T, rnd, xs);
    const xv = try ConstView(T).init(xs, 1, n, 3);
    const y = try al.alloc(T, n);
    defer al.free(y);
    fill(T, rnd, y);
    var d_re: f64 = 0;
    var d_im: f64 = 0;
    var s_re: f64 = 0;
    var s_im: f64 = 0;
    var sq: f64 = 0;
    var mag: f64 = 0;
    for (0..n) |i| {
        const x = xv.at(0, i);
        if (comptime isComplex(T)) {
            d_re += x.re * y[i].re + x.im * y[i].im;
            d_im += x.re * y[i].im - x.im * y[i].re;
            s_re += x.re;
            s_im += x.im;
            sq += x.re * x.re + x.im * x.im;
            mag += (@abs(x.re) + @abs(x.im)) * (@abs(y[i].re) + @abs(y[i].im));
        } else {
            d_re += @as(f64, x) * y[i];
            s_re += x;
            sq += @as(f64, x) * x;
            mag += @abs(@as(f64, x) * y[i]);
        }
    }
    const eps: f64 = std.math.floatEps(Real(T));
    const tol = 2 * @as(f64, @floatFromInt(n + 1)) * eps;
    const d = try dot(test_io, xv, y);
    const s = try sum(test_io, xv);
    const nrm: f64 = try norm(test_io, xv);
    try testing.expect(@abs(nrm - @sqrt(sq)) <= tol * @sqrt(sq));
    if (comptime isComplex(T)) {
        try testing.expect(@abs(d.re - d_re) <= tol * mag and @abs(d.im - d_im) <= tol * mag);
        try testing.expect(@abs(s.re - s_re) <= tol * @as(f64, @floatFromInt(n)) and @abs(s.im - s_im) <= tol * @as(f64, @floatFromInt(n)));
    } else {
        try testing.expect(@abs(@as(f64, d) - d_re) <= tol * mag);
        try testing.expect(@abs(@as(f64, s) - s_re) <= tol * @as(f64, @floatFromInt(n)));
        var hi = xv.at(0, 0);
        var lo = hi;
        for (1..n) |i| {
            hi = @max(hi, xv.at(0, i));
            lo = @min(lo, xv.at(0, i));
        }
        try testing.expectEqual(hi, try max(test_io, xv));
        try testing.expectEqual(lo, try min(test_io, xv));
    }
}

test "dot, norm, sum, max and min, inline and on the BLAS, on rows and columns" {
    var prng: std.Random.DefaultPrng = .init(0xd07);
    inline for (.{ f32, f64, C64 }) |T| for ([_]usize{ 1, 3, 8, matrix.inline_limit, matrix.inline_limit + 1, 100 }) |n| try expectVectorOps(T, n, prng.random());
}

test "up to the inline limit, vector reductions keep linalg.vector's bits" {
    var prng: std.Random.DefaultPrng = .init(0x5a5e);
    var a: [6]f64 = undefined;
    var b: [6]f64 = undefined;
    fill(f64, prng.random(), &a);
    fill(f64, prng.random(), &b);
    try testing.expectEqual(vector.dot(a, b), try dot(test_io, &a, &b));
    try testing.expectEqual(vector.length(a), try norm(test_io, &a));
    try testing.expectEqual(vector.sum(a), try sum(test_io, &a));
    var z: [5]C64 = undefined;
    fill(C64, prng.random(), &z);
    try testing.expectEqual(vector.dot(z, z), try dot(test_io, &z, &z));
}

test "the float-slice reductions on both paths: max and min return the first NaN, unlike @reduce" {
    inline for (.{ f32, f64 }) |F| {
        const U = @Int(.unsigned, @bitSizeOf(F));
        // two quiet NaNs told apart by their payloads
        const nan1: F = @bitCast(@as(U, @bitCast(std.math.nan(F))) | 1);
        const nan2: F = @bitCast(@as(U, @bitCast(std.math.nan(F))) | 2);
        for ([_]usize{ 5, 100 }) |n| {
            var buf: [100]F = undefined;
            for (buf[0..n], 0..) |*v, i| v.* = @floatFromInt(i % 7);
            buf[1] = nan1;
            buf[n - 1] = nan2;
            const x: []const F = buf[0..n];
            const want: U = @bitCast(nan1);
            try testing.expectEqual(want, @as(U, @bitCast(try max(test_io, x))));
            try testing.expectEqual(want, @as(U, @bitCast(try min(test_io, x))));
            try testing.expect(std.math.isNan(try sum(test_io, x)));
            // no NaN: the largest, the smallest and the sum of the slice
            buf[1] = 0;
            buf[n - 1] = 0;
            var total: F = 0;
            var hi: F = x[0];
            var lo: F = x[0];
            for (x) |v| {
                total += v;
                hi = @max(hi, v);
                lo = @min(lo, v);
            }
            try testing.expectEqual(hi, try max(test_io, x));
            try testing.expectEqual(lo, try min(test_io, x));
            try testing.expectApproxEqRel(total, try sum(test_io, x), 1e-6);
        }
    }
}

test "transposes, axpy and scale" {
    const al = testing.allocator;
    var a = matrix.fromRows(matrix.StaticMat(2, 3, f64), .{ .{ 1, 2, 3 }, .{ 4, 5, 6 } });
    const t = try transpose(al, test_io, &a);
    defer t.deinit(al);
    try testing.expectEqualSlices(f64, &.{ 1, 2, 3, 4, 5, 6 }, t.data);
    var too_small: [4]f64 = undefined;
    try testing.expectError(error.DimensionMismatch, transposeInto(test_io, &too_small, &a));

    var z: matrix.StaticMat(1, 2, C64) = .{ .data = .{ .{.init(1, 2)}, .{.init(3, -4)} } };
    const zh = try ctranspose(al, test_io, &z);
    defer zh.deinit(al);
    try testing.expectEqual(C64.init(1, -2), zh.data[0]);
    try testing.expectEqual(C64.init(3, 4), zh.data[1]);

    var y = [_]f64{ 1, 1, 1 };
    try axpy(test_io, 2, &[_]f64{ 1, 2, 3 }, &y);
    try testing.expectEqualSlices(f64, &.{ 3, 5, 7 }, &y);
    try testing.expectError(error.DimensionMismatch, axpy(test_io, 2, &[_]f64{ 1, 2 }, &y));
    scale(-1, &a);
    try testing.expectEqual(@as(f64, -6), a.data[2][1]);
    try testing.expectEqual(matrix.scale(matrix.fromCols(matrix.M2f32, .{ .{ 0.1, 0.2 }, .{ 0.3, 0.4 } }), 0.7).data[1][0], blk: {
        var s = matrix.fromCols(matrix.M2f32, .{ .{ 0.1, 0.2 }, .{ 0.3, 0.4 } });
        scale(0.7, &s);
        break :blk s.data[1][0];
    });
}

test "shape errors" {
    const al = testing.allocator;
    var a: matrix.StaticMat(2, 3, f64) = .{ .data = @splat(@splat(1)) };
    try testing.expectError(error.DimensionMismatch, dot(test_io, &a, &a));
    try testing.expectError(error.DimensionMismatch, dot(test_io, &[_]f64{1}, &[_]f64{ 1, 2 }));
    try testing.expectError(error.DimensionMismatch, mul(al, test_io, &a, &a));
    var out: matrix.M3f64 = undefined;
    try testing.expectError(error.DimensionMismatch, mulInto(al, test_io, &out, &a, &a));
    try testing.expectEqual(@as(f64, 0), try sum(test_io, &[_]f64{}));
    try testing.expectError(error.NotSquare, det(al, test_io, &a));
    try testing.expectError(error.NotSquare, inv(al, test_io, &a));
    try testing.expectError(error.NotSquare, chol(al, test_io, &a));
    try testing.expectError(error.DimensionMismatch, solve(al, test_io, &out, &a));
    try testing.expectError(error.DimensionMismatch, invInto(al, test_io, &a, &out));
    var x: [2]f64 = undefined;
    try testing.expectError(error.DimensionMismatch, solveInto(al, test_io, &x, &out, &[_]f64{ 1, 2, 3 }));
}

/// `det`, `inv` and `solve` on a `Mat` against `linalg.matrix`'s on the
/// `StaticMat` of the same numbers.
pub fn expectFixedBits(comptime T: type, comptime c: usize, rnd: std.Random) !void {
    const al = testing.allocator;
    for (0..10) |_| {
        var a: matrix.StaticMat(c, c, T) = undefined;
        var b: matrix.StaticMat(c, 3, T) = undefined;
        fill(T, rnd, a.view().data);
        fill(T, rnd, b.view().data);
        const m = try copyOf(T, al, a.constView());
        defer m.deinit(al);
        try testing.expect(sameBitsOrX87(T, &[_]T{matrix.det(a)}, &[_]T{try det(al, test_io, m)}));
        var want_inv = try matrix.inv(a);
        const got_inv = try inv(al, test_io, m);
        defer got_inv.deinit(al);
        try testing.expect(sameBitsOrX87(T, want_inv.view().data, got_inv.data));
        var want_x = try matrix.solve(a, b);
        const got_x = try solve(al, test_io, m, &b);
        defer got_x.deinit(al);
        try testing.expect(sameBitsOrX87(T, want_x.view().data, got_x.data));
    }
}

test "up to the inline limit, det, inv and solve keep a StaticMat's bits" {
    var prng: std.Random.DefaultPrng = .init(0x50f7);
    inline for (.{ f32, f64, C64 }) |T| inline for (.{ 1, 2, 3, 4, 6, 8, matrix.inline_limit }) |c| try expectFixedBits(T, c, prng.random());
}
