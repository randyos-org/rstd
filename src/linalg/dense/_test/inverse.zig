//! Tests: det, inv and solve.

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

const test_products = @import("products.zig");
const fill = test_products.fill;
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

pub fn absOf(x: anytype) f64 {
    return if (comptime isComplex(@TypeOf(x))) @abs(x.re) + @abs(x.im) else @abs(@as(f64, x));
}

pub fn maxAbs(comptime T: type, a: ConstView(T)) f64 {
    var r: f64 = 0;
    for (0..a.cols) |j| for (a.col(j)) |x| {
        r = @max(r, absOf(x));
    };
    return r;
}

/// `max |a(i, j) - b(i, j)|`, `|re| + |im|` for complex.
pub fn maxDiff(comptime T: type, a: ConstView(T), b: ConstView(T)) f64 {
    var r: f64 = 0;
    for (0..a.cols) |j| for (0..a.rows) |i| {
        const x = a.at(i, j);
        const y = b.at(i, j);
        const d = if (comptime isComplex(T)) @abs(x.re - y.re) + @abs(x.im - y.im) else @abs(@as(f64, x) - y);
        r = @max(r, d);
    };
    return r;
}

/// 100 `n` roundings of `T`.
pub fn tolOf(comptime T: type, n: usize) f64 {
    return 100 * @as(f64, @floatFromInt(@max(n, 1))) * std.math.floatEps(Real(T));
}

pub fn randomMat(comptime T: type, rnd: std.Random, m: usize, n: usize) !Mat(T) {
    const r = try Mat(T).init(testing.allocator, m, n);
    fill(T, rnd, r.data);
    return r;
}

/// `a(i, i) += s`.
pub fn shiftDiag(comptime T: type, a: Mat(T), s: f64) void {
    for (0..@min(a.rows, a.cols)) |i| {
        const d = &a.data[i + i * a.rows];
        d.* = if (comptime isComplex(T)) T.init(d.re + s, d.im) else d.* + @as(T, @floatCast(s));
    }
}

/// `Mᴴ M + n I`, exactly Hermitian.
pub fn hermitianPd(comptime T: type, rnd: std.Random, n: usize) !Mat(T) {
    const al = testing.allocator;
    const m = try randomMat(T, rnd, n, n);
    defer m.deinit(al);
    const mh = try ctranspose(al, test_io, m);
    defer mh.deinit(al);
    const a = try mul(al, test_io, mh, m);
    for (0..n) |j| {
        for (j + 1..n) |i| a.data[i + j * n] = ops.conj(T, a.data[j + i * n]);
        if (comptime isComplex(T)) a.data[j + j * n].im = 0;
    }
    shiftDiag(T, a, @floatFromInt(n));
    return a;
}

/// `‖a x - b‖ <= tol (n ‖a‖ ‖x‖ + ‖b‖)`, max norms.
pub fn expectSolves(comptime T: type, a: ConstView(T), x: ConstView(T), b: ConstView(T)) !void {
    const al = testing.allocator;
    const ax = try mul(al, test_io, a, x);
    defer ax.deinit(al);
    const n: f64 = @floatFromInt(a.cols);
    const bound = tolOf(T, a.cols) * (n * maxAbs(T, a) * maxAbs(T, x) + maxAbs(T, b));
    try testing.expect(maxDiff(T, ax.constView(), b) <= bound);
}

pub fn expectSolveRoutes(comptime T: type, n: usize, rnd: std.Random) !void {
    const al = testing.allocator;
    const b = try randomMat(T, rnd, n, 3);
    defer b.deinit(al);
    const general = try randomMat(T, rnd, n, n);
    defer general.deinit(al);
    shiftDiag(T, general, 2);
    const upper = try randomMat(T, rnd, n, n);
    defer upper.deinit(al);
    const lower = try randomMat(T, rnd, n, n);
    defer lower.deinit(al);
    for (0..n) |j| for (0..n) |i| {
        if (i > j) upper.data[i + j * n] = ops.zero(T);
        if (i < j) lower.data[i + j * n] = ops.zero(T);
    };
    shiftDiag(T, upper, 4);
    shiftDiag(T, lower, 4);
    const pd = try hermitianPd(T, rnd, n);
    defer pd.deinit(al);
    // Hermitian, positive diagonal, indefinite: `potrf` fails, LU solves
    const indefinite = try Mat(T).zeros(al, n, n);
    defer indefinite.deinit(al);
    shiftDiag(T, indefinite, 1);
    indefinite.data[1] = if (comptime isComplex(T)) T.init(2, 0) else 2;
    indefinite.data[n] = indefinite.data[1];
    try testing.expectError(error.NotPositiveDefinite, chol(al, test_io, indefinite));
    for ([_]Mat(T){ general, upper, lower, pd, indefinite }) |a| {
        const x = try solve(al, test_io, a, b);
        defer x.deinit(al);
        try expectSolves(T, a.constView(), x.constView(), b.constView());
        // one right-hand side: `trsv` and `gesv`'s single column
        const x1 = try Mat(T).init(al, n, 1);
        defer x1.deinit(al);
        try solveInto(al, test_io, x1, a, b.constView().col(0));
        try expectSolves(T, a.constView(), x1.constView(), (try b.constView().sub(0, 0, n, 1)));
    }
}

test "past the inline limit, solve routes by structure, each to a solution" {
    var prng: std.Random.DefaultPrng = .init(0x501e);
    inline for (.{ f32, f64, C64 }) |T| for ([_]usize{ matrix.inline_limit + 1, 30 }) |n| try expectSolveRoutes(T, n, prng.random());
    // a small square with many right-hand sides is past the limit too
    const al = testing.allocator;
    const a = try randomMat(f64, prng.random(), 3, 3);
    defer a.deinit(al);
    shiftDiag(f64, a, 2);
    const b = try randomMat(f64, prng.random(), 3, 20);
    defer b.deinit(al);
    const x = try solve(al, test_io, a, b);
    defer x.deinit(al);
    try expectSolves(f64, a.constView(), x.constView(), b.constView());
}

test "solve and inv: error.Singular leaves out untouched" {
    const al = testing.allocator;
    inline for (.{ 3, matrix.inline_limit + 1 }) |n| {
        // ones: Hermitian, not positive definite, LU's second pivot exactly 0
        const ones = try Mat(f64).init(al, n, n);
        defer ones.deinit(al);
        @memset(ones.data, 1);
        // zero: triangular, `trtrs`
        const zero = try Mat(f64).zeros(al, n, n);
        defer zero.deinit(al);
        for ([_]Mat(f64){ ones, zero }) |a| {
            var out: [n * n]f64 = @splat(7);
            const v = try View(f64).init(&out, n, n, n);
            try testing.expectError(error.Singular, solveInto(al, test_io, v, a, a));
            try testing.expectError(error.Singular, invInto(al, test_io, v, a));
            try testing.expect(std.mem.allEqual(f64, &out, 7));
            try testing.expectEqual(@as(f64, 0), try det(al, test_io, a));
        }
    }
}

pub fn expectLeastSquares(comptime T: type, m: usize, n: usize, rnd: std.Random) !void {
    const al = testing.allocator;
    const a = try randomMat(T, rnd, m, n);
    defer a.deinit(al);
    const x0 = try randomMat(T, rnd, n, 2);
    defer x0.deinit(al);
    const b = try mul(al, test_io, a, x0);
    defer b.deinit(al);
    const x = try solve(al, test_io, a, b);
    defer x.deinit(al);
    try expectSolves(T, a.constView(), x.constView(), b.constView());
    if (m > n) {
        // full column rank: the only solution
        try testing.expect(maxDiff(T, x.constView(), x0.constView()) <= tolOf(T, m) * 10 * maxAbs(T, x0.constView()));
        return;
    }
    // least norm: `x = aᴴ (a aᴴ)⁻¹ b`
    const ah = try ctranspose(al, test_io, a);
    defer ah.deinit(al);
    const g = try mul(al, test_io, a, ah);
    defer g.deinit(al);
    const y = try solve(al, test_io, g, b);
    defer y.deinit(al);
    const want = try mul(al, test_io, ah, y);
    defer want.deinit(al);
    try testing.expect(maxDiff(T, x.constView(), want.constView()) <= tolOf(T, n) * 10 * maxAbs(T, want.constView()));
}

test "rectangular solve is gelsy's least-norm least squares" {
    var prng: std.Random.DefaultPrng = .init(0x9e15);
    inline for (.{ f32, f64, C64 }) |T| {
        try expectLeastSquares(T, 15, 9, prng.random());
        try expectLeastSquares(T, 4, 12, prng.random());
        try expectLeastSquares(T, 3, 2, prng.random());
    }
    const al = testing.allocator;
    const empty = try Mat(f64).init(al, 0, 4);
    defer empty.deinit(al);
    const x = try solve(al, test_io, empty, try empty.constView().sub(0, 0, 0, 2));
    defer x.deinit(al);
    try testing.expect(x.rows == 4 and x.cols == 2 and std.mem.allEqual(f64, x.data, 0));
}

test "det past the inline limit: exact cases, scaled past the range, against inv" {
    const al = testing.allocator;
    // one past the limit
    const n = matrix.inline_limit + 1;
    // `diag(1e200, 1e200, -3e-100, 1, ...)` (design.md, M1's E1)
    const e1 = try Mat(f64).zeros(al, n, n);
    defer e1.deinit(al);
    shiftDiag(f64, e1, 1);
    e1.data[0] = 1e200;
    e1.data[n + 1] = 1e200;
    e1.data[2 * (n + 1)] = -3e-100;
    const d = try det(al, test_io, e1);
    try testing.expect(@abs(d + 3e300) <= 1e-14 * 3e300);
    // `2 I` and a row swap, exact
    const two = try Mat(f64).zeros(al, n, n);
    defer two.deinit(al);
    shiftDiag(f64, two, 2);
    try testing.expectEqual(@as(f64, 1 << n), try det(al, test_io, two));
    const swap = try Mat(f64).zeros(al, n, n);
    defer swap.deinit(al);
    shiftDiag(f64, swap, 1);
    swap.data[0] = 0;
    swap.data[n + 1] = 0;
    swap.data[1] = 1;
    swap.data[n] = 1;
    try testing.expectEqual(@as(f64, -1), try det(al, test_io, swap));
    // `det(i I) = i^n`
    var zi: [n * n]C64 = @splat(.init(0, 0));
    for (0..n) |i| zi[i * (n + 1)] = .init(0, 1);
    const i_pow = [4]C64{ .init(1, 0), .init(0, 1), .init(-1, 0), .init(0, -1) };
    try testing.expectEqual(i_pow[n % 4], try det(al, test_io, try ConstView(C64).init(&zi, n, n, n)));
    try testing.expectEqual(@as(f64, 1), try det(al, test_io, try ConstView(f64).init(&.{}, 0, 0, 1)));

    var prng: std.Random.DefaultPrng = .init(0xde7);
    inline for (.{ f64, C64 }) |T| {
        const a = try randomMat(T, prng.random(), n, n);
        defer a.deinit(al);
        shiftDiag(T, a, 1);
        const ai = try inv(al, test_io, a);
        defer ai.deinit(al);
        const p = eMul(try det(al, test_io, a), try det(al, test_io, ai));
        try testing.expect(absOf(if (comptime isComplex(T)) T.init(p.re - 1, p.im) else p - 1) <= 1e-10);
    }
}

test "inv past the inline limit: a a⁻¹ is I" {
    const al = testing.allocator;
    var prng: std.Random.DefaultPrng = .init(0x1a7);
    inline for (.{ f32, f64, C64 }) |T| for ([_]usize{ matrix.inline_limit + 1, 25 }) |n| {
        const a = try randomMat(T, prng.random(), n, n);
        defer a.deinit(al);
        shiftDiag(T, a, 2);
        const ai = try inv(al, test_io, a);
        defer ai.deinit(al);
        const id = try Mat(T).zeros(al, n, n);
        defer id.deinit(al);
        shiftDiag(T, id, 1);
        try expectSolves(T, a.constView(), ai.constView(), id.constView());
    };
}
