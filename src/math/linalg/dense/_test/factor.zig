//! Tests: factorizations, norms, SVD and eigenvalues.

const sysinfo = @import("builtin");
const std = @import("std");
const mem = @import("../../../../mem.zig");
const Allocator = std.mem.Allocator;
const Io = std.Io;
const blas = @import("../../../kernels/blas/__root__.zig");
const ext = @import("../../../kernels/blas/ext/interface.zig");
const l3 = @import("../../../kernels/blas/blas3/core.zig");
const l2 = @import("../../../kernels/blas/blas2/core.zig");
const l1 = @import("../../../kernels/blas/blas1/interface.zig");
const lapack = @import("../../../kernels/lapack/__root__.zig");
const vector = @import("../../vector.zig");
const matrix = @import("../../matrix.zig");
const ops = @import("../../../kernels/lapack/ops.zig");

const View = matrix.View;
const ConstView = matrix.ConstView;
const Mat = matrix.Mat;
const isComplex = vector.isComplex;
const Real = vector.Real;
const eMul = matrix.eMul;
const eMulAdd = matrix.eMulAdd;

const zigstd = @import("../../../../zigstd.zig");
const bitwise = @import("../../../kernels/blas/_test/bitwise.zig");
const testing = std.testing;
const sameBitsOrX87 = bitwise.sameBitsOrX87;
const C64 = std.math.Complex(f64);
const test_io = Io.failing;

const test_inverse = @import("inverse.zig");
const hermitianPd = test_inverse.hermitianPd;
const hermitianZeroDiag = test_inverse.hermitianZeroDiag;
const maxAbs = test_inverse.maxAbs;
const maxDiff = test_inverse.maxDiff;
const randomMat = test_inverse.randomMat;
const shiftDiag = test_inverse.shiftDiag;
const tolOf = test_inverse.tolOf;
const eigen = @import("../eigen.zig");
const eig = eigen.eig;
const eigGeneralized = eigen.eigGeneralized;
const factor = @import("../factor.zig");
const Lu = factor.Lu;
const QrMode = factor.QrMode;
const chol = factor.chol;
const ldl = factor.ldl;
const lu = factor.lu;
const qr = factor.qr;
const products = @import("../products.zig");
const ctranspose = products.ctranspose;
const mul = products.mul;
const singular = @import("../singular.zig");
const SvdMode = singular.SvdMode;
const cond = singular.cond;
const matrixNorm = singular.matrixNorm;
const pinv = singular.pinv;
const rank = singular.rank;
const svd = singular.svd;
const util = @import("../util.zig");
const copyOf = util.copyOf;

/// `P A = L U`, rows swapped by `ipiv` in order.
pub fn expectLu(comptime T: type, a: Mat(T)) !Lu(T) {
    const al = testing.allocator;
    const f = try lu(al, test_io, a);
    errdefer f.deinit(al);
    const l = try f.l(al);
    defer l.deinit(al);
    const u = try f.u(al);
    defer u.deinit(al);
    const p = try mul(al, test_io, l, u);
    defer p.deinit(al);
    const pa = try copyOf(T, al, a.constView());
    defer pa.deinit(al);
    for (f.ipiv, 0..) |r, i| {
        if (r == i) continue;
        for (0..pa.cols) |j| mem.swap(T, &pa.data[i + j * pa.rows], &pa.data[r + j * pa.rows]);
    }
    try testing.expect(maxDiff(T, pa.constView(), p.constView()) <= tolOf(T, a.cols) * maxAbs(T, a.constView()));
    return f;
}

test "lu: P A = L U, singular too" {
    const al = testing.allocator;
    var prng: std.Random.DefaultPrng = .init(0x10);
    inline for (.{ f32, f64, C64 }) |T| for ([_][2]usize{ .{ 10, 7 }, .{ 7, 10 }, .{ 12, 12 }, .{ 3, 2 } }) |s| {
        const a = try randomMat(T, prng.random(), s[0], s[1]);
        defer a.deinit(al);
        const f = try expectLu(T, a);
        defer f.deinit(al);
        try testing.expectEqual(@as(?usize, null), f.zero_pivot);
    };
    // a zero first column: `U(0, 0)` is 0, and the rest still factors
    const a = try randomMat(f64, prng.random(), 9, 9);
    defer a.deinit(al);
    @memset(a.data[0..9], 0);
    const f = try expectLu(f64, a);
    defer f.deinit(al);
    try testing.expectEqual(@as(?usize, 0), f.zero_pivot);
}

test "chol: Rᴴ R is A, R upper" {
    const al = testing.allocator;
    var prng: std.Random.DefaultPrng = .init(0xc401);
    inline for (.{ f32, f64, C64 }) |T| for ([_]usize{ 3, 20 }) |n| {
        const a = try hermitianPd(T, prng.random(), n);
        defer a.deinit(al);
        const r = try chol(al, test_io, a);
        defer r.deinit(al);
        for (0..n) |j| for (j + 1..n) |i| try testing.expect(ops.isZero(T, r.at(i, j)));
        const rh = try ctranspose(al, test_io, r);
        defer rh.deinit(al);
        const p = try mul(al, test_io, rh, r);
        defer p.deinit(al);
        try testing.expect(maxDiff(T, p.constView(), a.constView()) <= tolOf(T, n) * maxAbs(T, a.constView()));
    };
}

test "ldl: Pᵀ A P = L D Lᴴ, L unit lower, D Hermitian block diagonal" {
    const al = testing.allocator;
    var prng: std.Random.DefaultPrng = .init(0x1d1);
    inline for (.{ f32, f64, C64 }) |T| for ([_]usize{ 1, 2, 5, 20 }) |n| for ([_]bool{ false, true }) |zd| {
        const a = if (zd) try hermitianZeroDiag(T, prng.random(), n) else try hermitianPd(T, prng.random(), n);
        defer a.deinit(al);
        // the upper triangle is not read
        for (0..n) |j| for (0..j) |i| {
            a.data[i + j * n] = ops.zero(T);
        };
        const f = try ldl(al, test_io, a);
        defer f.deinit(al);
        for (0..n) |j| for (0..j + 1) |i| {
            const want = if (i == j) ops.one(T) else ops.zero(T);
            try testing.expectEqual(want, f.l.at(i, j));
        };
        // `D`: blocks of 1 or 2 on the diagonal, Hermitian
        for (0..n) |j| for (0..n) |i| {
            const d = f.d.at(i, j);
            if (i > j + 1 or j > i + 1) try testing.expect(ops.isZero(T, d));
            try testing.expectEqual(ops.conj(T, d), f.d.at(j, i));
        };
        const seen = try al.alloc(bool, n);
        defer al.free(seen);
        @memset(seen, false);
        for (f.perm) |p| seen[p] = true;
        for (seen) |s| try testing.expect(s);
        // `Pᵀ A P`, `A` Hermitian from its lower triangle
        const pap = try Mat(T).init(al, n, n);
        defer pap.deinit(al);
        for (0..n) |j| for (0..n) |i| {
            const r = f.perm[i];
            const c = f.perm[j];
            pap.data[i + j * n] = if (r >= c) a.at(r, c) else ops.conj(T, a.at(c, r));
        };
        const ld = try mul(al, test_io, f.l, f.d);
        defer ld.deinit(al);
        const lh = try ctranspose(al, test_io, f.l);
        defer lh.deinit(al);
        const ldlh = try mul(al, test_io, ld, lh);
        defer ldlh.deinit(al);
        try testing.expect(maxDiff(T, ldlh.constView(), pap.constView()) <= tolOf(T, n) * @max(maxAbs(T, pap.constView()), 1));
    };
    const wide = try Mat(f64).zeros(al, 2, 3);
    defer wide.deinit(al);
    try testing.expectError(error.NotSquare, ldl(al, test_io, wide));
}

pub fn expectQr(comptime T: type, a: Mat(T), mode: QrMode) !void {
    const al = testing.allocator;
    const f = try qr(al, test_io, a, mode);
    defer f.deinit(al);
    const m = a.rows;
    const qc = if (mode == .full) m else @min(m, a.cols);
    try testing.expect(f.q.rows == m and f.q.cols == qc and f.r.rows == qc and f.r.cols == a.cols);
    for (0..a.cols) |j| for (@min(j + 1, qc)..qc) |i| try testing.expect(ops.isZero(T, f.r.at(i, j)));
    const p = try mul(al, test_io, f.q, f.r);
    defer p.deinit(al);
    try testing.expect(maxDiff(T, p.constView(), a.constView()) <= tolOf(T, m) * @max(maxAbs(T, a.constView()), 1));
    const qh = try ctranspose(al, test_io, f.q);
    defer qh.deinit(al);
    const qhq = try mul(al, test_io, qh, f.q);
    defer qhq.deinit(al);
    const id = try Mat(T).zeros(al, qc, qc);
    defer id.deinit(al);
    shiftDiag(T, id, 1);
    try testing.expect(maxDiff(T, qhq.constView(), id.constView()) <= tolOf(T, m));
}

test "qr: Q R is A, Q orthonormal, both modes" {
    const al = testing.allocator;
    var prng: std.Random.DefaultPrng = .init(0x9a);
    inline for (.{ f32, f64, C64 }) |T| for ([_][2]usize{ .{ 12, 7 }, .{ 7, 12 }, .{ 5, 5 }, .{ 5, 0 } }) |s| {
        const a = try randomMat(T, prng.random(), s[0], s[1]);
        defer a.deinit(al);
        try expectQr(T, a, .economy);
        try expectQr(T, a, .full);
    };
}

test "matrixNorm: one, inf, fro and two" {
    const al = testing.allocator;
    const a = matrix.fromRows(matrix.M2f64, .{ .{ 1, -2 }, .{ 3, 4 } });
    try testing.expectEqual(@as(f64, 6), try matrixNorm(al, test_io, &a, .one));
    try testing.expectEqual(@as(f64, 7), try matrixNorm(al, test_io, &a, .inf));
    try testing.expectApproxEqRel(@sqrt(30.0), try matrixNorm(al, test_io, &a, .fro), 1e-15);
    try testing.expectApproxEqRel(@sqrt(15 + @sqrt(125.0)), try matrixNorm(al, test_io, &a, .two), 1e-14);
    const z: matrix.StaticMat(1, 2, C64) = .{ .data = .{ .{.init(3, 4)}, .{.init(0, 0)} } };
    try testing.expectApproxEqRel(@as(f64, 5), try matrixNorm(al, test_io, &z, .two), 1e-15);
    try testing.expectEqual(@as(f64, 0), try matrixNorm(al, test_io, try ConstView(f64).init(&.{}, 0, 3, 1), .two));
}

/// `x(:, j) *= d[j]`, real `d`.
pub fn scaleCols(comptime T: type, x: Mat(T), d: []const Real(T)) void {
    for (d, 0..) |dj, j| for (x.data[j * x.rows ..][0..x.rows]) |*e| {
        e.* = if (comptime isComplex(T)) T.init(e.re * dj, e.im * dj) else e.* * dj;
    };
}

/// `qᴴ q` is `I`.
pub fn expectOrthonormal(comptime T: type, q: Mat(T)) !void {
    const al = testing.allocator;
    const qh = try ctranspose(al, test_io, q);
    defer qh.deinit(al);
    const qhq = try mul(al, test_io, qh, q);
    defer qhq.deinit(al);
    const id = try Mat(T).zeros(al, q.cols, q.cols);
    defer id.deinit(al);
    shiftDiag(T, id, 1);
    try testing.expect(maxDiff(T, qhq.constView(), id.constView()) <= tolOf(T, q.rows));
}

pub fn expectSvd(comptime T: type, a: Mat(T), mode: SvdMode) !void {
    const al = testing.allocator;
    const m = a.rows;
    const n = a.cols;
    const k = @min(m, n);
    const f = try svd(al, test_io, a, mode);
    defer f.deinit(al);
    try testing.expectEqual(k, f.s.len);
    for (f.s, 0..) |x, i| try testing.expect(x >= 0 and (i == 0 or x <= f.s[i - 1]));
    const values = try svd(al, test_io, a, .values);
    defer values.deinit(al);
    for (f.s, values.s) |x, y| try testing.expect(@abs(x - y) <= tolOf(T, @max(m, n)) * @max(values.s[0], 1));
    if (mode == .values) return testing.expect(f.u == null and f.v == null);
    const u = f.u.?;
    const v = f.v.?;
    const uc = if (mode == .full) m else k;
    const vc = if (mode == .full) n else k;
    try testing.expect(u.rows == m and u.cols == uc and v.rows == n and v.cols == vc);
    try expectOrthonormal(T, u);
    try expectOrthonormal(T, v);
    // `U(:, 1:k) Σ V(:, 1:k)ᴴ`
    const us = try copyOf(T, al, try u.constView().sub(0, 0, m, k));
    defer us.deinit(al);
    scaleCols(T, us, f.s);
    const vh = try ctranspose(al, test_io, try v.constView().sub(0, 0, n, k));
    defer vh.deinit(al);
    const p = try mul(al, test_io, us, vh);
    defer p.deinit(al);
    try testing.expect(maxDiff(T, p.constView(), a.constView()) <= tolOf(T, @max(m, n)) * @max(maxAbs(T, a.constView()), 1));
}

test "svd: U Σ Vᴴ is A, U and V orthonormal, every mode" {
    const al = testing.allocator;
    var prng: std.Random.DefaultPrng = .init(0x5fd);
    inline for (.{ f32, f64, C64 }) |T| for ([_][2]usize{ .{ 12, 7 }, .{ 7, 12 }, .{ 5, 5 }, .{ 5, 0 }, .{ 0, 3 }, .{ 1, 1 } }) |s| {
        const a = try randomMat(T, prng.random(), s[0], s[1]);
        defer a.deinit(al);
        inline for (.{ .values, .economy, .full }) |mode| try expectSvd(T, a, mode);
    };
    // past `gesdd`'s divide and conquer threshold (128)
    const big = try randomMat(f64, prng.random(), 140, 130);
    defer big.deinit(al);
    try expectSvd(f64, big, .economy);
    // known values: `diag(3, -4)` padded with a zero row
    const d = matrix.fromRows(matrix.StaticMat(3, 2, f64), .{ .{ 3, 0 }, .{ 0, -4 }, .{ 0, 0 } });
    const f = try svd(al, test_io, &d, .full);
    defer f.deinit(al);
    try testing.expectApproxEqRel(@as(f64, 4), f.s[0], 1e-15);
    try testing.expectApproxEqRel(@as(f64, 3), f.s[1], 1e-15);
}

/// `a V = V Λ`, `V` orthonormal, `Λ` nondecreasing.
pub fn expectHermitianEig(comptime T: type, a: Mat(T)) !void {
    const al = testing.allocator;
    const n = a.rows;
    const e = try eig(al, test_io, a, true);
    defer e.deinit(al);
    const h = e.hermitian;
    const v = h.vectors.?;
    for (h.values[1..], h.values[0 .. n - 1]) |x, y| try testing.expect(x >= y);
    try expectOrthonormal(T, v);
    const av = try mul(al, test_io, a, v);
    defer av.deinit(al);
    const vl = try copyOf(T, al, v.constView());
    defer vl.deinit(al);
    scaleCols(T, vl, h.values);
    const bound = tolOf(T, n) * @as(f64, @floatFromInt(n)) * maxAbs(T, a.constView());
    try testing.expect(maxDiff(T, av.constView(), vl.constView()) <= bound);
    const values = try eig(al, test_io, a, false);
    defer values.deinit(al);
    try testing.expect(values.hermitian.vectors == null);
    for (h.values, values.hermitian.values) |x, y| try testing.expect(@abs(x - y) <= bound);
}

/// `a v = λ v` per column, each `v` of unit norm.
pub fn expectGeneralEig(comptime T: type, a: Mat(T)) !void {
    const C = std.math.Complex(Real(T));
    const al = testing.allocator;
    const n = a.rows;
    const e = try eig(al, test_io, a, true);
    defer e.deinit(al);
    const g = e.general;
    const v = g.vectors.?;
    const ac = try Mat(C).init(al, n, n);
    defer ac.deinit(al);
    for (ac.data, a.data) |*d, x| d.* = if (comptime isComplex(T)) x else C.init(x, 0);
    const av = try mul(al, test_io, ac, v);
    defer av.deinit(al);
    const vl = try copyOf(C, al, v.constView());
    defer vl.deinit(al);
    for (g.values, 0..) |l, j| for (vl.data[j * n ..][0..n]) |*x| {
        x.* = x.mul(l);
    };
    const bound = tolOf(T, n) * @as(f64, @floatFromInt(n)) * maxAbs(T, a.constView());
    try testing.expect(maxDiff(C, av.constView(), vl.constView()) <= bound);
    for (0..n) |j| {
        var s: f64 = 0;
        for (v.constView().col(j)) |x| s += x.re * x.re + x.im * x.im;
        try testing.expect(@abs(s - 1) <= tolOf(T, n));
    }
    const values = try eig(al, test_io, a, false);
    defer values.deinit(al);
    try testing.expect(values.general.vectors == null);
    for (g.values, values.general.values) |x, y| try testing.expect(x.re == y.re and x.im == y.im);
}

test "eig: Hermitian input real, the rest complex" {
    const al = testing.allocator;
    var prng: std.Random.DefaultPrng = .init(0xe16);
    inline for (.{ f32, f64, C64 }) |T| for ([_]usize{ 1, 5, 20 }) |n| {
        const h = try hermitianPd(T, prng.random(), n);
        defer h.deinit(al);
        try expectHermitianEig(T, h);
        const a = try randomMat(T, prng.random(), n, n);
        defer a.deinit(al);
        // exactly Hermitian at random is impossible past 1 by 1
        if (n > 1) try expectGeneralEig(T, a);
    };
    // a rotation: `±i`, the positive imaginary part first, vectors conjugate
    const r = matrix.fromRows(matrix.M2f64, .{ .{ 0, -1 }, .{ 1, 0 } });
    const e = try eig(al, test_io, &r, true);
    defer e.deinit(al);
    const w = e.general.values;
    try testing.expect(@abs(w[0].re) + @abs(w[0].im - 1) <= 1e-15 and @abs(w[1].re) + @abs(w[1].im + 1) <= 1e-15);
    const v = e.general.vectors.?;
    for (v.constView().col(0), v.constView().col(1)) |x, y| try testing.expect(x.re == y.re and x.im == -y.im);
    // symmetric: real, ascending
    const s = matrix.fromRows(matrix.M2f64, .{ .{ 2, 1 }, .{ 1, 2 } });
    const es = try eig(al, test_io, &s, false);
    defer es.deinit(al);
    try testing.expectApproxEqRel(@as(f64, 1), es.hermitian.values[0], 1e-15);
    try testing.expectApproxEqRel(@as(f64, 3), es.hermitian.values[1], 1e-15);
    // empty is Hermitian; not square fails
    const z = try eig(al, test_io, try ConstView(f64).init(&.{}, 0, 0, 1), true);
    defer z.deinit(al);
    try testing.expect(z.hermitian.values.len == 0 and z.hermitian.vectors.?.data.len == 0);
    try testing.expectError(error.NotSquare, eig(al, test_io, try ConstView(f64).init(&.{ 1, 2 }, 2, 1, 2), false));
}

test "a NaN gives error.NoConvergence" {
    const al = testing.allocator;
    const nan = std.math.nan(f64);
    const g = matrix.fromRows(matrix.M2f64, .{ .{ 1, nan }, .{ 3, 4 } });
    const s = matrix.fromRows(matrix.M2f64, .{ .{ nan, 1 }, .{ 1, 4 } });
    inline for (.{ &g, &s }) |a| {
        try testing.expectError(error.NoConvergence, eig(al, test_io, a, true));
        try testing.expectError(error.NoConvergence, svd(al, test_io, a, .values));
        try testing.expectError(error.NoConvergence, svd(al, test_io, a, .economy));
        try testing.expectError(error.NoConvergence, pinv(al, test_io, a, null));
        try testing.expectError(error.NoConvergence, rank(al, test_io, a, null));
        try testing.expectError(error.NoConvergence, cond(al, test_io, a));
    }
}

/// The four Moore-Penrose conditions: `a x a = a`, `x a x = x`, `a x` and
/// `x a` Hermitian.
pub fn expectPinv(comptime T: type, a: Mat(T)) !void {
    const al = testing.allocator;
    const x = try pinv(al, test_io, a, null);
    defer x.deinit(al);
    try testing.expect(x.rows == a.cols and x.cols == a.rows);
    const ax = try mul(al, test_io, a, x);
    defer ax.deinit(al);
    const xa = try mul(al, test_io, x, a);
    defer xa.deinit(al);
    const axa = try mul(al, test_io, ax, a);
    defer axa.deinit(al);
    const xax = try mul(al, test_io, xa, x);
    defer xax.deinit(al);
    const axh = try ctranspose(al, test_io, ax);
    defer axh.deinit(al);
    const xah = try ctranspose(al, test_io, xa);
    defer xah.deinit(al);
    const c = try cond(al, test_io, a);
    const t = tolOf(T, @max(a.rows, a.cols)) * @max(1, if (std.math.isInf(c)) 1 else c);
    const na = @max(maxAbs(T, a.constView()), 1);
    const nx = @max(maxAbs(T, x.constView()), 1);
    try testing.expect(maxDiff(T, axa.constView(), a.constView()) <= t * na);
    try testing.expect(maxDiff(T, xax.constView(), x.constView()) <= t * nx);
    try testing.expect(maxDiff(T, axh.constView(), ax.constView()) <= t);
    try testing.expect(maxDiff(T, xah.constView(), xa.constView()) <= t);
}

test "pinv, rank and cond" {
    const al = testing.allocator;
    var prng: std.Random.DefaultPrng = .init(0x919);
    inline for (.{ f32, f64, C64 }) |T| {
        for ([_][2]usize{ .{ 12, 7 }, .{ 7, 12 }, .{ 9, 9 } }) |s| {
            const a = try randomMat(T, prng.random(), s[0], s[1]);
            defer a.deinit(al);
            try expectPinv(T, a);
            try testing.expectEqual(@min(s[0], s[1]), try rank(al, test_io, a, null));
        }
        // rank 3: `x yᴴ`, 10 by 8
        const x = try randomMat(T, prng.random(), 10, 3);
        defer x.deinit(al);
        const y = try randomMat(T, prng.random(), 3, 8);
        defer y.deinit(al);
        const a = try mul(al, test_io, x, y);
        defer a.deinit(al);
        try testing.expectEqual(@as(usize, 3), try rank(al, test_io, a, null));
        try expectPinv(T, a);
        try testing.expect(try cond(al, test_io, a) > 1e-3 / std.math.floatEps(Real(T)));
    }
    // `rank`'s tolerance, and `cond` the ratio of the extremes
    const d = matrix.fromRows(matrix.M3f64, .{ .{ 4, 0, 0 }, .{ 0, 1e-3, 0 }, .{ 0, 0, -1e-9 } });
    try testing.expectEqual(@as(usize, 3), try rank(al, test_io, &d, null));
    try testing.expectEqual(@as(usize, 2), try rank(al, test_io, &d, 1e-6));
    try testing.expectApproxEqRel(@as(f64, 4e9), try cond(al, test_io, &d), 1e-14);
    // `pinv(diag(4, 1e-3, -1e-9), 1e-6)` drops the last value
    const p = try pinv(al, test_io, &d, 1e-6);
    defer p.deinit(al);
    try testing.expectApproxEqRel(@as(f64, 0.25), p.at(0, 0), 1e-15);
    try testing.expectApproxEqRel(@as(f64, 1e3), p.at(1, 1), 1e-14);
    try testing.expectEqual(@as(f64, 0), p.at(2, 2));
    // zero and empty
    const z: matrix.M3f64 = .{ .data = @splat(@splat(0)) };
    try testing.expectEqual(@as(usize, 0), try rank(al, test_io, &z, null));
    try testing.expect(std.math.isInf(try cond(al, test_io, &z)));
    const pz = try pinv(al, test_io, &z, null);
    defer pz.deinit(al);
    try testing.expect(maxAbs(f64, pz.constView()) == 0);
    const e = try ConstView(f64).init(&.{}, 0, 3, 1);
    try testing.expectEqual(@as(usize, 0), try rank(al, test_io, e, null));
    try testing.expectEqual(@as(f64, 0), try cond(al, test_io, e));
    const pe = try pinv(al, test_io, e, null);
    defer pe.deinit(al);
    try testing.expect(pe.rows == 3 and pe.cols == 0);
}

test "eigGeneralized: A V = B V Λ and Vᴴ B V = I for a Hermitian definite pencil; errors" {
    const al = testing.allocator;
    var prng: std.Random.DefaultPrng = .init(0xe16b);
    inline for (.{ f64, C64 }) |T| for ([_]usize{ 1, 5, 20 }) |n| {
        const a = try hermitianZeroDiag(T, prng.random(), n);
        defer a.deinit(al);
        const b = try hermitianPd(T, prng.random(), n);
        defer b.deinit(al);
        const e = try eigGeneralized(al, test_io, a, b, true);
        defer e.deinit(al);
        const h = e.hermitian;
        for (1..n) |i| try testing.expect(h.values[i - 1] <= h.values[i]);
        const v = h.vectors.?;
        const av = try mul(al, test_io, a, v);
        defer av.deinit(al);
        const bv = try mul(al, test_io, b, v);
        defer bv.deinit(al);
        for (0..n) |j| for (bv.data[j * n ..][0..n]) |*x| {
            x.* = ops.mulReal(T, x.*, h.values[j]);
        };
        const scale = (maxAbs(T, a.constView()) + maxAbs(T, b.constView()) * @abs(h.values[n - 1]) + 1) * @max(maxAbs(T, v.constView()), 1);
        try testing.expect(maxDiff(T, av.constView(), bv.constView()) <= tolOf(T, n) * scale);
        const vh = try ctranspose(al, test_io, v);
        defer vh.deinit(al);
        const bv2 = try mul(al, test_io, b, v);
        defer bv2.deinit(al);
        const vbv = try mul(al, test_io, vh, bv2);
        defer vbv.deinit(al);
        const id = try Mat(T).zeros(al, n, n);
        defer id.deinit(al);
        shiftDiag(T, id, 1);
        try testing.expect(maxDiff(T, vbv.constView(), id.constView()) <= tolOf(T, n) * maxAbs(T, b.constView()) * @max(maxAbs(T, v.constView()), 1));
        const vals = try eigGeneralized(al, test_io, a, b, false);
        defer vals.deinit(al);
        try testing.expect(vals.hermitian.vectors == null);
    };
    const g = try randomMat(f64, prng.random(), 3, 3);
    defer g.deinit(al);
    const pd = try hermitianPd(f64, prng.random(), 3);
    defer pd.deinit(al);
    try testing.expectError(error.NotHermitianPencil, eigGeneralized(al, test_io, g, pd, false));
    const neg = try Mat(f64).zeros(al, 3, 3);
    defer neg.deinit(al);
    shiftDiag(f64, neg, -1);
    try testing.expectError(error.NotPositiveDefinite, eigGeneralized(al, test_io, pd, neg, false));
}
