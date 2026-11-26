//! Tests: eigenvalues, standard and generalized.

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
const test_factor = @import("factor.zig");
const expectOrthonormal = test_factor.expectOrthonormal;
const scaleCols = test_factor.scaleCols;
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
