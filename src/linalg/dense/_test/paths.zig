//! Tests: each path at every size against a wide reference.

const sysinfo = @import("builtin");
const std = @import("std");
const mem = @import("../../../mem.zig");
const Allocator = std.mem.Allocator;
const Io = std.Io;
const blas = @import("../../blas/__root__.zig");
const ext = @import("../../blas/ext/interface.zig");
const l3 = @import("../../blas/blas3/core.zig");
const l2 = @import("../../blas/blas2/core.zig");
const l1 = @import("../../blas/blas1/interface.zig");
const lapack = @import("../../lapack/__root__.zig");
const vector = @import("../../vector.zig");
const matrix = @import("../../matrix.zig");
const ops = @import("../../lapack/ops.zig");

const View = matrix.View;
const ConstView = matrix.ConstView;
const Mat = matrix.Mat;
const isComplex = vector.isComplex;
const Real = vector.Real;
const eMul = matrix.eMul;
const eMulAdd = matrix.eMulAdd;

const zigstd = @import("../../../zigstd.zig");
const bitwise = @import("../../blas/_test/bitwise.zig");
const testing = std.testing;
const sameBitsOrX87 = bitwise.sameBitsOrX87;
const C64 = std.math.Complex(f64);
const test_io = Io.failing;

const test_inverse = @import("../_test/inverse.zig");
const expectSolves = test_inverse.expectSolves;
const maxAbs = test_inverse.maxAbs;
const maxDiff = test_inverse.maxDiff;
const randomMat = test_inverse.randomMat;
const shiftDiag = test_inverse.shiftDiag;
const tolOf = test_inverse.tolOf;
const inverse = @import("../inverse.zig");
const detWith = inverse.detWith;
const invWith = inverse.invWith;
const solveWith = inverse.solveWith;
const products = @import("../products.zig");
const mul = products.mul;
const mulWith = products.mulWith;
const reductions = @import("../reductions.zig");
const dotWith = reductions.dotWith;
const extremeWith = reductions.extremeWith;
const normWith = reductions.normWith;
const sumWith = reductions.sumWith;
const util = @import("../util.zig");
const Path = util.Path;
const constOf = util.constOf;
const strided = util.strided;

/// A wider `T` for references: `f64` for `f32`, `f128` for `f64`.
pub fn Wide(comptime T: type) type {
    const W = if (Real(T) == f32) f64 else f128;
    return if (isComplex(T)) std.math.Complex(W) else W;
}

pub fn widen(x: anytype) Wide(@TypeOf(x)) {
    return if (comptime isComplex(@TypeOf(x))) .init(x.re, x.im) else x;
}

/// `|re| + |im|` of a wide value, as `f64`.
pub fn wideAbs(x: anytype) f64 {
    return @floatCast(if (comptime isComplex(@TypeOf(x))) @abs(x.re) + @abs(x.im) else @abs(x));
}

pub fn wideDiff(x: anytype, ref: Wide(@TypeOf(x))) f64 {
    const w = widen(x);
    return wideAbs(if (comptime isComplex(@TypeOf(x))) w.sub(ref) else w - ref);
}

/// `acc + x y`, unfused.
pub fn wideMulAdd(acc: anytype, x: @TypeOf(acc), y: @TypeOf(acc)) @TypeOf(acc) {
    return if (comptime isComplex(@TypeOf(acc))) acc.add(x.mul(y)) else acc + x * y;
}

/// `det(a)` by LU with partial pivoting in the wide type.
pub fn wideDet(comptime T: type, a: ConstView(T)) !Wide(T) {
    const W = Wide(T);
    const n = a.rows;
    const f = try testing.allocator.alloc(W, n * n);
    defer testing.allocator.free(f);
    for (0..n) |j| for (0..n) |i| {
        f[i + j * n] = widen(a.at(i, j));
    };
    var d = ops.one(W);
    for (0..n) |k| {
        var p = k;
        for (k + 1..n) |i| {
            if (wideAbs(f[i + k * n]) > wideAbs(f[p + k * n])) p = i;
        }
        if (p != k) {
            d = ops.neg(W, d);
            for (0..n) |j| mem.swap(W, &f[k + j * n], &f[p + j * n]);
        }
        const piv = f[k + k * n];
        d = if (comptime isComplex(W)) d.mul(piv) else d * piv;
        for (k + 1..n) |i| {
            const l = ops.neg(W, if (comptime isComplex(W)) f[i + k * n].div(piv) else f[i + k * n] / piv);
            for (k + 1..n) |j| f[i + j * n] = wideMulAdd(f[i + j * n], l, f[k + j * n]);
        }
    }
    return d;
}

pub fn expectPaths(comptime T: type, n: usize, rnd: std.Random) !void {
    const al = testing.allocator;
    const W = Wide(T);
    const a = try randomMat(T, rnd, n, n);
    defer a.deinit(al);
    shiftDiag(T, a, 2);
    const b = try randomMat(T, rnd, n, 3);
    defer b.deinit(al);
    const x = b.constView().col(0);
    const y = b.constView().col(1);
    const sx = try strided(constOf(x));
    const sy = try strided(constOf(y));
    // the references, and the sums of magnitudes that bound their errors
    var dot_ref = ops.zero(W);
    var sum_ref = ops.zero(W);
    var sq: Real(W) = 0;
    var dot_mag: f64 = 0;
    var sum_mag: f64 = 0;
    for (x, y) |xi, yi| {
        dot_ref = wideMulAdd(dot_ref, ops.conj(W, widen(xi)), widen(yi));
        sum_ref = wideMulAdd(sum_ref, widen(xi), ops.one(W));
        const w = widen(xi);
        sq += if (comptime isComplex(W)) w.re * w.re + w.im * w.im else w * w;
        dot_mag += wideAbs(widen(xi)) * wideAbs(widen(yi));
        sum_mag += wideAbs(widen(xi));
    }
    const prod_ref = try al.alloc(W, n * 3);
    defer al.free(prod_ref);
    for (0..3) |j| for (0..n) |i| {
        var acc = ops.zero(W);
        for (0..n) |p| acc = wideMulAdd(acc, widen(a.at(i, p)), widen(b.at(p, j)));
        prod_ref[i + j * n] = acc;
    };
    const det_ref = try wideDet(T, a.constView());
    const tol = tolOf(T, n);
    const nf: f64 = @floatFromInt(n);
    for ([_]Path{ .fixed, .blas }) |path| {
        try testing.expect(wideDiff(try dotWith(path, test_io, sx, sy), dot_ref) <= tol * dot_mag);
        try testing.expect(wideDiff(try sumWith(path, test_io, sx), sum_ref) <= tol * sum_mag);
        try testing.expect(wideDiff(try normWith(path, test_io, sx), @sqrt(sq)) <= tol * @as(f64, @floatCast(@sqrt(sq))));
        if (comptime !isComplex(T)) {
            try testing.expectEqual(std.mem.max(T, x), try extremeWith(.max, path, test_io, sx));
            try testing.expectEqual(std.mem.min(T, x), try extremeWith(.min, path, test_io, sx));
        }
        // `gemm` and, one column, `gemv`
        const o = try Mat(T).init(al, n, 3);
        defer o.deinit(al);
        const mul_bound = tol * nf * maxAbs(T, a.constView()) * maxAbs(T, b.constView());
        for ([_]usize{ 3, 1 }) |k| {
            try mulWith(path, al, test_io, try o.view().sub(0, 0, n, k), a.constView(), try b.constView().sub(0, 0, n, k));
            for (0..k) |j| for (0..n) |i| try testing.expect(wideDiff(o.at(i, j), prod_ref[i + j * n]) <= mul_bound);
        }
        try testing.expect(wideDiff(try detWith(path, al, test_io, a.constView()), det_ref) <= tol * nf * wideAbs(det_ref));
        const ai = try Mat(T).init(al, n, n);
        defer ai.deinit(al);
        try invWith(path, al, test_io, ai.view(), a.constView());
        const p = try mul(al, test_io, a, ai);
        defer p.deinit(al);
        const id = try Mat(T).zeros(al, n, n);
        defer id.deinit(al);
        shiftDiag(T, id, 1);
        try testing.expect(maxDiff(T, p.constView(), id.constView()) <= tol * nf * maxAbs(T, a.constView()) * maxAbs(T, ai.constView()));
        try solveWith(path, al, test_io, o.view(), a.constView(), b.constView());
        try expectSolves(T, a.constView(), o.constView(), b.constView());
    }
}

test "each path at every size 1 to 16, against a wide reference" {
    var prng: std.Random.DefaultPrng = .init(0xd15);
    inline for (.{ f32, f64, C64 }) |T| for (1..17) |n| try expectPaths(T, n, prng.random());
}
