//! Tests: `sylvester` and `lyap`.

const std = @import("std");
const matrix = @import("../../matrix.zig");
const ops = @import("../../../kernels/lapack/ops.zig");
const products = @import("../products.zig");
const elementwise = @import("../elementwise.zig");
const sylvester_mod = @import("../sylvester.zig");
const test_inverse = @import("inverse.zig");

const testing = std.testing;
const Mat = matrix.Mat;
const C64 = std.math.Complex(f64);
const test_io = std.Io.failing;
const mul = products.mul;
const ctranspose = products.ctranspose;
const add = elementwise.add;
const maxAbs = test_inverse.maxAbs;
const maxDiff = test_inverse.maxDiff;
const randomMat = test_inverse.randomMat;
const shiftDiag = test_inverse.shiftDiag;
const tolOf = test_inverse.tolOf;
const sylvester = sylvester_mod.sylvester;
const lyap = sylvester_mod.lyap;

test "sylvester: A X + X B = C, real and complex, square and not" {
    const al = testing.allocator;
    var prng: std.Random.DefaultPrng = .init(0x5715);
    inline for (.{ f32, f64, C64 }) |T| for ([_][2]usize{ .{ 1, 1 }, .{ 4, 3 }, .{ 3, 7 }, .{ 12, 12 } }) |s| {
        const a = try randomMat(T, prng.random(), s[0], s[0]);
        defer a.deinit(al);
        const b = try randomMat(T, prng.random(), s[1], s[1]);
        defer b.deinit(al);
        // eigenvalues of A near 3, of -B near -3: well apart
        shiftDiag(T, a, 3);
        shiftDiag(T, b, 3);
        const c = try randomMat(T, prng.random(), s[0], s[1]);
        defer c.deinit(al);
        const x = try sylvester(al, test_io, a, b, c);
        defer x.deinit(al);
        const ax = try mul(al, test_io, a, x);
        defer ax.deinit(al);
        const xb = try mul(al, test_io, x, b);
        defer xb.deinit(al);
        const lhs = try add(al, ax, xb);
        defer lhs.deinit(al);
        const n: f64 = @floatFromInt(@max(s[0], s[1]));
        const bound = tolOf(T, s[0] + s[1]) * (n * (maxAbs(T, a.constView()) + maxAbs(T, b.constView())) * maxAbs(T, x.constView()) + maxAbs(T, c.constView()));
        try testing.expect(maxDiff(T, lhs.constView(), c.constView()) <= bound);
    };
}

test "lyap: A X + X Aᴴ + Q = 0, X Hermitian for a Hermitian Q" {
    const al = testing.allocator;
    var prng: std.Random.DefaultPrng = .init(0x1a9);
    inline for (.{ f64, C64 }) |T| for ([_]usize{ 2, 9 }) |n| {
        const a = try randomMat(T, prng.random(), n, n);
        defer a.deinit(al);
        shiftDiag(T, a, -4);
        const m = try randomMat(T, prng.random(), n, n);
        defer m.deinit(al);
        const mh = try ctranspose(al, test_io, m);
        defer mh.deinit(al);
        const q = try add(al, m, mh);
        defer q.deinit(al);
        const x = try lyap(al, test_io, a, q);
        defer x.deinit(al);
        const ax = try mul(al, test_io, a, x);
        defer ax.deinit(al);
        const ah = try ctranspose(al, test_io, a);
        defer ah.deinit(al);
        const xah = try mul(al, test_io, x, ah);
        defer xah.deinit(al);
        const s1 = try add(al, ax, xah);
        defer s1.deinit(al);
        const r = try add(al, s1, q);
        defer r.deinit(al);
        const zero = try Mat(T).zeros(al, n, n);
        defer zero.deinit(al);
        const bound = tolOf(T, 2 * n) * (@as(f64, @floatFromInt(n)) * 2 * maxAbs(T, a.constView()) * maxAbs(T, x.constView()) + maxAbs(T, q.constView()));
        try testing.expect(maxDiff(T, r.constView(), zero.constView()) <= bound);
        const xh = try ctranspose(al, test_io, x);
        defer xh.deinit(al);
        try testing.expect(maxDiff(T, x.constView(), xh.constView()) <= bound);
    };
}

test "sylvester: common eigenvalues are Singular; shape errors" {
    const al = testing.allocator;
    const one = try Mat(f64).zeros(al, 1, 1);
    defer one.deinit(al);
    one.data[0] = 1;
    const minus = try Mat(f64).zeros(al, 1, 1);
    defer minus.deinit(al);
    minus.data[0] = -1;
    // 1 x + x (-1) = 1 has no solution
    try testing.expectError(error.Singular, sylvester(al, test_io, one, minus, one));
    const wide = try Mat(f64).zeros(al, 1, 2);
    defer wide.deinit(al);
    try testing.expectError(error.NotSquare, sylvester(al, test_io, wide, one, one));
    try testing.expectError(error.DimensionMismatch, sylvester(al, test_io, one, one, wide));
    try testing.expectError(error.DimensionMismatch, lyap(al, test_io, one, wide));
}
