//! Tests: `sylvester` and `lyap`.

const std = @import("std");
const matrix = @import("../../matrix.zig");
const ops = @import("../../../kernels/lapack/ops.zig");
const products = @import("../products.zig");
const elementwise = @import("../elementwise.zig");
const sylvester_mod = @import("../sylvester.zig");
const lyap_generalized = @import("../lyap_generalized.zig");
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
const lyapSylvester = sylvester_mod.lyapSylvester;
const lyapGeneralized = lyap_generalized.lyapGeneralized;

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

test "lyapSylvester: A X + X B + C = 0" {
    const al = testing.allocator;
    var prng: std.Random.DefaultPrng = .init(0x1a95);
    inline for (.{ f64, C64 }) |T| for ([_][2]usize{ .{ 1, 1 }, .{ 5, 3 } }) |s| {
        const a = try randomMat(T, prng.random(), s[0], s[0]);
        defer a.deinit(al);
        const b = try randomMat(T, prng.random(), s[1], s[1]);
        defer b.deinit(al);
        shiftDiag(T, a, 3);
        shiftDiag(T, b, 3);
        const c = try randomMat(T, prng.random(), s[0], s[1]);
        defer c.deinit(al);
        const x = try lyapSylvester(al, test_io, a, b, c);
        defer x.deinit(al);
        const ax = try mul(al, test_io, a, x);
        defer ax.deinit(al);
        const xb = try mul(al, test_io, x, b);
        defer xb.deinit(al);
        const s1 = try add(al, ax, xb);
        defer s1.deinit(al);
        const r = try add(al, s1, c);
        defer r.deinit(al);
        const zero = try Mat(T).zeros(al, s[0], s[1]);
        defer zero.deinit(al);
        const n: f64 = @floatFromInt(@max(s[0], s[1]));
        const bound = tolOf(T, s[0] + s[1]) * (n * (maxAbs(T, a.constView()) + maxAbs(T, b.constView())) * maxAbs(T, x.constView()) + maxAbs(T, c.constView()));
        try testing.expect(maxDiff(T, r.constView(), zero.constView()) <= bound);
    };
}

/// `A X Eᴴ + E X Aᴴ + Q`.
fn generalizedResidual(comptime T: type, a: Mat(T), x: Mat(T), e: Mat(T), q: Mat(T)) !Mat(T) {
    const al = testing.allocator;
    const eh = try ctranspose(al, test_io, e);
    defer eh.deinit(al);
    const ah = try ctranspose(al, test_io, a);
    defer ah.deinit(al);
    const ax = try mul(al, test_io, a, x);
    defer ax.deinit(al);
    const axe = try mul(al, test_io, ax, eh);
    defer axe.deinit(al);
    const ex = try mul(al, test_io, e, x);
    defer ex.deinit(al);
    const exa = try mul(al, test_io, ex, ah);
    defer exa.deinit(al);
    const s1 = try add(al, axe, exa);
    defer s1.deinit(al);
    return add(al, s1, q);
}

test "lyapGeneralized: A X Eᴴ + E X Aᴴ + Q = 0; E = I is lyap" {
    const al = testing.allocator;
    var prng: std.Random.DefaultPrng = .init(0x9e1a);
    inline for (.{ f64, C64 }) |T| for ([_]usize{ 1, 3, 8 }) |n| {
        const a = try randomMat(T, prng.random(), n, n);
        defer a.deinit(al);
        shiftDiag(T, a, -4);
        const e = try randomMat(T, prng.random(), n, n);
        defer e.deinit(al);
        shiftDiag(T, e, 3);
        const m = try randomMat(T, prng.random(), n, n);
        defer m.deinit(al);
        const mh = try ctranspose(al, test_io, m);
        defer mh.deinit(al);
        const q = try add(al, m, mh);
        defer q.deinit(al);
        const x = try lyapGeneralized(al, test_io, a, q, e);
        defer x.deinit(al);
        const r = try generalizedResidual(T, a, x, e, q);
        defer r.deinit(al);
        const zero = try Mat(T).zeros(al, n, n);
        defer zero.deinit(al);
        const nf: f64 = @floatFromInt(n);
        const bound = tolOf(T, 4 * n) * (nf * nf * 2 * maxAbs(T, a.constView()) * maxAbs(T, e.constView()) * maxAbs(T, x.constView()) + maxAbs(T, q.constView()));
        try testing.expect(maxDiff(T, r.constView(), zero.constView()) <= bound);

        // E = I: the standard equation
        const eye = try Mat(T).zeros(al, n, n);
        defer eye.deinit(al);
        shiftDiag(T, eye, 1);
        const xg = try lyapGeneralized(al, test_io, a, q, eye);
        defer xg.deinit(al);
        const xs = try lyap(al, test_io, a, q);
        defer xs.deinit(al);
        try testing.expect(maxDiff(T, xg.constView(), xs.constView()) <= tolOf(T, 4 * n) * nf * maxAbs(T, xs.constView()) * maxAbs(T, a.constView()));
    };
}

test "lyapGeneralized: λ_i + conj(λ_j) = 0 is Singular; shape errors" {
    const al = testing.allocator;
    const a = try Mat(f64).zeros(al, 2, 2);
    defer a.deinit(al);
    a.data[0] = 1;
    a.data[3] = -1;
    const eye = try Mat(f64).zeros(al, 2, 2);
    defer eye.deinit(al);
    shiftDiag(f64, eye, 1);
    try testing.expectError(error.Singular, lyapGeneralized(al, test_io, a, eye, eye));
    const wide = try Mat(f64).zeros(al, 2, 3);
    defer wide.deinit(al);
    try testing.expectError(error.NotSquare, lyapGeneralized(al, test_io, wide, eye, eye));
    try testing.expectError(error.DimensionMismatch, lyapGeneralized(al, test_io, a, wide, eye));
    const one = try Mat(f64).zeros(al, 1, 1);
    defer one.deinit(al);
    try testing.expectError(error.DimensionMismatch, lyapGeneralized(al, test_io, a, eye, one));
}
