//! Tests: `cholupdate`.

const std = @import("std");
const matrix = @import("../../matrix.zig");
const vector = @import("../../vector.zig");
const products = @import("../products.zig");
const factor = @import("../factor.zig");
const cholupdate_ = @import("../cholupdate.zig");
const test_inverse = @import("inverse.zig");

const testing = std.testing;
const Mat = matrix.Mat;
const isComplex = vector.isComplex;
const C64 = std.math.Complex(f64);
const test_io = std.Io.failing;
const mul = products.mul;
const ctranspose = products.ctranspose;
const cholupdate = cholupdate_.cholupdate;
const maxAbs = test_inverse.maxAbs;
const maxDiff = test_inverse.maxDiff;
const randomMat = test_inverse.randomMat;
const shiftDiag = test_inverse.shiftDiag;
const tolOf = test_inverse.tolOf;

/// `Rᴴ R`.
fn gram(comptime T: type, r: Mat(T)) !Mat(T) {
    const al = testing.allocator;
    const rh = try ctranspose(al, test_io, r);
    defer rh.deinit(al);
    return mul(al, test_io, rh, r);
}

/// `A ± x xᴴ`, `x` `n x 1`.
fn plusOuter(comptime T: type, a: Mat(T), x: Mat(T), sign: f64) !Mat(T) {
    const al = testing.allocator;
    const n = a.rows;
    const out = try Mat(T).init(al, n, n);
    for (0..n) |j| for (0..n) |i| {
        const xx = if (comptime isComplex(T)) x.data[i].mul(x.data[j].conjugate()) else x.data[i] * x.data[j];
        out.data[i + j * n] = if (comptime isComplex(T)) a.data[i + j * n].add(.init(sign * xx.re, sign * xx.im)) else a.data[i + j * n] + @as(T, @floatCast(sign)) * xx;
    };
    return out;
}

fn expectUpperPositive(comptime T: type, r: Mat(T)) !void {
    for (0..r.rows) |j| {
        const d = r.at(j, j);
        if (comptime isComplex(T)) {
            try testing.expect(d.im == 0 and d.re > 0);
        } else try testing.expect(d > 0);
        for (j + 1..r.rows) |i| {
            const z = r.at(i, j);
            if (comptime isComplex(T)) try testing.expect(z.re == 0 and z.im == 0) else try testing.expectEqual(0, z);
        }
    }
}

test "cholupdate: R1ᴴ R1 = Rᴴ R ± x xᴴ, real and complex" {
    const al = testing.allocator;
    var prng: std.Random.DefaultPrng = .init(0xc401);
    inline for (.{ f32, f64, C64 }) |T| for ([_]usize{ 1, 2, 5, 17 }) |n| {
        // `A = Bᴴ B + n I`, its factor, and `x`
        const b = try randomMat(T, prng.random(), n, n);
        defer b.deinit(al);
        const a = try gram(T, b);
        defer a.deinit(al);
        shiftDiag(T, a, @floatFromInt(n));
        const r = try factor.chol(al, test_io, a);
        defer r.deinit(al);
        const x = try randomMat(T, prng.random(), n, 1);
        defer x.deinit(al);
        const tol = tolOf(T, 4 * n) * @as(f64, @floatFromInt(n)) * (maxAbs(T, a.constView()) + maxAbs(T, x.constView()) * maxAbs(T, x.constView()));

        // update
        const up = try cholupdate(al, r, x, .update);
        defer up.deinit(al);
        try expectUpperPositive(T, up);
        const want = try plusOuter(T, a, x, 1);
        defer want.deinit(al);
        const got = try gram(T, up);
        defer got.deinit(al);
        try testing.expect(maxDiff(T, got.constView(), want.constView()) <= tol);

        // and back down, `x` as a row this time
        const xr = try ctranspose(al, test_io, x);
        defer xr.deinit(al);
        if (comptime isComplex(T)) for (xr.data) |*v| {
            v.* = v.conjugate();
        };
        const down = try cholupdate(al, up, xr, .downdate);
        defer down.deinit(al);
        try expectUpperPositive(T, down);
        const back = try gram(T, down);
        defer back.deinit(al);
        try testing.expect(maxDiff(T, back.constView(), a.constView()) <= tol);
    };
}

test "cholupdate: a downdate past definiteness, and bad shapes" {
    const al = testing.allocator;
    const eye = try Mat(f64).zeros(al, 2, 2);
    defer eye.deinit(al);
    shiftDiag(f64, eye, 1);
    const x = try Mat(f64).zeros(al, 2, 1);
    defer x.deinit(al);
    x.data[0] = 2;
    // `I - x xᴴ` has the eigenvalue -3
    try testing.expectError(error.NotPositiveDefinite, cholupdate(al, eye, x, .downdate));
    const r = try cholupdate(al, eye, x, .update);
    defer r.deinit(al);
    try testing.expectApproxEqAbs(@sqrt(5.0), r.at(0, 0), 1e-15);
    const wide = try Mat(f64).zeros(al, 2, 3);
    defer wide.deinit(al);
    try testing.expectError(error.NotSquare, cholupdate(al, wide, x, .update));
    try testing.expectError(error.DimensionMismatch, cholupdate(al, eye, wide, .update));
    const three = try Mat(f64).zeros(al, 3, 1);
    defer three.deinit(al);
    try testing.expectError(error.DimensionMismatch, cholupdate(al, eye, three, .update));
    // a zero pivot is no Cholesky factor
    const zero = try Mat(f64).zeros(al, 2, 2);
    defer zero.deinit(al);
    try testing.expectError(error.NotPositiveDefinite, cholupdate(al, zero, x, .update));
}
