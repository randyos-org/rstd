//! Tests: `ordqz`.

const std = @import("std");
const matrix = @import("../../matrix.zig");
const vector = @import("../../vector.zig");
const products = @import("../products.zig");
const schur_form = @import("../schur_form.zig");
const ordqz_ = @import("../ordqz.zig");
const test_inverse = @import("inverse.zig");

const testing = std.testing;
const Mat = matrix.Mat;
const isComplex = vector.isComplex;
const Real = vector.Real;
const C64 = std.math.Complex(f64);
const test_io = std.Io.failing;
const mul = products.mul;
const ordqz = ordqz_.ordqz;
const maxAbs = test_inverse.maxAbs;
const maxDiff = test_inverse.maxDiff;
const randomMat = test_inverse.randomMat;
const tolOf = test_inverse.tolOf;

/// `‖Q M Z - X‖` over the entries.
fn residual(comptime T: type, q: Mat(T), m: Mat(T), z: Mat(T), x: Mat(T)) !f64 {
    const al = testing.allocator;
    const qm = try mul(al, test_io, q, m);
    defer qm.deinit(al);
    const qmz = try mul(al, test_io, qm, z);
    defer qmz.deinit(al);
    return maxDiff(T, qmz.constView(), x.constView());
}

/// The form's eigenvalues, a 1 x 1 block's `AA(k, k) / BB(k, k)`, a real
/// pair's by its quadratic.
fn eigs(comptime T: type, f: schur_form.Qz(T), out: []C64) void {
    const n = f.aa.rows;
    var k: usize = 0;
    while (k < n) : (k += 1) {
        const a = f.aa.data;
        const b = f.bb.data;
        if (comptime isComplex(T)) {
            const l = a[k + k * n].div(b[k + k * n]);
            out[k] = C64.init(l.re, l.im);
        } else if (k + 1 < n and a[k + 1 + k * n] != 0) {
            const db = b[k + k * n] * b[k + 1 + (k + 1) * n] - b[k + (k + 1) * n] * b[k + 1 + k * n];
            const s = a[k + k * n] * b[k + 1 + (k + 1) * n] + a[k + 1 + (k + 1) * n] * b[k + k * n] - a[k + (k + 1) * n] * b[k + 1 + k * n] - a[k + 1 + k * n] * b[k + (k + 1) * n];
            const da = a[k + k * n] * a[k + 1 + (k + 1) * n] - a[k + (k + 1) * n] * a[k + 1 + k * n];
            const re = s / (2 * db);
            const im = @sqrt(@max(da / db - re * re, 0));
            out[k] = C64.init(re, im);
            out[k + 1] = C64.init(re, -im);
            k += 1;
        } else out[k] = C64.init(a[k + k * n] / b[k + k * n], 0);
    }
}

fn expectForm(comptime T: type, a: Mat(T), b: Mat(T), f: schur_form.Qz(T)) !void {
    const n = a.rows;
    const tol = tolOf(T, 4 * n) * @as(f64, @floatFromInt(n * n)) * (maxAbs(T, a.constView()) + maxAbs(T, b.constView()));
    try testing.expect(try residual(T, f.q, a, f.z, f.aa) <= tol);
    try testing.expect(try residual(T, f.q, b, f.z, f.bb) <= tol);
    for (0..n) |j| for (j + 1..n) |i| {
        const bz = if (comptime isComplex(T)) f.bb.at(i, j).re == 0 and f.bb.at(i, j).im == 0 else f.bb.at(i, j) == 0;
        try testing.expect(bz);
        if (comptime isComplex(T)) try testing.expect(f.aa.at(i, j).re == 0 and f.aa.at(i, j).im == 0) else if (i > j + 1) try testing.expectEqual(0, f.aa.at(i, j));
    };
}

test "ordqz: a mask and the unit disk lead, Q A Z = AA kept" {
    const al = testing.allocator;
    var prng: std.Random.DefaultPrng = .init(0x0dc2);
    inline for (.{ f64, C64 }) |T| for ([_]usize{ 3, 8 }) |n| {
        const a = try randomMat(T, prng.random(), n, n);
        defer a.deinit(al);
        const b = try randomMat(T, prng.random(), n, n);
        defer b.deinit(al);
        const f = try schur_form.qz(al, test_io, a, b);
        defer f.deinit(al);
        var before: [16]C64 = undefined;
        eigs(T, f, before[0..n]);
        // the unit disk's eigenvalues lead
        const o = try ordqz(al, test_io, f.aa, f.bb, f.q, f.z, .{ .region = .udi });
        defer o.deinit(al);
        try expectForm(T, a, b, o);
        var after: [16]C64 = undefined;
        eigs(T, o, after[0..n]);
        var inside: usize = 0;
        for (before[0..n]) |e| inside += @intFromBool(e.magnitude() < 1);
        for (after[0..n], 0..) |e, k| try testing.expectEqual(k < inside, e.magnitude() < 1);
        // a mask picking the last eigenvalue (a real pair's both) leads it
        var pick: [16]bool = @splat(false);
        pick[n - 1] = true;
        const m = try ordqz(al, test_io, f.aa, f.bb, f.q, f.z, .{ .mask = pick[0..n] });
        defer m.deinit(al);
        try expectForm(T, a, b, m);
        var lead: [16]C64 = undefined;
        eigs(T, m, lead[0..n]);
        const want = before[n - 1];
        const got = if (want.im * lead[0].im < 0) lead[0].conjugate() else lead[0];
        try testing.expect(got.sub(want).magnitude() <= 1e-9 * (1 + want.magnitude()));
    };
}

test "ordqz: clusters, the highest first; an infinite eigenvalue is outside the disk" {
    const al = testing.allocator;
    // triangular with eigenvalues 1, 2, 3 and 4 / 0 (infinite), Q = Z = I
    const n = 4;
    const a = try Mat(f64).init(al, n, n);
    defer a.deinit(al);
    const b = try Mat(f64).init(al, n, n);
    defer b.deinit(al);
    const eye = try Mat(f64).zeros(al, n, n);
    defer eye.deinit(al);
    for (0..n) |j| for (0..n) |i| {
        a.data[i + j * n] = if (i > j) 0 else if (i == j) @floatFromInt(i + 1) else 0.5;
        b.data[i + j * n] = if (i > j) 0 else if (i == j) (if (i == 3) 0 else 1) else 0.25;
        if (i == j) eye.data[i + j * n] = 1;
    };
    const c = try ordqz(al, test_io, a, b, eye, eye, .{ .clusters = &.{ 1, 3, 1, 2 } });
    defer c.deinit(al);
    try expectForm(f64, a, b, c);
    // 2 (cluster 3), then the infinite one (2), then 1 and 3 (1)
    try testing.expectApproxEqAbs(2, c.aa.at(0, 0) / c.bb.at(0, 0), 1e-12);
    try testing.expect(@abs(c.bb.at(1, 1)) <= 1e-14 * maxAbs(f64, b.constView()));
    try testing.expectApproxEqAbs(1, c.aa.at(2, 2) / c.bb.at(2, 2), 1e-12);
    try testing.expectApproxEqAbs(3, c.aa.at(3, 3) / c.bb.at(3, 3), 1e-12);

    const u = try ordqz(al, test_io, a, b, eye, eye, .{ .region = .udo });
    defer u.deinit(al);
    // 1 is on the circle (`|λ| >= 1`): every eigenvalue is outside, order kept
    try testing.expectApproxEqAbs(1, u.aa.at(0, 0) / u.bb.at(0, 0), 1e-12);

    const wide = try Mat(f64).init(al, 2, 3);
    defer wide.deinit(al);
    try testing.expectError(error.NotSquare, ordqz(al, test_io, wide, b, eye, eye, .{ .region = .lhp }));
    try testing.expectError(error.DimensionMismatch, ordqz(al, test_io, a, b, wide, eye, .{ .region = .lhp }));
    try testing.expectError(error.DimensionMismatch, ordqz(al, test_io, a, b, eye, eye, .{ .mask = &.{true} }));
}
