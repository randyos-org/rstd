//! Tests: the real QZ (`qz_real.zig`).

const std = @import("std");
const check = @import("../../blas/_test/check.zig");
const pencil = @import("pencil.zig");
const gghrd_reduce = @import("../gghrd.zig");
const qz_real = @import("../qz_real.zig");

const testing = std.testing;
const SchurVectors = qz_real.SchurVectors;
const hgeqzReal = qz_real.hgeqzReal;

/// As `qz_complex`'s `qzCase`, for a real pencil: the real Schur form, its
/// blocks, and each complex pair conjugate and an eigenvalue of its block.
fn qzCase(comptime T: type, n: usize, zero_diag: []const usize, seed: u64, direct: bool) !void {
    const al = testing.allocator;
    const a0 = try al.alloc(T, n * n);
    defer al.free(a0);
    const b0 = try al.alloc(T, n * n);
    defer al.free(b0);
    check.fill(T, a0, seed);
    check.fill(T, b0, seed + 1);
    for (0..n) |j| for (j + 1..n) |i| {
        b0[i + j * n] = 0;
    };
    for (zero_diag) |k| b0[k + k * n] = 0;
    if (direct) for (0..n) |i| for (0..i -| 1) |j| {
        a0[i + j * n] = 0;
    };
    const h = try al.dupe(T, a0);
    defer al.free(h);
    const t = try al.dupe(T, b0);
    defer al.free(t);
    const q = try al.alloc(T, n * n);
    defer al.free(q);
    const z = try al.alloc(T, n * n);
    defer al.free(z);
    const w = try al.alloc(T, 3 * n);
    defer al.free(w);
    const ar = w[0..n];
    const ai = w[n .. 2 * n];
    const be = w[2 * n ..];
    if (!direct) try gghrd_reduce.gghrd(T, .from_identity, .from_identity, n, 0, n, h, n, t, n, q, n, z, n);
    const v: SchurVectors = if (direct) .from_identity else .accumulate;
    try hgeqzReal(check.io, T, .schur_form, v, v, n, 0, n, h, n, t, n, ar, ai, be, q, n, z, n);
    const eps = std.math.floatEps(T);
    const tol = 60 * @as(T, @floatFromInt(n * n)) * eps;
    var j: usize = 0;
    while (j < n) {
        for (j + 1..n) |i| try testing.expectEqual(0, t[i + j * n]);
        try testing.expect(t[j + j * n] >= 0);
        const pair = ai[j] != 0;
        for ((if (pair) j + 2 else j + 1)..n) |i| try testing.expectEqual(0, h[i + j * n]);
        if (!pair) {
            try testing.expectEqual(h[j + j * n], ar[j]);
            try testing.expectEqual(t[j + j * n], be[j]);
            j += 1;
            continue;
        }
        // a pair: P's block diagonal, conjugate eigenvalues, each one of the block's
        try testing.expect(j + 1 < n and ai[j] > 0 and ai[j + 1] < 0);
        try testing.expectEqual(0, t[j + (j + 1) * n]);
        const C = std.math.Complex(T);
        const l1 = C.init(ar[j] / be[j], ai[j] / be[j]);
        const l2 = C.init(ar[j + 1] / be[j + 1], ai[j + 1] / be[j + 1]);
        const mag = std.math.hypot(l1.re, l1.im);
        try testing.expect(@abs(l1.re - l2.re) + @abs(l1.im + l2.im) <= 1e3 * tol * (1 + mag));
        const d = C.init(h[j + j * n] - l1.re * t[j + j * n], -l1.im * t[j + j * n])
            .mul(C.init(h[j + 1 + (j + 1) * n] - l1.re * t[j + 1 + (j + 1) * n], -l1.im * t[j + 1 + (j + 1) * n]))
            .sub(C.init(h[j + (j + 1) * n], 0).mul(C.init(h[j + 1 + j * n], 0)));
        const scale = (@abs(h[j + j * n]) + @abs(h[j + 1 + (j + 1) * n]) + @abs(h[j + (j + 1) * n]) + @abs(h[j + 1 + j * n]) + mag * (t[j + j * n] + t[j + 1 + (j + 1) * n]));
        try testing.expect(std.math.hypot(d.re, d.im) <= 1e3 * tol * scale * scale);
        j += 2;
    }
    try testing.expect(pencil.unitary(T, n, q) <= tol);
    try testing.expect(pencil.unitary(T, n, z) <= tol);
    try testing.expect(pencil.residual(T, n, q, h, z, a0) <= tol * 4);
    try testing.expect(pencil.residual(T, n, q, t, z, b0) <= tol * 4);
    if (zero_diag.len > 0) {
        var infinite: usize = 0;
        for (be) |x| {
            if (@abs(x) <= tol) infinite += 1;
        }
        try testing.expect(infinite >= 1);
    }
}

test "hgeqzReal: the real generalized Schur form, pairs conjugate, infinite eigenvalues" {
    inline for (.{ f32, f64 }) |T| {
        for ([_]usize{ 1, 2, 3, 4, 7, 16, 30 }) |n| try qzCase(T, n, &.{}, n * 17, false);
        for ([_]bool{ false, true }) |direct| {
            try qzCase(T, 6, &.{2}, 91, direct);
            try qzCase(T, 9, &.{ 0, 4, 8 }, 92, direct);
            try qzCase(T, 5, &.{ 0, 1, 2, 3, 4 }, 93, direct);
            try qzCase(T, 7, &.{0}, 94, direct);
        }
    }
}

test "hgeqzReal: a rotation pencil's pair; errors" {
    // (R, I) with R a rotation by 0.5: eigenvalues exp(± 0.5 i)
    const c = @cos(@as(f64, 0.5));
    const sn = @sin(@as(f64, 0.5));
    var h = [_]f64{ c, sn, -sn, c };
    var t = [_]f64{ 1, 0, 0, 1 };
    var w: [6]f64 = undefined;
    var none: [1]f64 = undefined;
    try hgeqzReal(check.io, f64, .eigenvalues, .none, .none, 2, 0, 2, &h, 2, &t, 2, w[0..2], w[2..4], w[4..6], &none, 1, &none, 1);
    for (0..2) |k| {
        try testing.expectApproxEqAbs(c, w[k] / w[4 + k], 1e-14);
        try testing.expectApproxEqAbs(if (k == 0) sn else -sn, w[2 + k] / w[4 + k], 1e-14);
    }
    var bad = [_]f64{ std.math.nan(f64), 1, 1, 1 };
    var tt = [_]f64{ 1, 0, 0, 1 };
    try testing.expectError(error.NoConvergence, hgeqzReal(check.io, f64, .eigenvalues, .none, .none, 2, 0, 2, &bad, 2, &tt, 2, w[0..2], w[2..4], w[4..6], &none, 1, &none, 1));
    try testing.expectError(error.SliceTooShort, hgeqzReal(check.io, f64, .eigenvalues, .none, .none, 2, 0, 2, &h, 2, &t, 2, w[0..1], w[2..4], w[4..6], &none, 1, &none, 1));
}
