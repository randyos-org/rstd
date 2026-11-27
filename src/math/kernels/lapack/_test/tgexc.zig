//! Tests: `tgexc`, `tgsen` and `lagv2`, on `gges`'s forms of random pencils.

const std = @import("std");
const types = @import("../../blas/types.zig");
const check = @import("../../blas/_test/check.zig");
const ops = @import("../ops.zig");
const pencil = @import("pencil.zig");
const gges = @import("../gges.zig");
const tgexc_mod = @import("../tgexc.zig");
const tgsen_mod = @import("../tgsen.zig");
const lagv2_mod = @import("../lagv2.zig");
const qz_2x2 = @import("../qz_2x2.zig");
const schur = @import("../schur.zig");
const tgsyl_tests = @import("tgsyl.zig");

const testing = std.testing;
const Float = types.Float;
const isComplex = types.isComplex;
const tgexc = tgexc_mod.tgexc;
const tgsen = tgsen_mod.tgsen;

/// A random pencil, its generalized Schur form and `Q`, `Z`.
fn Form(comptime T: type) type {
    return struct {
        n: usize,
        a0: []T,
        b0: []T,
        s: []T,
        t: []T,
        q: []T,
        z: []T,

        fn init(n: usize, seed: u64) !@This() {
            const al = testing.allocator;
            const f: @This() = .{ .n = n, .a0 = try al.alloc(T, n * n), .b0 = try al.alloc(T, n * n), .s = try al.alloc(T, n * n), .t = try al.alloc(T, n * n), .q = try al.alloc(T, n * n), .z = try al.alloc(T, n * n) };
            check.fill(T, f.a0, seed);
            check.fill(T, f.b0, seed + 1);
            @memcpy(f.s, f.a0);
            @memcpy(f.t, f.b0);
            const w = try al.alloc(T, 3 * n);
            defer al.free(w);
            if (comptime isComplex(T)) {
                try gges.ggesComplex(al, check.io, T, .compute, .compute, n, f.s, n, f.t, n, w[0..n], w[n .. 2 * n], f.q, n, f.z, n);
            } else try gges.ggesReal(al, check.io, T, .compute, .compute, n, f.s, n, f.t, n, w[0..n], w[n .. 2 * n], w[2 * n ..], f.q, n, f.z, n);
            return f;
        }

        fn deinit(f: @This()) void {
            const al = testing.allocator;
            for ([_][]T{ f.a0, f.b0, f.s, f.t, f.q, f.z }) |x| al.free(x);
        }

        /// `Q S Zᴴ = A₀`, `Q T Zᴴ = B₀`, `Q`, `Z` unitary, the forms as promised.
        fn expectValid(f: @This()) !void {
            const n = f.n;
            const tol = 200 * @as(Float(T), @floatFromInt(n * n)) * std.math.floatEps(Float(T));
            try testing.expect(pencil.unitary(T, n, f.q) <= tol);
            try testing.expect(pencil.unitary(T, n, f.z) <= tol);
            try testing.expect(pencil.residual(T, n, f.q, f.s, f.z, f.a0) <= tol * 4);
            try testing.expect(pencil.residual(T, n, f.q, f.t, f.z, f.b0) <= tol * 4);
            for (0..n) |j| for (j + 1..n) |i| {
                try testing.expect(ops.isZero(T, f.t[i + j * n]));
                const pair_row = !isComplex(T) and i == j + 1;
                if (!pair_row) try testing.expect(ops.isZero(T, f.s[i + j * n]));
            };
            // no two 2 x 2 blocks overlap
            if (!isComplex(T)) for (0..n -| 2) |j| try testing.expect(ops.isZero(T, f.s[j + 1 + j * n]) or ops.isZero(T, f.s[j + 2 + (j + 1) * n]));
        }
    };
}

/// The eigenvalues `S(k, k) / T(k, k)`, a real pair's by `xLAG2`.
fn lambdas(comptime T: type, n: usize, s: []const T, t: []const T, out: []std.math.Complex(Float(T))) void {
    const C = std.math.Complex(Float(T));
    var k: usize = 0;
    while (k < n) : (k += 1) {
        if (comptime isComplex(T)) {
            out[k] = s[k + k * n].div(t[k + k * n]);
        } else if (k + 1 < n and s[k + 1 + k * n] != 0) {
            const e = qz_2x2.eig2x2(T, .{ s[k + k * n], s[k + (k + 1) * n], s[k + 1 + k * n], s[k + 1 + (k + 1) * n] }, .{ t[k + k * n], t[k + (k + 1) * n], 0, t[k + 1 + (k + 1) * n] });
            out[k] = C.init(e.wr1 / e.scale1, e.wi / e.scale1);
            out[k + 1] = C.init(e.wr2 / e.scale2, -e.wi / e.scale2);
            k += 1;
        } else out[k] = C.init(s[k + k * n] / t[k + k * n], 0);
    }
}

/// `x` and `y` hold the same values, within `tol` relative.
fn sameSet(comptime C: type, x: []const C, y: []const C, tol: anytype) !void {
    var used: [32]bool = @splat(false);
    for (x) |u| {
        var found = false;
        for (y, 0..) |v, j| {
            if (used[j]) continue;
            if (u.sub(v).magnitude() <= tol * (1 + u.magnitude())) {
                used[j] = true;
                found = true;
                break;
            }
        }
        try testing.expect(found);
    }
}

fn moveCase(comptime T: type, n: usize, seed: u64, ifst0: usize, ilst0: usize) !void {
    const C = std.math.Complex(Float(T));
    const f = try Form(T).init(n, seed);
    defer f.deinit();
    var before: [32]C = undefined;
    lambdas(T, n, f.s, f.t, before[0..n]);
    var ifst = ifst0;
    var ilst = ilst0;
    // a real move starts at its block's first row
    if (comptime !isComplex(T)) {
        if (ifst > 0 and f.s[ifst + (ifst - 1) * n] != 0) ifst -= 1;
    }
    const moved = before[ifst];
    try tgexc(testing.allocator, check.io, T, .update, .update, n, f.s, n, f.t, n, f.q, n, f.z, n, &ifst, &ilst);
    try f.expectValid();
    var after: [32]C = undefined;
    lambdas(T, n, f.s, f.t, after[0..n]);
    const tol = 1e3 * @as(Float(T), @floatFromInt(n)) * std.math.floatEps(Float(T));
    try sameSet(C, before[0..n], after[0..n], tol);
    // the block now at `ilst` is the one moved (a pair in either order)
    var at_ilst = [_]C{after[ilst]};
    const pair = !isComplex(T) and moved.im != 0;
    if (pair) at_ilst[0] = if (after[ilst].im * moved.im < 0) after[ilst].conjugate() else after[ilst];
    try sameSet(C, &.{moved}, &at_ilst, tol);
}

test "tgexc: blocks move forward and back, Q S Zᴴ kept" {
    // the real forms hold 2 x 2 blocks, so the pair swaps run
    for ([_][2]usize{ .{ 6, 41 }, .{ 12, 97 } }) |c| {
        const f = try Form(f64).init(c[0], c[1]);
        defer f.deinit();
        var pairs: usize = 0;
        for (0..c[0] - 1) |k| pairs += @intFromBool(f.s[k + 1 + k * c[0]] != 0);
        try testing.expect(pairs > 0);
    }
    inline for (.{ f32, f64, types.Complex32, types.Complex64 }) |T| {
        for ([_][2]usize{ .{ 0, 5 }, .{ 5, 0 }, .{ 2, 4 }, .{ 4, 1 }, .{ 3, 3 } }) |m| try moveCase(T, 6, 41 + m[0] * 7 + m[1], m[0], m[1]);
        for ([_][2]usize{ .{ 0, 11 }, .{ 11, 0 }, .{ 6, 2 } }) |m| try moveCase(T, 12, 97 + m[0] + m[1], m[0], m[1]);
    }
}

fn senCase(comptime T: type, n: usize, seed: u64) !void {
    const C = std.math.Complex(Float(T));
    const f = try Form(T).init(n, seed);
    defer f.deinit();
    var before: [32]C = undefined;
    lambdas(T, n, f.s, f.t, before[0..n]);
    // select the eigenvalues inside the unit disk
    var select: [32]bool = undefined;
    for (0..n) |k| select[k] = before[k].magnitude() < 1;
    var want: [32]C = undefined;
    var m: usize = 0;
    for (0..n) |k| if (select[k]) {
        want[m] = before[k];
        m += 1;
    };
    const w = try testing.allocator.alloc(T, 3 * n);
    defer testing.allocator.free(w);
    const alpha: schur.Eigenvalues(T) = if (comptime isComplex(T)) w[0..n] else .{ .re = w[0..n], .im = w[n .. 2 * n] };
    const beta = w[2 * n ..];
    const got = (try tgsen(testing.allocator, check.io, T, .none, .update, .update, select[0..n], n, f.s, n, f.t, n, alpha, beta, f.q, n, f.z, n)).selected;
    try f.expectValid();
    var after: [32]C = undefined;
    lambdas(T, n, f.s, f.t, after[0..n]);
    const tol = 1e3 * @as(Float(T), @floatFromInt(n)) * std.math.floatEps(Float(T));
    // a real pair picked by one row counts whole: only check the count covers the picks
    try testing.expect(got >= m);
    if (got == m) try sameSet(C, want[0..m], after[0..m], tol);
    try sameSet(C, before[0..n], after[0..n], tol);
    // `alpha / beta` the form's eigenvalues, `B`'s diagonal nonnegative
    for (0..n) |k| {
        const lam = if (comptime isComplex(T)) alpha[k].div(beta[k]) else C.init(alpha.re[k] / beta[k], alpha.im[k] / beta[k]);
        try testing.expect(lam.sub(after[k]).magnitude() <= tol * (1 + after[k].magnitude()));
        if (comptime isComplex(T)) {
            try testing.expect(f.t[k + k * n].im == 0 and f.t[k + k * n].re >= 0);
        } else if (alpha.im[k] == 0) try testing.expect(f.t[k + k * n] >= 0);
    }
}

test "tgsen: the selected eigenvalues lead, B's diagonal nonnegative" {
    inline for (.{ f32, f64, types.Complex32, types.Complex64 }) |T| for ([_]usize{ 1, 2, 5, 9, 16 }) |n| try senCase(T, n, n * 13 + 5);
}

test "tgexc, tgsen: shape and row errors" {
    const al = testing.allocator;
    var a = [_]f64{ 1, 0, 0, 2 };
    var b = [_]f64{ 1, 0, 0, 1 };
    var none: [1]f64 = undefined;
    var ifst: usize = 2;
    var ilst: usize = 0;
    try testing.expectError(error.RowOutOfRange, tgexc(al, check.io, f64, .none, .none, 2, &a, 2, &b, 2, &none, 1, &none, 1, &ifst, &ilst));
    ifst = 0;
    try testing.expectError(error.SliceTooShort, tgexc(al, check.io, f64, .update, .none, 2, &a, 2, &b, 2, &none, 2, &none, 1, &ifst, &ilst));
    try testing.expectError(error.LeadingDimensionTooSmall, tgexc(al, check.io, f64, .none, .none, 2, &a, 1, &b, 2, &none, 1, &none, 1, &ifst, &ilst));
    var w: [6]f64 = undefined;
    try testing.expectError(error.SliceTooShort, tgsen(al, check.io, f64, .none, .none, .none, &.{true}, 2, &a, 2, &b, 2, .{ .re = w[0..2], .im = w[2..4] }, w[4..6], &none, 1, &none, 1));
    // a diagonal pencil: swapping 1 and 2/1 puts 2 first
    ifst = 1;
    ilst = 0;
    try tgexc(al, check.io, f64, .none, .none, 2, &a, 2, &b, 2, &none, 1, &none, 1, &ifst, &ilst);
    try testing.expectApproxEqAbs(@as(f64, 2), a[0] / b[0], 1e-15);
}

test "lagv2: real eigenvalues give triangles, complex ones a diagonal B" {
    // `[1 2; 3 4]` over `I` has real eigenvalues; a rotation pencil complex ones
    var a = [_]f64{ 1, 3, 2, 4 };
    var b = [_]f64{ 2, 0, 1, 3 };
    const fa = std.math.hypot(std.math.hypot(@as(f64, 1), 3), std.math.hypot(@as(f64, 2), 4));
    const r = lagv2_mod.lagv2(f64, &a, 2, &b, 2);
    try testing.expectEqual(0, a[1]);
    try testing.expectEqual(0, b[1]);
    try testing.expectApproxEqAbs(fa, std.math.hypot(std.math.hypot(a[0], a[1]), std.math.hypot(a[2], a[3])), 1e-14);
    try testing.expectEqual(0, r.alphai[0]);
    var c = [_]f64{ 0, -1, 1, 0 };
    var d = [_]f64{ 1, 0, 0.5, 2 };
    const s = lagv2_mod.lagv2(f64, &c, 2, &d, 2);
    try testing.expectEqual(0, d[1]);
    try testing.expectEqual(0, d[2]);
    try testing.expect(s.alphai[0] > 0 and s.alphai[1] == -s.alphai[0]);
    // `det(A - λ B) = 0` at `λ = alpha / beta`, on the original pencil
    const lr = s.alphar[0] / s.beta[0];
    const li = s.alphai[0] / s.beta[0];
    // `[-λ, 1 - 0.5 λ; -1, -2 λ]`: `2 λ² + 1 - 0.5 λ = 0`
    const re = 2 * (lr * lr - li * li) + 1 - 0.5 * lr;
    const im = 2 * (2 * lr * li) - 0.5 * li;
    try testing.expect(@abs(re) + @abs(im) < 1e-13);
}

fn conditionCase(comptime T: type, n: usize, seed: u64, job: tgsen_mod.TgsenJob) !void {
    const al = testing.allocator;
    const f = try Form(T).init(n, seed);
    defer f.deinit();
    var select: [32]bool = undefined;
    for (0..n) |k| select[k] = k % 3 == 1;
    const w = try al.alloc(T, 3 * n);
    defer al.free(w);
    const alpha: schur.Eigenvalues(T) = if (comptime isComplex(T)) w[0..n] else .{ .re = w[0..n], .im = w[n .. 2 * n] };
    const r = try tgsen(al, check.io, T, job, .update, .update, select[0..n], n, f.s, n, f.t, n, alpha, w[2 * n ..], f.q, n, f.z, n);
    try f.expectValid();
    const m = r.selected;
    try testing.expect(m > 0 and m < n);
    const pq = r.projections orelse return error.TestUnexpectedResult;
    for (pq) |x| try testing.expect(x > 0 and x <= 1);
    const sep = r.separations orelse return error.TestUnexpectedResult;
    // `Dif_u` and `Dif_l` against `σ_min` of their Kronecker forms
    const m2 = n - m;
    const s = f.s;
    const t = f.t;
    const sig = [2]f64{
        try tgsyl_tests.kroneckerSigmaMin(T, al, m, m2, s, n, s[m + m * n ..], n, t, n, t[m + m * n ..], n),
        try tgsyl_tests.kroneckerSigmaMin(T, al, m2, m, s[m + m * n ..], n, s, n, t[m + m * n ..], n, t, n),
    };
    const lim = 4 * @sqrt(2 * @as(f64, @floatFromInt(m * m2)));
    for (sep, sig) |d, sg| {
        const dd: f64 = @floatCast(d);
        if (!(dd >= sg * (1 - 1e-3) / lim and dd <= sg * lim)) {
            check.say("tgsen {s} n={d} {t}: dif {e}, sigma_min {e}\n", .{ @typeName(T), n, job, dd, sg });
            return error.TestUnexpectedResult;
        }
    }
}

test "tgsen: projections in (0, 1], separations near σ_min of the Kronecker forms" {
    inline for (.{ f32, f64, types.Complex32, types.Complex64 }) |T| for ([_]usize{ 4, 7 }) |n| {
        try conditionCase(T, n, n * 7 + 3, .both_frobenius);
        try conditionCase(T, n, n * 7 + 3, .both_one_norm);
    };
}

test "tgsen: nothing moved, projections 1 and separations the pencil's norm" {
    const al = testing.allocator;
    const n = 5;
    const f = try Form(f64).init(n, 77);
    defer f.deinit();
    var fro: f64 = 0;
    for (f.s, f.t) |x, y| fro += x * x + y * y;
    var w: [3 * n]f64 = undefined;
    var none: [1]f64 = undefined;
    const r = try tgsen(al, check.io, f64, .both_one_norm, .none, .none, &@as([n]bool, @splat(false)), n, f.s, n, f.t, n, .{ .re = w[0..n], .im = w[n .. 2 * n] }, w[2 * n ..], &none, 1, &none, 1);
    try testing.expectEqual(0, r.selected);
    try testing.expectEqual([2]f64{ 1, 1 }, r.projections.?);
    try testing.expectApproxEqRel(@sqrt(fro), r.separations.?[0], 1e-14);
    try testing.expectEqual(r.separations.?[0], r.separations.?[1]);
}
