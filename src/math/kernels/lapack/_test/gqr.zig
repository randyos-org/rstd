//! Tests: `ggqrf` and `ggrqf`.

const std = @import("std");
const types = @import("../../blas/types.zig");
const check = @import("../../blas/_test/check.zig");
const lcheck = @import("check.zig");
const ops = @import("../ops.zig");
const qr = @import("../qr.zig");
const rq = @import("../rq.zig");
const gqr = @import("../gqr.zig");

const testing = std.testing;
const Allocator = std.mem.Allocator;
const Mat = ops.Mat;
const Float = types.Float;

fn maxAbs(comptime T: type, a: Mat(T)) Float(T) {
    var big: Float(T) = 0;
    for (0..a.cols) |j| for (0..a.rows) |i| {
        big = @max(big, ops.abs(T, a.get(i, j)));
    };
    return big;
}

/// `Q Qᴴ = I` for the square `q`.
fn expectUnitary(comptime T: type, al: Allocator, what: []const u8, q: Mat(T)) !void {
    const n = q.rows;
    const g = try ops.alloc(T, al, n, n);
    defer ops.free(T, al, g);
    lcheck.mul(T, false, q, true, q, g);
    const id = try ops.alloc(T, al, n, n);
    defer ops.free(T, al, id);
    id.identity();
    try lcheck.expectNear(T, what, g, id, lcheck.thresh * @as(Float(T), @floatFromInt(@max(n, 1))) * std.math.floatEps(Float(T)));
}

/// The `r x c` factor in `f`: on and above the diagonal (`upper`), or the
/// trapezoid ending at `(r - 1, c - 1)`; zero elsewhere.
fn trapezoid(comptime T: type, f: Mat(T), upper: bool, out: Mat(T)) void {
    for (0..f.cols) |j| for (0..f.rows) |i| {
        const in = if (upper) i <= j else j + f.rows >= i + f.cols;
        out.at(i, j).* = if (in) f.get(i, j) else ops.zero(T);
    };
}

/// `Q` (`n x n`) from `geqrf`'s `k` reflections in `f`'s columns.
fn qOf(comptime T: type, al: Allocator, f: Mat(T), k: usize, tau: []const T) !Mat(T) {
    const n = f.rows;
    const q = try ops.alloc(T, al, n, n);
    q.fill(ops.zero(T));
    q.sub(0, 0, n, k).copyFrom(f.sub(0, 0, n, k));
    try qr.orgqr(al, check.io, T, n, n, k, q.slice(), n, tau);
    return q;
}

/// `Q` (`n x n`) from `gerqf`'s `k` reflections in `f`'s last rows.
fn rqOf(comptime T: type, al: Allocator, f: Mat(T), k: usize, tau: []const T) !Mat(T) {
    const n = f.cols;
    const q = try ops.alloc(T, al, n, n);
    q.fill(ops.zero(T));
    q.sub(n - k, 0, k, n).copyFrom(f.sub(f.rows - k, 0, k, n));
    try rq.orgrq(al, check.io, T, n, n, k, q.slice(), n, tau);
    return q;
}

/// `dgqrts.f`'s ratios: `‖R - Qᴴ A‖`, `‖T - Qᴴ B Zᴴ‖` and the
/// orthogonality of `Q` and `Z`.
fn ggqrfCase(comptime T: type, al: Allocator, n: usize, m: usize, p: usize, seed: u64) !void {
    const F = Float(T);
    const eps = std.math.floatEps(F);
    const a0 = try ops.alloc(T, al, n, m);
    defer ops.free(T, al, a0);
    check.fill(T, a0.slice(), seed);
    const b0 = try ops.alloc(T, al, n, p);
    defer ops.free(T, al, b0);
    check.fill(T, b0.slice(), seed + 1);
    const a = try ops.alloc(T, al, n, m);
    defer ops.free(T, al, a);
    a.copyFrom(a0);
    const b = try ops.alloc(T, al, n, p);
    defer ops.free(T, al, b);
    b.copyFrom(b0);
    const taua = try al.alloc(T, @max(@min(n, m), 1));
    defer al.free(taua);
    const taub = try al.alloc(T, @max(@min(n, p), 1));
    defer al.free(taub);
    try gqr.ggqrf(al, check.io, T, n, m, p, a.slice(), a.ld, taua, b.slice(), b.ld, taub);
    const q = try qOf(T, al, a, @min(n, m), taua);
    defer ops.free(T, al, q);
    const z = try rqOf(T, al, b, @min(n, p), taub);
    defer ops.free(T, al, z);
    try expectUnitary(T, al, "ggqrf Q", q);
    try expectUnitary(T, al, "ggqrf Z", z);
    const big: F = @floatFromInt(@max(n, m, p, 1));
    // R = Qᴴ A
    const r = try ops.alloc(T, al, n, m);
    defer ops.free(T, al, r);
    lcheck.mul(T, true, q, false, a0, r);
    const want_r = try ops.alloc(T, al, n, m);
    defer ops.free(T, al, want_r);
    trapezoid(T, a, true, want_r);
    try lcheck.expectNear(T, "ggqrf R", r, want_r, lcheck.thresh * big * eps * @max(maxAbs(T, a0), 1));
    // T = Qᴴ B Zᴴ
    const qb = try ops.alloc(T, al, n, p);
    defer ops.free(T, al, qb);
    lcheck.mul(T, true, q, false, b0, qb);
    const t = try ops.alloc(T, al, n, p);
    defer ops.free(T, al, t);
    lcheck.mul(T, false, qb, true, z, t);
    const want_t = try ops.alloc(T, al, n, p);
    defer ops.free(T, al, want_t);
    trapezoid(T, b, false, want_t);
    try lcheck.expectNear(T, "ggqrf T", t, want_t, lcheck.thresh * big * eps * @max(maxAbs(T, b0), 1));
}

/// `dgrqts.f`'s ratios: `‖R - A Qᴴ‖`, `‖T - Zᴴ B Qᴴ‖` and the
/// orthogonality of `Q` and `Z`.
fn ggrqfCase(comptime T: type, al: Allocator, m: usize, p: usize, n: usize, seed: u64) !void {
    const F = Float(T);
    const eps = std.math.floatEps(F);
    const a0 = try ops.alloc(T, al, m, n);
    defer ops.free(T, al, a0);
    check.fill(T, a0.slice(), seed);
    const b0 = try ops.alloc(T, al, p, n);
    defer ops.free(T, al, b0);
    check.fill(T, b0.slice(), seed + 1);
    const a = try ops.alloc(T, al, m, n);
    defer ops.free(T, al, a);
    a.copyFrom(a0);
    const b = try ops.alloc(T, al, p, n);
    defer ops.free(T, al, b);
    b.copyFrom(b0);
    const taua = try al.alloc(T, @max(@min(m, n), 1));
    defer al.free(taua);
    const taub = try al.alloc(T, @max(@min(p, n), 1));
    defer al.free(taub);
    try gqr.ggrqf(al, check.io, T, m, p, n, a.slice(), a.ld, taua, b.slice(), b.ld, taub);
    const q = try rqOf(T, al, a, @min(m, n), taua);
    defer ops.free(T, al, q);
    const z = try qOf(T, al, b, @min(p, n), taub);
    defer ops.free(T, al, z);
    try expectUnitary(T, al, "ggrqf Q", q);
    try expectUnitary(T, al, "ggrqf Z", z);
    const big: F = @floatFromInt(@max(n, m, p, 1));
    // R = A Qᴴ
    const r = try ops.alloc(T, al, m, n);
    defer ops.free(T, al, r);
    lcheck.mul(T, false, a0, true, q, r);
    const want_r = try ops.alloc(T, al, m, n);
    defer ops.free(T, al, want_r);
    trapezoid(T, a, false, want_r);
    try lcheck.expectNear(T, "ggrqf R", r, want_r, lcheck.thresh * big * eps * @max(maxAbs(T, a0), 1));
    // T = Zᴴ B Qᴴ
    const zb = try ops.alloc(T, al, p, n);
    defer ops.free(T, al, zb);
    lcheck.mul(T, true, z, false, b0, zb);
    const t = try ops.alloc(T, al, p, n);
    defer ops.free(T, al, t);
    lcheck.mul(T, false, zb, true, q, t);
    const want_t = try ops.alloc(T, al, p, n);
    defer ops.free(T, al, want_t);
    trapezoid(T, b, true, want_t);
    try lcheck.expectNear(T, "ggrqf T", t, want_t, lcheck.thresh * big * eps * @max(maxAbs(T, b0), 1));
}

const shapes = [_][3]usize{ .{ 1, 1, 1 }, .{ 5, 3, 4 }, .{ 3, 5, 4 }, .{ 4, 4, 6 }, .{ 6, 4, 2 }, .{ 4, 2, 6 }, .{ 40, 30, 70 }, .{ 70, 40, 30 } };

test "lapack ggqrf: R = Qᴴ A, T = Qᴴ B Zᴴ, Q and Z unitary, every shape" {
    const al = testing.allocator;
    inline for (lcheck.test_types) |T| for (shapes, 0..) |s, i| {
        try ggqrfCase(T, al, s[0], s[1], s[2], 31 * i + 3);
    };
}

test "lapack ggrqf: R = A Qᴴ, T = Zᴴ B Qᴴ, Q and Z unitary, every shape" {
    const al = testing.allocator;
    inline for (lcheck.test_types) |T| for (shapes, 0..) |s, i| {
        try ggrqfCase(T, al, s[0], s[1], s[2], 37 * i + 5);
    };
}

test "lapack ggqrf, ggrqf: empty; bad arguments" {
    const al = testing.allocator;
    var a: [6]f64 = @splat(1);
    var b: [6]f64 = @splat(1);
    var ta: [2]f64 = undefined;
    var tb: [2]f64 = undefined;
    try gqr.ggqrf(al, check.io, f64, 0, 2, 2, &a, 1, &ta, &b, 1, &tb);
    try gqr.ggrqf(al, check.io, f64, 2, 2, 0, &a, 2, &ta, &b, 2, &tb);
    try gqr.ggrqf(al, check.io, f64, 0, 2, 3, &a, 1, &ta, &b, 2, &tb);
    try testing.expectError(error.LeadingDimensionTooSmall, gqr.ggqrf(al, check.io, f64, 2, 2, 2, &a, 1, &ta, &b, 2, &tb));
    try testing.expectError(error.SliceTooShort, gqr.ggqrf(al, check.io, f64, 2, 2, 4, &a, 2, &ta, &b, 2, &tb));
    try testing.expectError(error.SliceTooShort, gqr.ggqrf(al, check.io, f64, 2, 2, 2, &a, 2, ta[0..1], &b, 2, &tb));
    try testing.expectError(error.SliceTooShort, gqr.ggrqf(al, check.io, f64, 2, 2, 2, &a, 2, &ta, &b, 2, tb[0..1]));
    try testing.expectError(error.LeadingDimensionTooSmall, gqr.ggrqf(al, check.io, f64, 2, 3, 2, &a, 2, &ta, &b, 2, &tb));
}
