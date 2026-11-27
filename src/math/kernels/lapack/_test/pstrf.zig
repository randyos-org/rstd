//! Tests: `pstrf`.

const std = @import("std");
const types = @import("../../blas/types.zig");
const check = @import("../../blas/_test/check.zig");
const ops = @import("../ops.zig");
const pstrf_mod = @import("../pstrf.zig");

const testing = std.testing;
const Float = types.Float;
const pstrf = pstrf_mod.pstrf;

/// `A = B Bᴴ` (`n x n`, rank `min(r, n)`), packed.
fn semidefinite(comptime T: type, al: std.mem.Allocator, n: usize, r: usize, seed: u64) ![]T {
    const b = try al.alloc(T, n * r);
    defer al.free(b);
    check.fill(T, b, seed);
    const a = try al.alloc(T, n * n);
    for (0..n) |j| for (0..n) |i| {
        var s = ops.zero(T);
        for (0..r) |l| s = ops.add(T, s, ops.mul(T, b[i + l * n], ops.conj(T, b[j + l * n])));
        a[i + j * n] = s;
    };
    return a;
}

/// `max |(Pᵀ A P)(i, j) - (Fᴴ F)(i, j)|` over the factor's first `rank`
/// rows (`.u`) or columns (`.l`), `F` taken from `f`'s `uplo` triangle.
fn residual(comptime T: type, n: usize, a: []const T, f: []const T, upper: bool, piv: []const usize, rank: usize) Float(T) {
    var worst: Float(T) = 0;
    for (0..n) |i| for (0..n) |j| {
        var s = ops.zero(T);
        for (0..@min(rank, i + 1, j + 1)) |l| {
            // `U(l, i)`, as `conj(L(i, l))` when lower
            const ui = if (upper) f[l + i * n] else ops.conj(T, f[i + l * n]);
            const uj = if (upper) f[l + j * n] else ops.conj(T, f[j + l * n]);
            s = ops.add(T, s, ops.mul(T, ops.conj(T, ui), uj));
        }
        worst = @max(worst, ops.abs(T, ops.sub(T, s, a[piv[i] + piv[j] * n])));
    };
    return worst;
}

fn caseOf(comptime T: type, n: usize, r: usize, upper: bool, seed: u64) !void {
    const al = testing.allocator;
    const F = Float(T);
    const a = try semidefinite(T, al, n, r, seed);
    defer al.free(a);
    const f = try al.dupe(T, a);
    defer al.free(f);
    const piv = try al.alloc(usize, n);
    defer al.free(piv);
    const rank = try pstrf(al, check.io, T, if (upper) .u else .l, n, f, @max(n, 1), piv, -1);
    try testing.expectEqual(@min(r, n), rank);
    // a permutation, and nonincreasing pivots
    var seen = try al.alloc(bool, n);
    defer al.free(seen);
    @memset(seen, false);
    for (piv) |p| {
        try testing.expect(!seen[p]);
        seen[p] = true;
    }
    for (1..rank) |k| try testing.expect(ops.re(T, f[k + k * n]) <= ops.re(T, f[k - 1 + (k - 1) * n]) * (1 + 1e-3));
    var big: F = 0;
    for (a) |x| big = @max(big, ops.abs(T, x));
    const tol = 100 * @as(F, @floatFromInt(n + 1)) * std.math.floatEps(F) * @max(big, 1);
    const res = residual(T, n, a, f, upper, piv, rank);
    if (!(res <= tol)) {
        check.say("pstrf {s} n={d} r={d} upper={}: residual {e}, tol {e}\n", .{ @typeName(T), n, r, upper, res, tol });
        return error.TestUnexpectedResult;
    }
}

test "pstrf: Pᵀ A P = Uᴴ U or L Lᴴ, full rank and deficient, around the block" {
    inline for (.{ f64, f32, types.Complex64, types.Complex32 }) |T| for ([_]bool{ true, false }) |upper| {
        for ([_]usize{ 1, 2, 5, 63, 64, 65, 130 }) |n| try caseOf(T, n, n + 3, upper, n * 7 + 1);
        for ([_][2]usize{ .{ 6, 3 }, .{ 70, 20 }, .{ 130, 100 } }) |c| try caseOf(T, c[0], c[1], upper, c[0] + c[1]);
    };
}

test "pstrf: a zero matrix has rank 0; tol stops early; bad shapes" {
    const al = testing.allocator;
    var z = [_]f64{ 0, 0, 0, 0 };
    var piv: [2]usize = undefined;
    try testing.expectEqual(0, try pstrf(al, check.io, f64, .u, 2, &z, 2, &piv, -1));
    // diag(4, 1): the second pivot 1 is at tol 2, so the rank is 1
    var d = [_]f64{ 1, 0, 0, 4 };
    try testing.expectEqual(1, try pstrf(al, check.io, f64, .l, 2, &d, 2, &piv, 2));
    try testing.expectEqual(1, piv[0]);
    try testing.expectEqual(2, d[0]);
    try testing.expectError(error.LeadingDimensionTooSmall, pstrf(al, check.io, f64, .u, 2, &z, 1, &piv, -1));
    try testing.expectError(error.SliceTooShort, pstrf(al, check.io, f64, .u, 2, &z, 2, piv[0..1], -1));
    try testing.expectEqual(0, try pstrf(al, check.io, f64, .u, 0, z[0..0], 1, piv[0..0], -1));
}
