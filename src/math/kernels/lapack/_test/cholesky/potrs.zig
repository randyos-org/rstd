//! Tests of `potrs` on `potrf`'s factor.

const std = @import("std");
const types = @import("../../../blas/types.zig");
const Vtable = @import("../../../blas/Vtable.zig");
const cholesky = @import("../../cholesky.zig");
const check = @import("../../../blas/_test/check.zig");
const guard = @import("../../../blas/_test/guard.zig");
const reference = @import("../../../blas/_test/reference.zig");
const lcheck = @import("../check.zig");
const fixtures = @import("../fixtures.zig");
const ch = @import("../cholesky.zig");

const Allocator = std.mem.Allocator;
const Uplo = types.Uplo;
const Float = types.Float;
const potrs = cholesky.potrs;
const testing = std.testing;
const test_types = lcheck.test_types;
const W = lcheck.W;
const poisoned = lcheck.poisoned;
const isNanT = lcheck.isNanT;
const nan = std.math.nan(f64);
const own_sizes = ch.own_sizes;
const shrunk_sizes = ch.shrunk_sizes;

/// `potrs` on `potrf`'s factor of `p.a`, `nrhs` right-hand sides with
/// `ldb = n + 1` (the padding row NaN), then the residual `A X - B` in `W(T)`
/// within `8 (n + 2)` ulps of `|L| |L|ᴴ |X|`, and the padding untouched. With
/// `nan_at`, `B(i, j)` is NaN: column `j` must come out all NaN, the others
/// as without it, bit for bit.
pub fn solveCase(comptime T: type, al: Allocator, p: *fixtures.CholProblem(T), nrhs: usize, seed: u64, nan_at: ?[2]usize) !void {
    const n = p.n;
    const ldb = n + 1;
    const b = try al.alloc(T, ldb * nrhs);
    defer al.free(b);
    const b0 = try al.alloc(T, n * nrhs);
    defer al.free(b0);
    check.fill(T, b0, seed);
    for (0..nrhs) |j| {
        @memcpy(b[j * ldb ..][0..n], b0[j * n ..][0..n]);
        b[j * ldb + n] = poisoned(T);
    }
    var poison: guard.Poison = .{ .child = al };
    const uplo: Uplo = if (p.upper) .u else .l;
    try potrs(poison.allocator(), check.io, T, uplo, n, nrhs, p.a, p.lda, b, ldb);
    const F = Float(T);
    const eps: f128 = std.math.floatEps(F);
    const ulps: f128 = @floatFromInt(8 * (n + 2));
    const t = try al.alloc(f128, n);
    defer al.free(t);
    for (0..nrhs) |j| {
        const x = b[j * ldb ..][0..n];
        try testing.expect(isNanT(T, b[j * ldb + n]));
        // `|L|ᴴ |x|`, then `|L|` times it
        for (0..n) |k| {
            t[k] = 0;
            for (k..n) |i| t[k] += reference.abs1(reference.wide(T, p.got(i, k))) * reference.abs1(reference.wide(T, x[i]));
        }
        for (0..n) |i| {
            var s: W(T) = reference.wide(T, b0[i + j * n]);
            var mag: f128 = 0;
            for (0..n) |k| s = reference.subW(s, reference.mulW(p.aw[i + k * n], reference.wide(T, x[k])));
            for (0..i + 1) |k| mag += reference.abs1(reference.wide(T, p.got(i, k))) * t[k];
            const d = reference.abs1(s);
            if (!(d <= ulps * eps * mag)) {
                check.say("potrs {s} {s} n={d} nrhs={d}: residual ({d},{d}) {e}, allowed {e}\n", .{ @typeName(T), @tagName(uplo), n, nrhs, i, j, @as(f64, @floatCast(d)), @as(f64, @floatCast(ulps * eps * mag)) });
                return error.TestUnexpectedResult;
            }
        }
    }
    const at = nan_at orelse return;
    const bn = try al.alloc(T, ldb * nrhs);
    defer al.free(bn);
    for (0..nrhs) |j| @memcpy(bn[j * ldb ..][0..n], b0[j * n ..][0..n]);
    bn[at[0] + at[1] * ldb] = poisoned(T);
    try potrs(poison.allocator(), check.io, T, uplo, n, nrhs, p.a, p.lda, bn, ldb);
    try lcheck.expectNanColumn(T, "potrs", @tagName(uplo), n, nrhs, ldb, at, bn, b);
}

pub fn solveCases(vt: *const Vtable) anyerror!void {
    return lcheck.eachBlocking(vt, lcheck.shrink, test_types, solveOn);
}

pub fn solveOn(comptime T: type, v: *const Vtable, vi: usize) anyerror!void {
    _ = v;
    const al = testing.allocator;
    for ([_]bool{ true, false }) |upper| {
        const sizes: []const usize = if (vi == 0) &shrunk_sizes else &own_sizes;
        for (sizes) |n| {
            if (check.rig() and n > 70) continue;
            var p = try fixtures.CholProblem(T).init(al, n, upper, n * 17 + @as(u64, @intFromBool(upper)));
            defer p.deinit(al);
            _ = try p.factorIt(al);
            for ([_]usize{ 1, 3 }) |nrhs| {
                const nan_at: [2]usize = .{ (n * 7 + nrhs) % n, nrhs / 2 };
                solveCase(T, al, &p, nrhs, n + nrhs, nan_at) catch |e| {
                    check.say("{s}, blocking {s}\n", .{ @typeName(T), if (vi == 0) "shrunk" else "own" });
                    return e;
                };
            }
        }
    }
}

test "potrs: A X = B on potrf's factor, both triangles, NaN in one column of B, on every linked core" {
    try check.eachCore(solveCases);
}

test "potrs: a reversed uplo, short slices, small leading dimensions; nothing read when empty" {
    const fa = testing.failing_allocator;
    var a: [4]f64 = @splat(nan);
    var b: [4]f64 = @splat(nan);
    try testing.expectError(error.UnsupportedUplo, potrs(fa, check.io, f64, .v, 2, 2, &a, 2, &b, 2));
    try testing.expectError(error.LeadingDimensionTooSmall, potrs(fa, check.io, f64, .l, 2, 1, &a, 1, &b, 2));
    try testing.expectError(error.LeadingDimensionTooSmall, potrs(fa, check.io, f64, .l, 2, 1, &a, 2, &b, 1));
    try testing.expectError(error.LeadingDimensionTooSmall, potrs(fa, check.io, f64, .u, 0, 1, &a, 1, &b, 0));
    try testing.expectError(error.SliceTooShort, potrs(fa, check.io, f64, .u, 2, 3, &a, 2, &b, 2));
    try testing.expectError(error.SliceTooShort, potrs(fa, check.io, f64, .u, 2, 2, &a, 2, &b, std.math.maxInt(usize) / 2));
    try potrs(fa, check.io, f64, .u, 0, 2, &a, 1, &b, 1);
    try potrs(fa, check.io, f64, .l, 2, 0, &a, 2, &b, 2);
    for (b) |v| try testing.expect(std.math.isNan(v));
}
