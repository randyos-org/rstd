//! Tests of `../cholesky.zig`.

const std = @import("std");
const sysinfo = @import("builtin");
const blas = @import("../../blas/__root__.zig");
const types = @import("../../blas/types.zig");
const accel = @import("../../blas/apple_accel.zig");
const Vtable = @import("../../blas/Vtable.zig");
const ops = @import("../ops.zig");
const cholesky = @import("../cholesky.zig");
const check = @import("../../blas/_test/check.zig");
const guard = @import("../../blas/_test/guard.zig");
const lcheck = @import("check.zig");
const fixtures = @import("fixtures.zig");

const Uplo = types.Uplo;
const isComplex = types.isComplex;
const scalar = ops.scalar;
const re = ops.re;
const potf2 = cholesky.potf2;
const potrf = cholesky.potrf;
const potrfChecked = cholesky.potrfChecked;
const testing = std.testing;
const test_types = lcheck.test_types;
const fromF64 = lcheck.fromF64;
const nan = std.math.nan(f64);

/// Sizes around every block edge: on the shrunk blocking, and on a core's
/// own (`dtb_entries / 2` is 16 to 128 by core, so 130 is blocked on every
/// core, twice over where it is 32).
pub const shrunk_sizes = [_]usize{ 1, 2, 3, 4, 5, 7, 12, 13, 14, 17, 26 };
pub const own_sizes = [_]usize{ 1, 31, 33, 70, 130 };

pub fn factorCases(vt: *const Vtable) anyerror!void {
    return lcheck.eachBlocking(vt, lcheck.shrink, test_types, factorOn);
}

pub fn factorOn(comptime T: type, v: *const Vtable, vi: usize) anyerror!void {
    _ = v;
    const al = testing.allocator;
    for ([_]bool{ true, false }) |upper| {
        const sizes: []const usize = if (vi == 0) &shrunk_sizes else &own_sizes;
        for (sizes) |n| {
            if (check.rig() and n > 70) continue;
            var p = try fixtures.CholProblem(T).init(al, n, upper, n * 31 + @as(u64, @intFromBool(upper)));
            defer p.deinit(al);
            var poison: guard.Poison = .{ .child = al };
            const failed = p.factorIt(poison.allocator()) catch |e| {
                check.say("potrf {s} {s} n={d}\n", .{ @typeName(T), if (upper) "u" else "l", n });
                return e;
            };
            try testing.expectEqual(@as(?usize, null), failed);
            p.expectFactor() catch |e| {
                check.say("{s}, blocking {s}\n", .{ @typeName(T), if (vi == 0) "shrunk" else "own" });
                return e;
            };
        }
    }
}

test "potrf: A = L Lᴴ from a known L, both triangles, on every linked core" {
    try check.eachCore(factorCases);
}

/// `potrf` on `p` must fail at pivot `want` (0-based).
pub fn expectFailed(comptime T: type, p: *fixtures.CholProblem(T), want: usize, what: []const u8) !void {
    var poison: guard.Poison = .{ .child = testing.allocator };
    const got = try p.factorIt(poison.allocator());
    if (got != want) {
        check.say("{s} {s} {s} n={d}: failed at pivot {?d}, want {d}\n", .{ what, @typeName(T), if (p.upper) "u" else "l", p.n, got, want });
        return error.TestUnexpectedResult;
    }
}

/// Positions of a bad element: the first, last and either side of block
/// edges on both blockings.
pub const bad_at = [_]usize{ 0, 1, 2, 3, 6, 12, 13, 25, 32, 33, 69 };

pub fn failCases(vt: *const Vtable) anyerror!void {
    return lcheck.eachBlocking(vt, lcheck.shrink, test_types, failOn);
}

pub fn failOn(comptime T: type, v: *const Vtable, vi: usize) anyerror!void {
    _ = v;
    const al = testing.allocator;
    const inf = std.math.inf(f64);
    for ([_]bool{ true, false }) |upper| {
        const n: usize = if (vi == 0) 26 else 70;
        var p = try fixtures.CholProblem(T).init(al, n, upper, 7 + @as(u64, @intFromBool(upper)));
        defer p.deinit(al);
        for (bad_at) |k| {
            if (k >= n) continue;
            // a pivot made negative, NaN, or a NaN or infinity left of it
            for ([_]f64{ -100, nan }) |x| {
                p.store();
                p.a[k * (p.lda + 1)] = fromF64(T, x);
                try expectFailed(T, &p, k, "diagonal");
            }
            if (k == 0) continue;
            for ([_]f64{ nan, inf, -inf }) |x| for ([_]usize{ 0, k / 2, k - 1 }) |c| {
                p.store();
                // `(k, c)` in the lower triangle, `(c, k)` in the upper
                p.a[if (upper) c + k * p.lda else k + c * p.lda] = fromF64(T, x);
                try expectFailed(T, &p, k, "off the diagonal");
            };
        }
    }
}

test "potrf: the first minor that is not positive definite, by a pivot, NaN or infinity, on every linked core" {
    try check.eachCore(failCases);
}

test "potrf: a failed pivot is stored, real; potrfChecked returns its index" {
    // rstd's own (Accelerate's keeps a complex pivot's imaginary part), after
    // init, which turns routing on
    _ = blas.vtable();
    const was = accel.routing();
    _ = accel.route(false);
    defer _ = accel.route(was);
    inline for (.{ f64, types.Complex64 }) |T| for ([_]Uplo{ .u, .l }) |uplo| {
        // the identity with -3 at (2, 2): exact pivots
        var a: [16]T = @splat(scalar(T, 0));
        for (0..4) |i| a[i * 5] = scalar(T, 1);
        a[10] = if (comptime isComplex(T)) .init(-3, 5) else -3;
        const al = testing.allocator;
        try testing.expectError(error.NotPositiveDefinite, potrf(al, check.io, T, uplo, 4, &a, 4));
        try testing.expectEqual(scalar(T, -3), a[10]);
        try testing.expectEqual(scalar(T, 1), a[5]);
        // a zero pivot fails too
        a[10] = scalar(T, 0);
        try testing.expectEqual(@as(?usize, 2), try potrfChecked(al, check.io, T, uplo, 4, &a, 4));
        a[10] = scalar(T, 4);
        try testing.expectEqual(@as(?usize, null), try potrfChecked(al, check.io, T, uplo, 4, &a, 4));
        try testing.expectEqual(scalar(T, 2), a[10]);
    };
}

test "potrf: a reversed uplo, a short slice, a small leading dimension; nothing read for n = 0" {
    const fa = testing.failing_allocator;
    var a: [4]f64 = @splat(nan);
    try testing.expectError(error.UnsupportedUplo, potrf(fa, check.io, f64, .v, 2, &a, 2));
    try testing.expectError(error.UnsupportedUplo, potrf(fa, check.io, f64, .m, 2, &a, 2));
    try testing.expectError(error.LeadingDimensionTooSmall, potrf(fa, check.io, f64, .l, 2, &a, 1));
    try testing.expectError(error.LeadingDimensionTooSmall, potrf(fa, check.io, f64, .u, 0, &a, 0));
    try testing.expectError(error.SliceTooShort, potrf(fa, check.io, f64, .u, 3, &a, 3));
    try testing.expectError(error.SliceTooShort, potrf(fa, check.io, f64, .u, 2, &a, std.math.maxInt(usize) / 2));
    try potrf(fa, check.io, f64, .u, 0, &a, 1);
    for (a) |v| try testing.expect(std.math.isNan(v));
}

test "potf2: the unblocked factor, a failing pivot, bad arguments" {
    const al = testing.allocator;
    inline for (test_types) |T| for ([_]bool{ true, false }) |upper| {
        const n = 26;
        var p = try fixtures.CholProblem(T).init(al, n, upper, 5);
        defer p.deinit(al);
        try potf2(T, if (upper) .u else .l, n, p.a, p.lda);
        try p.expectFactor();
        p.store();
        p.a[7 * (p.lda + 1)] = fromF64(T, nan);
        try testing.expectError(error.NotPositiveDefinite, potf2(T, if (upper) .u else .l, n, p.a, p.lda));
        // the factor's diagonal positive before the NaN, the failed pivot at it
        for (0..8) |i| try testing.expectEqual(i < 7, re(T, p.a[i * (p.lda + 1)]) > 0);
    };
    var a: [4]f64 = @splat(nan);
    try testing.expectError(error.UnsupportedUplo, potf2(f64, .v, 2, &a, 2));
    try testing.expectError(error.LeadingDimensionTooSmall, potf2(f64, .l, 2, &a, 1));
    try testing.expectError(error.SliceTooShort, potf2(f64, .l, 3, &a, 3));
    try potf2(f64, .u, 0, &a, 1);
    for (a) |v| try testing.expect(std.math.isNan(v));
}

test "potrf goes to Accelerate where it routes, and agrees with rstd's" {
    _ = blas.coreName();
    if (!accel.routing()) return error.SkipZigTest;
    const al = testing.allocator;
    inline for (.{ f32, f64, types.Complex32, types.Complex64 }) |T| for ([_]bool{ true, false }) |upper| {
        for ([_]usize{ 1, 7, 70, 130 }) |n| {
            var p = try fixtures.CholProblem(T).init(al, n, upper, n * 3 + @as(u64, @intFromBool(upper)));
            defer p.deinit(al);
            const before = accel.calls.load(.monotonic);
            try testing.expectEqual(@as(?usize, null), try p.factorIt(al));
            // x86_64's complex `potrf` stays on rstd (`apple_accel.zig`)
            const routed = !(sysinfo.target.cpu.arch == .x86_64 and isComplex(T));
            try testing.expectEqual(routed, accel.calls.load(.monotonic) > before);
            p.expectFactor() catch |e| {
                check.say("Accelerate's {s}potrf_\n", .{@typeName(T)});
                return e;
            };
            if (n < 7) continue;
            p.store();
            p.a[(n / 2) * (p.lda + 1)] = fromF64(T, nan);
            try expectFailed(T, &p, n / 2, "Accelerate's NaN pivot");
            // finite values in the other triangle come back as they were (issue #537)
            p.store();
            for (0..n) |j| for (0..n) |i| {
                if (if (upper) i > j else i < j) p.a[i + j * p.lda] = fromF64(T, 777);
            };
            try testing.expectEqual(@as(?usize, null), try p.factorIt(al));
            for (0..n) |j| for (0..n) |i| {
                if (if (upper) i > j else i < j) try testing.expectEqual(fromF64(T, 777), p.a[i + j * p.lda]);
            };
        }
    };
}
