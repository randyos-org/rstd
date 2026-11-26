//! Tests of `potri` on `potrf`'s factor.

const std = @import("std");
const types = @import("../../../blas/types.zig");
const Vtable = @import("../../../blas/Vtable.zig");
const ops = @import("../../ops.zig");
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
const scalar = ops.scalar;
const conj = ops.conj;
const potri = cholesky.potri;
const testing = std.testing;
const test_types = lcheck.test_types;
const W = lcheck.W;
const isNanT = lcheck.isNanT;
const nan = std.math.nan(f64);
const own_sizes = ch.own_sizes;
const shrunk_sizes = ch.shrunk_sizes;

/// `potri` on `potrf`'s factor of `p`: `A X = I` within `32 (n + 2)` ulps
/// of `|A| |X|`, `X` Hermitian from its stored triangle, nothing outside
/// it written.
pub fn inverseCase(comptime T: type, al: Allocator, p: *fixtures.CholProblem(T)) !void {
    const n = p.n;
    _ = try p.factorIt(al);
    var poison: guard.Poison = .{ .child = al };
    try potri(poison.allocator(), check.io, T, if (p.upper) .u else .l, n, p.a, p.lda);
    const eps: f128 = std.math.floatEps(Float(T));
    const ulps: f128 = @floatFromInt(32 * (n + 2));
    for (0..n) |j| for (0..n) |i| {
        var s: W(T) = reference.wide(T, scalar(T, if (i == j) -1 else 0));
        var mag: f128 = 0;
        for (0..n) |k| {
            const inside = if (p.upper) k <= j else k >= j;
            const x = reference.wide(T, if (inside) p.a[k + j * p.lda] else conj(T, p.a[j + k * p.lda]));
            s = reference.addW(s, reference.mulW(p.aw[i + k * n], x));
            mag += reference.abs1(p.aw[i + k * n]) * reference.abs1(x);
        }
        const d = reference.abs1(s);
        if (!(d <= ulps * eps * mag)) {
            check.say("potri {s} {s} n={d}: (A X - I)({d},{d}) = {e}, allowed {e}\n", .{ @typeName(T), if (p.upper) "u" else "l", n, i, j, @as(f64, @floatCast(d)), @as(f64, @floatCast(ulps * eps * mag)) });
            return error.TestUnexpectedResult;
        }
    };
    for (0..n) |j| for (0..p.lda) |i| {
        const stored = i < n and (if (p.upper) i <= j else i >= j);
        if (!stored and !isNanT(T, p.a[i + j * p.lda])) {
            check.say("potri {s} n={d}: wrote ({d},{d}), outside its triangle\n", .{ @typeName(T), n, i, j });
            return error.TestUnexpectedResult;
        }
    };
}

pub fn inverseCases(vt: *const Vtable) anyerror!void {
    return lcheck.eachBlocking(vt, lcheck.shrink, test_types, inverseOn);
}

pub fn inverseOn(comptime T: type, v: *const Vtable, vi: usize) anyerror!void {
    _ = v;
    const al = testing.allocator;
    for ([_]bool{ true, false }) |upper| {
        const sizes: []const usize = if (vi == 0) &shrunk_sizes else &own_sizes;
        for (sizes) |n| {
            if (check.rig() and n > 70) continue;
            var p = try fixtures.CholProblem(T).init(al, n, upper, n * 11 + @as(u64, @intFromBool(upper)));
            defer p.deinit(al);
            inverseCase(T, al, &p) catch |e| {
                check.say("{s}, blocking {s}\n", .{ @typeName(T), if (vi == 0) "shrunk" else "own" });
                return e;
            };
        }
    }
}

test "potri: A X = I on potrf's factor, both triangles, on every linked core" {
    try check.eachCore(inverseCases);
}

test "potri: a zero on the factor's diagonal before A is touched; bad arguments" {
    const fa = testing.failing_allocator;
    inline for (.{ f64, types.Complex64 }) |T| for ([_]Uplo{ .u, .l }) |uplo| {
        var a: [9]T = @splat(scalar(T, 1));
        a[8] = scalar(T, 0);
        const before = a;
        try testing.expectError(error.Singular, potri(fa, check.io, T, uplo, 3, &a, 3));
        try testing.expectEqualSlices(T, &before, &a);
    };
    var a: [4]f64 = @splat(nan);
    try testing.expectError(error.UnsupportedUplo, potri(fa, check.io, f64, .m, 2, &a, 2));
    try testing.expectError(error.LeadingDimensionTooSmall, potri(fa, check.io, f64, .l, 2, &a, 1));
    try testing.expectError(error.SliceTooShort, potri(fa, check.io, f64, .l, 3, &a, 3));
    try potri(fa, check.io, f64, .u, 0, &a, 1);
    for (a) |v| try testing.expect(std.math.isNan(v));
}
