//! Tests of ReLAPACK's `pbtrf`.

const std = @import("std");
const types = @import("../../../blas/types.zig");
const check = @import("../../../blas/_test/check.zig");
const reference = @import("../../../blas/_test/reference.zig");
const lcheck = @import("../check.zig");
const ops = @import("../../ops.zig");
const relapack = @import("../../relapack.zig");

const testing = std.testing;
const Allocator = std.mem.Allocator;
const Float = types.Float;
const scalar = ops.scalar;
const test_types = lcheck.test_types;
const isNanT = lcheck.isNanT;
const W = lcheck.W;

// ---- pbtrf tests ------------------------------------------------------------------

/// A Hermitian (real: symmetric) positive-definite band matrix `A = L Lᴴ`,
/// `L` lower with `kd` subdiagonals, its diagonal in [2, 4) and the rest
/// within `1 / (kd + 1)`; stored in `upper`'s band form at `ldab = kd + 2`,
/// every unused place NaN.
pub fn Band(comptime T: type) type {
    return struct {
        n: usize,
        kd: usize,
        ldab: usize,
        upper: bool,
        ab: []T,
        aw: []W(T),

        const Self = @This();

        pub fn init(al: Allocator, n: usize, kd: usize, upper: bool, seed: u64) !Self {
            const v = try al.alloc(T, n * n);
            defer al.free(v);
            check.fill(T, v, seed);
            const scale: Float(T) = @floatCast(1 / @as(f64, @floatFromInt(kd + 1)));
            const l = try al.alloc(W(T), n * n);
            defer al.free(l);
            for (0..n) |j| for (0..n) |i| {
                const x = v[i + j * n];
                l[i + j * n] = reference.wide(T, if (i < j or i > j + kd) scalar(T, 0) else if (i == j) scalar(T, 3 + ops.re(T, x) / 2) else ops.mulReal(T, x, scale));
            };
            const aw = try al.alloc(W(T), n * n);
            errdefer al.free(aw);
            for (0..n) |j| for (0..n) |i| {
                var s = reference.wide(T, scalar(T, 0));
                for (0..@min(i, j) + 1) |k| s = reference.addW(s, reference.mulW(l[i + k * n], reference.conjW(l[j + k * n])));
                aw[i + j * n] = s;
            };
            const ldab = kd + 2;
            const ab = try al.alloc(T, @max(1, ldab * n));
            var self: Self = .{ .n = n, .kd = kd, .ldab = ldab, .upper = upper, .ab = ab, .aw = aw };
            self.store();
            return self;
        }

        pub fn deinit(self: *Self, al: Allocator) void {
            al.free(self.ab);
            al.free(self.aw);
        }

        /// The storage place of `(i, j)`, `i >= j` (`l`) or `i <= j` (`u`),
        /// inside the band, else null.
        pub fn place(self: *const Self, i: usize, j: usize) ?usize {
            if (self.upper) {
                if (i > j or j - i > self.kd) return null;
                return self.kd + i - j + j * self.ldab;
            }
            if (i < j or i - j > self.kd) return null;
            return i - j + j * self.ldab;
        }

        pub fn store(self: *Self) void {
            @memset(self.ab, lcheck.poisoned(T));
            for (0..self.n) |j| for (0..self.n) |i| if (self.place(i, j)) |at| {
                self.ab[at] = lcheck.toT(T, self.aw[i + j * self.n]);
            };
        }

        /// The factor's `(i, j)` as lower, `i >= j`, zero outside the band.
        pub fn got(self: *const Self, i: usize, j: usize) W(T) {
            const at = (if (self.upper) self.place(j, i) else self.place(i, j)) orelse return reference.wide(T, scalar(T, 0));
            const v = reference.wide(T, self.ab[at]);
            return if (self.upper) reference.conjW(v) else v;
        }

        /// `F Fᴴ` against `A` within `2 (n + 2)` ulps of `|F| |F|ᴴ`, the
        /// unused places untouched.
        pub fn expectFactor(self: *const Self) !void {
            const n = self.n;
            const eps: f128 = std.math.floatEps(Float(T));
            const ulps: f128 = @floatFromInt(2 * (n + 2));
            for (0..n) |j| for (j..n) |i| {
                var s = reference.wide(T, scalar(T, 0));
                var mag: f128 = 0;
                const lo = @min(if (i > self.kd) i - self.kd else 0, j + 1);
                for (lo..j + 1) |k| {
                    const x = self.got(i, k);
                    const y = reference.conjW(self.got(j, k));
                    s = reference.addW(s, reference.mulW(x, y));
                    mag += reference.abs1(x) * reference.abs1(y);
                }
                const d = reference.abs1(reference.subW(s, self.aw[i + j * n]));
                if (!(d <= ulps * eps * mag + std.math.floatMin(Float(T)))) {
                    check.say("pbtrf {s} {s} n={d} kd={d}: (F Fᴴ)({d},{d}) off by {e}\n", .{ @typeName(T), if (self.upper) "u" else "l", n, self.kd, i, j, @as(f64, @floatCast(d)) });
                    return error.TestUnexpectedResult;
                }
            };
            for (0..n) |j| for (0..self.ldab) |r| {
                const used = r <= self.kd and (if (self.upper) r + j >= self.kd else r + j < n);
                if (!used and !isNanT(T, self.ab[r + j * self.ldab])) {
                    check.say("pbtrf {s} n={d} kd={d}: wrote unused ({d},{d})\n", .{ @typeName(T), n, self.kd, r, j });
                    return error.TestUnexpectedResult;
                }
            };
        }
    };
}

pub const band_cases = [_][2]usize{ .{ 1, 0 }, .{ 5, 0 }, .{ 30, 1 }, .{ 30, 5 }, .{ 26, 25 }, .{ 50, 10 }, .{ 60, 30 }, .{ 40, 60 } };

test "relapack pbtrf: A = F Fᴴ in band storage, both triangles, every type" {
    const al = testing.allocator;
    inline for (test_types) |T| for ([_]bool{ true, false }) |upper| for (band_cases) |c| {
        var p = try Band(T).init(al, c[0], c[1], upper, c[0] * 3 + c[1]);
        defer p.deinit(al);
        relapack.pbtrf(al, check.io, T, lcheck.uploOf(upper), c[0], c[1], p.ab, p.ldab) catch |e| {
            check.say("pbtrf {s} n={d} kd={d}\n", .{ @typeName(T), c[0], c[1] });
            return e;
        };
        try p.expectFactor();
    };
}

test "relapack pbtrf: a rest wider than the band past the block (ReLAPACK's dense case), f64 and c64" {
    if (check.rig()) return error.SkipZigTest;
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64 }) |T| for ([_]bool{ true, false }) |upper| for ([_][2]usize{ .{ 100, 50 }, .{ 100, 60 }, .{ 120, 30 } }) |c| {
        var p = try Band(T).init(al, c[0], c[1], upper, 11);
        defer p.deinit(al);
        try relapack.pbtrf(al, check.io, T, lcheck.uploOf(upper), c[0], c[1], p.ab, p.ldab);
        try p.expectFactor();
    };
}

test "relapack pbtrf: the first pivot that is not positive, NaN too, is where it stops" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex32 }) |T| for ([_]bool{ true, false }) |upper| {
        var p = try Band(T).init(al, 60, 10, upper, 2);
        defer p.deinit(al);
        for ([_]usize{ 0, 23, 30, 59 }) |k| for ([_]f64{ -100, std.math.nan(f64) }) |x| {
            p.store();
            p.ab[p.place(k, k).?] = lcheck.fromF64(T, x);
            try testing.expectError(error.NotPositiveDefinite, relapack.pbtrf(al, check.io, T, lcheck.uploOf(upper), 60, 10, p.ab, p.ldab));
            // the factor's diagonal positive before `k`, the failed pivot at it
            for (0..k + 1) |i| {
                const d = ops.re(T, p.ab[p.place(i, i).?]);
                try testing.expectEqual(i < k, d > 0);
            }
        };
    };
}

test "relapack pbtrf: a reversed uplo, a short slice, a small leading dimension; nothing read for n = 0" {
    const al = testing.allocator;
    var ab: [6]f64 = @splat(lcheck.poisoned(f64));
    try testing.expectError(error.UnsupportedUplo, relapack.pbtrf(al, check.io, f64, .v, 2, 1, &ab, 2));
    try testing.expectError(error.LeadingDimensionTooSmall, relapack.pbtrf(al, check.io, f64, .l, 2, 1, &ab, 1));
    try testing.expectError(error.SliceTooShort, relapack.pbtrf(al, check.io, f64, .l, 4, 1, &ab, 2));
    try relapack.pbtrf(al, check.io, f64, .l, 0, 1, &ab, 2);
    for (ab) |x| try testing.expect(std.math.isNan(x));
}
