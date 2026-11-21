//! Tests of `../tgsyl.zig`.

const std = @import("std");
const types = @import("../../blas/types.zig");
const ops = @import("../ops.zig");
const sylvester = @import("../sylvester.zig");
const tgsyl_solve = @import("../tgsyl.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const isComplex = types.isComplex;
const Trans = types.Trans;
const scalar = ops.scalar;
const Dif = sylvester.Dif;
const tgsy2 = sylvester.tgsy2;
const tgsyl = tgsyl_solve.tgsyl;

const testing = std.testing;
const check = @import("../../blas/_test/check.zig");
const reference = @import("../../blas/_test/reference.zig");
const lcheck = @import("check.zig");
const fixtures = @import("fixtures.zig");
const faer_svd = @import("../svd.zig");

const schur = fixtures.schur;

/// An `n x n` upper triangular matrix like `schur`'s, without its 2 x 2
/// blocks: `D` and `E` of a generalized Schur form.
fn upperTri(comptime T: type, al: Allocator, n: usize, lo: f64, seed: u64) ![]T {
    const a = try schur(T, al, n, lo, n + 1, seed);
    if (comptime !isComplex(T)) for (0..n -| 1) |i| {
        a[i + 1 + i * n] = 0;
    };
    return a;
}

/// `tgsyl`'s operands, each at its own padded leading dimension so a mixed
/// up one shows: `(A, D)` with eigenvalues in `[0.8, 8]`, `(B, E)` in
/// `[-8, -0.8]`.
fn Gen(comptime T: type) type {
    return struct {
        m: usize,
        n: usize,
        a: []T,
        b: []T,
        d: []T,
        e: []T,
        c0: []T,
        f0: []T,
        c: []T,
        f: []T,
        const ld = .{ .a = 1, .b = 2, .c = 3, .d = 2, .e = 1, .f = 5 };

        const Self = @This();

        fn init(al: Allocator, m: usize, n: usize, seed: u64) !Self {
            var g: Self = undefined;
            g.m = m;
            g.n = n;
            g.a = try padded(al, m, ld.a, try schur(T, al, m, 2, 3, seed));
            g.d = try padded(al, m, ld.d, try upperTri(T, al, m, 0.5, seed + 1));
            g.b = try padded(al, n, ld.b, try schur(T, al, n, -4, 4, seed + 2));
            g.e = try padded(al, n, ld.e, try upperTri(T, al, n, 0.5, seed + 3));
            g.c0 = try al.alloc(T, (m + ld.c) * n);
            check.fill(T, g.c0, seed + 4);
            g.f0 = try al.alloc(T, (m + ld.f) * n);
            check.fill(T, g.f0, seed + 5);
            g.c = try al.dupe(T, g.c0);
            g.f = try al.dupe(T, g.f0);
            return g;
        }

        /// `x` (`k x k`, packed) at leading dimension `k + pad`, freed.
        fn padded(al: Allocator, k: usize, pad: usize, x: []T) ![]T {
            defer al.free(x);
            const y = try al.alloc(T, (k + pad) * k);
            @memset(y, lcheck.poisoned(T));
            for (0..k) |j| @memcpy(y[j * (k + pad) ..][0..k], x[j * k ..][0..k]);
            return y;
        }

        fn deinit(g: Self, al: Allocator) void {
            for ([_][]T{ g.a, g.b, g.d, g.e, g.c0, g.f0, g.c, g.f }) |x| al.free(x);
        }

        fn reset(g: Self) void {
            @memcpy(g.c, g.c0);
            @memcpy(g.f, g.f0);
        }

        /// `tgsyl`, unperturbed, its `scale` and any `dif` written out.
        fn run(g: Self, al: Allocator, io: Io, trans: Trans, ijob: usize, scale: *Float(T), dif: *Float(T)) !void {
            const m = g.m;
            const n = g.n;
            const r = try tgsyl(al, io, T, trans, ijob, m, n, g.a, m + ld.a, g.b, n + ld.b, g.c, m + ld.c, g.d, m + ld.d, g.e, n + ld.e, g.f, m + ld.f);
            if (r.perturbed) return error.TestUnexpectedResult;
            scale.* = r.scale;
            if (r.separation) |x| dif.* = x;
        }

        fn at(x: []const T, ldx: usize, i: usize, j: usize) Wt {
            return reference.wide(T, x[i + j * ldx]);
        }
        const Wt = @TypeOf(reference.wide(T, scalar(T, 0)));

        /// Both equations' residuals within `32 (m + n + 2)` ulps of their
        /// terms' magnitudes; padding untouched.
        fn expectSolved(g: Self, tr: bool, scale: Float(T)) !void {
            const m = g.m;
            const n = g.n;
            const lda = m + ld.a;
            const ldb = n + ld.b;
            const ldc = m + ld.c;
            const ldd = m + ld.d;
            const lde = n + ld.e;
            const ldf = m + ld.f;
            const eps: f128 = std.math.floatEps(Float(T));
            const ulps: f128 = @floatFromInt(32 * (m + n + 2));
            const cjw = struct {
                fn f(x: Wt, c: bool) Wt {
                    return if (c) reference.conjW(x) else x;
                }
            }.f;
            for (0..n) |j| for (0..m) |i| for (0..2) |eq| {
                const rhs = if (eq == 0) at(g.c0, ldc, i, j) else at(g.f0, ldf, i, j);
                // notrans: A R - L B, D R - L E; trans: Aᴴ R + Dᴴ L, R Bᴴ + L Eᴴ
                var r = reference.scaleW(if (tr and eq == 1) @as(f128, scale) else -@as(f128, scale), rhs);
                var mag: f128 = scale * reference.abs1(rhs);
                for (0..@max(m, n)) |k| {
                    var terms: [2][2]Wt = undefined;
                    var used: [2]bool = .{ false, false };
                    if (!tr) {
                        if (k < m) {
                            terms[0] = .{ at(if (eq == 0) g.a else g.d, if (eq == 0) lda else ldd, i, k), at(g.c, ldc, k, j) };
                            used[0] = true;
                        }
                        if (k < n) {
                            terms[1] = .{ reference.scaleW(-1, at(g.f, ldf, i, k)), at(if (eq == 0) g.b else g.e, if (eq == 0) ldb else lde, k, j) };
                            used[1] = true;
                        }
                    } else if (eq == 0) {
                        if (k < m) {
                            terms[0] = .{ cjw(at(g.a, lda, k, i), true), at(g.c, ldc, k, j) };
                            terms[1] = .{ cjw(at(g.d, ldd, k, i), true), at(g.f, ldf, k, j) };
                            used = .{ true, true };
                        }
                    } else if (k < n) {
                        terms[0] = .{ at(g.c, ldc, i, k), cjw(at(g.b, ldb, j, k), true) };
                        terms[1] = .{ at(g.f, ldf, i, k), cjw(at(g.e, lde, j, k), true) };
                        used = .{ true, true };
                    }
                    for (terms, used) |t, u| if (u) {
                        r = reference.addW(r, reference.mulW(t[0], t[1]));
                        mag += reference.abs1(t[0]) * reference.abs1(t[1]);
                    };
                }
                const dv = reference.abs1(r);
                if (!(dv <= ulps * eps * mag + std.math.floatMin(Float(T)))) {
                    check.say("tgsyl {s} {d}x{d} trans {any}: equation {d} ({d},{d}) residual {e}, allowed {e}\n", .{ @typeName(T), m, n, tr, eq + 1, i, j, @as(f64, @floatCast(dv)), @as(f64, @floatCast(ulps * eps * mag)) });
                    return error.TestUnexpectedResult;
                }
            };
        }
    };
}

fn transOf(comptime T: type, tr: bool) Trans {
    return if (!tr) .n else if (comptime isComplex(T)) .c else .t;
}

const tgsyl_shapes = [_][2]usize{ .{ 1, 1 }, .{ 2, 3 }, .{ 5, 4 }, .{ 24, 24 }, .{ 25, 10 }, .{ 10, 40 }, .{ 60, 50 }, .{ 33, 70 } };

test "relapack tgsyl: A R - L B = scale C, D R - L E = scale F and the transposed system, every type" {
    const al = testing.allocator;
    inline for (lcheck.test_types) |T| for (tgsyl_shapes) |s| {
        if (check.rig() and s[0] * s[1] > 1000) continue;
        const g = try Gen(T).init(al, s[0], s[1], s[0] * 7 + s[1]);
        defer g.deinit(al);
        for ([_]bool{ false, true }) |tr| {
            g.reset();
            var scale: Float(T) = 0;
            var dif: Float(T) = -1;
            try g.run(al, check.io, transOf(T, tr), 0, &scale, &dif);
            try testing.expect(scale > 0 and scale <= 1 and dif == -1);
            g.expectSolved(tr, scale) catch |e| {
                check.say("tgsyl {s} {d}x{d}\n", .{ @typeName(T), s[0], s[1] });
                return e;
            };
        }
        // ijob 1 and 2 solve as 0 does, bit for bit, and estimate
        const x = try al.dupe(T, g.c);
        defer al.free(x);
        const y = try al.dupe(T, g.f);
        defer al.free(y);
        g.reset();
        var scale: Float(T) = 0;
        var dif: Float(T) = -1;
        try g.run(al, check.io, .n, 0, &scale, &dif);
        @memcpy(x, g.c);
        @memcpy(y, g.f);
        for ([_]usize{ 1, 2 }) |job| {
            g.reset();
            try g.run(al, check.io, .n, job, &scale, &dif);
            try testing.expect(dif > 0 and std.math.isFinite(dif));
            // by value: f80's padding bytes are not its bits
            for (x, g.c) |u, v| try testing.expect(lcheck.sameBits(T, u, v));
            for (y, g.f) |u, v| try testing.expect(lcheck.sameBits(T, u, v));
        }
    };
}

/// `σ_min(Z)` of `dtgsyl.f`'s (2), from `gesvd` in `f64` or `c64`.
fn sigmaMin(comptime T: type, al: Allocator, g: Gen(T)) !f64 {
    const C = if (comptime isComplex(T)) types.Complex64 else f64;
    const m = g.m;
    const n = g.n;
    const k = m * n;
    const z = try al.alloc(C, 4 * k * k);
    defer al.free(z);
    @memset(z, scalar(C, 0));
    const lda = m + Gen(T).ld.a;
    const ldb = n + Gen(T).ld.b;
    const ldd = m + Gen(T).ld.d;
    const lde = n + Gen(T).ld.e;
    const to = struct {
        fn f(v: T) C {
            return if (comptime isComplex(T)) .init(@floatCast(v.re), @floatCast(v.im)) else @floatCast(v);
        }
    }.f;
    for (0..n) |j| for (0..m) |i| {
        const r = i + m * j;
        for (0..m) |kk| {
            z[r + 2 * k * (kk + m * j)] = to(g.a[i + kk * lda]);
            z[k + r + 2 * k * (kk + m * j)] = to(g.d[i + kk * ldd]);
        }
        for (0..n) |kk| {
            z[r + 2 * k * (k + i + m * kk)] = ops.neg(C, to(g.b[kk + j * ldb]));
            z[k + r + 2 * k * (k + i + m * kk)] = ops.neg(C, to(g.e[kk + j * lde]));
        }
    };
    const s = try al.alloc(f64, 2 * k);
    defer al.free(s);
    try faer_svd.gesvd(al, check.io, C, .none, .none, 2 * k, 2 * k, z, 2 * k, s, z[0..0], 1, z[0..0], 1, null);
    return s[2 * k - 1];
}

test "relapack tgsyl: dif estimates σ_min of the Kronecker form; ijob 1 and 3, 2 and 4 agree" {
    const al = testing.allocator;
    inline for (lcheck.test_types) |T| for ([_][2]usize{ .{ 1, 1 }, .{ 2, 2 }, .{ 3, 5 }, .{ 6, 4 } }) |s| {
        const g = try Gen(T).init(al, s[0], s[1], s[0] * 11 + s[1]);
        defer g.deinit(al);
        const sig = try sigmaMin(T, al, g);
        var difs: [5]Float(T) = undefined;
        for (1..5) |job| {
            g.reset();
            var scale: Float(T) = 0;
            difs[job] = -1;
            try g.run(al, check.io, .n, @intCast(job), &scale, &difs[job]);
            // an upper bound on σ_min, within sqrt(2 m n) of it as the
            // estimate's right-hand sides are about that long
            const d: f64 = @floatCast(difs[job]);
            const lim = 4 * @sqrt(2 * @as(f64, @floatFromInt(s[0] * s[1])));
            if (!(d >= sig * (1 - 1e-3) / lim and d <= sig * lim)) {
                check.say("tgsyl {s} {d}x{d} ijob {d}: dif {e}, sigma_min {e}\n", .{ @typeName(T), s[0], s[1], job, d, sig });
                return error.TestUnexpectedResult;
            }
        }
        try testing.expectEqual(difs[1], difs[3]);
        try testing.expectEqual(difs[2], difs[4]);
    };
}

test "relapack tgsyl: the recursion's dif is the unblocked xTGSY2's, every base case's blocks counted" {
    // ReLAPACK keeps only the last base case's `PQ`: its `ijob` 2 and 4
    // estimate would come out about sqrt(base cases) too small
    const al = testing.allocator;
    inline for (.{ f64, types.Complex32 }) |T| {
        const F = Float(T);
        const g = try Gen(T).init(al, 40, 30, 3);
        defer g.deinit(al);
        for ([_]usize{ 3, 4 }) |job| {
            g.reset();
            var scale: F = 0;
            var dif: F = 0;
            try g.run(al, check.io, .n, job, &scale, &dif);
            // the whole problem as one `xTGSY2`
            g.reset();
            for (0..g.n) |j| for (0..g.m) |i| {
                g.c[i + j * (g.m + 3)] = scalar(T, 0);
                g.f[i + j * (g.m + 5)] = scalar(T, 0);
            };
            var acc: Dif(F) = .{};
            _ = try tgsy2(al, check.io, T, false, @intCast(job - 2), g.m, g.n, g.a.ptr, g.m + 1, g.b.ptr, g.n + 2, g.c.ptr, g.m + 3, g.d.ptr, g.m + 2, g.e.ptr, g.n + 1, g.f.ptr, g.m + 5, &acc);
            const num: F = @floatFromInt(if (job == 3) 2 * g.m * g.n else acc.pq);
            const want = @sqrt(num) / (acc.scale * @sqrt(acc.sum));
            try testing.expect(@abs(dif - want) <= 1e3 * std.math.floatEps(F) * want);
        }
    }
}

test "relapack tgsyl: common eigenvalues are perturbed, written and finite; scale below 1; bad arguments" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64 }) |T| {
        const F = Float(T);
        // (A, D) = (B, E): Z is singular; a huge right-hand side scales
        const a = [_]T{ scalar(T, 2), scalar(T, 0), scalar(T, 1), scalar(T, 2) };
        const d = [_]T{ scalar(T, 1), scalar(T, 0), scalar(T, 0.5), scalar(T, 1) };
        for ([_]F{ 1, std.math.floatMax(F) / 4 }) |big| {
            var c: [4]T = @splat(scalar(T, big));
            var f: [4]T = @splat(scalar(T, big));
            const r = try tgsyl(al, check.io, T, .n, 0, 2, 2, &a, 2, &a, 2, &c, 2, &d, 2, &d, 2, &f, 2);
            try testing.expect(r.perturbed and r.separation == null);
            try testing.expect(r.scale > 0 and r.scale <= 1);
            if (big > 1) try testing.expect(r.scale < 1);
            for (c ++ f) |x| try testing.expect(if (comptime isComplex(T)) std.math.isFinite(x.re) and std.math.isFinite(x.im) else std.math.isFinite(x));
        }
        var c: [4]T = @splat(lcheck.poisoned(T));
        const t: Trans = if (comptime isComplex(T)) .t else .c;
        const z = &c;
        try testing.expectError(error.UnsupportedTranspose, tgsyl(al, check.io, T, t, 0, 2, 2, z, 2, z, 2, z, 2, z, 2, z, 2, z, 2));
        try testing.expectError(error.UnsupportedJob, tgsyl(al, check.io, T, .n, 5, 2, 2, z, 2, z, 2, z, 2, z, 2, z, 2, z, 2));
        try testing.expectError(error.LeadingDimensionTooSmall, tgsyl(al, check.io, T, .n, 0, 2, 2, z, 1, z, 2, z, 2, z, 2, z, 2, z, 2));
        try testing.expectError(error.LeadingDimensionTooSmall, tgsyl(al, check.io, T, .n, 0, 2, 2, z, 2, z, 1, z, 2, z, 2, z, 2, z, 2));
        try testing.expectError(error.LeadingDimensionTooSmall, tgsyl(al, check.io, T, .n, 0, 2, 2, z, 2, z, 2, z, 1, z, 2, z, 2, z, 2));
        try testing.expectError(error.LeadingDimensionTooSmall, tgsyl(al, check.io, T, .n, 0, 2, 2, z, 2, z, 2, z, 2, z, 1, z, 2, z, 2));
        try testing.expectError(error.LeadingDimensionTooSmall, tgsyl(al, check.io, T, .n, 0, 2, 2, z, 2, z, 2, z, 2, z, 2, z, 1, z, 2));
        try testing.expectError(error.LeadingDimensionTooSmall, tgsyl(al, check.io, T, .n, 0, 2, 2, z, 2, z, 2, z, 2, z, 2, z, 2, z, 1));
        try testing.expectError(error.SliceTooShort, tgsyl(al, check.io, T, .n, 0, 2, 2, z, 2, z, 2, z, 2, z, 2, z[0..3], 2, z, 2));
        // empty, which the reference refuses: nothing read or written
        for ([_][2]usize{ .{ 0, 2 }, .{ 2, 0 } }) |s| {
            const r = try tgsyl(al, check.io, T, transOf(T, true), 9, s[0], s[1], z, 2, z, 2, z, 2, z, 2, z, 2, z, 2);
            try testing.expect(r.scale == 1 and r.separation == null and !r.perturbed);
        }
        for (c) |x| try testing.expect(lcheck.isNanT(T, x));
    }
}

const BigGsyl = struct {
    g: Gen(f64),

    fn run(self: BigGsyl, io: Io, out: []f64) anyerror!void {
        self.g.reset();
        var scale: f64 = 0;
        var dif: f64 = 0;
        try self.g.run(testing.allocator, io, .n, 1, &scale, &dif);
        @memcpy(out[0..self.g.c.len], self.g.c);
        @memcpy(out[self.g.c.len..][0..self.g.f.len], self.g.f);
        out[out.len - 1] = dif;
    }
};

test "relapack tgsyl: the same bits inline and threaded" {
    const al = testing.allocator;
    const g = try Gen(f64).init(al, lcheck.big_n, lcheck.big_n, 9);
    defer g.deinit(al);
    try lcheck.expectIoBits(f64, BigGsyl{ .g = g }, BigGsyl.run, g.c.len + g.f.len + 1, "tgsyl");
}
