//! Tests of `../cholesky.zig`.

const std = @import("std");
const sysinfo = @import("builtin");
const blas = @import("../../blas/__root__.zig");
const types = @import("../../blas/types.zig");
const driver = @import("../../blas/blas3/driver.zig");
const accel = @import("../../blas/apple_accel.zig");
const Vtable = @import("../../blas/Vtable.zig");
const ops = @import("../ops.zig");
const interface = @import("../interface.zig");
const cholesky = @import("../cholesky.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Uplo = types.Uplo;
const Float = types.Float;
const isComplex = types.isComplex;
const scalar = ops.scalar;
const re = ops.re;
const conj = ops.conj;
const factorParallel = cholesky.factorParallel;
const potf2 = cholesky.potf2;
const potrf = cholesky.potrf;
const potrfChecked = cholesky.potrfChecked;
const potri = cholesky.potri;
const potrs = cholesky.potrs;

const testing = std.testing;
const check = @import("../../blas/_test/check.zig");
const test_io = @import("../../blas/_test/io.zig");
const guard = @import("../../blas/_test/guard.zig");
const reference = @import("../../blas/_test/reference.zig");
const lcheck = @import("check.zig");
const fixtures = @import("fixtures.zig");

const test_types = lcheck.test_types;
const W = lcheck.W;
const fromF64 = lcheck.fromF64;
const poisoned = lcheck.poisoned;
const isNanT = lcheck.isNanT;
const expectSameBits = lcheck.expectSameBits;

const nan = std.math.nan(f64);

/// Sizes around every block edge: on the shrunk blocking, and on a core's
/// own (`dtb_entries / 2` is 16 to 128 by core, so 130 is blocked on every
/// core, twice over where it is 32).
const shrunk_sizes = [_]usize{ 1, 2, 3, 4, 5, 7, 12, 13, 14, 17, 26 };
const own_sizes = [_]usize{ 1, 31, 33, 70, 130 };

fn factorCases(vt: *const Vtable) anyerror!void {
    return lcheck.eachBlocking(vt, lcheck.shrink, test_types, factorOn);
}

fn factorOn(comptime T: type, v: *const Vtable, vi: usize) anyerror!void {
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
fn expectFailed(comptime T: type, p: *fixtures.CholProblem(T), want: usize, what: []const u8) !void {
    var poison: guard.Poison = .{ .child = testing.allocator };
    const got = try p.factorIt(poison.allocator());
    if (got != want) {
        check.say("{s} {s} {s} n={d}: failed at pivot {?d}, want {d}\n", .{ what, @typeName(T), if (p.upper) "u" else "l", p.n, got, want });
        return error.TestUnexpectedResult;
    }
}

/// Positions of a bad element: the first, last and either side of block
/// edges on both blockings.
const bad_at = [_]usize{ 0, 1, 2, 3, 6, 12, 13, 25, 32, 33, 69 };

fn failCases(vt: *const Vtable) anyerror!void {
    return lcheck.eachBlocking(vt, lcheck.shrink, test_types, failOn);
}

fn failOn(comptime T: type, v: *const Vtable, vi: usize) anyerror!void {
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

/// `potrs` on `potrf`'s factor of `p.a`, `nrhs` right-hand sides with
/// `ldb = n + 1` (the padding row NaN), then the residual `A X - B` in `W(T)`
/// within `8 (n + 2)` ulps of `|L| |L|ᴴ |X|`, and the padding untouched. With
/// `nan_at`, `B(i, j)` is NaN: column `j` must come out all NaN, the others
/// as without it, bit for bit.
fn solveCase(comptime T: type, al: Allocator, p: *fixtures.CholProblem(T), nrhs: usize, seed: u64, nan_at: ?[2]usize) !void {
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

fn solveCases(vt: *const Vtable) anyerror!void {
    return lcheck.eachBlocking(vt, lcheck.shrink, test_types, solveOn);
}

fn solveOn(comptime T: type, v: *const Vtable, vi: usize) anyerror!void {
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

fn parallelCases(vt: *const Vtable) anyerror!void {
    return lcheck.eachBlocking(vt, lcheck.shrink, test_types, parallelOn);
}

fn parallelOn(comptime T: type, v: *const Vtable, vi: usize) anyerror!void {
    const al = testing.allocator;
    for ([_]bool{ true, false }) |upper| {
        const blk = driver.Blocking.of(T, v);
        const sizes: []const usize = if (vi == 0) &shrunk_sizes else &own_sizes;
        for (sizes) |n| {
            if (check.rig() and n > 70) continue;
            var p = try fixtures.CholProblem(T).init(al, n, upper, n * 29 + @as(u64, @intFromBool(upper)));
            defer p.deinit(al);
            var poison: guard.Poison = .{ .child = al };
            const k = try factorParallel(poison.allocator(), check.io, T, interface.table(T), upper, blk, n, p.a.ptr, p.lda, 1, 1, 1);
            try testing.expectEqual(@as(?usize, null), k);
            p.expectFactor() catch |e| {
                check.say("{s}, blocking {s}, parallel algorithm\n", .{ @typeName(T), if (vi == 0) "shrunk" else "own" });
                return e;
            };
            if (n < 13) continue;
            // a failure past the first block
            p.store();
            p.a[(n - 2) * (p.lda + 1)] = fromF64(T, -100);
            try testing.expectEqual(@as(?usize, n - 2), try factorParallel(al, check.io, T, interface.table(T), upper, blk, n, p.a.ptr, p.lda, 1, 1, 1));
        }
    }
}

test "potrf's parallel algorithm: A = L Lᴴ, and the failing minor, on every linked core" {
    try check.eachCore(parallelCases);
}

/// Inputs of a threaded-bits case: `A`, copied into each run's `out`.
fn Threaded(comptime T: type) type {
    return struct {
        p: *const fixtures.CholProblem(T),

        pub fn run(c: @This(), io: Io, want: u32, out: []T) anyerror!void {
            const vt = interface.table(T);
            @memcpy(out, c.p.a);
            _ = try factorParallel(testing.allocator, io, T, vt, c.p.upper, driver.Blocking.of(T, vt), c.p.n, out.ptr, c.p.lda, want, 1, 1);
        }

        /// The public `potrf`, its workers `blas.config`'s and its teams by
        /// `interface.team_grain`.
        pub fn grained(c: @This(), io: Io, want: u32, out: []T) anyerror!void {
            const saved = blas.config.cpu_count;
            defer blas.config.cpu_count = saved;
            blas.config.cpu_count = want;
            @memcpy(out, c.p.a);
            try potrf(testing.allocator, io, T, if (c.p.upper) .u else .l, c.p.n, out, c.p.lda);
        }

        /// The public `potri`, threading by `blas.config`'s workers.
        pub fn inverse(c: @This(), io: Io, want: u32, out: []T) anyerror!void {
            const saved = blas.config.cpu_count;
            defer blas.config.cpu_count = saved;
            blas.config.cpu_count = want;
            @memcpy(out, c.p.a);
            try potri(testing.allocator, io, T, if (c.p.upper) .u else .l, c.p.n, out, c.p.lda);
        }
    };
}

test "potrf's parallel algorithm and potri give one worker's bits at any worker count" {
    _ = blas.coreName();
    const was = accel.routing();
    _ = accel.route(false);
    defer _ = accel.route(was);
    const al = testing.allocator;
    inline for (.{ f32, f64, types.Complex32, types.Complex64 }) |T| for ([_]bool{ true, false }) |upper| {
        // past every type's `potri` threshold
        const n: usize = if (check.rig()) 70 else 210;
        var p = try fixtures.CholProblem(T).init(al, n, upper, 41);
        defer p.deinit(al);
        const c: Threaded(T) = .{ .p = &p };
        try expectSameBits(T, c, Threaded(T).run, p.a.len, "potrf");
        _ = try p.factorIt(al);
        try expectSameBits(T, c, Threaded(T).inverse, p.a.len, "potri");
    };
}

test "potrf gives one worker's bits at any worker count where its grain allows teams" {
    _ = blas.coreName();
    const was = accel.routing();
    _ = accel.route(false);
    defer _ = accel.route(was);
    const al = testing.allocator;
    inline for (.{ f32, f64, types.Complex32, types.Complex64 }) |T| for ([_]bool{ true, false }) |upper| {
        // the first `trsm` and `herk` take teams for every type, on Windows' grain too
        const n: usize = if (check.rig()) 70 else 420;
        var p = try fixtures.CholProblem(T).init(al, n, upper, 43);
        defer p.deinit(al);
        const c: Threaded(T) = .{ .p = &p };
        try expectSameBits(T, c, Threaded(T).grained, p.a.len, "potrf, grained");
    };
}

test "crew: one worker per grain, a complex product four, never more than asked" {
    const grain = interface.block_grain;
    try testing.expectEqual(@as(u32, 1), interface.crew(f64, 12, grain, 64, 64, 64));
    try testing.expectEqual(@as(u32, 1), interface.crew(f64, 12, grain, 100, 100, 199));
    try testing.expectEqual(@as(u32, 2), interface.crew(f64, 12, grain, 100, 100, 200));
    try testing.expectEqual(@as(u32, 4), interface.crew(types.Complex64, 12, grain, 100, 100, 100));
    try testing.expectEqual(@as(u32, 12), interface.crew(f32, 12, grain, 512, 512, 512));
    try testing.expectEqual(@as(u32, 1), interface.crew(f32, 1, grain, 512, 512, 512));
    try testing.expectEqual(@as(u32, 3), interface.crew(f32, 3, 0, 2, 2, 2));
    try testing.expectEqual(@as(u32, 5), interface.crew(f64, 5, grain, std.math.maxInt(u64), 2, 2));
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

/// `potri` on `potrf`'s factor of `p`: `A X = I` within `32 (n + 2)` ulps
/// of `|A| |X|`, `X` Hermitian from its stored triangle, nothing outside
/// it written.
fn inverseCase(comptime T: type, al: Allocator, p: *fixtures.CholProblem(T)) !void {
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

fn inverseCases(vt: *const Vtable) anyerror!void {
    return lcheck.eachBlocking(vt, lcheck.shrink, test_types, inverseOn);
}

fn inverseOn(comptime T: type, v: *const Vtable, vi: usize) anyerror!void {
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

test "potrf on a Team gives one worker's bits, every type, as the M6 bench found it did not" {
    if (comptime !test_io.hosted) return error.SkipZigTest;
    _ = blas.coreName();
    const was = accel.routing();
    _ = accel.route(false);
    defer _ = accel.route(was);
    const al = testing.allocator;
    var th: test_io.TestIo = .init(al, .{});
    defer th.deinit();
    var team: blas.Team = undefined;
    team.init(th.io(), 4);
    defer team.deinit();
    inline for (.{ f32, f64, types.Complex32, types.Complex64 }) |T| for ([_]bool{ true, false }) |upper| {
        const n: usize = 420;
        var p = try fixtures.CholProblem(T).init(al, n, upper, 47);
        defer p.deinit(al);
        const c: Threaded(T) = .{ .p = &p };
        const want = try al.alloc(T, p.a.len);
        defer al.free(want);
        const got = try al.alloc(T, p.a.len);
        defer al.free(got);
        try c.grained(Io.failing, 1, want);
        try c.grained(team.io(), 4, got);
        for (want, got, 0..) |x, y, i| if (!lcheck.sameBits(T, x, y)) {
            check.say("potrf {s} upper {}: a team differs from one worker at element {d}: {any} against {any}\n", .{ @typeName(T), upper, i, y, x });
            return error.TestUnexpectedResult;
        };
    };
}
