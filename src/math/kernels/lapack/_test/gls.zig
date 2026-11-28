//! Tests: `gglse` and `ggglm`, `dlsets.f`'s and `dglmts.f`'s residuals and
//! scalings, and the optimality conditions against `geqrf`'s null spaces.

const std = @import("std");
const types = @import("../../blas/types.zig");
const check = @import("../../blas/_test/check.zig");
const lcheck = @import("check.zig");
const ops = @import("../ops.zig");
const qr = @import("../qr.zig");
const gls = @import("../gls.zig");

const testing = std.testing;
const Allocator = std.mem.Allocator;
const Mat = ops.Mat;
const Float = types.Float;

fn asum(comptime T: type, v: []const T) Float(T) {
    var s: Float(T) = 0;
    for (v) |e| s += ops.abs(T, e);
    return s;
}

/// `‖A‖₁`, the largest column sum.
fn norm1(comptime T: type, a: Mat(T)) Float(T) {
    var big: Float(T) = 0;
    for (0..a.cols) |j| big = @max(big, asum(T, a.col(j)[0..a.rows]));
    return big;
}

fn vmax(comptime T: type, v: []const T) Float(T) {
    var big: Float(T) = 0;
    for (v) |e| big = @max(big, ops.abs(T, e));
    return big;
}

fn maxAbs(comptime T: type, a: Mat(T)) Float(T) {
    var big: Float(T) = 0;
    for (0..a.cols) |j| for (a.col(j)[0..a.rows]) |e| {
        big = @max(big, ops.abs(T, e));
    };
    return big;
}

/// `y := y - A x`.
fn subMul(comptime T: type, a: Mat(T), x: []const T, y: []T) void {
    for (0..a.cols) |j| for (0..a.rows) |i| {
        y[i] = ops.sub(T, y[i], ops.mul(T, a.get(i, j), x[j]));
    };
}

/// The `r x r` unitary `Q` of `w`'s QR (`r x c`, `r >= c`): its last
/// `r - c` columns span `w`'s complement.
fn fullQ(comptime T: type, al: Allocator, w: Mat(T)) !Mat(T) {
    const r = w.rows;
    const k = w.cols;
    const q = try ops.alloc(T, al, r, r);
    errdefer ops.free(T, al, q);
    q.fill(ops.zero(T));
    q.sub(0, 0, r, k).copyFrom(w);
    const tau = try al.alloc(T, @max(k, 1));
    defer al.free(tau);
    try qr.geqrf(al, check.io, T, r, k, q.slice(), q.ld, tau[0..k]);
    try qr.orgqr(al, check.io, T, r, r, k, q.slice(), q.ld, tau[0..k]);
    return q;
}

/// `‖Q(:, from..)ᴴ v‖∞`.
fn outside(comptime T: type, q: Mat(T), from: usize, v: []const T) Float(T) {
    var big: Float(T) = 0;
    for (from..q.cols) |j| {
        var s = ops.zero(T);
        for (0..q.rows) |i| s = ops.add(T, s, ops.mul(T, ops.conj(T, q.get(i, j)), v[i]));
        big = @max(big, ops.abs(T, s));
    }
    return big;
}

/// `dlsets.f`'s and `dglmts.f`'s scalings: common large, common tiny,
/// tiny right-hand sides only, tiny `A` side, tiny `B` side.
fn scales(comptime F: type, tnrm: F) [5]struct { name: []const u8, a: F, b: F, c: F, d: F } {
    const e = std.math.frexp(tnrm).exponent;
    const maxexp = std.math.frexp(std.math.floatMax(F)).exponent - 2;
    const tiny = std.math.ldexp(@as(F, 1), std.math.frexp(std.math.floatMin(F) / std.math.floatEps(F)).exponent - 4 - e);
    const big = std.math.ldexp(@as(F, 1), maxexp - e);
    return .{
        .{ .name = "common large", .a = big, .b = big, .c = big, .d = big },
        .{ .name = "common tiny", .a = tiny, .b = tiny, .c = tiny, .d = tiny },
        .{ .name = "tiny c, d", .a = 1, .b = 1, .c = tiny, .d = tiny },
        .{ .name = "tiny A, c", .a = tiny, .b = 1, .c = tiny, .d = 1 },
        .{ .name = "tiny B, d", .a = 1, .b = tiny, .c = 1, .d = tiny },
    };
}

fn scaleBy(comptime T: type, v: []T, s: Float(T)) void {
    for (v) |*e| e.* = ops.mulReal(T, e.*, s);
}

/// `gglse` on copies of `A`, `B`, `c`, `d` scaled by `sc`; `x` out, `c`'s
/// residual part into `cres` when given.
fn runLse(comptime T: type, al: Allocator, a0: Mat(T), b0: Mat(T), c0: []const T, d0: []const T, sc: anytype, x: []T, cres: ?[]T) !void {
    const m = a0.rows;
    const n = a0.cols;
    const p = b0.rows;
    const a = try ops.alloc(T, al, m, n);
    defer ops.free(T, al, a);
    a.copyFrom(a0);
    scaleBy(T, a.slice(), sc.a);
    const b = try ops.alloc(T, al, p, n);
    defer ops.free(T, al, b);
    b.copyFrom(b0);
    scaleBy(T, b.slice(), sc.b);
    const c = try al.dupe(T, c0);
    defer al.free(c);
    scaleBy(T, c, sc.c);
    const d = try al.dupe(T, d0);
    defer al.free(d);
    scaleBy(T, d, sc.d);
    try gls.gglse(al, check.io, T, m, n, p, a.slice(), a.ld, b.slice(), b.ld, c, d, x);
    if (cres) |r| @memcpy(r, c[n - p ..]);
}

/// `x` at the input scale: `B x = d` (`dget02.f`'s ratio) and `Aᴴ (c - A
/// x)` outside `null(B)` only, `null(B)` from `geqrf` of `Bᴴ`.
fn expectLse(comptime T: type, al: Allocator, a0: Mat(T), b0: Mat(T), bq: Mat(T), c0: []const T, d0: []const T, x: []const T, what: []const u8) !void {
    const F = Float(T);
    const eps = std.math.floatEps(F);
    const m = a0.rows;
    const n = a0.cols;
    const p = b0.rows;
    const xn = asum(T, x);
    if (p > 0) {
        const rb = try al.dupe(T, d0);
        defer al.free(rb);
        subMul(T, b0, x, rb);
        const ratio = asum(T, rb) / norm1(T, b0) / @max(xn, std.math.floatMin(F)) / eps;
        if (!(ratio <= lcheck.thresh)) {
            check.say("gglse {s} {s} {d} x {d}, p {d}: B x - d ratio {e}\n", .{ what, @typeName(T), m, n, p, @as(f64, @floatCast(ratio)) });
            return error.TestUnexpectedResult;
        }
    }
    if (n == p) return;
    // g = Aᴴ (c - A x)
    const r = try al.dupe(T, c0);
    defer al.free(r);
    subMul(T, a0, x, r);
    const g = try al.alloc(T, n);
    defer al.free(g);
    for (g, 0..) |*gj, j| {
        var s = ops.zero(T);
        for (0..m) |i| s = ops.add(T, s, ops.mul(T, ops.conj(T, a0.get(i, j)), r[i]));
        gj.* = s;
    }
    const an = norm1(T, a0);
    const tol = lcheck.thresh * @as(F, @floatFromInt(m + n)) * eps * an * (an * xn + asum(T, c0));
    const off = outside(T, bq, p, g);
    if (!(off <= tol)) {
        check.say("gglse {s} {s} {d} x {d}, p {d}: Aᴴ r in null(B) {e}, tolerance {e}\n", .{ what, @typeName(T), m, n, p, @as(f64, @floatCast(off)), @as(f64, @floatCast(tol)) });
        return error.TestUnexpectedResult;
    }
}

fn lseCase(comptime T: type, al: Allocator, m: usize, n: usize, p: usize, seed: u64) !void {
    const F = Float(T);
    const eps = std.math.floatEps(F);
    const a0 = try ops.alloc(T, al, m, n);
    defer ops.free(T, al, a0);
    check.fill(T, a0.slice(), seed);
    const b0 = try ops.alloc(T, al, p, n);
    defer ops.free(T, al, b0);
    check.fill(T, b0.slice(), seed + 1);
    const c0 = try al.alloc(T, m);
    defer al.free(c0);
    check.fill(T, c0, seed + 2);
    const d0 = try al.alloc(T, p);
    defer al.free(d0);
    check.fill(T, d0, seed + 3);
    // Bᴴ's full Q
    const bh = try ops.alloc(T, al, n, p);
    defer ops.free(T, al, bh);
    for (0..p) |j| for (0..n) |i| {
        bh.at(i, j).* = ops.conj(T, b0.get(j, i));
    };
    const bq = try fullQ(T, al, bh);
    defer ops.free(T, al, bq);
    const x = try al.alloc(T, n);
    defer al.free(x);
    const cres = try al.alloc(T, m + p - n);
    defer al.free(cres);
    try runLse(T, al, a0, b0, c0, d0, .{ .a = @as(F, 1), .b = @as(F, 1), .c = @as(F, 1), .d = @as(F, 1) }, x, cres);
    try expectLse(T, al, a0, b0, bq, c0, d0, x, "unscaled");
    // ‖c[n - p..m]‖ = ‖c - A x‖
    const r = try al.dupe(T, c0);
    defer al.free(r);
    subMul(T, a0, x, r);
    var want: F = 0;
    for (r) |e| want += ops.abs(T, e) * ops.abs(T, e);
    var got: F = 0;
    for (cres) |e| got += ops.abs(T, e) * ops.abs(T, e);
    const tol = lcheck.thresh * @as(F, @floatFromInt(m + n)) * eps * (norm1(T, a0) * asum(T, x) + asum(T, c0));
    if (!(@abs(@sqrt(got) - @sqrt(want)) <= tol)) {
        check.say("gglse {s} {d} x {d}, p {d}: residual {e}, c's part {e}\n", .{ @typeName(T), m, n, p, @as(f64, @floatCast(@sqrt(want))), @as(f64, @floatCast(@sqrt(got))) });
        return error.TestUnexpectedResult;
    }
    // not at an extended type, whose `nrm2` is `nrm2.S`'s unscaled sum of
    // squares and overflows there
    if (comptime types.isExtended(F)) return;
    const tnrm = @max(maxAbs(T, a0), maxAbs(T, b0), vmax(T, c0), vmax(T, d0));
    for (scales(F, tnrm), 0..) |sc, i| {
        runLse(T, al, a0, b0, c0, d0, sc, x, null) catch |e| {
            check.say("gglse {s} {s} {d} x {d}, p {d}: {t}\n", .{ sc.name, @typeName(T), m, n, p, e });
            return e;
        };
        // the right-hand sides alone scale x
        if (i == 2) scaleBy(T, x, 1 / sc.c);
        try expectLse(T, al, a0, b0, bq, c0, d0, x, sc.name);
    }
}

/// `dglmts.f`'s ratio, `‖d - A x - B y‖ / ((‖A‖ + ‖B‖)(‖x‖ + ‖y‖) eps)`,
/// and `y` in `Bᴴ null(Aᴴ)`, the minimum's condition.
fn expectGlm(comptime T: type, al: Allocator, a0: Mat(T), b0: Mat(T), w: ?Mat(T), d0: []const T, x: []const T, y: []const T, what: []const u8) !void {
    const F = Float(T);
    const eps = std.math.floatEps(F);
    const n = a0.rows;
    const m = a0.cols;
    const p = b0.cols;
    const r = try al.dupe(T, d0);
    defer al.free(r);
    subMul(T, a0, x, r);
    subMul(T, b0, y, r);
    const an = @max(norm1(T, a0), std.math.floatMin(F));
    const bn = @max(norm1(T, b0), std.math.floatMin(F));
    const xn = asum(T, x) + asum(T, y);
    const ratio = if (xn > 0) asum(T, r) / (an + bn) / xn / eps else 0;
    if (!(ratio <= lcheck.thresh)) {
        check.say("ggglm {s} {s} {d} x {d}, p {d}: ratio {e}\n", .{ what, @typeName(T), n, m, p, @as(f64, @floatCast(ratio)) });
        return error.TestUnexpectedResult;
    }
    const u = w orelse return;
    const tol = lcheck.thresh * @as(F, @floatFromInt(n + p)) * eps * (asum(T, y) + (asum(T, d0) + an * asum(T, x)) / bn);
    const off = outside(T, u, n - m, y);
    if (!(off <= tol)) {
        check.say("ggglm {s} {s} {d} x {d}, p {d}: y off Bᴴ null(Aᴴ) by {e}, tolerance {e}\n", .{ what, @typeName(T), n, m, p, @as(f64, @floatCast(off)), @as(f64, @floatCast(tol)) });
        return error.TestUnexpectedResult;
    }
}

fn runGlm(comptime T: type, al: Allocator, a0: Mat(T), b0: Mat(T), d0: []const T, sc: anytype, x: []T, y: []T) !void {
    const n = a0.rows;
    const a = try ops.alloc(T, al, n, a0.cols);
    defer ops.free(T, al, a);
    a.copyFrom(a0);
    scaleBy(T, a.slice(), sc.a);
    const b = try ops.alloc(T, al, n, b0.cols);
    defer ops.free(T, al, b);
    b.copyFrom(b0);
    scaleBy(T, b.slice(), sc.b);
    const d = try al.dupe(T, d0);
    defer al.free(d);
    scaleBy(T, d, sc.d);
    try gls.ggglm(al, check.io, T, n, a0.cols, b0.cols, a.slice(), a.ld, b.slice(), b.ld, d, x, y);
    scaleBy(T, x, sc.a / sc.d);
    scaleBy(T, y, sc.b / sc.d);
}

fn glmCase(comptime T: type, al: Allocator, n: usize, m: usize, p: usize, seed: u64) !void {
    const F = Float(T);
    const a0 = try ops.alloc(T, al, n, m);
    defer ops.free(T, al, a0);
    check.fill(T, a0.slice(), seed);
    const b0 = try ops.alloc(T, al, n, p);
    defer ops.free(T, al, b0);
    check.fill(T, b0.slice(), seed + 1);
    const d0 = try al.alloc(T, n);
    defer al.free(d0);
    check.fill(T, d0, seed + 2);
    // W = Bᴴ Q2, Q2 spanning null(Aᴴ); U's last p - (n - m) columns its
    // complement
    var u: ?Mat(T) = null;
    defer if (u) |q| ops.free(T, al, q);
    if (n > m) {
        const aq = try fullQ(T, al, a0);
        defer ops.free(T, al, aq);
        const w = try ops.alloc(T, al, p, n - m);
        defer ops.free(T, al, w);
        for (0..n - m) |j| for (0..p) |i| {
            var s = ops.zero(T);
            for (0..n) |k| s = ops.add(T, s, ops.mul(T, ops.conj(T, b0.get(k, i)), aq.get(k, m + j)));
            w.at(i, j).* = s;
        };
        u = try fullQ(T, al, w);
    }
    const x = try al.alloc(T, m);
    defer al.free(x);
    const y = try al.alloc(T, p);
    defer al.free(y);
    const one: F = 1;
    try runGlm(T, al, a0, b0, d0, .{ .a = one, .b = one, .d = one }, x, y);
    try expectGlm(T, al, a0, b0, u, d0, x, y, "unscaled");
    // not at an extended type, as `lseCase`
    if (comptime types.isExtended(F)) return;
    const tnrm = @max(maxAbs(T, a0), maxAbs(T, b0), vmax(T, d0));
    for (scales(F, tnrm)) |sc| {
        // `dglmts.f` scales `d` in every case
        runGlm(T, al, a0, b0, d0, .{ .a = sc.a, .b = sc.b, .d = @min(sc.a, sc.b, sc.c, sc.d) }, x, y) catch |e| {
            check.say("ggglm {s} {s} {d} x {d}, p {d}: {t}\n", .{ sc.name, @typeName(T), n, m, p, e });
            return e;
        };
        try expectGlm(T, al, a0, b0, u, d0, x, y, sc.name);
    }
}

test "lapack gglse: B x = d, the least-squares condition and the residual, scaled to the range's ends" {
    const al = testing.allocator;
    // (m, n, p), p <= n <= m + p
    const shapes = [_][3]usize{ .{ 1, 1, 1 }, .{ 6, 4, 2 }, .{ 4, 4, 4 }, .{ 3, 5, 2 }, .{ 2, 5, 3 }, .{ 5, 5, 0 }, .{ 10, 6, 3 }, .{ 40, 30, 10 } };
    inline for (lcheck.test_types) |T| for (shapes, 0..) |s, i| {
        try lseCase(T, al, s[0], s[1], s[2], 41 * i + 7);
    };
}

test "lapack ggglm: d = A x + B y with y minimal, scaled to the range's ends" {
    const al = testing.allocator;
    // (n, m, p), m <= n <= m + p
    const shapes = [_][3]usize{ .{ 1, 1, 1 }, .{ 4, 2, 3 }, .{ 4, 4, 0 }, .{ 5, 3, 2 }, .{ 6, 3, 8 }, .{ 3, 0, 3 }, .{ 30, 20, 40 } };
    inline for (lcheck.test_types) |T| for (shapes, 0..) |s, i| {
        try glmCase(T, al, s[0], s[1], s[2], 43 * i + 11);
    };
}

test "lapack gglse, ggglm: rank deficiency, empty, bad arguments" {
    const al = testing.allocator;
    var x: [3]f64 = @splat(9);
    var y: [3]f64 = @splat(9);
    // exact rank deficiencies, so the factors' zeros are exact
    // B's rows equal: rank(B) < p
    {
        var a = [_]f64{ 1, 0, 0, 1 };
        var b = [_]f64{ 0, 0, 1, 1 };
        var c = [_]f64{ 1, 1 };
        var d = [_]f64{ 1, 1 };
        try testing.expectError(error.BRankDeficient, gls.gglse(al, check.io, f64, 2, 2, 2, &a, 2, &b, 2, &c, &d, &x));
    }
    // A's and B's first columns zero: rank([A; B]) < n
    {
        var a = [_]f64{ 0, 0, 1, 1 };
        var b = [_]f64{ 0, 1 };
        var c = [_]f64{ 1, 1 };
        var d = [_]f64{1};
        try testing.expectError(error.ABRankDeficient, gls.gglse(al, check.io, f64, 2, 2, 1, &a, 2, &b, 1, &c, &d, &x));
    }
    // [A B] of rank 1 < n
    {
        var a = [_]f64{ 1, 0 };
        var b = [_]f64{ 1, 0 };
        var d = [_]f64{ 1, 1 };
        try testing.expectError(error.ABRankDeficient, gls.ggglm(al, check.io, f64, 2, 1, 1, &a, 2, &b, 2, &d, &x, &y));
    }
    // A's columns equal: rank(A) < m
    {
        var a = [_]f64{ 1, 0, 1, 0 };
        var b = [_]f64{ 0, 0 };
        var d = [_]f64{ 1, 1 };
        try testing.expectError(error.ARankDeficient, gls.ggglm(al, check.io, f64, 2, 2, 1, &a, 2, &b, 2, &d, &x, &y));
    }
    var a: [4]f64 = @splat(1);
    var b: [4]f64 = @splat(1);
    var v: [4]f64 = @splat(1);
    try gls.gglse(al, check.io, f64, 2, 0, 0, &a, 2, &b, 1, &v, &v, &x);
    try gls.ggglm(al, check.io, f64, 0, 0, 2, &a, 1, &b, 1, &v, &x, &y);
    try testing.expectEqualSlices(f64, &.{ 0, 0 }, y[0..2]);
    try testing.expectError(error.TooManyConstraints, gls.gglse(al, check.io, f64, 2, 1, 2, &a, 2, &b, 2, &v, &v, &x));
    try testing.expectError(error.TooFewRows, gls.gglse(al, check.io, f64, 1, 3, 1, &a, 1, &b, 1, &v, &v, &x));
    try testing.expectError(error.TooManyColumns, gls.ggglm(al, check.io, f64, 1, 2, 1, &a, 1, &b, 1, &v, &x, &y));
    try testing.expectError(error.TooFewColumns, gls.ggglm(al, check.io, f64, 3, 1, 1, &a, 3, &b, 3, &v, &x, &y));
    try testing.expectError(error.SliceTooShort, gls.gglse(al, check.io, f64, 2, 2, 1, &a, 2, &b, 1, &v, &v, x[0..1]));
    try testing.expectError(error.LeadingDimensionTooSmall, gls.ggglm(al, check.io, f64, 2, 1, 1, &a, 1, &b, 2, &v, &x, &y));
}
