//! Tests of `../bidiag_svd.zig`, `../bdsqr.zig`, `../bidiag_secular.zig` and `../bidiag_deflate.zig`.

const std = @import("std");
const ops = @import("../ops.zig");
const hh = @import("../householder.zig");
const bidiag_deflate = @import("../bidiag_deflate.zig");
const bidiag_svd = @import("../bidiag_svd.zig");
const bdsqr = @import("../bdsqr.zig");

const Mat = ops.Mat;
const deflate = bidiag_deflate.deflate;
const divideAndConquer = bidiag_svd.divideAndConquer;
const qrAlgorithm = bdsqr.qrAlgorithm;
const sortAscending = bidiag_svd.sortAscending;

const testing = std.testing;
const check = @import("../../blas/_test/check.zig");

/// `U [S; 0] Vᵀ` against the `(n + 1) x n` (or with `square`, `n x n`)
/// lower bidiagonal matrix, within faer's `test_divide_and_conquer` (and
/// `test_qr_algorithm`) tolerance: `128 eps` times `max(|d|, |e|) sqrt n`,
/// times 10 for divide and conquer; and `Uᵀ U = I`, `Vᵀ V = I` within it.
fn expectBidiagSvd(comptime F: type, diag: []const F, sub: []const F, d: []const F, u: Mat(F), v: Mat(F), square: bool, factor: F) !void {
    const n = diag.len;
    const rows = if (square) n else n + 1;
    const scale = @max(ops.absMax(F, diag), ops.absMax(F, sub));
    const tol = 128 * std.math.floatEps(F) * @max(scale, 1) * @sqrt(@as(F, @floatFromInt(n))) * factor;
    for (0..n) |j| for (0..rows) |i| {
        var x: F = 0;
        for (0..n) |kk| x += u.get(i, kk) * d[kk] * v.get(j, kk);
        const want: F = if (i == j) diag[j] else if (i == j + 1) sub[j] else 0;
        if (!(@abs(x - want) <= tol)) {
            check.say("bidiag svd n={d}: ({d},{d}) = {e}, want {e}\n", .{ n, i, j, x, want });
            return error.TestUnexpectedResult;
        }
    };
    for (0..n) |j| for (0..n) |i| {
        var x: F = 0;
        for (0..n) |kk| x += v.get(kk, i) * v.get(kk, j);
        try testing.expect(@abs(x - @as(F, if (i == j) 1 else 0)) <= tol);
    };
    for (0..rows) |j| for (0..rows) |i| {
        var x: F = 0;
        for (0..rows) |kk| x += u.get(kk, i) * u.get(kk, j);
        try testing.expect(@abs(x - @as(F, if (i == j) 1 else 0)) <= tol);
    };
    for (1..n) |i| try testing.expect(d[i - 1] >= d[i]);
    for (d) |x| try testing.expect(x >= 0);
}

/// A bidiagonal with entries spread over several magnitudes, and clusters.
fn testBidiag(comptime F: type, diag: []F, sub: []F, seed: u64) void {
    var prng = std.Random.DefaultPrng.init(seed);
    const r = prng.random();
    for (diag, sub, 0..) |*d, *e, i| {
        d.* = (r.float(F) * 2 - 1) * std.math.pow(F, 10, @floatFromInt(@as(i32, @intCast(i % 5)) - 2));
        e.* = (r.float(F) * 2 - 1);
        if (i % 7 == 3) d.* = diag[i - 1];
    }
}

test "faer bidiag svd: qr_algorithm and divide_and_conquer reconstruct B, both thresholds" {
    const al = testing.allocator;
    inline for (.{ f64, f32 }) |F| {
        for ([_]usize{ 1, 2, 3, 7, 20, 64, 129 }) |n| {
            const diag = try al.alloc(F, n);
            defer al.free(diag);
            const sub = try al.alloc(F, n);
            defer al.free(sub);
            testBidiag(F, diag, sub, n);
            const d = try al.alloc(F, n);
            defer al.free(d);
            const e = try al.alloc(F, n);
            defer al.free(e);
            // qr_algorithm on the square part
            {
                @memcpy(d, diag);
                @memcpy(e, sub);
                e[n - 1] = 0;
                const sq = try al.dupe(F, e);
                defer al.free(sq);
                const u = try ops.alloc(F, al, n, n);
                defer ops.free(F, al, u);
                const v = try ops.alloc(F, al, n, n);
                defer ops.free(F, al, v);
                u.identity();
                v.identity();
                try qrAlgorithm(F, d, e, u, v);
                for (e) |x| try testing.expectEqual(@as(F, 0), x);
                try expectBidiagSvd(F, diag, sq, d, u, v, true, 1);
            }
            // divide and conquer, faer's threshold of 4 and 8, and 128
            for ([_]usize{ 4, 8, 128 }) |thr| {
                @memcpy(d, diag);
                @memcpy(e, sub);
                const u = try ops.alloc(F, al, n + 1, n + 1);
                defer ops.free(F, al, u);
                const v = try ops.alloc(F, al, n, n);
                defer ops.free(F, al, v);
                try divideAndConquer(al, check.io, F, d, e, .{ .full = u }, v, thr);
                expectBidiagSvd(F, diag, sub, d, u, v, false, 10) catch |err| {
                    check.say("{s} n={d} threshold {d}\n", .{ @typeName(F), n, thr });
                    return err;
                };
                // the values alone, and without U, agree
                const d2 = try al.dupe(F, diag);
                defer al.free(d2);
                const e2 = try al.dupe(F, sub);
                defer al.free(e2);
                const u_bot = try ops.alloc(F, al, 2, n + 1);
                defer ops.free(F, al, u_bot);
                try divideAndConquer(al, check.io, F, d2, e2, .{ .two_rows_storage = u_bot }, null, thr);
                for (d, d2) |x, y| try testing.expect(@abs(x - y) <= 1024 * std.math.floatEps(F) * @max(d[0], 1));
            }
        }
    }
}

test "faer bidiag svd: a zero block deflates, with or without U; zink's bidiagonal keeps its smallest value" {
    const al = testing.allocator;
    // the zero block is the second half of the first half, so its `TwoRows`
    // last row (faer: 0, the identity: eₙ) becomes the coupling `lambda`
    const n = 24;
    var diag: [n]f64 = undefined;
    var sub: [n]f64 = undefined;
    for (0..n) |i| {
        const zero = i >= 7 and i < 12;
        diag[i] = if (zero) 0 else 1 + @as(f64, @floatFromInt(i));
        sub[i] = if (zero) 0 else 0.5;
    }
    var d = diag;
    var e = sub;
    const u = try ops.alloc(f64, al, n + 1, n + 1);
    defer ops.free(f64, al, u);
    const v = try ops.alloc(f64, al, n, n);
    defer ops.free(f64, al, v);
    try divideAndConquer(al, check.io, f64, &d, &e, .{ .full = u }, v, 4);
    try expectBidiagSvd(f64, &diag, &sub, &d, u, v, false, 10);
    var d2 = diag;
    var e2 = sub;
    const u_bot = try ops.alloc(f64, al, 2, n + 1);
    defer ops.free(f64, al, u_bot);
    try divideAndConquer(al, check.io, f64, &d2, &e2, .{ .two_rows_storage = u_bot }, null, 4);
    // within a few ulps; faer's 0 is off by 1.5e-13 at d[5]
    for (d, d2) |x, y| try testing.expectApproxEqAbs(x, y, 8 * std.math.floatEps(f64) * d[0]);
    // faer's `test_zink` (svd/mod.rs)
    var zd = [_]f64{ -9.931833701529301, -10.920807536026027, -52.33647796311243, 2.3685025127736967, 2.421701994236093, -0.5051763005624579, -0.04808263896606017, -0.003875251886338955, -0.0006413264967716465, -0.003381944152463707, 2.981152313236375e-5, 5.4290648208388795e-6, -6.329275972084404e-7, -6.879142344209158e-7, -5.265228263479126e-9, -2.941999902335516e-9, -1.3060984997930294e-10, 7.07516117218088e-12, 1.8657003929029376e-12, -6.216080089659131e-14 };
    var zs = [_]f64{ -57.8029649868477, 17.67263066467847, 8.884153814270894, -9.01998231080713, -1.028638150814966, 0.22247719217200435, 0.016389886745811315, -0.004090989452162578, 0.00036818904090536926, -0.0031394146217732367, -7.571300829706796e-6, 3.0045718718618155e-6, 2.1329796886727743e-6, 9.259701025627789e-8, 2.2291214755992877e-9, -2.3017207713252894e-9, 6.807967994979358e-11, 2.1677299575405587e-12, -3.07282771050034e-13, 0.0 };
    const zu = try ops.alloc(f64, al, 2, zd.len + 1);
    defer ops.free(f64, al, zu);
    try divideAndConquer(al, check.io, f64, &zd, &zs, .{ .two_rows_storage = zu }, null, 8);
    try testing.expect(zd[zd.len - 1] != 0);
}

test "faer bidiag svd: deflation 4.3 and 4.4, faer's cases" {
    // `test_deflation43`, `test_deflation44`, `test_both_deflation`
    const al = testing.allocator;
    const cases = [_]struct { diag: [7]f32, want: [2]usize }{
        .{ .diag = .{ 1, 5, 3, 1e-7, 4, 2, 2e-7 }, .want = .{ 2, 0 } },
        .{ .diag = .{ 1, 5, 3, 1, 4, 2, 1 }, .want = .{ 0, 1 } },
        .{ .diag = .{ 1, 5, 3, 2, 4, 2, 0 }, .want = .{ 1, 1 } },
    };
    for (cases) |c| {
        var diag = c.diag;
        var col: [7]f32 = @splat(1);
        var coeff: [7]hh.Rotation(f32) = undefined;
        var idx: [7]usize = undefined;
        var perm: [7]usize = undefined;
        var tr: [7]usize = .{ 0, 1, 2, 3, 4, 5, 6 };
        const got = try deflate(al, f32, &diag, &col, &coeff, &idx, &tr, &perm, 3);
        try testing.expectEqual(c.want, got);
        // M' in the permuted places, then the rotations undone, is M: the
        // diagonal and a first column of ones
        var m: [7][7]f32 = @splat(@splat(0));
        for (0..7) |i| {
            m[perm[i]][perm[i]] = diag[i];
            m[perm[i]][perm[0]] = col[i];
        }
        var t = got[1];
        while (t > 0) {
            t -= 1;
            const rot = coeff[got[0] + t];
            const pi = perm[idx[got[0] + t] - 1];
            const pj = perm[idx[got[0] + t]];
            rot.left(7, &m[pj], 1, &m[pi], 1);
            rot.left(7, @as([*]f32, @ptrCast(&m)) + pj, 7, @as([*]f32, @ptrCast(&m)) + pi, 7);
        }
        t = got[0];
        while (t > 0) {
            t -= 1;
            coeff[t].left(7, &m[idx[t]], 1, &m[0], 1);
        }
        for (0..7) |i| for (0..7) |j| {
            const want: f32 = if (j == 0) 1 else if (i == j) c.diag[i] else 0;
            try testing.expect(@abs(m[i][j] - want) <= 128 * std.math.floatEps(f32) * @max(1, @abs(want)));
        };
    }
}

test "faer bidiag svd: the merged values sort fully, several out of place, ties kept in order (T130)" {
    // one bubble pass, faer's, would leave 1 ahead of 0.5 here
    var s = [_]f64{ 3, 1, 0.5, 2, 2, 0 };
    var perm = [_]usize{ 0, 1, 2, 3, 4, 5 };
    sortAscending(f64, &s, &perm);
    try testing.expectEqualSlices(f64, &.{ 0, 0.5, 1, 2, 2, 3 }, &s);
    try testing.expectEqualSlices(usize, &.{ 5, 2, 1, 3, 4, 0 }, &perm);
    // already sorted, and faer's case (one value a step ahead): as its pass leaves them
    var t = [_]f64{ 0, 1, 2 };
    var q = [_]usize{ 0, 1, 2 };
    sortAscending(f64, &t, &q);
    try testing.expectEqualSlices(usize, &.{ 0, 1, 2 }, &q);
    var u = [_]f64{ 0.5, 1, 3, 2 };
    var r = [_]usize{ 0, 1, 2, 3 };
    sortAscending(f64, &u, &r);
    try testing.expectEqualSlices(usize, &.{ 0, 1, 3, 2 }, &r);
}
