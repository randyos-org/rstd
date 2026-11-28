//! Tests of `../svd.zig`.

const std = @import("std");
const ops = @import("../ops.zig");
const blas_calls = @import("../blas_calls.zig");
const svd_drivers = @import("../svd.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = ops.Mat;
const Float = ops.Float;
const SingularVectorJob = svd_drivers.SingularVectorJob;
const gesdd = svd_drivers.gesdd;
const gesvd = svd_drivers.gesvd;

const testing = std.testing;
const check = @import("../../blas/_test/check.zig");
const lcheck = @import("check.zig");
const types = @import("../../blas/types.zig");
const blas = @import("../../blas/__root__.zig");
const apple_accel = @import("../../blas/apple_accel.zig");
const test_io = @import("../../blas/_test/io.zig");

const test_types = .{ f32, f64, types.Complex32, types.Complex64 };

/// Faer's `test_svd` tolerance: `ApproxEq::eps()` (`128 eps`) times
/// `sqrt(8 max(m, n))`, here against `max(1, max |A|)` per element.
fn tolerance(comptime F: type, m: usize, n: usize, amax: F) F {
    return 128 * std.math.floatEps(F) * @sqrt(8 * @as(F, @floatFromInt(@max(m, n)))) * @max(amax, 1);
}

fn maxAbs(comptime T: type, a: Mat(T)) Float(T) {
    var x: Float(T) = 0;
    for (0..a.cols) |j| for (0..a.rows) |i| {
        x = @max(x, ops.abs(T, a.get(i, j)));
    };
    return x;
}

/// `‖A - U Σ Vᴴ‖` per element, `‖Uᴴ U - I‖` and `‖V Vᴴ - I‖` (over `Vᴴ`'s
/// rows) within `tolerance`; `s` nonincreasing and nonnegative.
fn expectSvd(comptime T: type, al: Allocator, a0: Mat(T), s: []const Float(T), u: Mat(T), vt: Mat(T), what: []const u8) !void {
    const F = Float(T);
    const m = a0.rows;
    const n = a0.cols;
    const k = s.len;
    const tol = tolerance(F, m, n, maxAbs(T, a0));
    for (0..k) |i| {
        if (!(s[i] >= 0) or (i > 0 and !(s[i - 1] >= s[i]))) {
            check.say("{s} {s} {d}x{d}: s[{d}] = {e} out of order\n", .{ what, @typeName(T), m, n, i, s[i] });
            return error.TestUnexpectedResult;
        }
    }
    // U Σ, then times Vᴴ
    const us = try ops.alloc(T, al, m, k);
    defer ops.free(T, al, us);
    for (0..k) |j| for (0..m) |i| {
        us.at(i, j).* = ops.mulReal(T, u.get(i, j), s[j]);
    };
    const r = try ops.alloc(T, al, m, n);
    defer ops.free(T, al, r);
    r.copyFrom(a0);
    try blas_calls.gemm(al, check.io, T, .n, .n, ops.neg(T, ops.one(T)), us, vt.sub(0, 0, k, n), ops.one(T), r);
    if (maxAbs(T, r) > tol) {
        check.say("{s} {s} {d}x{d}: A - U S Vh = {e}, tol {e}\n", .{ what, @typeName(T), m, n, maxAbs(T, r), tol });
        return error.TestUnexpectedResult;
    }
    const otol = tolerance(F, m, n, 1);
    const g = try ops.alloc(T, al, @max(u.cols, vt.rows), @max(u.cols, vt.rows));
    defer ops.free(T, al, g);
    const gu = g.sub(0, 0, u.cols, u.cols);
    try blas_calls.gemm(al, check.io, T, .c, .n, ops.one(T), u, u, ops.zero(T), gu);
    const gv = g.sub(0, 0, vt.rows, vt.rows);
    for ([_]Mat(T){ gu, gv }, 0..) |x, w| {
        if (w == 1) try blas_calls.gemm(al, check.io, T, .n, .c, ops.one(T), vt, vt, ops.zero(T), gv);
        for (0..x.cols) |j| for (0..x.rows) |i| {
            const want = if (i == j) ops.one(T) else ops.zero(T);
            if (!(ops.abs(T, ops.sub(T, x.get(i, j), want)) <= otol)) {
                check.say("{s} {s} {d}x{d}: {s} not orthonormal at ({d},{d})\n", .{ what, @typeName(T), m, n, if (w == 0) "U" else "Vh", i, j });
                return error.TestUnexpectedResult;
            }
        };
    }
}

const Routine = enum { gesvd, gesdd };

/// The routine on a copy of `a0` with jobs `ju`, `jv` (`gesdd`: `ju`), the
/// vectors gathered from wherever the jobs put them.
fn decompose(comptime T: type, al: Allocator, io: Io, r: Routine, ju: SingularVectorJob, jv: SingularVectorJob, a0: Mat(T), s: []Float(T), u: Mat(T), vt: Mat(T)) !void {
    const m = a0.rows;
    const n = a0.cols;
    const lda = m + 1;
    const a = try ops.alloc(T, al, lda, n);
    defer ops.free(T, al, a);
    a.fill(lcheck.poisoned(T));
    a.sub(0, 0, m, n).copyFrom(a0);
    switch (r) {
        .gesvd => try gesvd(al, io, T, ju, jv, m, n, a.slice(), lda, s, u.slice(), u.ld, vt.slice(), vt.ld, null),
        .gesdd => try gesdd(al, io, T, ju, m, n, a.slice(), lda, s, u.slice(), u.ld, vt.slice(), vt.ld),
    }
    const k = @min(m, n);
    const uo = if (r == .gesdd) ju == .overwrite and m >= n else ju == .overwrite;
    const vo = if (r == .gesdd) ju == .overwrite and m < n else jv == .overwrite;
    if (uo) u.sub(0, 0, m, k).copyFrom(a.sub(0, 0, m, k));
    if (vo) vt.sub(0, 0, k, n).copyFrom(a.sub(0, 0, k, n));
    // the padding row
    for (0..n) |j| try testing.expect(lcheck.isNanT(T, a.get(m, j)));
}

fn svdCases(comptime T: type, al: Allocator, r: Routine, m: usize, n: usize, a0: Mat(T)) !void {
    const k = @min(m, n);
    const s = try al.alloc(Float(T), @max(k, 1));
    defer al.free(s);
    const s2 = try al.alloc(Float(T), @max(k, 1));
    defer al.free(s2);
    const u = try ops.alloc(T, al, m, m);
    defer ops.free(T, al, u);
    const vt = try ops.alloc(T, al, n, n);
    defer ops.free(T, al, vt);
    const jobs = [_][2]SingularVectorJob{ .{ .all, .all }, .{ .leading, .leading }, .{ .overwrite, .leading }, .{ .leading, .overwrite }, .{ .all, .none }, .{ .none, .all } };
    const name = @tagName(r);
    for (jobs) |j| {
        if (r == .gesdd and j[0] != j[1] and !(j[0] == .overwrite)) continue;
        u.fill(lcheck.poisoned(T));
        vt.fill(lcheck.poisoned(T));
        try decompose(T, al, check.io, r, j[0], j[1], a0, s[0..k], u, vt);
        const gesdd_o = r == .gesdd and j[0] == .overwrite;
        const ucols: usize = if (gesdd_o) (if (m >= n) k else m) else switch (j[0]) {
            .all => m,
            .leading, .overwrite => k,
            .none => 0,
        };
        const vrows: usize = if (gesdd_o) (if (m >= n) n else k) else switch (if (r == .gesdd) j[0] else j[1]) {
            .all => n,
            .leading, .overwrite => k,
            .none => 0,
        };
        if (ucols > 0 and vrows > 0) try expectSvd(T, al, a0, s[0..k], u.sub(0, 0, m, ucols), vt.sub(0, 0, vrows, n), name);
    }
    // the values alone agree
    u.fill(lcheck.poisoned(T));
    vt.fill(lcheck.poisoned(T));
    try decompose(T, al, check.io, r, .none, .none, a0, s2[0..k], u, vt);
    try decompose(T, al, check.io, r, .all, .all, a0, s[0..k], u, vt);
    const tol = tolerance(Float(T), m, n, maxAbs(T, a0));
    for (s[0..k], s2[0..k]) |x, y| try testing.expect(@abs(x - y) <= tol);
}

/// Random `m x n`, faer's `test_real` and `test_cplx` shapes and more:
/// tall and wide past `11 / 6` (the QR path), and past 128 (divide and
/// conquer).
const shapes = [_][2]usize{ .{ 1, 1 }, .{ 2, 2 }, .{ 3, 2 }, .{ 2, 3 }, .{ 4, 4 }, .{ 15, 10 }, .{ 10, 15 }, .{ 17, 17 }, .{ 50, 50 }, .{ 150, 20 }, .{ 20, 150 }, .{ 140, 135 } };

test "faer svd: gesvd and gesdd, A = U S Vh with U and V orthonormal, every job, square, tall and wide" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        for (shapes) |mn| {
            const m = mn[0];
            const n = mn[1];
            const a0 = try ops.alloc(T, al, m, n);
            defer ops.free(T, al, a0);
            check.fill(T, a0.ptr[0 .. m * n], 13 * m + n);
            for ([_]Routine{ .gesvd, .gesdd }) |r| try svdCases(T, al, r, m, n, a0);
        }
    }
}

test "faer svd: special matrices, faer's test_special: zero, ones, identity; rank deficient; tiny and huge scales" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        for ([_][2]usize{ .{ 3, 2 }, .{ 2, 2 }, .{ 15, 10 }, .{ 10, 15 }, .{ 50, 50 }, .{ 150, 20 }, .{ 20, 150 }, .{ 130, 130 } }) |mn| {
            const m = mn[0];
            const n = mn[1];
            const a0 = try ops.alloc(T, al, m, n);
            defer ops.free(T, al, a0);
            for (0..7) |kind| {
                switch (kind) {
                    0 => a0.fill(ops.zero(T)),
                    1 => a0.fill(ops.one(T)),
                    2 => a0.identity(),
                    3 => try lcheck.rankProduct(T, al, a0, 3, 1),
                    4, 5, 6 => {
                        check.fill(T, a0.ptr[0 .. m * n], m + n);
                        // past the scaling bounds both ways, and a denormal
                        const sc: F = switch (kind) {
                            4 => std.math.floatMin(F) * 64,
                            5 => std.math.floatMax(F) / 64,
                            else => std.math.floatTrueMin(F) * 8,
                        };
                        if (kind == 6) {
                            a0.fill(ops.zero(T));
                            a0.at(m / 2, n / 2).* = ops.scalar(T, sc);
                        } else for (0..n) |j| for (a0.col(j)[0..m]) |*x| {
                            x.* = ops.mulReal(T, x.*, sc);
                        };
                    },
                    else => unreachable,
                }
                for ([_]Routine{ .gesvd, .gesdd }) |r| {
                    const k = @min(m, n);
                    const s = try al.alloc(F, k);
                    defer al.free(s);
                    const u = try ops.alloc(T, al, m, m);
                    defer ops.free(T, al, u);
                    const vt = try ops.alloc(T, al, n, n);
                    defer ops.free(T, al, vt);
                    decompose(T, al, check.io, r, .all, .all, a0, s, u, vt) catch |e| {
                        check.say("{s} {s} {d}x{d} kind {d}: {t}\n", .{ @tagName(r), @typeName(T), m, n, kind, e });
                        return e;
                    };
                    if (kind >= 4) {
                        // relative to the scale: compare the scaled-up values
                        // with those of the unscaled matrix
                        if (kind == 6) {
                            try testing.expectEqual(maxAbs(T, a0), s[0]);
                            for (s[1..]) |x| try testing.expectEqual(@as(F, 0), x);
                        } else {
                            const b = try ops.alloc(T, al, m, n);
                            defer ops.free(T, al, b);
                            check.fill(T, b.ptr[0 .. m * n], m + n);
                            const sb = try al.alloc(F, k);
                            defer al.free(sb);
                            try decompose(T, al, check.io, r, .none, .none, b, sb, u, vt);
                            const sc: F = if (kind == 4) std.math.floatMin(F) * 64 else std.math.floatMax(F) / 64;
                            const tol = tolerance(F, m, n, maxAbs(T, b));
                            for (s, sb) |x, y| try testing.expect(@abs(x / sc - y) <= tol);
                        }
                        continue;
                    }
                    try expectSvd(T, al, a0, s, u, vt, @tagName(r));
                    if (kind == 3) for (s[@min(3, k)..]) |x| try testing.expect(x <= tolerance(F, m, n, maxAbs(T, a0)));
                }
            }
        }
    }
}

test "faer svd: NaN and infinity: gesdd's NanInA (INFO = -4) for NaN, else NoConvergence" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex32 }) |T| {
        const F = Float(T);
        for ([_]F{ std.math.nan(F), std.math.inf(F), -std.math.inf(F) }) |bad| {
            var a: [12]T = undefined;
            check.fill(T, &a, 3);
            a[7] = ops.scalar(T, bad);
            var s: [3]F = undefined;
            var u: [16]T = undefined;
            var vt: [9]T = undefined;
            const want_dd: anyerror = if (std.math.isNan(bad)) error.NanInA else error.NoConvergence;
            var b = a;
            var count: usize = 99;
            try testing.expectError(error.NoConvergence, gesvd(al, check.io, T, .all, .all, 4, 3, &b, 4, &s, &u, 4, &vt, 3, &count));
            try testing.expectEqual(2, count);
            b = a;
            try testing.expectError(want_dd, gesdd(al, check.io, T, .all, 4, 3, &b, 4, &s, &u, 4, &vt, 3));
        }
    }
}

test "faer svd: gesvd's INFO > 0, the unconverged superdiagonals, QR iteration capped" {
    const al = testing.allocator;
    defer ops.test_iter_limit = null;
    inline for (test_types) |T| {
        const F = Float(T);
        // diag(A1, A2), 3 x 3 and 4 x 4: the bidiagonal splits between them
        // with an exact zero, so one sweep (on the bottom block) leaves
        // 2 + 3 of the 6 superdiagonals; uncapped, everything converges
        const n = 7;
        var a0: [n * n]T = undefined;
        check.fill(T, &a0, 17);
        for (0..n) |j| for (0..n) |i| {
            if ((i < 3) != (j < 3)) a0[i + j * n] = ops.zero(T);
        };
        var s: [n]F = undefined;
        var u: [n * n]T = undefined;
        var vt: [n * n]T = undefined;
        for ([_]?usize{ 1, null }) |lim| {
            ops.test_iter_limit = lim;
            var a = a0;
            var count: usize = 99;
            const r = gesvd(al, check.io, T, .all, .all, n, n, &a, n, &s, &u, n, &vt, n, &count);
            if (lim == null) {
                try r;
                try testing.expectEqual(99, count);
            } else {
                try testing.expectError(error.NoConvergence, r);
                try testing.expectEqual(5, count);
                // no out-parameter: the same error
                a = a0;
                try testing.expectError(error.NoConvergence, gesvd(al, check.io, T, .all, .none, n, n, &a, n, &s, &u, n, vt[0..0], 1, null));
                // values alone go to dqds, which the cap does not touch
                a = a0;
                try gesvd(al, check.io, T, .none, .none, n, n, &a, n, &s, u[0..0], 1, vt[0..0], 1, null);
            }
        }
    }
}

test "faer svd: bad arguments, and nothing read when empty" {
    const fa = testing.failing_allocator;
    const nan = std.math.nan(f64);
    var a: [4]f64 = @splat(nan);
    var s: [2]f64 = @splat(nan);
    var u: [4]f64 = @splat(nan);
    var vt: [4]f64 = @splat(nan);
    try testing.expectError(error.ConflictingJobs, gesvd(fa, check.io, f64, .overwrite, .overwrite, 2, 2, &a, 2, &s, &u, 2, &vt, 2, null));
    try testing.expectError(error.LeadingDimensionTooSmall, gesvd(fa, check.io, f64, .all, .all, 2, 2, &a, 1, &s, &u, 2, &vt, 2, null));
    try testing.expectError(error.SliceTooShort, gesvd(fa, check.io, f64, .none, .none, 2, 2, a[0..3], 2, &s, &u, 0, &vt, 0, null));
    try testing.expectError(error.SliceTooShort, gesvd(fa, check.io, f64, .none, .none, 2, 2, &a, 2, s[0..1], &u, 0, &vt, 0, null));
    try testing.expectError(error.LeadingDimensionTooSmall, gesvd(fa, check.io, f64, .all, .all, 2, 2, &a, 2, &s, &u, 1, &vt, 2, null));
    try testing.expectError(error.SliceTooShort, gesvd(fa, check.io, f64, .leading, .none, 2, 2, &a, 2, &s, u[0..3], 2, &vt, 0, null));
    try testing.expectError(error.LeadingDimensionTooSmall, gesvd(fa, check.io, f64, .none, .all, 2, 2, &a, 2, &s, &u, 0, &vt, 1, null));
    try testing.expectError(error.LeadingDimensionTooSmall, gesvd(fa, check.io, f64, .none, .leading, 2, 2, &a, 2, &s, &u, 0, &vt, 1, null));
    try testing.expectError(error.LeadingDimensionTooSmall, gesdd(fa, check.io, f64, .all, 2, 2, &a, 1, &s, &u, 2, &vt, 2));
    try testing.expectError(error.SliceTooShort, gesdd(fa, check.io, f64, .none, 2, 2, a[0..3], 2, &s, &u, 0, &vt, 0));
    try testing.expectError(error.LeadingDimensionTooSmall, gesdd(fa, check.io, f64, .all, 2, 2, &a, 2, &s, &u, 1, &vt, 2));
    try testing.expectError(error.LeadingDimensionTooSmall, gesdd(fa, check.io, f64, .overwrite, 1, 2, &a, 1, &s, &u, 0, &vt, 0));
    try testing.expectError(error.LeadingDimensionTooSmall, gesdd(fa, check.io, f64, .overwrite, 2, 2, &a, 2, &s, &u, 0, &vt, 1));
    try testing.expectError(error.LeadingDimensionTooSmall, gesdd(fa, check.io, f64, .leading, 2, 2, &a, 2, &s, &u, 2, &vt, 1));
    try testing.expectError(error.SliceTooShort, gesdd(fa, check.io, f64, .leading, 2, 2, &a, 2, &s, &u, 2, vt[0..3], 2));
    // empty: nothing read or written
    try gesvd(fa, check.io, f64, .all, .all, 0, 2, &a, 1, &s, &u, 1, &vt, 2, null);
    try gesdd(fa, check.io, f64, .all, 2, 0, &a, 2, &s, &u, 2, &vt, 1);
    for ([_][]const f64{ &a, &s, &u, &vt }) |x| for (x) |e| try testing.expect(std.math.isNan(e));
}

test "faer svd: 1 x 1 and a known 2 x 2, exact where the reference is" {
    const al = testing.allocator;
    var a = [_]f64{-3};
    var s: [1]f64 = undefined;
    var u: [1]f64 = undefined;
    var vt: [1]f64 = undefined;
    try gesvd(al, check.io, f64, .all, .all, 1, 1, &a, 1, &s, &u, 1, &vt, 1, null);
    try testing.expectEqual(@as(f64, 3), s[0]);
    try testing.expectEqual(@as(f64, -3), u[0] * s[0] * vt[0]);
    // [3 0; 4 5] has singular values sqrt(45) and sqrt(5)
    for ([_]Routine{ .gesvd, .gesdd }) |r| {
        var b = [_]f64{ 3, 4, 0, 5 };
        var s2: [2]f64 = undefined;
        var u_two: [4]f64 = undefined;
        var v2: [4]f64 = undefined;
        switch (r) {
            .gesvd => try gesvd(al, check.io, f64, .none, .none, 2, 2, &b, 2, &s2, &u_two, 1, &v2, 1, null),
            .gesdd => try gesdd(al, check.io, f64, .none, 2, 2, &b, 2, &s2, &u_two, 1, &v2, 1),
        }
        try testing.expectApproxEqRel(@sqrt(45.0), s2[0], 4e-16);
        try testing.expectApproxEqRel(@sqrt(5.0), s2[1], 4e-16);
    }
}

test "faer svd: gesdd gives one worker's bits at any worker count" {
    if (comptime !test_io.hosted) return error.SkipZigTest;
    _ = blas.coreName();
    // rstd's kernels: one worker's run is one chunk and the others many, and
    // Accelerate's bits follow a part's shape (its own split keeps the chunks)
    const was = apple_accel.routing();
    _ = apple_accel.route(false);
    defer _ = apple_accel.route(was);
    const al = testing.allocator;
    var runs: test_io.WorkerRuns = .init(al);
    defer runs.deinit();
    inline for (.{ f64, types.Complex32 }) |T| {
        // tall past 11 / 6, square past 128
        for ([_][2]usize{ .{ 300, 90 }, .{ 140, 140 } }) |mn| {
            const m = mn[0];
            const n = mn[1];
            const a0 = try ops.alloc(T, al, m, n);
            defer ops.free(T, al, a0);
            check.fill(T, a0.ptr[0 .. m * n], 99);
            const k = @min(m, n);
            var outs: [3][]T = undefined;
            var svals: [3][]Float(T) = undefined;
            for (0..3) |c| {
                outs[c] = try al.alloc(T, m * m + n * n);
                svals[c] = try al.alloc(Float(T), k);
            }
            defer for (0..3) |c| {
                al.free(outs[c]);
                al.free(svals[c]);
            };
            for (0..test_io.worker_runs.len) |c| {
                const io = runs.io(c);
                const u: Mat(T) = .{ .ptr = outs[c].ptr, .rows = m, .cols = m, .ld = m };
                const vt: Mat(T) = .{ .ptr = outs[c].ptr + m * m, .rows = n, .cols = n, .ld = n };
                try decompose(T, al, io, .gesdd, .all, .all, a0, svals[c], u, vt);
            }
            try test_io.expectSameBytes(outs);
            try test_io.expectSameBytes(svals);
        }
    }
}

test "gesvd, gesdd: a graded bidiagonal keeps its small singular values to relative accuracy (T126, M1's E4)" {
    const al = testing.allocator;
    inline for (.{ f64, f32 }) |T| {
        const F = Float(T);
        const n = 6;
        for (0..3) |kind| {
            // E4's: d(i) = 1e-4^i, e(i) = d(i) / 2 (f32: 1e-2^i); reversed; random grading
            var a: [n * n]T = @splat(0);
            var a128: [n * n]f128 = @splat(0);
            var prng = std.Random.DefaultPrng.init(kind + 1);
            const r = prng.random();
            const ratio: f64 = if (F == f64) 1e-4 else 1e-2;
            var g: f64 = 1;
            for (0..n) |i| {
                const k: f64 = @floatFromInt(if (kind == 1) n - 1 - i else i);
                const scale = if (kind == 2) g else std.math.pow(f64, ratio, k);
                const dv = scale * (if (kind == 2) 0.5 + r.float(f64) else 1);
                const ev = scale * (if (kind == 2) r.float(f64) else 0.5);
                a[i + i * n] = @floatCast(dv);
                a128[i + i * n] = @floatCast(@as(F, @floatCast(dv)));
                if (i + 1 < n) {
                    a[i + (i + 1) * n] = @floatCast(ev);
                    a128[i + (i + 1) * n] = @floatCast(@as(F, @floatCast(ev)));
                }
                if (kind == 2) g *= ratio * (0.5 + r.float(f64));
            }
            var s128: [n]f128 = undefined;
            var w128: [n * n]f128 = undefined;
            try gesvd(al, check.io, f128, .none, .none, n, n, &a128, n, &s128, &w128, 1, &w128, 1, null);
            for ([_]Routine{ .gesvd, .gesdd }) |rt| for ([_]bool{ false, true }) |vectors| {
                var b = a;
                var s: [n]F = undefined;
                var u: [n * n]T = undefined;
                var vt: [n * n]T = undefined;
                const job: SingularVectorJob = if (vectors) .all else .none;
                const ldv: usize = if (vectors) n else 1;
                switch (rt) {
                    .gesvd => try gesvd(al, check.io, T, job, job, n, n, &b, n, &s, &u, ldv, &vt, ldv, null),
                    .gesdd => try gesdd(al, check.io, T, job, n, n, &b, n, &s, &u, ldv, &vt, ldv),
                }
                for (s, s128, 0..) |x, ref, i| {
                    const err: f128 = @abs(@as(f128, x) - ref) / ref;
                    if (!(err <= 256 * @as(f128, n) * std.math.floatEps(F))) {
                        check.say("{s} {s} kind {d} vectors {}: s[{d}] = {e} against {e}, relative error {e}\n", .{ @tagName(rt), @typeName(T), kind, vectors, i, x, @as(f64, @floatCast(ref)), @as(f64, @floatCast(err)) });
                        return error.TestUnexpectedResult;
                    }
                }
            };
        }
    }
}
