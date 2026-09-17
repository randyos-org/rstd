//! SVD-based conveniences (`COND`, `RANK`, `PINV`) and `ORTH`, split out of
//! `api.zig` to keep files under this project's line-count guideline.

const std = @import("std");
const zigstd = @import("../../zigstd.zig");
const blas1 = @import("blas1.zig");
const svd_mod = @import("svd.zig");
const qr_mod = @import("qr.zig");

const Cx = blas1.Cx;
const ix = blas1.ix;

fn singularValues(allocator: zigstd.Allocator, ar: []const f64, ai: []const f64, n: usize, p: usize) ![]f64 {
    const work_r = try allocator.alloc(f64, n * p);
    defer allocator.free(work_r);
    const work_i = try allocator.alloc(f64, n * p);
    defer allocator.free(work_i);
    @memcpy(work_r, ar[0 .. n * p]);
    @memcpy(work_i, ai[0 .. n * p]);

    const mm = @min(n + 1, p);
    const sr = try allocator.alloc(f64, mm);
    errdefer allocator.free(sr);
    const si = try allocator.alloc(f64, mm);
    defer allocator.free(si);
    const er = try allocator.alloc(f64, p);
    defer allocator.free(er);
    const ei = try allocator.alloc(f64, p);
    defer allocator.free(ei);
    const workr = try allocator.alloc(f64, n);
    defer allocator.free(workr);
    const worki = try allocator.alloc(f64, n);
    defer allocator.free(worki);
    var dummy: [1]f64 = .{0.0};

    _ = svd_mod.wsvdc(work_r, work_i, n, n, p, sr, si, er, ei, &dummy, &dummy, 1, &dummy, &dummy, 1, workr, worki, 0);
    return sr;
}

/// COND: ratio of the largest to smallest singular value of the `n`x`p`
/// matrix `ar`/`ai` (2-norm condition number). `A` is not mutated.
pub fn cond(allocator: zigstd.Allocator, ar: []const f64, ai: []const f64, n: usize, p: usize) !f64 {
    const s = try singularValues(allocator, ar, ai, n, p);
    defer allocator.free(s);
    const smallest = s[@min(n, p) - 1];
    if (smallest == 0.0) return std.math.inf(f64);
    return s[0] / smallest;
}

/// RANK: number of singular values of the `n`x`p` matrix `ar`/`ai` that
/// exceed `tol` (pass `null` to use the same default the original used:
/// `max(n, p) * eps * largest_singular_value`). `A` is not mutated.
pub fn rank(allocator: zigstd.Allocator, ar: []const f64, ai: []const f64, n: usize, p: usize, tol: ?f64) !usize {
    const s = try singularValues(allocator, ar, ai, n, p);
    defer allocator.free(s);
    const mn = @min(n, p);
    const threshold = tol orelse @as(f64, @floatFromInt(@max(n, p))) * std.math.floatEps(f64) * s[0];
    var k: usize = 0;
    for (s[0..mn]) |sv| {
        if (sv <= threshold) break;
        k += 1;
    }
    return k;
}

/// PINV: Moore-Penrose pseudoinverse of the `n`x`p` matrix `ar`/`ai`, via
/// SVD (`pinv(A) = V * pinv(S) * ctrans(U)`, treating singular values at or
/// below `tol` as zero - `null` uses the same default as `rank`). Writes
/// the `p`x`n` result into caller-supplied `result_r`/`result_i`
/// (column-major, leading dimension `p`). `A` is not mutated.
pub fn pinv(
    allocator: zigstd.Allocator,
    ar: []const f64,
    ai: []const f64,
    n: usize,
    p: usize,
    tol: ?f64,
    result_r: []f64,
    result_i: []f64,
) !void {
    const work_r = try allocator.alloc(f64, n * p);
    defer allocator.free(work_r);
    const work_i = try allocator.alloc(f64, n * p);
    defer allocator.free(work_i);
    @memcpy(work_r, ar[0 .. n * p]);
    @memcpy(work_i, ai[0 .. n * p]);

    const mn = @min(n, p);
    const mm = @min(n + 1, p);
    const sr = try allocator.alloc(f64, mm);
    defer allocator.free(sr);
    const si = try allocator.alloc(f64, mm);
    defer allocator.free(si);
    const er = try allocator.alloc(f64, p);
    defer allocator.free(er);
    const ei = try allocator.alloc(f64, p);
    defer allocator.free(ei);
    const ur = try allocator.alloc(f64, n * n);
    defer allocator.free(ur);
    const ui = try allocator.alloc(f64, n * n);
    defer allocator.free(ui);
    const vr = try allocator.alloc(f64, p * p);
    defer allocator.free(vr);
    const vi = try allocator.alloc(f64, p * p);
    defer allocator.free(vi);
    const workr = try allocator.alloc(f64, n);
    defer allocator.free(workr);
    const worki = try allocator.alloc(f64, n);
    defer allocator.free(worki);

    _ = svd_mod.wsvdc(work_r, work_i, n, n, p, sr, si, er, ei, ur, ui, n, vr, vi, p, workr, worki, 11);

    const threshold = tol orelse @as(f64, @floatFromInt(@max(n, p))) * std.math.floatEps(f64) * sr[0];
    var k: usize = 0;
    for (sr[0..mn]) |sv| {
        if (sv <= threshold) break;
        k += 1;
    }

    var i: usize = 1;
    while (i <= p) : (i += 1) {
        var j: usize = 1;
        while (j <= n) : (j += 1) {
            var re: f64 = 0.0;
            var im: f64 = 0.0;
            var l: usize = 1;
            while (l <= k) : (l += 1) {
                const v = Cx.init(vr[ix(p, i, l)], vi[ix(p, i, l)]);
                const u_conj = Cx.init(ur[ix(n, j, l)], -ui[ix(n, j, l)]);
                const term = blas1.wmul(v, u_conj);
                re += term.re / sr[l - 1];
                im += term.im / sr[l - 1];
            }
            result_r[ix(p, i, j)] = re;
            result_i[ix(p, i, j)] = im;
        }
    }
}

/// ORTH: `Q` from a QR decomposition of the `m`x`n` matrix `ar`/`ai`
/// (`n <= m`, not mutated) - an orthonormal basis (as columns) for its
/// column space. Caller owns the returned `m`x`n` slices.
pub fn orth(allocator: zigstd.Allocator, ar: []const f64, ai: []const f64, lda: usize, m: usize, n: usize) !struct { r: []f64, i: []f64 } {
    const work_r = try allocator.alloc(f64, m * n);
    defer allocator.free(work_r);
    const work_i = try allocator.alloc(f64, m * n);
    defer allocator.free(work_i);
    var j: usize = 1;
    while (j <= n) : (j += 1) {
        @memcpy(work_r[(j - 1) * m ..][0..m], ar[ix(lda, 1, j)..][0..m]);
        @memcpy(work_i[(j - 1) * m ..][0..m], ai[ix(lda, 1, j)..][0..m]);
    }

    const qrauxr = try allocator.alloc(f64, n);
    defer allocator.free(qrauxr);
    const qrauxi = try allocator.alloc(f64, n);
    defer allocator.free(qrauxi);
    const workr = try allocator.alloc(f64, n);
    defer allocator.free(workr);
    const worki = try allocator.alloc(f64, n);
    defer allocator.free(worki);
    const jpvt = try allocator.alloc(isize, n);
    defer allocator.free(jpvt);
    @memset(jpvt, 0);

    qr_mod.wqrdc(work_r, work_i, m, m, n, qrauxr, qrauxi, jpvt, workr, worki, 0);

    const result_r = try allocator.alloc(f64, m * n);
    errdefer allocator.free(result_r);
    const result_i = try allocator.alloc(f64, m * n);
    errdefer allocator.free(result_i);

    var dummy: [1]f64 = .{0.0};
    const basis_r = try allocator.alloc(f64, m);
    defer allocator.free(basis_r);
    const basis_i = try allocator.alloc(f64, m);
    defer allocator.free(basis_i);
    const qy_r = try allocator.alloc(f64, m);
    defer allocator.free(qy_r);
    const qy_i = try allocator.alloc(f64, m);
    defer allocator.free(qy_i);

    var col: usize = 1;
    while (col <= m) : (col += 1) {
        @memset(basis_r, 0.0);
        @memset(basis_i, 0.0);
        basis_r[col - 1] = 1.0;

        _ = qr_mod.wqrsl(
            work_r,
            work_i,
            m,
            m,
            n,
            qrauxr,
            qrauxi,
            basis_r,
            basis_i,
            qy_r,
            qy_i,
            &dummy,
            &dummy,
            &dummy,
            &dummy,
            &dummy,
            &dummy,
            &dummy,
            &dummy,
            10000,
        );
        if (col <= n) {
            @memcpy(result_r[(col - 1) * m ..][0..m], qy_r);
            @memcpy(result_i[(col - 1) * m ..][0..m], qy_i);
        }
    }
    return .{ .r = result_r, .i = result_i };
}


test "cond, rank and pinv agree for a known rank-deficient matrix" {
    // A = [[1,2],[2,4]] (rank 1, columns parallel).
    const ar = [_]f64{ 1.0, 2.0, 2.0, 4.0 };
    const ai: [4]f64 = @splat(0.0);

    const c = try cond(zigstd.testAllocator, &ar, &ai, 2, 2);
    try zigstd.testExpectTrue(c > 1e10 or std.math.isInf(c));

    const r = try rank(zigstd.testAllocator, &ar, &ai, 2, 2, null);
    try zigstd.testExpectEqual(@as(usize, 1), r);

    var pr: [4]f64 = @splat(0.0);
    var pi: [4]f64 = @splat(0.0);
    try pinv(zigstd.testAllocator, &ar, &ai, 2, 2, null, &pr, &pi);
    // A * pinv(A) * A should reconstruct A (pseudoinverse identity).
    var apa_full: [4]f64 = @splat(0.0);
    for (0..2) |row| {
        for (0..2) |col| {
            var s: f64 = 0.0;
            for (0..2) |m| {
                var ap: f64 = 0.0;
                for (0..2) |l| ap += ar[ix(2, row + 1, l + 1)] * pr[ix(2, l + 1, m + 1)];
                s += ap * ar[ix(2, m + 1, col + 1)];
            }
            apa_full[ix(2, row + 1, col + 1)] = s;
        }
    }
    for (0..4) |idx| try zigstd.testExpectApproxEqAbs(ar[idx], apa_full[idx], 1e-6);
}

test "orth returns an orthonormal basis spanning the column space" {
    // A = [[1,1],[0,1],[1,0]] (3x2, full column rank).
    const ar = [_]f64{ 1.0, 0.0, 1.0, 1.0, 1.0, 0.0 };
    const ai: [6]f64 = @splat(0.0);
    const q = try orth(zigstd.testAllocator, &ar, &ai, 3, 3, 2);
    defer zigstd.testAllocator.free(q.r);
    defer zigstd.testAllocator.free(q.i);

    for (0..2) |a_col| {
        for (0..2) |b_col| {
            var dot: f64 = 0.0;
            for (0..3) |row| dot += q.r[ix(3, row + 1, a_col + 1)] * q.r[ix(3, row + 1, b_col + 1)];
            const expected: f64 = if (a_col == b_col) 1.0 else 0.0;
            try zigstd.testExpectApproxEqAbs(expected, dot, 1e-9);
        }
    }
}
