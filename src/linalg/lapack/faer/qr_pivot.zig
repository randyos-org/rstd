//! QR with column pivoting, from faer v0.24.4
//! (<https://codeberg.org/sarah-quinones/faer>, tag `faer-v0.24.4`, MIT):
//! `faer/src/linalg/qr/col_pivoting/factor.rs`'s `qr_in_place_unblocked`,
//! behind Reference-LAPACK's `?GEQP3` (`dgeqp3.f`, `zgeqp3.f`, 3.12.0 as
//! OpenBLAS v0.3.34-317 ships it: arguments, `JPVT`'s fixed columns, `tau`).
//!
//! `A` is scaled by its largest column norm; each step takes the column of
//! largest remaining norm, reflects it (`householder.make`) and updates the
//! next row at once and the rows below one step later, when every norm is
//! recomputed in full (faer's path without its SIMD-fused update, the same
//! arithmetic otherwise). The reference's `xLAQPS` downdates norms instead.
//! Columns `JPVT` fixes lead and are factored unpivoted, as `dgeqp3.f`'s
//! `xGEQRF` and `xORMQR` do. `tau` is LAPACK's, `1 / faer's` (0 where faer
//! skips a reflection); a complex `R`'s diagonal is not made real. One fix:
//! faer scales a zero matrix by an infinite reciprocal (NaN); here it is
//! left as is. Column norms are rstd's `nrm2`, at an extended type
//! `nrm2.S`'s unscaled sum of squares, so there a matrix near the underflow
//! threshold reads as zero.

const std = @import("std");
const types = @import("../../blas/types.zig");
const util = @import("../../blas/util.zig");
const ops = @import("ops.zig");
const hh = @import("householder.zig");
const chol = @import("../cholesky.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = ops.Mat;
const Float = types.Float;
const Error = ops.Error;

/// What `geqp3` returns.
pub const Geqp3Error = error{ InvalidM, InvalidN, InvalidLda } || Allocator.Error || Io.Cancelable;

/// `?GEQP3`: `A P = Q R` for the `m x n` `A`: `R` on and above the
/// diagonal, `Q = H(0) ... H(k-1)` (`k = min(m, n)`) as reflections below
/// it and their `tau`s in `tau` (`k`). On entry a nonzero `jpvt[j]` fixes
/// column `j` to the front; on exit `jpvt[j] = i` (1-based) when column `j`
/// of `A P` was column `i` of `A`.
pub fn geqp3(allocator: Allocator, io: Io, comptime T: type, m: isize, n: isize, a: [*]T, lda: isize, jpvt: [*]isize, tau: [*]T) Geqp3Error!void {
    comptime chol.floatType(T, "geqp3");
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, m)) return error.InvalidLda;
    if (n == 0) return;
    if (!util.fitsMatrix(m, n, lda, @sizeOf(T))) return error.InvalidLda;
    const am: Mat(T) = .{ .ptr = a, .rows = @intCast(m), .cols = @intCast(n), .ld = @intCast(lda) };
    try factor(allocator, io, T, am, jpvt, tau);
}

/// `geqp3` after its checks.
pub fn factor(al: Allocator, io: Io, comptime T: type, a: Mat(T), jpvt: [*]isize, tau: [*]T) Error!void {
    const n = a.cols;
    // `dgeqp3.f`: the fixed columns to the front
    var nfxd: usize = 0;
    for (0..n) |j| {
        if (jpvt[j] != 0) {
            if (j != nfxd) {
                swapCols(T, a, j, nfxd);
                jpvt[j] = jpvt[nfxd];
                jpvt[nfxd] = @intCast(j + 1);
            } else {
                jpvt[j] = @intCast(j + 1);
            }
            nfxd += 1;
        } else {
            jpvt[j] = @intCast(j + 1);
        }
    }
    if (a.rows == 0) return;
    try unblocked(al, io, T, a, nfxd, jpvt, tau);
}

fn swapCols(comptime T: type, a: Mat(T), i: usize, j: usize) void {
    for (0..a.rows) |r| std.mem.swap(T, a.at(r, i), a.at(r, j));
}

/// `qr_in_place_unblocked`, columns before `nfxd` unpivoted.
fn unblocked(al: Allocator, io: Io, comptime T: type, a: Mat(T), nfxd: usize, jpvt: [*]isize, tau: [*]T) Error!void {
    const F = Float(T);
    const m = a.rows;
    const n = a.cols;
    const size = @min(m, n);
    const dot = try al.alloc(T, n);
    defer al.free(dot);
    @memset(dot, ops.zero(T));
    const norm = try al.alloc(F, n);
    defer al.free(norm);
    var best: F = 0;
    for (0..n) |j| {
        norm[j] = ops.nrm2(io, T, m, a.col(j), 1);
        if (norm[j] > best) best = norm[j];
    }
    // scale by the largest norm (a zero matrix stays)
    const scale_fwd = best;
    const scale_bwd: F = if (best > 0) 1 / best else 1;
    for (0..n) |j| for (a.col(j)[0..m]) |*x| {
        x.* = ops.mulReal(T, x.*, scale_bwd);
    };
    for (norm) |*x| x.* *= scale_bwd;
    for (0..size) |k| {
        var best_col = k;
        if (k >= nfxd) {
            var new_best: F = 0;
            for (k..n) |j| if (norm[j] > new_best) {
                new_best = norm[j];
                best_col = j;
            };
        }
        if (k > 0) {
            // the previous reflection's update of the rows below, then the norms
            const a11 = a.sub(k, k, m - k, n - k);
            try ops.ger(al, io, T, ops.one(T), a.col(k - 1) + k, 1, dot.ptr + k, 1, a11);
            best = 0;
            for (k..n) |j| {
                norm[j] = ops.nrm2(io, T, m - k, a.col(j) + k, 1);
                if (k >= nfxd and norm[j] > best) {
                    best = norm[j];
                    best_col = j;
                }
            }
        }
        if (best_col != k) {
            std.mem.swap(isize, &jpvt[best_col], &jpvt[k]);
            swapCols(T, a, best_col, k);
            std.mem.swap(T, &dot[best_col], &dot[k]);
            std.mem.swap(F, &norm[best_col], &norm[k]);
        }
        const tail = m - k - 1;
        const v = a.col(k) + k + 1;
        const info = hh.make(io, T, a.at(k, k), tail, v, 1, null, 0);
        const tau_inv = 1 / info.tau;
        tau[k] = ops.scalar(T, tau_inv);
        if (k + 1 == size) break;
        // dot := A01 + vᴴ A11, then -dot / tau into the next row
        const rest = n - k - 1;
        for (0..rest) |j| dot[k + 1 + j] = a.get(k, k + 1 + j);
        try ops.gemv(al, io, T, .u, ops.one(T), a.sub(k + 1, k + 1, tail, rest), v, 1, ops.one(T), dot.ptr + k + 1, 1);
        for (k + 1..n) |j| {
            dot[j] = ops.neg(T, ops.mulReal(T, dot[j], tau_inv));
            const p = a.at(k, j);
            p.* = ops.add(T, p.*, dot[j]);
        }
    }
    // `R` back to scale
    if (scale_fwd > 0) for (0..n) |j| for (0..@min(j + 1, m)) |i| {
        a.at(i, j).* = ops.mulReal(T, a.get(i, j), scale_fwd);
    };
}

// ---- tests -----------------------------------------------------------------------

const testing = std.testing;
const check = @import("../../blas/_test/blas1/check.zig");
const qr = @import("qr.zig");

const test_types = chol.test_types;

/// `THRESH` of the reference's linear-equation tests (`TESTING/dtest.in`).
const thresh = 30;

/// `dqpt01.f`'s `‖A P - Q R‖ / (m ‖A‖ eps)` and `dqrt11.f`'s `‖I - Qᴴ Q‖ /
/// (m eps)` below the threshold, `|R(k, k)|` not increasing past the fixed
/// columns, and the fixed columns first in order.
fn qpCase(comptime T: type, al: Allocator, m: usize, n: usize, rank: usize, fixed: []const usize, seed: u64) !void {
    const F = Float(T);
    const eps = std.math.floatEps(F);
    const k = @min(m, n);
    const a0 = try ops.alloc(T, al, m, n);
    defer ops.free(T, al, a0);
    // rank `rank`: a product of random m x rank and rank x n
    {
        const x = try al.alloc(T, m * rank);
        defer al.free(x);
        const y = try al.alloc(T, rank * n);
        defer al.free(y);
        check.fill(T, x, seed);
        check.fill(T, y, seed + 1);
        for (0..n) |j| for (0..m) |i| {
            var s = ops.zero(T);
            for (0..rank) |p| s = ops.add(T, s, ops.mul(T, x[i + p * m], y[p + j * rank]));
            a0.at(i, j).* = s;
        };
    }
    const a = try ops.alloc(T, al, m, n);
    defer ops.free(T, al, a);
    a.copyFrom(a0);
    const jpvt = try al.alloc(isize, n);
    defer al.free(jpvt);
    @memset(jpvt, 0);
    for (fixed) |f| jpvt[f] = 1;
    const tau = try al.alloc(T, @max(k, 1));
    defer al.free(tau);
    try geqp3(al, check.io, T, @intCast(m), @intCast(n), a.ptr, @intCast(a.ld), jpvt.ptr, tau.ptr);
    // the fixed columns lead, in their order
    var nf: usize = 0;
    for (0..n) |j| if (std.mem.indexOfScalar(usize, fixed, j) != null) {
        try testing.expectEqual(@as(isize, @intCast(j + 1)), jpvt[nf]);
        nf += 1;
    };
    // jpvt is a permutation
    for (1..n + 1) |v| try testing.expect(std.mem.indexOfScalar(isize, jpvt, @intCast(v)) != null);
    // Q explicitly, then Q R against A P
    const q = try ops.alloc(T, al, m, m);
    defer ops.free(T, al, q);
    q.identity();
    if (k > 0) {
        const v = a.sub(0, 0, m, k);
        const h = try qr.factorsOf(al, check.io, T, v, tau.ptr);
        defer ops.free(T, al, h);
        try hh.applySequence(al, check.io, T, v, h, false, q);
    }
    var anorm: F = 0;
    for (a0.ptr[0 .. m * n]) |x| anorm = @max(anorm, ops.abs(T, x));
    const mf: F = @floatFromInt(@max(m, 1));
    for (0..n) |j| for (0..m) |i| {
        var s = ops.zero(T);
        for (0..@min(j + 1, m)) |p| s = ops.add(T, s, ops.mul(T, q.get(i, p), a.get(p, j)));
        const d = ops.abs(T, ops.sub(T, s, a0.get(i, @intCast(jpvt[j] - 1))));
        if (!(d <= thresh * mf * eps * @max(anorm, std.math.floatMin(F)))) {
            check.say("geqp3 {s} {d} x {d} rank {d}: A P - Q R at ({d}, {d}) is {e}\n", .{ @typeName(T), m, n, rank, i, j, @as(f64, @floatCast(d)) });
            return error.TestUnexpectedResult;
        }
    };
    for (0..m) |j| for (0..m) |i| {
        var s = ops.zero(T);
        for (0..m) |p| s = ops.add(T, s, ops.mul(T, ops.conj(T, q.get(p, i)), q.get(p, j)));
        if (i == j) s = ops.sub(T, s, ops.one(T));
        try testing.expect(ops.abs(T, s) <= thresh * mf * eps);
    };
    // pivoting: the diagonal does not grow
    for (nf + 1..k) |i| {
        const prev = ops.abs(T, a.get(i - 1, i - 1));
        try testing.expect(ops.abs(T, a.get(i, i)) <= prev * (1 + 16 * mf * eps) + thresh * mf * eps * anorm);
    }
}

test "lapack geqp3: A P = Q R, Q unitary, a non-increasing diagonal, fixed columns, full rank and rank deficient" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        for ([_][3]usize{ .{ 1, 1, 1 }, .{ 4, 4, 4 }, .{ 8, 5, 5 }, .{ 5, 8, 5 }, .{ 20, 12, 7 }, .{ 12, 20, 3 }, .{ 30, 30, 30 } }) |c| {
            try qpCase(T, al, c[0], c[1], c[2], &.{}, c[0] * 13 + c[2]);
            if (c[1] >= 4) try qpCase(T, al, c[0], c[1], c[2], &.{ 1, 3 }, c[0] * 17 + c[2]);
        }
        // zero
        try qpCase(T, al, 6, 4, 0, &.{}, 3);
    }
}

test "lapack geqp3: empty and bad arguments" {
    const al = testing.allocator;
    var a: [4]f64 = @splat(1);
    var jpvt = [_]isize{ 5, 5 };
    var tau: [2]f64 = undefined;
    try geqp3(al, check.io, f64, 0, 2, &a, 1, &jpvt, &tau);
    try testing.expectEqualSlices(isize, &.{ 1, 2 }, &jpvt);
    try geqp3(al, check.io, f64, 2, 0, &a, 2, &jpvt, &tau);
    try testing.expectError(error.InvalidM, geqp3(al, check.io, f64, -1, 0, &a, 1, &jpvt, &tau));
    try testing.expectError(error.InvalidN, geqp3(al, check.io, f64, 1, -1, &a, 1, &jpvt, &tau));
    try testing.expectError(error.InvalidLda, geqp3(al, check.io, f64, 2, 2, &a, 1, &jpvt, &tau));
}
