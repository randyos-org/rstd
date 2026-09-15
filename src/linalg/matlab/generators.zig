//! Row reduction and matrix generation, ported from `matlab88.f90`'s RREF
//! and HILBER. See `blas1.zig` for the shared matrix/vector conventions
//! this port uses.

const blas1 = @import("blas1.zig");
const Cx = blas1.Cx;
const ix = blas1.ix;

/// RREF: reduce the M×N complex matrix `ar`/`ai` to reduced row-echelon
/// form in place, treating any entry no larger than
/// `eps * 2*max(M,N) * max_column_norm(A)` as a numerical zero.
pub fn rref(ar: []f64, ai: []f64, lda: usize, m: usize, n: usize, eps: f64) void {
    var tol: f64 = 0.0;
    var j: usize = 1;
    while (j <= n) : (j += 1) {
        tol = @max(tol, blas1.wasum(m, ar[ix(lda, 1, j)..], ai[ix(lda, 1, j)..], 1));
    }
    tol *= eps * @as(f64, @floatFromInt(2 * @max(m, n)));

    const lda_stride: isize = @intCast(lda);
    var k: usize = 1;
    var l: usize = 1;
    while (k <= m and l <= n) {
        const piv = blas1.iwamax(m - k + 1, ar[ix(lda, k, l)..], ai[ix(lda, k, l)..], 1) + k - 1;
        if (blas1.cabs1(ar[ix(lda, piv, l)], ai[ix(lda, piv, l)]) <= tol) {
            blas1.wset(m - k + 1, Cx.init(0.0, 0.0), ar[ix(lda, k, l)..], ai[ix(lda, k, l)..], 1);
            l += 1;
            continue;
        }
        blas1.wswap(n - l + 1, ar[ix(lda, piv, l)..], ai[ix(lda, piv, l)..], lda_stride, ar[ix(lda, k, l)..], ai[ix(lda, k, l)..], lda_stride);
        const inv = blas1.wdiv(Cx.init(1.0, 0.0), Cx.init(ar[ix(lda, k, l)], ai[ix(lda, k, l)]));
        blas1.wscal(n - l + 1, inv, ar[ix(lda, k, l)..], ai[ix(lda, k, l)..], lda_stride);
        ar[ix(lda, k, l)] = 1.0;
        ai[ix(lda, k, l)] = 0.0;

        var row: usize = 1;
        while (row <= m) : (row += 1) {
            if (row != k) {
                const t = Cx.init(-ar[ix(lda, row, l)], -ai[ix(lda, row, l)]);
                blas1.waxpy(n - l + 1, t, ar[ix(lda, k, l)..], ai[ix(lda, k, l)..], lda_stride, ar[ix(lda, row, l)..], ai[ix(lda, row, l)..], lda_stride);
            }
        }
        k += 1;
        l += 1;
    }
}

/// HILBER: generate the exact inverse of the N×N Hilbert matrix directly
/// into the real array `a` (column-major, leading dimension `lda`), using
/// the closed-form recurrence rather than computing it via LU inversion -
/// the Hilbert matrix is notoriously ill-conditioned, so an explicit
/// formula is far more accurate than inverting it numerically. `a` is
/// real-only; zero the imaginary part yourself if a complex result is
/// needed.
pub fn hilber(a: []f64, lda: usize, n: usize) void {
    var p: f64 = @floatFromInt(n);
    var i: usize = 1;
    while (i <= n) : (i += 1) {
        if (i != 1) {
            const im1: f64 = @floatFromInt(i - 1);
            p = @as(f64, @floatFromInt(n - i + 1)) * p * @as(f64, @floatFromInt(n + i - 1)) / (im1 * im1);
        }
        var r = p * p;
        a[ix(lda, i, i)] = r / @as(f64, @floatFromInt(2 * i - 1));
        if (i == n) continue;

        var j: usize = i + 1;
        while (j <= n) : (j += 1) {
            const jm1: f64 = @floatFromInt(j - 1);
            r = -(@as(f64, @floatFromInt(n - j + 1)) * r * @as(f64, @floatFromInt(n + j - 1))) / (jm1 * jm1);
            const v = r / @as(f64, @floatFromInt(i + j - 1));
            a[ix(lda, i, j)] = v;
            a[ix(lda, j, i)] = v;
        }
    }
}

const std = @import("std");
const testing = std.testing;

test "hilber matches the known inverse of the 3x3 Hilbert matrix" {
    // inv(hilb(3)) = [[9,-36,30],[-36,192,-180],[30,-180,180]].
    var a: [9]f64 = @splat(0.0);
    hilber(&a, 3, 3);
    try testing.expectApproxEqAbs(9.0, a[ix(3, 1, 1)], 1e-6);
    try testing.expectApproxEqAbs(-36.0, a[ix(3, 1, 2)], 1e-6);
    try testing.expectApproxEqAbs(30.0, a[ix(3, 1, 3)], 1e-6);
    try testing.expectApproxEqAbs(192.0, a[ix(3, 2, 2)], 1e-6);
    try testing.expectApproxEqAbs(-180.0, a[ix(3, 2, 3)], 1e-6);
    try testing.expectApproxEqAbs(180.0, a[ix(3, 3, 3)], 1e-6);
}

test "rref reduces a simple rank-deficient system" {
    // [[1, 2, 3], [2, 4, 6]] -> rank 1, row-reduces to [[1,2,3],[0,0,0]].
    var ar = [_]f64{ 1.0, 2.0, 2.0, 4.0, 3.0, 6.0 }; // column-major, 2x3
    var ai: [6]f64 = @splat(0.0);
    rref(&ar, &ai, 2, 2, 3, 1e-12);
    try testing.expectApproxEqAbs(1.0, ar[ix(2, 1, 1)], 1e-9);
    try testing.expectApproxEqAbs(2.0, ar[ix(2, 1, 2)], 1e-9);
    try testing.expectApproxEqAbs(3.0, ar[ix(2, 1, 3)], 1e-9);
    try testing.expectApproxEqAbs(0.0, ar[ix(2, 2, 1)], 1e-9);
    try testing.expectApproxEqAbs(0.0, ar[ix(2, 2, 2)], 1e-9);
    try testing.expectApproxEqAbs(0.0, ar[ix(2, 2, 3)], 1e-9);
}
