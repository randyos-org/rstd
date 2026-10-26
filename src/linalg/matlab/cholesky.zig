//! Complex Cholesky factorization, ported from `matlab88.f90`'s WPOFA
//! (LINPACK's ZPOFA). See `blas1.zig` for the shared matrix/vector
//! conventions this port uses, and `CxMat.zig` for the `self`-shaped
//! calling convention matrix-level functions like this one now use.

const blas1 = @import("blas1.zig");
const CxMat = @import("CxMat.zig");
const Cx = blas1.Cx;
const ix = blas1.ix;

/// WPOFA: Cholesky-factor a Hermitian positive-definite matrix, using only
/// the diagonal and (strict) upper triangle of `self` on input - the lower
/// triangle is assumed to be their conjugate transpose and is never read.
/// On success, the upper triangle of `self` holds `R` such that
/// `ctrans(R)*R == A`; the strict lower triangle is left unmodified (zero
/// it yourself if a clean triangular matrix is required).
///
/// Returns 0 on success, or the 1-based order of the leading minor that
/// first failed to be positive definite.
pub fn wpofa(self: CxMat) usize {
    const ar = self.r;
    const ai = self.i;
    const lda = self.lda;
    const n = self.n;

    var j: usize = 1;
    while (j <= n) : (j += 1) {
        var s: f64 = 0.0;
        if (j >= 2) {
            var k: usize = 1;
            while (k <= j - 1) : (k += 1) {
                const tr = ar[ix(lda, k, j)] - blas1.wdotcr(k - 1, ar[ix(lda, 1, k)..], ai[ix(lda, 1, k)..], 1, ar[ix(lda, 1, j)..], ai[ix(lda, 1, j)..], 1);
                const ti = ai[ix(lda, k, j)] - blas1.wdotci(k - 1, ar[ix(lda, 1, k)..], ai[ix(lda, 1, k)..], 1, ar[ix(lda, 1, j)..], ai[ix(lda, 1, j)..], 1);
                const d = blas1.wdiv(Cx.init(tr, ti), Cx.init(ar[ix(lda, k, k)], ai[ix(lda, k, k)]));
                ar[ix(lda, k, j)] = d.re;
                ai[ix(lda, k, j)] = d.im;
                s += d.squaredMagnitude();
            }
        }
        s = ar[ix(lda, j, j)] - s;
        if (s <= 0.0 or ai[ix(lda, j, j)] != 0.0) return j;
        ar[ix(lda, j, j)] = @sqrt(s);
    }
    return 0;
}

const std = @import("std");

const zigstd = @import("../../zigstd.zig");

test "wpofa factors a known real SPD matrix" {
    // A = [[4, 2], [2, 3]] (upper triangle supplied, real).
    var ar = [_]f64{ 4.0, 0.0, 2.0, 3.0 }; // column-major; lower-left unused
    var ai: [4]f64 = @splat(0.0);
    const mat = CxMat.initPacked(&ar, &ai, 2, 2);
    const info = wpofa(mat);
    try zigstd.testExpectEqual(@as(usize, 0), info);

    // R = [[2, 1], [0, sqrt(2)]]; check ctrans(R)*R == A.
    const r11 = ar[ix(2, 1, 1)];
    const r12 = ar[ix(2, 1, 2)];
    const r22 = ar[ix(2, 2, 2)];
    try zigstd.testExpectApproxEqAbs(4.0, r11 * r11, 1e-10);
    try zigstd.testExpectApproxEqAbs(2.0, r11 * r12, 1e-10);
    try zigstd.testExpectApproxEqAbs(3.0, r12 * r12 + r22 * r22, 1e-10);
}

test "wpofa rejects a non positive-definite matrix" {
    var ar = [_]f64{ 1.0, 0.0, 5.0, 1.0 };
    var ai: [4]f64 = @splat(0.0);
    const mat = CxMat.initPacked(&ar, &ai, 2, 2);
    const info = wpofa(mat);
    try zigstd.testExpectEqual(@as(usize, 2), info);
}
