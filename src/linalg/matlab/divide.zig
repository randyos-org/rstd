//! Matrix left/right division (REPL `\`/`/`), kept as its own file rather
//! than folded into `CxMat.zig` to keep files under this project's
//! line-count guideline.

const std = @import("std");
const zigstd = @import("../../zigstd.zig");
const blas1 = @import("blas1.zig");
const CxMat = @import("CxMat.zig");
const lu_mod = @import("lu.zig");

const Cx = blas1.Cx;
const ix = blas1.ix;

/// Right matrix division for a square divisor: solve `X*B = A` for `X`,
/// via `B`'s LU factorization and a conjugate-transpose solve applied to
/// each row of `A` (`X*B=A` iff `ctrans(B)*ctrans(X) = ctrans(A)`, solved
/// one row/column at a time so no explicit transpose is ever materialized)
/// - REPL semantics for `A/B`. Takes two co-equal matrices rather than
/// mutating a single `self` the way single-matrix operations do: `a`
/// (`m`x`n`) is overwritten with `X`; `b` (`n`x`n`) is destroyed by
/// factorization. `ipvt`, `zr`/`zi`, and `rowr`/`rowi` are length-`n`
/// scratch. Returns `B`'s `rcond` so callers can judge accuracy.
pub fn mrdivide(
    a: CxMat,
    b: CxMat,
    ipvt: []usize,
    zr: []f64,
    zi: []f64,
    rowr: []f64,
    rowi: []f64,
) f64 {
    const ar = a.r;
    const ai = a.i;
    const lda = a.lda;
    const m = a.m;
    const n = a.n;

    const rc = lu_mod.wgeco(b, ipvt, zr, zi);
    var i: usize = 1;
    while (i <= m) : (i += 1) {
        var j: usize = 1;
        while (j <= n) : (j += 1) {
            rowr[j - 1] = ar[ix(lda, i, j)];
            rowi[j - 1] = -ai[ix(lda, i, j)];
        }
        lu_mod.wgesl(b, ipvt, rowr, rowi, 1);
        j = 1;
        while (j <= n) : (j += 1) {
            ar[ix(lda, i, j)] = rowr[j - 1];
            ai[ix(lda, i, j)] = -rowi[j - 1];
        }
    }
    return rc;
}

/// Left matrix division for a square `self`: solve `self*X = b` for `X` -
/// REPL semantics for `A\B` (square `A` only; a non-square least-squares
/// solve can be built directly from `wqrdc`/`wqrsl`, matching `\`'s own
/// documented fallback). `self` is destroyed by factorization; `b`
/// (`self.n`x`p`) is overwritten with `X`. `ipvt`, `zr`/`zi` (length
/// `self.n`) are scratch. Returns `self`'s `rcond` so callers can judge
/// accuracy.
pub fn mldivide(self: CxMat, b: CxMat, ipvt: []usize, zr: []f64, zi: []f64) f64 {
    const rc = lu_mod.wgeco(self, ipvt, zr, zi);
    var j: usize = 1;
    while (j <= b.n) : (j += 1) {
        const c = b.col(j);
        lu_mod.wgesl(self, ipvt, c.r, c.i, 0);
    }
    return rc;
}


test "mldivide solves A*X=B" {
    // A = [[2,0],[1,2]], B = [[1,2],[3,4]] (both square, real).
    var ar = [_]f64{ 2.0, 1.0, 0.0, 2.0 }; // column-major
    var ai: [4]f64 = @splat(0.0);
    var br = [_]f64{ 1.0, 3.0, 2.0, 4.0 };
    var bi: [4]f64 = @splat(0.0);
    var ipvt: [2]usize = @splat(0);
    var zr: [2]f64 = @splat(0.0);
    var zi: [2]f64 = @splat(0.0);

    const b_orig = br;
    const a_mat = CxMat.initPacked(&ar, &ai, 2, 2);
    const b_mat = CxMat.initPacked(&br, &bi, 2, 2);
    const a_for_check = [_]f64{ 2.0, 1.0, 0.0, 2.0 };
    const rc = mldivide(a_mat, b_mat, &ipvt, &zr, &zi);
    try zigstd.testExpectTrue(rc > 0.0);

    // Check A*X == original B.
    for (0..2) |row| {
        for (0..2) |col| {
            var s: f64 = 0.0;
            for (0..2) |k| s += a_for_check[ix(2, row + 1, k + 1)] * br[ix(2, k + 1, col + 1)];
            try zigstd.testExpectApproxEqAbs(b_orig[ix(2, row + 1, col + 1)], s, 1e-9);
        }
    }
}

test "mrdivide solves X*B=A" {
    // A = [[1,2],[3,4]], B = [[2,0],[1,2]] (both square, real).
    var ar = [_]f64{ 1.0, 3.0, 2.0, 4.0 }; // column-major
    var ai: [4]f64 = @splat(0.0);
    var br = [_]f64{ 2.0, 1.0, 0.0, 2.0 };
    var bi: [4]f64 = @splat(0.0);
    var ipvt: [2]usize = @splat(0);
    var zr: [2]f64 = @splat(0.0);
    var zi: [2]f64 = @splat(0.0);
    var rowr: [2]f64 = @splat(0.0);
    var rowi: [2]f64 = @splat(0.0);

    const b_for_check = br;
    const a_mat = CxMat.initPacked(&ar, &ai, 2, 2);
    const b_mat = CxMat.initPacked(&br, &bi, 2, 2);
    const rc = mrdivide(a_mat, b_mat, &ipvt, &zr, &zi, &rowr, &rowi);
    try zigstd.testExpectTrue(rc > 0.0);

    // Check X*B == original A.
    const a_orig = [_]f64{ 1.0, 3.0, 2.0, 4.0 };
    for (0..2) |row| {
        for (0..2) |col| {
            var s: f64 = 0.0;
            for (0..2) |k| s += ar[ix(2, row + 1, k + 1)] * b_for_check[ix(2, k + 1, col + 1)];
            try zigstd.testExpectApproxEqAbs(a_orig[ix(2, row + 1, col + 1)], s, 1e-9);
        }
    }
}
