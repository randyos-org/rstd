//! Friendly, MATLAB-command-named bindings on top of the low-level
//! Fortran-derived kernels in this directory - e.g. `eig` here is the
//! dispatch logic the original REPL's `EIG` command used to pick between
//! the Hermitian (`htridi`+`imtql2`+`htribk`) and general
//! (`corth`+`comqr3`) eigensolvers, not a new algorithm. Pure aliases
//! (`lu`, `svd`, `qr`, `hilb`) just repoint a REPL name at the underlying
//! kernel; see that kernel's doc comment for the actual contract.
//!
//! Functions here that need scratch space sized at runtime (`eig`,
//! `cond`, `rank`, `pinv`, `orth`, `roots`, and the `*m` matrix-function
//! family) take a `std.mem.Allocator` and return caller-owned slices;
//! everything else works directly on caller-supplied buffers like the
//! low-level kernels do, so this module stays usable without an allocator
//! where the original REPL command didn't need one either.
//!
//! This file is itself just the aliases plus the handful of LU-based
//! operations (`det`/`inv`/`chol`/`rcond`/`sum`/`prod`); the rest live in
//! sibling files kept under this project's line-count guideline and
//! re-exported here so every REPL name still resolves as `api.<name>`:
//! `eig.zig` (`eig`), `divide.zig` (`mrdivide`), `svdapi.zig`
//! (`cond`/`rank`/`pinv`/`orth`), `matrixfun.zig` (`poly`/`roots`, the
//! elementwise `sin`/`cos`/... family, and the matrix-valued
//! `expm`/`logm`/... family), `matmul.zig` (`matmul`, the `*` operator's
//! math with the REPL stack-machine bookkeeping left out), `schur.zig`
//! (`hess`/`schur`), `norm.zig` (`vecNorm2`/`vecNormP`/`vecNormInf`/
//! `matNorm2`/`matNorm1`/`matNormInf`/`matNormFro`), `elementary.zig`
//! (`wpow`, re-exported as `pow`), and `random.zig` (`Urand`, re-exported
//! as `Rand`).

const std = @import("std");
const zigstd = @import("../../zigstd.zig");
const blas1 = @import("blas1.zig");
const lu_mod = @import("lu.zig");
const cholesky = @import("cholesky.zig");
const generators = @import("generators.zig");
const svd_mod = @import("svd.zig");
const qr_mod = @import("qr.zig");
const eig_mod = @import("eig.zig");
const divide_mod = @import("divide.zig");
const svdapi = @import("svdapi.zig");
const matrixfun = @import("matrixfun.zig");
const matmul_mod = @import("matmul.zig");
const schur_mod = @import("schur.zig");
const norm_mod = @import("norm.zig");
const elementary = @import("elementary.zig");
const random = @import("random.zig");

const Cx = blas1.Cx;
const ix = blas1.ix;

/// LU: factor `A` in place (packed LINPACK/LAPACK-style: unit multipliers
/// below the diagonal, `U` on and above it) - REPL name for `wgefa`.
pub const lu = lu_mod.wgefa;

/// HILB: REPL name for `hilber` (the closed-form inverse Hilbert matrix).
pub const hilb = generators.hilber;

/// MAGIC: REPL name for `generators.magic`.
pub const magic = generators.magic;

/// EYE: REPL name for `generators.eye`.
pub const eye = generators.eye;

/// ONES: REPL name for `generators.ones`.
pub const ones = generators.ones;

/// TRIL: REPL name for `generators.tril`.
pub const tril = generators.tril;

/// TRIU: REPL name for `generators.triu`.
pub const triu = generators.triu;

/// KRON: REPL name for `generators.kron`.
pub const kron = generators.kron;

/// DIAG (vector -> matrix form): REPL name for `generators.diagFromVector`.
pub const diagFromVector = generators.diagFromVector;

/// DIAG (matrix -> vector form): REPL name for `generators.diagToVector`.
pub const diagToVector = generators.diagToVector;

/// BASE: REPL name for `base.base`.
pub const base = @import("base.zig").base;

/// RAT: REPL name for `rat.rat`.
pub const rat = @import("rat.zig").rat;

/// SVD: REPL name for `wsvdc`.
pub const svd = svd_mod.wsvdc;

/// QR: REPL name for `wqrdc`. Pair with `qrsl` (== `wqrsl`) to extract `Q`,
/// solve least squares, etc. - see `qr.zig`'s tests for a worked example.
pub const qr = qr_mod.wqrdc;
pub const qrsl = qr_mod.wqrsl;

/// RREF: REPL name for `generators.rref`.
pub const rref = generators.rref;

/// EIG: see `eig.zig`.
pub const eig = eig_mod.eig;
pub const EigResult = eig_mod.EigResult;
pub const EigError = eig_mod.EigError;

/// Right matrix division (`A/B`): see `divide.zig`.
pub const mrdivide = divide_mod.mrdivide;

/// COND/RANK/PINV/ORTH: see `svdapi.zig`.
pub const cond = svdapi.cond;
pub const rank = svdapi.rank;
pub const pinv = svdapi.pinv;
pub const orth = svdapi.orth;

/// POLY/ROOTS and the elementwise/matrix-valued elementary functions: see
/// `matrixfun.zig`.
pub const poly = matrixfun.poly;
pub const roots = matrixfun.roots;
pub const sin = matrixfun.sin;
pub const cos = matrixfun.cos;
pub const exp = matrixfun.exp;
pub const sqrt = matrixfun.sqrt;
pub const log = matrixfun.log;
pub const atan = matrixfun.atan;
pub const expm = matrixfun.expm;
pub const logm = matrixfun.logm;
pub const sqrtm = matrixfun.sqrtm;
pub const sinm = matrixfun.sinm;
pub const cosm = matrixfun.cosm;
pub const atanm = matrixfun.atanm;

/// `X*Y`: matrix multiplication. See `matmul.zig`.
pub const matmul = matmul_mod.matmul;

/// HESS/SCHUR: see `schur.zig`.
pub const hess = schur_mod.hess;
pub const HessResult = schur_mod.HessResult;
pub const schur = schur_mod.schur;
pub const SchurResult = schur_mod.SchurResult;

/// NORM: see `norm.zig`. Split by vector-vs-matrix and norm-kind, matching
/// this file's existing `sin`/`sinm`-style shape-explicit naming rather
/// than one function that dispatches on argument shape.
pub const vecNorm2 = norm_mod.vecNorm2;
pub const vecNormP = norm_mod.vecNormP;
pub const vecNormInf = norm_mod.vecNormInf;
pub const matNorm2 = norm_mod.matNorm2;
pub const matNorm1 = norm_mod.matNorm1;
pub const matNormInf = norm_mod.matNormInf;
pub const matNormFro = norm_mod.matNormFro;

/// POWER (`X**Y`): REPL name for `elementary.wpow`.
pub const pow = elementary.wpow;

/// RAND: REPL name for `random.Urand` - seed with `.init(seed)`, then
/// `.next()`/`.nextNormal()`/`.fillUniform()`/`.fillNormal()`.
pub const Rand = random.Urand;

/// DET: determinant of the `n`x`n` matrix `ar`/`ai` (column-major, leading
/// dimension `lda`), via LU factorization. `A` is destroyed; `ipvt`
/// (length `n`) is scratch.
pub fn det(ar: []f64, ai: []f64, lda: usize, n: usize, ipvt: []usize) lu_mod.Determinant {
    _ = lu_mod.wgefa(ar, ai, lda, n, ipvt);
    var result: lu_mod.Determinant = undefined;
    var no_workspace: [0]f64 = .{};
    lu_mod.wgedi(ar, ai, lda, n, ipvt, &no_workspace, &no_workspace, 10, &result);
    return result;
}

/// INV: invert the `n`x`n` matrix `ar`/`ai` in place, via LU factorization
/// with condition estimation. `ipvt` (length `n`) and `zr`/`zi` (length
/// `n`) are scratch. Returns the reciprocal condition number (`wgeco`'s
/// `rcond`) so callers can judge accuracy - a value so small that
/// `1.0 + rcond == 1.0` means `A` is numerically singular and the
/// "inverse" is not to be trusted.
pub fn inv(ar: []f64, ai: []f64, lda: usize, n: usize, ipvt: []usize, zr: []f64, zi: []f64) f64 {
    const rc = lu_mod.wgeco(ar, ai, lda, n, ipvt, zr, zi);
    // wgeco's own use of zr/zi (an approximate null vector) is done by the
    // time we get here, so it doubles as wgedi's scratch workspace.
    lu_mod.wgedi(ar, ai, lda, n, ipvt, zr, zi, 1, null);
    return rc;
}

/// CHOL: Cholesky-factor the Hermitian positive-definite matrix `ar`/`ai`
/// in place into upper-triangular `R` (`ctrans(R)*R == A`), zeroing the
/// strict lower triangle so the result is a clean triangular matrix
/// (unlike raw `wpofa`, which leaves it untouched). Returns 0 on success,
/// or the 1-based order of the leading minor that failed to be positive
/// definite (the matrix is left only partially factored in that case).
pub fn chol(ar: []f64, ai: []f64, lda: usize, n: usize) usize {
    const info = cholesky.wpofa(ar, ai, lda, n);
    if (info != 0) return info;
    var j: usize = 1;
    while (j <= n) : (j += 1) {
        blas1.wset(n - j, Cx.init(0.0, 0.0), ar[ix(lda, j + 1, j)..], ai[ix(lda, j + 1, j)..], 1);
    }
    return 0;
}

/// RCOND: REPL name for `lu.wgeco`'s reciprocal condition number estimate
/// (LU-based, distinct from `cond`'s SVD-based 2-norm ratio - the original
/// exposed both since they trade accuracy for speed differently).
pub const rcond = lu_mod.wgeco;

/// SUM: sum of all entries of a complex vector/matrix (shape doesn't
/// matter - `ar`/`ai` are read as a flat array).
pub fn sum(ar: []const f64, ai: []const f64) Cx {
    var sr: f64 = 0.0;
    var si: f64 = 0.0;
    for (ar, ai) |vr, vi| {
        sr += vr;
        si += vi;
    }
    return Cx.init(sr, si);
}

/// PROD: product of all entries of a complex vector/matrix (shape doesn't
/// matter - `ar`/`ai` are read as a flat array).
pub fn prod(ar: []const f64, ai: []const f64) Cx {
    var acc = Cx.init(1.0, 0.0);
    for (ar, ai) |vr, vi| acc = blas1.wmul(acc, Cx.init(vr, vi));
    return acc;
}


test "det, inv and chol agree on a small real matrix" {
    var ar = [_]f64{ 4.0, 2.0, 2.0, 3.0 }; // [[4,2],[2,3]], det=8, SPD
    var ai: [4]f64 = @splat(0.0);
    var ipvt: [2]usize = @splat(0);

    var ar_for_det = ar;
    var ai_for_det = ai;
    const d = det(&ar_for_det, &ai_for_det, 2, 2, &ipvt);
    try zigstd.testExpectApproxEqAbs(8.0, d.mantissa.re * std.math.pow(f64, 10.0, d.exponent10), 1e-9);

    var ar_for_inv = ar;
    var ai_for_inv = ai;
    var zr: [2]f64 = @splat(0.0);
    var zi: [2]f64 = @splat(0.0);
    const rc = inv(&ar_for_inv, &ai_for_inv, 2, 2, &ipvt, &zr, &zi);
    try zigstd.testExpectTrue(rc > 0.0);
    // inv([[4,2],[2,3]]) = [[3,-2],[-2,4]] / 8.
    try zigstd.testExpectApproxEqAbs(3.0 / 8.0, ar_for_inv[ix(2, 1, 1)], 1e-9);
    try zigstd.testExpectApproxEqAbs(-2.0 / 8.0, ar_for_inv[ix(2, 1, 2)], 1e-9);
    try zigstd.testExpectApproxEqAbs(4.0 / 8.0, ar_for_inv[ix(2, 2, 2)], 1e-9);

    const info = chol(&ar, &ai, 2, 2);
    try zigstd.testExpectEqual(@as(usize, 0), info);
    try zigstd.testExpectApproxEqAbs(0.0, ar[ix(2, 2, 1)], 1e-12); // lower triangle zeroed
}

test "sum and prod match hand computation" {
    const ar = [_]f64{ 1.0, 2.0, 3.0 };
    const ai = [_]f64{ 0.0, 1.0, 0.0 };
    const s = sum(&ar, &ai);
    try zigstd.testExpectApproxEqAbs(6.0, s.re, 1e-12);
    try zigstd.testExpectApproxEqAbs(1.0, s.im, 1e-12);

    const p = prod(&ar, &ai);
    // 1 * (2+i) * 3 = 6+3i
    try zigstd.testExpectApproxEqAbs(6.0, p.re, 1e-12);
    try zigstd.testExpectApproxEqAbs(3.0, p.im, 1e-12);
}
