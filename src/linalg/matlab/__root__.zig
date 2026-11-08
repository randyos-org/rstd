//! `matlab88.f90`'s function library (`ref/`) minus the REPL: `CxMat` and its
//! methods, and the commands below that are not matrix-shaped (`readme.md`).

const scalar = @import("scalar.zig");
const elementary = @import("elementary.zig");
const generators = @import("generators.zig");
const rat_mod = @import("rat.zig");
const base_mod = @import("base.zig");
const random = @import("random.zig");
const matrixfun = @import("matrixfun.zig");
const norm = @import("norm.zig");

pub const CxMat = @import("CxMat.zig");
pub const Cx = scalar.Cx;

/// HILB: REPL name for `hilber` (the closed-form inverse Hilbert matrix).
pub const hilb = generators.hilber;

/// MAGIC: REPL name for `generators.magic`.
pub const magic = generators.magic;

/// BASE: REPL name for `base.base`.
pub const base = base_mod.base;

/// RAT: REPL name for `rat.rat`.
pub const rat = rat_mod.rat;

/// POLY/ROOTS and the elementwise elementary functions: see
/// `matrixfun.zig`.
pub const poly = matrixfun.poly;
pub const roots = matrixfun.roots;
pub const sin = matrixfun.sin;
pub const cos = matrixfun.cos;
pub const exp = matrixfun.exp;
pub const sqrt = matrixfun.sqrt;
pub const log = matrixfun.log;
pub const atan = matrixfun.atan;

/// NORM (vector form): see `norm.zig`. The matrix-norm forms
/// (`matNorm2`/`matNorm1`/`matNormInf`/`matNormFro`) are `CxMat` methods.
pub const vecNorm2 = norm.vecNorm2;
pub const vecNormP = norm.vecNormP;
pub const vecNormInf = norm.vecNormInf;

/// POWER (`X**Y`): REPL name for `elementary.wpow`.
pub const pow = elementary.wpow;

/// RAND: REPL name for `random.Urand`; seed with `.init(seed)`, then
/// `.next()`/`.nextNormal()`/`.fillUniform()`/`.fillNormal()`.
pub const Rand = random.Urand;

/// SUM: sum of all entries of a complex vector/matrix (shape doesn't
/// matter; `ar`/`ai` are read as a flat array).
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
/// matter; `ar`/`ai` are read as a flat array).
pub fn prod(ar: []const f64, ai: []const f64) Cx {
    var acc = Cx.init(1.0, 0.0);
    for (ar, ai) |vr, vi| acc = acc.mul(Cx.init(vr, vi));
    return acc;
}

const zigstd = @import("../../zigstd.zig");
const lu = @import("lu.zig");
const cholesky = @import("cholesky.zig");
const eig = @import("eig.zig");
const svd = @import("svd.zig");
const divide = @import("divide.zig");
const schur = @import("schur.zig");
const matmul = @import("matmul.zig");
const elementwise = @import("elementwise.zig");

test {
    inline for (.{ CxMat, scalar, elementary, lu, divide, cholesky, generators, rat_mod, base_mod, random, eig, schur, svd, matrixfun, matmul, elementwise, norm }) |m| zigstd.testRefAllDecls(m);
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
