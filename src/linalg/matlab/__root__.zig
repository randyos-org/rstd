//! Zig port of `matlab88.f90`'s function library - the 1988 Argonne
//! National Laboratory "MATLAB" Fortran source (see `ref/matlab88_doc.txt`
//! / `ref/matlab88_help.txt` / `ref/readme.txt` for the original
//! documentation) - everything callable as a MATLAB function, minus the
//! REPL itself (parser, tokenizer, `<>` stack machine, `FOR`/`WHILE`/`IF`,
//! `SAVE`/`LOAD`/`DIARY`/`PLOT` I/O, `HELP` text) since those have no
//! meaning outside an interactive session. See `README.md` for a usage
//! guide and a full function table.
//!
//! The public surface is `CxMat` (bundles a complex matrix and carries
//! every self-shaped REPL command as a method - `mat.det(...)`,
//! `mat.eig(allocator, true)`, and so on) plus the handful of REPL
//! commands below that were never matrix-self-shaped to begin with
//! (`hilb`/`magic`, `sum`/`prod`, `pow`, the elementwise `sin`/`cos`/...
//! family, `poly`/`roots`, vector `norm`, `base`/`rat`/`Rand`). Everything
//! else (`lu`, `cholesky`, `svdapi`, `blas1`, `generators`, ...) is an internal
//! implementation detail behind those two surfaces, not part of the public
//! API - see `CxMat.zig`'s "Domain-file methods" section for how its
//! methods delegate to them.

const blas1 = @import("blas1.zig");
const elementary = @import("elementary.zig");
const generators = @import("generators.zig");
const rat_mod = @import("rat.zig");
const base_mod = @import("base.zig");
const random = @import("random.zig");
const matrixfun = @import("matrixfun.zig");
const norm = @import("norm.zig");

pub const CxMat = @import("CxMat.zig");
pub const Cx = blas1.Cx;

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
/// (`matNorm2`/`matNorm1`/`matNormInf`/`matNormFro`) are `CxMat` methods -
/// see `CxMat.zig`.
pub const vecNorm2 = norm.vecNorm2;
pub const vecNormP = norm.vecNormP;
pub const vecNormInf = norm.vecNormInf;

/// POWER (`X**Y`): REPL name for `elementary.wpow`.
pub const pow = elementary.wpow;

/// RAND: REPL name for `random.Urand` - seed with `.init(seed)`, then
/// `.next()`/`.nextNormal()`/`.fillUniform()`/`.fillNormal()`.
pub const Rand = random.Urand;

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

const zigstd = @import("../../zigstd.zig");
const lu = @import("lu.zig");
const cholesky = @import("cholesky.zig");
const eig = @import("eig.zig");
const divide = @import("divide.zig");
const svdapi = @import("svdapi.zig");
const matmul = @import("matmul.zig");
const schur = @import("schur.zig");
const elementwise = @import("elementwise.zig");

test {
    zigstd.testRefAllDecls(CxMat);
    zigstd.testRefAllDecls(blas1);
    zigstd.testRefAllDecls(elementary);
    zigstd.testRefAllDecls(lu);
    zigstd.testRefAllDecls(cholesky);
    zigstd.testRefAllDecls(generators);
    zigstd.testRefAllDecls(rat_mod);
    zigstd.testRefAllDecls(base_mod);
    zigstd.testRefAllDecls(random);
    zigstd.testRefAllDecls(eig);
    zigstd.testRefAllDecls(divide);
    zigstd.testRefAllDecls(svdapi);
    zigstd.testRefAllDecls(matrixfun);
    zigstd.testRefAllDecls(matmul);
    zigstd.testRefAllDecls(schur);
    zigstd.testRefAllDecls(norm);
    zigstd.testRefAllDecls(elementwise);
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
