//! Bundles a complex matrix's real/imaginary storage, leading dimension,
//! and shape into one value, instead of threading
//! `ar: []f64, ai: []f64, lda: usize, m: usize, n: usize` through every
//! matrix-level function call individually.
//!
//! Vector-level (BLAS-1) kernels in `blas1.zig` are unchanged by this -
//! they operate on a 1D strided view (a row, a column, or an arbitrary
//! vector with an explicit `incx`), not a whole matrix, matching real
//! BLAS's own Level-1-vs-Level-2/3 boundary. `CxMat` starts exactly where
//! "operates on a whole matrix" starts.
//!
//! Storage stays exactly what it always was: flat, column-major, `lda`
//! elements per column (may exceed `m`, enabling a submatrix view of a
//! larger allocated block without copying - see `idx`'s doc comment).
//! This struct changes nothing about layout or performance, only how many
//! parameters get threaded through call sites - every method here is a
//! thin wrapper that unpacks the same fields kernels already took as loose
//! arguments.

const std = @import("std");
const zigstd = @import("../../zigstd.zig");
const blas1 = @import("blas1.zig");
const Cx = blas1.Cx;

r: []f64,
i: []f64,
lda: usize,
m: usize,
n: usize,

const Self = @This();

/// Wrap existing storage. `r`/`i` must each have at least `lda*n`
/// elements (matching every kernel's own contract) - not copied.
pub fn init(r: []f64, i: []f64, lda: usize, m: usize, n: usize) Self {
    return .{ .r = r, .i = i, .lda = lda, .m = m, .n = n };
}

/// Wrap a tightly-packed `m`x`n` matrix (`lda == m`, no submatrix view) -
/// the common case.
pub fn initPacked(r: []f64, i: []f64, m: usize, n: usize) Self {
    return init(r, i, m, m, n);
}

/// 0-based flat offset of 1-based `(row, col)` - same convention as
/// `blas1.ix`, which this just forwards to.
pub inline fn idx(self: Self, row: usize, c: usize) usize {
    return blas1.ix(self.lda, row, c);
}

pub inline fn get(self: Self, row: usize, c: usize) Cx {
    const k = self.idx(row, c);
    return Cx.init(self.r[k], self.i[k]);
}

pub inline fn set(self: Self, row: usize, c: usize, v: Cx) void {
    const k = self.idx(row, c);
    self.r[k] = v.re;
    self.i[k] = v.im;
}

/// A contiguous length-`m` view of one column's real/imaginary parts -
/// unit stride, the fast case every BLAS1 kernel vectorizes. Returned by
/// `col`.
pub const ColView = struct { r: []f64, i: []f64 };

/// Column `c` (1-based) as a `ColView`.
pub inline fn col(self: Self, c: usize) ColView {
    const start = self.idx(1, c);
    return .{ .r = self.r[start..][0..self.m], .i = self.i[start..][0..self.m] };
}

/// A strided view of one row's real/imaginary parts, `stride` apart in
/// storage (not contiguous). Returned by `rowView`.
pub const RowView = struct { r: []f64, i: []f64, stride: isize };

/// Row `row` (1-based) as a `RowView`, strided by `lda` - matching every
/// existing row-wise access in this port (see e.g. `norm.matNormInf`).
/// Deliberately not materialized into a contiguous buffer: for a one-off,
/// non-reused row access this is consistently faster than gather-then-
/// vectorize (see `matmul.zig`'s doc comment for the benchmarked
/// reasoning) - only algorithms that actually reuse a gathered row across
/// many operations (like `matmul`) should pay that cost, not row access in
/// general.
pub inline fn rowView(self: Self, row: usize) RowView {
    const start = self.idx(row, 1);
    return .{ .r = self.r[start..], .i = self.i[start..], .stride = @intCast(self.lda) };
}

/// Copies `self` into freshly-allocated tightly-packed (`lda == m`)
/// storage. Caller frees the returned `CxMat`'s `.r`/`.i`.
pub fn clone(self: Self, allocator: zigstd.Allocator) !Self {
    const out_r = try allocator.alloc(f64, self.m * self.n);
    errdefer allocator.free(out_r);
    const out_i = try allocator.alloc(f64, self.m * self.n);
    errdefer allocator.free(out_i);
    const out = initPacked(out_r, out_i, self.m, self.n);
    var col_idx: usize = 1;
    while (col_idx <= self.n) : (col_idx += 1) {
        const src = self.col(col_idx);
        const dst = out.col(col_idx);
        @memcpy(dst.r, src.r);
        @memcpy(dst.i, src.i);
    }
    return out;
}

pub fn deinit(self: Self, allocator: zigstd.Allocator) void {
    allocator.free(self.r);
    allocator.free(self.i);
}

/// `X'`: conjugate transpose - REPL semantics for the `'` operator
/// (always conjugating, unlike plain `transpose`). Writes `out`
/// (`self.n`x`self.m`) with `out(j,i) = conj(self(i,j))`; `out` must not
/// alias `self`.
pub fn ctranspose(self: Self, out: Self) void {
    var j: usize = 1;
    while (j <= self.n) : (j += 1) {
        var i: usize = 1;
        while (i <= self.m) : (i += 1) {
            const v = self.get(i, j);
            out.set(j, i, Cx.init(v.re, -v.im));
        }
    }
}

/// Plain (non-conjugating) transpose: `out(j,i) = self(i,j)`. The REPL's
/// `'` operator is always the conjugate transpose (`ctranspose`); this is
/// for real-valued matrices or callers who explicitly don't want
/// conjugation.
pub fn transpose(self: Self, out: Self) void {
    var j: usize = 1;
    while (j <= self.n) : (j += 1) {
        var i: usize = 1;
        while (i <= self.m) : (i += 1) {
            out.set(j, i, self.get(i, j));
        }
    }
}

/// Scales every entry of `self` in place by the scalar `s` - REPL
/// semantics for `s*X`/`X*s` when `s` is a 1x1 matrix (e.g. `3*EYE`).
/// Loops per column rather than treating storage as one flat run, since
/// `self` may be a submatrix view (`self.lda > self.m`).
pub fn scale(self: Self, s: Cx) void {
    var j: usize = 1;
    while (j <= self.n) : (j += 1) {
        const c = self.col(j);
        blas1.wscal(self.m, s, c.r, c.i, 1);
    }
}

// --- Domain-file methods --------------------------------------------------
// Each alias below just gives dot-call syntax (`mat.wgefa(...)`) to a
// function that already takes `self: Self` as its first parameter but
// lives in a sibling file (kept there to stay under this project's
// line-count guideline, and because e.g. `eigen_general.zig`'s `comqr3` is
// a large, numerically-sensitive algorithm body not worth relocating).
// `a.foo(b, ...)`-shaped multi-matrix operations (`matmul`, `mrdivide`,
// `kron`) alias here too - Zig's method sugar only requires the first
// parameter be `Self`, not that every parameter be.

const cholesky = @import("cholesky.zig");
const lu = @import("lu.zig");
const hessenberg = @import("hessenberg.zig");
const tridiagonal = @import("tridiagonal.zig");
const eigen_general = @import("eigen_general.zig");
const generators = @import("generators.zig");
const divide = @import("divide.zig");
const matmul_mod = @import("matmul.zig");
const eig_mod = @import("eig.zig");
const schur_mod = @import("schur.zig");
const matrixfun = @import("matrixfun.zig");
const svd_mod = @import("svd.zig");
const qr_mod = @import("qr.zig");
const svdapi = @import("svdapi.zig");
const norm_mod = @import("norm.zig");
const elementwise = @import("elementwise.zig");

pub const wpofa = cholesky.wpofa;

pub const wgefa = lu.wgefa;
pub const wgesl = lu.wgesl;
pub const wgeco = lu.wgeco;
pub const wgedi = lu.wgedi;

pub const corth = hessenberg.corth;

pub const htridi = tridiagonal.htridi;
pub const htribk = tridiagonal.htribk;

pub const comqr3 = eigen_general.comqr3;

pub const rref = generators.rref;
pub const eye = generators.eye;
pub const ones = generators.ones;
pub const tril = generators.tril;
pub const triu = generators.triu;
pub const diagFromVector = generators.diagFromVector;
pub const diagToVector = generators.diagToVector;
pub const kron = generators.kron;

pub const wsvdc = svd_mod.wsvdc;

pub const wqrdc = qr_mod.wqrdc;
pub const wqrsl = qr_mod.wqrsl;

pub const cond = svdapi.cond;
pub const rank = svdapi.rank;
pub const pinv = svdapi.pinv;
pub const orth = svdapi.orth;

pub const matNorm2 = norm_mod.matNorm2;
pub const matNorm1 = norm_mod.matNorm1;
pub const matNormInf = norm_mod.matNormInf;
pub const matNormFro = norm_mod.matNormFro;

pub const mrdivide = divide.mrdivide;
pub const mldivide = divide.mldivide;

pub const matmul = matmul_mod.matmul;

pub const add = elementwise.add;
pub const sub = elementwise.sub;
pub const dotMul = elementwise.dotMul;
pub const dotDiv = elementwise.dotDiv;
pub const dotRdiv = elementwise.dotRdiv;

pub const eig = eig_mod.eig;
pub const EigResult = eig_mod.EigResult;
pub const EigError = eig_mod.EigError;

pub const hess = schur_mod.hess;
pub const HessResult = schur_mod.HessResult;
pub const schur = schur_mod.schur;
pub const SchurResult = schur_mod.SchurResult;

pub const expm = matrixfun.expm;
pub const logm = matrixfun.logm;
pub const sqrtm = matrixfun.sqrtm;
pub const sinm = matrixfun.sinm;
pub const cosm = matrixfun.cosm;
pub const atanm = matrixfun.atanm;
pub const powm = matrixfun.powm;

const lu_Determinant = lu.Determinant;

/// DET: determinant of `self`, via LU factorization. `self` is destroyed;
/// `ipvt` (length `self.n`) is scratch.
pub fn det(self: Self, ipvt: []usize) lu_Determinant {
    _ = lu.wgefa(self, ipvt);
    var result: lu_Determinant = undefined;
    var no_workspace: [0]f64 = .{};
    lu.wgedi(self, ipvt, &no_workspace, &no_workspace, 10, &result);
    return result;
}

/// INV: invert `self` in place, via LU factorization with condition
/// estimation. `ipvt` (length `self.n`) and `zr`/`zi` (length `self.n`)
/// are scratch. Returns the reciprocal condition number (`wgeco`'s
/// `rcond`) so callers can judge accuracy - a value so small that
/// `1.0 + rcond == 1.0` means `A` is numerically singular and the
/// "inverse" is not to be trusted.
pub fn inv(self: Self, ipvt: []usize, zr: []f64, zi: []f64) f64 {
    const rc = lu.wgeco(self, ipvt, zr, zi);
    // wgeco's own use of zr/zi (an approximate null vector) is done by the
    // time we get here, so it doubles as wgedi's scratch workspace.
    lu.wgedi(self, ipvt, zr, zi, 1, null);
    return rc;
}

/// CHOL: Cholesky-factor the Hermitian positive-definite matrix `self` in
/// place into upper-triangular `R` (`ctrans(R)*R == A`), zeroing the
/// strict lower triangle so the result is a clean triangular matrix
/// (unlike raw `wpofa`, which leaves it untouched). Returns 0 on success,
/// or the 1-based order of the leading minor that failed to be positive
/// definite (the matrix is left only partially factored in that case).
pub fn chol(self: Self) usize {
    const info = cholesky.wpofa(self);
    if (info != 0) return info;
    var j: usize = 1;
    while (j <= self.n) : (j += 1) {
        blas1.wset(self.n - j, Cx.init(0.0, 0.0), self.r[self.idx(j + 1, j)..], self.i[self.idx(j + 1, j)..], 1);
    }
    return 0;
}

test "ctranspose conjugates and transposes a non-square matrix" {
    // A = [[1+i, 2, 3-i]] (1x3) -> A' (3x1) = [[1-i],[2],[3+i]].
    var ar = [_]f64{ 1.0, 2.0, 3.0 };
    var ai = [_]f64{ 1.0, 0.0, -1.0 };
    var out_r: [3]f64 = undefined;
    var out_i: [3]f64 = undefined;
    const a = Self.initPacked(&ar, &ai, 1, 3);
    const out = Self.initPacked(&out_r, &out_i, 3, 1);

    a.ctranspose(out);
    try zigstd.testExpectEqualSlices(f64, &[_]f64{ 1.0, 2.0, 3.0 }, &out_r);
    try zigstd.testExpectEqualSlices(f64, &[_]f64{ -1.0, 0.0, 1.0 }, &out_i);
}

test "transpose does not conjugate" {
    var ar = [_]f64{ 1.0, 2.0 };
    var ai = [_]f64{ 1.0, -1.0 };
    var out_r: [2]f64 = undefined;
    var out_i: [2]f64 = undefined;
    const a = Self.initPacked(&ar, &ai, 1, 2);
    const out = Self.initPacked(&out_r, &out_i, 2, 1);

    a.transpose(out);
    try zigstd.testExpectEqualSlices(f64, &[_]f64{ 1.0, 2.0 }, &out_r);
    try zigstd.testExpectEqualSlices(f64, &[_]f64{ 1.0, -1.0 }, &out_i);
}

test "scale multiplies every entry by a complex scalar" {
    var ar = [_]f64{ 1.0, 2.0, 3.0, 4.0 };
    var ai: [4]f64 = @splat(0.0);
    const a = Self.initPacked(&ar, &ai, 2, 2);

    a.scale(Cx.init(0.0, 1.0)); // multiply by i
    try zigstd.testExpectEqualSlices(f64, &[_]f64{ 0.0, 0.0, 0.0, 0.0 }, &ar);
    try zigstd.testExpectEqualSlices(f64, &[_]f64{ 1.0, 2.0, 3.0, 4.0 }, &ai);
}

test "det, inv and chol agree on a small real matrix" {
    var ar = [_]f64{ 4.0, 2.0, 2.0, 3.0 }; // [[4,2],[2,3]], det=8, SPD
    var ai: [4]f64 = @splat(0.0);
    var ipvt: [2]usize = @splat(0);

    var ar_for_det = ar;
    var ai_for_det = ai;
    const d = Self.initPacked(&ar_for_det, &ai_for_det, 2, 2).det(&ipvt);
    try zigstd.testExpectApproxEqAbs(8.0, d.mantissa.re * std.math.pow(f64, 10.0, d.exponent10), 1e-9);

    var ar_for_inv = ar;
    var ai_for_inv = ai;
    var zr: [2]f64 = @splat(0.0);
    var zi: [2]f64 = @splat(0.0);
    const rc = Self.initPacked(&ar_for_inv, &ai_for_inv, 2, 2).inv(&ipvt, &zr, &zi);
    try zigstd.testExpectTrue(rc > 0.0);
    // inv([[4,2],[2,3]]) = [[3,-2],[-2,4]] / 8.
    try zigstd.testExpectApproxEqAbs(3.0 / 8.0, ar_for_inv[blas1.ix(2, 1, 1)], 1e-9);
    try zigstd.testExpectApproxEqAbs(-2.0 / 8.0, ar_for_inv[blas1.ix(2, 1, 2)], 1e-9);
    try zigstd.testExpectApproxEqAbs(4.0 / 8.0, ar_for_inv[blas1.ix(2, 2, 2)], 1e-9);

    const info = Self.initPacked(&ar, &ai, 2, 2).chol();
    try zigstd.testExpectEqual(@as(usize, 0), info);
    try zigstd.testExpectApproxEqAbs(0.0, ar[blas1.ix(2, 2, 1)], 1e-12); // lower triangle zeroed
}

test "init/get/set/idx agree with blas1.ix" {
    var r = [_]f64{ 0, 0, 0, 0 };
    var i = [_]f64{ 0, 0, 0, 0 };
    const mat = Self.initPacked(&r, &i, 2, 2);
    try zigstd.testExpectEqual(blas1.ix(2, 2, 1), mat.idx(2, 1));

    mat.set(2, 1, Cx.init(3.0, -4.0));
    try zigstd.testExpectEqual(@as(f64, 3.0), r[mat.idx(2, 1)]);
    const got = mat.get(2, 1);
    try zigstd.testExpectEqual(@as(f64, 3.0), got.re);
    try zigstd.testExpectEqual(@as(f64, -4.0), got.im);
}

test "col returns a contiguous unit-stride view" {
    // Column-major 2x2: col1 = [1,2], col2 = [3,4].
    var r = [_]f64{ 1, 2, 3, 4 };
    var i: [4]f64 = @splat(0.0);
    const mat = Self.initPacked(&r, &i, 2, 2);
    const c2 = mat.col(2);
    try zigstd.testExpectEqualSlices(f64, &[_]f64{ 3, 4 }, c2.r);
}

test "rowView strides by lda, not contiguously" {
    // Column-major 2x2: row1 = [1,3] (stride 2 apart in storage).
    var r = [_]f64{ 1, 2, 3, 4 };
    var i: [4]f64 = @splat(0.0);
    const mat = Self.initPacked(&r, &i, 2, 2);
    const row1 = mat.rowView(1);
    try zigstd.testExpectEqual(@as(isize, 2), row1.stride);
    try zigstd.testExpectEqual(@as(f64, 1), row1.r[0]);
    try zigstd.testExpectEqual(@as(f64, 3), row1.r[@intCast(row1.stride)]);
}

test "clone copies into independent tightly-packed storage" {
    var r = [_]f64{ 1, 2, 3, 4 };
    var i: [4]f64 = @splat(0.0);
    const mat = Self.initPacked(&r, &i, 2, 2);
    var copy = try mat.clone(zigstd.testAllocator);
    defer copy.deinit(zigstd.testAllocator);

    try zigstd.testExpectEqualSlices(f64, mat.r, copy.r);
    copy.r[0] = 99.0;
    try zigstd.testExpectEqual(@as(f64, 1.0), mat.r[0]); // independent storage
}
