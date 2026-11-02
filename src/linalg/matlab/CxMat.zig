//! Bundles a complex matrix's real/imaginary storage, leading dimension,
//! and shape into one value, instead of threading
//! `ar: []f64, ai: []f64, lda: usize, m: usize, n: usize` through every
//! matrix-level function call individually.
//!
//! Vector-level (BLAS-1) kernels, rstd's BLAS since M2, are unchanged by this -
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

const zigstd = @import("../../zigstd.zig");
const scalar = @import("scalar.zig");
const blas = @import("../blas/__root__.zig");
const Cx = scalar.Cx;

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
/// `scalar.ix`, which this just forwards to.
pub inline fn idx(self: Self, row: usize, c: usize) usize {
    return scalar.ix(self.lda, row, c);
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
/// existing row-wise access in this port (see e.g. `norm.matNormInf`),
/// not gathered into a contiguous buffer.
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

/// The interleaved element of the LAPACK boundary (M2).
pub const Complex = Cx;

/// Whether every imaginary part is zero: the real LAPACK routines then run
/// in place on `r`.
pub fn isReal(self: Self) bool {
    for (1..self.n + 1) |j| {
        for (self.col(j).i) |x| if (x != 0) return false;
    }
    return true;
}

/// `self` as a packed (`m` leading dimension) copy for LAPACK: interleaved
/// for `Complex`, `r` alone for `f64`. Caller frees.
pub fn toPacked(self: Self, comptime T: type, allocator: zigstd.Allocator) ![]T {
    const z = try allocator.alloc(T, self.m * self.n);
    for (1..self.n + 1) |j| {
        const c = self.col(j);
        const d = z[(j - 1) * self.m ..][0..self.m];
        if (T == f64) {
            @memcpy(d, c.r);
        } else for (d, c.r, c.i) |*e, re, im| e.* = Cx.init(re, im);
    }
    return z;
}

/// Writes `toPacked`'s layout `z` back into `self`; `f64` zeros `i`.
pub fn fromPacked(self: Self, comptime T: type, z: []const T) void {
    for (1..self.n + 1) |j| {
        const c = self.col(j);
        const s = z[(j - 1) * self.m ..][0..self.m];
        if (T == f64) {
            @memcpy(c.r, s);
            @memset(c.i, 0);
        } else for (s, c.r, c.i) |e, *re, *im| {
            re.* = e.re;
            im.* = e.im;
        }
    }
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
/// Per column, `drot` of its real and imaginary parts by `(re s, -im s)`:
/// `zscal`'s arithmetic over the split layout.
pub fn scale(self: Self, io: std.Io, s: Cx) std.Io.Cancelable!void {
    for (1..self.n + 1) |j| {
        const c = self.col(j);
        try scalar.level1(blas.rot(io, f64, @intCast(self.m), c.r.ptr, 1, c.i.ptr, 1, s.re, -s.im));
    }
}

// --- Domain-file methods --------------------------------------------------
// Each alias below gives dot-call syntax (`mat.det(...)`) to a function
// that takes `self: Self` first but lives in a sibling file (kept there to
// stay under this project's line-count guideline); those taking an
// allocator first (`eig`, `cond`, ...) are reached as `CxMat.eig(...)`.
// `a.foo(b, ...)`-shaped multi-matrix operations (`matmul`, `mrdivide`,
// `kron`) alias here too - Zig's method sugar only requires the first
// parameter be `Self`, not that every parameter be.

const cholesky = @import("cholesky.zig");
const lu = @import("lu.zig");
const generators = @import("generators.zig");
const divide = @import("divide.zig");
const matmul_mod = @import("matmul.zig");
const eig_mod = @import("eig.zig");
const schur_mod = @import("schur.zig");
const matrixfun = @import("matrixfun.zig");
const svdapi = @import("svdapi.zig");
const norm_mod = @import("norm.zig");
const elementwise = @import("elementwise.zig");

pub const chol = cholesky.chol;

pub const det = lu.det;
pub const inv = lu.inv;
pub const Determinant = lu.Determinant;

pub const rref = generators.rref;
pub const eye = generators.eye;
pub const ones = generators.ones;
pub const tril = generators.tril;
pub const triu = generators.triu;
pub const diagFromVector = generators.diagFromVector;
pub const diagToVector = generators.diagToVector;
pub const kron = generators.kron;

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

const std = @import("std");

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

    try a.scale(std.testing.io, Cx.init(0.0, 1.0)); // multiply by i
    try zigstd.testExpectEqualSlices(f64, &[_]f64{ 0.0, 0.0, 0.0, 0.0 }, &ar);
    try zigstd.testExpectEqualSlices(f64, &[_]f64{ 1.0, 2.0, 3.0, 4.0 }, &ai);
}

test "det, inv and chol agree on a small real matrix" {
    var ar = [_]f64{ 4.0, 2.0, 2.0, 3.0 }; // [[4,2],[2,3]], det=8, SPD
    var ai: [4]f64 = @splat(0.0);

    var ar_for_det = ar;
    var ai_for_det = ai;
    const d = try Self.initPacked(&ar_for_det, &ai_for_det, 2, 2).det(zigstd.testAllocator, std.testing.io);
    try zigstd.testExpectApproxEqAbs(8.0, d.mantissa.re * std.math.pow(f64, 10.0, d.exponent10), 1e-9);

    var ar_for_inv = ar;
    var ai_for_inv = ai;
    const rc = try Self.initPacked(&ar_for_inv, &ai_for_inv, 2, 2).inv(zigstd.testAllocator, std.testing.io);
    try zigstd.testExpectTrue(rc > 0.0);
    // inv([[4,2],[2,3]]) = [[3,-2],[-2,4]] / 8.
    try zigstd.testExpectApproxEqAbs(3.0 / 8.0, ar_for_inv[scalar.ix(2, 1, 1)], 1e-9);
    try zigstd.testExpectApproxEqAbs(-2.0 / 8.0, ar_for_inv[scalar.ix(2, 1, 2)], 1e-9);
    try zigstd.testExpectApproxEqAbs(4.0 / 8.0, ar_for_inv[scalar.ix(2, 2, 2)], 1e-9);

    try Self.initPacked(&ar, &ai, 2, 2).chol(zigstd.testAllocator, std.testing.io, null);
    try zigstd.testExpectApproxEqAbs(0.0, ar[scalar.ix(2, 2, 1)], 1e-12); // lower triangle zeroed
}

test "init/get/set/idx agree with scalar.ix" {
    var r = [_]f64{ 0, 0, 0, 0 };
    var i = [_]f64{ 0, 0, 0, 0 };
    const mat = Self.initPacked(&r, &i, 2, 2);
    try zigstd.testExpectEqual(scalar.ix(2, 2, 1), mat.idx(2, 1));

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
