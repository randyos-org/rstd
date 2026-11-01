//! `POLY`/`ROOTS`, and elementwise vs. matrix-valued elementary functions,
//! kept as its own file (rather than folded into `__root__.zig`/`CxMat.zig`,
//! which re-export these) to keep files under this project's line-count
//! guideline.

const std = @import("std");
const zigstd = @import("../../zigstd.zig");
const blas1 = @import("blas1.zig");
const CxMat = @import("CxMat.zig");
const eig_mod = @import("eig.zig");
const divide = @import("divide.zig");
const elementary = @import("elementary.zig");

const Cx = blas1.Cx;
const ix = blas1.ix;

/// POLY: characteristic polynomial coefficients (length `n+1`, highest
/// degree first, leading coefficient `1`) for the given `n` roots -
/// `prod(x - roots[k])` expanded.
pub fn poly(rootsr: []const f64, rootsi: []const f64, n: usize, cr: []f64, ci: []f64) void {
    @memset(cr[0 .. n + 1], 0.0);
    @memset(ci[0 .. n + 1], 0.0);
    cr[0] = 1.0;
    for (0..n) |idx| {
        const neg_root = Cx.init(-rootsr[idx], -rootsi[idx]);
        var k: usize = idx + 1;
        while (true) {
            const c = Cx.init(cr[k - 1], ci[k - 1]);
            const term = blas1.wmul(neg_root, c);
            cr[k] += term.re;
            ci[k] += term.im;
            if (k == 1) break;
            k -= 1;
        }
    }
}

/// ROOTS: the `n` roots of the degree-`n` polynomial with coefficients
/// `coeffr`/`coeffi` (length `>= n+1`, highest degree first; leading
/// zero coefficients are skipped), via the eigenvalues of its companion
/// matrix by `eig` (`geev`, which balances it, as MATLAB's `roots`). Writes
/// into `wr`/`wi` and returns the degree found (`0` for the zero
/// polynomial).
pub fn roots(allocator: zigstd.Allocator, io: std.Io, coeffr: []const f64, coeffi: []const f64, wr: []f64, wi: []f64) eig_mod.EigError!usize {
    var start: usize = 0;
    while (start < coeffr.len and blas1.cabs1(coeffr[start], coeffi[start]) == 0.0) start += 1;
    if (start >= coeffr.len) return 0;
    const n = coeffr.len - start - 1;
    if (n == 0) return 0;
    const c0 = Cx.init(coeffr[start], coeffi[start]);

    const compr = try allocator.alloc(f64, n * n);
    defer allocator.free(compr);
    const compi = try allocator.alloc(f64, n * n);
    defer allocator.free(compi);
    @memset(compr, 0.0);
    @memset(compi, 0.0);
    var j: usize = 1;
    while (j <= n) : (j += 1) {
        if (j < n) compr[ix(n, j + 1, j)] = 1.0;
        const cj = Cx.init(coeffr[start + j], coeffi[start + j]);
        const v = blas1.wdiv(Cx.init(-cj.re, -cj.im), c0);
        compr[ix(n, 1, j)] = v.re;
        compi[ix(n, 1, j)] = v.im;
    }

    const result = try eig_mod.eig(allocator, io, CxMat.initPacked(compr, compi, n, n), false);
    defer result.deinit(allocator);
    @memcpy(wr[0..n], result.values_r);
    @memcpy(wi[0..n], result.values_i);
    return n;
}

fn applyElementwise(comptime f: anytype, ar: []f64, ai: []f64) void {
    for (ar, ai) |*vr, *vi| {
        const r = f(Cx.init(vr.*, vi.*));
        vr.* = r.re;
        vi.* = r.im;
    }
}

/// SIN/COS/EXP/SQRT/LOG/ATAN (elementwise): apply the named elementary
/// function to every entry of `ar`/`ai` in place (any shape - vector or
/// matrix, read as a flat array). For the *matrix* function instead (e.g.
/// the matrix exponential, not entrywise `exp`), see `expm` and friends
/// below - the original conflated the two based on argument shape
/// (`sin([1,2,3])` elementwise vs. `sin(square_matrix)` matrix-valued);
/// this port keeps them as separate, explicitly named functions instead.
pub fn sin(ar: []f64, ai: []f64) void {
    applyElementwise(std.math.complex.sin, ar, ai);
}
pub fn cos(ar: []f64, ai: []f64) void {
    applyElementwise(std.math.complex.cos, ar, ai);
}
pub fn exp(ar: []f64, ai: []f64) void {
    applyElementwise(std.math.complex.exp, ar, ai);
}
pub fn sqrt(ar: []f64, ai: []f64) void {
    applyElementwise(std.math.complex.sqrt, ar, ai);
}
pub fn log(ar: []f64, ai: []f64) void {
    applyElementwise(blas1.wlog, ar, ai);
}
pub fn atan(ar: []f64, ai: []f64) void {
    applyElementwise(blas1.watan, ar, ai);
}

/// ABS/ROUND/REAL/IMAG/CONJ (elementwise): the remaining `FUN`-style
/// commands, each trivial enough not to need a loop helper - apply
/// directly with a `for` where you need them:
/// `for (ar, ai) |r, i| out[k] = blas1.pythag(r, i);` (ABS),
/// `for (ar) |*r| r.* = @round(r.*);` (ROUND, real part only - the
/// original ignored the imaginary part), `ar`/`ai` themselves are
/// REAL/IMAG, and `Cx.init(r, -i)` is CONJG.
/// Applies `f.call(lambda)` to each eigenvalue of the `n`x`n` matrix `a`
/// (not mutated) and reconstructs `f(A) = V * diag(f(lambda)) * inv(V)` -
/// the matrix function, not the entrywise one (see `sin`/`cos`/... above
/// for that). Only correct when `A` is diagonalizable; a defective matrix
/// (repeated eigenvalue without a full set of eigenvectors) makes `V`
/// singular and this will fail or blow up, matching the original, which
/// used the same eigendecomposition approach. `f` is a small value (often
/// zero-sized) with a `call(self, Cx) Cx` method rather than a bare
/// function, so callers needing runtime parameters (e.g. `powm`'s
/// exponent) can close over them without duplicating this function.
fn applyMatrixFunction(allocator: zigstd.Allocator, io: std.Io, f: anytype, a: CxMat, out: CxMat) !void {
    const n = a.n;
    const result = try eig_mod.eig(allocator, io, a, true);
    defer result.deinit(allocator);
    const vr = result.vectors_r.?;
    const vi = result.vectors_i.?;

    var j: usize = 1;
    while (j <= n) : (j += 1) {
        const s = f.call(Cx.init(result.values_r[j - 1], result.values_i[j - 1]));
        var i: usize = 1;
        while (i <= n) : (i += 1) {
            const t = blas1.wmul(Cx.init(vr[ix(n, i, j)], vi[ix(n, i, j)]), s);
            out.r[ix(n, i, j)] = t.re;
            out.i[ix(n, i, j)] = t.im;
        }
    }

    const v_copy_r = try allocator.alloc(f64, n * n);
    defer allocator.free(v_copy_r);
    const v_copy_i = try allocator.alloc(f64, n * n);
    defer allocator.free(v_copy_i);
    @memcpy(v_copy_r, vr);
    @memcpy(v_copy_i, vi);
    const v_copy_mat = CxMat.initPacked(v_copy_r, v_copy_i, n, n);
    _ = try divide.mrdivide(out, v_copy_mat, allocator, io);
}

const ExpCtx = struct {
    pub fn call(_: @This(), z: Cx) Cx {
        return std.math.complex.exp(z);
    }
};
const LogCtx = struct {
    pub fn call(_: @This(), z: Cx) Cx {
        return blas1.wlog(z);
    }
};
const SqrtCtx = struct {
    pub fn call(_: @This(), z: Cx) Cx {
        return std.math.complex.sqrt(z);
    }
};
const SinCtx = struct {
    pub fn call(_: @This(), z: Cx) Cx {
        return std.math.complex.sin(z);
    }
};
const CosCtx = struct {
    pub fn call(_: @This(), z: Cx) Cx {
        return std.math.complex.cos(z);
    }
};
const AtanCtx = struct {
    pub fn call(_: @This(), z: Cx) Cx {
        return blas1.watan(z);
    }
};
const PowCtx = struct {
    p: Cx,
    pub fn call(self: @This(), z: Cx) Cx {
        return elementary.wpow(z, self.p);
    }
};

pub fn expm(allocator: zigstd.Allocator, io: std.Io, a: CxMat, out: CxMat) !void {
    return applyMatrixFunction(allocator, io, ExpCtx{}, a, out);
}
pub fn logm(allocator: zigstd.Allocator, io: std.Io, a: CxMat, out: CxMat) !void {
    return applyMatrixFunction(allocator, io, LogCtx{}, a, out);
}
pub fn sqrtm(allocator: zigstd.Allocator, io: std.Io, a: CxMat, out: CxMat) !void {
    return applyMatrixFunction(allocator, io, SqrtCtx{}, a, out);
}
pub fn sinm(allocator: zigstd.Allocator, io: std.Io, a: CxMat, out: CxMat) !void {
    return applyMatrixFunction(allocator, io, SinCtx{}, a, out);
}
pub fn cosm(allocator: zigstd.Allocator, io: std.Io, a: CxMat, out: CxMat) !void {
    return applyMatrixFunction(allocator, io, CosCtx{}, a, out);
}
pub fn atanm(allocator: zigstd.Allocator, io: std.Io, a: CxMat, out: CxMat) !void {
    return applyMatrixFunction(allocator, io, AtanCtx{}, a, out);
}

/// POWER (`X**p`, matrix case): `a` raised to the complex power `p`, via
/// the same eigendecomposition approach as `expm`/`logm`/... above. The
/// original special-cased positive integer `p` as repeated multiplication
/// instead; this port always goes through the eigendecomposition, so a
/// caller wanting that accuracy for a known small integer power can chain
/// `matmul` calls directly.
pub fn powm(allocator: zigstd.Allocator, io: std.Io, a: CxMat, p: Cx, out: CxMat) !void {
    return applyMatrixFunction(allocator, io, PowCtx{ .p = p }, a, out);
}

test "poly and roots are inverse to each other for known roots" {
    // roots 1, 2, 3 -> (x-1)(x-2)(x-3) = x^3 - 6x^2 + 11x - 6.
    const rootsr = [_]f64{ 1.0, 2.0, 3.0 };
    const rootsi: [3]f64 = @splat(0.0);
    var cr: [4]f64 = @splat(0.0);
    var ci: [4]f64 = @splat(0.0);
    poly(&rootsr, &rootsi, 3, &cr, &ci);
    try zigstd.testExpectApproxEqAbs(1.0, cr[0], 1e-9);
    try zigstd.testExpectApproxEqAbs(-6.0, cr[1], 1e-9);
    try zigstd.testExpectApproxEqAbs(11.0, cr[2], 1e-9);
    try zigstd.testExpectApproxEqAbs(-6.0, cr[3], 1e-9);

    var wr: [3]f64 = @splat(0.0);
    var wi: [3]f64 = @splat(0.0);
    const n = try roots(zigstd.testAllocator, std.testing.io, &cr, &ci, &wr, &wi);
    try zigstd.testExpectEqual(@as(usize, 3), n);
    var sorted = wr;
    std.mem.sort(f64, &sorted, {}, std.sort.asc(f64));
    try zigstd.testExpectApproxEqAbs(1.0, sorted[0], 1e-7);
    try zigstd.testExpectApproxEqAbs(2.0, sorted[1], 1e-7);
    try zigstd.testExpectApproxEqAbs(3.0, sorted[2], 1e-7);
}

test "sin/cos/exp/sqrt/log/atan apply elementwise" {
    var ar = [_]f64{ 0.0, std.math.pi / 2.0 };
    var ai: [2]f64 = @splat(0.0);
    sin(&ar, &ai);
    try zigstd.testExpectApproxEqAbs(0.0, ar[0], 1e-9);
    try zigstd.testExpectApproxEqAbs(1.0, ar[1], 1e-9);

    var br = [_]f64{ 4.0, 9.0 };
    var bi: [2]f64 = @splat(0.0);
    sqrt(&br, &bi);
    try zigstd.testExpectApproxEqAbs(2.0, br[0], 1e-9);
    try zigstd.testExpectApproxEqAbs(3.0, br[1], 1e-9);
}

test "expm matches a hand-verifiable diagonalizable case" {
    // A = diag(0, ln(2)) -> expm(A) = diag(1, 2).
    const ar = [_]f64{ 0.0, 0.0, 0.0, @log(2.0) };
    const ai: [4]f64 = @splat(0.0);
    var out_r: [4]f64 = @splat(0.0);
    var out_i: [4]f64 = @splat(0.0);
    var ar_mut = ar;
    var ai_mut = ai;
    try expm(zigstd.testAllocator, std.testing.io, CxMat.initPacked(&ar_mut, &ai_mut, 2, 2), CxMat.initPacked(&out_r, &out_i, 2, 2));
    try zigstd.testExpectApproxEqAbs(1.0, out_r[ix(2, 1, 1)], 1e-8);
    try zigstd.testExpectApproxEqAbs(2.0, out_r[ix(2, 2, 2)], 1e-8);
    try zigstd.testExpectApproxEqAbs(0.0, out_r[ix(2, 1, 2)], 1e-8);
    try zigstd.testExpectApproxEqAbs(0.0, out_r[ix(2, 2, 1)], 1e-8);
}

test "powm matches a hand-verifiable diagonalizable case" {
    // A = diag(2, 3) -> A**3 = diag(8, 27).
    var ar = [_]f64{ 2.0, 0.0, 0.0, 3.0 };
    var ai: [4]f64 = @splat(0.0);
    var out_r: [4]f64 = @splat(0.0);
    var out_i: [4]f64 = @splat(0.0);
    try powm(zigstd.testAllocator, std.testing.io, CxMat.initPacked(&ar, &ai, 2, 2), Cx.init(3.0, 0.0), CxMat.initPacked(&out_r, &out_i, 2, 2));
    try zigstd.testExpectApproxEqAbs(8.0, out_r[ix(2, 1, 1)], 1e-8);
    try zigstd.testExpectApproxEqAbs(27.0, out_r[ix(2, 2, 2)], 1e-8);
    try zigstd.testExpectApproxEqAbs(0.0, out_r[ix(2, 1, 2)], 1e-8);
    try zigstd.testExpectApproxEqAbs(0.0, out_r[ix(2, 2, 1)], 1e-8);
}
