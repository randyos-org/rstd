//! `POLY`/`ROOTS`, and elementwise vs. matrix-valued elementary functions,
//! split out of `api.zig` to keep files under this project's line-count
//! guideline.

const std = @import("std");
const zigstd = @import("../../zigstd.zig");
const blas1 = @import("blas1.zig");
const eigen_general = @import("eigen_general.zig");
const eig_mod = @import("eig.zig");
const divide = @import("divide.zig");

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
/// matrix - which is already upper Hessenberg by construction, so this
/// skips `corth` and calls `comqr3` directly. Writes into `wr`/`wi` and
/// returns the degree found (`0` for the zero polynomial).
pub fn roots(allocator: zigstd.Allocator, coeffr: []const f64, coeffi: []const f64, wr: []f64, wi: []f64) !usize {
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

    const ortr = try allocator.alloc(f64, n);
    defer allocator.free(ortr);
    const orti = try allocator.alloc(f64, n);
    defer allocator.free(orti);
    @memset(ortr, 0.0);
    @memset(orti, 0.0);
    var dummy: [1]f64 = .{0.0};

    const ierr = eigen_general.comqr3(n, n, 1, n, ortr, orti, compr, compi, wr[0..n], wi[0..n], &dummy, &dummy, 0);
    if (ierr != 0) return error.EigenvalueDidNotConverge;
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
/// Applies `f` to the eigenvalues of the `n`x`n` matrix `ar`/`ai` (not
/// mutated) and reconstructs `f(A) = V * diag(f(lambda)) * inv(V)` - the
/// matrix function, not the entrywise one (see `sin`/`cos`/... above for
/// that). Only correct when `A` is diagonalizable; a defective matrix
/// (repeated eigenvalue without a full set of eigenvectors) makes `V`
/// singular and this will fail or blow up, matching the original, which
/// used the same eigendecomposition approach.
fn applyMatrixFunction(
    allocator: zigstd.Allocator,
    comptime f: anytype,
    ar: []const f64,
    ai: []const f64,
    n: usize,
    out_r: []f64,
    out_i: []f64,
) !void {
    const result = try eig_mod.eig(allocator, ar, ai, n, true);
    defer result.deinit(allocator);
    const vr = result.vectors_r.?;
    const vi = result.vectors_i.?;

    var j: usize = 1;
    while (j <= n) : (j += 1) {
        const s = f(Cx.init(result.values_r[j - 1], result.values_i[j - 1]));
        var i: usize = 1;
        while (i <= n) : (i += 1) {
            const t = blas1.wmul(Cx.init(vr[ix(n, i, j)], vi[ix(n, i, j)]), s);
            out_r[ix(n, i, j)] = t.re;
            out_i[ix(n, i, j)] = t.im;
        }
    }

    const v_copy_r = try allocator.alloc(f64, n * n);
    defer allocator.free(v_copy_r);
    const v_copy_i = try allocator.alloc(f64, n * n);
    defer allocator.free(v_copy_i);
    @memcpy(v_copy_r, vr);
    @memcpy(v_copy_i, vi);
    const ipvt = try allocator.alloc(usize, n);
    defer allocator.free(ipvt);
    const zr = try allocator.alloc(f64, n);
    defer allocator.free(zr);
    const zi = try allocator.alloc(f64, n);
    defer allocator.free(zi);
    const rowr = try allocator.alloc(f64, n);
    defer allocator.free(rowr);
    const rowi = try allocator.alloc(f64, n);
    defer allocator.free(rowi);

    _ = divide.mrdivide(out_r, out_i, n, n, n, v_copy_r, v_copy_i, n, ipvt, zr, zi, rowr, rowi);
}

pub fn expm(allocator: zigstd.Allocator, ar: []const f64, ai: []const f64, n: usize, out_r: []f64, out_i: []f64) !void {
    return applyMatrixFunction(allocator, std.math.complex.exp, ar, ai, n, out_r, out_i);
}
pub fn logm(allocator: zigstd.Allocator, ar: []const f64, ai: []const f64, n: usize, out_r: []f64, out_i: []f64) !void {
    return applyMatrixFunction(allocator, blas1.wlog, ar, ai, n, out_r, out_i);
}
pub fn sqrtm(allocator: zigstd.Allocator, ar: []const f64, ai: []const f64, n: usize, out_r: []f64, out_i: []f64) !void {
    return applyMatrixFunction(allocator, std.math.complex.sqrt, ar, ai, n, out_r, out_i);
}
pub fn sinm(allocator: zigstd.Allocator, ar: []const f64, ai: []const f64, n: usize, out_r: []f64, out_i: []f64) !void {
    return applyMatrixFunction(allocator, std.math.complex.sin, ar, ai, n, out_r, out_i);
}
pub fn cosm(allocator: zigstd.Allocator, ar: []const f64, ai: []const f64, n: usize, out_r: []f64, out_i: []f64) !void {
    return applyMatrixFunction(allocator, std.math.complex.cos, ar, ai, n, out_r, out_i);
}
pub fn atanm(allocator: zigstd.Allocator, ar: []const f64, ai: []const f64, n: usize, out_r: []f64, out_i: []f64) !void {
    return applyMatrixFunction(allocator, blas1.watan, ar, ai, n, out_r, out_i);
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
    const n = try roots(zigstd.testAllocator, &cr, &ci, &wr, &wi);
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
    try expm(zigstd.testAllocator, &ar, &ai, 2, &out_r, &out_i);
    try zigstd.testExpectApproxEqAbs(1.0, out_r[ix(2, 1, 1)], 1e-8);
    try zigstd.testExpectApproxEqAbs(2.0, out_r[ix(2, 2, 2)], 1e-8);
    try zigstd.testExpectApproxEqAbs(0.0, out_r[ix(2, 1, 2)], 1e-8);
    try zigstd.testExpectApproxEqAbs(0.0, out_r[ix(2, 2, 1)], 1e-8);
}
