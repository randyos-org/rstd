//! The singular value decomposition `A = U Σ Vᴴ`: LAPACK's `?GESVD` and
//! `?GESDD` (Reference-LAPACK's `SRC/dgesvd.f`, `dgesdd.f` and their `s`,
//! `c`, `z` forms, as OpenBLAS v0.3.34-317 ships them in `lapack-netlib/`:
//! arguments, `xerbla` order, the jobs and the scaling of `A`) over faer
//! v0.24.4's SVD (<https://codeberg.org/sarah-quinones/faer>, tag
//! `faer-v0.24.4`, MIT): `faer/src/linalg/svd/mod.rs`.
//!
//! Faer's SVD: a matrix more than `11 / 6` times as tall as wide is first
//! factored `Q R` (`qr.zig`), the `R` decomposed; the rest is bidiagonalized
//! (`bidiag.zig`), and the bidiagonal (taken as its conjugate transpose,
//! lower bidiagonal, with complex phases scaled out) decomposed
//! (`bidiag_svd.zig`); the reflections then build `U` and `V` with `gemm`.
//! A wide matrix is decomposed transposed. `gesvd` decomposes the bidiagonal
//! by QR iteration, as the reference's `xBDSQR` does; `gesdd` by divide and
//! conquer, as its `xBDSDC` does, with faer's threshold of 128 below which
//! it is QR iteration too.
//!
//! Departures, all documented in design.md's LAPACK section:
//! - no `WORK`, `LWORK`, `RWORK` or `IWORK`: workspace comes from the
//!   allocator, as the other routines here;
//! - failure is an error: `NoConvergence` for the reference's `INFO > 0`
//!   (without its count of unconverged superdiagonals, which faer does not
//!   keep), and also where `A` holds an infinity, or a NaN for `gesvd`,
//!   since faer's solver refuses a bidiagonal that is not finite
//!   (`compute_bidiag_real_svd`); `gesdd` returns `InvalidA` for a NaN, the
//!   reference's `INFO = -4` (`dgesdd.f:595`);
//! - four fixes to faer (issues.md): its tall path leaves `U`'s top right
//!   block as it finds it, here zeroed; reciprocals of denormal norms
//!   overflow, here divided by instead (`normalized`, `bidiag.zig`,
//!   `bidiag_svd.zig`); the secular solver's secant can stop far from a
//!   root by a pole, here bisected; a zero subproblem's `TwoRows` lacks its
//!   last 1 (`bidiag_svd.zig`).

const std = @import("std");
const ops = @import("ops.zig");
const hh = @import("householder.zig");
const qr = @import("qr.zig");
const bidiag = @import("bidiag.zig");
const bsvd = @import("bidiag_svd.zig");
const util = @import("../../blas/util.zig");
const chol = @import("../cholesky.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = ops.Mat;
const Float = ops.Float;
const isComplex = ops.isComplex;

/// `JOBU`, `JOBVT`, `JOBZ`: all of `U` (`m x m`) or `Vᴴ` (`n x n`), the first
/// `min(m, n)` columns of `U` or rows of `Vᴴ` (`s`), those overwritten on
/// `A` (`o`), or none (`n`).
pub const Job = enum { a, s, o, n };

/// What `gesvd` returns: `InvalidJobvt` when both jobs are `o`.
pub const GesvdError = error{ InvalidJobvt, InvalidM, InvalidN, InvalidLda, InvalidLdu, InvalidLdvt, NoConvergence } || Allocator.Error || Io.Cancelable;

/// What `gesdd` returns: `InvalidA` when `A` holds a NaN.
pub const GesddError = error{ InvalidM, InvalidN, InvalidA, InvalidLda, InvalidLdu, InvalidLdvt, NoConvergence } || Allocator.Error || Io.Cancelable;

/// `SvdParams::auto`'s `recursion_threshold`: below it, divide and conquer
/// is QR iteration.
pub const dc_threshold = 128;

/// `SvdParams::auto`'s `qr_ratio_threshold`.
const qr_ratio = 11.0 / 6.0;

/// `?GESVD`: `A = U Σ Vᴴ` for the `m x n` `A`; the singular values into `s`
/// (`min(m, n)`, nonincreasing, nonnegative), `U` by `jobu` into `u` (`ldu`)
/// or `A`, `Vᴴ` by `jobvt` into `vt` (`ldvt`) or `A`. `A` is destroyed
/// unless a job puts its vectors there. `u` and `vt` are not read where
/// their job does not use them. The bidiagonal by QR iteration.
pub fn gesvd(allocator: Allocator, io: Io, comptime T: type, jobu: Job, jobvt: Job, m: isize, n: isize, a: [*]T, lda: isize, s: [*]Float(T), u: [*]T, ldu: isize, vt: [*]T, ldvt: isize) GesvdError!void {
    comptime chol.floatType(T, "gesvd");
    const minmn = @min(m, n);
    if (jobu == .o and jobvt == .o) return error.InvalidJobvt;
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, m)) return error.InvalidLda;
    if (ldu < 1 or ((jobu == .a or jobu == .s) and ldu < m)) return error.InvalidLdu;
    if (ldvt < 1 or (jobvt == .a and ldvt < n) or (jobvt == .s and ldvt < minmn)) return error.InvalidLdvt;
    if (m == 0 or n == 0) return;
    if (!util.fitsMatrix(m, n, lda, @sizeOf(T))) return error.InvalidLda;
    if (!util.fitsMatrix(m, if (jobu == .a) m else if (jobu == .s) minmn else 0, ldu, @sizeOf(T))) return error.InvalidLdu;
    if (!util.fitsMatrix(if (jobvt == .a) n else if (jobvt == .s) minmn else 0, n, ldvt, @sizeOf(T))) return error.InvalidLdvt;
    run(allocator, io, T, jobu, jobvt, @intCast(m), @intCast(n), a, @intCast(lda), s, u, @intCast(ldu), vt, @intCast(ldvt), std.math.maxInt(usize)) catch |e| return switch (e) {
        error.NanInA => error.NoConvergence,
        else => |x| x,
    };
}

/// `?GESDD`: as `gesvd`, with one job for both: `a`, `s` and `n` as
/// `gesvd`'s for `U` and `Vᴴ` alike; `o` puts `U` on `A` and all of `Vᴴ` in
/// `vt` when `m >= n`, else all of `U` in `u` and `Vᴴ` on `A`. The
/// bidiagonal by divide and conquer.
pub fn gesdd(allocator: Allocator, io: Io, comptime T: type, jobz: Job, m: isize, n: isize, a: [*]T, lda: isize, s: [*]Float(T), u: [*]T, ldu: isize, vt: [*]T, ldvt: isize) GesddError!void {
    comptime chol.floatType(T, "gesdd");
    const minmn = @min(m, n);
    const as = jobz == .a or jobz == .s;
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, m)) return error.InvalidLda;
    if (ldu < 1 or (as and ldu < m) or (jobz == .o and m < n and ldu < m)) return error.InvalidLdu;
    if (ldvt < 1 or (jobz == .a and ldvt < n) or (jobz == .s and ldvt < minmn) or (jobz == .o and m >= n and ldvt < n)) return error.InvalidLdvt;
    if (m == 0 or n == 0) return;
    if (!util.fitsMatrix(m, n, lda, @sizeOf(T))) return error.InvalidLda;
    const ucols = switch (jobz) {
        .a => m,
        .s => minmn,
        .o => if (m < n) m else 0,
        .n => 0,
    };
    const vtrows = switch (jobz) {
        .a => n,
        .s => minmn,
        .o => if (m >= n) n else 0,
        .n => 0,
    };
    if (!util.fitsMatrix(m, ucols, ldu, @sizeOf(T))) return error.InvalidLdu;
    if (!util.fitsMatrix(vtrows, n, ldvt, @sizeOf(T))) return error.InvalidLdvt;
    const jobu: Job, const jobvt: Job = switch (jobz) {
        .o => if (m >= n) .{ .o, .a } else .{ .a, .o },
        .a, .s, .n => .{ jobz, jobz },
    };
    run(allocator, io, T, jobu, jobvt, @intCast(m), @intCast(n), a, @intCast(lda), s, u, @intCast(ldu), vt, @intCast(ldvt), dc_threshold) catch |e| return switch (e) {
        error.NanInA => error.InvalidA,
        else => |x| x,
    };
}

const RunError = error{ NanInA, NoConvergence } || ops.Error;

/// Both drivers after their checks: the scaling of `dgesvd.f:652-667` and
/// `dgesdd.f:585-606`, faer's SVD, and the vectors where the jobs put them.
fn run(al: Allocator, io: Io, comptime T: type, jobu: Job, jobvt: Job, m: usize, n: usize, a: [*]T, lda: usize, s: [*]Float(T), u: [*]T, ldu: usize, vt: [*]T, ldvt: usize, threshold: usize) RunError!void {
    const F = Float(T);
    const minmn = @min(m, n);
    const am: Mat(T) = .{ .ptr = a, .rows = m, .cols = n, .ld = lda };
    // `DLANGE('M')`, NaN first
    var anrm: F = 0;
    for (0..n) |j| for (0..m) |i| {
        const x = ops.abs(T, am.get(i, j));
        if (std.math.isNan(x)) return error.NanInA;
        anrm = @max(anrm, x);
    };
    if (std.math.isInf(anrm)) return error.NoConvergence;
    const smlnum = @sqrt(std.math.floatMin(F)) / std.math.floatEps(F);
    const bignum = 1 / smlnum;
    const to: ?F = if (anrm > 0 and anrm < smlnum) smlnum else if (anrm > bignum) bignum else null;
    if (to) |t| lascl(T, anrm, t, am);
    // U: into `u`, or for `o` into a copy that goes to `A`
    const ucols: usize = switch (jobu) {
        .a => m,
        .s, .o => minmn,
        .n => 0,
    };
    const want_u = jobu != .n;
    const want_v = jobvt != .n;
    const uo: ?Mat(T) = if (jobu == .o) try ops.alloc(T, al, m, minmn) else null;
    defer if (uo) |x| ops.free(T, al, x);
    const um: ?Mat(T) = if (!want_u) null else uo orelse Mat(T){ .ptr = u, .rows = m, .cols = ucols, .ld = ldu };
    const vcols: usize = if (jobvt == .a) n else minmn;
    const vm: ?Mat(T) = if (want_v) try ops.alloc(T, al, n, vcols) else null;
    defer if (vm) |x| ops.free(T, al, x);
    const sv = s[0..minmn];
    try svd(al, io, T, am, sv, um, vm, threshold);
    if (to) |t| {
        // `DLASCL` on S, back
        lascl(F, t, anrm, .{ .ptr = sv.ptr, .rows = minmn, .cols = 1, .ld = @max(minmn, 1) });
    }
    if (uo) |x| am.sub(0, 0, m, minmn).copyFrom(x);
    if (vm) |x| {
        // Vᴴ into `vt`, or for `o` into A's first rows
        const dst: Mat(T) = if (jobvt == .o) am.sub(0, 0, minmn, n) else .{ .ptr = vt, .rows = vcols, .cols = n, .ld = ldvt };
        for (0..n) |j| for (0..vcols) |i| {
            dst.at(i, j).* = ops.conj(T, x.get(j, i));
        };
    }
}

/// `DLASCL` type `G`: `A := A * (to / from)`, in steps that neither
/// overflow nor underflow (`lapack-netlib/SRC/dlascl.f`). `from` is finite
/// and nonzero.
fn lascl(comptime T: type, from: Float(T), to: Float(T), a: Mat(T)) void {
    const F = Float(T);
    const smlnum = std.math.floatMin(F);
    const bignum = 1 / smlnum;
    var cfromc = from;
    var ctoc = to;
    while (true) {
        const cfrom1 = cfromc * smlnum;
        var mul: F = undefined;
        var done = true;
        if (cfrom1 == cfromc) {
            mul = ctoc / cfromc;
        } else {
            const cto1 = ctoc / bignum;
            if (cto1 == ctoc) {
                mul = ctoc;
                cfromc = 1;
            } else if (@abs(cfrom1) > @abs(ctoc) and ctoc != 0) {
                mul = smlnum;
                done = false;
                cfromc = cfrom1;
            } else if (@abs(cto1) > @abs(cfromc)) {
                mul = bignum;
                done = false;
                ctoc = cto1;
            } else {
                mul = ctoc / cfromc;
                if (mul == 1) return;
            }
        }
        for (0..a.cols) |j| for (a.col(j)[0..a.rows]) |*x| {
            x.* = ops.mulReal(T, x.*, mul);
        };
        if (done) return;
    }
}

/// `svd`: `s` (`min(m, n)`), `u` (`m x m` or `m x min(m, n)`) and `v`
/// (`n x n` or `n x min(m, n)`) from `a`, which is only read; `threshold`
/// as `bidiagSvd`'s. A wide `a` is decomposed transposed, `U` and `V`
/// swapped and conjugated after.
pub fn svd(al: Allocator, io: Io, comptime T: type, a: Mat(T), s: []Float(T), u: ?Mat(T), v: ?Mat(T), threshold: usize) bsvd.Error!void {
    const m = a.rows;
    const n = a.cols;
    if (m >= n) return svdTall(al, io, T, a, s, u, v, threshold);
    const at = try ops.alloc(T, al, n, m);
    defer ops.free(T, al, at);
    for (0..m) |j| for (0..n) |i| {
        at.at(i, j).* = a.get(j, i);
    };
    try svdTall(al, io, T, at, s, v, u, threshold);
    if (comptime isComplex(T)) for ([_]?Mat(T){ u, v }) |mo| if (mo) |x| {
        for (0..x.cols) |j| for (x.col(j)[0..x.rows]) |*e| {
            e.* = ops.conj(T, e.*);
        };
    };
}

/// `svd` for `m >= n`.
fn svdTall(al: Allocator, io: Io, comptime T: type, a: Mat(T), s: []Float(T), u: ?Mat(T), v: ?Mat(T), threshold: usize) bsvd.Error!void {
    const m = a.rows;
    const n = a.cols;
    if (n == 0) {
        if (u) |x| x.identity();
        return;
    }
    if (@as(f64, @floatFromInt(m)) / @as(f64, @floatFromInt(n)) <= qr_ratio) return squareish(al, io, T, a, s, u, v, threshold);
    const bs = qr.blockSize(m, n);
    const q = try ops.alloc(T, al, m, n);
    defer ops.free(T, al, q);
    q.copyFrom(a);
    const h = try ops.alloc(T, al, bs, n);
    defer ops.free(T, al, h);
    _ = try qr.factor(al, io, T, q, h);
    {
        const r = try ops.alloc(T, al, n, n);
        defer ops.free(T, al, r);
        for (0..n) |j| for (0..n) |i| {
            r.at(i, j).* = if (i <= j) q.get(i, j) else ops.zero(T);
        };
        try squareish(al, io, T, r, s, if (u) |x| x.sub(0, 0, n, n) else null, v, threshold);
    }
    if (u) |x| {
        x.sub(n, 0, m - n, x.cols).fill(ops.zero(T));
        if (x.cols > n) x.sub(0, n, n, x.cols - n).fill(ops.zero(T));
        for (n..x.cols) |i| x.at(i, i).* = ops.one(T);
        try hh.applySequence(al, io, T, q, h, false, x);
    }
}

/// `svd_imp`: bidiagonalize, decompose `Bᴴ`, rebuild.
fn squareish(al: Allocator, io: Io, comptime T: type, a: Mat(T), s: []Float(T), u: ?Mat(T), v: ?Mat(T), threshold: usize) bsvd.Error!void {
    const m = a.rows;
    const n = a.cols;
    const bs = qr.blockSize(m, n);
    const bid = try ops.alloc(T, al, m, n);
    defer ops.free(T, al, bid);
    bid.copyFrom(a);
    const hl = try ops.alloc(T, al, bs, n);
    defer ops.free(T, al, hl);
    const hr = try ops.alloc(T, al, bs, n - 1);
    defer ops.free(T, al, hr);
    try bidiag.bidiagonalize(al, io, T, bid, hl, hr);
    const de = try al.alloc(T, 2 * n);
    defer al.free(de);
    const diag = de[0..n];
    const subdiag = de[n..];
    for (0..n) |i| {
        diag[i] = ops.conj(T, bid.get(i, i));
        subdiag[i] = if (i + 1 < n) ops.conj(T, bid.get(i, i + 1)) else ops.zero(T);
    }
    // Bᴴ = Ub S Vbᵀ, so U is Vb and V is Ub
    const ub: ?Mat(T) = if (v != null) try ops.alloc(T, al, n + 1, n + 1) else null;
    defer if (ub) |x| ops.free(T, al, x);
    const vb: ?Mat(T) = if (u != null) try ops.alloc(T, al, n, n) else null;
    defer if (vb) |x| ops.free(T, al, x);
    try bidiagSvd(al, io, T, diag, subdiag, ub, vb, threshold);
    for (s, diag) |*x, d| x.* = ops.re(T, d);
    if (u) |x| {
        const ncols = x.cols;
        x.sub(0, 0, n, n).copyFrom(vb.?);
        x.sub(n, 0, m - n, ncols).fill(ops.zero(T));
        x.sub(0, n, n, ncols - n).fill(ops.zero(T));
        for (n..ncols) |i| x.at(i, i).* = ops.one(T);
        try hh.applySequence(al, io, T, bid, hl, false, x);
    }
    if (v) |x| {
        x.copyFrom(ub.?.sub(0, 0, n, n));
        for (1..n) |j| for (0..j) |i| {
            bid.at(j, i).* = bid.get(i, j);
        };
        if (n > 1) try hh.applySequence(al, io, T, bid.sub(1, 0, n - 1, n - 1), hr, true, x.sub(1, 0, n - 1, x.cols));
    }
}

/// `compute_bidiag_real_svd` or, complex, `compute_bidiag_cplx_svd`: the
/// lower bidiagonal's SVD, `u` (`(n + 1) x (n + 1)`) and `v` (`n x n`)
/// where given; divide and conquer from `threshold` on, else QR iteration.
fn bidiagSvd(al: Allocator, io: Io, comptime T: type, diag: []T, subdiag: []T, u: ?Mat(T), v: ?Mat(T), threshold: usize) bsvd.Error!void {
    if (comptime !isComplex(T)) return realBidiagSvd(al, io, T, diag, subdiag, u, v, threshold);
    const F = Float(T);
    const n = diag.len;
    const reals = try al.alloc(F, 2 * n);
    defer al.free(reals);
    const dr = reals[0..n];
    const sr = reals[n..];
    const muls = try al.alloc(T, 2 * n);
    defer al.free(muls);
    const col_mul = muls[0..n];
    const row_mul = muls[n..];
    var cn = normalized(T, ops.conj(T, diag[0]));
    col_mul[0] = cn;
    dr[0] = ops.abs(T, diag[0]);
    sr[n - 1] = 0;
    for (1..n) |i| {
        const rn = normalized(T, ops.conj(T, ops.mul(T, subdiag[i - 1], cn)));
        sr[i - 1] = ops.abs(T, subdiag[i - 1]);
        row_mul[i - 1] = ops.conj(T, rn);
        cn = normalized(T, ops.conj(T, ops.mul(T, diag[i], rn)));
        dr[i] = ops.abs(T, diag[i]);
        col_mul[i] = cn;
    }
    const ur: ?Mat(F) = if (u != null) try ops.alloc(F, al, n + 1, n + 1) else null;
    defer if (ur) |x| ops.free(F, al, x);
    const vr: ?Mat(F) = if (v != null) try ops.alloc(F, al, n, n) else null;
    defer if (vr) |x| ops.free(F, al, x);
    try realBidiagSvd(al, io, F, dr, sr, ur, vr, threshold);
    for (diag, dr) |*d, x| d.* = ops.scalar(T, x);
    if (u) |x| {
        const r = ur.?;
        for (0..x.cols) |j| {
            x.at(0, j).* = ops.scalar(T, r.get(0, j));
            x.at(n, j).* = ops.scalar(T, r.get(n, j));
            for (1..n) |i| x.at(i, j).* = ops.mulReal(T, row_mul[i - 1], r.get(i, j));
        }
    }
    if (v) |x| {
        const r = vr.?;
        for (0..x.cols) |j| for (0..n) |i| {
            x.at(i, j).* = ops.mulReal(T, col_mul[i], r.get(i, j));
        };
    }
}

/// `x / |x|` (1 for 0), scaled first so `|x|` cannot overflow. Faer
/// multiplies by `1 / norm1`, infinite for a denormal `x`; that divides.
fn normalized(comptime T: type, x: T) T {
    if (ops.isZero(T, x)) return ops.one(T);
    const norm1 = @max(@abs(x.re), @abs(x.im));
    const y: T = if (norm1 >= std.math.floatMin(Float(T))) ops.mulReal(T, x, 1 / norm1) else .init(x.re / norm1, x.im / norm1);
    return ops.mulReal(T, y, 1 / ops.abs(T, y));
}

fn realBidiagSvd(al: Allocator, io: Io, comptime F: type, diag: []F, subdiag: []F, u: ?Mat(F), v: ?Mat(F), threshold: usize) bsvd.Error!void {
    const n = diag.len;
    for (diag, subdiag) |d, e| if (!(std.math.isFinite(d) and std.math.isFinite(e))) return error.NoConvergence;
    if (n < threshold) {
        if (u) |x| x.identity();
        if (v) |x| x.identity();
        return bsvd.qrAlgorithm(F, diag, subdiag, if (u) |x| x.sub(0, 0, n, n) else null, v);
    }
    if (u) |x| return bsvd.divideAndConquer(al, io, F, diag, subdiag, .{ .full = x }, v, threshold);
    const u_two = try ops.alloc(F, al, 2, n + 1);
    defer ops.free(F, al, u_two);
    return bsvd.divideAndConquer(al, io, F, diag, subdiag, .{ .two_rows_storage = u_two }, v, threshold);
}

// ---- tests -----------------------------------------------------------------------

const testing = std.testing;
const check = @import("../../blas/_test/blas1/check.zig");
const types = @import("../../blas/types.zig");
const blas = @import("../../blas/__root__.zig");
const threading = @import("../../blas/threading.zig");

const test_types = .{ f32, f64, types.Complex32, types.Complex64 };

/// Faer's `test_svd` tolerance: `ApproxEq::eps()` (`128 eps`) times
/// `sqrt(8 max(m, n))`, here against `max(1, max |A|)` per element.
fn tolerance(comptime F: type, m: usize, n: usize, amax: F) F {
    return 128 * std.math.floatEps(F) * @sqrt(8 * @as(F, @floatFromInt(@max(m, n)))) * @max(amax, 1);
}

fn maxAbs(comptime T: type, a: Mat(T)) Float(T) {
    var x: Float(T) = 0;
    for (0..a.cols) |j| for (0..a.rows) |i| {
        x = @max(x, ops.abs(T, a.get(i, j)));
    };
    return x;
}

/// `‖A - U Σ Vᴴ‖` per element, `‖Uᴴ U - I‖` and `‖V Vᴴ - I‖` (over `Vᴴ`'s
/// rows) within `tolerance`; `s` nonincreasing and nonnegative.
fn expectSvd(comptime T: type, al: Allocator, a0: Mat(T), s: []const Float(T), u: Mat(T), vt: Mat(T), what: []const u8) !void {
    const F = Float(T);
    const m = a0.rows;
    const n = a0.cols;
    const k = s.len;
    const tol = tolerance(F, m, n, maxAbs(T, a0));
    for (0..k) |i| {
        if (!(s[i] >= 0) or (i > 0 and !(s[i - 1] >= s[i]))) {
            check.say("{s} {s} {d}x{d}: s[{d}] = {e} out of order\n", .{ what, @typeName(T), m, n, i, s[i] });
            return error.TestUnexpectedResult;
        }
    }
    // U Σ, then times Vᴴ
    const us = try ops.alloc(T, al, m, k);
    defer ops.free(T, al, us);
    for (0..k) |j| for (0..m) |i| {
        us.at(i, j).* = ops.mulReal(T, u.get(i, j), s[j]);
    };
    const r = try ops.alloc(T, al, m, n);
    defer ops.free(T, al, r);
    r.copyFrom(a0);
    try ops.gemm(al, check.io, T, .n, .n, ops.neg(T, ops.one(T)), us, vt.sub(0, 0, k, n), ops.one(T), r);
    if (maxAbs(T, r) > tol) {
        check.say("{s} {s} {d}x{d}: A - U S Vh = {e}, tol {e}\n", .{ what, @typeName(T), m, n, maxAbs(T, r), tol });
        return error.TestUnexpectedResult;
    }
    const otol = tolerance(F, m, n, 1);
    const g = try ops.alloc(T, al, @max(u.cols, vt.rows), @max(u.cols, vt.rows));
    defer ops.free(T, al, g);
    const gu = g.sub(0, 0, u.cols, u.cols);
    try ops.gemm(al, check.io, T, .c, .n, ops.one(T), u, u, ops.zero(T), gu);
    const gv = g.sub(0, 0, vt.rows, vt.rows);
    for ([_]Mat(T){ gu, gv }, 0..) |x, w| {
        if (w == 1) try ops.gemm(al, check.io, T, .n, .c, ops.one(T), vt, vt, ops.zero(T), gv);
        for (0..x.cols) |j| for (0..x.rows) |i| {
            const want = if (i == j) ops.one(T) else ops.zero(T);
            if (!(ops.abs(T, ops.sub(T, x.get(i, j), want)) <= otol)) {
                check.say("{s} {s} {d}x{d}: {s} not orthonormal at ({d},{d})\n", .{ what, @typeName(T), m, n, if (w == 0) "U" else "Vh", i, j });
                return error.TestUnexpectedResult;
            }
        };
    }
}

const Routine = enum { gesvd, gesdd };

/// The routine on a copy of `a0` with jobs `ju`, `jv` (`gesdd`: `ju`), the
/// vectors gathered from wherever the jobs put them.
fn decompose(comptime T: type, al: Allocator, io: Io, r: Routine, ju: Job, jv: Job, a0: Mat(T), s: []Float(T), u: Mat(T), vt: Mat(T)) !void {
    const m = a0.rows;
    const n = a0.cols;
    const lda = m + 1;
    const a = try ops.alloc(T, al, lda, n);
    defer ops.free(T, al, a);
    a.fill(chol.poisoned(T));
    a.sub(0, 0, m, n).copyFrom(a0);
    switch (r) {
        .gesvd => try gesvd(al, io, T, ju, jv, @intCast(m), @intCast(n), a.ptr, @intCast(lda), s.ptr, u.ptr, @intCast(u.ld), vt.ptr, @intCast(vt.ld)),
        .gesdd => try gesdd(al, io, T, ju, @intCast(m), @intCast(n), a.ptr, @intCast(lda), s.ptr, u.ptr, @intCast(u.ld), vt.ptr, @intCast(vt.ld)),
    }
    const k = @min(m, n);
    const uo = if (r == .gesdd) ju == .o and m >= n else ju == .o;
    const vo = if (r == .gesdd) ju == .o and m < n else jv == .o;
    if (uo) u.sub(0, 0, m, k).copyFrom(a.sub(0, 0, m, k));
    if (vo) vt.sub(0, 0, k, n).copyFrom(a.sub(0, 0, k, n));
    // the padding row
    for (0..n) |j| try testing.expect(chol.isNanT(T, a.get(m, j)));
}

fn svdCases(comptime T: type, al: Allocator, r: Routine, m: usize, n: usize, a0: Mat(T)) !void {
    const k = @min(m, n);
    const s = try al.alloc(Float(T), @max(k, 1));
    defer al.free(s);
    const s2 = try al.alloc(Float(T), @max(k, 1));
    defer al.free(s2);
    const u = try ops.alloc(T, al, m, m);
    defer ops.free(T, al, u);
    const vt = try ops.alloc(T, al, n, n);
    defer ops.free(T, al, vt);
    const jobs = [_][2]Job{ .{ .a, .a }, .{ .s, .s }, .{ .o, .s }, .{ .s, .o }, .{ .a, .n }, .{ .n, .a } };
    const name = @tagName(r);
    for (jobs) |j| {
        if (r == .gesdd and j[0] != j[1] and !(j[0] == .o)) continue;
        u.fill(chol.poisoned(T));
        vt.fill(chol.poisoned(T));
        try decompose(T, al, check.io, r, j[0], j[1], a0, s[0..k], u, vt);
        const gesdd_o = r == .gesdd and j[0] == .o;
        const ucols: usize = if (gesdd_o) (if (m >= n) k else m) else switch (j[0]) {
            .a => m,
            .s, .o => k,
            .n => 0,
        };
        const vrows: usize = if (gesdd_o) (if (m >= n) n else k) else switch (if (r == .gesdd) j[0] else j[1]) {
            .a => n,
            .s, .o => k,
            .n => 0,
        };
        if (ucols > 0 and vrows > 0) try expectSvd(T, al, a0, s[0..k], u.sub(0, 0, m, ucols), vt.sub(0, 0, vrows, n), name);
    }
    // the values alone agree
    u.fill(chol.poisoned(T));
    vt.fill(chol.poisoned(T));
    try decompose(T, al, check.io, r, .n, .n, a0, s2[0..k], u, vt);
    try decompose(T, al, check.io, r, .a, .a, a0, s[0..k], u, vt);
    const tol = tolerance(Float(T), m, n, maxAbs(T, a0));
    for (s[0..k], s2[0..k]) |x, y| try testing.expect(@abs(x - y) <= tol);
}

/// Random `m x n`, faer's `test_real` and `test_cplx` shapes and more:
/// tall and wide past `11 / 6` (the QR path), and past 128 (divide and
/// conquer).
const shapes = [_][2]usize{ .{ 1, 1 }, .{ 2, 2 }, .{ 3, 2 }, .{ 2, 3 }, .{ 4, 4 }, .{ 15, 10 }, .{ 10, 15 }, .{ 17, 17 }, .{ 50, 50 }, .{ 150, 20 }, .{ 20, 150 }, .{ 140, 135 } };

test "faer svd: gesvd and gesdd, A = U S Vh with U and V orthonormal, every job, square, tall and wide" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        for (shapes) |mn| {
            const m = mn[0];
            const n = mn[1];
            const a0 = try ops.alloc(T, al, m, n);
            defer ops.free(T, al, a0);
            check.fill(T, a0.ptr[0 .. m * n], 13 * m + n);
            for ([_]Routine{ .gesvd, .gesdd }) |r| try svdCases(T, al, r, m, n, a0);
        }
    }
}

test "faer svd: special matrices, faer's test_special: zero, ones, identity; rank deficient; tiny and huge scales" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        for ([_][2]usize{ .{ 3, 2 }, .{ 2, 2 }, .{ 15, 10 }, .{ 10, 15 }, .{ 50, 50 }, .{ 150, 20 }, .{ 20, 150 }, .{ 130, 130 } }) |mn| {
            const m = mn[0];
            const n = mn[1];
            const a0 = try ops.alloc(T, al, m, n);
            defer ops.free(T, al, a0);
            for (0..7) |kind| {
                switch (kind) {
                    0 => a0.fill(ops.zero(T)),
                    1 => a0.fill(ops.one(T)),
                    2 => a0.identity(),
                    3 => {
                        // rank 3
                        const x = try ops.alloc(T, al, m, 3);
                        defer ops.free(T, al, x);
                        const y = try ops.alloc(T, al, 3, n);
                        defer ops.free(T, al, y);
                        check.fill(T, x.ptr[0 .. m * 3], 1);
                        check.fill(T, y.ptr[0 .. 3 * n], 2);
                        try ops.gemm(al, check.io, T, .n, .n, ops.one(T), x, y, ops.zero(T), a0);
                    },
                    4, 5, 6 => {
                        check.fill(T, a0.ptr[0 .. m * n], m + n);
                        // past the scaling bounds both ways, and a denormal
                        const sc: F = switch (kind) {
                            4 => std.math.floatMin(F) * 64,
                            5 => std.math.floatMax(F) / 64,
                            else => std.math.floatTrueMin(F) * 8,
                        };
                        if (kind == 6) {
                            a0.fill(ops.zero(T));
                            a0.at(m / 2, n / 2).* = ops.scalar(T, sc);
                        } else for (0..n) |j| for (a0.col(j)[0..m]) |*x| {
                            x.* = ops.mulReal(T, x.*, sc);
                        };
                    },
                    else => unreachable,
                }
                for ([_]Routine{ .gesvd, .gesdd }) |r| {
                    const k = @min(m, n);
                    const s = try al.alloc(F, k);
                    defer al.free(s);
                    const u = try ops.alloc(T, al, m, m);
                    defer ops.free(T, al, u);
                    const vt = try ops.alloc(T, al, n, n);
                    defer ops.free(T, al, vt);
                    decompose(T, al, check.io, r, .a, .a, a0, s, u, vt) catch |e| {
                        check.say("{s} {s} {d}x{d} kind {d}: {t}\n", .{ @tagName(r), @typeName(T), m, n, kind, e });
                        return e;
                    };
                    if (kind >= 4) {
                        // relative to the scale: compare the scaled-up values
                        // with those of the unscaled matrix
                        if (kind == 6) {
                            try testing.expectEqual(maxAbs(T, a0), s[0]);
                            for (s[1..]) |x| try testing.expectEqual(@as(F, 0), x);
                        } else {
                            const b = try ops.alloc(T, al, m, n);
                            defer ops.free(T, al, b);
                            check.fill(T, b.ptr[0 .. m * n], m + n);
                            const sb = try al.alloc(F, k);
                            defer al.free(sb);
                            try decompose(T, al, check.io, r, .n, .n, b, sb, u, vt);
                            const sc: F = if (kind == 4) std.math.floatMin(F) * 64 else std.math.floatMax(F) / 64;
                            const tol = tolerance(F, m, n, maxAbs(T, b));
                            for (s, sb) |x, y| try testing.expect(@abs(x / sc - y) <= tol);
                        }
                        continue;
                    }
                    try expectSvd(T, al, a0, s, u, vt, @tagName(r));
                    if (kind == 3) for (s[@min(3, k)..]) |x| try testing.expect(x <= tolerance(F, m, n, maxAbs(T, a0)));
                }
            }
        }
    }
}

test "faer svd: NaN and infinity: gesdd's INFO = -4 for NaN, else NoConvergence" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex32 }) |T| {
        const F = Float(T);
        for ([_]F{ std.math.nan(F), std.math.inf(F), -std.math.inf(F) }) |bad| {
            var a: [12]T = undefined;
            check.fill(T, &a, 3);
            a[7] = ops.scalar(T, bad);
            var s: [3]F = undefined;
            var u: [16]T = undefined;
            var vt: [9]T = undefined;
            const want_dd: anyerror = if (std.math.isNan(bad)) error.InvalidA else error.NoConvergence;
            var b = a;
            try testing.expectError(error.NoConvergence, gesvd(al, check.io, T, .a, .a, 4, 3, &b, 4, &s, &u, 4, &vt, 3));
            b = a;
            try testing.expectError(want_dd, gesdd(al, check.io, T, .a, 4, 3, &b, 4, &s, &u, 4, &vt, 3));
        }
    }
}

test "faer svd: reports the first bad argument and reads nothing when empty" {
    const fa = testing.failing_allocator;
    const nan = std.math.nan(f64);
    var a: [4]f64 = @splat(nan);
    var s: [2]f64 = @splat(nan);
    var u: [4]f64 = @splat(nan);
    var vt: [4]f64 = @splat(nan);
    try testing.expectError(error.InvalidJobvt, gesvd(fa, check.io, f64, .o, .o, -1, -1, &a, 0, &s, &u, 0, &vt, 0));
    try testing.expectError(error.InvalidM, gesvd(fa, check.io, f64, .a, .a, -1, -1, &a, 0, &s, &u, 0, &vt, 0));
    try testing.expectError(error.InvalidN, gesvd(fa, check.io, f64, .a, .a, 2, -1, &a, 0, &s, &u, 0, &vt, 0));
    try testing.expectError(error.InvalidLda, gesvd(fa, check.io, f64, .a, .a, 2, 2, &a, 1, &s, &u, 0, &vt, 0));
    try testing.expectError(error.InvalidLdu, gesvd(fa, check.io, f64, .a, .a, 2, 2, &a, 2, &s, &u, 1, &vt, 0));
    try testing.expectError(error.InvalidLdu, gesvd(fa, check.io, f64, .n, .a, 2, 2, &a, 2, &s, &u, 0, &vt, 0));
    try testing.expectError(error.InvalidLdvt, gesvd(fa, check.io, f64, .n, .a, 2, 2, &a, 2, &s, &u, 1, &vt, 1));
    try testing.expectError(error.InvalidLdvt, gesvd(fa, check.io, f64, .n, .s, 3, 2, &a, 3, &s, &u, 1, &vt, 1));
    try testing.expectError(error.InvalidM, gesdd(fa, check.io, f64, .a, -1, -1, &a, 0, &s, &u, 0, &vt, 0));
    try testing.expectError(error.InvalidN, gesdd(fa, check.io, f64, .a, 2, -1, &a, 0, &s, &u, 0, &vt, 0));
    try testing.expectError(error.InvalidLda, gesdd(fa, check.io, f64, .a, 2, 2, &a, 1, &s, &u, 0, &vt, 0));
    try testing.expectError(error.InvalidLdu, gesdd(fa, check.io, f64, .a, 2, 2, &a, 2, &s, &u, 1, &vt, 0));
    try testing.expectError(error.InvalidLdu, gesdd(fa, check.io, f64, .o, 1, 2, &a, 1, &s, &u, 0, &vt, 0));
    try testing.expectError(error.InvalidLdvt, gesdd(fa, check.io, f64, .o, 2, 2, &a, 2, &s, &u, 1, &vt, 1));
    try testing.expectError(error.InvalidLdvt, gesdd(fa, check.io, f64, .s, 3, 2, &a, 3, &s, &u, 3, &vt, 1));
    // empty: nothing read or written
    try gesvd(fa, check.io, f64, .a, .a, 0, 2, &a, 1, &s, &u, 1, &vt, 2);
    try gesdd(fa, check.io, f64, .a, 2, 0, &a, 2, &s, &u, 2, &vt, 1);
    for ([_][]const f64{ &a, &s, &u, &vt }) |x| for (x) |e| try testing.expect(std.math.isNan(e));
}

test "faer svd: 1 x 1 and a known 2 x 2, exact where the reference is" {
    const al = testing.allocator;
    var a = [_]f64{-3};
    var s: [1]f64 = undefined;
    var u: [1]f64 = undefined;
    var vt: [1]f64 = undefined;
    try gesvd(al, check.io, f64, .a, .a, 1, 1, &a, 1, &s, &u, 1, &vt, 1);
    try testing.expectEqual(@as(f64, 3), s[0]);
    try testing.expectEqual(@as(f64, -3), u[0] * s[0] * vt[0]);
    // [3 0; 4 5] has singular values sqrt(45) and sqrt(5)
    for ([_]Routine{ .gesvd, .gesdd }) |r| {
        var b = [_]f64{ 3, 4, 0, 5 };
        var s2: [2]f64 = undefined;
        var u_two: [4]f64 = undefined;
        var v2: [4]f64 = undefined;
        switch (r) {
            .gesvd => try gesvd(al, check.io, f64, .n, .n, 2, 2, &b, 2, &s2, &u_two, 1, &v2, 1),
            .gesdd => try gesdd(al, check.io, f64, .n, 2, 2, &b, 2, &s2, &u_two, 1, &v2, 1),
        }
        try testing.expectApproxEqRel(@sqrt(45.0), s2[0], 4e-16);
        try testing.expectApproxEqRel(@sqrt(5.0), s2[1], 4e-16);
    }
}

test "faer svd: gesdd gives one worker's bits at any worker count" {
    if (comptime !threading.hosted) return error.SkipZigTest;
    _ = blas.coreName();
    const al = testing.allocator;
    const cpus = blas.config.cpu_count;
    defer blas.config.cpu_count = cpus;
    defer threading.test_serial_max = null;
    var threaded: Io.Threaded = .init(al, .{ .async_limit = .limited(8) });
    defer threaded.deinit();
    inline for (.{ f64, types.Complex32 }) |T| {
        // tall past 11 / 6, square past 128
        for ([_][2]usize{ .{ 300, 90 }, .{ 140, 140 } }) |mn| {
            const m = mn[0];
            const n = mn[1];
            const a0 = try ops.alloc(T, al, m, n);
            defer ops.free(T, al, a0);
            check.fill(T, a0.ptr[0 .. m * n], 99);
            const k = @min(m, n);
            var outs: [3][]T = undefined;
            var svals: [3][]Float(T) = undefined;
            for (0..3) |c| {
                outs[c] = try al.alloc(T, m * m + n * n);
                svals[c] = try al.alloc(Float(T), k);
            }
            defer for (0..3) |c| {
                al.free(outs[c]);
                al.free(svals[c]);
            };
            for ([_]u32{ 1, 3, 8 }, 0..) |w, c| {
                blas.config.cpu_count = w;
                threading.test_serial_max = if (w == 1) null else 4096;
                const u: Mat(T) = .{ .ptr = outs[c].ptr, .rows = m, .cols = m, .ld = m };
                const vt: Mat(T) = .{ .ptr = outs[c].ptr + m * m, .rows = n, .cols = n, .ld = n };
                try decompose(T, al, if (w == 1) Io.failing else threaded.io(), .gesdd, .a, .a, a0, svals[c], u, vt);
            }
            for (1..3) |c| {
                try testing.expectEqualSlices(u8, std.mem.sliceAsBytes(outs[0]), std.mem.sliceAsBytes(outs[c]));
                try testing.expectEqualSlices(u8, std.mem.sliceAsBytes(svals[0]), std.mem.sliceAsBytes(svals[c]));
            }
        }
    }
}
