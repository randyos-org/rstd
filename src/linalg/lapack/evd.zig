//! Eigendecompositions: LAPACK's `?SYEVD`, `?HEEVD`, `?SYEV`, `?HEEV` and
//! `?GEEV` (Reference-LAPACK's `dsyevd.f`, `dsyev.f`, `dgeev.f` and their
//! `s`, `c`, `z` forms,
//! as OpenBLAS v0.3.34-317 ships them in `lapack-netlib/`: arguments,
//! `xerbla` order, the scaling of `A`) over faer v0.24.4's
//! (<https://codeberg.org/sarah-quinones/faer>, tag `faer-v0.24.4`, MIT)
//! `faer/src/linalg/evd/mod.rs`.
//!
//! Self-adjoint: tridiagonalize (`tridiag.zig`); take the real tridiagonal
//! (a complex one's subdiagonal by its moduli, the phases moved into the
//! vectors); decompose it (`tridiag_evd.zig`) by divide and conquer from
//! 128 up for `syevd` (the reference's `xSTEDC`), by QR iteration for
//! `syev` (its `xSTEQR`) and for values alone; apply the reflections.
//!
//! General (`geev`): reduce to Hessenberg form (`hessenberg.zig`), then to
//! Schur form by QR iteration (`real_schur.zig`, `complex_schur.zig`),
//! take the triangular factor's eigenvectors and multiply by the Schur
//! vectors; then `dgeev.f`'s normalization. Two parts come from
//! Reference-LAPACK, where faer has none: balancing (`dgeev.f`'s `xGEBAL`
//! and `xGEBAK`, `balance.zig`; faer's solvers then run over the whole
//! balanced matrix, its isolated parts already triangular), and the
//! triangular factor's eigenvectors (`xTREVC3`, `trevc.zig`), whose solves
//! scale against overflow where faer's back substitution
//! (`evd_from_real_schur_imp`) does not.
//!
//! Departures, all in design.md's LAPACK section: no `WORK`, `RWORK` or
//! `IWORK` (the allocator gives workspace); `NoConvergence` for the
//! reference's `INFO > 0`, and for a NaN or infinity in `A`, which faer's
//! solvers cannot take; with `jobz` `n`, `A` is left as it was.

const std = @import("std");
const ops = @import("ops.zig");
const hh = @import("householder.zig");
const qr = @import("qr.zig");
const tridiag = @import("tridiag.zig");
const tevd = @import("tridiag_evd.zig");
const hessenberg = @import("hessenberg.zig");
const real_schur = @import("real_schur.zig");
const complex_schur = @import("complex_schur.zig");
const balance = @import("balance.zig");
const schur = @import("schur.zig");
const trevc = @import("trevc.zig");
const util = @import("../blas/util.zig");
const chol = @import("cholesky.zig");
const types = @import("../blas/types.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = ops.Mat;
const Float = ops.Float;
const isComplex = ops.isComplex;
const Uplo = types.Uplo;

/// `JOBZ`, `JOBVL`, `JOBVR`: vectors (`v`) or not (`n`).
pub const Jobz = enum { n, v };

/// What `syevd`, `heevd`, `syev` and `heev` return.
pub const SyevError = error{ InvalidUplo, InvalidN, InvalidLda, NoConvergence } || Allocator.Error || Io.Cancelable;

/// `SelfAdjointEvdParams::auto`'s `recursion_threshold`: below it, divide
/// and conquer is QR iteration.
pub const dc_threshold = 128;

/// `?SYEVD`, `?HEEVD`: the eigenvalues of the `n x n` symmetric (Hermitian)
/// `A`, from its `uplo` triangle, into `w` (`n`, nondecreasing); with
/// `jobz` `v` the orthonormal eigenvectors over `A`. Divide and conquer.
pub fn syevd(allocator: Allocator, io: Io, comptime T: type, jobz: Jobz, uplo: Uplo, n: isize, a: [*]T, lda: isize, w: [*]Float(T)) SyevError!void {
    return driver(allocator, io, T, jobz, uplo, n, a, lda, w, dc_threshold);
}

/// `?SYEV`, `?HEEV`: as `syevd`, the tridiagonal by QR iteration.
pub fn syev(allocator: Allocator, io: Io, comptime T: type, jobz: Jobz, uplo: Uplo, n: isize, a: [*]T, lda: isize, w: [*]Float(T)) SyevError!void {
    return driver(allocator, io, T, jobz, uplo, n, a, lda, w, std.math.maxInt(usize));
}

fn driver(al: Allocator, io: Io, comptime T: type, jobz: Jobz, uplo: Uplo, n: isize, a: [*]T, lda: isize, w: [*]Float(T), threshold: usize) SyevError!void {
    comptime chol.floatType(T, "syevd");
    if (uplo.isRev()) return error.InvalidUplo;
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, n)) return error.InvalidLda;
    if (n == 0) return;
    if (!util.fitsMatrix(n, n, lda, @sizeOf(T))) return error.InvalidLda;
    const un: usize = @intCast(n);
    const am: Mat(T) = .{ .ptr = a, .rows = un, .cols = un, .ld = @intCast(lda) };
    return symmetric(al, io, T, jobz == .v, uplo.isLower(), am, w[0..un], threshold);
}

/// The drivers after their checks: `dsyevd.f:292-307`'s scaling into
/// `[sqrt(SMLNUM), sqrt(BIGNUM)]`, faer's `self_adjoint_evd`, the values
/// scaled back.
pub fn symmetric(al: Allocator, io: Io, comptime T: type, vectors: bool, lower: bool, a: Mat(T), w: []Float(T), threshold: usize) (ops.Error || error{NoConvergence})!void {
    const F = Float(T);
    const n = a.rows;
    // the triangle as a lower one (`dlansy('M')` on the way); the upper as
    // `J A J`'s lower (`J` the reversal), so faer's top-down reduction runs
    // from the bottom, as `xSYTRD` reduces `.u`, and a matrix graded large
    // at the bottom right keeps its small eigenvalues (T125); the vectors'
    // rows are reversed after
    const trid = try ops.alloc(T, al, n, n);
    defer ops.free(T, al, trid);
    var anrm: F = 0;
    for (0..n) |j| for (j..n) |i| {
        const x = if (lower) a.get(i, j) else a.get(n - 1 - i, n - 1 - j);
        trid.at(i, j).* = x;
        const v = ops.abs(T, x);
        if (!std.math.isFinite(v)) return error.NoConvergence;
        anrm = @max(anrm, v);
    };
    const smlnum = std.math.floatMin(F) / std.math.floatEps(F);
    const rmin = @sqrt(smlnum);
    const rmax = @sqrt(1 / smlnum);
    const sigma: ?F = if (anrm > 0 and anrm < rmin) rmin / anrm else if (anrm > rmax) rmax / anrm else null;
    if (sigma) |s| ops.lascl(T, 1, s, trid);
    try selfAdjointEvd(al, io, T, trid, w, if (vectors) a else null, threshold);
    if (vectors and !lower) for (0..n) |j| {
        const c = a.col(j);
        std.mem.reverse(T, c[0..n]);
    };
    if (sigma) |s| for (w) |*x| {
        x.* *= 1 / s;
    };
}

/// `self_adjoint_evd`: `trid`'s lower triangle (overwritten) to the
/// eigenvalues in `s`, nondecreasing, and the vectors in `u` if given.
pub fn selfAdjointEvd(al: Allocator, io: Io, comptime T: type, trid: Mat(T), s: []Float(T), u: ?Mat(T), threshold: usize) (ops.Error || error{NoConvergence})!void {
    const F = Float(T);
    const n = trid.rows;
    if (n == 0) return;
    const bs = qr.blockSize(n, n);
    const h = try ops.alloc(T, al, bs, @max(n, 2) - 1);
    defer ops.free(T, al, h);
    try tridiag.tridiagonalize(al, io, T, trid, h.sub(0, 0, bs, n - 1));
    const real = try al.alloc(F, 2 * n);
    defer al.free(real);
    const diag = real[0..n];
    const offdiag = real[n..];
    for (0..n) |i| {
        diag[i] = ops.re(T, trid.get(i, i));
        offdiag[i] = if (i + 1 < n) (if (comptime isComplex(T)) ops.abs(T, trid.get(i + 1, i)) else trid.get(i + 1, i)) else 0;
    }
    const um = u orelse {
        try tevd.qrAlgorithm(F, diag, offdiag, null);
        @memcpy(s, diag);
        return;
    };
    const ur: Mat(F) = if (comptime isComplex(T)) try ops.alloc(F, al, n, n) else um;
    defer if (comptime isComplex(T)) ops.free(F, al, ur);
    if (n < threshold) {
        try tevd.qrAlgorithm(F, diag, offdiag, ur);
    } else {
        try tevd.divideAndConquer(al, io, F, diag, offdiag, ur, threshold);
    }
    if (comptime isComplex(T)) {
        // U = diag(phases) Ur, the phases the subdiagonal's
        var x = ops.one(T);
        for (0..n) |i| {
            if (i > 0) x = normalized(T, ops.mul(T, trid.get(i, i - 1), x));
            for (0..n) |j| um.at(i, j).* = ops.mulReal(T, x, ur.get(i, j));
        }
    }
    if (n > 1) try hh.applySequence(al, io, T, trid.sub(1, 0, n - 1, n - 1), h.sub(0, 0, bs, n - 1), false, um.sub(1, 0, n - 1, n));
    @memcpy(s, diag);
}

/// `x / |x|` (1 for 0). Faer multiplies by `1 / |x|`, infinite for a
/// denormal `x`; that divides (as `svd.zig`'s).
fn normalized(comptime T: type, x: T) T {
    if (ops.isZero(T, x)) return ops.one(T);
    const norm1 = @max(@abs(x.re), @abs(x.im));
    const y: T = if (norm1 >= std.math.floatMin(Float(T))) ops.mulReal(T, x, 1 / norm1) else .init(x.re / norm1, x.im / norm1);
    return ops.mulReal(T, y, 1 / ops.abs(T, y));
}

// ---- the general eigenproblem ------------------------------------------------------

/// What `geev` returns: `InvalidA` when `A` holds a NaN (the reference's
/// `INFO = -4`, `dgeev.f:353`), `NoConvergence` for its `INFO > 0` and for
/// an infinity in `A`.
pub const GeevError = error{ InvalidN, InvalidA, InvalidLda, InvalidLdvl, InvalidLdvr, NoConvergence } || Allocator.Error || Io.Cancelable;

/// `geev`'s eigenvalues: real `wr` and `wi` (`dgeev`'s), complex `w`
/// (`zgeev`'s).
pub fn Eigenvalues(comptime T: type) type {
    return if (isComplex(T)) [*]T else struct { re: [*]T, im: [*]T };
}

/// `?GEEV`: the eigenvalues of the general `n x n` `A` into `w`, and with
/// `jobvl`, `jobvr` `v` its left (`uᴴ A = λ uᴴ`) and right (`A v = λ v`)
/// eigenvectors into `vl`, `vr`, each of unit 2-norm with its largest
/// entry real. Real: complex pairs come together, positive imaginary part
/// first, their vectors as `v(:, j) ± i v(:, j + 1)`. `A` is overwritten.
/// Balanced (`dgebal.f`, `balance.zig`), then the Hessenberg form's Schur
/// form by QR iteration and its eigenvectors by back substitution (faer's
/// `evd_real`, `evd_cplx`, over the whole balanced matrix).
///
/// On `error.NoConvergence` from the QR iteration, `unconverged` (when
/// given) gets the reference's `INFO = i > 0`: the first `i` eigenvalues
/// did not converge, `w[i..n]` hold those that did, and no eigenvectors are
/// computed (`vl`, `vr` undefined), as `dgeev.f` documents. For an infinity
/// in `A` it gets `n`.
pub fn geev(allocator: Allocator, io: Io, comptime T: type, jobvl: Jobz, jobvr: Jobz, n: isize, a: [*]T, lda: isize, w: Eigenvalues(T), vl: [*]T, ldvl: isize, vr: [*]T, ldvr: isize, unconverged: ?*usize) GeevError!void {
    comptime chol.floatType(T, "geev");
    const F = Float(T);
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, n)) return error.InvalidLda;
    if (ldvl < 1 or (jobvl == .v and ldvl < n)) return error.InvalidLdvl;
    if (ldvr < 1 or (jobvr == .v and ldvr < n)) return error.InvalidLdvr;
    if (n == 0) return;
    if (!util.fitsMatrix(n, n, lda, @sizeOf(T))) return error.InvalidLda;
    if (jobvl == .v and !util.fitsMatrix(n, n, ldvl, @sizeOf(T))) return error.InvalidLdvl;
    if (jobvr == .v and !util.fitsMatrix(n, n, ldvr, @sizeOf(T))) return error.InvalidLdvr;
    const un: usize = @intCast(n);
    const am: Mat(T) = .{ .ptr = a, .rows = un, .cols = un, .ld = @intCast(lda) };
    // `dgeev.f:341-359`: NaN, then the scale into [SMLNUM, BIGNUM]
    const anrm = ops.normMax(T, am);
    if (std.math.isNan(anrm)) return error.InvalidA;
    if (std.math.isInf(anrm)) {
        if (unconverged) |c| c.* = un;
        return error.NoConvergence;
    }
    const smlnum = @sqrt(std.math.floatMin(F)) / std.math.floatEps(F);
    const to = ops.scaleTo(F, anrm, smlnum);
    if (to) |t| ops.lascl(T, anrm, t, am);
    // `xGEBAL` (`balance.zig`), job `B`
    const scale = try allocator.alloc(F, un);
    defer allocator.free(scale);
    const bal = try balance.balance(io, T, am, scale);
    const ul: ?Mat(T) = if (jobvl == .v) Mat(T){ .ptr = vl, .rows = un, .cols = un, .ld = @intCast(ldvl) } else null;
    const ur: ?Mat(T) = if (jobvr == .v) Mat(T){ .ptr = vr, .rows = un, .cols = un, .ld = @intCast(ldvr) } else null;
    // `dgeev.f:398-404`: after a failure only `w[info..n]` is scaled back
    var info: usize = 0;
    const res = if (comptime isComplex(T))
        general(allocator, io, T, am, w[0..un], null, ul, ur, &info)
    else
        general(allocator, io, T, am, w.re[0..un], w.im[0..un], ul, ur, &info);
    res catch |e| if (e != error.NoConvergence) return e;
    if (to) |t| {
        const ws = if (comptime isComplex(T)) [_][*]T{w} else [_][*]T{ w.re, w.im };
        for (ws) |x| ops.lascl(T, t, anrm, .{ .ptr = x + info, .rows = un - info, .cols = 1, .ld = un });
    }
    res catch |e| {
        if (unconverged) |c| c.* = info;
        return e;
    };
    if (ul) |x| {
        try balance.back(io, T, bal, scale, true, x);
        normalizeVectors(io, T, x, if (comptime isComplex(T)) null else w.im);
    }
    if (ur) |x| {
        try balance.back(io, T, bal, scale, false, x);
        normalizeVectors(io, T, x, if (comptime isComplex(T)) null else w.im);
    }
}

/// `hseqr`'s job: the eigenvalues alone (`e`), or the Schur form too (`s`).
pub const HseqrJob = enum { e, s };
/// `hseqr`'s `Z`: none (`n`), from the identity (`i`), or accumulated into
/// the given one (`v`).
pub const Compz = enum { n, i, v };

/// What `hseqr` returns: `InvalidIlo` unless `1 <= ilo <= max(1, n)`,
/// `InvalidIhi` unless `min(ilo, n) <= ihi <= n`, `NoConvergence` for the
/// reference's `INFO > 0`.
pub const HseqrError = error{ InvalidN, InvalidIlo, InvalidIhi, InvalidLdh, InvalidLdz, NoConvergence };

/// `?HSEQR` (`dhseqr.f`, `zhseqr.f`): the eigenvalues of the upper
/// Hessenberg `n x n` `H` into `w`, as `geev`'s, by faer's `lahqr`
/// (`real_schur.zig`, `complex_schur.zig`; faer's multishift QR for large
/// sizes is not ported yet) on rows and columns `ilo..ihi` (1-based; the
/// rest triangular already, its eigenvalues read off the diagonal). With
/// `job` `s`, `H` becomes the Schur form `T` (real: quasi-triangular, its
/// 2 x 2 blocks standardized), else its contents are unspecified; entries
/// below the subdiagonal are zeroed first, as the reference clears them
/// after. With `compz` `i` or `v`, `Z` (`n x n`) is the identity or the
/// caller's, times the Schur vectors. `NoConvergence` writes the
/// reference's `INFO` to `unconverged`: `w[0..ilo - 1]` and
/// `w[INFO..n]` (0-based) converged.
pub fn hseqr(io: Io, comptime T: type, job: HseqrJob, compz: Compz, n: isize, ilo: isize, ihi: isize, h: [*]T, ldh: isize, w: Eigenvalues(T), z: [*]T, ldz: isize, unconverged: ?*usize) HseqrError!void {
    comptime chol.floatType(T, "hseqr");
    if (n < 0) return error.InvalidN;
    if (ilo < 1 or ilo > @max(1, n)) return error.InvalidIlo;
    if (ihi < @min(ilo, n) or ihi > n) return error.InvalidIhi;
    if (ldh < @max(1, n)) return error.InvalidLdh;
    if (ldz < 1 or (compz != .n and ldz < n)) return error.InvalidLdz;
    if (n == 0) return;
    if (!util.fitsMatrix(n, n, ldh, @sizeOf(T))) return error.InvalidLdh;
    if (compz != .n and !util.fitsMatrix(n, n, ldz, @sizeOf(T))) return error.InvalidLdz;
    const un: usize = @intCast(n);
    const lo: usize = @intCast(ilo - 1);
    const hi: usize = @intCast(ihi);
    const hm: Mat(T) = .{ .ptr = h, .rows = un, .cols = un, .ld = @intCast(ldh) };
    const zm: ?Mat(T) = if (compz == .n) null else .{ .ptr = z, .rows = un, .cols = un, .ld = @intCast(ldz) };
    if (compz == .i) zm.?.identity();
    for (0..un) |j| for (@min(j + 2, un)..un) |i| {
        hm.at(i, j).* = ops.zero(T);
    };
    // outside the block the matrix is triangular already
    for (0..un) |i| if (i < lo or i >= hi) {
        if (comptime isComplex(T)) w[i] = hm.get(i, i) else {
            w.re[i] = hm.get(i, i);
            w.im[i] = 0;
        }
    };
    const want = job == .s;
    const info = if (comptime isComplex(T))
        complex_schur.lahqr(T, want, hm, zm, w[0..un], lo, hi)
    else
        real_schur.lahqr(io, T, want, hm, zm, w.re[0..un], w.im[0..un], lo, hi);
    if (info != 0) {
        if (unconverged) |u| u.* = info;
        return error.NoConvergence;
    }
}

/// `gees`'s sort: which eigenvalues go first; real `(re, im)`, complex `λ`.
pub fn GeesSelect(comptime T: type) type {
    return if (isComplex(T)) *const fn (T) bool else *const fn (Float(T), Float(T)) bool;
}

/// What `gees` returns: `InvalidA` for a NaN in `A`, `NoConvergence` for the
/// reference's `INFO` 1 to `n`, `SelectChanged` for its `n + 1` (rounding
/// in the reordering changed which eigenvalues `select` picks),
/// `SwapRejected` for its `n + 2` (a swap too ill-conditioned to make).
pub const GeesError = error{ InvalidN, InvalidLda, InvalidLdvs, InvalidA, NoConvergence, SelectChanged, SwapRejected } || Allocator.Error || Io.Cancelable;

/// `?GEES` (`dgees.f`, `zgees.f`): the Schur factorization `A = Z T Zᴴ` of
/// the `n x n` `A`: `T` in place of `A` (real: quasi-triangular, its 2 x 2
/// blocks standardized), its eigenvalues in `w` (as `geev`'s), and with
/// `jobvs` `v` the Schur vectors `Z` in `vs`. As the reference: scaled into
/// `[SMLNUM, BIGNUM]`, balanced by permutation only (`gebal` job `p`),
/// `gehrd`, `orghr`, `hseqr`; with `select`, `trsen` then moves the
/// eigenvalues it picks to the top-left, `sdim` getting their count (a real
/// pair counts 2, and moves if either is picked); `Z` takes `gebak` back.
/// `NoConvergence` writes the reference's `INFO` to `unconverged`.
pub fn gees(allocator: Allocator, io: Io, comptime T: type, jobvs: Jobz, select: ?GeesSelect(T), n: isize, a: [*]T, lda: isize, sdim: *isize, w: Eigenvalues(T), vs: [*]T, ldvs: isize, unconverged: ?*usize) GeesError!void {
    comptime chol.floatType(T, "gees");
    const F = Float(T);
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, n)) return error.InvalidLda;
    if (ldvs < 1 or (jobvs == .v and ldvs < n)) return error.InvalidLdvs;
    sdim.* = 0;
    if (n == 0) return;
    if (!util.fitsMatrix(n, n, lda, @sizeOf(T))) return error.InvalidLda;
    if (jobvs == .v and !util.fitsMatrix(n, n, ldvs, @sizeOf(T))) return error.InvalidLdvs;
    const un: usize = @intCast(n);
    const am: Mat(T) = .{ .ptr = a, .rows = un, .cols = un, .ld = @intCast(lda) };
    // `dgees.f`: NaN, then the scale into [SMLNUM, BIGNUM], as `geev`
    const anrm = ops.normMax(T, am);
    if (std.math.isNan(anrm)) return error.InvalidA;
    if (std.math.isInf(anrm)) {
        if (unconverged) |c| c.* = un;
        return error.NoConvergence;
    }
    const smlnum = @sqrt(std.math.floatMin(F)) / std.math.floatEps(F);
    const to = ops.scaleTo(F, anrm, smlnum);
    if (to) |t| ops.lascl(T, anrm, t, am);
    const scale = try allocator.alloc(F, un);
    defer allocator.free(scale);
    var ilo: isize = 1;
    var ihi: isize = n;
    balance.gebal(allocator, io, T, .p, n, a, lda, &ilo, &ihi, scale.ptr) catch |e| switch (e) {
        error.Canceled => |x| return x,
        // checked above: no NaN, and the arguments
        else => unreachable,
    };
    const tau = try allocator.alloc(T, un);
    defer allocator.free(tau);
    hessenberg.gehrd(allocator, io, T, n, ilo, ihi, a, lda, tau.ptr) catch |e| return narrow(e);
    if (jobvs == .v) {
        const vm: Mat(T) = .{ .ptr = vs, .rows = un, .cols = un, .ld = @intCast(ldvs) };
        vm.copyFrom(am);
        hessenberg.orghr(allocator, io, T, n, ilo, ihi, vs, ldvs, tau.ptr) catch |e| return narrow(e);
    }
    const compz: Compz = if (jobvs == .v) .v else .n;
    hseqr(io, T, .s, compz, n, ilo, ihi, a, lda, w, vs, if (jobvs == .v) ldvs else 1, unconverged) catch |e| switch (e) {
        error.NoConvergence => |x| return x,
        else => unreachable,
    };
    if (select) |pick| {
        const bwork = try allocator.alloc(bool, un);
        defer allocator.free(bwork);
        for (bwork, 0..) |*b, i| b.* = picked(T, pick, w, i);
        var m: isize = 0;
        var s: F = 0;
        var sep: F = 0;
        schur.trsen(allocator, io, T, .n, if (jobvs == .v) .v else .n, bwork.ptr, n, a, lda, vs, if (jobvs == .v) ldvs else 1, w, &m, &s, &sep) catch |e| switch (e) {
            error.SwapRejected, error.OutOfMemory, error.Canceled => |x| return x,
            else => unreachable,
        };
        sdim.* = m;
        // `dgees.f`'s `INFO = N + 1`: the reordered eigenvalues, re-picked,
        // a real pair picked if either is
        const sm: usize = @intCast(m);
        var i: usize = 0;
        while (i < un) : (i += 1) {
            const pair = !isComplex(T) and i + 1 < un and (if (comptime isComplex(T)) false else w.im[i] > 0);
            const p = picked(T, pick, w, i) or (pair and picked(T, pick, w, i + 1));
            if (p != (i < sm)) return error.SelectChanged;
            if (pair) {
                if (p != (i + 1 < sm)) return error.SelectChanged;
                i += 1;
            }
        }
    }
    if (jobvs == .v) balance.gebak(allocator, io, T, .p, .r, n, ilo, ihi, scale.ptr, n, vs, ldvs) catch |e| switch (e) {
        error.Canceled => |x| return x,
        else => unreachable,
    };
    if (to) |t| {
        ops.lascl(T, t, anrm, am);
        const ws = if (comptime isComplex(T)) [_][*]T{w} else [_][*]T{ w.re, w.im };
        for (ws) |x| ops.lascl(T, t, anrm, .{ .ptr = x, .rows = un, .cols = 1, .ld = un });
    }
}

/// Whether `select` picks eigenvalue `i`, by its own `(re, im)` for real,
/// as `dgees.f` asks each.
fn picked(comptime T: type, pick: GeesSelect(T), w: Eigenvalues(T), i: usize) bool {
    return if (comptime isComplex(T)) pick(w[i]) else pick(w.re[i], w.im[i]);
}

/// `gehrd`'s and `orghr`'s errors here: their arguments are `gees`'s, checked.
fn narrow(e: hessenberg.GehrdError) GeesError {
    return switch (e) {
        error.OutOfMemory, error.Canceled => |x| x,
        else => unreachable,
    };
}

/// `evd_imp`: `a` (overwritten) to its eigenvalues (`s`, with `s_im` for
/// real), left vectors into `ul`, right into `ur` (unnormalized, real ones
/// packed as `geev`'s).
/// On `NoConvergence`, `info` gets the QR iteration's `INFO`: `s[info..]`
/// converged.
pub fn general(al: Allocator, io: Io, comptime T: type, a: Mat(T), s: []T, s_im: ?[]T, ul: ?Mat(T), ur: ?Mat(T), info: *usize) (ops.Error || error{NoConvergence})!void {
    const n = a.rows;
    const want = ul != null or ur != null;
    const bs = qr.blockSize(n - 1, n - 1);
    const h = try ops.alloc(T, al, bs, @max(n, 2) - 1);
    defer ops.free(T, al, h);
    const hseq = h.sub(0, 0, bs, n - 1);
    try hessenberg.hessenberg(al, io, T, a, hseq);
    const z: ?Mat(T) = if (want) try ops.alloc(T, al, n, n) else null;
    defer if (z) |x| ops.free(T, al, x);
    if (z) |x| {
        x.identity();
        if (n > 1) try hh.applySequenceRight(al, io, T, a.sub(1, 0, n - 1, n - 1), hseq, false, false, x.sub(1, 1, n - 1, n - 1));
    }
    for (0..n) |j| for (@min(j + 2, n)..n) |i| {
        a.at(i, j).* = ops.zero(T);
    };
    info.* = if (comptime isComplex(T)) complex_schur.lahqr(T, want, a, z, s, 0, n) else real_schur.lahqr(io, T, want, a, z, s, s_im.?, 0, n);
    if (info.* != 0) return error.NoConvergence;
    const zm = z orelse return;
    // `xTREVC3`'s vectors of the Schur form, then `Z X` (its blocked
    // back-transform)
    const x = try ops.alloc(T, al, n, n);
    defer ops.free(T, al, x);
    const rwork = try al.alloc(Float(T), n);
    defer al.free(rwork);
    const diag = try al.alloc(T, if (comptime isComplex(T)) n else 0);
    defer al.free(diag);
    for ([_]?Mat(T){ ur, ul }, [_]bool{ true, false }) |dst, right| if (dst) |u| {
        if (comptime isComplex(T)) {
            try trevc.complex(al, io, T, right, a, x, rwork, diag);
        } else trevc.real(T, right, a, x, rwork);
        try ops.gemm(al, io, T, .n, .n, ops.one(T), zm, x, ops.zero(T), u);
    };
}

/// `dgeev.f:430-487`, `zgeev.f:420-460`: each vector (a real pair as one)
/// to unit 2-norm, then turned so its largest entry (the first largest
/// `|v(k)|²`) is real.
fn normalizeVectors(io: Io, comptime T: type, v: Mat(T), wi: ?[*]const T) void {
    const n = v.rows;
    if (comptime isComplex(T)) {
        for (0..n) |j| {
            const col = v.col(j);
            const scl = 1 / ops.nrm2(io, T, n, col, 1);
            for (col[0..n]) |*x| x.* = ops.mulReal(T, x.*, scl);
            const k = argmaxAbs2(T, col, null, n);
            const t = ops.mulReal(T, ops.conj(T, col[k]), 1 / @sqrt(ops.abs2(T, col[k])));
            for (col[0..n]) |*x| x.* = ops.mul(T, x.*, t);
            col[k].im = 0;
        }
        return;
    }
    var j: usize = 0;
    while (j < n) {
        const col = v.col(j);
        if (wi.?[j] == 0) {
            const scl = 1 / ops.nrm2(io, T, n, col, 1);
            for (col[0..n]) |*x| x.* *= scl;
            j += 1;
            continue;
        }
        const c2 = v.col(j + 1);
        const scl = 1 / std.math.hypot(ops.nrm2(io, T, n, col, 1), ops.nrm2(io, T, n, c2, 1));
        for (col[0..n], c2[0..n]) |*x, *y| {
            x.* *= scl;
            y.* *= scl;
        }
        const k = argmaxAbs2(T, col, c2, n);
        // `DLARTG(f, g)`: c f + s g = r, -s f + c g = 0; then `DROT`
        const f = col[k];
        const g = c2[k];
        const d = std.math.hypot(f, g);
        const r = if (f < 0) -d else d;
        const cs: T = if (f == 0) 0 else @abs(f) / d;
        const sn: T = if (f == 0) (if (g < 0) -1 else 1) else g / r;
        for (col[0..n], c2[0..n]) |*x, *y| {
            const a = x.*;
            const b = y.*;
            x.* = cs * a + sn * b;
            y.* = cs * b - sn * a;
        }
        c2[k] = 0;
        j += 2;
    }
}

/// `IDAMAX` of `|x(k)|²` (with `y`: `x(k)² + y(k)²`): the first largest.
fn argmaxAbs2(comptime T: type, x: [*]const T, y: ?[*]const T, n: usize) usize {
    var best: usize = 0;
    var m: Float(T) = -1;
    for (0..n) |k| {
        const s = if (comptime isComplex(T)) ops.abs2(T, x[k]) else if (y) |yy| x[k] * x[k] + yy[k] * yy[k] else x[k] * x[k];
        if (s > m) {
            m = s;
            best = k;
        }
    }
    return best;
}

// ---- tests -----------------------------------------------------------------------

const testing = std.testing;
const check = @import("../blas/_test/blas1/check.zig");
const blas = @import("../blas/__root__.zig");
const apple_accel = @import("../blas/apple_accel.zig");
const threading = @import("../blas/threading.zig");

const test_types = .{ f32, f64, types.Complex32, types.Complex64 };

/// Faer's `test_self_adjoint_evd` tolerance: `ApproxEq::eps()` (`128 eps`)
/// times `sqrt(8 n)`, against `max(1, ‖A‖max)` per element.
fn tolerance(comptime F: type, n: usize, amax: F) F {
    return 128 * std.math.floatEps(F) * @sqrt(8 * @as(F, @floatFromInt(n))) * @max(amax, 1);
}

/// `‖A V - V Λ‖`, `‖Vᴴ V - I‖` per element within `tolerance`, `w`
/// nondecreasing.
fn expectSym(comptime T: type, al: Allocator, a0: Mat(T), w: []const Float(T), v: Mat(T), what: []const u8) !void {
    const F = Float(T);
    const n = a0.rows;
    for (1..n) |i| if (!(w[i - 1] <= w[i])) {
        check.say("{s} {s} n={d}: w out of order at {d}\n", .{ what, @typeName(T), n, i });
        return error.TestUnexpectedResult;
    };
    const r = try ops.alloc(T, al, n, n);
    defer ops.free(T, al, r);
    for (0..n) |j| for (0..n) |i| {
        r.at(i, j).* = ops.mulReal(T, v.get(i, j), w[j]);
    };
    try ops.gemm(al, check.io, T, .n, .n, ops.neg(T, ops.one(T)), a0, v, ops.one(T), r);
    const tol = tolerance(F, n, ops.normMax(T, a0));
    if (!(ops.normMax(T, r) <= tol)) {
        check.say("{s} {s} n={d}: A V - V L = {e}, tol {e}\n", .{ what, @typeName(T), n, ops.normMax(T, r), tol });
        return error.TestUnexpectedResult;
    }
    try ops.gemm(al, check.io, T, .c, .n, ops.one(T), v, v, ops.zero(T), r);
    for (0..n) |i| r.at(i, i).* = ops.sub(T, r.get(i, i), ops.one(T));
    if (!(ops.normMax(T, r) <= tolerance(F, n, 1))) {
        check.say("{s} {s} n={d}: Vh V - I = {e}\n", .{ what, @typeName(T), n, ops.normMax(T, r) });
        return error.TestUnexpectedResult;
    }
}

/// A random Hermitian `A + Aᴴ`.
fn hermitian(comptime T: type, a: Mat(T), seed: u64) void {
    const n = a.rows;
    check.fill(T, a.ptr[0 .. a.ld * n], seed);
    for (0..n) |j| for (0..j + 1) |i| {
        const s = ops.add(T, a.get(i, j), ops.conj(T, a.get(j, i)));
        a.at(i, j).* = s;
        a.at(j, i).* = ops.conj(T, s);
    };
}

const Routine = enum { syevd, syev };

/// The routine on a copy of `a0` (padded `lda`), both triangles hidden
/// behind NaN but the one `uplo` names; the vectors into `v`.
fn decompose(comptime T: type, al: Allocator, io: Io, r: Routine, jobz: Jobz, uplo: Uplo, a0: Mat(T), w: []Float(T), v: ?Mat(T)) !void {
    const n = a0.rows;
    const lda = n + 1;
    const a = try ops.alloc(T, al, lda, n);
    defer ops.free(T, al, a);
    a.fill(chol.poisoned(T));
    for (0..n) |j| for (0..n) |i| {
        if (if (uplo == .l) i >= j else i <= j) a.at(i, j).* = a0.get(i, j);
    };
    switch (r) {
        .syevd => try syevd(al, io, T, jobz, uplo, @intCast(n), a.ptr, @intCast(lda), w.ptr),
        .syev => try syev(al, io, T, jobz, uplo, @intCast(n), a.ptr, @intCast(lda), w.ptr),
    }
    for (0..n) |j| try testing.expect(chol.isNanT(T, a.get(n, j)));
    if (v) |x| x.copyFrom(a.sub(0, 0, n, n));
}

test "faer evd: syevd and syev, A V = V L with V orthonormal, both triangles, faer's sizes and past 128" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        for ([_]usize{ 1, 2, 4, 10, 15, 20, 50, 130, 150 }) |n| {
            const a0 = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, a0);
            hermitian(T, a0, 7 * n + 1);
            const w = try al.alloc(F, n);
            defer al.free(w);
            const w2 = try al.alloc(F, n);
            defer al.free(w2);
            const v = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, v);
            for ([_]Routine{ .syevd, .syev }) |r| for ([_]Uplo{ .l, .u }) |uplo| {
                try decompose(T, al, check.io, r, .v, uplo, a0, w, v);
                try expectSym(T, al, a0, w, v, @tagName(r));
                try decompose(T, al, check.io, r, .n, uplo, a0, w2, null);
                for (w, w2) |x, y| try testing.expect(@abs(x - y) <= tolerance(F, n, ops.normMax(T, a0)));
            };
        }
    }
}

test "faer evd: special matrices, faer's test_special: zero, ones, identity; rank one, clustered, graded" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        for ([_]usize{ 1, 2, 4, 15, 50, 150 }) |n| for (0..6) |kind| {
            const a0 = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, a0);
            switch (kind) {
                0 => a0.fill(ops.zero(T)),
                1 => a0.fill(ops.one(T)),
                2 => a0.identity(),
                // clustered: I + eps-sized noise
                3 => {
                    hermitian(T, a0, n);
                    for (0..n) |j| for (a0.col(j)[0..n]) |*x| {
                        x.* = ops.mulReal(T, x.*, std.math.floatEps(F));
                    };
                    for (0..n) |i| a0.at(i, i).* = ops.add(T, a0.get(i, i), ops.one(T));
                },
                // graded: D H D, D from 1 down by 10 per two rows
                4 => {
                    hermitian(T, a0, n + 1);
                    for (0..n) |j| for (0..n) |i| {
                        const s = std.math.pow(F, 10, -@as(F, @floatFromInt(((i + j) / 4) % 8)));
                        a0.at(i, j).* = ops.mulReal(T, a0.get(i, j), s);
                    };
                },
                // tridiagonal Wilkinson W+ (a close pair at the top)
                else => {
                    a0.fill(ops.zero(T));
                    for (0..n) |i| {
                        const c = @as(F, @floatFromInt(n / 2)) - @as(F, @floatFromInt(i));
                        a0.at(i, i).* = ops.scalar(T, @abs(c));
                        if (i + 1 < n) {
                            a0.at(i + 1, i).* = ops.one(T);
                            a0.at(i, i + 1).* = ops.one(T);
                        }
                    }
                },
            }
            const w = try al.alloc(F, n);
            defer al.free(w);
            const v = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, v);
            for ([_]Routine{ .syevd, .syev }) |r| {
                decompose(T, al, check.io, r, .v, .l, a0, w, v) catch |e| {
                    check.say("{s} {s} n={d} kind {d}: {t}\n", .{ @tagName(r), @typeName(T), n, kind, e });
                    return e;
                };
                try expectSym(T, al, a0, w, v, @tagName(r));
            }
        };
    }
}

test "faer evd: OpenBLAS's test_laed3_nan.c: dsyevd on I + 10 v vᵀ, nearly degenerate, no NaN" {
    // OpenBLAS `test/test_laed3_nan.c` (71f22601): n = 8, v = (1, 1 + 1e-8,
    // 1, ...), jobz 'V', uplo 'U'; here through `syevd` (n below 128: QR
    // iteration), through divide and conquer at faer's test threshold 8 and
    // at 4, and at n = 200, where `syevd` divides
    const al = testing.allocator;
    for ([_]usize{ 8, 200 }) |n| {
        const a0 = try ops.alloc(f64, al, n, n);
        defer ops.free(f64, al, a0);
        const v = try al.alloc(f64, n);
        defer al.free(v);
        @memset(v, 1);
        v[1] = 1 + 1e-8;
        for (0..n) |j| for (0..n) |i| {
            a0.at(i, j).* = (if (i == j) @as(f64, 1) else 0) + 10 * v[i] * v[j];
        };
        const w = try al.alloc(f64, n);
        defer al.free(w);
        const a = try ops.alloc(f64, al, n, n);
        defer ops.free(f64, al, a);
        a.copyFrom(a0);
        try syevd(al, check.io, f64, .v, .u, @intCast(n), a.ptr, @intCast(n), w.ptr);
        for (w) |x| try testing.expect(!std.math.isNan(x));
        for (a.ptr[0 .. n * n]) |x| try testing.expect(!std.math.isNan(x));
        try expectSym(f64, al, a0, w, a, "laed3");
        for ([_]usize{ 4, 8 }) |thr| {
            a.copyFrom(a0);
            try symmetric(al, check.io, f64, true, false, a, w, thr);
            for (a.ptr[0 .. n * n]) |x| try testing.expect(!std.math.isNan(x));
            try expectSym(f64, al, a0, w, a, "laed3 dc");
        }
    }
}

test "faer evd: NaN and infinity are NoConvergence; scaled past the bounds; 1 x 1 and empty" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex32 }) |T| {
        const F = Float(T);
        for ([_]F{ std.math.nan(F), std.math.inf(F) }) |bad| {
            var a: [9]T = @splat(ops.one(T));
            a[1] = ops.scalar(T, bad);
            var w: [3]F = undefined;
            try testing.expectError(error.NoConvergence, syevd(al, check.io, T, .v, .l, 3, &a, 3, &w));
            // the other triangle is not read
            a[1] = ops.one(T);
            a[3] = ops.scalar(T, bad);
            try syevd(al, check.io, T, .n, .l, 3, &a, 3, &w);
        }
        const n = 20;
        for ([_]F{ std.math.floatMin(F) * 4, std.math.floatMax(F) / 4 }) |sc| {
            const a0 = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, a0);
            hermitian(T, a0, 3);
            var w0: [n]F = undefined;
            var w1: [n]F = undefined;
            const a = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, a);
            a.copyFrom(a0);
            try syevd(al, check.io, T, .n, .l, n, a.ptr, n, &w0);
            for (0..n) |j| for (a0.col(j)[0..n]) |*x| {
                x.* = ops.mulReal(T, x.*, sc / 8);
            };
            try syevd(al, check.io, T, .n, .l, n, a0.ptr, n, &w1);
            const tol = tolerance(F, n, 1) * (sc / 8) * 8;
            for (w0, w1) |x, y| try testing.expect(@abs(x * (sc / 8) - y) <= tol);
        }
        var one = [_]T{ops.scalar(T, -3)};
        var w: [1]F = undefined;
        try syevd(al, check.io, T, .v, .u, 1, &one, 1, &w);
        try testing.expectEqual(@as(F, -3), w[0]);
        try testing.expectEqual(@as(F, 1), ops.re(T, one[0]));
    }
    const fa = testing.failing_allocator;
    var a: [4]f64 = @splat(std.math.nan(f64));
    var w: [2]f64 = @splat(std.math.nan(f64));
    try testing.expectError(error.InvalidUplo, syevd(fa, check.io, f64, .v, .v, -1, &a, 0, &w));
    try testing.expectError(error.InvalidN, syevd(fa, check.io, f64, .v, .l, -1, &a, 0, &w));
    try testing.expectError(error.InvalidLda, syev(fa, check.io, f64, .v, .u, 2, &a, 1, &w));
    try syevd(fa, check.io, f64, .v, .l, 0, &a, 1, &w);
    for (a ++ w) |x| try testing.expect(std.math.isNan(x));
}

test "faer evd: syevd gives one worker's bits at any worker count" {
    if (comptime !threading.hosted) return error.SkipZigTest;
    _ = blas.coreName();
    const al = testing.allocator;
    const cpus = blas.config.cpu_count;
    defer blas.config.cpu_count = cpus;
    defer threading.test_serial_max = null;
    var threaded: Io.Threaded = .init(al, .{ .async_limit = .limited(8) });
    defer threaded.deinit();
    inline for (.{ f64, types.Complex32 }) |T| {
        const n = 260;
        const a0 = try ops.alloc(T, al, n, n);
        defer ops.free(T, al, a0);
        hermitian(T, a0, 5);
        var outs: [3][]T = undefined;
        var ws: [3][]Float(T) = undefined;
        for (0..3) |c| {
            outs[c] = try al.alloc(T, n * n);
            ws[c] = try al.alloc(Float(T), n);
        }
        defer for (0..3) |c| {
            al.free(outs[c]);
            al.free(ws[c]);
        };
        for ([_]u32{ 1, 3, 8 }, 0..) |wk, c| {
            blas.config.cpu_count = wk;
            threading.test_serial_max = if (wk == 1) null else 4096;
            @memcpy(outs[c], a0.ptr[0 .. n * n]);
            try syevd(al, if (wk == 1) Io.failing else threaded.io(), T, .v, .l, n, outs[c].ptr, n, ws[c].ptr);
        }
        for (1..3) |c| {
            try testing.expectEqualSlices(u8, std.mem.sliceAsBytes(outs[0]), std.mem.sliceAsBytes(outs[c]));
            try testing.expectEqualSlices(u8, std.mem.sliceAsBytes(ws[0]), std.mem.sliceAsBytes(ws[c]));
        }
    }
}

/// Column `j` of `geev`'s vectors as complex, a real pair unpacked.
fn vectorOf(comptime T: type, comptime C: type, wi: ?[]const Float(T), v: Mat(T), j: usize, out: []C) void {
    for (out, 0..) |*o, i| {
        if (comptime isComplex(T)) {
            o.* = v.get(i, j);
        } else if (wi.?[j] == 0) {
            o.* = .init(v.get(i, j), 0);
        } else if (wi.?[j] > 0) {
            o.* = .init(v.get(i, j), v.get(i, j + 1));
        } else {
            o.* = .init(v.get(i, j - 1), -v.get(i, j));
        }
    }
}

/// `‖A v - λ v‖` and `‖uᴴ A - λ uᴴ‖` per vector within `tolerance`, each
/// of unit norm with its largest entry real.
fn expectGeev(comptime T: type, al: Allocator, a0: Mat(T), wr: []const Float(T), wi: ?[]const Float(T), w: ?[]const T, vl: Mat(T), vr: Mat(T)) !void {
    const F = Float(T);
    const C = std.math.Complex(F);
    const n = a0.rows;
    const tol = tolerance(F, n, ops.normMax(T, a0));
    const x = try al.alloc(C, n);
    defer al.free(x);
    for (0..n) |j| {
        const lam: C = if (comptime isComplex(T)) w.?[j] else .init(wr[j], wi.?[j]);
        for ([_]Mat(T){ vr, vl }, 0..) |v, side| {
            vectorOf(T, C, wi, v, j, x);
            // unit norm, largest entry real
            var nrm: F = 0;
            var big: F = -1;
            var kbig: usize = 0;
            for (x, 0..) |e, k| {
                const e2 = e.re * e.re + e.im * e.im;
                nrm += e2;
                if (e2 > big) {
                    big = e2;
                    kbig = k;
                }
            }
            try testing.expect(@abs(@sqrt(nrm) - 1) <= tolerance(F, n, 1));
            try testing.expect(@abs(x[kbig].im) <= tolerance(F, n, 1));
            // right: (A x)(i) - λ x(i); left: (xᴴ A)(i) - λ conj(x(i))
            var worst: F = 0;
            for (0..n) |i| {
                var s: C = .init(0, 0);
                for (0..n) |k| {
                    const e = if (side == 0) a0.get(i, k) else a0.get(k, i);
                    const aik: C = if (comptime isComplex(T)) e else .init(e, 0);
                    const xk = if (side == 0) x[k] else x[k].conjugate();
                    s = s.add(aik.mul(xk));
                }
                const xi = if (side == 0) x[i] else x[i].conjugate();
                const r = s.sub(lam.mul(xi));
                worst = @max(worst, @sqrt(r.re * r.re + r.im * r.im));
            }
            if (!(worst <= tol)) {
                check.say("geev {s} n={d}: {s} vector {d} residual {e}, tol {e}\n", .{ @typeName(T), n, if (side == 0) "right" else "left", j, worst, tol });
                return error.TestUnexpectedResult;
            }
        }
    }
}

/// `geev` on a copy of `a0` with both vector sets, then values alone, which
/// must agree.
fn geevCase(comptime T: type, al: Allocator, io: Io, a0: Mat(T)) !void {
    const F = Float(T);
    const n = a0.rows;
    const a = try ops.alloc(T, al, n + 1, n);
    defer ops.free(T, al, a);
    const vl = try ops.alloc(T, al, n, n);
    defer ops.free(T, al, vl);
    const vr = try ops.alloc(T, al, n, n);
    defer ops.free(T, al, vr);
    const ev = try al.alloc(F, 4 * n);
    defer al.free(ev);
    const cw = try al.alloc(T, 2 * n);
    defer al.free(cw);
    for (0..2) |pass| {
        a.fill(chol.poisoned(T));
        a.sub(0, 0, n, n).copyFrom(a0);
        const jv: Jobz = if (pass == 0) .v else .n;
        const off = pass * 2 * n;
        const w: Eigenvalues(T) = if (comptime isComplex(T)) cw.ptr + pass * n else .{ .re = ev.ptr + off, .im = ev.ptr + off + n };
        try geev(al, io, T, jv, jv, @intCast(n), a.ptr, @intCast(n + 1), w, vl.ptr, @intCast(n), vr.ptr, @intCast(n), null);
        for (0..n) |j| try testing.expect(chol.isNanT(T, a.get(n, j)));
    }
    const etol = tolerance(F, n, ops.normMax(T, a0));
    if (comptime isComplex(T)) {
        try expectGeev(T, al, a0, &.{}, null, cw[0..n], vl, vr);
        for (cw[0..n], cw[n..]) |x, y| try testing.expect(ops.abs(T, ops.sub(T, x, y)) <= etol);
    } else {
        try expectGeev(T, al, a0, ev[0..n], ev[n .. 2 * n], null, vl, vr);
        for (ev[0 .. 2 * n], ev[2 * n ..]) |x, y| try testing.expect(@abs(x - y) <= etol);
        // pairs: positive imaginary part first
        var j: usize = 0;
        while (j < n) : (j += 1) if (ev[n + j] != 0) {
            try testing.expect(ev[n + j] > 0 and ev[n + j + 1] == -ev[n + j] and ev[j] == ev[j + 1]);
            j += 1;
        };
    }
}

test "faer evd: geev, A v = λ v and uᴴ A = λ uᴴ, unit vectors, faer's sizes" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        for ([_]usize{ 1, 2, 3, 4, 10, 15, 20, 50, 100 }) |n| {
            const a0 = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, a0);
            check.fill(T, a0.ptr[0 .. n * n], 3 * n + 2);
            geevCase(T, al, check.io, a0) catch |e| {
                check.say("geev {s} n={d}\n", .{ @typeName(T), n });
                return e;
            };
        }
    }
}

test "faer evd: geev in the extended types, real and complex (T124)" {
    const al = testing.allocator;
    inline for (check.x_real ++ check.x_complex ++ check.f80_unslotted ++ .{ f128, std.math.Complex(f128) }) |T| {
        for ([_]usize{ 1, 2, 3, 10 }) |n| {
            const a0 = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, a0);
            check.fill(T, a0.ptr[0 .. n * n], 3 * n + 2);
            geevCase(T, al, check.io, a0) catch |e| {
                check.say("geev {s} n={d}\n", .{ @typeName(T), n });
                return e;
            };
        }
    }
}

test "faer evd: the complex square root in every float type agrees with f64's (T124)" {
    const C64 = std.math.Complex(f64);
    const cases = [_]C64{ .init(5, 3), .init(-5, 3), .init(-5, -3), .init(0, 2), .init(-4, 0), .init(1e300, 1e300), .init(-1e-300, 1e-310) };
    inline for (.{ f80, f128 }) |F| {
        const C = std.math.Complex(F);
        for (cases) |z| {
            const want = std.math.complex.sqrt(z);
            const got = ops.sqrt(C, .init(z.re, z.im));
            try testing.expectApproxEqRel(want.re, @as(f64, @floatCast(got.re)), 4 * std.math.floatEps(f64));
            try testing.expectApproxEqRel(want.im, @as(f64, @floatCast(got.im)), 4 * std.math.floatEps(f64));
        }
        const inf = std.math.inf(F);
        const r = ops.sqrt(C, .init(-inf, 1));
        try testing.expect(r.re == 0 and r.im == inf);
        try testing.expect(std.math.isNan(ops.sqrt(C, .init(std.math.nan(F), 1)).re));
    }
}

test "faer evd: geev on special matrices: zero, identity, Jordan, rotation pairs, companion, triangular" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        for ([_]usize{ 1, 2, 5, 12 }) |n| for (0..6) |kind| {
            // a companion's roots are ill-conditioned past n = 5 or so; a
            // Jordan block's vectors are only as good as eps^(1/n)
            if ((kind == 4 and n > 5) or (kind == 2 and n > 2)) continue;
            const a0 = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, a0);
            a0.fill(ops.zero(T));
            switch (kind) {
                0 => {},
                1 => a0.identity(),
                2 => for (0..n) |i| {
                    a0.at(i, i).* = ops.scalar(T, 2);
                    if (i + 1 < n) a0.at(i, i + 1).* = ops.one(T);
                },
                // rotations: eigenvalues ±i, and 1 for odd n
                3 => {
                    var i: usize = 0;
                    while (i + 1 < n) : (i += 2) {
                        a0.at(i + 1, i).* = ops.one(T);
                        a0.at(i, i + 1).* = ops.scalar(T, -1);
                    }
                    if (n % 2 == 1) a0.at(n - 1, n - 1).* = ops.one(T);
                },
                // the companion matrix of (x - 1)(x - 2)...(x - n)
                4 => {
                    var c: [13]F = @splat(0);
                    c[0] = 1;
                    for (1..n + 1) |r| {
                        var k = r;
                        while (k > 0) : (k -= 1) c[k] = c[k] - @as(F, @floatFromInt(r)) * c[k - 1];
                    }
                    for (0..n) |j| a0.at(0, j).* = ops.scalar(T, -c[j + 1]);
                    for (1..n) |i| a0.at(i, i - 1).* = ops.one(T);
                },
                else => {
                    check.fill(T, a0.ptr[0 .. n * n], n);
                    for (0..n) |j| for (j + 1..n) |i| {
                        a0.at(i, j).* = ops.zero(T);
                    };
                },
            }
            geevCase(T, al, check.io, a0) catch |e| {
                check.say("geev {s} n={d} kind {d}\n", .{ @typeName(T), n, kind });
                return e;
            };
        };
    }
}

test "faer evd: geev's eigenvectors that overflow unscaled (xTREVC3's scaling)" {
    // upper triangular, T(i, i) = i h, ones above: back substitution for
    // the last (first, for left vectors) eigenvalue grows like
    // 1 / (k! h^k), past the largest float by about 2^500 for f64 (2^60
    // for f32), so only a scaled solve stays finite
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        const n = 20;
        const e = (3 * std.math.floatExponentMax(F)) / (2 * (n - 1)) + 1;
        const h = std.math.ldexp(@as(F, 1), -e);
        const a0 = try ops.alloc(T, al, n, n);
        defer ops.free(T, al, a0);
        a0.fill(ops.zero(T));
        for (0..n) |i| {
            a0.at(i, i).* = ops.scalar(T, @as(F, @floatFromInt(i)) * h);
            if (i + 1 < n) a0.at(i, i + 1).* = ops.one(T);
        }
        geevCase(T, al, check.io, a0) catch |err| {
            check.say("geev overflow {s}\n", .{@typeName(T)});
            return err;
        };
    }
}

test "faer evd: geev's NaN is InvalidA, infinity NoConvergence; scaled past the bounds; arguments" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex32 }) |T| {
        const F = Float(T);
        var a: [9]T = undefined;
        check.fill(T, &a, 1);
        var e: [6]F = undefined;
        var ce: [3]T = undefined;
        const w: Eigenvalues(T) = if (comptime isComplex(T)) &ce else .{ .re = &e, .im = e[3..].ptr };
        var v: [9]T = undefined;
        a[4] = ops.scalar(T, std.math.nan(F));
        try testing.expectError(error.InvalidA, geev(al, check.io, T, .v, .n, 3, &a, 3, w, &v, 3, &v, 1, null));
        a[4] = ops.scalar(T, std.math.inf(F));
        try testing.expectError(error.NoConvergence, geev(al, check.io, T, .n, .n, 3, &a, 3, w, &v, 1, &v, 1, null));
        // scaled past SMLNUM and BIGNUM: the unscaled matrix's vectors
        // (residuals on it, as they would overflow on the scaled one) and
        // its values times the scale
        // powers of 2, so the scaled matrix is exact: normal, past SMLNUM;
        // and below the overflow threshold, past BIGNUM
        const exps = .{ std.math.floatExponentMin(F) + 64, std.math.floatExponentMax(F) - 8 };
        inline for (exps) |ex| {
            const sc = 8 * std.math.ldexp(@as(F, 1), ex);
            const n = 12;
            const a0 = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, a0);
            check.fill(T, a0.ptr[0 .. n * n], 9);
            const as = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, as);
            for (0..n) |j| for (0..n) |i| {
                as.at(i, j).* = ops.mulReal(T, a0.get(i, j), sc / 8);
            };
            const vl = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, vl);
            const vr = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, vr);
            var ev: [2 * n]F = undefined;
            var cw: [n]T = undefined;
            const ws: Eigenvalues(T) = if (comptime isComplex(T)) &cw else .{ .re = &ev, .im = ev[n..].ptr };
            try geev(al, check.io, T, .v, .v, n, as.ptr, n, ws, vl.ptr, n, vr.ptr, n, null);
            for (&ev) |*x| x.* /= sc / 8;
            for (&cw) |*x| x.* = ops.mulReal(T, x.*, 1 / (sc / 8));
            if (comptime isComplex(T)) {
                try expectGeev(T, al, a0, &.{}, null, &cw, vl, vr);
            } else {
                try expectGeev(T, al, a0, ev[0..n], ev[n..], null, vl, vr);
            }
        }
    }
    const fa = testing.failing_allocator;
    var a: [4]f64 = @splat(std.math.nan(f64));
    var wr: [2]f64 = @splat(std.math.nan(f64));
    var wi: [2]f64 = @splat(std.math.nan(f64));
    const w: Eigenvalues(f64) = .{ .re = &wr, .im = &wi };
    try testing.expectError(error.InvalidN, geev(fa, check.io, f64, .v, .v, -1, &a, 0, w, &a, 0, &a, 0, null));
    try testing.expectError(error.InvalidLda, geev(fa, check.io, f64, .v, .v, 2, &a, 1, w, &a, 0, &a, 0, null));
    try testing.expectError(error.InvalidLdvl, geev(fa, check.io, f64, .v, .v, 2, &a, 2, w, &a, 1, &a, 0, null));
    try testing.expectError(error.InvalidLdvl, geev(fa, check.io, f64, .n, .v, 2, &a, 2, w, &a, 0, &a, 2, null));
    try testing.expectError(error.InvalidLdvr, geev(fa, check.io, f64, .n, .v, 2, &a, 2, w, &a, 1, &a, 1, null));
    try geev(fa, check.io, f64, .v, .v, 0, &a, 1, w, &a, 1, &a, 1, null);
    for (a ++ wr ++ wi) |x| try testing.expect(std.math.isNan(x));
}

test "faer evd: geev's INFO > 0: w[i..n] converged and scaled back, the count in unconverged" {
    const al = testing.allocator;
    defer ops.test_iter_limit = null;
    inline for (test_types) |T| {
        const F = Float(T);
        const n = 12;
        var a0: [n * n]T = undefined;
        check.fill(T, &a0, 31);
        var v: [1]T = undefined;
        // the true eigenvalues, then capped runs: unscaled and past SMLNUM
        var ev: [2 * n]F = undefined;
        var cw: [n]T = undefined;
        var a = a0;
        const w0: Eigenvalues(T) = if (comptime isComplex(T)) &cw else .{ .re = &ev, .im = ev[n..].ptr };
        try geev(al, check.io, T, .n, .n, n, &a, n, w0, &v, 1, &v, 1, null);
        for ([_]F{ 1, std.math.ldexp(@as(F, 1), std.math.floatExponentMin(F) + 8) }) |sc| {
            ops.test_iter_limit = 10;
            defer ops.test_iter_limit = null;
            for (&a, a0) |*x, y| x.* = ops.mulReal(T, y, sc);
            var ev1: [2 * n]F = undefined;
            var cw1: [n]T = undefined;
            const w1: Eigenvalues(T) = if (comptime isComplex(T)) &cw1 else .{ .re = &ev1, .im = ev1[n..].ptr };
            var info: usize = 0;
            var vr: [n * n]T = undefined;
            try testing.expectError(error.NoConvergence, geev(al, check.io, T, .n, .v, n, &a, n, w1, &v, 1, &vr, n, &info));
            try testing.expect(info > 0 and info < n);
            const tol = tolerance(F, n, ops.normMax(T, .{ .ptr = &a0, .rows = n, .cols = n, .ld = n }));
            for (info..n) |i| {
                var best = std.math.inf(F);
                for (0..n) |j| {
                    const d = if (comptime isComplex(T)) ops.abs(T, ops.sub(T, ops.mulReal(T, cw1[i], 1 / sc), cw[j])) else std.math.hypot(ev1[i] / sc - ev[j], ev1[n + i] / sc - ev[n + j]);
                    best = @min(best, d);
                }
                try testing.expect(best <= tol);
            }
        }
        // an infinity: none converged
        a = a0;
        a[5] = ops.scalar(T, std.math.inf(F));
        var info: usize = 0;
        try testing.expectError(error.NoConvergence, geev(al, check.io, T, .n, .n, n, &a, n, w0, &v, 1, &v, 1, &info));
        try testing.expectEqual(n, info);
    }
}

test "faer evd: geev gives one worker's bits at any worker count" {
    if (comptime !threading.hosted) return error.SkipZigTest;
    _ = blas.coreName();
    // rstd's kernels: one worker's run is one chunk and the others many, and
    // Accelerate's bits follow a part's shape (its own split keeps the chunks)
    const was = apple_accel.routing();
    _ = apple_accel.route(false);
    defer _ = apple_accel.route(was);
    const al = testing.allocator;
    const cpus = blas.config.cpu_count;
    defer blas.config.cpu_count = cpus;
    defer threading.test_serial_max = null;
    var threaded: Io.Threaded = .init(al, .{ .async_limit = .limited(8) });
    defer threaded.deinit();
    inline for (.{ f64, types.Complex32 }) |T| {
        const n = 150;
        const a0 = try ops.alloc(T, al, n, n);
        defer ops.free(T, al, a0);
        check.fill(T, a0.ptr[0 .. n * n], 11);
        var outs: [3][]T = undefined;
        var evs: [3][]T = undefined;
        for (0..3) |c| {
            outs[c] = try al.alloc(T, 3 * n * n);
            evs[c] = try al.alloc(T, 2 * n);
        }
        defer for (0..3) |c| {
            al.free(outs[c]);
            al.free(evs[c]);
        };
        for ([_]u32{ 1, 3, 8 }, 0..) |wk, c| {
            blas.config.cpu_count = wk;
            threading.test_serial_max = if (wk == 1) null else 4096;
            const o = outs[c];
            @memcpy(o[0 .. n * n], a0.ptr[0 .. n * n]);
            const w: Eigenvalues(T) = if (comptime isComplex(T)) evs[c].ptr else .{ .re = evs[c].ptr, .im = evs[c].ptr + n };
            try geev(al, if (wk == 1) Io.failing else threaded.io(), T, .v, .v, n, o.ptr, n, w, o.ptr + n * n, n, o.ptr + 2 * n * n, n, null);
        }
        for (1..3) |c| {
            try testing.expectEqualSlices(u8, std.mem.sliceAsBytes(outs[0]), std.mem.sliceAsBytes(outs[c]));
            try testing.expectEqualSlices(u8, std.mem.sliceAsBytes(evs[0]), std.mem.sliceAsBytes(evs[c]));
        }
    }
}

test "faer evd: geev balances, so D M D⁻¹ keeps M's eigenvalues to M's scale (dgebal.f)" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        const n = 16;
        const m = try ops.alloc(T, al, n, n);
        defer ops.free(T, al, m);
        check.fill(T, m.ptr[0 .. n * n], 21);
        const a = try ops.alloc(T, al, n, n);
        defer ops.free(T, al, a);
        // D = diag(2^(3i)): entries from 2^-45 to 2^45 times M's
        for (0..n) |j| for (0..n) |i| {
            const d = std.math.pow(F, 2, 3 * (@as(F, @floatFromInt(i)) - @as(F, @floatFromInt(j))));
            a.at(i, j).* = ops.mulReal(T, m.get(i, j), d);
        };
        var ev: [4 * n]F = undefined;
        var cw: [2 * n]T = undefined;
        var v: [1]T = undefined;
        for ([_]Mat(T){ m, a }, 0..) |x, p| {
            const w: Eigenvalues(T) = if (comptime isComplex(T)) cw[p * n ..].ptr else .{ .re = ev[2 * p * n ..].ptr, .im = ev[2 * p * n + n ..].ptr };
            try geev(al, check.io, T, .n, .n, n, x.ptr, n, w, &v, 1, &v, 1, null);
        }
        // each of A's within tolerance of one of M's
        const tol = tolerance(F, n, 4);
        for (0..n) |i| {
            var best = std.math.inf(F);
            for (0..n) |j| {
                const d = if (comptime isComplex(T)) ops.abs(T, ops.sub(T, cw[n + i], cw[j])) else std.math.hypot(ev[2 * n + i] - ev[j], ev[3 * n + i] - ev[n + j]);
                best = @min(best, d);
            }
            try testing.expect(best <= tol);
        }
    }
}

test "hseqr: Zᴴ A Z = T after gehrd and orghr, the ilo..ihi block, both jobs, real and complex" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64, f32, types.Complex32 }) |T| {
        const F = Float(T);
        for ([_][3]usize{ .{ 1, 1, 1 }, .{ 2, 1, 2 }, .{ 6, 1, 6 }, .{ 9, 3, 7 }, .{ 30, 4, 27 } }) |c| {
            const n = c[0];
            const lo = c[1] - 1;
            const hi = c[2];
            const a0 = try al.alloc(T, n * n);
            defer al.free(a0);
            check.fill(T, a0, 11 * n + lo);
            for (0..n) |j| for (j + 1..n) |i| {
                if (j < lo or i >= hi) a0[i + j * n] = ops.zero(T);
            };
            const t = try al.dupe(T, a0);
            defer al.free(t);
            const tau = try al.alloc(T, @max(n, 2) - 1);
            defer al.free(tau);
            const ni: isize = @intCast(n);
            try hessenberg.gehrd(al, check.io, T, ni, @intCast(lo + 1), @intCast(hi), t.ptr, ni, tau.ptr);
            const z = try al.dupe(T, t);
            defer al.free(z);
            try hessenberg.orghr(al, check.io, T, ni, @intCast(lo + 1), @intCast(hi), z.ptr, ni, tau.ptr);
            const h2 = try al.dupe(T, t);
            defer al.free(h2);
            const ev = try al.alloc(F, 4 * n);
            defer al.free(ev);
            const cw = try al.alloc(T, 2 * n);
            defer al.free(cw);
            const w: Eigenvalues(T) = if (comptime isComplex(T)) cw.ptr else .{ .re = ev.ptr, .im = ev.ptr + n };
            const w2: Eigenvalues(T) = if (comptime isComplex(T)) cw.ptr + n else .{ .re = ev.ptr + 2 * n, .im = ev.ptr + 3 * n };
            try hseqr(check.io, T, .s, .v, ni, @intCast(lo + 1), @intCast(hi), t.ptr, ni, w, z.ptr, ni, null);
            var dummy: [1]T = undefined;
            try hseqr(check.io, T, .e, .n, ni, @intCast(lo + 1), @intCast(hi), h2.ptr, ni, w2, &dummy, 1, null);
            const nrm = @max(ops.normMax(T, .{ .ptr = a0.ptr, .rows = n, .cols = n, .ld = n }), 1);
            const tol = 128 * std.math.floatEps(F) * @as(F, @floatFromInt(n + 1)) * nrm;
            for (0..n) |j| for (0..n) |i| {
                var s = ops.zero(T);
                for (0..n) |k| for (0..n) |l| {
                    s = ops.add(T, s, ops.mul(T, ops.mul(T, ops.conj(T, z[k + i * n]), a0[k + l * n]), z[l + j * n]));
                };
                // T upper triangular; real: a subdiagonal entry only inside a complex pair's block
                const pair = !isComplex(T) and i == j + 1 and (if (comptime isComplex(T)) false else ev[n + j] > 0);
                const want = if (i <= j or pair) t[i + j * n] else ops.zero(T);
                if (i > j and !pair and !ops.isZero(T, t[i + j * n])) return error.TestUnexpectedResult;
                var g = ops.zero(T);
                for (0..n) |k| g = ops.add(T, g, ops.mul(T, ops.conj(T, z[k + i * n]), z[k + j * n]));
                const err = ops.abs(T, ops.sub(T, s, want));
                const gerr = ops.abs(T, ops.sub(T, g, if (i == j) ops.one(T) else ops.zero(T)));
                if (!(err <= tol) or !(gerr <= tol / nrm)) {
                    check.say("hseqr {s} n={d} ilo={d} ihi={d}: ({d},{d}) T off by {e}, ZᴴZ by {e}\n", .{ @typeName(T), n, lo + 1, hi, i, j, err, gerr });
                    return error.TestUnexpectedResult;
                }
            };
            // the eigenvalues alone agree with the Schur form's
            if (comptime isComplex(T)) {
                for (cw[0..n], cw[n..]) |x, y| try testing.expect(ops.abs(T, ops.sub(T, x, y)) <= tol);
            } else for (0..n) |i| {
                try testing.expect(@abs(ev[i] - ev[2 * n + i]) <= tol and @abs(ev[n + i] - ev[3 * n + i]) <= tol);
            }
        }
    }
    var a: [4]f64 = @splat(0);
    var wr: [2]f64 = undefined;
    var wi: [2]f64 = undefined;
    const w: Eigenvalues(f64) = .{ .re = &wr, .im = &wi };
    try testing.expectError(error.InvalidN, hseqr(check.io, f64, .e, .n, -1, 1, 0, &a, 1, w, &a, 1, null));
    try testing.expectError(error.InvalidIlo, hseqr(check.io, f64, .e, .n, 2, 3, 2, &a, 2, w, &a, 1, null));
    try testing.expectError(error.InvalidIhi, hseqr(check.io, f64, .e, .n, 2, 2, 1, &a, 2, w, &a, 1, null));
    try testing.expectError(error.InvalidLdh, hseqr(check.io, f64, .e, .n, 2, 1, 2, &a, 1, w, &a, 1, null));
    try testing.expectError(error.InvalidLdz, hseqr(check.io, f64, .s, .i, 2, 1, 2, &a, 2, w, &a, 1, null));
    try hseqr(check.io, f64, .e, .n, 0, 1, 0, &a, 1, w, &a, 1, null);
}

/// `‖A v - λ v‖` (right) or `‖vᴴ A - λ vᴴ‖` (left) for column `j` of `v`
/// (a real pair's two columns as one complex vector), and its largest
/// `|re| + |im|`.
fn eigResidual(comptime T: type, n: usize, a: []const T, v: []const T, j: usize, pair: bool, lam: std.math.Complex(Float(T)), right: bool) [2]Float(T) {
    const F = Float(T);
    const C = std.math.Complex(F);
    var worst: F = 0;
    var big: F = 0;
    for (0..n) |k| {
        const x: C = if (comptime isComplex(T)) v[k + j * n] else .init(v[k + j * n], if (pair) v[k + (j + 1) * n] else 0);
        big = @max(big, @abs(x.re) + @abs(x.im));
    }
    for (0..n) |i| {
        var s: C = .init(0, 0);
        for (0..n) |k| {
            const e = if (right) a[i + k * n] else a[k + i * n];
            const aik: C = if (comptime isComplex(T)) e else .init(e, 0);
            const xk_: C = if (comptime isComplex(T)) v[k + j * n] else .init(v[k + j * n], if (pair) v[k + (j + 1) * n] else 0);
            s = s.add(aik.mul(if (right) xk_ else xk_.conjugate()));
        }
        const xi_: C = if (comptime isComplex(T)) v[i + j * n] else .init(v[i + j * n], if (pair) v[i + (j + 1) * n] else 0);
        const r = s.sub(if (right) lam.mul(xi_) else lam.mul(xi_.conjugate()));
        worst = @max(worst, @sqrt(r.re * r.re + r.im * r.im));
    }
    return .{ worst, big };
}

test "trevc3: A's eigenvectors after gehrd, orghr and hseqr, back-transformed and selected, real and complex" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64, f32, types.Complex32 }) |T| {
        const F = Float(T);
        const C = std.math.Complex(F);
        for ([_]usize{ 1, 2, 5, 12, 31 }) |n| {
            const ni: isize = @intCast(n);
            const a0 = try al.alloc(T, n * n);
            defer al.free(a0);
            check.fill(T, a0, 13 * n + 1);
            const t = try al.dupe(T, a0);
            defer al.free(t);
            const tau = try al.alloc(T, @max(n, 2) - 1);
            defer al.free(tau);
            try hessenberg.gehrd(al, check.io, T, ni, 1, ni, t.ptr, ni, tau.ptr);
            const z = try al.dupe(T, t);
            defer al.free(z);
            try hessenberg.orghr(al, check.io, T, ni, 1, ni, z.ptr, ni, tau.ptr);
            const ev = try al.alloc(F, 2 * n);
            defer al.free(ev);
            const cw = try al.alloc(T, n);
            defer al.free(cw);
            const w: Eigenvalues(T) = if (comptime isComplex(T)) cw.ptr else .{ .re = ev.ptr, .im = ev.ptr + n };
            try hseqr(check.io, T, .s, .v, ni, 1, ni, t.ptr, ni, w, z.ptr, ni, null);
            const vl = try al.dupe(T, z);
            defer al.free(vl);
            const vr = try al.dupe(T, z);
            defer al.free(vr);
            var m: isize = 0;
            try trevc.trevc3(al, check.io, T, .b, .b, null, ni, t.ptr, ni, vl.ptr, ni, vr.ptr, ni, ni, &m);
            try testing.expectEqual(ni, m);
            const tol = 256 * std.math.floatEps(F) * @as(F, @floatFromInt(n + 1)) * @max(ops.normMax(T, .{ .ptr = a0.ptr, .rows = n, .cols = n, .ld = n }), 1);
            var j: usize = 0;
            while (j < n) : (j += 1) {
                const pair = !isComplex(T) and ev[n + j] > 0;
                const lam: C = if (comptime isComplex(T)) cw[j] else .init(ev[j], ev[n + j]);
                for ([_]bool{ true, false }) |right| {
                    const r = eigResidual(T, n, a0, if (right) vr else vl, j, pair, lam, right);
                    if (!(r[0] <= tol) or !(@abs(r[1] - 1) <= tol)) {
                        check.say("trevc3 {s} n={d} vector {d} {s}: residual {e}, largest {e}\n", .{ @typeName(T), n, j, if (right) "right" else "left", r[0], r[1] });
                        return error.TestUnexpectedResult;
                    }
                }
                if (pair) j += 1;
            }
            // selected: every other eigenvalue, a pair picked by its second entry
            const sel = try al.alloc(bool, n);
            defer al.free(sel);
            for (sel, 0..) |*s, i| s.* = i % 2 == 1;
            const all = try al.alloc(T, n * n);
            defer al.free(all);
            const some = try al.alloc(T, n * n);
            defer al.free(some);
            try trevc.trevc3(al, check.io, T, .r, .a, null, ni, t.ptr, ni, vl.ptr, ni, all.ptr, ni, ni, &m);
            try trevc.trevc3(al, check.io, T, .r, .s, sel.ptr, ni, t.ptr, ni, vl.ptr, ni, some.ptr, ni, ni, &m);
            var c: usize = 0;
            j = 0;
            while (j < n) : (j += 1) {
                const pair = !isComplex(T) and ev[n + j] > 0;
                if (sel[j]) {
                    for (0..if (pair) 2 else 1) |o| {
                        try testing.expectEqualSlices(T, all[(j + o) * n ..][0..n], some[(c + o) * n ..][0..n]);
                    }
                    c += if (pair) 2 else 1;
                }
                if (pair) {
                    try testing.expect(!sel[j + 1]);
                    j += 1;
                }
            }
            try testing.expectEqual(@as(isize, @intCast(c)), m);
        }
    }
    var a: [4]f64 = @splat(0);
    var m: isize = 0;
    try testing.expectError(error.InvalidN, trevc.trevc3(al, check.io, f64, .r, .a, null, -1, &a, 1, &a, 1, &a, 1, 0, &m));
    try testing.expectError(error.InvalidSelect, trevc.trevc3(al, check.io, f64, .r, .s, null, 2, &a, 2, &a, 1, &a, 2, 2, &m));
    try testing.expectError(error.InvalidLdvr, trevc.trevc3(al, check.io, f64, .r, .a, null, 2, &a, 2, &a, 1, &a, 1, 2, &m));
    try testing.expectError(error.InvalidMm, trevc.trevc3(al, check.io, f64, .r, .a, null, 2, &a, 2, &a, 1, &a, 2, 1, &m));
}

/// `gees`'s test sort: the left half plane first.
fn LeftHalf(comptime T: type) type {
    return struct {
        fn real(re: Float(T), im: Float(T)) bool {
            _ = im;
            return re < 0;
        }
        fn cplx(x: T) bool {
            return x.re < 0;
        }
    };
}

test "gees: A = Z T Zᴴ, unsorted and the left half plane first, tiny scale, NaN; real and complex" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64, f32, types.Complex32 }) |T| {
        const F = Float(T);
        const pick: GeesSelect(T) = if (comptime isComplex(T)) &LeftHalf(T).cplx else &LeftHalf(T).real;
        for ([_]usize{ 1, 2, 6, 17 }) |n| for ([_]bool{ false, true }) |sort| for ([_]F{ 1, 1e-30 }) |sc| {
            const ni: isize = @intCast(n);
            const a0 = try al.alloc(T, n * n);
            defer al.free(a0);
            check.fill(T, a0, 17 * n + @intFromBool(sort));
            for (a0) |*e| e.* = ops.mulReal(T, e.*, sc);
            const t = try al.dupe(T, a0);
            defer al.free(t);
            const z = try al.alloc(T, n * n);
            defer al.free(z);
            const ev = try al.alloc(F, 2 * n);
            defer al.free(ev);
            const cw = try al.alloc(T, n);
            defer al.free(cw);
            const w: Eigenvalues(T) = if (comptime isComplex(T)) cw.ptr else .{ .re = ev.ptr, .im = ev.ptr + n };
            var sdim: isize = -1;
            try gees(al, check.io, T, .v, if (sort) pick else null, ni, t.ptr, ni, &sdim, w, z.ptr, ni, null);
            const nrm = @max(ops.normMax(T, .{ .ptr = a0.ptr, .rows = n, .cols = n, .ld = n }), std.math.floatMin(F));
            const tol = 256 * std.math.floatEps(F) * @as(F, @floatFromInt(n + 1)) * nrm;
            for (0..n) |j| for (0..n) |i| {
                var s = ops.zero(T);
                for (0..n) |k| for (0..n) |l| {
                    s = ops.add(T, s, ops.mul(T, ops.mul(T, z[i + k * n], t[k + l * n]), ops.conj(T, z[j + l * n])));
                };
                const pair = !isComplex(T) and i == j + 1 and (if (comptime isComplex(T)) false else ev[n + j] > 0);
                if (i > j and !pair and !ops.isZero(T, t[i + j * n])) return error.TestUnexpectedResult;
                const err = ops.abs(T, ops.sub(T, s, a0[i + j * n]));
                if (!(err <= tol)) {
                    check.say("gees {s} n={d} sort={} scale={e}: ({d},{d}) Z T Zᴴ off by {e}\n", .{ @typeName(T), n, sort, sc, i, j, err });
                    return error.TestUnexpectedResult;
                }
            };
            if (!sort) {
                try testing.expectEqual(@as(isize, 0), sdim);
                continue;
            }
            var neg: isize = 0;
            for (0..n) |i| {
                const re = if (comptime isComplex(T)) cw[i].re else ev[i];
                if (re < 0) neg += 1;
                try testing.expectEqual(re < 0, @as(isize, @intCast(i)) < sdim);
            }
            try testing.expectEqual(neg, sdim);
        };
    }
    var a: [4]f64 = .{ 1, std.math.nan(f64), 0, 1 };
    var z: [4]f64 = undefined;
    var wr: [2]f64 = undefined;
    var wi: [2]f64 = undefined;
    var sdim: isize = 0;
    const w: Eigenvalues(f64) = .{ .re = &wr, .im = &wi };
    try testing.expectError(error.InvalidA, gees(al, check.io, f64, .v, null, 2, &a, 2, &sdim, w, &z, 2, null));
    try testing.expectError(error.InvalidN, gees(al, check.io, f64, .v, null, -1, &a, 1, &sdim, w, &z, 1, null));
    try testing.expectError(error.InvalidLda, gees(al, check.io, f64, .n, null, 2, &a, 1, &sdim, w, &z, 1, null));
    try testing.expectError(error.InvalidLdvs, gees(al, check.io, f64, .v, null, 2, &a, 2, &sdim, w, &z, 1, null));
    try gees(al, check.io, f64, .n, null, 0, &a, 1, &sdim, w, &z, 1, null);
}

test "syevd, syev: graded large at the bottom right given as .u, and at the top left as .l, keep the small eigenvalues (T125)" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64 }) |T| {
        const F = Float(T);
        const W = if (comptime isComplex(T)) std.math.Complex(f128) else f128;
        for ([_]usize{ 4, 6, 10 }) |n| for ([_]bool{ true, false }) |increasing| for ([_]bool{ false, true }) |use_qr| {
            // M1's E2: D B D, B Hermitian with a dominant diagonal, D geometric
            const ratio: f64 = if (n == 4) 1e-4 else 1e-2;
            var prng = std.Random.DefaultPrng.init(n * 7 + @intFromBool(increasing));
            const r = prng.random();
            const a = try al.alloc(T, n * n);
            defer al.free(a);
            const a128 = try al.alloc(W, n * n);
            defer al.free(a128);
            for (0..n) |j| for (0..j + 1) |i| {
                const di = std.math.pow(f64, ratio, @floatFromInt(if (increasing) n - 1 - i else i));
                const dj = std.math.pow(f64, ratio, @floatFromInt(if (increasing) n - 1 - j else j));
                const br = if (i == j) 2 * @as(f64, @floatFromInt(n)) + r.float(f64) else r.float(f64) - 0.5;
                const bi = if (i == j or !isComplex(T)) 0 else r.float(f64) - 0.5;
                const re = di * br * dj;
                const im = di * bi * dj;
                if (comptime isComplex(T)) {
                    a[i + j * n] = .init(re, im);
                    a[j + i * n] = .init(re, -im);
                    a128[i + j * n] = .init(re, im);
                    a128[j + i * n] = .init(re, -im);
                } else {
                    a[i + j * n] = re;
                    a[j + i * n] = re;
                    a128[i + j * n] = re;
                    a128[j + i * n] = re;
                }
            };
            const uplo: Uplo = if (increasing) .u else .l;
            const w = try al.alloc(F, n);
            defer al.free(w);
            const w128 = try al.alloc(f128, n);
            defer al.free(w128);
            const ni: isize = @intCast(n);
            if (use_qr) try syev(al, check.io, T, .n, uplo, ni, a.ptr, ni, w.ptr) else try syevd(al, check.io, T, .n, uplo, ni, a.ptr, ni, w.ptr);
            try syevd(al, check.io, W, .n, uplo, ni, a128.ptr, ni, w128.ptr);
            // graded small eigenvalues hold to a few hundred eps (matlab88's
            // bottom-up reduction reaches 1.7e-13 here, M1's E2); reduced
            // from the wrong end they were off by 0.3 to 3e8
            for (w, w128, 0..) |x, ref, i| {
                const err: f128 = @abs(@as(f128, x) - ref) / @abs(ref);
                if (!(err <= 1024 * @as(f128, @floatFromInt(n)) * std.math.floatEps(F))) {
                    check.say("{s} {s} n={d} increasing={}: eigenvalue {d} relative error {e}\n", .{ if (use_qr) "syev" else "syevd", @typeName(T), n, increasing, i, @as(f64, @floatCast(err)) });
                    return error.TestUnexpectedResult;
                }
            }
        };
    }
}
