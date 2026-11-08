//! `?SYEVD`, `?HEEVD`, `?SYEV`, `?HEEV`, `?GEEV`, `?HSEQR` and `?GEES`:
//! Reference-LAPACK's arguments and scaling (`lapack-netlib/SRC/dsyevd.f`,
//! `dsyev.f`, `dgeev.f`, `dgees.f` and their forms, OpenBLAS v0.3.34-317) over
//! faer v0.24.4's methods (<https://codeberg.org/sarah-quinones/faer>, tag
//! `faer-v0.24.4`, MIT, `faer/src/linalg/evd/mod.rs`). Departures: design.md's
//! LAPACK section.

const std = @import("std");
const ops = @import("ops.zig");
const blas_calls = @import("blas_calls.zig");
const auxiliary = @import("auxiliary.zig");
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
    comptime ops.requireFloat(T, "syevd");
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
    if (sigma) |s| auxiliary.lascl(T, false, 1, s, trid);
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
            if (i > 0) x = ops.normalized(T, ops.mul(T, trid.get(i, i - 1), x));
            for (0..n) |j| um.at(i, j).* = ops.mulReal(T, x, ur.get(i, j));
        }
    }
    if (n > 1) try hh.applySequence(al, io, T, trid.sub(1, 0, n - 1, n - 1), h.sub(0, 0, bs, n - 1), false, um.sub(1, 0, n - 1, n));
    @memcpy(s, diag);
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
    comptime ops.requireFloat(T, "geev");
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
    const scaled = auxiliary.scaleIn(T, am) catch |e| switch (e) {
        error.NanInA => return error.InvalidA,
        error.InfInA => {
            if (unconverged) |c| c.* = un;
            return error.NoConvergence;
        },
    };
    const anrm = scaled.anrm;
    const to = scaled.to;
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
        for (ws) |x| auxiliary.lascl(T, false, t, anrm, .{ .ptr = x + info, .rows = un - info, .cols = 1, .ld = un });
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
    comptime ops.requireFloat(T, "hseqr");
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
    comptime ops.requireFloat(T, "gees");
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
    const scaled = auxiliary.scaleIn(T, am) catch |e| switch (e) {
        error.NanInA => return error.InvalidA,
        error.InfInA => {
            if (unconverged) |c| c.* = un;
            return error.NoConvergence;
        },
    };
    const anrm = scaled.anrm;
    const to = scaled.to;
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
    hessenberg.gehrd(allocator, io, T, n, ilo, ihi, a, lda, tau.ptr) catch |e| return ops.narrow(ops.Error, e);
    if (jobvs == .v) {
        const vm: Mat(T) = .{ .ptr = vs, .rows = un, .cols = un, .ld = @intCast(ldvs) };
        vm.copyFrom(am);
        hessenberg.orghr(allocator, io, T, n, ilo, ihi, vs, ldvs, tau.ptr) catch |e| return ops.narrow(ops.Error, e);
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
        auxiliary.lascl(T, false, t, anrm, am);
        const ws = if (comptime isComplex(T)) [_][*]T{w} else [_][*]T{ w.re, w.im };
        for (ws) |x| auxiliary.lascl(T, false, t, anrm, .{ .ptr = x, .rows = un, .cols = 1, .ld = un });
    }
}

/// Whether `select` picks eigenvalue `i`, by its own `(re, im)` for real,
/// as `dgees.f` asks each.
fn picked(comptime T: type, pick: GeesSelect(T), w: Eigenvalues(T), i: usize) bool {
    return if (comptime isComplex(T)) pick(w[i]) else pick(w.re[i], w.im[i]);
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
        try blas_calls.gemm(al, io, T, .n, .n, ops.one(T), zm, x, ops.zero(T), u);
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
            const scl = 1 / blas_calls.nrm2(io, T, n, col, 1);
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
            const scl = 1 / blas_calls.nrm2(io, T, n, col, 1);
            for (col[0..n]) |*x| x.* *= scl;
            j += 1;
            continue;
        }
        const c2 = v.col(j + 1);
        const scl = 1 / std.math.hypot(blas_calls.nrm2(io, T, n, col, 1), blas_calls.nrm2(io, T, n, c2, 1));
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
