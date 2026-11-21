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
const stebr_dc = @import("stebr.zig");
const sytrd = @import("sytrd.zig");
const hessenberg = @import("hessenberg.zig");
const real_schur = @import("real_schur.zig");
const complex_schur = @import("complex_schur.zig");
const balance = @import("balance.zig");
const schur = @import("schur.zig");
const trevc = @import("trevc.zig");
const util = @import("../blas/util.zig");
const types = @import("../blas/types.zig");
const interface = @import("interface.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = ops.Mat;
const Float = ops.Float;
const isComplex = ops.isComplex;
const Eigenvalues = schur.Eigenvalues;
const Uplo = types.Uplo;

pub const VectorJob = interface.VectorJob;

/// What `syevd`, `heevd`, `syev` and `heev` return: a reversed `uplo`, a
/// slice too short (`a`, or `w` under `n`) or a leading dimension too
/// small; `NoConvergence`, also for a NaN or infinity in `A`.
pub const SyevError = interface.TriangleError || error{NoConvergence} || Allocator.Error || Io.Cancelable;

/// `SelfAdjointEvdParams::auto`'s `recursion_threshold`: below it, divide
/// and conquer is QR iteration.
pub const dc_threshold = 128;

/// `?SYEVD`, `?HEEVD`: the eigenvalues of the `n x n` symmetric (Hermitian)
/// `A`, from its `uplo` triangle, into `w` (`n`, nondecreasing); with
/// `jobz` `.compute` the orthonormal eigenvectors over `A`. Divide and conquer.
pub fn syevd(allocator: Allocator, io: Io, comptime T: type, jobz: VectorJob, uplo: Uplo, n: usize, a: []T, lda: usize, w: []Float(T)) SyevError!void {
    return driver(allocator, io, T, jobz, uplo, n, a, lda, w, dc_threshold);
}

/// The order from which the values alone come from `stebr.zig`'s divide and
/// conquer rather than QR iteration (both drivers, real and complex). On
/// Zen 3, one thread, f64, two runs (T186, 2026-11-19), `stebr` took of
/// `qrAlgorithm`'s time 0.90 to 0.94 at n = 26 on random and 1-2-1
/// matrices but 1.12 clustered, 0.72 to 0.91 at 40 on all three, 0.42 to
/// 0.48 at 100. Below 26 `stebr` is `qrAlgorithm`.
pub const stebr_crossover = 40;

/// `?SYEV`, `?HEEV`: as `syevd`, the tridiagonal's vectors by QR iteration.
pub fn syev(allocator: Allocator, io: Io, comptime T: type, jobz: VectorJob, uplo: Uplo, n: usize, a: []T, lda: usize, w: []Float(T)) SyevError!void {
    return driver(allocator, io, T, jobz, uplo, n, a, lda, w, std.math.maxInt(usize));
}

fn driver(al: Allocator, io: Io, comptime T: type, jobz: VectorJob, uplo: Uplo, n: usize, a: []T, lda: usize, w: []Float(T), threshold: usize) SyevError!void {
    comptime ops.requireFloat(T, "syevd");
    try interface.holdsTriangle(uplo, n, a.len, lda);
    if (w.len < n) return error.SliceTooShort;
    if (n == 0) return;
    const am: Mat(T) = .{ .ptr = a.ptr, .rows = n, .cols = n, .ld = lda };
    return symmetric(al, io, T, jobz == .compute, uplo.isLower(), am, w[0..n], threshold);
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

/// The eigenvalues alone of the self-adjoint `trid` (lower triangle,
/// overwritten) into `s`, nondecreasing: LAPACK's blocked `?SYTRD`/`?HETRD`
/// (`sytrd.zig`; faer's unblocked reduction ran at 0.13 to 0.77 of
/// OpenBLAS's `syevd` rate, issue #552), then `stebr` or QR iteration.
fn values(al: Allocator, io: Io, comptime T: type, trid: Mat(T), s: []Float(T)) (ops.Error || error{NoConvergence})!void {
    const F = Float(T);
    const n = trid.rows;
    const real = try al.alloc(F, 2 * n);
    defer al.free(real);
    const diag = real[0..n];
    const offdiag = real[n..];
    @memset(offdiag, 0);
    const tau = try al.alloc(T, @max(n, 2) - 1);
    defer al.free(tau);
    try sytrd.sytrd(al, io, T, trid, diag, offdiag[0 .. n - 1], tau);
    if (n < stebr_crossover) {
        try tevd.qrAlgorithm(F, diag, offdiag, null);
    } else try stebr_dc.stebr(al, F, diag, offdiag[0 .. n - 1]);
    @memcpy(s, diag);
}

/// `self_adjoint_evd`: `trid`'s lower triangle (overwritten) to the
/// eigenvalues in `s`, nondecreasing, and the vectors in `u` if given.
pub fn selfAdjointEvd(al: Allocator, io: Io, comptime T: type, trid: Mat(T), s: []Float(T), u: ?Mat(T), threshold: usize) (ops.Error || error{NoConvergence})!void {
    const F = Float(T);
    const n = trid.rows;
    if (n == 0) return;
    if (u == null) return values(al, io, T, trid, s);
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
    const um = u.?;
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

/// What `geev` returns: a slice too short or a leading dimension too small
/// (`vl`, `vr` checked only when computed), `NanInA` for a NaN in `A` (the
/// reference's `INFO = -4`, `dgeev.f:353`), `NoConvergence` for its
/// `INFO > 0` and for an infinity in `A`.
pub const GeevError = interface.ShapeError || interface.NanError || error{NoConvergence} || Allocator.Error || Io.Cancelable;

/// `?GEEV`: the eigenvalues of the general `n x n` `A` into `w`, and with
/// `jobvl`, `jobvr` `.compute` its left (`uᴴ A = λ uᴴ`) and right (`A v = λ v`)
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
pub fn geev(allocator: Allocator, io: Io, comptime T: type, jobvl: VectorJob, jobvr: VectorJob, n: usize, a: []T, lda: usize, w: Eigenvalues(T), vl: []T, ldvl: usize, vr: []T, ldvr: usize, unconverged: ?*usize) GeevError!void {
    comptime ops.requireFloat(T, "geev");
    const F = Float(T);
    try util.holds(a.len, n, n, lda);
    try schur.eigenvaluesHold(T, w, n);
    if (jobvl == .compute) try util.holds(vl.len, n, n, ldvl);
    if (jobvr == .compute) try util.holds(vr.len, n, n, ldvr);
    if (n == 0) return;
    const un = n;
    const am: Mat(T) = .{ .ptr = a.ptr, .rows = un, .cols = un, .ld = lda };
    // `dgeev.f:341-359`: NaN, then the scale into [SMLNUM, BIGNUM]
    const scaled = auxiliary.scaleIn(T, am) catch |e| switch (e) {
        error.NanInA => return error.NanInA,
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
    const ul: ?Mat(T) = if (jobvl == .compute) Mat(T){ .ptr = vl.ptr, .rows = un, .cols = un, .ld = ldvl } else null;
    const ur: ?Mat(T) = if (jobvr == .compute) Mat(T){ .ptr = vr.ptr, .rows = un, .cols = un, .ld = ldvr } else null;
    // `dgeev.f:398-404`: after a failure only `w[info..n]` is scaled back
    var info: usize = 0;
    const res = if (comptime isComplex(T))
        general(allocator, io, T, am, w[0..un], null, ul, ur, &info)
    else
        general(allocator, io, T, am, w.re[0..un], w.im[0..un], ul, ur, &info);
    res catch |e| if (e != error.NoConvergence) return e;
    if (to) |t| {
        const ws = if (comptime isComplex(T)) [_][*]T{w.ptr} else [_][*]T{ w.re.ptr, w.im.ptr };
        for (ws) |x| auxiliary.lascl(T, false, t, anrm, .{ .ptr = x + info, .rows = un - info, .cols = 1, .ld = un });
    }
    res catch |e| {
        if (unconverged) |c| c.* = info;
        return e;
    };
    if (ul) |x| {
        try balance.back(io, T, bal, scale, true, x);
        normalizeVectors(io, T, x, if (comptime isComplex(T)) null else w.im.ptr);
    }
    if (ur) |x| {
        try balance.back(io, T, bal, scale, false, x);
        normalizeVectors(io, T, x, if (comptime isComplex(T)) null else w.im.ptr);
    }
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
