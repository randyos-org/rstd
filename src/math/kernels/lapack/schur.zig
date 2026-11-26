//! The Schur form: `hseqr` and `gees` over the QR iterations of
//! `real_schur.zig` and `complex_schur.zig`, and the names of its reordering,
//! `trexc` and `trsen` (`trexc.zig`, `trsen.zig`, over `laexc.zig`'s 2 x 2
//! kernels), from Reference-LAPACK 3.12.0 as OpenBLAS v0.3.34-317 ships it
//! (`lapack-netlib/SRC/`, modified BSD).

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const interface = @import("interface.zig");
const ops = @import("ops.zig");
const auxiliary = @import("auxiliary.zig");
const hessenberg = @import("hessenberg.zig");
const real_schur = @import("real_schur.zig");
const complex_schur = @import("complex_schur.zig");
const balance = @import("balance.zig");
const strided = @import("../blas/strided.zig");
const laexc_swap = @import("laexc.zig");
const trexc_move = @import("trexc.zig");
const trsen_reorder = @import("trsen.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const isComplex = types.isComplex;
const Mat = ops.Mat;
const VectorJob = interface.VectorJob;

pub const lanv2 = laexc_swap.lanv2;
pub const TrexcError = trexc_move.TrexcError;
pub const trexc = trexc_move.trexc;
pub const ConditionJob = trsen_reorder.ConditionJob;
pub const TrsenError = trsen_reorder.TrsenError;
pub const ReorderedSchur = trsen_reorder.ReorderedSchur;
pub const trsen = trsen_reorder.trsen;

/// LAPACK's `COMPQ`: update `Q` (`V`) or not (`N`).
pub const ReorderVectors = enum { none, update };

/// `geev`'s eigenvalues: real `wr` and `wi` (`dgeev`'s), complex `w`
/// (`zgeev`'s).
pub fn Eigenvalues(comptime T: type) type {
    return if (isComplex(T)) []T else struct { re: []T, im: []T };
}

/// `w` holds `n` eigenvalues.
pub fn eigenvaluesHold(comptime T: type, w: Eigenvalues(T), n: usize) strided.SliceError!void {
    const short = if (comptime isComplex(T)) w.len < n else w.re.len < n or w.im.len < n;
    if (short) return error.SliceTooShort;
}

// ---- xHSEQR and xGEES (moved from evd.zig, issue #87) -----------------------------

/// `hseqr`'s `JOB`: the eigenvalues alone (`E`), or the Schur form too (`S`).
pub const SchurJob = enum { eigenvalues, schur_form };
/// `hseqr`'s `COMPZ`: no `Z` (`N`), from the identity (`I`), or accumulated
/// into the given one (`V`).
pub const SchurVectors = enum { none, from_identity, accumulate };

/// What `hseqr` returns: `BlockOutOfRange` unless `lo <= hi <= n`, a slice
/// too short or a leading dimension too small (`z` checked only with
/// `compz` other than `.none`), `NoConvergence` for the reference's `INFO > 0`.
pub const HseqrError = interface.BlockError || interface.ShapeError || error{NoConvergence};

/// `?HSEQR` (`dhseqr.f`, `zhseqr.f`): the eigenvalues of the upper
/// Hessenberg `n x n` `H` into `w`, as `geev`'s, by faer's `lahqr`
/// (`real_schur.zig`, `complex_schur.zig`; faer's multishift QR for large
/// sizes is not ported yet) on rows and columns `lo..hi` (0-based, `hi`
/// past the block's end, as `gebal` returns them; the rest triangular
/// already, its eigenvalues read off the diagonal). With
/// `job` `.schur_form`, `H` becomes the Schur form `T` (real: quasi-triangular, its
/// 2 x 2 blocks standardized), else its contents are unspecified; entries
/// below the subdiagonal are zeroed first, as the reference clears them
/// after. With `compz` `.from_identity` or `.accumulate`, `Z` (`n x n`) is the identity or the
/// caller's, times the Schur vectors. `NoConvergence` writes the
/// reference's `INFO` to `unconverged`: `w[0..lo]` and `w[INFO..n]`
/// converged.
pub fn hseqr(io: Io, comptime T: type, job: SchurJob, compz: SchurVectors, n: usize, lo: usize, hi: usize, h: []T, ldh: usize, w: Eigenvalues(T), z: []T, ldz: usize, unconverged: ?*usize) HseqrError!void {
    comptime ops.requireFloat(T, "hseqr");
    try interface.blockIn(lo, hi, n);
    try util.holds(h.len, n, n, ldh);
    try eigenvaluesHold(T, w, n);
    if (compz != .none) try util.holds(z.len, n, n, ldz);
    if (n == 0) return;
    const un = n;
    const hm: Mat(T) = .{ .ptr = h.ptr, .rows = un, .cols = un, .ld = ldh };
    const zm: ?Mat(T) = if (compz == .none) null else .{ .ptr = z.ptr, .rows = un, .cols = un, .ld = ldz };
    if (compz == .from_identity) zm.?.identity();
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
    const want = job == .schur_form;
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
pub fn EigenvalueSelector(comptime T: type) type {
    return if (isComplex(T)) *const fn (T) bool else *const fn (Float(T), Float(T)) bool;
}

/// What `gees` returns: a slice too short or a leading dimension too small
/// (`vs` checked only with `jobvs` `.compute`), `NanInA` for a NaN in `A`,
/// `NoConvergence` for the reference's `INFO` 1 to `n`, `SelectChanged` for
/// its `n + 1` (rounding in the reordering changed which eigenvalues
/// `select` picks), `SwapRejected` for its `n + 2` (a swap too
/// ill-conditioned to make).
pub const GeesError = interface.ShapeError || interface.NanError || error{ NoConvergence, SelectChanged, SwapRejected } || Allocator.Error || Io.Cancelable;

/// `?GEES` (`dgees.f`, `zgees.f`): the Schur factorization `A = Z T Zᴴ` of
/// the `n x n` `A`: `T` in place of `A` (real: quasi-triangular, its 2 x 2
/// blocks standardized), its eigenvalues in `w` (as `geev`'s), and with
/// `jobvs` `.compute` the Schur vectors `Z` in `vs`. As the reference: scaled into
/// `[SMLNUM, BIGNUM]`, balanced by permutation only (`gebal` job `.permute`),
/// `gehrd`, `orghr`, `hseqr`; with `select`, `trsen` then moves the
/// eigenvalues it picks to the top-left; `Z` takes `gebak` back. Returns
/// the reference's `sdim`: the count `select` picked (a real pair counts 2,
/// and moves if either is picked), 0 without `select`. `NoConvergence`
/// writes the reference's `INFO` to `unconverged`.
pub fn gees(allocator: Allocator, io: Io, comptime T: type, jobvs: VectorJob, select: ?EigenvalueSelector(T), n: usize, a: []T, lda: usize, w: Eigenvalues(T), vs: []T, ldvs: usize, unconverged: ?*usize) GeesError!usize {
    comptime ops.requireFloat(T, "gees");
    const F = Float(T);
    try util.holds(a.len, n, n, lda);
    try eigenvaluesHold(T, w, n);
    if (jobvs == .compute) try util.holds(vs.len, n, n, ldvs);
    if (n == 0) return 0;
    const un = n;
    const am: Mat(T) = .{ .ptr = a.ptr, .rows = un, .cols = un, .ld = lda };
    // `dgees.f`: NaN, then the scale into [SMLNUM, BIGNUM], as `geev`
    const scaled = auxiliary.scaleIn(T, am) catch |e| switch (e) {
        error.NanInA => return error.NanInA,
        error.InfInA => {
            if (unconverged) |c| c.* = un;
            return error.NoConvergence;
        },
    };
    const anrm = scaled.anrm;
    const to = scaled.to;
    const scale = try allocator.alloc(F, un);
    defer allocator.free(scale);
    const bal = balance.gebal(io, T, .permute, n, a, lda, scale) catch |e| switch (e) {
        error.Canceled => |x| return x,
        // checked above: no NaN, and the arguments
        else => unreachable,
    };
    const tau = try allocator.alloc(T, un);
    defer allocator.free(tau);
    hessenberg.gehrd(allocator, io, T, n, bal.lo, bal.hi, a, lda, tau) catch |e| return ops.narrow(ops.Error, e);
    if (jobvs == .compute) {
        const vm: Mat(T) = .{ .ptr = vs.ptr, .rows = un, .cols = un, .ld = ldvs };
        vm.copyFrom(am);
        hessenberg.orghr(allocator, io, T, n, bal.lo, bal.hi, vs, ldvs, tau) catch |e| return ops.narrow(ops.Error, e);
    }
    const compz: SchurVectors = if (jobvs == .compute) .accumulate else .none;
    hseqr(io, T, .schur_form, compz, n, bal.lo, bal.hi, a, lda, w, vs, ldvs, unconverged) catch |e| switch (e) {
        error.NoConvergence => |x| return x,
        else => unreachable,
    };
    var sdim: usize = 0;
    if (select) |pick| {
        const bwork = try allocator.alloc(bool, un);
        defer allocator.free(bwork);
        for (bwork, 0..) |*b, i| b.* = picked(T, pick, w, i);
        const r = trsen(allocator, io, T, .none, if (jobvs == .compute) .update else .none, bwork, n, a, lda, vs, ldvs, w) catch |e| switch (e) {
            error.SwapRejected, error.OutOfMemory, error.Canceled => |x| return x,
            else => unreachable,
        };
        sdim = r.selected;
        // `dgees.f`'s `INFO = N + 1`: the reordered eigenvalues, re-picked,
        // a real pair picked if either is
        const sm = r.selected;
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
    if (jobvs == .compute) balance.gebak(io, T, .permute, .r, n, bal.lo, bal.hi, scale, n, vs, ldvs) catch |e| switch (e) {
        error.Canceled => |x| return x,
        else => unreachable,
    };
    if (to) |t| {
        auxiliary.lascl(T, false, t, anrm, am);
        const ws = if (comptime isComplex(T)) [_][*]T{w.ptr} else [_][*]T{ w.re.ptr, w.im.ptr };
        for (ws) |x| auxiliary.lascl(T, false, t, anrm, .{ .ptr = x, .rows = un, .cols = 1, .ld = un });
    }
    return sdim;
}

/// Whether `select` picks eigenvalue `i`, by its own `(re, im)` for real,
/// as `dgees.f` asks each.
fn picked(comptime T: type, pick: EigenvalueSelector(T), w: Eigenvalues(T), i: usize) bool {
    return if (comptime isComplex(T)) pick(w[i]) else pick(w.re[i], w.im[i]);
}
