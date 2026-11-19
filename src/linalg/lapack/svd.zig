//! `?GESVD` and `?GESDD`: Reference-LAPACK's arguments, jobs and scaling
//! (`lapack-netlib/SRC/dgesvd.f`, `dgesdd.f` and their `s`, `c`, `z` forms,
//! OpenBLAS v0.3.34-317) over faer v0.24.4's SVD
//! (<https://codeberg.org/sarah-quinones/faer>, tag `faer-v0.24.4`, MIT,
//! `faer/src/linalg/svd/mod.rs`). Departures and fixes: design.md's LAPACK
//! section.

const std = @import("std");
const ops = @import("ops.zig");
const auxiliary = @import("auxiliary.zig");
const hh = @import("householder.zig");
const qr = @import("qr.zig");
const bidiag = @import("bidiag.zig");
const bsvd = @import("bidiag_svd.zig");
const util = @import("../blas/util.zig");
const interface = @import("interface.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = ops.Mat;
const Float = ops.Float;
const isComplex = ops.isComplex;

/// `JOBU`, `JOBVT`, `JOBZ`: all of `U` (`m x m`) or `Vᴴ` (`n x n`), the first
/// `min(m, n)` columns of `U` or rows of `Vᴴ` (`s`), those overwritten on
/// `A` (`o`), or none (`n`).
pub const Job = enum { a, s, o, n };

/// What `gesvd` returns: `ConflictingJobs` when both jobs are `o`, a slice
/// too short (`a`, `s` under `min(m, n)`, `u` and `vt` where their jobs
/// write them) or a leading dimension too small, `NoConvergence`.
pub const GesvdError = error{ConflictingJobs} || interface.ShapeError || error{NoConvergence} || Allocator.Error || Io.Cancelable;

/// What `gesdd` returns: `gesvd`'s shape errors, `NanInA` when `A` holds a
/// NaN, `NoConvergence`.
pub const GesddError = interface.ShapeError || interface.NanError || error{NoConvergence} || Allocator.Error || Io.Cancelable;

/// `SvdParams::auto`'s `recursion_threshold`: below it, divide and conquer
/// is QR iteration.
pub const dc_threshold = 128;

/// `SvdParams::auto`'s `qr_ratio_threshold`.
const qr_ratio = 11.0 / 6.0;

/// `?GESVD`: `A = U Σ Vᴴ` for the `m x n` `A`; the singular values into `s`
/// (`min(m, n)`, nonincreasing, nonnegative), `U` by `jobu` into `u` (`ldu`)
/// or `A`, `Vᴴ` by `jobvt` into `vt` (`ldvt`) or `A`. `A` is destroyed
/// unless a job puts its vectors there. `u` and `vt` are neither read nor
/// checked where their job does not use them. The bidiagonal by QR
/// iteration.
///
/// On `error.NoConvergence`, `unconverged` (when given) gets the reference's
/// `INFO > 0`: how many superdiagonals of the bidiagonal did not converge to
/// zero (`min(m, n) - 1`, all of them, for an infinity or NaN in `A`); `s`,
/// `U` and `Vᴴ` are then undefined. The reference also leaves that
/// bidiagonal in `S` and `WORK`; rstd has no `WORK`.
pub fn gesvd(allocator: Allocator, io: Io, comptime T: type, jobu: Job, jobvt: Job, m: usize, n: usize, a: []T, lda: usize, s: []Float(T), u: []T, ldu: usize, vt: []T, ldvt: usize, unconverged: ?*usize) GesvdError!void {
    comptime ops.requireFloat(T, "gesvd");
    const minmn = @min(m, n);
    if (jobu == .o and jobvt == .o) return error.ConflictingJobs;
    try util.holds(a.len, m, n, lda);
    if (s.len < minmn) return error.SliceTooShort;
    if (jobu == .a or jobu == .s) try util.holds(u.len, m, if (jobu == .a) m else minmn, ldu);
    if (jobvt == .a or jobvt == .s) try util.holds(vt.len, if (jobvt == .a) n else minmn, n, ldvt);
    if (m == 0 or n == 0) return;
    run(allocator, io, T, jobu, jobvt, m, n, a.ptr, lda, s.ptr, u.ptr, ldu, vt.ptr, ldvt, std.math.maxInt(usize), unconverged) catch |e| return switch (e) {
        error.NanInA => error.NoConvergence,
        else => |x| x,
    };
}

/// `?GESDD`: as `gesvd`, with one job for both: `a`, `s` and `n` as
/// `gesvd`'s for `U` and `Vᴴ` alike; `o` puts `U` on `A` and all of `Vᴴ` in
/// `vt` when `m >= n`, else all of `U` in `u` and `Vᴴ` on `A`. The
/// bidiagonal by divide and conquer.
pub fn gesdd(allocator: Allocator, io: Io, comptime T: type, jobz: Job, m: usize, n: usize, a: []T, lda: usize, s: []Float(T), u: []T, ldu: usize, vt: []T, ldvt: usize) GesddError!void {
    comptime ops.requireFloat(T, "gesdd");
    const minmn = @min(m, n);
    try util.holds(a.len, m, n, lda);
    if (s.len < minmn) return error.SliceTooShort;
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
    if (ucols > 0) try util.holds(u.len, m, ucols, ldu);
    if (vtrows > 0) try util.holds(vt.len, vtrows, n, ldvt);
    if (m == 0 or n == 0) return;
    const jobu: Job, const jobvt: Job = switch (jobz) {
        .o => if (m >= n) .{ .o, .a } else .{ .a, .o },
        .a, .s, .n => .{ jobz, jobz },
    };
    try run(allocator, io, T, jobu, jobvt, m, n, a.ptr, lda, s.ptr, u.ptr, ldu, vt.ptr, ldvt, dc_threshold, null);
}

const RunError = error{ NanInA, NoConvergence } || ops.Error;

/// Both drivers after their checks: the scaling of `dgesvd.f:652-667` and
/// `dgesdd.f:585-606`, faer's SVD, and the vectors where the jobs put them.
fn run(al: Allocator, io: Io, comptime T: type, jobu: Job, jobvt: Job, m: usize, n: usize, a: [*]T, lda: usize, s: [*]Float(T), u: [*]T, ldu: usize, vt: [*]T, ldvt: usize, threshold: usize, unconverged: ?*usize) RunError!void {
    const F = Float(T);
    const minmn = @min(m, n);
    const am: Mat(T) = .{ .ptr = a, .rows = m, .cols = n, .ld = lda };
    const scaled = auxiliary.scaleIn(T, am) catch |e| {
        if (unconverged) |c| c.* = minmn - 1;
        return switch (e) {
            error.NanInA => error.NanInA,
            error.InfInA => error.NoConvergence,
        };
    };
    const anrm = scaled.anrm;
    const to = scaled.to;
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
    try svd(al, io, T, am, sv, um, vm, threshold, unconverged);
    if (to) |t| {
        // `DLASCL` on S, back
        auxiliary.lascl(F, false, t, anrm, .{ .ptr = sv.ptr, .rows = minmn, .cols = 1, .ld = @max(minmn, 1) });
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

/// `svd`: `s` (`min(m, n)`), `u` (`m x m` or `m x min(m, n)`) and `v`
/// (`n x n` or `n x min(m, n)`) from `a`, which is only read; `threshold`
/// as `bidiagSvd`'s. A wide `a` is decomposed transposed, `U` and `V`
/// swapped and conjugated after.
pub fn svd(al: Allocator, io: Io, comptime T: type, a: Mat(T), s: []Float(T), u: ?Mat(T), v: ?Mat(T), threshold: usize, unconverged: ?*usize) bsvd.Error!void {
    const m = a.rows;
    const n = a.cols;
    if (m >= n) return svdTall(al, io, T, a, s, u, v, threshold, unconverged);
    const at = try ops.alloc(T, al, n, m);
    defer ops.free(T, al, at);
    for (0..m) |j| for (0..n) |i| {
        at.at(i, j).* = a.get(j, i);
    };
    try svdTall(al, io, T, at, s, v, u, threshold, unconverged);
    if (comptime isComplex(T)) for ([_]?Mat(T){ u, v }) |mo| if (mo) |x| {
        for (0..x.cols) |j| for (x.col(j)[0..x.rows]) |*e| {
            e.* = ops.conj(T, e.*);
        };
    };
}

/// `svd` for `m >= n`.
fn svdTall(al: Allocator, io: Io, comptime T: type, a: Mat(T), s: []Float(T), u: ?Mat(T), v: ?Mat(T), threshold: usize, unconverged: ?*usize) bsvd.Error!void {
    const m = a.rows;
    const n = a.cols;
    if (n == 0) {
        if (u) |x| x.identity();
        return;
    }
    if (@as(f64, @floatFromInt(m)) / @as(f64, @floatFromInt(n)) <= qr_ratio) return squareish(al, io, T, a, s, u, v, threshold, unconverged);
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
        try squareish(al, io, T, r, s, if (u) |x| x.sub(0, 0, n, n) else null, v, threshold, unconverged);
    }
    if (u) |x| {
        x.sub(n, 0, m - n, x.cols).fill(ops.zero(T));
        if (x.cols > n) x.sub(0, n, n, x.cols - n).fill(ops.zero(T));
        for (n..x.cols) |i| x.at(i, i).* = ops.one(T);
        try hh.applySequence(al, io, T, q, h, false, x);
    }
}

/// `svd_imp`: bidiagonalize, decompose `Bᴴ`, rebuild.
fn squareish(al: Allocator, io: Io, comptime T: type, a: Mat(T), s: []Float(T), u: ?Mat(T), v: ?Mat(T), threshold: usize, unconverged: ?*usize) bsvd.Error!void {
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
    try bidiagSvd(al, io, T, diag, subdiag, ub, vb, threshold, unconverged);
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
fn bidiagSvd(al: Allocator, io: Io, comptime T: type, diag: []T, subdiag: []T, u: ?Mat(T), v: ?Mat(T), threshold: usize, unconverged: ?*usize) bsvd.Error!void {
    if (comptime !isComplex(T)) return realBidiagSvd(al, io, T, diag, subdiag, u, v, threshold, unconverged);
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
    var cn = ops.normalized(T, ops.conj(T, diag[0]));
    col_mul[0] = cn;
    dr[0] = ops.abs(T, diag[0]);
    sr[n - 1] = 0;
    for (1..n) |i| {
        const rn = ops.normalized(T, ops.conj(T, ops.mul(T, subdiag[i - 1], cn)));
        sr[i - 1] = ops.abs(T, subdiag[i - 1]);
        row_mul[i - 1] = ops.conj(T, rn);
        cn = ops.normalized(T, ops.conj(T, ops.mul(T, diag[i], rn)));
        dr[i] = ops.abs(T, diag[i]);
        col_mul[i] = cn;
    }
    const ur: ?Mat(F) = if (u != null) try ops.alloc(F, al, n + 1, n + 1) else null;
    defer if (ur) |x| ops.free(F, al, x);
    const vr: ?Mat(F) = if (v != null) try ops.alloc(F, al, n, n) else null;
    defer if (vr) |x| ops.free(F, al, x);
    try realBidiagSvd(al, io, F, dr, sr, ur, vr, threshold, unconverged);
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

fn realBidiagSvd(al: Allocator, io: Io, comptime F: type, diag: []F, subdiag: []F, u: ?Mat(F), v: ?Mat(F), threshold: usize, unconverged: ?*usize) bsvd.Error!void {
    const n = diag.len;
    for (diag, subdiag) |d, e| if (!(std.math.isFinite(d) and std.math.isFinite(e))) {
        // every superdiagonal counts as unconverged
        if (unconverged) |c| c.* = n - 1;
        return error.NoConvergence;
    };
    if (n < threshold) {
        if (u) |x| x.identity();
        if (v) |x| x.identity();
        return bsvd.qrAlgorithm(F, diag, subdiag, if (u) |x| x.sub(0, 0, n, n) else null, v) catch |e| {
            // `xBDSQR`'s `INFO`: the superdiagonals not yet zero
            if (unconverged) |c| {
                c.* = 0;
                for (subdiag[0 .. n - 1]) |x| c.* += @intFromBool(x != 0);
            }
            return e;
        };
    }
    if (u) |x| return bsvd.divideAndConquer(al, io, F, diag, subdiag, .{ .full = x }, v, threshold);
    const u_two = try ops.alloc(F, al, 2, n + 1);
    defer ops.free(F, al, u_two);
    return bsvd.divideAndConquer(al, io, F, diag, subdiag, .{ .two_rows_storage = u_two }, v, threshold);
}
