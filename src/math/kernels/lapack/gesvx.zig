//! `xGESVX`: the expert LU driver, equilibration, refinement and error bounds,
//! Reference-LAPACK 3.12.0 (`lapack-netlib/SRC/dgesvx.f`, `zgesvx.f`, OpenBLAS
//! v0.3.34-317; modified BSD).

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const interface = @import("interface.zig");
const ops = @import("ops.zig");
const norm = @import("norm.zig");
const norm_structured = @import("norm_structured.zig");
const equilibrate = @import("equilibrate.zig");
const refine = @import("gerfs.zig");
const lu = @import("lu.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const Trans = types.Trans;
const NormType = norm.NormType;
const Equilibration = interface.Equilibration;
pub const FactorizationInput = interface.FactorizationInput;

/// What `gesvx` reports besides `X`: `rcond`, `gecon`'s estimate of the
/// reciprocal condition number of `A` after equilibration (0 when
/// singular); `pivot_growth`, `max |A| / max |U|` (1 for a zero `U`), over
/// the leading `singular + 1` columns when singular. `singular` is the
/// first `k` (0-based) whose `U(k, k)` is exactly zero, the reference's
/// `INFO = k + 1`: `X`, `ferr` and `berr` are then not computed.
/// `ill_conditioned` is the reference's `INFO = n + 1`, `rcond` below the
/// unit roundoff (`DLAMCH('Epsilon')`): `X` is computed, but `A` is singular to working
/// precision.
pub fn LuSolveReport(comptime F: type) type {
    return struct { rcond: F, pivot_growth: F, singular: ?usize = null, ill_conditioned: bool = false };
}

/// What `gesvx` returns: a `trans` other than `n`, `t` or `c`; for
/// `.factored`, `ScaleFactorNotPositive` for a factor of `r` or `c` that
/// `equed` applies at or below zero (the reference's `INFO = -11`, `-12`) and
/// `PivotOutOfRange` for a pivot that is not a row; a slice too short or a
/// leading dimension too small; `gecon`'s `NormOutOfRange` for an `A` with a
/// NaN or infinite entry and `EstimateFailed`; the allocator's and
/// cancellation.
pub const GesvxError = error{ ScaleFactorNotPositive, NormOutOfRange, EstimateFailed } || interface.TransError || interface.PivotError || interface.ShapeError || Allocator.Error || Io.Cancelable;

/// `?GESVX` (`dgesvx.f`, `zgesvx.f`): solves `op(A) X = B` for the
/// `n x nrhs` `X`, `trans` `n`, `t` or `c` (real `c` is `t`). By `fact`:
/// `.factor` copies `A` to `af` and factors it there (`getrf`, pivots in
/// `ipiv`); `.equilibrate_and_factor` first scales `A` in place to
/// `diag(r) A diag(c)` where `geequ` and `laqge` find it worthwhile,
/// `equed` saying what was scaled; `.factored` takes `af` and `ipiv` as a
/// factor of `A` as given, scaled as `equed`, `r` and `c` say. `B` is
/// scaled to match and kept so on return; `X` is the solution of the
/// original system. Refined by `gerfs`, with each column's forward bound
/// `ferr` and backward error `berr`. `r` and `c` hold `n` factors each.
pub fn gesvx(allocator: Allocator, io: Io, comptime T: type, fact: FactorizationInput, trans: Trans, n: usize, nrhs: usize, a: []T, lda: usize, af: []T, ldaf: usize, ipiv: []usize, equed: *Equilibration, r: []Float(T), c: []Float(T), b: []T, ldb: usize, x: []T, ldx: usize, ferr: []Float(T), berr: []Float(T)) GesvxError!LuSolveReport(Float(T)) {
    comptime ops.requireFloat(T, "gesvx");
    const F = Float(T);
    const notran = switch (trans) {
        .n => true,
        .t, .c => false,
        .r, .o, .u, .s, .d => return error.UnsupportedTranspose,
    };
    try util.holds(a.len, n, n, lda);
    try util.holds(af.len, n, n, ldaf);
    try util.holds(b.len, n, nrhs, ldb);
    try util.holds(x.len, n, nrhs, ldx);
    if (ipiv.len < n or r.len < n or c.len < n or ferr.len < nrhs or berr.len < nrhs) return error.SliceTooShort;
    var rowcnd: F = 1;
    var colcnd: F = 1;
    if (fact == .factored) {
        try interface.pivotsIn(ipiv, n);
        if (scalesRows(equed.*)) rowcnd = try ratio(F, r[0..n]);
        if (scalesColumns(equed.*)) colcnd = try ratio(F, c[0..n]);
    } else equed.* = .none;
    // as the reference's calls give for `n = 0`: `gecon`'s 1, no growth, `gerfs`'s zeros
    if (n == 0) {
        @memset(ferr[0..nrhs], 0);
        @memset(berr[0..nrhs], 0);
        return .{ .rcond = 1, .pivot_growth = 1 };
    }

    if (fact == .equilibrate_and_factor) {
        // the reference's `INFEQU > 0`, a zero row or column, leaves `A` as it is
        if (equilibrate.geequ(T, n, n, a, lda, r, c)) |s| {
            equed.* = try equilibrate.laqge(T, n, n, a, lda, r, c, s);
            rowcnd = s.rowcnd;
            colcnd = s.colcnd;
        } else |e| switch (e) {
            error.ZeroRow, error.ZeroColumn => {},
            else => |err| return err,
        }
    }
    const rows = scalesRows(equed.*);
    const cols = scalesColumns(equed.*);
    if (if (notran) rows else cols) scale(T, n, nrhs, b, ldb, if (notran) r else c);

    if (fact != .factored) {
        for (0..n) |j| @memcpy(af[j * ldaf ..][0..n], a[j * lda ..][0..n]);
        if (try lu.getrf(allocator, io, T, n, n, af, ldaf, ipiv)) |k| {
            return .{ .rcond = 0, .pivot_growth = try growth(allocator, io, T, n, k + 1, a, lda, af, ldaf), .singular = k };
        }
    }
    const kind: NormType = if (notran) .one else .infinity;
    const anorm = try norm.lange(allocator, io, T, kind, n, n, a, lda);
    const pivot_growth = try growth(allocator, io, T, n, n, a, lda, af, ldaf);
    const rcond = norm.gecon(allocator, io, T, kind, n, af, ldaf, anorm) catch |e| switch (e) {
        // `kind` is one gecon takes
        error.UnsupportedNorm => unreachable,
        else => |err| return err,
    };
    for (0..nrhs) |j| @memcpy(x[j * ldx ..][0..n], b[j * ldb ..][0..n]);
    try lu.getrs(allocator, io, T, trans, n, nrhs, af, ldaf, ipiv, x, ldx);
    try refine.gerfs(allocator, io, T, trans, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr);
    // back to the original system's `X`
    if (if (notran) cols else rows) {
        scale(T, n, nrhs, x, ldx, if (notran) c else r);
        const cnd = if (notran) colcnd else rowcnd;
        for (ferr[0..nrhs]) |*e| e.* /= cnd;
    }
    // `DLAMCH('Epsilon')`, half of `eps`
    return .{ .rcond = rcond, .pivot_growth = pivot_growth, .ill_conditioned = rcond < std.math.floatEps(F) / 2 };
}

fn scalesRows(e: Equilibration) bool {
    return e == .rows or e == .both;
}

fn scalesColumns(e: Equilibration) bool {
    return e == .columns or e == .both;
}

/// The smallest factor over the largest, each clamped to `[smlnum, bignum]`
/// (`dgesvx.f`'s `ROWCND`, `COLCND` for given factors); 1 for none.
fn ratio(comptime F: type, s: []const F) error{ScaleFactorNotPositive}!F {
    const smlnum = std.math.floatMin(F);
    const bignum = 1 / smlnum;
    var lo = bignum;
    var hi: F = 0;
    for (s) |v| {
        lo = @min(lo, v);
        hi = @max(hi, v);
    }
    if (lo <= 0) return error.ScaleFactorNotPositive;
    if (s.len == 0) return 1;
    return @max(lo, smlnum) / @min(hi, bignum);
}

/// Row `i` of the `n x nrhs` `m` times `s[i]`.
fn scale(comptime T: type, n: usize, nrhs: usize, m: []T, ld: usize, s: []const Float(T)) void {
    for (0..nrhs) |j| for (m[j * ld ..][0..n], s[0..n]) |*e, f| {
        e.* = ops.mulReal(T, e.*, f);
    };
}

/// The reciprocal pivot growth over the leading `k` columns: the largest
/// `|A(i, j)|` over the largest `|U(i, j)|`, 1 for a zero `U` (`dgesvx.f`'s
/// `RPVGRW`, by `lange` and `lantr`'s `'M'`).
fn growth(al: Allocator, io: Io, comptime T: type, n: usize, k: usize, a: []const T, lda: usize, af: []const T, ldaf: usize) (interface.ShapeError || Allocator.Error)!Float(T) {
    const umax = norm_structured.lantr(al, io, T, .max_abs, .u, .n, k, k, af, ldaf) catch |e| switch (e) {
        // `.u` is a triangle lantr takes
        error.UnsupportedUplo => unreachable,
        else => |err| return err,
    };
    if (umax == 0) return 1;
    return try norm.lange(al, io, T, .max_abs, n, k, a, lda) / umax;
}
