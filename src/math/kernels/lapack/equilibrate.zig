//! `xGEEQU` and `xLAQGE`: row and column equilibration of a general matrix,
//! Reference-LAPACK 3.12.0 (`lapack-netlib/SRC/dgeequ.f`, `zgeequ.f`,
//! `dlaqge.f`, `zlaqge.f`, OpenBLAS v0.3.34-317; modified BSD).

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const interface = @import("interface.zig");
const ops = @import("ops.zig");

const Float = types.Float;
const isComplex = types.isComplex;
pub const Equilibration = interface.Equilibration;

/// `geequ`'s measures besides the factors: the smallest over the largest
/// row factor's reciprocal, the same for columns, and `A`'s largest entry.
pub fn Scaling(comptime F: type) type {
    return struct { rowcnd: F, colcnd: F, amax: F };
}

/// What `geequ` returns: a slice too short or a leading dimension too small;
/// `ZeroRow` or `ZeroColumn` for the reference's `INFO > 0`, with `r` and `c`
/// then partly written.
pub const GeequError = error{ ZeroRow, ZeroColumn } || interface.ShapeError;

/// `?GEEQU` (`dgeequ.f`, `zgeequ.f`): row factors `r` (`m`) and column
/// factors `c` (`n`) that bring the largest entry of each row and column of
/// `diag(r) A diag(c)` (`m x n`) near 1, each kept in `[smlnum, bignum]`
/// (`smlnum` the safe minimum). Entries are measured by `|re| + |im|` for
/// complex, as the reference's `CABS1`. For `m = 0` or `n = 0` returns
/// `rowcnd = colcnd = 1` and `amax = 0`.
pub fn geequ(comptime T: type, m: usize, n: usize, a: []const T, lda: usize, r: []Float(T), c: []Float(T)) GeequError!Scaling(Float(T)) {
    comptime ops.requireFloat(T, "geequ");
    const F = Float(T);
    try util.holds(a.len, m, n, lda);
    if (r.len < m or c.len < n) return error.SliceTooShort;
    if (m == 0 or n == 0) return .{ .rowcnd = 1, .colcnd = 1, .amax = 0 };
    const smlnum = std.math.floatMin(F);
    const bignum = 1 / smlnum;
    @memset(r[0..m], 0);
    for (0..n) |j| for (0..m) |i| {
        r[i] = @max(r[i], ops.abs1(T, a[i + j * lda]));
    };
    var rcmin = bignum;
    var rcmax: F = 0;
    for (r[0..m]) |v| {
        rcmax = @max(rcmax, v);
        rcmin = @min(rcmin, v);
    }
    const amax = rcmax;
    if (rcmin == 0) return error.ZeroRow;
    for (r[0..m]) |*v| v.* = 1 / @min(@max(v.*, smlnum), bignum);
    const rowcnd = @max(rcmin, smlnum) / @min(rcmax, bignum);
    @memset(c[0..n], 0);
    for (0..n) |j| for (0..m) |i| {
        c[j] = @max(c[j], ops.abs1(T, a[i + j * lda]) * r[i]);
    };
    rcmin = bignum;
    rcmax = 0;
    for (c[0..n]) |v| {
        rcmin = @min(rcmin, v);
        rcmax = @max(rcmax, v);
    }
    if (rcmin == 0) return error.ZeroColumn;
    for (c[0..n]) |*v| v.* = 1 / @min(@max(v.*, smlnum), bignum);
    return .{ .rowcnd = rowcnd, .colcnd = @max(rcmin, smlnum) / @min(rcmax, bignum), .amax = amax };
}

/// What `laqge` returns: a slice too short or a leading dimension too small.
pub const LaqgeError = interface.ShapeError;

/// `dlaqge.f`'s `THRESH`: a ratio below it calls for scaling.
const thresh = 0.1;

/// `?LAQGE` (`dlaqge.f`, `zlaqge.f`): scales `A` (`m x n`) by `geequ`'s
/// `r` and `c` where they matter: rows when `scaling.rowcnd < 0.1` or
/// `amax` is outside `[small, 1 / small]` (`small` the safe minimum over
/// the precision), columns when `colcnd < 0.1`. Returns what it scaled.
pub fn laqge(comptime T: type, m: usize, n: usize, a: []T, lda: usize, r: []const Float(T), c: []const Float(T), scaling: Scaling(Float(T))) LaqgeError!Equilibration {
    comptime ops.requireFloat(T, "laqge");
    const F = Float(T);
    try util.holds(a.len, m, n, lda);
    if (r.len < m or c.len < n) return error.SliceTooShort;
    if (m == 0 or n == 0) return .none;
    const small = std.math.floatMin(F) / std.math.floatEps(F);
    const large = 1 / small;
    const rows = !(scaling.rowcnd >= thresh and scaling.amax >= small and scaling.amax <= large);
    const cols = scaling.colcnd < thresh;
    if (!rows and !cols) return .none;
    for (0..n) |j| for (0..m) |i| {
        const s: F = (if (cols) c[j] else 1) * (if (rows) r[i] else 1);
        const x = &a[i + j * lda];
        x.* = if (comptime isComplex(T)) ops.mulReal(T, x.*, s) else s * x.*;
    };
    return if (rows and cols) .both else if (rows) .rows else .columns;
}
