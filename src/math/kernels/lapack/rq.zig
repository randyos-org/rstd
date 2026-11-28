//! RQ: `?GERQF`, `?ORGRQ` and `?ORMRQ` behind Reference-LAPACK's arguments
//! (`dgerqf.f`, `dorgrq.f`, `dormrq.f` and their `z` forms, 3.12.0 as
//! OpenBLAS v0.3.34-794 ships it; modified BSD) over `lq.zig`'s LQ: with
//! `J` the reversal, `J A J = (J R J)(J Q J)` is an LQ, whose reflections,
//! reversed, sit where `dgerqf.f` puts its own (`conj(v)` in each row).

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const interface = @import("interface.zig");
const ops = @import("ops.zig");
const lq = @import("lq.zig");
const reversal = @import("reversal.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = ops.Mat;
const Side = types.Side;
const Trans = types.Trans;
const flip = reversal.flip;
const reverse = reversal.reverse;
const reverseAlong = reversal.reverseAlong;

/// What `gerqf` returns: a slice too short (`a`, or `tau` under
/// `min(m, n)`) or a leading dimension too small, and the allocator's and
/// cancellation.
pub const GerqfError = lq.GelqfError;
/// What `orgrq` returns: `TooManyRows` for `m > n`, `TooManyReflections`
/// for `k > m`, and `GerqfError`'s.
pub const OrgrqError = lq.OrglqError;
/// What `ormrq` returns: a `trans` other than `n`, `t` or `c` (complex: `n`
/// or `c`), `TooManyReflections` for `k` past the side's order, and
/// `GerqfError`'s.
pub const OrmrqError = lq.OrmlqError;

/// `?GERQF`: `A = R Q` for the `m x n` `A`: `R` in the upper trapezoid
/// ending at `A(m - 1, n - 1)`, `Q = H(0)ᴴ ... H(k-1)ᴴ` (`k = min(m, n)`)
/// as reflections left of it (`conj(v)` in row `m - k + i`, left of column
/// `n - k + i`) and their `tau`s in `tau` (`k`). As `gelqf`'s, a complex
/// `R`'s diagonal stays complex.
pub fn gerqf(allocator: Allocator, io: Io, comptime T: type, m: usize, n: usize, a: []T, lda: usize, tau: []T) GerqfError!void {
    comptime ops.requireFloat(T, "gerqf");
    try util.holds(a.len, m, n, lda);
    const k = @min(m, n);
    if (tau.len < k) return error.SliceTooShort;
    if (k == 0) return;
    const am: Mat(T) = .{ .ptr = a.ptr, .rows = m, .cols = n, .ld = lda };
    // the QR of (J A J)ᴴ
    const b = try ops.alloc(T, allocator, n, m);
    defer ops.free(T, allocator, b);
    flip(T, true, am, b);
    const th = try allocator.alloc(T, k);
    defer allocator.free(th);
    try lq.factor(allocator, io, T, b, th.ptr);
    flip(T, true, b, am);
    reverse(T, th, tau[0..k]);
}

/// `?ORGRQ`, `?UNGRQ`: the `m x n` `Q` with orthonormal rows, the last `m`
/// of `H(0)ᴴ ... H(k-1)ᴴ`, over the reflections `gerqf` left in `A`'s last
/// `k` rows.
pub fn orgrq(allocator: Allocator, io: Io, comptime T: type, m: usize, n: usize, k: usize, a: []T, lda: usize, tau: []const T) OrgrqError!void {
    comptime ops.requireFloat(T, "orgrq");
    if (m > n) return error.TooManyRows;
    if (k > m) return error.TooManyReflections;
    try util.holds(a.len, m, n, lda);
    if (tau.len < k) return error.SliceTooShort;
    if (m == 0) return;
    const am: Mat(T) = .{ .ptr = a.ptr, .rows = m, .cols = n, .ld = lda };
    const w = try ops.alloc(T, allocator, m, n);
    defer ops.free(T, allocator, w);
    flip(T, false, am, w);
    const th = try allocator.alloc(T, k);
    defer allocator.free(th);
    reverse(T, tau[0..k], th);
    try lq.orglq(allocator, io, T, m, n, k, w.slice(), w.ld, th);
    flip(T, false, w, am);
}

/// `?ORMRQ`, `?UNMRQ`: `C := op(Q) C` (side `l`) or `C op(Q)` (side `r`),
/// `op` by `trans` (`t` and `c` both `Qᴴ`), for the `m x n` `C` and the
/// `Q` of `k` reflections in `A`'s `k` rows (`m` columns for `l`, `n` for
/// `r`), as `gerqf` left them in its last `k`. Only the part of `A` left
/// of each reflection's unit is read.
pub fn ormrq(allocator: Allocator, io: Io, comptime T: type, side: Side, trans: Trans, m: usize, n: usize, k: usize, a: []const T, lda: usize, tau: []const T, c: []T, ldc: usize) OrmrqError!void {
    comptime ops.requireFloat(T, "ormrq");
    try interface.takesAdjoint(T, trans);
    const nq = if (side == .l) m else n;
    if (k > nq) return error.TooManyReflections;
    try util.holds(a.len, k, nq, lda);
    try util.holds(c.len, m, n, ldc);
    if (tau.len < k) return error.SliceTooShort;
    if (m == 0 or n == 0 or k == 0) return;
    const am: Mat(T) = .{ .ptr = @constCast(a.ptr), .rows = k, .cols = nq, .ld = lda };
    const cm: Mat(T) = .{ .ptr = c.ptr, .rows = m, .cols = n, .ld = ldc };
    const w = try ops.alloc(T, allocator, k, nq);
    defer ops.free(T, allocator, w);
    flip(T, false, am, w);
    const th = try allocator.alloc(T, k);
    defer allocator.free(th);
    reverse(T, tau[0..k], th);
    // op(J Q̂ J) C = J op(Q̂) (J C); C op(J Q̂ J) = (C J) op(Q̂) J
    reverseAlong(T, side, cm);
    try lq.ormlq(allocator, io, T, side, trans, m, n, k, w.slice(), w.ld, th, c, ldc);
    reverseAlong(T, side, cm);
}
