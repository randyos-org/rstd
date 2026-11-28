//! `?ORMQL`, behind Reference-LAPACK's arguments (`dormql.f`, `zunmql.f`,
//! 3.12.0 as OpenBLAS v0.3.34-794-gdbf2c24c9 ships it; modified BSD), for
//! `ormtr`'s upper triangle, over `qr.zig`'s `ormqr`: with `J` the reversal,
//! a QL's `Q = H(k-1) ... H(0)` is `J Q̂ J` for the QR-stored `Q̂` of its
//! reflections flipped (`reversal.zig`), their `tau`s reversed, as `rq.zig`
//! does for RQ over LQ.

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const interface = @import("interface.zig");
const ops = @import("ops.zig");
const qr = @import("qr.zig");
const reversal = @import("reversal.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = ops.Mat;
const Side = types.Side;
const Trans = types.Trans;

/// What `ormql` returns: `ormqr`'s.
pub const OrmqlError = qr.OrmqrError;

/// `?ORMQL`, `?UNMQL`: `C := op(Q) C` (side `l`) or `C op(Q)` (side `r`),
/// `op` by `trans` (`t` and `c` both `Qᴴ`), for the `m x n` `C` and the
/// `Q = H(k-1) ... H(0)` of `k` reflections in `A`'s `k` columns (`m` rows
/// for `l`, `n` for `r`), as `geqlf` leaves them: reflection `i`'s unit at
/// row `nq - k + i`, its vector above it in column `i`. Only that part of
/// `A` is read.
pub fn ormql(allocator: Allocator, io: Io, comptime T: type, side: Side, trans: Trans, m: usize, n: usize, k: usize, a: []const T, lda: usize, tau: []const T, c: []T, ldc: usize) OrmqlError!void {
    comptime ops.requireFloat(T, "ormql");
    try interface.takesAdjoint(T, trans);
    const nq = if (side == .l) m else n;
    if (k > nq) return error.TooManyReflections;
    try util.holds(a.len, nq, k, lda);
    try util.holds(c.len, m, n, ldc);
    if (tau.len < k) return error.SliceTooShort;
    if (m == 0 or n == 0 or k == 0) return;
    const am: Mat(T) = .{ .ptr = @constCast(a.ptr), .rows = nq, .cols = k, .ld = lda };
    const cm: Mat(T) = .{ .ptr = c.ptr, .rows = m, .cols = n, .ld = ldc };
    const w = try ops.alloc(T, allocator, nq, k);
    defer ops.free(T, allocator, w);
    reversal.flip(T, false, am, w);
    const th = try allocator.alloc(T, k);
    defer allocator.free(th);
    reversal.reverse(T, tau[0..k], th);
    // op(J Q̂ J) C = J op(Q̂) (J C); C op(J Q̂ J) = (C J) op(Q̂) J
    reversal.reverseAlong(T, side, cm);
    try qr.ormqr(allocator, io, T, side, trans, m, n, k, w.slice(), w.ld, th, c, ldc);
    reversal.reverseAlong(T, side, cm);
}
