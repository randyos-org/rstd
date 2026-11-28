//! `?ORMTR`, `?UNMTR`: multiplying by `sytrd`'s `Q`, Reference-LAPACK
//! 3.12.0 (`lapack-netlib/SRC/dormtr.f`, `zunmtr.f`, OpenBLAS
//! v0.3.34-794-gdbf2c24c9; modified BSD), over `ql.zig`'s `ormql` for the
//! upper triangle and `qr.zig`'s `ormqr` for the lower, as the reference.

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const interface = @import("interface.zig");
const ops = @import("ops.zig");
const qr = @import("qr.zig");
const ql = @import("ql.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Side = types.Side;
const Trans = types.Trans;
const Uplo = types.Uplo;

/// What `ormtr` returns: a reversed `uplo`, a `trans` other than `n`, `t`
/// or `c` (complex: `n` or `c`), a slice too short or a leading dimension
/// too small, and the allocator's and cancellation.
pub const OrmtrError = interface.UploError || qr.OrmqrError;

/// `?ORMTR`, `?UNMTR` (`dormtr.f`, `zunmtr.f`): `C := op(Q) C` (side `l`)
/// or `C op(Q)` (side `r`), `op` by `trans` (`t` and `c` both `Qᴴ`), for
/// the `m x n` `C` and the `Q` of `sytrd`'s reduction of an order `nq`
/// matrix (`m` for `l`, `n` for `r`): its reflections in `a` (`nq x nq`) on
/// `uplo`'s side of the diagonal, their factors in `tau` (`nq - 1`).
pub fn ormtr(allocator: Allocator, io: Io, comptime T: type, side: Side, uplo: Uplo, trans: Trans, m: usize, n: usize, a: []const T, lda: usize, tau: []const T, c: []T, ldc: usize) OrmtrError!void {
    comptime ops.requireFloat(T, "ormtr");
    if (uplo.isRev()) return error.UnsupportedUplo;
    try interface.takesAdjoint(T, trans);
    const nq = if (side == .l) m else n;
    try util.holds(a.len, nq, nq, lda);
    try util.holds(c.len, m, n, ldc);
    if (tau.len < nq -| 1) return error.SliceTooShort;
    if (m == 0 or n == 0 or nq == 1) return;
    // `Q` acts on the last `nq - 1` rows (columns) for `.l`, the first for `.u`
    const mi = if (side == .l) m - 1 else m;
    const ni = if (side == .l) n else n - 1;
    switch (uplo) {
        // `dormtr.f`: `DORMQL` on `A(1, 2)` and `C(1, 1)`
        .u => try ql.ormql(allocator, io, T, side, trans, mi, ni, nq - 1, a[lda..], lda, tau, c, ldc),
        // `DORMQR` on `A(2, 1)` and `C(2, 1)` or `C(1, 2)`
        .l => try qr.ormqr(allocator, io, T, side, trans, mi, ni, nq - 1, a[1..], lda, tau, c[if (side == .l) 1 else ldc..], ldc),
        .v, .m => unreachable,
    }
}
