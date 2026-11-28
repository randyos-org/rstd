//! The generalized QR and RQ factorizations, `?GGQRF` and `?GGRQF`:
//! Reference-LAPACK 3.12.0 (`lapack-netlib/SRC/dggqrf.f`, `zggqrf.f`,
//! `dggrqf.f`, `zggrqf.f`, OpenBLAS v0.3.34-794; modified BSD), over rstd's
//! `geqrf`, `ormqr`, `gerqf` and `ormrq`.

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const interface = @import("interface.zig");
const ops = @import("ops.zig");
const qr = @import("qr.zig");
const rq = @import("rq.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;

/// What `ggqrf` and `ggrqf` return: a slice too short (`a`, `b`, or a
/// `tau` under its count) or a leading dimension too small, and the
/// allocator's and cancellation.
pub const GgqrfError = interface.ShapeError || Allocator.Error || Io.Cancelable;
pub const GgrqfError = GgqrfError;

/// `?GGQRF`: for the `n x m` `A` and `n x p` `B`, `A = Q R` and `B = Q T
/// Z`: `A` gets `geqrf`'s factor (`taua`, `min(n, m)`), then `B := Qᴴ B`
/// gets `gerqf`'s (`taub`, `min(n, p)`).
pub fn ggqrf(allocator: Allocator, io: Io, comptime T: type, n: usize, m: usize, p: usize, a: []T, lda: usize, taua: []T, b: []T, ldb: usize, taub: []T) GgqrfError!void {
    comptime ops.requireFloat(T, "ggqrf");
    try util.holds(a.len, n, m, lda);
    try util.holds(b.len, n, p, ldb);
    if (taua.len < @min(n, m) or taub.len < @min(n, p)) return error.SliceTooShort;
    try qr.geqrf(allocator, io, T, n, m, a, lda, taua);
    qr.ormqr(allocator, io, T, .l, .c, n, p, @min(n, m), a, lda, taua, b, ldb) catch |e| return ops.narrow(GgqrfError, e);
    try rq.gerqf(allocator, io, T, n, p, b, ldb, taub);
}

/// `?GGRQF`: for the `m x n` `A` and `p x n` `B`, `A = R Q` and `B = Z T
/// Q`: `A` gets `gerqf`'s factor (`taua`, `min(m, n)`), then `B := B Qᴴ`
/// gets `geqrf`'s (`taub`, `min(p, n)`).
pub fn ggrqf(allocator: Allocator, io: Io, comptime T: type, m: usize, p: usize, n: usize, a: []T, lda: usize, taua: []T, b: []T, ldb: usize, taub: []T) GgrqfError!void {
    comptime ops.requireFloat(T, "ggrqf");
    try util.holds(a.len, m, n, lda);
    try util.holds(b.len, p, n, ldb);
    const k = @min(m, n);
    if (taua.len < k or taub.len < @min(p, n)) return error.SliceTooShort;
    try rq.gerqf(allocator, io, T, m, n, a, lda, taua);
    // the reflections are `A`'s last `k` rows
    if (k > 0 and p > 0) rq.ormrq(allocator, io, T, .r, .c, p, n, k, a[m - k ..], lda, taua, b, ldb) catch |e| return ops.narrow(GgrqfError, e);
    try qr.geqrf(allocator, io, T, p, n, b, ldb, taub);
}
