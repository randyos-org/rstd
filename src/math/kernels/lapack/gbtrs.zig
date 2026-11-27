//! `?GBTRS`: solves with `gbtrf`'s band factor, Reference-LAPACK 3.12.0
//! (`lapack-netlib/SRC/dgbtrs.f`, `zgbtrs.f`, OpenBLAS v0.3.34-317;
//! modified BSD).

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const interface = @import("interface.zig");
const ops = @import("ops.zig");
const blas_calls = @import("blas_calls.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Trans = types.Trans;
const Mat = ops.Mat;

/// What `gbtrs` returns: a `trans` it does not take, a pivot outside
/// `0..n - 1` (checked before `B` is touched; the reference checks none), a
/// slice too short or a leading dimension too small (`ldab` under `2 kl +
/// ku + 1`), and the allocator's and cancellation.
pub const GbtrsError = interface.TransError || interface.PivotError || interface.ShapeError || Allocator.Error || Io.Cancelable;

/// `?GBTRS` (`dgbtrs.f`, `zgbtrs.f`): solves `op(A) X = B` for the `n x
/// nrhs` `X`, in place of `B`, given `gbtrf`'s `P L U` of the band `A` (`kl`
/// sub-, `ku` superdiagonals) in `ab` and its `ipiv` (0-based). `trans` as
/// `getrs`'s: `n`, `t`, `r` (`conj(A)`, rstd's, by conjugating `B` around
/// the plain solve) or `c`. No singularity check, as in the reference.
pub fn gbtrs(allocator: Allocator, io: Io, comptime T: type, trans: Trans, n: usize, kl: usize, ku: usize, nrhs: usize, ab: []const T, ldab: usize, ipiv: []const usize, b: []T, ldb: usize) GbtrsError!void {
    comptime ops.requireFloat(T, "gbtrs");
    try interface.takes(trans);
    const rows = 2 * kl + ku + 1;
    if (ldab < rows) return error.LeadingDimensionTooSmall;
    try util.holds(ab.len, rows, n, ldab);
    try util.holds(b.len, n, nrhs, ldb);
    try interface.pivotsIn(ipiv, n);
    if (n == 0 or nrhs == 0) return;
    const cx = types.isComplex(T) and trans.conjA();
    const kd = kl + ku;
    const a = ab.ptr;
    const bm: Mat(T) = .{ .ptr = b.ptr, .rows = n, .cols = nrhs, .ld = ldb };
    const minus = ops.scalar(T, -1);
    if (!trans.transposed()) {
        if (cx) conjugate(T, bm);
        // `L⁻¹ P B`, one column of `L` at a time
        if (kl > 0) for (0..n - 1) |j| {
            const lm = @min(kl, n - 1 - j);
            const l = ipiv[j];
            if (l != j) try blas_calls.swap(io, T, nrhs, b.ptr + l, ldb, b.ptr + j, ldb);
            try blas_calls.ger(allocator, io, T, minus, @constCast(a + kd + 1 + j * ldab), 1, b.ptr + j, ldb, bm.sub(j + 1, 0, lm, nrhs));
        };
        for (0..nrhs) |i| try blas_calls.tbsv(allocator, io, T, true, .n, false, n, kd, a, ldab, bm.col(i), 1);
        if (cx) conjugate(T, bm);
        return;
    }
    // `op(U)⁻¹`, then `op(L)⁻¹` and the interchanges backward
    const op: Trans = if (cx) .c else .t;
    for (0..nrhs) |i| try blas_calls.tbsv(allocator, io, T, true, op, false, n, kd, a, ldab, bm.col(i), 1);
    if (kl > 0) {
        var j = n - 1;
        while (j > 0) {
            j -= 1;
            const lm = @min(kl, n - 1 - j);
            if (cx) conjugateRow(T, bm, j);
            try blas_calls.gemv(allocator, io, T, op, minus, bm.sub(j + 1, 0, lm, nrhs), @constCast(a + kd + 1 + j * ldab), 1, ops.one(T), b.ptr + j, ldb);
            if (cx) conjugateRow(T, bm, j);
            const l = ipiv[j];
            if (l != j) try blas_calls.swap(io, T, nrhs, b.ptr + l, ldb, b.ptr + j, ldb);
        }
    }
}

/// `B := conj(B)`.
fn conjugate(comptime T: type, m: Mat(T)) void {
    if (comptime !types.isComplex(T)) return;
    for (0..m.cols) |j| for (m.col(j)[0..m.rows]) |*x| {
        x.* = x.conjugate();
    };
}

/// `?LACGV` on row `i`.
fn conjugateRow(comptime T: type, m: Mat(T), i: usize) void {
    if (comptime !types.isComplex(T)) return;
    for (0..m.cols) |j| m.at(i, j).* = m.get(i, j).conjugate();
}
