//! `?PPTRF` and `?PPTRS`: Cholesky in packed storage, Reference-LAPACK
//! 3.12.0 (`lapack-netlib/SRC/dpptrf.f`, `zpptrf.f`, `dpptrs.f`,
//! `zpptrs.f`, OpenBLAS v0.3.34-794; modified BSD).

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const interface = @import("interface.zig");
const ops = @import("ops.zig");
const blas_calls = @import("blas_calls.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Trans = types.Trans;
const Uplo = types.Uplo;
const Float = types.Float;

/// What `pptrf` returns: `NotPositiveDefinite` when a leading minor is not
/// (a pivot not positive, or NaN), a `uplo` that is neither triangle, a
/// slice too short, and the allocator's and cancellation.
pub const PptrfError = error{NotPositiveDefinite} || interface.UploError || interface.ShapeError || Allocator.Error || Io.Cancelable;

/// `?PPTRF` (`dpptrf.f`, `zpptrf.f`): `A = Uᴴ U` (`uplo` `.u`) or `L Lᴴ`
/// for the Hermitian (real: symmetric) positive definite `n x n` `A`,
/// packed by columns in `ap` (`n (n + 1) / 2`), the factor over it, its
/// diagonal real. On failure, `ap` holds the factor's columns before the
/// pivot and the failed pivot. NaN fails too, as `potrf`'s.
pub fn pptrf(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, ap: []T) PptrfError!void {
    comptime ops.requireFloat(T, "pptrf");
    if (uplo.isRev()) return error.UnsupportedUplo;
    try util.holdsPacked(ap.len, n);
    const re = ops.re;
    if (uplo == .u) {
        var jj: usize = 0;
        for (0..n) |j| {
            // column `j` starts at `jc`, its diagonal at `jj`
            const jc = jj;
            jj += j;
            const col = ap.ptr + jc;
            if (j > 0) try blas_calls.tpsv(allocator, io, T, true, if (comptime types.isComplex(T)) .c else .t, false, j, ap.ptr, col);
            const d = try blas_calls.dot(io, T, true, j, col, 1, col, 1);
            const ajj = re(T, ap[jj]) - re(T, d);
            if (!(ajj > 0)) return fail(T, ap, jj, ajj);
            ap[jj] = ops.scalar(T, @sqrt(ajj));
            jj += 1;
        }
        return;
    }
    var jj: usize = 0;
    for (0..n) |j| {
        const ajj = re(T, ap[jj]);
        if (!(ajj > 0)) return fail(T, ap, jj, ajj);
        const s = @sqrt(ajj);
        ap[jj] = ops.scalar(T, s);
        const rest = n - 1 - j;
        if (rest > 0) {
            try blas_calls.scal(io, T, rest, 1 / s, ap.ptr + jj + 1, 1);
            try blas_calls.packedRank1(allocator, io, T, false, rest, -1, ap.ptr + jj + 1, ap.ptr + jj + rest + 1);
        }
        jj += rest + 1;
    }
}

/// The failed pivot stored, as the reference leaves it.
fn fail(comptime T: type, ap: []T, jj: usize, ajj: Float(T)) PptrfError {
    ap[jj] = ops.scalar(T, ajj);
    return error.NotPositiveDefinite;
}

/// What `pptrs` returns: a `uplo` that is neither triangle, a slice too
/// short or a leading dimension too small, and the allocator's and
/// cancellation.
pub const PptrsError = interface.UploError || interface.ShapeError || Allocator.Error || Io.Cancelable;

/// `?PPTRS` (`dpptrs.f`, `zpptrs.f`): solves `A X = B` for the `n x nrhs`
/// `X`, in place of `B`, given `pptrf`'s packed factor in `ap`: two `tpsv`
/// solves per column of `B`.
pub fn pptrs(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, nrhs: usize, ap: []const T, b: []T, ldb: usize) PptrsError!void {
    comptime ops.requireFloat(T, "pptrs");
    if (uplo.isRev()) return error.UnsupportedUplo;
    try util.holdsPacked(ap.len, n);
    try util.holds(b.len, n, nrhs, ldb);
    if (n == 0 or nrhs == 0) return;
    const upper = uplo == .u;
    const adjoint: Trans = if (comptime types.isComplex(T)) .c else .t;
    for (0..nrhs) |j| {
        const x = b.ptr + j * ldb;
        // `Uᴴ Y = B`, then `U X = Y`; or `L Y = B`, then `Lᴴ X = Y`
        try blas_calls.tpsv(allocator, io, T, upper, if (upper) adjoint else .n, false, n, ap.ptr, x);
        try blas_calls.tpsv(allocator, io, T, upper, if (upper) .n else adjoint, false, n, ap.ptr, x);
    }
}
