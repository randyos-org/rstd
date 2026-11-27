//! `?PSTRF`: the Cholesky factorization with complete pivoting of a
//! positive semidefinite matrix, Reference-LAPACK 3.12.0
//! (`lapack-netlib/SRC/dpstrf.f`, `zpstrf.f`, `dpstf2.f`, `zpstf2.f`,
//! OpenBLAS v0.3.34-317; modified BSD). One routine for both: `?PSTF2` is
//! `?PSTRF`'s arithmetic in a single block, which `?PSTRF` itself runs
//! when its block covers the matrix.

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const interface = @import("interface.zig");
const ops = @import("ops.zig");
const blas_calls = @import("blas_calls.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const Mat = ops.Mat;

/// What `pstrf` returns: a slice too short or a leading dimension too
/// small. A rank below `n` is no error (the reference's `INFO = 1`).
pub const PstrfError = interface.ShapeError || Allocator.Error || Io.Cancelable;

/// `ILAENV(1, 'DPOTRF')`, `?PSTRF`'s block (`ilaenv.f`: 64).
const nb = 64;

/// `?PSTRF` (`dpstrf.f`, `zpstrf.f`): `Pᵀ A P = Uᴴ U` (`uplo` `.u`) or `L
/// Lᴴ` (`.l`) for the `n x n` positive semidefinite (Hermitian) `A` in
/// `a`'s `uplo` triangle, the other not read, the factor over it. Returns
/// the computed rank `r`: the factor's leading `r` columns (rows) are
/// complete and the trailing block holds what was left, as the reference
/// leaves it. `piv` (`n`) is the permutation, 0-based: column `k` of `P`
/// is `e(piv[k])`. A diagonal pivot at or below `tol` stops it; a negative
/// `tol` is `n eps max(diag(A))`, `DLAMCH('Epsilon')` being half of `eps`.
pub fn pstrf(allocator: Allocator, io: Io, comptime T: type, uplo: types.Uplo, n: usize, a: []T, lda: usize, piv: []usize, tol: Float(T)) PstrfError!usize {
    comptime ops.requireFloat(T, "pstrf");
    try util.holds(a.len, n, n, lda);
    if (piv.len < n) return error.SliceTooShort;
    if (n == 0) return 0;
    const work = try allocator.alloc(Float(T), 2 * n);
    defer allocator.free(work);
    const m: Mat(T) = .{ .ptr = a.ptr, .rows = n, .cols = n, .ld = lda };
    return factor(allocator, io, T, uplo == .u, m, piv[0..n], tol, work, if (nb >= n) n else nb);
}

/// `pstrf` after its checks, in blocks of `block` columns.
fn factor(al: Allocator, io: Io, comptime T: type, upper: bool, a: Mat(T), piv: []usize, tol: Float(T), work: []Float(T), block: usize) PstrfError!usize {
    const F = Float(T);
    const n = a.rows;
    const re = ops.re;
    const conj = ops.conj;
    for (piv, 0..) |*p, i| p.* = i;
    // the largest diagonal entry, the first pivot
    var pvt: usize = 0;
    var ajj = re(T, a.get(0, 0));
    for (1..n) |i| if (re(T, a.get(i, i)) > ajj) {
        pvt = i;
        ajj = re(T, a.get(i, i));
    };
    if (!(ajj > 0)) return 0;
    const dstop = if (tol < 0) @as(F, @floatFromInt(n)) * (std.math.floatEps(F) / 2) * ajj else tol;
    // `A(i, j)` of the factor's side: `U(i, j)` is `A(i, j)`, `L(i, j)` is `A(j, i)` stored lower
    const at = struct {
        fn f(m: Mat(T), up: bool, i: usize, j: usize) *T {
            return if (up) m.at(i, j) else m.at(j, i);
        }
    }.f;
    var k: usize = 0;
    while (k < n) {
        const jb = @min(block, n - k);
        @memset(work[k..n], 0);
        for (k..k + jb) |j| {
            // the partial pivots `A(i, i) - Σ |U(l, i)|²` over this block
            for (j..n) |i| {
                if (j > k) {
                    const u = at(a, upper, j - 1, i).*;
                    work[i] += re(T, ops.mul(T, conj(T, u), u));
                }
                work[n + i] = re(T, a.get(i, i)) - work[i];
            }
            if (j > 0) {
                // `MAXLOC`: the first largest
                pvt = j;
                ajj = work[n + j];
                for (j + 1..n) |i| if (work[n + i] > ajj) {
                    pvt = i;
                    ajj = work[n + i];
                };
                if (ajj <= dstop or std.math.isNan(ajj)) {
                    a.at(j, j).* = ops.scalar(T, ajj);
                    return j;
                }
            }
            if (j != pvt) {
                // swap row and column `j` with `pvt`, the triangle's way
                a.at(pvt, pvt).* = a.get(j, j);
                const ld = a.ld;
                if (upper) {
                    try blas_calls.swap(io, T, j, a.col(j), 1, a.col(pvt), 1);
                    if (pvt + 1 < n) try blas_calls.swap(io, T, n - 1 - pvt, a.ptr + j + (pvt + 1) * ld, ld, a.ptr + pvt + (pvt + 1) * ld, ld);
                } else {
                    try blas_calls.swap(io, T, j, a.ptr + j, ld, a.ptr + pvt, ld);
                    if (pvt + 1 < n) try blas_calls.swap(io, T, n - 1 - pvt, a.col(j) + pvt + 1, 1, a.col(pvt) + pvt + 1, 1);
                }
                for (j + 1..pvt) |i| {
                    const x = at(a, upper, j, i);
                    const y = at(a, upper, i, pvt);
                    const t = conj(T, x.*);
                    x.* = conj(T, y.*);
                    y.* = t;
                }
                const jp = at(a, upper, j, pvt);
                jp.* = conj(T, jp.*);
                std.mem.swap(F, &work[j], &work[pvt]);
                std.mem.swap(usize, &piv[j], &piv[pvt]);
            }
            ajj = @sqrt(ajj);
            a.at(j, j).* = ops.scalar(T, ajj);
            if (j + 1 < n) {
                // row (column) `j` of the factor: `(A(j, j+1:n) - U(k:j, j)ᴴ U(k:j, j+1:n)) / ajj`
                const rest = n - j - 1;
                if (upper) {
                    conjugate(T, a.col(j), 1, j);
                    try blas_calls.gemv(al, io, T, .t, ops.scalar(T, -1), a.sub(k, j + 1, j - k, rest), a.col(j) + k, 1, ops.one(T), a.ptr + j + (j + 1) * a.ld, a.ld);
                    conjugate(T, a.col(j), 1, j);
                    try blas_calls.scal(io, T, rest, 1 / ajj, a.ptr + j + (j + 1) * a.ld, a.ld);
                } else {
                    conjugate(T, a.ptr + j, a.ld, j);
                    try blas_calls.gemv(al, io, T, .n, ops.scalar(T, -1), a.sub(j + 1, k, rest, j - k), a.ptr + j + k * a.ld, a.ld, ops.one(T), a.col(j) + j + 1, 1);
                    conjugate(T, a.ptr + j, a.ld, j);
                    try blas_calls.scal(io, T, rest, 1 / ajj, a.col(j) + j + 1, 1);
                }
            }
        }
        // the trailing block's update, `?SYRK` (`?HERK`)
        const j = k + jb;
        if (j < n) {
            if (upper) {
                try blas_calls.gram(al, io, T, true, -1, a.sub(k, j, jb, n - j), 1, a.sub(j, j, n - j, n - j));
            } else try blas_calls.gram(al, io, T, false, -1, a.sub(j, k, n - j, jb), 1, a.sub(j, j, n - j, n - j));
        }
        k += jb;
    }
    return n;
}

/// `?LACGV`: `x := conj(x)`, `n` elements `inc` apart; nothing for real.
fn conjugate(comptime T: type, x: [*]T, inc: usize, n: usize) void {
    if (comptime !types.isComplex(T)) return;
    for (0..n) |i| x[i * inc] = x[i * inc].conjugate();
}
