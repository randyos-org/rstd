//! ReLAPACK's recursive `pbtrf` (see `relapack.zig` for the source).

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const driver = @import("../blas/blas3/driver.zig");
const ops = @import("ops.zig");
const interface = @import("interface.zig");
const relapack = @import("relapack.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Uplo = types.Uplo;
const isComplex = types.isComplex;
const scalar = ops.scalar;
const Error = driver.Error;
const crossover = relapack.crossover;
const ctrans = relapack.ctrans;
const gemm = relapack.gemm;
const potrfRec = relapack.potrfRec;
const split = relapack.split;
const Tview = relapack.Tview;
const syrk = relapack.syrk;
const trsm = relapack.trsm;

/// What `pbtrf` returns: bad arguments in the reference's order,
/// `NotPositiveDefinite` (with the minor, as `potrf`'s), `error.OutOfMemory`
/// or `error.Canceled`.
/// What `pbtrf` returns: `NotPositiveDefinite` as `potrf`, a reversed
/// `uplo`, a slice too short or a leading dimension under `kd + 1`, and the
/// allocator's and cancellation.
pub const PbtrfError = interface.TriangleError || Allocator.Error || Io.Cancelable || error{NotPositiveDefinite};

/// `?PBTRF` by ReLAPACK's recursion (`relapack/src/dpbtrf.c`, `zpbtrf.c`
/// and the rest; OpenBLAS v0.3.34-317): the Cholesky factor of the Hermitian
/// (real: symmetric) positive-definite band matrix with `kd` off-diagonals
/// in LAPACK's band storage `ab` (`ldab >= kd + 1`; `u`: `A(i, j)` at row
/// `kd + i - j` of column `j`, `l`: at row `i - j`), into the same storage.
/// Blocks of up to `kd` columns are `potrf`'s (this file's), the part below
/// them a `trsm` and `herk`, the triangle past the band's edge through a
/// workspace, as ReLAPACK; the base case is Reference-LAPACK's `xPBTF2`
/// (`pbtf2`). A failure is `potrf`'s: a leading minor that is not positive
/// definite (NaN included) is `error.NotPositiveDefinite`.
///
/// The reference wins where ReLAPACK factors the rest as a dense matrix
/// whenever `kd` exceeds the block (`if (*kd > n1)`), which reads and
/// writes outside the band once the rest is wider than `kd + 1`; there it
/// recurses on the band. The workspace is `n1 x n1` for the top level's
/// block `n1`, which every level's `n22 x n1` fits (`dpbtrf.c`'s `n22`;
/// `zpbtrf.c` bounds it by `kd`, past the block).
pub fn pbtrf(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, kd: usize, ab: []T, ldab: usize) PbtrfError!void {
    comptime ops.requireFloat(T, "pbtrf");
    if (uplo.isRev()) return error.UnsupportedUplo;
    try util.holds(ab.len, kd + 1, n, ldab);
    if (n == 0) return;
    const n1 = @min(split(T, n), kd);
    // zero, as `?laset` leaves the triangle the copies never write
    const work = try allocator.alloc(T, n1 * n1);
    defer allocator.free(work);
    @memset(work, scalar(T, 0));
    if (try pbtrfRec(allocator, io, T, uplo, n, kd, ab.ptr, ldab, work.ptr, @max(1, n1)) != null) return error.NotPositiveDefinite;
}

/// `RELAPACK_?pbtrf_rec`: the first pivot that is not positive (0-based),
/// or null.
fn pbtrfRec(al: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, kd: usize, ab: [*]T, ldab: usize, work: [*]T, ldw: usize) Error!?usize {
    if (n <= crossover or ldab == 1) return pbtf2(T, uplo == .u, n, kd, ab, ldab);
    const one = scalar(T, 1);
    const minus = scalar(T, -1);
    const ct = ctrans(T);
    // the band as a dense matrix of leading dimension `ldab - 1`
    const ld = ldab - 1;
    const a = ab + (if (uplo == .l) 0 else kd);
    const n1 = @min(split(T, n), kd);
    const n2 = n - n1;
    if (try potrfRec(al, io, T, uplo, n1, a, ld)) |k| return k;
    const n21 = @min(n2, kd - n1);
    const n22 = @min(n2 - n21, n1);
    const a_br = a + ld * n1 + n1;
    if (uplo == .l) {
        const blt = a + n1;
        const blb = blt + n21;
        // A_BLt = A_BLt / A_TL'; A_BRtl -= A_BLt A_BLt'
        try trsm(al, io, T, true, .r, .l, ct, .n, n21, n1, one, a, ld, blt, ld);
        try syrk(al, io, T, true, .l, .n, n21, n1, minus, blt, ld, null, 0, one, a_br, ld);
        // the triangle of A_BLb in the band, through the workspace
        copyTri(T, true, n22, n1, blb, ld, work, ldw);
        try trsm(al, io, T, true, .r, .l, ct, .n, n22, n1, one, a, ld, work, ldw);
        try gemm(al, io, T, .n, ct, n22, n21, n1, minus, work, ldw, blt, ld, one, a_br + n21, ld);
        try syrk(al, io, T, true, .l, .n, n22, n1, minus, work, ldw, null, 0, one, a_br + ld * n21 + n21, ld);
        copyTri(T, true, n22, n1, work, ldw, blb, ld);
    } else {
        const trl = a + ld * n1;
        const trr = trl + ld * n21;
        // A_TRl = A_TL' \ A_TRl; A_BRtl -= A_TRl' A_TRl
        try trsm(al, io, T, true, .l, .u, ct, .n, n1, n21, one, a, ld, trl, ld);
        try syrk(al, io, T, true, .u, ct, n21, n1, minus, trl, ld, null, 0, one, a_br, ld);
        copyTri(T, false, n1, n22, trr, ld, work, ldw);
        try trsm(al, io, T, true, .l, .u, ct, .n, n1, n22, one, a, ld, work, ldw);
        try gemm(al, io, T, ct, .n, n21, n22, n1, minus, trl, ld, work, ldw, one, a_br + ld * n21, ld);
        try syrk(al, io, T, true, .u, ct, n22, n1, minus, work, ldw, null, 0, one, a_br + ld * n21 + n21, ld);
        copyTri(T, false, n1, n22, work, ldw, trr, ld);
    }
    // the rest, dense where it lies inside the band, else banded
    const k2 = if (kd > n1 and n2 <= kd + 1)
        try potrfRec(al, io, T, uplo, n2, a_br, ld)
    else
        try pbtrfRec(al, io, T, uplo, n2, kd, ab + ldab * n1, ldab, work, ldw);
    return if (k2) |k| k + n1 else null;
}

/// `?LACPY` of the upper (`upper`) or lower triangle of `m x n`.
fn copyTri(comptime T: type, upper: bool, m: usize, n: usize, src: [*]const T, lds: usize, dst: [*]T, ldd: usize) void {
    for (0..n) |j| {
        const lo = if (upper) 0 else @min(j, m);
        const hi = if (upper) @min(j + 1, m) else m;
        for (lo..hi) |i| dst[i + j * ldd] = src[i + j * lds];
    }
}

/// Reference-LAPACK's `xPBTF2` (`lapack-netlib/SRC/dpbtf2.f`, `zpbtf2.f`;
/// OpenBLAS v0.3.34-317): column by column, `xSYR` or `xHER` on the
/// `kd`-wide window below (`u`: right of) the pivot, in the reference
/// BLAS's order. `upper` runs the lower algorithm on the transpose (the
/// reference conjugates the row it takes). A pivot that is not positive
/// fails, NaN too, as `potf2`'s (the reference's `AJJ.LE.ZERO` passes NaN).
fn pbtf2(comptime T: type, upper: bool, n: usize, kd: usize, ab: [*]T, ldab: usize) ?usize {
    const herm = comptime isComplex(T);
    const a: Tview(T) = .{ .p = ab + (if (upper) kd else 0), .ld = ldab - 1, .t = upper };
    for (0..n) |j| {
        const ajj0 = ops.re(T, a.at(j, j).*);
        if (!(ajj0 > 0)) {
            a.at(j, j).* = scalar(T, ajj0);
            return j;
        }
        const ajj = @sqrt(ajj0);
        a.at(j, j).* = scalar(T, ajj);
        const kn = @min(kd, n - j - 1);
        if (kn == 0) continue;
        const r = 1 / ajj;
        for (j + 1..j + 1 + kn) |i| a.at(i, j).* = ops.mulReal(T, a.at(i, j).*, r);
        // `xSYR`, `xHER` with -1 on the window
        for (j + 1..j + 1 + kn) |c| {
            const x = a.at(c, j).*;
            if (ops.isZero(T, x)) {
                if (comptime herm) a.at(c, c).* = scalar(T, a.at(c, c).re);
                continue;
            }
            const t = ops.neg(T, ops.conj(T, x));
            if (comptime herm) a.at(c, c).* = scalar(T, a.at(c, c).re + ops.mul(T, x, t).re);
            for (if (herm) c + 1 else c..j + 1 + kn) |i| a.at(i, c).* = ops.add(T, a.at(i, c).*, ops.mul(T, a.at(i, j).*, t));
        }
    }
    return null;
}
