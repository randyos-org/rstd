//! Cholesky factor, solve and inverse: OpenBLAS v0.3.34-317's `potrf`,
//! `potf2` and `potri` (`interface/lapack/potrf.c`, `potf2.c`, `potri.c` and
//! their `z` forms, `lapack/potrf/`, `lapack/potf2/`) with Reference-LAPACK's
//! `xPOTRF2` semantics (`lapack-netlib/SRC/zpotrf2.f`), and `xPOTRS`
//! (`zpotrs.f`). Departures and threading: design.md's LAPACK section.

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const driver = @import("../blas/blas3/driver.zig");
const accel = @import("../blas/apple_accel.zig");
const Vtable = @import("../blas/Vtable.zig");
const ops = @import("ops.zig");
const interface = @import("interface.zig");
const tri = @import("triangular.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Uplo = types.Uplo;
const Float = types.Float;
const isComplex = types.isComplex;
const scalar = ops.scalar;
const re = ops.re;
const abs2 = ops.abs2;
const add = ops.add;
const sub = ops.sub;
const mulConj = ops.mulConj;
const divReal = ops.divReal;

/// What `potrs` returns: a reversed `uplo`, a slice too short or a leading
/// dimension too small, and the allocator's and cancellation.
pub const PotrsError = interface.TriangleError || Allocator.Error || Io.Cancelable;

/// What `potrfChecked` returns.
pub const CheckedError = interface.TriangleError || Allocator.Error || Io.Cancelable;

/// What `potrf` returns: `NotPositiveDefinite` when a leading minor is not
/// positive definite, and `CheckedError`'s.
pub const PotrfError = CheckedError || error{NotPositiveDefinite};

/// What `potf2` returns, which allocates nothing.
pub const Potf2Error = interface.TriangleError || error{NotPositiveDefinite};

/// What `potri` returns: `Singular` when the factor has a zero on its
/// diagonal.
pub const PotriError = interface.TriangleError || Allocator.Error || Io.Cancelable || error{Singular};

// ---- potrf, potf2 ---------------------------------------------------------------

/// `?POTRF`: the Cholesky factor of the Hermitian (real: symmetric)
/// positive-definite `n x n` `A`, from and into its `uplo` triangle:
/// `A = Uᴴ U` for `u`, `A = L Lᴴ` for `l`. The other triangle and the
/// imaginary parts of the diagonal are never read; the factor's diagonal
/// comes out real and positive.
///
/// When a leading minor is not positive definite (a pivot not positive, or
/// NaN), returns `error.NotPositiveDefinite`, `A` partly factored: the
/// factor's columns (`u`: rows) before that pivot, the failed pivot on the
/// diagonal, partial updates elsewhere. The order of that minor is not
/// reported (issue #314); `potrfChecked` returns it.
///
/// Routed to Accelerate (macOS), a failure leaves its partial factor in `A`.
pub fn potrf(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, a: []T, lda: usize) PotrfError!void {
    if (try potrfChecked(allocator, io, T, uplo, n, a, lda) != null) return error.NotPositiveDefinite;
}

/// `accel.potrf`, or null where rstd runs it. A real upper one keeps the
/// strictly lower triangle across the call: Accelerate uses it as scratch on
/// Apple Silicon (issue #537).
fn potrfAccel(allocator: Allocator, comptime T: type, uplo: Uplo, n: usize, a: []T, lda: usize) Allocator.Error!?isize {
    const scratched = comptime accel.potrf_scratches_lower and !isComplex(T);
    if (!scratched or uplo != .u or n < 2 or !accel.routing()) return accel.potrf(T, uplo, @intCast(n), a.ptr, @intCast(lda));
    const kept = try allocator.alloc(T, n * (n - 1) / 2);
    defer allocator.free(kept);
    var at: usize = 0;
    for (0..n - 1) |j| {
        const len = n - 1 - j;
        @memcpy(kept[at..][0..len], a[j * lda + j + 1 ..][0..len]);
        at += len;
    }
    const info = accel.potrf(T, uplo, @intCast(n), a.ptr, @intCast(lda));
    at = 0;
    for (0..n - 1) |j| {
        const len = n - 1 - j;
        @memcpy(a[j * lda + j + 1 ..][0..len], kept[at..][0..len]);
        at += len;
    }
    return info;
}

/// `potrf`, a failure returned as the index (0-based) of the first pivot
/// that is not positive, or null: the leading minor one larger is not
/// positive definite, `A` partly factored as `potrf` leaves it. Not in
/// `rstd.linalg.lapack`: MATLAB's `[R, p] = chol(A)` `p`, for a caller wanting it.
pub fn potrfChecked(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, a: []T, lda: usize) CheckedError!?usize {
    comptime ops.requireFloat(T, "potrf");
    try interface.holdsTriangle(uplo, n, a.len, lda);
    if (n == 0) return null;
    // Accelerate's `INFO`, 1-based at its ABI
    if (try potrfAccel(allocator, T, uplo, n, a, lda)) |info| return if (info == 0) null else @intCast(info - 1);
    const vt = interface.table(T);
    const blk = driver.Blocking.of(T, vt);
    const upper = uplo == .u;
    return if (interface.threads(.potrf, T, n)) |w|
        try factorParallel(allocator, io, T, vt, upper, blk, n, a.ptr, lda, w)
    else
        try factor(allocator, io, T, vt, upper, interface.dtbOf(vt), blk.q, n, a.ptr, lda);
}

/// `?POTF2`: `potrf` unblocked (`potf2_U.c` and the rest, as loops), with
/// its errors.
pub fn potf2(comptime T: type, uplo: Uplo, n: usize, a: []T, lda: usize) Potf2Error!void {
    comptime ops.requireFloat(T, "potf2");
    try interface.holdsTriangle(uplo, n, a.len, lda);
    if (n == 0) return;
    if (unblocked(T, uplo == .u, n, a.ptr, lda) != null) return error.NotPositiveDefinite;
}

/// `potrf_U_single.c` and `potrf_L_single.c`: the index of the first pivot
/// that is not positive, or null. Upstream fuses the panel's `trsm`
/// with the update's `herk` column block by column block; here each is one
/// driver call, the same arithmetic.
fn factor(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, upper: bool, dtb: usize, q: usize, n: usize, a: [*]T, lda: usize) driver.Error!?usize {
    if (n <= dtb / 2) return unblocked(T, upper, n, a, lda);
    // `(n + 3) / 4` and `n / 4` as upstream has them
    const blocking = if (n > 4 * q) q else if (upper) (n + 3) / 4 else @max(1, n / 4);
    var j: usize = 0;
    while (j < n) : (j += blocking) {
        const bk = @min(n - j, blocking);
        const d = a + j * (lda + 1);
        if (try factor(allocator, io, T, vt, upper, dtb, q, bk, d, lda)) |k| return k + j;
        try trailing(allocator, io, T, vt, upper, n - j - bk, bk, d, lda, 1);
    }
    return null;
}

/// `potrf_U_parallel.c` and `potrf_L_parallel.c`: blocks of half the order
/// rounded up to `unroll_n`, at most `Q`, each factored recursively, until
/// `4 * unroll_n` is left for `factor`; the drivers take up to `want`
/// workers. Upstream drops to `factor` for one thread too; here only the
/// size does (`threads`).
pub fn factorParallel(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, upper: bool, blk: driver.Blocking, n: usize, a: [*]T, lda: usize, want: u32) driver.Error!?usize {
    if (n <= 4 * blk.un) return factor(allocator, io, T, vt, upper, interface.dtbOf(vt), blk.q, n, a, lda);
    const blocking = @min((n / 2 + blk.un - 1) / blk.un * blk.un, blk.q);
    var i: usize = 0;
    while (i < n) : (i += blocking) {
        const bk = @min(n - i, blocking);
        const d = a + i * (lda + 1);
        if (try factorParallel(allocator, io, T, vt, upper, blk, bk, d, lda, want)) |k| return k + i;
        try trailing(allocator, io, T, vt, upper, n - i - bk, bk, d, lda, want);
    }
    return null;
}

/// The `rest` columns (`u`: rows) past the `bk` just factored at `d`:
/// `A12 := U11^-H A12; A22 -= A12^H A12`, or `A21 := A21 L11^-H;
/// A22 -= A21 A21^H` (`TRSM_LCUN`, `HERK_UC`; `TRSM_RCLN`, `HERK_LN`).
fn trailing(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, upper: bool, rest: usize, bk: usize, d: [*]T, lda: usize, want: u32) driver.Error!void {
    if (rest == 0) return;
    const complex = comptime isComplex(T);
    const one = scalar(T, 1);
    const minus = scalar(T, -1);
    const p = if (upper) d + bk * lda else d + bk;
    const c = d + bk * (lda + 1);
    if (upper) {
        try driver.trsm(allocator, io, T, vt, true, true, true, complex, false, bk, rest, one, d, lda, p, lda, want);
    } else {
        try driver.trsm(allocator, io, T, vt, false, false, true, complex, false, rest, bk, one, d, lda, p, lda, want);
    }
    try driver.syrk(allocator, io, T, vt, complex, upper, upper, rest, bk, minus, p, lda, null, 0, one, c, lda, want);
}

/// `?POTRS`: solves `A X = B` for the `n x nrhs` `X`, in place of `B`,
/// given `potrf`'s factor of `A` in its `uplo` triangle: `Uᴴ U X = B` or
/// `L Lᴴ X = B`, as two `trsm`s. The other triangle is never read.
pub fn potrs(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, nrhs: usize, a: []const T, lda: usize, b: []T, ldb: usize) PotrsError!void {
    comptime ops.requireFloat(T, "potrs");
    try interface.holdsTriangle(uplo, n, a.len, lda);
    try util.holds(b.len, n, nrhs, ldb);
    if (n == 0 or nrhs == 0) return;
    const vt = interface.table(T);
    const complex = comptime isComplex(T);
    const one = scalar(T, 1);
    const upper = uplo == .u;
    // `u`: `Uᴴ Y = B`, then `U X = Y`; `l`: `L Y = B`, then `Lᴴ X = Y`
    try driver.trsm(allocator, io, T, vt, true, upper, upper, upper and complex, false, n, nrhs, one, a.ptr, lda, b.ptr, ldb, 1);
    try driver.trsm(allocator, io, T, vt, true, upper, !upper, !upper and complex, false, n, nrhs, one, a.ptr, lda, b.ptr, ldb, 1);
}

/// `potf2_U.c`, `potf2_L.c` and their `z` forms: the unblocked factor, with
/// upstream's `dot`, `gemv` and `scal` as loops. A pivot that is not
/// positive is stored, and its index (0-based) returned; else null.
pub fn unblocked(comptime T: type, upper: bool, n: usize, a: [*]T, lda: usize) ?usize {
    const F = Float(T);
    for (0..n) |j| {
        const djj = &a[j * (lda + 1)];
        // `upper`: column `j` above the diagonal; else row `j` left of it
        const s: usize = if (upper) 1 else lda;
        const v = if (upper) a + j * lda else a + j;
        var sum: F = 0;
        for (0..j) |i| sum += abs2(T, v[i * s]);
        const ajj = re(T, djj.*) - sum;
        // `!(ajj > 0)`: NaN fails too, as `xPOTRF2`'s `DISNAN`
        if (!(ajj > 0)) {
            djj.* = scalar(T, ajj);
            return j;
        }
        const d = @sqrt(ajj);
        djj.* = scalar(T, d);
        if (upper) {
            // row `j` right of the diagonal: `(A(j, k) - U(:j, j)^H U(:j, k)) / d`
            for (j + 1..n) |k| {
                const c = a + k * lda;
                var t: T = scalar(T, 0);
                for (0..j) |i| t = add(T, t, mulConj(T, v[i], c[i]));
                c[j] = divReal(T, sub(T, c[j], t), d);
            }
        } else {
            // column `j` below the diagonal: `(A(k, j) - L(k, :j) L(j, :j)^H) / d`
            const c = a + j * lda;
            for (0..j) |i| {
                const l = a + i * lda;
                const x = v[i * lda];
                for (j + 1..n) |k| c[k] = sub(T, c[k], mulConj(T, x, l[k]));
            }
            for (j + 1..n) |k| c[k] = divReal(T, c[k], d);
        }
    }
    return null;
}

// ---- potri ---------------------------------------------------------------------

/// `?POTRI`: `A^-1` of the Hermitian (real: symmetric) positive-definite
/// `A` from `potrf`'s factor in its `uplo` triangle, into that triangle:
/// `U^-1 U^-H` or `L^-H L^-1`, as `trtri` (non-unit) then `lauum`. The other
/// triangle is never read.
///
/// A zero on the factor's diagonal returns `error.Singular` before `A` is
/// touched, as the reference's `xTRTRI` fails; upstream's `potri.c` calls `trtri`'s algorithm without that
/// check and divides by the zero.
pub fn potri(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, a: []T, lda: usize) PotriError!void {
    comptime ops.requireFloat(T, "potri");
    try interface.holdsTriangle(uplo, n, a.len, lda);
    if (n == 0) return;
    if (tri.firstZero(T, n, a.ptr, lda) != null) return error.Singular;
    // upstream runs both parallel algorithms, or both single ones
    const w = interface.threads(.potri, T, n);
    try tri.inverse(allocator, io, T, uplo == .u, false, n, a.ptr, lda, w);
    try tri.product(allocator, io, T, uplo == .u, n, a.ptr, lda, w);
}
