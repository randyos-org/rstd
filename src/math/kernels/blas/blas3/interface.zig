//! The level-3 interface: the checks and quick returns of `interface/gemm.c`,
//! `gemmt.c`, `symm.c`, `syrk.c`, `syr2k.c` and `trsm.c`, in `xerbla`'s order;
//! leading dimensions against `max(1, rows)`, as the reference BLAS does.

const std = @import("std");
const blas = @import("../__root__.zig");
const driver = @import("driver.zig");
const driver3m = @import("driver3m.zig");
const threaded = @import("threaded.zig");
const types = @import("../types.zig");
const util = @import("../util.zig");
const accel = @import("../apple_accel.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const isComplex = types.isComplex;
const Float = types.Float;
const Trans = types.Trans;
const Uplo = types.Uplo;
const Diag = types.Diag;
const Side = types.Side;

fn conjugated(comptime T: type, t: Trans) bool {
    return isComplex(T) and t.conjA();
}

/// The reference's `xerbla` argument, by name.
pub const ArgError = error{
    InvalidUplo,
    InvalidTrans,
    InvalidM,
    InvalidN,
    InvalidK,
    InvalidLda,
    InvalidLdb,
    InvalidLdc,
};
/// What a level-3 routine returns.
pub const Error = ArgError || Allocator.Error || Io.Cancelable;

const table = blas.tableOf;

/// `?GEMM`: `C := alpha*op(A)*op(B) + beta*C`, `op(A)` being `m x k`,
/// `op(B)` `k x n`. `C` is scaled by `beta` first (a zero `beta` never reads
/// it), and `A` and `B` are not read for a zero `alpha` or `k`.
pub fn gemm(allocator: Allocator, io: Io, comptime T: type, transa: Trans, transb: Trans, m: isize, n: isize, k: isize, alpha: T, a: [*]const T, lda: isize, b: [*]const T, ldb: isize, beta: T, c: [*]T, ldc: isize) Error!void {
    comptime util.requireFloat(T, "gemm");
    return gemmImpl(allocator, io, T, false, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

/// `?GEMM3M` (`interface/gemm.c` with `GEMM3M`): `gemm` by three real
/// products where `gemm` takes four (`driver3m.zig`), taking its
/// arguments as `gemm` does. One part of each element comes from a
/// difference of products, accurate only against the operands' norms
/// (`comparisons.md`, 7). Not routed to Accelerate. Complex only.
pub fn gemm3m(allocator: Allocator, io: Io, comptime T: type, transa: Trans, transb: Trans, m: isize, n: isize, k: isize, alpha: T, a: [*]const T, lda: isize, b: [*]const T, ldb: isize, beta: T, c: [*]T, ldc: isize) Error!void {
    comptime util.requireComplex(T, "gemm3m");
    return gemmImpl(allocator, io, T, true, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

/// `gemm`'s argument checks, in `xerbla`'s order.
pub fn gemmArgs(transa: Trans, transb: Trans, m: isize, n: isize, k: isize, lda: isize, ldb: isize, ldc: isize) ArgError!void {
    if (transa.conjX() or transb.conjX()) return error.InvalidTrans;
    const ta = transa.transposed();
    const tb = transb.transposed();
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    if (k < 0) return error.InvalidK;
    if (lda < @max(1, if (ta) k else m)) return error.InvalidLda;
    if (ldb < @max(1, if (tb) n else k)) return error.InvalidLdb;
    if (ldc < @max(1, m)) return error.InvalidLdc;
}

/// Refuses `gemm` operands that could not exist in memory (`util.fitsMatrix`),
/// `A` and `B` of `size_in` bytes an element, `C` of `size_out`.
pub fn gemmFits(transa: Trans, transb: Trans, m: isize, n: isize, k: isize, lda: isize, ldb: isize, ldc: isize, size_in: usize, size_out: usize) ArgError!void {
    const ta = transa.transposed();
    const tb = transb.transposed();
    if (!util.fitsMatrix(if (ta) k else m, if (ta) m else k, lda, size_in)) return error.InvalidLda;
    if (!util.fitsMatrix(if (tb) n else k, if (tb) k else n, ldb, size_in)) return error.InvalidLdb;
    if (!util.fitsMatrix(m, n, ldc, size_out)) return error.InvalidLdc;
}

fn gemmImpl(allocator: Allocator, io: Io, comptime T: type, comptime three: bool, transa: Trans, transb: Trans, m: isize, n: isize, k: isize, alpha: T, a: [*]const T, lda: isize, b: [*]const T, ldb: isize, beta: T, c: [*]T, ldc: isize) Error!void {
    try gemmArgs(transa, transb, m, n, k, lda, ldb, ldc);
    const ta = transa.transposed();
    const tb = transb.transposed();
    if (m == 0 or n == 0) return;
    // upstream, CBLAS row-major `sgemm` may take `sgemm_direct` or
    // `sgemm_direct_alpha_beta` here; rstd has no row-major layer
    try gemmFits(transa, transb, m, n, k, lda, ldb, ldc, @sizeOf(T), @sizeOf(T));
    if (!three and try accel.gemm(io, T, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)) return;
    const conj: driver.Conj = .of(conjugated(T, transa), conjugated(T, transb));
    const w = threaded.gemmWorkers(T, @intCast(m), @intCast(n), @intCast(k));
    const f = if (three) driver3m.gemm3m else driver.gemm;
    try f(allocator, io, T, table(T), ta, tb, conj, @intCast(m), @intCast(n), @intCast(k), alpha, a, @intCast(lda), b, @intCast(ldb), beta, c, @intCast(ldc), w);
}

/// `?GEMMT` (`interface/gemmt.c`): `C := alpha*op(A)*op(B) + beta*C` on the
/// `uplo` triangle of the `m x m` `C`, `op(A)` being `m x k` and `op(B)`
/// `k x m`, `op` as `gemm`'s. The other triangle is never read or written;
/// the triangle is scaled by `beta` first (a zero `beta` never reads it), and
/// `A` and `B` are not read for a zero `alpha` or `k`. `B` is only read,
/// where upstream conjugates it in place for `r` and `c`.
pub fn gemmt(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, transa: Trans, transb: Trans, m: isize, k: isize, alpha: T, a: [*]const T, lda: isize, b: [*]const T, ldb: isize, beta: T, c: [*]T, ldc: isize) Error!void {
    comptime util.requireFloat(T, "gemmt");
    try gemmtArgs(uplo, transa, transb, m, k, lda, ldb, ldc);
    const ta = transa.transposed();
    const tb = transb.transposed();
    if (m == 0) return;
    try gemmtFits(transa, transb, m, k, lda, ldb, ldc, @sizeOf(T), @sizeOf(T));
    const conj: driver.Conj = .of(conjugated(T, transa), conjugated(T, transb));
    const w = threaded.gemmtWorkers(T, @intCast(m), @intCast(k));
    try driver.gemmt(allocator, io, T, table(T), uplo == .u, ta, tb, conj, @intCast(m), @intCast(k), alpha, a, @intCast(lda), b, @intCast(ldb), beta, c, @intCast(ldc), w);
}

/// `gemmt`'s argument checks, in `xerbla`'s order.
pub fn gemmtArgs(uplo: Uplo, transa: Trans, transb: Trans, m: isize, k: isize, lda: isize, ldb: isize, ldc: isize) ArgError!void {
    if (uplo.isRev()) return error.InvalidUplo;
    if (transa.conjX() or transb.conjX()) return error.InvalidTrans;
    const ta = transa.transposed();
    const tb = transb.transposed();
    if (m < 0) return error.InvalidM;
    if (k < 0) return error.InvalidK;
    if (lda < @max(1, if (ta) k else m)) return error.InvalidLda;
    if (ldb < @max(1, if (tb) m else k)) return error.InvalidLdb;
    if (ldc < @max(1, m)) return error.InvalidLdc;
}

/// `gemmFits` for `gemmt`'s operands.
pub fn gemmtFits(transa: Trans, transb: Trans, m: isize, k: isize, lda: isize, ldb: isize, ldc: isize, size_in: usize, size_out: usize) ArgError!void {
    const ta = transa.transposed();
    const tb = transb.transposed();
    if (!util.fitsMatrix(if (ta) k else m, if (ta) m else k, lda, size_in)) return error.InvalidLda;
    if (!util.fitsMatrix(if (tb) m else k, if (tb) k else m, ldb, size_in)) return error.InvalidLdb;
    if (!util.fitsMatrix(m, m, ldc, size_out)) return error.InvalidLdc;
}

/// `?SYMM` (`interface/symm.c`): `C := alpha*S*B + beta*C` (`side` `l`) or
/// `alpha*B*S + beta*C`, `S` the symmetric matrix whose `uplo` triangle `a`
/// stores, `B` and `C` `m x n`. `C` is scaled by `beta` first (a zero `beta`
/// never reads it), and nothing is read for a zero `alpha`.
pub fn symm(allocator: Allocator, io: Io, comptime T: type, side: Side, uplo: Uplo, m: isize, n: isize, alpha: T, a: [*]const T, lda: isize, b: [*]const T, ldb: isize, beta: T, c: [*]T, ldc: isize) Error!void {
    comptime util.requireFloat(T, "symm");
    return symmImpl(allocator, io, T, false, false, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
}

/// `?SYMM3M` (`interface/symm.c` with `GEMM3M`, `driver/level3/symm3m_k.c`):
/// `symm` by three real products, as `gemm3m`. Complex only.
pub fn symm3m(allocator: Allocator, io: Io, comptime T: type, side: Side, uplo: Uplo, m: isize, n: isize, alpha: T, a: [*]const T, lda: isize, b: [*]const T, ldb: isize, beta: T, c: [*]T, ldc: isize) Error!void {
    comptime util.requireComplex(T, "symm3m");
    return symmImpl(allocator, io, T, false, true, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
}

/// `?HEMM` (`interface/symm.c` as `HEMM`): the same for the Hermitian `S`,
/// whose diagonal's imaginary parts are never read. Complex only.
pub fn hemm(allocator: Allocator, io: Io, comptime T: type, side: Side, uplo: Uplo, m: isize, n: isize, alpha: T, a: [*]const T, lda: isize, b: [*]const T, ldb: isize, beta: T, c: [*]T, ldc: isize) Error!void {
    comptime util.requireComplex(T, "hemm");
    return symmImpl(allocator, io, T, true, false, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
}

/// `?HEMM3M` (`driver/level3/hemm3m_k.c`): `hemm` by three real products,
/// as `gemm3m`. Complex only.
pub fn hemm3m(allocator: Allocator, io: Io, comptime T: type, side: Side, uplo: Uplo, m: isize, n: isize, alpha: T, a: [*]const T, lda: isize, b: [*]const T, ldb: isize, beta: T, c: [*]T, ldc: isize) Error!void {
    comptime util.requireComplex(T, "hemm3m");
    return symmImpl(allocator, io, T, true, true, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
}

fn symmImpl(allocator: Allocator, io: Io, comptime T: type, comptime herm: bool, comptime three: bool, side: Side, uplo: Uplo, m: isize, n: isize, alpha: T, a: [*]const T, lda: isize, b: [*]const T, ldb: isize, beta: T, c: [*]T, ldc: isize) Error!void {
    const ka = if (side == .l) m else n;
    if (uplo.isRev()) return error.InvalidUplo;
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, ka)) return error.InvalidLda;
    if (ldb < @max(1, m)) return error.InvalidLdb;
    if (ldc < @max(1, m)) return error.InvalidLdc;
    if (m == 0 or n == 0) return;
    // upstream, CBLAS row-major `ssymm` may take
    // `ssymm_direct_alpha_beta*` here; rstd has no row-major layer
    if (!util.fitsMatrix(ka, ka, lda, @sizeOf(T))) return error.InvalidLda;
    if (!util.fitsMatrix(m, n, ldb, @sizeOf(T))) return error.InvalidLdb;
    if (!util.fitsMatrix(m, n, ldc, @sizeOf(T))) return error.InvalidLdc;
    if (!three and try accel.symm(io, T, herm, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)) return;
    const w = threaded.symmWorkers(T, @intCast(m), @intCast(n), @intCast(ka));
    const f = if (three) driver3m.symm3m else driver.symm;
    try f(allocator, io, T, table(T), herm, side == .l, uplo == .u, @intCast(m), @intCast(n), alpha, a, @intCast(lda), b, @intCast(ldb), beta, c, @intCast(ldc), w);
}

/// `?SYRK` (`interface/syrk.c`): `C := alpha*op(A)*op(A)^T + beta*C` on the
/// `uplo` triangle of the `n x n` `C`, `op(A)` being `n x k`: `A` for
/// `trans` `n`, `A^T` for `t` (and real `c`). The other triangle is never
/// read or written; `C` is scaled by `beta` first (a zero `beta` never reads
/// it), and `A` is not read for a zero `alpha` or `k`.
pub fn syrk(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, n: isize, k: isize, alpha: T, a: [*]const T, lda: isize, beta: T, c: [*]T, ldc: isize) Error!void {
    comptime util.requireFloat(T, "syrk");
    return rankImpl(allocator, io, T, false, uplo, trans, n, k, alpha, a, lda, null, 0, beta, c, ldc);
}

/// `?HERK` (`interface/syrk.c` as `HERK`): `C := alpha*op(A)*op(A)^H +
/// beta*C`, `op(A)` being `A` (`trans` `n`) or `A^H` (`c`), with real
/// `alpha` and `beta`. The diagonal comes out real: its imaginary parts are
/// read only where `beta` is 1 and nothing is added. Complex only.
pub fn herk(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, n: isize, k: isize, alpha: Float(T), a: [*]const T, lda: isize, beta: Float(T), c: [*]T, ldc: isize) Error!void {
    comptime util.requireComplex(T, "herk");
    return rankImpl(allocator, io, T, true, uplo, trans, n, k, .init(alpha, 0), a, lda, null, 0, .init(beta, 0), c, ldc);
}

/// `?SYR2K` (`interface/syr2k.c`): `C := alpha*op(A)*op(B)^T +
/// alpha*op(B)*op(A)^T + beta*C` on the `uplo` triangle, `op` as `syrk`'s.
pub fn syr2k(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, n: isize, k: isize, alpha: T, a: [*]const T, lda: isize, b: [*]const T, ldb: isize, beta: T, c: [*]T, ldc: isize) Error!void {
    comptime util.requireFloat(T, "syr2k");
    return rankImpl(allocator, io, T, false, uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

/// `?HER2K` (`interface/syr2k.c` as `HER2K`): `C := alpha*op(A)*op(B)^H +
/// conj(alpha)*op(B)*op(A)^H + beta*C`, `op` as `herk`'s, real `beta`.
/// Complex only.
pub fn her2k(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, n: isize, k: isize, alpha: T, a: [*]const T, lda: isize, b: [*]const T, ldb: isize, beta: Float(T), c: [*]T, ldc: isize) Error!void {
    comptime util.requireComplex(T, "her2k");
    return rankImpl(allocator, io, T, true, uplo, trans, n, k, alpha, a, lda, b, ldb, .init(beta, 0), c, ldc);
}

fn rankImpl(allocator: Allocator, io: Io, comptime T: type, comptime herm: bool, uplo: Uplo, trans: Trans, n: isize, k: isize, alpha: T, a: [*]const T, lda: isize, b: ?[*]const T, ldb: isize, beta: T, c: [*]T, ldc: isize) Error!void {
    if (uplo.isRev()) return error.InvalidUplo;
    // real `c` is `t`; complex `syrk` takes `t`, `herk` `c`, and none `r`
    const t = switch (trans) {
        .n => false,
        .t => if (herm) return error.InvalidTrans else true,
        .c => if (isComplex(T) and !herm) return error.InvalidTrans else true,
        .r, .o, .u, .s, .d => return error.InvalidTrans,
    };
    if (n < 0) return error.InvalidN;
    if (k < 0) return error.InvalidK;
    const rows = if (t) k else n;
    const cols = if (t) n else k;
    if (lda < @max(1, rows)) return error.InvalidLda;
    if (b != null and ldb < @max(1, rows)) return error.InvalidLdb;
    if (ldc < @max(1, n)) return error.InvalidLdc;
    if (n == 0) return;
    // upstream, CBLAS row-major `ssyrk` and `ssyr2k` may take their
    // `*_direct_alpha_beta*` kernels here; rstd has no row-major layer
    if (!util.fitsMatrix(rows, cols, lda, @sizeOf(T))) return error.InvalidLda;
    if (b != null and !util.fitsMatrix(rows, cols, ldb, @sizeOf(T))) return error.InvalidLdb;
    if (!util.fitsMatrix(n, n, ldc, @sizeOf(T))) return error.InvalidLdc;
    if (accel.syrk(T, herm, uplo, t, n, k, alpha, a, lda, b, ldb, beta, c, ldc)) return;
    const w = threaded.syrkWorkers(T, b != null, @intCast(n), @intCast(k));
    try driver.syrk(allocator, io, T, table(T), herm, uplo == .u, t, @intCast(n), @intCast(k), alpha, a, @intCast(lda), b, @intCast(ldb), beta, c, @intCast(ldc), w);
}

/// `?TRMM` (`interface/trsm.c` as `TRMM`): `B := alpha*op(A)*B` (`side`
/// `l`) or `alpha*B*op(A)`, `A` triangular (`uplo`; a `unit` diagonal is
/// never read, nor the other triangle), `op` as `transa` says, `B` `m x n`.
/// `B` is scaled by `alpha` first, and neither `A` nor `B` is read for a zero
/// one.
pub fn trmm(allocator: Allocator, io: Io, comptime T: type, side: Side, uplo: Uplo, transa: Trans, diag: Diag, m: isize, n: isize, alpha: T, a: [*]const T, lda: isize, b: [*]T, ldb: isize) Error!void {
    comptime util.requireFloat(T, "trmm");
    return triImpl(allocator, io, T, false, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
}

/// `?TRSM` (`interface/trsm.c`): `B := alpha*inv(op(A))*B` (`side` `l`) or
/// `alpha*B*inv(op(A))`, taking its arguments as `trmm` does. A zero on a
/// diagonal that is not unit divides by zero, as the reference does.
pub fn trsm(allocator: Allocator, io: Io, comptime T: type, side: Side, uplo: Uplo, transa: Trans, diag: Diag, m: isize, n: isize, alpha: T, a: [*]const T, lda: isize, b: [*]T, ldb: isize) Error!void {
    comptime util.requireFloat(T, "trsm");
    return triImpl(allocator, io, T, true, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
}

fn triImpl(allocator: Allocator, io: Io, comptime T: type, comptime solve: bool, side: Side, uplo: Uplo, transa: Trans, diag: Diag, m: isize, n: isize, alpha: T, a: [*]const T, lda: isize, b: [*]T, ldb: isize) Error!void {
    const ka = if (side == .l) m else n;
    if (uplo.isRev()) return error.InvalidUplo;
    if (transa.conjX()) return error.InvalidTrans;
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, ka)) return error.InvalidLda;
    if (ldb < @max(1, m)) return error.InvalidLdb;
    if (m == 0 or n == 0) return;
    // upstream, CBLAS row-major `strmm` may take `strmm_direct_*`
    // here; rstd has no row-major layer
    if (!util.fitsMatrix(ka, ka, lda, @sizeOf(T))) return error.InvalidLda;
    if (!util.fitsMatrix(m, n, ldb, @sizeOf(T))) return error.InvalidLdb;
    if (try accel.trmm(io, T, solve, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)) return;
    const f = if (solve) driver.trsm else driver.trmm;
    const w = threaded.triWorkers(T, @intCast(m), @intCast(n), @intCast(ka));
    try f(allocator, io, T, table(T), side == .l, uplo == .u, transa.transposed(), conjugated(T, transa), diag == .u, @intCast(m), @intCast(n), alpha, a, @intCast(lda), b, @intCast(ldb), w);
}
