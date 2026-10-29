//! The level-3 interface: the checks and quick returns of OpenBLAS's
//! `interface/gemm.c`, `gemmt.c`, `symm.c`, `syrk.c`, `syr2k.c` and `trsm.c`
//! before the drivers (`driver.zig`). Each takes the element type `T` explicitly, any
//! float type or its `Complex` (the vtable's kernels for those it has slots
//! for, the portable ones built for `T` otherwise), and sizes and leading
//! dimensions as `isize`; matrices are column-major.
//! A matrix that would not fit the address space is its leading
//! dimension's error. Bad arguments return an `ArgError` in the
//! order the reference numbers them for `xerbla`, so the first bad one wins.
//! Leading dimensions are checked against `max(1, rows)`, as the reference
//! BLAS does (OpenBLAS lets `0` through for empty matrices). Packing space
//! comes from `allocator`; small problems use the stack. Every routine but
//! `gemmt` (serial, as upstream) threads on `io` above its threshold
//! (`threading.md`).

const std = @import("std");
const blas = @import("../__root__.zig");
const driver = @import("driver.zig");
const threaded = @import("threaded.zig");
const threading = @import("../threading.zig");
const types = @import("../types.zig");
const portable = @import("../portable.zig");
const util = @import("../util.zig");
const Vtable = @import("../Vtable.zig");
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

fn floatType(comptime T: type, comptime name: []const u8) void {
    if (@typeInfo(Float(T)) != .float) @compileError(name ++ ": " ++ @typeName(T) ++ " is not a float type or its Complex");
}

fn complexType(comptime T: type, comptime name: []const u8) void {
    if (!isComplex(T)) @compileError(name ++ ": " ++ @typeName(T) ++ " is not complex");
}

/// The vtable for `T`: a type without slots runs the portable kernels and
/// never needs it.
fn table(comptime T: type) *const Vtable {
    return if (comptime portable.letter(T) != null) blas.vtable() else &no_table;
}

/// `?GEMM`: `C := alpha*op(A)*op(B) + beta*C`, `op(A)` being `m x k`,
/// `op(B)` `k x n`. `C` is scaled by `beta` first (a zero `beta` never reads
/// it), and `A` and `B` are not read for a zero `alpha` or `k`.
pub fn gemm(allocator: Allocator, io: Io, comptime T: type, transa: Trans, transb: Trans, m: isize, n: isize, k: isize, alpha: T, a: [*]const T, lda: isize, b: [*]const T, ldb: isize, beta: T, c: [*]T, ldc: isize) Error!void {
    comptime floatType(T, "gemm");
    if (transa.conjX() or transb.conjX()) return error.InvalidTrans;
    const ta = transa.transposed();
    const tb = transb.transposed();
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    if (k < 0) return error.InvalidK;
    if (lda < @max(1, if (ta) k else m)) return error.InvalidLda;
    if (ldb < @max(1, if (tb) n else k)) return error.InvalidLdb;
    if (ldc < @max(1, m)) return error.InvalidLdc;
    if (m == 0 or n == 0) return;
    // upstream, CBLAS row-major `sgemm` may take `sgemm_direct` or
    // `sgemm_direct_alpha_beta` here; rstd has no row-major layer
    if (!util.fitsMatrix(if (ta) k else m, if (ta) m else k, lda, @sizeOf(T))) return error.InvalidLda;
    if (!util.fitsMatrix(if (tb) n else k, if (tb) k else n, ldb, @sizeOf(T))) return error.InvalidLdb;
    if (!util.fitsMatrix(m, n, ldc, @sizeOf(T))) return error.InvalidLdc;
    if (accel.gemm(T, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)) return;
    const conj: driver.Conj = .of(conjugated(T, transa), conjugated(T, transb));
    const w = threaded.workers(.gemm, T, threading.product(&.{ @intCast(m), @intCast(n), @intCast(k) }));
    try driver.gemm(allocator, io, T, table(T), ta, tb, conj, @intCast(m), @intCast(n), @intCast(k), alpha, a, @intCast(lda), b, @intCast(ldb), beta, c, @intCast(ldc), w);
}

/// `?GEMMT` (`interface/gemmt.c`): `C := alpha*op(A)*op(B) + beta*C` on the
/// `uplo` triangle of the `m x m` `C`, `op(A)` being `m x k` and `op(B)`
/// `k x m`, `op` as `gemm`'s. The other triangle is never read or written;
/// the triangle is scaled by `beta` first (a zero `beta` never reads it), and
/// `A` and `B` are not read for a zero `alpha` or `k`. `B` is only read,
/// where upstream conjugates it in place for `r` and `c`.
pub fn gemmt(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, transa: Trans, transb: Trans, m: isize, k: isize, alpha: T, a: [*]const T, lda: isize, b: [*]const T, ldb: isize, beta: T, c: [*]T, ldc: isize) Error!void {
    comptime floatType(T, "gemmt");
    if (uplo.isRev()) return error.InvalidUplo;
    if (transa.conjX() or transb.conjX()) return error.InvalidTrans;
    const ta = transa.transposed();
    const tb = transb.transposed();
    if (m < 0) return error.InvalidM;
    if (k < 0) return error.InvalidK;
    if (lda < @max(1, if (ta) k else m)) return error.InvalidLda;
    if (ldb < @max(1, if (tb) m else k)) return error.InvalidLdb;
    if (ldc < @max(1, m)) return error.InvalidLdc;
    if (m == 0) return;
    if (!util.fitsMatrix(if (ta) k else m, if (ta) m else k, lda, @sizeOf(T))) return error.InvalidLda;
    if (!util.fitsMatrix(if (tb) m else k, if (tb) k else m, ldb, @sizeOf(T))) return error.InvalidLdb;
    if (!util.fitsMatrix(m, m, ldc, @sizeOf(T))) return error.InvalidLdc;
    const conj: driver.Conj = .of(conjugated(T, transa), conjugated(T, transb));
    try driver.gemmt(allocator, io, T, table(T), uplo == .u, ta, tb, conj, @intCast(m), @intCast(k), alpha, a, @intCast(lda), b, @intCast(ldb), beta, c, @intCast(ldc));
}

/// `?SYMM` (`interface/symm.c`): `C := alpha*S*B + beta*C` (`side` `l`) or
/// `alpha*B*S + beta*C`, `S` the symmetric matrix whose `uplo` triangle `a`
/// stores, `B` and `C` `m x n`. `C` is scaled by `beta` first (a zero `beta`
/// never reads it), and nothing is read for a zero `alpha`.
pub fn symm(allocator: Allocator, io: Io, comptime T: type, side: Side, uplo: Uplo, m: isize, n: isize, alpha: T, a: [*]const T, lda: isize, b: [*]const T, ldb: isize, beta: T, c: [*]T, ldc: isize) Error!void {
    comptime floatType(T, "symm");
    return symmImpl(allocator, io, T, false, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
}

/// `?HEMM` (`interface/symm.c` as `HEMM`): the same for the Hermitian `S`,
/// whose diagonal's imaginary parts are never read. Complex only.
pub fn hemm(allocator: Allocator, io: Io, comptime T: type, side: Side, uplo: Uplo, m: isize, n: isize, alpha: T, a: [*]const T, lda: isize, b: [*]const T, ldb: isize, beta: T, c: [*]T, ldc: isize) Error!void {
    comptime complexType(T, "hemm");
    return symmImpl(allocator, io, T, true, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
}

fn symmImpl(allocator: Allocator, io: Io, comptime T: type, comptime herm: bool, side: Side, uplo: Uplo, m: isize, n: isize, alpha: T, a: [*]const T, lda: isize, b: [*]const T, ldb: isize, beta: T, c: [*]T, ldc: isize) Error!void {
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
    if (accel.symm(T, herm, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)) return;
    // `symm.c` counts `2 * m * m * n` on either side
    const w = threaded.workers(.symm, T, threading.product(&.{ 2, @intCast(m), @intCast(m), @intCast(n) }));
    try driver.symm(allocator, io, T, table(T), herm, side == .l, uplo == .u, @intCast(m), @intCast(n), alpha, a, @intCast(lda), b, @intCast(ldb), beta, c, @intCast(ldc), w);
}

/// `?SYRK` (`interface/syrk.c`): `C := alpha*op(A)*op(A)^T + beta*C` on the
/// `uplo` triangle of the `n x n` `C`, `op(A)` being `n x k`: `A` for
/// `trans` `n`, `A^T` for `t` (and real `c`). The other triangle is never
/// read or written; `C` is scaled by `beta` first (a zero `beta` never reads
/// it), and `A` is not read for a zero `alpha` or `k`.
pub fn syrk(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, n: isize, k: isize, alpha: T, a: [*]const T, lda: isize, beta: T, c: [*]T, ldc: isize) Error!void {
    comptime floatType(T, "syrk");
    return rankImpl(allocator, io, T, false, uplo, trans, n, k, alpha, a, lda, null, 0, beta, c, ldc);
}

/// `?HERK` (`interface/syrk.c` as `HERK`): `C := alpha*op(A)*op(A)^H +
/// beta*C`, `op(A)` being `A` (`trans` `n`) or `A^H` (`c`), with real
/// `alpha` and `beta`. The diagonal comes out real: its imaginary parts are
/// read only where `beta` is 1 and nothing is added. Complex only.
pub fn herk(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, n: isize, k: isize, alpha: Float(T), a: [*]const T, lda: isize, beta: Float(T), c: [*]T, ldc: isize) Error!void {
    comptime complexType(T, "herk");
    return rankImpl(allocator, io, T, true, uplo, trans, n, k, .init(alpha, 0), a, lda, null, 0, .init(beta, 0), c, ldc);
}

/// `?SYR2K` (`interface/syr2k.c`): `C := alpha*op(A)*op(B)^T +
/// alpha*op(B)*op(A)^T + beta*C` on the `uplo` triangle, `op` as `syrk`'s.
pub fn syr2k(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, n: isize, k: isize, alpha: T, a: [*]const T, lda: isize, b: [*]const T, ldb: isize, beta: T, c: [*]T, ldc: isize) Error!void {
    comptime floatType(T, "syr2k");
    return rankImpl(allocator, io, T, false, uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

/// `?HER2K` (`interface/syr2k.c` as `HER2K`): `C := alpha*op(A)*op(B)^H +
/// conj(alpha)*op(B)*op(A)^H + beta*C`, `op` as `herk`'s, real `beta`.
/// Complex only.
pub fn her2k(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, n: isize, k: isize, alpha: T, a: [*]const T, lda: isize, b: [*]const T, ldb: isize, beta: Float(T), c: [*]T, ldc: isize) Error!void {
    comptime complexType(T, "her2k");
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
    const w = if (b == null)
        threaded.workers(.syrk, T, threading.product(&.{ @intCast(n + 1), @intCast(n), @intCast(k) }))
    else
        threaded.workers(.syr2k, T, threading.product(&.{ @intCast(n), @intCast(k) }));
    try driver.syrk(allocator, io, T, table(T), herm, uplo == .u, t, @intCast(n), @intCast(k), alpha, a, @intCast(lda), b, @intCast(ldb), beta, c, @intCast(ldc), w);
}

/// `?TRMM` (`interface/trsm.c` as `TRMM`): `B := alpha*op(A)*B` (`side`
/// `l`) or `alpha*B*op(A)`, `A` triangular (`uplo`; a `unit` diagonal is
/// never read, nor the other triangle), `op` as `transa` says, `B` `m x n`.
/// `B` is scaled by `alpha` first, and neither `A` nor `B` is read for a zero
/// one.
pub fn trmm(allocator: Allocator, io: Io, comptime T: type, side: Side, uplo: Uplo, transa: Trans, diag: Diag, m: isize, n: isize, alpha: T, a: [*]const T, lda: isize, b: [*]T, ldb: isize) Error!void {
    comptime floatType(T, "trmm");
    return triImpl(allocator, io, T, false, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb);
}

/// `?TRSM` (`interface/trsm.c`): `B := alpha*inv(op(A))*B` (`side` `l`) or
/// `alpha*B*inv(op(A))`, taking its arguments as `trmm` does. A zero on a
/// diagonal that is not unit divides by zero, as the reference does.
pub fn trsm(allocator: Allocator, io: Io, comptime T: type, side: Side, uplo: Uplo, transa: Trans, diag: Diag, m: isize, n: isize, alpha: T, a: [*]const T, lda: isize, b: [*]T, ldb: isize) Error!void {
    comptime floatType(T, "trsm");
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
    if (accel.trmm(T, solve, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)) return;
    const f = if (solve) driver.trsm else driver.trmm;
    const w = threaded.workers(.trsm, T, threading.product(&.{ @intCast(m), @intCast(n) }));
    try f(allocator, io, T, table(T), side == .l, uplo == .u, transa.transposed(), conjugated(T, transa), diag == .u, @intCast(m), @intCast(n), alpha, a, @intCast(lda), b, @intCast(ldb), w);
}

const no_table: Vtable = .{};

// ---- tests -----------------------------------------------------------------

const testing = std.testing;
const check = @import("../_test/blas3/check.zig");

test "interface: gemm reports the first bad argument and does nothing for empty problems" {
    const fa = testing.failing_allocator;
    var a: [16]f64 = @splat(1);
    var c: [16]f64 = @splat(std.math.nan(f64));
    const nan = std.math.nan(f64);
    try testing.expectError(error.InvalidTrans, gemm(fa, check.io, f64, .n, .o, -1, -1, -1, 1, &a, 0, &a, 0, 0, &c, 0));
    try testing.expectError(error.InvalidM, gemm(fa, check.io, f64, .n, .n, -1, -1, -1, 1, &a, 0, &a, 0, 0, &c, 0));
    try testing.expectError(error.InvalidN, gemm(fa, check.io, f64, .n, .n, 0, -1, -1, 1, &a, 0, &a, 0, 0, &c, 0));
    try testing.expectError(error.InvalidK, gemm(fa, check.io, f64, .n, .n, 0, 0, -1, 1, &a, 0, &a, 0, 0, &c, 0));
    try testing.expectError(error.InvalidLda, gemm(fa, check.io, f64, .n, .n, 3, 2, 4, 1, &a, 2, &a, 1, 0, &c, 1));
    try testing.expectError(error.InvalidLda, gemm(fa, check.io, f64, .t, .n, 3, 2, 4, 1, &a, 3, &a, 1, 0, &c, 1));
    try testing.expectError(error.InvalidLdb, gemm(fa, check.io, f64, .n, .n, 3, 2, 4, 1, &a, 3, &a, 3, 0, &c, 1));
    try testing.expectError(error.InvalidLdb, gemm(fa, check.io, f64, .n, .c, 3, 2, 4, 1, &a, 3, &a, 1, 0, &c, 1));
    try testing.expectError(error.InvalidLdc, gemm(fa, check.io, f64, .n, .n, 3, 2, 4, 1, &a, 3, &a, 4, 0, &c, 2));
    try testing.expectError(error.InvalidLdc, gemm(fa, check.io, f64, .n, .n, 0, 2, 0, 1, &a, 1, &a, 1, 0, &c, 0));
    // empty: nothing read or written, nothing allocated
    try gemm(fa, check.io, f64, .n, .n, 0, 3, 2, nan, &a, 1, &a, 2, 0, &c, 1);
    try gemm(fa, check.io, f64, .n, .n, 3, 0, 2, nan, &a, 3, &a, 2, 0, &c, 3);
    try testing.expect(std.math.isNan(c[0]));
    // k == 0 only scales C, and allocates nothing
    try gemm(fa, check.io, f64, .n, .n, 2, 2, 0, nan, &a, 2, &a, 1, 0, &c, 2);
    try testing.expectEqualSlices(f64, &.{ 0, 0, 0, 0 }, c[0..4]);
}

const all = .{ f32, f64, types.Complex32, types.Complex64 };
const unslotted = .{ f16, f128, std.math.Complex(f16), std.math.Complex(f128) } ++ check.f80_unslotted;

test "interface: gemm matches the wide reference through the global vtable" {
    inline for (all) |T| try check.interfaceCases(T, .gemm);
}

test "interface: gemm over types the vtable has no slots for runs the portable kernels" {
    inline for (unslotted) |T| try check.interfaceCases(T, .gemm);
}

test "interface: a matrix past the address space is its leading dimension's error" {
    const fa = testing.failing_allocator;
    var a: [4]f64 = @splat(1);
    const big = std.math.maxInt(isize) / 2;
    try testing.expectError(error.InvalidLda, gemm(fa, check.io, f64, .n, .n, 2, 2, 3, 1, &a, big, &a, 3, 0, &a, 2));
    try testing.expectError(error.InvalidLdb, gemm(fa, check.io, f64, .n, .n, 2, 2, 3, 1, &a, 2, &a, big, 0, &a, 2));
    try testing.expectError(error.InvalidLdc, gemm(fa, check.io, f64, .n, .n, 2, 2, 3, 1, &a, 2, &a, 3, 0, &a, big));
}

test "interface: gemmt matches the wide reference" {
    inline for (all ++ unslotted) |T| try check.interfaceCases(T, .gemmt);
}

test "interface: gemmt reports the first bad argument, and returns for empty problems" {
    const fa = testing.failing_allocator;
    const nan = std.math.nan(f64);
    var a: [16]f64 = @splat(nan);
    var c: [16]f64 = @splat(nan);
    try testing.expectError(error.InvalidUplo, gemmt(fa, check.io, f64, .v, .o, .n, -1, -1, 1, &a, 0, &a, 0, 0, &c, 0));
    try testing.expectError(error.InvalidTrans, gemmt(fa, check.io, f64, .u, .o, .n, -1, -1, 1, &a, 0, &a, 0, 0, &c, 0));
    try testing.expectError(error.InvalidTrans, gemmt(fa, check.io, f64, .u, .n, .d, -1, -1, 1, &a, 0, &a, 0, 0, &c, 0));
    try testing.expectError(error.InvalidM, gemmt(fa, check.io, f64, .u, .n, .n, -1, -1, 1, &a, 0, &a, 0, 0, &c, 0));
    try testing.expectError(error.InvalidK, gemmt(fa, check.io, f64, .l, .n, .n, 0, -1, 1, &a, 0, &a, 0, 0, &c, 0));
    try testing.expectError(error.InvalidLda, gemmt(fa, check.io, f64, .l, .n, .n, 3, 2, 1, &a, 2, &a, 2, 0, &c, 3));
    try testing.expectError(error.InvalidLda, gemmt(fa, check.io, f64, .l, .t, .n, 3, 4, 1, &a, 3, &a, 4, 0, &c, 3));
    try testing.expectError(error.InvalidLdb, gemmt(fa, check.io, f64, .u, .n, .n, 3, 4, 1, &a, 3, &a, 3, 0, &c, 3));
    try testing.expectError(error.InvalidLdb, gemmt(fa, check.io, f64, .u, .n, .c, 3, 4, 1, &a, 3, &a, 2, 0, &c, 3));
    try testing.expectError(error.InvalidLdc, gemmt(fa, check.io, f64, .u, .n, .n, 3, 4, 1, &a, 3, &a, 4, 0, &c, 2));
    try testing.expectError(error.InvalidLdc, gemmt(fa, check.io, f64, .u, .n, .n, 0, 4, 1, &a, 1, &a, 4, 0, &c, 0));
    try gemmt(fa, check.io, f64, .u, .n, .n, 0, 3, nan, &a, 1, &a, 3, nan, &c, 1);
    try testing.expect(std.math.isNan(c[0]));
    // `k` 0 scales the triangle alone, and allocates nothing
    try gemmt(fa, check.io, f64, .l, .n, .n, 2, 0, nan, &a, 2, &a, 1, 0, &c, 2);
    try testing.expectEqual(0, c[0]);
    try testing.expectEqual(0, c[1]);
    try testing.expect(std.math.isNan(c[2]));
    try testing.expectEqual(0, c[3]);
}

test "interface: symm and hemm match the wide reference" {
    inline for (all ++ unslotted) |T| {
        try check.interfaceCases(T, .symm);
        if (comptime isComplex(T)) try check.interfaceCases(T, .hemm);
    }
}

test "interface: symm and hemm report the first bad argument, and return for empty problems" {
    const fa = testing.failing_allocator;
    var a: [16]f64 = @splat(std.math.nan(f64));
    var c: [16]f64 = @splat(std.math.nan(f64));
    try testing.expectError(error.InvalidUplo, symm(fa, check.io, f64, .l, .v, -1, -1, 1, &a, 0, &a, 0, 0, &c, 0));
    try testing.expectError(error.InvalidM, symm(fa, check.io, f64, .l, .u, -1, -1, 1, &a, 0, &a, 0, 0, &c, 0));
    try testing.expectError(error.InvalidN, symm(fa, check.io, f64, .l, .u, 0, -1, 1, &a, 0, &a, 0, 0, &c, 0));
    // the order of `S` is `m` on the left, `n` on the right
    try testing.expectError(error.InvalidLda, symm(fa, check.io, f64, .l, .u, 3, 0, 1, &a, 2, &a, 3, 0, &c, 3));
    try testing.expectError(error.InvalidLda, symm(fa, check.io, f64, .r, .l, 1, 3, 1, &a, 2, &a, 1, 0, &c, 1));
    try testing.expectError(error.InvalidLdb, symm(fa, check.io, f64, .r, .u, 2, 1, 1, &a, 1, &a, 1, 0, &c, 2));
    try testing.expectError(error.InvalidLdc, symm(fa, check.io, f64, .l, .u, 2, 1, 1, &a, 2, &a, 2, 0, &c, 1));
    try symm(fa, check.io, f64, .l, .u, 0, 3, std.math.nan(f64), &a, 1, &a, 1, 0, &c, 1);
    try symm(fa, check.io, f64, .r, .l, 3, 0, std.math.nan(f64), &a, 1, &a, 3, 0, &c, 3);
    try testing.expect(std.math.isNan(c[0]));
    var z: [4]types.Complex64 = @splat(.init(0, 0));
    try testing.expectError(error.InvalidLda, hemm(fa, check.io, types.Complex64, .r, .u, 1, 2, .init(1, 0), &z, 1, &z, 1, .init(0, 0), &z, 1));
}

test "interface: syrk, herk, syr2k and her2k match the wide reference" {
    inline for (all ++ unslotted) |T| {
        try check.interfaceCases(T, .syrk);
        try check.interfaceCases(T, .syr2k);
        if (comptime isComplex(T)) {
            try check.interfaceCases(T, .herk);
            try check.interfaceCases(T, .her2k);
        }
    }
}

test "interface: the rank-k updates report the first bad argument, and return for empty problems" {
    const fa = testing.failing_allocator;
    const Z = types.Complex64;
    var a: [16]f64 = @splat(std.math.nan(f64));
    var c: [16]f64 = @splat(std.math.nan(f64));
    var z: [4]Z = @splat(.init(0, 0));
    try testing.expectError(error.InvalidUplo, syrk(fa, check.io, f64, .m, .r, -1, 0, 1, &a, 1, 0, &c, 1));
    try testing.expectError(error.InvalidTrans, syrk(fa, check.io, f64, .u, .r, -1, 0, 1, &a, 1, 0, &c, 1));
    try testing.expectError(error.InvalidTrans, syr2k(fa, check.io, f64, .u, .u, -1, 0, 1, &a, 1, &a, 1, 0, &c, 1));
    try testing.expectError(error.InvalidTrans, syrk(fa, check.io, Z, .u, .c, -1, 0, .init(1, 0), &z, 1, .init(0, 0), &z, 1));
    try testing.expectError(error.InvalidTrans, herk(fa, check.io, Z, .l, .t, -1, 0, 1, &z, 1, 0, &z, 1));
    try testing.expectError(error.InvalidN, syrk(fa, check.io, f64, .u, .c, -1, -1, 1, &a, 1, 0, &c, 1));
    try testing.expectError(error.InvalidK, syrk(fa, check.io, f64, .u, .n, 0, -1, 1, &a, 1, 0, &c, 1));
    try testing.expectError(error.InvalidLda, syrk(fa, check.io, f64, .u, .t, 1, 3, 1, &a, 2, 0, &c, 1));
    try testing.expectError(error.InvalidLdb, syr2k(fa, check.io, f64, .l, .n, 3, 1, 1, &a, 3, &a, 2, 0, &c, 3));
    try testing.expectError(error.InvalidLdc, syr2k(fa, check.io, f64, .l, .n, 3, 1, 1, &a, 3, &a, 3, 0, &c, 2));
    try testing.expectError(error.InvalidLdc, her2k(fa, check.io, Z, .u, .c, 2, 0, .init(1, 0), &z, 1, &z, 1, 0, &z, 1));
    try syrk(fa, check.io, f64, .u, .n, 0, 3, std.math.nan(f64), &a, 1, 0, &c, 1);
    try testing.expect(std.math.isNan(c[0]));
    // nothing to add and `beta` 1: a Hermitian diagonal keeps its imaginary part
    z[0] = .init(2, 5);
    try herk(fa, check.io, Z, .u, .n, 1, 0, 1, &z, 1, 1, &z, 1);
    try testing.expectEqual(Z.init(2, 5), z[0]);
    try herk(fa, check.io, Z, .u, .n, 1, 0, 1, &z, 1, 2, &z, 1);
    try testing.expectEqual(Z.init(4, 0), z[0]);
}

test "interface: trmm and trsm match the wide reference" {
    inline for (all) |T| {
        try check.interfaceCases(T, .trmm);
        try check.interfaceCases(T, .trsm);
    }
}

test "interface: unslotted trmm and trsm match the wide reference" {
    inline for (unslotted) |T| {
        try check.interfaceCases(T, .trmm);
        try check.interfaceCases(T, .trsm);
    }
}

test "interface: trmm and trsm report the first bad argument, and return for empty problems" {
    const fa = testing.failing_allocator;
    var a: [16]f64 = @splat(std.math.nan(f64));
    var b: [16]f64 = @splat(std.math.nan(f64));
    inline for (.{ trmm, trsm }) |f| {
        try testing.expectError(error.InvalidUplo, f(fa, check.io, f64, .l, .v, .s, .n, -1, -1, 1, &a, 0, &b, 0));
        try testing.expectError(error.InvalidTrans, f(fa, check.io, f64, .l, .u, .s, .n, -1, -1, 1, &a, 0, &b, 0));
        try testing.expectError(error.InvalidM, f(fa, check.io, f64, .l, .u, .n, .n, -1, -1, 1, &a, 0, &b, 0));
        try testing.expectError(error.InvalidN, f(fa, check.io, f64, .r, .l, .t, .u, 0, -1, 1, &a, 0, &b, 0));
        // the order of `A` is `m` on the left, `n` on the right
        try testing.expectError(error.InvalidLda, f(fa, check.io, f64, .l, .u, .n, .n, 2, 0, 1, &a, 1, &b, 2));
        try testing.expectError(error.InvalidLda, f(fa, check.io, f64, .r, .u, .c, .n, 0, 2, 1, &a, 1, &b, 1));
        try testing.expectError(error.InvalidLdb, f(fa, check.io, f64, .r, .l, .n, .n, 2, 1, 1, &a, 1, &b, 1));
        try f(fa, check.io, f64, .l, .u, .n, .n, 0, 3, std.math.nan(f64), &a, 1, &b, 1);
        try f(fa, check.io, f64, .r, .u, .n, .n, 3, 0, std.math.nan(f64), &a, 1, &b, 3);
        try testing.expect(std.math.isNan(b[0]));
        // a zero `alpha` zeroes `B` without reading either
        try f(fa, check.io, f64, .l, .l, .t, .u, 2, 2, 0, &a, 2, &b, 3);
        try testing.expectEqualSlices(f64, &.{ 0, 0 }, b[0..2]);
        try testing.expectEqualSlices(f64, &.{ 0, 0 }, b[3..5]);
        try testing.expect(std.math.isNan(b[2]));
        @memset(&b, std.math.nan(f64));
    }
}

test "the routines through Accelerate: special scalars, and special values in each operand" {
    if (comptime !accel.covers(f64)) return error.SkipZigTest;
    _ = blas.vtable();
    const was = accel.routing();
    defer _ = accel.route(was);
    if (!accel.route(true)) return error.SkipZigTest;
    var failed: ?anyerror = null;
    inline for (.{ f32, f64, types.Complex32, types.Complex64 }) |T| {
        inline for (.{ .gemm, .symm, .hemm, .syrk, .herk, .syr2k, .her2k, .trmm, .trsm }) |r| {
            if (comptime isComplex(T) or !(r == .hemm or r == .herk or r == .her2k)) {
                check.accelerateCases(T, r) catch |e| {
                    failed = e;
                };
            }
        }
    }
    if (failed) |e| return e;
}
