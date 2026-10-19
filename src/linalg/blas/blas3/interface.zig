//! The level-3 interface: the checks and quick returns of OpenBLAS's
//! `interface/gemm.c`, `symm.c`, `syrk.c`, `syr2k.c` and `trsm.c` before the
//! drivers (`level3.zig`). Each takes the element type `T` explicitly, any
//! float type or its `Complex` (the vtable's kernels for those it has slots
//! for, the portable ones built for `T` otherwise), and sizes and leading
//! dimensions as `isize`; matrices are column-major.
//! A matrix that would not fit the address space is its leading
//! dimension's error. Bad arguments return an `ArgError` in the
//! order the reference numbers them for `xerbla`, so the first bad one wins.
//! Leading dimensions are checked against `max(1, rows)`, as the reference
//! BLAS does (OpenBLAS lets `0` through for empty matrices). Packing space
//! comes from `allocator`; small problems use the stack. Single-threaded.

const std = @import("std");
const blas = @import("../__root__.zig");
const level3 = @import("level3.zig");
const types = @import("../types.zig");
const portable = @import("../portable.zig");
const util = @import("../util.zig");
const Vtable = @import("../Vtable.zig");

const Allocator = std.mem.Allocator;
const isComplex = types.isComplex;
const Float = types.Float;

/// `TRANSA`, `TRANSB`: `n` is the matrix, `t` its transpose, `c` its
/// conjugate transpose and `r` its conjugate. Real `r` is `n` and `c` is `t`.
pub const Trans = enum { n, t, r, c };

fn transposed(t: Trans) bool {
    return switch (t) {
        .t, .c => true,
        .n, .r => false,
    };
}

fn conjugated(comptime T: type, t: Trans) bool {
    return isComplex(T) and (t == .r or t == .c);
}

/// The reference's `xerbla` argument, by name.
pub const ArgError = error{
    InvalidTrans,
    InvalidM,
    InvalidN,
    InvalidK,
    InvalidLda,
    InvalidLdb,
    InvalidLdc,
};
pub const Error = ArgError || Allocator.Error;

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
pub fn gemm(allocator: Allocator, comptime T: type, transa: Trans, transb: Trans, m: isize, n: isize, k: isize, alpha: T, a: [*]const T, lda: isize, b: [*]const T, ldb: isize, beta: T, c: [*]T, ldc: isize) Error!void {
    comptime floatType(T, "gemm");
    const ta = transposed(transa);
    const tb = transposed(transb);
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    if (k < 0) return error.InvalidK;
    if (lda < @max(1, if (ta) k else m)) return error.InvalidLda;
    if (ldb < @max(1, if (tb) n else k)) return error.InvalidLdb;
    if (ldc < @max(1, m)) return error.InvalidLdc;
    if (m == 0 or n == 0) return;
    if (!util.fitsMatrix(if (ta) k else m, if (ta) m else k, lda, @sizeOf(T))) return error.InvalidLda;
    if (!util.fitsMatrix(if (tb) n else k, if (tb) k else n, ldb, @sizeOf(T))) return error.InvalidLdb;
    if (!util.fitsMatrix(m, n, ldc, @sizeOf(T))) return error.InvalidLdc;
    const conj: level3.Conj = .of(conjugated(T, transa), conjugated(T, transb));
    try level3.gemm(allocator, T, table(T), ta, tb, conj, @intCast(m), @intCast(n), @intCast(k), alpha, a, @intCast(lda), b, @intCast(ldb), beta, c, @intCast(ldc));
}

/// `SIDE`: whether the symmetric or triangular matrix multiplies from the
/// left (`l`) or the right (`r`).
pub const Side = enum { l, r };

/// `UPLO`: the triangle a matrix stores, or a rank-k update computes.
pub const Uplo = enum { u, l };

/// `DIAG`: a unit (`u`) diagonal is taken as ones and never read.
pub const Diag = enum { u, n };

/// `?SYMM` (`interface/symm.c`): `C := alpha*S*B + beta*C` (`side` `l`) or
/// `alpha*B*S + beta*C`, `S` the symmetric matrix whose `uplo` triangle `a`
/// stores, `B` and `C` `m x n`. `C` is scaled by `beta` first (a zero `beta`
/// never reads it), and nothing is read for a zero `alpha`.
pub fn symm(allocator: Allocator, comptime T: type, side: Side, uplo: Uplo, m: isize, n: isize, alpha: T, a: [*]const T, lda: isize, b: [*]const T, ldb: isize, beta: T, c: [*]T, ldc: isize) Error!void {
    comptime floatType(T, "symm");
    return symmImpl(allocator, T, false, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
}

/// `?HEMM` (`interface/symm.c` as `HEMM`): the same for the Hermitian `S`,
/// whose diagonal's imaginary parts are never read. Complex only.
pub fn hemm(allocator: Allocator, comptime T: type, side: Side, uplo: Uplo, m: isize, n: isize, alpha: T, a: [*]const T, lda: isize, b: [*]const T, ldb: isize, beta: T, c: [*]T, ldc: isize) Error!void {
    comptime complexType(T, "hemm");
    return symmImpl(allocator, T, true, side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);
}

fn symmImpl(allocator: Allocator, comptime T: type, comptime herm: bool, side: Side, uplo: Uplo, m: isize, n: isize, alpha: T, a: [*]const T, lda: isize, b: [*]const T, ldb: isize, beta: T, c: [*]T, ldc: isize) Error!void {
    const ka = if (side == .l) m else n;
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, ka)) return error.InvalidLda;
    if (ldb < @max(1, m)) return error.InvalidLdb;
    if (ldc < @max(1, m)) return error.InvalidLdc;
    if (m == 0 or n == 0) return;
    if (!util.fitsMatrix(ka, ka, lda, @sizeOf(T))) return error.InvalidLda;
    if (!util.fitsMatrix(m, n, ldb, @sizeOf(T))) return error.InvalidLdb;
    if (!util.fitsMatrix(m, n, ldc, @sizeOf(T))) return error.InvalidLdc;
    try level3.symm(allocator, T, table(T), herm, side == .l, uplo == .u, @intCast(m), @intCast(n), alpha, a, @intCast(lda), b, @intCast(ldb), beta, c, @intCast(ldc));
}

/// `?SYRK` (`interface/syrk.c`): `C := alpha*op(A)*op(A)^T + beta*C` on the
/// `uplo` triangle of the `n x n` `C`, `op(A)` being `n x k`: `A` for
/// `trans` `n`, `A^T` for `t` (and real `c`). The other triangle is never
/// read or written; `C` is scaled by `beta` first (a zero `beta` never reads
/// it), and `A` is not read for a zero `alpha` or `k`.
pub fn syrk(allocator: Allocator, comptime T: type, uplo: Uplo, trans: Trans, n: isize, k: isize, alpha: T, a: [*]const T, lda: isize, beta: T, c: [*]T, ldc: isize) Error!void {
    comptime floatType(T, "syrk");
    return rankImpl(allocator, T, false, uplo, trans, n, k, alpha, a, lda, null, 0, beta, c, ldc);
}

/// `?HERK` (`interface/syrk.c` as `HERK`): `C := alpha*op(A)*op(A)^H +
/// beta*C`, `op(A)` being `A` (`trans` `n`) or `A^H` (`c`), with real
/// `alpha` and `beta`. The diagonal comes out real: its imaginary parts are
/// read only where `beta` is 1 and nothing is added. Complex only.
pub fn herk(allocator: Allocator, comptime T: type, uplo: Uplo, trans: Trans, n: isize, k: isize, alpha: Float(T), a: [*]const T, lda: isize, beta: Float(T), c: [*]T, ldc: isize) Error!void {
    comptime complexType(T, "herk");
    return rankImpl(allocator, T, true, uplo, trans, n, k, .init(alpha, 0), a, lda, null, 0, .init(beta, 0), c, ldc);
}

/// `?SYR2K` (`interface/syr2k.c`): `C := alpha*op(A)*op(B)^T +
/// alpha*op(B)*op(A)^T + beta*C` on the `uplo` triangle, `op` as `syrk`'s.
pub fn syr2k(allocator: Allocator, comptime T: type, uplo: Uplo, trans: Trans, n: isize, k: isize, alpha: T, a: [*]const T, lda: isize, b: [*]const T, ldb: isize, beta: T, c: [*]T, ldc: isize) Error!void {
    comptime floatType(T, "syr2k");
    return rankImpl(allocator, T, false, uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
}

/// `?HER2K` (`interface/syr2k.c` as `HER2K`): `C := alpha*op(A)*op(B)^H +
/// conj(alpha)*op(B)*op(A)^H + beta*C`, `op` as `herk`'s, real `beta`.
/// Complex only.
pub fn her2k(allocator: Allocator, comptime T: type, uplo: Uplo, trans: Trans, n: isize, k: isize, alpha: T, a: [*]const T, lda: isize, b: [*]const T, ldb: isize, beta: Float(T), c: [*]T, ldc: isize) Error!void {
    comptime complexType(T, "her2k");
    return rankImpl(allocator, T, true, uplo, trans, n, k, alpha, a, lda, b, ldb, .init(beta, 0), c, ldc);
}

fn rankImpl(allocator: Allocator, comptime T: type, comptime herm: bool, uplo: Uplo, trans: Trans, n: isize, k: isize, alpha: T, a: [*]const T, lda: isize, b: ?[*]const T, ldb: isize, beta: T, c: [*]T, ldc: isize) Error!void {
    // real `c` is `t`; complex `syrk` takes `t`, `herk` `c`, and none `r`
    const t = switch (trans) {
        .n => false,
        .t => if (herm) return error.InvalidTrans else true,
        .c => if (isComplex(T) and !herm) return error.InvalidTrans else true,
        .r => return error.InvalidTrans,
    };
    if (n < 0) return error.InvalidN;
    if (k < 0) return error.InvalidK;
    const rows = if (t) k else n;
    const cols = if (t) n else k;
    if (lda < @max(1, rows)) return error.InvalidLda;
    if (b != null and ldb < @max(1, rows)) return error.InvalidLdb;
    if (ldc < @max(1, n)) return error.InvalidLdc;
    if (n == 0) return;
    if (!util.fitsMatrix(rows, cols, lda, @sizeOf(T))) return error.InvalidLda;
    if (b != null and !util.fitsMatrix(rows, cols, ldb, @sizeOf(T))) return error.InvalidLdb;
    if (!util.fitsMatrix(n, n, ldc, @sizeOf(T))) return error.InvalidLdc;
    try level3.syrk(allocator, T, table(T), herm, uplo == .u, t, @intCast(n), @intCast(k), alpha, a, @intCast(lda), b, @intCast(ldb), beta, c, @intCast(ldc));
}

const no_table: Vtable = .{};

// ---- tests -----------------------------------------------------------------

const testing = std.testing;
const check = @import("check.zig");

test "interface: gemm reports the first bad argument and does nothing for empty problems" {
    const fa = testing.failing_allocator;
    var a: [16]f64 = @splat(1);
    var c: [16]f64 = @splat(std.math.nan(f64));
    const nan = std.math.nan(f64);
    try testing.expectError(error.InvalidM, gemm(fa, f64, .n, .n, -1, -1, -1, 1, &a, 0, &a, 0, 0, &c, 0));
    try testing.expectError(error.InvalidN, gemm(fa, f64, .n, .n, 0, -1, -1, 1, &a, 0, &a, 0, 0, &c, 0));
    try testing.expectError(error.InvalidK, gemm(fa, f64, .n, .n, 0, 0, -1, 1, &a, 0, &a, 0, 0, &c, 0));
    try testing.expectError(error.InvalidLda, gemm(fa, f64, .n, .n, 3, 2, 4, 1, &a, 2, &a, 1, 0, &c, 1));
    try testing.expectError(error.InvalidLda, gemm(fa, f64, .t, .n, 3, 2, 4, 1, &a, 3, &a, 1, 0, &c, 1));
    try testing.expectError(error.InvalidLdb, gemm(fa, f64, .n, .n, 3, 2, 4, 1, &a, 3, &a, 3, 0, &c, 1));
    try testing.expectError(error.InvalidLdb, gemm(fa, f64, .n, .c, 3, 2, 4, 1, &a, 3, &a, 1, 0, &c, 1));
    try testing.expectError(error.InvalidLdc, gemm(fa, f64, .n, .n, 3, 2, 4, 1, &a, 3, &a, 4, 0, &c, 2));
    try testing.expectError(error.InvalidLdc, gemm(fa, f64, .n, .n, 0, 2, 0, 1, &a, 1, &a, 1, 0, &c, 0));
    // empty: nothing read or written, nothing allocated
    try gemm(fa, f64, .n, .n, 0, 3, 2, nan, &a, 1, &a, 2, 0, &c, 1);
    try gemm(fa, f64, .n, .n, 3, 0, 2, nan, &a, 3, &a, 2, 0, &c, 3);
    try testing.expect(std.math.isNan(c[0]));
    // k == 0 only scales C, and allocates nothing
    try gemm(fa, f64, .n, .n, 2, 2, 0, nan, &a, 2, &a, 1, 0, &c, 2);
    try testing.expectEqualSlices(f64, &.{ 0, 0, 0, 0 }, c[0..4]);
}

const Iface = @This();
const all = .{ f32, f64, types.Complex32, types.Complex64 };
const unslotted = .{ f16, f128, std.math.Complex(f16), std.math.Complex(f128) };

test "interface: gemm matches the wide reference through the global vtable" {
    inline for (all) |T| try check.interfaceCases(T, Iface, .gemm);
}

test "interface: gemm over types the vtable has no slots for runs the portable kernels" {
    inline for (unslotted) |T| try check.interfaceCases(T, Iface, .gemm);
}

test "interface: a matrix past the address space is its leading dimension's error" {
    const fa = testing.failing_allocator;
    var a: [4]f64 = @splat(1);
    const big = std.math.maxInt(isize) / 2;
    try testing.expectError(error.InvalidLda, gemm(fa, f64, .n, .n, 2, 2, 3, 1, &a, big, &a, 3, 0, &a, 2));
    try testing.expectError(error.InvalidLdb, gemm(fa, f64, .n, .n, 2, 2, 3, 1, &a, 2, &a, big, 0, &a, 2));
    try testing.expectError(error.InvalidLdc, gemm(fa, f64, .n, .n, 2, 2, 3, 1, &a, 2, &a, 3, 0, &a, big));
}

test "interface: symm and hemm match the wide reference" {
    inline for (all ++ unslotted) |T| {
        try check.interfaceCases(T, Iface, .symm);
        if (comptime isComplex(T)) try check.interfaceCases(T, Iface, .hemm);
    }
}

test "interface: symm and hemm report the first bad argument, and return for empty problems" {
    const fa = testing.failing_allocator;
    var a: [16]f64 = @splat(std.math.nan(f64));
    var c: [16]f64 = @splat(std.math.nan(f64));
    try testing.expectError(error.InvalidM, symm(fa, f64, .l, .u, -1, -1, 1, &a, 0, &a, 0, 0, &c, 0));
    try testing.expectError(error.InvalidN, symm(fa, f64, .l, .u, 0, -1, 1, &a, 0, &a, 0, 0, &c, 0));
    // the order of `S` is `m` on the left, `n` on the right
    try testing.expectError(error.InvalidLda, symm(fa, f64, .l, .u, 3, 0, 1, &a, 2, &a, 3, 0, &c, 3));
    try testing.expectError(error.InvalidLda, symm(fa, f64, .r, .l, 1, 3, 1, &a, 2, &a, 1, 0, &c, 1));
    try testing.expectError(error.InvalidLdb, symm(fa, f64, .r, .u, 2, 1, 1, &a, 1, &a, 1, 0, &c, 2));
    try testing.expectError(error.InvalidLdc, symm(fa, f64, .l, .u, 2, 1, 1, &a, 2, &a, 2, 0, &c, 1));
    try symm(fa, f64, .l, .u, 0, 3, std.math.nan(f64), &a, 1, &a, 1, 0, &c, 1);
    try symm(fa, f64, .r, .l, 3, 0, std.math.nan(f64), &a, 1, &a, 3, 0, &c, 3);
    try testing.expect(std.math.isNan(c[0]));
    var z: [4]types.Complex64 = @splat(.init(0, 0));
    try testing.expectError(error.InvalidLda, hemm(fa, types.Complex64, .r, .u, 1, 2, .init(1, 0), &z, 1, &z, 1, .init(0, 0), &z, 1));
}

test "interface: syrk, herk, syr2k and her2k match the wide reference" {
    inline for (all ++ unslotted) |T| {
        try check.interfaceCases(T, Iface, .syrk);
        try check.interfaceCases(T, Iface, .syr2k);
        if (comptime isComplex(T)) {
            try check.interfaceCases(T, Iface, .herk);
            try check.interfaceCases(T, Iface, .her2k);
        }
    }
}

test "interface: the rank-k updates report the first bad argument, and return for empty problems" {
    const fa = testing.failing_allocator;
    const Z = types.Complex64;
    var a: [16]f64 = @splat(std.math.nan(f64));
    var c: [16]f64 = @splat(std.math.nan(f64));
    var z: [4]Z = @splat(.init(0, 0));
    try testing.expectError(error.InvalidTrans, syrk(fa, f64, .u, .r, -1, 0, 1, &a, 1, 0, &c, 1));
    try testing.expectError(error.InvalidTrans, syrk(fa, Z, .u, .c, -1, 0, .init(1, 0), &z, 1, .init(0, 0), &z, 1));
    try testing.expectError(error.InvalidTrans, herk(fa, Z, .l, .t, -1, 0, 1, &z, 1, 0, &z, 1));
    try testing.expectError(error.InvalidN, syrk(fa, f64, .u, .c, -1, -1, 1, &a, 1, 0, &c, 1));
    try testing.expectError(error.InvalidK, syrk(fa, f64, .u, .n, 0, -1, 1, &a, 1, 0, &c, 1));
    try testing.expectError(error.InvalidLda, syrk(fa, f64, .u, .t, 1, 3, 1, &a, 2, 0, &c, 1));
    try testing.expectError(error.InvalidLdb, syr2k(fa, f64, .l, .n, 3, 1, 1, &a, 3, &a, 2, 0, &c, 3));
    try testing.expectError(error.InvalidLdc, syr2k(fa, f64, .l, .n, 3, 1, 1, &a, 3, &a, 3, 0, &c, 2));
    try testing.expectError(error.InvalidLdc, her2k(fa, Z, .u, .c, 2, 0, .init(1, 0), &z, 1, &z, 1, 0, &z, 1));
    try syrk(fa, f64, .u, .n, 0, 3, std.math.nan(f64), &a, 1, 0, &c, 1);
    try testing.expect(std.math.isNan(c[0]));
    // nothing to add and `beta` 1: a Hermitian diagonal keeps its imaginary part
    z[0] = .init(2, 5);
    try herk(fa, Z, .u, .n, 1, 0, 1, &z, 1, 1, &z, 1);
    try testing.expectEqual(Z.init(2, 5), z[0]);
    try herk(fa, Z, .u, .n, 1, 0, 1, &z, 1, 2, &z, 1);
    try testing.expectEqual(Z.init(4, 0), z[0]);
}
