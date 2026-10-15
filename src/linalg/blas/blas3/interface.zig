//! The level-3 interface: the checks and quick returns of OpenBLAS's
//! `interface/gemm.c` before the driver (`level3.zig`). Each takes the
//! element type `T` explicitly, any real float type (the vtable's kernels
//! for those it has slots for, the portable ones built for `T` otherwise),
//! and sizes and leading dimensions as `isize`; matrices are column-major.
//! A matrix that would not fit the address space is its leading
//! dimension's error. Bad arguments return an `ArgError` in the
//! order the reference numbers them for `xerbla`, so the first bad one wins.
//! Leading dimensions are checked against `max(1, rows)`, as the reference
//! BLAS does (OpenBLAS lets `0` through for empty matrices). Packing space
//! comes from `allocator`; small problems use the stack. Single-threaded.

const std = @import("std");
const blas = @import("../__root__.zig");
const level3 = @import("level3.zig");
const portable = @import("../portable.zig");
const util = @import("../util.zig");
const Vtable = @import("../Vtable.zig");

const Allocator = std.mem.Allocator;

/// `TRANSA`, `TRANSB`: `n` is the matrix, `t` its transpose, `c` its
/// conjugate transpose and `r` its conjugate. Real `r` is `n` and `c` is `t`.
pub const Trans = enum { n, t, r, c };

fn transposed(t: Trans) bool {
    return switch (t) {
        .t, .c => true,
        .n, .r => false,
    };
}

/// The reference's `xerbla` argument, by name.
pub const ArgError = error{
    InvalidM,
    InvalidN,
    InvalidK,
    InvalidLda,
    InvalidLdb,
    InvalidLdc,
};
pub const Error = ArgError || Allocator.Error;

/// `?GEMM`: `C := alpha*op(A)*op(B) + beta*C`, `op(A)` being `m x k`,
/// `op(B)` `k x n`. `C` is scaled by `beta` first (a zero `beta` never reads
/// it), and `A` and `B` are not read for a zero `alpha` or `k`.
pub fn gemm(allocator: Allocator, comptime T: type, transa: Trans, transb: Trans, m: isize, n: isize, k: isize, alpha: T, a: [*]const T, lda: isize, b: [*]const T, ldb: isize, beta: T, c: [*]T, ldc: isize) Error!void {
    if (@typeInfo(T) != .float) @compileError("gemm takes a real float type; complex is not ported yet");
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
    // a type without slots runs the portable kernels and never needs the table
    const vt: *const Vtable = if (comptime portable.letter(T) != null) blas.vtable() else &no_table;
    try level3.gemm(allocator, T, vt, ta, tb, @intCast(m), @intCast(n), @intCast(k), alpha, a, @intCast(lda), b, @intCast(ldb), beta, c, @intCast(ldc));
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

test "interface: gemm matches the wide reference through the global vtable" {
    try check.interfaceCases(f32, Trans, gemm);
    try check.interfaceCases(f64, Trans, gemm);
}

test "interface: gemm over types the vtable has no slots for runs the portable kernels" {
    try check.interfaceCases(f16, Trans, gemm);
    try check.interfaceCases(f128, Trans, gemm);
}

test "interface: a matrix past the address space is its leading dimension's error" {
    const fa = testing.failing_allocator;
    var a: [4]f64 = @splat(1);
    const big = std.math.maxInt(isize) / 2;
    try testing.expectError(error.InvalidLda, gemm(fa, f64, .n, .n, 2, 2, 3, 1, &a, big, &a, 3, 0, &a, 2));
    try testing.expectError(error.InvalidLdb, gemm(fa, f64, .n, .n, 2, 2, 3, 1, &a, 2, &a, big, 0, &a, 2));
    try testing.expectError(error.InvalidLdc, gemm(fa, f64, .n, .n, 2, 2, 3, 1, &a, 2, &a, 3, 0, &a, big));
}
