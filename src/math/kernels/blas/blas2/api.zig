//! The public level-2 routines (issue #310): vectors as strided views
//! (`strided.zig`), matrices as a slice and its leading dimension,
//! column-major. The sizes come from the views: `m` and `n` from `y` and `x`
//! (by `trans`), `n` from `x` where the matrix is square. Each checks that the
//! views agree and the slice holds the matrix, then delegates to the
//! reference-shaped core (`interface.zig`, `triangular.zig`, `banded.zig`,
//! `rank.zig`), which the LAPACK routines and the ported tests call.

const std = @import("std");
const types = @import("../types.zig");
const strided = @import("../strided.zig");
const iface = @import("interface.zig");
const tri = @import("triangular.zig");
const band = @import("banded.zig");
const rank = @import("rank.zig");
const util = @import("../util.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const Trans = types.Trans;
const Uplo = types.Uplo;
const Diag = types.Diag;
const Strided = strided.Strided;
const ConstStrided = strided.ConstStrided;
const holds = util.holds;
const holdsPacked = util.holdsPacked;

pub const HalfPrecision = iface.HalfPrecision;

/// A level-2 routine's bad input: counts that disagree, a slice too short for
/// its matrix, a leading dimension below the matrix's rows; and
/// `error.OutOfMemory` and `error.Canceled`.
pub const Error = strided.LengthError || strided.SliceError || strided.LeadingDimensionError || Allocator.Error || Io.Cancelable;
/// `gemvHalf`'s, which allocates nothing.
pub const HalfError = strided.LengthError || strided.SliceError || strided.LeadingDimensionError || Io.Cancelable;

fn len(n: usize) isize {
    return @intCast(n);
}

/// The core's stride: any, for at most one element, where it has no meaning.
fn inc(v: anytype) isize {
    return if (v.n <= 1) 1 else v.stride;
}

fn ptr(comptime T: type, v: anytype) [*]T {
    return @constCast(v.data.ptr);
}

fn mat(comptime T: type, a: []const T) [*]T {
    return @constCast(a.ptr);
}

/// `op(A)`'s rows and columns from `y` and `x`: `y.n` by `x.n`.
fn dims(trans: Trans, x_n: usize, y_n: usize) [2]usize {
    return if (trans.transposed()) .{ x_n, y_n } else .{ y_n, x_n };
}

/// The core's errors past these checks: the allocator's, cancellation, and
/// `InvalidLda` should the core's rule be stricter than `holds`.
fn mapped(e: iface.Error) Error {
    return switch (e) {
        error.OutOfMemory => error.OutOfMemory,
        error.Canceled => error.Canceled,
        error.InvalidLda => error.LeadingDimensionTooSmall,
        error.InvalidUplo, error.InvalidTrans, error.InvalidDiag, error.InvalidM, error.InvalidN, error.InvalidK, error.InvalidIncx, error.InvalidIncy => unreachable,
    };
}

fn mappedHalf(e: iface.HalfError) HalfError {
    return switch (e) {
        error.Canceled => error.Canceled,
        error.InvalidLda => error.LeadingDimensionTooSmall,
        error.InvalidUplo, error.InvalidTrans, error.InvalidDiag, error.InvalidM, error.InvalidN, error.InvalidK, error.InvalidIncx, error.InvalidIncy => unreachable,
    };
}

// ---- general -----------------------------------------------------------------

/// `y := alpha op(A) x + beta y`, `A` being `m` by `n`: `op(A)` is `y.n` by
/// `x.n`.
pub fn gemv(allocator: Allocator, io: Io, comptime T: type, trans: Trans, alpha: T, a: []const T, lda: usize, x: ConstStrided(T), beta: T, y: Strided(T)) Error!void {
    const mn = dims(trans, x.n, y.n);
    try holds(a.len, mn[0], mn[1], lda);
    iface.gemv(allocator, io, T, trans, len(mn[0]), len(mn[1]), alpha, mat(T, a), len(lda), ptr(T, x), inc(x), beta, ptr(T, y), inc(y)) catch |e| return mapped(e);
}

/// `gemv` in half precision: `kind` picks bfloat16 or `f16` in, `f32` or
/// bfloat16 out.
pub fn gemvHalf(io: Io, comptime kind: HalfPrecision, trans: Trans, alpha: kind.Out(), a: []const kind.In(), lda: usize, x: ConstStrided(kind.In()), beta: kind.Out(), y: Strided(kind.Out())) HalfError!void {
    const mn = dims(trans, x.n, y.n);
    try holds(a.len, mn[0], mn[1], lda);
    iface.gemvHalf(io, kind, trans, len(mn[0]), len(mn[1]), alpha, mat(kind.In(), a), len(lda), ptr(kind.In(), x), inc(x), beta, ptr(kind.Out(), y), inc(y)) catch |e| return mappedHalf(e);
}

/// `A := alpha x yᵀ + A`, `A` being `x.n` by `y.n`.
pub fn ger(allocator: Allocator, io: Io, comptime T: type, alpha: T, x: ConstStrided(T), y: ConstStrided(T), a: []T, lda: usize) Error!void {
    try holds(a.len, x.n, y.n, lda);
    iface.ger(allocator, io, T, len(x.n), len(y.n), alpha, ptr(T, x), inc(x), ptr(T, y), inc(y), a.ptr, len(lda)) catch |e| return mapped(e);
}

/// `A := alpha x yᴴ + A`. Complex only.
pub fn gerc(allocator: Allocator, io: Io, comptime T: type, alpha: T, x: ConstStrided(T), y: ConstStrided(T), a: []T, lda: usize) Error!void {
    try holds(a.len, x.n, y.n, lda);
    iface.gerc(allocator, io, T, len(x.n), len(y.n), alpha, ptr(T, x), inc(x), ptr(T, y), inc(y), a.ptr, len(lda)) catch |e| return mapped(e);
}

/// `gemv` on a band matrix of `kl` subdiagonals and `ku` superdiagonals, in
/// band storage (`kl + ku + 1` rows).
pub fn gbmv(allocator: Allocator, io: Io, comptime T: type, trans: Trans, kl: usize, ku: usize, alpha: T, a: []const T, lda: usize, x: ConstStrided(T), beta: T, y: Strided(T)) Error!void {
    const mn = dims(trans, x.n, y.n);
    try holds(a.len, kl + ku + 1, mn[1], lda);
    band.gbmv(allocator, io, T, trans, len(mn[0]), len(mn[1]), len(kl), len(ku), alpha, mat(T, a), len(lda), ptr(T, x), inc(x), beta, ptr(T, y), inc(y)) catch |e| return mapped(e);
}

// ---- symmetric and Hermitian -------------------------------------------------

/// `y := alpha A x + beta y`, `A` symmetric (`uplo`'s triangle), order `x.n`.
pub fn symv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, alpha: T, a: []const T, lda: usize, x: ConstStrided(T), beta: T, y: Strided(T)) Error!void {
    try y.sameLength(x);
    try holds(a.len, x.n, x.n, lda);
    iface.symv(allocator, io, T, uplo, len(x.n), alpha, mat(T, a), len(lda), ptr(T, x), inc(x), beta, ptr(T, y), inc(y)) catch |e| return mapped(e);
}

/// `symv` with `A` Hermitian. Complex only.
pub fn hemv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, alpha: T, a: []const T, lda: usize, x: ConstStrided(T), beta: T, y: Strided(T)) Error!void {
    try y.sameLength(x);
    try holds(a.len, x.n, x.n, lda);
    iface.hemv(allocator, io, T, uplo, len(x.n), alpha, mat(T, a), len(lda), ptr(T, x), inc(x), beta, ptr(T, y), inc(y)) catch |e| return mapped(e);
}

/// `symv` on a band of `k` off-diagonals.
pub fn sbmv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, k: usize, alpha: T, a: []const T, lda: usize, x: ConstStrided(T), beta: T, y: Strided(T)) Error!void {
    try y.sameLength(x);
    try holds(a.len, k + 1, x.n, lda);
    band.sbmv(allocator, io, T, uplo, len(x.n), len(k), alpha, mat(T, a), len(lda), ptr(T, x), inc(x), beta, ptr(T, y), inc(y)) catch |e| return mapped(e);
}

/// `hemv` on a band of `k` off-diagonals. Complex only.
pub fn hbmv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, k: usize, alpha: T, a: []const T, lda: usize, x: ConstStrided(T), beta: T, y: Strided(T)) Error!void {
    try y.sameLength(x);
    try holds(a.len, k + 1, x.n, lda);
    band.hbmv(allocator, io, T, uplo, len(x.n), len(k), alpha, mat(T, a), len(lda), ptr(T, x), inc(x), beta, ptr(T, y), inc(y)) catch |e| return mapped(e);
}

/// `symv` on a packed triangle.
pub fn spmv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, alpha: T, ap: []const T, x: ConstStrided(T), beta: T, y: Strided(T)) Error!void {
    try y.sameLength(x);
    try holdsPacked(ap.len, x.n);
    band.spmv(allocator, io, T, uplo, len(x.n), alpha, mat(T, ap), ptr(T, x), inc(x), beta, ptr(T, y), inc(y)) catch |e| return mapped(e);
}

/// `hemv` on a packed triangle. Complex only.
pub fn hpmv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, alpha: T, ap: []const T, x: ConstStrided(T), beta: T, y: Strided(T)) Error!void {
    try y.sameLength(x);
    try holdsPacked(ap.len, x.n);
    band.hpmv(allocator, io, T, uplo, len(x.n), alpha, mat(T, ap), ptr(T, x), inc(x), beta, ptr(T, y), inc(y)) catch |e| return mapped(e);
}

// ---- triangular ----------------------------------------------------------------

/// `x := op(A) x`, `A` triangular, order `x.n`.
pub fn trmv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, diag: Diag, a: []const T, lda: usize, x: Strided(T)) Error!void {
    try holds(a.len, x.n, x.n, lda);
    tri.trmv(allocator, io, T, uplo, trans, diag, len(x.n), mat(T, a), len(lda), ptr(T, x), inc(x)) catch |e| return mapped(e);
}

/// `x := op(A)⁻¹ x`, `A` triangular, order `x.n`.
pub fn trsv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, diag: Diag, a: []const T, lda: usize, x: Strided(T)) Error!void {
    try holds(a.len, x.n, x.n, lda);
    tri.trsv(allocator, io, T, uplo, trans, diag, len(x.n), mat(T, a), len(lda), ptr(T, x), inc(x)) catch |e| return mapped(e);
}

/// `trmv` on a band of `k` off-diagonals.
pub fn tbmv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, diag: Diag, k: usize, a: []const T, lda: usize, x: Strided(T)) Error!void {
    try holds(a.len, k + 1, x.n, lda);
    tri.tbmv(allocator, io, T, uplo, trans, diag, len(x.n), len(k), mat(T, a), len(lda), ptr(T, x), inc(x)) catch |e| return mapped(e);
}

/// `trsv` on a band of `k` off-diagonals.
pub fn tbsv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, diag: Diag, k: usize, a: []const T, lda: usize, x: Strided(T)) Error!void {
    try holds(a.len, k + 1, x.n, lda);
    tri.tbsv(allocator, io, T, uplo, trans, diag, len(x.n), len(k), mat(T, a), len(lda), ptr(T, x), inc(x)) catch |e| return mapped(e);
}

/// `trmv` on a packed triangle.
pub fn tpmv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, diag: Diag, ap: []const T, x: Strided(T)) Error!void {
    try holdsPacked(ap.len, x.n);
    tri.tpmv(allocator, io, T, uplo, trans, diag, len(x.n), mat(T, ap), ptr(T, x), inc(x)) catch |e| return mapped(e);
}

/// `trsv` on a packed triangle.
pub fn tpsv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, diag: Diag, ap: []const T, x: Strided(T)) Error!void {
    try holdsPacked(ap.len, x.n);
    tri.tpsv(allocator, io, T, uplo, trans, diag, len(x.n), mat(T, ap), ptr(T, x), inc(x)) catch |e| return mapped(e);
}

// ---- rank updates --------------------------------------------------------------

/// `A := alpha x xᵀ + A`, symmetric (`uplo`'s triangle), order `x.n`.
pub fn syr(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, alpha: T, x: ConstStrided(T), a: []T, lda: usize) Error!void {
    try holds(a.len, x.n, x.n, lda);
    rank.syr(allocator, io, T, uplo, len(x.n), alpha, ptr(T, x), inc(x), a.ptr, len(lda)) catch |e| return mapped(e);
}

/// `A := alpha x xᴴ + A`, Hermitian, real `alpha`. Complex only.
pub fn her(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, alpha: Float(T), x: ConstStrided(T), a: []T, lda: usize) Error!void {
    try holds(a.len, x.n, x.n, lda);
    rank.her(allocator, io, T, uplo, len(x.n), alpha, ptr(T, x), inc(x), a.ptr, len(lda)) catch |e| return mapped(e);
}

/// `syr` on a packed triangle.
pub fn spr(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, alpha: T, x: ConstStrided(T), ap: []T) Error!void {
    try holdsPacked(ap.len, x.n);
    rank.spr(allocator, io, T, uplo, len(x.n), alpha, ptr(T, x), inc(x), ap.ptr) catch |e| return mapped(e);
}

/// `her` on a packed triangle. Complex only.
pub fn hpr(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, alpha: Float(T), x: ConstStrided(T), ap: []T) Error!void {
    try holdsPacked(ap.len, x.n);
    rank.hpr(allocator, io, T, uplo, len(x.n), alpha, ptr(T, x), inc(x), ap.ptr) catch |e| return mapped(e);
}

/// `A := alpha x yᵀ + alpha y xᵀ + A`, symmetric, order `x.n`.
pub fn syr2(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, alpha: T, x: ConstStrided(T), y: ConstStrided(T), a: []T, lda: usize) Error!void {
    try x.sameLength(y);
    try holds(a.len, x.n, x.n, lda);
    rank.syr2(allocator, io, T, uplo, len(x.n), alpha, ptr(T, x), inc(x), ptr(T, y), inc(y), a.ptr, len(lda)) catch |e| return mapped(e);
}

/// `A := alpha x yᴴ + conj(alpha) y xᴴ + A`, Hermitian. Complex only.
pub fn her2(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, alpha: T, x: ConstStrided(T), y: ConstStrided(T), a: []T, lda: usize) Error!void {
    try x.sameLength(y);
    try holds(a.len, x.n, x.n, lda);
    rank.her2(allocator, io, T, uplo, len(x.n), alpha, ptr(T, x), inc(x), ptr(T, y), inc(y), a.ptr, len(lda)) catch |e| return mapped(e);
}

/// `syr2` on a packed triangle.
pub fn spr2(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, alpha: T, x: ConstStrided(T), y: ConstStrided(T), ap: []T) Error!void {
    try x.sameLength(y);
    try holdsPacked(ap.len, x.n);
    rank.spr2(allocator, io, T, uplo, len(x.n), alpha, ptr(T, x), inc(x), ptr(T, y), inc(y), ap.ptr) catch |e| return mapped(e);
}

/// `her2` on a packed triangle. Complex only.
pub fn hpr2(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, alpha: T, x: ConstStrided(T), y: ConstStrided(T), ap: []T) Error!void {
    try x.sameLength(y);
    try holdsPacked(ap.len, x.n);
    rank.hpr2(allocator, io, T, uplo, len(x.n), alpha, ptr(T, x), inc(x), ptr(T, y), inc(y), ap.ptr) catch |e| return mapped(e);
}

// ---- tests -------------------------------------------------------------------

const check = @import("../_test/check.zig");
const testing = std.testing;

test "api: gemv on views gives the core's bits, transposed and with a backward x" {
    const al = testing.allocator;
    var a: [12]f64 = undefined;
    check.fill(f64, &a, 1);
    var xs: [7]f64 = undefined;
    check.fill(f64, &xs, 2);
    var y1: [4]f64 = undefined;
    check.fill(f64, &y1, 3);
    var y2 = y1;
    // A is 4 x 3 at lda 4; x is 3 elements, stride -3
    try gemv(al, check.io, f64, .n, 1.5, &a, 4, try .init(&xs, 3, -3), 0.5, .of(&y1));
    try iface.gemv(al, check.io, f64, .n, 4, 3, 1.5, &a, 4, &xs, -3, 0.5, &y2, 1);
    try testing.expectEqualSlices(f64, &y2, &y1);
    var z1: [3]f64 = @splat(1);
    var z2 = z1;
    try gemv(al, check.io, f64, .t, 1, &a, 4, .of(y1[0..4]), 0, .of(&z1));
    try iface.gemv(al, check.io, f64, .t, 4, 3, 1, &a, 4, &y1, 1, 0, &z2, 1);
    try testing.expectEqualSlices(f64, &z2, &z1);
}

test "api: a short slice, a small leading dimension and counts that disagree are errors" {
    const al = testing.allocator;
    var a: [12]f64 = @splat(1);
    const x = [_]f64{ 1, 2, 3 };
    var y: [4]f64 = @splat(0);
    // 4 x 3 needs 12 elements at lda 4
    try testing.expectError(error.SliceTooShort, gemv(al, check.io, f64, .n, 1, a[0..11], 4, .of(&x), 0, .of(&y)));
    try testing.expectError(error.LeadingDimensionTooSmall, gemv(al, check.io, f64, .n, 1, &a, 3, .of(&x), 0, .of(&y)));
    try testing.expectError(error.LengthMismatch, symv(al, check.io, f64, .u, 1, &a, 3, .of(&x), 0, .of(&y)));
    try testing.expectError(error.SliceTooShort, spmv(al, check.io, f64, .u, 1, a[0..5], .of(&x), 0, .of(y[0..3])));
    try testing.expectError(error.LeadingDimensionTooSmall, sbmv(al, check.io, f64, .u, 2, 1, &a, 2, .of(&x), 0, .of(y[0..3])));
    try spmv(al, check.io, f64, .u, 1, a[0..6], .of(&x), 0, .of(y[0..3]));
}

test "api: the rank updates and triangular solves take their order from x" {
    const al = testing.allocator;
    var a1: [9]f64 = undefined;
    check.fill(f64, &a1, 4);
    var a2 = a1;
    const x = [_]f64{ 1, -2, 0.5 };
    try syr(al, check.io, f64, .l, 2, .of(&x), &a1, 3);
    try rank.syr(al, check.io, f64, .l, 3, 2, @constCast(&x), 1, &a2, 3);
    try testing.expectEqualSlices(f64, &a2, &a1);
    for (0..3) |i| a1[i * 4] = 4;
    var b1 = x;
    var b2 = x;
    try trsv(al, check.io, f64, .l, .n, .n, &a1, 3, .of(&b1));
    try tri.trsv(al, check.io, f64, .l, .n, .n, 3, &a1, 3, &b2, 1);
    try testing.expectEqualSlices(f64, &b2, &b1);
}
