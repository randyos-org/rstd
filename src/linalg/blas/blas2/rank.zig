//! The level-2 rank updates: `?SPR`, `?SPR2`, `?HPR`, `?HPR2`, `?SYR`,
//! `?SYR2`, `?HER`, `?HER2` (the complex symmetric ones OpenBLAS
//! extensions), each `interface/*.c` over `driver/level2/*_k.c`'s columns,
//! threaded above OpenBLAS's thresholds as `syr_thread.c` divides them.

const std = @import("std");
const blas = @import("../__root__.zig");
const types = @import("../types.zig");
const util = @import("../util.zig");
const accel = @import("../apple_accel.zig");
const threading = @import("../threading.zig");
const banded = @import("banded.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const BlasLong = types.BlasLong;
const Uplo = types.Uplo;
const Float = types.Float;
const Error = banded.Error;
const off = util.off;
const at = util.at;
const first = util.first;
const isZero = util.isZero;
const Scratch = banded.Scratch;
const axpyK = banded.axpyK;
const mulRounded = banded.mulRounded;
const scaleF = banded.scaleF;
const conj = banded.conj;
const columnAt = banded.columnAt;
const spans = banded.spans;

/// `run(ctx, c0, c1)` over columns `c0..c1` of a triangle of order `n`: all
/// of them at once below `family`'s threshold, else in chunks balanced by
/// their lengths, as `syr_thread.c` and the rest divide them.
fn byColumns(io: Io, comptime T: type, comptime family: threading.Family, comptime lower: bool, n: BlasLong, ctx: anytype, comptime run: fn (@TypeOf(ctx), BlasLong, BlasLong) void) Io.Cancelable!void {
    const nu: usize = @intCast(n);
    // multiples of 8 columns, as their `mask`
    const p = threading.plan(family, T, nu, nu, if (lower) .falling else .rising, 8);
    if (p.workers <= 1) return run(ctx, 0, n);
    const S = struct {
        fn one(c: @TypeOf(ctx), ch: threading.Chunk) void {
            run(c, @intCast(ch.start), @intCast(ch.end));
        }
    };
    return threading.forEach(io, &p, ctx, S.one);
}

/// `driver/level2/{spr,syr,zspr,zsyr,zhpr,zher}_k.c`'s columns `c0..c1`:
/// `A += alpha*x*x^T`, or `alpha*x*x^H` for `herm` (real `alpha`), on one
/// triangle in full (`lda`) or packed storage, over a unit-stride `X`. A
/// column whose `x` is zero is passed over, as the reference does;
/// OpenBLAS's `zher_k.c` and `zhpr_k.c` do not.
fn rankColumns(comptime T: type, comptime herm: bool, comptime lower: bool, comptime packed_: bool, n: BlasLong, alpha: if (herm) Float(T) else T, X: [*]T, a: [*]T, lda: BlasLong, c0: BlasLong, c1: BlasLong) void {
    var col = off(T, a, columnAt(lower, packed_, n, lda, c0));
    var i: BlasLong = c0;
    while (i < c1) : (i += 1) {
        const len = if (lower) n - i else i + 1;
        const xs = if (lower) off(T, X, i) else X;
        const xi = at(T, X, i).*;
        if (!isZero(T, xi)) axpyK(T, false, len, if (herm) scaleF(T, conj(T, xi), alpha) else mulRounded(T, alpha, xi), xs, col);
        if (herm) at(T, col, if (lower) 0 else i).im = 0;
        col = off(T, col, if (!packed_) (if (lower) lda + 1 else lda) else len);
    }
}

/// `driver/level2/{spr2,syr2,zspr2,zsyr2,zhpr2,zher2}_k.c`'s columns
/// `c0..c1`: `A += alpha*x*y^T + alpha*y*x^T`, or `alpha*x*y^H +
/// conj(alpha)*y*x^H` for `herm`, over unit-stride `X` and `Y`. A column
/// whose `x` and `y` are both zero is passed over, as the reference does;
/// OpenBLAS's drivers do not.
fn rank2Columns(comptime T: type, comptime herm: bool, comptime lower: bool, comptime packed_: bool, n: BlasLong, alpha: T, X: [*]T, Y: [*]T, a: [*]T, lda: BlasLong, c0: BlasLong, c1: BlasLong) void {
    var col = off(T, a, columnAt(lower, packed_, n, lda, c0));
    var i: BlasLong = c0;
    while (i < c1) : (i += 1) {
        const len = if (lower) n - i else i + 1;
        const xs = if (lower) off(T, X, i) else X;
        const ys = if (lower) off(T, Y, i) else Y;
        const xi = at(T, X, i).*;
        const yi = at(T, Y, i).*;
        if (!isZero(T, xi) or !isZero(T, yi)) {
            axpyK(T, false, len, if (herm) conj(T, mulRounded(T, alpha, xi)) else mulRounded(T, alpha, xi), ys, col);
            axpyK(T, false, len, if (herm) mulRounded(T, alpha, conj(T, yi)) else mulRounded(T, alpha, yi), xs, col);
        }
        if (herm) at(T, col, if (lower) 0 else i).im = 0;
        col = off(T, col, if (!packed_) (if (lower) lda + 1 else lda) else len);
    }
}

/// The rank updates' columns over unit-stride copies of strided vectors; a
/// unit-stride vector is taken as it is, so OpenBLAS's small-`n` paths
/// (`interface/syr.c` and the rest) are these same columns.
fn rankInterface(allocator: Allocator, io: Io, comptime T: type, comptime herm: bool, comptime packed_: bool, lower: bool, n: isize, alpha: if (herm) Float(T) else T, x: [*]T, incx: isize, a: [*]T, lda: isize) Error!void {
    if (n == 0) return;
    try spans(T, packed_, n, lda);
    if (!util.fits(n, incx, @sizeOf(T))) return error.InvalidIncx;
    if (comptime herm) {
        if (alpha == 0) return;
    } else if (isZero(T, alpha)) return;
    _ = blas.vtable();
    if (accel.syr(T, herm, packed_, if (lower) .l else .u, n, alpha, x, incx, a, lda)) return;
    var sx: Scratch(T) = .{};
    defer sx.deinit(allocator);
    const C = struct {
        alpha: if (herm) Float(T) else T,
        X: [*]T,
        a: [*]T,
        n: BlasLong,
        lda: BlasLong,
    };
    const c: C = .{ .alpha = alpha, .X = try sx.unit(allocator, n, first(T, x, n, incx), incx), .a = a, .n = n, .lda = lda };
    switch (lower) {
        inline else => |lo| {
            const S = struct {
                fn run(q: C, c0: BlasLong, c1: BlasLong) void {
                    rankColumns(T, herm, lo, packed_, q.n, q.alpha, q.X, q.a, q.lda, c0, c1);
                }
            };
            try byColumns(io, T, if (herm) .her else .syr, lo, n, c, S.run);
        },
    }
}

fn rank2Interface(allocator: Allocator, io: Io, comptime T: type, comptime herm: bool, comptime packed_: bool, lower: bool, n: isize, alpha: T, x: [*]T, incx: isize, y: [*]T, incy: isize, a: [*]T, lda: isize) Error!void {
    if (n == 0 or isZero(T, alpha)) return;
    try spans(T, packed_, n, lda);
    if (!util.fits(n, incx, @sizeOf(T))) return error.InvalidIncx;
    if (!util.fits(n, incy, @sizeOf(T))) return error.InvalidIncy;
    _ = blas.vtable();
    if (accel.syr2(T, herm, packed_, if (lower) .l else .u, n, alpha, x, incx, y, incy, a, lda)) return;
    var sx: Scratch(T) = .{};
    var sy: Scratch(T) = .{};
    defer sx.deinit(allocator);
    defer sy.deinit(allocator);
    const C = struct { alpha: T, X: [*]T, Y: [*]T, a: [*]T, n: BlasLong, lda: BlasLong };
    const c: C = .{
        .alpha = alpha,
        .X = try sx.unit(allocator, n, first(T, x, n, incx), incx),
        .Y = try sy.unit(allocator, n, first(T, y, n, incy), incy),
        .a = a,
        .n = n,
        .lda = lda,
    };
    switch (lower) {
        inline else => |lo| {
            const S = struct {
                fn run(q: C, c0: BlasLong, c1: BlasLong) void {
                    rank2Columns(T, herm, lo, packed_, q.n, q.alpha, q.X, q.Y, q.a, q.lda, c0, c1);
                }
            };
            try byColumns(io, T, if (herm) .her else .syr, lo, n, c, S.run);
        },
    }
}

/// `?SPR` (`interface/spr.c`, `zspr.c`): `A += alpha*x*x^T` on the packed
/// triangle. Columns where `x[i]` is zero are skipped.
pub fn spr(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, alpha: T, x: [*]T, incx: isize, a: [*]T) Error!void {
    if (uplo.isRev()) return error.InvalidUplo;
    if (n < 0) return error.InvalidN;
    if (incx == 0) return error.InvalidIncx;
    try rankInterface(allocator, io, T, false, true, uplo.isLower(), n, alpha, x, incx, a, 0);
}

/// `?SYR` (`interface/syr.c`, `zsyr.c`): `spr` for full storage.
pub fn syr(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, alpha: T, x: [*]T, incx: isize, a: [*]T, lda: isize) Error!void {
    if (uplo.isRev()) return error.InvalidUplo;
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, n)) return error.InvalidLda;
    if (incx == 0) return error.InvalidIncx;
    try rankInterface(allocator, io, T, false, false, uplo.isLower(), n, alpha, x, incx, a, lda);
}

/// `?HPR` (`interface/zhpr.c`): `A += alpha*x*x^H` for real `alpha`; the
/// diagonal's imaginary part becomes zero.
pub fn hpr(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, alpha: Float(T), x: [*]T, incx: isize, a: [*]T) Error!void {
    comptime util.requireComplex(T, "hpr");
    if (uplo.isRev()) return error.InvalidUplo;
    if (n < 0) return error.InvalidN;
    if (incx == 0) return error.InvalidIncx;
    try rankInterface(allocator, io, T, true, true, uplo.isLower(), n, alpha, x, incx, a, 0);
}

/// `?HER` (`interface/zher.c`): `hpr` for full storage.
pub fn her(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, alpha: Float(T), x: [*]T, incx: isize, a: [*]T, lda: isize) Error!void {
    comptime util.requireComplex(T, "her");
    if (uplo.isRev()) return error.InvalidUplo;
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, n)) return error.InvalidLda;
    if (incx == 0) return error.InvalidIncx;
    try rankInterface(allocator, io, T, true, false, uplo.isLower(), n, alpha, x, incx, a, lda);
}

/// `?SPR2` (`interface/spr2.c`, `zspr2.c`): `A += alpha*(x*y^T + y*x^T)` on
/// the packed triangle.
pub fn spr2(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, alpha: T, x: [*]T, incx: isize, y: [*]T, incy: isize, a: [*]T) Error!void {
    if (uplo.isRev()) return error.InvalidUplo;
    if (n < 0) return error.InvalidN;
    if (incx == 0) return error.InvalidIncx;
    if (incy == 0) return error.InvalidIncy;
    try rank2Interface(allocator, io, T, false, true, uplo.isLower(), n, alpha, x, incx, y, incy, a, 0);
}

/// `?SYR2` (`interface/syr2.c`, `zsyr2.c`): `spr2` for full storage.
pub fn syr2(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, alpha: T, x: [*]T, incx: isize, y: [*]T, incy: isize, a: [*]T, lda: isize) Error!void {
    if (uplo.isRev()) return error.InvalidUplo;
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, n)) return error.InvalidLda;
    if (incx == 0) return error.InvalidIncx;
    if (incy == 0) return error.InvalidIncy;
    try rank2Interface(allocator, io, T, false, false, uplo.isLower(), n, alpha, x, incx, y, incy, a, lda);
}

/// `?HPR2` (`interface/zhpr2.c`): `A += alpha*x*y^H + conj(alpha)*y*x^H`; the
/// diagonal's imaginary part becomes zero.
pub fn hpr2(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, alpha: T, x: [*]T, incx: isize, y: [*]T, incy: isize, a: [*]T) Error!void {
    comptime util.requireComplex(T, "hpr2");
    if (uplo.isRev()) return error.InvalidUplo;
    if (n < 0) return error.InvalidN;
    if (incx == 0) return error.InvalidIncx;
    if (incy == 0) return error.InvalidIncy;
    try rank2Interface(allocator, io, T, true, true, uplo.isLower(), n, alpha, x, incx, y, incy, a, 0);
}

/// `?HER2` (`interface/zher2.c`): `hpr2` for full storage.
pub fn her2(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, alpha: T, x: [*]T, incx: isize, y: [*]T, incy: isize, a: [*]T, lda: isize) Error!void {
    comptime util.requireComplex(T, "her2");
    if (uplo.isRev()) return error.InvalidUplo;
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, n)) return error.InvalidLda;
    if (incx == 0) return error.InvalidIncx;
    if (incy == 0) return error.InvalidIncy;
    try rank2Interface(allocator, io, T, true, false, uplo.isLower(), n, alpha, x, incx, y, incy, a, lda);
}
