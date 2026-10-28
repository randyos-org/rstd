//! The level-2 banded, packed and rank-update interface: `?GBMV`, `?SBMV`,
//! `?HBMV`, `?SPMV`, `?HPMV`, `?SPR`, `?SPR2`, `?HPR`, `?HPR2`, `?SYR`,
//! `?SYR2`, `?HER`, `?HER2`. Each is `interface/*.c` (validation, `beta`
//! scaling, negative strides) over `driver/level2/*_k.c` (one column at a
//! time through the level-1 kernels: copy, scal, dot, axpy, axpyc).
//!
//! The complex symmetric forms (`ZSBMV`, `ZSPMV`, `ZSPR`, `ZSPR2`, `ZSYR`,
//! `ZSYR2`) are OpenBLAS extensions. `V` and `M` in the Hermitian `uplo` are
//! its `HEMVREV` forms (upper and lower, the transpose of the Hermitian
//! matrix); the Fortran interface accepts them.
//!
//! Element type `T` is any float type or its `Complex`; sizes and strides
//! are `isize`. Vectors with a negative stride start at the
//! pointer given, as in the Fortran interface. No vtable slot exists for
//! these routines; every kernel used is a level-1 one (`portable.kernel`),
//! and `gbmv` uses no gemv kernel.
//!
//! There is no xerbla here: bad arguments return `interface.ArgError`
//! (`info` is the first failing argument in the reference's check order).
//! Each routine takes an allocator, for copies of strided vectors beyond
//! `stack_bytes`, and `io`, for the threading to come; every one is serial.

const std = @import("std");
const sysinfo = @import("builtin");
const blas = @import("../__root__.zig");
const types = @import("../types.zig");
const util = @import("../util.zig");
const interface = @import("interface.zig");
const portable = @import("../portable.zig");
const accel = @import("../apple_accel.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const BlasLong = types.BlasLong;
const Trans = types.Trans;
const Uplo = types.Uplo;

pub const Error = interface.Error;

/// Bytes of each vector copy kept on the stack (OpenBLAS's `STACK_ALLOC`).
pub const stack_bytes = interface.stack_bytes;

// ---- element helpers ---------------------------------------------------------

const isComplex = types.isComplex;
const Float = types.Float;

fn raw(comptime T: type, p: [*]T) [*]Float(T) {
    return @ptrCast(p);
}

/// `p + i` for any sign of `i`.
fn at(comptime T: type, p: [*]T, i: BlasLong) [*]T {
    return @ptrCast(util.at(T, p, i));
}

fn el(comptime T: type, p: [*]T, i: BlasLong) *T {
    return util.at(T, p, i);
}

/// Where a vector of `n` with stride `inc` starts: its last element when
/// `inc < 0`.
fn first(comptime T: type, p: [*]T, n: BlasLong, inc: BlasLong) [*]T {
    return if (inc < 0) at(T, p, -(n - 1) * inc) else p;
}

fn zero(comptime T: type) T {
    return if (comptime isComplex(T)) .init(0, 0) else 0;
}

fn mul(comptime T: type, a: T, b: T) T {
    return if (comptime isComplex(T)) a.mul(b) else a * b;
}

fn add(comptime T: type, a: T, b: T) T {
    return if (comptime isComplex(T)) a.add(b) else a + b;
}

fn conj(comptime T: type, a: T) T {
    return if (comptime isComplex(T)) a.conjugate() else a;
}

fn scaleF(comptime T: type, a: T, f: Float(T)) T {
    return if (comptime isComplex(T)) T.init(a.re * f, a.im * f) else a * f;
}

fn isZero(comptime T: type, a: T) bool {
    return if (comptime isComplex(T)) a.re == 0 and a.im == 0 else a == 0;
}

fn isOne(comptime T: type, a: T) bool {
    return if (comptime isComplex(T)) a.re == 1 and a.im == 0 else a == 1;
}

// ---- level-1 kernels, unit stride --------------------------------------------

/// `y += alpha*x`, or `alpha*conj(x)`; the kernels take `n >= 1`.
fn axpyK(comptime T: type, comptime cj: bool, n: BlasLong, alpha: T, x: [*]T, y: [*]T) void {
    if (n <= 0) return;
    if (comptime isComplex(T)) {
        _ = portable.fetch(T, if (cj) .axpyc_k else .axpy_k, blas.table)(n, 0, 0, alpha.re, alpha.im, raw(T, x), 1, raw(T, y), 1, raw(T, y), 0);
    } else {
        _ = portable.fetch(T, .axpy_k, blas.table)(n, 0, 0, alpha, x, 1, y, 1, y, 0);
    }
}

/// `x . y`, or `conj(x) . y`; `cj` is ignored for real `T`.
fn dotK(comptime T: type, comptime cj: bool, n: BlasLong, x: [*]T, y: [*]T) T {
    if (n <= 0) return zero(T);
    if (comptime isComplex(T)) {
        return types.complexOf(portable.fetch(T, if (cj) .dotc_k else .dotu_k, blas.table)(n, raw(T, x), 1, raw(T, y), 1));
    }
    return portable.fetch(T, .dot_k, blas.table)(n, x, 1, y, 1);
}

fn copyK(comptime T: type, n: BlasLong, x: [*]T, incx: BlasLong, y: [*]T, incy: BlasLong) void {
    if (n <= 0) return;
    _ = portable.fetch(T, .copy_k, blas.table)(n, raw(T, x), incx, raw(T, y), incy);
}

/// `y := beta*y`, with the interface's flag 0: a zero `beta` stores zeros.
fn scalK(comptime T: type, n: BlasLong, beta: T, y: [*]T, incy: BlasLong) void {
    if (n <= 0) return;
    const p = raw(T, y);
    if (comptime isComplex(T)) {
        _ = portable.fetch(T, .scal_k, blas.table)(n, 0, 0, beta.re, beta.im, p, incy, p, 0, p, 0);
    } else {
        _ = portable.fetch(T, .scal_k, blas.table)(n, 0, 0, beta, p, incy, p, 0, p, 0);
    }
}

// ---- unit-stride copies ------------------------------------------------------

/// Room for one vector copy: on the stack when it fits, else from the
/// allocator.
fn Scratch(comptime T: type) type {
    return struct {
        stack: [stack_bytes / @sizeOf(T)]T = undefined,
        heap: ?[]T = null,

        /// `v` (length `n`, stride `inc`) as a unit-stride vector: `v` itself
        /// at stride 1, else a copy.
        fn unit(self: *@This(), allocator: Allocator, n: BlasLong, v: [*]T, inc: BlasLong) Allocator.Error![*]T {
            if (inc == 1) return v;
            const len: usize = @intCast(n);
            const buf: []T = if (len <= self.stack.len) self.stack[0..len] else blk: {
                self.heap = try allocator.alloc(T, len);
                break :blk self.heap.?;
            };
            // poisoned in tests, as `interface.zig`'s scratch is
            if (sysinfo.is_test) @memset(std.mem.sliceAsBytes(buf), 0xff);
            copyK(T, n, v, inc, buf.ptr, 1);
            return buf.ptr;
        }

        fn deinit(self: *@This(), allocator: Allocator) void {
            if (self.heap) |h| allocator.free(h);
        }
    };
}

// ---- drivers -----------------------------------------------------------------

/// `driver/level2/gbmv_k.c`, `zgbmv_k.c`: `y += alpha*op(A)*x` for an `m` by
/// `n` band matrix. The flags are `TRANS`, `CONJ` (of `A`) and `XCONJ` (of `x`).
fn gbmvDriver(comptime T: type, comptime tr: bool, comptime cj: bool, comptime xc: bool, allocator: Allocator, m: BlasLong, n: BlasLong, ku: BlasLong, kl: BlasLong, alpha: T, a: [*]T, lda: BlasLong, x: [*]T, incx: BlasLong, y: [*]T, incy: BlasLong) Allocator.Error!void {
    const leny = if (tr) n else m;
    const lenx = if (tr) m else n;
    var sy: Scratch(T) = .{};
    var sx: Scratch(T) = .{};
    defer sy.deinit(allocator);
    defer sx.deinit(allocator);
    const Y = try sy.unit(allocator, leny, y, incy);
    const X = try sx.unit(allocator, lenx, x, incx);
    var offset_u = ku;
    var offset_l = ku + m;
    var col = a;
    var i: BlasLong = 0;
    while (i < @min(n, m + ku)) : (i += 1) {
        const start = @max(offset_u, 0);
        const end = @min(offset_l, ku + kl + 1);
        const length = end - start;
        if (length > 0) {
            const ap = at(T, col, start);
            const xs = at(T, X, start - offset_u);
            if (comptime !tr) {
                const xi = el(T, X, i).*;
                axpyK(T, cj, length, mul(T, alpha, if (xc) conj(T, xi) else xi), ap, at(T, Y, start - offset_u));
            } else {
                const t = if (comptime !xc) dotK(T, cj, length, ap, xs) else dotK(T, !cj, length, xs, ap);
                const v = mul(T, alpha, if (comptime xc and cj) conj(T, t) else t);
                el(T, Y, i).* = add(T, el(T, Y, i).*, v);
            }
        }
        offset_u -= 1;
        offset_l -= 1;
        col = at(T, col, lda);
    }
    if (incy != 1) copyK(T, leny, Y, 1, y, incy);
}

/// Which matrix a stored triangle stands for.
const Kind = enum { symmetric, hermitian, hermitian_rev };

/// `driver/level2/{sbmv,spmv,zsbmv,zspmv,zhbmv,zhpmv}_k.c`: `y += alpha*S*x`
/// for `S` in band storage (`k` bands, `lda`) or packed storage.
fn mvDriver(comptime T: type, comptime kind: Kind, comptime lower: bool, comptime packed_: bool, allocator: Allocator, n: BlasLong, k: BlasLong, alpha: T, a: [*]T, lda: BlasLong, x: [*]T, incx: BlasLong, y: [*]T, incy: BlasLong) Allocator.Error!void {
    var sy: Scratch(T) = .{};
    var sx: Scratch(T) = .{};
    defer sy.deinit(allocator);
    defer sx.deinit(allocator);
    const Y = try sy.unit(allocator, n, y, incy);
    const X = try sx.unit(allocator, n, x, incx);
    const rev = kind == .hermitian_rev;
    var cs: BlasLong = 0; // start of column i in `a`
    var i: BlasLong = 0;
    while (i < n) : (i += 1) {
        const room: BlasLong = if (lower) n - 1 - i else i;
        const noff = if (packed_) room else @min(room, k);
        const diag = cs + (if (lower) 0 else if (packed_) i else k);
        const dp = at(T, a, diag);
        const off_a = if (lower) at(T, a, diag + 1) else at(T, a, diag - noff);
        const off_i = if (lower) i + 1 else i - noff; // the off-diagonals' rows
        const xi = el(T, X, i).*;
        const yi = el(T, Y, i);
        const ax = mul(T, alpha, xi);
        switch (kind) {
            .symmetric => if (packed_ and lower) {
                yi.* = add(T, yi.*, mul(T, alpha, dotK(T, false, n - i, dp, at(T, X, i))));
                axpyK(T, false, noff, ax, off_a, at(T, Y, off_i));
            } else {
                axpyK(T, false, noff + 1, ax, if (lower) dp else off_a, at(T, Y, if (lower) i else off_i));
                if (noff > 0) yi.* = add(T, yi.*, mul(T, alpha, dotK(T, false, noff, off_a, at(T, X, off_i))));
            },
            .hermitian, .hermitian_rev => {
                axpyK(T, rev, noff, ax, off_a, at(T, Y, off_i));
                const d = if (comptime isComplex(T)) dp[0].re else dp[0];
                yi.* = add(T, yi.*, mul(T, alpha, scaleF(T, xi, d)));
                if (noff > 0) yi.* = add(T, yi.*, mul(T, alpha, dotK(T, !rev, noff, off_a, at(T, X, off_i))));
            },
        }
        cs += if (!packed_) lda else if (lower) n - i else i + 1;
    }
    if (incy != 1) copyK(T, n, Y, 1, y, incy);
}

/// `driver/level2/{spr,syr,zspr,zsyr,zhpr,zher}_k.c`: `A += alpha*x*x^T`, or
/// `alpha*x*x^H` for the Hermitian kinds (real `alpha`), on one triangle in
/// full (`lda`) or packed storage. A column whose `x` is zero is passed over,
/// as the reference does; OpenBLAS's `zher_k.c` and `zhpr_k.c` do not.
fn rankDriver(comptime T: type, comptime kind: Kind, comptime lower: bool, comptime packed_: bool, allocator: Allocator, n: BlasLong, alpha: if (kind == .symmetric) T else Float(T), x: [*]T, incx: BlasLong, a: [*]T, lda: BlasLong) Allocator.Error!void {
    var sx: Scratch(T) = .{};
    defer sx.deinit(allocator);
    rankColumns(T, kind, lower, packed_, n, alpha, try sx.unit(allocator, n, x, incx), a, lda);
}

/// `rankDriver`'s columns over a unit-stride `X`.
fn rankColumns(comptime T: type, comptime kind: Kind, comptime lower: bool, comptime packed_: bool, n: BlasLong, alpha: if (kind == .symmetric) T else Float(T), X: [*]T, a: [*]T, lda: BlasLong) void {
    var col = a;
    var i: BlasLong = 0;
    while (i < n) : (i += 1) {
        const len = if (lower) n - i else i + 1;
        const xs = if (lower) at(T, X, i) else X;
        const xi = el(T, X, i).*;
        if (!isZero(T, xi)) switch (kind) {
            .symmetric => axpyK(T, false, len, mul(T, alpha, xi), xs, col),
            .hermitian => axpyK(T, false, len, scaleF(T, conj(T, xi), alpha), xs, col),
            .hermitian_rev => axpyK(T, true, len, scaleF(T, xi, alpha), xs, col),
        };
        if (kind != .symmetric) el(T, col, if (lower) 0 else i).im = 0;
        col = at(T, col, if (!packed_) (if (lower) lda + 1 else lda) else len);
    }
}

/// `driver/level2/{spr2,syr2,zspr2,zsyr2,zhpr2,zher2}_k.c`: `A += alpha*x*y^T
/// + alpha*y*x^T`, or `alpha*x*y^H + conj(alpha)*y*x^H` for the Hermitian
/// kinds. A column whose `x` and `y` are both zero is passed over, as the
/// reference does; OpenBLAS's drivers do not.
fn rank2Driver(comptime T: type, comptime kind: Kind, comptime lower: bool, comptime packed_: bool, allocator: Allocator, n: BlasLong, alpha: T, x: [*]T, incx: BlasLong, y: [*]T, incy: BlasLong, a: [*]T, lda: BlasLong) Allocator.Error!void {
    var sx: Scratch(T) = .{};
    var sy: Scratch(T) = .{};
    defer sx.deinit(allocator);
    defer sy.deinit(allocator);
    const X = try sx.unit(allocator, n, x, incx);
    const Y = try sy.unit(allocator, n, y, incy);
    rank2Columns(T, kind, lower, packed_, n, alpha, X, Y, a, lda);
}

/// `rank2Driver`'s columns over unit-stride `X` and `Y`.
fn rank2Columns(comptime T: type, comptime kind: Kind, comptime lower: bool, comptime packed_: bool, n: BlasLong, alpha: T, X: [*]T, Y: [*]T, a: [*]T, lda: BlasLong) void {
    var col = a;
    var i: BlasLong = 0;
    while (i < n) : (i += 1) {
        const len = if (lower) n - i else i + 1;
        const xs = if (lower) at(T, X, i) else X;
        const ys = if (lower) at(T, Y, i) else Y;
        const xi = el(T, X, i).*;
        const yi = el(T, Y, i).*;
        if (!isZero(T, xi) or !isZero(T, yi)) switch (kind) {
            .symmetric => {
                axpyK(T, false, len, mul(T, alpha, xi), ys, col);
                axpyK(T, false, len, mul(T, alpha, yi), xs, col);
            },
            .hermitian => {
                axpyK(T, false, len, conj(T, mul(T, alpha, xi)), ys, col);
                axpyK(T, false, len, mul(T, alpha, conj(T, yi)), xs, col);
            },
            .hermitian_rev => {
                axpyK(T, true, len, mul(T, alpha, xi), ys, col);
                axpyK(T, true, len, mul(T, conj(T, alpha), yi), xs, col);
            },
        };
        if (kind != .symmetric) el(T, col, if (lower) 0 else i).im = 0;
        col = at(T, col, if (!packed_) (if (lower) lda + 1 else lda) else len);
    }
}

// ---- interface ---------------------------------------------------------------

/// Refuses a matrix of order `n` that could not exist in memory
/// (`util.fitsMatrix`): `lda` columns of band or full storage, or the packed
/// triangle, bounded by `n` columns of `n`.
fn spans(comptime T: type, comptime packed_: bool, n: isize, lda: isize) interface.ArgError!void {
    if (packed_) {
        if (!util.fitsMatrix(n, n, n, @sizeOf(T))) return error.InvalidN;
    } else if (!util.fitsMatrix(lda, n, lda, @sizeOf(T))) return error.InvalidLda;
}

/// `?GBMV` (`interface/gbmv.c`, `zgbmv.c`): `y := alpha*op(A)*x + beta*y` for
/// an `m` by `n` matrix with `kl` sub- and `ku` superdiagonals in band
/// storage. `beta == 0` stores `alpha*op(A)*x` without reading `y`.
pub fn gbmv(allocator: Allocator, io: Io, comptime T: type, trans: Trans, m: isize, n: isize, kl: isize, ku: isize, alpha: T, a: [*]T, lda: isize, x: [*]T, incx: isize, beta: T, y: [*]T, incy: isize) Error!void {
    _ = io;
    if (!isComplex(T) and trans.conjX()) return error.InvalidTrans;
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    if (kl < 0 or ku < 0) return error.InvalidK;
    if (lda < kl + ku + 1) return error.InvalidLda;
    if (incx == 0) return error.InvalidIncx;
    if (incy == 0) return error.InvalidIncy;
    if (m == 0 or n == 0) return;
    _ = blas.vtable();
    const tr = trans.transposed();
    const lenx = if (tr) m else n;
    const leny = if (tr) n else m;
    if (!util.fitsMatrix(lda, n, lda, @sizeOf(T))) return error.InvalidLda;
    if (!util.fits(lenx, incx, @sizeOf(T))) return error.InvalidIncx;
    if (!util.fits(leny, incy, @sizeOf(T))) return error.InvalidIncy;
    if (accel.gbmv(T, trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)) return;
    if (!isOne(T, beta)) scalK(T, leny, beta, y, @intCast(@abs(incy)));
    if (isZero(T, alpha)) return;
    const px = first(T, x, lenx, incx);
    const py = first(T, y, leny, incy);
    if (comptime isComplex(T)) {
        switch (trans) {
            inline else => |t| try gbmvDriver(T, comptime t.transposed(), comptime t.conjA(), comptime t.conjX(), allocator, m, n, ku, kl, alpha, a, lda, px, incx, py, incy),
        }
    } else switch (tr) {
        inline else => |t| try gbmvDriver(T, t, false, false, allocator, m, n, ku, kl, alpha, a, lda, px, incx, py, incy),
    }
}

/// Shared tail of the matrix-vector products: `beta` scaling, `alpha == 0`,
/// start-of-vector pointers, and the uplo dispatch.
fn mvInterface(allocator: Allocator, io: Io, comptime T: type, comptime kind: Kind, comptime packed_: bool, lower: bool, n: isize, k: isize, alpha: T, a: [*]T, lda: isize, x: [*]T, incx: isize, beta: T, y: [*]T, incy: isize) Error!void {
    _ = io;
    if (n == 0) return;
    try spans(T, packed_, n, lda);
    if (!util.fits(n, incx, @sizeOf(T))) return error.InvalidIncx;
    if (!util.fits(n, incy, @sizeOf(T))) return error.InvalidIncy;
    _ = blas.vtable();
    // CBLAS has no `HEMVREV` triangles
    if (kind != .hermitian_rev and accel.sbmv(T, kind == .hermitian, packed_, if (lower) .l else .u, n, k, alpha, a, lda, x, incx, beta, y, incy)) return;
    const kn = n;
    if (!isOne(T, beta)) scalK(T, kn, beta, y, @intCast(@abs(incy)));
    if (isZero(T, alpha)) return;
    const px = first(T, x, kn, incx);
    const py = first(T, y, kn, incy);
    switch (lower) {
        inline else => |lo| try mvDriver(T, kind, lo, packed_, allocator, kn, k, alpha, a, lda, px, incx, py, incy),
    }
}

/// `?SBMV` (`interface/sbmv.c`, `zsbmv.c`): `y := alpha*S*x + beta*y` for the
/// symmetric `n` by `n` matrix of `k` bands in band storage.
pub fn sbmv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, k: isize, alpha: T, a: [*]T, lda: isize, x: [*]T, incx: isize, beta: T, y: [*]T, incy: isize) Error!void {
    if (uplo.isRev()) return error.InvalidUplo;
    if (n < 0) return error.InvalidN;
    if (k < 0) return error.InvalidK;
    if (lda < k + 1) return error.InvalidLda;
    if (incx == 0) return error.InvalidIncx;
    if (incy == 0) return error.InvalidIncy;
    try mvInterface(allocator, io, T, .symmetric, false, uplo.isLower(), n, k, alpha, a, lda, x, incx, beta, y, incy);
}

/// `?HBMV` (`interface/zhbmv.c`): `sbmv` for a Hermitian matrix; the
/// diagonal's imaginary part is not read.
pub fn hbmv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, k: isize, alpha: T, a: [*]T, lda: isize, x: [*]T, incx: isize, beta: T, y: [*]T, incy: isize) Error!void {
    if (comptime !isComplex(T)) @compileError("hbmv: " ++ @typeName(T) ++ " is not complex");
    if (n < 0) return error.InvalidN;
    if (k < 0) return error.InvalidK;
    if (lda < k + 1) return error.InvalidLda;
    if (incx == 0) return error.InvalidIncx;
    if (incy == 0) return error.InvalidIncy;
    switch (uplo.isRev()) {
        inline else => |rev| try mvInterface(allocator, io, T, if (rev) .hermitian_rev else .hermitian, false, uplo.isLower(), n, k, alpha, a, lda, x, incx, beta, y, incy),
    }
}

/// `?SPMV` (`interface/spmv.c`, `zspmv.c`): `sbmv` for packed storage.
pub fn spmv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, alpha: T, a: [*]T, x: [*]T, incx: isize, beta: T, y: [*]T, incy: isize) Error!void {
    if (uplo.isRev()) return error.InvalidUplo;
    if (n < 0) return error.InvalidN;
    if (incx == 0) return error.InvalidIncx;
    if (incy == 0) return error.InvalidIncy;
    try mvInterface(allocator, io, T, .symmetric, true, uplo.isLower(), n, 0, alpha, a, 0, x, incx, beta, y, incy);
}

/// `?HPMV` (`interface/zhpmv.c`): `hbmv` for packed storage.
pub fn hpmv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, alpha: T, a: [*]T, x: [*]T, incx: isize, beta: T, y: [*]T, incy: isize) Error!void {
    if (comptime !isComplex(T)) @compileError("hpmv: " ++ @typeName(T) ++ " is not complex");
    if (n < 0) return error.InvalidN;
    if (incx == 0) return error.InvalidIncx;
    if (incy == 0) return error.InvalidIncy;
    switch (uplo.isRev()) {
        inline else => |rev| try mvInterface(allocator, io, T, if (rev) .hermitian_rev else .hermitian, true, uplo.isLower(), n, 0, alpha, a, 0, x, incx, beta, y, incy),
    }
}

/// OpenBLAS's small-`n` paths (`interface/syr.c`, `spr.c`, `syr2.c`, `spr2.c`,
/// `zsyr.c`): below this `n`, a unit-stride update runs its columns in the
/// interface, with no buffer; 0 where upstream has none.
fn smallN(comptime T: type, comptime kind: Kind, comptime packed_: bool, comptime two: bool) BlasLong {
    if (kind != .symmetric) return 0;
    if (isComplex(T)) return if (!two and !packed_) 50 else 0;
    return if (two and packed_) 50 else 100;
}

fn rankInterface(allocator: Allocator, io: Io, comptime T: type, comptime kind: Kind, comptime packed_: bool, lower: bool, n: isize, alpha: if (kind == .symmetric) T else Float(T), x: [*]T, incx: isize, a: [*]T, lda: isize) Error!void {
    _ = io;
    if (n == 0) return;
    try spans(T, packed_, n, lda);
    if (!util.fits(n, incx, @sizeOf(T))) return error.InvalidIncx;
    if (comptime kind == .symmetric) {
        if (isZero(T, alpha)) return;
    } else if (alpha == 0) return;
    _ = blas.vtable();
    if (kind != .hermitian_rev and accel.syr(T, kind == .hermitian, packed_, if (lower) .l else .u, n, alpha, x, incx, a, lda)) return;
    const small = incx == 1 and n < comptime smallN(T, kind, packed_, false);
    const px = first(T, x, n, incx);
    switch (lower) {
        inline else => |lo| if (small)
            rankColumns(T, kind, lo, packed_, n, alpha, x, a, lda)
        else
            try rankDriver(T, kind, lo, packed_, allocator, n, alpha, px, incx, a, lda),
    }
}

fn rank2Interface(allocator: Allocator, io: Io, comptime T: type, comptime kind: Kind, comptime packed_: bool, lower: bool, n: isize, alpha: T, x: [*]T, incx: isize, y: [*]T, incy: isize, a: [*]T, lda: isize) Error!void {
    _ = io;
    if (n == 0 or isZero(T, alpha)) return;
    try spans(T, packed_, n, lda);
    if (!util.fits(n, incx, @sizeOf(T))) return error.InvalidIncx;
    if (!util.fits(n, incy, @sizeOf(T))) return error.InvalidIncy;
    _ = blas.vtable();
    if (kind != .hermitian_rev and accel.syr2(T, kind == .hermitian, packed_, if (lower) .l else .u, n, alpha, x, incx, y, incy, a, lda)) return;
    const small = incx == 1 and incy == 1 and n < comptime smallN(T, kind, packed_, true);
    const px = first(T, x, n, incx);
    const py = first(T, y, n, incy);
    switch (lower) {
        inline else => |lo| if (small)
            rank2Columns(T, kind, lo, packed_, n, alpha, x, y, a, lda)
        else
            try rank2Driver(T, kind, lo, packed_, allocator, n, alpha, px, incx, py, incy, a, lda),
    }
}

/// `?SPR` (`interface/spr.c`, `zspr.c`): `A += alpha*x*x^T` on the packed
/// triangle. Columns where `x[i]` is zero are skipped.
pub fn spr(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, alpha: T, x: [*]T, incx: isize, a: [*]T) Error!void {
    if (uplo.isRev()) return error.InvalidUplo;
    if (n < 0) return error.InvalidN;
    if (incx == 0) return error.InvalidIncx;
    try rankInterface(allocator, io, T, .symmetric, true, uplo.isLower(), n, alpha, x, incx, a, 0);
}

/// `?SYR` (`interface/syr.c`, `zsyr.c`): `spr` for full storage.
pub fn syr(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, alpha: T, x: [*]T, incx: isize, a: [*]T, lda: isize) Error!void {
    if (uplo.isRev()) return error.InvalidUplo;
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, n)) return error.InvalidLda;
    if (incx == 0) return error.InvalidIncx;
    try rankInterface(allocator, io, T, .symmetric, false, uplo.isLower(), n, alpha, x, incx, a, lda);
}

/// `?HPR` (`interface/zhpr.c`): `A += alpha*x*x^H` for real `alpha`; the
/// diagonal's imaginary part becomes zero.
pub fn hpr(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, alpha: Float(T), x: [*]T, incx: isize, a: [*]T) Error!void {
    if (comptime !isComplex(T)) @compileError("hpr: " ++ @typeName(T) ++ " is not complex");
    if (n < 0) return error.InvalidN;
    if (incx == 0) return error.InvalidIncx;
    switch (uplo.isRev()) {
        inline else => |rev| try rankInterface(allocator, io, T, if (rev) .hermitian_rev else .hermitian, true, uplo.isLower(), n, alpha, x, incx, a, 0),
    }
}

/// `?HER` (`interface/zher.c`): `hpr` for full storage.
pub fn her(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, alpha: Float(T), x: [*]T, incx: isize, a: [*]T, lda: isize) Error!void {
    if (comptime !isComplex(T)) @compileError("her: " ++ @typeName(T) ++ " is not complex");
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, n)) return error.InvalidLda;
    if (incx == 0) return error.InvalidIncx;
    switch (uplo.isRev()) {
        inline else => |rev| try rankInterface(allocator, io, T, if (rev) .hermitian_rev else .hermitian, false, uplo.isLower(), n, alpha, x, incx, a, lda),
    }
}

/// `?SPR2` (`interface/spr2.c`, `zspr2.c`): `A += alpha*(x*y^T + y*x^T)` on
/// the packed triangle.
pub fn spr2(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, alpha: T, x: [*]T, incx: isize, y: [*]T, incy: isize, a: [*]T) Error!void {
    if (uplo.isRev()) return error.InvalidUplo;
    if (n < 0) return error.InvalidN;
    if (incx == 0) return error.InvalidIncx;
    if (incy == 0) return error.InvalidIncy;
    try rank2Interface(allocator, io, T, .symmetric, true, uplo.isLower(), n, alpha, x, incx, y, incy, a, 0);
}

/// `?SYR2` (`interface/syr2.c`, `zsyr2.c`): `spr2` for full storage.
pub fn syr2(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, alpha: T, x: [*]T, incx: isize, y: [*]T, incy: isize, a: [*]T, lda: isize) Error!void {
    if (uplo.isRev()) return error.InvalidUplo;
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, n)) return error.InvalidLda;
    if (incx == 0) return error.InvalidIncx;
    if (incy == 0) return error.InvalidIncy;
    try rank2Interface(allocator, io, T, .symmetric, false, uplo.isLower(), n, alpha, x, incx, y, incy, a, lda);
}

/// `?HPR2` (`interface/zhpr2.c`): `A += alpha*x*y^H + conj(alpha)*y*x^H`; the
/// diagonal's imaginary part becomes zero.
pub fn hpr2(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, alpha: T, x: [*]T, incx: isize, y: [*]T, incy: isize, a: [*]T) Error!void {
    if (comptime !isComplex(T)) @compileError("hpr2: " ++ @typeName(T) ++ " is not complex");
    if (n < 0) return error.InvalidN;
    if (incx == 0) return error.InvalidIncx;
    if (incy == 0) return error.InvalidIncy;
    switch (uplo.isRev()) {
        inline else => |rev| try rank2Interface(allocator, io, T, if (rev) .hermitian_rev else .hermitian, true, uplo.isLower(), n, alpha, x, incx, y, incy, a, 0),
    }
}

/// `?HER2` (`interface/zher2.c`): `hpr2` for full storage.
pub fn her2(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, alpha: T, x: [*]T, incx: isize, y: [*]T, incy: isize, a: [*]T, lda: isize) Error!void {
    if (comptime !isComplex(T)) @compileError("her2: " ++ @typeName(T) ++ " is not complex");
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, n)) return error.InvalidLda;
    if (incx == 0) return error.InvalidIncx;
    if (incy == 0) return error.InvalidIncy;
    switch (uplo.isRev()) {
        inline else => |rev| try rank2Interface(allocator, io, T, if (rev) .hermitian_rev else .hermitian, false, uplo.isLower(), n, alpha, x, incx, y, incy, a, lda),
    }
}

// ---- tests -------------------------------------------------------------------

const Vtable = @import("../Vtable.zig");
const check = @import("../_test/blas1/check.zig");
const Complex32 = types.Complex32;
const Complex64 = types.Complex64;

const all = .{ f32, f64, Complex32, Complex64 } ++ check.x_real ++ check.x_complex;
const talloc = std.testing.allocator;
const max_n = 9;
const max_mem = 256;

/// Where logical element `i` of a vector of `n` with stride `inc` lives.
fn vecAt(i: usize, n: usize, inc: isize) usize {
    const s: usize = @abs(inc);
    return if (inc > 0) i * s else (n - 1 - i) * s;
}

const Store = enum { full, packed_, band };

/// Memory index of the stored entry `(i, j)`, or null when it is not in the
/// stored triangle or band.
fn storeIdx(comptime st: Store, lower: bool, n: usize, ld: usize, k: usize, i: usize, j: usize) ?usize {
    if (lower) {
        if (i < j) return null;
    } else if (i > j) return null;
    const d = if (lower) i - j else j - i;
    return switch (st) {
        .full => i + j * ld,
        .packed_ => if (lower) j * (2 * n - j + 1) / 2 + d else j * (j + 1) / 2 + i,
        .band => if (d > k) null else (if (lower) d else k - d) + j * ld,
    };
}

/// Entry `(i, j)` of the matrix a stored triangle stands for.
fn dense(comptime T: type, comptime st: Store, comptime kind: Kind, lower: bool, n: usize, ld: usize, k: usize, a: []const T, i: usize, j: usize) T {
    if (i == j) {
        const d = a[storeIdx(st, lower, n, ld, k, i, i).?];
        return if (kind == .symmetric) d else if (comptime isComplex(T)) T.init(d.re, 0) else d;
    }
    if (storeIdx(st, lower, n, ld, k, i, j)) |p| {
        return if (kind == .hermitian_rev) conj(T, a[p]) else a[p];
    }
    if (storeIdx(st, lower, n, ld, k, j, i)) |p| {
        return if (kind == .hermitian) conj(T, a[p]) else a[p];
    }
    return zero(T);
}

fn tol(comptime T: type, terms: usize) Float(T) {
    return check.tolerance(Float(T), 8 * terms);
}

const cases = [_]struct { n: usize, k: usize }{ .{ .n = 1, .k = 0 }, .{ .n = 2, .k = 1 }, .{ .n = 3, .k = 0 }, .{ .n = 5, .k = 2 }, .{ .n = 9, .k = 1 }, .{ .n = 9, .k = 4 }, .{ .n = 9, .k = 12 } };
const inc_pairs = [_][2]isize{ .{ 1, 1 }, .{ 2, 3 }, .{ -1, 1 }, .{ 1, -2 }, .{ -3, -2 } };

fn scalarOf(comptime T: type, re: f64, im: f64) T {
    return if (comptime isComplex(T)) T.init(@floatCast(re), @floatCast(im)) else @floatCast(re);
}

/// Checks `y := alpha*S*x + beta*y` against dense arithmetic.
fn checkMv(comptime T: type, comptime st: Store, comptime kind: Kind, lower: bool) !void {
    var seed: u64 = 1;
    for (cases) |c| for (inc_pairs) |inc| for ([_]f64{ 1, 0, -0.5 }) |bre| {
        seed += 1;
        const n = c.n;
        const k = if (st == .band) c.k else 0;
        const ld: usize = k + 1 + @as(usize, @intCast(seed & 1));
        var a: [max_mem]T = undefined;
        var xs: [max_mem]T = undefined;
        var ys: [max_mem]T = undefined;
        check.fill(T, &a, seed);
        check.fill(T, &xs, seed + 100);
        check.fill(T, &ys, seed + 200);
        const alpha = scalarOf(T, 1.25, -0.5);
        const beta = scalarOf(T, bre, if (bre == -0.5) 0.75 else 0);
        var want = ys;
        for (0..n) |i| {
            var acc = zero(T);
            for (0..n) |j| acc = add(T, acc, mul(T, dense(T, st, kind, lower, n, ld, k, &a, i, j), xs[vecAt(j, n, inc[0])]));
            const yp = vecAt(i, n, inc[1]);
            want[yp] = add(T, if (bre == 0) zero(T) else mul(T, beta, ys[yp]), mul(T, alpha, acc));
        }
        if (bre == 0) for (0..n) |i| {
            ys[vecAt(i, n, inc[1])] = scalarOf(T, std.math.nan(f64), std.math.nan(f64));
        };
        const un: Uplo = if (lower) .l else .u;
        const hu: Uplo = switch (kind) {
            .hermitian_rev => if (lower) .m else .v,
            else => un,
        };
        const nn: isize = @intCast(n);
        const lda: isize = @intCast(ld);
        switch (st) {
            .band => if (kind == .symmetric)
                try sbmv(talloc, check.io, T, un, nn, @intCast(k), alpha, &a, lda, &xs, inc[0], beta, &ys, inc[1])
            else
                try hbmv(talloc, check.io, T, hu, nn, @intCast(k), alpha, &a, lda, &xs, inc[0], beta, &ys, inc[1]),
            .packed_ => if (kind == .symmetric)
                try spmv(talloc, check.io, T, un, nn, alpha, &a, &xs, inc[0], beta, &ys, inc[1])
            else
                try hpmv(talloc, check.io, T, hu, nn, alpha, &a, &xs, inc[0], beta, &ys, inc[1]),
            .full => unreachable,
        }
        check.expectCloseSlices(T, &want, &ys, tol(T, n)) catch |e| {
            check.say("mv {s} {s} n={d} k={d} inc={d},{d} lower={}: ", .{ @tagName(st), @tagName(kind), n, k, inc[0], inc[1], lower });
            return e;
        };
    };
}

fn kindsOf(comptime T: type) []const Kind {
    return if (comptime isComplex(T)) &.{ .symmetric, .hermitian, .hermitian_rev } else &.{.symmetric};
}

/// Checks the rank-1 and rank-2 updates against dense arithmetic.
fn checkRank(comptime T: type, comptime st: Store, comptime kind: Kind, lower: bool) !void {
    var seed: u64 = 1000;
    for (cases) |c| for (inc_pairs) |inc| for ([_]bool{ false, true }) |two| {
        seed += 1;
        const n = c.n;
        const ld: usize = n + @as(usize, @intCast(seed & 1));
        var a: [max_mem]T = undefined;
        var xs: [max_mem]T = undefined;
        var ys: [max_mem]T = undefined;
        check.fill(T, &a, seed);
        check.fill(T, &xs, seed + 100);
        check.fill(T, &ys, seed + 200);
        if (c.n > 2) xs[vecAt(1, n, inc[0])] = zero(T); // the zero-skip path
        const alpha = scalarOf(T, 1.25, -0.5);
        const alpha_f: Float(T) = 1.25;
        const a0 = alpha;
        var want = a;
        for (0..n) |c_| for (0..n) |r| {
            const p = storeIdx(st, lower, n, ld, 0, r, c_) orelse continue;
            const xr = xs[vecAt(r, n, inc[0])];
            const xc = xs[vecAt(c_, n, inc[0])];
            const yr = ys[vecAt(r, n, inc[1])];
            const yc = ys[vecAt(c_, n, inc[1])];
            var v: T = undefined;
            switch (kind) {
                .symmetric => v = if (two)
                    mul(T, alpha, add(T, mul(T, xr, yc), mul(T, yr, xc)))
                else
                    mul(T, alpha, mul(T, xr, xc)),
                .hermitian, .hermitian_rev => {
                    v = if (two)
                        add(T, mul(T, alpha, mul(T, xr, conj(T, yc))), mul(T, conj(T, alpha), mul(T, yr, conj(T, xc))))
                    else
                        scaleF(T, mul(T, xr, conj(T, xc)), alpha_f);
                    if (kind == .hermitian_rev) v = conj(T, v);
                },
            }
            want[p] = add(T, want[p], v);
            if (r == c_ and kind != .symmetric) {
                if (comptime isComplex(T)) want[p].im = 0;
            }
        };
        const un: Uplo = if (lower) .l else .u;
        const hu: Uplo = switch (kind) {
            .hermitian_rev => if (lower) .m else .v,
            else => un,
        };
        const nn: isize = @intCast(n);
        const lda: isize = @intCast(ld);
        switch (st) {
            .packed_ => if (two) {
                if (kind == .symmetric)
                    try spr2(talloc, check.io, T, un, nn, a0, &xs, inc[0], &ys, inc[1], &a)
                else
                    try hpr2(talloc, check.io, T, hu, nn, a0, &xs, inc[0], &ys, inc[1], &a);
            } else {
                if (kind == .symmetric)
                    try spr(talloc, check.io, T, un, nn, a0, &xs, inc[0], &a)
                else
                    try hpr(talloc, check.io, T, hu, nn, alpha_f, &xs, inc[0], &a);
            },
            .full => if (two) {
                if (kind == .symmetric)
                    try syr2(talloc, check.io, T, un, nn, a0, &xs, inc[0], &ys, inc[1], &a, lda)
                else
                    try her2(talloc, check.io, T, hu, nn, a0, &xs, inc[0], &ys, inc[1], &a, lda);
            } else {
                if (kind == .symmetric)
                    try syr(talloc, check.io, T, un, nn, a0, &xs, inc[0], &a, lda)
                else
                    try her(talloc, check.io, T, hu, nn, alpha_f, &xs, inc[0], &a, lda);
            },
            .band => unreachable,
        }
        check.expectCloseSlices(T, &want, &a, tol(T, 4)) catch |e| {
            check.say("rank{s} {s} {s} n={d} inc={d},{d} lower={}: ", .{ if (two) "2" else "1", @tagName(st), @tagName(kind), n, inc[0], inc[1], lower });
            return e;
        };
    };
}

/// `gbmv` against dense arithmetic over every `trans`.
fn checkGbmv(comptime T: type) !void {
    const shapes = [_][4]usize{ .{ 1, 1, 0, 0 }, .{ 4, 4, 1, 1 }, .{ 5, 3, 2, 0 }, .{ 3, 6, 0, 2 }, .{ 7, 5, 3, 2 }, .{ 2, 9, 4, 4 } };
    var seed: u64 = 5000;
    inline for (std.enums.values(Trans)) |tr| {
        if (comptime !isComplex(T) and (tr.conjA() or tr.conjX())) continue;
        for (shapes) |s| for (inc_pairs) |inc| for ([_]f64{ 1, 0, 0.5 }) |bre| {
            seed += 1;
            const m, const n, const kl, const ku = s;
            const ld: usize = kl + ku + 1 + @as(usize, @intCast(seed & 1));
            const lenx = if (tr.transposed()) m else n;
            const leny = if (tr.transposed()) n else m;
            var a: [max_mem]T = undefined;
            var xs: [max_mem]T = undefined;
            var ys: [max_mem]T = undefined;
            check.fill(T, &a, seed);
            check.fill(T, &xs, seed + 100);
            check.fill(T, &ys, seed + 200);
            const alpha = scalarOf(T, 0.75, 1.5);
            const beta = scalarOf(T, bre, if (bre == 0.5) -1 else 0);
            var want = ys;
            for (0..leny) |o| {
                var acc = zero(T);
                for (0..lenx) |p| {
                    // A(i, j) with i the row, j the column of the m by n matrix
                    const i = if (tr.transposed()) p else o;
                    const j = if (tr.transposed()) o else p;
                    if (i + ku < j or i > j + kl) continue;
                    var e = a[ku + i - j + j * ld];
                    if (tr.conjA()) e = conj(T, e);
                    var xv = xs[vecAt(p, lenx, inc[0])];
                    if (tr.conjX()) xv = conj(T, xv);
                    acc = add(T, acc, mul(T, e, xv));
                }
                const yp = vecAt(o, leny, inc[1]);
                want[yp] = add(T, if (bre == 0) zero(T) else mul(T, beta, ys[yp]), mul(T, alpha, acc));
            }
            if (bre == 0) for (0..leny) |o| {
                ys[vecAt(o, leny, inc[1])] = scalarOf(T, std.math.nan(f64), std.math.nan(f64));
            };
            try gbmv(talloc, check.io, T, tr, @intCast(m), @intCast(n), @intCast(kl), @intCast(ku), alpha, &a, @intCast(ld), &xs, inc[0], beta, &ys, inc[1]);
            check.expectCloseSlices(T, &want, &ys, tol(T, lenx)) catch |e| {
                check.say("gbmv {s} {s} m={d} n={d} kl={d} ku={d} inc={d},{d}: ", .{ @typeName(T), @tagName(tr), m, n, kl, ku, inc[0], inc[1] });
                return e;
            };
        };
    }
}

fn runAll(comptime list: anytype) !void {
    inline for (list) |T| {
        try checkGbmv(T);
        inline for (.{ Store.band, Store.packed_ }) |st| {
            inline for (comptime kindsOf(T)) |kind| for ([_]bool{ false, true }) |lower| try checkMv(T, st, kind, lower);
        }
        inline for (.{ Store.full, Store.packed_ }) |st| {
            inline for (comptime kindsOf(T)) |kind| for ([_]bool{ false, true }) |lower| try checkRank(T, st, kind, lower);
        }
    }
}

/// Runs the cases with `vt`'s kernels in the global vtable.
fn withCore(vt: *const Vtable) anyerror!void {
    const saved = blas.testUse(vt);
    defer blas.testRestore(saved);
    try runAll(all);
}

test "banded, packed and rank updates match dense arithmetic (portable kernels)" {
    try runAll(all);
}

test "f80 without slots: banded, packed and rank updates match dense arithmetic (portable kernels)" {
    if (comptime check.f80_unslotted.len == 0) return error.SkipZigTest;
    try runAll(check.f80_unslotted);
}

test "banded, packed and rank updates match dense arithmetic on every core" {
    try check.eachCore(withCore);
}

/// The small-`n` path against the driver, bit for bit, for every `n` below
/// the threshold, with zeros in `x` and `y` for the skip.
fn checkSmall(comptime T: type, comptime packed_: bool, comptime two: bool) !void {
    const limit = comptime smallN(T, .symmetric, packed_, two);
    if (limit == 0) return;
    const len: usize = @intCast(limit * limit);
    const a = try talloc.alloc(T, len);
    defer talloc.free(a);
    const b = try talloc.alloc(T, len);
    defer talloc.free(b);
    var xs: [limit]T = undefined;
    var ys: [limit]T = undefined;
    check.fill(T, &xs, 7);
    check.fill(T, &ys, 8);
    for (0..limit) |i| {
        if (i % 3 == 1) xs[i] = zero(T);
        if (i % 6 == 1 or i % 5 == 3) ys[i] = zero(T);
    }
    const alpha = scalarOf(T, 1.25, -0.5);
    for (1..limit) |n| for ([_]bool{ false, true }) |lower| {
        const nn: isize = @intCast(n);
        const lda: isize = if (packed_) 0 else nn + @as(isize, @intCast(n & 1));
        check.fill(T, a, n);
        @memcpy(b, a);
        const u: Uplo = if (lower) .l else .u;
        if (two) {
            if (packed_) try spr2(talloc, check.io, T, u, nn, alpha, &xs, 1, &ys, 1, a.ptr) else try syr2(talloc, check.io, T, u, nn, alpha, &xs, 1, &ys, 1, a.ptr, lda);
        } else {
            if (packed_) try spr(talloc, check.io, T, u, nn, alpha, &xs, 1, a.ptr) else try syr(talloc, check.io, T, u, nn, alpha, &xs, 1, a.ptr, lda);
        }
        switch (lower) {
            inline else => |lo| if (two)
                try rank2Driver(T, .symmetric, lo, packed_, talloc, nn, alpha, &xs, 1, &ys, 1, b.ptr, lda)
            else
                try rankDriver(T, .symmetric, lo, packed_, talloc, nn, alpha, &xs, 1, b.ptr, lda),
        }
        for (a, b, 0..) |p, q, i| if (!check.same(T, p, q, true)) {
            check.say("small rank{s} {s} {s} n={d} lower={} at {d}: ", .{ if (two) "2" else "1", @typeName(T), if (packed_) "packed" else "full", n, lower, i });
            return error.TestExpectedEqual;
        };
    };
}

fn runSmall() !void {
    inline for (.{ f16, f32, f64, f80, f128 }) |F| {
        inline for (.{ F, std.math.Complex(F) }) |T| {
            inline for (.{ false, true }) |packed_| inline for (.{ false, true }) |two| try checkSmall(T, packed_, two);
        }
    }
}

fn withCoreSmall(vt: *const Vtable) anyerror!void {
    const saved = blas.testUse(vt);
    defer blas.testRestore(saved);
    try runSmall();
}

test "the small-n rank updates match the driver bit for bit" {
    // rstd's own paths
    const was = blas.useAccelerate(false);
    defer _ = blas.useAccelerate(was);
    try runSmall();
    try check.eachCore(withCoreSmall);
}

test "strided vectors beyond the stack buffer are copied from the allocator" {
    const n = 1200; // > stack_bytes / 8
    const k = 3;
    const ld = k + 1;
    const a = try talloc.alloc(f64, ld * n);
    defer talloc.free(a);
    const xs = try talloc.alloc(f64, 2 * n);
    defer talloc.free(xs);
    const ys = try talloc.alloc(f64, 3 * n);
    defer talloc.free(ys);
    check.fill(f64, a, 1);
    check.fill(f64, xs, 2);
    check.fill(f64, ys, 3);
    const want = try talloc.dupe(f64, ys);
    defer talloc.free(want);
    for (0..n) |i| {
        var acc: f64 = 0;
        for (0..n) |j| {
            const d = if (j > i) j - i else i - j;
            if (d > k) continue;
            acc += a[k - d + @max(i, j) * ld] * xs[vecAt(j, n, -2)];
        }
        want[vecAt(i, n, 3)] += 2 * acc;
    }
    try sbmv(talloc, check.io, f64, .u, n, k, 2, a.ptr, ld, xs.ptr, -2, 1, ys.ptr, 3);
    try check.expectCloseSlices(f64, want, ys, tol(f64, k));
}

test "an allocation failure comes back as OutOfMemory" {
    // rstd's own path, which allocates
    const was = blas.useAccelerate(false);
    defer _ = blas.useAccelerate(was);
    var failing = std.testing.FailingAllocator.init(talloc, .{ .fail_index = 0 });
    const n = 1200;
    const a = try talloc.alloc(f32, n);
    defer talloc.free(a);
    const xs = try talloc.alloc(f32, 2 * n);
    defer talloc.free(xs);
    @memset(a, 1);
    @memset(xs, 1);
    try std.testing.expectError(error.OutOfMemory, syr(failing.allocator(), check.io, f32, .l, n, 1, xs.ptr, 2, a.ptr, n));
}

test "invalid arguments and quick returns" {
    var a = [_]f64{ 1, 2, 3, 4 };
    var x = [_]f64{ 1, 2 };
    var y = [_]f64{ 5, 6 };
    try std.testing.expectError(error.InvalidLda, sbmv(talloc, check.io, f64, .u, 2, 1, 1, &a, 1, &x, 1, 1, &y, 1));
    try std.testing.expectError(error.InvalidIncx, sbmv(talloc, check.io, f64, .u, 2, 0, 1, &a, 1, &x, 0, 1, &y, 1));
    try std.testing.expectError(error.InvalidN, spmv(talloc, check.io, f64, .u, -1, 1, &a, &x, 1, 1, &y, 1));
    try std.testing.expectError(error.InvalidLda, syr(talloc, check.io, f64, .l, 2, 1, &x, 1, &a, 1));
    try std.testing.expectError(error.InvalidLda, gbmv(talloc, check.io, f64, .n, 2, 2, 1, 1, 1, &a, 2, &x, 1, 1, &y, 1));
    // n == 0 and alpha == 0 leave A alone; beta still scales y
    try syr(talloc, check.io, f64, .u, 0, 1, &x, 1, &a, 1);
    try syr(talloc, check.io, f64, .u, 2, 0, &x, 1, &a, 2);
    try std.testing.expectEqualSlices(f64, &.{ 1, 2, 3, 4 }, &a);
    try spmv(talloc, check.io, f64, .u, 2, 0, &a, &x, 1, 2, &y, 1);
    try std.testing.expectEqualSlices(f64, &.{ 10, 12 }, &y);
    // the symmetric routines take no Hermitian `uplo`, real `gbmv` no `o` to `d`
    try std.testing.expectError(error.InvalidUplo, sbmv(talloc, check.io, f64, .v, -1, -1, 1, &a, 0, &x, 0, 1, &y, 0));
    try std.testing.expectError(error.InvalidUplo, syr2(talloc, check.io, f64, .m, -1, 1, &x, 0, &y, 0, &a, 0));
    try std.testing.expectError(error.InvalidTrans, gbmv(talloc, check.io, f64, .d, -1, 2, 1, 1, 1, &a, 2, &x, 1, 1, &y, 1));
    // real `r` and `c` are `n` and `t`
    var yn = [_]f64{ 0, 0 };
    var yr = yn;
    try gbmv(talloc, check.io, f64, .t, 2, 2, 1, 0, 1, &a, 2, &x, 1, 0, &yn, 1);
    try gbmv(talloc, check.io, f64, .c, 2, 2, 1, 0, 1, &a, 2, &x, 1, 0, &yr, 1);
    try std.testing.expectEqualSlices(f64, &yn, &yr);
}
