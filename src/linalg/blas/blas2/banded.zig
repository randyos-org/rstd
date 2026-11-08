//! `?GBMV`, `?SBMV`, `?HBMV`, `?SPMV`, `?HPMV`: `interface/*.c` over
//! `driver/level2/*_k.c`, a column at a time through the level-1 kernels.
//! `?HBMV` takes `HEMVREV`'s `V` and `M`, as the Fortran interface does.
//! Transposed `gbmv` threads by columns; the others through partial vectors
//! folded in chunk order.

const std = @import("std");
const blas = @import("../__root__.zig");
const types = @import("../types.zig");
const util = @import("../util.zig");
const interface = @import("interface.zig");
const portable = @import("../portable.zig");
const accel = @import("../apple_accel.zig");
const threading = @import("../threading.zig");

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

const raw = util.raw;
const off = util.off;
const at = util.at;
const first = util.first;
const zero = util.zero;
const isZero = util.isZero;
const isOne = util.isOne;

fn mul(comptime T: type, a: T, b: T) T {
    return if (comptime isComplex(T)) a.mul(b) else a * b;
}

/// `mul`, each step rounded on the x87 (`util.cmul`): the rank updates'
/// scalars.
pub fn mulRounded(comptime T: type, a: T, b: T) T {
    if (comptime !isComplex(T)) return a * b;
    const p = util.cmul(Float(T), a.re, a.im, b.re, b.im);
    return .init(p[0], p[1]);
}

fn add(comptime T: type, a: T, b: T) T {
    return if (comptime isComplex(T)) a.add(b) else a + b;
}

pub fn conj(comptime T: type, a: T) T {
    return if (comptime isComplex(T)) a.conjugate() else a;
}

pub fn scaleF(comptime T: type, a: T, f: Float(T)) T {
    return if (comptime isComplex(T)) T.init(a.re * f, a.im * f) else a * f;
}

// ---- level-1 kernels, unit stride --------------------------------------------

/// `y += alpha*x`, or `alpha*conj(x)`; the kernels take `n >= 1`.
pub fn axpyK(comptime T: type, comptime cj: bool, n: BlasLong, alpha: T, x: [*]T, y: [*]T) void {
    if (n <= 0) return;
    const op: portable.Op = comptime if (cj and isComplex(T)) .axpyc_k else .axpy_k;
    _ = @call(.auto, portable.fetch(T, op, blas.table), .{ n, 0, 0 } ++ util.scalarArgs(T, alpha) ++ .{ raw(T, x), 1, raw(T, y), 1, raw(T, y), 0 });
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
    _ = @call(.auto, portable.fetch(T, .scal_k, blas.table), .{ n, 0, 0 } ++ util.scalarArgs(T, beta) ++ .{ p, incy, p, 0, p, 0 });
}

// ---- unit-stride copies ------------------------------------------------------

/// Room for one vector copy: on the stack when it fits, else from the
/// allocator.
pub fn Scratch(comptime T: type) type {
    return struct {
        stack: [stack_bytes]u8 align(@alignOf(T)) = undefined,
        buf: ?util.Scratch(T, @alignOf(T), @alignOf(T)) = null,

        /// `v` (length `n`, stride `inc`) as a unit-stride vector: `v` itself
        /// at stride 1, else a copy.
        pub fn unit(self: *@This(), allocator: Allocator, n: BlasLong, v: [*]T, inc: BlasLong) Allocator.Error![*]T {
            if (inc == 1) return v;
            const buf = try util.Scratch(T, @alignOf(T), @alignOf(T)).init(allocator, &self.stack, @intCast(n));
            self.buf = buf;
            copyK(T, n, v, inc, buf.ptr, 1);
            return buf.ptr;
        }

        pub fn deinit(self: *@This(), allocator: Allocator) void {
            if (self.buf) |b| b.deinit(allocator);
        }
    };
}

// ---- drivers -----------------------------------------------------------------

/// `driver/level2/gbmv_k.c`, `zgbmv_k.c`: `y += alpha*op(A)*x` for an `m` by
/// `n` band matrix. The flags are `TRANS`, `CONJ` (of `A`) and `XCONJ` (of `x`).
/// Above the threshold, as `gbmv_thread.c`, the columns split: with a
/// transpose each element of `y` is one column's dot; without, each chunk's
/// columns go into a partial vector of their own, folded into `y` in chunk
/// order (`threading.Partials`).
fn gbmvDriver(comptime T: type, comptime tr: bool, comptime cj: bool, comptime xc: bool, allocator: Allocator, io: Io, m: BlasLong, n: BlasLong, ku: BlasLong, kl: BlasLong, alpha: T, a: [*]T, lda: BlasLong, x: [*]T, incx: BlasLong, y: [*]T, incy: BlasLong) (Allocator.Error || Io.Cancelable)!void {
    const leny = if (tr) n else m;
    const lenx = if (tr) m else n;
    var sy: Scratch(T) = .{};
    var sx: Scratch(T) = .{};
    defer sy.deinit(allocator);
    defer sx.deinit(allocator);
    var g: Gbmv(T, tr, cj, xc) = .{
        .m = m,
        .ku = ku,
        .kl = kl,
        .alpha = alpha,
        .a = a,
        .lda = lda,
        .X = try sx.unit(allocator, lenx, x, incx),
        .Y = undefined,
    };
    const cols = @min(n, m + ku);
    // `gbmv.c`, `zgbmv.c`: `m*n` under 250000 (125000 complex) or `kl + ku`
    // under 15 runs serially
    const p = threading.plan(.gbmv, T, threading.product(&.{ @intCast(m), @intCast(n) }), @intCast(cols), .even, 8);
    const split = kl + ku >= 15 and p.chunks.count > 1;
    if (!tr and split) {
        const parts: threading.Partials(T) = try .init(allocator, &p, &g, @TypeOf(g).rows);
        defer parts.deinit(allocator);
        const C = struct {
            g: *const @TypeOf(g),
            parts: *const threading.Partials(T),

            fn one(c: @This(), ch: threading.Chunk) void {
                var h = c.g.*;
                h.Y = c.parts.take(ch.index).ptr;
                h.yb = @intCast(c.parts.rows[ch.index].lo);
                h.columns(@intCast(ch.start), @intCast(ch.end));
            }
        };
        try threading.forEach(io, &p, C{ .g = &g, .parts = &parts }, C.one);
        return parts.fold(io, y, incy);
    }
    g.Y = try sy.unit(allocator, leny, y, incy);
    if (tr and split and p.workers > 1) {
        try threading.forEach(io, &p, &g, Gbmv(T, tr, cj, xc).chunk);
    } else g.columns(0, cols);
    if (incy != 1) copyK(T, leny, g.Y, 1, y, incy);
}

/// `gbmvDriver`'s columns, over unit-stride `X` and `Y`; element 0 of `Y` is
/// row `yb`.
fn Gbmv(comptime T: type, comptime tr: bool, comptime cj: bool, comptime xc: bool) type {
    return struct {
        m: BlasLong,
        ku: BlasLong,
        kl: BlasLong,
        alpha: T,
        a: [*]T,
        lda: BlasLong,
        X: [*]T,
        Y: [*]T,
        yb: BlasLong = 0,

        /// The rows of `y` columns `c.start..c.end` reach without a transpose.
        fn rows(g: *const @This(), c: threading.Chunk) threading.Rows {
            const mu: usize = @intCast(g.m);
            return .{ .lo = c.start -| @as(usize, @intCast(g.ku)), .hi = @min(mu, c.end + @as(usize, @intCast(g.kl))) };
        }

        fn chunk(g: *const @This(), c: threading.Chunk) void {
            g.columns(@intCast(c.start), @intCast(c.end));
        }

        fn columns(g: *const @This(), c0: BlasLong, c1: BlasLong) void {
            var offset_u = g.ku - c0;
            var offset_l = g.ku + g.m - c0;
            var col = off(T, g.a, c0 * g.lda);
            var i: BlasLong = c0;
            while (i < c1) : (i += 1) {
                const start = @max(offset_u, 0);
                const end = @min(offset_l, g.ku + g.kl + 1);
                const length = end - start;
                if (length > 0) {
                    const ap = off(T, col, start);
                    const xs = off(T, g.X, start - offset_u);
                    if (comptime !tr) {
                        const xi = at(T, g.X, i).*;
                        axpyK(T, cj, length, mul(T, g.alpha, if (xc) conj(T, xi) else xi), ap, off(T, g.Y, start - offset_u - g.yb));
                    } else {
                        const t = if (comptime !xc) dotK(T, cj, length, ap, xs) else dotK(T, !cj, length, xs, ap);
                        const v = mul(T, g.alpha, if (comptime xc and cj) conj(T, t) else t);
                        at(T, g.Y, i).* = add(T, at(T, g.Y, i).*, v);
                    }
                }
                offset_u -= 1;
                offset_l -= 1;
                col = off(T, col, g.lda);
            }
        }
    };
}

/// Which matrix a stored triangle stands for.
const Kind = @import("symmetric.zig").Kind;

/// `driver/level2/{sbmv,spmv,zsbmv,zspmv,zhbmv,zhpmv}_k.c`: `y += alpha*S*x`
/// for `S` in band storage (`k` bands, `lda`) or packed storage. Over more
/// than one chunk, as `sbmv_thread.c` and `spmv_thread.c` split the columns,
/// each chunk's columns go into a partial vector of their own, folded into `y`
/// in chunk order (`threading.Partials`).
fn mvDriver(comptime T: type, comptime kind: Kind, comptime lower: bool, comptime packed_: bool, allocator: Allocator, io: Io, n: BlasLong, k: BlasLong, alpha: T, a: [*]T, lda: BlasLong, x: [*]T, incx: BlasLong, y: [*]T, incy: BlasLong) (Allocator.Error || Io.Cancelable)!void {
    var sx: Scratch(T) = .{};
    defer sx.deinit(allocator);
    const X = try sx.unit(allocator, n, x, incx);
    const nu: usize = @intCast(n);
    const g: Mv(T, kind, lower, packed_) = .{ .n = n, .k = k, .alpha = alpha, .a = a, .lda = lda, .X = X };
    // `spmv.c`, `zspmv.c` and `zhpmv.c` as their full siblings `symv.c` and
    // `zsymv.c`; `sbmv.c`, `zsbmv.c` and `zhbmv.c` thread every call, here
    // `gbmv.c`'s rule for the band matrix `S` is (`kl = ku = k`)
    const p = if (packed_)
        threading.plan(.symv, T, nu, nu, if (lower) .falling else .rising, 8)
    else
        threading.plan(.gbmv, T, threading.product(&.{ nu, nu }), nu, .even, 8);
    if (p.chunks.count > 1 and (packed_ or 2 * k >= 15)) {
        const parts: threading.Partials(T) = try .init(allocator, &p, &g, @TypeOf(g).rows);
        defer parts.deinit(allocator);
        const C = struct {
            g: *const @TypeOf(g),
            parts: *const threading.Partials(T),

            fn one(c: @This(), ch: threading.Chunk) void {
                const part = c.parts.take(ch.index);
                c.g.columns(part.ptr, @intCast(c.parts.rows[ch.index].lo), @intCast(ch.start), @intCast(ch.end));
            }
        };
        try threading.forEach(io, &p, C{ .g = &g, .parts = &parts }, C.one);
        return parts.fold(io, y, incy);
    }
    var sy: Scratch(T) = .{};
    defer sy.deinit(allocator);
    const Y = try sy.unit(allocator, n, y, incy);
    g.columns(Y, 0, 0, n);
    if (incy != 1) copyK(T, n, Y, 1, y, incy);
}

/// `mvDriver`'s columns, over a unit-stride `X`.
fn Mv(comptime T: type, comptime kind: Kind, comptime lower: bool, comptime packed_: bool) type {
    return struct {
        n: BlasLong,
        k: BlasLong,
        alpha: T,
        a: [*]T,
        lda: BlasLong,
        X: [*]T,

        /// The rows columns `c.start..c.end` reach.
        fn rows(g: *const @This(), c: threading.Chunk) threading.Rows {
            const nu: usize = @intCast(g.n);
            const k: usize = if (packed_) nu else @intCast(g.k);
            return if (lower) .{ .lo = c.start, .hi = @min(nu, c.end + k) } else .{ .lo = c.start -| k, .hi = c.end };
        }

        /// Columns `c0..c1` added to `Y`, whose element 0 is row `yb`.
        fn columns(g: *const @This(), Y: [*]T, yb: BlasLong, c0: BlasLong, c1: BlasLong) void {
            const n = g.n;
            const alpha = g.alpha;
            const X = g.X;
            const rev = kind == .hermitian_rev;
            // start of column i in `a`
            var cs: BlasLong = if (packed_) columnAt(lower, true, n, 0, c0) else c0 * g.lda;
            var i: BlasLong = c0;
            while (i < c1) : (i += 1) {
                const room: BlasLong = if (lower) n - 1 - i else i;
                const noff = if (packed_) room else @min(room, g.k);
                const diag = cs + (if (lower) 0 else if (packed_) i else g.k);
                const dp = off(T, g.a, diag);
                const off_a = if (lower) off(T, g.a, diag + 1) else off(T, g.a, diag - noff);
                const off_i = if (lower) i + 1 else i - noff; // the off-diagonals' rows
                const xi = at(T, X, i).*;
                const yi = at(T, Y, i - yb);
                const ax = mul(T, alpha, xi);
                switch (kind) {
                    .symmetric => if (packed_ and lower) {
                        yi.* = add(T, yi.*, mul(T, alpha, dotK(T, false, n - i, dp, off(T, X, i))));
                        axpyK(T, false, noff, ax, off_a, off(T, Y, off_i - yb));
                    } else {
                        axpyK(T, false, noff + 1, ax, if (lower) dp else off_a, off(T, Y, (if (lower) i else off_i) - yb));
                        if (noff > 0) yi.* = add(T, yi.*, mul(T, alpha, dotK(T, false, noff, off_a, off(T, X, off_i))));
                    },
                    .hermitian, .hermitian_rev => {
                        axpyK(T, rev, noff, ax, off_a, off(T, Y, off_i - yb));
                        const d = if (comptime isComplex(T)) dp[0].re else dp[0];
                        yi.* = add(T, yi.*, mul(T, alpha, scaleF(T, xi, d)));
                        if (noff > 0) yi.* = add(T, yi.*, mul(T, alpha, dotK(T, !rev, noff, off_a, off(T, X, off_i))));
                    },
                }
                cs += if (!packed_) g.lda else if (lower) n - i else i + 1;
            }
        }
    };
}

/// Where column `i` of a triangle of order `n` starts, in full (`lda`) or
/// packed storage.
pub fn columnAt(comptime lower: bool, comptime packed_: bool, n: BlasLong, lda: BlasLong, i: BlasLong) BlasLong {
    if (!packed_) return i * (if (lower) lda + 1 else lda);
    return if (lower) i * n - @divExact(i * (i - 1), 2) else @divExact(i * (i + 1), 2);
}

// ---- interface ---------------------------------------------------------------

/// Refuses a matrix of order `n` that could not exist in memory
/// (`util.fitsMatrix`): `lda` columns of band or full storage, or the packed
/// triangle, bounded by `n` columns of `n`.
pub fn spans(comptime T: type, comptime packed_: bool, n: isize, lda: isize) interface.ArgError!void {
    if (packed_) {
        if (!util.fitsMatrix(n, n, n, @sizeOf(T))) return error.InvalidN;
    } else if (!util.fitsMatrix(lda, n, lda, @sizeOf(T))) return error.InvalidLda;
}

/// `?GBMV` (`interface/gbmv.c`, `zgbmv.c`): `y := alpha*op(A)*x + beta*y` for
/// an `m` by `n` matrix with `kl` sub- and `ku` superdiagonals in band
/// storage. `beta == 0` stores `alpha*op(A)*x` without reading `y`.
pub fn gbmv(allocator: Allocator, io: Io, comptime T: type, trans: Trans, m: isize, n: isize, kl: isize, ku: isize, alpha: T, a: [*]T, lda: isize, x: [*]T, incx: isize, beta: T, y: [*]T, incy: isize) Error!void {
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
            inline else => |t| try gbmvDriver(T, comptime t.transposed(), comptime t.conjA(), comptime t.conjX(), allocator, io, m, n, ku, kl, alpha, a, lda, px, incx, py, incy),
        }
    } else switch (tr) {
        inline else => |t| try gbmvDriver(T, t, false, false, allocator, io, m, n, ku, kl, alpha, a, lda, px, incx, py, incy),
    }
}

/// Shared tail of the matrix-vector products: `beta` scaling, `alpha == 0`,
/// start-of-vector pointers, and the uplo dispatch.
fn mvInterface(allocator: Allocator, io: Io, comptime T: type, comptime kind: Kind, comptime packed_: bool, lower: bool, n: isize, k: isize, alpha: T, a: [*]T, lda: isize, x: [*]T, incx: isize, beta: T, y: [*]T, incy: isize) Error!void {
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
        inline else => |lo| try mvDriver(T, kind, lo, packed_, allocator, io, kn, k, alpha, a, lda, px, incx, py, incy),
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
/// diagonal's imaginary part is not read. `v` and `m` take `conj(S)`.
pub fn hbmv(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, k: isize, alpha: T, a: [*]T, lda: isize, x: [*]T, incx: isize, beta: T, y: [*]T, incy: isize) Error!void {
    comptime util.requireComplex(T, "hbmv");
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
    comptime util.requireComplex(T, "hpmv");
    if (uplo.isRev()) return error.InvalidUplo;
    if (n < 0) return error.InvalidN;
    if (incx == 0) return error.InvalidIncx;
    if (incy == 0) return error.InvalidIncy;
    try mvInterface(allocator, io, T, .hermitian, true, uplo.isLower(), n, 0, alpha, a, 0, x, incx, beta, y, incy);
}

// ---- tests -------------------------------------------------------------------

const rank = @import("rank.zig");
const spr = rank.spr;
const syr = rank.syr;
const hpr = rank.hpr;
const her = rank.her;
const spr2 = rank.spr2;
const syr2 = rank.syr2;
const hpr2 = rank.hpr2;
const her2 = rank.her2;

const Vtable = @import("../Vtable.zig");
const check = @import("../_test/check.zig");
const Complex32 = types.Complex32;
const Complex64 = types.Complex64;

const all = check.all_types;
const talloc = std.testing.allocator;
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

const scalarOf = check.val;

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
                try hpmv(talloc, check.io, T, un, nn, alpha, &a, &xs, inc[0], beta, &ys, inc[1]),
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
fn checkRank(comptime T: type, comptime st: Store, comptime herm: bool, lower: bool) !void {
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
            const v = if (herm and two)
                add(T, mul(T, alpha, mul(T, xr, conj(T, yc))), mul(T, conj(T, alpha), mul(T, yr, conj(T, xc))))
            else if (herm)
                scaleF(T, mul(T, xr, conj(T, xc)), alpha_f)
            else if (two)
                mul(T, alpha, add(T, mul(T, xr, yc), mul(T, yr, xc)))
            else
                mul(T, alpha, mul(T, xr, xc));
            want[p] = add(T, want[p], v);
            if (herm) {
                if (r == c_) want[p].im = 0;
            }
        };
        const un: Uplo = if (lower) .l else .u;
        const nn: isize = @intCast(n);
        const lda: isize = @intCast(ld);
        switch (st) {
            .packed_ => if (two) {
                if (!herm)
                    try spr2(talloc, check.io, T, un, nn, a0, &xs, inc[0], &ys, inc[1], &a)
                else
                    try hpr2(talloc, check.io, T, un, nn, a0, &xs, inc[0], &ys, inc[1], &a);
            } else {
                if (!herm)
                    try spr(talloc, check.io, T, un, nn, a0, &xs, inc[0], &a)
                else
                    try hpr(talloc, check.io, T, un, nn, alpha_f, &xs, inc[0], &a);
            },
            .full => if (two) {
                if (!herm)
                    try syr2(talloc, check.io, T, un, nn, a0, &xs, inc[0], &ys, inc[1], &a, lda)
                else
                    try her2(talloc, check.io, T, un, nn, a0, &xs, inc[0], &ys, inc[1], &a, lda);
            } else {
                if (!herm)
                    try syr(talloc, check.io, T, un, nn, a0, &xs, inc[0], &a, lda)
                else
                    try her(talloc, check.io, T, un, nn, alpha_f, &xs, inc[0], &a, lda);
            },
            .band => unreachable,
        }
        check.expectCloseSlices(T, &want, &a, tol(T, 4)) catch |e| {
            check.say("rank{s} {s} herm={} n={d} inc={d},{d} lower={}: ", .{ if (two) "2" else "1", @tagName(st), herm, n, inc[0], inc[1], lower });
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
            inline for (comptime kindsOf(T)) |kind| {
                // `hpmv` takes no `HEMVREV` triangle
                if (st == .packed_ and kind == .hermitian_rev) continue;
                for ([_]bool{ false, true }) |lower| try checkMv(T, st, kind, lower);
            }
        }
        inline for (.{ Store.full, Store.packed_ }) |st| {
            inline for (if (comptime isComplex(T)) .{ false, true } else .{false}) |herm| for ([_]bool{ false, true }) |lower| try checkRank(T, st, herm, lower);
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
    _ = blas.vtable();
    const was = accel.routing();
    _ = accel.route(false);
    defer _ = accel.route(was);
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
    // only `hbmv` of these takes `v` and `m`; upstream's others get them from
    // row-major CBLAS alone
    var ca: [4]Complex64 = @splat(.init(1, 0));
    var cx: [2]Complex64 = @splat(.init(1, 0));
    var cy = cx;
    const c1: Complex64 = .init(1, 0);
    for ([_]Uplo{ .v, .m }) |u| {
        try std.testing.expectError(error.InvalidUplo, hpmv(talloc, check.io, Complex64, u, -1, c1, &ca, &cx, 0, c1, &cy, 0));
        try std.testing.expectError(error.InvalidUplo, her(talloc, check.io, Complex64, u, -1, 1, &cx, 0, &ca, 0));
        try std.testing.expectError(error.InvalidUplo, hpr(talloc, check.io, Complex64, u, -1, 1, &cx, 0, &ca));
        try std.testing.expectError(error.InvalidUplo, her2(talloc, check.io, Complex64, u, -1, c1, &cx, 0, &cy, 0, &ca, 0));
        try std.testing.expectError(error.InvalidUplo, hpr2(talloc, check.io, Complex64, u, -1, c1, &cx, 0, &cy, 0, &ca));
        try hbmv(talloc, check.io, Complex64, u, 2, 1, c1, &ca, 2, &cx, 1, c1, &cy, 1);
    }
    try std.testing.expectError(error.InvalidTrans, gbmv(talloc, check.io, f64, .d, -1, 2, 1, 1, 1, &a, 2, &x, 1, 1, &y, 1));
    // real `r` and `c` are `n` and `t`
    var yn = [_]f64{ 0, 0 };
    var yr = yn;
    try gbmv(talloc, check.io, f64, .t, 2, 2, 1, 0, 1, &a, 2, &x, 1, 0, &yn, 1);
    try gbmv(talloc, check.io, f64, .c, 2, 2, 1, 0, 1, &a, 2, &x, 1, 0, &yr, 1);
    try std.testing.expectEqualSlices(f64, &yn, &yr);
}
