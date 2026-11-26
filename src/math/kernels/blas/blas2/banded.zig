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
const ops = @import("ops.zig");
const Num = ops.Num;
const Ops = ops.Ops;

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
            Ops(T).init().copy(n, v, inc, buf.ptr, 1);
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
    if (incy != 1) Ops(T).init().copy(leny, g.Y, 1, y, incy);
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
            const N = Num(T);
            const o: Ops(T) = .init();
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
                        o.axpy(cj, length, N.mul(g.alpha, if (xc) N.conj(xi, true) else xi, false), ap, off(T, g.Y, start - offset_u - g.yb));
                    } else {
                        const t = if (comptime !xc) o.dot(cj, length, ap, xs) else o.dot(!cj, length, xs, ap);
                        const v = N.mul(g.alpha, if (comptime xc and cj) N.conj(t, true) else t, false);
                        at(T, g.Y, i).* = N.add(at(T, g.Y, i).*, v);
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
    if (incy != 1) Ops(T).init().copy(n, Y, 1, y, incy);
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
            const N = Num(T);
            const o: Ops(T) = .init();
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
                const ax = N.mul(alpha, xi, false);
                switch (kind) {
                    .symmetric => if (packed_ and lower) {
                        yi.* = N.add(yi.*, N.mul(alpha, o.dot(false, n - i, dp, off(T, X, i)), false));
                        o.axpy(false, noff, ax, off_a, off(T, Y, off_i - yb));
                    } else {
                        o.axpy(false, noff + 1, ax, if (lower) dp else off_a, off(T, Y, (if (lower) i else off_i) - yb));
                        if (noff > 0) yi.* = N.add(yi.*, N.mul(alpha, o.dot(false, noff, off_a, off(T, X, off_i)), false));
                    },
                    .hermitian, .hermitian_rev => {
                        o.axpy(rev, noff, ax, off_a, off(T, Y, off_i - yb));
                        const d = if (comptime isComplex(T)) dp[0].re else dp[0];
                        yi.* = N.add(yi.*, N.mul(alpha, N.scale(xi, d), false));
                        if (noff > 0) yi.* = N.add(yi.*, N.mul(alpha, o.dot(!rev, noff, off_a, off(T, X, off_i)), false));
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
    if (!isOne(T, beta)) Ops(T).init().scal(leny, beta, y, @intCast(@abs(incy)));
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
    if (!isOne(T, beta)) Ops(T).init().scal(kn, beta, y, @intCast(@abs(incy)));
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
