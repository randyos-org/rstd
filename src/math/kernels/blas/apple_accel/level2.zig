//! Accelerate's level-2 routines (`apple_accel.zig`).

const std = @import("std");
const types = @import("../types.zig");
const threading = @import("../threading.zig");
const util = @import("../util.zig");
const accel = @import("../apple_accel.zig");
const binding = @import("binding.zig");
const routines = @import("routines.zig");
const parts = @import("split.zig");
const quirks = @import("quirks.zig");
const speed = @import("speed.zig");

const Io = std.Io;
const isComplex = types.isComplex;
const Float = types.Float;
const isZero = util.isZero;
const x86_64 = binding.x86_64;
const Blas = binding.Blas;
const routine = routines.routine;
const Transpose = binding.Transpose;
const covers = accel.covers;
const use = accel.use;
const call = accel.call;
const span = accel.span;
const mat = accel.mat;
const matrix_quantum = parts.matrix_quantum;
const split2 = parts.split2;
const fits = parts.fits;
const mag = parts.mag;
const sub = parts.sub;
const column = parts.column;
const units = parts.units;
const scalar = parts.scalar;
const transpose = parts.transpose;
const uplo = parts.uplo;
const diag = parts.diag;
const col = parts.col;
const realDiagonal = quirks.realDiagonal;
const routes = speed.routes;
const zeroPart = quirks.zeroPart;
const nonFinite = quirks.nonFinite;
const zeroColumn = quirks.zeroColumn;

/// `gemv`'s parts: rows of `y` without a transpose, columns with one, as
/// `gemv_thread.c` splits them.
fn Gemv(comptime T: type) type {
    return struct {
        const Self = @This();
        const name = routine(.gemv, T);

        b: *const Blas,
        t: Transpose,
        m: isize,
        n: isize,
        alpha: T,
        a: [*]T,
        lda: isize,
        x: [*]T,
        incx: isize,
        beta: T,
        y: [*]T,
        incy: isize,

        fn part(g: *const Self, ch: threading.Chunk) void {
            const tr = g.t != .no_trans;
            const m = if (tr) g.m else units(ch);
            const n = if (tr) units(ch) else g.n;
            const a = if (tr) column(g.a, ch.start, g.lda) else g.a + ch.start;
            const y = sub(g.y, if (tr) g.n else g.m, g.incy, ch);
            _ = call(name, g.b, &.{}, .{ col, g.t, m, n, scalar(T, &g.alpha), a, g.lda, g.x, g.incx, scalar(T, &g.beta), y, g.incy }).?;
        }
    };
}

pub inline fn gemv(io: Io, comptime T: type, trans: types.Trans, m: isize, n: isize, alpha: T, a: [*]T, lda: isize, x: [*]T, incx: isize, beta: T, y: [*]T, incy: isize) Io.Cancelable!bool {
    if (comptime !covers(T) or !routes(.gemv, T)) return false else {
        const b = use(T) orelse return false;
        const t = transpose(T, trans) orelse return false;
        const lenx = if (trans.transposed()) m else n;
        const leny = if (trans.transposed()) n else m;
        // arm64 only
        if (!x86_64 and zeroPart(T, beta) and !isZero(T, beta) and nonFinite(T, leny, y, incy)) return false;
        const G = Gemv(T);
        if (!fits(G.name, b, &.{ m, n, lda, mag(incx), mag(incy), mat(m, n, lda), span(lenx, incx), span(leny, incy) })) return false;
        const g: G = .{ .b = b, .t = t, .m = m, .n = n, .alpha = alpha, .a = a, .lda = lda, .x = x, .incx = incx, .beta = beta, .y = y, .incy = incy };
        try split2(io, .gemv, T, threading.product(&.{ @intCast(m), @intCast(n) }), @intCast(leny), matrix_quantum, threading.product(&.{ @intCast(m), @intCast(n), @sizeOf(T) }), &g, G.part);
        return true;
    }
}

pub inline fn gbmv(comptime T: type, trans: types.Trans, m: isize, n: isize, kl: isize, ku: isize, alpha: T, a: [*]T, lda: isize, x: [*]T, incx: isize, beta: T, y: [*]T, incy: isize) bool {
    if (comptime !covers(T) or !routes(.gbmv, T)) return false else {
        const b = use(T) orelse return false;
        const t = transpose(T, trans) orelse return false;
        const lenx = if (trans.transposed()) m else n;
        const leny = if (trans.transposed()) n else m;
        const ext = &.{ mat(lda, n, lda), span(lenx, incx), span(leny, incy) };
        return call(routine(.gbmv, T), b, ext, .{ col, t, m, n, kl, ku, scalar(T, &alpha), a, lda, x, incx, scalar(T, &beta), y, incy }) != null;
    }
}

/// `ger`'s parts: columns of `A`, as `ger_thread.c` splits them.
fn Ger(comptime T: type, comptime conj: bool) type {
    return struct {
        const Self = @This();
        const name = routine(if (!isComplex(T)) .ger else if (conj) .gerc else .geru, T);

        b: *const Blas,
        m: isize,
        n: isize,
        alpha: T,
        x: [*]T,
        incx: isize,
        y: [*]T,
        incy: isize,
        a: [*]T,
        lda: isize,

        fn part(g: *const Self, ch: threading.Chunk) void {
            _ = call(name, g.b, &.{}, .{ col, g.m, units(ch), scalar(T, &g.alpha), g.x, g.incx, sub(g.y, g.n, g.incy, ch), g.incy, column(g.a, ch.start, g.lda), g.lda }).?;
        }
    };
}

/// `ger`, `geru`, and `gerc` for `conj`.
pub inline fn ger(io: Io, comptime T: type, comptime conj: bool, m: isize, n: isize, alpha: T, x: [*]T, incx: isize, y: [*]T, incy: isize, a: [*]T, lda: isize) Io.Cancelable!bool {
    // `gerc` follows `ger`, the one timed
    if (comptime !covers(T) or !routes(.ger, T)) return false else {
        const b = use(T) orelse return false;
        if (nonFinite(T, m, x, incx)) return false;
        const G = Ger(T, conj);
        if (!fits(G.name, b, &.{ m, n, lda, mag(incx), mag(incy), span(m, incx), span(n, incy), mat(m, n, lda) })) return false;
        const g: G = .{ .b = b, .m = m, .n = n, .alpha = alpha, .x = x, .incx = incx, .y = y, .incy = incy, .a = a, .lda = lda };
        try split2(io, .ger, T, threading.product(&.{ @intCast(m), @intCast(n) }), @intCast(n), matrix_quantum, threading.product(&.{ @intCast(m), @intCast(n), @sizeOf(T) }), &g, G.part);
        return true;
    }
}

/// `symv`, or `hemv` for `herm`.
pub inline fn symv(comptime T: type, comptime herm: bool, ul: types.Uplo, n: isize, alpha: T, a: [*]T, lda: isize, x: [*]T, incx: isize, beta: T, y: [*]T, incy: isize) bool {
    if (comptime !covers(T) or herm != isComplex(T) or !routes(if (herm) .hemv else .symv, T)) return false else {
        const name = comptime routine(if (herm) .hemv else .symv, T);
        const b = use(T) orelse return false;
        const u = uplo(ul) orelse return false;
        // real `symv`'s shows on arm64 only (issue #410)
        if (!x86_64 and !isComplex(T) and nonFinite(T, n, x, incx)) return false;
        return call(name, b, &.{ mat(n, n, lda), span(n, incx), span(n, incy) }, .{ col, u, n, scalar(T, &alpha), a, lda, x, incx, scalar(T, &beta), y, incy }) != null;
    }
}

/// `sbmv` or `hbmv`, packed (`spmv`, `hpmv`) for `packed_`.
pub inline fn sbmv(comptime T: type, comptime herm: bool, comptime packed_: bool, ul: types.Uplo, n: isize, k: isize, alpha: T, a: [*]T, lda: isize, x: [*]T, incx: isize, beta: T, y: [*]T, incy: isize) bool {
    if (comptime !covers(T) or herm != isComplex(T) or !routes(if (herm) (if (packed_) .hpmv else .hbmv) else (if (packed_) .spmv else .sbmv), T)) return false else {
        const name = comptime routine(if (herm) (if (packed_) .hpmv else .hbmv) else (if (packed_) .spmv else .sbmv), T);
        const b = use(T) orelse return false;
        const u = uplo(ul) orelse return false;
        const al = scalar(T, &alpha);
        const be = scalar(T, &beta);
        if (packed_) return call(name, b, &.{ @divExact(n * (n + 1), 2), span(n, incx), span(n, incy) }, .{ col, u, n, al, a, x, incx, be, y, incy }) != null;
        return call(name, b, &.{ mat(lda, n, lda), span(n, incx), span(n, incy) }, .{ col, u, n, k, al, a, lda, x, incx, be, y, incy }) != null;
    }
}

pub const TriKind = enum { trmv, trsv, tbmv, tbsv, tpmv, tpsv };

pub inline fn tri(comptime T: type, comptime kind: TriKind, ul: types.Uplo, trans: types.Trans, dg: types.Diag, n: isize, k: isize, a: [*]T, lda: isize, x: [*]T, incx: isize) bool {
    if (comptime !covers(T) or !routes(switch (kind) {
        .trmv => .trmv,
        .trsv => .trsv,
        .tbmv => .tbmv,
        .tbsv => .tbsv,
        .tpmv => .tpmv,
        .tpsv => .tpsv,
    }, T)) return false else {
        const b = use(T) orelse return false;
        const u = uplo(ul) orelse return false;
        const t = transpose(T, trans) orelse return false;
        // real `trmv`'s shows on arm64 only
        if (!(x86_64 and kind == .trmv and !isComplex(T)) and nonFinite(T, n, x, incx)) return false;
        const name = comptime routine(switch (kind) {
            .trmv => .trmv,
            .trsv => .trsv,
            .tbmv => .tbmv,
            .tbsv => .tbsv,
            .tpmv => .tpmv,
            .tpsv => .tpsv,
        }, T);
        return switch (kind) {
            .trmv, .trsv => call(name, b, &.{ mat(n, n, lda), span(n, incx) }, .{ col, u, t, diag(dg), n, a, lda, x, incx }),
            .tbmv, .tbsv => call(name, b, &.{ mat(lda, n, lda), span(n, incx) }, .{ col, u, t, diag(dg), n, k, a, lda, x, incx }),
            .tpmv, .tpsv => call(name, b, &.{ @divExact(n * (n + 1), 2), span(n, incx) }, .{ col, u, t, diag(dg), n, a, x, incx }),
        } != null;
    }
}

/// `syr` or `her` (real `alpha`), packed (`spr`, `hpr`) for `packed_`.
pub inline fn syr(comptime T: type, comptime herm: bool, comptime packed_: bool, ul: types.Uplo, n: isize, alpha: if (herm) Float(T) else T, x: [*]T, incx: isize, a: [*]T, lda: isize) bool {
    if (comptime !covers(T) or herm != isComplex(T) or !routes(if (herm) (if (packed_) .hpr else .her) else (if (packed_) .spr else .syr), T)) return false else {
        const name = comptime routine(if (herm) (if (packed_) .hpr else .her) else (if (packed_) .spr else .syr), T);
        const b = use(T) orelse return false;
        const u = uplo(ul) orelse return false;
        if (zeroColumn(T, n, x, incx, null, 0) and nonFinite(T, n, x, incx)) return false;
        const ok = if (packed_)
            call(name, b, &.{ span(n, incx), @divExact(n * (n + 1), 2) }, .{ col, u, n, alpha, x, incx, a }) != null
        else
            call(name, b, &.{ span(n, incx), mat(n, n, lda) }, .{ col, u, n, alpha, x, incx, a, lda }) != null;
        if (ok and herm) realDiagonal(T, n, a, if (packed_) 0 else lda, ul == .l);
        return ok;
    }
}

/// `syr2` or `her2`, packed (`spr2`, `hpr2`) for `packed_`.
pub inline fn syr2(comptime T: type, comptime herm: bool, comptime packed_: bool, ul: types.Uplo, n: isize, alpha: T, x: [*]T, incx: isize, y: [*]T, incy: isize, a: [*]T, lda: isize) bool {
    if (comptime !covers(T) or herm != isComplex(T) or !routes(if (herm) (if (packed_) .hpr2 else .her2) else (if (packed_) .spr2 else .syr2), T)) return false else {
        const name = comptime routine(if (herm) (if (packed_) .hpr2 else .her2) else (if (packed_) .spr2 else .syr2), T);
        const b = use(T) orelse return false;
        const u = uplo(ul) orelse return false;
        if (zeroColumn(T, n, x, incx, y, incy) and (nonFinite(T, n, x, incx) or nonFinite(T, n, y, incy))) return false;
        const al = scalar(T, &alpha);
        const ok = if (packed_)
            call(name, b, &.{ span(n, incx), span(n, incy), @divExact(n * (n + 1), 2) }, .{ col, u, n, al, x, incx, y, incy, a }) != null
        else
            call(name, b, &.{ span(n, incx), span(n, incy), mat(n, n, lda) }, .{ col, u, n, al, x, incx, y, incy, a, lda }) != null;
        if (ok and herm) realDiagonal(T, n, a, if (packed_) 0 else lda, ul == .l);
        return ok;
    }
}
