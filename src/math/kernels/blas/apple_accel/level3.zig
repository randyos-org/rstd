//! Accelerate's level-3 routines (`apple_accel.zig`).

const std = @import("std");
const types = @import("../types.zig");
const threading = @import("../threading.zig");
const accel = @import("../apple_accel.zig");
const binding = @import("binding.zig");
const routines = @import("routines.zig");
const parts = @import("split.zig");
const quirks = @import("quirks.zig");
const speed = @import("speed.zig");

const Io = std.Io;
const isComplex = types.isComplex;
const Blas = binding.Blas;
const routine = routines.routine;
const Transpose = binding.Transpose;
const Uplo = binding.Uplo;
const Diag = binding.Diag;
const Side = binding.Side;
const covers = accel.covers;
const rosetta = accel.rosetta;
const use = accel.use;
const call = accel.call;
const mat = accel.mat;
const matrix_quantum = parts.matrix_quantum;
const split = parts.split;
const fits = parts.fits;
const column = parts.column;
const units = parts.units;
const scalar = parts.scalar;
const transpose = parts.transpose;
const uplo = parts.uplo;
const diag = parts.diag;
const side = parts.side;
const col = parts.col;
const realDiagonal = quirks.realDiagonal;
const nonFiniteTri = quirks.nonFiniteTri;
const tinyDiagonal = quirks.tinyDiagonal;
const nonFiniteMat = quirks.nonFiniteMat;
const routes = speed.routes;

/// `gemm`'s parts: columns of `C` if it has as many as rows, else rows;
/// each takes the matching part of `op(B)` or `op(A)` and all of the other.
fn Gemm(comptime T: type) type {
    return struct {
        const Self = @This();
        const name = routine(.gemm, T);

        b: *const Blas,
        ta: Transpose,
        tb: Transpose,
        by_columns: bool,
        m: isize,
        n: isize,
        k: isize,
        alpha: T,
        a: [*]const T,
        lda: isize,
        b_: [*]const T,
        ldb: isize,
        beta: T,
        c: [*]T,
        ldc: isize,

        fn part(g: *const Self, ch: threading.Chunk) void {
            var m = g.m;
            var n = g.n;
            var a = g.a;
            var b_ = g.b_;
            var c = g.c;
            if (g.by_columns) {
                n = units(ch);
                b_ = if (g.tb == .no_trans) column(b_, ch.start, g.ldb) else b_ + ch.start;
                c = column(c, ch.start, g.ldc);
            } else {
                m = units(ch);
                a = if (g.ta == .no_trans) a + ch.start else column(a, ch.start, g.lda);
                c += ch.start;
            }
            _ = call(name, g.b, &.{}, .{ col, g.ta, g.tb, m, n, g.k, scalar(T, &g.alpha), a, g.lda, b_, g.ldb, scalar(T, &g.beta), c, g.ldc }).?;
        }
    };
}

pub inline fn gemm(io: Io, comptime T: type, transa: types.Trans, transb: types.Trans, m: isize, n: isize, k: isize, alpha: T, a: [*]const T, lda: isize, b_: [*]const T, ldb: isize, beta: T, c: [*]T, ldc: isize) Io.Cancelable!bool {
    if (comptime !covers(T) or !routes(.gemm, T)) return false else {
        const b = use(T) orelse return false;
        const ta = transpose(T, transa) orelse return false;
        const tb = transpose(T, transb) orelse return false;
        if (rosetta() and (if (transa.transposed()) nonFiniteMat(T, k, m, a, lda) else nonFiniteMat(T, m, k, a, lda))) return false;
        const ea = if (transa.transposed()) mat(k, m, lda) else mat(m, k, lda);
        const eb = if (transb.transposed()) mat(n, k, ldb) else mat(k, n, ldb);
        const G = Gemm(T);
        if (!fits(G.name, b, &.{ m, n, k, lda, ldb, ldc, ea, eb, mat(m, n, ldc) })) return false;
        const g: G = .{ .b = b, .ta = ta, .tb = tb, .by_columns = n >= m, .m = m, .n = n, .k = k, .alpha = alpha, .a = a, .lda = lda, .b_ = b_, .ldb = ldb, .beta = beta, .c = c, .ldc = ldc };
        const work = threading.product(&.{ @intCast(m), @intCast(n), @intCast(k) });
        try split(io, .gemm, T, work, @intCast(@max(m, n)), matrix_quantum, &g, G.part);
        return true;
    }
}

/// `symm`'s parts: columns of `B` and `C` for `S` on the left, rows on the
/// right, each with all of `S`.
fn Symm(comptime T: type, comptime herm: bool) type {
    return struct {
        const Self = @This();
        const name = routine(if (herm) .hemm else .symm, T);

        b: *const Blas,
        sd: Side,
        u: Uplo,
        m: isize,
        n: isize,
        alpha: T,
        a: [*]const T,
        lda: isize,
        b_: [*]const T,
        ldb: isize,
        beta: T,
        c: [*]T,
        ldc: isize,

        fn part(g: *const Self, ch: threading.Chunk) void {
            const left = g.sd == .left;
            const m = if (left) g.m else units(ch);
            const n = if (left) units(ch) else g.n;
            const b_ = if (left) column(g.b_, ch.start, g.ldb) else g.b_ + ch.start;
            const c = if (left) column(g.c, ch.start, g.ldc) else g.c + ch.start;
            _ = call(name, g.b, &.{}, .{ col, g.sd, g.u, m, n, scalar(T, &g.alpha), g.a, g.lda, b_, g.ldb, scalar(T, &g.beta), c, g.ldc }).?;
        }
    };
}

/// `symm`, or `hemm` for `herm`.
pub inline fn symm(io: Io, comptime T: type, comptime herm: bool, sd: types.Side, ul: types.Uplo, m: isize, n: isize, alpha: T, a: [*]const T, lda: isize, b_: [*]const T, ldb: isize, beta: T, c: [*]T, ldc: isize) Io.Cancelable!bool {
    if (comptime !covers(T) or (herm and !isComplex(T)) or !routes(if (herm) .hemm else .symm, T)) return false else {
        const b = use(T) orelse return false;
        const u = uplo(ul) orelse return false;
        const ka = if (sd == .l) m else n;
        const S = Symm(T, herm);
        if (!fits(S.name, b, &.{ m, n, lda, ldb, ldc, mat(ka, ka, lda), mat(m, n, ldb), mat(m, n, ldc) })) return false;
        const g: S = .{ .b = b, .sd = side(sd), .u = u, .m = m, .n = n, .alpha = alpha, .a = a, .lda = lda, .b_ = b_, .ldb = ldb, .beta = beta, .c = c, .ldc = ldc };
        // the work `blas3/interface.zig` thresholds
        const work = threading.product(&.{ 2, @intCast(m), @intCast(m), @intCast(n) });
        try split(io, .symm, T, work, @intCast(if (sd == .l) n else m), matrix_quantum, &g, S.part);
        return true;
    }
}

/// `syrk`, `herk` (real `alpha` and `beta`), `syr2k` or `her2k` (real
/// `beta`), as `herm` and `b_` say; `t` is a transposed `op`.
pub inline fn syrk(comptime T: type, comptime herm: bool, ul: types.Uplo, t: bool, n: isize, k: isize, alpha: T, a: [*]const T, lda: isize, b_: ?[*]const T, ldb: isize, beta: T, c: [*]T, ldc: isize) bool {
    if (comptime !covers(T) or (herm and !isComplex(T))) return false else {
        if (b_ == null) {
            if (comptime !routes(if (herm) .herk else .syrk, T)) return false;
        } else if (comptime !routes(if (herm) .her2k else .syr2k, T)) return false;
        const b = use(T) orelse return false;
        const u = uplo(ul) orelse return false;
        const tr: Transpose = if (!t) .no_trans else if (herm) .conj_trans else .trans;
        if (comptime herm) {
            if (beta.re == 1 and (k == 0 or (alpha.re == 0 and alpha.im == 0))) return false;
            if (beta.re != 0 and nonFiniteTri(T, n, c, ldc, ul == .l)) return false;
        }
        const ea = if (t) mat(k, n, lda) else mat(n, k, lda);
        const ec = mat(n, n, ldc);
        const ok = if (b_) |bb| blk: {
            const eb = if (t) mat(k, n, ldb) else mat(n, k, ldb);
            const ext = &.{ ea, eb, ec };
            if (comptime herm) break :blk call(routine(.her2k, T), b, ext, .{ col, u, tr, n, k, scalar(T, &alpha), a, lda, bb, ldb, beta.re, c, ldc }) != null;
            break :blk call(routine(.syr2k, T), b, ext, .{ col, u, tr, n, k, scalar(T, &alpha), a, lda, bb, ldb, scalar(T, &beta), c, ldc }) != null;
        } else if (comptime herm)
            call(routine(.herk, T), b, &.{ ea, ec }, .{ col, u, tr, n, k, alpha.re, a, lda, beta.re, c, ldc }) != null
        else
            call(routine(.syrk, T), b, &.{ ea, ec }, .{ col, u, tr, n, k, scalar(T, &alpha), a, lda, scalar(T, &beta), c, ldc }) != null;
        if (ok and herm) realDiagonal(T, n, c, ldc, ul == .l);
        return ok;
    }
}

/// `trmm`'s parts: columns of `B` for `A` on the left, rows on the right,
/// as `gemm_thread_n` and `gemm_thread_m` split them.
fn Trmm(comptime T: type, comptime solve: bool) type {
    return struct {
        const Self = @This();
        const name = routine(if (solve) .trsm else .trmm, T);

        b: *const Blas,
        sd: Side,
        u: Uplo,
        t: Transpose,
        dg: Diag,
        m: isize,
        n: isize,
        alpha: T,
        a: [*]const T,
        lda: isize,
        b_: [*]T,
        ldb: isize,

        fn part(g: *const Self, ch: threading.Chunk) void {
            const left = g.sd == .left;
            const m = if (left) g.m else units(ch);
            const n = if (left) units(ch) else g.n;
            const b_ = if (left) column(g.b_, ch.start, g.ldb) else g.b_ + ch.start;
            _ = call(name, g.b, &.{}, .{ col, g.sd, g.u, g.t, g.dg, m, n, scalar(T, &g.alpha), g.a, g.lda, b_, g.ldb }).?;
        }
    };
}

/// `trmm`, or `trsm` for `solve`.
pub inline fn trmm(io: Io, comptime T: type, comptime solve: bool, sd: types.Side, ul: types.Uplo, transa: types.Trans, dg: types.Diag, m: isize, n: isize, alpha: T, a: [*]const T, lda: isize, b_: [*]T, ldb: isize) Io.Cancelable!bool {
    if (comptime !covers(T) or !routes(if (solve) .trsm else .trmm, T)) return false else {
        const b = use(T) orelse return false;
        const u = uplo(ul) orelse return false;
        const t = transpose(T, transa) orelse return false;
        const ka = if (sd == .l) m else n;
        if (nonFiniteMat(T, m, n, b_, ldb)) return false;
        if (solve and dg == .n and tinyDiagonal(T, ka, a, lda)) return false;
        const R = Trmm(T, solve);
        if (!fits(R.name, b, &.{ m, n, lda, ldb, mat(ka, ka, lda), mat(m, n, ldb) })) return false;
        const g: R = .{ .b = b, .sd = side(sd), .u = u, .t = t, .dg = diag(dg), .m = m, .n = n, .alpha = alpha, .a = a, .lda = lda, .b_ = b_, .ldb = ldb };
        const work = threading.product(&.{ @intCast(m), @intCast(n) });
        try split(io, .trsm, T, work, @intCast(if (sd == .l) n else m), matrix_quantum, &g, R.part);
        return true;
    }
}
