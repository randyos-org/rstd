//! The level-3 entries: `gemm`, `symm`, `syrk`, `trmm`, `trsm`.

const rstd = @import("rstd");
const bb = @import("../blas_bench.zig");
const libs = @import("../libs.zig");
const xr = @import("../ext.zig");
const pr = @import("pair.zig");
const wk = @import("work.zig");
const rc = @import("race.zig");

const blas = rstd.math.kernels.blas;
const Prec = bb.Prec;
const cast = libs.cast;
const Num = libs.Num;
const rstd_alloc = libs.rstd_alloc;
const Buf = xr.Buf;
const Fill = xr.Fill;
const Restore = xr.Restore;
const Shape = xr.Shape;
const Ch = pr.Ch;
const In = pr.In;
const Pair = pr.Pair;
const Sc = pr.Sc;
const cadd = pr.cadd;
const cdist = pr.cdist;
const cmag = pr.cmag;
const cmul = pr.cmul;
const fc = pr.fc;
const num = pr.num;
const sc = pr.sc;
const Work = wk.Work;
const symAt = rc.symAt;

/// `c0 + A B` at (i, j) against `c`, `A` read by `a(w, i, k)`.
pub fn productAt(comptime p: Prec, w: *const Work(p), i: usize, j: usize, comptime a: fn (*const Work(p), usize, usize) Pair) [2]f64 {
    var s = w.at(.c0, i + j * w.n);
    var mag = cmag(s);
    for (0..w.n) |k| {
        const x = a(w, i, k);
        const y = w.at(.b0, k + j * w.n);
        s = cadd(s, cmul(x, y));
        mag += cmag(x) * cmag(y);
    }
    return .{ cdist(s, w.at(.c, i + j * w.n)), mag };
}

/// `trmm`'s and `trsm`'s Fortran signature: side, uplo, transa, diag, m,
/// n, alpha, a, lda, b, ldb.
pub fn SigTri(comptime p: Prec) type {
    return fn (Ch, Ch, Ch, Ch, In, In, *const Sc(p), [*]p.F, In, [*]p.F, In, usize, usize, usize, usize) callconv(.c) void;
}

/// `trmm`'s and `trsm`'s call: left, upper, no transpose, non-unit, `alpha` 1.
pub fn libTri(comptime p: Prec, f: *const SigTri(p), w: *Work(p)) void {
    const n = w.ni();
    const one = sc(p, 1);
    f(&fc.L, &fc.U, &fc.N, &fc.N, &n, &n, &one, w.a.ptr, &n, w.b.ptr, &n, 1, 1, 1, 1);
}

/// Row i of upper triangular `A` (`a0`) times column j of `x`, and the
/// size of its terms.
pub fn triRow(comptime p: Prec, w: *const Work(p), comptime x: Buf, i: usize, j: usize) [2]Pair {
    const n = w.n;
    var s: Pair = .{ 0, 0 };
    var mag: f64 = 0;
    for (i..n) |k| {
        const a = w.at(.a0, i + k * n);
        const v = w.at(x, k + j * n);
        s = cadd(s, cmul(a, v));
        mag += cmag(a) * cmag(v);
    }
    return .{ s, .{ mag, 0 } };
}

pub const Gemm3 = struct {
    pub const name = "gemm";
    pub const shape: Shape = .mm;
    pub const out: Buf = .c;
    pub fn flops(n: f64) f64 {
        return 2 * n * n * n;
    }
    pub fn Sig(comptime p: Prec) type {
        return fn (Ch, Ch, In, In, In, *const Sc(p), [*]p.F, In, [*]p.F, In, *const Sc(p), [*]p.F, In, usize, usize) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        const one = sc(p, 1);
        f(&fc.N, &fc.N, &n, &n, &n, &one, w.a.ptr, &n, w.b.ptr, &n, &one, w.c.ptr, &n, 1, 1);
    }
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        const n = w.in();
        const u: usize = @intCast(n);
        try blas.gemm(rstd_alloc, libs.rstd_io, Num(p), .n, .n, u, u, u, num(p, 1), libs.matrix(Num(p), cast(p, w.a.ptr), n, n, n), u, libs.matrix(Num(p), cast(p, w.b.ptr), n, n, n), u, num(p, 1), libs.matrix(Num(p), cast(p, w.c.ptr), n, n, n), u);
    }
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        return productAt(p, w, i, j, struct {
            fn a(v: *const Work(p), r: usize, k: usize) Pair {
                return v.at(.a0, r + k * v.n);
            }
        }.a);
    }
};

pub const Symm = struct {
    pub const name = "symm";
    pub const shape: Shape = .mm;
    pub const out: Buf = .c;
    pub fn flops(n: f64) f64 {
        return 2 * n * n * n;
    }
    pub fn Sig(comptime p: Prec) type {
        return fn (Ch, Ch, In, In, *const Sc(p), [*]p.F, In, [*]p.F, In, *const Sc(p), [*]p.F, In, usize, usize) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        const one = sc(p, 1);
        f(&fc.L, &fc.U, &n, &n, &one, w.a.ptr, &n, w.b.ptr, &n, &one, w.c.ptr, &n, 1, 1);
    }
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        const n = w.in();
        const u: usize = @intCast(n);
        try blas.symm(rstd_alloc, libs.rstd_io, Num(p), .l, .u, u, u, num(p, 1), libs.matrix(Num(p), cast(p, w.a.ptr), n, n, n), u, libs.matrix(Num(p), cast(p, w.b.ptr), n, n, n), u, num(p, 1), libs.matrix(Num(p), cast(p, w.c.ptr), n, n, n), u);
    }
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        return productAt(p, w, i, j, struct {
            fn a(v: *const Work(p), r: usize, k: usize) Pair {
                return symAt(p, v, .a0, false, r, k);
            }
        }.a);
    }
};

pub const Syrk = struct {
    pub const name = "syrk";
    pub const shape: Shape = .mm;
    pub const out: Buf = .c;
    pub fn flops(n: f64) f64 {
        return n * n * (n + 1);
    }
    pub fn Sig(comptime p: Prec) type {
        return fn (Ch, Ch, In, In, *const Sc(p), [*]p.F, In, *const Sc(p), [*]p.F, In, usize, usize) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        const one = sc(p, 1);
        f(&fc.U, &fc.N, &n, &n, &one, w.a.ptr, &n, &one, w.c.ptr, &n, 1, 1);
    }
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        const n = w.in();
        const u: usize = @intCast(n);
        try blas.syrk(rstd_alloc, libs.rstd_io, Num(p), .u, .n, u, u, num(p, 1), libs.matrix(Num(p), cast(p, w.a.ptr), n, n, n), u, num(p, 1), libs.matrix(Num(p), cast(p, w.c.ptr), n, n, n), u);
    }
    /// `c0 + A Aᵀ`, upper triangle.
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        if (i > j) return null;
        const n = w.n;
        var s = w.at(.c0, i + j * n);
        var mag = cmag(s);
        for (0..n) |k| {
            const x = w.at(.a0, i + k * n);
            const y = w.at(.a0, j + k * n);
            s = cadd(s, cmul(x, y));
            mag += cmag(x) * cmag(y);
        }
        return .{ cdist(s, w.at(.c, i + j * n)), mag };
    }
};

pub const Trmm = struct {
    pub const name = "trmm";
    pub const shape: Shape = .mm;
    pub const restore: Restore = .{ .b = true };
    pub const out: Buf = .b;
    pub fn flops(n: f64) f64 {
        return n * n * n;
    }
    pub const Sig = SigTri;
    pub const lib = libTri;
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        const n = w.in();
        const u: usize = @intCast(n);
        try blas.trmm(rstd_alloc, libs.rstd_io, Num(p), .l, .u, .n, .n, u, u, num(p, 1), libs.matrix(Num(p), cast(p, w.a.ptr), n, n, n), u, libs.matrix(Num(p), cast(p, w.b.ptr), n, n, n), u);
    }
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        const r = triRow(p, w, .b0, i, j);
        return .{ cdist(r[0], w.at(.b, i + j * w.n)), r[1][0] };
    }
};

pub const Trsm = struct {
    pub const name = "trsm";
    pub const shape: Shape = .mm;
    pub const fill: Fill = .dominant;
    pub const restore: Restore = .{ .b = true };
    pub const out: Buf = .b;
    pub fn flops(n: f64) f64 {
        return n * n * n;
    }
    pub const Sig = SigTri;
    pub const lib = libTri;
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        const n = w.in();
        const u: usize = @intCast(n);
        try blas.trsm(rstd_alloc, libs.rstd_io, Num(p), .l, .u, .n, .n, u, u, num(p, 1), libs.matrix(Num(p), cast(p, w.a.ptr), n, n, n), u, libs.matrix(Num(p), cast(p, w.b.ptr), n, n, n), u);
    }
    /// `A X` against the right-hand sides.
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        const r = triRow(p, w, .b, i, j);
        const b = w.at(.b0, i + j * w.n);
        return .{ cdist(r[0], b), r[1][0] + cmag(b) };
    }
};
