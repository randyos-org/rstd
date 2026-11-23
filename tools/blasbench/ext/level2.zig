//! The level-2 entries: `ger`, `symv`, `trsv`, `gbmv`.

const rstd = @import("rstd");
const bb = @import("../blas_bench.zig");
const kernels = @import("../kernels.zig");
const libs = @import("../libs.zig");
const xr = @import("../ext.zig");
const pr = @import("pair.zig");
const wk = @import("work.zig");
const rc = @import("race.zig");

const blas = rstd.math.kernels.blas;
const Prec = bb.Prec;
const axpy_alpha = kernels.axpy_alpha;
const gemv_alpha = kernels.gemv_alpha;
const cast = libs.cast;
const Num = libs.Num;
const rstd_alloc = libs.rstd_alloc;
const Buf = xr.Buf;
const Fill = xr.Fill;
const Restore = xr.Restore;
const Shape = xr.Shape;
const band_k = xr.band_k;
const Ch = pr.Ch;
const In = pr.In;
const Pair = pr.Pair;
const Sc = pr.Sc;
const bandRows = pr.bandRows;
const band_kc = pr.band_kc;
const cadd = pr.cadd;
const cdist = pr.cdist;
const cmag = pr.cmag;
const cmul = pr.cmul;
const cnum = pr.cnum;
const fc = pr.fc;
const num = pr.num;
const one_c = pr.one_c;
const sc = pr.sc;
const Work = wk.Work;
const symAt = rc.symAt;

// ---- level 2 entries

pub const Ger = struct {
    pub const name = "ger";
    pub const shape: Shape = .mv;
    pub const out: Buf = .a;
    pub fn flops(n: f64) f64 {
        return 2 * n * n;
    }
    pub fn sym(comptime p: Prec) []const u8 {
        return p.name ++ (if (p.complex) "geru" else "ger");
    }
    pub fn Sig(comptime p: Prec) type {
        return fn (In, In, *const Sc(p), [*]p.F, In, [*]p.F, In, [*]p.F, In) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        const al = sc(p, axpy_alpha);
        f(&n, &n, &al, w.b.ptr, &one_c, w.c.ptr, &one_c, w.a.ptr, &n);
    }
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        const n = w.in();
        try blas.ger(rstd_alloc, libs.rstd_io, Num(p), num(p, axpy_alpha), libs.view(Num(p), cast(p, w.b.ptr), n, 1).toConst(), libs.view(Num(p), cast(p, w.c.ptr), n, 1).toConst(), libs.matrix(Num(p), cast(p, w.a.ptr), n, n, n), @intCast(n));
    }
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        const k = i + j * w.n;
        const al = cnum(p, axpy_alpha);
        const x = w.at(.b0, i);
        const y = w.at(.c0, j);
        return .{ cdist(cadd(w.at(.a0, k), cmul(al, cmul(x, y))), w.at(.a, k)), cmag(w.at(.a0, k)) + cmag(al) * cmag(x) * cmag(y) };
    }
};

pub const Symv = struct {
    pub const name = "symv";
    pub const shape: Shape = .mv;
    pub const out: Buf = .c;
    pub fn flops(n: f64) f64 {
        return 2 * n * n;
    }
    pub fn Sig(comptime p: Prec) type {
        return fn (Ch, In, *const Sc(p), [*]p.F, In, [*]p.F, In, *const Sc(p), [*]p.F, In, usize) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        const al = sc(p, gemv_alpha);
        const one = sc(p, 1);
        f(&fc.U, &n, &al, w.a.ptr, &n, w.b.ptr, &one_c, &one, w.c.ptr, &one_c, 1);
    }
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        const n = w.in();
        try blas.symv(rstd_alloc, libs.rstd_io, Num(p), .u, num(p, gemv_alpha), libs.matrix(Num(p), cast(p, w.a.ptr), n, n, n), @intCast(n), libs.view(Num(p), cast(p, w.b.ptr), n, 1).toConst(), num(p, 1), libs.view(Num(p), cast(p, w.c.ptr), n, 1));
    }
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, _: usize) ?[2]f64 {
        var s: Pair = .{ 0, 0 };
        var mag: f64 = 0;
        for (0..w.n) |k| {
            const a = symAt(p, w, .a0, false, i, k);
            const x = w.at(.b0, k);
            s = cadd(s, cmul(a, x));
            mag += cmag(a) * cmag(x);
        }
        const al = cnum(p, gemv_alpha);
        const y0 = w.at(.c0, i);
        return .{ cdist(cadd(y0, cmul(al, s)), w.at(.c, i)), cmag(y0) + cmag(al) * mag };
    }
};

pub const Trsv = struct {
    pub const name = "trsv";
    pub const shape: Shape = .mv;
    pub const fill: Fill = .dominant;
    pub const restore: Restore = .{ .b = true };
    pub const out: Buf = .b;
    pub fn flops(n: f64) f64 {
        return n * n;
    }
    pub fn Sig(comptime p: Prec) type {
        return fn (Ch, Ch, Ch, In, [*]p.F, In, [*]p.F, In, usize, usize, usize) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        f(&fc.U, &fc.N, &fc.N, &n, w.a.ptr, &n, w.b.ptr, &one_c, 1, 1, 1);
    }
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        const n = w.in();
        try blas.trsv(rstd_alloc, libs.rstd_io, Num(p), .u, .n, .n, libs.matrix(Num(p), cast(p, w.a.ptr), n, n, n), @intCast(n), libs.view(Num(p), cast(p, w.b.ptr), n, 1));
    }
    /// `A x` against the right-hand side.
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, _: usize) ?[2]f64 {
        var s: Pair = .{ 0, 0 };
        var mag: f64 = 0;
        for (i..w.n) |k| {
            const a = w.at(.a0, i + k * w.n);
            const x = w.at(.b, k);
            s = cadd(s, cmul(a, x));
            mag += cmag(a) * cmag(x);
        }
        const r = w.at(.b0, i);
        return .{ cdist(s, r), mag + cmag(r) };
    }
};

pub const Gbmv = struct {
    pub const name = "gbmv";
    pub const shape: Shape = .band;
    pub const out: Buf = .c;
    pub fn flops(n: f64) f64 {
        return 2 * n * @as(f64, @floatFromInt(bandRows(@This())));
    }
    pub fn Sig(comptime p: Prec) type {
        return fn (Ch, In, In, In, In, *const Sc(p), [*]p.F, In, [*]p.F, In, *const Sc(p), [*]p.F, In, usize) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        const lda: c_int = @intCast(bandRows(@This()));
        const al = sc(p, gemv_alpha);
        const one = sc(p, 1);
        f(&fc.N, &n, &n, &band_kc, &band_kc, &al, w.a.ptr, &lda, w.b.ptr, &one_c, &one, w.c.ptr, &one_c, 1);
    }
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        const n = w.in();
        const rows: isize = @intCast(bandRows(@This()));
        try blas.gbmv(rstd_alloc, libs.rstd_io, Num(p), .n, band_k, band_k, num(p, gemv_alpha), libs.matrix(Num(p), cast(p, w.a.ptr), rows, n, rows), @intCast(rows), libs.view(Num(p), cast(p, w.b.ptr), n, 1).toConst(), num(p, 1), libs.view(Num(p), cast(p, w.c.ptr), n, 1));
    }
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, _: usize) ?[2]f64 {
        const lda = bandRows(@This());
        var s: Pair = .{ 0, 0 };
        var mag: f64 = 0;
        for (i -| band_k..@min(w.n, i + band_k + 1)) |k| {
            const a = w.at(.a0, band_k + i - k + k * lda);
            const x = w.at(.b0, k);
            s = cadd(s, cmul(a, x));
            mag += cmag(a) * cmag(x);
        }
        const al = cnum(p, gemv_alpha);
        const y0 = w.at(.c0, i);
        return .{ cdist(cadd(y0, cmul(al, s)), w.at(.c, i)), cmag(y0) + cmag(al) * mag };
    }
};
