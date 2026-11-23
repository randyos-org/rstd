//! The LAPACK entries.

const rstd = @import("rstd");
const bb = @import("../blas_bench.zig");
const libs = @import("../libs.zig");
const xr = @import("../ext.zig");
const pr = @import("pair.zig");
const wk = @import("work.zig");
const rc = @import("race.zig");

const lapack = rstd.math.kernels.lapack;
const Prec = bb.Prec;
const cast = libs.cast;
const Num = libs.Num;
const rstd_alloc = libs.rstd_alloc;
const Buf = xr.Buf;
const Cols = xr.Cols;
const Fill = xr.Fill;
const Level = xr.Level;
const Prep = xr.Prep;
const Restore = xr.Restore;
const Shape = xr.Shape;
const band_k = xr.band_k;
const Ch = pr.Ch;
const In = pr.In;
const Pair = pr.Pair;
const band_kc = pr.band_kc;
const cadd = pr.cadd;
const cconj = pr.cconj;
const cdist = pr.cdist;
const cmag = pr.cmag;
const cmul = pr.cmul;
const fc = pr.fc;
const one_c = pr.one_c;
const Work = wk.Work;
const symAt = rc.symAt;

/// `orig X` at (i, j), over k in `lo..hi`, against the identity; `X` read
/// by `x(w, k, j)`: the residual of an inverse.
pub fn inverseAt(comptime p: Prec, w: *const Work(p), i: usize, j: usize, lo: usize, hi: usize, comptime x: fn (*const Work(p), usize, usize) Pair) [2]f64 {
    var s: Pair = .{ 0, 0 };
    var mag: f64 = 0;
    for (lo..hi) |k| {
        const a = w.at(.orig, i + k * w.n);
        const v = x(w, k, j);
        s = cadd(s, cmul(a, v));
        mag += cmag(a) * cmag(v);
    }
    const id: f64 = if (i == j) 1 else 0;
    return .{ cdist(s, .{ id, 0 }), mag + id };
}

/// `orig X` at (i, j) against the right-hand sides `b0`; the solution is
/// in `b`.
pub fn solveAt(comptime p: Prec, w: *const Work(p), i: usize, j: usize) [2]f64 {
    const n = w.n;
    var s: Pair = .{ 0, 0 };
    var mag: f64 = 0;
    for (0..n) |k| {
        const a = w.at(.orig, i + k * n);
        const x = w.at(.b, k + j * n);
        s = cadd(s, cmul(a, x));
        mag += cmag(a) * cmag(x);
    }
    const r = w.at(.b0, i + j * n);
    return .{ cdist(s, r), mag + cmag(r) };
}

/// `uplo, n, a, lda, info`: `potrf`, `potri`.
pub fn SigUplo(comptime p: Prec) type {
    return fn (Ch, In, [*]p.F, In, *c_int, usize) callconv(.c) void;
}

pub fn libUplo(comptime p: Prec, f: *const SigUplo(p), w: *Work(p)) void {
    const n = w.ni();
    f(&fc.U, &n, w.a.ptr, &n, &w.info, 1);
}

pub const Potrf = struct {
    pub const name = "potrf";
    pub const shape: Shape = .mm;
    pub const fill: Fill = .spd;
    pub const restore: Restore = .{ .a = true };
    pub const out: Buf = .a;
    pub fn flops(n: f64) f64 {
        return n * n * n / 3;
    }
    pub const Sig = SigUplo;
    pub const lib = libUplo;
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        const n = w.in();
        try lapack.potrf(rstd_alloc, libs.rstd_io, Num(p), .u, w.n, libs.matrix(Num(p), cast(p, w.a.ptr), n, n, n), w.n);
    }
    /// `Uᴴ U` against `A`, upper triangle.
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        if (i > j) return null;
        const n = w.n;
        var s: Pair = .{ 0, 0 };
        var mag: f64 = 0;
        for (0..i + 1) |k| {
            const x = cconj(w.at(.a, k + i * n));
            const y = w.at(.a, k + j * n);
            s = cadd(s, cmul(x, y));
            mag += cmag(x) * cmag(y);
        }
        return .{ cdist(s, w.at(.orig, i + j * n)), mag };
    }
};

pub const Getrf = struct {
    pub const name = "getrf";
    pub const shape: Shape = .mm;
    pub const fill: Fill = .dominant;
    pub const restore: Restore = .{ .a = true };
    pub const out: Buf = .a;
    pub fn flops(n: f64) f64 {
        return 2 * n * n * n / 3;
    }
    pub fn Sig(comptime p: Prec) type {
        return fn (In, In, [*]p.F, In, [*]c_int, *c_int) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        f(&n, &n, w.a.ptr, &n, w.ipiv.ptr, &w.info);
    }
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        const n = w.in();
        _ = try lapack.getrf(rstd_alloc, libs.rstd_io, Num(p), w.n, w.n, libs.matrix(Num(p), cast(p, w.a.ptr), n, n, n), w.n, w.rpiv);
        // the libraries' 1-based form, which `residual` reads
        for (w.ipiv[0..w.n], w.rpiv[0..w.n]) |*d, s| d.* = @intCast(s + 1);
    }
    /// `L U` against `P A`.
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        const n = w.n;
        // the row of `A` the exchanges, made in order, bring to row i
        var row = i;
        var k = n;
        while (k > 0) {
            k -= 1;
            const q: usize = @intCast(w.ipiv[k] - 1);
            if (row == k) row = q else if (row == q) row = k;
        }
        var s: Pair = .{ 0, 0 };
        var mag: f64 = 0;
        for (0..@min(i, j) + 1) |l| {
            const x: Pair = if (l == i) .{ 1, 0 } else w.at(.a, i + l * n);
            const y = w.at(.a, l + j * n);
            s = cadd(s, cmul(x, y));
            mag += cmag(x) * cmag(y);
        }
        return .{ cdist(s, w.at(.orig, row + j * n)), mag };
    }
};

pub const Getrs = struct {
    pub const name = "getrs";
    pub const shape: Shape = .mm;
    pub const fill: Fill = .dominant;
    pub const prep: Prep = .getrf;
    pub const restore: Restore = .{ .b = true };
    pub const out: Buf = .b;
    pub fn flops(n: f64) f64 {
        return 2 * n * n * n;
    }
    pub fn Sig(comptime p: Prec) type {
        return fn (Ch, In, In, [*]p.F, In, [*]c_int, [*]p.F, In, *c_int, usize) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        f(&fc.N, &n, &n, w.a.ptr, &n, w.ipiv.ptr, w.b.ptr, &n, &w.info, 1);
    }
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        const n = w.in();
        try lapack.getrs(rstd_alloc, libs.rstd_io, Num(p), .n, w.n, w.n, libs.matrix(Num(p), cast(p, w.a.ptr), n, n, n), w.n, w.rpiv, libs.matrix(Num(p), cast(p, w.b.ptr), n, n, n), w.n);
    }
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        return solveAt(p, w, i, j);
    }
};

pub const Gesv = struct {
    pub const name = "gesv";
    pub const shape: Shape = .mm;
    pub const fill: Fill = .dominant;
    pub const restore: Restore = .{ .a = true, .b = true };
    pub const out: Buf = .b;
    pub const cols: Cols = .one;
    pub fn flops(n: f64) f64 {
        return 2 * n * n * n / 3 + 2 * n * n;
    }
    pub fn Sig(comptime p: Prec) type {
        return fn (In, In, [*]p.F, In, [*]c_int, [*]p.F, In, *c_int) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        f(&n, &one_c, w.a.ptr, &n, w.ipiv.ptr, w.b.ptr, &n, &w.info);
    }
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        const n = w.in();
        try lapack.gesv(rstd_alloc, libs.rstd_io, Num(p), w.n, 1, libs.matrix(Num(p), cast(p, w.a.ptr), n, n, n), w.n, w.rpiv, libs.matrix(Num(p), cast(p, w.b.ptr), n, 1, n), w.n);
    }
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        return solveAt(p, w, i, j);
    }
};

pub const Trtri = struct {
    pub const name = "trtri";
    pub const shape: Shape = .mm;
    pub const fill: Fill = .dominant;
    pub const restore: Restore = .{ .a = true };
    pub const out: Buf = .a;
    pub fn flops(n: f64) f64 {
        return n * n * n / 3;
    }
    pub fn Sig(comptime p: Prec) type {
        return fn (Ch, Ch, In, [*]p.F, In, *c_int, usize, usize) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        f(&fc.U, &fc.N, &n, w.a.ptr, &n, &w.info, 1, 1);
    }
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        const n = w.in();
        try lapack.trtri(rstd_alloc, libs.rstd_io, Num(p), .u, .n, w.n, libs.matrix(Num(p), cast(p, w.a.ptr), n, n, n), w.n);
    }
    /// `A X` against `I`, upper triangle.
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        if (i > j) return null;
        return inverseAt(p, w, i, j, i, j + 1, struct {
            fn x(v: *const Work(p), k: usize, c: usize) Pair {
                return v.at(.a, k + c * v.n);
            }
        }.x);
    }
};

pub const Potri = struct {
    pub const name = "potri";
    pub const shape: Shape = .mm;
    pub const fill: Fill = .spd;
    pub const prep: Prep = .potrf;
    pub const restore: Restore = .{ .a = true };
    pub const out: Buf = .a;
    pub fn flops(n: f64) f64 {
        return 2 * n * n * n / 3;
    }
    pub const Sig = SigUplo;
    pub const lib = libUplo;
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        const n = w.in();
        try lapack.potri(rstd_alloc, libs.rstd_io, Num(p), .u, w.n, libs.matrix(Num(p), cast(p, w.a.ptr), n, n, n), w.n);
    }
    /// `A X` against `I`, `X` Hermitian from its upper triangle.
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        return inverseAt(p, w, i, j, 0, w.n, struct {
            fn x(v: *const Work(p), k: usize, c: usize) Pair {
                return symAt(p, v, .a, true, k, c);
            }
        }.x);
    }
};

/// `orig` against `U diag(s) Vᵀ`.
pub fn svdAt(comptime p: Prec, w: *const Work(p), i: usize, j: usize) [2]f64 {
    const n = w.n;
    var s: Pair = .{ 0, 0 };
    var mag: f64 = 0;
    for (0..n) |k| {
        const u = w.at(.u, i + k * n);
        const v = w.at(.vt, k + j * n);
        const sk = w.at(.s, k)[0];
        const x = cmul(u, v);
        s = cadd(s, .{ x[0] * sk, x[1] * sk });
        mag += cmag(u) * @abs(sk) * cmag(v);
    }
    return .{ cdist(s, w.at(.orig, i + j * n)), mag };
}

/// faer: thin `U` and `Vᵀ` (`S`).
pub const Gesvd = struct {
    pub const name = "gesvd";
    pub const shape: Shape = .mm;
    pub const restore: Restore = .{ .a = true };
    pub const out: Buf = .s;
    pub fn flops(n: f64) f64 {
        return 22 * n * n * n;
    }
    pub fn Sig(comptime p: Prec) type {
        const F = p.F;
        return if (p.complex)
            fn (Ch, Ch, In, In, [*]F, In, [*]F, [*]F, In, [*]F, In, [*]F, In, [*]F, *c_int, usize, usize) callconv(.c) void
        else
            fn (Ch, Ch, In, In, [*]F, In, [*]F, [*]F, In, [*]F, In, [*]F, In, *c_int, usize, usize) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        const lw = w.lw();
        if (p.complex)
            f(&fc.S, &fc.S, &n, &n, w.a.ptr, &n, w.s.ptr, w.u.ptr, &n, w.vt.ptr, &n, w.work.ptr, &lw, w.rwork.ptr, &w.info, 1, 1)
        else
            f(&fc.S, &fc.S, &n, &n, w.a.ptr, &n, w.s.ptr, w.u.ptr, &n, w.vt.ptr, &n, w.work.ptr, &lw, &w.info, 1, 1);
    }
    pub const has_rstd = @hasDecl(lapack, "gesvd");
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        const n = w.in();
        try lapack.gesvd(rstd_alloc, libs.rstd_io, Num(p), .leading, .leading, w.n, w.n, libs.matrix(Num(p), cast(p, w.a.ptr), n, n, n), w.n, w.s[0..w.n], libs.matrix(Num(p), cast(p, w.u.ptr), n, n, n), w.n, libs.matrix(Num(p), cast(p, w.vt.ptr), n, n, n), w.n, null);
    }
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        return svdAt(p, w, i, j);
    }
};

/// faer: `jobz` `S`.
pub const Gesdd = struct {
    pub const name = "gesdd";
    pub const shape: Shape = .mm;
    pub const restore: Restore = .{ .a = true };
    pub const out: Buf = .s;
    pub fn flops(n: f64) f64 {
        return 22 * n * n * n;
    }
    pub fn Sig(comptime p: Prec) type {
        const F = p.F;
        return if (p.complex)
            fn (Ch, In, In, [*]F, In, [*]F, [*]F, In, [*]F, In, [*]F, In, [*]F, [*]c_int, *c_int, usize) callconv(.c) void
        else
            fn (Ch, In, In, [*]F, In, [*]F, [*]F, In, [*]F, In, [*]F, In, [*]c_int, *c_int, usize) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        const lw = w.lw();
        if (p.complex)
            f(&fc.S, &n, &n, w.a.ptr, &n, w.s.ptr, w.u.ptr, &n, w.vt.ptr, &n, w.work.ptr, &lw, w.rwork.ptr, w.iwork.ptr, &w.info, 1)
        else
            f(&fc.S, &n, &n, w.a.ptr, &n, w.s.ptr, w.u.ptr, &n, w.vt.ptr, &n, w.work.ptr, &lw, w.iwork.ptr, &w.info, 1);
    }
    pub const has_rstd = @hasDecl(lapack, "gesdd");
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        const n = w.in();
        try lapack.gesdd(rstd_alloc, libs.rstd_io, Num(p), .leading, w.n, w.n, libs.matrix(Num(p), cast(p, w.a.ptr), n, n, n), w.n, w.s[0..w.n], libs.matrix(Num(p), cast(p, w.u.ptr), n, n, n), w.n, libs.matrix(Num(p), cast(p, w.vt.ptr), n, n, n), w.n);
    }
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        return svdAt(p, w, i, j);
    }
};

/// faer's QR; rstd's call waits on issue #120. `tau` in `s`.
pub const Geqrf = struct {
    pub const name = "geqrf";
    pub const shape: Shape = .mm;
    pub const restore: Restore = .{ .a = true };
    pub const out: Buf = .a;
    pub fn flops(n: f64) f64 {
        return 4 * n * n * n / 3;
    }
    pub fn Sig(comptime p: Prec) type {
        return fn (In, In, [*]p.F, In, [*]p.F, [*]p.F, In, *c_int) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        const lw = w.lw();
        f(&n, &n, w.a.ptr, &n, w.s.ptr, w.work.ptr, &lw, &w.info);
    }
};

/// `syevd`, `heevd` complex, with vectors; rstd's call waits on issue #120.
pub const Syevd = struct {
    pub const name = "syevd";
    pub const shape: Shape = .mm;
    pub const fill: Fill = .spd;
    pub const restore: Restore = .{ .a = true };
    pub const out: Buf = .s;
    pub fn flops(n: f64) f64 {
        return 9 * n * n * n;
    }
    pub fn sym(comptime p: Prec) []const u8 {
        return p.name ++ (if (p.complex) "heevd" else "syevd");
    }
    pub fn Sig(comptime p: Prec) type {
        const F = p.F;
        return if (p.complex)
            fn (Ch, Ch, In, [*]F, In, [*]F, [*]F, In, [*]F, In, [*]c_int, In, *c_int, usize, usize) callconv(.c) void
        else
            fn (Ch, Ch, In, [*]F, In, [*]F, [*]F, In, [*]c_int, In, *c_int, usize, usize) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        const lw = w.lw();
        const liw = w.liw();
        if (p.complex) {
            const lrw = w.lrw();
            f(&fc.V, &fc.U, &n, w.a.ptr, &n, w.s.ptr, w.work.ptr, &lw, w.rwork.ptr, &lrw, w.iwork.ptr, &liw, &w.info, 1, 1);
        } else f(&fc.V, &fc.U, &n, w.a.ptr, &n, w.s.ptr, w.work.ptr, &lw, w.iwork.ptr, &liw, &w.info, 1, 1);
    }
    /// `A Z` against `Z diag(w)`.
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        const n = w.n;
        var s: Pair = .{ 0, 0 };
        var mag: f64 = 0;
        for (0..n) |k| {
            const a = w.at(.orig, i + k * n);
            const z = w.at(.a, k + j * n);
            s = cadd(s, cmul(a, z));
            mag += cmag(a) * cmag(z);
        }
        const ev = w.at(.s, j)[0];
        const z = w.at(.a, i + j * n);
        return .{ cdist(s, .{ z[0] * ev, z[1] * ev }), mag + @abs(ev) * cmag(z) };
    }
};

/// `syevd`, `heevd` with `jobz` `N`: the values alone, the reduction and
/// the tridiagonal's eigenvalues (rstd's `stebr.zig` from
/// `evd.stebr_crossover`); a nominal 4n³/3, the reduction's. No residual
/// without vectors: only the agreement with rstd.
pub const SyevdN = struct {
    pub const name = "syevdN";
    pub const shape: Shape = .mm;
    pub const fill: Fill = .spd;
    pub const restore: Restore = .{ .a = true };
    pub const out: Buf = .s;
    pub fn flops(n: f64) f64 {
        return 4 * n * n * n / 3;
    }
    pub fn sym(comptime p: Prec) []const u8 {
        return p.name ++ (if (p.complex) "heevd" else "syevd");
    }
    pub const Sig = Syevd.Sig;
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        const lw = w.lw();
        const liw = w.liw();
        if (p.complex) {
            const lrw = w.lrw();
            f(&fc.N, &fc.U, &n, w.a.ptr, &n, w.s.ptr, w.work.ptr, &lw, w.rwork.ptr, &lrw, w.iwork.ptr, &liw, &w.info, 1, 1);
        } else f(&fc.N, &fc.U, &n, w.a.ptr, &n, w.s.ptr, w.work.ptr, &lw, w.iwork.ptr, &liw, &w.info, 1, 1);
    }
    pub const has_rstd = @hasDecl(lapack, "syevd");
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        const n = w.in();
        try lapack.syevd(rstd_alloc, libs.rstd_io, Num(p), .none, .u, w.n, libs.matrix(Num(p), cast(p, w.a.ptr), n, n, n), w.n, w.s[0..w.n]);
    }
};

/// ReLAPACK: Bunch-Kaufman, symmetric (not Hermitian) complex too;
/// no residual, which would need the pivots' blocks.
pub const Sytrf = struct {
    pub const name = "sytrf";
    pub const shape: Shape = .mm;
    pub const fill: Fill = .spd;
    pub const restore: Restore = .{ .a = true };
    pub const out: Buf = .a;
    pub fn flops(n: f64) f64 {
        return n * n * n / 3;
    }
    pub fn Sig(comptime p: Prec) type {
        return fn (Ch, In, [*]p.F, In, [*]c_int, [*]p.F, In, *c_int, usize) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        const lw = w.lw();
        f(&fc.U, &n, w.a.ptr, &n, w.ipiv.ptr, w.work.ptr, &lw, &w.info, 1);
    }
    pub const has_rstd = @hasDecl(lapack, "sytrf");
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        const n = w.in();
        _ = try lapack.sytrf(rstd_alloc, libs.rstd_io, Num(p), .u, w.n, libs.matrix(Num(p), cast(p, w.a.ptr), n, n, n), w.n, w.spiv);
    }
};

/// ReLAPACK: kl = ku = `band_k` in LAPACK's band storage, kl rows
/// for the fill, at `gemv`'s band sizes; no residual.
pub const Gbtrf = struct {
    pub const name = "gbtrf";
    pub const level: Level = .l2;
    pub const shape: Shape = .band;
    pub const band_rows: usize = 3 * band_k + 1;
    pub const restore: Restore = .{ .a = true };
    pub const out: Buf = .a;
    pub fn flops(n: f64) f64 {
        return 4 * n * band_k * band_k;
    }
    pub fn Sig(comptime p: Prec) type {
        return fn (In, In, In, In, [*]p.F, In, [*]c_int, *c_int) callconv(.c) void;
    }
    pub fn lib(comptime p: Prec, f: *const Sig(p), w: *Work(p)) void {
        const n = w.ni();
        const ld: c_int = band_rows;
        f(&n, &n, &band_kc, &band_kc, w.a.ptr, &ld, w.ipiv.ptr, &w.info);
    }
    pub const has_rstd = @hasDecl(lapack, "gbtrf");
    pub fn rstd(comptime p: Prec, w: *Work(p)) !void {
        _ = try lapack.gbtrf(Num(p), w.n, w.n, band_k, band_k, libs.matrix(Num(p), cast(p, w.a.ptr), band_rows, w.in(), band_rows), band_rows, w.rpiv);
    }
};
