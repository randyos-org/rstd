//! The level-2 programs `test/[sdcz]blat2.f` with their `.dat` tables and
//! their CBLAS `ctest/` versions, OpenBLAS v0.3.34-317; `gbmv`'s negative `kl`
//! or `ku` is `InvalidK`, and an error on a valid call fails at once.

const std = @import("std");
const ob = @import("__root__.zig");
const blas = @import("../../__root__.zig");
const check = @import("../check.zig");
const guard = @import("../guard.zig");

const Float = ob.Float;
const isComplex = ob.isComplex;
const testing = std.testing;

/// `NMAX` and `INCMAX`: the largest order and stride the programs take.
const nmax = 65;
const incmax = 2;

/// A `.dat` file's sweep: values of `N`, `K`, `INCX` and `INCY`, `ALPHA`
/// and `BETA` (complex as re, im), and the threshold of the test ratio.
pub const Params = struct {
    thresh: f64 = 16,
    n: []const usize,
    k: []const usize = &.{ 0, 1, 2, 4 },
    inc: []const isize = &.{ 1, 2, -1, -2 },
    alpha: []const [2]f64,
    beta: []const [2]f64,
};

/// `test/[sd]blat2.dat`.
pub const real_dat: Params = .{
    .n = &.{ 0, 1, 2, 3, 7, 31, 63 },
    .alpha = &.{ .{ 0, 0 }, .{ 1, 0 }, .{ 0.7, 0 } },
    .beta = &.{ .{ 0, 0 }, .{ 1, 0 }, .{ 0.9, 0 } },
};

/// `test/[cz]blat2.dat`.
pub const complex_dat: Params = .{
    .n = &.{ 0, 1, 2, 3, 7, 31, 63 },
    .alpha = &.{ .{ 0, 0 }, .{ 1, 0 }, .{ 0.7, -0.9 } },
    .beta = &.{ .{ 0, 0 }, .{ 1, 0 }, .{ 1.3, -1.1 } },
};

/// `ctest/[sdcz]in2`: the same but for `N`.
const ctest_n = [_]usize{ 0, 1, 2, 3, 5, 9, 63 };

/// The rigs' `N`: the example data file in each program's header.
pub const rig_n = [_]usize{ 0, 1, 2, 3, 5, 9 };

fn dat(comptime T: type) Params {
    return if (comptime isComplex(T)) complex_dat else real_dat;
}

fn withN(p: Params, n: []const usize) Params {
    var q = p;
    q.n = n;
    return q;
}

/// The routines, by their real names; complex runs `HEMV`, `HBMV`, `HPMV`,
/// `GERU`, `HER`, `HPR`, `HER2`, `HPR2` in their place.
const Routine = enum { gemv, gbmv, symv, sbmv, spmv, trmv, tbmv, tpmv, trsv, tbsv, tpsv, ger, gerc, syr, spr, syr2, spr2 };

fn name(comptime T: type, r: Routine) []const u8 {
    const c = comptime isComplex(T);
    const s = switch (r) {
        .symv => if (c) "HEMV" else "SYMV",
        .sbmv => if (c) "HBMV" else "SBMV",
        .spmv => if (c) "HPMV" else "SPMV",
        .ger => if (c) "GERU" else "GER",
        .syr => if (c) "HER" else "SYR",
        .spr => if (c) "HPR" else "SPR",
        .syr2 => if (c) "HER2" else "SYR2",
        .spr2 => if (c) "HPR2" else "SPR2",
        inline else => |t| comptime upper(@tagName(t)),
    };
    return s;
}

fn upper(comptime s: []const u8) []const u8 {
    var out: [s.len]u8 = undefined;
    for (s, 0..) |ch, i| out[i] = std.ascii.toUpper(ch);
    const final = out;
    return &final;
}

fn letter(comptime T: type) u8 {
    return switch (Float(T)) {
        f32 => if (isComplex(T)) 'C' else 'S',
        else => if (isComplex(T)) 'Z' else 'D',
    };
}

// ---- SMAKE, LSERES -----------------------------------------------------

const Kind = enum { ge, gb, sy, sb, sp, tr, tb, tp };

fn rogue(comptime T: type) T {
    return ob.val(T, -1e10, 1e10);
}

/// `SMAKE`/`CMAKE`: an `m` by `n` matrix within the band `kl`, `ku`, into
/// `a` (leading dimension `ld`) and stored as `kind` in `aa`, unused
/// elements rogue. `sy`, `sb`, `sp` are Hermitian for complex `T`.
fn make(comptime T: type, kind: Kind, uplo: u8, diag: u8, m: usize, n: usize, a: []T, ld: usize, aa: []T, lda: usize, kl: isize, ku: isize, gen: *ob.Gen, transl: f64) void {
    const gen_ = kind == .ge or kind == .gb;
    const sym = kind == .sy or kind == .sb or kind == .sp;
    const tri_ = kind == .tr or kind == .tb or kind == .tp;
    const up = (sym or tri_) and uplo == 'U';
    const lo = (sym or tri_) and uplo == 'L';
    const unit = tri_ and diag == 'U';
    const herm = sym and comptime isComplex(T);
    const z = ob.zero(T);
    const im: isize = @intCast(m);
    const in: isize = @intCast(n);
    const A = struct {
        fn at(buf: []T, l: usize, i: isize, j: isize) *T {
            return &buf[@as(usize, @intCast(i - 1)) + @as(usize, @intCast(j - 1)) * l];
        }
    }.at;
    var j: isize = 1;
    while (j <= in) : (j += 1) {
        var i: isize = 1;
        while (i <= im) : (i += 1) {
            if (gen_ or (up and i <= j) or (lo and i >= j)) {
                A(a, ld, i, j).* = if ((i <= j and j - i <= ku) or (i >= j and i - j <= kl))
                    ob.add(gen.next(T), ob.val(T, transl, 0))
                else
                    z;
                if (i != j) {
                    if (sym) {
                        A(a, ld, j, i).* = ob.conj(A(a, ld, i, j).*);
                    } else if (tri_) A(a, ld, j, i).* = z;
                }
            }
        }
        if (herm) A(a, ld, j, j).* = ob.val(T, re(A(a, ld, j, j).*), 0);
        if (tri_) A(a, ld, j, j).* = ob.add(A(a, ld, j, j).*, ob.val(T, 1, 0));
        if (unit) A(a, ld, j, j).* = ob.val(T, 1, 0);
    }
    const ilda: isize = @intCast(lda);
    const AA = struct {
        fn at(buf: []T, k: isize) *T {
            return &buf[@intCast(k - 1)];
        }
    }.at;
    switch (kind) {
        .ge => {
            j = 1;
            while (j <= in) : (j += 1) {
                var i: isize = 1;
                while (i <= ilda) : (i += 1) AA(aa, i + (j - 1) * ilda).* = if (i <= im) A(a, ld, i, j).* else rogue(T);
            }
        },
        .gb => {
            j = 1;
            while (j <= in) : (j += 1) {
                const first = @max(1, ku + 2 - j);
                const last = @min(kl + ku + 1, ku + 1 + im - j);
                var i: isize = 1;
                while (i <= ilda) : (i += 1) {
                    AA(aa, i + (j - 1) * ilda).* = if (i >= first and i <= last) A(a, ld, i + j - ku - 1, j).* else rogue(T);
                }
            }
        },
        .sy, .tr => {
            j = 1;
            while (j <= in) : (j += 1) {
                const ibeg: isize = if (up) 1 else if (unit) j + 1 else j;
                const iend: isize = if (up) (if (unit) j - 1 else j) else in;
                var i: isize = 1;
                while (i <= ilda) : (i += 1) {
                    AA(aa, i + (j - 1) * ilda).* = if (i >= ibeg and i <= iend) A(a, ld, i, j).* else rogue(T);
                }
                if (herm) setIm(T, AA(aa, j + (j - 1) * ilda), -1e10);
            }
        },
        .sb, .tb => {
            j = 1;
            while (j <= in) : (j += 1) {
                const kk: isize = if (up) kl + 1 else 1;
                const ibeg: isize = if (up) @max(1, kl + 2 - j) else if (unit) 2 else 1;
                const iend: isize = if (up) (if (unit) kl else kl + 1) else @min(kl + 1, 1 + im - j);
                var i: isize = 1;
                while (i <= ilda) : (i += 1) {
                    AA(aa, i + (j - 1) * ilda).* = if (i >= ibeg and i <= iend) A(a, ld, i + j - kk, j).* else rogue(T);
                }
                if (herm) setIm(T, AA(aa, kk + (j - 1) * ilda), -1e10);
            }
        },
        .sp, .tp => {
            var ioff: isize = 0;
            j = 1;
            while (j <= in) : (j += 1) {
                const ibeg: isize = if (up) 1 else j;
                const iend: isize = if (up) j else in;
                var i = ibeg;
                while (i <= iend) : (i += 1) {
                    ioff += 1;
                    AA(aa, ioff).* = A(a, ld, i, j).*;
                    if (i == j) {
                        if (unit) AA(aa, ioff).* = rogue(T);
                        if (herm) setIm(T, AA(aa, ioff), -1e10);
                    }
                }
            }
        },
    }
}

fn re(v: anytype) Float(@TypeOf(v)) {
    return if (comptime isComplex(@TypeOf(v))) v.re else v;
}

fn setIm(comptime T: type, p: *T, v: Float(T)) void {
    if (comptime isComplex(T)) p.im = v;
}

/// `SMAKE('GE', ' ', ' ', 1, n, x, 1, xx, inc, ...)`: a vector, and its
/// storage at stride `inc` with rogue gaps.
fn makeVec(comptime T: type, n: usize, x: []T, xx: []T, inc: usize, gen: *ob.Gen, transl: f64) void {
    make(T, .ge, ' ', ' ', 1, n, x, 1, xx, inc, 0, @as(isize, @intCast(n)) - 1, gen, transl);
}

/// `makeVec`, and the middle element zero in both, as the Fortran sets it
/// for `n > 1`.
fn makeVecZeroMid(comptime T: type, n: usize, x: []T, xx: []T, inc: usize, gen: *ob.Gen, transl: f64) void {
    makeVec(T, n, x, xx, inc, gen, transl);
    if (n > 1) {
        x[n / 2 - 1] = ob.zero(T);
        xx[inc * (n / 2 - 1)] = ob.zero(T);
    }
}

/// `LSERES`: the elements of the `m` by `n` storage outside the result agree,
/// rows past `m` for `ge`, the other triangle for `sy`.
fn lseres(comptime T: type, kind: Kind, uplo: u8, m: usize, n: usize, aa: []const T, as: []const T, lda: usize) bool {
    for (0..n) |j| for (0..lda) |i| {
        const inside = switch (kind) {
            .ge => i < m,
            .sy => if (uplo == 'U') i <= j else i >= j and i < n,
            else => true,
        };
        if (!inside and !ob.eql(aa[i + j * lda], as[i + j * lda])) return false;
    };
    return true;
}

// ---- SMVCH ----------------------------------------------------------------

const Result = ob.Result;

/// `SMVCH`/`CMVCH`: `alpha*op(A)*x + beta*y` from the dense `a` (leading
/// dimension `ld`) and compacted `x`, `y` into `yt`, gauges into `g`, and the
/// test ratio of `yy` (at stride `|incy|`) against it.
fn mvch(comptime T: type, trans: u8, m: usize, n: usize, alpha: T, a: []const T, ld: usize, x: []const T, incx: isize, beta: T, y: []const T, incy: isize, yt: []ob.Ref(T), g: []Float(ob.Ref(T)), yy: []const T) Result {
    const R = ob.Ref(T);
    const RF = Float(R);
    const tran = trans == 'T' or trans == 'C';
    const ctran = trans == 'C';
    const ml = if (tran) n else m;
    const nl = if (tran) m else n;
    const ra = ob.cast(R, alpha);
    const rb = ob.cast(R, beta);
    for (0..ml) |i| {
        const iy = if (incy < 0) ml - 1 - i else i;
        var acc = ob.zero(R);
        var gauge: RF = 0;
        for (0..nl) |j| {
            const jx = if (incx < 0) nl - 1 - j else j;
            var e = ob.cast(R, if (tran) a[j + i * ld] else a[i + j * ld]);
            if (ctran) e = ob.conj(e);
            const xv = ob.cast(R, x[jx]);
            acc = ob.add(acc, ob.mul(e, xv));
            gauge += ob.abs1(e) * ob.abs1(xv);
        }
        const yv = ob.cast(R, y[iy]);
        yt[iy] = ob.add(ob.mul(ra, acc), ob.mul(rb, yv));
        g[iy] = ob.abs1(ra) * gauge + ob.abs1(rb) * ob.abs1(yv);
    }
    const eps: RF = std.math.floatEps(Float(T));
    var err: RF = 0;
    const ay: usize = @abs(incy);
    for (0..ml) |i| {
        var erri = ob.absv(ob.sub(yt[i], ob.cast(R, yy[i * ay]))) / eps;
        if (g[i] != 0) erri /= g[i];
        err = @max(err, erri);
        if (err * @sqrt(eps) >= 1) {
            ob.say("row {d}: expected {e}, computed {e}\n", .{ i + 1, re(yt[i]), re(yy[i * ay]) });
            return .{ .err = @floatCast(err), .fatal = true };
        }
    }
    return .{ .err = @floatCast(err), .fatal = false };
}

// ---- the calls ----------------------------------------------------------------

fn Work(comptime T: type) type {
    const R = ob.Ref(T);
    return struct {
        a: []T,
        aa: []T,
        as: []T,
        x: []T,
        xx: []T,
        xs: []T,
        y: []T,
        yy: []T,
        ys: []T,
        yt: []R,
        g: []Float(R),
        z: []T,

        fn init(al: std.mem.Allocator) !@This() {
            return .{
                .a = try al.alloc(T, nmax * nmax),
                .aa = try al.alloc(T, nmax * nmax),
                .as = try al.alloc(T, nmax * nmax),
                .x = try al.alloc(T, nmax),
                .xx = try al.alloc(T, nmax * incmax),
                .xs = try al.alloc(T, nmax * incmax),
                .y = try al.alloc(T, nmax),
                .yy = try al.alloc(T, nmax * incmax),
                .ys = try al.alloc(T, nmax * incmax),
                .yt = try al.alloc(R, nmax),
                .g = try al.alloc(Float(R), nmax),
                .z = try al.alloc(T, 2 * nmax),
            };
        }
    };
}

/// One call's arguments; `k` is `kl` too.
const Call = struct {
    trans: u8 = 'N',
    uplo: u8 = 'U',
    diag: u8 = 'N',
    m: usize = 0,
    n: usize = 0,
    k: isize = 0,
    ku: isize = 0,
    lda: usize = 0,
    incx: isize = 1,
    incy: isize = 1,
    alpha: [2]f64 = .{ 0, 0 },
    beta: [2]f64 = .{ 0, 0 },

    fn describe(c: Call, comptime T: type, r: Routine, nc: usize) void {
        ob.say("{c}{s} call {d}: trans={c} uplo={c} diag={c} m={d} n={d} k={d} ku={d} lda={d} incx={d} incy={d} alpha=({d},{d}) beta=({d},{d})\n", .{
            letter(T), name(T, r), nc, c.trans, c.uplo, c.diag, c.m, c.n, c.k, c.ku, c.lda, c.incx, c.incy, c.alpha[0], c.alpha[1], c.beta[0], c.beta[1],
        });
    }
};

const transOf = ob.transOf;
const uploOf = ob.uploOf;

/// Calls routine `r` on the storage in `w`.
fn run(comptime T: type, comptime r: Routine, c: Call, w: anytype) !void {
    var poison: guard.Poison = .{ .child = testing.allocator };
    const al = poison.allocator();
    const cx = comptime isComplex(T);
    const alpha = ob.val(T, c.alpha[0], c.alpha[1]);
    const beta = ob.val(T, c.beta[0], c.beta[1]);
    const m: isize = @intCast(c.m);
    const n: isize = @intCast(c.n);
    const lda: isize = @intCast(c.lda);
    const a = w.aa.ptr;
    const x = w.xx.ptr;
    const y = w.yy.ptr;
    switch (r) {
        .gemv => try blas.gemv(al, ob.io, T, transOf(c.trans), m, n, alpha, a, lda, x, c.incx, beta, y, c.incy),
        .gbmv => try blas.gbmv(al, ob.io, T, transOf(c.trans), m, n, c.k, c.ku, alpha, a, lda, x, c.incx, beta, y, c.incy),
        .symv => if (cx)
            try blas.hemv(al, ob.io, T, uploOf(c.uplo), n, alpha, a, lda, x, c.incx, beta, y, c.incy)
        else
            try blas.symv(al, ob.io, T, uploOf(c.uplo), n, alpha, a, lda, x, c.incx, beta, y, c.incy),
        .sbmv => if (cx)
            try blas.hbmv(al, ob.io, T, uploOf(c.uplo), n, c.k, alpha, a, lda, x, c.incx, beta, y, c.incy)
        else
            try blas.sbmv(al, ob.io, T, uploOf(c.uplo), n, c.k, alpha, a, lda, x, c.incx, beta, y, c.incy),
        .spmv => if (cx)
            try blas.hpmv(al, ob.io, T, uploOf(c.uplo), n, alpha, a, x, c.incx, beta, y, c.incy)
        else
            try blas.spmv(al, ob.io, T, uploOf(c.uplo), n, alpha, a, x, c.incx, beta, y, c.incy),
        .trmv, .tbmv, .tpmv, .trsv, .tbsv, .tpsv => {
            const u = ob.letter(blas.Uplo, c.uplo).?;
            const t = ob.letter(blas.Trans, c.trans).?;
            const d = ob.letter(blas.Diag, c.diag).?;
            switch (r) {
                .trmv => try blas.trmv(al, ob.io, T, u, t, d, n, a, lda, x, c.incx),
                .tbmv => try blas.tbmv(al, ob.io, T, u, t, d, n, c.k, a, lda, x, c.incx),
                .tpmv => try blas.tpmv(al, ob.io, T, u, t, d, n, a, x, c.incx),
                .trsv => try blas.trsv(al, ob.io, T, u, t, d, n, a, lda, x, c.incx),
                .tbsv => try blas.tbsv(al, ob.io, T, u, t, d, n, c.k, a, lda, x, c.incx),
                .tpsv => try blas.tpsv(al, ob.io, T, u, t, d, n, a, x, c.incx),
                else => unreachable,
            }
        },
        .ger => try blas.ger(al, ob.io, T, m, n, alpha, x, c.incx, y, c.incy, a, lda),
        .gerc => if (cx) try blas.gerc(al, ob.io, T, m, n, alpha, x, c.incx, y, c.incy, a, lda) else unreachable,
        .syr => if (cx)
            try blas.her(al, ob.io, T, uploOf(c.uplo), n, alpha.re, x, c.incx, a, lda)
        else
            try blas.syr(al, ob.io, T, uploOf(c.uplo), n, alpha, x, c.incx, a, lda),
        .spr => if (cx)
            try blas.hpr(al, ob.io, T, uploOf(c.uplo), n, alpha.re, x, c.incx, a)
        else
            try blas.spr(al, ob.io, T, uploOf(c.uplo), n, alpha, x, c.incx, a),
        .syr2 => if (cx)
            try blas.her2(al, ob.io, T, uploOf(c.uplo), n, alpha, x, c.incx, y, c.incy, a, lda)
        else
            try blas.syr2(al, ob.io, T, uploOf(c.uplo), n, alpha, x, c.incx, y, c.incy, a, lda),
        .spr2 => if (cx)
            try blas.hpr2(al, ob.io, T, uploOf(c.uplo), n, alpha, x, c.incx, y, c.incy, a)
        else
            try blas.spr2(al, ob.io, T, uploOf(c.uplo), n, alpha, x, c.incx, y, c.incy, a),
    }
}

/// The running state of one routine's sweep: calls made and the largest
/// test ratio.
const Sweep = struct {
    nc: usize = 0,
    errmax: f64 = 0,

    /// Makes the call, failing on an error exit (the Fortran's `.NOT.OK`).
    fn call(s: *Sweep, comptime T: type, comptime r: Routine, c: Call, w: anytype) !void {
        s.nc += 1;
        run(T, r, c, w) catch |e| {
            ob.say("error exit taken on a valid call: ", .{});
            c.describe(T, r, s.nc);
            return e;
        };
    }

    /// Fails if `ok` is false: argument `what` was changed.
    fn same(s: *Sweep, comptime T: type, r: Routine, c: Call, what: []const u8, ok: bool) !void {
        if (ok) return;
        ob.say("{s} was changed incorrectly: ", .{what});
        c.describe(T, r, s.nc);
        return error.TestUnexpectedResult;
    }

    fn result(s: *Sweep, comptime T: type, r: Routine, c: Call, res: Result) !void {
        s.errmax = @max(s.errmax, res.err);
        if (!res.fatal) return;
        ob.say("computed result is less than half accurate: ", .{});
        c.describe(T, r, s.nc);
        return error.TestUnexpectedResult;
    }

    fn finish(s: *Sweep, comptime T: type, r: Routine, p: Params) !void {
        if (s.errmax < p.thresh) return;
        ob.say("{c}{s} completed {d} calls with maximum test ratio {d}\n", .{ letter(T), name(T, r), s.nc, s.errmax });
        return error.TestUnexpectedResult;
    }
};

// ---- SCHK1 to SCHK6 -----------------------------------------------------------

/// `SCHK1`: `gemv` and `gbmv`.
fn chk1(comptime T: type, comptime r: Routine, p: Params) !void {
    const band = r == .gbmv;
    var arena: std.heap.ArenaAllocator = .init(testing.allocator);
    defer arena.deinit();
    const w = try Work(T).init(arena.allocator());
    var gen: ob.Gen = .{};
    var s: Sweep = .{};
    for (p.n) |n| {
        const nd = n / 2 + 1;
        next_m: for (0..2) |im| {
            const m = if (im == 0) n -| nd else @min(n + nd, nmax);
            const nk = if (band) p.k.len else 1;
            for (0..nk) |iku| {
                const ku: isize = if (band) @intCast(p.k[iku]) else @as(isize, @intCast(n)) - 1;
                const kl: isize = if (band) @max(ku - 1, 0) else @as(isize, @intCast(m)) - 1;
                const lda = ob.bump(if (band) @intCast(kl + ku + 1) else m, nmax) orelse continue;
                const laa = lda * n;
                const null_ = n == 0 or m == 0;
                make(T, if (band) .gb else .ge, ' ', ' ', m, n, w.a, nmax, w.aa, lda, kl, ku, &gen, 0);
                for ("NTC") |trans| {
                    const tran = trans != 'N';
                    const ml = if (tran) n else m;
                    const nl = if (tran) m else n;
                    for (p.inc) |incx| {
                        const ax: usize = @abs(incx);
                        const lx = ax * nl;
                        makeVecZeroMid(T, nl, w.x, w.xx, ax, &gen, 0.5);
                        for (p.inc) |incy| {
                            const ay: usize = @abs(incy);
                            const ly = ay * ml;
                            for (p.alpha) |alpha| for (p.beta) |beta| {
                                makeVec(T, ml, w.y, w.yy, ay, &gen, 0);
                                @memcpy(w.as[0..laa], w.aa[0..laa]);
                                @memcpy(w.xs[0..lx], w.xx[0..lx]);
                                @memcpy(w.ys[0..ly], w.yy[0..ly]);
                                const c: Call = .{ .trans = trans, .m = m, .n = n, .k = kl, .ku = ku, .lda = lda, .incx = incx, .incy = incy, .alpha = alpha, .beta = beta };
                                try s.call(T, r, c, w);
                                try s.same(T, r, c, "A", ob.lse(T, w.as[0..laa], w.aa[0..laa]));
                                try s.same(T, r, c, "X", ob.lse(T, w.xs[0..lx], w.xx[0..lx]));
                                try s.same(T, r, c, "Y", if (null_) ob.lse(T, w.ys[0..ly], w.yy[0..ly]) else lseres(T, .ge, ' ', 1, ml, w.yy, w.ys, ay));
                                if (null_) continue :next_m;
                                try s.result(T, r, c, mvch(T, trans, m, n, ob.val(T, alpha[0], alpha[1]), w.a, nmax, w.x, incx, ob.val(T, beta[0], beta[1]), w.y, incy, w.yt, w.g, w.yy));
                            };
                        }
                    }
                }
            }
        }
    }
    try s.finish(T, r, p);
}

/// `SCHK2`: `symv`, `sbmv`, `spmv` (complex: `hemv`, `hbmv`, `hpmv`).
fn chk2(comptime T: type, comptime r: Routine, p: Params) !void {
    const band = r == .sbmv;
    const pack = r == .spmv;
    var arena: std.heap.ArenaAllocator = .init(testing.allocator);
    defer arena.deinit();
    const w = try Work(T).init(arena.allocator());
    var gen: ob.Gen = .{};
    var s: Sweep = .{};
    next_n: for (p.n) |n| {
        const nk = if (band) p.k.len else 1;
        for (0..nk) |ik| {
            const k: isize = if (band) @intCast(p.k[ik]) else @as(isize, @intCast(n)) - 1;
            const lda = ob.bump(if (band) @intCast(k + 1) else n, nmax) orelse continue;
            const laa = if (pack) n * (n + 1) / 2 else lda * n;
            const null_ = n == 0;
            for ("UL") |uplo| {
                make(T, if (band) .sb else if (pack) .sp else .sy, uplo, ' ', n, n, w.a, nmax, w.aa, lda, k, k, &gen, 0);
                for (p.inc) |incx| {
                    const ax: usize = @abs(incx);
                    const lx = ax * n;
                    makeVecZeroMid(T, n, w.x, w.xx, ax, &gen, 0.5);
                    for (p.inc) |incy| {
                        const ay: usize = @abs(incy);
                        const ly = ay * n;
                        for (p.alpha) |alpha| for (p.beta) |beta| {
                            makeVec(T, n, w.y, w.yy, ay, &gen, 0);
                            @memcpy(w.as[0..laa], w.aa[0..laa]);
                            @memcpy(w.xs[0..lx], w.xx[0..lx]);
                            @memcpy(w.ys[0..ly], w.yy[0..ly]);
                            const c: Call = .{ .uplo = uplo, .n = n, .k = k, .lda = lda, .incx = incx, .incy = incy, .alpha = alpha, .beta = beta };
                            try s.call(T, r, c, w);
                            try s.same(T, r, c, "A", ob.lse(T, w.as[0..laa], w.aa[0..laa]));
                            try s.same(T, r, c, "X", ob.lse(T, w.xs[0..lx], w.xx[0..lx]));
                            try s.same(T, r, c, "Y", if (null_) ob.lse(T, w.ys[0..ly], w.yy[0..ly]) else lseres(T, .ge, ' ', 1, n, w.yy, w.ys, ay));
                            if (null_) continue :next_n;
                            try s.result(T, r, c, mvch(T, 'N', n, n, ob.val(T, alpha[0], alpha[1]), w.a, nmax, w.x, incx, ob.val(T, beta[0], beta[1]), w.y, incy, w.yt, w.g, w.yy));
                        };
                    }
                }
            }
        }
    }
    try s.finish(T, r, p);
}

/// `SCHK3`: `trmv`, `tbmv`, `tpmv`, `trsv`, `tbsv`, `tpsv`.
fn chk3(comptime T: type, comptime r: Routine, p: Params) !void {
    const band = r == .tbmv or r == .tbsv;
    const pack = r == .tpmv or r == .tpsv;
    const solve = r == .trsv or r == .tbsv or r == .tpsv;
    var arena: std.heap.ArenaAllocator = .init(testing.allocator);
    defer arena.deinit();
    const w = try Work(T).init(arena.allocator());
    @memset(w.z, ob.zero(T));
    var gen: ob.Gen = .{};
    var s: Sweep = .{};
    const one = ob.val(T, 1, 0);
    next_n: for (p.n) |n| {
        const nk = if (band) p.k.len else 1;
        for (0..nk) |ik| {
            const k: isize = if (band) @intCast(p.k[ik]) else @as(isize, @intCast(n)) - 1;
            const lda = ob.bump(if (band) @intCast(k + 1) else n, nmax) orelse continue;
            const laa = if (pack) n * (n + 1) / 2 else lda * n;
            const null_ = n == 0;
            for ("UL") |uplo| for ("NTC") |trans| for ("UN") |diag| {
                make(T, if (band) .tb else if (pack) .tp else .tr, uplo, diag, n, n, w.a, nmax, w.aa, lda, k, k, &gen, 0);
                for (p.inc) |incx| {
                    const ax: usize = @abs(incx);
                    const lx = ax * n;
                    makeVecZeroMid(T, n, w.x, w.xx, ax, &gen, 0.5);
                    @memcpy(w.as[0..laa], w.aa[0..laa]);
                    @memcpy(w.xs[0..lx], w.xx[0..lx]);
                    const c: Call = .{ .uplo = uplo, .trans = trans, .diag = diag, .n = n, .k = k, .lda = lda, .incx = incx };
                    try s.call(T, r, c, w);
                    try s.same(T, r, c, "A", ob.lse(T, w.as[0..laa], w.aa[0..laa]));
                    try s.same(T, r, c, "X", if (null_) ob.lse(T, w.xs[0..lx], w.xx[0..lx]) else lseres(T, .ge, ' ', 1, n, w.xx, w.xs, ax));
                    if (null_) continue :next_n;
                    if (!solve) {
                        try s.result(T, r, c, mvch(T, trans, n, n, one, w.a, nmax, w.x, incx, ob.zero(T), w.z, incx, w.yt, w.g, w.xx));
                    } else {
                        // the approximation to the original vector
                        for (0..n) |i| {
                            w.z[i] = w.xx[i * ax];
                            w.xx[i * ax] = w.x[i];
                        }
                        try s.result(T, r, c, mvch(T, trans, n, n, one, w.a, nmax, w.z, incx, ob.zero(T), w.x, incx, w.yt, w.g, w.xx));
                    }
                }
            };
        }
    }
    try s.finish(T, r, p);
}

/// `SCHK4`: `ger` (complex: `geru`, `gerc`).
fn chk4(comptime T: type, comptime r: Routine, p: Params) !void {
    var arena: std.heap.ArenaAllocator = .init(testing.allocator);
    defer arena.deinit();
    const w = try Work(T).init(arena.allocator());
    var gen: ob.Gen = .{};
    var s: Sweep = .{};
    const one = ob.val(T, 1, 0);
    for (p.n) |n| {
        const nd = n / 2 + 1;
        next_m: for (0..2) |im| {
            const m = if (im == 0) n -| nd else @min(n + nd, nmax);
            const lda = ob.bump(m, nmax) orelse continue;
            const laa = lda * n;
            const null_ = n == 0 or m == 0;
            for (p.inc) |incx| {
                const ax: usize = @abs(incx);
                const lx = ax * m;
                makeVecZeroMid(T, m, w.x, w.xx, ax, &gen, 0.5);
                for (p.inc) |incy| {
                    const ay: usize = @abs(incy);
                    const ly = ay * n;
                    makeVecZeroMid(T, n, w.y, w.yy, ay, &gen, 0);
                    for (p.alpha) |alpha| {
                        make(T, .ge, ' ', ' ', m, n, w.a, nmax, w.aa, lda, @as(isize, @intCast(m)) - 1, @as(isize, @intCast(n)) - 1, &gen, 0);
                        @memcpy(w.as[0..laa], w.aa[0..laa]);
                        @memcpy(w.xs[0..lx], w.xx[0..lx]);
                        @memcpy(w.ys[0..ly], w.yy[0..ly]);
                        const c: Call = .{ .m = m, .n = n, .lda = lda, .incx = incx, .incy = incy, .alpha = alpha };
                        try s.call(T, r, c, w);
                        try s.same(T, r, c, "X", ob.lse(T, w.xs[0..lx], w.xx[0..lx]));
                        try s.same(T, r, c, "Y", ob.lse(T, w.ys[0..ly], w.yy[0..ly]));
                        try s.same(T, r, c, "A", if (null_) ob.lse(T, w.as[0..laa], w.aa[0..laa]) else lseres(T, .ge, ' ', m, n, w.aa, w.as, lda));
                        if (null_) continue :next_m;
                        const z = w.z[0..m];
                        for (z, 0..) |*e, i| e.* = w.x[if (incx > 0) i else m - 1 - i];
                        for (0..n) |j| {
                            var wv = [1]T{w.y[if (incy > 0) j else n - 1 - j]};
                            if (r == .gerc) wv[0] = ob.conj(wv[0]);
                            try s.result(T, r, c, mvch(T, 'N', m, 1, ob.val(T, alpha[0], alpha[1]), z, nmax, &wv, 1, one, w.a[j * nmax ..], 1, w.yt, w.g, w.aa[j * lda ..]));
                        }
                    }
                }
            }
        }
    }
    try s.finish(T, r, p);
}

/// `SCHK5`: `syr`, `spr` (complex: `her`, `hpr`, with the real part of
/// alpha).
fn chk5(comptime T: type, comptime r: Routine, p: Params) !void {
    const pack = r == .spr;
    var arena: std.heap.ArenaAllocator = .init(testing.allocator);
    defer arena.deinit();
    const w = try Work(T).init(arena.allocator());
    var gen: ob.Gen = .{};
    var s: Sweep = .{};
    const one = ob.val(T, 1, 0);
    next_n: for (p.n) |n| {
        const lda = ob.bump(n, nmax) orelse continue;
        const laa = if (pack) n * (n + 1) / 2 else lda * n;
        for ("UL") |uplo| {
            const up = uplo == 'U';
            for (p.inc) |incx| {
                const ax: usize = @abs(incx);
                const lx = ax * n;
                makeVecZeroMid(T, n, w.x, w.xx, ax, &gen, 0.5);
                for (p.alpha) |alpha_c| {
                    const alpha: [2]f64 = .{ alpha_c[0], 0 };
                    const null_ = n == 0 or alpha[0] == 0;
                    make(T, if (pack) .sp else .sy, uplo, ' ', n, n, w.a, nmax, w.aa, lda, @as(isize, @intCast(n)) - 1, @as(isize, @intCast(n)) - 1, &gen, 0);
                    @memcpy(w.as[0..laa], w.aa[0..laa]);
                    @memcpy(w.xs[0..lx], w.xx[0..lx]);
                    const c: Call = .{ .uplo = uplo, .n = n, .lda = lda, .incx = incx, .alpha = alpha };
                    try s.call(T, r, c, w);
                    try s.same(T, r, c, "X", ob.lse(T, w.xs[0..lx], w.xx[0..lx]));
                    try s.same(T, r, c, "A", if (null_) ob.lse(T, w.as[0..laa], w.aa[0..laa]) else pack or lseres(T, .sy, uplo, n, n, w.aa, w.as, lda));
                    if (null_) {
                        if (n == 0) continue :next_n;
                        continue;
                    }
                    const z = w.z[0..n];
                    for (z, 0..) |*e, i| e.* = w.x[if (incx > 0) i else n - 1 - i];
                    var ja: usize = 0;
                    for (0..n) |j| {
                        const wv = [1]T{ob.conj(z[j])};
                        const jj = if (up) 0 else j;
                        const lj = if (up) j + 1 else n - j;
                        try s.result(T, r, c, mvch(T, 'N', lj, 1, ob.val(T, alpha[0], alpha[1]), z[jj..], lj, &wv, 1, one, w.a[jj + j * nmax ..], 1, w.yt, w.g, w.aa[ja..]));
                        ja += if (pack) lj else if (up) lda else lda + 1;
                    }
                }
            }
        }
    }
    try s.finish(T, r, p);
}

/// `SCHK6`: `syr2`, `spr2` (complex: `her2`, `hpr2`).
fn chk6(comptime T: type, comptime r: Routine, p: Params) !void {
    const pack = r == .spr2;
    var arena: std.heap.ArenaAllocator = .init(testing.allocator);
    defer arena.deinit();
    const w = try Work(T).init(arena.allocator());
    var gen: ob.Gen = .{};
    var s: Sweep = .{};
    const one = ob.val(T, 1, 0);
    next_n: for (p.n) |n| {
        const lda = ob.bump(n, nmax) orelse continue;
        const laa = if (pack) n * (n + 1) / 2 else lda * n;
        for ("UL") |uplo| {
            const up = uplo == 'U';
            for (p.inc) |incx| {
                const ax: usize = @abs(incx);
                const lx = ax * n;
                makeVecZeroMid(T, n, w.x, w.xx, ax, &gen, 0.5);
                for (p.inc) |incy| {
                    const ay: usize = @abs(incy);
                    const ly = ay * n;
                    makeVecZeroMid(T, n, w.y, w.yy, ay, &gen, 0);
                    for (p.alpha) |alpha| {
                        const av = ob.val(T, alpha[0], alpha[1]);
                        const null_ = n == 0 or ob.isZero(av);
                        make(T, if (pack) .sp else .sy, uplo, ' ', n, n, w.a, nmax, w.aa, lda, @as(isize, @intCast(n)) - 1, @as(isize, @intCast(n)) - 1, &gen, 0);
                        @memcpy(w.as[0..laa], w.aa[0..laa]);
                        @memcpy(w.xs[0..lx], w.xx[0..lx]);
                        @memcpy(w.ys[0..ly], w.yy[0..ly]);
                        const c: Call = .{ .uplo = uplo, .n = n, .lda = lda, .incx = incx, .incy = incy, .alpha = alpha };
                        try s.call(T, r, c, w);
                        try s.same(T, r, c, "X", ob.lse(T, w.xs[0..lx], w.xx[0..lx]));
                        try s.same(T, r, c, "Y", ob.lse(T, w.ys[0..ly], w.yy[0..ly]));
                        try s.same(T, r, c, "A", if (null_) ob.lse(T, w.as[0..laa], w.aa[0..laa]) else pack or lseres(T, .sy, uplo, n, n, w.aa, w.as, lda));
                        if (null_) {
                            if (n == 0) continue :next_n;
                            continue;
                        }
                        // Z(:, 1) is x, Z(:, 2) is y, both forwards
                        const z1 = w.z[0..n];
                        const z2 = w.z[nmax..][0..n];
                        for (0..n) |i| {
                            z1[i] = w.x[if (incx > 0) i else n - 1 - i];
                            z2[i] = w.y[if (incy > 0) i else n - 1 - i];
                        }
                        var ja: usize = 0;
                        for (0..n) |j| {
                            const wv: [2]T = if (comptime isComplex(T))
                                .{ ob.mul(av, ob.conj(z2[j])), ob.mul(ob.conj(av), ob.conj(z1[j])) }
                            else
                                .{ z2[j], z1[j] };
                            const ma = if (comptime isComplex(T)) one else av;
                            const jj = if (up) 0 else j;
                            const lj = if (up) j + 1 else n - j;
                            try s.result(T, r, c, mvch(T, 'N', lj, 2, ma, w.z[jj..], nmax, &wv, 1, one, w.a[jj + j * nmax ..], 1, w.yt, w.g, w.aa[ja..]));
                            ja += if (pack) lj else if (up) lda else lda + 1;
                        }
                    }
                }
            }
        }
    }
    try s.finish(T, r, p);
}

fn chk(comptime T: type, comptime r: Routine, p: Params) !void {
    const f = switch (r) {
        .gemv, .gbmv => chk1,
        .symv, .sbmv, .spmv => chk2,
        .trmv, .tbmv, .tpmv, .trsv, .tbsv, .tpsv => chk3,
        .ger, .gerc => chk4,
        .syr, .spr => chk5,
        .syr2, .spr2 => chk6,
    };
    try f(T, r, p);
}

const all_types = check.sd_types;

/// Every type `r` exists for, with the `.dat` table natively, ctest's too,
/// and the header's on the rigs.
fn sweep(comptime r: Routine) !void {
    inline for (all_types) |T| {
        if (r == .gerc and !comptime isComplex(T)) continue;
        if (ob.reduced()) {
            try chk(T, r, withN(dat(T), &rig_n));
        } else {
            try chk(T, r, dat(T));
            try chk(T, r, withN(dat(T), &ctest_n));
        }
    }
}

/// Every routine at `T` on the `N` of `n`, after the error exits, for
/// `check.eachCore`.
pub fn allOf(comptime T: type, n: []const usize) !void {
    try errorExits(T);
    inline for (comptime std.enums.values(Routine)) |r| {
        if (r != .gerc or comptime isComplex(T)) try chk(T, r, withN(dat(T), n));
    }
}

test "blat2: SMVCH self-check on exact data" {
    inline for (all_types) |T| {
        const R = ob.Ref(T);
        const n = 32;
        var a: [n * n]T = undefined;
        var x: [n]T = undefined;
        var y: [n]T = undefined;
        var yy: [n]T = undefined;
        var yt: [n]R = undefined;
        var g: [n]Float(R) = undefined;
        for (0..n) |jj| {
            const j: i64 = @intCast(jj + 1);
            for (0..n) |ii| {
                const i: i64 = @intCast(ii + 1);
                a[ii + jj * n] = ob.val(T, @floatFromInt(@max(i - j + 1, 0)), 0);
            }
            x[jj] = ob.val(T, @floatFromInt(j), 0);
            y[jj] = ob.zero(T);
            yy[jj] = ob.val(T, @floatFromInt(@divTrunc(j * ((j + 1) * j), 2) - @divTrunc((j + 1) * j * (j - 1), 3)), 0);
        }
        for ([_]struct { u8, isize }{ .{ 'N', 1 }, .{ 'T', -1 } }) |c| {
            const res = mvch(T, c[0], n, n, ob.val(T, 1, 0), &a, n, &x, c[1], ob.zero(T), &y, c[1], &yt, &g, &yy);
            for (yt, yy) |t, e| try ob.expect(ob.eql(t, ob.cast(R, e)), "SMVCH: wrong");
            try ob.expect(res.err == 0 and !res.fatal, "SMVCH: ratio not 0");
        }
    }
}

test "blat2: xGEMV" {
    try sweep(.gemv);
}

test "blat2: xGBMV" {
    try sweep(.gbmv);
}

test "blat2: xSYMV, xHEMV" {
    try sweep(.symv);
}

test "blat2: xSBMV, xHBMV" {
    try sweep(.sbmv);
}

test "blat2: xSPMV, xHPMV" {
    try sweep(.spmv);
}

test "blat2: xTRMV" {
    try sweep(.trmv);
}

test "blat2: xTBMV" {
    try sweep(.tbmv);
}

test "blat2: xTPMV" {
    try sweep(.tpmv);
}

test "blat2: xTRSV" {
    try sweep(.trsv);
}

test "blat2: xTBSV" {
    try sweep(.tbsv);
}

test "blat2: xTPSV" {
    try sweep(.tpsv);
}

test "blat2: xGER, xGERU" {
    try sweep(.ger);
}

test "blat2: xGERC" {
    try sweep(.gerc);
}

test "blat2: xSYR, xHER" {
    try sweep(.syr);
}

test "blat2: xSPR, xHPR" {
    try sweep(.spr);
}

test "blat2: xSYR2, xHER2" {
    try sweep(.syr2);
}

test "blat2: xSPR2, xHPR2" {
    try sweep(.spr2);
}

// ---- SCHKE ------------------------------------------------------------------

/// `SCHKE`/`CCHKE`: each call has one bad argument, and must report it.
fn errorExits(comptime T: type) !void {
    const al = testing.allocator;
    const cx = comptime isComplex(T);
    var a: [2]T = @splat(ob.zero(T));
    var x: [2]T = @splat(ob.zero(T));
    var y: [2]T = @splat(ob.zero(T));
    const al0 = ob.zero(T);
    const ra: Float(T) = 0;
    const E = ob.expectError;
    // gemv: a bad letter is one real types do not take
    if (!cx) try E(error.InvalidTrans, blas.gemv(al, ob.io, T, .o, 0, 0, al0, &a, 1, &x, 1, al0, &y, 1), "GEMV INFOT=1");
    try E(error.InvalidM, blas.gemv(al, ob.io, T, .n, -1, 0, al0, &a, 1, &x, 1, al0, &y, 1), "GEMV INFOT=2");
    try E(error.InvalidN, blas.gemv(al, ob.io, T, .n, 0, -1, al0, &a, 1, &x, 1, al0, &y, 1), "GEMV INFOT=3");
    try E(error.InvalidLda, blas.gemv(al, ob.io, T, .n, 2, 0, al0, &a, 1, &x, 1, al0, &y, 1), "GEMV INFOT=6");
    try E(error.InvalidIncx, blas.gemv(al, ob.io, T, .n, 0, 0, al0, &a, 1, &x, 0, al0, &y, 1), "GEMV INFOT=8");
    try E(error.InvalidIncy, blas.gemv(al, ob.io, T, .n, 0, 0, al0, &a, 1, &x, 1, al0, &y, 0), "GEMV INFOT=11");
    // gbmv
    if (!cx) try E(error.InvalidTrans, blas.gbmv(al, ob.io, T, .o, 0, 0, 0, 0, al0, &a, 1, &x, 1, al0, &y, 1), "GBMV INFOT=1");
    try E(error.InvalidM, blas.gbmv(al, ob.io, T, .n, -1, 0, 0, 0, al0, &a, 1, &x, 1, al0, &y, 1), "GBMV INFOT=2");
    try E(error.InvalidN, blas.gbmv(al, ob.io, T, .n, 0, -1, 0, 0, al0, &a, 1, &x, 1, al0, &y, 1), "GBMV INFOT=3");
    try E(error.InvalidK, blas.gbmv(al, ob.io, T, .n, 0, 0, -1, 0, al0, &a, 1, &x, 1, al0, &y, 1), "GBMV INFOT=4");
    try E(error.InvalidK, blas.gbmv(al, ob.io, T, .n, 2, 0, 0, -1, al0, &a, 1, &x, 1, al0, &y, 1), "GBMV INFOT=5");
    try E(error.InvalidLda, blas.gbmv(al, ob.io, T, .n, 0, 0, 1, 0, al0, &a, 1, &x, 1, al0, &y, 1), "GBMV INFOT=8");
    try E(error.InvalidIncx, blas.gbmv(al, ob.io, T, .n, 0, 0, 0, 0, al0, &a, 1, &x, 0, al0, &y, 1), "GBMV INFOT=10");
    try E(error.InvalidIncy, blas.gbmv(al, ob.io, T, .n, 0, 0, 0, 0, al0, &a, 1, &x, 1, al0, &y, 0), "GBMV INFOT=13");
    // symv, hemv
    if (cx) {
        try E(error.InvalidN, blas.hemv(al, ob.io, T, .u, -1, al0, &a, 1, &x, 1, al0, &y, 1), "HEMV INFOT=2");
        try E(error.InvalidLda, blas.hemv(al, ob.io, T, .u, 2, al0, &a, 1, &x, 1, al0, &y, 1), "HEMV INFOT=5");
        try E(error.InvalidIncx, blas.hemv(al, ob.io, T, .u, 0, al0, &a, 1, &x, 0, al0, &y, 1), "HEMV INFOT=7");
        try E(error.InvalidIncy, blas.hemv(al, ob.io, T, .u, 0, al0, &a, 1, &x, 1, al0, &y, 0), "HEMV INFOT=10");
    } else {
        try E(error.InvalidUplo, blas.symv(al, ob.io, T, .v, 0, al0, &a, 1, &x, 1, al0, &y, 1), "SYMV INFOT=1");
        try E(error.InvalidN, blas.symv(al, ob.io, T, .u, -1, al0, &a, 1, &x, 1, al0, &y, 1), "SYMV INFOT=2");
        try E(error.InvalidLda, blas.symv(al, ob.io, T, .u, 2, al0, &a, 1, &x, 1, al0, &y, 1), "SYMV INFOT=5");
        try E(error.InvalidIncx, blas.symv(al, ob.io, T, .u, 0, al0, &a, 1, &x, 0, al0, &y, 1), "SYMV INFOT=7");
        try E(error.InvalidIncy, blas.symv(al, ob.io, T, .u, 0, al0, &a, 1, &x, 1, al0, &y, 0), "SYMV INFOT=10");
    }
    // sbmv, hbmv
    if (!cx) try E(error.InvalidUplo, blas.sbmv(al, ob.io, T, .v, 0, 0, al0, &a, 1, &x, 1, al0, &y, 1), "SBMV INFOT=1");
    if (cx) {
        try E(error.InvalidN, blas.hbmv(al, ob.io, T, .u, -1, 0, al0, &a, 1, &x, 1, al0, &y, 1), "HBMV INFOT=2");
        try E(error.InvalidK, blas.hbmv(al, ob.io, T, .u, 0, -1, al0, &a, 1, &x, 1, al0, &y, 1), "HBMV INFOT=3");
        try E(error.InvalidLda, blas.hbmv(al, ob.io, T, .u, 0, 1, al0, &a, 1, &x, 1, al0, &y, 1), "HBMV INFOT=6");
        try E(error.InvalidIncx, blas.hbmv(al, ob.io, T, .u, 0, 0, al0, &a, 1, &x, 0, al0, &y, 1), "HBMV INFOT=8");
        try E(error.InvalidIncy, blas.hbmv(al, ob.io, T, .u, 0, 0, al0, &a, 1, &x, 1, al0, &y, 0), "HBMV INFOT=11");
        // `HEMVREV` letters reach `hpmv` and the four below only from
        // row-major CBLAS upstream
        try E(error.InvalidUplo, blas.hpmv(al, ob.io, T, .v, 0, al0, &a, &x, 1, al0, &y, 1), "HPMV INFOT=1");
        try E(error.InvalidN, blas.hpmv(al, ob.io, T, .u, -1, al0, &a, &x, 1, al0, &y, 1), "HPMV INFOT=2");
        try E(error.InvalidIncx, blas.hpmv(al, ob.io, T, .u, 0, al0, &a, &x, 0, al0, &y, 1), "HPMV INFOT=6");
        try E(error.InvalidIncy, blas.hpmv(al, ob.io, T, .u, 0, al0, &a, &x, 1, al0, &y, 0), "HPMV INFOT=9");
    } else {
        try E(error.InvalidN, blas.sbmv(al, ob.io, T, .u, -1, 0, al0, &a, 1, &x, 1, al0, &y, 1), "SBMV INFOT=2");
        try E(error.InvalidK, blas.sbmv(al, ob.io, T, .u, 0, -1, al0, &a, 1, &x, 1, al0, &y, 1), "SBMV INFOT=3");
        try E(error.InvalidLda, blas.sbmv(al, ob.io, T, .u, 0, 1, al0, &a, 1, &x, 1, al0, &y, 1), "SBMV INFOT=6");
        try E(error.InvalidIncx, blas.sbmv(al, ob.io, T, .u, 0, 0, al0, &a, 1, &x, 0, al0, &y, 1), "SBMV INFOT=8");
        try E(error.InvalidIncy, blas.sbmv(al, ob.io, T, .u, 0, 0, al0, &a, 1, &x, 1, al0, &y, 0), "SBMV INFOT=11");
        try E(error.InvalidN, blas.spmv(al, ob.io, T, .u, -1, al0, &a, &x, 1, al0, &y, 1), "SPMV INFOT=2");
        try E(error.InvalidIncx, blas.spmv(al, ob.io, T, .u, 0, al0, &a, &x, 0, al0, &y, 1), "SPMV INFOT=6");
        try E(error.InvalidIncy, blas.spmv(al, ob.io, T, .u, 0, al0, &a, &x, 1, al0, &y, 0), "SPMV INFOT=9");
    }
    // the triangular six: a Hermitian `uplo`, an `x`-conjugating `trans`;
    // every `Diag` is valid
    try E(error.InvalidUplo, blas.trmv(al, ob.io, T, .v, .n, .n, 0, &a, 1, &x, 1), "TRMV INFOT=1");
    try E(error.InvalidTrans, blas.trmv(al, ob.io, T, .u, .o, .n, 0, &a, 1, &x, 1), "TRMV INFOT=2");
    inline for (.{ blas.trmv, blas.trsv }, .{ "TRMV", "TRSV" }) |f, nm| {
        try E(error.InvalidN, f(al, ob.io, T, .u, .n, .n, -1, &a, 1, &x, 1), nm ++ " INFOT=4");
        try E(error.InvalidLda, f(al, ob.io, T, .u, .n, .n, 2, &a, 1, &x, 1), nm ++ " INFOT=6");
        try E(error.InvalidIncx, f(al, ob.io, T, .u, .n, .n, 0, &a, 1, &x, 0), nm ++ " INFOT=8");
    }
    inline for (.{ blas.tbmv, blas.tbsv }, .{ "TBMV", "TBSV" }) |f, nm| {
        try E(error.InvalidN, f(al, ob.io, T, .u, .n, .n, -1, 0, &a, 1, &x, 1), nm ++ " INFOT=4");
        try E(error.InvalidK, f(al, ob.io, T, .u, .n, .n, 0, -1, &a, 1, &x, 1), nm ++ " INFOT=5");
        try E(error.InvalidLda, f(al, ob.io, T, .u, .n, .n, 0, 1, &a, 1, &x, 1), nm ++ " INFOT=7");
        try E(error.InvalidIncx, f(al, ob.io, T, .u, .n, .n, 0, 0, &a, 1, &x, 0), nm ++ " INFOT=9");
    }
    inline for (.{ blas.tpmv, blas.tpsv }, .{ "TPMV", "TPSV" }) |f, nm| {
        try E(error.InvalidN, f(al, ob.io, T, .u, .n, .n, -1, &a, &x, 1), nm ++ " INFOT=4");
        try E(error.InvalidIncx, f(al, ob.io, T, .u, .n, .n, 0, &a, &x, 0), nm ++ " INFOT=7");
    }
    // ger, geru, gerc
    inline for (if (cx) .{ blas.ger, blas.gerc } else .{blas.ger}) |f| {
        try E(error.InvalidM, f(al, ob.io, T, -1, 0, al0, &x, 1, &y, 1, &a, 1), "GER INFOT=1");
        try E(error.InvalidN, f(al, ob.io, T, 0, -1, al0, &x, 1, &y, 1, &a, 1), "GER INFOT=2");
        try E(error.InvalidIncx, f(al, ob.io, T, 0, 0, al0, &x, 0, &y, 1, &a, 1), "GER INFOT=5");
        try E(error.InvalidIncy, f(al, ob.io, T, 0, 0, al0, &x, 1, &y, 0, &a, 1), "GER INFOT=7");
        try E(error.InvalidLda, f(al, ob.io, T, 2, 0, al0, &x, 1, &y, 1, &a, 1), "GER INFOT=9");
    }
    // syr, spr, syr2, spr2; her, hpr, her2, hpr2
    if (cx) {
        try E(error.InvalidUplo, blas.her(al, ob.io, T, .m, 0, ra, &x, 1, &a, 1), "HER INFOT=1");
        try E(error.InvalidN, blas.her(al, ob.io, T, .u, -1, ra, &x, 1, &a, 1), "HER INFOT=2");
        try E(error.InvalidIncx, blas.her(al, ob.io, T, .u, 0, ra, &x, 0, &a, 1), "HER INFOT=5");
        try E(error.InvalidLda, blas.her(al, ob.io, T, .u, 2, ra, &x, 1, &a, 1), "HER INFOT=7");
        try E(error.InvalidUplo, blas.hpr(al, ob.io, T, .v, 0, ra, &x, 1, &a), "HPR INFOT=1");
        try E(error.InvalidN, blas.hpr(al, ob.io, T, .u, -1, ra, &x, 1, &a), "HPR INFOT=2");
        try E(error.InvalidIncx, blas.hpr(al, ob.io, T, .u, 0, ra, &x, 0, &a), "HPR INFOT=5");
        try E(error.InvalidUplo, blas.her2(al, ob.io, T, .m, 0, al0, &x, 1, &y, 1, &a, 1), "HER2 INFOT=1");
        try E(error.InvalidN, blas.her2(al, ob.io, T, .u, -1, al0, &x, 1, &y, 1, &a, 1), "HER2 INFOT=2");
        try E(error.InvalidIncx, blas.her2(al, ob.io, T, .u, 0, al0, &x, 0, &y, 1, &a, 1), "HER2 INFOT=5");
        try E(error.InvalidIncy, blas.her2(al, ob.io, T, .u, 0, al0, &x, 1, &y, 0, &a, 1), "HER2 INFOT=7");
        try E(error.InvalidLda, blas.her2(al, ob.io, T, .u, 2, al0, &x, 1, &y, 1, &a, 1), "HER2 INFOT=9");
        try E(error.InvalidUplo, blas.hpr2(al, ob.io, T, .v, 0, al0, &x, 1, &y, 1, &a), "HPR2 INFOT=1");
        try E(error.InvalidN, blas.hpr2(al, ob.io, T, .u, -1, al0, &x, 1, &y, 1, &a), "HPR2 INFOT=2");
        try E(error.InvalidIncx, blas.hpr2(al, ob.io, T, .u, 0, al0, &x, 0, &y, 1, &a), "HPR2 INFOT=5");
        try E(error.InvalidIncy, blas.hpr2(al, ob.io, T, .u, 0, al0, &x, 1, &y, 0, &a), "HPR2 INFOT=7");
    } else {
        try E(error.InvalidN, blas.syr(al, ob.io, T, .u, -1, al0, &x, 1, &a, 1), "SYR INFOT=2");
        try E(error.InvalidIncx, blas.syr(al, ob.io, T, .u, 0, al0, &x, 0, &a, 1), "SYR INFOT=5");
        try E(error.InvalidLda, blas.syr(al, ob.io, T, .u, 2, al0, &x, 1, &a, 1), "SYR INFOT=7");
        try E(error.InvalidN, blas.spr(al, ob.io, T, .u, -1, al0, &x, 1, &a), "SPR INFOT=2");
        try E(error.InvalidIncx, blas.spr(al, ob.io, T, .u, 0, al0, &x, 0, &a), "SPR INFOT=5");
        try E(error.InvalidN, blas.syr2(al, ob.io, T, .u, -1, al0, &x, 1, &y, 1, &a, 1), "SYR2 INFOT=2");
        try E(error.InvalidIncx, blas.syr2(al, ob.io, T, .u, 0, al0, &x, 0, &y, 1, &a, 1), "SYR2 INFOT=5");
        try E(error.InvalidIncy, blas.syr2(al, ob.io, T, .u, 0, al0, &x, 1, &y, 0, &a, 1), "SYR2 INFOT=7");
        try E(error.InvalidLda, blas.syr2(al, ob.io, T, .u, 2, al0, &x, 1, &y, 1, &a, 1), "SYR2 INFOT=9");
        try E(error.InvalidN, blas.spr2(al, ob.io, T, .u, -1, al0, &x, 1, &y, 1, &a), "SPR2 INFOT=2");
        try E(error.InvalidIncx, blas.spr2(al, ob.io, T, .u, 0, al0, &x, 0, &y, 1, &a), "SPR2 INFOT=5");
        try E(error.InvalidIncy, blas.spr2(al, ob.io, T, .u, 0, al0, &x, 1, &y, 0, &a), "SPR2 INFOT=7");
    }
}

test "blat2: error exits" {
    inline for (all_types) |T| errorExits(T) catch |e| {
        ob.say("{c}\n", .{letter(T)});
        return e;
    };
}
