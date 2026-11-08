//! The level-3 programs `test/[sdcz]blat3.f`, `[cz]blat3_3m.f` and their CBLAS
//! `ctest/` versions with their `.dat` tables, OpenBLAS v0.3.34-317; `SMMCH`
//! keeps the largest ratio over all columns and gauges 3M at three times the bound.

const std = @import("std");
const ob = @import("__root__.zig");
const blas = @import("../../__root__.zig");
const check = @import("../check.zig");
const guard = @import("../guard.zig");

const Float = ob.Float;
const isComplex = ob.isComplex;
const testing = std.testing;
const C32 = std.math.Complex(f32);
const C64 = std.math.Complex(f64);

const nmax = 65;

// ---- parameters ----------------------------------------------------------------

/// A `.dat` file's sweep: values of `N`, `ALPHA`, `BETA` and the threshold.
pub const Params = struct {
    thresh: f64 = 16,
    n: []const usize,
    alpha: []const [2]f64,
    beta: []const [2]f64,
};

/// `test/[sd]blat3.dat`: the first 6 of its `N` (`NIDIM` is 6).
pub const real_dat: Params = .{
    .n = &.{ 0, 1, 2, 3, 7, 31 },
    .alpha = &.{ .{ 0, 0 }, .{ 1, 0 }, .{ 0.7, 0 } },
    .beta = &.{ .{ 0, 0 }, .{ 1, 0 }, .{ 1.3, 0 } },
};

/// `test/[cz]blat3.dat`, `[cz]blat3_3m.dat`.
pub const complex_dat: Params = .{
    .n = &.{ 0, 1, 2, 3, 7, 31 },
    .alpha = &.{ .{ 0, 0 }, .{ 1, 0 }, .{ 0.7, -0.9 } },
    .beta = &.{ .{ 0, 0 }, .{ 1, 0 }, .{ 1.3, -1.1 } },
};

/// `ctest/[sdcz]in3` (and `[cz]in3_3m`, the same as `[cz]in3`): `N` by type.
fn ctestN(comptime T: type) []const usize {
    return switch (T) {
        f32 => &.{ 0, 1, 2, 3, 5, 9, 35 },
        f64 => &.{ 1, 2, 3, 5, 7, 9, 35 },
        C32 => &.{ 0, 1, 2, 3, 5, 9 },
        else => &.{ 0, 1, 2, 3, 5, 9, 35 },
    };
}

/// The rigs' `N`: the example data file in each program's header.
const rig_n = [_]usize{ 0, 1, 2, 3, 5, 9 };

fn dat(comptime T: type) Params {
    return if (comptime isComplex(T)) complex_dat else real_dat;
}

/// The `.dat` table natively, ctest's after it; the header's on the rigs.
fn tables(comptime T: type) [2]?Params {
    var p = dat(T);
    if (ob.reduced()) {
        p.n = &rig_n;
        return .{ p, null };
    }
    var c = p;
    c.n = ctestN(T);
    return .{ p, c };
}

const Routine = enum { gemm, gemm3m, symm, hemm, trmm, trsm, syrk, herk, syr2k, her2k };

fn letter(comptime T: type) u8 {
    return switch (Float(T)) {
        f32 => if (isComplex(T)) 'C' else 'S',
        else => if (isComplex(T)) 'Z' else 'D',
    };
}

// ---- SMAKE, LSERES, SMMCH ------------------------------------------------------

const Kind = enum { ge, sy, he, tr };

fn rogue(comptime T: type) T {
    return ob.val(T, -1e10, 1e10);
}

fn re(v: anytype) Float(@TypeOf(v)) {
    return if (comptime isComplex(@TypeOf(v))) v.re else v;
}

/// `SMAKE`/`CMAKE` (level 3): an `m` by `n` matrix into `a` (leading
/// dimension `ld`) and stored as `kind` in `aa`, unused elements rogue; off
/// the diagonal, column `n/2` is zero when `n > 3`.
fn make(comptime T: type, kind: Kind, uplo: u8, diag: u8, m: usize, n: usize, a: []T, ld: usize, aa: []T, lda: usize, gen: *ob.Gen) void {
    const gen_ = kind == .ge;
    const her = kind == .he;
    const sym = kind == .sy;
    const tri = kind == .tr;
    const up = (her or sym or tri) and uplo == 'U';
    const lo = (her or sym or tri) and uplo == 'L';
    const unit = tri and diag == 'U';
    const z = ob.zero(T);
    for (1..n + 1) |j| {
        for (1..m + 1) |i| {
            if (gen_ or (up and i <= j) or (lo and i >= j)) {
                a[(i - 1) + (j - 1) * ld] = gen.next(T);
                if (i != j) {
                    if (n > 3 and j == n / 2) a[(i - 1) + (j - 1) * ld] = z;
                    const v = a[(i - 1) + (j - 1) * ld];
                    if (her) {
                        a[(j - 1) + (i - 1) * ld] = ob.conj(v);
                    } else if (sym) {
                        a[(j - 1) + (i - 1) * ld] = v;
                    } else if (tri) a[(j - 1) + (i - 1) * ld] = z;
                }
            }
        }
        const d = &a[(j - 1) + (j - 1) * ld];
        if (her) d.* = ob.val(T, re(d.*), 0);
        if (tri) d.* = ob.add(d.*, ob.val(T, 1, 0));
        if (unit) d.* = ob.val(T, 1, 0);
    }
    for (1..n + 1) |j| {
        const ibeg: usize = if (gen_ or up) 1 else if (unit) j + 1 else j;
        const iend: usize = if (gen_) m else if (up) (if (unit) j - 1 else j) else n;
        for (1..lda + 1) |i| {
            aa[(i - 1) + (j - 1) * lda] = if (i >= ibeg and i <= iend) a[(i - 1) + (j - 1) * ld] else rogue(T);
        }
        if (her) {
            const p = &aa[(j - 1) + (j - 1) * lda];
            if (comptime isComplex(T)) p.im = -1e10;
        }
    }
}

/// `LSERES`: the storage outside the result is unchanged: rows past `m` for
/// `ge`, the other triangle for `sy` and `he`.
fn lseres(comptime T: type, kind: Kind, uplo: u8, m: usize, n: usize, aa: []const T, as: []const T, lda: usize) bool {
    for (0..n) |j| for (0..lda) |i| {
        const inside = switch (kind) {
            .ge => i < m,
            else => if (uplo == 'U') i <= j else i >= j and i < n,
        };
        if (!inside and !ob.eql(aa[i + j * lda], as[i + j * lda])) return false;
    };
    return true;
}

const Result = ob.Result;

/// `SMMCH`/`CMMCH`: `alpha*op(A)*op(B) + beta*C` for `m` by `n` and `kk`
/// from the dense `a`, `b`, `c`, against `cc`; gauges are sums of `ABS1`.
fn mmch(comptime T: type, transa: u8, transb: u8, m: usize, n: usize, kk: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: T, c: []const T, ldc: usize, ct: []ob.Ref(T), g: []Float(ob.Ref(T)), cc: []const T, ldcc: usize) Result {
    return mmchLoose(T, .none, transa, transb, m, n, kk, alpha, a, lda, b, ldb, beta, c, ldc, ct, g, cc, ldcc);
}

/// Which part of a result is held to three times the bound.
const Loose = enum { none, re, im };

/// `mmch` for `gemm3m`: the part it forms from three products, a difference
/// of them, gauged at three times the bound (Higham, SIAM J. Matrix Anal.
/// Appl. 13(3), 1992): the imaginary one, or the real one where exactly one
/// operand is conjugated.
fn mmch3m(comptime T: type, transa: u8, transb: u8, m: usize, n: usize, kk: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: T, c: []const T, ldc: usize, ct: []ob.Ref(T), g: []Float(ob.Ref(T)), cc: []const T, ldcc: usize) Result {
    const loose: Loose = if ((transa == 'C') == (transb == 'C')) .im else .re;
    return mmchLoose(T, loose, transa, transb, m, n, kk, alpha, a, lda, b, ldb, beta, c, ldc, ct, g, cc, ldcc);
}

/// An error `d` with its `loose` part a third.
fn loosen(d: anytype, loose: Loose) @TypeOf(d) {
    if (comptime !isComplex(@TypeOf(d))) return d;
    var v = d;
    switch (loose) {
        .none => {},
        .re => v.re /= 3,
        .im => v.im /= 3,
    }
    return v;
}

fn mmchLoose(comptime T: type, loose: Loose, transa: u8, transb: u8, m: usize, n: usize, kk: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: T, c: []const T, ldc: usize, ct: []ob.Ref(T), g: []Float(ob.Ref(T)), cc: []const T, ldcc: usize) Result {
    const R = ob.Ref(T);
    const RF = Float(R);
    const tra = transa == 'T' or transa == 'C';
    const trb = transb == 'T' or transb == 'C';
    const ra = ob.cast(R, alpha);
    const rb = ob.cast(R, beta);
    const eps: RF = std.math.floatEps(Float(T));
    var err: RF = 0;
    for (0..n) |j| {
        for (0..m) |i| {
            var acc = ob.zero(R);
            var gauge: RF = 0;
            for (0..kk) |k| {
                var x = ob.cast(R, if (tra) a[k + i * lda] else a[i + k * lda]);
                var y = ob.cast(R, if (trb) b[j + k * ldb] else b[k + j * ldb]);
                if (transa == 'C') x = ob.conj(x);
                if (transb == 'C') y = ob.conj(y);
                acc = ob.add(acc, ob.mul(x, y));
                gauge += ob.abs1(x) * ob.abs1(y);
            }
            const cv = ob.cast(R, c[i + j * ldc]);
            ct[i] = ob.add(ob.mul(ra, acc), ob.mul(rb, cv));
            g[i] = ob.abs1(ra) * gauge + ob.abs1(rb) * ob.abs1(cv);
        }
        for (0..m) |i| {
            var erri = ob.abs1(loosen(ob.sub(ct[i], ob.cast(R, cc[i + j * ldcc])), loose)) / eps;
            if (g[i] != 0) erri /= g[i];
            err = @max(err, erri);
            if (err * @sqrt(eps) >= 1) {
                ob.say("({d},{d}): expected {e}, computed {e}\n", .{ i + 1, j + 1, re(ct[i]), re(cc[i + j * ldcc]) });
                return .{ .err = @floatCast(err), .fatal = true };
            }
        }
    }
    return .{ .err = @floatCast(err), .fatal = false };
}

// ---- the calls -----------------------------------------------------------------

fn Work(comptime T: type) type {
    const R = ob.Ref(T);
    return struct {
        a: []T,
        aa: []T,
        as: []T,
        b: []T,
        bb: []T,
        bs: []T,
        c: []T,
        cc: []T,
        cs: []T,
        ab: []T,
        w: []T,
        ct: []R,
        g: []Float(R),

        fn init(al: std.mem.Allocator) !@This() {
            const nn = nmax * nmax;
            return .{
                .a = try al.alloc(T, nn),
                .aa = try al.alloc(T, nn),
                .as = try al.alloc(T, nn),
                .b = try al.alloc(T, nn),
                .bb = try al.alloc(T, nn),
                .bs = try al.alloc(T, nn),
                .c = try al.alloc(T, nn),
                .cc = try al.alloc(T, nn),
                .cs = try al.alloc(T, nn),
                .ab = try al.alloc(T, 2 * nn),
                .w = try al.alloc(T, 2 * nmax),
                .ct = try al.alloc(R, nmax),
                .g = try al.alloc(Float(R), nmax),
            };
        }
    };
}

const Call = struct {
    transa: u8 = 'N',
    transb: u8 = 'N',
    side: u8 = 'L',
    uplo: u8 = 'U',
    diag: u8 = 'N',
    m: usize = 0,
    n: usize = 0,
    k: usize = 0,
    lda: usize = 0,
    ldb: usize = 0,
    ldc: usize = 0,
    alpha: [2]f64 = .{ 0, 0 },
    beta: [2]f64 = .{ 0, 0 },

    fn describe(c: Call, comptime T: type, r: Routine, nc: usize) void {
        ob.say("{c}{s} call {d}: transa={c} transb={c} side={c} uplo={c} diag={c} m={d} n={d} k={d} lda={d} ldb={d} ldc={d} alpha=({d},{d}) beta=({d},{d})\n", .{
            letter(T), @tagName(r), nc, c.transa, c.transb, c.side, c.uplo, c.diag, c.m, c.n, c.k, c.lda, c.ldb, c.ldc, c.alpha[0], c.alpha[1], c.beta[0], c.beta[1],
        });
    }
};

const asTrans = ob.transOf;
const asSide = ob.sideOf;
const asUplo = ob.uploOf;
const asDiag = ob.diagOf;

fn run(comptime T: type, comptime r: Routine, c: Call, w: anytype) !void {
    var poison: guard.Poison = .{ .child = testing.allocator };
    const al = poison.allocator();
    const alpha = ob.val(T, c.alpha[0], c.alpha[1]);
    const beta = ob.val(T, c.beta[0], c.beta[1]);
    const m: isize = @intCast(c.m);
    const n: isize = @intCast(c.n);
    const k: isize = @intCast(c.k);
    const lda: isize = @intCast(c.lda);
    const ldb: isize = @intCast(c.ldb);
    const ldc: isize = @intCast(c.ldc);
    const a = w.aa.ptr;
    const b = w.bb.ptr;
    const cp = w.cc.ptr;
    switch (r) {
        .gemm, .gemm3m => try @field(blas, @tagName(r))(al, ob.io, T, asTrans(c.transa), asTrans(c.transb), m, n, k, alpha, a, lda, b, ldb, beta, cp, ldc),
        .symm, .hemm => try @field(blas, @tagName(r))(al, ob.io, T, asSide(c.side), asUplo(c.uplo), m, n, alpha, a, lda, b, ldb, beta, cp, ldc),
        .trmm, .trsm => try @field(blas, @tagName(r))(al, ob.io, T, asSide(c.side), asUplo(c.uplo), asTrans(c.transa), asDiag(c.diag), m, n, alpha, a, lda, b, ldb),
        .syrk => try blas.syrk(al, ob.io, T, asUplo(c.uplo), asTrans(c.transa), n, k, alpha, a, lda, beta, cp, ldc),
        .herk => try blas.herk(al, ob.io, T, asUplo(c.uplo), asTrans(c.transa), n, k, re(alpha), a, lda, re(beta), cp, ldc),
        .syr2k => try blas.syr2k(al, ob.io, T, asUplo(c.uplo), asTrans(c.transa), n, k, alpha, a, lda, b, ldb, beta, cp, ldc),
        .her2k => try blas.her2k(al, ob.io, T, asUplo(c.uplo), asTrans(c.transa), n, k, alpha, a, lda, b, ldb, re(beta), cp, ldc),
    }
}

/// One routine's sweep: calls made, and the largest test ratio.
const Sweep = struct {
    nc: usize = 0,
    errmax: f64 = 0,

    fn call(s: *Sweep, comptime T: type, comptime r: Routine, c: Call, w: anytype) !void {
        s.nc += 1;
        run(T, r, c, w) catch |e| {
            if (e == error.SkipZigTest) return e;
            ob.say("error exit taken on a valid call: ", .{});
            c.describe(T, r, s.nc);
            return e;
        };
    }

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
        ob.say("{c}{s} completed {d} calls with maximum test ratio {d}\n", .{ letter(T), @tagName(r), s.nc, s.errmax });
        return error.TestUnexpectedResult;
    }
};

fn bump(ld: usize) ?usize {
    return ob.bump(ld, nmax);
}

fn val(comptime T: type, v: [2]f64) T {
    return ob.val(T, v[0], v[1]);
}

// ---- SCHK1 to SCHK5 -------------------------------------------------------------

/// `SCHK1`: `gemm` (and `gemm3m`).
fn chk1(comptime T: type, comptime r: Routine, p: Params) !void {
    var arena: std.heap.ArenaAllocator = .init(testing.allocator);
    defer arena.deinit();
    const w = try Work(T).init(arena.allocator());
    var gen: ob.Gen = .{};
    var s: Sweep = .{};
    for (p.n) |m| for (p.n) |n| {
        const ldc = bump(m) orelse continue;
        const lcc = ldc * n;
        const null_ = n == 0 or m == 0;
        for (p.n) |k| for ("NTC") |ta| {
            const tra = ta != 'N';
            const ma = if (tra) k else m;
            const na = if (tra) m else k;
            const lda = bump(ma) orelse continue;
            const laa = lda * na;
            make(T, .ge, ' ', ' ', ma, na, w.a, nmax, w.aa, lda, &gen);
            for ("NTC") |tb| {
                const trb = tb != 'N';
                const mb = if (trb) n else k;
                const nb = if (trb) k else n;
                const ldb = bump(mb) orelse continue;
                const lbb = ldb * nb;
                make(T, .ge, ' ', ' ', mb, nb, w.b, nmax, w.bb, ldb, &gen);
                for (p.alpha) |alpha| for (p.beta) |beta| {
                    make(T, .ge, ' ', ' ', m, n, w.c, nmax, w.cc, ldc, &gen);
                    @memcpy(w.as[0..laa], w.aa[0..laa]);
                    @memcpy(w.bs[0..lbb], w.bb[0..lbb]);
                    @memcpy(w.cs[0..lcc], w.cc[0..lcc]);
                    const c: Call = .{ .transa = ta, .transb = tb, .m = m, .n = n, .k = k, .lda = lda, .ldb = ldb, .ldc = ldc, .alpha = alpha, .beta = beta };
                    try s.call(T, r, c, w);
                    try s.same(T, r, c, "A", ob.lse(T, w.as[0..laa], w.aa[0..laa]));
                    try s.same(T, r, c, "B", ob.lse(T, w.bs[0..lbb], w.bb[0..lbb]));
                    try s.same(T, r, c, "C", if (null_) ob.lse(T, w.cs[0..lcc], w.cc[0..lcc]) else lseres(T, .ge, ' ', m, n, w.cc, w.cs, ldc));
                    if (null_) continue;
                    const mm = if (r == .gemm3m) mmch3m else mmch;
                    try s.result(T, r, c, mm(T, ta, tb, m, n, k, val(T, alpha), w.a, nmax, w.b, nmax, val(T, beta), w.c, nmax, w.ct, w.g, w.cc, ldc));
                };
            }
        };
    };
    try s.finish(T, r, p);
}

/// `SCHK2`: `symm` (complex: `symm`, `hemm`).
fn chk2(comptime T: type, comptime r: Routine, p: Params) !void {
    var arena: std.heap.ArenaAllocator = .init(testing.allocator);
    defer arena.deinit();
    const w = try Work(T).init(arena.allocator());
    var gen: ob.Gen = .{};
    var s: Sweep = .{};
    for (p.n) |m| for (p.n) |n| {
        const ldc = bump(m) orelse continue;
        const lcc = ldc * n;
        const null_ = n == 0 or m == 0;
        const ldb = bump(m) orelse continue;
        const lbb = ldb * n;
        make(T, .ge, ' ', ' ', m, n, w.b, nmax, w.bb, ldb, &gen);
        for ("LR") |sd| {
            const left = sd == 'L';
            const na = if (left) m else n;
            const lda = bump(na) orelse continue;
            const laa = lda * na;
            for ("UL") |ul| {
                make(T, if (r == .hemm) .he else .sy, ul, ' ', na, na, w.a, nmax, w.aa, lda, &gen);
                for (p.alpha) |alpha| for (p.beta) |beta| {
                    make(T, .ge, ' ', ' ', m, n, w.c, nmax, w.cc, ldc, &gen);
                    @memcpy(w.as[0..laa], w.aa[0..laa]);
                    @memcpy(w.bs[0..lbb], w.bb[0..lbb]);
                    @memcpy(w.cs[0..lcc], w.cc[0..lcc]);
                    const c: Call = .{ .side = sd, .uplo = ul, .m = m, .n = n, .lda = lda, .ldb = ldb, .ldc = ldc, .alpha = alpha, .beta = beta };
                    try s.call(T, r, c, w);
                    try s.same(T, r, c, "A", ob.lse(T, w.as[0..laa], w.aa[0..laa]));
                    try s.same(T, r, c, "B", ob.lse(T, w.bs[0..lbb], w.bb[0..lbb]));
                    try s.same(T, r, c, "C", if (null_) ob.lse(T, w.cs[0..lcc], w.cc[0..lcc]) else lseres(T, .ge, ' ', m, n, w.cc, w.cs, ldc));
                    if (null_) continue;
                    try s.result(T, r, c, if (left)
                        mmch(T, 'N', 'N', m, n, m, val(T, alpha), w.a, nmax, w.b, nmax, val(T, beta), w.c, nmax, w.ct, w.g, w.cc, ldc)
                    else
                        mmch(T, 'N', 'N', m, n, n, val(T, alpha), w.b, nmax, w.a, nmax, val(T, beta), w.c, nmax, w.ct, w.g, w.cc, ldc));
                };
            }
        }
    };
    try s.finish(T, r, p);
}

/// `SCHK3`: `trmm`, `trsm`.
fn chk3(comptime T: type, comptime r: Routine, p: Params) !void {
    var arena: std.heap.ArenaAllocator = .init(testing.allocator);
    defer arena.deinit();
    const w = try Work(T).init(arena.allocator());
    @memset(w.c, ob.zero(T));
    var gen: ob.Gen = .{};
    var s: Sweep = .{};
    const one = ob.val(T, 1, 0);
    const z = ob.zero(T);
    for (p.n) |m| next_n: for (p.n) |n| {
        const ldb = bump(m) orelse continue;
        const lbb = ldb * n;
        const null_ = m == 0 or n == 0;
        for ("LR") |sd| {
            const left = sd == 'L';
            const na = if (left) m else n;
            const lda = bump(na) orelse continue :next_n;
            const laa = lda * na;
            for ("UL") |ul| for ("NTC") |ta| for ("UN") |dg| for (p.alpha) |alpha| {
                make(T, .tr, ul, dg, na, na, w.a, nmax, w.aa, lda, &gen);
                make(T, .ge, ' ', ' ', m, n, w.b, nmax, w.bb, ldb, &gen);
                @memcpy(w.as[0..laa], w.aa[0..laa]);
                @memcpy(w.bs[0..lbb], w.bb[0..lbb]);
                const c: Call = .{ .side = sd, .uplo = ul, .transa = ta, .diag = dg, .m = m, .n = n, .lda = lda, .ldb = ldb, .alpha = alpha };
                try s.call(T, r, c, w);
                try s.same(T, r, c, "A", ob.lse(T, w.as[0..laa], w.aa[0..laa]));
                try s.same(T, r, c, "B", if (null_) ob.lse(T, w.bs[0..lbb], w.bb[0..lbb]) else lseres(T, .ge, ' ', m, n, w.bb, w.bs, ldb));
                if (null_) continue;
                const av = val(T, alpha);
                if (r == .trmm) {
                    try s.result(T, r, c, if (left)
                        mmch(T, ta, 'N', m, n, m, av, w.a, nmax, w.b, nmax, z, w.c, nmax, w.ct, w.g, w.bb, ldb)
                    else
                        mmch(T, 'N', ta, m, n, n, av, w.b, nmax, w.a, nmax, z, w.c, nmax, w.ct, w.g, w.bb, ldb));
                } else {
                    // the computed solution into C, alpha*B where it was
                    for (0..n) |j| for (0..m) |i| {
                        w.c[i + j * nmax] = w.bb[i + j * ldb];
                        w.bb[i + j * ldb] = ob.mul(av, w.b[i + j * nmax]);
                    };
                    try s.result(T, r, c, if (left)
                        mmch(T, ta, 'N', m, n, m, one, w.a, nmax, w.c, nmax, z, w.b, nmax, w.ct, w.g, w.bb, ldb)
                    else
                        mmch(T, 'N', ta, m, n, n, one, w.c, nmax, w.a, nmax, z, w.b, nmax, w.ct, w.g, w.bb, ldb));
                }
            };
        }
    };
    try s.finish(T, r, p);
}

/// `SCHK4`: `syrk` (complex: `syrk` and `herk`, with real alpha and beta).
fn chk4(comptime T: type, comptime r: Routine, p: Params) !void {
    const herm = r == .herk;
    var arena: std.heap.ArenaAllocator = .init(testing.allocator);
    defer arena.deinit();
    const w = try Work(T).init(arena.allocator());
    var gen: ob.Gen = .{};
    var s: Sweep = .{};
    const transes = if (comptime isComplex(T)) "NC" else "NTC";
    for (p.n) |n| {
        const ldc = bump(n) orelse continue;
        const lcc = ldc * n;
        for (p.n) |k| for (transes) |t0| {
            const tran = t0 != 'N';
            const tc: u8 = if (tran and comptime isComplex(T) and !herm) 'T' else t0;
            const ma = if (tran) k else n;
            const na = if (tran) n else k;
            const lda = bump(ma) orelse continue;
            const laa = lda * na;
            make(T, .ge, ' ', ' ', ma, na, w.a, nmax, w.aa, lda, &gen);
            for ("UL") |ul| {
                const up = ul == 'U';
                for (p.alpha) |alpha0| for (p.beta) |beta0| {
                    const alpha: [2]f64 = if (herm) .{ alpha0[0], 0 } else alpha0;
                    const beta: [2]f64 = if (herm) .{ beta0[0], 0 } else beta0;
                    var null_ = n == 0;
                    if (herm) null_ = null_ or ((k == 0 or alpha[0] == 0) and beta[0] == 1);
                    make(T, if (herm) .he else .sy, ul, ' ', n, n, w.c, nmax, w.cc, ldc, &gen);
                    @memcpy(w.as[0..laa], w.aa[0..laa]);
                    @memcpy(w.cs[0..lcc], w.cc[0..lcc]);
                    const c: Call = .{ .uplo = ul, .transa = tc, .n = n, .k = k, .lda = lda, .ldc = ldc, .alpha = alpha, .beta = beta };
                    try s.call(T, r, c, w);
                    try s.same(T, r, c, "A", ob.lse(T, w.as[0..laa], w.aa[0..laa]));
                    try s.same(T, r, c, "C", if (null_) ob.lse(T, w.cs[0..lcc], w.cc[0..lcc]) else lseres(T, .sy, ul, n, n, w.cc, w.cs, ldc));
                    if (null_) continue;
                    const transt: u8 = if (herm) 'C' else 'T';
                    var jc: usize = 0;
                    for (0..n) |j| {
                        const jj = if (up) 0 else j;
                        const lj = if (up) j + 1 else n - j;
                        try s.result(T, r, c, if (tran)
                            mmch(T, transt, 'N', lj, 1, k, val(T, alpha), w.a[jj * nmax ..], nmax, w.a[j * nmax ..], nmax, val(T, beta), w.c[jj + j * nmax ..], nmax, w.ct, w.g, w.cc[jc..], ldc)
                        else
                            mmch(T, 'N', transt, lj, 1, k, val(T, alpha), w.a[jj..], nmax, w.a[j..], nmax, val(T, beta), w.c[jj + j * nmax ..], nmax, w.ct, w.g, w.cc[jc..], ldc));
                        jc += if (up) ldc else ldc + 1;
                    }
                };
            }
        };
    }
    try s.finish(T, r, p);
}

/// `SCHK5`: `syr2k` (complex: `syr2k` and `her2k`, with real beta). `ab`
/// holds `A` and `B` side by side (leading dimension `nmax`), or one over
/// the other (leading dimension `2*nmax`) when transposed.
fn chk5(comptime T: type, comptime r: Routine, p: Params) !void {
    const herm = r == .her2k;
    const cx = comptime isComplex(T);
    var arena: std.heap.ArenaAllocator = .init(testing.allocator);
    defer arena.deinit();
    const w = try Work(T).init(arena.allocator());
    var gen: ob.Gen = .{};
    var s: Sweep = .{};
    const one = ob.val(T, 1, 0);
    const transes = if (cx) "NC" else "NTC";
    for (p.n) |n| {
        const ldc = bump(n) orelse continue;
        const lcc = ldc * n;
        for (p.n) |k| for (transes) |t0| {
            const tran = t0 != 'N';
            const tc: u8 = if (tran and cx and !herm) 'T' else t0;
            const ma = if (tran) k else n;
            const na = if (tran) n else k;
            const lda = bump(ma) orelse continue;
            const laa = lda * na;
            const ld_ab: usize = if (tran) 2 * nmax else nmax;
            make(T, .ge, ' ', ' ', ma, na, w.ab, ld_ab, w.aa, lda, &gen);
            const ldb = lda;
            const lbb = laa;
            make(T, .ge, ' ', ' ', ma, na, w.ab[if (tran) k else k * nmax..], ld_ab, w.bb, ldb, &gen);
            for ("UL") |ul| {
                const up = ul == 'U';
                for (p.alpha) |alpha| for (p.beta) |beta0| {
                    const beta: [2]f64 = if (herm) .{ beta0[0], 0 } else beta0;
                    const av = val(T, alpha);
                    var null_ = n == 0;
                    if (herm) null_ = null_ or ((k == 0 or ob.isZero(av)) and beta[0] == 1);
                    make(T, if (herm) .he else .sy, ul, ' ', n, n, w.c, nmax, w.cc, ldc, &gen);
                    @memcpy(w.as[0..laa], w.aa[0..laa]);
                    @memcpy(w.bs[0..lbb], w.bb[0..lbb]);
                    @memcpy(w.cs[0..lcc], w.cc[0..lcc]);
                    const c: Call = .{ .uplo = ul, .transa = tc, .n = n, .k = k, .lda = lda, .ldb = ldb, .ldc = ldc, .alpha = alpha, .beta = beta };
                    try s.call(T, r, c, w);
                    try s.same(T, r, c, "A", ob.lse(T, w.as[0..laa], w.aa[0..laa]));
                    try s.same(T, r, c, "B", ob.lse(T, w.bs[0..lbb], w.bb[0..lbb]));
                    try s.same(T, r, c, "C", if (null_) ob.lse(T, w.cs[0..lcc], w.cc[0..lcc]) else lseres(T, .sy, ul, n, n, w.cc, w.cs, ldc));
                    if (null_) continue;
                    const transt: u8 = if (herm) 'C' else 'T';
                    // real `SCHK5` keeps alpha in the check; complex folds it into W
                    const ma_ = if (cx) one else av;
                    var jjab: usize = 0;
                    var jc: usize = 0;
                    for (0..n) |j| {
                        const jj = if (up) 0 else j;
                        const lj = if (up) j + 1 else n - j;
                        for (0..k) |i| {
                            const bv = if (tran) w.ab[j * 2 * nmax + k + i] else w.ab[(k + i) * nmax + j];
                            const avv = if (tran) w.ab[j * 2 * nmax + i] else w.ab[i * nmax + j];
                            if (!cx) {
                                w.w[i] = bv;
                                w.w[k + i] = avv;
                            } else if (tran) {
                                w.w[i] = ob.mul(av, bv);
                                w.w[k + i] = ob.mul(if (herm) ob.conj(av) else av, avv);
                            } else if (herm) {
                                w.w[i] = ob.mul(av, ob.conj(bv));
                                w.w[k + i] = ob.conj(ob.mul(av, avv));
                            } else {
                                w.w[i] = ob.mul(av, bv);
                                w.w[k + i] = ob.mul(av, avv);
                            }
                        }
                        try s.result(T, r, c, if (tran)
                            mmch(T, transt, 'N', lj, 1, 2 * k, ma_, w.ab[jjab..], 2 * nmax, w.w, 2 * nmax, val(T, beta), w.c[jj + j * nmax ..], nmax, w.ct, w.g, w.cc[jc..], ldc)
                        else
                            mmch(T, 'N', 'N', lj, 1, 2 * k, ma_, w.ab[jj..], nmax, w.w, 2 * nmax, val(T, beta), w.c[jj + j * nmax ..], nmax, w.ct, w.g, w.cc[jc..], ldc));
                        if (up) {
                            jc += ldc;
                        } else {
                            jc += ldc + 1;
                            if (tran) jjab += 2 * nmax;
                        }
                    }
                };
            }
        };
    }
    try s.finish(T, r, p);
}

fn chk(comptime T: type, comptime r: Routine, p: Params) !void {
    const f = switch (r) {
        .gemm, .gemm3m => chk1,
        .symm, .hemm => chk2,
        .trmm, .trsm => chk3,
        .syrk, .herk => chk4,
        .syr2k, .her2k => chk5,
    };
    try f(T, r, p);
}

/// `r` over `types`, each on its tables, after its error exits.
fn sweep(comptime r: Routine, comptime types: anytype) !void {
    inline for (types) |T| {
        try errorExits(T, r);
        for (tables(T)) |t| if (t) |p| try chk(T, r, p);
    }
}

/// Real `gemm` on the `N` of `n`, after its error exits, for `check.eachCore`.
pub fn all(n: []const usize) !void {
    inline for (.{ f32, f64 }) |T| {
        try errorExits(T, .gemm);
        var p = dat(T);
        p.n = n;
        try chk(T, .gemm, p);
    }
}

// ---- SCHKE ----------------------------------------------------------------------

const E = ob.expectError;

/// `SCHKE`/`CCHKE` for `r`: one bad argument per call, the first a bad letter
/// (not translatable) or one `r` does not take.
fn errorExits(comptime T: type, comptime r: Routine) !void {
    const al = testing.failing_allocator;
    var a: [4]T = @splat(ob.zero(T));
    var b: [4]T = @splat(ob.zero(T));
    var c: [4]T = @splat(ob.zero(T));
    const z = ob.zero(T);
    const rz: Float(T) = 0;
    const nm = @tagName(r);
    switch (r) {
        .gemm, .gemm3m => {
            const G = struct { ta: blas.Trans, tb: blas.Trans, m: isize, n: isize, k: isize, lda: isize, ldb: isize, ldc: isize, e: anyerror };
            const cases = [_]G{
                .{ .ta = .n, .tb = .n, .m = -1, .n = 0, .k = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidM },
                .{ .ta = .n, .tb = .t, .m = -1, .n = 0, .k = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidM },
                .{ .ta = .t, .tb = .n, .m = -1, .n = 0, .k = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidM },
                .{ .ta = .t, .tb = .t, .m = -1, .n = 0, .k = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidM },
                .{ .ta = .n, .tb = .n, .m = 0, .n = -1, .k = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidN },
                .{ .ta = .n, .tb = .t, .m = 0, .n = -1, .k = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidN },
                .{ .ta = .t, .tb = .n, .m = 0, .n = -1, .k = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidN },
                .{ .ta = .t, .tb = .t, .m = 0, .n = -1, .k = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidN },
                .{ .ta = .n, .tb = .n, .m = 0, .n = 0, .k = -1, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidK },
                .{ .ta = .n, .tb = .t, .m = 0, .n = 0, .k = -1, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidK },
                .{ .ta = .t, .tb = .n, .m = 0, .n = 0, .k = -1, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidK },
                .{ .ta = .t, .tb = .t, .m = 0, .n = 0, .k = -1, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidK },
                .{ .ta = .n, .tb = .n, .m = 2, .n = 0, .k = 0, .lda = 1, .ldb = 1, .ldc = 2, .e = error.InvalidLda },
                .{ .ta = .n, .tb = .t, .m = 2, .n = 0, .k = 0, .lda = 1, .ldb = 1, .ldc = 2, .e = error.InvalidLda },
                .{ .ta = .t, .tb = .n, .m = 0, .n = 0, .k = 2, .lda = 1, .ldb = 2, .ldc = 1, .e = error.InvalidLda },
                .{ .ta = .t, .tb = .t, .m = 0, .n = 0, .k = 2, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidLda },
                .{ .ta = .n, .tb = .n, .m = 0, .n = 0, .k = 2, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidLdb },
                .{ .ta = .t, .tb = .n, .m = 0, .n = 0, .k = 2, .lda = 2, .ldb = 1, .ldc = 1, .e = error.InvalidLdb },
                .{ .ta = .n, .tb = .t, .m = 0, .n = 2, .k = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidLdb },
                .{ .ta = .t, .tb = .t, .m = 0, .n = 2, .k = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidLdb },
                .{ .ta = .n, .tb = .n, .m = 2, .n = 0, .k = 0, .lda = 2, .ldb = 1, .ldc = 1, .e = error.InvalidLdc },
                .{ .ta = .n, .tb = .t, .m = 2, .n = 0, .k = 0, .lda = 2, .ldb = 1, .ldc = 1, .e = error.InvalidLdc },
                .{ .ta = .t, .tb = .n, .m = 2, .n = 0, .k = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidLdc },
                .{ .ta = .t, .tb = .t, .m = 2, .n = 0, .k = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidLdc },
            };
            for (cases) |g| try E(g.e, @field(blas, nm)(al, ob.io, T, g.ta, g.tb, g.m, g.n, g.k, z, &a, g.lda, &b, g.ldb, z, &c, g.ldc), nm);
        },
        .symm, .hemm => {
            const S = struct { s: u8, u: u8, m: isize, n: isize, lda: isize, ldb: isize, ldc: isize, e: anyerror };
            const cases = [_]S{
                .{ .s = 'L', .u = 'U', .m = -1, .n = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidM },
                .{ .s = 'R', .u = 'U', .m = -1, .n = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidM },
                .{ .s = 'L', .u = 'L', .m = -1, .n = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidM },
                .{ .s = 'R', .u = 'L', .m = -1, .n = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidM },
                .{ .s = 'L', .u = 'U', .m = 0, .n = -1, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidN },
                .{ .s = 'R', .u = 'U', .m = 0, .n = -1, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidN },
                .{ .s = 'L', .u = 'L', .m = 0, .n = -1, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidN },
                .{ .s = 'R', .u = 'L', .m = 0, .n = -1, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidN },
                .{ .s = 'L', .u = 'U', .m = 2, .n = 0, .lda = 1, .ldb = 2, .ldc = 2, .e = error.InvalidLda },
                .{ .s = 'R', .u = 'U', .m = 0, .n = 2, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidLda },
                .{ .s = 'L', .u = 'L', .m = 2, .n = 0, .lda = 1, .ldb = 2, .ldc = 2, .e = error.InvalidLda },
                .{ .s = 'R', .u = 'L', .m = 0, .n = 2, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidLda },
                .{ .s = 'L', .u = 'U', .m = 2, .n = 0, .lda = 2, .ldb = 1, .ldc = 2, .e = error.InvalidLdb },
                .{ .s = 'R', .u = 'U', .m = 2, .n = 0, .lda = 1, .ldb = 1, .ldc = 2, .e = error.InvalidLdb },
                .{ .s = 'L', .u = 'L', .m = 2, .n = 0, .lda = 2, .ldb = 1, .ldc = 2, .e = error.InvalidLdb },
                .{ .s = 'R', .u = 'L', .m = 2, .n = 0, .lda = 1, .ldb = 1, .ldc = 2, .e = error.InvalidLdb },
                .{ .s = 'L', .u = 'U', .m = 2, .n = 0, .lda = 2, .ldb = 2, .ldc = 1, .e = error.InvalidLdc },
                .{ .s = 'R', .u = 'U', .m = 2, .n = 0, .lda = 1, .ldb = 2, .ldc = 1, .e = error.InvalidLdc },
                .{ .s = 'L', .u = 'L', .m = 2, .n = 0, .lda = 2, .ldb = 2, .ldc = 1, .e = error.InvalidLdc },
                .{ .s = 'R', .u = 'L', .m = 2, .n = 0, .lda = 1, .ldb = 2, .ldc = 1, .e = error.InvalidLdc },
            };
            for (cases) |g| try E(g.e, @field(blas, nm)(al, ob.io, T, asSide(g.s), asUplo(g.u), g.m, g.n, z, &a, g.lda, &b, g.ldb, z, &c, g.ldc), nm);
        },
        .trmm, .trsm => {
            const S = struct { s: u8, u: u8, t: u8, m: isize, n: isize, lda: isize, ldb: isize, e: anyerror };
            var cases: [32]S = undefined;
            var i: usize = 0;
            // INFOT 5, 6, 9, 11 over every side, triangle and trans N, T
            for ("UL") |u| for ("LR") |s| for ("NT") |t| {
                const left = s == 'L';
                cases[i] = .{ .s = s, .u = u, .t = t, .m = -1, .n = 0, .lda = 1, .ldb = 1, .e = error.InvalidM };
                cases[i + 1] = .{ .s = s, .u = u, .t = t, .m = 0, .n = -1, .lda = 1, .ldb = 1, .e = error.InvalidN };
                cases[i + 2] = if (left)
                    .{ .s = s, .u = u, .t = t, .m = 2, .n = 0, .lda = 1, .ldb = 2, .e = error.InvalidLda }
                else
                    .{ .s = s, .u = u, .t = t, .m = 0, .n = 2, .lda = 1, .ldb = 1, .e = error.InvalidLda };
                cases[i + 3] = .{ .s = s, .u = u, .t = t, .m = 2, .n = 0, .lda = if (left) 2 else 1, .ldb = 1, .e = error.InvalidLdb };
                i += 4;
            };
            for (cases) |g| try E(g.e, @field(blas, nm)(al, ob.io, T, asSide(g.s), asUplo(g.u), asTrans(g.t), .n, g.m, g.n, z, &a, g.lda, &b, g.ldb), nm);
        },
        .syrk, .herk, .syr2k, .her2k => {
            const two = r == .syr2k or r == .her2k;
            // the second letter of each pair is the transpose this routine takes
            const tt: u8 = if (r == .herk or r == .her2k) 'C' else 'T';
            // a transpose it does not take: `herk`'s T, complex `syrk`'s C
            if (comptime isComplex(T)) try E(error.InvalidTrans, call2(T, r, al, 'U', if (tt == 'C') 'T' else 'C', 0, 0, 1, 1, 1, &a, &b, &c, z, rz), nm);
            const S = struct { u: u8, t: u8, n: isize, k: isize, lda: isize, ldb: isize, ldc: isize, e: anyerror };
            var cases: [20]S = undefined;
            var i: usize = 0;
            for ("UL") |u| for ([_]u8{ 'N', tt }) |t| {
                const nt = t == 'N';
                cases[i] = .{ .u = u, .t = t, .n = -1, .k = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidN };
                cases[i + 1] = .{ .u = u, .t = t, .n = 0, .k = -1, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidK };
                cases[i + 2] = if (nt)
                    .{ .u = u, .t = t, .n = 2, .k = 0, .lda = 1, .ldb = 1, .ldc = 2, .e = error.InvalidLda }
                else
                    .{ .u = u, .t = t, .n = 0, .k = 2, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidLda };
                cases[i + 3] = if (nt)
                    .{ .u = u, .t = t, .n = 2, .k = 0, .lda = 2, .ldb = 1, .ldc = 2, .e = error.InvalidLdb }
                else
                    .{ .u = u, .t = t, .n = 0, .k = 2, .lda = 2, .ldb = 1, .ldc = 1, .e = error.InvalidLdb };
                cases[i + 4] = if (nt)
                    .{ .u = u, .t = t, .n = 2, .k = 0, .lda = 2, .ldb = 2, .ldc = 1, .e = error.InvalidLdc }
                else
                    .{ .u = u, .t = t, .n = 2, .k = 0, .lda = 1, .ldb = 1, .ldc = 1, .e = error.InvalidLdc };
                i += 5;
            };
            for (cases) |g| {
                if (!two and g.e == error.InvalidLdb) continue;
                try E(g.e, call2(T, r, al, g.u, g.t, g.n, g.k, g.lda, g.ldb, g.ldc, &a, &b, &c, z, rz), nm);
            }
        },
    }
}

/// The rank-k updates with `alpha` and `beta` zero.
fn call2(comptime T: type, comptime r: Routine, al: std.mem.Allocator, u: u8, t: u8, n: isize, k: isize, lda: isize, ldb: isize, ldc: isize, a: [*]T, b: [*]T, c: [*]T, z: T, rz: Float(T)) blas.Level3Error!void {
    return switch (r) {
        .syrk => blas.syrk(al, ob.io, T, asUplo(u), asTrans(t), n, k, z, a, lda, z, c, ldc),
        .herk => blas.herk(al, ob.io, T, asUplo(u), asTrans(t), n, k, rz, a, lda, rz, c, ldc),
        .syr2k => blas.syr2k(al, ob.io, T, asUplo(u), asTrans(t), n, k, z, a, lda, b, ldb, z, c, ldc),
        .her2k => blas.her2k(al, ob.io, T, asUplo(u), asTrans(t), n, k, z, a, lda, b, ldb, rz, c, ldc),
        else => unreachable,
    };
}

// ---- tests ----------------------------------------------------------------------

test "blat3: SMMCH self-check on exact data" {
    inline for (.{ f32, f64, C32, C64 }) |T| {
        const R = ob.Ref(T);
        const n = 32;
        var ab: [n * (n + 1)]T = undefined;
        var c: [n]T = @splat(ob.zero(T));
        var cc: [n]T = undefined;
        var ct: [n]R = undefined;
        var g: [n]Float(R) = undefined;
        for (1..n + 1) |j| for (1..n + 1) |i| {
            ab[(i - 1) + (j - 1) * n] = ob.val(T, @floatFromInt(@as(i64, @intCast(i)) - @as(i64, @intCast(j)) + 1), 0);
            if (i < j) ab[(i - 1) + (j - 1) * n] = ob.zero(T);
        };
        const bcol = ab[n * n ..];
        for ([_][2]u8{ "NN".*, "NT".*, "TN".*, "TT".* }) |t| {
            const rev = t[0] == 'T';
            for (1..n + 1) |jj| {
                const j: i64 = @intCast(jj);
                bcol[jj - 1] = ob.val(T, @floatFromInt(if (rev) n - j + 1 else j), 0);
                const v = @divTrunc(j * ((j + 1) * j), 2) - @divTrunc((j + 1) * j * (j - 1), 3);
                cc[if (rev) n - jj else jj - 1] = ob.val(T, @floatFromInt(v), 0);
            }
            // B is `bcol` as a column (N) or a row with leading dimension 1 (T)
            const res = mmch(T, t[0], t[1], n, 1, n, ob.val(T, 1, 0), &ab, n, bcol, if (t[1] == 'T') 1 else n, ob.zero(T), &c, n, &ct, &g, &cc, n);
            for (ct, cc) |x, y| try ob.expect(ob.eql(x, ob.cast(R, y)), "SMMCH: wrong");
            try ob.expect(res.err == 0 and !res.fatal, "SMMCH: ratio not 0");
        }
    }
}

test "blat3: SGEMM, DGEMM" {
    try sweep(.gemm, .{ f32, f64 });
}

test "blat3: CGEMM, ZGEMM" {
    try sweep(.gemm, .{ C32, C64 });
}

test "blat3_3m: CGEMM3M, ZGEMM3M" {
    try sweep(.gemm3m, .{ C32, C64 });
}

test "blat3: xSYMM" {
    try sweep(.symm, .{ f32, f64, C32, C64 });
}

test "blat3: CHEMM, ZHEMM" {
    try sweep(.hemm, .{ C32, C64 });
}

test "blat3: xTRMM" {
    try sweep(.trmm, .{ f32, f64, C32, C64 });
}

test "blat3: xTRSM" {
    try sweep(.trsm, .{ f32, f64, C32, C64 });
}

test "blat3: xSYRK" {
    try sweep(.syrk, .{ f32, f64, C32, C64 });
}

test "blat3: CHERK, ZHERK" {
    try sweep(.herk, .{ C32, C64 });
}

test "blat3: xSYR2K" {
    try sweep(.syr2k, .{ f32, f64, C32, C64 });
}

test "blat3: CHER2K, ZHER2K" {
    try sweep(.her2k, .{ C32, C64 });
}
