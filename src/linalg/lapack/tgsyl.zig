//! `tgsyl`: ReLAPACK's recursion (`relapack/src/dtgsyl.c` and the rest,
//! OpenBLAS v0.3.34-317, MIT) over `sylvester.zig`'s `tgsy2`.

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const driver = @import("../blas/blas3/driver.zig");
const ops = @import("ops.zig");
const relapack = @import("relapack.zig");
const sylvester = @import("sylvester.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const isComplex = types.isComplex;
const Trans = types.Trans;
const scalar = ops.scalar;
const Solved = sylvester.Solved;
const scaleBlock = sylvester.scaleBlock;
const Dif = sylvester.Dif;
const tgsy2 = sylvester.tgsy2;

const Error = driver.Error;

/// What `tgsyl` returns: bad arguments in the reference's order (`m` or `n`
/// below 1 is invalid, as there); `Perturbed` (the reference's `INFO > 0`:
/// `(A, D)` and `(B, E)` have common or close eigenvalues, so perturbed
/// values were used), with `C`, `F`, `scale` and `dif` written;
/// `error.OutOfMemory` or `error.Canceled`.
pub const TgsylError = error{ InvalidTrans, InvalidIjob, InvalidM, InvalidN, InvalidLda, InvalidLdb, InvalidLdc, InvalidLdd, InvalidLde, InvalidLdf, Perturbed } || Allocator.Error || Io.Cancelable;

/// `?TGSYL` by ReLAPACK's recursion (`relapack/src/dtgsyl.c` and the rest):
/// `A R - L B = scale C`, `D R - L E = scale F` for `R` and `L` (`m x n`),
/// into `C` and `F`, `(A, D)` and `(B, E)` in generalized Schur form;
/// `trans` `t` (real) or `c` (complex) solves `Aᴴ R + Dᴴ L = scale C`,
/// `R Bᴴ + L Eᴴ = -scale F`. `scale <= 1` keeps `R` and `L` finite. `ijob`
/// (read only for `trans` `n`): 0 solves; 3 and 4 only estimate
/// `Dif[(A, D), (B, E)]` into `dif` (by `xLATDF`'s look-ahead or its `xGECON`
/// null vectors), leaving that estimate's `R` and `L` in `C` and `F`; 1 and 2
/// do both. `dif` is not written for `ijob` 0 or `trans`.
pub fn tgsyl(allocator: Allocator, io: Io, comptime T: type, trans: Trans, ijob: isize, m: isize, n: isize, a: [*]const T, lda: isize, b: [*]const T, ldb: isize, c: [*]T, ldc: isize, d: [*]const T, ldd: isize, e: [*]const T, lde: isize, f: [*]T, ldf: isize, scale: *Float(T), dif: *Float(T)) TgsylError!void {
    comptime ops.requireFloat(T, "tgsyl");
    const F = Float(T);
    // `N`, or `T` for real and `C` for complex, as `dtgsyl.f` and `ztgsyl.f`
    const tr = switch (trans) {
        .n => false,
        .t => if (comptime isComplex(T)) return error.InvalidTrans else true,
        .c => if (comptime isComplex(T)) true else return error.InvalidTrans,
        .r, .o, .u, .s, .d => return error.InvalidTrans,
    };
    if (!tr and (ijob < 0 or ijob > 4)) return error.InvalidIjob;
    if (m <= 0) return error.InvalidM;
    if (n <= 0) return error.InvalidN;
    if (lda < m) return error.InvalidLda;
    if (ldb < n) return error.InvalidLdb;
    if (ldc < m) return error.InvalidLdc;
    if (ldd < m) return error.InvalidLdd;
    if (lde < n) return error.InvalidLde;
    if (ldf < m) return error.InvalidLdf;
    inline for (.{ .{ m, lda, error.InvalidLda }, .{ n, ldb, error.InvalidLdb }, .{ m, ldd, error.InvalidLdd }, .{ n, lde, error.InvalidLde } }) |x| {
        if (!util.fitsMatrix(x[0], x[0], x[1], @sizeOf(T))) return x[2];
    }
    if (!util.fitsMatrix(m, n, ldc, @sizeOf(T))) return error.InvalidLdc;
    if (!util.fitsMatrix(m, n, ldf, @sizeOf(T))) return error.InvalidLdf;
    const um: usize = @intCast(m);
    const un: usize = @intCast(n);
    const p: Gsyl(T) = .{ .tr = tr, .a = a, .lda = @intCast(lda), .b = b, .ldb = @intCast(ldb), .c = c, .ldc = @intCast(ldc), .d = d, .ldd = @intCast(ldd), .e = e, .lde = @intCast(lde), .f = f, .ldf = @intCast(ldf) };
    const job: usize = if (tr) 0 else @intCast(ijob);
    // `ijob` 1 and 2 solve, then estimate on zeros; 3 and 4 only estimate
    var ifunc: u2 = if (job >= 3) @intCast(job - 2) else 0;
    const isolve: usize = if (job == 1 or job == 2) 2 else 1;
    if (job >= 3) zeroBoth(T, p, um, un);
    const saved = try allocator.alloc(T, if (isolve == 2) 2 * um * un else 0);
    defer allocator.free(saved);
    var perturbed = false;
    var scale2: F = 1;
    for (0..isolve) |round| {
        var acc: Dif(F) = .{};
        const r = try tgsylRec(allocator, io, T, ifunc, um, un, p, &acc);
        scale.* = r.scale;
        if (r.perturbed) perturbed = true;
        if (acc.scale != 0) {
            const num: F = @floatFromInt(if (job == 1 or job == 3) 2 * um * un else acc.pq);
            dif.* = @sqrt(num) / (acc.scale * @sqrt(acc.sum));
        }
        if (isolve == 2 and round == 0) {
            ifunc = @intCast(job);
            scale2 = scale.*;
            copyBlock(T, um, un, c, p.ldc, saved.ptr, um);
            copyBlock(T, um, un, f, p.ldf, saved.ptr + um * un, um);
            zeroBoth(T, p, um, un);
        } else if (isolve == 2) {
            copyBlock(T, um, un, saved.ptr, um, c, p.ldc);
            copyBlock(T, um, un, saved.ptr + um * un, um, f, p.ldf);
            scale.* = scale2;
        }
    }
    if (perturbed) return error.Perturbed;
}

/// `tgsyl`'s six matrices, as one recursion level sees them.
fn Gsyl(comptime T: type) type {
    return struct {
        tr: bool,
        a: [*]const T,
        lda: usize,
        b: [*]const T,
        ldb: usize,
        c: [*]T,
        ldc: usize,
        d: [*]const T,
        ldd: usize,
        e: [*]const T,
        lde: usize,
        f: [*]T,
        ldf: usize,

        /// Rows `i..` of `C`, `F` and the `A`, `D` blocks from `(i, i)`.
        fn rows(s: @This(), i: usize) @This() {
            var x = s;
            x.a += i * s.lda + i;
            x.d += i * s.ldd + i;
            x.c += i;
            x.f += i;
            return x;
        }

        /// Columns `j..` of `C`, `F` and the `B`, `E` blocks from `(j, j)`.
        fn cols(s: @This(), j: usize) @This() {
            var x = s;
            x.b += j * s.ldb + j;
            x.e += j * s.lde + j;
            x.c += j * s.ldc;
            x.f += j * s.ldf;
            return x;
        }
    };
}

fn zeroBoth(comptime T: type, p: Gsyl(T), m: usize, n: usize) void {
    for (0..n) |j| for (0..m) |i| {
        p.c[i + j * p.ldc] = scalar(T, 0);
        p.f[i + j * p.ldf] = scalar(T, 0);
    };
}

/// `?LACPY('F')`.
fn copyBlock(comptime T: type, m: usize, n: usize, src: [*]const T, lds: usize, dst: [*]T, ldd: usize) void {
    for (0..n) |j| @memcpy(dst[j * ldd ..][0..m], src[j * lds ..][0..m]);
}

/// `RELAPACK_?tgsyl_rec`: the halves of the larger side, `dif`'s sums
/// carried through.
fn tgsylRec(al: Allocator, io: Io, comptime T: type, ifunc: u2, m: usize, n: usize, p: Gsyl(T), acc: *Dif(Float(T))) Error!Solved(Float(T)) {
    if (m <= relapack.crossover and n <= relapack.crossover) return tgsy2(al, io, T, p.tr, ifunc, m, n, p.a, p.lda, p.b, p.ldb, p.c, p.ldc, p.d, p.ldd, p.e, p.lde, p.f, p.ldf, acc);
    const one = scalar(T, 1);
    const minus = scalar(T, -1);
    const ct = relapack.ctrans(T);
    var r1: Solved(Float(T)) = undefined;
    var r2: Solved(Float(T)) = undefined;
    if (m > n) {
        var m1 = relapack.split(T, m);
        // never inside a 2 x 2 block
        if (comptime !isComplex(T)) if (p.a[m1 + p.lda * (m1 - 1)] != 0) {
            m1 += 1;
        };
        const m2 = m - m1;
        const bot = p.rows(m1);
        const a_tr = p.a + p.lda * m1;
        const d_tr = p.d + p.ldd * m1;
        if (!p.tr) {
            // C_T = scale1 C_T - A_TR C_B, F_T = scale1 F_T - D_TR C_B
            r1 = try tgsylRec(al, io, T, ifunc, m2, n, bot, acc);
            const s1 = scalar(T, r1.scale);
            try relapack.gemm(al, io, T, .n, .n, m1, n, m2, minus, a_tr, p.lda, bot.c, p.ldc, s1, p.c, p.ldc);
            try relapack.gemm(al, io, T, .n, .n, m1, n, m2, minus, d_tr, p.ldd, bot.c, p.ldc, s1, p.f, p.ldf);
            r2 = try tgsylRec(al, io, T, ifunc, m1, n, p, acc);
            if (r2.scale != 1) {
                scaleBlock(T, m2, n, r2.scale, bot.c, p.ldc);
                scaleBlock(T, m2, n, r2.scale, bot.f, p.ldf);
            }
        } else {
            // C_B = scale1 C_B - A_TRᴴ C_T - D_TRᴴ F_T, F_B = scale1 F_B
            r1 = try tgsylRec(al, io, T, ifunc, m1, n, p, acc);
            if (r1.scale != 1) scaleBlock(T, m2, n, r1.scale, bot.f, p.ldf);
            try relapack.gemm(al, io, T, ct, .n, m2, n, m1, minus, a_tr, p.lda, p.c, p.ldc, scalar(T, r1.scale), bot.c, p.ldc);
            try relapack.gemm(al, io, T, ct, .n, m2, n, m1, minus, d_tr, p.ldd, p.f, p.ldf, one, bot.c, p.ldc);
            r2 = try tgsylRec(al, io, T, ifunc, m2, n, bot, acc);
            if (r2.scale != 1) {
                scaleBlock(T, m1, n, r2.scale, p.c, p.ldc);
                scaleBlock(T, m1, n, r2.scale, p.f, p.ldf);
            }
        }
    } else {
        var n1 = relapack.split(T, n);
        if (comptime !isComplex(T)) if (p.b[n1 + p.ldb * (n1 - 1)] != 0) {
            n1 += 1;
        };
        const n2 = n - n1;
        const right = p.cols(n1);
        const b_tr = p.b + p.ldb * n1;
        const e_tr = p.e + p.lde * n1;
        if (!p.tr) {
            // C_R = scale1 C_R + F_L B_TR, F_R = scale1 F_R + F_L E_TR
            r1 = try tgsylRec(al, io, T, ifunc, m, n1, p, acc);
            const s1 = scalar(T, r1.scale);
            try relapack.gemm(al, io, T, .n, .n, m, n2, n1, one, p.f, p.ldf, b_tr, p.ldb, s1, right.c, p.ldc);
            try relapack.gemm(al, io, T, .n, .n, m, n2, n1, one, p.f, p.ldf, e_tr, p.lde, s1, right.f, p.ldf);
            r2 = try tgsylRec(al, io, T, ifunc, m, n2, right, acc);
            if (r2.scale != 1) {
                scaleBlock(T, m, n1, r2.scale, p.c, p.ldc);
                scaleBlock(T, m, n1, r2.scale, p.f, p.ldf);
            }
        } else {
            // F_L = scale1 F_L + C_R B_TRᴴ + F_R E_TRᴴ, C_L = scale1 C_L
            r1 = try tgsylRec(al, io, T, ifunc, m, n2, right, acc);
            if (r1.scale != 1) scaleBlock(T, m, n1, r1.scale, p.c, p.ldc);
            try relapack.gemm(al, io, T, .n, ct, m, n1, n2, one, right.c, p.ldc, b_tr, p.ldb, scalar(T, r1.scale), p.f, p.ldf);
            try relapack.gemm(al, io, T, .n, ct, m, n1, n2, one, right.f, p.ldf, e_tr, p.lde, one, p.f, p.ldf);
            r2 = try tgsylRec(al, io, T, ifunc, m, n1, p, acc);
            if (r2.scale != 1) {
                scaleBlock(T, m, n2, r2.scale, right.c, p.ldc);
                scaleBlock(T, m, n2, r2.scale, right.f, p.ldf);
            }
        }
    }
    return .{ .scale = r1.scale * r2.scale, .perturbed = r1.perturbed or r2.perturbed };
}

// ---- tests ---------------------------------------------------------------------

const testing = std.testing;
const check = @import("../blas/_test/check.zig");
const reference = @import("../blas/_test/reference.zig");
const lcheck = @import("_test/check.zig");
const fixtures = @import("_test/fixtures.zig");
const faer_svd = @import("svd.zig");

const schur = fixtures.schur;

/// An `n x n` upper triangular matrix like `schur`'s, without its 2 x 2
/// blocks: `D` and `E` of a generalized Schur form.
fn upperTri(comptime T: type, al: Allocator, n: usize, lo: f64, seed: u64) ![]T {
    const a = try schur(T, al, n, lo, n + 1, seed);
    if (comptime !isComplex(T)) for (0..n -| 1) |i| {
        a[i + 1 + i * n] = 0;
    };
    return a;
}

/// `tgsyl`'s operands, each at its own padded leading dimension so a mixed
/// up one shows: `(A, D)` with eigenvalues in `[0.8, 8]`, `(B, E)` in
/// `[-8, -0.8]`.
fn Gen(comptime T: type) type {
    return struct {
        m: usize,
        n: usize,
        a: []T,
        b: []T,
        d: []T,
        e: []T,
        c0: []T,
        f0: []T,
        c: []T,
        f: []T,
        const ld = .{ .a = 1, .b = 2, .c = 3, .d = 2, .e = 1, .f = 5 };

        const Self = @This();

        fn init(al: Allocator, m: usize, n: usize, seed: u64) !Self {
            var g: Self = undefined;
            g.m = m;
            g.n = n;
            g.a = try padded(al, m, ld.a, try schur(T, al, m, 2, 3, seed));
            g.d = try padded(al, m, ld.d, try upperTri(T, al, m, 0.5, seed + 1));
            g.b = try padded(al, n, ld.b, try schur(T, al, n, -4, 4, seed + 2));
            g.e = try padded(al, n, ld.e, try upperTri(T, al, n, 0.5, seed + 3));
            g.c0 = try al.alloc(T, (m + ld.c) * n);
            check.fill(T, g.c0, seed + 4);
            g.f0 = try al.alloc(T, (m + ld.f) * n);
            check.fill(T, g.f0, seed + 5);
            g.c = try al.dupe(T, g.c0);
            g.f = try al.dupe(T, g.f0);
            return g;
        }

        /// `x` (`k x k`, packed) at leading dimension `k + pad`, freed.
        fn padded(al: Allocator, k: usize, pad: usize, x: []T) ![]T {
            defer al.free(x);
            const y = try al.alloc(T, (k + pad) * k);
            @memset(y, lcheck.poisoned(T));
            for (0..k) |j| @memcpy(y[j * (k + pad) ..][0..k], x[j * k ..][0..k]);
            return y;
        }

        fn deinit(g: Self, al: Allocator) void {
            for ([_][]T{ g.a, g.b, g.d, g.e, g.c0, g.f0, g.c, g.f }) |x| al.free(x);
        }

        fn reset(g: Self) void {
            @memcpy(g.c, g.c0);
            @memcpy(g.f, g.f0);
        }

        fn run(g: Self, al: Allocator, io: Io, trans: Trans, ijob: isize, scale: *Float(T), dif: *Float(T)) TgsylError!void {
            const m = g.m;
            const n = g.n;
            return tgsyl(al, io, T, trans, ijob, @intCast(m), @intCast(n), g.a.ptr, @intCast(m + ld.a), g.b.ptr, @intCast(n + ld.b), g.c.ptr, @intCast(m + ld.c), g.d.ptr, @intCast(m + ld.d), g.e.ptr, @intCast(n + ld.e), g.f.ptr, @intCast(m + ld.f), scale, dif);
        }

        fn at(x: []const T, ldx: usize, i: usize, j: usize) Wt {
            return reference.wide(T, x[i + j * ldx]);
        }
        const Wt = @TypeOf(reference.wide(T, scalar(T, 0)));

        /// Both equations' residuals within `32 (m + n + 2)` ulps of their
        /// terms' magnitudes; padding untouched.
        fn expectSolved(g: Self, tr: bool, scale: Float(T)) !void {
            const m = g.m;
            const n = g.n;
            const lda = m + ld.a;
            const ldb = n + ld.b;
            const ldc = m + ld.c;
            const ldd = m + ld.d;
            const lde = n + ld.e;
            const ldf = m + ld.f;
            const eps: f128 = std.math.floatEps(Float(T));
            const ulps: f128 = @floatFromInt(32 * (m + n + 2));
            const cjw = struct {
                fn f(x: Wt, c: bool) Wt {
                    return if (c) reference.conjW(x) else x;
                }
            }.f;
            for (0..n) |j| for (0..m) |i| for (0..2) |eq| {
                const rhs = if (eq == 0) at(g.c0, ldc, i, j) else at(g.f0, ldf, i, j);
                // notrans: A R - L B, D R - L E; trans: Aᴴ R + Dᴴ L, R Bᴴ + L Eᴴ
                var r = reference.scaleW(if (tr and eq == 1) @as(f128, scale) else -@as(f128, scale), rhs);
                var mag: f128 = scale * reference.abs1(rhs);
                for (0..@max(m, n)) |k| {
                    var terms: [2][2]Wt = undefined;
                    var used: [2]bool = .{ false, false };
                    if (!tr) {
                        if (k < m) {
                            terms[0] = .{ at(if (eq == 0) g.a else g.d, if (eq == 0) lda else ldd, i, k), at(g.c, ldc, k, j) };
                            used[0] = true;
                        }
                        if (k < n) {
                            terms[1] = .{ reference.scaleW(-1, at(g.f, ldf, i, k)), at(if (eq == 0) g.b else g.e, if (eq == 0) ldb else lde, k, j) };
                            used[1] = true;
                        }
                    } else if (eq == 0) {
                        if (k < m) {
                            terms[0] = .{ cjw(at(g.a, lda, k, i), true), at(g.c, ldc, k, j) };
                            terms[1] = .{ cjw(at(g.d, ldd, k, i), true), at(g.f, ldf, k, j) };
                            used = .{ true, true };
                        }
                    } else if (k < n) {
                        terms[0] = .{ at(g.c, ldc, i, k), cjw(at(g.b, ldb, j, k), true) };
                        terms[1] = .{ at(g.f, ldf, i, k), cjw(at(g.e, lde, j, k), true) };
                        used = .{ true, true };
                    }
                    for (terms, used) |t, u| if (u) {
                        r = reference.addW(r, reference.mulW(t[0], t[1]));
                        mag += reference.abs1(t[0]) * reference.abs1(t[1]);
                    };
                }
                const dv = reference.abs1(r);
                if (!(dv <= ulps * eps * mag + std.math.floatMin(Float(T)))) {
                    check.say("tgsyl {s} {d}x{d} trans {any}: equation {d} ({d},{d}) residual {e}, allowed {e}\n", .{ @typeName(T), m, n, tr, eq + 1, i, j, @as(f64, @floatCast(dv)), @as(f64, @floatCast(ulps * eps * mag)) });
                    return error.TestUnexpectedResult;
                }
            };
        }
    };
}

fn transOf(comptime T: type, tr: bool) Trans {
    return if (!tr) .n else if (comptime isComplex(T)) .c else .t;
}

const tgsyl_shapes = [_][2]usize{ .{ 1, 1 }, .{ 2, 3 }, .{ 5, 4 }, .{ 24, 24 }, .{ 25, 10 }, .{ 10, 40 }, .{ 60, 50 }, .{ 33, 70 } };

test "relapack tgsyl: A R - L B = scale C, D R - L E = scale F and the transposed system, every type" {
    const al = testing.allocator;
    inline for (lcheck.test_types) |T| for (tgsyl_shapes) |s| {
        if (check.rig() and s[0] * s[1] > 1000) continue;
        const g = try Gen(T).init(al, s[0], s[1], s[0] * 7 + s[1]);
        defer g.deinit(al);
        for ([_]bool{ false, true }) |tr| {
            g.reset();
            var scale: Float(T) = 0;
            var dif: Float(T) = -1;
            try g.run(al, check.io, transOf(T, tr), 0, &scale, &dif);
            try testing.expect(scale > 0 and scale <= 1 and dif == -1);
            g.expectSolved(tr, scale) catch |e| {
                check.say("tgsyl {s} {d}x{d}\n", .{ @typeName(T), s[0], s[1] });
                return e;
            };
        }
        // ijob 1 and 2 solve as 0 does, bit for bit, and estimate
        const x = try al.dupe(T, g.c);
        defer al.free(x);
        const y = try al.dupe(T, g.f);
        defer al.free(y);
        g.reset();
        var scale: Float(T) = 0;
        var dif: Float(T) = -1;
        try g.run(al, check.io, .n, 0, &scale, &dif);
        @memcpy(x, g.c);
        @memcpy(y, g.f);
        for ([_]isize{ 1, 2 }) |job| {
            g.reset();
            try g.run(al, check.io, .n, job, &scale, &dif);
            try testing.expect(dif > 0 and std.math.isFinite(dif));
            // by value: f80's padding bytes are not its bits
            for (x, g.c) |u, v| try testing.expect(lcheck.sameBits(T, u, v));
            for (y, g.f) |u, v| try testing.expect(lcheck.sameBits(T, u, v));
        }
    };
}

/// `σ_min(Z)` of `dtgsyl.f`'s (2), from `gesvd` in `f64` or `c64`.
fn sigmaMin(comptime T: type, al: Allocator, g: Gen(T)) !f64 {
    const C = if (comptime isComplex(T)) types.Complex64 else f64;
    const m = g.m;
    const n = g.n;
    const k = m * n;
    const z = try al.alloc(C, 4 * k * k);
    defer al.free(z);
    @memset(z, scalar(C, 0));
    const lda = m + Gen(T).ld.a;
    const ldb = n + Gen(T).ld.b;
    const ldd = m + Gen(T).ld.d;
    const lde = n + Gen(T).ld.e;
    const to = struct {
        fn f(v: T) C {
            return if (comptime isComplex(T)) .init(@floatCast(v.re), @floatCast(v.im)) else @floatCast(v);
        }
    }.f;
    for (0..n) |j| for (0..m) |i| {
        const r = i + m * j;
        for (0..m) |kk| {
            z[r + 2 * k * (kk + m * j)] = to(g.a[i + kk * lda]);
            z[k + r + 2 * k * (kk + m * j)] = to(g.d[i + kk * ldd]);
        }
        for (0..n) |kk| {
            z[r + 2 * k * (k + i + m * kk)] = ops.neg(C, to(g.b[kk + j * ldb]));
            z[k + r + 2 * k * (k + i + m * kk)] = ops.neg(C, to(g.e[kk + j * lde]));
        }
    };
    const s = try al.alloc(f64, 2 * k);
    defer al.free(s);
    var dummy: [1]C = undefined;
    try faer_svd.gesvd(al, check.io, C, .n, .n, @intCast(2 * k), @intCast(2 * k), z.ptr, @intCast(2 * k), s.ptr, &dummy, 1, &dummy, 1, null);
    return s[2 * k - 1];
}

test "relapack tgsyl: dif estimates σ_min of the Kronecker form; ijob 1 and 3, 2 and 4 agree" {
    const al = testing.allocator;
    inline for (lcheck.test_types) |T| for ([_][2]usize{ .{ 1, 1 }, .{ 2, 2 }, .{ 3, 5 }, .{ 6, 4 } }) |s| {
        const g = try Gen(T).init(al, s[0], s[1], s[0] * 11 + s[1]);
        defer g.deinit(al);
        const sig = try sigmaMin(T, al, g);
        var difs: [5]Float(T) = undefined;
        for (1..5) |job| {
            g.reset();
            var scale: Float(T) = 0;
            difs[job] = -1;
            try g.run(al, check.io, .n, @intCast(job), &scale, &difs[job]);
            // an upper bound on σ_min, within sqrt(2 m n) of it as the
            // estimate's right-hand sides are about that long
            const d: f64 = @floatCast(difs[job]);
            const lim = 4 * @sqrt(2 * @as(f64, @floatFromInt(s[0] * s[1])));
            if (!(d >= sig * (1 - 1e-3) / lim and d <= sig * lim)) {
                check.say("tgsyl {s} {d}x{d} ijob {d}: dif {e}, sigma_min {e}\n", .{ @typeName(T), s[0], s[1], job, d, sig });
                return error.TestUnexpectedResult;
            }
        }
        try testing.expectEqual(difs[1], difs[3]);
        try testing.expectEqual(difs[2], difs[4]);
    };
}

test "relapack tgsyl: the recursion's dif is the unblocked xTGSY2's, every base case's blocks counted" {
    // ReLAPACK keeps only the last base case's `PQ`: its `ijob` 2 and 4
    // estimate would come out about sqrt(base cases) too small
    const al = testing.allocator;
    inline for (.{ f64, types.Complex32 }) |T| {
        const F = Float(T);
        const g = try Gen(T).init(al, 40, 30, 3);
        defer g.deinit(al);
        for ([_]isize{ 3, 4 }) |job| {
            g.reset();
            var scale: F = 0;
            var dif: F = 0;
            try g.run(al, check.io, .n, job, &scale, &dif);
            // the whole problem as one `xTGSY2`
            g.reset();
            for (0..g.n) |j| for (0..g.m) |i| {
                g.c[i + j * (g.m + 3)] = scalar(T, 0);
                g.f[i + j * (g.m + 5)] = scalar(T, 0);
            };
            var acc: Dif(F) = .{};
            _ = try tgsy2(al, check.io, T, false, @intCast(job - 2), g.m, g.n, g.a.ptr, g.m + 1, g.b.ptr, g.n + 2, g.c.ptr, g.m + 3, g.d.ptr, g.m + 2, g.e.ptr, g.n + 1, g.f.ptr, g.m + 5, &acc);
            const num: F = @floatFromInt(if (job == 3) 2 * g.m * g.n else acc.pq);
            const want = @sqrt(num) / (acc.scale * @sqrt(acc.sum));
            try testing.expect(@abs(dif - want) <= 1e3 * std.math.floatEps(F) * want);
        }
    }
}

test "relapack tgsyl: common eigenvalues are perturbed, written and finite; scale below 1; the first bad argument" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64 }) |T| {
        const F = Float(T);
        // (A, D) = (B, E): Z is singular; a huge right-hand side scales
        const a = [_]T{ scalar(T, 2), scalar(T, 0), scalar(T, 1), scalar(T, 2) };
        const d = [_]T{ scalar(T, 1), scalar(T, 0), scalar(T, 0.5), scalar(T, 1) };
        for ([_]F{ 1, std.math.floatMax(F) / 4 }) |big| {
            var c: [4]T = @splat(scalar(T, big));
            var f: [4]T = @splat(scalar(T, big));
            var scale: F = 0;
            var dif: F = -1;
            try testing.expectError(error.Perturbed, tgsyl(al, check.io, T, .n, 0, 2, 2, &a, 2, &a, 2, &c, 2, &d, 2, &d, 2, &f, 2, &scale, &dif));
            try testing.expect(scale > 0 and scale <= 1);
            if (big > 1) try testing.expect(scale < 1);
            for (c ++ f) |x| try testing.expect(if (comptime isComplex(T)) std.math.isFinite(x.re) and std.math.isFinite(x.im) else std.math.isFinite(x));
        }
        var c: [4]T = @splat(lcheck.poisoned(T));
        var scale: F = 0;
        var dif: F = 0;
        const t: Trans = if (comptime isComplex(T)) .t else .c;
        const z = &c;
        try testing.expectError(error.InvalidTrans, tgsyl(al, check.io, T, t, 0, 2, 2, z, 2, z, 2, z, 2, z, 2, z, 2, z, 2, &scale, &dif));
        try testing.expectError(error.InvalidIjob, tgsyl(al, check.io, T, .n, 5, 2, 2, z, 2, z, 2, z, 2, z, 2, z, 2, z, 2, &scale, &dif));
        try testing.expectError(error.InvalidM, tgsyl(al, check.io, T, transOf(T, true), 9, 0, 2, z, 2, z, 2, z, 2, z, 2, z, 2, z, 2, &scale, &dif));
        try testing.expectError(error.InvalidN, tgsyl(al, check.io, T, .n, 0, 2, 0, z, 2, z, 2, z, 2, z, 2, z, 2, z, 2, &scale, &dif));
        try testing.expectError(error.InvalidLda, tgsyl(al, check.io, T, .n, 0, 2, 2, z, 1, z, 2, z, 2, z, 2, z, 2, z, 2, &scale, &dif));
        try testing.expectError(error.InvalidLdb, tgsyl(al, check.io, T, .n, 0, 2, 2, z, 2, z, 1, z, 2, z, 2, z, 2, z, 2, &scale, &dif));
        try testing.expectError(error.InvalidLdc, tgsyl(al, check.io, T, .n, 0, 2, 2, z, 2, z, 2, z, 1, z, 2, z, 2, z, 2, &scale, &dif));
        try testing.expectError(error.InvalidLdd, tgsyl(al, check.io, T, .n, 0, 2, 2, z, 2, z, 2, z, 2, z, 1, z, 2, z, 2, &scale, &dif));
        try testing.expectError(error.InvalidLde, tgsyl(al, check.io, T, .n, 0, 2, 2, z, 2, z, 2, z, 2, z, 2, z, 1, z, 2, &scale, &dif));
        try testing.expectError(error.InvalidLdf, tgsyl(al, check.io, T, .n, 0, 2, 2, z, 2, z, 2, z, 2, z, 2, z, 2, z, 1, &scale, &dif));
        for (c) |x| try testing.expect(lcheck.isNanT(T, x));
    }
}

const BigGsyl = struct {
    g: Gen(f64),

    fn run(self: BigGsyl, io: Io, out: []f64) anyerror!void {
        self.g.reset();
        var scale: f64 = 0;
        var dif: f64 = 0;
        try self.g.run(testing.allocator, io, .n, 1, &scale, &dif);
        @memcpy(out[0..self.g.c.len], self.g.c);
        @memcpy(out[self.g.c.len..][0..self.g.f.len], self.g.f);
        out[out.len - 1] = dif;
    }
};

test "relapack tgsyl: the same bits inline and threaded" {
    const al = testing.allocator;
    const g = try Gen(f64).init(al, lcheck.big_n, lcheck.big_n, 9);
    defer g.deinit(al);
    try lcheck.expectIoBits(f64, BigGsyl{ .g = g }, BigGsyl.run, g.c.len + g.f.len + 1, "tgsyl");
}
