//! `tgsyl`: ReLAPACK's recursion (`relapack/src/dtgsyl.c` and the rest,
//! OpenBLAS v0.3.34-317, MIT) over `sylvester.zig`'s `tgsy2`.

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const driver = @import("../blas/blas3/driver.zig");
const interface = @import("interface.zig");
const ops = @import("ops.zig");
const relapack = @import("relapack.zig");
const sylvester = @import("sylvester.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const isComplex = types.isComplex;
const Trans = types.Trans;
const scalar = ops.scalar;
const SylvesterSolution = sylvester.SylvesterSolution;
const scaleBlock = sylvester.scaleBlock;
const Dif = sylvester.Dif;
const tgsy2 = sylvester.tgsy2;

const Error = driver.Error;

/// What `tgsyl` returns: a `trans` it does not take, `UnsupportedJob` for an
/// `ijob` past 4 (with `trans` `n`), a slice too short or a leading
/// dimension too small; `error.OutOfMemory` or `error.Canceled`.
pub const TgsylError = interface.TransError || error{UnsupportedJob} || interface.ShapeError || Allocator.Error || Io.Cancelable;

/// What `tgsyl` found: `scale`, the `Dif` estimate (null for `ijob` 0 or a
/// transposed solve), and `perturbed` for the reference's `INFO > 0`
/// (`(A, D)` and `(B, E)` have common or close eigenvalues, so perturbed
/// values were used); `C` and `F` are written either way.
pub fn GeneralizedSylvesterSolution(comptime F: type) type {
    return struct { scale: F, separation: ?F = null, perturbed: bool };
}

/// `?TGSYL` by ReLAPACK's recursion (`relapack/src/dtgsyl.c` and the rest):
/// `A R - L B = scale C`, `D R - L E = scale F` for `R` and `L` (`m x n`),
/// into `C` and `F`, `(A, D)` and `(B, E)` in generalized Schur form;
/// `trans` `t` (real) or `c` (complex) solves `Aᴴ R + Dᴴ L = scale C`,
/// `R Bᴴ + L Eᴴ = -scale F`. The returned `scale <= 1` keeps `R` and `L`
/// finite. `ijob` (read only for `trans` `n`): 0 solves; 3 and 4 only
/// estimate `Dif[(A, D), (B, E)]` (by `xLATDF`'s look-ahead or its `xGECON`
/// null vectors), leaving that estimate's `R` and `L` in `C` and `F`; 1 and
/// 2 do both. An empty `m` or `n`, which the reference refuses, solves
/// nothing.
pub fn tgsyl(allocator: Allocator, io: Io, comptime T: type, trans: Trans, ijob: usize, m: usize, n: usize, a: []const T, lda: usize, b: []const T, ldb: usize, c: []T, ldc: usize, d: []const T, ldd: usize, e: []const T, lde: usize, f: []T, ldf: usize) TgsylError!GeneralizedSylvesterSolution(Float(T)) {
    comptime ops.requireFloat(T, "tgsyl");
    const F = Float(T);
    // `N`, or `T` for real and `C` for complex, as `dtgsyl.f` and `ztgsyl.f`
    const tr = switch (trans) {
        .n => false,
        .t => if (comptime isComplex(T)) return error.UnsupportedTranspose else true,
        .c => if (comptime isComplex(T)) true else return error.UnsupportedTranspose,
        .r, .o, .u, .s, .d => return error.UnsupportedTranspose,
    };
    if (!tr and ijob > 4) return error.UnsupportedJob;
    try util.holds(a.len, m, m, lda);
    try util.holds(b.len, n, n, ldb);
    try util.holds(c.len, m, n, ldc);
    try util.holds(d.len, m, m, ldd);
    try util.holds(e.len, n, n, lde);
    try util.holds(f.len, m, n, ldf);
    if (m == 0 or n == 0) return .{ .scale = 1, .perturbed = false };
    const um = m;
    const un = n;
    const p: Gsyl(T) = .{ .tr = tr, .a = a.ptr, .lda = lda, .b = b.ptr, .ldb = ldb, .c = c.ptr, .ldc = ldc, .d = d.ptr, .ldd = ldd, .e = e.ptr, .lde = lde, .f = f.ptr, .ldf = ldf };
    const job: usize = if (tr) 0 else ijob;
    // `ijob` 1 and 2 solve, then estimate on zeros; 3 and 4 only estimate
    var ifunc: u2 = if (job >= 3) @intCast(job - 2) else 0;
    const isolve: usize = if (job == 1 or job == 2) 2 else 1;
    if (job >= 3) zeroBoth(T, p, um, un);
    const saved = try allocator.alloc(T, if (isolve == 2) 2 * um * un else 0);
    defer allocator.free(saved);
    var perturbed = false;
    var scale: F = 1;
    var dif: ?F = null;
    var scale2: F = 1;
    for (0..isolve) |round| {
        var acc: Dif(F) = .{};
        const r = try tgsylRec(allocator, io, T, ifunc, um, un, p, &acc);
        scale = r.scale;
        if (r.perturbed) perturbed = true;
        if (acc.scale != 0) {
            const num: F = @floatFromInt(if (job == 1 or job == 3) 2 * um * un else acc.pq);
            dif = @sqrt(num) / (acc.scale * @sqrt(acc.sum));
        }
        if (isolve == 2 and round == 0) {
            ifunc = @intCast(job);
            scale2 = scale;
            copyBlock(T, um, un, c.ptr, p.ldc, saved.ptr, um);
            copyBlock(T, um, un, f.ptr, p.ldf, saved.ptr + um * un, um);
            zeroBoth(T, p, um, un);
        } else if (isolve == 2) {
            copyBlock(T, um, un, saved.ptr, um, c.ptr, p.ldc);
            copyBlock(T, um, un, saved.ptr + um * un, um, f.ptr, p.ldf);
            scale = scale2;
        }
    }
    return .{ .scale = scale, .separation = dif, .perturbed = perturbed };
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
fn tgsylRec(al: Allocator, io: Io, comptime T: type, ifunc: u2, m: usize, n: usize, p: Gsyl(T), acc: *Dif(Float(T))) Error!SylvesterSolution(Float(T)) {
    if (m <= relapack.crossover and n <= relapack.crossover) return tgsy2(al, io, T, p.tr, ifunc, m, n, p.a, p.lda, p.b, p.ldb, p.c, p.ldc, p.d, p.ldd, p.e, p.lde, p.f, p.ldf, acc);
    const one = scalar(T, 1);
    const minus = scalar(T, -1);
    const ct = relapack.ctrans(T);
    var r1: SylvesterSolution(Float(T)) = undefined;
    var r2: SylvesterSolution(Float(T)) = undefined;
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
