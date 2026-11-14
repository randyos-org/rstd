//! The table-driven races of `lib2`, `lib3` and `lapack`.

const std = @import("std");
const Io = std.Io;
const rstd = @import("rstd");
const bb = @import("blas_bench.zig");
const kernels = @import("kernels.zig");
const libs = @import("libs.zig");

const blas = rstd.linalg.blas;
const lapack = rstd.linalg.lapack;

const clock = bb.clock;
const l1_timing = bb.l1_timing;
const l3_timing = bb.l3_timing;
const Prec = bb.Prec;
const precs = bb.precs;

const axpy_alpha = kernels.axpy_alpha;
const gemm_sizes = kernels.gemm_sizes;
const gemv_alpha = kernels.gemv_alpha;
const l2_bytes = kernels.l2_bytes;

const cast = libs.cast;
const Col = libs.Col;
const fail = libs.fail;
const lib_cols = libs.lib_cols;
const lookupF = libs.lookupF;
const Num = libs.Num;
const printHeader = libs.printHeader;
const rstd_alloc = libs.rstd_alloc;
const rstd_col = libs.rstd_col;
const rstdFirst = libs.rstdFirst;
const sectionMark = libs.sectionMark;
const sumBound = libs.sumBound;

// ---- table-driven races ------------------------------------------------------
//
// The routines beyond `gemv` and `gemm` (`ext_l2`, `ext_l3`) and LAPACK
// (`lapack_ops`), one entry each, a struct of declarations:
//
// - `name`; `shape` (`Shape`); `fill` (`A`'s `Fill`, default `general`);
//   `prep` (`Prep`, default `none`); `restore` (what a call overwrites and
//   must find fresh again; default none, it accumulates); `out` (the `Buf`
//   checked, compared and hashed); `cols` (`out`'s columns where it is `b`
//   or `c` of an `mm` shape, default n); `level` (default its table's);
//   `band_rows` (a `band` shape's, default `2 band_k + 1`);
// - `flops(n)`, real; a complex call counts 4 times as many;
// - `sym(p)` (default `p.name ++ name`), `Sig(p)` and `lib(p, f, w)`: the
//   Fortran-ABI call (`lookupF`); a library without it shows `-`;
// - `rstd(p, w)`, absent until rstd has the routine; `has_rstd` gates one
//   written ahead of it;
// - `residual(p, w, i, j)`: `|want - got|` at output entry (i, j) and the
//   size of its terms, null off the stored part; without it only the
//   agreement with rstd is shown.
//
// Adding a routine is adding an entry to its table.

/// Where an entry's orders come from: `l2_bytes`, `gemm_sizes`,
/// `lapack_sizes`.
const Level = enum { l2, l3, lapack };

/// LAPACK's orders: `gemm_sizes` less the largest, where reference LAPACK
/// and `prepare` take many seconds a call.
const lapack_sizes = gemm_sizes[0 .. gemm_sizes.len - 1].*;

const Shape = enum {
    /// `A` n by n, `x` in `b` and `y` in `c`, n long; n from `l2_bytes`
    mv,
    /// `A` in band storage, `band_rows` by n; `x` and `y` as `mv`; n from
    /// `l2_bytes` over the rows
    band,
    /// `A`, `B` and `C` n by n
    mm,
};

/// `gbmv`'s and `gbtrf`'s `kl` and `ku`.
const band_k = 32;

const Fill = enum {
    /// every float in [-1, 1]
    general,
    /// `general` with 2n added to the diagonal's real parts: strictly
    /// diagonally dominant by rows and columns, so LU with partial pivoting
    /// exchanges no rows and every triangle is well conditioned
    dominant,
    /// `dominant` made Hermitian (symmetric if real), its diagonal real:
    /// positive definite
    spd,
};

const Prep = enum {
    none,
    /// `A` replaced by its Cholesky factor `U`, `A = Uᴴ U` (`prepare`)
    potrf,
    /// `A` replaced by `L` and `U` of `A = L U`, pivots 1 to n (`dominant`
    /// takes no row exchange; `prepare`)
    getrf,
};

/// `orig` is `A` before `prep`; `a0`, `b0` and `c0` what each call starts
/// from; `s` singular values, eigenvalues or `tau`.
const Buf = enum { a, a0, orig, b, b0, c, c0, s, u, vt };

const Restore = struct { a: bool = false, b: bool = false, c: bool = false };

const Cols = enum { one, n };

/// The table-driven races' contestants: the libraries' Fortran-ABI routines
/// and rstd's public ones. `max_n` holds in level 3 and LAPACK only.
pub const ext_cols = lib_cols[0..4].* ++ [_]Col{
    .{ .label = "reference", .short = "ref", .lib = .reference, .max_n = 1024 },
    .{ .label = "eigen", .short = "eigen", .lib = .eigen },
    .{ .label = "aocl-libflame", .short = "aocl", .lib = .aoclflame },
};

pub const ext_l2 = .{ Ger, Symv, Trsv, Gbmv };
pub const ext_l3 = .{ Gemm3, Symm, Syrk, Trmm, Trsm };
/// rstd's LAPACK, then the hooks: T34's SVD (faer), T111's QR and
/// eigensolver, T110's ReLAPACK set; each races the libraries now and
/// takes rstd's column once its routine lands.
const lapack_ops = .{ Potrf, Getrf, Getrs, Gesv, Trtri, Potri, Gesvd, Gesdd, Geqrf, Syevd, Sytrf, Gbtrf };

pub fn libLapack(io: Io, gpa: std.mem.Allocator) !void {
    sectionMark("lapack");
    try clock(io, "lapack start");
    try extTable(io, gpa, .lapack, lapack_ops);
}

pub fn extTable(io: Io, gpa: std.mem.Allocator, comptime lv: Level, comptime ops: anytype) !void {
    const t = if (lv == .l2) l1_timing else l3_timing;
    std.debug.print("\n{s} GFLOP/s, best of {d} batches of at least {d} ms; ratios are over OpenBLAS. {s} {s}\n\n", .{ ext_leads[@backingInt(lv)], t.batches, @divTrunc(t.min_ns, std.time.ns_per_ms), libs.io_note, ext_note });
    printHeader("routine | p | n", "------- | - | -:", &ext_cols);
    inline for (ops) |E| {
        inline for (precs) |p| try raceExt(io, gpa, E, lv, p);
    }
}

const ext_leads = [_][]const u8{
    "Level 2 beyond `gemv`, each library's Fortran-ABI routine, unit strides, n as `gemv`'s (`A` of 32 KiB, 4 MiB, 64 MiB): `ger` `A += 1e-6 x yᵀ` (`geru` complex), `symv` `y += 1e-3 S x` (upper), `trsv` upper, no transpose, non-unit on a diagonally dominant `A`, `x` restored before each call (counted), `gbmv` `y += 1e-3 A x` with kl = ku = 32 (`A`'s band of the same bytes); rates by 2n², 2n², n², 2n (kl + ku + 1), four times as many complex.",
    "Level 3 beyond `gemm`, each library's Fortran-ABI routine, square n, `alpha` and `beta` 1: `gemm` NN, `symm` left upper, `syrk` upper NN, `trmm` and `trsm` left upper NN non-unit (`trsm` on a diagonally dominant `A`), `B` restored before each `trmm` and `trsm` call (an n² copy, counted); rates by 2n³, 2n³, n²(n + 1), n³, n³, four times as many complex; `reference` to n = 1024.",
    "LAPACK, square n, upper where it matters: rates by LAPACK's counts, `potrf` n³/3, `getrf` 2n³/3, `getrs` 2n³ (n right-hand sides), `gesv` 2n³/3 + 2n² (one), `trtri` n³/3, `potri` 2n³/3, and nominal ones: `gesvd` and `gesdd` 22n³ (thin `U` and `Vᵀ`, Golub and Van Loan's Golub-Reinsch count), `geqrf` 4n³/3, `syevd` (`heevd` complex) 9n³ with vectors, `sytrf` n³/3, `gbtrf` 4nk² (kl = ku = k = 32, at `gemv`'s band sizes); four times as many complex. `A` is diagonally dominant for `getrf`, `getrs`, `gesv` and `trtri` (no row exchanges, so pivots agree), Hermitian positive definite for `potrf`, `potri`, `syevd` and `sytrf`; `getrs` and `potri` start from a factorization made here in f64. Each call restores what it overwrites first (an n² copy, counted). `gesvd`, `gesdd`, `geqrf`, `syevd`, `sytrf` and `gbtrf` race the libraries until rstd has them.",
};

const ext_note = "`rstd err` is rstd's largest residual at the checked entries (all, or 256 sampled), in eps of the type relative to the size of its terms; `rstd-<lib>` is rstd's largest difference from that library's output in eps relative to the output's largest entry, `!` past 2 (terms + 8). `rstd bits` hashes rstd's output. Eigen runs one thread; `-` where a library lacks the routine";

/// A complex value, or a real one with 0 imaginary part, in f64.
const Pair = [2]f64;

fn cadd(a: Pair, b: Pair) Pair {
    return .{ a[0] + b[0], a[1] + b[1] };
}
fn csub(a: Pair, b: Pair) Pair {
    return .{ a[0] - b[0], a[1] - b[1] };
}
fn cmul(a: Pair, b: Pair) Pair {
    return .{ a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0] };
}
fn cdiv(a: Pair, b: Pair) Pair {
    const d = b[0] * b[0] + b[1] * b[1];
    return .{ (a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d };
}
fn cconj(a: Pair) Pair {
    return .{ a[0], -a[1] };
}
/// `|re| + |im|`, the size `checkGemv` gives a complex.
fn cmag(a: Pair) f64 {
    return @abs(a[0]) + @abs(a[1]);
}
fn cdist(a: Pair, b: Pair) f64 {
    return @max(@abs(a[0] - b[0]), @abs(a[1] - b[1]));
}
/// `v` rounded to `p`'s type, as the routines get it.
fn cnum(comptime p: Prec, v: f64) Pair {
    return .{ @as(p.F, @floatCast(v)), 0 };
}

/// Fortran's `CHARACTER*1` arguments, by reference; each call also passes
/// their hidden lengths (1) last, which a library built without them
/// ignores.
const fc = struct {
    const U: u8 = 'U';
    const L: u8 = 'L';
    const N: u8 = 'N';
    const S: u8 = 'S';
    const V: u8 = 'V';
};
const Ch = *const u8;
const In = *const c_int;
const one_c: c_int = 1;
const band_kc: c_int = band_k;

/// A Fortran scalar of `p` by reference: real, or (re, im).
fn Sc(comptime p: Prec) type {
    return if (p.complex) [2]p.F else p.F;
}
fn sc(comptime p: Prec, v: f64) Sc(p) {
    return if (p.complex) .{ @floatCast(v), 0 } else @floatCast(v);
}
/// The same scalar as rstd takes it.
fn num(comptime p: Prec, v: f64) Num(p) {
    return if (p.complex) .init(@floatCast(v), 0) else @floatCast(v);
}

fn decl(comptime E: type, comptime name: []const u8, comptime default: anytype) @TypeOf(default) {
    return if (@hasDecl(E, name)) @field(E, name) else default;
}
fn hasRstd(comptime E: type) bool {
    return decl(E, "has_rstd", @hasDecl(E, "rstd"));
}
fn bandRows(comptime E: type) usize {
    return decl(E, "band_rows", @as(usize, 2 * band_k + 1));
}
fn symOf(comptime E: type, comptime p: Prec) []const u8 {
    return if (@hasDecl(E, "sym")) E.sym(p) else p.name ++ E.name;
}

/// An entry's orders at level `lv`.
fn orders(comptime E: type, comptime lv: Level, comptime p: Prec) []const usize {
    return comptime switch (lv) {
        .l2 => blk: {
            var ns: [l2_bytes.len]usize = undefined;
            for (l2_bytes, &ns) |bytes, *n| n.* = switch (E.shape) {
                .band => bytes / (p.es() * bandRows(E)),
                .mv, .mm => std.math.sqrt(bytes / p.es()),
            };
            const out = ns;
            break :blk &out;
        },
        .l3 => &gemm_sizes,
        .lapack => &lapack_sizes,
    };
}

/// Elements of `A`, and of `B` and `C`, at order n.
fn lenA(comptime E: type, n: usize) usize {
    return switch (E.shape) {
        .mv, .mm => n * n,
        .band => bandRows(E) * n,
    };
}
fn lenV(comptime E: type, n: usize) usize {
    return switch (E.shape) {
        .mv, .band => n,
        .mm => n * n,
    };
}

/// Rows and columns of `E`'s output at order n.
fn outDims(comptime E: type, n: usize) [2]usize {
    return switch (E.out) {
        .a => .{ if (E.shape == .band) bandRows(E) else n, n },
        .b, .c => if (E.shape == .mm) .{ n, if (decl(E, "cols", Cols.n) == .one) 1 else n } else .{ n, 1 },
        .s => .{ n, 1 },
        else => @compileError("an output is a, b, c or s"),
    };
}

/// Floats of `E`'s output at order n; `s` is real.
fn outLen(comptime E: type, comptime p: Prec, n: usize) usize {
    const d = outDims(E, n);
    return d[0] * d[1] * (if (E.out == .s) 1 else p.cs());
}

/// Terms of one output entry, for the bounds.
fn terms(comptime E: type, n: usize) usize {
    return if (E.shape == .band) bandRows(E) else n;
}

/// A residual past this many eps is wrong, not rounding.
fn wrongLimit(comptime E: type, n: usize) f64 {
    return 64 * @as(f64, @floatFromInt(terms(E, n) + 8));
}

/// An entry's operands at its largest order, shared by its contestants.
fn Work(comptime p: Prec) type {
    return struct {
        const W = @This();
        const F = p.F;
        const cs = p.cs();

        n: usize = 0,
        a: []align(64) F,
        a0: []align(64) F,
        orig: []align(64) F,
        b: []align(64) F,
        b0: []align(64) F,
        c: []align(64) F,
        c0: []align(64) F,
        s: []align(64) F,
        u: []align(64) F,
        vt: []align(64) F,
        work: []align(64) F,
        rwork: []align(64) F,
        ipiv: []c_int,
        ipiv0: []c_int,
        rpiv: []isize,
        rpiv0: []isize,
        iwork: []c_int,
        info: c_int = 0,

        fn init(gpa: std.mem.Allocator, comptime E: type, comptime lv: Level, n: usize) !W {
            const la = lenA(E, n) * cs;
            const lb = lenV(E, n) * cs;
            const big = lv == .lapack;
            var w: W = undefined;
            w.n = 0;
            w.info = 0;
            inline for (.{ .{ "a", la }, .{ "a0", la }, .{ "orig", la }, .{ "b", lb }, .{ "b0", lb }, .{ "c", lb }, .{ "c0", lb } }) |f|
                @field(w, f[0]) = try gpa.alignedAlloc(F, .@"64", f[1]);
            w.s = try gpa.alignedAlloc(F, .@"64", n * cs);
            const sq = if (big) n * n * cs else 0;
            w.u = try gpa.alignedAlloc(F, .@"64", sq);
            w.vt = try gpa.alignedAlloc(F, .@"64", sq);
            w.work = try gpa.alignedAlloc(F, .@"64", if (big) lwork(n) * cs else 0);
            w.rwork = try gpa.alignedAlloc(F, .@"64", if (big) lrwork(n) else 0);
            w.iwork = try gpa.alloc(c_int, if (big) liwork(n) else 0);
            w.ipiv = try gpa.alloc(c_int, n);
            w.ipiv0 = try gpa.alloc(c_int, n);
            w.rpiv = try gpa.alloc(isize, n);
            w.rpiv0 = try gpa.alloc(isize, n);
            return w;
        }

        fn deinit(w: *W, gpa: std.mem.Allocator) void {
            inline for (.{ "a", "a0", "orig", "b", "b0", "c", "c0", "s", "u", "vt", "work", "rwork", "iwork", "ipiv", "ipiv0", "rpiv", "rpiv0" }) |f| gpa.free(@field(w, f));
        }

        /// Workspace for every routine's blocked path: past each one's
        /// documented need (`gesdd`'s 4n² + 7n, `syevd`'s 2n² + 6n + 1,
        /// `n nb` with `nb` up to 64).
        fn lwork(n: usize) usize {
            return 5 * n * n + 128 * n;
        }
        fn lrwork(n: usize) usize {
            return 5 * n * n + 64 * n;
        }
        fn liwork(n: usize) usize {
            return 8 * n + 64;
        }

        fn ni(w: *const W) c_int {
            return @intCast(w.n);
        }
        fn in(w: *const W) isize {
            return @intCast(w.n);
        }
        fn lw(w: *const W) c_int {
            return @intCast(lwork(w.n));
        }
        fn lrw(w: *const W) c_int {
            return @intCast(lrwork(w.n));
        }
        fn liw(w: *const W) c_int {
            return @intCast(liwork(w.n));
        }

        fn buf(w: *const W, comptime b: Buf) []const F {
            return @field(w, @tagName(b));
        }

        /// Element k of `b`, as (re, im); `s` is real.
        fn at(w: *const W, comptime b: Buf, k: usize) Pair {
            const v = w.buf(b);
            if (b == .s or !p.complex) return .{ v[k], 0 };
            return .{ v[2 * k], v[2 * k + 1] };
        }

        fn set(v: []F, k: usize, x: Pair) void {
            if (p.complex) {
                v[2 * k] = @floatCast(x[0]);
                v[2 * k + 1] = @floatCast(x[1]);
            } else v[k] = @floatCast(x[0]);
        }

        /// `E`'s output at the current order.
        fn out(w: *const W, comptime E: type) []const F {
            return w.buf(E.out)[0..outLen(E, p, w.n)];
        }

        /// Every operand fresh, before a check.
        fn reset(w: *W, comptime E: type) void {
            const la = lenA(E, w.n) * cs;
            const lb = lenV(E, w.n) * cs;
            @memcpy(w.a[0..la], w.a0[0..la]);
            @memcpy(w.b[0..lb], w.b0[0..lb]);
            @memcpy(w.c[0..lb], w.c0[0..lb]);
            @memcpy(w.ipiv[0..w.n], w.ipiv0[0..w.n]);
            @memcpy(w.rpiv[0..w.n], w.rpiv0[0..w.n]);
            w.info = 0;
        }

        /// What a call of `E` overwrites, fresh; inside its timing.
        fn restore(w: *W, comptime E: type) void {
            const r = decl(E, "restore", Restore{});
            if (r.a) {
                const la = lenA(E, w.n) * cs;
                @memcpy(w.a[0..la], w.a0[0..la]);
            }
            const lb = lenV(E, w.n) * cs;
            if (r.b) @memcpy(w.b[0..lb], w.b0[0..lb]);
            if (r.c) @memcpy(w.c[0..lb], w.c0[0..lb]);
        }

        /// The operands at order n: random, as `E`'s fill, and prepared.
        fn setup(w: *W, gpa: std.mem.Allocator, comptime E: type, n: usize) !void {
            w.n = n;
            var rng = std.Random.DefaultPrng.init(5);
            const r = rng.random();
            const la = lenA(E, n) * cs;
            const lb = lenV(E, n) * cs;
            for (w.orig[0..la]) |*v| v.* = @floatCast(r.float(f64) * 2 - 1);
            for (w.b0[0..lb]) |*v| v.* = @floatCast(r.float(f64) * 2 - 1);
            for (w.c0[0..lb]) |*v| v.* = @floatCast(r.float(f64) * 2 - 1);
            const fill = comptime decl(E, "fill", Fill.general);
            if (fill != .general) {
                if (E.shape == .band) @compileError("a band is filled general");
                if (fill == .spd) {
                    for (0..n) |j| for (0..j) |i| set(w.orig, j + i * n, cconj(w.at(.orig, i + j * n)));
                    for (0..n) |i| set(w.orig, i + i * n, .{ w.at(.orig, i + i * n)[0], 0 });
                }
                const two_n: F = @floatFromInt(2 * n);
                for (0..n) |i| w.orig[(i + i * n) * cs] += two_n;
            }
            @memcpy(w.a0[0..la], w.orig[0..la]);
            for (w.ipiv0[0..n], w.rpiv0[0..n], 1..) |*d, *e, k| {
                d.* = @intCast(k);
                e.* = @intCast(k);
            }
            const prep = comptime decl(E, "prep", Prep.none);
            if (prep != .none) try prepare(gpa, p, w, prep);
        }
    };
}

/// `prep`'s factorization of `orig` into `a0`, unblocked in f64: `potrf`
/// writes `U` over the upper triangle, `getrf` `L` and `U` over all of it.
fn prepare(gpa: std.mem.Allocator, comptime p: Prec, w: *Work(p), comptime what: Prep) !void {
    const n = w.n;
    const m = try gpa.alloc(Pair, n * n);
    defer gpa.free(m);
    for (m, 0..) |*v, k| v.* = w.at(.orig, k);
    switch (what) {
        .none => return,
        .potrf => for (0..n) |j| {
            for (0..j + 1) |i| {
                var s = m[i + j * n];
                for (0..i) |k| s = csub(s, cmul(cconj(m[k + i * n]), m[k + j * n]));
                m[i + j * n] = if (i < j) cdiv(s, m[i + i * n]) else .{ @sqrt(s[0]), 0 };
            }
        },
        .getrf => for (0..n) |k| {
            const piv = m[k + k * n];
            for (k + 1..n) |i| m[i + k * n] = cdiv(m[i + k * n], piv);
            for (k + 1..n) |j| {
                const u = m[k + j * n];
                for (k + 1..n) |i| m[i + j * n] = csub(m[i + j * n], cmul(m[i + k * n], u));
            }
        },
    }
    for (0..n) |j| for (0..n) |i| {
        if (what == .potrf and i > j) continue;
        Work(p).set(w.a0, i + j * n, m[i + j * n]);
    };
}

/// A contestant's largest residual over `E`'s output in eps of the type
/// relative to the size of its terms: every entry of a small output, else
/// 256 sampled; null without `residual`.
fn checkExt(comptime E: type, comptime p: Prec, w: *const Work(p)) ?f64 {
    if (!@hasDecl(E, "residual")) return null;
    const eps: f64 = std.math.floatEps(p.F);
    const d = outDims(E, w.n);
    var worst: f64 = 0;
    if (d[0] * d[1] <= 4096) {
        for (0..d[1]) |j| for (0..d[0]) |i| {
            const r = E.residual(p, w, i, j) orelse continue;
            worst = @max(worst, r[0] / (eps * @max(r[1], 1e-300)));
        };
    } else {
        var rng = std.Random.DefaultPrng.init(7);
        const rand = rng.random();
        var hits: usize = 0;
        var tries: usize = 0;
        while (hits < 256 and tries < 4096) : (tries += 1) {
            const r = E.residual(p, w, rand.uintLessThan(usize, d[0]), rand.uintLessThan(usize, d[1])) orelse continue;
            hits += 1;
            worst = @max(worst, r[0] / (eps * @max(r[1], 1e-300)));
        }
    }
    return if (std.math.isNan(worst)) std.math.inf(f64) else worst;
}

/// rstd's largest difference from `got` in eps relative to its own largest
/// entry.
fn diffNorm(comptime F: type, mine: []const F, got: []const F) f64 {
    var big: f64 = 0;
    var worst: f64 = 0;
    for (mine, got) |a, b| {
        big = @max(big, @abs(@as(f64, a)));
        worst = @max(worst, @abs(@as(f64, a) - b));
    }
    const r = worst / (std.math.floatEps(F) * @max(big, 1e-300));
    return if (std.math.isNan(r)) std.math.inf(f64) else r;
}

fn ExtCall(comptime E: type, comptime p: Prec) type {
    return struct {
        /// null: rstd's
        f: ?*const E.Sig(p),
        w: *Work(p),
        fn call(s: *const @This()) void {
            s.w.restore(E);
            if (s.f) |f| return E.lib(p, f, s.w);
            if (comptime hasRstd(E)) E.rstd(p, s.w) catch |e| fail(e);
        }
        pub fn run(o: *anyopaque) void {
            const s: *@This() = @ptrCast(@alignCast(o));
            s.call();
        }
    };
}

fn raceExt(io: Io, gpa: std.mem.Allocator, comptime E: type, comptime table_level: Level, comptime p: Prec) !void {
    const lv: Level = comptime decl(E, "level", table_level);
    const ns = comptime orders(E, lv, p);
    var w = try Work(p).init(gpa, E, lv, ns[ns.len - 1]);
    defer w.deinit(gpa);
    const kout = try gpa.alloc(p.F, outLen(E, p, ns[ns.len - 1]));
    defer gpa.free(kout);
    const Sig = *const E.Sig(p);
    const sym = comptime symOf(E, p);
    const fns = [ext_cols.len]?Sig{
        lookupF(.openblas, Sig, sym),
        null,
        lookupF(.blis, Sig, sym),
        lookupF(.accelerate, Sig, sym),
        lookupF(.reference, Sig, sym),
        lookupF(.eigen, Sig, sym),
        lookupF(.aoclflame, Sig, sym),
    };
    const EC = ExtCall(E, p);
    const timing = if (lv == .l2) l1_timing else l3_timing;
    for (ns) |n| {
        if (n > libs.n_max) continue;
        try w.setup(gpa, E, n);
        const len = outLen(E, p, n);
        var es: libs.Entrants(EC, ext_cols.len) = .{};
        var ok = true;
        var have_mine = false;
        for (0..ext_cols.len) |o| {
            const i = rstdFirst(o);
            const col = ext_cols[i];
            if (!col.shown() or (lv != .l2 and n > col.max_n)) continue;
            if (i == rstd_col) {
                if (!comptime hasRstd(E)) continue;
            } else if (fns[i] == null) continue;
            const ctx = es.add(i, .{ .f = fns[i], .w = &w });
            w.reset(E);
            ctx.call();
            if (libs.lib_err) |err| return err;
            const e = checkExt(E, p, &w);
            if (w.info != 0 or (e != null and e.? > wrongLimit(E, n))) {
                std.debug.print("WRONG: {s} {s} {s} n={d}: info {d}, residual {?d:.1} eps\n", .{ E.name, p.name, col.label, n, w.info, e });
                ok = false;
            }
            const got = w.out(E);
            if (i == rstd_col) {
                @memcpy(kout[0..len], got);
                es.mine = .{ .err = e, .bits = std.hash.Wyhash.hash(0, std.mem.sliceAsBytes(got)) };
                have_mine = true;
            } else if (have_mine and col.lib != null) es.diffs[i] = .{ .eps = diffNorm(p.F, kout[0..len], got), .bound = sumBound(terms(E, n)) };
        }
        if (!ok or es.count == 0) continue;
        w.reset(E);
        const rates = try es.race(io, timing, E.flops(@floatFromInt(n)) * @as(f64, if (p.complex) 4 else 1));
        std.debug.print("| {s} | {s} | {d} | ", .{ E.name, p.name, n });
        es.print(&ext_cols, &rates);
    }
}

/// `A[i, k]` of a symmetric (`herm` false) or Hermitian matrix whose upper
/// triangle `b` stores, order n.
fn symAt(comptime p: Prec, w: *const Work(p), comptime b: Buf, comptime herm: bool, i: usize, k: usize) Pair {
    const n = w.n;
    if (i <= k) return w.at(b, i + k * n);
    const v = w.at(b, k + i * n);
    return if (herm) cconj(v) else v;
}

// ---- level 2 entries

const Ger = struct {
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
        try blas.ger(rstd_alloc, libs.rstd_io, Num(p), w.in(), w.in(), num(p, axpy_alpha), cast(p, w.b.ptr), 1, cast(p, w.c.ptr), 1, cast(p, w.a.ptr), w.in());
    }
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        const k = i + j * w.n;
        const al = cnum(p, axpy_alpha);
        const x = w.at(.b0, i);
        const y = w.at(.c0, j);
        return .{ cdist(cadd(w.at(.a0, k), cmul(al, cmul(x, y))), w.at(.a, k)), cmag(w.at(.a0, k)) + cmag(al) * cmag(x) * cmag(y) };
    }
};

const Symv = struct {
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
        try blas.symv(rstd_alloc, libs.rstd_io, Num(p), .u, w.in(), num(p, gemv_alpha), cast(p, w.a.ptr), w.in(), cast(p, w.b.ptr), 1, num(p, 1), cast(p, w.c.ptr), 1);
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

const Trsv = struct {
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
        try blas.trsv(rstd_alloc, libs.rstd_io, Num(p), .u, .n, .n, w.in(), cast(p, w.a.ptr), w.in(), cast(p, w.b.ptr), 1);
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

const Gbmv = struct {
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
        try blas.gbmv(rstd_alloc, libs.rstd_io, Num(p), .n, w.in(), w.in(), band_k, band_k, num(p, gemv_alpha), cast(p, w.a.ptr), @intCast(bandRows(@This())), cast(p, w.b.ptr), 1, num(p, 1), cast(p, w.c.ptr), 1);
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

// ---- level 3 entries

/// `c0 + A B` at (i, j) against `c`, `A` read by `a(w, i, k)`.
fn productAt(comptime p: Prec, w: *const Work(p), i: usize, j: usize, comptime a: fn (*const Work(p), usize, usize) Pair) [2]f64 {
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
fn SigTri(comptime p: Prec) type {
    return fn (Ch, Ch, Ch, Ch, In, In, *const Sc(p), [*]p.F, In, [*]p.F, In, usize, usize, usize, usize) callconv(.c) void;
}

/// `trmm`'s and `trsm`'s call: left, upper, no transpose, non-unit, `alpha` 1.
fn libTri(comptime p: Prec, f: *const SigTri(p), w: *Work(p)) void {
    const n = w.ni();
    const one = sc(p, 1);
    f(&fc.L, &fc.U, &fc.N, &fc.N, &n, &n, &one, w.a.ptr, &n, w.b.ptr, &n, 1, 1, 1, 1);
}

/// Row i of upper triangular `A` (`a0`) times column j of `x`, and the
/// size of its terms.
fn triRow(comptime p: Prec, w: *const Work(p), comptime x: Buf, i: usize, j: usize) [2]Pair {
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

const Gemm3 = struct {
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
        try blas.gemm(rstd_alloc, libs.rstd_io, Num(p), .n, .n, n, n, n, num(p, 1), cast(p, w.a.ptr), n, cast(p, w.b.ptr), n, num(p, 1), cast(p, w.c.ptr), n);
    }
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        return productAt(p, w, i, j, struct {
            fn a(v: *const Work(p), r: usize, k: usize) Pair {
                return v.at(.a0, r + k * v.n);
            }
        }.a);
    }
};

const Symm = struct {
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
        try blas.symm(rstd_alloc, libs.rstd_io, Num(p), .l, .u, n, n, num(p, 1), cast(p, w.a.ptr), n, cast(p, w.b.ptr), n, num(p, 1), cast(p, w.c.ptr), n);
    }
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        return productAt(p, w, i, j, struct {
            fn a(v: *const Work(p), r: usize, k: usize) Pair {
                return symAt(p, v, .a0, false, r, k);
            }
        }.a);
    }
};

const Syrk = struct {
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
        try blas.syrk(rstd_alloc, libs.rstd_io, Num(p), .u, .n, n, n, num(p, 1), cast(p, w.a.ptr), n, num(p, 1), cast(p, w.c.ptr), n);
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

const Trmm = struct {
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
        try blas.trmm(rstd_alloc, libs.rstd_io, Num(p), .l, .u, .n, .n, n, n, num(p, 1), cast(p, w.a.ptr), n, cast(p, w.b.ptr), n);
    }
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        const r = triRow(p, w, .b0, i, j);
        return .{ cdist(r[0], w.at(.b, i + j * w.n)), r[1][0] };
    }
};

const Trsm = struct {
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
        try blas.trsm(rstd_alloc, libs.rstd_io, Num(p), .l, .u, .n, .n, n, n, num(p, 1), cast(p, w.a.ptr), n, cast(p, w.b.ptr), n);
    }
    /// `A X` against the right-hand sides.
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        const r = triRow(p, w, .b, i, j);
        const b = w.at(.b0, i + j * w.n);
        return .{ cdist(r[0], b), r[1][0] + cmag(b) };
    }
};

// ---- LAPACK entries

/// `orig X` at (i, j), over k in `lo..hi`, against the identity; `X` read
/// by `x(w, k, j)`: the residual of an inverse.
fn inverseAt(comptime p: Prec, w: *const Work(p), i: usize, j: usize, lo: usize, hi: usize, comptime x: fn (*const Work(p), usize, usize) Pair) [2]f64 {
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
fn solveAt(comptime p: Prec, w: *const Work(p), i: usize, j: usize) [2]f64 {
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
fn SigUplo(comptime p: Prec) type {
    return fn (Ch, In, [*]p.F, In, *c_int, usize) callconv(.c) void;
}

fn libUplo(comptime p: Prec, f: *const SigUplo(p), w: *Work(p)) void {
    const n = w.ni();
    f(&fc.U, &n, w.a.ptr, &n, &w.info, 1);
}

const Potrf = struct {
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
        try lapack.potrf(rstd_alloc, libs.rstd_io, Num(p), .u, w.in(), cast(p, w.a.ptr), w.in(), null);
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

const Getrf = struct {
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
        try lapack.getrf(rstd_alloc, libs.rstd_io, Num(p), w.in(), w.in(), cast(p, w.a.ptr), w.in(), w.rpiv.ptr, null);
        for (w.ipiv[0..w.n], w.rpiv[0..w.n]) |*d, s| d.* = @intCast(s);
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

const Getrs = struct {
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
        try lapack.getrs(rstd_alloc, libs.rstd_io, Num(p), .n, n, n, cast(p, w.a.ptr), n, w.rpiv.ptr, cast(p, w.b.ptr), n);
    }
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        return solveAt(p, w, i, j);
    }
};

const Gesv = struct {
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
        try lapack.gesv(rstd_alloc, libs.rstd_io, Num(p), n, 1, cast(p, w.a.ptr), n, w.rpiv.ptr, cast(p, w.b.ptr), n, null);
    }
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        return solveAt(p, w, i, j);
    }
};

const Trtri = struct {
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
        try lapack.trtri(rstd_alloc, libs.rstd_io, Num(p), .u, .n, w.in(), cast(p, w.a.ptr), w.in(), null);
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

const Potri = struct {
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
        try lapack.potri(rstd_alloc, libs.rstd_io, Num(p), .u, w.in(), cast(p, w.a.ptr), w.in(), null);
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
fn svdAt(comptime p: Prec, w: *const Work(p), i: usize, j: usize) [2]f64 {
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

/// T34 (faer): thin `U` and `Vᵀ` (`S`).
const Gesvd = struct {
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
        try lapack.gesvd(rstd_alloc, libs.rstd_io, Num(p), .s, .s, n, n, cast(p, w.a.ptr), n, w.s.ptr, cast(p, w.u.ptr), n, cast(p, w.vt.ptr), n, null);
    }
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        return svdAt(p, w, i, j);
    }
};

/// T34 (faer): `jobz` `S`.
const Gesdd = struct {
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
        try lapack.gesdd(rstd_alloc, libs.rstd_io, Num(p), .s, n, n, cast(p, w.a.ptr), n, w.s.ptr, cast(p, w.u.ptr), n, cast(p, w.vt.ptr), n);
    }
    pub fn residual(comptime p: Prec, w: *const Work(p), i: usize, j: usize) ?[2]f64 {
        return svdAt(p, w, i, j);
    }
};

/// T111 (faer): add `rstd` when its `geqrf` lands; `tau` in `s`.
const Geqrf = struct {
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

/// T111 (faer): `syevd`, `heevd` complex, with vectors; add `rstd` when its
/// eigensolver lands.
const Syevd = struct {
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

/// T110 (ReLAPACK): Bunch-Kaufman, symmetric (not Hermitian) complex too;
/// no residual, which would need the pivots' blocks.
const Sytrf = struct {
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
        try lapack.sytrf(rstd_alloc, libs.rstd_io, Num(p), .u, w.in(), cast(p, w.a.ptr), w.in(), w.rpiv.ptr, null);
    }
};

/// T110 (ReLAPACK): kl = ku = `band_k` in LAPACK's band storage, kl rows
/// for the fill, at `gemv`'s band sizes; no residual.
const Gbtrf = struct {
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
        try lapack.gbtrf(rstd_alloc, libs.rstd_io, Num(p), w.in(), w.in(), band_k, band_k, cast(p, w.a.ptr), band_rows, w.rpiv.ptr, null);
    }
};
