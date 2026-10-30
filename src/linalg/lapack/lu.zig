//! LU factorization with partial pivoting, and solve. `getrf` is OpenBLAS's
//! `interface/lapack/getrf.c` and `zgetrf.c` over
//! `lapack/getrf/getrf_single.c`, whose base case is `lapack/getf2/getf2_k.c`
//! and `zgetf2_k.c` (OpenBLAS v0.3.34-317). The semantics are
//! Reference-LAPACK's `xGETRF2` (`lapack-netlib/SRC/zgetrf2.f` in that
//! checkout).
//!
//! Where they differ, the reference wins:
//! - a nonzero pivot below the safe minimum (denormal, or NaN) is swapped in
//!   and divided by, where upstream's `getf2` neither swaps nor scales yet
//!   records the swap in `ipiv`, so `P L U` is not `A`;
//! - complex: that test is on `|p|`, as `zgetrf2.f`'s `ABS`, where upstream
//!   tests each part.
//!
//! The updates are the level-3 drivers' `trsm` and `gemm`, called directly as
//! upstream calls its kernels, so Accelerate never runs them; upstream fuses
//! the interchanges, `trsm` and `gemm` per `GEMM_R` columns, the same
//! arithmetic. Its interchanges are `laswp_plus` then `gemm_oncopy` (the
//! `#if 1`), not `laswp_ncopy`, and nothing calls `neg_tcopy`; nor here.
//! `getf2`'s `dot`, `gemv`, `swap` and `scal` are loops, its `iamax` the
//! core's kernel: the pivot is the first NaN, else the first largest
//! `|re| + |im|`, as `izamax`.
//!
//! Threaded, `getrf` is `getrf_parallel.c`'s lookahead: the caller updates
//! and factors the next block while the other workers update the columns
//! past it. Upstream sizes its blocks by its thread count (`FORMULA1`,
//! `FORMULA2`); here they are `getrf_single.c`'s, so the bits are the serial
//! ones at any worker count. `getrf_parallel.c`'s `laswp_ncopy` path is
//! dead (`if (0 && ...)`); its closing interchanges, threaded there
//! (`blas_level1_thread`), are serial here.
//!
//! `getrs` is `lapack/getrs/getrs_single.c` and `zgetrs_single.c`, as
//! Reference-LAPACK's `xGETRS`: the interchanges and two triangular solves
//! through the drivers, `trsv` for one real right-hand side as upstream.
//! `getrs_parallel.c` runs that per column block on each thread; here the
//! drivers' `trsm` splits the columns alike. `gesv` is
//! `interface/lapack/gesv.c`: `getrf`, then `getrs` unless `U` is singular,
//! as `xGESV`.

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const portable = @import("../blas/portable.zig");
const threading = @import("../blas/threading.zig");
const driver = @import("../blas/blas3/driver.zig");
const triangular = @import("../blas/blas2/triangular.zig");
const laswp = @import("../blas/ext/laswp/kernel.zig");
const Vtable = @import("../blas/Vtable.zig");
const chol = @import("cholesky.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Trans = types.Trans;
const Float = types.Float;
const isComplex = types.isComplex;
const scalar = chol.scalar;
const isZero = chol.isZero;
const add = chol.add;
const sub = chol.sub;
const mul = chol.mul;
const Work = driver.Work;

/// What `getrf` returns: `Singular` when `U` has a zero on its diagonal.
pub const GetrfError = error{ InvalidM, InvalidN, InvalidLda, Singular } || Allocator.Error || Io.Cancelable;

/// `?GETRF`: `A = P L U` for the `m x n` `A`, in place: `L` unit lower
/// trapezoidal below the diagonal, `U` upper trapezoidal on and above it.
/// Row `k` (1-based) was interchanged with row `ipiv[k - 1]`, for the
/// `min(m, n)` entries of `ipiv`.
///
/// When `U(i, i)` is exactly zero, the factorization still completes, as the
/// reference's does, then returns `error.Singular` with `pivot.*` set to the
/// first such `i` (the reference's `INFO`); `A` and `ipiv` hold the whole
/// factorization. `pivot` is written only then. A NaN pivot is not zero:
/// its column below fills with NaN, and nothing fails.
pub fn getrf(allocator: Allocator, io: Io, comptime T: type, m: isize, n: isize, a: [*]T, lda: isize, ipiv: [*]isize, pivot: ?*usize) GetrfError!void {
    comptime chol.floatType(T, "getrf");
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, m)) return error.InvalidLda;
    if (m == 0 or n == 0) return;
    if (!util.fitsMatrix(m, n, lda, @sizeOf(T))) return error.InvalidLda;
    const want = chol.threads(.getrf, T, @as(u64, @intCast(m)) * @as(u64, @intCast(n))) orelse 1;
    const k = try factorWith(allocator, io, T, @intCast(m), @intCast(n), a, @intCast(lda), ipiv, want);
    return singular(k, pivot);
}

/// `error.Singular` for a nonzero `k`, written to `pivot`.
fn singular(k: usize, pivot: ?*usize) error{Singular}!void {
    if (k == 0) return;
    if (pivot) |p| p.* = k;
    return error.Singular;
}

/// The factorization on up to `want` workers: the first zero pivot, or 0.
fn factorWith(allocator: Allocator, io: Io, comptime T: type, m: usize, n: usize, a: [*]T, lda: usize, ipiv: [*]isize, want: u32) driver.Error!usize {
    const vt = chol.table(T);
    const blk = driver.Blocking.of(T, vt);
    if (want > 1) return factorAhead(allocator, io, T, vt, blk, m, n, a, lda, ipiv, want);
    return factor(allocator, io, T, vt, blk.un, blk.q, m, n, a, lda, ipiv);
}

/// `getrf_single.c`'s block width for `min(m, n) = mn`: half of it rounded up
/// to `unroll_n`, at most `Q`; at most `2 * unroll_n` is `getf2`'s
/// (upstream's POWER8 build stops at `unroll_n`).
fn blockWidth(mn: usize, un: usize, q: usize) ?usize {
    const blocking = @min((mn / 2 + un - 1) / un * un, q);
    return if (blocking <= 2 * un) null else blocking;
}

/// `getrf_single.c`: the first zero pivot's index, or 0. Pivots are 1-based
/// from `a`'s first row. Blocks of `blockWidth`, each factored recursively.
fn factor(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, un: usize, q: usize, m: usize, n: usize, a: [*]T, lda: usize, ipiv: [*]isize) driver.Error!usize {
    const mn = @min(m, n);
    const blocking = blockWidth(mn, un, q) orelse return getf2(T, vt, m, n, a, lda, ipiv);
    const one = scalar(T, 1);
    const minus = scalar(T, -1);
    var info: usize = 0;
    var j: usize = 0;
    while (j < mn) : (j += blocking) {
        const jb = @min(mn - j, blocking);
        const d = a + j * (lda + 1);
        const k = try factor(allocator, io, T, vt, un, q, m - j, jb, d, lda, ipiv + j);
        if (k != 0 and info == 0) info = k + j;
        for (ipiv[j..][0..jb]) |*p| p.* += @intCast(j);
        const rest = n - j - jb;
        if (rest == 0) continue;
        // A12 := L11^-1 P A12; A22 -= A21 A12
        const c = a + (j + jb) * lda;
        laswp.swapRows(T, rest, c, lda, j + 1, j + jb, ipiv + j, 1);
        try driver.trsm(allocator, io, T, vt, true, false, false, false, true, jb, rest, one, d, lda, c + j, lda, 1);
        if (m > j + jb) try driver.gemm(allocator, io, T, vt, false, false, .n, m - j - jb, rest, jb, minus, d + jb, lda, c + j, lda, one, c + j + jb, lda, 1);
    }
    swapBack(T, mn, blocking, a, lda, ipiv);
    return info;
}

/// The later interchanges, into each block's columns.
fn swapBack(comptime T: type, mn: usize, blocking: usize, a: [*]T, lda: usize, ipiv: [*]const isize) void {
    var j: usize = 0;
    while (j < mn) : (j += blocking) {
        const jb = @min(mn - j, blocking);
        laswp.swapRows(T, jb, a + j * lda, lda, j + jb + 1, mn, ipiv + j + jb, 1);
    }
}

/// `getrf_parallel.c` over `factor`'s blocks: the first block recursively,
/// then for each block its interchanges, `trsm` and `gemm` on the columns
/// right of it, which the caller does for the next block's columns and then
/// factors that block (`GETRF_SINGLE`), while the other workers do the rest
/// (`inner_basic_thread`). Every column gets `factor`'s kernel calls, cut at
/// whole `unroll_n` panels, so the bits are `factor`'s.
fn factorAhead(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, blk: driver.Blocking, m: usize, n: usize, a: [*]T, lda: usize, ipiv: [*]isize, want: u32) driver.Error!usize {
    const mn = @min(m, n);
    const blocking = blockWidth(mn, blk.un, blk.q) orelse return getf2(T, vt, m, n, a, lda, ipiv);
    var jb = @min(mn, blocking);
    var info = try factorAhead(allocator, io, T, vt, blk, m, jb, a, lda, ipiv, want);
    var j: usize = 0;
    while (j < mn) {
        for (ipiv[j..][0..jb]) |*p| p.* += @intCast(j);
        const next = @min(mn - j - jb, blocking);
        if (n > j + jb) {
            var step: Step(T) = .{ .vt = vt, .blk = blk, .m = m, .a = a, .lda = lda, .ipiv = ipiv, .j = j, .jb = jb, .next = next, .allocator = allocator };
            const k = try step.run(io, n - j - jb, want);
            if (k != 0 and info == 0) info = k + j + jb;
        }
        j += jb;
        jb = next;
    }
    swapBack(T, mn, blocking, a, lda, ipiv);
    return info;
}

/// One step of `factorAhead`: block `j` (`jb` wide) applied to the columns
/// right of it, and the next block (`next` wide) factored. Part 0 is the
/// caller's: the next block's columns, rounded up to whole `unroll_n`
/// panels; the rest is cut into chunks at panel edges, which the other
/// workers claim, and the caller too once its block is done. Each part has
/// packing space of its own, allocated before any task starts.
fn Step(comptime T: type) type {
    return struct {
        const Self = @This();

        vt: *const Vtable,
        blk: driver.Blocking,
        m: usize,
        a: [*]T,
        lda: usize,
        ipiv: [*]isize,
        j: usize,
        jb: usize,
        next: usize,
        allocator: Allocator,
        /// Part 0's width, the chunks of the rest, and where each part's
        /// packing space starts in `buf` (part `i + 1` is chunk `i`).
        w0: usize = 0,
        chunks: threading.Chunks = undefined,
        at: [threading.chunks_max + 2]usize = undefined,
        buf: [*]T = undefined,
        claimed: std.atomic.Value(usize) = .init(0),
        /// The next block's first zero pivot (relative), and the caller's
        /// allocation failure there.
        info: usize = 0,
        failed: bool = false,

        fn run(s: *Self, io: Io, rest: usize, want: u32) driver.Error!usize {
            const un = s.blk.un;
            s.w0 = @min(rest, (s.next + un - 1) / un * un);
            s.chunks = .init(rest - s.w0, want - 1, .even, un);
            const parts = 1 + if (rest > s.w0) s.chunks.count else 0;
            var total: usize = 0;
            for (0..parts) |i| {
                s.at[i] = total;
                total += Work(T).up(s.space(s.width(i)));
            }
            const buf = try s.allocator.alignedAlloc(T, .fromByteUnits(Work(T).alignment), @max(total, 1));
            defer s.allocator.free(buf);
            s.buf = buf.ptr;
            if (parts == 1) {
                try s.work(io, 0, null);
            } else {
                try threading.team(io, s, want);
            }
            if (s.failed) return error.OutOfMemory;
            return s.info;
        }

        /// Part `i`'s first column and width, from the first column right of
        /// block `j`.
        fn first(s: *const Self, i: usize) usize {
            return if (i == 0) 0 else s.w0 + s.chunks.get(i - 1).start;
        }

        fn width(s: *const Self, i: usize) usize {
            if (i == 0) return s.w0;
            const c = s.chunks.get(i - 1);
            return c.end - c.start;
        }

        /// Elements of packing space for `cols` columns: the `trsm`'s, then
        /// the `gemm`'s, one after the other.
        fn space(s: *const Self, cols: usize) usize {
            if (cols == 0) return 0;
            const t = Work(T).total(driver.triLens(s.blk, s.jb, cols, s.jb));
            const g = if (s.m > s.j + s.jb) Work(T).total(Work(T).lens(s.blk, s.m - s.j - s.jb, cols, s.jb, false)) else 0;
            return @max(t, g);
        }

        /// `team`'s plan: every worker that started takes part.
        pub fn plan(_: *Self, started: u32) u32 {
            return started;
        }

        /// Worker 0, the caller, does part 0 and the next block, then claims
        /// chunks as the others do.
        pub fn work(s: *Self, io: Io, i: u32, lent: ?*threading.Lent) Io.Cancelable!void {
            _ = lent;
            if (i == 0) {
                s.update(0);
                if (s.next > 0) {
                    const off = s.j + s.jb;
                    const d = s.a + off * (s.lda + 1);
                    s.info = factor(s.allocator, io, T, s.vt, s.blk.un, s.blk.q, s.m - off, s.next, d, s.lda, s.ipiv + off) catch |e| switch (e) {
                        error.Canceled => return error.Canceled,
                        error.OutOfMemory => {
                            s.failed = true;
                            return;
                        },
                    };
                }
            }
            const count = if (s.chunks.bounds[s.chunks.count] > 0) s.chunks.count else 0;
            while (true) {
                try threading.checkCancel(io);
                const c = s.claimed.fetchAdd(1, .monotonic);
                if (c >= count) return;
                s.update(c + 1);
            }
        }

        /// `inner_basic_thread` on part `i`: block `j`'s interchanges, then
        /// `A12 := L11^-1 A12` and `A22 -= A21 A12` on its columns, with
        /// `factor`'s driver loops on the part's own packing space.
        fn update(s: *const Self, i: usize) void {
            const cols = s.width(i);
            if (cols == 0) return;
            const j = s.j;
            const jb = s.jb;
            const lda = s.lda;
            const d = s.a + j * (lda + 1);
            const b = s.a + (j + jb + s.first(i)) * lda;
            const base = s.buf + s.at[i];
            laswp.swapRows(T, cols, b, lda, j + 1, j + jb, s.ipiv + j, 1);
            const o: driver.Tri = .{ .left = true, .upper = false, .trans = false, .conj = false, .unit = true };
            driver.triDrive(T, s.vt, true, o, s.blk, jb, cols, d, lda, null, b + j, lda, Work(T).at(base, driver.triLens(s.blk, jb, cols, jb)));
            if (s.m <= j + jb) return;
            const rows = s.m - j - jb;
            const work_g = Work(T).at(base, Work(T).lens(s.blk, rows, cols, jb, false));
            driver.drive(T, s.vt, .{ .rows = .n, .cols = .n }, s.blk, rows, cols, jb, scalar(T, -1), .{ .p = d + jb, .ld = lda }, .{ .p = b + j, .ld = lda }, null, b + j + jb, lda, work_g);
        }
    };
}

/// What `getrs` returns: `InvalidIpiv` for a pivot that is not a row of `B`.
pub const GetrsError = error{ InvalidTrans, InvalidN, InvalidNrhs, InvalidLda, InvalidLdb, InvalidIpiv } || Allocator.Error || Io.Cancelable;

/// `?GETRS`: solves `op(A) X = B` for the `n x nrhs` `X`, in place of `B`,
/// given `getrf`'s `P L U` of `A` and its `ipiv`. `trans` is `n`, `t`, `r`
/// (`conj(A)`) or `c` (`Aᴴ`); real `r` and `c` are `n` and `t`, as in the
/// level-3 routines, where the reference takes no `r`. No singularity check,
/// as in the reference. Every pivot must lie in `1..n`, else
/// `error.InvalidIpiv` before `B` is touched; the reference checks none.
pub fn getrs(allocator: Allocator, io: Io, comptime T: type, trans: Trans, n: isize, nrhs: isize, a: [*]const T, lda: isize, ipiv: [*]const isize, b: [*]T, ldb: isize) GetrsError!void {
    comptime chol.floatType(T, "getrs");
    if (trans.conjX()) return error.InvalidTrans;
    if (n < 0) return error.InvalidN;
    if (nrhs < 0) return error.InvalidNrhs;
    if (lda < @max(1, n)) return error.InvalidLda;
    if (ldb < @max(1, n)) return error.InvalidLdb;
    if (n == 0 or nrhs == 0) return;
    if (!util.fitsMatrix(n, n, lda, @sizeOf(T))) return error.InvalidLda;
    if (!util.fitsMatrix(n, nrhs, ldb, @sizeOf(T))) return error.InvalidLdb;
    for (ipiv[0..@intCast(n)]) |p| if (p < 1 or p > n) return error.InvalidIpiv;
    const m: usize = @intCast(n);
    const want = chol.threads(.getrs, T, m) orelse 1;
    try solve(allocator, io, T, trans, m, @intCast(nrhs), a, @intCast(lda), ipiv, b, @intCast(ldb), want);
}

/// `getrs_single.c`, `zgetrs_single.c`, the `trsm`s on up to `want` workers
/// (`getrs_parallel.c`).
fn solve(allocator: Allocator, io: Io, comptime T: type, trans: Trans, m: usize, r: usize, a: [*]const T, la: usize, ipiv: [*]const isize, b: [*]T, lb: usize, want: u32) driver.Error!void {
    const vt = chol.table(T);
    const t = trans.transposed();
    // `L Y = P B`, then `U X = Y`; transposed `op(U) Y = B`, then
    // `op(L) Z = Y` and `X = Pᵀ Z`: the first solve is upper when transposed,
    // and `L` is the unit one
    if (!t) laswp.swapRows(T, r, b, lb, 1, m, ipiv, 1);
    for ([_]bool{ t, !t }) |upper| {
        if (!isComplex(T) and r == 1) {
            try triangular.trsvTable(allocator, io, T, vt, if (upper) .u else .l, trans, if (upper) .n else .u, m, a, la, b);
        } else {
            try driver.trsm(allocator, io, T, vt, true, upper, t, trans.conjA(), !upper, m, r, scalar(T, 1), a, la, b, lb, want);
        }
    }
    if (t) laswp.swapRows(T, r, b, lb, 1, m, ipiv, -1);
}

/// What `gesv` returns: `Singular` when `U` has a zero on its diagonal.
pub const GesvError = error{ InvalidN, InvalidNrhs, InvalidLda, InvalidLdb, Singular } || Allocator.Error || Io.Cancelable;

/// `?GESV`: solves `A X = B` for the `n x nrhs` `X`, in place of `B`, by
/// `getrf` (into `A` and `ipiv`) and `getrs`. When `U(i, i)` is exactly
/// zero, returns `error.Singular` with `i` in `pivot`, the factorization
/// complete and `B` untouched, as the reference's `INFO`.
pub fn gesv(allocator: Allocator, io: Io, comptime T: type, n: isize, nrhs: isize, a: [*]T, lda: isize, ipiv: [*]isize, b: [*]T, ldb: isize, pivot: ?*usize) GesvError!void {
    comptime chol.floatType(T, "gesv");
    if (n < 0) return error.InvalidN;
    if (nrhs < 0) return error.InvalidNrhs;
    if (lda < @max(1, n)) return error.InvalidLda;
    if (ldb < @max(1, n)) return error.InvalidLdb;
    if (n == 0) return;
    if (!util.fitsMatrix(n, n, lda, @sizeOf(T))) return error.InvalidLda;
    if (!util.fitsMatrix(n, nrhs, ldb, @sizeOf(T))) return error.InvalidLdb;
    const m: usize = @intCast(n);
    const la: usize = @intCast(lda);
    const want = chol.threads(.gesv, T, @as(u64, m) * m) orelse 1;
    try singular(try factorWith(allocator, io, T, m, m, a, la, ipiv, want), pivot);
    if (nrhs == 0) return;
    try solve(allocator, io, T, .n, m, @intCast(nrhs), a, la, ipiv, b, @intCast(ldb), want);
}

/// `getf2_k.c`, `zgetf2_k.c`: left-looking, a column at a time. Column `j`
/// takes the earlier interchanges, then `U(:j, j)` by forward substitution,
/// then `A(j:, j) -= L(j:, :j) U(:j, j)`, its pivot, and `L(j + 1:, j)`.
fn getf2(comptime T: type, vt: *const Vtable, m: usize, n: usize, a: [*]T, lda: usize, ipiv: [*]isize) usize {
    const iamax = portable.kernel(T, vt, .iamax_k);
    var info: usize = 0;
    for (0..n) |j| {
        const b = a + j * lda;
        const len = @min(j, m);
        for (0..len) |i| {
            const ip: usize = @intCast(ipiv[i] - 1);
            if (ip != i) std.mem.swap(T, &b[i], &b[ip]);
        }
        for (1..@max(1, len)) |i| {
            var s = scalar(T, 0);
            for (0..i) |k| s = add(T, s, mul(T, a[i + k * lda], b[k]));
            b[i] = sub(T, b[i], s);
        }
        if (j >= m) continue;
        for (0..j) |k| {
            const x = b[k];
            const l = a + k * lda;
            for (j..m) |i| b[i] = sub(T, b[i], mul(T, l[i], x));
        }
        const p = j + @as(usize, @intCast(iamax(@intCast(m - j), @ptrCast(b + j), 1))) - 1;
        ipiv[j] = @intCast(p + 1);
        const piv = b[p];
        if (isZero(T, piv)) {
            if (info == 0) info = j + 1;
            continue;
        }
        if (p != j) for (0..j + 1) |k| std.mem.swap(T, &a[j + k * lda], &a[p + k * lda]);
        scaleBelow(T, b[j + 1 .. m], piv);
    }
    return info;
}

/// `x /= p`, as `xGETRF2`'s base case: times `1 / p` where `|p|` is at least
/// the safe minimum, else divided by `p`, whose reciprocal may overflow.
/// Complex quotients by Smith's algorithm, as `util.cdiv` and upstream.
fn scaleBelow(comptime T: type, x: []T, p: T) void {
    const F = Float(T);
    const sfmin = std.math.floatMin(F);
    if (comptime isComplex(T)) {
        if (std.math.hypot(p.re, p.im) >= sfmin) {
            const r = util.cdiv(F, 1, 0, p.re, p.im);
            for (x) |*v| {
                const z = util.cmul(F, r[0], r[1], v.re, v.im);
                v.* = .init(z[0], z[1]);
            }
        } else for (x) |*v| {
            const z = util.cdiv(F, v.re, v.im, p.re, p.im);
            v.* = .init(z[0], z[1]);
        }
    } else if (@abs(p) >= sfmin) {
        const r = 1 / p;
        for (x) |*v| v.* *= r;
    } else for (x) |*v| {
        v.* /= p;
    }
}

// ---- tests ---------------------------------------------------------------------

const testing = std.testing;
const blas = @import("../blas/__root__.zig");
const check = @import("../blas/_test/blas1/check.zig");
const guard = @import("../blas/_test/blas2/guard.zig");

const W = chol.W;
const test_types = chol.test_types;
const poisoned = chol.poisoned;
const isNanT = chol.isNanT;
const fromF64 = chol.fromF64;

const nan = std.math.nan(f64);
const inf = std.math.inf(f64);

/// A random `m x n` `A` in [-2, 2), stored with `lda = m + 2`, the padding
/// rows NaN; `aw` keeps it in `W(T)`.
fn Problem(comptime T: type) type {
    return struct {
        m: usize,
        n: usize,
        lda: usize,
        a: []T,
        aw: []W(T),
        ipiv: []isize,

        const Self = @This();

        fn init(al: Allocator, m: usize, n: usize, seed: u64) !Self {
            const lda = m + 2;
            const a = try al.alloc(T, @max(1, lda * n));
            errdefer al.free(a);
            const aw = try al.alloc(W(T), m * n);
            errdefer al.free(aw);
            const ipiv = try al.alloc(isize, @max(1, @min(m, n)));
            const v = try al.alloc(T, m * n);
            defer al.free(v);
            check.fill(T, v, seed);
            for (aw, v) |*w, x| w.* = check.wide(T, x);
            var self: Self = .{ .m = m, .n = n, .lda = lda, .a = a, .aw = aw, .ipiv = ipiv };
            self.store();
            return self;
        }

        /// `aw` into `a`, the padding poisoned; `ipiv` poisoned too.
        fn store(self: *Self) void {
            @memset(self.a, poisoned(T));
            for (0..self.n) |j| for (0..self.m) |i| {
                self.a[i + j * self.lda] = chol.toT(T, self.aw[i + j * self.m]);
            };
            @memset(self.ipiv, -1);
        }

        fn deinit(self: *Self, al: Allocator) void {
            al.free(self.a);
            al.free(self.aw);
            al.free(self.ipiv);
        }

        fn factorIt(self: *Self, al: Allocator, pivot: ?*usize) GetrfError!void {
            return getrf(al, check.io, T, @intCast(self.m), @intCast(self.n), self.a.ptr, @intCast(self.lda), self.ipiv.ptr, pivot);
        }

        fn l(self: *const Self, i: usize, k: usize) W(T) {
            return if (i == k) check.wide(T, scalar(T, 1)) else check.wide(T, self.a[i + k * self.lda]);
        }

        fn u(self: *const Self, k: usize, j: usize) W(T) {
            return check.wide(T, self.a[k + j * self.lda]);
        }

        /// `A` with `ipiv`'s interchanges applied, in `W(T)`.
        fn permuted(self: *const Self, al: Allocator) ![]W(T) {
            const pa = try al.dupe(W(T), self.aw);
            for (0..@min(self.m, self.n)) |k| {
                const p: usize = @intCast(self.ipiv[k] - 1);
                for (0..self.n) |j| std.mem.swap(W(T), &pa[k + j * self.m], &pa[p + j * self.m]);
            }
            return pa;
        }

        /// Pivots inside their rows, `|L|` at most 1 (complex: `sqrt 2`, the
        /// pivot being largest in `|re| + |im|`), `L U` against `P A` within
        /// `2 (mn + 2)` ulps of `|L| |U|`, and the padding untouched.
        fn expectFactor(self: *const Self, al: Allocator) !void {
            const m = self.m;
            const n = self.n;
            const mn = @min(m, n);
            const F = Float(T);
            const eps: f128 = std.math.floatEps(F);
            for (self.ipiv[0..mn], 0..) |p, k| if (p < k + 1 or p > m) {
                check.say("getrf {d}x{d}: ipiv[{d}] = {d}\n", .{ m, n, k, p });
                return error.TestUnexpectedResult;
            };
            const bound: f128 = (if (comptime isComplex(T)) 2 else 1) * (1 + 8 * eps);
            for (0..mn) |k| for (k + 1..m) |i| {
                const x = self.l(i, k);
                const s2 = if (comptime isComplex(T)) x.re * x.re + x.im * x.im else x * x;
                if (!(s2 <= bound)) {
                    check.say("getrf {d}x{d}: |L({d},{d})|^2 = {e}\n", .{ m, n, i, k, @as(f64, @floatCast(s2)) });
                    return error.TestUnexpectedResult;
                }
            };
            const pa = try self.permuted(al);
            defer al.free(pa);
            const ulps: f128 = @floatFromInt(2 * (mn + 2));
            for (0..n) |j| for (0..m) |i| {
                var s: W(T) = check.wide(T, scalar(T, 0));
                var mag: f128 = 0;
                for (0..@min(i, j, mn - 1) + 1) |k| {
                    const x = self.l(i, k);
                    const y = self.u(k, j);
                    s = check.addW(s, check.mulW(x, y));
                    mag += check.abs1(x) * check.abs1(y);
                }
                const d = check.abs1(check.subW(s, pa[i + j * m]));
                if (!(d <= ulps * eps * mag + std.math.floatMin(F))) {
                    check.say("getrf {d}x{d}: (L U)({d},{d}) off by {e}, allowed {e}\n", .{ m, n, i, j, @as(f64, @floatCast(d)), @as(f64, @floatCast(ulps * eps * mag)) });
                    return error.TestUnexpectedResult;
                }
            };
            for (0..n) |j| for (m..self.lda) |i| if (!isNanT(T, self.a[i + j * self.lda])) {
                check.say("getrf {d}x{d}: wrote padding ({d},{d})\n", .{ m, n, i, j });
                return error.TestUnexpectedResult;
            };
        }
    };
}

/// `vt` with every precision's `q` at `2 * unroll_n + 1`, so the blocks are
/// capped at an odd width just past the base case's; static, as a vtable is
/// large for a rig's stack.
var shrunk: Vtable = undefined;

fn shrink(vt: *const Vtable) *const Vtable {
    shrunk = vt.*;
    inline for (.{ "s", "d", "q", "c", "z", "x" }) |p| {
        const un = @max(1, @field(shrunk, p ++ "gemm_unroll_n"));
        @field(shrunk, p ++ "gemm_q") = 2 * un + 1;
    }
    return &shrunk;
}

/// `m x n` shapes: square, tall, wide, a row and a column, across both
/// blockings' block edges.
const shapes = [_][2]usize{
    .{ 1, 1 },   .{ 2, 2 },   .{ 3, 3 },   .{ 7, 7 },   .{ 1, 9 },   .{ 9, 1 },   .{ 17, 17 },   .{ 33, 33 },
    .{ 40, 25 }, .{ 25, 40 }, .{ 50, 50 }, .{ 70, 70 }, .{ 90, 61 }, .{ 61, 90 }, .{ 130, 130 },
};

fn factorCases(vt: *const Vtable) anyerror!void {
    const al = testing.allocator;
    inline for (test_types) |T| for ([_]*const Vtable{ shrink(vt), vt }, 0..) |v, vi| {
        const saved = blas.testUse(v);
        defer blas.testRestore(saved);
        for (shapes) |s| {
            if (check.rig() and s[0] * s[1] > 70 * 70) continue;
            var p = try Problem(T).init(al, s[0], s[1], s[0] * 31 + s[1]);
            defer p.deinit(al);
            var poison: guard.Poison = .{ .child = al };
            var pivot: usize = 99;
            p.factorIt(poison.allocator(), &pivot) catch |e| {
                check.say("getrf {s} {d}x{d}\n", .{ @typeName(T), s[0], s[1] });
                return e;
            };
            try testing.expectEqual(@as(usize, 99), pivot);
            p.expectFactor(al) catch |e| {
                check.say("{s}, blocking {s}\n", .{ @typeName(T), if (vi == 0) "shrunk" else "own" });
                return e;
            };
        }
    };
}

test "getrf: P L U = A, square, tall and wide, on every linked core" {
    try check.eachCore(factorCases);
}

/// `getrf` on `p` must fail at pivot `want`, its factorization complete.
fn expectSingular(comptime T: type, al: Allocator, p: *Problem(T), want: usize, what: []const u8) !void {
    var pivot: usize = 0;
    var poison: guard.Poison = .{ .child = al };
    if (p.factorIt(poison.allocator(), &pivot)) |_| {
        check.say("{s} {s} {d}x{d}: no error, want pivot {d}\n", .{ what, @typeName(T), p.m, p.n, want });
        return error.TestExpectedError;
    } else |e| if (e != error.Singular or pivot != want) {
        check.say("{s} {s} {d}x{d}: {s} at pivot {d}, want {d}\n", .{ what, @typeName(T), p.m, p.n, @errorName(e), pivot, want });
        return error.TestUnexpectedResult;
    }
    p.expectFactor(al) catch |e| {
        check.say("{s} {s} {d}x{d}, singular at {d}\n", .{ what, @typeName(T), p.m, p.n, want });
        return e;
    };
}

/// Columns made zero: the first and last, and either side of block edges on
/// both blockings.
const zero_at = [_]usize{ 0, 1, 2, 5, 12, 13, 17, 25, 32, 33, 69 };

fn singularCases(vt: *const Vtable) anyerror!void {
    const al = testing.allocator;
    inline for (test_types) |T| for ([_]*const Vtable{ shrink(vt), vt }, 0..) |v, vi| {
        const saved = blas.testUse(v);
        defer blas.testRestore(saved);
        const n: usize = if (vi == 0) 40 else 70;
        for ([_]usize{ n, n + 9 }) |m| {
            var p = try Problem(T).init(al, m, n, 5 + m);
            defer p.deinit(al);
            const aw = try al.dupe(W(T), p.aw);
            defer al.free(aw);
            for (zero_at) |k| {
                if (k >= n) continue;
                // a zero column `k`, then a second one after it
                @memcpy(p.aw, aw);
                for (0..m) |i| p.aw[i + k * m] = check.wide(T, scalar(T, 0));
                if (k + 3 < n) for (0..m) |i| {
                    p.aw[i + (k + 3) * m] = check.wide(T, scalar(T, 0));
                };
                p.store();
                try expectSingular(T, al, &p, k + 1, "zero column");
            }
            // rank one, `x yᵀ` with `x` of powers of two and `y(0) = 1`: every
            // step is exact, so `U(2, 2)` is zero
            for (0..m) |i| for (0..n) |j| {
                const xi = std.math.ldexp(@as(f64, if (i % 3 == 1) -1 else 1), -@as(i32, @intCast((i * 7) % 5)));
                const yj: f64 = if (j == 0) 1 else @floatFromInt(@as(i32, @intCast(j % 7)) - 3);
                p.aw[i + j * m] = check.wide(T, fromF64(T, xi * yj));
            };
            p.store();
            try expectSingular(T, al, &p, 2, "rank one");
        }
    };
}

test "getrf: the first zero pivot, the factorization completed, on every linked core" {
    try check.eachCore(singularCases);
}

/// A NaN anywhere is a pivot's or reaches one, and is never zero: no error,
/// even where it reaches a zero column, as the reference's `INFO`.
fn nanCases(vt: *const Vtable) anyerror!void {
    const al = testing.allocator;
    inline for (test_types) |T| for ([_]*const Vtable{ shrink(vt), vt }) |v| {
        const saved = blas.testUse(v);
        defer blas.testRestore(saved);
        const n = 40;
        var p = try Problem(T).init(al, n, n, 11);
        defer p.deinit(al);
        for ([_][2]usize{ .{ 0, 0 }, .{ 7, 3 }, .{ 39, 20 }, .{ 0, 39 } }) |at| {
            p.store();
            p.a[at[0] + at[1] * p.lda] = fromF64(T, nan);
            // a zero column past the NaN's
            const z = @min(at[1] + 2, n - 1);
            if (z != at[1]) for (0..n) |i| {
                p.a[i + z * p.lda] = scalar(T, 0);
            };
            var pivot: usize = 99;
            var poison: guard.Poison = .{ .child = al };
            p.factorIt(poison.allocator(), &pivot) catch |e| {
                check.say("getrf {s}: NaN at ({d},{d}): {s} at {d}\n", .{ @typeName(T), at[0], at[1], @errorName(e), pivot });
                return e;
            };
            // the last pivot is NaN
            try testing.expect(isNanT(T, p.a[(n - 1) * (p.lda + 1)]));
        }
    };
}

test "getrf: NaN is never a zero pivot, on every linked core" {
    try check.eachCore(nanCases);
}

/// The `m x 1` column `col` factored: its pivot (1-based) and `L`.
fn column(comptime T: type, comptime m: usize, col: [m]T) !struct { usize, [m]T } {
    var a = col;
    var ipiv: [1]isize = .{-1};
    try getrf(testing.failing_allocator, check.io, T, m, 1, &a, m, &ipiv, null);
    return .{ @intCast(ipiv[0]), a };
}

test "getrf: the first largest |re| + |im| pivots, the first NaN before it" {
    inline for (test_types) |T| {
        const c = struct {
            fn f(re: f64, im: f64) T {
                return if (comptime isComplex(T)) .init(@floatCast(re), @floatCast(im)) else @floatCast(re);
            }
        }.f;
        // ties keep the first; no swap for the diagonal's
        try testing.expectEqual(@as(usize, 1), (try column(T, 3, .{ c(1, 0), c(-1, 0), c(1, 0) }))[0]);
        try testing.expectEqual(@as(usize, 2), (try column(T, 4, .{ c(-1, 0), c(2, 0), c(-2, 0), c(2, 0) }))[0]);
        if (comptime isComplex(T)) {
            // 7, 5, 7, 7: the first; and 3 - 4i over 5, the same modulus
            try testing.expectEqual(@as(usize, 1), (try column(T, 4, .{ c(3, 4), c(5, 0), c(-4, 3), c(0, -7) }))[0]);
            try testing.expectEqual(@as(usize, 3), (try column(T, 3, .{ c(1, 0), c(5, 0), c(3, -4) }))[0]);
        }
        // the first NaN wins, and fills L
        const r = try column(T, 4, .{ c(1, 0), c(nan, 0), c(5, 0), c(nan, 0) });
        try testing.expectEqual(@as(usize, 2), r[0]);
        for (r[1][1..]) |x| try testing.expect(if (comptime isComplex(T)) std.math.isNan(x.re) else std.math.isNan(x));
        // an infinite pivot: finite elements over it are zero, infinite ones NaN
        const s = try column(T, 3, .{ c(1, 0), c(inf, 0), c(-inf, 0) });
        try testing.expectEqual(@as(usize, 2), s[0]);
        try testing.expectEqual(c(0, 0), s[1][1]);
        try testing.expect(if (comptime isComplex(T)) std.math.isNan(s[1][2].re) else std.math.isNan(s[1][2]));
    }
}

test "getrf: a denormal pivot is swapped in and divided by" {
    inline for (test_types) |T| {
        const F = Float(T);
        const tm = std.math.floatTrueMin(F);
        // [[4tm, 1], [8tm, 3]]: rows swapped, L(2, 1) = 1/2, U(2, 2) = -1/2
        var a: [4]T = .{ scalar(T, 4 * tm), scalar(T, 8 * tm), scalar(T, 1), scalar(T, 3) };
        var ipiv: [2]isize = @splat(-1);
        try getrf(testing.failing_allocator, check.io, T, 2, 2, &a, 2, &ipiv, null);
        try testing.expectEqual([2]isize{ 2, 2 }, ipiv);
        try testing.expectEqual([4]T{ scalar(T, 8 * tm), scalar(T, 0.5), scalar(T, 3), scalar(T, -0.5) }, a);
        if (comptime isComplex(T)) {
            // 4tm over 8tm i: -i/2
            const r = try column(T, 2, .{ .init(4 * tm, 0), .init(0, 8 * tm) });
            try testing.expectEqual(@as(usize, 2), r[0]);
            try testing.expectEqual(T.init(0, -0.5), r[1][1]);
        }
    }
}

test "getrf: an exactly singular U; pivot is written only on failure" {
    inline for (.{ f64, types.Complex64 }) |T| {
        const al = testing.failing_allocator;
        var pivot: usize = 99;
        var ipiv: [3]isize = @splat(-1);
        // a zero matrix fails at 1, a zero row's `U(1, 1)` too
        var z: [9]T = @splat(scalar(T, 0));
        try testing.expectError(error.Singular, getrf(al, check.io, T, 3, 3, &z, 3, &ipiv, &pivot));
        try testing.expectEqual(@as(usize, 1), pivot);
        try testing.expectEqual([3]isize{ 1, 2, 3 }, ipiv);
        var row: [3]T = .{ scalar(T, 0), scalar(T, 2), scalar(T, 3) };
        try testing.expectError(error.Singular, getrf(al, check.io, T, 1, 3, &row, 1, &ipiv, &pivot));
        try testing.expectEqual(@as(usize, 1), pivot);
        // [[1, 2], [2, 4]]: rows swapped, then U(2, 2) = 2 - 1/2 4 = 0
        var a: [4]T = .{ scalar(T, 1), scalar(T, 2), scalar(T, 2), scalar(T, 4) };
        try testing.expectError(error.Singular, getrf(al, check.io, T, 2, 2, &a, 2, &ipiv, &pivot));
        try testing.expectEqual(@as(usize, 2), pivot);
        try testing.expectEqual([4]T{ scalar(T, 2), scalar(T, 0.5), scalar(T, 4), scalar(T, 0) }, a);
        a = .{ scalar(T, 1), scalar(T, 2), scalar(T, 2), scalar(T, 5) };
        pivot = 99;
        try getrf(al, check.io, T, 2, 2, &a, 2, &ipiv, &pivot);
        try testing.expectEqual(@as(usize, 99), pivot);
        try testing.expectEqual([4]T{ scalar(T, 2), scalar(T, 0.5), scalar(T, 5), scalar(T, -0.5) }, a);
    }
}

test "getrf: reports the first bad argument and touches nothing when empty" {
    const fa = testing.failing_allocator;
    var a: [4]f64 = @splat(nan);
    var ipiv: [2]isize = @splat(-1);
    try testing.expectError(error.InvalidM, getrf(fa, check.io, f64, -1, -1, &a, 0, &ipiv, null));
    try testing.expectError(error.InvalidN, getrf(fa, check.io, f64, 2, -1, &a, 0, &ipiv, null));
    try testing.expectError(error.InvalidLda, getrf(fa, check.io, f64, 2, 2, &a, 1, &ipiv, null));
    try testing.expectError(error.InvalidLda, getrf(fa, check.io, f64, 0, 2, &a, 0, &ipiv, null));
    try testing.expectError(error.InvalidLda, getrf(fa, check.io, f64, 2, 2, &a, std.math.maxInt(isize) / 2, &ipiv, null));
    try getrf(fa, check.io, f64, 0, 2, &a, 1, &ipiv, null);
    try getrf(fa, check.io, f64, 2, 0, &a, 2, &ipiv, null);
    for (a) |v| try testing.expect(std.math.isNan(v));
    try testing.expectEqual([2]isize{ -1, -1 }, ipiv);
}

/// `getrs` on `p`'s factor, `nrhs` right-hand sides with `ldb = n + 1` (the
/// padding row NaN), then the residual in `W(T)` within `8 (n + 2)` ulps of
/// `L` and `U` times `X` in magnitude: `P (A X - B)` against
/// `|L| |U| |X|`, or for `op` transposed `op(A) X - B` against
/// `(|L| |U|)ᵀ P |X|`. With `nan_at`, `B(i, j)` is NaN: column `j` must come
/// out all NaN, the others as without it, bit for bit.
fn solveCase(comptime T: type, al: Allocator, p: *const Problem(T), trans: Trans, nrhs: usize, seed: u64, nan_at: ?[2]usize) !void {
    const n = p.n;
    const ldb = n + 1;
    const b = try al.alloc(T, ldb * nrhs);
    defer al.free(b);
    const b0 = try al.alloc(T, n * nrhs);
    defer al.free(b0);
    check.fill(T, b0, seed);
    for (0..nrhs) |j| {
        @memcpy(b[j * ldb ..][0..n], b0[j * n ..][0..n]);
        b[j * ldb + n] = poisoned(T);
    }
    var poison: guard.Poison = .{ .child = al };
    try getrs(poison.allocator(), check.io, T, trans, @intCast(n), @intCast(nrhs), p.a.ptr, @intCast(p.lda), p.ipiv.ptr, b.ptr, @intCast(ldb));
    // `|L| |U|`
    const mag = try al.alloc(f128, n * n);
    defer al.free(mag);
    for (0..n) |j| for (0..n) |i| {
        var s: f128 = 0;
        for (0..@min(i, j) + 1) |k| s += check.abs1(p.l(i, k)) * check.abs1(p.u(k, j));
        mag[i + j * n] = s;
    };
    const t = trans.transposed();
    const cj = isComplex(T) and trans.conjA();
    const F = Float(T);
    const eps: f128 = std.math.floatEps(F);
    const ulps: f128 = @floatFromInt(8 * (n + 2));
    const res = try al.alloc(W(T), n);
    defer al.free(res);
    const xm = try al.alloc(f128, n);
    defer al.free(xm);
    for (0..nrhs) |j| {
        const x = b[j * ldb ..][0..n];
        try testing.expect(isNanT(T, b[j * ldb + n]));
        for (0..n) |i| {
            var s: W(T) = check.wide(T, b0[i + j * n]);
            for (0..n) |k| {
                const v = if (t) p.aw[k + i * n] else p.aw[i + k * n];
                s = check.subW(s, check.mulW(if (cj) check.conjW(v) else v, check.wide(T, x[k])));
            }
            res[i] = s;
            xm[i] = check.abs1(check.wide(T, x[i]));
        }
        // `P` on the residual, or on `|X|`
        for (0..n) |k| {
            const q: usize = @intCast(p.ipiv[k] - 1);
            if (t) std.mem.swap(f128, &xm[k], &xm[q]) else std.mem.swap(W(T), &res[k], &res[q]);
        }
        for (0..n) |i| {
            var bound: f128 = 0;
            for (0..n) |k| bound += (if (t) mag[k + i * n] else mag[i + k * n]) * xm[k];
            const d = check.abs1(res[i]);
            if (!(d <= ulps * eps * bound + std.math.floatMin(F))) {
                check.say("getrs {s} {s} n={d} nrhs={d}: residual ({d},{d}) {e}, allowed {e}\n", .{ @typeName(T), @tagName(trans), n, nrhs, i, j, @as(f64, @floatCast(d)), @as(f64, @floatCast(ulps * eps * bound)) });
                return error.TestUnexpectedResult;
            }
        }
    }
    const at = nan_at orelse return;
    const bn = try al.alloc(T, ldb * nrhs);
    defer al.free(bn);
    for (0..nrhs) |j| @memcpy(bn[j * ldb ..][0..n], b0[j * n ..][0..n]);
    bn[at[0] + at[1] * ldb] = poisoned(T);
    try getrs(poison.allocator(), check.io, T, trans, @intCast(n), @intCast(nrhs), p.a.ptr, @intCast(p.lda), p.ipiv.ptr, bn.ptr, @intCast(ldb));
    for (0..nrhs) |j| for (0..n) |i| {
        const got = bn[i + j * ldb];
        const ok = if (j == at[1]) isNanT(T, got) else chol.sameBits(T, got, b[i + j * ldb]);
        if (!ok) {
            check.say("getrs {s} {s} n={d}: NaN at B({d},{d}), X({d},{d}) is {any}\n", .{ @typeName(T), @tagName(trans), n, at[0], at[1], i, j, got });
            return error.TestUnexpectedResult;
        }
    };
}

const solve_sizes = [_]usize{ 1, 2, 3, 7, 17, 33, 50, 70, 130 };

/// The real or the complex test types: a test each, so sparc64-debug's images
/// fit the firmware's window.
fn solveCases(comptime complex: bool, vt: *const Vtable) anyerror!void {
    const al = testing.allocator;
    inline for (test_types) |T| if (comptime isComplex(T) == complex) for ([_]*const Vtable{ shrink(vt), vt }, 0..) |v, vi| {
        const saved = blas.testUse(v);
        defer blas.testRestore(saved);
        for (solve_sizes) |n| {
            if (check.rig() and n > 70) continue;
            var p = try Problem(T).init(al, n, n, n * 17 + vi);
            defer p.deinit(al);
            try p.factorIt(al, null);
            for ([_]Trans{ .n, .t, .r, .c }) |trans| for ([_]usize{ 1, 3 }) |nrhs| {
                const nan_at: [2]usize = .{ (n * 7 + nrhs) % n, nrhs / 2 };
                solveCase(T, al, &p, trans, nrhs, n + nrhs, nan_at) catch |e| {
                    check.say("{s}, blocking {s}\n", .{ @typeName(T), if (vi == 0) "shrunk" else "own" });
                    return e;
                };
            };
        }
    };
}

fn solveReal(vt: *const Vtable) anyerror!void {
    return solveCases(false, vt);
}

fn solveComplex(vt: *const Vtable) anyerror!void {
    return solveCases(true, vt);
}

test "getrs: op(A) X = B on getrf's factor, every op, NaN in one column of B, real, on every linked core" {
    try check.eachCore(solveReal);
}

test "getrs: op(A) X = B on getrf's factor, every op, NaN in one column of B, complex, on every linked core" {
    try check.eachCore(solveComplex);
}

test "getrs: a known 2 x 2 solve, A and its transpose" {
    inline for (.{ f64, types.Complex64 }) |T| {
        // [[1, 2], [3, 4]] x = [5, 6]: x = [-4, 9/2]; transposed, [-1, 2]
        for ([_]struct { Trans, [2]f64 }{ .{ .n, .{ -4, 4.5 } }, .{ .t, .{ -1, 2 } } }) |c| {
            var a: [4]T = .{ scalar(T, 1), scalar(T, 3), scalar(T, 2), scalar(T, 4) };
            var ipiv: [2]isize = undefined;
            try getrf(testing.failing_allocator, check.io, T, 2, 2, &a, 2, &ipiv, null);
            var b: [2]T = .{ scalar(T, 5), scalar(T, 6) };
            try getrs(testing.allocator, check.io, T, c[0], 2, 1, &a, 2, &ipiv, &b, 2);
            for (b, c[1]) |x, w| try check.expectClose(T, scalar(T, @floatCast(w)), x, 1e-14);
        }
    }
}

test "getrs: reports the first bad argument and a pivot outside A, and touches nothing when empty" {
    const fa = testing.failing_allocator;
    var a: [4]f64 = @splat(nan);
    var b: [4]f64 = @splat(nan);
    const ipiv: [2]isize = .{ 1, 2 };
    try testing.expectError(error.InvalidTrans, getrs(fa, check.io, f64, .o, -1, -1, &a, 0, &ipiv, &b, 0));
    try testing.expectError(error.InvalidN, getrs(fa, check.io, f64, .n, -1, -1, &a, 0, &ipiv, &b, 0));
    try testing.expectError(error.InvalidNrhs, getrs(fa, check.io, f64, .t, 0, -1, &a, 0, &ipiv, &b, 0));
    try testing.expectError(error.InvalidLda, getrs(fa, check.io, f64, .n, 2, 1, &a, 1, &ipiv, &b, 0));
    try testing.expectError(error.InvalidLdb, getrs(fa, check.io, f64, .n, 2, 1, &a, 2, &ipiv, &b, 1));
    try testing.expectError(error.InvalidLdb, getrs(fa, check.io, f64, .n, 0, 1, &a, 1, &ipiv, &b, 0));
    try testing.expectError(error.InvalidLda, getrs(fa, check.io, f64, .n, 2, 2, &a, std.math.maxInt(isize) / 2, &ipiv, &b, 2));
    try testing.expectError(error.InvalidLdb, getrs(fa, check.io, f64, .t, 2, 2, &a, 2, &ipiv, &b, std.math.maxInt(isize) / 2));
    for ([_][2]isize{ .{ 0, 2 }, .{ 1, 3 }, .{ -1, 1 } }) |bad| {
        try testing.expectError(error.InvalidIpiv, getrs(fa, check.io, f64, .t, 2, 2, &a, 2, &bad, &b, 2));
    }
    try getrs(fa, check.io, f64, .n, 0, 2, &a, 1, &ipiv, &b, 1);
    try getrs(fa, check.io, f64, .n, 2, 0, &a, 2, &ipiv, &b, 2);
    for (b) |x| try testing.expect(std.math.isNan(x));
}

/// Inputs of a threaded-bits case: `A`, copied into each run's `out`, then
/// `ipiv` and the first zero pivot after it.
fn Threaded(comptime T: type) type {
    return struct {
        p: *const Problem(T),

        fn len(c: @This()) usize {
            return c.p.a.len + @min(c.p.m, c.p.n) + 1;
        }

        /// `factor` at one worker, else `factorAhead`.
        fn run(c: @This(), io: Io, want: u32, out: []T) anyerror!void {
            const p = c.p;
            const al = testing.allocator;
            const ipiv = try al.alloc(isize, @min(p.m, p.n));
            defer al.free(ipiv);
            @memcpy(out[0..p.a.len], p.a);
            const vt = chol.table(T);
            const blk = driver.Blocking.of(T, vt);
            const k = if (want == 1)
                try factor(al, io, T, vt, blk.un, blk.q, p.m, p.n, out.ptr, p.lda, ipiv.ptr)
            else
                try factorAhead(al, io, T, vt, blk, p.m, p.n, out.ptr, p.lda, ipiv.ptr, want);
            for (ipiv, out[p.a.len..][0..ipiv.len]) |v, *o| o.* = fromF64(T, @floatFromInt(v));
            out[out.len - 1] = fromF64(T, @floatFromInt(k));
        }

        /// The public `getrs` and `gesv` (their solutions after their
        /// factors), threading by `blas.config`'s workers.
        fn solve(c: @This(), io: Io, want: u32, out: []T) anyerror!void {
            const saved = blas.config.cpu_count;
            defer blas.config.cpu_count = saved;
            blas.config.cpu_count = want;
            const p = c.p;
            const n = p.n;
            const al = testing.allocator;
            const ipiv = try al.alloc(isize, n);
            defer al.free(ipiv);
            const b = out[p.a.len..];
            const nrhs = b.len / (2 * n);
            check.fill(T, b, 4);
            @memcpy(out[0..p.a.len], p.a);
            try gesv(al, io, T, @intCast(n), @intCast(nrhs), out.ptr, @intCast(p.lda), ipiv.ptr, b.ptr, @intCast(n), null);
            try getrs(al, io, T, .c, @intCast(n), @intCast(nrhs), out.ptr, @intCast(p.lda), ipiv.ptr, b.ptr + n * nrhs, @intCast(n));
        }
    };
}

fn aheadCases(vt: *const Vtable) anyerror!void {
    _ = blas.coreName();
    const al = testing.allocator;
    inline for (.{ f32, f64, types.Complex32, types.Complex64 }) |T| for ([_]*const Vtable{ shrink(vt), vt }) |v| {
        const saved = blas.testUse(v);
        defer blas.testRestore(saved);
        for ([_][2]usize{ .{ 130, 130 }, .{ 150, 90 }, .{ 90, 150 }, .{ 41, 41 } }) |s| {
            if (check.rig() and s[0] * s[1] > 70 * 70) continue;
            var p = try Problem(T).init(al, s[0], s[1], s[0] + 7 * s[1]);
            defer p.deinit(al);
            try chol.expectSameBits(T, Threaded(T){ .p = &p }, Threaded(T).run, (Threaded(T){ .p = &p }).len(), "getrf");
            // two zero columns: the first one's pivot, the factorization completed
            for (0..s[0]) |i| for ([_]usize{ 17, 30 }) |j| {
                p.a[i + j * p.lda] = scalar(T, 0);
            };
            try chol.expectSameBits(T, Threaded(T){ .p = &p }, Threaded(T).run, (Threaded(T){ .p = &p }).len(), "getrf, singular");
        }
    };
}

test "getrf's lookahead gives the serial bits at any worker count, singular too, on every linked core" {
    try check.eachCore(aheadCases);
}

test "gesv and getrs give one worker's bits at any worker count" {
    _ = blas.coreName();
    const al = testing.allocator;
    inline for (.{ f32, f64, types.Complex32, types.Complex64 }) |T| {
        // past every type's `gesv` threshold
        const n: usize = if (check.rig()) 70 else 210;
        var p = try Problem(T).init(al, n, n, 77);
        defer p.deinit(al);
        try chol.expectSameBits(T, Threaded(T){ .p = &p }, Threaded(T).solve, p.a.len + 2 * n * 5, "gesv, getrs");
    }
}

test "gesv: getrf then getrs, bit for bit; a singular U before B is touched; the first bad argument" {
    const al = testing.allocator;
    inline for (test_types) |T| for ([_]usize{ 1, 7, 33, 70, 130 }) |n| {
        if (check.rig() and n > 70) continue;
        var p = try Problem(T).init(al, n, n, n + 3);
        defer p.deinit(al);
        const nrhs = 3;
        const b = try al.alloc(T, n * nrhs);
        defer al.free(b);
        check.fill(T, b, n);
        const b2 = try al.dupe(T, b);
        defer al.free(b2);
        const a2 = try al.dupe(T, p.a);
        defer al.free(a2);
        const ipiv2 = try al.alloc(isize, n);
        defer al.free(ipiv2);
        try p.factorIt(al, null);
        try getrs(al, check.io, T, .n, @intCast(n), nrhs, p.a.ptr, @intCast(p.lda), p.ipiv.ptr, b.ptr, @intCast(n));
        var poison: guard.Poison = .{ .child = al };
        try gesv(poison.allocator(), check.io, T, @intCast(n), nrhs, a2.ptr, @intCast(p.lda), ipiv2.ptr, b2.ptr, @intCast(n), null);
        try testing.expectEqualSlices(isize, p.ipiv, ipiv2);
        for (p.a, a2) |x, y| try testing.expect(chol.sameBits(T, x, y));
        for (b, b2) |x, y| try testing.expect(chol.sameBits(T, x, y));
    };
    inline for (.{ f64, types.Complex64 }) |T| {
        // [[1, 2], [2, 4]]: rows swapped, then U(2, 2) = 0
        var a: [4]T = .{ scalar(T, 1), scalar(T, 2), scalar(T, 2), scalar(T, 4) };
        var b: [2]T = .{ scalar(T, 5), scalar(T, 6) };
        var ipiv: [2]isize = @splat(-1);
        var pivot: usize = 99;
        try testing.expectError(error.Singular, gesv(al, check.io, T, 2, 1, &a, 2, &ipiv, &b, 2, &pivot));
        try testing.expectEqual(@as(usize, 2), pivot);
        try testing.expectEqual([4]T{ scalar(T, 2), scalar(T, 0.5), scalar(T, 4), scalar(T, 0) }, a);
        try testing.expectEqual([2]T{ scalar(T, 5), scalar(T, 6) }, b);
        // no right-hand side: factored all the same
        a = .{ scalar(T, 1), scalar(T, 3), scalar(T, 2), scalar(T, 4) };
        try gesv(al, check.io, T, 2, 0, &a, 2, &ipiv, &b, 2, null);
        try testing.expectEqual([2]isize{ 2, 2 }, ipiv);
    }
    const fa = testing.failing_allocator;
    var a: [4]f64 = @splat(nan);
    var b: [4]f64 = @splat(nan);
    var ipiv: [2]isize = @splat(-1);
    try testing.expectError(error.InvalidN, gesv(fa, check.io, f64, -1, -1, &a, 0, &ipiv, &b, 0, null));
    try testing.expectError(error.InvalidNrhs, gesv(fa, check.io, f64, 2, -1, &a, 0, &ipiv, &b, 0, null));
    try testing.expectError(error.InvalidLda, gesv(fa, check.io, f64, 2, 1, &a, 1, &ipiv, &b, 0, null));
    try testing.expectError(error.InvalidLdb, gesv(fa, check.io, f64, 2, 1, &a, 2, &ipiv, &b, 1, null));
    try testing.expectError(error.InvalidLdb, gesv(fa, check.io, f64, 2, 2, &a, 2, &ipiv, &b, std.math.maxInt(isize) / 2, null));
    try gesv(fa, check.io, f64, 0, 2, &a, 1, &ipiv, &b, 1, null);
    for (a) |v| try testing.expect(std.math.isNan(v));
    try testing.expectEqual([2]isize{ -1, -1 }, ipiv);
}
