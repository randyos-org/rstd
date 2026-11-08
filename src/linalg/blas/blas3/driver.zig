//! The level-3 drivers: OpenBLAS's `driver/level3/level3.c`, Goto's loops
//! (ACM TOMS 34(3), 2008, <https://doi.org/10.1145/1356052.1356053>) over a
//! vtable's packed kernels, for every routine of two packed operands (`gemm/notes.md`).

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("../types.zig");
const Vtable = @import("../Vtable.zig");
const portable = @import("../portable.zig");
const util = @import("../util.zig");
const gemm_kernels = @import("gemm/kernel.zig");
const small = @import("gemm/small.zig");
const threaded = @import("threaded.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const BlasLong = types.BlasLong;
const Float = types.Float;
const isComplex = types.isComplex;
const Op = portable.Op;
pub const Conj = gemm_kernels.Conj;

/// `P`, `Q`, `R` and the unrolls the loops round to, from the vtable (a
/// vtable no core filled, or a type it has no slots for, gets the generic
/// target's sizes and the portable kernels' unrolls).
pub const Blocking = struct {
    p: usize,
    q: usize,
    r: usize,
    um: usize,
    un: usize,

    /// The sizes a vtable without them gets: `Q` the generic table's, `P`
    /// and `R` ours.
    pub const default_p = 224;
    pub const default_q = 128;
    pub const default_r = 4096;

    pub fn of(comptime T: type, vt: *const Vtable) Blocking {
        if (comptime portable.letter(T) == null) {
            const K = gemm_kernels.Kernels(T, gemm_kernels.geometry(T));
            const um: usize = K.unroll_m;
            const un: usize = K.unroll_n;
            return .{ .p = @max(um, default_p / um * um), .q = default_q, .r = @max(un, default_r), .um = um, .un = un };
        }
        return fromSlots(vt, comptime portable.letter(T).? ++ "gemm_");
    }

    /// From `vt`'s slots `<pre>p`, `q`, `r`, `unroll_m` and `unroll_n`; a
    /// size no rule set, or a rule left negative (`params.derivedR` for a
    /// wide element), takes its default.
    pub fn fromSlots(vt: *const Vtable, comptime pre: []const u8) Blocking {
        const um: usize = @intCast(@max(1, @field(vt, pre ++ "unroll_m")));
        const un: usize = @intCast(@max(1, @field(vt, pre ++ "unroll_n")));
        const p: usize = @intCast(positive(@field(vt, pre ++ "p"), default_p));
        const q: usize = @intCast(positive(@field(vt, pre ++ "q"), default_q));
        const r: usize = @intCast(positive(@field(vt, pre ++ "r"), default_r));
        return .{
            // `p` a whole number of row panels
            .p = @max(um, p / um * um),
            .q = q,
            .r = @max(un, r),
            .um = um,
            .un = un,
        };
    }

    fn positive(v: i32, fallback: i32) i32 {
        return if (v > 0) v else fallback;
    }

    /// `level3.c`'s split of `rest` rows (or depth) for a block of `max`:
    /// all of it up to `max`, half of it (rounded to `unroll`) below twice
    /// `max`, else `max`.
    pub fn split(rest: usize, max: usize, unroll: usize) usize {
        if (rest >= 2 * max) return max;
        if (rest > max) return @min(rest, (rest / 2 + unroll - 1) / unroll * unroll);
        return rest;
    }

    /// The column chunk `B` is packed and multiplied in: `level3.c`'s, whose
    /// `2 * GEMM_UNROLL_N` step is commented out upstream. Its `6 *
    /// GEMM_UNROLL_N` for the Skylake-X family serves those cores' assembly,
    /// which is not converted; the Zig kernels' `unroll_n` is already wide there.
    pub fn chunk(b: Blocking, rest: usize) usize {
        if (rest >= 3 * b.un) return 3 * b.un;
        if (rest > b.un) return b.un;
        return rest;
    }
};

/// The packing space of one call, in elements: `sa` holds a row block of
/// `op(A)` (`split` never gives more than `p` rows, while a depth slice can
/// pass `q` by its rounding), `sb` a slice of `op(B)` (one chunk when a
/// single row block takes every row, so each chunk is used once), and `tmp`
/// one tile, where a triangle of `C` is kept.
pub fn Work(comptime T: type) type {
    return struct {
        const Self = @This();
        /// Every buffer starts on a cache line; more than any kernel loads need.
        pub const alignment = 64;
        sa: [*]T,
        sb: [*]T,
        tmp: [*]T,

        pub fn lens(blk: Blocking, m: usize, n: usize, k: usize, tile: bool) [3]usize {
            const rows = @min(m, blk.p);
            const depth = @min(k, blk.q + blk.um);
            const cols = if (m <= blk.p) @min(n, 3 * blk.un) else @min(n, blk.r);
            return .{ rows * depth, depth * cols, if (tile) blk.um * blk.un else 0 };
        }

        /// Elements per step of an offset: the fewest that span a whole
        /// number of cache lines (3 for a 12-byte `f80` on x86).
        const step = std.math.lcm(alignment, @sizeOf(T)) / @sizeOf(T);

        pub fn up(n: usize) usize {
            return (n + step - 1) / step * step;
        }

        /// Elements to allocate for `lens`, each buffer on a cache line.
        pub fn total(l: [3]usize) usize {
            return up(l[0]) + up(l[1]) + l[2];
        }

        pub fn at(base: [*]T, l: [3]usize) Self {
            return .{ .sa = base, .sb = base + up(l[0]), .tmp = base + up(l[0]) + up(l[1]) };
        }
    };
}

/// Bytes of packing space on the stack before the allocator is asked.
pub const stack_bytes = util.stack_bytes;

pub const Error = Allocator.Error || Io.Cancelable;

const raw = util.raw;

const isZero = util.isZero;

const isOne = util.isOne;

// ---- kernel calls ------------------------------------------------------------

/// `C += alpha * A * B` on packed `A` (`m x k`) and `B` (`k x n`), complex
/// conjugating as `conj` says.
pub fn kernel(comptime T: type, vt: *const Vtable, conj: Conj, m: usize, n: usize, k: usize, alpha: T, sa: [*]T, sb: [*]T, c: [*]T, ldc: usize) void {
    if (comptime isComplex(T)) {
        const f = switch (conj) {
            .n => portable.kernel(T, vt, .gemm_kernel_n),
            .l => portable.kernel(T, vt, .gemm_kernel_l),
            .r => portable.kernel(T, vt, .gemm_kernel_r),
            .b => portable.kernel(T, vt, .gemm_kernel_b),
        };
        _ = f(@intCast(m), @intCast(n), @intCast(k), alpha.re, alpha.im, raw(T, sa), raw(T, sb), raw(T, c), @intCast(ldc));
    } else {
        _ = portable.kernel(T, vt, .gemm_kernel)(@intCast(m), @intCast(n), @intCast(k), alpha, sa, sb, c, @intCast(ldc));
    }
}

/// `C := beta*C` on `m x n`, never reading `C` for a zero `beta`.
pub fn scale(comptime T: type, vt: *const Vtable, m: usize, n: usize, beta: T, c: [*]T, ldc: usize) void {
    if (comptime isComplex(T)) {
        const f = raw(T, c);
        _ = portable.kernel(T, vt, .gemm_beta)(@intCast(m), @intCast(n), 0, beta.re, beta.im, f, 0, f, 0, f, @intCast(ldc));
    } else {
        _ = portable.kernel(T, vt, .gemm_beta)(@intCast(m), @intCast(n), 0, beta, c, 0, c, 0, c, @intCast(ldc));
    }
}

/// A gemm packing copy (`gemm_incopy` and the like).
inline fn copy(comptime T: type, vt: *const Vtable, comptime op: Op, m: usize, n: usize, a: [*]T, lda: usize, dst: [*]T) void {
    if (comptime isComplex(T)) {
        _ = portable.kernel(T, vt, op)(@intCast(m), @intCast(n), raw(T, a), @intCast(lda), raw(T, dst));
    } else {
        _ = portable.kernel(T, vt, op)(@intCast(m), @intCast(n), a, @intCast(lda), dst);
    }
}

/// A copy that places a block itself (`symm_iutcopy` and the like): `x` and
/// `y` are the block's first line and depth in the whole matrix at `a`.
pub inline fn placed(comptime T: type, vt: *const Vtable, comptime op: Op, m: usize, n: usize, a: [*]T, lda: usize, x: usize, y: usize, dst: [*]T) void {
    if (comptime isComplex(T)) {
        _ = portable.kernel(T, vt, op)(@intCast(m), @intCast(n), raw(T, a), @intCast(lda), @intCast(x), @intCast(y), raw(T, dst));
    } else {
        _ = portable.kernel(T, vt, op)(@intCast(m), @intCast(n), a, @intCast(lda), @intCast(x), @intCast(y), dst);
    }
}

// ---- the loops -----------------------------------------------------------------

/// How an operand is packed: as it is (`n`), transposed (`t`), or whole from
/// its stored `upper` or `lower` triangle (`symm`, or `hemm` with `herm`).
pub const Pack = enum { n, t, upper, lower };

/// Which part of `C` a routine computes; a triangle includes the diagonal.
pub const Region = enum { full, upper, lower };

/// One routine's loops. Values, not compile-time choices: the loops only
/// call kernels, and one instantiation per type keeps debug builds small.
pub const Desc = struct {
    /// How `op(A)` (the rows) and `op(B)` (the columns) are packed.
    rows: Pack,
    cols: Pack,
    /// `upper` and `lower` pack Hermitian, not symmetric.
    herm: bool = false,
    conj: Conj = .n,
    region: Region = .full,
    /// Diagonal of `C` real (`herk`, `her2k`).
    real_diag: bool = false,
};

/// An operand: a matrix and its leading dimension.
pub fn Matrix(comptime T: type) type {
    return struct { p: [*]T, ld: usize };
}

/// A symmetric (or `herm`, Hermitian) copy: `io` the side, `upper` the triangle.
fn symCopy(comptime T: type, vt: *const Vtable, comptime io: []const u8, herm: bool, upper: bool, m: usize, n: usize, a: [*]T, lda: usize, x: usize, y: usize, dst: [*]T) void {
    if (comptime isComplex(T)) {
        if (herm) return if (upper)
            placed(T, vt, @field(Op, "hemm_" ++ io ++ "utcopy"), m, n, a, lda, x, y, dst)
        else
            placed(T, vt, @field(Op, "hemm_" ++ io ++ "ltcopy"), m, n, a, lda, x, y, dst);
    }
    if (upper) placed(T, vt, @field(Op, "symm_" ++ io ++ "utcopy"), m, n, a, lda, x, y, dst) else placed(T, vt, @field(Op, "symm_" ++ io ++ "ltcopy"), m, n, a, lda, x, y, dst);
}

/// Rows of `op(A)` (`gemm_incopy` for `A^T`, else `gemm_itcopy`).
pub fn rowsCopy(comptime T: type, vt: *const Vtable, trans: bool, m: usize, n: usize, a: [*]T, lda: usize, dst: [*]T) void {
    if (trans) copy(T, vt, .gemm_incopy, m, n, a, lda, dst) else copy(T, vt, .gemm_itcopy, m, n, a, lda, dst);
}

/// Columns of `op(B)` (`gemm_otcopy` for `B^T`, else `gemm_oncopy`).
pub fn colsCopy(comptime T: type, vt: *const Vtable, trans: bool, m: usize, n: usize, a: [*]T, lda: usize, dst: [*]T) void {
    if (trans) copy(T, vt, .gemm_otcopy, m, n, a, lda, dst) else copy(T, vt, .gemm_oncopy, m, n, a, lda, dst);
}

/// Packs rows `is..` of `op(A)` over depth `ls..` (`ICOPY_OPERATION`).
pub fn packRows(comptime T: type, vt: *const Vtable, d: Desc, x: Matrix(T), min_l: usize, min_i: usize, ls: usize, is: usize, dst: [*]T) void {
    switch (d.rows) {
        .n => rowsCopy(T, vt, false, min_l, min_i, x.p + is + ls * x.ld, x.ld, dst),
        .t => rowsCopy(T, vt, true, min_l, min_i, x.p + ls + is * x.ld, x.ld, dst),
        .upper, .lower => symCopy(T, vt, "i", d.herm, d.rows == .upper, min_l, min_i, x.p, x.ld, is, ls, dst),
    }
}

/// Packs columns `jjs..` of `op(B)` over depth `ls..` (`OCOPY_OPERATION`).
pub fn packCols(comptime T: type, vt: *const Vtable, d: Desc, y: Matrix(T), min_l: usize, min_jj: usize, ls: usize, jjs: usize, dst: [*]T) void {
    switch (d.cols) {
        .n => colsCopy(T, vt, false, min_l, min_jj, y.p + ls + jjs * y.ld, y.ld, dst),
        .t => colsCopy(T, vt, true, min_l, min_jj, y.p + jjs + ls * y.ld, y.ld, dst),
        .upper, .lower => symCopy(T, vt, "o", d.herm, d.cols == .upper, min_l, min_jj, y.p, y.ld, jjs, ls, dst),
    }
}

/// `C := beta*C + alpha * op(A) * op(B)` over `C`'s `m x n` (`region` of
/// it), depth `k`, in the js / ls / is loops of `level3.c`; a null `beta`
/// is 1.
pub fn drive(comptime T: type, vt: *const Vtable, d: Desc, blk: Blocking, m: usize, n: usize, k: usize, alpha: T, x: Matrix(T), y: Matrix(T), beta: ?T, c: [*]T, ldc: usize, work: Work(T)) void {
    var js: usize = 0;
    while (js < n) : (js += blk.r) {
        const min_j = @min(blk.r, n - js);
        // the rows the column block meets
        const r_lo = if (d.region == .lower) js else 0;
        const r_hi = if (d.region == .upper) @min(m, js + min_j) else m;
        var ls: usize = 0;
        var min_l: usize = 0;
        while (ls < k) : (ls += min_l) {
            min_l = Blocking.split(k - ls, blk.q, blk.um);
            const s = if (ls == 0) beta else null;
            // the first row block packs `B` as it goes
            var min_i = Blocking.split(r_hi - r_lo, blk.p, blk.um);
            const once = min_i == r_hi - r_lo;
            packRows(T, vt, d, x, min_l, min_i, ls, r_lo, work.sa);
            var jjs = js;
            while (jjs < js + min_j) {
                const min_jj = blk.chunk(js + min_j - jjs);
                const pb = if (once) work.sb else work.sb + min_l * (jjs - js);
                packCols(T, vt, d, y, min_l, min_jj, ls, jjs, pb);
                block(T, vt, d, blk, r_lo, min_i, jjs, min_jj, min_l, alpha, s, work.sa, pb, c, ldc, work.tmp);
                jjs += min_jj;
            }
            var is = r_lo + min_i;
            while (is < r_hi) : (is += min_i) {
                min_i = Blocking.split(r_hi - is, blk.p, blk.um);
                packRows(T, vt, d, x, min_l, min_i, ls, is, work.sa);
                jjs = js;
                while (jjs < js + min_j) {
                    const min_jj = blk.chunk(js + min_j - jjs);
                    block(T, vt, d, blk, is, min_i, jjs, min_jj, min_l, alpha, s, work.sa, work.sb + min_l * (jjs - js), c, ldc, work.tmp);
                    jjs += min_jj;
                }
            }
        }
    }
}

/// Whether rows `[ib, ie)` and columns `[jb, je)` lie wholly in `region`,
/// not touching the diagonal (strictly, so every diagonal element goes
/// through `tmp`, where `real_diag` sees it), or wholly outside it.
fn inside(region: Region, ib: usize, ie: usize, jb: usize, je: usize) bool {
    return switch (region) {
        .full => true,
        .upper => ie <= jb,
        .lower => ib >= je,
    };
}
fn outside(region: Region, ib: usize, ie: usize, jb: usize, je: usize) bool {
    return switch (region) {
        .full => false,
        .upper => ib >= je,
        .lower => ie <= jb,
    };
}

/// The kernel on packed rows `[ib, ib + mi)` (`pa`) and columns `[jb, jb +
/// nj)` (`pb`), within `d.region`: blocks inside it at once, then panel by
/// panel, and the tiles across the diagonal (`syrk_kernel.c`) into `tmp`,
/// added for the elements that belong; scaled by `beta` first, unless null.
pub fn block(comptime T: type, vt: *const Vtable, d: Desc, blk: Blocking, ib: usize, mi: usize, jb: usize, nj: usize, kl: usize, alpha: T, beta: ?T, pa: [*]T, pb: [*]T, c: [*]T, ldc: usize, tmp: [*]T) void {
    const reg = d.region;
    if (beta) |s| scaleBlock(T, vt, reg, d.real_diag, ib, ib + mi, jb, jb + nj, s, c, ldc);
    if (inside(reg, ib, ib + mi, jb, jb + nj)) return kernel(T, vt, d.conj, mi, nj, kl, alpha, pa, pb, c + ib + jb * ldc, ldc);
    if (outside(reg, ib, ib + mi, jb, jb + nj)) return;
    var r: usize = 0;
    while (r < mi) {
        const w = gemm_kernels.panelWidth(mi - r, blk.um);
        const ri = ib + r;
        if (inside(reg, ri, ri + w, jb, jb + nj)) {
            kernel(T, vt, d.conj, w, nj, kl, alpha, pa + r * kl, pb, c + ri + jb * ldc, ldc);
        } else if (!outside(reg, ri, ri + w, jb, jb + nj)) {
            var q: usize = 0;
            while (q < nj) {
                const wq = gemm_kernels.panelWidth(nj - q, blk.un);
                const qj = jb + q;
                if (inside(reg, ri, ri + w, qj, qj + wq)) {
                    kernel(T, vt, d.conj, w, wq, kl, alpha, pa + r * kl, pb + q * kl, c + ri + qj * ldc, ldc);
                } else if (!outside(reg, ri, ri + w, qj, qj + wq)) {
                    scale(T, vt, w, wq, zero(T), tmp, w);
                    kernel(T, vt, d.conj, w, wq, kl, alpha, pa + r * kl, pb + q * kl, tmp, w);
                    for (0..wq) |jj| for (0..w) |ii| {
                        const gi = ri + ii;
                        const gj = qj + jj;
                        if (if (reg == .upper) gi > gj else gi < gj) continue;
                        const t = tmp[ii + jj * w];
                        const p = &c[gi + gj * ldc];
                        if (comptime isComplex(T)) {
                            p.re += t.re;
                            p.im = if (d.real_diag and gi == gj) 0 else p.im + t.im;
                        } else p.* += t;
                    };
                }
                q += wq;
            }
        }
        r += w;
    }
}

const zero = util.zero;

/// Runs `d`'s loops with packing space from `allocator`, `beta` (null for
/// 1) folded into them; `twice` runs them again with the operands swapped
/// and `alpha2` (`syr2k`). More than one worker `want`ed, on `io`'s threads.
fn run(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, d: Desc, m: usize, n: usize, k: usize, alpha: T, x: Matrix(T), y: Matrix(T), twice: ?T, beta: ?T, c: [*]T, ldc: usize, want: u32) Error!void {
    if (want > 1) return threaded.run(allocator, io, T, vt, d, m, n, k, alpha, x, y, twice, beta, c, ldc, want);
    const blk = Blocking.of(T, vt);
    const lens = Work(T).lens(blk, m, n, k, d.region != .full);
    const len = Work(T).total(lens);
    var stack: [stack_bytes]u8 align(Work(T).alignment) = undefined;
    const sc = try util.Scratch(T, Work(T).alignment, Work(T).alignment).init(allocator, &stack, len);
    defer sc.deinit(allocator);
    const work = Work(T).at(sc.ptr, lens);
    drive(T, vt, d, blk, m, n, k, alpha, x, y, beta, c, ldc, work);
    if (twice) |alpha2| drive(T, vt, d, blk, m, n, k, alpha2, y, x, null, c, ldc, work);
}

// ---- the routines ----------------------------------------------------------------

/// `C := alpha*op(A)*op(B) + beta*C` with the kernels of `vt`; `ta`, `tb`
/// say whether `op` transposes, `conj` which it conjugates (complex).
/// Arguments are checked by the caller (`interface.zig`), and `m`, `n` are
/// positive. Up to `want` workers run it (`threaded.workers`).
pub fn gemm(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, ta: bool, tb: bool, conj: Conj, m: usize, n: usize, k: usize, alpha: T, a: [*]const T, lda: usize, b: [*]const T, ldb: usize, beta: T, c: [*]T, ldc: usize, want: u32) Error!void {
    if (k == 0 or isZero(T, alpha)) {
        if (!isOne(T, beta)) scale(T, vt, m, n, beta, c, ldc);
        return;
    }
    const cj = comptime isComplex(T);
    if (small.gemm(T, vt, ta, tb, cj and conj.ofA(), cj and conj.ofB(), m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)) return;
    const x: Matrix(T) = .{ .p = @constCast(a), .ld = lda };
    const y: Matrix(T) = .{ .p = @constCast(b), .ld = ldb };
    const d: Desc = .{ .rows = if (ta) .t else .n, .cols = if (tb) .t else .n, .conj = if (comptime isComplex(T)) conj else .n };
    try run(allocator, io, T, vt, d, m, n, k, alpha, x, y, null, unlessOne(T, beta), c, ldc, want);
}

/// `C := alpha*op(A)*op(B) + beta*C` on the `upper` or lower triangle of the
/// `m x m` `C`, `op(A)` being `m x k`: `gemm`'s operands over `syrk`'s
/// triangle. Upstream's `interface/gemmt.c` takes a `gemv` per column
/// instead, conjugating `B` in place for it; this reads `B` only. `m`
/// positive.
pub fn gemmt(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, upper: bool, ta: bool, tb: bool, conj: Conj, m: usize, k: usize, alpha: T, a: [*]const T, lda: usize, b: [*]const T, ldb: usize, beta: T, c: [*]T, ldc: usize) Error!void {
    const region: Region = if (upper) .upper else .lower;
    if (k == 0 or isZero(T, alpha)) {
        if (!isOne(T, beta)) scaleBlock(T, vt, region, false, 0, m, 0, m, beta, c, ldc);
        return;
    }
    const x: Matrix(T) = .{ .p = @constCast(a), .ld = lda };
    const y: Matrix(T) = .{ .p = @constCast(b), .ld = ldb };
    const d: Desc = .{ .rows = if (ta) .t else .n, .cols = if (tb) .t else .n, .conj = if (comptime isComplex(T)) conj else .n, .region = region };
    try run(allocator, io, T, vt, d, m, m, k, alpha, x, y, null, unlessOne(T, beta), c, ldc, 1);
}

/// `C := alpha*S*B + beta*C` (`left`) or `alpha*B*S + beta*C`, `S` the
/// symmetric (`herm`: Hermitian) matrix whose `upper` or lower triangle `a`
/// holds; `symm_k.c` and `zhemm_k.c` over `level3.c`. `m`, `n` positive.
/// Up to `want` workers run it.
pub fn symm(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, comptime herm: bool, left: bool, upper: bool, m: usize, n: usize, alpha: T, a: [*]const T, lda: usize, b: [*]const T, ldb: usize, beta: T, c: [*]T, ldc: usize, want: u32) Error!void {
    if (isZero(T, alpha)) {
        if (!isOne(T, beta)) scale(T, vt, m, n, beta, c, ldc);
        return;
    }
    const s: Matrix(T) = .{ .p = @constCast(a), .ld = lda };
    const g: Matrix(T) = .{ .p = @constCast(b), .ld = ldb };
    const t: Pack = if (upper) .upper else .lower;
    // `hemm`'s right side packs `S^T`: the kernel conjugates it (`NC`)
    const d: Desc = if (left) .{ .rows = t, .cols = .n, .herm = herm } else .{ .rows = .n, .cols = t, .herm = herm, .conj = if (herm) .r else .n };
    try run(allocator, io, T, vt, d, m, n, if (left) m else n, alpha, if (left) s else g, if (left) g else s, null, unlessOne(T, beta), c, ldc, want);
}

/// `C := beta*C` on rows `[ib, ie)` and columns `[jb, je)` of `C`, the part
/// of them in `region`; with `real_diag` `beta` is real and the diagonal
/// comes out real (`zherk_k.c`'s `syrk_beta`). A zero `beta` never reads `C`.
fn scaleBlock(comptime T: type, vt: *const Vtable, region: Region, real_diag: bool, ib: usize, ie: usize, jb: usize, je: usize, beta: T, c: [*]T, ldc: usize) void {
    if (region == .full) return scale(T, vt, ie - ib, je - jb, beta, c + ib + jb * ldc, ldc);
    for (jb..je) |j| {
        const lo = if (region == .lower) @max(ib, j) else ib;
        const hi = if (region == .upper) @min(ie, j + 1) else ie;
        if (lo >= hi) continue;
        const col = c + lo + j * ldc;
        if (comptime isComplex(T)) {
            if (real_diag) {
                // a real `beta` scales each part, as the reference's real times complex
                const rows = hi - lo;
                _ = portable.kernel(Float(T), vt, .gemm_beta)(@intCast(2 * rows), 1, 0, beta.re, raw(T, col), 0, raw(T, col), 0, raw(T, col), @intCast(2 * rows));
                if (lo <= j and j < hi) c[j + j * ldc].im = 0;
                continue;
            }
        }
        scale(T, vt, hi - lo, 1, beta, col, ldc);
    }
}

/// `beta` as `run` takes it: null for 1.
pub fn unlessOne(comptime T: type, beta: T) ?T {
    return if (isOne(T, beta)) null else beta;
}

/// `C := alpha*op(A)*op(A)' + beta*C` on the `upper` or lower triangle of
/// the `n x n` `C`, `op(A)` being `n x k`; `'` is the transpose, or with
/// `herm` the conjugate transpose (`alpha`, `beta` then real, given as
/// complex with zero imaginary parts). `trans` takes `op(A) = A'`. With `b`,
/// `syr2k`: `alpha*op(A)*op(B)' + alpha2*op(B)*op(A)' + beta*C`, `alpha2`
/// being `alpha` (or its conjugate for `herm`). `level3_syrk.c` and
/// `level3_syr2k.c`, over `syrk_k.c`'s triangle. `n` positive. Up to `want`
/// workers run it.
pub fn syrk(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, comptime herm: bool, upper: bool, trans: bool, n: usize, k: usize, alpha: T, a: [*]const T, lda: usize, b: ?[*]const T, ldb: usize, beta: T, c: [*]T, ldc: usize, want: u32) Error!void {
    const region: Region = if (upper) .upper else .lower;
    if (k == 0 or isZero(T, alpha)) {
        if (!isOne(T, beta)) scaleBlock(T, vt, region, herm, 0, n, 0, n, beta, c, ldc);
        return;
    }
    const x: Matrix(T) = .{ .p = @constCast(a), .ld = lda };
    const y: Matrix(T) = if (b) |p| .{ .p = @constCast(p), .ld = ldb } else x;
    const alpha2: ?T = if (b == null) null else if (comptime herm) alpha.conjugate() else alpha;
    // `herk`'s `op(A)'` conjugates the columns (N) or `op(A)` the rows (C)
    try run(allocator, io, T, vt, .{
        .rows = if (trans) .t else .n,
        .cols = if (trans) .n else .t,
        .conj = if (herm) (if (trans) .l else .r) else .n,
        .region = region,
        .real_diag = herm,
    }, n, n, k, alpha, x, y, alpha2, unlessOne(T, beta), c, ldc, want);
}

// ---- trmm and trsm: `triangular.zig` -------------------------------------------

const tri = @import("triangular.zig");

pub const Tri = tri.Tri;
pub const trmm = tri.trmm;
pub const trsm = tri.trsm;
pub const triLens = tri.triLens;
pub const triDrive = tri.triDrive;

// ---- tests ---------------------------------------------------------------------

const check = @import("../_test/check.zig");
const check3 = @import("../_test/blas3/check.zig");

const fromReal = util.fromReal;

const Routine = enum { gemm, symm, hemm, syrk, herk, syr2k, her2k, trmm, trsm };

/// The operands of one `foldedCase`.
fn Operands(comptime T: type) type {
    return struct { a: []T, b: []T, c: []T, ld: usize };
}

/// `r` with `s` as its `beta` (`trmm`'s and `trsm`'s `alpha`), options
/// `opt`, on `o`.
fn callWith(comptime T: type, vt: *const Vtable, r: Routine, opt: [3]bool, s: T, m: usize, n: usize, k: usize, o: Operands(T)) !void {
    const al = std.testing.allocator;
    const alpha: T = if (comptime isComplex(T)) .init(-0.75, if (r == .herk) 0 else 0.25) else -0.75;
    const ld = o.ld;
    const a = o.a.ptr;
    const b = o.b.ptr;
    const c = o.c.ptr;
    switch (r) {
        .gemm => try gemm(al, check.io, T, vt, opt[0], opt[1], if (opt[2] and isComplex(T)) .b else .n, m, n, k, alpha, a, ld, b, ld, s, c, ld, 1),
        .symm => try symm(al, check.io, T, vt, false, opt[0], opt[1], m, n, alpha, a, ld, b, ld, s, c, ld, 1),
        .syrk => try syrk(al, check.io, T, vt, false, opt[0], opt[1], n, k, alpha, a, ld, null, 0, s, c, ld, 1),
        .syr2k => try syrk(al, check.io, T, vt, false, opt[0], opt[1], n, k, alpha, a, ld, b, ld, s, c, ld, 1),
        .hemm, .herk, .her2k => if (comptime isComplex(T)) switch (r) {
            .hemm => try symm(al, check.io, T, vt, true, opt[0], opt[1], m, n, alpha, a, ld, b, ld, s, c, ld, 1),
            .herk => try syrk(al, check.io, T, vt, true, opt[0], opt[1], n, k, alpha, a, ld, null, 0, s, c, ld, 1),
            else => try syrk(al, check.io, T, vt, true, opt[0], opt[1], n, k, alpha, a, ld, b, ld, s, c, ld, 1),
        } else unreachable,
        .trmm => try trmm(al, check.io, T, vt, opt[0], opt[1], opt[2], false, false, m, n, s, a, ld, c, ld, 1),
        .trsm => try trsm(al, check.io, T, vt, opt[0], opt[1], opt[2], isComplex(T), false, m, n, s, a, ld, c, ld, 1),
    }
}

/// `r` with `s` folded in against a separate pass of `s` first and then 1,
/// on a blocking shrunk so every edge comes at small sizes: the same bits.
fn foldedCase(comptime T: type, vt: *const Vtable, r: Routine, opt: [3]bool, s: T) !void {
    const al = std.testing.allocator;
    const blk = Blocking.of(T, vt);
    const m = 2 * blk.p + blk.um + 1;
    const n = 2 * blk.r + 3;
    const k = 2 * blk.q + 3;
    const ld = @max(m, n, k) + 1;
    const a = try al.alloc(T, ld * ld);
    defer al.free(a);
    const b = try al.alloc(T, ld * ld);
    defer al.free(b);
    const c = try al.alloc(T, ld * n);
    defer al.free(c);
    const c2 = try al.alloc(T, ld * n);
    defer al.free(c2);
    check.fill(T, a, 1);
    check.fill(T, b, 2);
    check.fill(T, c, 3);
    // a tame triangle for `trsm`: small off the diagonal, at least 1 on it
    const f: Float(T) = 1 / @as(Float(T), @floatFromInt(ld));
    for (0..ld) |j| for (0..ld) |i| {
        const p = &a[i + j * ld];
        const g: Float(T) = if (i == j) 1 else f;
        const h: Float(T) = if (i == j) 2 else 0;
        p.* = if (comptime isComplex(T)) .init(p.re * g + h, p.im * g) else p.* * g + h;
    };
    // a zero `s` must not read `C`
    if (isZero(T, s)) @memset(c, if (comptime isComplex(T)) .init(std.math.nan(Float(T)), 0) else std.math.nan(T));
    @memcpy(c2, c);
    const o: Operands(T) = .{ .a = a, .b = b, .c = c, .ld = ld };
    try callWith(T, vt, r, opt, s, m, n, k, o);
    const rank = r == .syrk or r == .herk or r == .syr2k or r == .her2k;
    if (rank) {
        scaleBlock(T, vt, if (opt[0]) .upper else .lower, r == .herk or r == .her2k, 0, n, 0, n, s, c2.ptr, ld);
    } else scale(T, vt, m, n, s, c2.ptr, ld);
    try callWith(T, vt, r, opt, util.one(T), m, n, k, .{ .a = a, .b = b, .c = c2, .ld = ld });
    try std.testing.expectEqualSlices(u8, std.mem.sliceAsBytes(c2), std.mem.sliceAsBytes(c));
}

test "beta folded into the first slice gives the bits of a separate pass first" {
    inline for (check.sd_types) |T| {
        var vt: Vtable = .{};
        portable.install(&vt);
        const own = Blocking.of(T, &vt);
        check3.shrinkBlocking(T, &vt, false, own.um, own.un);
        const complex = comptime isComplex(T);
        const ss = [_]T{ fromReal(T, 0), if (complex) .init(0.5, -0.25) else 0.5, fromReal(T, -2) };
        for (std.enums.values(Routine)) |r| {
            if (!complex and (r == .hemm or r == .herk or r == .her2k)) continue;
            for (0..8) |bits| {
                const opt: [3]bool = .{ bits & 1 != 0, bits & 2 != 0, bits & 4 != 0 };
                for (ss) |s0| {
                    // `trmm`'s and `trsm`'s zero `alpha` returns before the drivers
                    if ((r == .trmm or r == .trsm) and isZero(T, s0)) continue;
                    // `herk`'s and `her2k`'s `beta` is real
                    const s = if (complex and (r == .herk or r == .her2k)) T.init(s0.re, 0) else s0;
                    foldedCase(T, &vt, r, opt, s) catch |e| {
                        std.debug.print("{s} {t} opt={any} s={any}\n", .{ @typeName(T), r, opt, s });
                        return e;
                    };
                }
            }
        }
    }
}
