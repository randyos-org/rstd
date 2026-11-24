//! The level-3 drivers: OpenBLAS's `driver/level3/level3.c`, Goto's loops
//! (ACM TOMS 34(3), 2008, <https://doi.org/10.1145/1356052.1356053>) over a
//! vtable's packed kernels, for every routine of two packed operands (`gemm/notes.md`).

const std = @import("std");
const types = @import("../types.zig");
const Vtable = @import("../Vtable.zig");
const portable = @import("../portable.zig");
const util = @import("../util.zig");
const gemm_kernels = @import("gemm/kernel.zig");
const threaded = @import("threaded.zig");
const loops = @import("loops.zig");
const drout = @import("driver/routines.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const isComplex = types.isComplex;
const Op = portable.Op;
pub const Conj = gemm_kernels.Conj;
pub const gemm = drout.gemm;
pub const gemmBlocked = drout.gemmBlocked;
pub const gemmt = drout.gemmt;
const scaleBlock = drout.scaleBlock;
pub const symm = drout.symm;
pub const syrk = drout.syrk;
pub const unlessOne = drout.unlessOne;

const driver = @This();

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
        return fromSlots(vt, comptime .of(T, false));
    }

    /// From `vt`'s `p`, `q`, `r`, `unroll_m` and `unroll_n` of `g`; a
    /// size no rule set, or a rule left negative (`params.derivedR` for a
    /// wide element), takes its default.
    pub fn fromSlots(vt: *const Vtable, comptime g: portable.Gemm) Blocking {
        const um: usize = @intCast(@max(1, portable.paramPtr(vt, g, .unroll_m).*));
        const un: usize = @intCast(@max(1, portable.paramPtr(vt, g, .unroll_n).*));
        const p: usize = @intCast(positive(portable.paramPtr(vt, g, .p).*, default_p));
        const q: usize = @intCast(positive(portable.paramPtr(vt, g, .q).*, default_q));
        const r: usize = @intCast(positive(portable.paramPtr(vt, g, .r).*, default_r));
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

pub const raw = util.raw;

pub const isZero = util.isZero;

pub const isOne = util.isOne;

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

/// Which operand a copy packs, as the reference names its copies: `i` the
/// inner (`A`'s panel), `o` the outer (`B`'s).
pub const Side = enum { i, o };

/// `<symm|hemm>_<io><u|l>tcopy`.
fn symCopyOp(comptime herm: bool, comptime io: Side, comptime upper: bool) Op {
    return switch (io) {
        .i => if (herm) (if (upper) .hemm_iutcopy else .hemm_iltcopy) else (if (upper) .symm_iutcopy else .symm_iltcopy),
        .o => if (herm) (if (upper) .hemm_outcopy else .hemm_oltcopy) else (if (upper) .symm_outcopy else .symm_oltcopy),
    };
}

/// A symmetric (or `herm`, Hermitian) copy: `io` the side, `upper` the triangle.
fn symCopy(comptime T: type, vt: *const Vtable, comptime io: Side, herm: bool, upper: bool, m: usize, n: usize, a: [*]T, lda: usize, x: usize, y: usize, dst: [*]T) void {
    if (comptime isComplex(T)) {
        if (herm) return if (upper)
            placed(T, vt, comptime symCopyOp(true, io, true), m, n, a, lda, x, y, dst)
        else
            placed(T, vt, comptime symCopyOp(true, io, false), m, n, a, lda, x, y, dst);
    }
    if (upper) placed(T, vt, comptime symCopyOp(false, io, true), m, n, a, lda, x, y, dst) else placed(T, vt, comptime symCopyOp(false, io, false), m, n, a, lda, x, y, dst);
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
        .upper, .lower => symCopy(T, vt, .i, d.herm, d.rows == .upper, min_l, min_i, x.p, x.ld, is, ls, dst),
    }
}

/// Packs columns `jjs..` of `op(B)` over depth `ls..` (`OCOPY_OPERATION`).
pub fn packCols(comptime T: type, vt: *const Vtable, d: Desc, y: Matrix(T), min_l: usize, min_jj: usize, ls: usize, jjs: usize, dst: [*]T) void {
    switch (d.cols) {
        .n => colsCopy(T, vt, false, min_l, min_jj, y.p + ls + jjs * y.ld, y.ld, dst),
        .t => colsCopy(T, vt, true, min_l, min_jj, y.p + jjs + ls * y.ld, y.ld, dst),
        .upper, .lower => symCopy(T, vt, .o, d.herm, d.cols == .upper, min_l, min_jj, y.p, y.ld, jjs, ls, dst),
    }
}

/// `C := beta*C + alpha * op(A) * op(B)` over `C`'s `m x n` (`region` of
/// it), depth `k`, in the js / ls / is loops of `level3.c`; a null `beta`
/// is 1.
pub fn drive(comptime T: type, vt: *const Vtable, d: Desc, blk: Blocking, m: usize, n: usize, k: usize, alpha: T, x: Matrix(T), y: Matrix(T), beta: ?T, c: [*]T, ldc: usize, work: Work(T)) void {
    const Ctx = struct {
        vt: *const Vtable,
        d: Desc,
        blk: Blocking,
        m: usize,
        alpha: T,
        x: Matrix(T),
        y: Matrix(T),
        beta: ?T,
        c: [*]T,
        ldc: usize,
        work: Work(T),

        pub fn depth(s: *const @This(), rest: usize) usize {
            return Blocking.split(rest, s.blk.q, s.blk.um);
        }
        pub fn chunk(s: *const @This(), rest: usize) usize {
            return s.blk.chunk(rest);
        }
        pub fn slice(s: *const @This(), js: usize, min_j: usize, ls: usize, min_l: usize) void {
            // the rows the column block meets
            const r_lo = if (s.d.region == .lower) js else 0;
            const r_hi = if (s.d.region == .upper) @min(s.m, js + min_j) else s.m;
            loops.slice(s, s.blk, r_lo, r_hi, js, min_j, ls, min_l, if (ls == 0) s.beta else null, s.work.sa, s.work.sb);
        }
        pub fn packRows(s: *const @This(), min_l: usize, min_i: usize, ls: usize, is: usize, dst: [*]T) void {
            driver.packRows(T, s.vt, s.d, s.x, min_l, min_i, ls, is, dst);
        }
        pub fn packCols(s: *const @This(), min_l: usize, min_jj: usize, ls: usize, jjs: usize, dst: [*]T) void {
            driver.packCols(T, s.vt, s.d, s.y, min_l, min_jj, ls, jjs, dst);
        }
        pub fn block(s: *const @This(), ib: usize, mi: usize, jb: usize, nj: usize, kl: usize, b: ?T, pa: [*]T, pb: [*]T) void {
            driver.block(T, s.vt, s.d, s.blk, ib, mi, jb, nj, kl, s.alpha, b, pa, pb, s.c, s.ldc, s.work.tmp);
        }
    };
    const ctx: Ctx = .{ .vt = vt, .d = d, .blk = blk, .m = m, .alpha = alpha, .x = x, .y = y, .beta = beta, .c = c, .ldc = ldc, .work = work };
    loops.columns(blk, 0, n, k, &ctx);
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
pub fn run(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, d: Desc, m: usize, n: usize, k: usize, alpha: T, x: Matrix(T), y: Matrix(T), twice: ?T, beta: ?T, c: [*]T, ldc: usize, want: u32) Error!void {
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

// ---- trmm and trsm: `triangular.zig` -------------------------------------------

const tri = @import("triangular.zig");

pub const Tri = tri.Tri;
pub const trmm = tri.trmm;
pub const trsm = tri.trsm;
pub const triLens = tri.triLens;
pub const triDrive = tri.triDrive;
