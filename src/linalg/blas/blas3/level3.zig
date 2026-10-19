//! The level-3 drivers: OpenBLAS's `driver/level3/level3.c`, Goto's loops
//! over the packed kernels of a vtable, single-threaded. `C` is split into
//! `R`-wide column blocks, the depth into `Q`-deep slices and the rows into
//! `P`-tall blocks; each slice of `op(B)` is packed once (`sb`), each block
//! of `op(A)` once per slice (`sa`), and the kernel runs on the pair.
//!
//! The same loops serve every routine that is a product of two packed
//! operands (`Desc`): `gemm`; `symm` and `hemm`, whose symmetric operand is
//! packed whole from its stored triangle (`symm_k.c`, `zhemm_k.c`); and the
//! rank-k updates, which keep one triangle of `C` (`level3_syrk.c`,
//! `level3_syr2k.c`).
//!
//! The packing space comes from the caller's allocator, sized to the problem
//! (`Work`), with a small stack fast path. Unlike the reference, which hands
//! the kernel all of a slice's `sb` at once for the later row blocks, the
//! kernel is called on the same `3 * unroll_n` column chunks it was packed
//! in: equal access order, and correct whatever unroll the vtable states.

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("../types.zig");
const Vtable = @import("../Vtable.zig");
const portable = @import("../portable.zig");
const gemm_kernels = @import("gemm/kernel.zig");

const Allocator = std.mem.Allocator;
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

    pub fn of(comptime T: type, vt: *const Vtable) Blocking {
        if (comptime portable.letter(T) == null) {
            const K = gemm_kernels.Kernels(T, gemm_kernels.geometry(T));
            const um: usize = K.unroll_m;
            const un: usize = K.unroll_n;
            return .{ .p = @max(um, 224 / um * um), .q = 128, .r = @max(un, 4096), .um = um, .un = un };
        }
        const pre = comptime portable.letter(T).? ++ "gemm_";
        // sizes no rule set, or a rule left negative (`params.derivedR` for a
        // wide element), take the generic target's
        const um: usize = @intCast(@max(1, @field(vt, pre ++ "unroll_m")));
        const un: usize = @intCast(@max(1, @field(vt, pre ++ "unroll_n")));
        const p: usize = @intCast(positive(@field(vt, pre ++ "p"), 224));
        const q: usize = @intCast(positive(@field(vt, pre ++ "q"), 128));
        const r: usize = @intCast(positive(@field(vt, pre ++ "r"), 4096));
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
    fn split(rest: usize, max: usize, unroll: usize) usize {
        if (rest >= 2 * max) return max;
        if (rest > max) return @min(rest, (rest / 2 + unroll - 1) / unroll * unroll);
        return rest;
    }

    /// The column chunk `B` is packed and multiplied in: `level3.c`'s, whose
    /// `2 * GEMM_UNROLL_N` step is commented out upstream. Its `6 *
    /// GEMM_UNROLL_N` for the Skylake-X family serves those cores' assembly,
    /// which is not converted; the Zig kernels' `unroll_n` is already wide there.
    fn chunk(b: Blocking, rest: usize) usize {
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

        fn up(n: usize) usize {
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
pub const stack_bytes = 2048;

pub const Error = Allocator.Error;

fn raw(comptime T: type, p: [*]T) [*]Float(T) {
    return @ptrCast(p);
}

fn isZero(comptime T: type, v: T) bool {
    return if (comptime isComplex(T)) v.re == 0 and v.im == 0 else v == 0;
}

fn isOne(comptime T: type, v: T) bool {
    return if (comptime isComplex(T)) v.re == 1 and v.im == 0 else v == 1;
}

// ---- kernel calls ------------------------------------------------------------

/// `C += alpha * A * B` on packed `A` (`m x k`) and `B` (`k x n`), complex
/// conjugating as `conj` says.
inline fn kernel(comptime T: type, vt: *const Vtable, comptime conj: Conj, m: usize, n: usize, k: usize, alpha: T, sa: [*]T, sb: [*]T, c: [*]T, ldc: usize) void {
    if (comptime isComplex(T)) {
        const op: Op = switch (conj) {
            .n => .gemm_kernel_n,
            .l => .gemm_kernel_l,
            .r => .gemm_kernel_r,
            .b => .gemm_kernel_b,
        };
        _ = portable.kernel(T, vt, op)(@intCast(m), @intCast(n), @intCast(k), alpha.re, alpha.im, raw(T, sa), raw(T, sb), raw(T, c), @intCast(ldc));
    } else {
        _ = portable.kernel(T, vt, .gemm_kernel)(@intCast(m), @intCast(n), @intCast(k), alpha, sa, sb, c, @intCast(ldc));
    }
}

/// `C := beta*C` on `m x n`, never reading `C` for a zero `beta`.
fn scale(comptime T: type, vt: *const Vtable, m: usize, n: usize, beta: T, c: [*]T, ldc: usize) void {
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
inline fn placed(comptime T: type, vt: *const Vtable, comptime op: Op, m: usize, n: usize, a: [*]T, lda: usize, x: usize, y: usize, dst: [*]T) void {
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

/// One routine's loops, all fixed at compile time.
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

/// Packs rows `is..` of `op(A)` over depth `ls..` (`ICOPY_OPERATION`).
inline fn packRows(comptime T: type, vt: *const Vtable, comptime d: Desc, x: Matrix(T), min_l: usize, min_i: usize, ls: usize, is: usize, dst: [*]T) void {
    const pre = if (d.herm) "hemm_" else "symm_";
    switch (d.rows) {
        .n => copy(T, vt, .gemm_itcopy, min_l, min_i, x.p + is + ls * x.ld, x.ld, dst),
        .t => copy(T, vt, .gemm_incopy, min_l, min_i, x.p + ls + is * x.ld, x.ld, dst),
        .upper => placed(T, vt, @field(Op, pre ++ "iutcopy"), min_l, min_i, x.p, x.ld, is, ls, dst),
        .lower => placed(T, vt, @field(Op, pre ++ "iltcopy"), min_l, min_i, x.p, x.ld, is, ls, dst),
    }
}

/// Packs columns `jjs..` of `op(B)` over depth `ls..` (`OCOPY_OPERATION`).
inline fn packCols(comptime T: type, vt: *const Vtable, comptime d: Desc, y: Matrix(T), min_l: usize, min_jj: usize, ls: usize, jjs: usize, dst: [*]T) void {
    const pre = if (d.herm) "hemm_" else "symm_";
    switch (d.cols) {
        .n => copy(T, vt, .gemm_oncopy, min_l, min_jj, y.p + ls + jjs * y.ld, y.ld, dst),
        .t => copy(T, vt, .gemm_otcopy, min_l, min_jj, y.p + jjs + ls * y.ld, y.ld, dst),
        .upper => placed(T, vt, @field(Op, pre ++ "outcopy"), min_l, min_jj, y.p, y.ld, jjs, ls, dst),
        .lower => placed(T, vt, @field(Op, pre ++ "oltcopy"), min_l, min_jj, y.p, y.ld, jjs, ls, dst),
    }
}

/// `C += alpha * op(A) * op(B)` over `C`'s `m x n` (`region` of it), depth
/// `k`, in the js / ls / is loops of `level3.c`.
fn drive(comptime T: type, vt: *const Vtable, comptime d: Desc, blk: Blocking, m: usize, n: usize, k: usize, alpha: T, x: Matrix(T), y: Matrix(T), c: [*]T, ldc: usize, work: Work(T)) void {
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
            // the first row block packs `B` as it goes
            var min_i = Blocking.split(r_hi - r_lo, blk.p, blk.um);
            const once = min_i == r_hi - r_lo;
            packRows(T, vt, d, x, min_l, min_i, ls, r_lo, work.sa);
            var jjs = js;
            while (jjs < js + min_j) {
                const min_jj = blk.chunk(js + min_j - jjs);
                const pb = if (once) work.sb else work.sb + min_l * (jjs - js);
                packCols(T, vt, d, y, min_l, min_jj, ls, jjs, pb);
                block(T, vt, d, blk, r_lo, min_i, jjs, min_jj, min_l, alpha, work.sa, pb, c, ldc, work.tmp);
                jjs += min_jj;
            }
            var is = r_lo + min_i;
            while (is < r_hi) : (is += min_i) {
                min_i = Blocking.split(r_hi - is, blk.p, blk.um);
                packRows(T, vt, d, x, min_l, min_i, ls, is, work.sa);
                jjs = js;
                while (jjs < js + min_j) {
                    const min_jj = blk.chunk(js + min_j - jjs);
                    block(T, vt, d, blk, is, min_i, jjs, min_jj, min_l, alpha, work.sa, work.sb + min_l * (jjs - js), c, ldc, work.tmp);
                    jjs += min_jj;
                }
            }
        }
    }
}

/// Whether rows `[ib, ie)` and columns `[jb, je)` lie wholly in `region`,
/// not touching the diagonal (strictly, so every diagonal element goes
/// through `tmp`, where `real_diag` sees it), or wholly outside it.
fn inside(comptime region: Region, ib: usize, ie: usize, jb: usize, je: usize) bool {
    return switch (region) {
        .full => true,
        .upper => ie <= jb,
        .lower => ib >= je,
    };
}
fn outside(comptime region: Region, ib: usize, ie: usize, jb: usize, je: usize) bool {
    return switch (region) {
        .full => false,
        .upper => ib >= je,
        .lower => ie <= jb,
    };
}

/// The kernel on packed rows `[ib, ib + mi)` (`pa`) and columns `[jb, jb +
/// nj)` (`pb`), within `d.region`: blocks inside it at once, then panel by
/// panel, and the tiles across the diagonal (`syrk_kernel.c`) into `tmp`,
/// added for the elements that belong.
fn block(comptime T: type, vt: *const Vtable, comptime d: Desc, blk: Blocking, ib: usize, mi: usize, jb: usize, nj: usize, kl: usize, alpha: T, pa: [*]T, pb: [*]T, c: [*]T, ldc: usize, tmp: [*]T) void {
    const reg = d.region;
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

fn zero(comptime T: type) T {
    return if (comptime isComplex(T)) .init(0, 0) else 0;
}

/// Runs `d`'s loops with packing space from `allocator`; `twice` runs them
/// again with the operands swapped and `alpha2` (`syr2k`).
fn run(allocator: Allocator, comptime T: type, vt: *const Vtable, comptime d: Desc, m: usize, n: usize, k: usize, alpha: T, x: Matrix(T), y: Matrix(T), twice: ?T, c: [*]T, ldc: usize) Error!void {
    const blk = Blocking.of(T, vt);
    const lens = Work(T).lens(blk, m, n, k, d.region != .full);
    const len = Work(T).total(lens);
    var stack: [stack_bytes]u8 align(Work(T).alignment) = undefined;
    var heap: ?[]align(Work(T).alignment) T = null;
    defer if (heap) |h| allocator.free(h);
    const base: [*]T = if (len * @sizeOf(T) <= stack_bytes) stk: {
        // tests poison it as `check.zig`'s allocator does the heap: a packing
        // buffer read before it is written shows as NaN
        if (sysinfo.is_test) @memset(&stack, 0xff);
        break :stk @ptrCast(&stack);
    } else hp: {
        heap = try allocator.alignedAlloc(T, .fromByteUnits(Work(T).alignment), len);
        break :hp heap.?.ptr;
    };
    const work = Work(T).at(base, lens);
    drive(T, vt, d, blk, m, n, k, alpha, x, y, c, ldc, work);
    if (twice) |alpha2| drive(T, vt, d, blk, m, n, k, alpha2, y, x, c, ldc, work);
}

// ---- the routines ----------------------------------------------------------------

/// `C := alpha*op(A)*op(B) + beta*C` with the kernels of `vt`; `ta`, `tb`
/// say whether `op` transposes, `conj` which it conjugates (complex).
/// Arguments are checked by the caller (`interface.zig`), and `m`, `n` are
/// positive.
pub fn gemm(allocator: Allocator, comptime T: type, vt: *const Vtable, ta: bool, tb: bool, conj: Conj, m: usize, n: usize, k: usize, alpha: T, a: [*]const T, lda: usize, b: [*]const T, ldb: usize, beta: T, c: [*]T, ldc: usize) Error!void {
    if (!isOne(T, beta)) scale(T, vt, m, n, beta, c, ldc);
    if (k == 0 or isZero(T, alpha)) return;
    const x: Matrix(T) = .{ .p = @constCast(a), .ld = lda };
    const y: Matrix(T) = .{ .p = @constCast(b), .ld = ldb };
    switch (ta) {
        inline else => |cta| switch (tb) {
            inline else => |ctb| switch (if (comptime isComplex(T)) conj else .n) {
                inline else => |cj| try run(allocator, T, vt, .{ .rows = if (cta) .t else .n, .cols = if (ctb) .t else .n, .conj = cj }, m, n, k, alpha, x, y, null, c, ldc),
            },
        },
    }
}

/// `C := alpha*S*B + beta*C` (`left`) or `alpha*B*S + beta*C`, `S` the
/// symmetric (`herm`: Hermitian) matrix whose `upper` or lower triangle `a`
/// holds; `symm_k.c` and `zhemm_k.c` over `level3.c`. `m`, `n` positive.
pub fn symm(allocator: Allocator, comptime T: type, vt: *const Vtable, comptime herm: bool, left: bool, upper: bool, m: usize, n: usize, alpha: T, a: [*]const T, lda: usize, b: [*]const T, ldb: usize, beta: T, c: [*]T, ldc: usize) Error!void {
    if (!isOne(T, beta)) scale(T, vt, m, n, beta, c, ldc);
    if (isZero(T, alpha)) return;
    const s: Matrix(T) = .{ .p = @constCast(a), .ld = lda };
    const g: Matrix(T) = .{ .p = @constCast(b), .ld = ldb };
    const tri: Pack = if (upper) .upper else .lower;
    switch (left) {
        inline else => |l| switch (tri) {
            inline .upper, .lower => |t| {
                // `hemm`'s right side packs `S^T`: the kernel conjugates it (`NC`)
                const d: Desc = if (l) .{ .rows = t, .cols = .n, .herm = herm } else .{ .rows = .n, .cols = t, .herm = herm, .conj = if (herm) .r else .n };
                try run(allocator, T, vt, d, m, n, if (l) m else n, alpha, if (l) s else g, if (l) g else s, null, c, ldc);
            },
            else => unreachable,
        },
    }
}

/// `C := beta*C` on the `upper` or lower triangle of the `n x n` `C`; with
/// `herm` `beta` is real and the diagonal comes out real (`zherk_k.c`'s
/// `syrk_beta`). A zero `beta` never reads `C`.
fn scaleTriangle(comptime T: type, vt: *const Vtable, comptime herm: bool, upper: bool, n: usize, beta: T, c: [*]T, ldc: usize) void {
    for (0..n) |j| {
        const ib = if (upper) 0 else j;
        const rows = if (upper) j + 1 else n - j;
        const col = c + ib + j * ldc;
        if (comptime herm) {
            // a real `beta` scales each part, as the reference's real times complex
            const F = Float(T);
            _ = portable.kernel(F, vt, .gemm_beta)(@intCast(2 * rows), 1, 0, beta.re, raw(T, col), 0, raw(T, col), 0, raw(T, col), @intCast(2 * rows));
            c[j + j * ldc].im = 0;
        } else scale(T, vt, rows, 1, beta, col, ldc);
    }
}

/// `C := alpha*op(A)*op(A)' + beta*C` on the `upper` or lower triangle of
/// the `n x n` `C`, `op(A)` being `n x k`; `'` is the transpose, or with
/// `herm` the conjugate transpose (`alpha`, `beta` then real, given as
/// complex with zero imaginary parts). `trans` takes `op(A) = A'`. With `b`,
/// `syr2k`: `alpha*op(A)*op(B)' + alpha2*op(B)*op(A)' + beta*C`, `alpha2`
/// being `alpha` (or its conjugate for `herm`). `level3_syrk.c` and
/// `level3_syr2k.c`, over `syrk_k.c`'s triangle. `n` positive.
pub fn syrk(allocator: Allocator, comptime T: type, vt: *const Vtable, comptime herm: bool, upper: bool, trans: bool, n: usize, k: usize, alpha: T, a: [*]const T, lda: usize, b: ?[*]const T, ldb: usize, beta: T, c: [*]T, ldc: usize) Error!void {
    if (!isOne(T, beta)) scaleTriangle(T, vt, herm, upper, n, beta, c, ldc);
    if (k == 0 or isZero(T, alpha)) return;
    const x: Matrix(T) = .{ .p = @constCast(a), .ld = lda };
    const y: Matrix(T) = if (b) |p| .{ .p = @constCast(p), .ld = ldb } else x;
    const alpha2: ?T = if (b == null) null else if (comptime herm) alpha.conjugate() else alpha;
    switch (upper) {
        inline else => |u| switch (trans) {
            // `herk`'s `op(A)'` conjugates the columns (N) or `op(A)` the rows (C)
            inline else => |t| try run(allocator, T, vt, .{
                .rows = if (t) .t else .n,
                .cols = if (t) .n else .t,
                .conj = if (herm) (if (t) .l else .r) else .n,
                .region = if (u) .upper else .lower,
                .real_diag = herm,
            }, n, n, k, alpha, x, y, alpha2, c, ldc),
        },
    }
}

// ---- trmm and trsm -------------------------------------------------------------

/// A copy of the triangle (`trsm_iunucopy` and the like); `offset` places
/// the block's lines against its depth.
inline fn triCopy(comptime T: type, vt: *const Vtable, comptime op: Op, m: usize, n: usize, a: [*]T, lda: usize, offset: isize, dst: [*]T) void {
    if (comptime isComplex(T)) {
        _ = portable.kernel(T, vt, op)(@intCast(m), @intCast(n), raw(T, a), @intCast(lda), offset, raw(T, dst));
    } else {
        _ = portable.kernel(T, vt, op)(@intCast(m), @intCast(n), a, @intCast(lda), offset, dst);
    }
}

/// A `trmm` or `trsm` kernel: `alpha` for `trmm` (`trsm`'s is a dummy).
inline fn triKernel(comptime T: type, vt: *const Vtable, comptime op: Op, m: usize, n: usize, k: usize, alpha: T, sa: [*]T, sb: [*]T, c: [*]T, ldc: usize, offset: isize) void {
    if (comptime isComplex(T)) {
        _ = portable.kernel(T, vt, op)(@intCast(m), @intCast(n), @intCast(k), alpha.re, alpha.im, raw(T, sa), raw(T, sb), raw(T, c), @intCast(ldc), offset);
    } else {
        _ = portable.kernel(T, vt, op)(@intCast(m), @intCast(n), @intCast(k), alpha, sa, sb, c, @intCast(ldc), offset);
    }
}

fn one(comptime T: type, v: Float(T)) T {
    return if (comptime isComplex(T)) .init(v, 0) else v;
}

/// The options of a `trmm` or `trsm`, fixed at compile time.
const Tri = struct {
    left: bool,
    upper: bool,
    trans: bool,
    conj: bool,
    unit: bool,

    /// `<routine>_<i|o><u|l><n|t><u|n>copy`.
    fn copyOp(comptime t: Tri, comptime routine: []const u8, comptime io: []const u8, comptime ul: []const u8, comptime nt: []const u8) Op {
        return @field(Op, routine ++ "_" ++ io ++ ul ++ nt ++ (if (t.unit) "u" else "n") ++ "copy");
    }
    /// `<routine>_kernel_<L|R><N|T>`, or `R` and `C` for `N` and `T` conjugated.
    fn kernelOp(comptime t: Tri, comptime routine: []const u8, comptime n: bool) Op {
        const k: []const u8 = if (t.conj) (if (n) "R" else "C") else (if (n) "N" else "T");
        return @field(Op, routine ++ "_kernel_" ++ (if (t.left) "L" else "R") ++ k);
    }
    /// The gemm kernel's conjugation: the triangle's side.
    fn gemmConj(comptime t: Tri) Conj {
        return if (!t.conj) .n else if (t.left) .l else .r;
    }
};

/// Packing space for `trmm` and `trsm`: a row block of `P` by a slice of `Q`,
/// and that slice by `R` columns (the reference's kernels take all of it).
fn triWork(allocator: Allocator, comptime T: type, blk: Blocking, m: usize, n: usize, ka: usize, stack: *align(Work(T).alignment) [stack_bytes]u8, heap: *?[]align(Work(T).alignment) T) Error!Work(T) {
    const depth = @min(ka, blk.q);
    const lens: [3]usize = .{ @min(m, blk.p) * depth, depth * @min(n, blk.r), 0 };
    const len = Work(T).total(lens);
    const base: [*]T = if (len * @sizeOf(T) <= stack_bytes) stk: {
        if (sysinfo.is_test) @memset(stack, 0xff);
        break :stk @ptrCast(stack);
    } else hp: {
        heap.* = try allocator.alignedAlloc(T, .fromByteUnits(Work(T).alignment), len);
        break :hp heap.*.?.ptr;
    };
    return Work(T).at(base, lens);
}

/// `B := alpha*op(A)*B` (`left`) or `alpha*B*op(A)`, `A` triangular
/// (`upper` or lower, `unit` diagonal), `op` transposing (`trans`) and
/// conjugating (`conj`, complex) as asked; `B` is `m x n`. `B` is scaled by
/// `alpha` first, never read for a zero one (`trmm_L.c`, `trmm_R.c`).
pub fn trmm(allocator: Allocator, comptime T: type, vt: *const Vtable, left: bool, upper: bool, trans: bool, conj: bool, unit: bool, m: usize, n: usize, alpha: T, a: [*]const T, lda: usize, b: [*]T, ldb: usize) Error!void {
    return triangular(allocator, T, vt, false, left, upper, trans, conj, unit, m, n, alpha, a, lda, b, ldb);
}

/// `B := alpha*inv(op(A))*B` (`left`) or `alpha*B*inv(op(A))`, as `trmm`
/// takes its options (`trsm_L.c`, `trsm_R.c`).
pub fn trsm(allocator: Allocator, comptime T: type, vt: *const Vtable, left: bool, upper: bool, trans: bool, conj: bool, unit: bool, m: usize, n: usize, alpha: T, a: [*]const T, lda: usize, b: [*]T, ldb: usize) Error!void {
    return triangular(allocator, T, vt, true, left, upper, trans, conj, unit, m, n, alpha, a, lda, b, ldb);
}

fn triangular(allocator: Allocator, comptime T: type, vt: *const Vtable, comptime solve: bool, left: bool, upper: bool, trans: bool, conj: bool, unit: bool, m: usize, n: usize, alpha: T, a: [*]const T, lda: usize, b: [*]T, ldb: usize) Error!void {
    // the reference passes `alpha` as `beta`
    if (!isOne(T, alpha)) scale(T, vt, m, n, alpha, b, ldb);
    if (isZero(T, alpha)) return;
    const blk = Blocking.of(T, vt);
    var stack: [stack_bytes]u8 align(Work(T).alignment) = undefined;
    var heap: ?[]align(Work(T).alignment) T = null;
    defer if (heap) |h| allocator.free(h);
    const work = try triWork(allocator, T, blk, m, n, if (left) m else n, &stack, &heap);
    const ap: [*]T = @constCast(a);
    const cj = isComplex(T) and conj;
    switch (left) {
        inline else => |l| switch (upper) {
            inline else => |u| switch (trans) {
                inline else => |t| switch (cj) {
                    inline else => |c| switch (unit) {
                        inline else => |un| {
                            const opt: Tri = .{ .left = l, .upper = u, .trans = t, .conj = c and comptime isComplex(T), .unit = un };
                            if (solve) {
                                if (l) trsmL(T, vt, opt, blk, m, n, ap, lda, b, ldb, work) else trsmR(T, vt, opt, blk, m, n, ap, lda, b, ldb, work);
                            } else {
                                if (l) trmmL(T, vt, opt, blk, m, n, ap, lda, b, ldb, work) else trmmR(T, vt, opt, blk, m, n, ap, lda, b, ldb, work);
                            }
                        },
                    },
                },
            },
        },
    }
}

/// `min_i` of the `trmm` drivers: at most `P`, whole row panels past one.
fn trmmRows(blk: Blocking, rest: usize) usize {
    const v = @min(rest, blk.p);
    return if (v > blk.um) v / blk.um * blk.um else v;
}

fn trmmL(comptime T: type, vt: *const Vtable, comptime o: Tri, blk: Blocking, m: usize, n: usize, a: [*]T, lda: usize, b: [*]T, ldb: usize, w: Work(T)) void {
    const p1 = one(T, 1);
    const gc = comptime o.gemmConj();
    const icopy: Op = if (o.trans) .gemm_incopy else .gemm_itcopy;
    // `op(A)` upper: the slices go down, each row taking the depths from it (`LN`)
    const down = o.upper != o.trans;
    const tcopy = comptime if (down) o.copyOp("trmm", "i", if (o.trans) "l" else "u", if (o.trans) "n" else "t") else o.copyOp("trmm", "i", if (o.trans) "u" else "l", if (o.trans) "n" else "t");
    const tkern = comptime o.kernelOp("trmm", down);
    var js: usize = 0;
    while (js < n) : (js += blk.r) {
        const min_j = @min(n - js, blk.r);
        if (down) {
            var min_l = @min(m, blk.q);
            var min_i = trmmRows(blk, min_l);
            placed(T, vt, tcopy, min_l, min_i, a, lda, 0, 0, w.sa);
            var jjs = js;
            while (jjs < js + min_j) {
                const min_jj = blk.chunk(js + min_j - jjs);
                const pb = w.sb + min_l * (jjs - js);
                copy(T, vt, .gemm_oncopy, min_l, min_jj, b + jjs * ldb, ldb, pb);
                triKernel(T, vt, tkern, min_i, min_jj, min_l, p1, w.sa, pb, b + jjs * ldb, ldb, 0);
                jjs += min_jj;
            }
            var is = min_i;
            while (is < min_l) : (is += min_i) {
                min_i = trmmRows(blk, min_l - is);
                placed(T, vt, tcopy, min_l, min_i, a, lda, 0, is, w.sa);
                triKernel(T, vt, tkern, min_i, min_j, min_l, p1, w.sa, w.sb, b + is + js * ldb, ldb, @intCast(is));
            }
            var ls = min_l;
            while (ls < m) : (ls += blk.q) {
                min_l = @min(m - ls, blk.q);
                min_i = trmmRows(blk, ls);
                copy(T, vt, icopy, min_l, min_i, if (o.trans) a + ls else a + ls * lda, lda, w.sa);
                jjs = js;
                while (jjs < js + min_j) {
                    const min_jj = blk.chunk(js + min_j - jjs);
                    const pb = w.sb + min_l * (jjs - js);
                    copy(T, vt, .gemm_oncopy, min_l, min_jj, b + ls + jjs * ldb, ldb, pb);
                    kernel(T, vt, gc, min_i, min_jj, min_l, p1, w.sa, pb, b + jjs * ldb, ldb);
                    jjs += min_jj;
                }
                is = min_i;
                while (is < ls) : (is += min_i) {
                    min_i = trmmRows(blk, ls - is);
                    copy(T, vt, icopy, min_l, min_i, if (o.trans) a + ls + is * lda else a + is + ls * lda, lda, w.sa);
                    kernel(T, vt, gc, min_i, min_j, min_l, p1, w.sa, w.sb, b + is + js * ldb, ldb);
                }
                is = ls;
                while (is < ls + min_l) : (is += min_i) {
                    min_i = trmmRows(blk, ls + min_l - is);
                    placed(T, vt, tcopy, min_l, min_i, a, lda, ls, is, w.sa);
                    triKernel(T, vt, tkern, min_i, min_j, min_l, p1, w.sa, w.sb, b + is + js * ldb, ldb, @intCast(is - ls));
                }
            }
        } else {
            // `op(A)` lower: the slices go up, each row taking the depths up to it (`LT`)
            var min_l = @min(m, blk.q);
            var min_i = trmmRows(blk, min_l);
            placed(T, vt, tcopy, min_l, min_i, a, lda, m - min_l, m - min_l, w.sa);
            var jjs = js;
            while (jjs < js + min_j) {
                const min_jj = blk.chunk(js + min_j - jjs);
                const pb = w.sb + min_l * (jjs - js);
                copy(T, vt, .gemm_oncopy, min_l, min_jj, b + (m - min_l) + jjs * ldb, ldb, pb);
                triKernel(T, vt, tkern, min_i, min_jj, min_l, p1, w.sa, pb, b + (m - min_l) + jjs * ldb, ldb, 0);
                jjs += min_jj;
            }
            var is = m - min_l + min_i;
            while (is < m) : (is += min_i) {
                min_i = trmmRows(blk, m - is);
                placed(T, vt, tcopy, min_l, min_i, a, lda, m - min_l, is, w.sa);
                triKernel(T, vt, tkern, min_i, min_j, min_l, p1, w.sa, w.sb, b + is + js * ldb, ldb, @intCast(is - (m - min_l)));
            }
            var ls = m - min_l;
            while (ls > 0) : (ls -= min_l) {
                min_l = @min(ls, blk.q);
                const l0 = ls - min_l;
                min_i = trmmRows(blk, min_l);
                placed(T, vt, tcopy, min_l, min_i, a, lda, l0, l0, w.sa);
                jjs = js;
                while (jjs < js + min_j) {
                    const min_jj = blk.chunk(js + min_j - jjs);
                    const pb = w.sb + min_l * (jjs - js);
                    copy(T, vt, .gemm_oncopy, min_l, min_jj, b + l0 + jjs * ldb, ldb, pb);
                    triKernel(T, vt, tkern, min_i, min_jj, min_l, p1, w.sa, pb, b + l0 + jjs * ldb, ldb, 0);
                    jjs += min_jj;
                }
                is = l0 + min_i;
                while (is < ls) : (is += min_i) {
                    min_i = trmmRows(blk, ls - is);
                    placed(T, vt, tcopy, min_l, min_i, a, lda, l0, is, w.sa);
                    triKernel(T, vt, tkern, min_i, min_j, min_l, p1, w.sa, w.sb, b + is + js * ldb, ldb, @intCast(is - l0));
                }
                is = ls;
                while (is < m) : (is += min_i) {
                    min_i = trmmRows(blk, m - is);
                    copy(T, vt, icopy, min_l, min_i, if (o.trans) a + l0 + is * lda else a + is + l0 * lda, lda, w.sa);
                    kernel(T, vt, gc, min_i, min_j, min_l, p1, w.sa, w.sb, b + is + js * ldb, ldb);
                }
            }
        }
    }
}

fn trmmR(comptime T: type, vt: *const Vtable, comptime o: Tri, blk: Blocking, m: usize, n: usize, a: [*]T, lda: usize, b: [*]T, ldb: usize, w: Work(T)) void {
    const p1 = one(T, 1);
    const gc = comptime o.gemmConj();
    const ocopy: Op = if (o.trans) .gemm_otcopy else .gemm_oncopy;
    // `op(A)` lower: the columns go right, each taking the depths from it (`RT`)
    const right = o.upper == o.trans;
    const tcopy = comptime if (right) o.copyOp("trmm", "o", if (o.trans) "u" else "l", if (o.trans) "t" else "n") else o.copyOp("trmm", "o", if (o.trans) "l" else "u", if (o.trans) "t" else "n");
    const tkern = comptime o.kernelOp("trmm", !right);
    // `op(A)` from row `l`, column `j`, as `ocopy` takes it
    const at = struct {
        inline fn f(p: [*]T, ld: usize, l: usize, j: usize) [*]T {
            return if (o.trans) p + j + l * ld else p + l + j * ld;
        }
    }.f;
    const min_i0 = @min(m, blk.p);
    if (right) {
        var js: usize = 0;
        while (js < n) : (js += blk.r) {
            const min_j = @min(n - js, blk.r);
            var ls = js;
            while (ls < js + min_j) : (ls += blk.q) {
                const min_l = @min(js + min_j - ls, blk.q);
                copy(T, vt, .gemm_itcopy, min_l, min_i0, b + ls * ldb, ldb, w.sa);
                var jjs: usize = 0;
                while (jjs < ls - js) {
                    const min_jj = blk.chunk(ls - js - jjs);
                    copy(T, vt, ocopy, min_l, min_jj, at(a, lda, ls, js + jjs), lda, w.sb + min_l * jjs);
                    kernel(T, vt, gc, min_i0, min_jj, min_l, p1, w.sa, w.sb + min_l * jjs, b + (js + jjs) * ldb, ldb);
                    jjs += min_jj;
                }
                jjs = 0;
                while (jjs < min_l) {
                    const min_jj = blk.chunk(min_l - jjs);
                    const pb = w.sb + min_l * (ls - js + jjs);
                    placed(T, vt, tcopy, min_l, min_jj, a, lda, ls, ls + jjs, pb);
                    triKernel(T, vt, tkern, min_i0, min_jj, min_l, p1, w.sa, pb, b + (ls + jjs) * ldb, ldb, -@as(isize, @intCast(jjs)));
                    jjs += min_jj;
                }
                var is = min_i0;
                while (is < m) : (is += blk.p) {
                    const min_i = @min(m - is, blk.p);
                    copy(T, vt, .gemm_itcopy, min_l, min_i, b + is + ls * ldb, ldb, w.sa);
                    kernel(T, vt, gc, min_i, ls - js, min_l, p1, w.sa, w.sb, b + is + js * ldb, ldb);
                    triKernel(T, vt, tkern, min_i, min_l, min_l, p1, w.sa, w.sb + (ls - js) * min_l, b + is + ls * ldb, ldb, 0);
                }
            }
            ls = js + min_j;
            while (ls < n) : (ls += blk.q) {
                const min_l = @min(n - ls, blk.q);
                copy(T, vt, .gemm_itcopy, min_l, min_i0, b + ls * ldb, ldb, w.sa);
                var jjs = js;
                while (jjs < js + min_j) {
                    const min_jj = blk.chunk(js + min_j - jjs);
                    const pb = w.sb + min_l * (jjs - js);
                    copy(T, vt, ocopy, min_l, min_jj, at(a, lda, ls, jjs), lda, pb);
                    kernel(T, vt, gc, min_i0, min_jj, min_l, p1, w.sa, pb, b + jjs * ldb, ldb);
                    jjs += min_jj;
                }
                var is = min_i0;
                while (is < m) : (is += blk.p) {
                    const min_i = @min(m - is, blk.p);
                    copy(T, vt, .gemm_itcopy, min_l, min_i, b + is + ls * ldb, ldb, w.sa);
                    kernel(T, vt, gc, min_i, min_j, min_l, p1, w.sa, w.sb, b + is + js * ldb, ldb);
                }
            }
        }
    } else {
        // `op(A)` upper: the columns go left, each taking the depths up to it (`RN`)
        var js = n;
        while (js > 0) : (js -= @min(js, blk.r)) {
            const min_j = @min(js, blk.r);
            const j0 = js - min_j;
            var start_ls = j0;
            while (start_ls + blk.q < js) start_ls += blk.q;
            var ls = start_ls;
            while (true) {
                const min_l = @min(js - ls, blk.q);
                copy(T, vt, .gemm_itcopy, min_l, min_i0, b + ls * ldb, ldb, w.sa);
                var jjs: usize = 0;
                while (jjs < min_l) {
                    const min_jj = blk.chunk(min_l - jjs);
                    placed(T, vt, tcopy, min_l, min_jj, a, lda, ls, ls + jjs, w.sb + min_l * jjs);
                    triKernel(T, vt, tkern, min_i0, min_jj, min_l, p1, w.sa, w.sb + min_l * jjs, b + (ls + jjs) * ldb, ldb, -@as(isize, @intCast(jjs)));
                    jjs += min_jj;
                }
                jjs = 0;
                while (jjs < js - ls - min_l) {
                    const min_jj = blk.chunk(js - ls - min_l - jjs);
                    const pb = w.sb + min_l * (min_l + jjs);
                    copy(T, vt, ocopy, min_l, min_jj, at(a, lda, ls, ls + min_l + jjs), lda, pb);
                    kernel(T, vt, gc, min_i0, min_jj, min_l, p1, w.sa, pb, b + (ls + min_l + jjs) * ldb, ldb);
                    jjs += min_jj;
                }
                var is = min_i0;
                while (is < m) : (is += blk.p) {
                    const min_i = @min(m - is, blk.p);
                    copy(T, vt, .gemm_itcopy, min_l, min_i, b + is + ls * ldb, ldb, w.sa);
                    triKernel(T, vt, tkern, min_i, min_l, min_l, p1, w.sa, w.sb, b + is + ls * ldb, ldb, 0);
                    if (js - ls - min_l > 0) kernel(T, vt, gc, min_i, js - ls - min_l, min_l, p1, w.sa, w.sb + min_l * min_l, b + is + (ls + min_l) * ldb, ldb);
                }
                if (ls < j0 + blk.q) break;
                ls -= blk.q;
            }
            ls = 0;
            while (ls < j0) : (ls += blk.q) {
                const min_l = @min(j0 - ls, blk.q);
                copy(T, vt, .gemm_itcopy, min_l, min_i0, b + ls * ldb, ldb, w.sa);
                var jjs = j0;
                while (jjs < js) {
                    const min_jj = blk.chunk(js - jjs);
                    const pb = w.sb + min_l * (jjs - j0);
                    copy(T, vt, ocopy, min_l, min_jj, at(a, lda, ls, jjs), lda, pb);
                    kernel(T, vt, gc, min_i0, min_jj, min_l, p1, w.sa, pb, b + jjs * ldb, ldb);
                    jjs += min_jj;
                }
                var is = min_i0;
                while (is < m) : (is += blk.p) {
                    const min_i = @min(m - is, blk.p);
                    copy(T, vt, .gemm_itcopy, min_l, min_i, b + is + ls * ldb, ldb, w.sa);
                    kernel(T, vt, gc, min_i, min_j, min_l, p1, w.sa, w.sb, b + is + j0 * ldb, ldb);
                }
            }
        }
    }
}

fn trsmL(comptime T: type, vt: *const Vtable, comptime o: Tri, blk: Blocking, m: usize, n: usize, a: [*]T, lda: usize, b: [*]T, ldb: usize, w: Work(T)) void {
    const m1 = one(T, -1);
    const gc = comptime o.gemmConj();
    const icopy: Op = if (o.trans) .gemm_incopy else .gemm_itcopy;
    // `op(A)` lower: forward from the top (`LT`), else back from the bottom (`LN`)
    const forward = o.upper == o.trans;
    const tcopy = comptime if (forward) o.copyOp("trsm", "i", if (o.trans) "u" else "l", if (o.trans) "n" else "t") else o.copyOp("trsm", "i", if (o.trans) "l" else "u", if (o.trans) "n" else "t");
    const tkern = comptime o.kernelOp("trsm", !forward);
    // `op(A)` from row `i`, column `l`, as the `i` copies take it
    const at = struct {
        inline fn f(p: [*]T, ld: usize, i: usize, l: usize) [*]T {
            return if (o.trans) p + l + i * ld else p + i + l * ld;
        }
    }.f;
    var js: usize = 0;
    while (js < n) : (js += blk.r) {
        const min_j = @min(n - js, blk.r);
        if (forward) {
            var ls: usize = 0;
            while (ls < m) : (ls += blk.q) {
                const min_l = @min(m - ls, blk.q);
                var min_i = @min(min_l, blk.p);
                triCopy(T, vt, tcopy, min_l, min_i, at(a, lda, ls, ls), lda, 0, w.sa);
                var jjs = js;
                while (jjs < js + min_j) {
                    const min_jj = blk.chunk(js + min_j - jjs);
                    const pb = w.sb + min_l * (jjs - js);
                    copy(T, vt, .gemm_oncopy, min_l, min_jj, b + ls + jjs * ldb, ldb, pb);
                    triKernel(T, vt, tkern, min_i, min_jj, min_l, m1, w.sa, pb, b + ls + jjs * ldb, ldb, 0);
                    jjs += min_jj;
                }
                var is = ls + min_i;
                while (is < ls + min_l) : (is += blk.p) {
                    min_i = @min(ls + min_l - is, blk.p);
                    triCopy(T, vt, tcopy, min_l, min_i, at(a, lda, is, ls), lda, @intCast(is - ls), w.sa);
                    triKernel(T, vt, tkern, min_i, min_j, min_l, m1, w.sa, w.sb, b + is + js * ldb, ldb, @intCast(is - ls));
                }
                is = ls + min_l;
                while (is < m) : (is += blk.p) {
                    min_i = @min(m - is, blk.p);
                    copy(T, vt, icopy, min_l, min_i, at(a, lda, is, ls), lda, w.sa);
                    kernel(T, vt, gc, min_i, min_j, min_l, m1, w.sa, w.sb, b + is + js * ldb, ldb);
                }
            }
        } else {
            var ls = m;
            while (ls > 0) : (ls -= @min(ls, blk.q)) {
                const min_l = @min(ls, blk.q);
                const l0 = ls - min_l;
                var start_is = l0;
                while (start_is + blk.p < ls) start_is += blk.p;
                var min_i = @min(ls - start_is, blk.p);
                triCopy(T, vt, tcopy, min_l, min_i, at(a, lda, start_is, l0), lda, @intCast(start_is - l0), w.sa);
                var jjs = js;
                while (jjs < js + min_j) {
                    const min_jj = blk.chunk(js + min_j - jjs);
                    const pb = w.sb + min_l * (jjs - js);
                    copy(T, vt, .gemm_oncopy, min_l, min_jj, b + l0 + jjs * ldb, ldb, pb);
                    triKernel(T, vt, tkern, min_i, min_jj, min_l, m1, w.sa, pb, b + start_is + jjs * ldb, ldb, @intCast(start_is - l0));
                    jjs += min_jj;
                }
                var is = start_is;
                while (is >= l0 + blk.p) {
                    is -= blk.p;
                    min_i = @min(ls - is, blk.p);
                    triCopy(T, vt, tcopy, min_l, min_i, at(a, lda, is, l0), lda, @intCast(is - l0), w.sa);
                    triKernel(T, vt, tkern, min_i, min_j, min_l, m1, w.sa, w.sb, b + is + js * ldb, ldb, @intCast(is - l0));
                }
                is = 0;
                while (is < l0) : (is += blk.p) {
                    min_i = @min(l0 - is, blk.p);
                    copy(T, vt, icopy, min_l, min_i, at(a, lda, is, l0), lda, w.sa);
                    kernel(T, vt, gc, min_i, min_j, min_l, m1, w.sa, w.sb, b + is + js * ldb, ldb);
                }
            }
        }
    }
}

fn trsmR(comptime T: type, vt: *const Vtable, comptime o: Tri, blk: Blocking, m: usize, n: usize, a: [*]T, lda: usize, b: [*]T, ldb: usize, w: Work(T)) void {
    const m1 = one(T, -1);
    const gc = comptime o.gemmConj();
    const ocopy: Op = if (o.trans) .gemm_otcopy else .gemm_oncopy;
    // `op(A)` upper: forward from the left (`RN`), else back from the right (`RT`)
    const forward = o.upper != o.trans;
    const tcopy = comptime if (forward) o.copyOp("trsm", "o", if (o.trans) "l" else "u", if (o.trans) "t" else "n") else o.copyOp("trsm", "o", if (o.trans) "u" else "l", if (o.trans) "t" else "n");
    const tkern = comptime o.kernelOp("trsm", forward);
    // `op(A)` from row `l`, column `j`, as `ocopy` takes it
    const at = struct {
        inline fn f(p: [*]T, ld: usize, l: usize, j: usize) [*]T {
            return if (o.trans) p + j + l * ld else p + l + j * ld;
        }
    }.f;
    const min_i0 = @min(m, blk.p);
    if (forward) {
        var js: usize = 0;
        while (js < n) : (js += blk.r) {
            const min_j = @min(n - js, blk.r);
            var ls: usize = 0;
            while (ls < js) : (ls += blk.q) {
                const min_l = @min(js - ls, blk.q);
                copy(T, vt, .gemm_itcopy, min_l, min_i0, b + ls * ldb, ldb, w.sa);
                var jjs = js;
                while (jjs < js + min_j) {
                    const min_jj = blk.chunk(js + min_j - jjs);
                    const pb = w.sb + min_l * (jjs - js);
                    copy(T, vt, ocopy, min_l, min_jj, at(a, lda, ls, jjs), lda, pb);
                    kernel(T, vt, gc, min_i0, min_jj, min_l, m1, w.sa, pb, b + jjs * ldb, ldb);
                    jjs += min_jj;
                }
                var is = min_i0;
                while (is < m) : (is += blk.p) {
                    const min_i = @min(m - is, blk.p);
                    copy(T, vt, .gemm_itcopy, min_l, min_i, b + is + ls * ldb, ldb, w.sa);
                    kernel(T, vt, gc, min_i, min_j, min_l, m1, w.sa, w.sb, b + is + js * ldb, ldb);
                }
            }
            ls = js;
            while (ls < js + min_j) : (ls += blk.q) {
                const min_l = @min(js + min_j - ls, blk.q);
                copy(T, vt, .gemm_itcopy, min_l, min_i0, b + ls * ldb, ldb, w.sa);
                triCopy(T, vt, tcopy, min_l, min_l, at(a, lda, ls, ls), lda, 0, w.sb);
                triKernel(T, vt, tkern, min_i0, min_l, min_l, m1, w.sa, w.sb, b + ls * ldb, ldb, 0);
                var jjs: usize = 0;
                while (jjs < js + min_j - ls - min_l) {
                    const min_jj = blk.chunk(js + min_j - ls - min_l - jjs);
                    const pb = w.sb + min_l * (min_l + jjs);
                    copy(T, vt, ocopy, min_l, min_jj, at(a, lda, ls, ls + min_l + jjs), lda, pb);
                    kernel(T, vt, gc, min_i0, min_jj, min_l, m1, w.sa, pb, b + (ls + min_l + jjs) * ldb, ldb);
                    jjs += min_jj;
                }
                var is = min_i0;
                while (is < m) : (is += blk.p) {
                    const min_i = @min(m - is, blk.p);
                    copy(T, vt, .gemm_itcopy, min_l, min_i, b + is + ls * ldb, ldb, w.sa);
                    triKernel(T, vt, tkern, min_i, min_l, min_l, m1, w.sa, w.sb, b + is + ls * ldb, ldb, 0);
                    kernel(T, vt, gc, min_i, js + min_j - ls - min_l, min_l, m1, w.sa, w.sb + min_l * min_l, b + is + (ls + min_l) * ldb, ldb);
                }
            }
        }
    } else {
        var js = n;
        while (js > 0) : (js -= @min(js, blk.r)) {
            const min_j = @min(js, blk.r);
            const j0 = js - min_j;
            var ls = js;
            while (ls < n) : (ls += blk.q) {
                const min_l = @min(n - ls, blk.q);
                copy(T, vt, .gemm_itcopy, min_l, min_i0, b + ls * ldb, ldb, w.sa);
                var jjs = j0;
                while (jjs < js) {
                    const min_jj = blk.chunk(js - jjs);
                    const pb = w.sb + min_l * (jjs - j0);
                    copy(T, vt, ocopy, min_l, min_jj, at(a, lda, ls, jjs), lda, pb);
                    kernel(T, vt, gc, min_i0, min_jj, min_l, m1, w.sa, pb, b + jjs * ldb, ldb);
                    jjs += min_jj;
                }
                var is = min_i0;
                while (is < m) : (is += blk.p) {
                    const min_i = @min(m - is, blk.p);
                    copy(T, vt, .gemm_itcopy, min_l, min_i, b + is + ls * ldb, ldb, w.sa);
                    kernel(T, vt, gc, min_i, min_j, min_l, m1, w.sa, w.sb, b + is + j0 * ldb, ldb);
                }
            }
            var start_ls = j0;
            while (start_ls + blk.q < js) start_ls += blk.q;
            ls = start_ls;
            while (true) {
                const min_l = @min(js - ls, blk.q);
                const tb = w.sb + min_l * (ls - j0);
                copy(T, vt, .gemm_itcopy, min_l, min_i0, b + ls * ldb, ldb, w.sa);
                triCopy(T, vt, tcopy, min_l, min_l, at(a, lda, ls, ls), lda, 0, tb);
                triKernel(T, vt, tkern, min_i0, min_l, min_l, m1, w.sa, tb, b + ls * ldb, ldb, 0);
                var jjs: usize = 0;
                while (jjs < ls - j0) {
                    const min_jj = blk.chunk(ls - j0 - jjs);
                    const pb = w.sb + min_l * jjs;
                    copy(T, vt, ocopy, min_l, min_jj, at(a, lda, ls, j0 + jjs), lda, pb);
                    kernel(T, vt, gc, min_i0, min_jj, min_l, m1, w.sa, pb, b + (j0 + jjs) * ldb, ldb);
                    jjs += min_jj;
                }
                var is = min_i0;
                while (is < m) : (is += blk.p) {
                    const min_i = @min(m - is, blk.p);
                    copy(T, vt, .gemm_itcopy, min_l, min_i, b + is + ls * ldb, ldb, w.sa);
                    triKernel(T, vt, tkern, min_i, min_l, min_l, m1, w.sa, tb, b + is + ls * ldb, ldb, 0);
                    kernel(T, vt, gc, min_i, ls - j0, min_l, m1, w.sa, w.sb, b + is + j0 * ldb, ldb);
                }
                if (ls < j0 + blk.q) break;
                ls -= blk.q;
            }
        }
    }
}
