//! The level-3 driver: OpenBLAS's `driver/level3/level3.c`, Goto's loops
//! over the packed kernels of a vtable, single-threaded. `C` is split into
//! `R`-wide column blocks, the depth into `Q`-deep slices and the rows into
//! `P`-tall blocks; each slice of `op(B)` is packed once (`sb`), each block
//! of `op(A)` once per slice (`sa`), and the kernel runs on the pair.
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
/// single row block takes every row, so each chunk is used once).
pub fn Work(comptime T: type) type {
    return struct {
        const Self = @This();
        /// Every buffer starts on a cache line; more than any kernel loads need.
        pub const alignment = 64;
        sa: [*]T,
        sb: [*]T,

        pub fn lens(blk: Blocking, m: usize, n: usize, k: usize) [2]usize {
            const rows = @min(m, blk.p);
            const depth = @min(k, blk.q + blk.um);
            const cols = if (m <= blk.p) @min(n, 3 * blk.un) else @min(n, blk.r);
            return .{ rows * depth, depth * cols };
        }

        /// Elements per step of `sb`'s offset: the fewest that span a whole
        /// number of cache lines (3 for a 12-byte `f80` on x86).
        const step = std.math.lcm(alignment, @sizeOf(T)) / @sizeOf(T);

        /// Elements to allocate for `lens`, with `sb` moved to a cache line.
        pub fn total(l: [2]usize) usize {
            return (l[0] + step - 1) / step * step + l[1];
        }

        pub fn at(base: [*]T, l: [2]usize) Self {
            return .{ .sa = base, .sb = base + (total(l) - l[1]) };
        }
    };
}

/// Bytes of packing space on the stack before the allocator is asked.
pub const stack_bytes = 2048;

pub const Error = Allocator.Error;

/// `C := alpha*op(A)*op(B) + beta*C` with the kernels of `vt`; `ta`, `tb`
/// say whether `op` transposes. Arguments are checked by the caller
/// (`interface.zig`), and `m`, `n` are positive.
pub fn gemm(allocator: Allocator, comptime T: type, vt: *const Vtable, ta: bool, tb: bool, m: usize, n: usize, k: usize, alpha: T, a: [*]const T, lda: usize, b: [*]const T, ldb: usize, beta: T, c: [*]T, ldc: usize) Error!void {
    if (beta != 1) _ = portable.kernel(T, vt, .gemm_beta)(@intCast(m), @intCast(n), 0, beta, c, 0, c, 0, c, @intCast(ldc));
    if (k == 0 or alpha == 0) return;

    const blk = Blocking.of(T, vt);
    const lens = Work(T).lens(blk, m, n, k);
    const len = Work(T).total(lens);
    var stack: [stack_bytes]u8 align(Work(T).alignment) = undefined;
    var heap: ?[]align(Work(T).alignment) T = null;
    defer if (heap) |h| allocator.free(h);
    const base: [*]T = if (len * @sizeOf(T) <= stack_bytes) stk: {
        // tests poison it as `check.zig`'s allocator does the heap: a packing
        // buffer read before it is written shows as NaN
        if (sysinfo.is_test) @memset(&stack, 0xff);
        break :stk @ptrCast(&stack);
    } else blk: {
        heap = try allocator.alignedAlloc(T, .fromByteUnits(Work(T).alignment), len);
        break :blk heap.?.ptr;
    };
    const work = Work(T).at(base, lens);

    switch (ta) {
        inline else => |cta| switch (tb) {
            inline else => |ctb| run(T, vt, cta, ctb, blk, m, n, k, alpha, @constCast(a), lda, @constCast(b), ldb, c, ldc, work),
        },
    }
}

/// The js / ls / is loops of `level3.c`.
fn run(comptime T: type, vt: *const Vtable, comptime ta: bool, comptime tb: bool, blk: Blocking, m: usize, n: usize, k: usize, alpha: T, a: [*]T, lda: usize, b: [*]T, ldb: usize, c: [*]T, ldc: usize, work: Work(T)) void {
    const kernel = portable.kernel(T, vt, .gemm_kernel);
    const icopy = portable.kernel(T, vt, if (ta) .gemm_incopy else .gemm_itcopy);
    const ocopy = portable.kernel(T, vt, if (tb) .gemm_otcopy else .gemm_oncopy);
    const L = BlasLong;
    const at = struct {
        // element (i, l) of op(A), (l, j) of op(B)
        inline fn opA(p: [*]T, ld: usize, i: usize, l: usize) [*]T {
            return if (ta) p + l + i * ld else p + i + l * ld;
        }
        inline fn opB(p: [*]T, ld: usize, l: usize, j: usize) [*]T {
            return if (tb) p + j + l * ld else p + l + j * ld;
        }
    };

    var js: usize = 0;
    while (js < n) : (js += blk.r) {
        const min_j = @min(blk.r, n - js);
        var ls: usize = 0;
        var min_l: usize = 0;
        while (ls < k) : (ls += min_l) {
            min_l = Blocking.split(k - ls, blk.q, blk.um);
            // the first row block packs `B` as it goes
            var min_i = Blocking.split(m, blk.p, blk.um);
            const once = min_i == m;
            _ = icopy(@intCast(min_l), @intCast(min_i), at.opA(a, lda, 0, ls), @intCast(lda), work.sa);
            var jjs = js;
            while (jjs < js + min_j) {
                const min_jj = blk.chunk(js + min_j - jjs);
                const pb = if (once) work.sb else work.sb + min_l * (jjs - js);
                _ = ocopy(@intCast(min_l), @intCast(min_jj), at.opB(b, ldb, ls, jjs), @intCast(ldb), pb);
                _ = kernel(@intCast(min_i), @intCast(min_jj), @intCast(min_l), alpha, work.sa, pb, c + jjs * ldc, @as(L, @intCast(ldc)));
                jjs += min_jj;
            }
            var is = min_i;
            while (is < m) : (is += min_i) {
                min_i = Blocking.split(m - is, blk.p, blk.um);
                _ = icopy(@intCast(min_l), @intCast(min_i), at.opA(a, lda, is, ls), @intCast(lda), work.sa);
                jjs = js;
                while (jjs < js + min_j) {
                    const min_jj = blk.chunk(js + min_j - jjs);
                    _ = kernel(@intCast(min_i), @intCast(min_jj), @intCast(min_l), alpha, work.sa, work.sb + min_l * (jjs - js), c + is + jjs * ldc, @as(L, @intCast(ldc)));
                    jjs += min_jj;
                }
            }
        }
    }
}
