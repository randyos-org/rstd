//! What the SME2 direct real `gemm` shares between the main library's driver
//! (`driver.zig`) and the `vortexm4` object's kernel (`kernel.zig`): the job,
//! the streaming vector's shape and the blocking model. Ported from OpenBLAS
//! `kernel/arm64/sme2_gemm_impl.h` and `sme2_gemm_tile.h` (`v0.3.34-794-gdbf2c24c9`,
//! BSD-3-Clause, the OpenBLAS Project), the design of Deng, Yang, Fang and
//! Dong, "Demystifying ARM SME to Optimize General Matrix Multiplications",
//! arXiv:2512.21473. Plain Zig, no `asm`: imports stay in `linalg/blas/`.

const std = @import("std");

/// The streaming vector length the kernels are written for, in bytes:
/// `s2_usable` takes 16 `f32` lanes only.
pub const svl_bytes = 64;

/// The tile shape of `F` (`S2_VL`, `S2_NT`, `S2_NR`, `S2_W2`).
pub fn Shape(comptime F: type) type {
    return struct {
        /// Lanes of one streaming vector.
        pub const vl: usize = svl_bytes / @sizeOf(F);
        /// ZA tiles, one per byte of the element: 4 of `f32`, 8 of `f64`.
        pub const nt: usize = @sizeOf(F);
        /// Columns of one row of tiles, the main kernel's width.
        pub const nr: usize = nt * vl;
        /// The half-width kernel's columns.
        pub const w2: usize = nt / 2 * vl;
    };
}

/// One row-major core problem `C = alpha A B + beta C`, `C` of `m x n` (row
/// stride `ldc`), depth `k` (`s2_job`), with its blocking and workspace.
/// `a_cols`: `A` has contiguous columns (`(i, l)` at `i + l * lda`), else
/// rows (`i * lda + l`); `b_cols` likewise for `B`'s `(l, j)`.
pub fn Job(comptime F: type) type {
    return extern struct {
        m: usize,
        n: usize,
        k: usize,
        alpha: F,
        beta: F,
        a: [*]const F,
        b: [*]const F,
        c: [*]F,
        lda: usize,
        ldb: usize,
        ldc: usize,
        a_cols: bool,
        b_cols: bool,
        blk: Blocking,
        /// Packed `A` (`work.a_len` elements), then packed `B`.
        ac: [*]F,
        bc: [*]F,
    };
}

/// Rows, columns and depth of one block (`s2_blocking`).
pub const Blocking = extern struct { mc: usize, nc: usize, kc: usize };

// The model's machine (`S2_L2_BYTES`, `S2_PAGE_BYTES`, `S2_TLB_ENTRIES`:
// Apple M4's SME load bandwidth holds to 8 MB, 16 KB pages).
const l2_bytes = 8 << 20;
const page_bytes = 16384;
const tlb_entries = 160;

fn roundUp(x: usize, m: usize) usize {
    return (x + m - 1) / m * m;
}

/// The deepest `kc` whose panels fit the TLB (`s2_model_kc_max`).
fn kcMax(comptime F: type, mr: usize, nr: usize) usize {
    var best: usize = 16;
    var kc: usize = 16;
    while (kc <= 1 << 16) : (kc += 16) {
        const ta = (mr * kc * @sizeOf(F) + page_bytes - 1) / page_bytes + 1;
        const tb = (nr * kc * @sizeOf(F) + page_bytes - 1) / page_bytes + 1;
        if (ta + 2 * tb + mr < tlb_entries) best = kc else break;
    }
    return best;
}

/// `block` evened out over `total` in multiples of `step` (`s2_balance`).
fn balance(total: usize, block: usize, step: usize) usize {
    if (block >= total) return roundUp(total, step);
    const nblk = (total + block - 1) / block;
    return @min(roundUp((total + nblk - 1) / nblk, step), block);
}

/// The blocking of an `m x n` by `k` job: the paper's model, eqs. 1 to 3
/// (`s2_model`), rounded as `s2_blocking_for` does. Upstream keeps each
/// thread's last eight shapes; this searches every call, tens of
/// microseconds at most against the multiply.
pub fn blocking(comptime F: type, m: usize, n: usize, k: usize) Blocking {
    const S = Shape(F);
    const mr = S.vl;
    const nr = S.nr;
    const budget: usize = l2_bytes / @sizeOf(F);
    const kcap = @min(kcMax(F, mr, nr), roundUp(k, 16));
    const mcap = roundUp(m, mr);
    const ncap = roundUp(n, nr);
    var best: f64 = -1;
    var r: Blocking = .{ .mc = mr, .nc = nr, .kc = 16 };
    var kc: usize = 16;
    while (kc <= kcap) : (kc += 16) {
        var mc: usize = mr;
        while (mc <= mcap) : (mc += mr) {
            if (mc * kc >= budget) break;
            const room = budget - mc * kc;
            var nc = room / (2 * kc + 2 * mc);
            nc = @min(nc / nr * nr, ncap);
            if (nc < nr) break;
            const fm: f64 = @floatFromInt(mc);
            const fn_: f64 = @floatFromInt(nc);
            const fk: f64 = @floatFromInt(kc);
            const cmr = 2.0 * fm * fn_ * fk / (fm * fk + fk * fn_ + 2.0 * fm * fn_);
            if (cmr > best * (1 + 1e-9)) {
                best = cmr;
                r = .{ .mc = mc, .nc = nc, .kc = kc };
            }
        }
    }
    r.kc = @max(balance(k, r.kc, 16), 1);
    r.mc = roundUp(@max(balance(m, r.mc, mr), mr), mr);
    r.nc = roundUp(@max(balance(n, r.nc, nr), nr), nr);
    return r;
}

/// Elements of packed `A` and of the whole workspace for `b` (`s2_work_bytes`,
/// with `A` rounded up to whole 128-byte lines).
pub const Work = struct {
    a_len: usize,
    len: usize,

    pub fn of(comptime F: type, b: Blocking) Work {
        const line = 128 / @sizeOf(F);
        const a_len = roundUp(roundUp(b.mc, 4 * Shape(F).vl) * b.kc, line);
        return .{ .a_len = a_len, .len = a_len + b.nc * b.kc };
    }
};

test "the shapes are upstream's" {
    try std.testing.expectEqual(16, Shape(f32).vl);
    try std.testing.expectEqual(4, Shape(f32).nt);
    try std.testing.expectEqual(64, Shape(f32).nr);
    try std.testing.expectEqual(32, Shape(f32).w2);
    try std.testing.expectEqual(8, Shape(f64).vl);
    try std.testing.expectEqual(8, Shape(f64).nt);
    try std.testing.expectEqual(64, Shape(f64).nr);
    try std.testing.expectEqual(32, Shape(f64).w2);
}

test "the blocking stays in the model's budget and covers small jobs whole" {
    inline for (.{ f32, f64 }) |F| {
        for ([_][3]usize{ .{ 1, 1, 1 }, .{ 64, 64, 64 }, .{ 2048, 2048, 2048 }, .{ 100000, 64, 64 }, .{ 37, 5000, 3 } }) |d| {
            const b = blocking(F, d[0], d[1], d[2]);
            try std.testing.expect(b.mc % Shape(F).vl == 0 and b.nc % Shape(F).nr == 0 and b.kc >= 1);
            try std.testing.expect(Work.of(F, b).len * @sizeOf(F) <= 2 * l2_bytes);
        }
        const s = blocking(F, 64, 64, 64);
        try std.testing.expect(s.mc >= 64 and s.nc >= 64 and s.kc >= 64);
    }
}
