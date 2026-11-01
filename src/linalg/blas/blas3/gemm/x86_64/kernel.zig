//! x86_64's own small-matrix gemm kernels (`small.zig`): Skylake-X's
//! `kernel/x86_64/sgemm_small_kernel_{nn,nt,tn,tt}_skylakex.c` and the
//! `dgemm_` ones, each also built `B0`, on `@Vector`s of 512 bits; Skylake-X,
//! Cooper Lake and Sapphire Rapids take them (`KERNEL.SKYLAKEX`, which the
//! other two include). The gemm tile's registry exports them in place of
//! `small_kernel.zig`'s on those cores (`small_kernel.forCore`). Cooper
//! Lake's `sbgemm` ones (`Cooperlake`) are a port of their own (`ports`).
//!
//! Every block sums each element over the depth in order, fused, then stores
//! `alpha*sum`, rounded, fused with `beta*C` as `vfmadd231p[sd]` does: the
//! generic kernels' bits, but for the dot blocks (`tn`, and `nn`'s last rows
//! when at most half a vector of them is left and the depth is at least two
//! vectors), whose lanes sum every `L`-th product and reduce as upstream's
//! `REDUCE_4` (pairs of neighbors) or `_mm512_reduce_add_p[sd]` (halves).
//! Upstream's `nn` copies those rows to a `malloc`'d buffer first; here they
//! are read in place, the same values. Runs only on AVX-512 hardware; the
//! tests run the same code split into narrower vectors (`small.zig`).

const std = @import("std");
const types = @import("../../../types.zig");
const coretypes = @import("../../../cpu/_coretypes.zig");
const fma = @import("../../../fma.zig");
const registry = @import("../../../asm_registry.zig");
const util = @import("../../../util.zig");
const bf16_kernel = @import("../../../ext/bf16/kernel.zig");
const gk = @import("../kernel.zig");
const small_kernel = @import("../small_kernel.zig");

const BlasLong = types.BlasLong;
const BF16 = types.BF16;
const Slot = util.Slot;

pub const Core = coretypes.X86CoreType;

/// The small kernels of `T`: Skylake-X's for real `T`, the generic ones for
/// complex, which upstream does not tune.
pub fn Kernels(comptime T: type, comptime geo: gk.Geometry) type {
    return if (comptime types.isComplex(T)) small_kernel.Kernels(T, geo) else Skylakex(T);
}

/// The two operand forms: `A` across rows (`nn`, `nt`), or broadcast (`nt`,
/// `tt`, for a row tail), and the dot blocks (`tn`, `nn`).
pub fn Skylakex(comptime T: type) type {
    const cc = types.kernelCc(T);
    // a zmm register's lanes
    const L: usize = 64 / @sizeOf(T);
    const V = @Vector(L, T);
    return struct {
        pub const lanes = L;

        pub fn gemm_small_kernel_nn(m: BlasLong, n: BlasLong, k: BlasLong, a: [*]T, lda: BlasLong, alpha: T, b: [*]T, ldb: BlasLong, beta: T, c: [*]T, ldc: BlasLong) callconv(cc) i32 {
            nn(.init(m, n, k, a, lda, alpha, b, ldb, beta, false, c, ldc));
            return 0;
        }
        pub fn gemm_small_kernel_nt(m: BlasLong, n: BlasLong, k: BlasLong, a: [*]T, lda: BlasLong, alpha: T, b: [*]T, ldb: BlasLong, beta: T, c: [*]T, ldc: BlasLong) callconv(cc) i32 {
            nt(.init(m, n, k, a, lda, alpha, b, ldb, beta, false, c, ldc));
            return 0;
        }
        pub fn gemm_small_kernel_tn(m: BlasLong, n: BlasLong, k: BlasLong, a: [*]T, lda: BlasLong, alpha: T, b: [*]T, ldb: BlasLong, beta: T, c: [*]T, ldc: BlasLong) callconv(cc) i32 {
            const o: Ops = .init(m, n, k, a, lda, alpha, b, ldb, beta, false, c, ldc);
            dotRows(true, o, 0);
            return 0;
        }
        pub fn gemm_small_kernel_tt(m: BlasLong, n: BlasLong, k: BlasLong, a: [*]T, lda: BlasLong, alpha: T, b: [*]T, ldb: BlasLong, beta: T, c: [*]T, ldc: BlasLong) callconv(cc) i32 {
            tt(.init(m, n, k, a, lda, alpha, b, ldb, beta, false, c, ldc));
            return 0;
        }
        pub fn gemm_small_kernel_b0_nn(m: BlasLong, n: BlasLong, k: BlasLong, a: [*]T, lda: BlasLong, alpha: T, b: [*]T, ldb: BlasLong, c: [*]T, ldc: BlasLong) callconv(cc) i32 {
            nn(.init(m, n, k, a, lda, alpha, b, ldb, 0, true, c, ldc));
            return 0;
        }
        pub fn gemm_small_kernel_b0_nt(m: BlasLong, n: BlasLong, k: BlasLong, a: [*]T, lda: BlasLong, alpha: T, b: [*]T, ldb: BlasLong, c: [*]T, ldc: BlasLong) callconv(cc) i32 {
            nt(.init(m, n, k, a, lda, alpha, b, ldb, 0, true, c, ldc));
            return 0;
        }
        pub fn gemm_small_kernel_b0_tn(m: BlasLong, n: BlasLong, k: BlasLong, a: [*]T, lda: BlasLong, alpha: T, b: [*]T, ldb: BlasLong, c: [*]T, ldc: BlasLong) callconv(cc) i32 {
            const o: Ops = .init(m, n, k, a, lda, alpha, b, ldb, 0, true, c, ldc);
            dotRows(true, o, 0);
            return 0;
        }
        pub fn gemm_small_kernel_b0_tt(m: BlasLong, n: BlasLong, k: BlasLong, a: [*]T, lda: BlasLong, alpha: T, b: [*]T, ldb: BlasLong, c: [*]T, ldc: BlasLong) callconv(cc) i32 {
            tt(.init(m, n, k, a, lda, alpha, b, ldb, 0, true, c, ldc));
            return 0;
        }

        /// The operands, unpacked; `b0`: `C` is written, not read.
        const Ops = struct {
            m: usize,
            n: usize,
            k: usize,
            a: [*]const T,
            lda: usize,
            b: [*]const T,
            ldb: usize,
            alpha: T,
            beta: T,
            b0: bool,
            c: [*]T,
            ldc: usize,

            fn init(m: BlasLong, n: BlasLong, k: BlasLong, a: [*]T, lda: BlasLong, alpha: T, b: [*]T, ldb: BlasLong, beta: T, b0: bool, c: [*]T, ldc: BlasLong) Ops {
                return .{ .m = @intCast(m), .n = @intCast(n), .k = @intCast(k), .a = a, .lda = @intCast(lda), .b = b, .ldb = @intCast(ldb), .alpha = alpha, .beta = beta, .b0 = b0, .c = c, .ldc = @intCast(ldc) };
            }
        };

        // ---- the loops, as each file nests them ----------------------------------

        /// `nn`: rows in blocks of 4, 2 and 1 vectors, then the rest by a
        /// masked vector, or by dot blocks for at most half a vector of rows
        /// over a depth of two vectors or more.
        fn nn(o: Ops) void {
            var i: usize = 0;
            while (i + 4 * L <= o.m) : (i += 4 * L) across(false, 4, &.{ 4, 2, 1 }, o, i, 4 * L);
            while (i + 2 * L <= o.m) : (i += 2 * L) across(false, 2, &.{ 6, 2, 1 }, o, i, 2 * L);
            while (i + L <= o.m) : (i += L) across(false, 1, &.{ 6, 2, 1 }, o, i, L);
            const mm = o.m - i;
            if (mm == 0) return;
            if (mm > L / 2 or o.k < 2 * L) across(false, 1, &.{ 6, 2, 1 }, o, i, mm) else dotRows(false, o, i);
        }

        /// `nt`: rows in blocks of 4, 2 and 1 vectors, the rest by a masked
        /// vector from three quarters of one, else by broadcast rows.
        fn nt(o: Ops) void {
            var i: usize = 0;
            while (i + 4 * L <= o.m) : (i += 4 * L) across(true, 4, &.{ 6, 2, 1 }, o, i, 4 * L);
            while (i + 2 * L <= o.m) : (i += 2 * L) across(true, 2, &.{ 8, 4, 2, 1 }, o, i, 2 * L);
            while (i + L <= o.m) : (i += L) across(true, 1, &.{ 8, 4, 2, 1 }, o, i, L);
            const mm = o.m - i;
            if (mm >= 3 * L / 4) {
                across(true, 1, &.{ 8, 4, 2, 1 }, o, i, mm);
            } else if (mm > 0) {
                inline for (.{ 4, 2, 1 }) |r| {
                    while (i + r <= o.m) : (i += r) along(false, r, &.{ 4, 2 }, o, i);
                }
            }
        }

        /// `tt`: broadcast rows in blocks of 8, 4, 2 and 1.
        fn tt(o: Ops) void {
            var i: usize = 0;
            while (i + 8 <= o.m) : (i += 8) along(true, 8, &.{2}, o, i);
            inline for (.{ 4, 2, 1 }) |r| {
                while (i + r <= o.m) : (i += r) along(true, r, &.{ 4, 2 }, o, i);
            }
        }

        /// Dot blocks from row `from`: rows by 4, 2 and 1, columns by 4, 2 and 1.
        fn dotRows(comptime ta: bool, o: Ops, from: usize) void {
            var i = from;
            inline for (.{ 4, 2, 1 }) |r| {
                while (i + r <= o.m) : (i += r) {
                    var j: usize = 0;
                    inline for (.{ 4, 2, 1 }) |w| {
                        while (j + w <= o.n) : (j += w) dot(ta, r, w, o, i, j);
                    }
                }
            }
        }

        // ---- the blocks ----------------------------------------------------------

        /// Rows `[i, i + rows)` (`nv` vectors, or one masked) of `op(A)`
        /// down columns, times broadcast `op(B)`, over all columns by
        /// `widths` (`MATMUL_512`, `STORE_512`, `MASK_STORE_512`).
        fn across(comptime tb: bool, comptime nv: usize, comptime widths: []const usize, o: Ops, i: usize, rows: usize) void {
            var j: usize = 0;
            inline for (widths) |w| {
                while (j + w <= o.n) : (j += w) {
                    if (nv == 1 and rows < L) acrossBlock(tb, 1, w, true, o, i, j, rows) else acrossBlock(tb, nv, w, false, o, i, j, rows);
                }
            }
        }

        fn acrossBlock(comptime tb: bool, comptime nv: usize, comptime wn: usize, comptime masked: bool, o: Ops, i: usize, j: usize, rows: usize) void {
            var acc: [wn][nv]V = @splat(@splat(@splat(0)));
            for (0..o.k) |l| {
                var av: [nv]V = undefined;
                inline for (0..nv) |v| av[v] = load(o.a + o.lda * l + i + v * L, if (masked) rows else L);
                inline for (0..wn) |q| {
                    const y = if (tb) o.b[o.ldb * l + j + q] else o.b[l + o.ldb * (j + q)];
                    inline for (0..nv) |v| acc[q][v] = madd(av[v], @splat(y), acc[q][v]);
                }
            }
            inline for (0..wn) |q| inline for (0..nv) |v| {
                const p = o.c + (j + q) * o.ldc + i + v * L;
                const len = if (masked) rows else L;
                var r = acc[q][v] * @as(V, @splat(o.alpha));
                if (!o.b0) r = madd(load(p, len), @splat(o.beta), r);
                store(p, r, len);
            };
        }

        /// `rows` broadcast rows of `op(A)` from `i` times `op(B)` along its
        /// rows: column blocks of each of `nvs` vectors, then by one masked
        /// vector (`BROADCAST_LOAD_A_512`, `SCATTER_STORE_512`, `tt`'s
        /// `REORDER_STORE`s, which store the same elements).
        fn along(comptime ta: bool, comptime rows: usize, comptime nvs: []const usize, o: Ops, i: usize) void {
            var j: usize = 0;
            inline for (nvs) |nv| {
                while (j + nv * L <= o.n) : (j += nv * L) alongBlock(ta, rows, nv, false, o, i, j, L);
            }
            while (j < o.n) : (j += L) {
                const cols = @min(L, o.n - j);
                if (cols == L) alongBlock(ta, rows, 1, false, o, i, j, L) else alongBlock(ta, rows, 1, true, o, i, j, cols);
            }
        }

        fn alongBlock(comptime ta: bool, comptime rows: usize, comptime nv: usize, comptime masked: bool, o: Ops, i: usize, j: usize, cols: usize) void {
            var acc: [rows][nv]V = @splat(@splat(@splat(0)));
            for (0..o.k) |l| {
                var bv: [nv]V = undefined;
                inline for (0..nv) |v| bv[v] = load(o.b + o.ldb * l + j + v * L, if (masked) cols else L);
                inline for (0..rows) |r| {
                    const x = if (ta) o.a[l + o.lda * (i + r)] else o.a[o.lda * l + i + r];
                    inline for (0..nv) |v| acc[r][v] = madd(@splat(x), bv[v], acc[r][v]);
                }
            }
            inline for (0..rows) |r| inline for (0..nv) |v| {
                const p = o.c + (j + v * L) * o.ldc + i + r;
                const len = if (masked) cols else L;
                var s = acc[r][v] * @as(V, @splat(o.alpha));
                if (!o.b0) s = madd(gather(L, p, o.ldc, len), @splat(o.beta), s);
                scatter(L, p, o.ldc, s, len);
            };
        }

        /// `rows x cols` elements from `(i, j)`, each lane summing every
        /// `L`-th product of the depth, the last part masked to zeros
        /// (`LOAD_KA_512`, `MASK_LOAD_KA_512`); then four rows (or four
        /// columns of fewer rows) reduce together (`REDUCE_4`), the rest one
        /// at a time (`STORE_REDUCE`).
        fn dot(comptime ta: bool, comptime rows: usize, comptime cols: usize, o: Ops, i: usize, j: usize) void {
            var acc: [rows][cols]V = @splat(@splat(@splat(0)));
            var l: usize = 0;
            while (l < o.k) : (l += L) {
                const len = @min(L, o.k - l);
                var av: [rows]V = undefined;
                inline for (0..rows) |r| av[r] = if (ta) load(o.a + (i + r) * o.lda + l, len) else gather(L, o.a + i + r + o.lda * l, o.lda, len);
                inline for (0..cols) |q| {
                    const bv = load(o.b + (j + q) * o.ldb + l, len);
                    inline for (0..rows) |r| acc[r][q] = madd(av[r], bv, acc[r][q]);
                }
            }
            if (rows == 4) {
                // `REDUCE_M4`: column `q`'s four rows, stored down it
                inline for (0..cols) |q| {
                    var s: @Vector(4, T) = undefined;
                    inline for (0..4) |r| s[r] = pairs(L, acc[r][q]);
                    four(o, s, o.c + (j + q) * o.ldc + i, 1);
                }
            } else if (cols == 4) {
                // `REDUCE_N4`: row `r`'s four columns, stored across them
                inline for (0..rows) |r| {
                    var s: @Vector(4, T) = undefined;
                    inline for (0..4) |q| s[q] = pairs(L, acc[r][q]);
                    four(o, s, o.c + j * o.ldc + i + r, o.ldc);
                }
            } else inline for (0..rows) |r| inline for (0..cols) |q| {
                const p = &o.c[(j + q) * o.ldc + i + r];
                const sum = halves(L, acc[r][q]);
                // `alpha * sum + beta * C`, contracted as the expression rule says
                p.* = if (o.b0) o.alpha * sum else fma.madd(fma.fuse, T, o.alpha, sum, o.beta * p.*);
            };
        }

        /// Four reduced sums to `C` from `p`, `step` apart: `alpha*s`, then
        /// fused with `beta*C` (`STORE_REDUCE_M4`, `_N4`).
        fn four(o: Ops, sums: @Vector(4, T), p: [*]T, step: usize) void {
            var s = sums * @as(@Vector(4, T), @splat(o.alpha));
            if (!o.b0) s = fma.madd(true, @Vector(4, T), gather(4, p, step, 4), @splat(o.beta), s);
            scatter(4, p, step, s, 4);
        }

        // ---- lanes ---------------------------------------------------------------

        inline fn madd(a: V, b: V, c: V) V {
            return fma.madd(true, V, a, b, c);
        }

        /// `len` elements from `p`, zeros past them (`maskz_loadu`).
        inline fn load(p: [*]const T, len: usize) V {
            if (len == L) return p[0..L].*;
            var v: [L]T = @splat(0);
            @memcpy(v[0..len], p[0..len]);
            return v;
        }

        inline fn store(p: [*]T, v: V, len: usize) void {
            if (len == L) {
                p[0..L].* = v;
                return;
            }
            const w: [L]T = v;
            @memcpy(p[0..len], w[0..len]);
        }

        /// `len` elements `step` apart from `p`, zeros past them.
        inline fn gather(comptime n: usize, p: [*]const T, step: usize, len: usize) @Vector(n, T) {
            var v: [n]T = @splat(0);
            for (v[0..len], 0..) |*e, x| e.* = p[x * step];
            return v;
        }

        inline fn scatter(comptime n: usize, p: [*]T, step: usize, v: @Vector(n, T), len: usize) void {
            const w: [n]T = v;
            for (w[0..len], 0..) |e, x| p[x * step] = e;
        }

        /// `REDUCE_4`'s sum of one vector: neighbors in pairs, then pairs
        /// of those, down to one.
        fn pairs(comptime n: usize, x: @Vector(n, T)) T {
            if (n == 1) return x[0] else {
                const h = n / 2;
                const even = comptime std.simd.iota(i32, h) * @as(@Vector(h, i32), @splat(2));
                const odd = comptime even + @as(@Vector(h, i32), @splat(1));
                return pairs(h, @shuffle(T, x, undefined, even) + @shuffle(T, x, undefined, odd));
            }
        }

        /// `_mm512_reduce_add_p[sd]`'s sum: the upper half onto the lower,
        /// down to one (gcc's `__MM512_REDUCE_OP`, LLVM's shuffle reduction).
        fn halves(comptime n: usize, x: @Vector(n, T)) T {
            if (n == 1) return x[0] else {
                const h = n / 2;
                const low = comptime std.simd.iota(i32, h);
                const high = comptime low + @as(@Vector(h, i32), @splat(h));
                return halves(h, @shuffle(T, x, undefined, low) + @shuffle(T, x, undefined, high));
            }
        }
    };
}

// ---- Cooper Lake's sbgemm --------------------------------------------------------

/// `kernel/x86_64/sbgemm_small_kernel_{nn,nt,tn,tt}_cooperlake.c`, each also
/// built `B0` (`sbgemm_small_kernel_template_cooperlake.c`): `C` scaled by
/// `beta`, or zeroed, then the depth in slices of `slice`, each adding
/// `alpha` times its sums to `C`, over blocks of 32 rows by 8 columns
/// (`sbgemm_blocking_kernel_*`). Upstream sums with `vdpbf16ps`, pairs of
/// products rounded as the instruction does, without denormals; here each
/// slice sums in `f32` in depth order, unfused, as the portable half gemm
/// (`ext/bf16/kernel.zig`), so only the slicing and the shape are upstream's.
/// Sapphire Rapids empties these variables (`KERNEL.SAPPHIRERAPIDS`), so it
/// takes the generic permit, which allows nothing.
pub const Cooperlake = struct {
    /// `BF16_BLOCK_THRES_K`, `_M` and `BF16_BLOCK_STEP_N`
    /// (`sbgemm_microk_cooperlake_template.c`).
    pub const slice = 1024;
    const block_m = 32;
    const step_n = 8;
    const L = 16;
    const V = @Vector(L, f32);

    pub fn entry(comptime slot: Slot) util.SlotFn(slot) {
        const name = @tagName(slot);
        const op = name[name.len - 2 ..];
        const K = Ops(op[0] == 't', op[1] == 't');
        return if (comptime std.mem.indexOf(u8, name, "_b0_") != null) &K.zero else &K.scaled;
    }

    fn Ops(comptime ta: bool, comptime tb: bool) type {
        return struct {
            fn scaled(m: BlasLong, n: BlasLong, k: BlasLong, a: [*]BF16, lda: BlasLong, alpha: f32, b: [*]BF16, ldb: BlasLong, beta: f32, c: [*]f32, ldc: BlasLong) callconv(.c) i32 {
                run(ta, tb, false, @intCast(m), @intCast(n), @intCast(k), a, @intCast(lda), alpha, b, @intCast(ldb), beta, c, @intCast(ldc));
                return 0;
            }
            fn zero(m: BlasLong, n: BlasLong, k: BlasLong, a: [*]BF16, lda: BlasLong, alpha: f32, b: [*]BF16, ldb: BlasLong, c: [*]f32, ldc: BlasLong) callconv(.c) i32 {
                run(ta, tb, true, @intCast(m), @intCast(n), @intCast(k), a, @intCast(lda), alpha, b, @intCast(ldb), 0, c, @intCast(ldc));
                return 0;
            }
        };
    }

    fn run(comptime ta: bool, comptime tb: bool, b0: bool, m: usize, n: usize, k: usize, a: [*]const BF16, lda: usize, alpha: f32, b: [*]const BF16, ldb: usize, beta: f32, c: [*]f32, ldc: usize) void {
        // `sbgemm_zero_operation`, `sbgemm_scal_operation`
        for (0..n) |j| for (c[j * ldc ..][0..m]) |*e| {
            e.* = if (b0) 0 else beta * e.*;
        };
        var l0: usize = 0;
        while (l0 < k) : (l0 += slice) {
            const kl = @min(slice, k - l0);
            var j: usize = 0;
            while (j + step_n <= n) : (j += step_n) rows(ta, tb, step_n, m, l0, kl, a, lda, alpha, b, ldb, c, ldc, j);
            while (j < n) : (j += 1) rows(ta, tb, 1, m, l0, kl, a, lda, alpha, b, ldb, c, ldc, j);
        }
    }

    /// Rows by `block_m`, then by a vector, then the rest masked.
    fn rows(comptime ta: bool, comptime tb: bool, comptime wn: usize, m: usize, l0: usize, kl: usize, a: [*]const BF16, lda: usize, alpha: f32, b: [*]const BF16, ldb: usize, c: [*]f32, ldc: usize, j: usize) void {
        var i: usize = 0;
        while (i + block_m <= m) : (i += block_m) block(ta, tb, block_m / L, wn, L, i, l0, kl, a, lda, alpha, b, ldb, c, ldc, j);
        while (i + L <= m) : (i += L) block(ta, tb, 1, wn, L, i, l0, kl, a, lda, alpha, b, ldb, c, ldc, j);
        if (i < m) block(ta, tb, 1, wn, m - i, i, l0, kl, a, lda, alpha, b, ldb, c, ldc, j);
    }

    /// `nv` vectors of rows (`len` of the last) by `wn` columns over one
    /// slice: `C += alpha*sum`, the product and the sum each rounded.
    fn block(comptime ta: bool, comptime tb: bool, comptime nv: usize, comptime wn: usize, len: usize, i: usize, l0: usize, kl: usize, a: [*]const BF16, lda: usize, alpha: f32, b: [*]const BF16, ldb: usize, c: [*]f32, ldc: usize, j: usize) void {
        var acc: [wn][nv]V = @splat(@splat(@splat(0)));
        for (l0..l0 + kl) |l| {
            var av: [nv]V = undefined;
            inline for (0..nv) |v| {
                const r0 = i + v * L;
                if (!ta and len == L) {
                    av[v] = bf16_kernel.widen(BF16, L, a + r0 + l * lda);
                } else {
                    var w: [L]f32 = @splat(0);
                    for (w[0..len], 0..) |*e, r| e.* = bf16_kernel.wide(BF16, if (ta) a[l + (r0 + r) * lda] else a[r0 + r + l * lda]);
                    av[v] = w;
                }
            }
            inline for (0..wn) |q| {
                const y = bf16_kernel.wide(BF16, if (tb) b[j + q + l * ldb] else b[l + (j + q) * ldb]);
                inline for (0..nv) |v| acc[q][v] = fma.madd(false, V, av[v], @splat(y), acc[q][v]);
            }
        }
        inline for (0..wn) |q| inline for (0..nv) |v| {
            const p = c + i + v * L + (j + q) * ldc;
            store(p, fma.madd(false, V, acc[q][v], @splat(alpha), load(p, len)), len);
        };
    }

    fn load(p: [*]const f32, len: usize) V {
        if (len == L) return p[0..L].*;
        var w: [L]f32 = @splat(0);
        @memcpy(w[0..len], p[0..len]);
        return w;
    }

    fn store(p: [*]f32, v: V, len: usize) void {
        const w: [L]f32 = v;
        @memcpy(p[0..len], w[0..len]);
    }
};

/// The ports and where they serve (`KERNEL.COOPERLAKE`'s `SBGEMM_SMALL_K_*`).
pub const ports = [_]registry.Port(Core){
    .{ .module = Cooperlake, .cores = &.{.cooperlake}, .slots = &.{ .sbgemm_small_kernel_nn, .sbgemm_small_kernel_nt, .sbgemm_small_kernel_tn, .sbgemm_small_kernel_tt, .sbgemm_small_kernel_b0_nn, .sbgemm_small_kernel_b0_nt, .sbgemm_small_kernel_b0_tn, .sbgemm_small_kernel_b0_tt }, .written = .vector },
};

const P = registry.Ports(Core, "x86_64", &ports);

pub const has = P.has;
pub const hasSlot = P.hasSlot;
pub const exportCore = P.exportCore;
pub const install = P.install;
