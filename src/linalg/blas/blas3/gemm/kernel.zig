//! The Zig gemm kernels for real `T` and its `Complex`: the packing copies,
//! `beta`, and the register-tile kernel `C += alpha * A * B` on packed panels
//! (the slots `<p>gemm_{incopy,itcopy,oncopy,otcopy,beta,kernel}`, complex
//! `kernel_{n,l,r,b}`). Compiled into the main library as the portable
//! kernels (`installGeneric`), and into every per-core object (`Registry`),
//! where the core's features set the vector width and register count, so one
//! source serves every core. The other level-3 routines' packing and kernels
//! (`sets`) go with them: they share the tile.
//!
//! Packed layout, shared by the copies and the kernel: `n / u` panels of
//! width `u`, then one of each power of two below `u` set in `n % u`. Panel
//! `p` at column `c0` with width `w` sits at `m * c0`, as `m` groups of `w`.
//! For a power-of-two `u` this is OpenBLAS's `kernel/generic` layout. A
//! complex element is its two parts, real first.
//!
//! The tile is `mr x nr`: `mv` vectors of `vl` lanes by `nr` columns, the
//! most accumulators the registers hold beside the `mv` loads of A and a
//! broadcast of B (`geometry`). Complex lanes hold whole elements, and each
//! has two accumulators, one per part of B (`kernel/x86_64/zgemm_kernel_4x2_haswell.S`'s
//! scheme); without vectors the parts are separate scalars. Imports stay
//! inside `linalg/blas/`: the per-core objects compile this file on their own.

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("../../types.zig");
const fma = @import("../../fma.zig");
const util = @import("../../util.zig");
const Vtable = @import("../../Vtable.zig");

const BlasLong = types.BlasLong;
const XDouble = types.XDouble;
const isComplex = types.isComplex;

/// A kernel set's tile: `mv` vectors of `vl` lanes (`mr` rows) by `nr`
/// columns; `mr` and `nr` are also the copies' unrolls. For complex, `vl`
/// counts elements, and `split` holds their parts in scalars.
pub const Geometry = struct {
    vl: u16,
    mv: u16,
    nr: u16,
    split: bool = false,

    pub fn mr(g: Geometry) u16 {
        return g.vl * g.mv;
    }
};

/// Bits of the widest vector of `T` the target computes on natively; 0 for
/// none. Not `std.simd.suggestVectorLength`: that counts SVE as 256-bit
/// fixed vectors (LLVM splits them into NEON pairs) and leaves AVX at 128.
fn vectorBits(comptime T: type, comptime cpu: std.Target.Cpu) u16 {
    const double = T == f64;
    if (T != f32 and T != f64) return 0;
    return switch (cpu.arch) {
        .x86_64, .x86 => if (cpu.has(.x86, .avx512f) and sysinfo.zig_backend != .stage2_x86_64)
            512
        else if (cpu.has(.x86, .avx))
            256
        else if (cpu.has(.x86, if (double) .sse2 else .sse))
            128
        else
            0,
        .aarch64, .aarch64_be => if (cpu.has(.aarch64, .neon)) 128 else 0,
        // ARMv7 NEON has no double lanes
        .arm, .armeb, .thumb, .thumbeb => if (!double and cpu.has(.arm, .neon)) 128 else 0,
        // AltiVec has no double lanes; VSX has
        .powerpc, .powerpcle, .powerpc64, .powerpc64le => if (cpu.has(.powerpc, .vsx) or (!double and cpu.has(.powerpc, .altivec))) 128 else 0,
        .riscv32, .riscv32be, .riscv64, .riscv64be => if (!cpu.has(.riscv, .v))
            0
        else if (cpu.has(.riscv, .zvl512b))
            512
        else if (cpu.has(.riscv, .zvl256b))
            256
        else
            128,
        .loongarch32, .loongarch64 => if (cpu.has(.loongarch, .lasx)) 256 else if (cpu.has(.loongarch, .lsx)) 128 else 0,
        .mips, .mipsel, .mips64, .mips64el => if (cpu.has(.mips, .msa)) 128 else 0,
        // z13 vectors hold doubles; singles need the z14 enhancements
        .s390x => if (cpu.has(.s390x, .vector) and (double or cpu.has(.s390x, .vector_enhancements_1))) 128 else 0,
        .wasm32, .wasm64 => if (cpu.has(.wasm, .simd128)) 128 else 0,
        else => 0,
    };
}

/// Registers that hold one such vector (or one scalar where `bits` is 0).
fn registers(comptime bits: u16, comptime cpu: std.Target.Cpu) u16 {
    return switch (cpu.arch) {
        .x86_64 => if (bits == 512) 32 else 16,
        .x86 => 8,
        .arm, .armeb, .thumb, .thumbeb => if (bits == 0 and cpu.has(.arm, .d32)) 32 else 16,
        .powerpc, .powerpcle, .powerpc64, .powerpc64le => if (cpu.has(.powerpc, .vsx)) 64 else 32,
        .s390x => if (bits == 0) 16 else 32,
        .wasm32, .wasm64 => 16,
        .mips, .mipsel, .sparc => 16,
        else => 32,
    };
}

/// Whether LLVM multiplies by a lane of a register (AArch64 `fmla .s[i]`,
/// NEON `vmla` by scalar): B then comes as whole vectors, `nr / vl` registers,
/// rather than one broadcast at a time.
fn byLane(comptime cpu: std.Target.Cpu, comptime bits: u16) bool {
    return bits > 0 and switch (cpu.arch) {
        .aarch64, .aarch64_be, .arm, .armeb, .thumb, .thumbeb => true,
        else => false,
    };
}

/// The tile for `T` on `cpu`: the most accumulators (`mv * nr`) that fit
/// beside `mv` A vectors, B (one broadcast, or `nr / vl` vectors where
/// `byLane`) and one spare register, with no more loads per step than three
/// quarters of the multiply-adds (relaxed where nothing meets that); ties
/// go to fewer loads, then to wider `nr`. Without the spare, LLVM spills
/// accumulators in the loop (`notes.md`). Complex counts two accumulators
/// per element and two broadcasts of B, and twice the registers where split.
pub fn geometryFor(comptime T: type, comptime cpu: std.Target.Cpu) Geometry {
    const complex = isComplex(T);
    const F = types.Float(T);
    const raw_bits = vectorBits(F, cpu);
    // a complex lane is two floats
    const split = complex and raw_bits < 2 * @bitSizeOf(F);
    const bits: u16 = if (split) 0 else raw_bits;
    const vl: u16 = if (bits == 0) 1 else bits / @bitSizeOf(F) / @as(u16, if (complex) 2 else 1);
    const regs = registers(bits, cpu);
    @setEvalBranchQuota(100_000);
    inline for (.{ true, false }) |strict| {
        var best: ?Geometry = null;
        var best_loads: usize = 0;
        for (1..5) |mv| for (1..17) |nr| {
            const g: Geometry = .{ .vl = vl, .mv = mv, .nr = nr, .split = split };
            // accumulators, registers, loads and multiply-adds per step
            const acc = if (complex) 2 * mv * nr else mv * nr;
            const a_regs = if (split) 2 * mv else mv;
            const b_regs = if (complex) 2 else if (byLane(cpu, bits)) (nr + vl - 1) / vl else 1;
            const loads = a_regs + if (complex) 2 * nr else nr;
            const madds = if (split) 4 * mv * nr else acc;
            if (acc + a_regs + b_regs + 1 > regs) continue;
            if (strict and 4 * loads > 3 * madds) continue;
            if (best) |b| {
                const bacc = @as(usize, b.mv) * b.nr;
                if (mv * nr < bacc) continue;
                if (mv * nr == bacc and loads > best_loads) continue;
                if (mv * nr == bacc and loads == best_loads and nr < b.nr) continue;
            }
            best = g;
            best_loads = loads;
        };
        if (best) |b| return b;
    }
    unreachable;
}

/// The compile target's tile.
pub fn geometry(comptime T: type) Geometry {
    @setEvalBranchQuota(100_000);
    return comptime geometryFor(T, sysinfo.cpu);
}

/// The powers of two below `u`, largest first: the remainder panels' widths.
pub fn tails(comptime u: usize) []const usize {
    comptime {
        var out: []const usize = &.{};
        if (u < 2) return out;
        var w = std.math.floorPowerOfTwo(usize, u - 1);
        while (w >= 1) : (w /= 2) out = out ++ [_]usize{w};
        return out;
    }
}

/// The width of the panel at the start of `rest` lines packed in panels of
/// `u`: `u`, or past the whole panels the highest bit of what is left.
pub fn panelWidth(rest: usize, u: usize) usize {
    return if (rest >= u) u else std.math.floorPowerOfTwo(usize, rest);
}

fn Lanes(comptime T: type, comptime n: usize) type {
    return if (n == 1) T else @Vector(n, T);
}

inline fn load(comptime T: type, comptime n: usize, p: [*]const T) Lanes(T, n) {
    return if (n == 1) p[0] else util.load(@Vector(n, T), p);
}

inline fn store(comptime T: type, comptime n: usize, p: [*]T, v: Lanes(T, n)) void {
    if (n == 1) p[0] = v else util.store(@Vector(n, T), p, v);
}

inline fn splat(comptime T: type, comptime n: usize, x: T) Lanes(T, n) {
    return if (n == 1) x else @splat(x);
}

/// `a*b + c` under the build's contraction policy (`fma.zig`).
pub inline fn madd(comptime T: type, comptime n: usize, a: Lanes(T, n), b: Lanes(T, n), c: Lanes(T, n)) Lanes(T, n) {
    const fused = (if (n == 1) fma.fuse else fma.fuse_vector) and fma.fusible(T);
    return if (fused) @mulAdd(Lanes(T, n), a, b, c) else a * b + c;
}

/// Packs `n` lines of depth `m` of elements `E` in panels of `u`. `across`:
/// line `j` is element `j` of each of `m` rows `lda` apart (`tcopy`); else it
/// is column `j`, `lda` apart (`ncopy`).
fn pack(comptime E: type, comptime u: usize, comptime across: bool, m: BlasLong, n: BlasLong, a: [*]const E, lda: BlasLong, b: [*]E) void {
    const depth: usize = @intCast(m);
    const lines: usize = @intCast(n);
    const ld: usize = @intCast(lda);
    var c0: usize = 0;
    while (c0 + u <= lines) : (c0 += u) panel(E, u, across, depth, a, ld, c0, b + depth * c0);
    inline for (comptime tails(u)) |w| {
        if ((lines - c0) & w != 0) {
            panel(E, w, across, depth, a, ld, c0, b + depth * c0);
            c0 += w;
        }
    }
}

inline fn panel(comptime E: type, comptime w: usize, comptime across: bool, depth: usize, a: [*]const E, ld: usize, c0: usize, dst: [*]E) void {
    if (across) {
        for (0..depth) |l| dst[l * w ..][0..w].* = a[l * ld + c0 ..][0..w].*;
    } else {
        var cols: [w][*]const E = undefined;
        inline for (0..w) |c| cols[c] = a + (c0 + c) * ld;
        for (0..depth) |l| {
            inline for (0..w) |c| dst[l * w + c] = cols[c][l];
        }
    }
}

/// Packs `lines` lines of `depth` as `pack` does, element `(j, l)` being
/// `src.at(j, l)`; where that is an optional, null leaves it unwritten. For
/// the other routines' copies, which transform what they pack.
pub inline fn packWith(comptime E: type, comptime u: usize, depth: usize, lines: usize, dst: [*]E, src: anytype) void {
    var c0: usize = 0;
    while (c0 + u <= lines) : (c0 += u) panelWith(E, u, depth, c0, dst + depth * c0, src);
    inline for (comptime tails(u)) |w| {
        if ((lines - c0) & w != 0) {
            panelWith(E, w, depth, c0, dst + depth * c0, src);
            c0 += w;
        }
    }
}

inline fn panelWith(comptime E: type, comptime w: usize, depth: usize, c0: usize, dst: [*]E, src: anytype) void {
    for (0..depth) |l| {
        inline for (0..w) |c| {
            const v = src.at(c0 + c, l);
            if (@TypeOf(v) == ?E) {
                if (v) |x| dst[l * w + c] = x;
            } else dst[l * w + c] = v;
        }
    }
}

/// Which operand a complex kernel conjugates: `n` neither, `l` A, `r` B,
/// `b` both (`GEMM_KERNEL_N`, `_L`, `_R`, `_B`).
pub const Conj = enum {
    n,
    l,
    r,
    b,

    pub fn of(a: bool, b_: bool) Conj {
        return if (a) (if (b_) .b else .l) else (if (b_) .r else .n);
    }
    pub fn ofA(c: Conj) bool {
        return c == .l or c == .b;
    }
    pub fn ofB(c: Conj) bool {
        return c == .r or c == .b;
    }
};

/// The kernel set of `T` with tile `geo`; its functions have the vtable's
/// signatures, named as the slots less the precision letter.
pub fn Kernels(comptime T: type, comptime geo: Geometry) type {
    return if (isComplex(T)) ComplexKernels(T, geo) else RealKernels(T, geo);
}

fn RealKernels(comptime T: type, comptime geo: Geometry) type {
    const vl = geo.vl;
    // C for the vtable's types, Zig's own for the rest
    const cc = types.kernelCc(T);
    // `f16` sums in `f32`; vectors are only ever `f32` or `f64`, summed as is
    const A = types.Acc(T);
    if (A != T and vl > 1) @compileError("a widened sum needs scalar lanes");
    const mr = geo.mr();
    const nr = geo.nr;
    return struct {
        pub const unroll_m: i32 = mr;
        pub const unroll_n: i32 = nr;

        /// `C := beta*C` on `m x n`; a zero `beta` writes zeros without
        /// reading `C`.
        pub fn gemm_beta(m: BlasLong, n: BlasLong, _: BlasLong, b: T, _: [*]T, _: BlasLong, _: [*]T, _: BlasLong, c: [*]T, ldc: BlasLong) callconv(cc) i32 {
            const height: usize = @intCast(m);
            const ld: usize = @intCast(ldc);
            for (0..@intCast(n)) |j| {
                const col = c[j * ld ..][0..height];
                if (b == 0) {
                    @memset(col, 0);
                } else {
                    for (col) |*v| v.* *= b;
                }
            }
            return 0;
        }

        /// Rows of `op(A) = A^T` (`m` deep, `n` rows), `mr` a panel.
        pub fn gemm_incopy(m: BlasLong, n: BlasLong, a: [*]T, lda: BlasLong, b: [*]T) callconv(cc) i32 {
            pack(T, mr, false, m, n, a, lda, b);
            return 0;
        }
        /// Rows of `op(A) = A`.
        pub fn gemm_itcopy(m: BlasLong, n: BlasLong, a: [*]T, lda: BlasLong, b: [*]T) callconv(cc) i32 {
            pack(T, mr, true, m, n, a, lda, b);
            return 0;
        }
        /// Columns of `op(B) = B`, `nr` a panel.
        pub fn gemm_oncopy(m: BlasLong, n: BlasLong, a: [*]T, lda: BlasLong, b: [*]T) callconv(cc) i32 {
            pack(T, nr, false, m, n, a, lda, b);
            return 0;
        }
        /// Columns of `op(B) = B^T`.
        pub fn gemm_otcopy(m: BlasLong, n: BlasLong, a: [*]T, lda: BlasLong, b: [*]T) callconv(cc) i32 {
            pack(T, nr, true, m, n, a, lda, b);
            return 0;
        }

        /// `C += alpha * A * B` for `A` packed `bm x bk` (rows in panels of
        /// `mr`) and `B` packed `bk x bn` (columns in panels of `nr`). Columns
        /// outermost, so a B panel stays in L1 while the A panels stream from L2.
        pub fn gemm_kernel(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha: T, ba: [*]T, bb: [*]T, c: [*]T, ldc: BlasLong) callconv(cc) i32 {
            const m: usize = @intCast(bm);
            const n: usize = @intCast(bn);
            const k: usize = @intCast(bk);
            const ld: usize = @intCast(ldc);
            var q0: usize = 0;
            while (q0 + nr <= n) : (q0 += nr) rows(nr, m, k, alpha, ba, bb + q0 * k, c + q0 * ld, ld);
            inline for (comptime tails(nr)) |w| {
                if ((n - q0) & w != 0) {
                    rows(w, m, k, alpha, ba, bb + q0 * k, c + q0 * ld, ld);
                    q0 += w;
                }
            }
            return 0;
        }

        /// Every row panel against one B panel of width `wn`: the full tile
        /// inline, the row tails of whole vectors as calls of their own, and
        /// the narrower ones one row at a time (`thin`), so the unrolled tiles
        /// stay few.
        fn rows(comptime wn: usize, m: usize, k: usize, alpha: T, ba: [*]const T, pb: [*]const T, c: [*]T, ld: usize) void {
            var r0: usize = 0;
            while (r0 + mr <= m) : (r0 += mr) @call(.always_inline, tile, .{ mr, wn, k, alpha, ba + r0 * k, pb, c + r0, ld });
            inline for (comptime tails(mr)) |w| {
                if ((m - r0) & w != 0) {
                    if (w >= vl) {
                        @call(.never_inline, tile, .{ w, wn, k, alpha, ba + r0 * k, pb, c + r0, ld });
                    } else {
                        @call(.never_inline, thin, .{ wn, w, k, alpha, ba + r0 * k, pb, c + r0, ld });
                    }
                    r0 += w;
                }
            }
        }

        /// A row panel narrower than a vector, `w` rows: each row against
        /// the `wn` columns as one vector, in two chains over `k`; one
        /// function for every row tail of a column width.
        fn thin(comptime wn: usize, w: usize, k: usize, alpha: T, pa: [*]const T, pb: [*]const T, c: [*]T, ld: usize) void {
            const V = @Vector(wn, T);
            const fused = fma.fuse_vector and fma.fusible(T);
            for (0..w) |r| {
                var acc: [2]V = .{ @splat(0), @splat(0) };
                var l: usize = 0;
                while (l + 2 <= k) : (l += 2) {
                    inline for (0..2) |h| {
                        const x: V = @splat(pa[(l + h) * w + r]);
                        const y = util.load(V, pb + (l + h) * wn);
                        acc[h] = if (fused) @mulAdd(V, x, y, acc[h]) else x * y + acc[h];
                    }
                }
                if (l < k) {
                    const x: V = @splat(pa[l * w + r]);
                    const y = util.load(V, pb + l * wn);
                    acc[0] = if (fused) @mulAdd(V, x, y, acc[0]) else x * y + acc[0];
                }
                const sum: [wn]T = acc[0] + acc[1];
                inline for (0..wn) |j| c[r + j * ld] = madd(T, 1, sum[j], alpha, c[r + j * ld]);
            }
        }

        /// One `wm x wn` tile: `wm / lanes` vectors per column in registers
        /// over all of `k`, then `C += alpha * acc`.
        fn tile(comptime wm: usize, comptime wn: usize, k: usize, alpha: T, pa: [*]const T, pb: [*]const T, c: [*]T, ld: usize) void {
            @setEvalBranchQuota(100_000);
            const lanes = @min(vl, wm);
            const nv = wm / lanes;
            const V = Lanes(A, lanes);
            var acc: [wn][nv]V = undefined;
            inline for (0..wn) |j| inline for (0..nv) |v| {
                acc[j][v] = splat(A, lanes, 0);
            };
            var a = pa;
            var b = pb;
            for (0..k) |_| {
                var av: [nv]V = undefined;
                inline for (0..nv) |v| av[v] = wide(lanes, a + v * lanes);
                inline for (0..wn) |j| {
                    const bv = splat(A, lanes, b[j]);
                    inline for (0..nv) |v| acc[j][v] = madd(A, lanes, av[v], bv, acc[j][v]);
                }
                a += wm;
                b += wn;
            }
            const al = splat(A, lanes, alpha);
            inline for (0..wn) |j| inline for (0..nv) |v| {
                const p = c + j * ld + v * lanes;
                const r = madd(A, lanes, acc[j][v], al, wide(lanes, p));
                if (A == T) store(T, lanes, p, r) else p[0] = @floatCast(r);
            };
        }

        /// `lanes` elements at `p` as `A`.
        inline fn wide(comptime lanes: usize, p: [*]const T) Lanes(A, lanes) {
            return if (A == T) load(T, lanes, p) else p[0];
        }
    };
}

fn ComplexKernels(comptime T: type, comptime geo: Geometry) type {
    const F = types.Float(T);
    const cc = types.kernelCc(F);
    const A = types.Acc(F);
    if (A != F and !geo.split) @compileError("a widened sum needs split parts");
    const vl = geo.vl;
    const mr = geo.mr();
    const nr = geo.nr;
    // one element's parts, as the copies move them
    const P = [2]F;
    return struct {
        pub const unroll_m: i32 = mr;
        pub const unroll_n: i32 = nr;

        /// `C := beta*C` on `m x n`; a zero `beta` writes zeros without
        /// reading `C` (`kernel/generic/zgemm_beta.c`).
        pub fn gemm_beta(m: BlasLong, n: BlasLong, _: BlasLong, br: F, bi: F, _: [*]F, _: BlasLong, _: [*]F, _: BlasLong, c: [*]F, ldc: BlasLong) callconv(cc) i32 {
            const height: usize = @intCast(m);
            const ld: usize = @intCast(ldc);
            for (0..@intCast(n)) |j| {
                const col = c[2 * j * ld ..][0 .. 2 * height];
                if (br == 0 and bi == 0) {
                    @memset(col, 0);
                } else {
                    for (0..height) |i| {
                        const x = col[2 * i];
                        const y = col[2 * i + 1];
                        col[2 * i] = br * x - bi * y;
                        col[2 * i + 1] = br * y + bi * x;
                    }
                }
            }
            return 0;
        }

        pub fn gemm_incopy(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, b: [*]F) callconv(cc) i32 {
            pack(P, mr, false, m, n, @ptrCast(a), lda, @ptrCast(b));
            return 0;
        }
        pub fn gemm_itcopy(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, b: [*]F) callconv(cc) i32 {
            pack(P, mr, true, m, n, @ptrCast(a), lda, @ptrCast(b));
            return 0;
        }
        pub fn gemm_oncopy(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, b: [*]F) callconv(cc) i32 {
            pack(P, nr, false, m, n, @ptrCast(a), lda, @ptrCast(b));
            return 0;
        }
        pub fn gemm_otcopy(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, b: [*]F) callconv(cc) i32 {
            pack(P, nr, true, m, n, @ptrCast(a), lda, @ptrCast(b));
            return 0;
        }

        pub fn gemm_kernel_n(bm: BlasLong, bn: BlasLong, bk: BlasLong, ar: F, ai: F, ba: [*]F, bb: [*]F, c: [*]F, ldc: BlasLong) callconv(cc) i32 {
            return run(.n, bm, bn, bk, ar, ai, ba, bb, c, ldc);
        }
        pub fn gemm_kernel_l(bm: BlasLong, bn: BlasLong, bk: BlasLong, ar: F, ai: F, ba: [*]F, bb: [*]F, c: [*]F, ldc: BlasLong) callconv(cc) i32 {
            return run(.l, bm, bn, bk, ar, ai, ba, bb, c, ldc);
        }
        pub fn gemm_kernel_r(bm: BlasLong, bn: BlasLong, bk: BlasLong, ar: F, ai: F, ba: [*]F, bb: [*]F, c: [*]F, ldc: BlasLong) callconv(cc) i32 {
            return run(.r, bm, bn, bk, ar, ai, ba, bb, c, ldc);
        }
        pub fn gemm_kernel_b(bm: BlasLong, bn: BlasLong, bk: BlasLong, ar: F, ai: F, ba: [*]F, bb: [*]F, c: [*]F, ldc: BlasLong) callconv(cc) i32 {
            return run(.b, bm, bn, bk, ar, ai, ba, bb, c, ldc);
        }

        /// `C += alpha * op(A) * op(B)`, `op` conjugating as `conj` says;
        /// the loops of the real kernel.
        fn run(comptime conj: Conj, bm: BlasLong, bn: BlasLong, bk: BlasLong, ar: F, ai: F, ba: [*]F, bb: [*]F, c: [*]F, ldc: BlasLong) i32 {
            const m: usize = @intCast(bm);
            const n: usize = @intCast(bn);
            const k: usize = @intCast(bk);
            const ld: usize = @intCast(ldc);
            var q0: usize = 0;
            while (q0 + nr <= n) : (q0 += nr) rows(conj, nr, m, k, ar, ai, ba, bb + 2 * q0 * k, c + 2 * q0 * ld, ld);
            inline for (comptime tails(nr)) |w| {
                if ((n - q0) & w != 0) {
                    rows(conj, w, m, k, ar, ai, ba, bb + 2 * q0 * k, c + 2 * q0 * ld, ld);
                    q0 += w;
                }
            }
            return 0;
        }

        fn rows(comptime conj: Conj, comptime wn: usize, m: usize, k: usize, ar: F, ai: F, ba: [*]const F, pb: [*]const F, c: [*]F, ld: usize) void {
            var r0: usize = 0;
            while (r0 + mr <= m) : (r0 += mr) @call(.always_inline, tile, .{ conj, mr, wn, k, ar, ai, ba + 2 * r0 * k, pb, c + 2 * r0, ld });
            inline for (comptime tails(mr)) |w| {
                if ((m - r0) & w != 0) {
                    @call(.never_inline, tile, .{ conj, w, wn, k, ar, ai, ba + 2 * r0 * k, pb, c + 2 * r0, ld });
                    r0 += w;
                }
            }
        }

        fn tile(comptime conj: Conj, comptime wm: usize, comptime wn: usize, k: usize, ar: F, ai: F, pa: [*]const F, pb: [*]const F, c: [*]F, ld: usize) void {
            if (geo.split) splitTile(conj, wm, wn, k, ar, ai, pa, pb, c, ld) else vectorTile(conj, wm, wn, k, ar, ai, pa, pb, c, ld);
        }

        /// `[1, 0, 3, 2, ...]`: each element's parts swapped.
        fn swapMask(comptime len: usize) @Vector(len, i32) {
            var s: [len]i32 = undefined;
            for (0..len) |i| s[i] = @intCast(i ^ 1);
            return s;
        }

        /// `x` in the real parts, `y` in the imaginary ones.
        fn alternate(comptime len: usize, x: F, y: F) @Vector(len, F) {
            var s: [len]F = undefined;
            inline for (0..len) |i| s[i] = if (i % 2 == 0) x else y;
            return s;
        }

        /// One `wm x wn` tile on vectors of whole elements: `acc1` sums A
        /// times the real parts of B, `acc2` A times the imaginary ones; the
        /// products' cross terms come from swapping `acc2`'s parts at the end.
        fn vectorTile(comptime conj: Conj, comptime wm: usize, comptime wn: usize, k: usize, ar: F, ai: F, pa: [*]const F, pb: [*]const F, c: [*]F, ld: usize) void {
            @setEvalBranchQuota(100_000);
            const lanes: usize = 2 * @as(usize, @min(vl, wm));
            const nv = 2 * wm / lanes;
            const V = @Vector(lanes, F);
            var acc1: [wn][nv]V = undefined;
            var acc2: [wn][nv]V = undefined;
            inline for (0..wn) |j| inline for (0..nv) |v| {
                acc1[j][v] = @splat(0);
                acc2[j][v] = @splat(0);
            };
            var a = pa;
            var b = pb;
            for (0..k) |_| {
                var av: [nv]V = undefined;
                inline for (0..nv) |v| av[v] = util.load(V, a + v * lanes);
                inline for (0..wn) |j| {
                    const br: V = @splat(b[2 * j]);
                    const bi: V = @splat(b[2 * j + 1]);
                    inline for (0..nv) |v| {
                        acc1[j][v] = madd(F, lanes, av[v], br, acc1[j][v]);
                        acc2[j][v] = madd(F, lanes, av[v], bi, acc2[j][v]);
                    }
                }
                a += 2 * wm;
                b += 2 * wn;
            }
            const swap = comptime swapMask(lanes);
            // `(x + iy)(u + iv)` with x, y from A and u, v from B: acc1 holds
            // (xu, yu), the swapped acc2 (yv, xv)
            const sb = comptime switch (conj) {
                .n => alternate(lanes, -1, 1),
                .l => alternate(lanes, 1, 1),
                .r => alternate(lanes, 1, -1),
                .b => alternate(lanes, -1, -1),
            };
            const sa = comptime alternate(lanes, 1, -1);
            const alr: V = @splat(ar);
            const ali = alternate(lanes, -ai, ai);
            inline for (0..wn) |j| inline for (0..nv) |v| {
                const s = @shuffle(F, acc2[j][v], undefined, swap);
                const first = if (comptime conj.ofA()) acc1[j][v] * sa else acc1[j][v];
                const res = madd(F, lanes, s, sb, first);
                const p = c + 2 * j * ld + v * lanes;
                const rs = @shuffle(F, res, undefined, swap);
                util.store(V, p, madd(F, lanes, res, alr, madd(F, lanes, rs, ali, util.load(V, p))));
            };
        }

        /// The same with each part a scalar, summed in `A`.
        fn splitTile(comptime conj: Conj, comptime wm: usize, comptime wn: usize, k: usize, ar: F, ai: F, pa: [*]const F, pb: [*]const F, c: [*]F, ld: usize) void {
            @setEvalBranchQuota(100_000);
            var re: [wn][wm]A = @splat(@splat(0));
            var im: [wn][wm]A = @splat(@splat(0));
            var a = pa;
            var b = pb;
            for (0..k) |_| {
                var x: [wm]A = undefined;
                var y: [wm]A = undefined;
                inline for (0..wm) |i| {
                    x[i] = a[2 * i];
                    y[i] = if (comptime conj.ofA()) -@as(A, a[2 * i + 1]) else a[2 * i + 1];
                }
                inline for (0..wn) |j| {
                    const u: A = b[2 * j];
                    const v: A = if (comptime conj.ofB()) -@as(A, b[2 * j + 1]) else b[2 * j + 1];
                    inline for (0..wm) |i| {
                        re[j][i] = madd(A, 1, -y[i], v, madd(A, 1, x[i], u, re[j][i]));
                        im[j][i] = madd(A, 1, x[i], v, madd(A, 1, y[i], u, im[j][i]));
                    }
                }
                a += 2 * wm;
                b += 2 * wn;
            }
            const war: A = ar;
            const wai: A = ai;
            inline for (0..wn) |j| inline for (0..wm) |i| {
                const p = c + 2 * (i + j * ld);
                const nre = madd(A, 1, war, re[j][i], madd(A, 1, -wai, im[j][i], p[0]));
                const nim = madd(A, 1, war, im[j][i], madd(A, 1, wai, re[j][i], p[1]));
                p[0] = @floatCast(nre);
                p[1] = @floatCast(nim);
            };
        }
    };
}

/// Every module with kernels over the gemm tile: each `Kernels(T, geo)` has
/// a public function per vtable slot it fills, named as the slot less its
/// precision letter.
pub const sets = .{ @This(), @import("../symm/kernel.zig"), @import("../trmm/kernel.zig"), @import("../trsm/kernel.zig") };

/// The precisions these kernels serve, by vtable prefix.
const precisions = .{ .{ "s", f32 }, .{ "d", f64 }, .{ "q", XDouble }, .{ "c", types.Complex32 }, .{ "z", types.Complex64 }, .{ "x", types.ComplexXD } };
/// Those a core's object may replace: `q` and `x` slots keep the portable kernels.
const core_precisions = .{ precisions[0], precisions[1], precisions[3], precisions[4] };

/// The slot that says whose kernels a precision has.
fn witness(comptime T: type) []const u8 {
    return if (isComplex(T)) "gemm_kernel_n" else "gemm_kernel";
}

/// The names of `K`'s functions that are slots with prefix `pre`.
fn slotNames(comptime pre: []const u8, comptime K: type) []const []const u8 {
    comptime {
        var out: []const []const u8 = &.{};
        for (std.meta.declarations(K)) |d| {
            if (@typeInfo(@TypeOf(@field(K, d))) == .@"fn" and @hasField(Vtable, pre ++ d)) out = out ++ [_][]const u8{d};
        }
        return out;
    }
}

/// The portable kernel of `T` named `name` (a slot less its letter), from
/// whichever of `sets` has it: `portable.kernel`'s source for level 3.
pub fn find(comptime T: type, comptime name: []const u8) @TypeOf(&@field(Owner(T, name).Kernels(T, geometry(T)), name)) {
    return &@field(Owner(T, name).Kernels(T, geometry(T)), name);
}

fn Owner(comptime T: type, comptime name: []const u8) type {
    inline for (sets) |s| {
        if (@hasDecl(s.Kernels(T, geometry(T)), name)) return s;
    }
    @compileError("no level-3 kernel " ++ name ++ " for " ++ @typeName(T));
}

/// Points `vt`'s level-3 slots at the kernels of each of `list`, and its
/// gemm unrolls at the tiles, through `source` (`function`, `unroll`).
fn fill(vt: *Vtable, comptime list: anytype, comptime source: anytype) void {
    inline for (list) |p| {
        inline for (sets) |s| {
            inline for (comptime slotNames(p[0], s.Kernels(p[1], geometry(p[1])))) |name| {
                @field(vt, p[0] ++ name) = source.function(p[0], p[1], s, name);
            }
        }
        const um = source.unroll(p[0], p[1], "m");
        const un = source.unroll(p[0], p[1], "n");
        @field(vt, p[0] ++ "gemm_unroll_m") = um;
        @field(vt, p[0] ++ "gemm_unroll_n") = un;
        @field(vt, p[0] ++ "gemm_unroll_mn") = @max(um, un);
    }
}

/// Puts back the portable kernels' unrolls where `vt`'s slots still hold
/// them: a core's blocking (`params.installCore`), set after them, writes
/// `param.h`'s, and no per-core object may come to replace both.
pub fn settle(vt: *Vtable) void {
    inline for (precisions) |p| {
        const K = Kernels(p[1], geometry(p[1]));
        if (@field(vt, p[0] ++ witness(p[1])) == &@field(K, witness(p[1]))) {
            @field(vt, p[0] ++ "gemm_unroll_m") = K.unroll_m;
            @field(vt, p[0] ++ "gemm_unroll_n") = K.unroll_n;
            @field(vt, p[0] ++ "gemm_unroll_mn") = @max(K.unroll_m, K.unroll_n);
        }
    }
}

/// The portable kernels: the main library's own build of these.
pub fn installGeneric(vt: *Vtable) void {
    fill(vt, precisions, struct {
        fn function(comptime pre: []const u8, comptime T: type, comptime s: type, comptime name: []const u8) @typeInfo(@FieldType(Vtable, pre ++ name)).optional.child {
            return &@field(s.Kernels(T, geometry(T)), name);
        }
        fn unroll(comptime _: []const u8, comptime T: type, comptime which: []const u8) i32 {
            return @field(Kernels(T, geometry(T)), "unroll_" ++ which);
        }
    });
}

/// The per-core kernels, for an architecture's registry (`kernels/<arch>.zig`):
/// every per-core object exports them, built for its own features, with
/// their unrolls as data (the tile depends on the object's features, which
/// the main library does not have). They never make a core need an object.
pub fn Registry(comptime Core: type) type {
    return struct {
        pub fn has(comptime core: Core) bool {
            _ = core;
            return false;
        }

        fn Unrolls(comptime T: type) type {
            return struct {
                const m: i32 = Kernels(T, geometry(T)).unroll_m;
                const n: i32 = Kernels(T, geometry(T)).unroll_n;
            };
        }

        pub fn exportCore(comptime core: Core) void {
            inline for (core_precisions) |p| {
                inline for (sets) |s| {
                    const K = s.Kernels(p[1], geometry(p[1]));
                    inline for (comptime slotNames(p[0], K)) |name| {
                        @export(&@field(K, name), .{ .name = util.kernelSymbol(p[0] ++ name, @tagName(core)), .visibility = .hidden });
                    }
                }
                inline for (.{ "m", "n" }) |w| {
                    @export(&@field(Unrolls(p[1]), w), .{ .name = util.kernelSymbol(p[0] ++ "gemm_unroll_" ++ w, @tagName(core)), .visibility = .hidden });
                }
            }
        }

        pub fn install(vt: *Vtable, comptime core: Core) void {
            fill(vt, core_precisions, struct {
                fn function(comptime pre: []const u8, comptime _: type, comptime _: type, comptime name: []const u8) @typeInfo(@FieldType(Vtable, pre ++ name)).optional.child {
                    const Fn = @typeInfo(@FieldType(Vtable, pre ++ name)).optional.child;
                    return @extern(Fn, .{ .name = util.kernelSymbol(pre ++ name, @tagName(core)), .visibility = .hidden });
                }
                fn unroll(comptime pre: []const u8, comptime _: type, comptime which: []const u8) i32 {
                    return @extern(*const i32, .{ .name = util.kernelSymbol(pre ++ "gemm_unroll_" ++ which, @tagName(core)), .visibility = .hidden }).*;
                }
            });
        }
    };
}

// ---- tests -----------------------------------------------------------------

fn fillSeq(comptime T: type, buf: []T, seed: usize) void {
    for (buf, 0..) |*v, i| v.* = @floatFromInt(@as(i64, @intCast((i * 7 + seed * 13) % 23)) - 11);
}

/// The packed kernel against plain arithmetic for every panel remainder, the
/// operands packed both ways (exact: small integers). Storage follows the
/// tile, so every size runs whatever the tile. Complex `T` runs every
/// conjugation; its parts are `fillSeq`'s floats in pairs.
fn checkLayout(comptime T: type, comptime geo: Geometry) !void {
    const K = Kernels(T, geo);
    const complex = comptime isComplex(T);
    const F = types.Float(T);
    const cs: usize = if (complex) 2 else 1;
    const al = std.testing.allocator;
    const sizes = [_]usize{ 1, 2, 3, 5, geo.mr() - 1, geo.mr(), geo.mr() + 1, 2 * geo.mr() + 3, geo.nr - 1, geo.nr + 1, 2 * geo.nr + 1 };
    const ld = std.mem.max(usize, &sizes);
    const kmax = 13;
    const bufs = try al.alloc(F, 9 * cs * ld * ld + 6 * cs * ld * kmax);
    defer al.free(bufs);
    var rest = bufs;
    const take = struct {
        fn f(r: *[]F, n: usize) []F {
            const s = r.*[0..n];
            r.* = r.*[n..];
            return s;
        }
    }.f;
    const a = take(&rest, cs * ld * kmax);
    const at = take(&rest, cs * kmax * ld);
    const b = take(&rest, cs * kmax * ld);
    const bt = take(&rest, cs * ld * kmax);
    const pa = take(&rest, cs * ld * kmax);
    const pb = take(&rest, cs * ld * kmax);
    const c = take(&rest, cs * ld * ld);
    const want = take(&rest, cs * ld * ld);
    const got = take(&rest, cs * ld * ld);
    fillSeq(F, a, 1);
    fillSeq(F, b, 2);
    fillSeq(F, c, 3);
    // `at` is A transposed (`k x m`, ld `kmax`), `bt` B transposed (`n x k`)
    for (0..ld) |i| for (0..kmax) |l| for (0..cs) |p| {
        at[cs * (l + i * kmax) + p] = a[cs * (i + l * ld) + p];
        bt[cs * (i + l * ld) + p] = b[cs * (l + i * kmax) + p];
    };
    const conjs: []const Conj = if (complex) &.{ .n, .l, .r, .b } else &.{.n};
    for (sizes) |m| for (sizes) |n| for ([_]usize{ 1, 4, kmax }) |k| for (conjs) |cj| {
        if (m == 0 or n == 0) continue;
        @memcpy(want, c);
        for (0..n) |j| for (0..m) |i| {
            var sr: F = 0;
            var si: F = 0;
            for (0..k) |l| {
                const x = a[cs * (i + l * ld)];
                const y = if (complex) a[cs * (i + l * ld) + 1] * @as(F, if (cj.ofA()) -1 else 1) else 0;
                const u = b[cs * (l + j * kmax)];
                const v = if (complex) b[cs * (l + j * kmax) + 1] * @as(F, if (cj.ofB()) -1 else 1) else 0;
                sr += x * u - y * v;
                si += x * v + y * u;
            }
            // alpha 2, or (2, -1) for complex
            want[cs * (i + j * ld)] += 2 * sr + if (complex) si else 0;
            if (complex) want[cs * (i + j * ld) + 1] += 2 * si - sr;
        };
        inline for (.{ false, true }) |trans| {
            @memcpy(got, c);
            const lm: BlasLong = @intCast(ld);
            if (trans) {
                _ = K.gemm_incopy(@intCast(k), @intCast(m), at.ptr, kmax, pa.ptr);
                _ = K.gemm_otcopy(@intCast(k), @intCast(n), bt.ptr, lm, pb.ptr);
            } else {
                _ = K.gemm_itcopy(@intCast(k), @intCast(m), a.ptr, lm, pa.ptr);
                _ = K.gemm_oncopy(@intCast(k), @intCast(n), b.ptr, kmax, pb.ptr);
            }
            if (complex) {
                const f = switch (cj) {
                    .n => &K.gemm_kernel_n,
                    .l => &K.gemm_kernel_l,
                    .r => &K.gemm_kernel_r,
                    .b => &K.gemm_kernel_b,
                };
                _ = f(@intCast(m), @intCast(n), @intCast(k), 2, -1, pa.ptr, pb.ptr, got.ptr, lm);
            } else {
                _ = K.gemm_kernel(@intCast(m), @intCast(n), @intCast(k), 2, pa.ptr, pb.ptr, got.ptr, lm);
            }
            if (!std.mem.eql(F, want, got)) return error.TestUnexpectedResult;
        }
    };
}

test "packed kernels match plain arithmetic for every panel remainder" {
    try checkLayout(f64, geometry(f64));
    try checkLayout(f32, geometry(f32));
    try checkLayout(f64, .{ .vl = 1, .mv = 3, .nr = 6 });
    try checkLayout(f32, .{ .vl = 4, .mv = 2, .nr = 5 });
    // AVX-512 and RVV zvl512b tiles, wider than the old fixed storage
    try checkLayout(f32, .{ .vl = 16, .mv = 2, .nr = 14 });
    try checkLayout(f64, .{ .vl = 8, .mv = 2, .nr = 14 });
}

test "complex packed kernels match plain arithmetic for every panel remainder and conjugation" {
    try checkLayout(types.Complex64, geometry(types.Complex64));
    try checkLayout(types.Complex32, geometry(types.Complex32));
    try checkLayout(types.Complex64, .{ .vl = 1, .mv = 3, .nr = 5, .split = true });
    try checkLayout(types.Complex32, .{ .vl = 2, .mv = 3, .nr = 3 });
    try checkLayout(types.Complex32, .{ .vl = 8, .mv = 2, .nr = 6 });
}

test "geometry: the most accumulators the registers hold" {
    const cpu = comptime blk: {
        var c = std.Target.Cpu.baseline(.x86_64, sysinfo.os);
        c.features.addFeature(@backingInt(std.Target.x86.Feature.avx2));
        c.features.populateDependencies(&std.Target.x86.all_features);
        break :blk c;
    };
    try std.testing.expectEqual(Geometry{ .vl = 8, .mv = 2, .nr = 6 }, comptime geometryFor(f32, cpu));
    try std.testing.expectEqual(Geometry{ .vl = 4, .mv = 2, .nr = 6 }, comptime geometryFor(f64, cpu));
    // the Haswell `zgemm` tile, 4 x 2
    try std.testing.expectEqual(Geometry{ .vl = 2, .mv = 2, .nr = 2 }, comptime geometryFor(types.Complex64, cpu));
    const arm64 = comptime std.Target.Cpu.baseline(.aarch64, sysinfo.os);
    try std.testing.expectEqual(Geometry{ .vl = 4, .mv = 4, .nr = 6 }, comptime geometryFor(f32, arm64));
    try std.testing.expectEqual(Geometry{ .vl = 2, .mv = 4, .nr = 6 }, comptime geometryFor(f64, arm64));
    try std.testing.expectEqualSlices(usize, &.{ 4, 2, 1 }, comptime tails(6));
    try std.testing.expectEqualSlices(usize, &.{ 8, 4, 2, 1 }, comptime tails(9));
}

test "beta zero does not read c; other values scale the block only" {
    const K = Kernels(f64, geometry(f64));
    var c = [_]f64{ std.math.nan(f64), 2, 3, 4, 5, 6, 7, 8 };
    _ = K.gemm_beta(2, 2, 0, 0, &c, 0, &c, 0, &c, 4);
    try std.testing.expectEqualSlices(f64, &.{ 0, 0, 3, 4, 0, 0, 7, 8 }, &c);
    c[4] = 5;
    _ = K.gemm_beta(2, 2, 0, 2, &c, 0, &c, 0, &c, 4);
    try std.testing.expectEqualSlices(f64, &.{ 0, 0, 3, 4, 10, 0, 7, 8 }, &c);
    const Z = Kernels(types.Complex64, geometry(types.Complex64));
    var z = [_]f64{ std.math.nan(f64), 1, 3, 4, 1, 2, 7, 8 };
    _ = Z.gemm_beta(1, 2, 0, 0, 0, &z, 0, &z, 0, &z, 2);
    try std.testing.expectEqualSlices(f64, &.{ 0, 0, 3, 4, 0, 0, 7, 8 }, &z);
    // (1 + 2i) i, then (-2 + i)(2 + i)
    z[0] = 1;
    z[1] = 2;
    _ = Z.gemm_beta(1, 1, 0, 0, 1, &z, 0, &z, 0, &z, 2);
    try std.testing.expectEqualSlices(f64, &.{ -2, 1, 3, 4, 0, 0, 7, 8 }, &z);
    _ = Z.gemm_beta(1, 1, 0, 2, 1, &z, 0, &z, 0, &z, 2);
    try std.testing.expectEqualSlices(f64, &.{ -5, 0, 3, 4, 0, 0, 7, 8 }, &z);
}
