//! The Zig gemm kernels for real `T` and its `Complex`: the packing copies,
//! `beta` and the register tile, in the main library and every per-core
//! object (`../kernels.zig`) alike; layout and tile in `notes.md`. Imports stay in `linalg/blas/`.

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("../../types.zig");
const fma = @import("../../fma.zig");
const util = @import("../../util.zig");
const loops = @import("loops.zig");

const BlasLong = types.BlasLong;
const isComplex = types.isComplex;

// s390x vector features: `vx`, `vxeh` on Zig 0.18; `vector`, `vector_enhancements_1` on 0.17 (issue #266)
const S390x = std.Target.s390x.Feature;
const s390x_vx: S390x = @field(S390x, if (@hasField(S390x, "vx")) "vx" else "vector");
const s390x_vxeh: S390x = @field(S390x, if (@hasField(S390x, "vxeh")) "vxeh" else "vector_enhancements_1");

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
/// none, or for lanes that flush denormals (`util.vectorIeeeOn`). Not
/// `std.simd.suggestVectorLength`: that counts SVE as 256-bit fixed vectors
/// (LLVM splits them into NEON pairs) and leaves AVX at 128.
fn vectorBits(comptime T: type, comptime cpu: std.Target.Cpu) u16 {
    const double = T == f64;
    if ((T != f32 and T != f64) or !util.vectorIeeeOn(T, cpu)) return 0;
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
        .s390x => if (cpu.has(.s390x, s390x_vx) and (double or cpu.has(.s390x, s390x_vxeh))) 128 else 0,
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
        // wasm has locals, not registers: ours, an x86_64 engine's 16 XMM
        .wasm32, .wasm64 => 16,
        // o32's 32-bit FPU registers pair up for doubles; MSA has 32 of 128 bits
        // (and Zig builds it with 64-bit FPU registers, `-mfp64`)
        .mips, .mipsel => if (cpu.has(.mips, .msa)) 32 else 16,
        .sparc => 16,
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
    return comptime geometryFor(T, sysinfo.target.cpu);
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

pub fn Lanes(comptime T: type, comptime n: usize) type {
    return if (n == 1) T else @Vector(n, T);
}

// whole arrays, not `util.load` and `util.store`: the self-hosted x86_64
// backend stores a row tail's two-lane `f32` vector as sixteen bytes
pub inline fn load(comptime T: type, comptime n: usize, p: [*]const T) Lanes(T, n) {
    return if (n == 1) p[0] else p[0..n].*;
}

pub inline fn store(comptime T: type, comptime n: usize, p: [*]T, v: Lanes(T, n)) void {
    if (n == 1) p[0] = v else p[0..n].* = v;
}

pub inline fn splat(comptime T: type, comptime n: usize, x: T) Lanes(T, n) {
    return if (n == 1) x else @splat(x);
}

/// `a*b + c` under the build's contraction policy (`fma.zig`). On the x87
/// each step is rounded to `T`: a full tile's sums stay in registers, wider
/// than a partial tile's in memory, and threads move columns between them.
pub inline fn madd(comptime T: type, comptime n: usize, a: Lanes(T, n), b: Lanes(T, n), c: Lanes(T, n)) Lanes(T, n) {
    const fused = comptime fma.fuseFor(T);
    if (comptime n == 1 and !fused and util.onX87(T)) return util.rounded(T, util.rounded(T, a * b) + c);
    return fma.madd(fused, Lanes(T, n), a, b, c);
}

/// Packs `n` lines of depth `m` of elements `E` in panels of `u`. `across`:
/// line `j` is element `j` of each of `m` rows `lda` apart (`tcopy`); else it
/// is column `j`, `lda` apart (`ncopy`).
pub fn pack(comptime E: type, comptime u: usize, comptime across: bool, m: BlasLong, n: BlasLong, a: [*]const E, lda: BlasLong, b: [*]E) void {
    const depth: usize = @intCast(m);
    const lines: usize = @intCast(n);
    const ld: usize = @intCast(lda);
    var c0: usize = 0;
    while (c0 + u <= lines) : (c0 += u) panel(E, u, across, depth, a, ld, c0, b + depth * c0);
    if (c0 < lines) packRest(E, across, depth, a, ld, c0, lines, b);
}

/// The panels past the whole ones, widths as values: they are few, and one
/// function for every width keeps debug builds small.
fn packRest(comptime E: type, comptime across: bool, depth: usize, a: [*]const E, ld: usize, from: usize, lines: usize, b: [*]E) void {
    var c0 = from;
    while (c0 < lines) {
        const w = std.math.floorPowerOfTwo(usize, lines - c0);
        const out = b + depth * c0;
        for (0..depth) |l| for (0..w) |c| {
            out[l * w + c] = if (across) a[l * ld + c0 + c] else a[(c0 + c) * ld + l];
        };
        c0 += w;
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

/// Packs `lines` lines of `depth` as `pack` does, in panels of `u`, element
/// `(j, l)` being `src.at(j, l)`; where that is an optional, null leaves it
/// unwritten. For the other routines' copies, which transform what they pack;
/// the widths are run-time values, so each copy is one small function.
pub fn packWith(comptime E: type, u: usize, depth: usize, lines: usize, dst: [*]E, src: anytype) void {
    var c0: usize = 0;
    while (c0 < lines) {
        const w = panelWidth(lines - c0, u);
        const out = dst + depth * c0;
        for (0..depth) |l| for (0..w) |c| {
            const v = src.at(c0 + c, l);
            if (@TypeOf(v) == ?E) {
                if (v) |x| out[l * w + c] = x;
            } else out[l * w + c] = v;
        };
        c0 += w;
    }
}

/// A copy in the reference's ABI on `packWith`: `m` deep, `n` lines into
/// `b`, from `src`; the trmm, trsm, symm and 3M copies.
pub fn copyFrom(comptime E: type, u: usize, m: BlasLong, n: BlasLong, b: anytype, src: anytype) i32 {
    packWith(E, u, @intCast(m), @intCast(n), @ptrCast(b), src);
    return 0;
}

/// `C += alpha * A * B` on packed panels through `G`'s gemm kernel
/// (`Kernels(T, geo)`), conjugating the triangle's side, `left` or right,
/// for `conj`: trmm's and trsm's kernels. `alpha` is `[re, im]`, or `[re]`.
pub fn sideKernel(comptime G: type, comptime F: type, comptime complex: bool, left: bool, conj: bool, m: usize, n: usize, k: usize, alpha: [*]const F, a: anytype, b: anytype, c: anytype, ldc: usize) void {
    const bm: BlasLong = @intCast(m);
    const bn: BlasLong = @intCast(n);
    const bk: BlasLong = @intCast(k);
    if (comptime complex) {
        const f = if (!conj) &G.gemm_kernel_n else if (left) &G.gemm_kernel_l else &G.gemm_kernel_r;
        _ = f(bm, bn, bk, alpha[0], alpha[1], @ptrCast(a), @ptrCast(b), @ptrCast(c), @intCast(ldc));
    } else {
        _ = G.gemm_kernel(bm, bn, bk, alpha[0], @ptrCast(a), @ptrCast(b), @ptrCast(c), @intCast(ldc));
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
                    for (col) |*v| v.* = util.scalar(T, v.* * b);
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
        /// `mr`) and `B` packed `bk x bn` (columns in panels of `nr`), through
        /// `loops.zig`.
        pub fn gemm_kernel(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha: T, ba: [*]T, bb: [*]T, c: [*]T, ldc: BlasLong) callconv(cc) i32 {
            L.run(@intCast(bm), @intCast(bn), @intCast(bk), alpha, ba, bb, c, @intCast(ldc));
            return 0;
        }

        const L = loops.Loops(T, geo, Out);

        /// `C += alpha * sum`, rounded back to `T`.
        const Out = struct {
            pub const C = T;
            pub const scale = 1;
            pub const Alpha = T;

            pub inline fn prep(comptime lanes: usize, alpha: T) Lanes(A, lanes) {
                return splat(A, lanes, alpha);
            }

            pub inline fn put(p: [*]T, sum: A, alpha: T) void {
                const r = madd(A, 1, sum, alpha, p[0]);
                p[0] = if (A == T) r else @floatCast(r);
            }

            pub inline fn putLanes(comptime lanes: usize, p: [*]T, sums: Lanes(A, lanes), al: Lanes(A, lanes)) void {
                const r = madd(A, lanes, sums, al, L.wide(lanes, p));
                if (A == T) store(T, lanes, p, r) else p[0] = @floatCast(r);
            }
        };
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
                        col[2 * i] = util.scalar(F, br * x - bi * y);
                        col[2 * i + 1] = util.scalar(F, br * y + bi * x);
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
        /// the loops of the real kernel. `conj` comptime: a run-time one
        /// selects between vectors, which LLVM cannot compile on MIPS64r6
        /// (issue #192).
        fn run(comptime conj: Conj, bm: BlasLong, bn: BlasLong, bk: BlasLong, ar: F, ai: F, ba: [*]F, bb: [*]F, c: [*]F, ldc: BlasLong) i32 {
            const m: usize = @intCast(bm);
            const n: usize = @intCast(bn);
            const k: usize = @intCast(bk);
            const ld: usize = @intCast(ldc);
            if (geo.split) {
                // one tile function, shaped by values: the scalar fallback's
                // tiles would otherwise be one per width pair
                var p0: usize = 0;
                while (p0 < n) {
                    const wn = panelWidth(n - p0, nr);
                    var r0: usize = 0;
                    while (r0 < m) {
                        const wm = panelWidth(m - r0, mr);
                        splitTile(conj, wm, wn, k, ar, ai, ba + 2 * r0 * k, bb + 2 * p0 * k, c + 2 * (r0 + p0 * ld), ld);
                        r0 += wm;
                    }
                    p0 += wn;
                }
                return 0;
            }
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
            if (geo.split) unreachable;
            vectorTile(conj, wm, wn, k, ar, ai, pa, pb, c, ld);
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
                inline for (0..nv) |v| av[v] = a[v * lanes ..][0..lanes].*;
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
            const sb = switch (conj) {
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
                const first = if (conj.ofA()) acc1[j][v] * sa else acc1[j][v];
                const res = madd(F, lanes, s, sb, first);
                const p = c + 2 * j * ld + v * lanes;
                const rs = @shuffle(F, res, undefined, swap);
                // whole arrays, not `util.store`: the self-hosted x86_64 backend
                // stores a two-lane `f32` vector as sixteen bytes
                const cv: V = p[0..lanes].*;
                p[0..lanes].* = madd(F, lanes, res, alr, madd(F, lanes, rs, ali, cv));
            };
        }

        /// The same with each part a scalar, summed in `A`.
        fn splitTile(conj: Conj, wm: usize, wn: usize, k: usize, ar: F, ai: F, pa: [*]const F, pb: [*]const F, c: [*]F, ld: usize) void {
            // the tile's most, of which `wm x wn` are used
            var re: [nr][mr]A = @splat(@splat(0));
            var im: [nr][mr]A = @splat(@splat(0));
            var a = pa;
            var b = pb;
            // a conjugated operand's imaginary parts change sign
            const ya: A = if (conj.ofA()) -1 else 1;
            const yb: A = if (conj.ofB()) -1 else 1;
            // loops over the tile, not unrolled: a debug build of the unrolled
            // tile overran its frame on MIPS release 6
            for (0..k) |_| {
                var x: [mr]A = undefined;
                var y: [mr]A = undefined;
                for (0..wm) |i| {
                    x[i] = a[2 * i];
                    y[i] = @as(A, a[2 * i + 1]) * ya;
                }
                for (0..wn) |j| {
                    const u: A = b[2 * j];
                    const v: A = @as(A, b[2 * j + 1]) * yb;
                    for (0..wm) |i| {
                        re[j][i] = madd(A, 1, -y[i], v, madd(A, 1, x[i], u, re[j][i]));
                        im[j][i] = madd(A, 1, x[i], v, madd(A, 1, y[i], u, im[j][i]));
                    }
                }
                a += 2 * wm;
                b += 2 * wn;
            }
            const war: A = ar;
            const wai: A = ai;
            for (0..wn) |j| for (0..wm) |i| {
                const p = c + 2 * (i + j * ld);
                const nre = madd(A, 1, war, re[j][i], madd(A, 1, -wai, im[j][i], p[0]));
                const nim = madd(A, 1, war, im[j][i], madd(A, 1, wai, re[j][i], p[1]));
                p[0] = @floatCast(nre);
                p[1] = @floatCast(nim);
            };
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
    return checkConj(T, geo, if (comptime isComplex(T)) &.{ .n, .l, .r, .b } else &.{.n});
}

/// `checkLayout` for the conjugations `conjs` only, which links only their
/// kernels: a wide complex tile's four are past the sparc64 rig's window.
fn checkConj(comptime T: type, comptime geo: Geometry, comptime conjs: []const Conj) !void {
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
                const f = inline for (conjs) |kc| {
                    if (kc == cj) break &@field(K, "gemm_kernel_" ++ @tagName(kc));
                } else unreachable;
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
}

test "complex packed kernels of an 8-lane tile match plain arithmetic, A plain" {
    try checkConj(types.Complex32, .{ .vl = 8, .mv = 2, .nr = 6 }, &.{ .n, .r });
}

test "complex packed kernels of an 8-lane tile match plain arithmetic, A conjugated" {
    try checkConj(types.Complex32, .{ .vl = 8, .mv = 2, .nr = 6 }, &.{ .l, .b });
}

test "geometry: the most accumulators the registers hold" {
    const cpu = comptime blk: {
        var c = std.Target.Cpu.baseline(.x86_64, sysinfo.target.os);
        c.features.addFeature(@backingInt(std.Target.x86.Feature.avx2));
        c.features.populateDependencies(&std.Target.x86.all_features);
        break :blk c;
    };
    try std.testing.expectEqual(Geometry{ .vl = 8, .mv = 2, .nr = 6 }, comptime geometryFor(f32, cpu));
    try std.testing.expectEqual(Geometry{ .vl = 4, .mv = 2, .nr = 6 }, comptime geometryFor(f64, cpu));
    // the Haswell `zgemm` tile, 4 x 2
    try std.testing.expectEqual(Geometry{ .vl = 2, .mv = 2, .nr = 2 }, comptime geometryFor(types.Complex64, cpu));
    const arm64 = comptime std.Target.Cpu.baseline(.aarch64, sysinfo.target.os);
    try std.testing.expectEqual(Geometry{ .vl = 4, .mv = 4, .nr = 6 }, comptime geometryFor(f32, arm64));
    try std.testing.expectEqual(Geometry{ .vl = 2, .mv = 4, .nr = 6 }, comptime geometryFor(f64, arm64));
    // 32-bit MIPS with MSA has the 64-bit cores' 32 vector registers
    const p5600 = comptime blk: {
        var c = std.Target.mips.cpu.p5600.toCpu(.mipsel);
        c.features.addFeature(@backingInt(std.Target.mips.Feature.msa));
        break :blk c;
    };
    // the I6400 with its SIMD unit: newer Zig's model lists it no more
    const r6 = comptime blk: {
        var c = std.Target.mips.cpu.i6400.toCpu(.mips64el);
        c.features.addFeature(@backingInt(std.Target.mips.Feature.msa));
        break :blk c;
    };
    inline for (.{ f32, f64, types.Complex32, types.Complex64 }) |T|
        try std.testing.expectEqual(comptime geometryFor(T, r6), comptime geometryFor(T, p5600));
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
