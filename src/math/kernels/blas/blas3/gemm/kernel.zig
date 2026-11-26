//! The Zig gemm kernels for real `T` and its `Complex`: the packing copies,
//! `beta` and the register tile, in the main library and every per-core
//! object (`../kernels.zig`) alike; layout and tile in `notes.md`. Imports stay in `linalg/blas/`.

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("../../types.zig");
const fma = @import("../../fma.zig");
const util = @import("../../util.zig");
const kpack = @import("kernel/pack.zig");
const kreal = @import("kernel/real.zig");
const kcplx = @import("kernel/complex.zig");

const BlasLong = types.BlasLong;
const isComplex = types.isComplex;
// s390x vector features: `vx`, `vxeh` on Zig 0.18; `vector`, `vector_enhancements_1` on 0.17 (issue #266)
const S390x = std.Target.s390x.Feature;
pub const copyFrom = kpack.copyFrom;
pub const pack = kpack.pack;
pub const packWith = kpack.packWith;
pub const sideKernel = kpack.sideKernel;
const RealKernels = kreal.RealKernels;
const ComplexKernels = kcplx.ComplexKernels;

const s390x_vx: S390x = @field(S390x, if (@hasField(S390x, "vx")) "vx" else "vector");
const s390x_vxeh: S390x = @field(S390x, if (@hasField(S390x, "vxeh")) "vxeh" else "vector_enhancements_1");

/// A kernel set's tile: `mv` vectors of `vl` lanes (`mr` rows) by `nr`
/// columns; `mr` and `nr` are also the copies' unrolls. For complex, `vl`
/// counts elements, and `split` holds their parts in scalars. `contracted`
/// complex tiles keep one accumulator per vector and column (`contractedTile`).
pub const Geometry = struct {
    vl: u16,
    mv: u16,
    nr: u16,
    split: bool = false,
    contracted: bool = false,

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
/// With FMA3 on 256-bit vectors, complex tiles are contracted, as OpenBLAS's
/// Haswell `cgemm` and `zgemm` (T21): one accumulator per vector and column,
/// one broadcast of B at a time, A then A swapped.
pub fn geometryFor(comptime T: type, comptime cpu: std.Target.Cpu) Geometry {
    const complex = isComplex(T);
    const F = types.Float(T);
    const raw_bits = vectorBits(F, cpu);
    // a complex lane is two floats
    const split = complex and raw_bits < 2 * @bitSizeOf(F);
    const bits: u16 = if (split) 0 else raw_bits;
    const vl: u16 = if (bits == 0) 1 else bits / @bitSizeOf(F) / @as(u16, if (complex) 2 else 1);
    const regs = util.registers(bits, cpu);
    const contracted = complex and !split and bits == 256 and (cpu.arch == .x86_64 or cpu.arch == .x86) and cpu.has(.x86, .fma);
    @setEvalBranchQuota(100_000);
    inline for (.{ true, false }) |strict| {
        var best: ?Geometry = null;
        var best_loads: usize = 0;
        for (1..5) |mv| for (1..17) |nr| {
            const g: Geometry = .{ .vl = vl, .mv = mv, .nr = nr, .split = split, .contracted = contracted };
            // accumulators, registers, loads and multiply-adds per step
            const acc = if (complex and !contracted) 2 * mv * nr else mv * nr;
            const a_regs = if (split) 2 * mv else mv;
            const b_regs = if (contracted) 1 else if (complex) 2 else if (byLane(cpu, bits)) (nr + vl - 1) / vl else 1;
            const loads = (if (contracted) 2 * a_regs else a_regs) + if (complex) 2 * nr else nr;
            const madds = if (split) 4 * mv * nr else if (contracted) 2 * acc else acc;
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
                    if (kc == cj) break switch (kc) {
                        .n => &K.gemm_kernel_n,
                        .l => &K.gemm_kernel_l,
                        .r => &K.gemm_kernel_r,
                        .b => &K.gemm_kernel_b,
                    };
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
    // Haswell's and Zen's contracted tiles (T21)
    try checkLayout(types.Complex32, .{ .vl = 4, .mv = 2, .nr = 6, .contracted = true });
    try checkLayout(types.Complex64, .{ .vl = 2, .mv = 2, .nr = 6, .contracted = true });
}

test "complex packed kernels of an 8-lane tile match plain arithmetic, A plain" {
    try checkConj(types.Complex32, .{ .vl = 8, .mv = 2, .nr = 6 }, &.{ .n, .r });
}

test "complex packed kernels of an 8-lane tile match plain arithmetic, A conjugated" {
    try checkConj(types.Complex32, .{ .vl = 8, .mv = 2, .nr = 6 }, &.{ .l, .b });
}

test "contracted complex tiles match plain arithmetic under every conjugation" {
    try checkConj(types.Complex32, .{ .vl = 4, .mv = 2, .nr = 6, .contracted = true }, &.{ .n, .l, .r, .b });
    try checkConj(types.Complex64, .{ .vl = 2, .mv = 2, .nr = 6, .contracted = true }, &.{ .n, .l, .r, .b });
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
    // with FMA3, OpenBLAS's contracted Haswell tiles, `cgemm` 8 x 6 and `zgemm` 4 x 6
    const fma3 = comptime blk: {
        var c = cpu;
        c.features.addFeature(@backingInt(std.Target.x86.Feature.fma));
        break :blk c;
    };
    try std.testing.expectEqual(Geometry{ .vl = 4, .mv = 2, .nr = 6, .contracted = true }, comptime geometryFor(types.Complex32, fma3));
    try std.testing.expectEqual(Geometry{ .vl = 2, .mv = 2, .nr = 6, .contracted = true }, comptime geometryFor(types.Complex64, fma3));
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
