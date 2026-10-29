//! x86_64 gemv kernels, one set per core, built into the per-core kernel
//! objects (`src/linalg/blas/_kernel.zig`). A core's kernel is a `.S` kernel
//! (`asm/`) or the C driver `{s,d,c,z}gemv_{n,t}_4.c` (`../kernel.zig`) around
//! the blocks `uk/index_n.zig` and `uk/index_t.zig` name: the reference's
//! inline-assembly microkernels, or for the AVX-512 cores' `dgemv_n`, whose
//! blocks are intrinsics, the same arithmetic on `@Vector`s. Skylake-X's
//! `sgemv` takes a second, intrinsics block besides Haswell's
//! (`HAVE_SGEMV_N_SKYLAKE_KERNEL`), which is not ported, so those cores keep
//! the portable `sgemv` (`uk_block.Kernels`). Tests live in `gemv/__root__.zig`.

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("../../../types.zig");
const coretypes = @import("../../../cpu/_coretypes.zig");
const fma = @import("../../../fma.zig");
const util = @import("../../../util.zig");
const kernel = @import("../kernel.zig");
const uk = @import("../../../blas1/uk_block.zig");
const index_n = @import("uk/index_n.zig");
const index_t = @import("uk/index_t.zig");

const BlasLong = types.BlasLong;
const Slot = util.Slot;

pub const Core = coretypes.X86CoreType;

/// The reference compiles the blocks with a 512-bit pass first when
/// `__AVX512CD__` is defined.
const has_avx512cd = std.Target.x86.featureSetHas(sysinfo.target.cpu.features, .avx512cd);

/// `dgemv_n_microk_skylakex-4.c`: `y += alpha * (a0*x0 + a1*x1 + ...)` per row,
/// as vector expressions the compiler contracts (`fuse`): the sum's first
/// product fused with the second, each later one into it, and `y` fused with
/// `sum * alpha`. 512-bit vectors first where `wide`, then 256-bit ones.
pub fn VectorBlockWith(comptime slot: Slot, comptime wide: bool, comptime fuse: bool) type {
    if (slot != .dgemv_n) @compileError(@tagName(slot) ++ " has no vector block");
    const F = f64;
    const V4 = @Vector(4, F);
    const V8 = @Vector(8, F);
    return struct {
        fn rows(comptime V: type, comptime cols: comptime_int, ap: [*]const [*]F, x: [*]F, y: [*]F, alpha: F, i: usize) void {
            var sum = fma.madd(fuse, V, util.load(V, ap[0] + i), @splat(x[0]), util.load(V, ap[1] + i) * @as(V, @splat(x[1])));
            inline for (2..cols) |c| sum = fma.madd(fuse, V, util.load(V, ap[c] + i), @splat(x[c]), sum);
            util.store(V, y + i, fma.madd(fuse, V, sum, @splat(alpha), util.load(V, y + i)));
        }

        fn run(comptime cols: comptime_int, n: BlasLong, ap: [*]const [*]F, x: [*]F, y: [*]F, alpha: F) void {
            const len: usize = @intCast(n);
            var i: usize = 0;
            if (wide and cols == 4) {
                const n5 = len & ~@as(usize, 7);
                while (i < n5) : (i += 8) rows(V8, cols, ap, x, y, alpha, i);
            }
            while (i < len) : (i += 4) rows(V4, cols, ap, x, y, alpha, i);
        }

        pub fn dgemv_kernel_4x4(n: BlasLong, ap: [*]const [*]F, x: [*]F, y: [*]F, alpha: [*]const F) void {
            run(4, n, ap, x, y, alpha[0]);
        }

        pub fn dgemv_kernel_4x2(n: BlasLong, ap: [*]const [*]F, x: [*]F, y: [*]F, alpha: [*]const F) void {
            run(2, n, ap, x, y, alpha[0]);
        }
    };
}

fn VectorBlock(comptime slot: Slot) type {
    return VectorBlockWith(slot, has_avx512cd, fma.fuse);
}

/// Both indexes as one.
const ukidx = struct {
    pub const entries = index_n.entries ++ index_t.entries;
    pub fn find(comptime slot: Slot, comptime core: []const u8) ?uk.Entry {
        return index_n.find(slot, core) orelse index_t.find(slot, core);
    }
};

const K = uk.Kernels(Core, kernel, @import("asm/__root__.zig"), ukidx, VectorBlock);

pub const hasSlot = K.hasSlot;
pub const skips = K.skips;
pub const has = K.has;
pub const exportCore = K.exportCore;
pub const install = K.install;
