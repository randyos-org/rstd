//! x86_64 real scal kernels, one set per core, built into the per-core kernel
//! objects (`src/linalg/blas/_kernel.zig`). A core's kernel is a `.S` kernel
//! (`asm/`) or the C driver `sscal.c`/`dscal.c` (`../kernel.zig`) around the
//! blocks `uk/__root__.zig` names: the reference's inline-assembly
//! microkernels, or for the AVX-512 cores, whose blocks are intrinsics, the
//! same loops on `@Vector`s. Tests live in `scal/__root__.zig`.

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("../../../types.zig");
const coretypes = @import("../../../cpu/_coretypes.zig");
const util = @import("../../../util.zig");
const kernel = @import("../kernel.zig");
const uk = @import("../../uk_block.zig");

const BlasLong = types.BlasLong;
const Slot = util.Slot;

pub const Core = coretypes.X86CoreType;

/// The reference compiles the blocks with 512-bit vectors when
/// `__AVX512CD__` is defined, else 256-bit ones.
const has_avx512cd = std.Target.x86.featureSetHas(sysinfo.target.cpu.features, .avx512cd);

/// `sscal_microk_skylakex-2.c` and `dscal_microk_skylakex-2.c`: `x *= alpha`
/// and its zeroing twin over a whole block, one vector at a time.
fn VectorBlock(comptime slot: Slot) type {
    const F = kernel.Elem(slot);
    const lanes = (if (has_avx512cd) 64 else 32) / @sizeOf(F);
    const V = @Vector(lanes, F);
    return struct {
        pub fn run(n: BlasLong, alpha: [*]const F, x: [*]F) void {
            const a: V = @splat(alpha[0]);
            const len: usize = @intCast(n);
            var i: usize = 0;
            while (i < len) : (i += lanes) util.store(V, x + i, a * util.load(V, x + i));
        }
        pub fn zero(n: BlasLong, _: [*]const F, x: [*]F) void {
            const len: usize = @intCast(n);
            var i: usize = 0;
            while (i < len) : (i += lanes) util.store(V, x + i, @splat(0));
        }
    };
}

const K = uk.Kernels(Core, kernel, @import("asm/__root__.zig"), @import("uk/__root__.zig"), VectorBlock);

pub const hasSlot = K.hasSlot;
pub const skips = K.skips;
pub const has = K.has;
pub const exportCore = K.exportCore;
pub const install = K.install;
