//! x86_64 real rot kernels, one set per core, built into the per-core kernel
//! objects (`src/linalg/blas/_kernel.zig`). A core's kernel is a `.S` kernel
//! (`asm/`) or the C driver `srot.c`/`drot.c` (`../kernel.zig`) around the
//! block `uk/__root__.zig` names, which upstream is intrinsics
//! (`<p>rot_microk_haswell-2.c`, `_skylakex-2.c`), here the same loops on
//! `@Vector`s. Tests live in `rot/__root__.zig`.

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

/// The Skylake-X block is 512 bits wide, Haswell's 256.
const has_avx512f = std.Target.x86.featureSetHas(sysinfo.target.cpu.features, .avx512f);

/// `<p>rot_kernel` of the microkernels: four vectors at a time, then one,
/// then a masked vector, which here is the same fused pair per element.
pub fn VectorBlock(comptime slot: Slot) type {
    const F = kernel.Elem(slot);
    const l = (if (has_avx512f) 64 else 32) / @sizeOf(F);
    const V = @Vector(l, F);
    return struct {
        pub fn run(n: BlasLong, x: [*]F, y: [*]F, c: F, s: F) void {
            const len: usize = @intCast(n);
            const cv: V = @splat(c);
            const sv: V = @splat(s);
            const x4 = len & ~@as(usize, 4 * l - 1);
            const x1 = len & ~@as(usize, l - 1);
            var i: usize = 0;
            while (i < x4) : (i += 4 * l) kernel.rotateVectors(V, 4, x, y, i, cv, sv);
            while (i < x1) : (i += l) kernel.rotateVectors(V, 1, x, y, i, cv, sv);
            while (i < len) : (i += 1) kernel.rotate(F, true, &x[i], &y[i], c, s);
        }
    };
}

const K = uk.Kernels(Core, kernel, @import("asm/__root__.zig"), @import("uk/__root__.zig"), VectorBlock);

pub const hasSlot = K.hasSlot;
pub const skips = K.skips;
pub const has = K.has;
pub const exportCore = K.exportCore;
pub const install = K.install;
