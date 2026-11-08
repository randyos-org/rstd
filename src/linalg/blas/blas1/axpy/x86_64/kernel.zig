//! x86_64 axpy kernels per core (`asm/`, or `../driver.zig` around a `uk/` block or `@Vector`s for AVX-512); tests in `../__root__.zig`.

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("../../../types.zig");
const coretypes = @import("../../../cpu/coretypes.zig");
const fma = @import("../../../fma.zig");
const util = @import("../../../util.zig");
const kernel = @import("../driver.zig");
const uk = @import("../../../uk_block.zig");

const BlasLong = types.BlasLong;
const Slot = util.Slot;

pub const Core = coretypes.X86CoreType;

const has_avx512f = std.Target.x86.featureSetHas(sysinfo.target.cpu.features, .avx512f);

/// `y += alpha * x` over a whole block, on vectors as wide as the core's
/// widest: what the reference's intrinsics block does (real slots only).
fn VectorBlock(comptime slot: Slot) type {
    const T = kernel.form(slot).F;
    const width = (if (has_avx512f) 64 else 32) / @sizeOf(T);
    const V = @Vector(width, T);
    return struct {
        pub fn run(n1: BlasLong, alpha: T, x: [*]T, y: [*]T) void {
            const av: V = @splat(alpha);
            const len: usize = @intCast(n1);
            var i: usize = 0;
            while (i + width <= len) : (i += width) {
                const xv: V = x[i..][0..width].*;
                const yv: V = y[i..][0..width].*;
                y[i..][0..width].* = fma.madd(fma.fuse, V, av, xv, yv);
            }
        }
    };
}

const K = uk.Kernels(Core, kernel, @import("asm/__root__.zig"), @import("uk/__root__.zig"), VectorBlock);

pub const hasSlot = K.hasSlot;
pub const has = K.has;
pub const exportCore = K.exportCore;
pub const install = K.install;
