//! x86_64 axpy kernels per core (`asm/`, or `../driver.zig` around a `uk/` block or `@Vector`s for AVX-512); tests in `../__root__.zig`.

const coretypes = @import("../../../cpu/coretypes.zig");
const util = @import("../../../util.zig");
const kernel = @import("../driver.zig");
const uk = @import("../../../uk_block.zig");

const Slot = util.Slot;

pub const Core = coretypes.X86CoreType;

const has_avx512f = util.x86Has(.avx512f);

/// `y += alpha * x` over a whole block, on vectors as wide as the core's
/// widest: what the reference's intrinsics block does (real slots only).
fn VectorBlock(comptime slot: Slot) type {
    const T = kernel.form(slot).F;
    return kernel.VectorBlock(T, (if (has_avx512f) 64 else 32) / @sizeOf(T), 1);
}

const K = uk.Kernels(Core, kernel, @import("asm/__root__.zig"), @import("uk/__root__.zig"), VectorBlock);

pub const hasSlot = K.hasSlot;
pub const has = K.has;
pub const exportCore = K.exportCore;
pub const install = K.install;
