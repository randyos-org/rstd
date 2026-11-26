//! Which kernel registry serves an architecture's family
//! (`coretypes.familyOf`). Only the registries and `cpu/coretypes.zig`:
//! the per-core object root (`_kernel.zig`) is compiled on its own and
//! cannot reach the rest of `cpu/`, where `cpu.GetCpuArch(arch).kernel`
//! makes the same choice for the library.

const std = @import("std");
const coretypes = @import("../cpu/coretypes.zig");

pub fn ForArch(comptime arch: std.Target.Cpu.Arch) type {
    return ForFamily(coretypes.familyOf(arch));
}

/// The registry of `family`.
pub fn ForFamily(comptime family: coretypes.Family) type {
    return switch (family) {
        .x86_64 => @import("x86_64.zig").Registry,
        .x86 => @import("x86.zig").Registry,
        .arm64 => @import("arm64.zig").Registry,
        .arm => @import("arm.zig").Registry,
        .power => @import("power.zig").Registry,
        .riscv => @import("riscv.zig").Registry,
        .s390x => @import("s390x.zig").Registry,
        .mips => @import("mips.zig").Registry,
        .loongarch => @import("loongarch.zig").Registry,
        .sparc => @import("sparc.zig").Registry,
        .wasm => @import("wasm.zig").Registry,
        .generic => @import("generic.zig").Registry,
    };
}
