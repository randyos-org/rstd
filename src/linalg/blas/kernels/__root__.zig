//! Which kernel registry serves an architecture. Import-free (only the
//! registries, which stay inside `linalg/blas/`) because the per-core object
//! root (`_kernel.zig`) is compiled on its own and cannot reach `cpu/`, which
//! is where the same choice is made for the rest of the library
//! (`cpu.GetCpuArch(arch).kernel`); keep the two in step.

const std = @import("std");

pub fn ForArch(comptime arch: std.Target.Cpu.Arch) type {
    return switch (arch) {
        .x86, .x86_64 => @import("x86.zig"),
        .aarch64, .aarch64_be => @import("arm64.zig"),
        .arm, .armeb, .thumb, .thumbeb => @import("arm.zig"),
        .powerpc, .powerpcle, .powerpc64, .powerpc64le => @import("power.zig"),
        .riscv32, .riscv32be, .riscv64, .riscv64be => @import("riscv.zig"),
        .s390x => @import("s390x.zig"),
        .mips, .mipsel, .mips64, .mips64el => @import("mips.zig"),
        .loongarch32, .loongarch64 => @import("loongarch.zig"),
        .sparc, .sparc64 => @import("sparc.zig"),
        .wasm32, .wasm64 => @import("wasm.zig"),
        else => @import("generic.zig"),
    };
}
