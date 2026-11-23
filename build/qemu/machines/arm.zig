//! The aarch64 and 32-bit arm rigs.

const std = @import("std");
const machines = @import("../machines.zig");

const Machine = machines.Machine;
const sve_parts = machines.sve_parts;

pub const list = [_]Machine{
    .{
        .name = "aarch64",
        .query = .{ .cpu_arch = .aarch64, .os_tag = .freestanding, .cpu_model = .{ .explicit = &std.Target.aarch64.cpu.cortex_a72 } },
        .qemu = "qemu-system-aarch64",
        .args = &.{ "-M", "virt", "-cpu", "cortex-a72" },
        // no semihosting: the boot code falls back to the PL011 and PSCI
        .hvf = &.{ "-accel", "hvf", "-M", "virt", "-cpu", "host" },
        .kvm = &.{ "-accel", "kvm", "-M", "virt,gic-version=host", "-cpu", "host" },
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/aarch64.ld",
        } },
    },
    .{
        // the threaded tests on four cores, started by PSCI
        .name = "aarch64-smp",
        .query = .{ .cpu_arch = .aarch64, .os_tag = .freestanding, .cpu_model = .{ .explicit = &std.Target.aarch64.cpu.cortex_a72 } },
        .qemu = "qemu-system-aarch64",
        .args = &.{ "-M", "virt", "-cpu", "cortex-a72" },
        .hvf = &.{ "-accel", "hvf", "-M", "virt", "-cpu", "host" },
        .kvm = &.{ "-accel", "kvm", "-M", "virt,gic-version=host", "-cpu", "host" },
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/aarch64.ld",
        } },
        .cpus = 4,
    },
    .{
        // TCG on a Mac too, which has no SVE
        .name = "aarch64-sve",
        .query = .{ .cpu_arch = .aarch64, .os_tag = .freestanding, .cpu_model = .{ .explicit = &std.Target.aarch64.cpu.neoverse_v1 } },
        .qemu = "qemu-system-aarch64",
        // max has SVE2, SME2 and BF16, so every core but Apple's runs
        .args = &.{ "-M", "virt", "-cpu", "max,pauth-impdef=on" },
        .split = &sve_parts,
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/aarch64.ld",
        } },
    },
    .{
        .name = "arm",
        // the A15's FPU and divide: Zig's model has neither (LLVM adds them
        // from its target parser), so the bare model is soft-float
        .query = .{
            .cpu_arch = .arm,
            .os_tag = .freestanding,
            .abi = .eabihf,
            .cpu_model = .{ .explicit = &std.Target.arm.cpu.cortex_a15 },
            .cpu_features_add = std.Target.arm.featureSet(&.{ .vfp4, .neon, .d32, .hwdiv, .hwdiv_arm }),
        },
        .qemu = "qemu-system-arm",
        .args = &.{ "-M", "virt", "-cpu", "cortex-a15" },
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/arm.ld",
        } },
    },
    .{
        // the soft-float ABI over the A15's FPU (softfp): floats in integer
        // registers, the kernels' `sf` variant; the features pinned as in
        // `arm`, since Zig 0.18's model has none and 0.17's had them
        .name = "arm-softfp",
        .query = .{
            .cpu_arch = .arm,
            .os_tag = .freestanding,
            .abi = .eabi,
            .cpu_model = .{ .explicit = &std.Target.arm.cpu.cortex_a15 },
            .cpu_features_add = std.Target.arm.featureSet(&.{ .vfp4, .neon, .d32, .hwdiv, .hwdiv_arm }),
        },
        .qemu = "qemu-system-arm",
        .args = &.{ "-M", "virt", "-cpu", "cortex-a15" },
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/arm.ld",
        } },
    },
};
