//! The UEFI rigs.

const std = @import("std");
const machines = @import("../machines.zig");

const Machine = machines.Machine;
const sve_parts = machines.sve_parts;

pub const list = [_]Machine{
    .{
        // RandyOS's loader's platforms: the tests as a UEFI application under
        // QEMU's OVMF and AAVMF, on rstd's UEFI `std.Io`
        .name = "x86_64-uefi",
        .query = .{ .cpu_arch = .x86_64, .os_tag = .uefi },
        .qemu = "qemu-system-x86_64",
        .args = &.{ "-M", "q35", "-cpu", "max" },
        .hvf = &.{ "-accel", "hvf", "-M", "q35", "-cpu", "host" },
        .kvm = &.{ "-accel", "kvm", "-M", "q35", "-cpu", "host" },
        .memory = "1G",
        .boot = .{ .kernel = .{
            .linker_script = "",
            .llvm = true,
            .firmware = "edk2-x86_64-code.fd",
        } },
    },
    .{
        .name = "aarch64-uefi",
        .query = .{ .cpu_arch = .aarch64, .os_tag = .uefi },
        .qemu = "qemu-system-aarch64",
        .args = &.{ "-M", "virt", "-cpu", "cortex-a72" },
        .hvf = &.{ "-accel", "hvf", "-M", "virt", "-cpu", "host" },
        .kvm = &.{ "-accel", "kvm", "-M", "virt,gic-version=host", "-cpu", "host" },
        .memory = "1G",
        .boot = .{ .kernel = .{
            .linker_script = "",
            .firmware = "edk2-aarch64-code.fd",
        } },
    },
    .{
        // the SVE and SME cores, linked only for an SVE target; TCG, as `aarch64-sve`
        .name = "aarch64-uefi-sve",
        .query = .{ .cpu_arch = .aarch64, .os_tag = .uefi, .cpu_model = .{ .explicit = &std.Target.aarch64.cpu.neoverse_v1 } },
        .qemu = "qemu-system-aarch64",
        .args = &.{ "-M", "virt", "-cpu", "max,pauth-impdef=on" },
        .memory = "1G",
        .split = &sve_parts,
        .boot = .{ .kernel = .{
            .linker_script = "",
            .firmware = "edk2-aarch64-code.fd",
        } },
    },
};
