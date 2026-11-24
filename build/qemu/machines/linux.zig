//! The Linux user-mode rigs.

const std = @import("std");
const machines = @import("../machines.zig");

const Machine = machines.Machine;
const sve_parts = machines.sve_parts;
const tsan_parts = machines.tsan_parts;

pub const list = [_]Machine{
    .{
        .name = "linux-aarch64",
        .query = .{ .cpu_arch = .aarch64, .os_tag = .linux },
        .qemu = "qemu-system-aarch64",
        // four CPUs, for the tests that race threads
        .args = &.{ "-M", "virt", "-cpu", "max", "-smp", "4" },
        .hvf = &.{ "-accel", "hvf", "-M", "virt", "-cpu", "host", "-smp", "4" },
        .hvf_core = "vortexm4",
        .memory = "1G",
        .tsan_memory = "8G",
        .tsan_parts = &tsan_parts,
        .boot = .{ .linux = .{ .nix = "aarch64-linux", .console = "ttyAMA0" } },
    },
    .{
        // the SVE and SME cores, linked only for an SVE target; TCG, as `aarch64-sve`
        .name = "linux-aarch64-sve",
        .query = .{ .cpu_arch = .aarch64, .os_tag = .linux, .cpu_model = .{ .explicit = &std.Target.aarch64.cpu.neoverse_v1 } },
        .qemu = "qemu-system-aarch64",
        .args = &.{ "-M", "virt", "-cpu", "max,pauth-impdef=on", "-smp", "4" },
        .core = "armv9sme",
        .memory = "1G",
        .split = &sve_parts,
        .boot = .{ .linux = .{ .nix = "aarch64-linux", .console = "ttyAMA0" } },
    },
    .{
        .name = "linux-x86_64",
        .query = .{ .cpu_arch = .x86_64, .os_tag = .linux },
        .qemu = "qemu-system-x86_64",
        .args = &.{ "-M", "pc", "-cpu", "max", "-smp", "4" },
        .hvf = &.{ "-accel", "hvf", "-M", "pc", "-cpu", "host", "-smp", "4" },
        .kvm = &.{ "-accel", "kvm", "-M", "pc", "-cpu", "host", "-smp", "4" },
        // the initramfs and its unpacked copy: the test, with every core's
        // kernels and their debug info, is about 700 MB, and `/` (where tests
        // write files) gets half the memory
        .memory = "3G",
        .tsan_memory = "12G",
        .tsan_parts = &tsan_parts,
        .boot = .{ .linux = .{ .nix = "x86_64-linux", .console = "ttyS0" } },
    },
    .{
        .name = "linux-arm",
        .query = .{ .cpu_arch = .arm, .os_tag = .linux },
        .qemu = "qemu-system-arm",
        .args = &.{ "-M", "virt", "-cpu", "cortex-a15", "-smp", "4" },
        .memory = "1G",
        .stub = "no armv7l-linux kernel in nixpkgs' binary cache; CI's qemu-user arm-linux leg covers it (issue #656)",
        .boot = .{ .linux = .{ .nix = "armv7l-linux", .console = "ttyAMA0" } },
    },
    .{
        .name = "linux-riscv64",
        .query = .{ .cpu_arch = .riscv64, .os_tag = .linux },
        .qemu = "qemu-system-riscv64",
        // the vector extension at the widest linked core's VLEN, as `riscv64`
        .args = &.{ "-M", "virt", "-cpu", "rv64,v=true,vlen=512", "-smp", "4" },
        .memory = "1G",
        .stub = "no riscv64-linux kernel in nixpkgs' binary cache; CI's qemu-user riscv64-linux leg covers it (issue #656)",
        .boot = .{ .linux = .{ .nix = "riscv64-linux", .console = "ttyS0" } },
    },
    .{
        .name = "linux-powerpc64le",
        // pwr8, `f128` by library calls: clear of the FastISel bug the other
        // POWER rigs build `.safe` for (issue #186)
        .query = .{ .cpu_arch = .powerpc64le, .os_tag = .linux },
        .qemu = "qemu-system-ppc64",
        // detection sees a POWER10 (the baseline build links up to POWER9's
        // kernels); the caps TCG lacks are off, so qemu warns of none
        .args = &.{ "-M", "pseries,cap-cfpc=broken,cap-sbbc=broken,cap-ibs=broken,cap-ccf-assist=off", "-cpu", "power10", "-smp", "4" },
        .memory = "1G",
        .stub = "no powerpc64le-linux kernel in nixpkgs' binary cache; CI's qemu-user powerpc64le-linux legs cover it (issue #656)",
        .boot = .{ .linux = .{ .nix = "powerpc64le-linux", .console = "hvc0" } },
    },
    .{
        // a POWER10 build, so power10's kernels link and run
        .name = "linux-powerpc64le-power10",
        // native `f128` from pwr9: FastISel's bug, as the freestanding rigs (issue #186)
        .optimize = .safe,
        .query = .{ .cpu_arch = .powerpc64le, .os_tag = .linux, .cpu_model = .{ .explicit = &std.Target.powerpc.cpu.pwr10 } },
        .qemu = "qemu-system-ppc64",
        .args = &.{ "-M", "pseries,cap-cfpc=broken,cap-sbbc=broken,cap-ibs=broken,cap-ccf-assist=off", "-cpu", "power10", "-smp", "4" },
        .core = "power10",
        .memory = "1G",
        .stub = "no powerpc64le-linux kernel in nixpkgs' binary cache; CI's qemu-user powerpc64le-linux legs cover it (issue #656)",
        .boot = .{ .linux = .{ .nix = "powerpc64le-linux", .console = "hvc0" } },
    },
    .{
        .name = "linux-i386",
        .query = .{ .cpu_arch = .x86, .os_tag = .linux },
        .qemu = "qemu-system-i386",
        .args = &.{ "-M", "pc", "-cpu", "max", "-smp", "4" },
        .memory = "1G",
        // boots and runs; filtered runs pass
        .stub = "TCG's fldt/fstpt misplace the significand bytes past a page boundary, which i386's 4-aligned f80 arrays cross (issue #221)",
        .boot = .{ .linux = .{ .nix = "i686-linux", .console = "ttyS0" } },
    },
};
