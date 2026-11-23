//! The i386 and x86_64 rigs.

const std = @import("std");
const machines = @import("../machines.zig");

const Machine = machines.Machine;

pub const list = [_]Machine{
    .{
        .name = "i386",
        .query = .{ .cpu_arch = .x86, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.x86.cpu.pentium3 } },
        .qemu = "qemu-system-i386",
        .args = &.{ "-M", "pc", "-cpu", "pentium3" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/x86.ld",
        } },
    },
    .{
        // the SSE2-and-up cores, which a `pentium3` cannot run; a CPU past the
        // baseline (`pentium4`) also builds each core's object for its own ISA
        .name = "i386-max",
        .query = .{ .cpu_arch = .x86, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.x86.cpu.nehalem } },
        .qemu = "qemu-system-i386",
        .args = &.{ "-M", "pc", "-cpu", "max" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/x86.ld",
        } },
    },
    .{
        .name = "x86_64",
        .query = .{ .cpu_arch = .x86_64, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.x86.cpu.haswell } },
        .qemu = "qemu-system-x86_64",
        .args = &.{ "-M", "pc", "-cpu", "max" },
        .hvf = &.{ "-accel", "hvf", "-M", "pc", "-cpu", "host" },
        .kvm = &.{ "-accel", "kvm", "-M", "pc", "-cpu", "host" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/x86_64.ld",
            .llvm = true,
        } },
    },
    .{
        // the threaded tests on four cores, started by INIT-SIPI-SIPI
        .name = "x86_64-smp",
        .query = .{ .cpu_arch = .x86_64, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.x86.cpu.haswell } },
        .qemu = "qemu-system-x86_64",
        .args = &.{ "-M", "pc", "-cpu", "max" },
        .hvf = &.{ "-accel", "hvf", "-M", "pc", "-cpu", "host" },
        .kvm = &.{ "-accel", "kvm", "-M", "pc", "-cpu", "host" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/x86_64.ld",
            .llvm = true,
        } },
        .cpus = 4,
    },
    .{
        // The host CPU under Windows' hypervisor: TCG has no AVX-512, FMA4 or
        // XOP (issue #220), so their cores run only here, on a host that has
        // them. WHPX keeps only the model's features the host has; `host`
        // and `max` stop in SeaBIOS.
        .name = "x86_64-whpx",
        .query = .{ .cpu_arch = .x86_64, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.x86.cpu.haswell } },
        .qemu = "qemu-system-x86_64",
        .args = &.{ "-accel", "whpx", "-M", "pc", "-cpu", "SapphireRapids,+fma4,+xop" },
        .whpx = true,
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/x86_64.ld",
            .llvm = true,
        } },
    },
};
