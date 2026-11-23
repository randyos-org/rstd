//! The PowerPC and s390x rigs.

const std = @import("std");
const machines = @import("../machines.zig");

const Machine = machines.Machine;

pub const list = [_]Machine{
    .{
        .name = "powerpc64le",
        // LLVM PowerPC FastISel (debug) uses FPR instructions on VSX registers (issue #186)
        .optimize = .safe,
        .query = .{ .cpu_arch = .powerpc64le, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.powerpc.cpu.pwr9 } },
        .qemu = "qemu-system-ppc64",
        .args = &.{ "-M", "pseries,x-vof=on", "-cpu", "power9" },
        .memory = "512M",
        .boot = .{
            .kernel = .{
                .linker_script = "test/qemu/ppc64.ld",
                // the firmware loader's 0x400000 shift, see the linker script
                .stderr_note = "kernel_addr changed",
            },
        },
    },
    .{
        .name = "powerpc64le-power10",
        // LLVM PowerPC FastISel (debug) uses FPR instructions on VSX registers (issue #186)
        .optimize = .safe,
        .query = .{ .cpu_arch = .powerpc64le, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.powerpc.cpu.pwr10 } },
        .qemu = "qemu-system-ppc64",
        .args = &.{ "-M", "pseries,x-vof=on", "-cpu", "power10" },
        .memory = "512M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/ppc64.ld",
            .stderr_note = "kernel_addr changed",
        } },
    },
    .{
        // big-endian, ELFv2 as Zig emits it; POWER10 runs every 64-bit core
        .name = "powerpc64",
        // LLVM PowerPC FastISel (debug) uses FPR instructions on VSX registers (issue #186)
        .optimize = .safe,
        .query = .{ .cpu_arch = .powerpc64, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.powerpc.cpu.pwr10 } },
        .qemu = "qemu-system-ppc64",
        .args = &.{ "-M", "pseries,x-vof=on", "-cpu", "power10" },
        .memory = "512M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/ppc64.ld",
            .stderr_note = "kernel_addr changed",
        } },
    },
    .{
        // a G4 under OpenBIOS: the ppcg4 core
        .name = "powerpc",
        .query = .{ .cpu_arch = .powerpc, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.powerpc.cpu.g4 } },
        .qemu = "qemu-system-ppc",
        .args = &.{ "-M", "mac99", "-cpu", "G4" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/ppc32.ld",
        } },
    },
    .{
        // RandyOS's iBook G3 target (a 750FX, no AltiVec) on `mac99`: soft
        // float, so no kernel objects; the portable kernels run. A 750, as
        // OpenBIOS hangs on QEMU's `750fx` (issue #225); the same ISA
        .name = "powerpc-750",
        .query = .{ .cpu_arch = .powerpc, .os_tag = .freestanding, .abi = .eabi, .cpu_model = .{ .explicit = &std.Target.powerpc.cpu.@"750" } },
        .qemu = "qemu-system-ppc",
        .args = &.{ "-M", "mac99", "-cpu", "750" },
        .memory = "256M",
        .blas_cores = .none,
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/ppc32.ld",
        } },
    },
    .{
        // the 750 with hard float: every 32-bit core linked, none runnable
        .name = "powerpc-750-hf",
        .query = .{ .cpu_arch = .powerpc, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.powerpc.cpu.@"750" } },
        .qemu = "qemu-system-ppc",
        .args = &.{ "-M", "mac99", "-cpu", "750" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/ppc32.ld",
        } },
    },
    .{
        // the G4 build in a POWER10's 32-bit mode under VOF: the server cores
        // through power6, built 32-bit
        .name = "powerpc-pseries",
        .query = .{ .cpu_arch = .powerpc, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.powerpc.cpu.g4 } },
        .qemu = "qemu-system-ppc64",
        .args = &.{ "-M", "pseries,x-vof=on", "-cpu", "power10" },
        .memory = "512M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/ppc64.ld",
            .stderr_note = "kernel_addr changed",
        } },
    },
    .{
        // a 32-bit VSX build on a POWER9 under VOF: the power8 and power9
        // cores' `@Vector` kernels as 32-bit code. The build drops `64bit`
        // and `fpcvt`, with which LLVM fails on compiler-rt in 32-bit mode;
        // ISA 3.0 fails it too, so power10's kernels cannot link (issue #187).
        .name = "powerpc-power9",
        .query = .{
            .cpu_arch = .powerpc,
            .os_tag = .freestanding,
            .abi = .none,
            .cpu_model = .{ .explicit = &std.Target.powerpc.cpu.pwr8 },
            .cpu_features_sub = std.Target.powerpc.featureSet(&.{ .@"64bit", .fpcvt }),
        },
        .qemu = "qemu-system-ppc64",
        .args = &.{ "-M", "pseries,x-vof=on", "-cpu", "power9" },
        .memory = "512M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/ppc64.ld",
            .stderr_note = "kernel_addr changed",
        } },
    },
    .{
        // Mach-O on PowerPC, the one big-endian 32-bit Mach-O target qemu can
        // boot (OpenBIOS on `mac99`); covers the object format the other rigs
        // do not.
        .name = "powerpc-macho",
        .query = .{ .cpu_arch = .powerpc, .os_tag = .freestanding, .abi = .none, .ofmt = .macho, .cpu_model = .{ .explicit = &std.Target.powerpc.cpu.g4 } },
        .qemu = "qemu-system-ppc",
        .args = &.{ "-M", "mac99", "-cpu", "G4" },
        .memory = "256M",
        .stub = "needs the powerpc-macos Zig port, a linker script, and a check that OpenBIOS boots the image",
        .boot = .{ .kernel = .{
            .linker_script = "",
        } },
    },
    .{
        .name = "s390x",
        .query = .{ .cpu_arch = .s390x, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.s390x.cpu.z14 } },
        .qemu = "qemu-system-s390x",
        .args = &.{ "-M", "s390-ccw-virtio", "-cpu", "max" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/s390x.ld",
        } },
    },
};
