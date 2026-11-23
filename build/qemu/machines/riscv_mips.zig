//! The riscv64, loongarch64 and MIPS rigs.

const std = @import("std");
const machines = @import("../machines.zig");

const Machine = machines.Machine;

pub const list = [_]Machine{
    .{
        .name = "riscv64",
        .query = .{ .cpu_arch = .riscv64, .os_tag = .freestanding, .cpu_model = .baseline },
        .qemu = "qemu-system-riscv64",
        // the vector extension, which the RVV kernels need, at the VLEN of
        // the widest linked core (x280's 512)
        .args = &.{ "-M", "virt", "-cpu", "rv64,v=true,vlen=512", "-smp", "1", "-bios", "none" },
        .boot = .{
            .kernel = .{
                .linker_script = "test/qemu/riscv64.ld",
                // The kernel sits at 0x80000000, out of reach of the default `medlow`.
                .code_model = .medany,
            },
        },
    },
    .{
        // `riscv64` at VLEN 128 and 256: `zvl128b`'s and `zvl256b`'s objects
        // on the hardware they are named for, the VLEN-agnostic `axpy` at both
        .name = "riscv64-vlen128",
        .query = .{ .cpu_arch = .riscv64, .os_tag = .freestanding, .cpu_model = .baseline },
        .qemu = "qemu-system-riscv64",
        .args = &.{ "-M", "virt", "-cpu", "rv64,v=true,vlen=128", "-smp", "1", "-bios", "none" },
        .boot = .{ .kernel = .{ .linker_script = "test/qemu/riscv64.ld", .code_model = .medany } },
    },
    .{
        .name = "riscv64-vlen256",
        .query = .{ .cpu_arch = .riscv64, .os_tag = .freestanding, .cpu_model = .baseline },
        .qemu = "qemu-system-riscv64",
        .args = &.{ "-M", "virt", "-cpu", "rv64,v=true,vlen=256", "-smp", "1", "-bios", "none" },
        .boot = .{ .kernel = .{ .linker_script = "test/qemu/riscv64.ld", .code_model = .medany } },
    },
    .{
        .name = "loongarch64",
        .query = .{ .cpu_arch = .loongarch64, .os_tag = .freestanding, .cpu_model = .baseline },
        .qemu = "qemu-system-loongarch64",
        .args = &.{ "-M", "virt", "-cpu", "la464" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/loongarch64.ld",
        } },
    },
    .{
        .name = "mips64r2el",
        .query = .{ .cpu_arch = .mips64el, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.mips.cpu.mips64r2 } },
        .qemu = "qemu-system-mips64el",
        // with MSA, so the Loongson 3A R4 core runs; the big-endian rig keeps
        // the 3A1000, a release 2 machine without it
        .args = &.{ "-M", "malta", "-cpu", "Loongson-3A4000" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/mips64.ld",
        } },
    },
    .{
        .name = "mips64r6el",
        .query = .{ .cpu_arch = .mips64el, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.mips.cpu.mips64r6 } },
        .qemu = "qemu-system-mips64el",
        .args = &.{ "-M", "malta", "-cpu", "I6400" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/mips64.ld",
        } },
    },
    .{
        .name = "mips32r2el",
        .query = .{ .cpu_arch = .mipsel, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.mips.cpu.mips32r2 } },
        .qemu = "qemu-system-mipsel",
        .args = &.{ "-M", "malta", "-cpu", "24Kf" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/mips32.ld",
        } },
    },
    .{
        .name = "mips32r5el",
        .query = .{ .cpu_arch = .mipsel, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.mips.cpu.p5600 } },
        .qemu = "qemu-system-mipsel",
        .args = &.{ "-M", "malta", "-cpu", "P5600" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/mips32.ld",
        } },
    },
    .{
        .name = "mips64r2",
        .query = .{ .cpu_arch = .mips64, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.mips.cpu.mips64r2 } },
        .qemu = "qemu-system-mips64",
        .args = &.{ "-M", "malta", "-cpu", "Loongson-3A1000" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/mips64.ld",
        } },
    },
    .{
        .name = "mips64r6",
        .query = .{ .cpu_arch = .mips64, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.mips.cpu.mips64r6 } },
        .qemu = "qemu-system-mips64",
        .args = &.{ "-M", "malta", "-cpu", "I6400" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/mips64.ld",
        } },
    },
    .{
        .name = "mips32r2",
        .query = .{ .cpu_arch = .mips, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.mips.cpu.mips32r2 } },
        .qemu = "qemu-system-mips",
        .args = &.{ "-M", "malta", "-cpu", "24Kf" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/mips32.ld",
        } },
    },
    .{
        .name = "mips32r5",
        .query = .{ .cpu_arch = .mips, .os_tag = .freestanding, .abi = .none, .cpu_model = .{ .explicit = &std.Target.mips.cpu.p5600 } },
        .qemu = "qemu-system-mips",
        .args = &.{ "-M", "malta", "-cpu", "P5600" },
        .memory = "256M",
        .boot = .{ .kernel = .{
            .linker_script = "test/qemu/mips32.ld",
        } },
    },
};
