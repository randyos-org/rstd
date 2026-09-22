//! MIPS / MIPS64 CPU detection, mirroring
//! `driver/others/dynamic_mips64.c`.
//!
//! The reference prefers Loongson's `cpucfg` instruction when
//! `AT_HWCAP` advertises it, and otherwise parses `/proc/cpuinfo` for the
//! model name -- the only detector in the tree that has to read a file.

const std = @import("std");
const coretypes = @import("_coretypes.zig");
const sysinfo = @import("builtin");
const auxv = @import("_auxv.zig");

pub const CoreType = coretypes.MipsCoreType;
pub const kernel = @import("../kernels/mips.zig");

pub const HWCAP_LOONGSON_CPUCFG: usize = 1 << 14;

const LOONGSON_CFG1: u32 = 0x1;
const LOONGSON_CONF1_MMI: u32 = 1 << 4;
const LOONGSON_CONF1_MSA1: u32 = 1 << 5;

pub fn cpucfg(word: u32) u32 {
    // Loongson's CPUCFG is encoded as a coprocessor-0 style instruction;
    // the reference emits it via `.word`. Kept behind the same hwcap gate
    // it uses, since it is not a baseline MIPS instruction.
    return asm volatile (
        \\.set push
        \\.set noreorder
        \\.insn
        \\.word 0xc8080118
        \\.set pop
        : [out] "={$8}" (-> u32),
        : [in] "{$4}" (word),
    );
}

/// Parses `/proc/cpuinfo` for the model name, the reference's fallback when
/// `cpucfg` is unavailable. This is the only detector in the tree that has
/// to read a file, so it goes straight to the raw syscalls rather than
/// through `std.fs` -- no allocator, no `Io` instance, and it runs at
/// library init where neither is necessarily set up yet.
fn detectFromCpuinfo() CoreType {
    if (sysinfo.os.tag != .linux) return .mips64_generic;

    const linux = std.os.linux;
    const rc = linux.open("/proc/cpuinfo", .{ .ACCMODE = .RDONLY }, 0);
    if (@as(isize, @bitCast(rc)) < 0) return .mips64_generic;
    const fd: linux.fd_t = @intCast(rc);
    defer _ = linux.close(fd);

    var buf: [2048]u8 = undefined;
    const got = linux.read(fd, &buf, buf.len);
    if (@as(isize, @bitCast(got)) <= 0) return .mips64_generic;
    const text = buf[0..got];

    if (std.mem.indexOf(u8, text, "Loongson-3A R4") != null) return .loongson3r4;
    if (std.mem.indexOf(u8, text, "Loongson-3A R3") != null) return .loongson3r3;
    if (std.mem.indexOf(u8, text, "Loongson-3") != null) return .loongson3r3;
    return .mips64_generic;
}

pub fn detect() CoreType {
    if (auxv.hwcap() & HWCAP_LOONGSON_CPUCFG == 0) return detectFromCpuinfo();

    const cfg = cpucfg(LOONGSON_CFG1);
    if (cfg & LOONGSON_CONF1_MSA1 != 0) return .loongson3r4;
    if (cfg & LOONGSON_CONF1_MMI != 0) return .loongson3r3;
    return .mips64_generic;
}

const Vtable = @import("../Vtable.zig");
const axpy = @import("../blas1/axpy/__root__.zig");
const blas = @import("../__root__.zig");
const gemm_params = @import("../blas3/gemm/params.zig");
pub fn vtableFromCoreType(core_type: CoreType) Vtable {
    var vt: Vtable = .{};
    axpy.install(&vt, core_type);
    gemm_params.installCore(&vt, (if (@sizeOf(usize) == 8) .mips64 else .mips32), core_type, .fromConfig(blas.config));
    return vt;
}

test "detect runs on MIPS targets" {
    if (comptime !sysinfo.cpu.arch.isMIPS()) return error.SkipZigTest;
    _ = detect();
}
