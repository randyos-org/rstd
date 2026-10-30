//! RISC-V CPU detection, mirroring `driver/others/dynamic_riscv64.c`.
//!
//! The reference invokes the `riscv_hwprobe` syscall directly rather than
//! through glibc, falls back to `AT_HWCAP`'s 'V' bit on kernels without it,
//! then checks the vector unit really is RVV 1.0 (T-Head's 0.7.1 also sets the
//! bit) and reads `vlenb` to tell VLEN 128 from 256 and up.

const std = @import("std");
const coretypes = @import("../_coretypes.zig");
const select = @import("../_select.zig");
const sysinfo = @import("builtin");
const auxv = @import("../_auxv.zig");

pub const CoreType = coretypes.RiscvCoreType;
pub const Decision = select.Decision(CoreType);

/// `AT_HWCAP` reports base ISA extensions as letter bits: bit 0 is 'A',
/// so 'V' is bit 21.
pub const HWCAP_ISA_V: usize = 1 << ('V' - 'A');

const RISCV_HWPROBE_KEY_IMA_EXT_0: i64 = 4;
const RISCV_HWPROBE_IMA_V: u64 = 1 << 2;

const Pair = extern struct { key: i64, value: u64 };

/// The syscall number is stable across RISC-V Linux; the reference
/// hardcodes it for the same reason we do: `__NR_arch_specific_syscall` (244,
/// `include/uapi/asm-generic/unistd.h`) + 14.
const NR_riscv_hwprobe: usize = 258;

fn hwprobe(pairs: []Pair, flags: usize) bool {
    if (sysinfo.target.os.tag != .linux) return false;
    const rc = std.os.linux.syscall5(
        @fromBackingInt(@intCast(NR_riscv_hwprobe)),
        @intFromPtr(pairs.ptr),
        pairs.len,
        0, // cpu set size
        0, // cpu set pointer
        flags,
    );
    return @as(isize, @bitCast(rc)) == 0;
}

/// `vlenb` (CSR 0xc22), in bytes. Only legal where a vector unit exists.
fn readVlenb() u32 {
    return asm volatile ("csrr %[out], 0xc22"
        : [out] "=r" (-> u32),
    );
}

/// Whether the vector unit takes RVV 1.0's `vsetvli ... ta, ma`: 0.7.1 parts
/// set `vill` in `vtype`, which reads back negative.
fn isRvv100() bool {
    const r = asm volatile (
        \\.option push
        \\.option arch, +v
        \\vsetvli t0, x0, e8, m1, ta, ma
        \\csrr %[out], vtype
        \\slt %[out], x0, %[out]
        \\.option pop
        : [out] "=r" (-> usize),
        :
        : .{ .x5 = true });
    return r != 0;
}

fn csr(comptime name: []const u8) usize {
    return asm volatile ("csrr %[out], " ++ name
        : [out] "=r" (-> usize),
    );
}

/// With no OS to ask, in M-mode (as the rig runs): `misa`'s V (the hwcap's
/// letter bit) and the unit on in `mstatus.VS`, before any vector instruction.
fn machineCaps() Caps {
    if (csr("misa") & HWCAP_ISA_V == 0 or csr("mstatus") & 0x600 == 0) return .{};
    if (!isRvv100()) return .{};
    return .{ .vector = true, .vlenb = readVlenb() };
}

pub const Caps = struct {
    /// RVV 1.0 present and usable
    vector: bool = false,
    /// bytes per vector register, 0 if unknown
    vlenb: u32 = 0,
};

pub fn probe() Caps {
    if (comptime !sysinfo.target.cpu.arch.isRiscv()) return .{};
    switch (comptime sysinfo.target.os.tag) {
        .linux => {},
        .freestanding => return machineCaps(),
        else => return .{},
    }
    var pairs = [_]Pair{.{ .key = RISCV_HWPROBE_KEY_IMA_EXT_0, .value = 0 }};
    if (hwprobe(&pairs, 0)) {
        if (pairs[0].value & RISCV_HWPROBE_IMA_V == 0) return .{};
    } else {
        // older kernels: the hwcap letter bit, then the 0.7.1 check
        if (auxv.hwcap() & HWCAP_ISA_V == 0) return .{};
        if (!isRvv100()) return .{};
    }
    return .{ .vector = true, .vlenb = readVlenb() };
}

/// Whether the kernels of `core` run: the vector ones want the vector unit and
/// the VLEN they are named for.
pub fn runnable(core: CoreType, c: Caps) bool {
    return switch (core) {
        .riscv64_generic, .c910v, .u74 => true,
        .riscv64_zvl128b => c.vector and c.vlenb >= 16,
        .riscv64_zvl256b => c.vector and c.vlenb >= 32,
        .x280 => c.vector and c.vlenb >= 64,
    };
}

/// VLEN 256 and up is ZVL256B, 128 ZVL128B, anything else generic.
pub fn decide(c: Caps) Decision {
    if (!c.vector or c.vlenb < 16) return .{ .core = .riscv64_generic };
    if (c.vlenb >= 32) return .{ .core = .riscv64_zvl256b };
    return .{ .core = .riscv64_zvl128b };
}

/// Whether the kernels of `core` run on this CPU.
pub fn supportsCore(core: CoreType) bool {
    return runnable(core, probe());
}

pub fn detect() CoreType {
    return decide(probe()).core;
}

pub fn detectNoted() Decision {
    return decide(probe());
}

test "decide reads the vector length" {
    try std.testing.expectEqual(CoreType.riscv64_generic, decide(.{}).core);
    try std.testing.expectEqual(CoreType.riscv64_generic, decide(.{ .vector = true, .vlenb = 0 }).core);
    try std.testing.expectEqual(CoreType.riscv64_generic, decide(.{ .vector = false, .vlenb = 32 }).core);
    try std.testing.expectEqual(CoreType.riscv64_zvl128b, decide(.{ .vector = true, .vlenb = 16 }).core);
    try std.testing.expectEqual(CoreType.riscv64_zvl256b, decide(.{ .vector = true, .vlenb = 32 }).core);
    try std.testing.expectEqual(CoreType.riscv64_zvl256b, decide(.{ .vector = true, .vlenb = 64 }).core);
}

test "the vector cores need their VLEN" {
    const v128: Caps = .{ .vector = true, .vlenb = 16 };
    try std.testing.expect(runnable(.riscv64_zvl128b, v128));
    try std.testing.expect(!runnable(.riscv64_zvl256b, v128));
    try std.testing.expect(!runnable(.x280, v128));
    try std.testing.expect(runnable(.x280, .{ .vector = true, .vlenb = 64 }));
    try std.testing.expect(runnable(.riscv64_generic, .{}));
    try std.testing.expect(!runnable(.riscv64_zvl128b, .{}));
}

test "detect runs on RISC-V targets" {
    if (comptime !sysinfo.target.cpu.arch.isRiscv()) return error.SkipZigTest;
    try std.testing.expect(supportsCore(detect()));
}
