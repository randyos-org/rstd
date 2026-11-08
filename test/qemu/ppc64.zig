//! 64-bit POWER boot, either endianness, for `pseries` under VOF; console `H_PUT_TERM_CHAR`, exit RTAS `power-off` (the summary carries the result).
//! `_start` flips MSR.LE by the usual `tdi` trampoline when little-endian, and sets MSR.SF: VOF enters in 32-bit mode, where carries come
//! from the low 32 bits. Zig emits ELFv2 either way, so no function descriptor. Load addresses: `ppc64.ld`.

const std = @import("std");

comptime {
    asm (
        \\.section .text.boot, "ax"
        \\.globl _start
        \\_start:
        \\  tdi 0, 0, 0x48
        \\  b . + 44
        \\  .byte 0x7d, 0x60, 0x00, 0xa6
        \\  .byte 0x69, 0x6b, 0x00, 0x01
        \\  .byte 0x39, 0x40, 0x00, 0x00
        \\  .byte 0x7d, 0x41, 0x01, 0x64
        \\  .byte 0x42, 0x9f, 0x00, 0x05
        \\  .byte 0x7d, 0x48, 0x02, 0xa6
        \\  .byte 0x39, 0x4a, 0x00, 0x14
        \\  .byte 0x7d, 0x5a, 0x03, 0xa6
        \\  .byte 0x7d, 0x7b, 0x03, 0xa6
        \\  .byte 0x4c, 0x00, 0x00, 0x24
        \\  bl 0f
        \\0:
        \\  mflr 12
        \\  addis 2, 12, .TOC.-0b@ha
        \\  addi 2, 2, .TOC.-0b@l
        \\  mfmsr 3
        \\  li 4, 1
        \\  sldi 4, 4, 63
        \\  oris 4, 4, 0x0280
        \\  ori 4, 4, 0x2000
        \\  or 3, 3, 4
        \\  mtmsrd 3
        \\  isync
        \\  lis 1, __stack_top@highest
        \\  ori 1, 1, __stack_top@higher
        \\  sldi 1, 1, 32
        \\  oris 1, 1, __stack_top@h
        \\  ori 1, 1, __stack_top@l
        \\  li 0, 0
        \\  stdu 0, -64(1)
        \\  lis 6, __bss_start@highest
        \\  ori 6, 6, __bss_start@higher
        \\  sldi 6, 6, 32
        \\  oris 6, 6, __bss_start@h
        \\  ori 6, 6, __bss_start@l
        \\  lis 7, __bss_end@highest
        \\  ori 7, 7, __bss_end@higher
        \\  sldi 7, 7, 32
        \\  oris 7, 7, __bss_end@h
        \\  ori 7, 7, __bss_end@l
        \\1:
        \\  cmpld 6, 7
        \\  bge 2f
        \\  std 0, 0(6)
        \\  addi 6, 6, 8
        \\  b 1b
        \\2:
        \\  bl qemutest_main
        \\  nop
        \\3:
        \\  b 3b
    );
}

/// PAPR hypercall: token in r3, arguments r4-r7, status back in r3.
inline fn hcall(token: usize, a1: usize, a2: usize, a3: usize, a4: usize) usize {
    return asm volatile ("sc 1"
        : [ret] "={r3}" (-> usize),
        : [tok] "{r3}" (token),
          [a1] "{r4}" (a1),
          [a2] "{r5}" (a2),
          [a3] "{r6}" (a3),
          [a4] "{r7}" (a4),
        : .{ .memory = true, .r0 = true, .r8 = true, .r9 = true, .r10 = true, .r11 = true, .r12 = true, .cr0 = true });
}

const H_PUT_TERM_CHAR = 0x58;
/// Register of the machine's default vty (`/vdevice/vty@71000000`).
const vty_reg = 0x71000000;

pub fn writeByte(c: u8) void {
    // The character sits in the top byte of the first data register.
    _ = hcall(H_PUT_TERM_CHAR, vty_reg, 1, @as(usize, c) << 56, 0);
}

const H_RTAS = 0xf000;
/// QEMU's RTAS token for `power-off` (`RTAS_TOKEN_BASE` 0x2000 + 3).
const RTAS_POWER_OFF = 0x2003;

/// QEMU services RTAS directly through `H_RTAS`, so no firmware RTAS blob is
/// needed. The argument buffer holds big-endian `u32` cells.
pub fn exit(code: usize) noreturn {
    _ = code;
    // `power-off` takes two (ignored) arguments and returns one status.
    const be = std.mem.nativeToBig;
    var args = [6]u32{ be(u32, RTAS_POWER_OFF), be(u32, 2), be(u32, 1), 0, 0, 0 };
    _ = hcall(H_RTAS, @intFromPtr(&args), 0, 0, 0);
    while (true) {}
}
