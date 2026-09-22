//! MIPS32 boot for QEMU's `malta` machine (`-kernel`) and UHI semihosting for
//! console and exit (`-semihosting`), for release 2 (`-cpu 24Kf`) and release
//! 5 (`-cpu P5600`). Release 5 has MSA, which `_start` turns on together with
//! the 64-bit floating-point registers it needs (`Status.FR`); release 2 has
//! the 32-bit FPU only.

const builtin = @import("builtin");

comptime {
    if (builtin.cpu.has(.mips, .mips32r5)) {
        asm (
            \\.set noreorder
            \\.set noat
            \\.section .text.boot
            \\.globl _start
            \\.globl __start
            \\_start:
            \\__start:
            \\  mfc0 $t0, $12
            \\  lui $t1, 0x2400
            \\  or $t0, $t0, $t1
            \\  mtc0 $t0, $12
            \\  ehb
            \\  mfc0 $t0, $16, 5
            \\  lui $t1, 0x0800
            \\  or $t0, $t0, $t1
            \\  mtc0 $t0, $16, 5
            \\  ehb
            \\  la $sp, __stack_top
            \\  la $t0, __bss_start
            \\  la $t1, __bss_end
            \\1:
            \\  sltu $t2, $t0, $t1
            \\  beqz $t2, 2f
            \\  nop
            \\  sw $zero, 0($t0)
            \\  b 1b
            \\  addiu $t0, $t0, 4
            \\2:
            \\  jal qemutest_main
            \\  nop
            \\3:
            \\  wait
            \\  b 3b
            \\  nop
        );
    } else {
        asm (
            \\.set noreorder
            \\.set noat
            \\.section .text.boot
            \\.globl _start
            \\.globl __start
            \\_start:
            \\__start:
            \\  mfc0 $t0, $12
            \\  lui $t1, 0x2000
            \\  or $t0, $t0, $t1
            \\  mtc0 $t0, $12
            \\  ehb
            \\  la $sp, __stack_top
            \\  la $t0, __bss_start
            \\  la $t1, __bss_end
            \\1:
            \\  sltu $t2, $t0, $t1
            \\  beqz $t2, 2f
            \\  nop
            \\  sw $zero, 0($t0)
            \\  b 1b
            \\  addiu $t0, $t0, 4
            \\2:
            \\  jal qemutest_main
            \\  nop
            \\3:
            \\  wait
            \\  b 3b
            \\  nop
        );
    }
}

const UHI_exit = 1;
const UHI_write = 5;

inline fn uhi3(op: usize, a0: usize, a1: usize, a2: usize) void {
    asm volatile ("sdbbp 1"
        :
        : [op] "{$25}" (op),
          [a0] "{$4}" (a0),
          [a1] "{$5}" (a1),
          [a2] "{$6}" (a2),
        : .{ .memory = true, .r2 = true, .r3 = true });
}

pub fn writeByte(c: u8) void {
    var ch = c;
    uhi3(UHI_write, 1, @intFromPtr(&ch), 1);
}

pub fn exit(code: usize) noreturn {
    uhi3(UHI_exit, code, 0, 0);
    while (true) {}
}
