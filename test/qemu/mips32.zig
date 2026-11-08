//! MIPS32 boot for `malta` with UHI semihosting console and exit, release 2 (24Kf, 32-bit FPU) and 5 (P5600), where `_start` turns on
//! MSA and the 64-bit FP registers it needs (`Status.FR`).

const sysinfo = @import("builtin");
const uhi = @import("uhi.zig");

comptime {
    if (sysinfo.target.cpu.has(.mips, .mips32r5)) {
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

pub const writeByte = uhi.writeByte;
pub const exit = uhi.exit;
