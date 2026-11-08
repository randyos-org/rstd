//! MIPS64 boot for `malta` with UHI semihosting console and exit, release 2 (Loongson 3A1000, 3A4000) and 6 (I6400). `_start` turns on
//! the FPU with 64-bit registers (`Status.FR`) and MSA where present (r6; `Config3.MSAP` on r2).

const sysinfo = @import("builtin");
const uhi = @import("uhi.zig");

comptime {
    if (sysinfo.target.cpu.has(.mips, .mips64r6)) {
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
            \\  dla $sp, __stack_top
            \\  dla $t0, __bss_start
            \\  dla $t1, __bss_end
            \\1:
            \\  sltu $t2, $t0, $t1
            \\  beqzc $t2, 2f
            \\  sd $zero, 0($t0)
            \\  daddiu $t0, $t0, 8
            \\  bc 1b
            \\2:
            \\  balc qemutest_main
            \\3:
            \\  wait
            \\  bc 3b
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
            \\  lui $t1, 0x2400
            \\  or $t0, $t0, $t1
            \\  mtc0 $t0, $12
            \\  ehb
            \\  mfc0 $t0, $16, 3
            \\  lui $t1, 0x1000
            \\  and $t0, $t0, $t1
            \\  beqz $t0, 4f
            \\  nop
            \\  mfc0 $t0, $16, 5
            \\  lui $t1, 0x0800
            \\  or $t0, $t0, $t1
            \\  mtc0 $t0, $16, 5
            \\  ehb
            \\4:
            \\  dla $sp, __stack_top
            \\  dla $t0, __bss_start
            \\  dla $t1, __bss_end
            \\1:
            \\  sltu $t2, $t0, $t1
            \\  beqz $t2, 2f
            \\  nop
            \\  sd $zero, 0($t0)
            \\  b 1b
            \\  daddiu $t0, $t0, 8
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
