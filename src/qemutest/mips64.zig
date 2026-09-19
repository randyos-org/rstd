//! MIPS64 (release 6) boot for QEMU's `malta` machine (`-cpu I6400 -kernel`)
//! and UHI semihosting for console and exit (`-semihosting`).
//!
//! Release 6 handles unaligned accesses in hardware and uses compact
//! branches (no delay slots), which keeps this boot code simple. `_start`
//! turns on the FPU (with 64-bit registers, `Status.FR`) and MSA.

comptime {
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
