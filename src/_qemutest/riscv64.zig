//! RISC-V 64 boot for QEMU's `virt` machine (`-bios none -kernel`, M-mode, one
//! hart, floating point and vector state on) and RISC-V semihosting for console and exit (`-semihosting`).

comptime {
    asm (
        \\.section .text.boot
        \\.globl _start
        \\_start:
        \\  csrr t0, mstatus
        \\  li t1, 0x6600
        \\  or t0, t0, t1
        \\  csrw mstatus, t0
        \\  la sp, __stack_top
        \\  la t0, __bss_start
        \\  la t1, __bss_end
        \\1:
        \\  bgeu t0, t1, 2f
        \\  sd zero, 0(t0)
        \\  addi t0, t0, 8
        \\  j 1b
        \\2:
        \\  call qemutest_main
        \\3:
        \\  wfi
        \\  j 3b
    );
}

const SYS_WRITEC = 0x03;
const SYS_EXIT_EXTENDED = 0x20;
const ADP_Stopped_ApplicationExit = 0x20026;

inline fn semihost(op: usize, arg: usize) void {
    asm volatile (
        \\.option push
        \\.option norvc
        \\slli zero, zero, 0x1f
        \\ebreak
        \\srai zero, zero, 7
        \\.option pop
        :
        : [op] "{a0}" (op),
          [arg] "{a1}" (arg),
        : .{ .memory = true });
}

pub fn writeByte(c: u8) void {
    var ch = c;
    semihost(SYS_WRITEC, @intFromPtr(&ch));
}

pub fn exit(code: usize) noreturn {
    var block = [2]usize{ ADP_Stopped_ApplicationExit, code };
    semihost(SYS_EXIT_EXTENDED, @intFromPtr(&block));
    while (true) {}
}
