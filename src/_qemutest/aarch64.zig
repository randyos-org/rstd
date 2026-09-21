//! AArch64 boot for QEMU's `virt` machine (`-kernel`, EL1) and ARM
//! semihosting for console and exit (`-semihosting`).
//!
//! With the MMU off every data access is Device memory, which faults on the
//! unaligned accesses compiled code makes, so `_start` builds a one-block
//! identity map first: 0-1 GiB as device, 1-2 GiB (RAM) as normal cached.

comptime {
    asm (
        \\.section .text.boot
        \\.globl _start
        \\_start:
        \\  mrs x0, cpacr_el1
        \\  orr x0, x0, #(3 << 20)
        \\  msr cpacr_el1, x0
        \\  isb
        \\  ldr x0, =0x00ff
        \\  msr mair_el1, x0
        \\  ldr x0, =0x200803519
        \\  msr tcr_el1, x0
        \\  adrp x0, qemutest_l1
        \\  msr ttbr0_el1, x0
        \\  isb
        \\  mrs x0, sctlr_el1
        \\  orr x0, x0, #1
        \\  orr x0, x0, #(1 << 2)
        \\  orr x0, x0, #(1 << 12)
        \\  bic x0, x0, #(1 << 1)
        \\  msr sctlr_el1, x0
        \\  isb
        \\  ldr x0, =__stack_top
        \\  mov sp, x0
        \\  ldr x0, =__bss_start
        \\  ldr x1, =__bss_end
        \\1:
        \\  cmp x0, x1
        \\  b.hs 2f
        \\  str xzr, [x0], #8
        \\  b 1b
        \\2:
        \\  bl qemutest_main
        \\3:
        \\  wfe
        \\  b 3b
        \\.section .data.pgtbl, "aw"
        \\.balign 4096
        \\qemutest_l1:
        \\  .quad 0x00000405
        \\  .quad 0x40000701
        \\  .fill 510, 8, 0
    );
}

const SYS_WRITEC = 0x03;
const SYS_EXIT_EXTENDED = 0x20;
const ADP_Stopped_ApplicationExit = 0x20026;

inline fn semihost(op: usize, arg: usize) void {
    asm volatile ("hlt #0xf000"
        :
        : [op] "{x0}" (op),
          [arg] "{x1}" (arg),
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
