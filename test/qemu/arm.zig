//! 32-bit ARM boot for `virt` (cortex-a15, SVC mode) with semihosting console and exit. `_start` turns on the FPU (`vmsr` by encoding: a
//! build without VFP rejects the mnemonic, exiting 0) and the MMU with a 1 MiB section identity map, since with it off every access
//! is strongly ordered and unaligned ones fault.

comptime {
    asm (
        \\.section .text.boot
        \\.globl _start
        \\_start:
        \\  mrc p15, 0, r0, c1, c0, 2
        \\  orr r0, r0, #(0xf << 20)
        \\  mcr p15, 0, r0, c1, c0, 2
        \\  isb
        \\  mov r0, #0x40000000
        \\  .inst 0xeee80a10 @ vmsr fpexc, r0
        \\  ldr r0, =qemutest_ttb
        \\  mcr p15, 0, r0, c2, c0, 0
        \\  mov r0, #0
        \\  mcr p15, 0, r0, c2, c0, 2
        \\  mov r0, #3
        \\  mcr p15, 0, r0, c3, c0, 0
        \\  mrc p15, 0, r0, c1, c0, 0
        \\  orr r0, r0, #1
        \\  orr r0, r0, #4
        \\  orr r0, r0, #0x1000
        \\  bic r0, r0, #2
        \\  mcr p15, 0, r0, c1, c0, 0
        \\  isb
        \\  ldr sp, =__stack_top
        \\  ldr r0, =__bss_start
        \\  ldr r1, =__bss_end
        \\  mov r2, #0
        \\1:
        \\  cmp r0, r1
        \\  bhs 2f
        \\  str r2, [r0], #4
        \\  b 1b
        \\2:
        \\  bl qemutest_main
        \\3:
        \\  wfe
        \\  b 3b
        \\.ltorg
        \\.section .data.ttb, "aw"
        \\.balign 16384
        \\qemutest_ttb:
        \\  .set i, 0
        \\  .rept 1024
        \\    .word (i << 20) | 0xC06
        \\    .set i, i + 1
        \\  .endr
        \\  .rept 128
        \\    .word (i << 20) | 0x1C0E
        \\    .set i, i + 1
        \\  .endr
        \\  .fill 4096 - 1152, 4, 0
    );
}

const SYS_WRITEC = 0x03;
const SYS_EXIT_EXTENDED = 0x20;
const ADP_Stopped_ApplicationExit = 0x20026;

inline fn semihost(op: usize, arg: usize) void {
    asm volatile ("svc #0x123456"
        :
        : [op] "{r0}" (op),
          [arg] "{r1}" (arg),
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
