//! 32-bit x86 boot for `pc` (Multiboot, protected mode, paging off). Console COM1; exit powers off through the ACPI PM control register
//! (the summary carries the result).

const pc = @import("pc.zig");

comptime {
    asm (
        \\.section .multiboot, "a"
        \\.balign 4
        \\.long 0x1BADB002
        \\.long 0
        \\.long -0x1BADB002
        \\.section .text.boot, "ax"
        \\.globl _start
        \\_start:
        \\  cli
        \\  movl $__stack_top, %esp
        \\  movl %cr0, %eax
        \\  andl $0xfffffffb, %eax
        \\  orl $2, %eax
        \\  movl %eax, %cr0
        \\  movl %cr4, %eax
        \\  orl $0x600, %eax
        \\  movl %eax, %cr4
        \\  fninit
        \\  movl $__bss_start, %edi
        \\  movl $__bss_end, %ecx
        \\  subl %edi, %ecx
        \\  shrl $2, %ecx
        \\  xorl %eax, %eax
        \\  rep stosl
        \\  call qemutest_main
        \\1:
        \\  hlt
        \\  jmp 1b
    );
}

pub const writeByte = pc.writeByte;

pub fn exit(code: usize) noreturn {
    _ = code;
    pc.powerOff();
}
