//! 32-bit x86 boot for QEMU's `pc` machine (`-kernel`, Multiboot, 32-bit
//! protected mode with paging off). Console is COM1; exit powers the machine
//! off through the ACPI PM control register, so the run step judges the
//! result from the printed summary rather than an exit code.

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

inline fn outb(port: u16, value: u8) void {
    asm volatile ("outb %[value], %[port]"
        :
        : [value] "{al}" (value),
          [port] "{dx}" (port),
    );
}

inline fn outw(port: u16, value: u16) void {
    asm volatile ("outw %[value], %[port]"
        :
        : [value] "{ax}" (value),
          [port] "{dx}" (port),
    );
}

inline fn inb(port: u16) u8 {
    return asm volatile ("inb %[port], %[value]"
        : [value] "={al}" (-> u8),
        : [port] "{dx}" (port),
    );
}

const com1: u16 = 0x3f8;

pub fn writeByte(c: u8) void {
    // the transmit holding register is empty
    while (inb(com1 + 5) & 0x20 == 0) {}
    outb(com1, c);
}

pub fn exit(code: usize) noreturn {
    _ = code;
    // S5, at the PIIX4 PM control register and at the newer machines' one
    outw(0xb004, 0x2000);
    outw(0x604, 0x2000);
    while (true) asm volatile ("hlt");
}
