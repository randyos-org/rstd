//! x86_64 boot for QEMU's `pc` machine (`-kernel`). QEMU takes an ELF64 kernel
//! only through the PVH entry note, which starts it in 32-bit protected mode
//! with paging off; `_start32` builds an identity map of the first gigabyte,
//! enters long mode and turns on the SSE, AVX and AVX-512 state the CPU has
//! (XCR0 from CPUID leaf 0xd). Console is COM1; exit powers the machine off
//! through the ACPI PM control register, so the run step judges the result
//! from the printed summary.

comptime {
    asm (
        \\.section .note.Xen, "a", @note
        \\.balign 4
        \\.long 4
        \\.long 4
        \\.long 18
        \\.asciz "Xen"
        \\.long _start32
        \\.section .rodata
        \\.balign 8
        \\gdt64:
        \\.quad 0
        \\.quad 0x00af9a000000ffff
        \\.quad 0x00cf92000000ffff
        \\gdt64_ptr:
        \\.word 23
        \\.long gdt64
        \\.section .text.boot32, "ax"
        \\.code32
        \\.globl _start32
        \\.globl _start
        \\_start32:
        \\_start:
        \\  cli
        \\  movl $__stack_top, %esp
        \\  movl $__pdpt, %eax
        \\  orl $3, %eax
        \\  movl %eax, __pml4
        \\  movl $__pd, %eax
        \\  orl $3, %eax
        \\  movl %eax, __pdpt
        \\  xorl %ecx, %ecx
        \\1:
        \\  movl %ecx, %eax
        \\  shll $21, %eax
        \\  orl $0x83, %eax
        \\  movl %eax, __pd(,%ecx,8)
        \\  incl %ecx
        \\  cmpl $512, %ecx
        \\  jne 1b
        \\  movl %cr4, %eax
        \\  orl $0x620, %eax
        \\  movl %eax, %cr4
        \\  movl $__pml4, %eax
        \\  movl %eax, %cr3
        \\  movl $0xc0000080, %ecx
        \\  rdmsr
        \\  orl $0x100, %eax
        \\  wrmsr
        \\  movl %cr0, %eax
        \\  andl $0xfffffffb, %eax
        \\  orl $0x80000003, %eax
        \\  movl %eax, %cr0
        \\  lgdt gdt64_ptr
        \\  ljmp $0x08, $long_mode
        \\.code64
        \\long_mode:
        \\  movw $0x10, %ax
        \\  movw %ax, %ds
        \\  movw %ax, %es
        \\  movw %ax, %ss
        \\  xorl %eax, %eax
        \\  movw %ax, %fs
        \\  movw %ax, %gs
        \\  movq $__stack_top, %rsp
        \\  call qemutest_cpu_init
        \\  jmp 8f
        \\// what every core sets up (`smp.zig`); keeps %rdi
        \\.globl qemutest_cpu_init
        \\qemutest_cpu_init:
        \\  movl $1, %eax
        \\  cpuid
        \\  btl $26, %ecx
        \\  jnc 2f
        \\  movq %cr4, %rax
        \\  orl $0x40000, %eax
        \\  movq %rax, %cr4
        \\  movl $0xd, %eax
        \\  xorl %ecx, %ecx
        \\  cpuid
        \\  andl $0xe7, %eax
        \\  xorl %edx, %edx
        \\  xorl %ecx, %ecx
        \\  xsetbv
        \\2:
        \\  fninit
        \\  ret
        \\8:
        \\  movq $__bss_start, %rdi
        \\  movq $__bss_end, %rcx
        \\  subq %rdi, %rcx
        \\  shrq $3, %rcx
        \\  xorl %eax, %eax
        \\  rep stosq
        \\  call qemutest_main
        \\3:
        \\  hlt
        \\  jmp 3b
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
