ENTRY(_start32)
SECTIONS
{
    . = 0x100000;
    .text : { *(.text.boot32) *(.text.boot) *(.text .text.*) }
    .note.Xen : { KEEP(*(.note.Xen)) }
    .rodata : { *(.rodata .rodata.*) }
    .data : { *(.data .data.*) }
    /* thread-locals: each core's copy is made at start (`smp.zig`) */
    . = ALIGN(64);
    .tdata : { *(.tdata .tdata.*) }
    .tbss : { *(.tbss .tbss.*) }
    __tls_start = ADDR(.tdata);
    __tls_image = SIZEOF(.tdata);
    __tls_size = ADDR(.tbss) + SIZEOF(.tbss) - ADDR(.tdata);
    __tls_align = MAX(ALIGNOF(.tdata), ALIGNOF(.tbss));
    . = ALIGN(4096);
    .bss : {
        __bss_start = .;
        *(.bss .bss.*) *(COMMON)
        . = ALIGN(4096);
        __bss_end = .;
    }
    . = ALIGN(4096);
    __pml4 = .;
    . += 4096;
    __pdpt = .;
    . += 4096;
    __pd = .;
    . += 4096;
    . = ALIGN(16);
    __stack_bottom = .;
    . += 0x100000;
    __stack_top = .;
}
