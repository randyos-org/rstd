ENTRY(_start32)
SECTIONS
{
    . = 0x100000;
    .text : { *(.text.boot32) *(.text.boot) *(.text .text.*) }
    .note.Xen : { KEEP(*(.note.Xen)) }
    .rodata : { *(.rodata .rodata.*) }
    .data : { *(.data .data.*) }
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
    . += 0x100000;
    __stack_top = .;
}
