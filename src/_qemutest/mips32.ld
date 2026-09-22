ENTRY(_start)
SECTIONS
{
    . = 0x80100000;
    .text : { *(.text.boot) *(.text .text.*) }
    .rodata : { *(.rodata .rodata.*) }
    .data : { *(.data .data.*) *(.sdata .sdata.*) }
    . = ALIGN(4);
    .bss : {
        __bss_start = .;
        *(.bss .bss.*) *(.sbss .sbss.*) *(COMMON)
        . = ALIGN(4);
        __bss_end = .;
    }
    . = ALIGN(16);
    . += 0x40000;
    __stack_top = .;
}
