ENTRY(_start)
SECTIONS
{
    . = 0xffffffff80100000;
    .text : { *(.text.boot) *(.text .text.*) }
    .rodata : { *(.rodata .rodata.*) }
    .data : { *(.data .data.*) *(.sdata .sdata.*) }
    . = ALIGN(8);
    .bss : {
        __bss_start = .;
        *(.bss .bss.*) *(.sbss .sbss.*) *(COMMON)
        . = ALIGN(8);
        __bss_end = .;
    }
    . = ALIGN(16);
    . += 0x40000;
    __stack_top = .;
}
