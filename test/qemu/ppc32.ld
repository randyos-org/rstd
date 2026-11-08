/* mac99 loads each segment at its physical address's low 28 bits + 0x1000000 and starts there, not at the entry: linked at 0x1000000,
 * physical 0, `_start` first. */
ENTRY(_start)
PHDRS
{
    image PT_LOAD AT(0) FLAGS(7);
}
SECTIONS
{
    . = 0x1000000;
    .text : { *(.text.boot) *(.text .text.*) } :image
    .rodata : { *(.rodata .rodata.*) } :image
    .data : { *(.data .data.*) *(.sdata .sdata.*) *(.got .got2) } :image
    . = ALIGN(8);
    .bss : {
        __bss_start = .;
        *(.bss .bss.*) *(.sbss .sbss.*) *(COMMON)
        . = ALIGN(8);
        __bss_end = .;
    } :image
    . = ALIGN(16);
    __stack_bottom = .;
    . += 0x100000;
    __stack_top = .;
}
