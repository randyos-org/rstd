/*
 * QEMU's mac99 loader puts each segment at its physical address's low 28
 * bits plus 0x1000000, and OpenBIOS starts the image there, not at the entry.
 * So the image is linked at 0x1000000 with a physical address of 0, and
 * `_start` comes first.
 */
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
