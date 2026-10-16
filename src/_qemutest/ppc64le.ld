/*
 * The firmware loader shifts every load address (and the entry point) up by
 * 0x400000, so the image is linked at 0xc00000 but given a physical address
 * of 0x800000, and the entry is expressed 0x400000 below `_start`. Net effect:
 * it runs where it was linked, and absolute addresses in code and data hold.
 */
ENTRY(qemutest_entry)
PHDRS
{
    image PT_LOAD AT(0x800000) FLAGS(7);
}
SECTIONS
{
    . = 0xc00000;
    .text : { *(.text.boot) *(.text .text.*) } :image
    .rodata : { *(.rodata .rodata.*) } :image
    .data : { *(.data .data.*) *(.got .toc .tocbss) } :image
    . = ALIGN(8);
    .bss : {
        __bss_start = .;
        *(.bss .bss.*) *(COMMON)
        . = ALIGN(8);
        __bss_end = .;
    } :image
    . = ALIGN(16);
    __stack_bottom = .;
    . += 0x100000;
    __stack_top = .;
    qemutest_entry = _start - 0x400000;
}
