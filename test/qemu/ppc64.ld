/* The loader shifts load addresses and entry up 0x400000: linked at 0xc00000, physical 0x800000, entry 0x400000 below `_start`. */
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
