ENTRY(_start)
SECTIONS
{
    . = 0x400000;
    .text : { *(.text.boot) *(.text .text.*) }
    .rodata : { *(.rodata .rodata.*) }
    /* The firmware maps the image as loaded, so bss and the stack are in the
       file too, as part of .data. */
    .data : {
        *(.data .data.*)
        . = ALIGN(8);
        __bss_start = .;
        *(.bss .bss.*) *(COMMON)
        . = ALIGN(8);
        __bss_end = .;
        . = ALIGN(16);
        __stack_bottom = .;
        /* 512 KiB: the heap comes from the firmware, so the image has room */
        . += 0x80000;
        __stack_top = .;
    }
}
