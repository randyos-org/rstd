//! LoongArch64 boot for QEMU's `virt` machine (`-kernel`, direct address
//! mode, no MMU needed). Console is the 16550 UART; exit powers the machine
//! off through the ACPI GED sleep register, so the run step judges the
//! result from the printed summary rather than an exit code.

comptime {
    asm (
        \\.section .text.boot
        \\.globl _start
        \\_start:
        \\  li.w $t0, 7
        \\  csrwr $t0, 0x2
        \\  la.local $sp, __stack_top
        \\  la.local $t0, __bss_start
        \\  la.local $t1, __bss_end
        \\1:
        \\  bgeu $t0, $t1, 2f
        \\  st.d $zero, $t0, 0
        \\  addi.d $t0, $t0, 8
        \\  b 1b
        \\2:
        \\  bl qemutest_main
        \\3:
        \\  idle 0
        \\  b 3b
    );
}

const uart: [*]volatile u8 = @ptrFromInt(0x1fe001e0);
const ged_sleep_ctl: *volatile u8 = @ptrFromInt(0x100e001c);

pub fn writeByte(c: u8) void {
    while (uart[5] & 0x20 == 0) {}
    uart[0] = c;
}

pub fn exit(code: usize) noreturn {
    _ = code;
    ged_sleep_ctl.* = 0x34; // S5
    while (true) {}
}
