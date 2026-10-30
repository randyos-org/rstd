//! AArch64 boot for QEMU's `virt` machine (`-kernel`, EL1) and ARM
//! semihosting for console and exit (`-semihosting`). Where `hlt` does not
//! semihost but traps (HVF), the console is the PL011 and exit PSCI's.
//!
//! With the MMU off every data access is Device memory, which faults on the
//! unaligned accesses compiled code makes, so `_start` builds a one-block
//! identity map first: 0-1 GiB as device, 1-2 GiB (RAM) as normal cached.
//!
//! SVE and SME, where the CPU has them, get vectors of up to 512 bits
//! (`ZCR_EL1`, `SMCR_EL1`; reset leaves 128), and SME2 gets ZT0 (`EZT0`).

comptime {
    asm (
        \\.section .text.boot
        \\.globl _start
        \\_start:
        \\  bl qemutest_cpu_init
        \\  b 8f
        \\// what every core sets up (`smp.zig`); uses x0 to x2 only
        \\.globl qemutest_cpu_init
        \\qemutest_cpu_init:
        \\  adr x0, qemutest_vectors
        \\  msr vbar_el1, x0
        \\  mrs x0, cpacr_el1
        \\  orr x0, x0, #(3 << 20)
        \\  orr x0, x0, #(3 << 16)
        \\  orr x0, x0, #(3 << 24)
        \\  msr cpacr_el1, x0
        \\  isb
        \\  mrs x1, id_aa64pfr0_el1
        \\  ubfx x1, x1, #32, #4
        \\  cbz x1, 4f
        \\  mov x1, #3
        \\  msr s3_0_c1_c2_0, x1
        \\4:
        \\  mrs x1, id_aa64pfr1_el1
        \\  ubfx x1, x1, #24, #4
        \\  cbz x1, 5f
        \\  mov x2, #3
        \\  cmp x1, #2
        \\  b.lo 6f
        \\  orr x2, x2, #(1 << 30)
        \\6:
        \\  msr s3_0_c1_c2_6, x2
        \\5:
        \\  isb
        \\  ldr x0, =0x00ff
        \\  msr mair_el1, x0
        \\  ldr x0, =0x200803519
        \\  msr tcr_el1, x0
        \\  adrp x0, qemutest_l1
        \\  msr ttbr0_el1, x0
        \\  isb
        \\  mrs x0, sctlr_el1
        \\  orr x0, x0, #1
        \\  orr x0, x0, #(1 << 2)
        \\  orr x0, x0, #(1 << 12)
        \\  bic x0, x0, #(1 << 1)
        \\  msr sctlr_el1, x0
        \\  isb
        \\  ret
        \\8:
        \\  ldr x0, =__stack_top
        \\  mov sp, x0
        \\  ldr x0, =__bss_start
        \\  ldr x1, =__bss_end
        \\1:
        \\  cmp x0, x1
        \\  b.hs 2f
        \\  str xzr, [x0], #8
        \\  b 1b
        \\2:
        \\  bl qemutest_main
        \\3:
        \\  wfe
        \\  b 3b
        \\// every exception: a trapped `hlt #0xf000` is noted and skipped; the
        \\// rest stop here
        \\.balign 2048
        \\qemutest_vectors:
        \\.rept 16
        \\  b qemutest_sync
        \\  .balign 128
        \\.endr
        \\qemutest_sync:
        \\  stp x2, x3, [sp, #-16]!
        \\  mrs x2, esr_el1
        \\  lsr x2, x2, #26
        \\  cbnz x2, 7f
        \\  mrs x2, elr_el1
        \\  ldr w2, [x2]
        \\  movz w3, #0xd45e, lsl #16
        \\  cmp w2, w3
        \\  b.ne 7f
        \\  adrp x2, qemutest_no_semihost
        \\  mov w3, #1
        \\  strb w3, [x2, :lo12:qemutest_no_semihost]
        \\  mrs x2, elr_el1
        \\  add x2, x2, #4
        \\  msr elr_el1, x2
        \\  ldp x2, x3, [sp], #16
        \\  eret
        \\7:
        \\  wfe
        \\  b 7b
        \\.section .data.pgtbl, "aw"
        \\.balign 4096
        \\qemutest_l1:
        \\  .quad 0x00000405
        \\  .quad 0x40000701
        \\  .fill 510, 8, 0
    );
}

const SYS_WRITEC = 0x03;
const SYS_EXIT_EXTENDED = 0x20;
const ADP_Stopped_ApplicationExit = 0x20026;

inline fn semihost(op: usize, arg: usize) void {
    asm volatile ("hlt #0xf000"
        :
        : [op] "{x0}" (op),
          [arg] "{x1}" (arg),
        : .{ .memory = true });
}

/// Set by the vector when `hlt` trapped rather than semihosted.
export var qemutest_no_semihost: u8 = 0;

fn semihosted() bool {
    return @as(*volatile u8, &qemutest_no_semihost).* == 0;
}

/// `virt`'s PL011.
const uart: [*]volatile u32 = @ptrFromInt(0x0900_0000);
var uart_on = false;

pub fn writeByte(c: u8) void {
    var ch = c;
    if (semihosted()) {
        semihost(SYS_WRITEC, @intFromPtr(&ch));
        if (semihosted()) return;
    }
    if (!uart_on) {
        // UARTEN, TXE, RXE
        uart[0x30 / 4] = 0x301;
        uart_on = true;
    }
    // FR.TXFF
    while (uart[0x18 / 4] & 0x20 != 0) {}
    uart[0] = c;
}

/// Without semihosting the status is lost; the printed summary carries it.
pub fn exit(code: usize) noreturn {
    var block = [2]usize{ ADP_Stopped_ApplicationExit, code };
    semihost(SYS_EXIT_EXTENDED, @intFromPtr(&block));
    // PSCI SYSTEM_OFF
    asm volatile ("hvc #0"
        :
        : [f] "{x0}" (0x84000008),
        : .{ .memory = true });
    while (true) {}
}
