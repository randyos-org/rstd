//! The aarch64 rigs' core bring-up for `smp.zig`: PSCI `CPU_ON`, the event stream, and
//! `wfe`/`sev` waits.

const std = @import("std");
const sysinfo = @import("builtin");
const smp = @import("../smp.zig");

const Atomic = std.atomic.Value;
const stack_size = smp.stack_size;

comptime {
    if (!sysinfo.single_threaded and sysinfo.target.cpu.arch == .aarch64) asm (std.fmt.comptimePrint(
            \\.text
            \\.globl qemutest_ap_entry
            \\qemutest_ap_entry:
            \\  mov x19, x0
            \\  bl qemutest_cpu_init
            \\  adrp x1, qemutest_ap_stacks
            \\  add x1, x1, :lo12:qemutest_ap_stacks
            \\  mov x2, #{[stack]d}
            \\  madd x1, x19, x2, x1
            \\  mov sp, x1
            \\  mov x0, x19
            \\  bl qemutest_ap_main
            \\1:
            \\  wfe
            \\  b 1b
        , .{ .stack = stack_size }));
}

extern fn qemutest_ap_entry() callconv(.naked) noreturn;

/// `virt`'s cores past the first, by PSCI `CPU_ON` (over `hvc`, as
/// `aarch64.zig`'s exit) until one is missing: affinity 0 is the index
/// up to 8 cores, passed on as the context.
pub fn startCores(max: usize) usize {
    var n: usize = 0;
    while (n < max) : (n += 1) {
        const i = n + 1;
        if (psci(0xc400_0003, i, @intFromPtr(&qemutest_ap_entry), i) != 0) break;
    }
    return n;
}

fn psci(f: u32, a: usize, b: usize, c: usize) isize {
    var x1: usize = undefined;
    var x2: usize = undefined;
    var x3: usize = undefined;
    return asm volatile ("hvc #0"
        : [r] "={x0}" (-> isize),
          [o1] "={x1}" (x1),
          [o2] "={x2}" (x2),
          [o3] "={x3}" (x3),
        : [f] "{x0}" (@as(usize, f)),
          [a] "{x1}" (a),
          [b] "{x2}" (b),
          [c] "{x3}" (c),
        : .{ .x4 = true, .x5 = true, .x6 = true, .x7 = true, .x8 = true, .x9 = true, .x10 = true, .x11 = true, .x12 = true, .x13 = true, .x14 = true, .x15 = true, .x16 = true, .x17 = true, .memory = true });
}

pub fn setThreadPointer(tp: usize) void {
    asm volatile ("msr tpidr_el0, %[tp]"
        :
        : [tp] "r" (tp),
        : .{ .memory = true });
}

/// The virtual counter's event stream (`CNTKCTL_EL1.EVNTEN`, bit 12's
/// rising edges, about 130 µs at QEMU's 62.5 MHz) ends every `wfe` in
/// time, as Linux's arch timer does: QEMU 11's TCG can lose a `sev`
/// to a core going to sleep (issue #467).
pub fn initCore() void {
    const old = asm volatile ("mrs %[v], cntkctl_el1"
        : [v] "=r" (-> u64),
    );
    // EVNTEN (bit 2) on, EVNTDIR (3) rising, EVNTI (7:4) 12
    const v = (old & ~@as(u64, 0xfc)) | 0x4 | (12 << 4);
    asm volatile ("msr cntkctl_el1, %[v]; isb"
        :
        : [v] "r" (v),
        : .{ .memory = true });
}

pub fn arm(_: *Atomic(u32)) void {}

/// `wfe`: a `signal` since `seen` was read left the event set, so it
/// returns at once.
pub fn sleep(doorbell: *Atomic(u32), seen: u32) void {
    if (doorbell.load(.acquire) == seen) asm volatile ("wfe" ::: .{ .memory = true });
}

pub fn signal() void {
    asm volatile ("dsb ish; sev" ::: .{ .memory = true });
}
