//! SPARC64 boot for `sun4u` (OpenBIOS); console and exit through the Open Firmware client from `%o4`, the only one the firmware leaves.
//! Exit is `reset-all`, which `-no-reboot` turns into qemu exiting (the summary carries the result).

const ofw = @import("ofw.zig");

// non-zero so it lives in .data, out of the bss `_start` clears after storing it
export var client_interface: usize = 1;

comptime {
    asm (
        \\.section .text.boot, "ax"
        \\.globl _start
        \\_start:
        \\  setx client_interface, %g2, %g1
        \\  stx %o4, [%g1]
        \\  setx __stack_top, %g2, %g1
        \\  sub %g1, 2047 + 512, %sp
        \\  rdpr %pstate, %g1
        \\  or %g1, 0x10, %g1
        \\  wrpr %g1, %pstate
        \\  wr %g0, 4, %fprs
        \\  setx __bss_start, %g2, %o0
        \\  setx __bss_end, %g2, %o1
        \\1:
        \\  cmp %o0, %o1
        \\  bgeu,pn %xcc, 2f
        \\  nop
        \\  stx %g0, [%o0]
        \\  ba 1b
        \\  add %o0, 8, %o0
        \\2:
        \\  call qemutest_main
        \\  nop
        \\3:
        \\  ba 3b
        \\  nop
    );
}

const of = ofw.Client(u64, u64, &client_interface);
const call = of.call;
const serviceName = of.serviceName;

/// The tests' heap, claimed from the firmware: mapped wherever it chooses, so
/// none of the 4 MB it maps of the image goes to it. Empty if refused.
pub fn heapRegion(comptime len: usize, comptime alignment: usize) []align(alignment) u8 {
    var claim = [_]u64{ serviceName("claim"), 3, 1, 0, len, alignment, 0 };
    call(&claim);
    if (claim[6] == 0 or claim[6] == ~@as(u64, 0)) return &.{};
    const p: [*]align(alignment) u8 = @ptrFromInt(claim[6]);
    return p[0..len];
}

var stdout_handle: u32 = 0;

fn openStdout() void {
    var find = [_]u64{ serviceName("finddevice"), 1, 1, serviceName("/chosen"), 0 };
    call(&find);
    var handle: [4]u8 align(4) = undefined;
    var prop = [_]u64{ serviceName("getprop"), 4, 1, find[4], serviceName("stdout"), @intFromPtr(&handle), 4, 0 };
    call(&prop);
    stdout_handle = (@as(u32, handle[0]) << 24) | (@as(u32, handle[1]) << 16) | (@as(u32, handle[2]) << 8) | handle[3];
}

pub fn writeByte(c: u8) void {
    if (stdout_handle == 0) openStdout();
    of.write(stdout_handle, c);
}

pub fn exit(code: usize) noreturn {
    _ = code;
    of.resetAll();
    while (true) {}
}
