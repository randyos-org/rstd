//! 32-bit PowerPC boot for `mac99` (OpenBIOS) and `pseries` under VOF (64-bit CPUs in 32-bit mode); both pass the Open Firmware client in r5.
//! Under OpenBIOS console and exit are its `write` and `reset-all` (`-no-reboot` exits qemu); on pseries, ppc64.zig's hypercalls.
//! Load addresses: `ppc32.ld` (mac99), `ppc64.ld` (pseries).

const sysinfo = @import("builtin");
const ofw = @import("ofw.zig");

// non-zero so it lives in .data, out of the bss `_start` clears after storing it
export var client_interface: usize = 1;

/// MSR.VEC, where the CPU has AltiVec, and MSR.VSX, where it has VSX.
const vec = if (sysinfo.target.cpu.has(.powerpc, .vsx))
    "  oris 3, 3, 0x0280\n"
else if (sysinfo.target.cpu.has(.powerpc, .altivec))
    "  oris 3, 3, 0x0200\n"
else
    "";

comptime {
    asm (
        \\.section .text.boot, "ax"
        \\.globl _start
        \\_start:
        \\  lis 6, client_interface@ha
        \\  stw 5, client_interface@l(6)
        \\  lis 1, __stack_top@ha
        \\  addi 1, 1, __stack_top@l
        \\  li 0, 0
        \\  stwu 0, -16(1)
        \\  mfmsr 3
        \\  ori 3, 3, 0x2000
        \\
    ++ vec ++
        \\  mtmsr 3
        \\  isync
        \\  lis 6, __bss_start@ha
        \\  addi 6, 6, __bss_start@l
        \\  lis 7, __bss_end@ha
        \\  addi 7, 7, __bss_end@l
        \\1:
        \\  cmplw 6, 7
        \\  bge 2f
        \\  stw 0, 0(6)
        \\  addi 6, 6, 4
        \\  b 1b
        \\2:
        \\  bl qemutest_main
        \\3:
        \\  b 3b
    );
}

const of = ofw.Client(u32, i32, &client_interface);
const call = of.call;
const serviceName = of.serviceName;

/// Whether this is pseries, where the console and exit are PAPR hypercalls,
/// as in ppc64.zig: its VOF has no `interpret`, and its `write` shows nothing.
var papr: ?bool = null;

fn isPapr() bool {
    if (papr == null) {
        var find = [_]u32{ serviceName("finddevice"), 1, 1, serviceName("/vdevice"), 0 };
        call(&find);
        papr = find[4] != 0xffffffff;
    }
    return papr.?;
}

var stdout_handle: u32 = 0;

fn openStdout() void {
    var find = [_]u32{ serviceName("finddevice"), 1, 1, serviceName("/chosen"), 0 };
    call(&find);
    var handle: u32 = 0;
    var prop = [_]u32{ serviceName("getprop"), 4, 1, find[4], serviceName("stdout"), @intFromPtr(&handle), 4, 0 };
    call(&prop);
    stdout_handle = handle;
}

/// A PAPR hypercall, as ppc64.zig's; the 64-bit CPU takes 64-bit registers.
inline fn hcall(token: usize, a1: usize, a2: usize, a3: usize) void {
    _ = asm volatile ("sc 1"
        : [ret] "={r3}" (-> usize),
        : [tok] "{r3}" (token),
          [a1] "{r4}" (a1),
          [a2] "{r5}" (a2),
          [a3] "{r6}" (a3),
        : .{ .memory = true, .r0 = true, .r7 = true, .r8 = true, .r9 = true, .r10 = true, .r11 = true, .r12 = true, .cr0 = true });
}

// PAPR's numbers, as ppc64.zig names them
const H_PUT_TERM_CHAR = 0x58;
const vty_reg = 0x71000000;
const H_RTAS = 0xf000;
const RTAS_POWER_OFF = 0x2003;

pub fn writeByte(c: u8) void {
    if (isPapr()) {
        // `H_PUT_TERM_CHAR` to the default vty, the character in the top
        // byte of r6's 64 bits, which only a 64-bit shift reaches
        _ = asm volatile (
            \\sldi 6, %[c], 56
            \\sc 1
            : [ret] "={r3}" (-> usize),
            : [tok] "{r3}" (H_PUT_TERM_CHAR),
              [reg] "{r4}" (vty_reg),
              [len] "{r5}" (1),
              [c] "r" (c),
            : .{ .memory = true, .r0 = true, .r6 = true, .r7 = true, .r8 = true, .r9 = true, .r10 = true, .r11 = true, .r12 = true, .cr0 = true });
        return;
    }
    if (stdout_handle == 0) openStdout();
    of.write(stdout_handle, c);
}

pub fn exit(code: usize) noreturn {
    _ = code;
    if (isPapr()) {
        // RTAS `power-off` through `H_RTAS`
        var args = [6]u32{ RTAS_POWER_OFF, 2, 1, 0, 0, 0 };
        hcall(H_RTAS, @intFromPtr(&args), 0, 0);
    } else {
        of.resetAll();
    }
    while (true) {}
}
