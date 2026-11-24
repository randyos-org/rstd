//! The x86_64 rigs' core bring-up for `smp.zig`: the real-mode trampoline, INIT-SIPI-SIPI
//! through the local APIC, and `monitor`/`mwait` waits.

const std = @import("std");
const sysinfo = @import("builtin");
const smp = @import("../smp.zig");
const pc = @import("../pc.zig");

const Atomic = std.atomic.Value;
const max_cores = smp.max_cores;
const stack_size = smp.stack_size;
const fail = smp.fail;

/// Where the real-mode trampoline is copied; the SIPI vector is its page.
const trampoline = 0x8000;

comptime {
    if (!sysinfo.single_threaded and sysinfo.target.cpu.arch == .x86_64) asm (std.fmt.comptimePrint(
            \\.section .text.qemutest_ap, "ax"
            \\.code16
            \\.globl qemutest_ap_tramp
            \\qemutest_ap_tramp:
            \\  cli
            \\  movw %cs, %ax
            \\  movw %ax, %ds
            \\  lgdtl qemutest_ap_gdt_ptr - qemutest_ap_tramp
            \\  movl %cr0, %eax
            \\  orl $1, %eax
            \\  movl %eax, %cr0
            \\  ljmpl $0x08, ${[t]d} + (qemutest_ap_32 - qemutest_ap_tramp)
            \\.code32
            \\qemutest_ap_32:
            \\  movw $0x10, %ax
            \\  movw %ax, %ds
            \\  movw %ax, %es
            \\  movw %ax, %ss
            \\  movl %cr4, %eax
            \\  orl $0x620, %eax
            \\  movl %eax, %cr4
            \\  movl $__pml4, %eax
            \\  movl %eax, %cr3
            \\  movl $0xc0000080, %ecx
            \\  rdmsr
            \\  orl $0x100, %eax
            \\  wrmsr
            \\  movl %cr0, %eax
            \\  andl $0xfffffffb, %eax
            \\  orl $0x80000003, %eax
            \\  movl %eax, %cr0
            \\  ljmpl $0x18, $qemutest_ap_64
            \\.balign 8
            \\qemutest_ap_gdt:
            \\  .quad 0
            \\  .quad 0x00cf9a000000ffff
            \\  .quad 0x00cf92000000ffff
            \\  .quad 0x00af9a000000ffff
            \\qemutest_ap_gdt_ptr:
            \\  .word 31
            \\  .long {[t]d} + (qemutest_ap_gdt - qemutest_ap_tramp)
            \\.globl qemutest_ap_tramp_end
            \\qemutest_ap_tramp_end:
            \\.code64
            \\qemutest_ap_64:
            \\  movw $0x10, %ax
            \\  movw %ax, %ds
            \\  movw %ax, %es
            \\  movw %ax, %ss
            \\  xorl %eax, %eax
            \\  movw %ax, %fs
            \\  movw %ax, %gs
            \\  movl $1, %edi
            \\  lock xaddl %edi, qemutest_ap_ticket(%rip)
            \\  cmpl ${[tickets]d}, %edi
            \\  jae 1f
            \\  incl %edi
            \\  movl %edi, %eax
            \\  imulq ${[stack]d}, %rax, %rax
            \\  leaq qemutest_ap_stacks(%rip), %rsp
            \\  addq %rax, %rsp
            \\  call qemutest_cpu_init
            \\  call qemutest_ap_main
            \\1:
            \\  hlt
            \\  jmp 1b
            \\.text
        , .{ .t = trampoline, .tickets = max_cores - 1, .stack = stack_size }));
}

/// Tickets the started cores take, in the order they come; index is ticket + 1.
export var qemutest_ap_ticket: u32 = 0;

extern const qemutest_ap_tramp: u8;
extern const qemutest_ap_tramp_end: u8;

/// Maps 3-4 GiB uncached (the local APIC's page) through `pd_hi`.
var pd_hi: [512]u64 align(4096) = @splat(0);

/// `monitor`/`mwait` (CPUID.1:ECX.MONITOR).
var has_monitor = false;

/// The rest of QEMU's CPUs (fw_cfg's count), all at once by broadcast.
pub fn startCores(max: usize) usize {
    has_monitor = cpuid(1)[2] & (1 << 3) != 0;
    const others = @min(cpuCount() -| 1, max);
    if (others == 0) return 0;
    const apic_base = rdmsr(0x1b) & 0xffff_f000;
    if (apic_base >> 30 != 3) fail("local APIC outside 3-4 GiB");
    for (&pd_hi, 0..) |*e, i| e.* = ((3 << 30) + (i << 21)) | 0x9b;
    const pdpt: *[512]u64 = @ptrCast(@alignCast(@extern(*u8, .{ .name = "__pdpt" })));
    pdpt[3] = @intFromPtr(&pd_hi) | 3;
    asm volatile ("movq %%cr3, %%rax; movq %%rax, %%cr3" ::: .{ .rax = true, .memory = true });
    const len = @intFromPtr(&qemutest_ap_tramp_end) - @intFromPtr(&qemutest_ap_tramp);
    const dst: [*]u8 = @ptrFromInt(trampoline);
    @memcpy(dst[0..len], @as([*]const u8, @ptrCast(&qemutest_ap_tramp))[0..len]);
    const apic: [*]volatile u32 = @ptrFromInt(apic_base);
    // software-enabled, spurious vector 0xff
    apic[0xf0 / 4] |= 0x1ff;
    // INIT, then SIPI twice, to all but this core
    ipi(apic, 0xc4500);
    delay(1 << 20);
    for (0..2) |_| {
        ipi(apic, 0xc4600 | (trampoline >> 12));
        delay(1 << 16);
    }
    return others;
}

fn ipi(apic: [*]volatile u32, icr: u32) void {
    apic[0x310 / 4] = 0;
    apic[0x300 / 4] = icr;
    // delivery status
    while (apic[0x300 / 4] & (1 << 12) != 0) std.atomic.spinLoopHint();
}

fn delay(n: usize) void {
    for (0..n) |_| std.atomic.spinLoopHint();
}

/// QEMU's CPU count from fw_cfg (`FW_CFG_NB_CPUS`), 1 if it has none.
fn cpuCount() usize {
    var sig: [4]u8 = undefined;
    fwCfgRead(0, &sig);
    if (!std.mem.eql(u8, &sig, "QEMU")) return 1;
    var n: [2]u8 = undefined;
    fwCfgRead(5, &n);
    return std.mem.readInt(u16, &n, .little);
}

fn fwCfgRead(key: u16, out: []u8) void {
    pc.outw(0x510, key);
    for (out) |*b| b.* = pc.inb(0x511);
}

fn cpuid(leaf: u32) [4]u32 {
    var a: u32 = undefined;
    var b: u32 = undefined;
    var c: u32 = undefined;
    var d: u32 = undefined;
    asm volatile ("cpuid"
        : [a] "={eax}" (a),
          [b] "={ebx}" (b),
          [c] "={ecx}" (c),
          [d] "={edx}" (d),
        : [l] "{eax}" (leaf),
          [s] "{ecx}" (@as(u32, 0)),
    );
    return .{ a, b, c, d };
}

fn rdmsr(msr: u32) u64 {
    var lo: u32 = undefined;
    var hi: u32 = undefined;
    asm volatile ("rdmsr"
        : [lo] "={eax}" (lo),
          [hi] "={edx}" (hi),
        : [m] "{ecx}" (msr),
    );
    return (@as(u64, hi) << 32) | lo;
}

fn wrmsr(msr: u32, v: u64) void {
    asm volatile ("wrmsr"
        :
        : [m] "{ecx}" (msr),
          [lo] "{eax}" (@as(u32, @truncate(v))),
          [hi] "{edx}" (@as(u32, @truncate(v >> 32))),
        : .{ .memory = true });
}

pub fn initCore() void {}

/// FS base, with the self pointer local-exec code loads from `%fs:0`.
pub fn setThreadPointer(tp: usize) void {
    @as(*usize, @ptrFromInt(tp)).* = tp;
    wrmsr(0xc000_0100, tp);
}

pub fn arm(doorbell: *Atomic(u32)) void {
    if (has_monitor) asm volatile ("monitor"
        :
        : [a] "{rax}" (doorbell),
          [c] "{ecx}" (@as(u32, 0)),
          [d] "{edx}" (@as(u32, 0)),
        : .{ .memory = true });
}

/// `mwait` on the armed doorbell (a pause under TCG with more than one
/// CPU), else a pause.
pub fn sleep(doorbell: *Atomic(u32), seen: u32) void {
    if (has_monitor and doorbell.load(.acquire) == seen) {
        asm volatile ("mwait"
            :
            : [a] "{eax}" (@as(u32, 0)),
              [c] "{ecx}" (@as(u32, 0)),
            : .{ .memory = true });
    } else std.atomic.spinLoopHint();
}

/// A doorbell's write ends `mwait`.
pub fn signal() void {}
