//! UEFI boot for the OVMF (x86_64) and AAVMF (aarch64) rigs: `main` takes heap and stack from boot services, turns on vector state,
//! masks interrupts and never calls the firmware again (its handlers save only 128-bit SIMD, clobbering SVE and SME). So the console
//! is the UART (COM1, PL011) and exit `ResetSystem`. `std.Io` is rstd's UEFI one, its clock unstarted (it ticks on interrupts).

const std = @import("std");
const sysinfo = @import("builtin");
const uefi = std.os.uefi;
const heap = @import("heap.zig");
const pc = @import("pc.zig");

/// rstd's UEFI `std.Io`, for `std_os_options.Io` and the runner's console Io.
pub const Io = @import("qemutest_root").uefi_io;

/// The tests' stack, as the bare rigs' linker scripts give it.
pub var stack: []align(page) u8 = &.{};
const stack_len = 1 << 20;
const page = 4096;

var heap_region: []align(heap.max_align) u8 = &.{};

extern fn qemutest_main() callconv(.c) noreturn;

pub fn main() noreturn {
    const bs = uefi.system_table.boot_services.?;
    // the watchdog BDS arms for a boot option (five minutes); the tests outlast it
    bs.setWatchdogTimer(0, 0, null) catch {};
    // for `addr2line`: a panic's addresses less this, plus the PE image base
    if (bs.handleProtocol(uefi.protocol.LoadedImage, uefi.handle) catch null) |image| {
        say("\nUEFI image at 0x");
        for (0..@sizeOf(usize) * 2) |i| writeByte("0123456789abcdef"[(@intFromPtr(image.image_base) >> @intCast((@sizeOf(usize) * 2 - 1 - i) * 4)) & 0xf]);
        say("\n");
    }
    heap_region = pages(bs, heap.max_len, heap.max_align);
    stack = pages(bs, stack_len, page);
    if (stack.len == 0) {
        say("no stack from the firmware\n");
        exit(1);
    }
    enableVectors();
    const top = @intFromPtr(stack.ptr) + stack.len;
    // a zero frame pointer ends the panic handler's walk
    switch (sysinfo.target.cpu.arch) {
        // the Microsoft ABI's 32 bytes of shadow space above the return address
        .x86_64 => asm volatile (
            \\cli
            \\leaq -32(%%rax), %%rsp
            \\xorl %%ebp, %%ebp
            \\callq *%%rcx
            \\ud2
            :
            : [top] "{rax}" (top),
              [f] "{rcx}" (&qemutest_main),
            : .{ .memory = true }),
        .aarch64 => asm volatile (
            \\msr daifset, #3
            \\mov sp, x0
            \\mov x29, xzr
            \\mov x30, xzr
            \\blr x1
            \\brk #0
            :
            : [top] "{x0}" (top),
              [f] "{x1}" (&qemutest_main),
            : .{ .memory = true }),
        else => @compileError("qemutest: no UEFI boot for " ++ @tagName(sysinfo.target.cpu.arch)),
    }
    unreachable;
}

/// `len` bytes at `alignment` from boot services; empty if they have none.
fn pages(bs: *uefi.tables.BootServices, comptime len: usize, comptime alignment: usize) []align(alignment) u8 {
    const got = bs.allocatePages(.any, .loader_data, (len + alignment - page) / page) catch return &.{};
    const start = std.mem.alignForward(usize, @intFromPtr(got.ptr), alignment);
    return @as([*]align(alignment) u8, @ptrFromInt(start))[0..len];
}

/// The heap `main` took, which must be there before interrupts are masked.
pub fn heapRegion(comptime len: usize, comptime alignment: usize) []align(alignment) u8 {
    comptime std.debug.assert(len == heap.max_len and alignment == heap.max_align);
    return heap_region;
}

/// SSE's successors on x86_64 (XCR0 from CPUID leaf 0xd), SVE and SME on
/// aarch64 at up to 512 bits, as `x86_64.zig` and `aarch64.zig` do; the
/// firmware leaves them off. AAVMF runs at EL1 here.
fn enableVectors() void {
    switch (sysinfo.target.cpu.arch) {
        .x86_64 => {
            const ecx = asm volatile ("cpuid"
                : [_] "={ecx}" (-> u32),
                : [_] "{eax}" (1),
                  [_] "{ecx}" (0),
                : .{ .ebx = true, .edx = true });
            // XSAVE
            if (ecx & (1 << 26) == 0) return;
            asm volatile (
                \\movq %%cr4, %%rax
                \\orl $0x40000, %%eax
                \\movq %%rax, %%cr4
                ::: .{ .rax = true, .cc = true });
            const xcr0 = asm volatile ("cpuid"
                : [_] "={eax}" (-> u32),
                : [_] "{eax}" (0xd),
                  [_] "{ecx}" (0),
                : .{ .ebx = true, .ecx = true, .edx = true });
            // x87, SSE, AVX and the three AVX-512 states
            asm volatile ("xsetbv"
                :
                : [_] "{eax}" (xcr0 & 0xe7),
                  [_] "{edx}" (0),
                  [_] "{ecx}" (0),
            );
        },
        .aarch64 => asm volatile (
            \\mrs x0, cpacr_el1
            \\orr x0, x0, #(3 << 20)
            \\orr x0, x0, #(3 << 16)
            \\orr x0, x0, #(3 << 24)
            \\msr cpacr_el1, x0
            \\isb
            \\mrs x1, id_aa64pfr0_el1
            \\ubfx x1, x1, #32, #4
            \\cbz x1, 1f
            \\mov x1, #3
            \\msr s3_0_c1_c2_0, x1
            \\1:
            \\mrs x1, id_aa64pfr1_el1
            \\ubfx x1, x1, #24, #4
            \\cbz x1, 2f
            \\mov x2, #3
            \\cmp x1, #2
            \\b.lo 3f
            \\orr x2, x2, #(1 << 30)
            \\3:
            \\msr s3_0_c1_c2_6, x2
            \\2:
            \\isb
            ::: .{ .x0 = true, .x1 = true, .x2 = true, .nzcv = true }),
        else => comptime unreachable,
    }
}

fn say(s: []const u8) void {
    for (s) |c| writeByte(c);
}

pub fn writeByte(c: u8) void {
    switch (sysinfo.target.cpu.arch) {
        .x86_64 => pc.writeByte(c),
        .aarch64 => {
            // `virt`'s PL011, which AAVMF set up and maps as device memory
            const uart: [*]volatile u32 = @ptrFromInt(0x0900_0000);
            // FR.TXFF
            while (uart[0x18 / 4] & 0x20 != 0) {}
            uart[0] = c;
        },
        else => comptime unreachable,
    }
}

/// Powers off; the status is lost under QEMU, the printed summary carries it.
pub fn exit(code: usize) noreturn {
    uefi.system_table.runtime_services.resetSystem(.shutdown, if (code == 0) .success else .aborted, null);
}
