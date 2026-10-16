//! Freestanding test runner for `qemu-system-*`. Prints results through the
//! architecture's semihosting-style console and exits with the failure count
//! (0 = all passed), which the run step maps to the process exit code.
//!
//! Use with `zig test -target <arch>-freestanding --test-runner runner.zig`
//! and the matching `<arch>.ld`; `<arch>.zig` supplies the boot code.

const std = @import("std");
const sysinfo = @import("builtin");

const arch = switch (sysinfo.cpu.arch) {
    .aarch64 => @import("aarch64.zig"),
    .riscv64 => @import("riscv64.zig"),
    .arm => @import("arm.zig"),
    .loongarch64 => @import("loongarch64.zig"),
    .mips64el => @import("mips64.zig"),
    .mipsel => @import("mips32.zig"),
    .x86 => @import("x86.zig"),
    .x86_64 => @import("x86_64.zig"),
    .sparc64 => @import("sparc64.zig"),
    .s390x => @import("s390x.zig"),
    .powerpc64le => @import("ppc64le.zig"),
    else => @compileError("qemutest: no boot code for " ++ @tagName(sysinfo.cpu.arch)),
};

/// The message on the console, then a failed exit; a silent trap would just
/// end the log mid-test.
pub const panic = std.debug.FullPanic(struct {
    fn call(msg: []const u8, first: ?usize) noreturn {
        print("\nPANIC ");
        print(msg);
        // return addresses for `addr2line -e <kernel>`, by the frame pointers
        print(" at 0x");
        printHex(first orelse @returnAddress());
        var fp = @frameAddress();
        for (0..12) |_| {
            if (fp == 0 or fp % @alignOf(usize) != 0) break;
            const frame: [*]const usize = @ptrFromInt(fp);
            if (frame[1] == 0) break;
            print(" 0x");
            printHex(frame[1]);
            if (frame[0] <= fp) break;
            fp = frame[0];
        }
        print("\n");
        arch.exit(1);
    }
}.call);

/// Freestanding targets have no page size to infer, and no page allocator;
/// `heap.zig` is one, under `std.testing.allocator`.
pub const std_os_options: std.Options.OperatingSystem = .{
    .heap = struct {
        pub const page_size_min = heap.page_size;
        pub const page_size_max = heap.page_size;
        pub const page_allocator = @import("heap.zig").allocator;
        pub fn queryPageSize() usize {
            return heap.page_size;
        }
    },
    // one CPU and no threads; `std.heap.SafeAllocator` asks for the count
    .Thread = struct {
        pub const Id = usize;
        pub const ThreadHandle = void;
        pub fn getCurrentId() Id {
            return 0;
        }
        pub fn getCpuCount() std.Thread.CpuCountError!usize {
            return 1;
        }
        pub fn yield() std.Thread.YieldError!void {}
        pub fn spawn(_: std.Thread.SpawnConfig, comptime _: anytype, _: anytype) std.Thread.SpawnError!@This() {
            return error.SystemResources;
        }
        pub fn getHandle(_: @This()) ThreadHandle {}
        pub fn detach(_: @This()) void {}
        pub fn join(_: @This()) void {}
    },
};

/// `std.testing` failure messages go through `std.debug.print`; there is no OS
/// to print to, and the default Io would drag in the threaded implementation.
pub const std_options_debug_io: std.Io = std.Io.failing;

comptime {
    _ = arch;
}

fn print(s: []const u8) void {
    for (s) |c| arch.writeByte(c);
}

/// The console, for tests with something to say; there is no stderr here.
pub fn testWrite(s: []const u8) void {
    print(s);
}

fn printInt(v: usize) void {
    var buf: [20]u8 = undefined;
    var i: usize = buf.len;
    var n = v;
    while (true) {
        i -= 1;
        buf[i] = '0' + @as(u8, @intCast(n % 10));
        n /= 10;
        if (n == 0) break;
    }
    print(buf[i..]);
}

fn printHex(v: usize) void {
    var buf: [@sizeOf(usize) * 2]u8 = undefined;
    for (&buf, 0..) |*c, i| c.* = "0123456789abcdef"[(v >> @intCast((buf.len - 1 - i) * 4)) & 0xf];
    print(&buf);
}

// The stack, from each linker script; it grows down from `__stack_top`.
extern const __stack_bottom: u8;
extern const __stack_top: u8;

/// What unused stack holds, so an overflow and the high-water mark show.
const paint: usize = @truncate(@as(u64, 0x5a5a_c3c3_5a5a_c3c3));
/// Words at the bottom that must keep `paint`; past them is the data below.
const canary_words = 64;

fn stackWords() []usize {
    const lo = std.mem.alignForward(usize, @intFromPtr(&__stack_bottom), @alignOf(usize));
    const p: [*]usize = @ptrFromInt(lo);
    return p[0 .. (@intFromPtr(&__stack_top) - lo) / @sizeOf(usize)];
}

/// Paints the stack below this frame, less a page of margin.
fn paintStack() void {
    const words = stackWords();
    const below = (@frameAddress() -| 4096 -| @intFromPtr(words.ptr)) / @sizeOf(usize);
    @memset(words[0..@min(below, words.len)], paint);
}

fn stackIntact() bool {
    for (stackWords()[0..canary_words]) |w| if (w != paint) return false;
    return true;
}

/// The most stack used since the last `paintStack`, in bytes.
fn stackHighWater() usize {
    const words = stackWords();
    for (words, 0..) |w, i| if (w != paint) return (words.len - i) * @sizeOf(usize);
    return 0;
}

/// A test's stack use, when it is big enough to matter on the smaller stacks.
fn printStack(used: usize) void {
    if (used < 64 << 10) return;
    print(" (stack ");
    printInt(used);
    print(")");
}

const heap = @import("heap.zig");

/// The heap, where the boot code cannot get one from the firmware
/// (`heapRegion`): sparc64's maps only 4 MB of the image.
var static_heap: [if (@hasDecl(arch, "heapRegion")) 0 else heap.max_len]u8 align(heap.max_align) = undefined;

export fn qemutest_main() callconv(.c) noreturn {
    print("qemutest ");
    print(@tagName(sysinfo.cpu.arch));
    print("\n");
    heap.use(if (@hasDecl(arch, "heapRegion")) arch.heapRegion(heap.max_len, heap.max_align) else &static_heap);
    paintStack();
    var pass: usize = 0;
    var fail: usize = 0;
    var skip: usize = 0;
    var peak: usize = 0;
    for (sysinfo.test_functions) |t| {
        // as std's runner does: a fresh testing allocator, leaks fail the test
        std.testing.allocator_instance = .init(std.heap.page_allocator, .{});
        print(t.name);
        print("... ");
        const result = t.func();
        if (!stackIntact()) {
            // what lies below the stack is overwritten; nothing after is trustworthy
            print("FAIL stack overflow (");
            printInt(stackWords().len * @sizeOf(usize));
            print(" byte stack)\n");
            arch.exit(1);
        }
        const used = stackHighWater();
        peak = @max(peak, used);
        paintStack();
        const leaks = std.testing.allocator_instance.deinit();
        if (result) |_| {
            if (leaks != 0) {
                fail += 1;
                print("FAIL leaked ");
                printInt(leaks);
                print("\n");
            } else {
                pass += 1;
                print("ok");
                printStack(used);
                print("\n");
            }
        } else |err| switch (err) {
            error.SkipZigTest => {
                skip += 1;
                print("skip\n");
            },
            else => {
                fail += 1;
                print("FAIL ");
                print(@errorName(err));
                print("\n");
            },
        }
    }
    printInt(pass);
    print(" passed, ");
    printInt(skip);
    print(" skipped, ");
    printInt(fail);
    print(" failed\n");
    print("stack high water ");
    printInt(peak);
    print(" of ");
    printInt(stackWords().len * @sizeOf(usize));
    print(" bytes\n");
    arch.exit(if (fail == 0) 0 else 1);
}
