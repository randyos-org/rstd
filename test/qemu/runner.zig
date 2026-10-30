//! Freestanding test runner for `qemu-system-*`. Prints results through the
//! architecture's semihosting-style console and exits with the failure count
//! (0 = all passed), which the run step maps to the process exit code.
//!
//! Use with `zig test -target <arch>-freestanding --test-runner runner.zig`
//! and the matching `<arch>.ld`; `<arch>.zig` supplies the boot code. For a
//! `<arch>-uefi` target, `uefi.zig` does, and no linker script is used.

const std = @import("std");
const sysinfo = @import("builtin");

const is_uefi = sysinfo.target.os.tag == .uefi;

const arch = if (is_uefi) @import("uefi.zig") else switch (sysinfo.target.cpu.arch) {
    .aarch64 => @import("aarch64.zig"),
    .riscv64 => @import("riscv64.zig"),
    .arm => @import("arm.zig"),
    .loongarch64 => @import("loongarch64.zig"),
    .mips64el, .mips64 => @import("mips64.zig"),
    .mipsel, .mips => @import("mips32.zig"),
    .x86 => @import("x86.zig"),
    .x86_64 => @import("x86_64.zig"),
    .sparc64 => @import("sparc64.zig"),
    .s390x => @import("s390x.zig"),
    .powerpc64le, .powerpc64 => @import("ppc64.zig"),
    .powerpc => @import("ppc32.zig"),
    else => @compileError("qemutest: no boot code for " ++ @tagName(sysinfo.target.cpu.arch)),
};

/// UEFI's entry, which std's `EfiMain` calls; the bare rigs' boot code calls
/// `qemutest_main` itself.
pub const main = if (is_uefi) arch.main else {};

/// The message on the console, then a failed exit; a silent trap would just
/// end the log mid-test.
pub const panic = std.debug.FullPanic(struct {
    fn call(msg: []const u8, first: ?usize) noreturn {
        print("\nPANIC ");
        print(msg);
        // return addresses for `addr2line -e <kernel>`, by the frame pointers
        print(" at 0x");
        printHex(first orelse @returnAddress());
        // 64-bit POWER keeps the saved LR two words up, the back chain at 0
        const ra = if (sysinfo.target.cpu.arch.isPowerpc64()) 2 else 1;
        var fp = @frameAddress();
        for (0..12) |_| {
            if (fp == 0 or fp % @alignOf(usize) != 0) break;
            const frame: [*]const usize = @ptrFromInt(fp);
            if (frame[ra] == 0) break;
            print(" 0x");
            printHex(frame[ra]);
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
    // the cores `qemutest_smp` started, or one, and no threads to spawn;
    // `std.heap.SafeAllocator` asks for the count
    .Thread = struct {
        pub const Id = usize;
        pub const ThreadHandle = void;
        pub fn getCurrentId() Id {
            return if (qemutest_smp) |smp| smp.coreIndex() else 0;
        }
        pub fn getCpuCount() std.Thread.CpuCountError!usize {
            return if (qemutest_smp) |smp| smp.count() else 1;
        }
        pub fn yield() std.Thread.YieldError!void {
            std.atomic.spinLoopHint();
        }
        pub fn spawn(_: std.Thread.SpawnConfig, comptime _: anytype, _: anytype) std.Thread.SpawnError!@This() {
            return error.SystemResources;
        }
        pub fn getHandle(_: @This()) ThreadHandle {}
        pub fn detach(_: @This()) void {}
        pub fn join(_: @This()) void {}
    },
    // UEFI: rstd's `std.Io`, as RandyOS's loader has it, for the tests too
    .Io = if (is_uefi) arch.Io else null,
    .testing = if (is_uefi) struct {
        pub const allocator = std.testing.allocator_instance.allocator();
        pub const io = std_options_debug_io;
    } else null,
};

/// `std.testing` failure messages go through `std.debug.print`: `Io.failing`
/// (the default Io would drag in the threaded implementation), or UEFI's, with
/// stderr on the console. `failing` alone reaches `unreachable` there, so
/// every failed `expectEqual` was a panic.
pub const std_options_debug_io: std.Io = .{ .userdata = null, .vtable = &console_vtable };

const console_vtable = blk: {
    var vt = (if (is_uefi) arch.Io.io_inst else std.Io.failing).vtable.*;
    vt.swapCancelProtection = struct {
        fn f(_: ?*anyopaque, _: std.Io.CancelProtection) std.Io.CancelProtection {
            return .unblocked;
        }
    }.f;
    const stderr = struct {
        fn lock(_: ?*anyopaque, _: ?std.Io.Terminal.Mode) std.Io.Cancelable!std.Io.LockedStderr {
            return .{ .file_writer = &console_writer, .terminal_mode = .no_color };
        }
        fn tryLock(_: ?*anyopaque, _: ?std.Io.Terminal.Mode) std.Io.Cancelable!?std.Io.LockedStderr {
            return .{ .file_writer = &console_writer, .terminal_mode = .no_color };
        }
    };
    vt.lockStderr = stderr.lock;
    vt.tryLockStderr = stderr.tryLock;
    vt.unlockStderr = struct {
        fn f(_: ?*anyopaque) void {
            console_writer.interface.flush() catch {};
        }
    }.f;
    break :blk vt;
};

var console_writer: std.Io.File.Writer = .{
    .io = undefined,
    .file = undefined,
    .interface = .{ .vtable = &.{ .drain = consoleDrain }, .buffer = &.{} },
};

fn consoleDrain(w: *std.Io.Writer, data: []const []const u8, splat: usize) std.Io.Writer.Error!usize {
    print(w.buffered());
    w.end = 0;
    var n: usize = 0;
    for (data[0 .. data.len - 1]) |d| {
        print(d);
        n += d.len;
    }
    for (0..splat) |_| print(data[data.len - 1]);
    return n + data[data.len - 1].len * splat;
}

comptime {
    _ = arch;
}

/// A multithreaded kernel's other cores and their `Io` (`smp.zig`), which
/// the threaded tests use.
pub const qemutest_smp: ?type = if (sysinfo.single_threaded) null else @import("smp.zig");

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

// The stack, from each linker script or UEFI's boot code; it grows down
// from the top.
extern const __stack_bottom: u8;
extern const __stack_top: u8;

/// What unused stack holds, so an overflow and the high-water mark show.
const paint: usize = @truncate(@as(u64, 0x5a5a_c3c3_5a5a_c3c3));
/// Words at the bottom that must keep `paint`; past them is the data below.
const canary_words = 64;

fn stackWords() []usize {
    if (@hasDecl(arch, "stack")) return @ptrCast(arch.stack);
    const lo = std.mem.alignForward(usize, @intFromPtr(&__stack_bottom), @alignOf(usize));
    const p: [*]usize = @ptrFromInt(lo);
    return p[0 .. (@intFromPtr(&__stack_top) - lo) / @sizeOf(usize)];
}

/// Paints the stack below this frame, less a page of margin. Never inlined:
/// inlined, the frame address is the caller's top, and a caller's frame past
/// the margin (`qemutest_main`'s in release builds) gets painted.
noinline fn paintStack() void {
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

/// A test's heap use, when it is a good part of the heap.
fn printHeap(used: usize) void {
    if (used < heap.max_len / 4) return;
    print(" (heap ");
    printInt(used);
    print(")");
}

/// The request that ran the heap out, where one did.
fn printHeapFailure(failed: @TypeOf(heap.failed)) void {
    const f = failed orelse return;
    print(" (heap: ");
    printInt(f.len);
    print(" bytes wanted, ");
    printInt(f.in_use);
    print(" in use, largest free ");
    printInt(f.largest);
    print(")");
}

const heap = @import("heap.zig");

/// The heap, where the boot code cannot get one from the firmware
/// (`heapRegion`): sparc64's maps only 4 MB of the image; UEFI's.
var static_heap: [if (@hasDecl(arch, "heapRegion")) 0 else heap.max_len]u8 align(heap.max_align) = undefined;

export fn qemutest_main() callconv(.c) noreturn {
    print("qemutest ");
    print(@tagName(sysinfo.target.cpu.arch));
    print("\n");
    if (qemutest_smp) |smp| smp.start();
    heap.use(if (@hasDecl(arch, "heapRegion")) arch.heapRegion(heap.max_len, heap.max_align) else &static_heap);
    paintStack();
    var pass: usize = 0;
    var fail: usize = 0;
    var skip: usize = 0;
    var peak: usize = 0;
    var heap_peak: usize = 0;
    for (sysinfo.test_functions) |t| {
        _ = heap.resetPeak();
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
        const heap_failed = heap.failed;
        const heap_used = heap.resetPeak();
        heap_peak = @max(heap_peak, heap_used);
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
                printHeap(heap_used);
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
                printHeapFailure(heap_failed);
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
    print("heap high water ");
    printInt(heap_peak);
    print(" of ");
    printInt(heap.len());
    print(" bytes\n");
    arch.exit(if (fail == 0) 0 else 1);
}
