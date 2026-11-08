//! The rigs' `std.heap.page_allocator`: a region the runner hands over (`use`; static, or from the firmware) in whole pages, tracked
//! by a bitmap so freed pages are reused, as a test cycling through every core needs.

const std = @import("std");
const sysinfo = @import("builtin");

/// The page size the rigs use, on every target.
pub const page_size = 4096;
const page = page_size;
/// Page indices align as addresses do, up to this.
pub const max_align = 1 << 16;
/// The largest region, and what the runner asks for: on x86_64 and aarch64,
/// room for a core's whole gemm packing buffer (T132; 20 MB on AVX) beside a
/// test's operands, threaded too.
pub const max_len = switch (sysinfo.target.cpu.arch) {
    .x86_64, .aarch64 => 64 << 20,
    else => 32 << 20,
};

var region: []align(max_align) u8 = &.{};
var used: std.bit_set.Static(max_len / page) = .empty;

/// Pages in use, and the most since `resetPeak`.
var in_use: usize = 0;
var peak: usize = 0;
/// The first request that failed since `resetPeak`: its bytes, the bytes in
/// use and the largest free run then.
pub var failed: ?struct { len: usize, in_use: usize, largest: usize } = null;

/// Bytes of the region.
pub fn len() usize {
    return pages() * page;
}

/// The most bytes in use since the last call, which starts a new count.
pub fn resetPeak() usize {
    defer {
        peak = in_use;
        failed = null;
    }
    return peak * page;
}

fn take(i: usize, n: usize, v: bool) void {
    used.setRangeValue(.{ .start = i, .end = i + n }, v);
    in_use = if (v) in_use + n else in_use - n;
    peak = @max(peak, in_use);
}

/// The longest run of free pages, in bytes.
fn largest() usize {
    var best: usize = 0;
    var run: usize = 0;
    for (0..pages()) |j| {
        run = if (used.isSet(j)) 0 else run + 1;
        best = @max(best, run);
    }
    return best * page;
}

/// Hands out `r` from now on; before, every allocation fails.
pub fn use(r: []align(max_align) u8) void {
    region = r[0..@min(r.len, max_len)];
}

/// Held through each call where cores share the heap (`smp.zig`).
var lock: std.atomic.Mutex = .unlocked;

fn acquire() void {
    if (!sysinfo.single_threaded) while (!lock.tryLock()) std.atomic.spinLoopHint();
}

fn release() void {
    if (!sysinfo.single_threaded) lock.unlock();
}

pub const allocator: std.mem.Allocator = .{ .ptr = undefined, .vtable = &.{
    .alloc = alloc,
    .resize = resize,
    .remap = remap,
    .free = free,
} };

fn pages() usize {
    return region.len / page;
}

fn count(bytes: usize) usize {
    return (bytes + page - 1) / page;
}

fn first(p: [*]u8) usize {
    return (@intFromPtr(p) - @intFromPtr(region.ptr)) / page;
}

/// The first `n` free pages starting at a multiple of `step`.
fn findRun(n: usize, step: usize) ?usize {
    var i: usize = 0;
    outer: while (i + n <= pages()) {
        for (i..i + n) |j| if (used.isSet(j)) {
            i = std.mem.alignForward(usize, j + 1, step);
            continue :outer;
        };
        return i;
    }
    return null;
}

fn alloc(_: *anyopaque, bytes: usize, alignment: std.mem.Alignment, _: usize) ?[*]u8 {
    acquire();
    defer release();
    if (alignment.toByteUnits() > max_align) return null;
    const n = count(bytes);
    const i = findRun(n, @max(1, alignment.toByteUnits() / page)) orelse {
        if (failed == null) failed = .{ .len = bytes, .in_use = in_use * page, .largest = largest() };
        return null;
    };
    take(i, n, true);
    return region[i * page ..].ptr;
}

fn resize(_: *anyopaque, memory: []u8, _: std.mem.Alignment, new_len: usize, _: usize) bool {
    acquire();
    defer release();
    const i = first(memory.ptr);
    const old = count(memory.len);
    const new = count(new_len);
    if (new <= old) {
        take(i + new, old - new, false);
        return true;
    }
    if (i + new > pages()) return false;
    for (i + old..i + new) |j| if (used.isSet(j)) return false;
    take(i + old, new - old, true);
    return true;
}

fn remap(ctx: *anyopaque, memory: []u8, alignment: std.mem.Alignment, new_len: usize, ret_addr: usize) ?[*]u8 {
    return if (resize(ctx, memory, alignment, new_len, ret_addr)) memory.ptr else null;
}

fn free(_: *anyopaque, memory: []u8, _: std.mem.Alignment, _: usize) void {
    acquire();
    defer release();
    const i = first(memory.ptr);
    take(i, count(memory.len), false);
}
