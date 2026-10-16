//! The rigs' page allocator, `std.heap.page_allocator` there: a static region
//! handed out in whole pages and tracked by a bitmap, so freed pages are
//! reused. One test cycling through every core frees as much as it takes.

const std = @import("std");
const sysinfo = @import("builtin");
const mem = @import("../mem.zig");

const page = mem.page_size;
/// sparc64's OpenBIOS maps only about 4 MB of the image.
const pages = (if (sysinfo.cpu.arch == .sparc64) 1 << 20 else 32 << 20) / page;

/// Page indices then align as addresses do, up to this.
const max_align = 1 << 16;
var region: [pages * page]u8 align(max_align) = undefined;
var used: std.StaticBitSet(pages) = .empty;

pub const allocator: std.mem.Allocator = .{ .ptr = undefined, .vtable = &.{
    .alloc = alloc,
    .resize = resize,
    .remap = remap,
    .free = free,
} };

fn count(len: usize) usize {
    return (len + page - 1) / page;
}

fn first(p: [*]u8) usize {
    return (@intFromPtr(p) - @intFromPtr(&region)) / page;
}

/// The first `n` free pages starting at a multiple of `step`.
fn findRun(n: usize, step: usize) ?usize {
    var i: usize = 0;
    outer: while (i + n <= pages) {
        for (i..i + n) |j| if (used.isSet(j)) {
            i = std.mem.alignForward(usize, j + 1, step);
            continue :outer;
        };
        return i;
    }
    return null;
}

fn alloc(_: *anyopaque, len: usize, alignment: std.mem.Alignment, _: usize) ?[*]u8 {
    if (alignment.toByteUnits() > max_align) return null;
    const n = count(len);
    const i = findRun(n, @max(1, alignment.toByteUnits() / page)) orelse return null;
    used.setRangeValue(.{ .start = i, .end = i + n }, true);
    return region[i * page ..].ptr;
}

fn resize(_: *anyopaque, memory: []u8, _: std.mem.Alignment, new_len: usize, _: usize) bool {
    const i = first(memory.ptr);
    const old = count(memory.len);
    const new = count(new_len);
    if (new <= old) {
        used.setRangeValue(.{ .start = i + new, .end = i + old }, false);
        return true;
    }
    if (i + new > pages) return false;
    for (i + old..i + new) |j| if (used.isSet(j)) return false;
    used.setRangeValue(.{ .start = i + old, .end = i + new }, true);
    return true;
}

fn remap(ctx: *anyopaque, memory: []u8, alignment: std.mem.Alignment, new_len: usize, ret_addr: usize) ?[*]u8 {
    return if (resize(ctx, memory, alignment, new_len, ret_addr)) memory.ptr else null;
}

fn free(_: *anyopaque, memory: []u8, _: std.mem.Alignment, _: usize) void {
    const i = first(memory.ptr);
    used.setRangeValue(.{ .start = i, .end = i + count(memory.len) }, false);
}
