//! The rigs' page allocator, `std.heap.page_allocator` there: a region
//! handed out in whole pages and tracked by a bitmap, so freed pages are
//! reused. One test cycling through every core frees as much as it takes.
//! The runner gives it the region (`use`): a static buffer, or one the boot
//! code gets from the firmware.

const std = @import("std");
/// The page size the rigs use, on every target.
pub const page_size = 4096;
const page = page_size;
/// Page indices align as addresses do, up to this.
pub const max_align = 1 << 16;
/// The largest region, and what the runner asks for.
pub const max_len = 32 << 20;

var region: []align(max_align) u8 = &.{};
var used: std.bit_set.Static(max_len / page) = .empty;

/// Hands out `r` from now on; before, every allocation fails.
pub fn use(r: []align(max_align) u8) void {
    region = r[0..@min(r.len, max_len)];
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

fn count(len: usize) usize {
    return (len + page - 1) / page;
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
    if (i + new > pages()) return false;
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
