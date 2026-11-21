//! An allocator that keeps what a routine's level-3 calls free for its next
//! ones, and frees it all at the end: upstream takes its packing space once
//! per interface call (`blas_memory_alloc` in `interface/lapack/potrf.c`),
//! where each of the drivers here asks its allocator. Fresh pages cost a
//! threaded `potrf` half its rate on Windows (issue #546, `threading.md`).
//!
//! One thread only: the drivers allocate and free on the calling thread,
//! before their workers start and after they end. Nothing resizes in place.

const std = @import("std");

const Allocator = std.mem.Allocator;
const Alignment = std.mem.Alignment;
const Kept = @This();

child: Allocator,
/// Blocks out on loan and blocks kept, by their whole length; `ptr` null
/// for an empty slot.
out: [slots]Block = @splat(.{}),
free: [slots]Block = @splat(.{}),

/// More than one call has live at once (`threaded.zig`: the packing space
/// and the barriers).
const slots = 4;

const Block = struct {
    ptr: ?[*]u8 = null,
    len: usize = 0,
    alignment: Alignment = .@"1",
};

pub fn init(child: Allocator) Kept {
    return .{ .child = child };
}

/// Frees what is kept; everything lent must be back.
pub fn deinit(k: *Kept) void {
    for (&k.free) |*b| if (b.ptr) |p| {
        k.child.rawFree(p[0..b.len], b.alignment, @returnAddress());
        b.* = .{};
    };
    for (k.out) |b| std.debug.assert(b.ptr == null);
}

pub fn allocator(k: *Kept) Allocator {
    return .{ .ptr = k, .vtable = &.{ .alloc = alloc, .resize = Allocator.noResize, .remap = Allocator.noRemap, .free = release } };
}

fn alloc(ctx: *anyopaque, len: usize, alignment: Alignment, ret: usize) ?[*]u8 {
    const k: *Kept = @ptrCast(@alignCast(ctx));
    const spot = for (&k.out) |*o| {
        if (o.ptr == null) break o;
    } else return k.child.rawAlloc(len, alignment, ret);
    // the smallest kept block that holds it
    var best: ?*Block = null;
    for (&k.free) |*b| {
        if (b.ptr == null or b.len < len or b.alignment.compare(.lt, alignment)) continue;
        if (best == null or b.len < best.?.len) best = b;
    }
    if (best) |b| {
        spot.* = b.*;
        b.* = .{};
        return spot.ptr;
    }
    const p = k.child.rawAlloc(len, alignment, ret) orelse return null;
    spot.* = .{ .ptr = p, .len = len, .alignment = alignment };
    return p;
}

fn release(ctx: *anyopaque, memory: []u8, alignment: Alignment, ret: usize) void {
    const k: *Kept = @ptrCast(@alignCast(ctx));
    const lent = for (&k.out) |*o| {
        if (o.ptr == memory.ptr) break o;
    } else return k.child.rawFree(memory, alignment, ret);
    const block = lent.*;
    lent.* = .{};
    // into an empty slot, else in place of the smallest if that is smaller
    var into: *Block = &k.free[0];
    for (&k.free) |*b| {
        if (b.ptr == null) {
            into = b;
            break;
        }
        if (b.len < into.len) into = b;
    }
    if (into.ptr) |p| {
        if (into.len >= block.len) return k.child.rawFree(block.ptr.?[0..block.len], block.alignment, ret);
        k.child.rawFree(p[0..into.len], into.alignment, ret);
    }
    into.* = block;
}

// ---- tests -------------------------------------------------------------------------

const testing = std.testing;

test "a freed block serves the next request it holds; deinit frees what is kept" {
    var k: Kept = .init(testing.allocator);
    defer k.deinit();
    const al = k.allocator();
    const a = try al.alignedAlloc(u8, .@"64", 4096);
    const first = a.ptr;
    al.free(a);
    // smaller, and a weaker alignment: the same block
    const b = try al.alloc(u8, 100);
    try testing.expectEqual(@as([*]u8, first), b.ptr);
    // a second live block is new
    const c = try al.alloc(u8, 100);
    try testing.expect(c.ptr != b.ptr);
    al.free(b);
    al.free(c);
    // larger than anything kept: new, and both kept blocks stay
    const d = try al.alloc(u8, 1 << 16);
    try testing.expect(d.ptr != first);
    al.free(d);
    // the smallest that holds it
    const e = try al.alloc(u8, 50);
    defer al.free(e);
    try testing.expectEqual(@as(usize, 100), k.out[0].len);
}

test "past its slots it passes requests through and drops the smallest kept block" {
    var k: Kept = .init(testing.allocator);
    defer k.deinit();
    const al = k.allocator();
    var held: [slots + 2][]u8 = undefined;
    for (&held, 0..) |*h, i| h.* = try al.alloc(u8, @as(usize, 64) << @intCast(i));
    for (k.out) |o| try testing.expect(o.ptr != null);
    for (held) |h| al.free(h);
    var kept: usize = 0;
    for (k.free) |b| kept += b.len;
    try testing.expectEqual(@as(usize, 64 + 128 + 256 + 512), kept);
    // six more, each larger than every kept block: the kept ones go
    for (&held, 0..) |*h, i| h.* = try al.alloc(u8, @as(usize, 4096) << @intCast(i));
    for (held) |h| al.free(h);
    for (k.free) |b| try testing.expect(b.len >= 4096);
}
