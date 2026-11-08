//! The global vtable a test points the composed kernels at; level 2's work
//! space, exactly the buffer its contract promises (`buffers.md`) between
//! canaries; static operands; and an allocator handing out NaN-filled memory.

const std = @import("std");
const blas = @import("../__root__.zig");
const Vtable = @import("../Vtable.zig");
const portable = @import("../portable.zig");
const check = @import("check.zig");

const page = 4096;
/// A complex `f80` gemv's buffer past the x87 4096-row block.
const max_bytes = 136 * 1024;
var space: [2 * page + max_bytes]u8 align(page) = undefined;

/// The vtable `use` replaced; static, since a vtable is large for a rig's stack.
var saved: blas.TestSaved = undefined;
var in_use = false;

/// Points the global vtable, which the composed kernels call through, at `vt`
/// until `restore`. Not nested.
pub fn use(vt: *const Vtable) void {
    std.debug.assert(!in_use);
    in_use = true;
    saved = blas.testUse(vt);
}

pub fn restore() void {
    std.debug.assert(in_use);
    blas.testRestore(saved);
    in_use = false;
}

var portable_vt: Vtable = .{};
var portable_ready = false;

/// Every level's portable kernels, in static storage.
pub fn portableVtable() *const Vtable {
    if (!portable_ready) {
        portable.install(&portable_vt);
        portable_ready = true;
    }
    return &portable_vt;
}

/// What `open` puts in the buffer: finite junk, for kernels whose callers must
/// keep NaN out of it (x86_64's SSE `symv`, which reads past its copies), or
/// NaN, for the rest.
pub const Fill = enum { junk, nan };

/// A buffer of exactly `bytes`, `off` past a page boundary, filled as `with`
/// says (never zeros: a kernel must not rely on them) between canaries.
pub fn open(comptime F: type, bytes: usize, off: usize, with: Fill) [*]F {
    std.debug.assert(bytes + off <= max_bytes);
    @memset(&space, 0xaa);
    @memset(space[page + off ..][0..bytes], switch (with) {
        .junk => 0x41,
        .nan => 0xff,
    });
    return @ptrCast(@alignCast(&space[page + off]));
}

/// Whether the kernel stayed inside the buffer `open` gave for `bytes`.
pub fn intact(bytes: usize, off: usize) !void {
    for (space[0 .. page + off], 0..) |b, i| if (b != 0xaa) {
        check.say("buffer underflow at {d}\n", .{i});
        return error.TestUnexpectedResult;
    };
    for (space[page + off + bytes ..], 0..) |b, i| if (b != 0xaa) {
        check.say("buffer overflow at +{d}\n", .{i});
        return error.TestUnexpectedResult;
    };
}

/// Elements a vector of `n` at stride `inc` spans.
pub fn span(n: usize, inc: usize) usize {
    return (n - 1) * inc + 1;
}

/// The largest element the tests use: a complex `f80` or `f128` `XDouble`.
const max_elem = 32;

/// Operands for shapes too big for the stack: a matrix and its expected
/// value, vectors, and copies of them.
var store_a: [5120 * max_elem]u8 align(16) = undefined;
var store_w: [5120 * max_elem]u8 align(16) = undefined;
var store_x: [400 * max_elem]u8 align(16) = undefined;
var store_y: [400 * max_elem]u8 align(16) = undefined;
var store_x0: [400 * max_elem]u8 align(16) = undefined;
var store_y0: [400 * max_elem]u8 align(16) = undefined;

pub const Operands = enum { a, want, x, y, x0, y0 };

/// `n` elements of `T` of operand `which`.
pub fn operand(comptime T: type, which: Operands, n: usize) []T {
    const store: []align(16) u8 = switch (which) {
        .a => &store_a,
        .want => &store_w,
        .x => &store_x,
        .y => &store_y,
        .x0 => &store_x0,
        .y0 => &store_y0,
    };
    // whole elements only: `f80` on x86 is 12 bytes
    return std.mem.bytesAsSlice(T, store[0 .. store.len / @sizeOf(T) * @sizeOf(T)])[0..n];
}

/// An allocator over `child` whose fresh memory is all ones bits (NaN), so a
/// caller that must clear a buffer is caught when it does not, as is memory a
/// resize adds. Freed memory is set to the same first, so a read after free
/// gives NaN rather than the old values; a shrink's cut-off tail is left as is,
/// since the child may already have released it.
pub const Poison = struct {
    child: std.mem.Allocator,

    pub fn allocator(self: *Poison) std.mem.Allocator {
        return .{ .ptr = self, .vtable = &.{ .alloc = alloc, .resize = resize, .remap = remap, .free = free } };
    }

    fn alloc(ctx: *anyopaque, len: usize, alignment: std.mem.Alignment, ret: usize) ?[*]u8 {
        const self: *Poison = @ptrCast(@alignCast(ctx));
        const p = self.child.rawAlloc(len, alignment, ret) orelse return null;
        @memset(p[0..len], 0xff);
        return p;
    }
    fn resize(ctx: *anyopaque, memory: []u8, alignment: std.mem.Alignment, new_len: usize, ret: usize) bool {
        const self: *Poison = @ptrCast(@alignCast(ctx));
        if (!self.child.rawResize(memory, alignment, new_len, ret)) return false;
        if (new_len > memory.len) @memset(memory.ptr[memory.len..new_len], 0xff);
        return true;
    }
    fn remap(ctx: *anyopaque, memory: []u8, alignment: std.mem.Alignment, new_len: usize, ret: usize) ?[*]u8 {
        const self: *Poison = @ptrCast(@alignCast(ctx));
        const p = self.child.rawRemap(memory, alignment, new_len, ret) orelse return null;
        if (new_len > memory.len) @memset(p[memory.len..new_len], 0xff);
        return p;
    }
    fn free(ctx: *anyopaque, memory: []u8, alignment: std.mem.Alignment, ret: usize) void {
        const self: *Poison = @ptrCast(@alignCast(ctx));
        @memset(memory, 0xff);
        self.child.rawFree(memory, alignment, ret);
    }
};
