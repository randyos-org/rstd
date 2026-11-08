const std = @import("std");
const zigstd = @import("../zigstd.zig");
const Allocator = zigstd.Allocator;

// A guillotine-style rectangle bin packer: starts from one free space (the
// whole canvas) and, for each rect placed, splits the space it landed in
// into up to two new free spaces along the leftover width/height; no
// rotation, no multiple bins, just placement + splitting. General-purpose
// (not font-specific); used here to lay out a glyph atlas bitmaps.

pub const Rect = struct {
    w: i32,
    h: i32,
    // x/y are overwritten during packing and are effectively return values.
    x: i32 = 0,
    y: i32 = 0,

    pub inline fn area(self: Rect) i32 {
        return self.w * self.h;
    }

    pub inline fn perimeter(self: Rect) i32 {
        return 2 * (self.w + self.h);
    }
};

/// Only useful if you don't want to track which output rect corresponds to
/// which input yourself.
pub const IdRect = struct {
    id: i32,
    rect: Rect,

    pub inline fn area(self: IdRect) i32 {
        return self.rect.area();
    }

    pub inline fn perimeter(self: IdRect) i32 {
        return self.rect.perimeter();
    }
};

pub const Context = struct {
    allocator: Allocator,
    list: std.ArrayList(Rect),
    w: i32,
    h: i32,

    /// spaces_to_prealloc should be set to 2x rects supplied to pack() until
    /// deinit if setting assume_capacity to true in pack(). Lower values
    /// will likely work, but are not guaranteed to.
    pub fn create(allocator: Allocator, w: i32, h: i32, opts: struct { spaces_to_prealloc: u32 = 0 }) Allocator.Error!Context {
        var list: std.ArrayList(Rect) = if (opts.spaces_to_prealloc > 0)
            try .initCapacity(allocator, opts.spaces_to_prealloc)
        else
            .empty;
        try list.append(allocator, .{ .w = w, .h = h, .x = 0, .y = 0 });
        return .{ .allocator = allocator, .list = list, .w = w, .h = h };
    }

    /// Capacity stays intact, so it's safe to use assume_capacity again if
    /// it was supplied in create().
    pub fn clear(self: *Context) Allocator.Error!void {
        self.list.clearRetainingCapacity();
        try self.list.append(self.allocator, .{ .w = self.w, .h = self.h, .x = 0, .y = 0 });
    }

    pub fn areaFree(self: Context) f32 {
        const total_space: f32 = @floatFromInt(self.w * self.h);
        var free_space: f32 = 0.0;
        for (self.list.items) |space| free_space += @floatFromInt(space.w * space.h);
        return if (free_space == 0.0) return 0.0 else free_space / total_space;
    }

    pub fn deinit(self: *Context) void {
        self.list.deinit(self.allocator);
    }
};

inline fn append(
    comptime assume_capacity: bool,
    allocator: Allocator,
    noalias list: *std.ArrayList(Rect),
    rect: Rect,
) Allocator.Error!void {
    if (assume_capacity)
        list.appendAssumeCapacity(rect)
    else
        try list.append(allocator, rect);
}

/// Sorting is highly recommended for most real-world inputs; it tends to
/// pack tighter and faster when the rects have high w/h entropy. Good
/// metrics to sort on (descending), roughly in order of usefulness:
/// max(w, h), then min(w, h), then area, then perimeter.
pub fn pack(
    comptime RectType: type,
    noalias ctx: *Context,
    noalias rects: []RectType,
    comptime opts: struct {
        assume_capacity: bool = false,
        sortLessThanFn: ?fn (context: void, lhs: RectType, rhs: RectType) bool = null,
    },
) (error{NoSpaceLeft} || Allocator.Error)!void {
    if (RectType != Rect and RectType != IdRect)
        @compileError("Invalid rect type. Valid ones are 'Rect' and 'IdRect'.");

    if (opts.sortLessThanFn) |lessThanFn| std.sort.pdq(RectType, rects, {}, lessThanFn);

    rectLoop: for (rects) |*any_rect| {
        var rect = if (RectType == IdRect) &any_rect.rect else any_rect;

        if (ctx.list.items.len == 0) return error.NoSpaceLeft;

        var iter = std.mem.reverseIterator(ctx.list.items);
        var i: usize = ctx.list.items.len -% 1;
        while (iter.next()) |space| : (i -%= 1) {
            const free_w = space.w - rect.w;
            const free_h = space.h - rect.h;
            if (free_w < 0 or free_h < 0) continue;

            defer {
                rect.x = space.x;
                rect.y = space.y;
                _ = ctx.list.orderedRemove(i);
            }

            if (free_w == 0 and free_h == 0) continue :rectLoop;

            if (free_w > 0 and free_h == 0) {
                try append(
                    opts.assume_capacity,
                    ctx.allocator,
                    &ctx.list,
                    .{ .x = space.x + rect.w, .y = space.y, .w = space.w - rect.w, .h = space.h },
                );
                continue :rectLoop;
            }

            if (free_w == 0 and free_h > 0) {
                try append(
                    opts.assume_capacity,
                    ctx.allocator,
                    &ctx.list,
                    .{ .x = space.x, .y = space.y + rect.h, .w = space.w, .h = space.h - rect.h },
                );
                continue :rectLoop;
            }

            if (free_w > free_h) {
                try append(
                    opts.assume_capacity,
                    ctx.allocator,
                    &ctx.list,
                    .{ .x = space.x + rect.w, .y = space.y, .w = free_w, .h = space.h },
                );
                try append(
                    opts.assume_capacity,
                    ctx.allocator,
                    &ctx.list,
                    .{ .x = space.x, .y = space.y + rect.h, .w = rect.w, .h = free_h },
                );
                continue :rectLoop;
            }

            try append(
                opts.assume_capacity,
                ctx.allocator,
                &ctx.list,
                .{ .x = space.x, .y = space.y + rect.h, .w = space.w, .h = free_h },
            );
            try append(
                opts.assume_capacity,
                ctx.allocator,
                &ctx.list,
                .{ .x = space.x + rect.w, .y = space.y, .w = free_w, .h = rect.h },
            );
            continue :rectLoop;
        }

        return error.NoSpaceLeft;
    }
}

fn rectsOverlap(a: Rect, b: Rect) bool {
    return a.x < b.x + b.w and a.x + a.w > b.x and a.y < b.y + b.h and a.y + a.h > b.y;
}

test "pack places non-overlapping rects that all fit within the canvas" {
    const alloc = zigstd.testAllocator;
    var ctx = try Context.create(alloc, 64, 64, .{});
    defer ctx.deinit();

    var rects = [_]Rect{
        .{ .w = 20, .h = 30 },
        .{ .w = 40, .h = 10 },
        .{ .w = 15, .h = 15 },
    };
    try pack(Rect, &ctx, &rects, .{ .sortLessThanFn = struct {
        fn lessThan(_: void, a: Rect, b: Rect) bool {
            return @max(a.w, a.h) > @max(b.w, b.h);
        }
    }.lessThan });

    for (rects) |r| {
        try zigstd.testExpectTrue(r.x >= 0 and r.y >= 0);
        try zigstd.testExpectTrue(r.x + r.w <= 64 and r.y + r.h <= 64);
    }
    for (0..rects.len) |i| {
        for (i + 1..rects.len) |j| {
            try zigstd.testExpectTrue(!rectsOverlap(rects[i], rects[j]));
        }
    }
}

test "pack reports NoSpaceLeft when a rect can't possibly fit" {
    const alloc = zigstd.testAllocator;
    var ctx = try Context.create(alloc, 16, 16, .{});
    defer ctx.deinit();

    var rects = [_]Rect{.{ .w = 32, .h = 32 }};
    try zigstd.testExpectError(error.NoSpaceLeft, pack(Rect, &ctx, &rects, .{}));
}
