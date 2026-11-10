const std = @import("std");
const zigstd = @import("../zigstd.zig");
const Allocator = zigstd.Allocator;
const geometry = @import("geometry.zig");

// A guillotine-style rectangle bin packer: starts from one free space (the
// whole canvas) and, for each rect placed, splits the space it landed in
// into up to two new free spaces along the leftover width/height; no
// rotation, no multiple bins, just placement + splitting. General-purpose
// (not font-specific); used here to lay out a glyph atlas bitmaps.

/// A rect to place: `pack` writes its `x` and `y`.
pub const Rect = geometry.Rect;

/// Only useful if you don't want to track which output rect corresponds to
/// which input yourself.
pub const IdRect = struct {
    id: i32,
    rect: Rect,

    pub inline fn area(self: IdRect) u64 {
        return self.rect.area();
    }

    pub inline fn perimeter(self: IdRect) u64 {
        return self.rect.perimeter();
    }
};

pub const Context = struct {
    allocator: Allocator,
    list: std.ArrayList(Rect),
    w: u32,
    h: u32,

    /// spaces_to_prealloc should be set to 2x rects supplied to pack() until
    /// deinit if setting assume_capacity to true in pack(). Lower values
    /// will likely work, but are not guaranteed to.
    pub fn init(allocator: Allocator, w: u32, h: u32, opts: struct { spaces_to_prealloc: u32 = 0 }) Allocator.Error!Context {
        var list: std.ArrayList(Rect) = if (opts.spaces_to_prealloc > 0)
            try .initCapacity(allocator, opts.spaces_to_prealloc)
        else
            .empty;
        try list.append(allocator, .{ .w = w, .h = h, .x = 0, .y = 0 });
        return .{ .allocator = allocator, .list = list, .w = w, .h = h };
    }

    /// Capacity stays intact, so it's safe to use assume_capacity again if
    /// it was supplied in init().
    pub fn clear(self: *Context) Allocator.Error!void {
        self.list.clearRetainingCapacity();
        try self.list.append(self.allocator, .{ .w = self.w, .h = self.h, .x = 0, .y = 0 });
    }

    pub fn areaFree(self: Context) f32 {
        const total_space: f32 = @floatFromInt(@as(u64, self.w) * self.h);
        var free_space: f32 = 0.0;
        for (self.list.items) |space| free_space += @floatFromInt(space.area());
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
            if (rect.w > space.w or rect.h > space.h) continue;
            const free_w = space.w - rect.w;
            const free_h = space.h - rect.h;
            // where the free spaces right of and below `rect` begin
            const right = space.x + @as(i32, @intCast(rect.w));
            const below = space.y + @as(i32, @intCast(rect.h));

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
                    .{ .x = right, .y = space.y, .w = free_w, .h = space.h },
                );
                continue :rectLoop;
            }

            if (free_w == 0 and free_h > 0) {
                try append(
                    opts.assume_capacity,
                    ctx.allocator,
                    &ctx.list,
                    .{ .x = space.x, .y = below, .w = space.w, .h = free_h },
                );
                continue :rectLoop;
            }

            if (free_w > free_h) {
                try append(
                    opts.assume_capacity,
                    ctx.allocator,
                    &ctx.list,
                    .{ .x = right, .y = space.y, .w = free_w, .h = space.h },
                );
                try append(
                    opts.assume_capacity,
                    ctx.allocator,
                    &ctx.list,
                    .{ .x = space.x, .y = below, .w = rect.w, .h = free_h },
                );
                continue :rectLoop;
            }

            try append(
                opts.assume_capacity,
                ctx.allocator,
                &ctx.list,
                .{ .x = space.x, .y = below, .w = space.w, .h = free_h },
            );
            try append(
                opts.assume_capacity,
                ctx.allocator,
                &ctx.list,
                .{ .x = right, .y = space.y, .w = free_w, .h = rect.h },
            );
            continue :rectLoop;
        }

        return error.NoSpaceLeft;
    }
}

fn rectsOverlap(a: Rect, b: Rect) bool {
    const a1 = a.bottomRight();
    const b1 = b.bottomRight();
    return a.x < b1.x and a1.x > b.x and a.y < b1.y and a1.y > b.y;
}

test "pack places non-overlapping rects that all fit within the canvas" {
    const alloc = zigstd.testAllocator;
    var ctx = try Context.init(alloc, 64, 64, .{});
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
        try zigstd.testExpectTrue(r.bottomRight().x <= 64 and r.bottomRight().y <= 64);
    }
    for (0..rects.len) |i| {
        for (i + 1..rects.len) |j| {
            try zigstd.testExpectTrue(!rectsOverlap(rects[i], rects[j]));
        }
    }
}

test "pack reports NoSpaceLeft when a rect can't possibly fit" {
    const alloc = zigstd.testAllocator;
    var ctx = try Context.init(alloc, 16, 16, .{});
    defer ctx.deinit();

    var rects = [_]Rect{.{ .w = 32, .h = 32 }};
    try zigstd.testExpectError(error.NoSpaceLeft, pack(Rect, &ctx, &rects, .{}));
}
