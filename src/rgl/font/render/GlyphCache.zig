const std = @import("std");
const zigstd = @import("../../../zigstd.zig");
const Allocator = zigstd.Allocator;

const Sfnt = @import("../sfnt/Sfnt.zig");
const types = @import("../sfnt/types.zig");
const GlyphIndex = types.GlyphIndex;
const pack = @import("../../pack.zig");
const Image = @import("../../image.zig").Image;
const test_font = @import("_test/font.zig");

// A persistent, incrementally-populated glyph atlas: unlike
// render/atlas.zig's buildTextAtlas (which rasterizes and packs every glyph
// from scratch on every call), a GlyphCache rasterizes a glyph once and
// reuses its placement for as long as the cache lives. One cache is meant
// to back one (font, pixel size) pairing, since a glyph's rasterized size
// depends on scale; mixing sizes through the same cache would require a
// compound key this doesn't bother with.
const Self = @This();

/// 8bpp coverage, same convention as Sfnt.glyphBitmap. Grows (see grow())
/// but is never shrunk or relocated mid-cell, so every previously-returned
/// PixelRect stays valid for the cache's whole lifetime.
atlas: Image(u8),
pack_ctx: pack.Context,
locations: std.AutoHashMapUnmanaged(GlyphIndex, PixelRect) = .empty,

pub const PixelRect = struct {
    x: u32,
    y: u32,
    w: u32,
    h: u32,
};

/// Sentinel cached for a glyph with no visible outline (space, ...), so
/// repeated lookups don't re-run rasterization just to learn that again.
const empty_rect: PixelRect = .{ .x = 0, .y = 0, .w = 0, .h = 0 };

pub fn init(gpa: Allocator, initial_size: u32) !Self {
    var atlas: Image(u8) = try .init(gpa, initial_size, initial_size, 0);
    errdefer atlas.deinit(gpa);

    const pack_ctx = try pack.Context.init(gpa, @intCast(initial_size), @intCast(initial_size), .{});

    return .{ .atlas = atlas, .pack_ctx = pack_ctx };
}

pub fn deinit(self: *Self, gpa: Allocator) void {
    self.atlas.deinit(gpa);
    self.pack_ctx.deinit();
    self.locations.deinit(gpa);
    self.* = undefined;
}

/// Drops every cached placement and re-clears the backing texture, keeping
/// the current canvas size and allocation. Useful when reusing a cache for
/// a different font/size pairing rather than allocating a fresh one.
pub fn reset(self: *Self) !void {
    self.locations.clearRetainingCapacity();
    try self.pack_ctx.clear();
    @memset(self.atlas.data, 0);
}

/// Returns this glyph's placement within `atlas`, lazily
/// rasterizing and packing it on first request; later requests for the same
/// glyph are a hashmap lookup. Returns null for a glyph with no visible
/// outline (space, ...), matching Sfnt.glyphBitmap's "empty" convention.
pub fn getOrAdd(self: *Self, gpa: Allocator, sfnt: Sfnt, glyph: GlyphIndex, scale: f32) !?PixelRect {
    if (self.locations.get(glyph)) |rect| {
        return if (rect.w == 0) null else rect;
    }

    var raster_pixels: std.ArrayList(u8) = .empty;
    defer raster_pixels.deinit(gpa);

    const bmp = try sfnt.glyphBitmap(gpa, &raster_pixels, glyph, scale, scale);
    if (bmp.width == 0 or bmp.height == 0) {
        try self.locations.put(gpa, glyph, empty_rect);
        return null;
    }

    const rect = try self.place(gpa, bmp.width, bmp.height);
    self.atlas.blit(rect.x, rect.y, bmp.width, raster_pixels.items);
    try self.locations.put(gpa, glyph, rect);
    return rect;
}

fn place(self: *Self, gpa: Allocator, w: u16, h: u16) !PixelRect {
    while (true) {
        var rects = [_]pack.Rect{.{ .w = w, .h = h }};
        if (pack.pack(pack.Rect, &self.pack_ctx, &rects, .{})) {
            return .{ .x = @intCast(rects[0].x), .y = @intCast(rects[0].y), .w = w, .h = h };
        } else |err| switch (err) {
            error.NoSpaceLeft => try self.grow(gpa),
            else => |e| return e,
        }
    }
}

/// Doubles whichever dimension keeps the canvas closer to square. Existing
/// pixel content and every already-cached PixelRect stay valid: growth only
/// ever extends the canvas outward from the already-used top-left corner,
/// so nothing already placed ever needs to move. Only the newly-added strip
/// is pushed onto the packer's free-space list; not a full re-pack.
fn grow(self: *Self, gpa: Allocator) !void {
    const width = self.atlas.width;
    const height = self.atlas.height;
    const grow_width = width <= height;
    const new_width: u32 = if (grow_width) width * 2 else width;
    const new_height: u32 = if (grow_width) height else height * 2;

    var grown: Image(u8) = try .init(gpa, new_width, new_height, 0);
    grown.blit(0, 0, width, self.atlas.data);
    self.atlas.deinit(gpa);
    self.atlas = grown;

    const new_area: pack.Rect = if (grow_width)
        .{ .x = @intCast(width), .y = 0, .w = @intCast(width), .h = @intCast(height) }
    else
        .{ .x = 0, .y = @intCast(height), .w = @intCast(new_width), .h = @intCast(height) };
    try self.pack_ctx.list.append(self.pack_ctx.allocator, new_area);
}

test "GlyphCache.getOrAdd caches a glyph's placement across repeated calls" {
    const alloc = zigstd.testAllocator;
    const font_bytes = try test_font.build(alloc);
    defer alloc.free(font_bytes);
    const sfnt = try Sfnt.init(alloc, font_bytes, .{ .check = true });

    var cache = try Self.init(alloc, 256);
    defer cache.deinit(alloc);

    const glyph: GlyphIndex = @fromBackingInt(@intCast(1));
    const first = (try cache.getOrAdd(alloc, sfnt, glyph, 1.0)).?;
    try zigstd.testExpectEqual(@as(u32, 200), first.w);
    try zigstd.testExpectEqual(@as(u32, 300), first.h);

    try zigstd.testExpectEqual(@as(usize, 1), cache.locations.count());
    const second = (try cache.getOrAdd(alloc, sfnt, glyph, 1.0)).?;
    try zigstd.testExpectEqual(first.x, second.x);
    try zigstd.testExpectEqual(first.y, second.y);
    try zigstd.testExpectEqual(@as(usize, 1), cache.locations.count());
}

test "GlyphCache.getOrAdd returns null and caches the miss for a glyph with no visible outline" {
    const alloc = zigstd.testAllocator;
    const font_bytes = try test_font.build(alloc);
    defer alloc.free(font_bytes);
    const sfnt = try Sfnt.init(alloc, font_bytes, .{ .check = true });

    var cache = try Self.init(alloc, 256);
    defer cache.deinit(alloc);

    const space: GlyphIndex = @fromBackingInt(@intCast(2));
    try zigstd.testExpectEqual(@as(?Self.PixelRect, null), try cache.getOrAdd(alloc, sfnt, space, 1.0));
    try zigstd.testExpectEqual(@as(usize, 1), cache.locations.count());
    try zigstd.testExpectEqual(@as(?Self.PixelRect, null), try cache.getOrAdd(alloc, sfnt, space, 1.0));
}

test "GlyphCache grows the canvas to fit a glyph larger than the initial size, preserving its pixels" {
    const alloc = zigstd.testAllocator;
    const font_bytes = try test_font.build(alloc);
    defer alloc.free(font_bytes);
    const sfnt = try Sfnt.init(alloc, font_bytes, .{ .check = true });

    // A 200x300px bitmap (scale 1.0) can't possibly fit an 8x8 canvas, so
    // this forces several grow() doublings before the first placement.
    var cache = try Self.init(alloc, 8);
    defer cache.deinit(alloc);

    const glyph: GlyphIndex = @fromBackingInt(@intCast(1));
    const rect = (try cache.getOrAdd(alloc, sfnt, glyph, 1.0)).?;
    try zigstd.testExpectEqual(@as(u32, 200), rect.w);
    try zigstd.testExpectEqual(@as(u32, 300), rect.h);
    try zigstd.testExpectTrue(cache.atlas.width >= rect.x + rect.w);
    try zigstd.testExpectTrue(cache.atlas.height >= rect.y + rect.h);

    // A fully-covered interior pixel of the rasterized box should have
    // landed intact at its placement in the (grown) atlas.
    const mid_x = rect.x + rect.w / 2;
    const mid_y = rect.y + rect.h / 2;
    try zigstd.testExpectTrue(cache.atlas.getPixel(mid_x, mid_y) > 0);
}
