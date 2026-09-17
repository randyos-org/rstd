const std = @import("std");
const zigstd = @import("../../../../zigstd.zig");
const Allocator = zigstd.Allocator;

const Sfnt = @import("../Sfnt.zig");
const types = @import("../stream_types.zig");
const GlyphIndex = types.GlyphIndex;
const pack = @import("pack.zig");

// A persistent, incrementally-populated glyph atlas: unlike
// render/atlas.zig's buildTextAtlas (which rasterizes and packs every glyph
// from scratch on every call), a GlyphCache rasterizes a glyph once and
// reuses its placement for as long as the cache lives. One cache is meant
// to back one (font, pixel size) pairing, since a glyph's rasterized size
// depends on scale - mixing sizes through the same cache would require a
// compound key this doesn't bother with.
const Self = @This();

/// Row-major 8bpp coverage, `width` bytes per row, same convention as
/// Sfnt.glyphBitmap. Grows (see grow()) but is never shrunk or relocated
/// mid-cell, so every previously-returned PixelRect stays valid for the
/// cache's whole lifetime.
pixels: []u8,
width: u32,
height: u32,
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
    const pixels = try gpa.alloc(u8, @as(usize, initial_size) * initial_size);
    @memset(pixels, 0);
    errdefer gpa.free(pixels);

    const pack_ctx = try pack.Context.create(gpa, @intCast(initial_size), @intCast(initial_size), .{});

    return .{
        .pixels = pixels,
        .width = initial_size,
        .height = initial_size,
        .pack_ctx = pack_ctx,
    };
}

pub fn deinit(self: *Self, gpa: Allocator) void {
    gpa.free(self.pixels);
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
    @memset(self.pixels, 0);
}

/// Returns this glyph's placement within `pixels`/`width`/`height`, lazily
/// rasterizing and packing it on first request; later requests for the same
/// glyph are a hashmap lookup. Returns null for a glyph with no visible
/// outline (space, ...), matching Sfnt.glyphBitmap's "empty" convention.
pub fn getOrAdd(self: *Self, gpa: Allocator, sfnt: Sfnt, glyph: GlyphIndex, scale: f32) !?PixelRect {
    if (self.locations.get(glyph)) |rect| {
        return if (rect.w == 0) null else rect;
    }

    var raster_pixels: std.ArrayListUnmanaged(u8) = .empty;
    defer raster_pixels.deinit(gpa);

    const bmp = try sfnt.glyphBitmap(gpa, &raster_pixels, glyph, scale, scale);
    if (bmp.width == 0 or bmp.height == 0) {
        try self.locations.put(gpa, glyph, empty_rect);
        return null;
    }

    const rect = try self.place(gpa, bmp.width, bmp.height);
    self.blit(rect, raster_pixels.items, bmp.width);
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
/// is pushed onto the packer's free-space list - not a full re-pack.
fn grow(self: *Self, gpa: Allocator) !void {
    const grow_width = self.width <= self.height;
    const new_width: u32 = if (grow_width) self.width * 2 else self.width;
    const new_height: u32 = if (grow_width) self.height else self.height * 2;

    const new_pixels = try gpa.alloc(u8, @as(usize, new_width) * new_height);
    @memset(new_pixels, 0);
    for (0..self.height) |row| {
        const src = self.pixels[row * self.width ..][0..self.width];
        @memcpy(new_pixels[row * new_width ..][0..self.width], src);
    }
    gpa.free(self.pixels);
    self.pixels = new_pixels;

    const new_area: pack.Rect = if (grow_width)
        .{ .x = @intCast(self.width), .y = 0, .w = @intCast(self.width), .h = @intCast(self.height) }
    else
        .{ .x = 0, .y = @intCast(self.height), .w = @intCast(new_width), .h = @intCast(self.height) };
    try self.pack_ctx.list.append(self.pack_ctx.allocator, new_area);

    self.width = new_width;
    self.height = new_height;
}

fn blit(self: *Self, rect: PixelRect, src: []const u8, src_stride: u16) void {
    for (0..rect.h) |row| {
        const src_row = src[row * src_stride ..][0..rect.w];
        const dst_off = (rect.y + row) * self.width + rect.x;
        @memcpy(self.pixels[dst_off..][0..rect.w], src_row);
    }
}

// ---
// Same synthetic-font approach as layout.zig's tests: units_per_em = 1000,
// glyph 1 is a visible 200x300-funit box mapped from every uppercase ASCII
// letter, glyph 2 is an invisible (no outline) space.

const bytes = @import("../../../../bytes.zig");
const general = @import("../tables/general.zig");
const TableTag = general.TableTag;

fn appendTable(buf: []u8, table_dir_off: usize, table_idx: usize, tag: TableTag, data: []const u8, off: usize) void {
    const rec_off = table_dir_off + table_idx * 16;
    bytes.setU32BE(buf[rec_off..], tag.toInt());
    bytes.setU32BE(buf[rec_off + 4 ..], 0);
    bytes.setU32BE(buf[rec_off + 8 ..], @intCast(off));
    bytes.setU32BE(buf[rec_off + 12 ..], @intCast(data.len));
    @memcpy(buf[off..][0..data.len], data);
}

fn buildTestFont(alloc: Allocator) ![]u8 {
    const num_tables = 7;
    const table_dir_off = 12;
    const tables_start = table_dir_off + num_tables * 16;

    var head_buf: [54]u8 = @splat(0);
    bytes.setU16BE(head_buf[18..], 1000);
    bytes.setI16BE(head_buf[50..], 1);

    var hhea_buf: [36]u8 = @splat(0);
    bytes.setI16BE(hhea_buf[4..], 800);
    bytes.setI16BE(hhea_buf[6..], -200);
    bytes.setI16BE(hhea_buf[8..], 0);
    bytes.setU16BE(hhea_buf[34..], 3);

    var maxp_buf: [6]u8 = @splat(0);
    bytes.setU32BE(maxp_buf[0..], 0x00005000);
    bytes.setU16BE(maxp_buf[4..], 3);

    var hmtx_buf: [12]u8 = @splat(0);
    bytes.setU16BE(hmtx_buf[4..], 250);
    bytes.setU16BE(hmtx_buf[8..], 150);

    var loca_buf: [16]u8 = @splat(0);
    bytes.setU32BE(loca_buf[0..], 0);
    bytes.setU32BE(loca_buf[4..], 0);
    bytes.setU32BE(loca_buf[8..], 34);
    bytes.setU32BE(loca_buf[12..], 34);

    var glyf_buf: [34]u8 = @splat(0);
    bytes.setI16BE(glyf_buf[0..], 1);
    bytes.setI16BE(glyf_buf[2..], 0);
    bytes.setI16BE(glyf_buf[4..], 0);
    bytes.setI16BE(glyf_buf[6..], 200);
    bytes.setI16BE(glyf_buf[8..], 300);
    bytes.setU16BE(glyf_buf[10..], 3);
    bytes.setU16BE(glyf_buf[12..], 0);
    const flags_start = 14;
    @memset(glyf_buf[flags_start..][0..4], 0x01);
    const x_start = flags_start + 4;
    bytes.setI16BE(glyf_buf[x_start + 0 ..], 0);
    bytes.setI16BE(glyf_buf[x_start + 2 ..], 200);
    bytes.setI16BE(glyf_buf[x_start + 4 ..], 0);
    bytes.setI16BE(glyf_buf[x_start + 6 ..], -200);
    const y_start = x_start + 8;
    bytes.setI16BE(glyf_buf[y_start + 0 ..], 0);
    bytes.setI16BE(glyf_buf[y_start + 2 ..], 0);
    bytes.setI16BE(glyf_buf[y_start + 4 ..], 300);
    bytes.setI16BE(glyf_buf[y_start + 6 ..], 0);

    var cmap_buf: [4 + 8 + 262]u8 = @splat(0);
    bytes.setU16BE(cmap_buf[0..], 0);
    bytes.setU16BE(cmap_buf[2..], 1);
    bytes.setU16BE(cmap_buf[4..], 3);
    bytes.setU16BE(cmap_buf[6..], 1);
    bytes.setU32BE(cmap_buf[8..], 12);
    const sub = cmap_buf[12..];
    bytes.setU16BE(sub, 0);
    bytes.setU16BE(sub[2..], 262);
    bytes.setU16BE(sub[4..], 0);
    var c: u8 = 'A';
    while (c <= 'Z') : (c += 1) sub[6 + @as(usize, c)] = 1;
    sub[6 + @as(usize, ' ')] = 2;

    var buf = try alloc.alloc(u8, tables_start + head_buf.len + hhea_buf.len + maxp_buf.len + hmtx_buf.len + loca_buf.len + cmap_buf.len + glyf_buf.len);
    errdefer alloc.free(buf);
    @memset(buf, 0);

    bytes.setU32BE(buf[0..], 0x00010000);
    bytes.setU16BE(buf[4..], num_tables);

    var off: usize = tables_start;
    appendTable(buf, table_dir_off, 0, .head, &head_buf, off);
    off += head_buf.len;
    appendTable(buf, table_dir_off, 1, .hhea, &hhea_buf, off);
    off += hhea_buf.len;
    appendTable(buf, table_dir_off, 2, .maxp, &maxp_buf, off);
    off += maxp_buf.len;
    appendTable(buf, table_dir_off, 3, .hmtx, &hmtx_buf, off);
    off += hmtx_buf.len;
    appendTable(buf, table_dir_off, 4, .loca, &loca_buf, off);
    off += loca_buf.len;
    appendTable(buf, table_dir_off, 5, .cmap, &cmap_buf, off);
    off += cmap_buf.len;
    appendTable(buf, table_dir_off, 6, .glyf, &glyf_buf, off);
    off += glyf_buf.len;

    return buf;
}

test "GlyphCache.getOrAdd caches a glyph's placement across repeated calls" {
    const alloc = zigstd.testAllocator;
    const font_bytes = try buildTestFont(alloc);
    defer alloc.free(font_bytes);
    const sfnt = try Sfnt.init(alloc, font_bytes, .{ .check = true });

    var cache = try Self.init(alloc, 256);
    defer cache.deinit(alloc);

    const glyph: GlyphIndex = @enumFromInt(1);
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
    const font_bytes = try buildTestFont(alloc);
    defer alloc.free(font_bytes);
    const sfnt = try Sfnt.init(alloc, font_bytes, .{ .check = true });

    var cache = try Self.init(alloc, 256);
    defer cache.deinit(alloc);

    const space: GlyphIndex = @enumFromInt(2);
    try zigstd.testExpectEqual(@as(?Self.PixelRect, null), try cache.getOrAdd(alloc, sfnt, space, 1.0));
    try zigstd.testExpectEqual(@as(usize, 1), cache.locations.count());
    try zigstd.testExpectEqual(@as(?Self.PixelRect, null), try cache.getOrAdd(alloc, sfnt, space, 1.0));
}

test "GlyphCache grows the canvas to fit a glyph larger than the initial size, preserving its pixels" {
    const alloc = zigstd.testAllocator;
    const font_bytes = try buildTestFont(alloc);
    defer alloc.free(font_bytes);
    const sfnt = try Sfnt.init(alloc, font_bytes, .{ .check = true });

    // A 200x300px bitmap (scale 1.0) can't possibly fit an 8x8 canvas, so
    // this forces several grow() doublings before the first placement.
    var cache = try Self.init(alloc, 8);
    defer cache.deinit(alloc);

    const glyph: GlyphIndex = @enumFromInt(1);
    const rect = (try cache.getOrAdd(alloc, sfnt, glyph, 1.0)).?;
    try zigstd.testExpectEqual(@as(u32, 200), rect.w);
    try zigstd.testExpectEqual(@as(u32, 300), rect.h);
    try zigstd.testExpectTrue(cache.width >= rect.x + rect.w);
    try zigstd.testExpectTrue(cache.height >= rect.y + rect.h);

    // A fully-covered interior pixel of the rasterized box should have
    // landed intact at its placement in the (grown) atlas.
    const mid_x = rect.x + rect.w / 2;
    const mid_y = rect.y + rect.h / 2;
    try zigstd.testExpectTrue(cache.pixels[mid_y * cache.width + mid_x] > 0);
}
