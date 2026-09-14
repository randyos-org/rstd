const std = @import("std");
const zigstd = @import("../../../zigstd.zig");
const Allocator = zigstd.Allocator;

const Sfnt = @import("../Sfnt.zig");
const types = @import("../stream_types.zig");
const pack = @import("pack.zig");

const GlyphIndex = types.GlyphIndex;

/// One glyph's placement within a TextAtlas: where to draw it (in the same
/// coordinate space as buildTextAtlas's start_x/baseline_y, already
/// including this glyph's own bearing and pen advance) and where its
/// bitmap sits within TextAtlas.pixels.
pub const GlyphPlacement = struct {
    pos: [2]f32,
    size: [2]f32,
    atlas_pos: [2]u32,
};

/// A run of text rasterized into one shared 8-bit coverage atlas, laid out
/// left to right with hmtx advances and kern-table/GPOS kerning between
/// consecutive pairs. Caller owns both slices.
pub const TextAtlas = struct {
    /// Row-major 8bpp coverage, `width` bytes per row. 0 transparent, 255 opaque.
    pixels: []u8,
    width: u32,
    height: u32,
    /// One entry per glyph that actually produced a visible bitmap; glyphs
    /// with no outline (space, ...) still advance the pen but are omitted
    /// here, matching GlyphBitmap's own "empty" convention.
    glyphs: []GlyphPlacement,

    pub fn deinit(self: *TextAtlas, gpa: Allocator) void {
        gpa.free(self.pixels);
        gpa.free(self.glyphs);
        self.* = undefined;
    }
};

const RasterGlyph = struct {
    /// Offset into the scratch raster buffer, or null for a glyph with no
    /// visible bitmap (width/height 0).
    bitmap_off: ?usize,
    width: u32,
    height: u32,
    bearing_x: f32,
    bearing_y: f32,
    advance: f32,
};

const Placement = struct { x: u32, y: u32 };

fn packSortLessThan(_: void, a: pack.IdRect, b: pack.IdRect) bool {
    return @max(a.rect.w, a.rect.h) > @max(b.rect.w, b.rect.h);
}

pub const AtlasType = enum {
    /// Plain 8-bit per-pixel coverage, same convention as glyphBitmap.
    coverage,
    /// Multi-channel signed distance field. Not implemented yet.
    msdf,
};

pub const Options = struct {
    type: AtlasType = .coverage,
};

/// Rasterizes every glyph in `text` (via glyphBitmapSubpixel - the real
/// analytic-coverage scanline rasterizer), packs the bitmaps into one
/// shared atlas (via the guillotine packer in pack.zig, growing the canvas
/// until everything fits), and lays them out left to right from
/// (start_x, baseline_y).
pub fn buildTextAtlas(
    self: Sfnt,
    gpa: Allocator,
    text: []const u8,
    pixel_size: f32,
    start_x: f32,
    baseline_y: f32,
    options: Options,
) !TextAtlas {
    switch (options.type) {
        .coverage => {},
        .msdf => @panic("TODO: MSDF atlas generation is not implemented yet"),
    }

    const head = self.head orelse return error.MissingHead;
    const scale = pixel_size / @as(f32, @floatFromInt(head.unitsPerEm()));

    var glyphs: std.ArrayList(GlyphIndex) = .empty;
    defer glyphs.deinit(gpa);
    var rastered: std.ArrayList(RasterGlyph) = .empty;
    defer rastered.deinit(gpa);

    var raster_pixels: std.ArrayListUnmanaged(u8) = .empty;
    defer raster_pixels.deinit(gpa);

    for (text) |ch| {
        const glyph = self.getGlyphIndex(ch) orelse continue;
        const advance = @as(f32, @floatFromInt(self.getGlyphMetrics(glyph).advance_width)) * scale;

        const off = raster_pixels.items.len;
        const bmp = try self.glyphBitmapSubpixel(gpa, &raster_pixels, glyph, scale, scale, 0, 0);

        try glyphs.append(gpa, glyph);
        try rastered.append(gpa, .{
            .bitmap_off = if (bmp.width == 0 or bmp.height == 0) null else off,
            .width = bmp.width,
            .height = bmp.height,
            .bearing_x = @floatFromInt(bmp.off_x),
            .bearing_y = @floatFromInt(bmp.off_y),
            .advance = advance,
        });
    }

    // Pack every rastered glyph's bitmap into one shared atlas via a
    // guillotine rect packer (pack.zig - general-purpose, not font-specific).
    var id_rects: std.ArrayList(pack.IdRect) = .empty;
    defer id_rects.deinit(gpa);
    for (rastered.items, 0..) |g, i| {
        if (g.bitmap_off == null) continue;
        try id_rects.append(gpa, .{ .id = @intCast(i), .rect = .{ .w = @intCast(g.width), .h = @intCast(g.height) } });
    }

    // Start from a modest canvas and double it until every glyph fits,
    // rather than assuming a fixed size up front.
    var canvas: i32 = 256;
    var pack_ctx: pack.Context = while (true) {
        var ctx = try pack.Context.create(gpa, canvas, canvas, .{});
        const result = pack.pack(pack.IdRect, &ctx, id_rects.items, .{ .sortLessThanFn = packSortLessThan });
        if (result) {
            break ctx;
        } else |err| switch (err) {
            error.NoSpaceLeft => {
                ctx.deinit();
                canvas *= 2;
            },
            else => |e| return e,
        }
    };
    defer pack_ctx.deinit();

    var atlas_w: u32 = 1;
    var atlas_h: u32 = 1;
    const placements = try gpa.alloc(Placement, rastered.items.len);
    defer gpa.free(placements);
    @memset(placements, .{ .x = 0, .y = 0 });
    for (id_rects.items) |ir| {
        const idx: usize = @intCast(ir.id);
        placements[idx] = .{ .x = @intCast(ir.rect.x), .y = @intCast(ir.rect.y) };
        atlas_w = @max(atlas_w, @as(u32, @intCast(ir.rect.x + ir.rect.w)));
        atlas_h = @max(atlas_h, @as(u32, @intCast(ir.rect.y + ir.rect.h)));
    }

    const atlas_pixels = try gpa.alloc(u8, @as(usize, atlas_w) * atlas_h);
    errdefer gpa.free(atlas_pixels);
    @memset(atlas_pixels, 0);
    for (rastered.items, placements) |g, p| {
        const off = g.bitmap_off orelse continue;
        const src = raster_pixels.items[off..];
        for (0..g.height) |row| {
            const src_row = src[row * g.width ..][0..g.width];
            const dst_off = (p.y + row) * atlas_w + p.x;
            @memcpy(atlas_pixels[dst_off..][0..g.width], src_row);
        }
    }

    // Second pass: walk the glyphs again in pen order, applying kerning
    // between consecutive pairs, to build the final placement list.
    var out_glyphs: std.ArrayList(GlyphPlacement) = .empty;
    errdefer out_glyphs.deinit(gpa);

    var pen_x: f32 = start_x;
    var prev_glyph: ?GlyphIndex = null;
    for (rastered.items, glyphs.items, placements) |g, gi, p| {
        if (prev_glyph) |pg| pen_x += @as(f32, @floatFromInt(self.glyphKernAdvance(pg, gi))) * scale;
        prev_glyph = gi;

        if (g.bitmap_off != null) {
            try out_glyphs.append(gpa, .{
                .pos = .{ pen_x + g.bearing_x, baseline_y + g.bearing_y },
                .size = .{ @floatFromInt(g.width), @floatFromInt(g.height) },
                .atlas_pos = .{ p.x, p.y },
            });
        }

        pen_x += g.advance;
    }

    return .{
        .pixels = atlas_pixels,
        .width = atlas_w,
        .height = atlas_h,
        .glyphs = try out_glyphs.toOwnedSlice(gpa),
    };
}
