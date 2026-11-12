const std = @import("std");
const zigstd = @import("../../../zigstd.zig");
const Allocator = zigstd.Allocator;

const Sfnt = @import("../sfnt/Sfnt.zig");
const types = @import("../sfnt/types.zig");
const pack = @import("../../pack.zig");
const layout = @import("layout.zig");
const GlyphCache = @import("GlyphCache.zig");
const rasterize = @import("rasterize.zig");
const msdf = @import("../msdf/Generator.zig");
const Image = @import("../../image.zig").Image;

const GlyphIndex = types.GlyphIndex;
const test_font = @import("_test/font.zig");

/// One glyph's placement within a TextAtlas: where to draw it (in the same
/// coordinate space as buildTextAtlas's start_x/baseline_y, already
/// including this glyph's own bearing and pen advance) and where its
/// bitmap sits within TextAtlas.pixels.
pub const GlyphPlacement = struct {
    pos: [2]f32,
    size: [2]f32,
    atlas_pos: [2]u32,
};

/// A run of text rasterized into one shared atlas texture, laid out left to
/// right with hmtx advances and kern-table/GPOS kerning between consecutive
/// pairs. Caller owns both slices.
pub const TextAtlas = struct {
    /// Row-major, `width * channels` bytes per row. For `.coverage`, one
    /// 8-bit coverage byte per pixel (0 transparent, 255 opaque). For
    /// `.msdf`, `Options.msdf.sdf_type.bytesPerPixel()` signed-distance bytes
    /// per pixel (msdf10: one packed `msdf.Msdf10Pixel`).
    pixels: []u8,
    width: u32,
    height: u32,
    channels: u8,
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
    /// Multi-channel signed distance field, via the msdf module's Generator.
    msdf,
};

/// Only consulted when `Options.type == .msdf`. Mirrors the handful of
/// msdf.Options fields relevant to a plain atlas build; the rest
/// (var_font_args, disable_concurrency, winding, ...) keep Generator's own
/// defaults since they don't apply here.
pub const MsdfOptions = struct {
    sdf_type: msdf.SdfType = .msdf,
    /// Width, in pixels, of the distance-field falloff around each edge.
    px_range: u16 = 8,
    coloring_method: msdf.ColoringMethod = .distance,
    corner_angle_threshold: f64 = 3.0,
    coloring_rng_seed: u64 = 0,
    error_correction_opts: ?msdf.error_correction.Options = .{},
    /// Fails the whole atlas build with error.InvalidShape instead of
    /// silently rendering a bad outline if a glyph's decoded contours don't
    /// form a closed loop.
    validate_shape: bool = true,
    /// Deconverges near-degenerate/converging curve corners before
    /// generating the SDF; fixes some corner-adjacent stray-texel
    /// artifacts that show up at small sizes.
    normalize_shape: bool = true,
    /// Geometrically determines each contour's correct winding direction,
    /// rather than relying on the single-sample-point `winding: .guess`
    /// heuristic; matters for multi-contour glyphs (most letters with a
    /// counter/hole).
    orient_contours: bool = true,
};

pub const Options = struct {
    type: AtlasType = .coverage,
    msdf: MsdfOptions = .{},
};

/// Rasterizes every glyph in `text`; via glyphBitmapSubpixel (the analytic
/// coverage rasterizer) for `.coverage`, or via the msdf module's Generator
/// for `.msdf`; packs the bitmaps into one shared atlas (via the guillotine
/// packer in pack.zig, growing the canvas until everything fits), and lays
/// them out left to right from (start_x, baseline_y).
pub fn buildTextAtlas(
    self: Sfnt,
    gpa: Allocator,
    text: []const u8,
    pixel_size: f32,
    start_x: f32,
    baseline_y: f32,
    options: Options,
) !TextAtlas {
    const channels: u8 = switch (options.type) {
        .coverage => 1,
        .msdf => options.msdf.sdf_type.bytesPerPixel(),
    };

    const head = self.head orelse return error.MissingHead;
    const scale = pixel_size / @as(f32, @floatFromInt(head.unitsPerEm()));

    var generator: msdf = .init(self);
    const msdf_opts: msdf.Options = .{
        .sdf_type = options.msdf.sdf_type,
        .px_size = @round(pixel_size),
        .px_range = options.msdf.px_range,
        .coloring_rng_seed = options.msdf.coloring_rng_seed,
        .coloring_method = options.msdf.coloring_method,
        .corner_angle_threshold = options.msdf.corner_angle_threshold,
        .error_correction_opts = options.msdf.error_correction_opts,
        .validate_shape = options.msdf.validate_shape,
        .normalize_shape = options.msdf.normalize_shape,
        .orient_contours = options.msdf.orient_contours,
    };
    // Half the falloff range, in pixels - the padding the msdf module's
    // Generator adds around a glyph's raw outline bounds on every side (see
    // its generateSingle:
    // tfm = {bounds.left - px_range/2, bounds.bottom - px_range/2}, w/h
    // include a full px_range). Needed below to recover where a glyph's
    // pixel buffer's top-left corner sits relative to its own outline
    // bounding box, which the MsdfOptions/GlyphMetricsEm bearing fields
    // alone don't capture.
    const half_px_range: f32 = @as(f32, @floatFromInt(options.msdf.px_range)) / 2.0;

    var glyphs: std.ArrayList(GlyphIndex) = .empty;
    defer glyphs.deinit(gpa);
    var rastered: std.ArrayList(RasterGlyph) = .empty;
    defer rastered.deinit(gpa);

    var raster_pixels: std.ArrayList(u8) = .empty;
    defer raster_pixels.deinit(gpa);

    const view = std.unicode.Utf8View.init(text) catch return error.InvalidUtf8;
    var code_points = view.iterator();
    while (code_points.nextCodepoint()) |ch| {
        const glyph = self.getGlyphIndex(ch) orelse continue;

        switch (options.type) {
            .coverage => {
                const advance = @as(f32, @floatFromInt(self.getGlyphMetrics(glyph).advance_width)) * scale;

                const off = raster_pixels.items.len;
                const bmp = try rasterize.glyphBitmapSubpixel(self, gpa, &raster_pixels, glyph, scale, scale, 0, 0);

                try glyphs.append(gpa, glyph);
                try rastered.append(gpa, .{
                    .bitmap_off = if (bmp.width == 0 or bmp.height == 0) null else off,
                    .width = bmp.width,
                    .height = bmp.height,
                    .bearing_x = @floatFromInt(bmp.off_x),
                    .bearing_y = @floatFromInt(bmp.off_y),
                    .advance = advance,
                });
            },
            .msdf => {
                // generateSingle resolves the glyph index itself from the
                // codepoint, and returns everything (including .width/
                // .height, already in pixels) em-scaled rather than
                // pixel-scaled - multiply by pixel_size to get pixels.
                var generated = try generator.generateSingle(gpa, ch, &msdf_opts);
                defer generated.deinit(gpa);

                const w = generated.metrics.width;
                const h = generated.metrics.height;
                const off = raster_pixels.items.len;
                if (w != 0 and h != 0) try raster_pixels.appendSlice(gpa, generated.pixels);

                const bearing_x_px: f32 = @as(f32, @floatCast(generated.metrics.bearing_x)) * pixel_size;

                // Anchor from y_min (the bottom, unaffected by the truncation),
                // then subtract the actual (truncated) height to get the top.
                // Returns null for a glyph with no outline at all (e.g. space);
                // generateSingle already reflects that above via w==0/h==0,
                // so bottom_offset is simply unused in that case.
                const bottom_offset: f32 = if (self.getGlyphBoundingBox(glyph)) |bbox| blk: {
                    const y_min_px: f32 = @as(f32, @floatFromInt(bbox.y_min)) / @as(f32, @floatFromInt(head.unitsPerEm())) * pixel_size;
                    break :blk half_px_range - y_min_px;
                } else 0;

                try glyphs.append(gpa, glyph);
                try rastered.append(gpa, .{
                    .bitmap_off = if (w == 0 or h == 0) null else off,
                    .width = w,
                    .height = h,
                    // Top-left corner of the padded SDF box, relative to the
                    // glyph origin (screen convention: x right, y down).
                    .bearing_x = bearing_x_px - half_px_range,
                    .bearing_y = bottom_offset - @as(f32, @floatFromInt(h)),
                    .advance = @as(f32, @floatCast(generated.metrics.advance)) * pixel_size,
                });
            },
        }
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
    var canvas: u32 = 256;
    var pack_ctx: pack.Context = while (true) {
        var ctx = try pack.Context.init(gpa, canvas, canvas, .{});
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
        atlas_w = @max(atlas_w, @as(u32, @intCast(ir.rect.bottomRight().x)));
        atlas_h = @max(atlas_h, @as(u32, @intCast(ir.rect.bottomRight().y)));
    }

    const atlas_pixels: []u8 = switch (channels) {
        inline 1, 3, 4 => |n| blk: {
            const Px = [n]u8;
            const atlas: Image(Px) = try .init(gpa, atlas_w, atlas_h, @splat(0));
            for (rastered.items, placements) |g, p| {
                const off = g.bitmap_off orelse continue;
                const src = raster_pixels.items[off..][0 .. @as(usize, g.width) * g.height * n];
                atlas.blit(p.x, p.y, g.width, std.mem.bytesAsSlice(Px, src));
            }
            break :blk std.mem.sliceAsBytes(atlas.data);
        },
        else => unreachable,
    };
    errdefer gpa.free(atlas_pixels);

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
        .channels = channels,
        .glyphs = try out_glyphs.toOwnedSlice(gpa),
    };
}

/// One glyph's placement when drawn against a persistent `GlyphCache`:
/// unlike GlyphPlacement, `atlas_pos` here refers into the *cache's* texture
/// (`cache.atlas`), which the caller already holds and which outlives any
/// single buildTextAtlasCached call.
pub const CachedGlyphPlacement = struct {
    pos: [2]f32,
    size: [2]f32,
    atlas_pos: [2]u32,
};

pub const CachedTextLayout = struct {
    glyphs: []CachedGlyphPlacement,
    width: u32,
    height: u32,

    pub fn deinit(self: *CachedTextLayout, gpa: Allocator) void {
        gpa.free(self.glyphs);
        self.* = undefined;
    }
};

/// Like buildTextAtlas, but through a persistent `GlyphCache`: glyphs it has
/// seen are reused, not re-rasterized; `layout.Options` adds word wrap.
/// buildTextAtlas stays separate: its bitmap-box rounding may differ from this
/// layout-box rounding, and it has no tests to check a shared path against.
///
/// `cache` must be for this same font and pixel_size (see GlyphCache).
pub fn buildTextAtlasCached(
    self: Sfnt,
    gpa: Allocator,
    cache: *GlyphCache,
    text: []const u8,
    pixel_size: f32,
    start_x: f32,
    baseline_y: f32,
    layout_options: layout.Options,
) !CachedTextLayout {
    const head = self.head orelse return error.MissingHead;
    const scale = pixel_size / @as(f32, @floatFromInt(head.unitsPerEm()));

    var text_layout = try layout.layoutText(gpa, self, text, pixel_size, layout_options);
    defer text_layout.deinit(gpa);

    var out_glyphs: std.ArrayList(CachedGlyphPlacement) = .empty;
    errdefer out_glyphs.deinit(gpa);

    for (text_layout.glyphs) |g| {
        const rect = try cache.getOrAdd(gpa, self, g.glyph, scale) orelse continue;
        try out_glyphs.append(gpa, .{
            // pixel_y2 is the glyph's topmost point above the baseline (Y-up,
            // from layout's funit bounding box); screen/texture space is
            // Y-down from baseline_y, so the bitmap's top edge lands at
            // baseline_y - pixel_y2.
            .pos = .{
                start_x + @as(f32, @floatFromInt(g.pixel_x1)),
                baseline_y - @as(f32, @floatFromInt(g.pixel_y2)),
            },
            .size = .{ @floatFromInt(rect.w), @floatFromInt(rect.h) },
            .atlas_pos = .{ rect.x, rect.y },
        });
    }

    return .{
        .glyphs = try out_glyphs.toOwnedSlice(gpa),
        .width = text_layout.width(),
        .height = text_layout.height(),
    };
}

test "buildTextAtlas with .msdf produces a correctly-shaped multi-channel atlas" {
    const alloc = zigstd.testAllocator;
    const font_bytes = try test_font.build(alloc);
    defer alloc.free(font_bytes);
    const sfnt = try Sfnt.init(alloc, font_bytes, .{ .check = true });

    var atlas = try buildTextAtlas(sfnt, alloc, "A", 32, 0, 0, .{ .type = .msdf });
    defer atlas.deinit(alloc);

    try zigstd.testExpectEqual(@as(u8, 3), atlas.channels); // default sdf_type is .msdf
    try zigstd.testExpectEqual(@as(usize, 1), atlas.glyphs.len);
    try zigstd.testExpectEqual(
        @as(usize, atlas.width) * atlas.height * atlas.channels,
        atlas.pixels.len,
    );
    try zigstd.testExpectTrue(atlas.glyphs[0].size[0] > 0);
    try zigstd.testExpectTrue(atlas.glyphs[0].size[1] > 0);
}

test "buildTextAtlas decodes UTF-8 and rejects invalid UTF-8" {
    const alloc = zigstd.testAllocator;
    const font_bytes = try test_font.build(alloc);
    defer alloc.free(font_bytes);
    const sfnt = try Sfnt.init(alloc, font_bytes, .{ .check = true });

    // "AÉ": byte by byte, 0xC3 and 0x89 map to no glyph
    var atlas = try buildTextAtlas(sfnt, alloc, "A\u{C9}", 32, 0, 0, .{});
    defer atlas.deinit(alloc);
    try zigstd.testExpectEqual(@as(usize, 2), atlas.glyphs.len);
    try zigstd.testExpectError(error.InvalidUtf8, buildTextAtlas(sfnt, alloc, "A\xC3", 32, 0, 0, .{}));
}

test "buildTextAtlas .msdf10 packs the same distances as .msdf" {
    const alloc = zigstd.testAllocator;
    const font_bytes = try test_font.build(alloc);
    defer alloc.free(font_bytes);
    const sfnt = try Sfnt.init(alloc, font_bytes, .{ .check = true });

    var m8 = try buildTextAtlas(sfnt, alloc, "A", 32, 0, 0, .{ .type = .msdf });
    defer m8.deinit(alloc);
    var m10 = try buildTextAtlas(sfnt, alloc, "A", 32, 0, 0, .{ .type = .msdf, .msdf = .{ .sdf_type = .msdf10 } });
    defer m10.deinit(alloc);

    try zigstd.testExpectEqual(@as(u8, 4), m10.channels);
    try zigstd.testExpectEqual(m8.width * m8.height, m10.width * m10.height);
    for (0..m8.width * m8.height) |i| {
        const p = std.mem.bytesToValue(msdf.Msdf10Pixel, m10.pixels[i * 4 ..][0..4]);
        for ([_]u10{ p.r, p.g, p.b }, m8.pixels[i * 3 ..][0..3]) |c10, c8| {
            const scaled = @as(f64, @floatFromInt(c10)) * 255.0 / 1023.0;
            try zigstd.testExpectTrue(@abs(scaled - @as(f64, @floatFromInt(c8))) <= 1.0);
        }
    }
}

test "buildTextAtlasCached reuses a glyph's atlas placement across separate calls" {
    const alloc = zigstd.testAllocator;
    const font_bytes = try test_font.build(alloc);
    defer alloc.free(font_bytes);
    const sfnt = try Sfnt.init(alloc, font_bytes, .{ .check = true });

    var cache = try GlyphCache.init(alloc, 256);
    defer cache.deinit(alloc);

    var first = try buildTextAtlasCached(sfnt, alloc, &cache, "AA", 1000, 0, 0, .{});
    defer first.deinit(alloc);
    try zigstd.testExpectEqual(@as(usize, 2), first.glyphs.len);
    try zigstd.testExpectEqual(first.glyphs[0].atlas_pos, first.glyphs[1].atlas_pos);

    var second = try buildTextAtlasCached(sfnt, alloc, &cache, "A", 1000, 0, 0, .{});
    defer second.deinit(alloc);
    try zigstd.testExpectEqual(first.glyphs[0].atlas_pos, second.glyphs[0].atlas_pos);
}

test "buildTextAtlasCached wraps text across lines when a wrap width is given" {
    const alloc = zigstd.testAllocator;
    const font_bytes = try test_font.build(alloc);
    defer alloc.free(font_bytes);
    const sfnt = try Sfnt.init(alloc, font_bytes, .{ .check = true });

    var cache = try GlyphCache.init(alloc, 256);
    defer cache.deinit(alloc);

    var result = try buildTextAtlasCached(sfnt, alloc, &cache, "AB CD", 1000, 0, 0, .{ .wrap_width_px = 500 });
    defer result.deinit(alloc);

    try zigstd.testExpectEqual(@as(usize, 4), result.glyphs.len);
    try zigstd.testExpectTrue(result.glyphs[2].pos[1] > result.glyphs[0].pos[1]);
}
