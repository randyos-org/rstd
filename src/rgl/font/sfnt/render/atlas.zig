const std = @import("std");
const zigstd = @import("../../../../zigstd.zig");
const Allocator = zigstd.Allocator;

const Sfnt = @import("../Sfnt.zig");
const types = @import("../stream_types.zig");
const pack = @import("pack.zig");
const layout = @import("layout.zig");
const GlyphCache = @import("GlyphCache.zig");
const msdf = @import("../msdf/__root__.zig");

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

/// A run of text rasterized into one shared atlas texture, laid out left to
/// right with hmtx advances and kern-table/GPOS kerning between consecutive
/// pairs. Caller owns both slices.
pub const TextAtlas = struct {
    /// Row-major, `width * channels` bytes per row. For `.coverage`, one
    /// 8-bit coverage byte per pixel (0 transparent, 255 opaque). For
    /// `.msdf`, `channels` signed-distance bytes per pixel per
    /// `Options.msdf.sdf_type.numChannels()`.
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
    error_correction_opts: ?msdf.ErrorCorrection.Options = .{},
    /// Fails the whole atlas build with error.InvalidShape instead of
    /// silently rendering a bad outline if a glyph's decoded contours don't
    /// form a closed loop.
    validate_shape: bool = true,
    /// Deconverges near-degenerate/converging curve corners before
    /// generating the SDF - fixes some corner-adjacent stray-texel
    /// artifacts that show up at small sizes.
    normalize_shape: bool = true,
    /// Geometrically determines each contour's correct winding direction,
    /// rather than relying on the single-sample-point `winding: .guess`
    /// heuristic - matters for multi-contour glyphs (most letters with a
    /// counter/hole).
    orient_contours: bool = true,
};

pub const Options = struct {
    type: AtlasType = .coverage,
    msdf: MsdfOptions = .{},
};

/// Rasterizes every glyph in `text` - via glyphBitmapSubpixel (the analytic
/// coverage rasterizer) for `.coverage`, or via the msdf module's Generator
/// for `.msdf` - packs the bitmaps into one shared atlas (via the guillotine
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
        .msdf => options.msdf.sdf_type.numChannels(),
    };

    const head = self.head orelse return error.MissingHead;
    const scale = pixel_size / @as(f32, @floatFromInt(head.unitsPerEm()));

    var generator: msdf = .create(self);
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

    for (text) |ch| {
        const glyph = self.getGlyphIndex(ch) orelse continue;

        switch (options.type) {
            .coverage => {
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

    const atlas_pixels = try gpa.alloc(u8, @as(usize, atlas_w) * atlas_h * channels);
    errdefer gpa.free(atlas_pixels);
    @memset(atlas_pixels, 0);
    for (rastered.items, placements) |g, p| {
        const off = g.bitmap_off orelse continue;
        const src = raster_pixels.items[off..];
        for (0..g.height) |row| {
            const src_row = src[row * g.width * channels ..][0 .. g.width * channels];
            const dst_off = ((p.y + row) * atlas_w + p.x) * channels;
            @memcpy(atlas_pixels[dst_off..][0 .. g.width * channels], src_row);
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
        .channels = channels,
        .glyphs = try out_glyphs.toOwnedSlice(gpa),
    };
}

/// One glyph's placement when drawn against a persistent `GlyphCache`:
/// unlike GlyphPlacement, `atlas_pos` here refers into the *cache's* texture
/// (`cache.pixels`/`width`/`height`), which the caller already holds and
/// which outlives any single buildTextAtlasCached call.
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

/// Like buildTextAtlas, but backed by a persistent, incrementally-populated
/// `GlyphCache` instead of rasterizing and packing every glyph from
/// scratch: a glyph already seen by `cache` (from this call or an earlier
/// one) is reused rather than re-rasterized. Also exposes `layout.Options`,
/// so callers can request word-wrapping - buildTextAtlas itself stays
/// single-line-only and untouched, since its bitmap-box rounding and this
/// function's layout-box rounding aren't guaranteed pixel-identical, and
/// buildTextAtlas has no test coverage to safely verify a shared code path
/// against.
///
/// `cache` must have been created for this same (font, pixel_size) pairing;
/// mixing sizes through one cache isn't supported (see GlyphCache's own
/// doc comment).
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

// ---
// Same synthetic-font approach as layout.zig's/GlyphCache.zig's tests.

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

test "buildTextAtlas with .msdf produces a correctly-shaped multi-channel atlas" {
    const alloc = zigstd.testAllocator;
    const font_bytes = try buildTestFont(alloc);
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

test "buildTextAtlasCached reuses a glyph's atlas placement across separate calls" {
    const alloc = zigstd.testAllocator;
    const font_bytes = try buildTestFont(alloc);
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
    const font_bytes = try buildTestFont(alloc);
    defer alloc.free(font_bytes);
    const sfnt = try Sfnt.init(alloc, font_bytes, .{ .check = true });

    var cache = try GlyphCache.init(alloc, 256);
    defer cache.deinit(alloc);

    var result = try buildTextAtlasCached(sfnt, alloc, &cache, "AB CD", 1000, 0, 0, .{ .wrap_width_px = 500 });
    defer result.deinit(alloc);

    try zigstd.testExpectEqual(@as(usize, 4), result.glyphs.len);
    try zigstd.testExpectTrue(result.glyphs[2].pos[1] > result.glyphs[0].pos[1]);
}
