const std = @import("std");
const zigstd = @import("../../../zigstd.zig");
const Allocator = zigstd.Allocator;

const Sfnt = @import("../sfnt/Sfnt.zig");
const aux_types = @import("../sfnt/aux_types.zig");
const pack = @import("../../pack.zig");

const coloring = @import("coloring.zig");
const EdgeColor = coloring.EdgeColor;
const EdgeSegment = @import("EdgeSegment.zig");
const math = @import("math.zig");
const Scanline = @import("Scanline.zig");
const Shape = @import("Shape.zig");
const Image = @import("../../image.zig").Image;

pub const error_correction = @import("error_correction.zig");

const Vertex = aux_types.Vertex;
const V2f64 = aux_types.V2f64;
const f64i = math.f64i;

const Generator = @This();

pub const FontMetrics = struct {
    line_height: f64,
    ascender: f64,
    descender: f64,
    underline_y: f64,
    underline_thickness: f64,
};

/// Em-scaled
pub const GlyphMetricsEm = struct {
    advance: f64,
    bearing_x: f64,
    bearing_y: f64,
    width: u16,
    height: u16,
};

pub const KerningPair = struct {
    codepoint_1: u21,
    codepoint_2: u21,
    x: f64,
    y: f64,
};

pub const GeneratedGlyph = struct {
    metrics: GlyphMetricsEm,
    pixels: []const u8,

    pub fn deinit(self: GeneratedGlyph, allocator: Allocator) void {
        allocator.free(self.pixels);
    }
};

pub const GeneratedAtlasGlyph = struct {
    glyph_data: GlyphMetricsEm,
    codepoint: u21,
    /// This is the glyph's unscaled, unnormalized
    /// bounding box on the atlas, with padding included.
    tex_bounds: pack.Rect,
};

pub const GeneratedAtlas = struct {
    glyphs: []const GeneratedAtlasGlyph,
    kernings: []const KerningPair,
    /// In case of `msdf10` you should reinterpret this as a `Msdf10Pixel` slice,
    /// like with 10-bit ABGR GPU formats (assuming a little endian host).
    pixels: []const u8,

    pub fn deinit(self: GeneratedAtlas, allocator: Allocator) void {
        allocator.free(self.glyphs);
        allocator.free(self.kernings);
        allocator.free(self.pixels);
    }
};

pub const SdfType = enum {
    sdf,
    psdf,
    msdf,
    mtsdf,
    /// Experimental: A packed BGR MSDF where each channel is 10-bit,
    /// with a 2-bit alpha channel that is ignored (set to u2 max).
    ///
    /// Can prove useful in place of MSDFs as native `R8G8B8_X` (and equivalent)
    /// format support is scarce. Additionally, 3-channel images are often padded
    /// to have an alignment of 4 bytes per pixel on a lot of hardware,
    /// which results in the final byte getting wasted on such formats.
    msdf10,

    pub fn numChannels(self: SdfType) u8 {
        return switch (self) {
            .sdf, .psdf => 1,
            .msdf, .msdf10 => 3,
            .mtsdf => 4,
        };
    }

    /// Bytes per output texel: msdf10 packs its three channels into four.
    pub fn bytesPerPixel(self: SdfType) u8 {
        return switch (self) {
            .msdf10 => 4,
            else => self.numChannels(),
        };
    }

    pub fn requiresColoring(self: SdfType) bool {
        return switch (self) {
            .msdf, .msdf10, .mtsdf => true,
            else => false,
        };
    }
};

pub const ColoringMethod = enum {
    simple,
    /// Only for use with ink trap fonts, as the coloring remains correct
    /// after removing the edges required for trapping ink.
    ink_trap,
    /// Performs the coloring based on edge distances.
    /// Somewhat slower than other methods, but it produces a better result most of the time.
    distance,
};

pub const Winding = enum {
    /// Attempts to figure out winding on its own, by checking
    /// the polarity of an OOB point's distance.
    guess,
    positive,
    negative,
};

pub const VarFontArgument = struct {
    name: []const u8,
    value: f64,
};

pub const Options = struct {
    sdf_type: SdfType,
    px_size: u16,
    px_range: u16,
    /// Has no effect if `sdf_type.requiresColoring()` is false.
    coloring_rng_seed: u64 = 0,
    /// The method with which to perform the MSDF 3-coloring.
    /// While the implementations are based on msdfgen, they're (intentionally)
    /// not equivalent, but should resolve corners similarly well.
    ///
    /// Has no effect if `sdf_type.requiresColoring()` is false.
    coloring_method: ColoringMethod = .distance,
    /// The angle which is considered to be a corner, in radians.
    corner_angle_threshold: f64 = 3.0,
    winding: Winding = .guess,
    /// Validates that the given (or generated) shapes' contours form a
    /// closed loop, with each edge connecting to each other properly.
    validate_shape: bool = false,
    normalize_shape: bool = false,
    orient_contours: bool = false,
    /// Requires `orient_contours` to be disabled.
    scanline_fill_rule: ?Scanline.FillRule = null,
    /// Only MSDFs (both their normal and their 10-bit versions) and MTSDFs can be error corrected.
    error_correction_opts: ?error_correction.Options = .{},
    /// The list of arguments to use if the given font has multiple masters.
    var_font_args: []const VarFontArgument = &.{},
    /// Whether to use async tasks over concurrent ones during atlas generation.
    /// Currently has no effect outside of atlas generation.
    disable_concurrency: bool = false,
};

pub const Msdf10Pixel = packed struct(u32) {
    r: u10 = 0,
    g: u10 = 0,
    b: u10 = 0,
    a: u2 = std.math.maxInt(u2),
};

/// A distance texel: one value per channel.
fn DistPixel(comptime ty: SdfType) type {
    return [ty.numChannels()]f64;
}

/// An output texel, as stored in GeneratedGlyph/GeneratedAtlas pixels.
fn OutPixel(comptime ty: SdfType) type {
    return [ty.bytesPerPixel()]u8;
}

sfnt: Sfnt,

/// `sfnt` should stay valid and available during the entire Generator
/// lifecycle: it's a borrowed view over the font's bytes, re-read on every
/// call, since Sfnt itself owns no allocated state (see Sfnt.init's doc
/// comment); there's nothing for a deinit() to release, and none exists.
pub fn init(sfnt: Sfnt) Generator {
    return .{ .sfnt = sfnt };
}

pub fn fontMetrics(self: *Generator) !FontMetrics {
    const head = self.sfnt.head orelse return error.MissingHead;
    const scale = 1.0 / f64i(head.unitsPerEm());
    const vm = self.sfnt.verticalMetrics();
    return .{
        .line_height = scale * f64i(vm.line_height),
        .ascender = scale * f64i(vm.ascender),
        .descender = scale * f64i(vm.descender),
        .underline_y = scale * f64i(vm.underline_position),
        .underline_thickness = scale * f64i(vm.underline_thickness),
    };
}

fn handleVarFont(self: *Generator, var_args: []const VarFontArgument) void {
    if (var_args.len == 0) return;

    if (self.sfnt.fvar == null) {
        std.log.warn("Var font args supplied, but the font has no fvar table", .{});
        return;
    }

    @panic("TODO: variable font axis coordinates are not yet applied to glyph outlines (gvar deformation is unimplemented)");
}

/// Walks a decomposed glyph outline (see Sfnt.glyphShape) into `shape`,
/// scaling design-unit coordinates into em space. Replaces what upstream
/// mist did via FT_Outline_Decompose + its four ftMoveTo/ftLineTo/ftConicTo/
/// ftCubicTo callbacks; same edge-building logic, just driven by a plain
/// loop instead of a C callback API.
fn appendVertices(shape: *Shape, allocator: Allocator, vertices: []const Vertex, scale: f64) !void {
    var pos: V2f64 = @splat(0.0);
    var contour: ?*Shape.Contour = null;
    for (vertices) |vtx| switch (vtx.type) {
        .vmove => {
            contour = try shape.contours.addOne(allocator);
            contour.?.* = .{};
            pos = scaledVertexVec(vtx.x, vtx.y, scale);
        },
        .vline => {
            const endpoint = scaledVertexVec(vtx.x, vtx.y, scale);
            if (!std.meta.eql(endpoint, pos)) {
                try contour.?.edges.append(allocator, .create(pos, endpoint, null, null, .all));
                pos = endpoint;
            }
        },
        .vcurve => {
            const endpoint = scaledVertexVec(vtx.x, vtx.y, scale);
            if (!std.meta.eql(endpoint, pos)) {
                try contour.?.edges.append(
                    allocator,
                    .create(pos, scaledVertexVec(vtx.cx, vtx.cy, scale), endpoint, null, .all),
                );
                pos = endpoint;
            }
        },
        .vcubic => {
            const endpoint = scaledVertexVec(vtx.x, vtx.y, scale);
            const c1 = scaledVertexVec(vtx.cx, vtx.cy, scale);
            const c2 = scaledVertexVec(vtx.cx1, vtx.cy1, scale);
            if (!std.meta.eql(endpoint, pos) or math.cross(c1 - endpoint, c2 - endpoint) != 0.0) {
                try contour.?.edges.append(allocator, .create(pos, c1, c2, endpoint, .all));
                pos = endpoint;
            }
        },
        _ => {},
    };
}

fn scaledVertexVec(x: i16, y: i16, scale: f64) V2f64 {
    return .{ f64i(x) * scale, f64i(y) * scale };
}

/// The result is under the caller's ownership (call `deinit()` or deallocate fields manually)
pub fn generateSingle(
    self: *Generator,
    allocator: Allocator,
    codepoint: u21,
    opts: *const Options,
) !GeneratedGlyph {
    self.handleVarFont(opts.var_font_args);

    const head = self.sfnt.head orelse return error.MissingHead;
    const scale = 1.0 / f64i(head.unitsPerEm());
    const glyph_index = self.sfnt.getGlyphIndex(codepoint) orelse return error.InvalidCodepoint;

    var shape: Shape = .{};
    defer shape.deinit(allocator);

    const outline_vertices = try self.sfnt.glyphShape(allocator, glyph_index);
    defer allocator.free(outline_vertices);
    try appendVertices(&shape, allocator, outline_vertices, scale);

    const metrics = self.sfnt.getGlyphMetrics(glyph_index);
    if (shape.contours.items.len == 0)
        return .{
            .metrics = .{
                .advance = scale * f64i(metrics.advance_width),
                .bearing_x = scale * f64i(metrics.bearing_x),
                .bearing_y = scale * f64i(metrics.bearing_y),
                .width = 0,
                .height = 0,
            },
            .pixels = &.{},
        };

    var contour_it = std.mem.reverseIterator(shape.contours.items);
    var i: isize = @intCast(shape.contours.items.len - 1);
    while (contour_it.next()) |contour| : (i -= 1)
        if (contour.edges.items.len == 0) {
            _ = shape.contours.swapRemove(@intCast(i));
        };

    if (opts.validate_shape and !shape.validate()) return error.InvalidShape;
    if (opts.orient_contours) try shape.orientContours(allocator);
    if (opts.normalize_shape) try shape.normalize(allocator);

    const px_size = f64i(opts.px_size);
    const px_range = f64i(opts.px_range) / px_size;

    var bounds = shape.calcBounds();
    if (bounds.left >= bounds.right or bounds.bottom >= bounds.top)
        bounds = .whole_frame;

    const bound_w = bounds.right - bounds.left;
    const bound_h = bounds.top - bounds.bottom;
    const w: u16 = @ceil((bound_w + px_range) * px_size);
    const h: u16 = @ceil((bound_h + px_range) * px_size);

    if (opts.winding == .negative or
        (opts.winding == .guess and findDistanceAt(
            .sdf,
            shape,
            .{
                bounds.left - px_range - bound_w - 1.0,
                bounds.bottom - px_range - bound_h - 1.0,
            },
            px_range,
        ) > 0)) for (shape.contours.items) |*contour| contour.reverse();

    return .{
        .metrics = .{
            .advance = scale * f64i(metrics.advance_width),
            .bearing_x = scale * f64i(metrics.bearing_x),
            .bearing_y = scale * f64i(metrics.bearing_y),
            .width = w,
            .height = h,
        },
        .pixels = try getSdfPixels(
            allocator,
            opts,
            w,
            h,
            &shape,
            .{
                bounds.left - px_range / 2.0,
                bounds.bottom - px_range / 2.0,
            },
        ),
    };
}

/// The result is under the caller's ownership (call `deinit()` or deallocate fields manually)
pub fn generateAtlas(
    self: *Generator,
    allocator: Allocator,
    io: std.Io,
    codepoints: []const u21,
    atlas_w: u16,
    atlas_h: u16,
    glyph_padding: u8,
    use_kerning: bool,
    opts: *const Options,
) !GeneratedAtlas {
    const glyphs = try allocator.alloc(GeneratedAtlasGlyph, codepoints.len);
    errdefer allocator.free(glyphs);

    var kernings: std.ArrayList(KerningPair) = .empty;
    errdefer kernings.deinit(allocator);
    kerning: {
        if (!use_kerning)
            break :kerning;

        if (!self.sfnt.hasKerning()) {
            std.log.warn("Kerning requested, but none were found in the font file.", .{});
            break :kerning;
        }

        const head = self.sfnt.head orelse break :kerning;
        const scale = 1.0 / f64i(head.unitsPerEm());
        for (codepoints) |codepoint_a| for (codepoints) |codepoint_b| {
            const idx_a = self.sfnt.getGlyphIndex(codepoint_a) orelse return error.InvalidCodepoint;
            const idx_b = self.sfnt.getGlyphIndex(codepoint_b) orelse return error.InvalidCodepoint;
            if (idx_a != idx_b) {
                const kern_x = self.sfnt.glyphKernAdvance(idx_a, idx_b);
                if (kern_x != 0)
                    try kernings.append(allocator, .{
                        .codepoint_1 = codepoint_a,
                        .codepoint_2 = codepoint_b,
                        .x = scale * f64i(kern_x),
                        .y = 0.0,
                    });
            }
        };
    }

    const id_rects = try allocator.alloc(pack.IdRect, codepoints.len);
    defer allocator.free(id_rects);

    const rect_pixels = try allocator.alloc([]const u8, codepoints.len);
    defer {
        for (rect_pixels) |px| allocator.free(px);
        allocator.free(rect_pixels);
    }
    @memset(rect_pixels, &.{});

    var process_group: std.Io.Group = .init;
    for (
        codepoints,
        glyphs,
        id_rects,
        rect_pixels,
        0..,
    ) |codepoint, *glyph, *id_rect, *rect_px, i| {
        id_rect.id = @intCast(i);
        const args = .{
            self,
            allocator,
            opts,
            codepoint,
            glyph,
            id_rect,
            rect_px,
            glyph_padding,
        };

        if (opts.disable_concurrency)
            process_group.async(io, processAtlasCodepoint, args)
        else
            try process_group.concurrent(io, processAtlasCodepoint, args);
    }
    try process_group.await(io);

    var pack_ctx: pack.Context = try .init(allocator, atlas_w, atlas_h, .{});
    defer pack_ctx.deinit();
    try pack.pack(
        pack.IdRect,
        &pack_ctx,
        id_rects,
        .{ .sortLessThanFn = struct {
            fn lessThan(_: void, a: pack.IdRect, b: pack.IdRect) bool {
                return @max(a.rect.w, a.rect.h) > @max(b.rect.w, b.rect.h);
            }
        }.lessThan },
    );

    const pixels: []u8 = switch (opts.sdf_type) {
        inline else => |ty| blk: {
            const Px = OutPixel(ty);
            const atlas: Image(Px) = try .init(allocator, atlas_w, atlas_h, @splat(0));
            for (id_rects) |id_rect| {
                const index: usize = @intCast(id_rect.id);
                const rect = id_rect.rect;
                glyphs[index].tex_bounds = rect;
                if (rect.w <= 0 or rect.h <= 0)
                    continue;

                atlas.blit(
                    @intCast(rect.x + glyph_padding),
                    @intCast(rect.y + glyph_padding),
                    @intCast(rect.w - glyph_padding * 2),
                    std.mem.bytesAsSlice(Px, rect_pixels[index]),
                );
            }
            break :blk std.mem.sliceAsBytes(atlas.data);
        },
    };
    errdefer allocator.free(pixels);

    return .{
        .glyphs = glyphs,
        .pixels = pixels,
        .kernings = if (use_kerning and kernings.items.len > 0)
            try kernings.toOwnedSlice(allocator)
        else
            &.{},
    };
}

fn processAtlasCodepoint(
    self: *Generator,
    allocator: Allocator,
    opts: *const Options,
    codepoint: u21,
    glyph: *GeneratedAtlasGlyph,
    id_rect: *pack.IdRect,
    rect_px: *[]const u8,
    padding: u8,
) std.Io.Cancelable!void {
    self.processAtlasCodepointInner(
        allocator,
        opts,
        glyph,
        id_rect,
        rect_px,
        padding,
        codepoint,
    ) catch {
        if (@errorReturnTrace()) |trace| std.debug.dumpErrorReturnTrace(trace);
        return std.Io.Cancelable.Canceled;
    };
}

fn processAtlasCodepointInner(
    self: *Generator,
    allocator: Allocator,
    opts: *const Options,
    glyph: *GeneratedAtlasGlyph,
    id_rect: *pack.IdRect,
    rect_px: *[]const u8,
    padding: u8,
    codepoint: u21,
) !void {
    const single_glyph = try self.generateSingle(allocator, codepoint, opts);
    const glyph_w = single_glyph.metrics.width;
    const glyph_h = single_glyph.metrics.height;
    if (glyph_w == 0 or glyph_h == 0) {
        glyph.* = .{
            .glyph_data = .{
                .advance = single_glyph.metrics.advance,
                .bearing_x = single_glyph.metrics.bearing_x,
                .bearing_y = single_glyph.metrics.bearing_y,
                .width = 0.0,
                .height = 0.0,
            },
            .codepoint = codepoint,
            .tex_bounds = .{ .w = 0, .h = 0 },
        };
        id_rect.rect = .{ .w = 0, .h = 0 };
        return;
    }

    rect_px.* = single_glyph.pixels;
    id_rect.rect = .{
        .w = glyph_w + padding * 2,
        .h = glyph_h + padding * 2,
    };
    glyph.* = .{
        .glyph_data = .{
            .advance = single_glyph.metrics.advance,
            .bearing_x = single_glyph.metrics.bearing_x,
            .bearing_y = single_glyph.metrics.bearing_y,
            .width = @intCast(id_rect.rect.w),
            .height = @intCast(id_rect.rect.h),
        },
        .codepoint = codepoint,
        // placed by the packer (`generate`)
        .tex_bounds = id_rect.rect,
    };
}

fn getSdfPixels(
    allocator: Allocator,
    opts: *const Options,
    w: u16,
    h: u16,
    shape: *Shape,
    tfm: V2f64,
) ![]const u8 {
    const px_size = f64i(opts.px_size);
    const px_range = f64i(opts.px_range) / px_size;

    switch (opts.sdf_type) {
        inline else => |ty| {
            const n = comptime ty.numChannels();
            var dist: Image(DistPixel(ty)) = try .init(allocator, w, h, null);
            defer dist.deinit(allocator);

            if (ty.requiresColoring()) switch (opts.coloring_method) {
                .simple => try coloring.colorSimple(allocator, opts.coloring_rng_seed, shape, opts.corner_angle_threshold),
                .ink_trap => try coloring.colorInkTrap(allocator, opts.coloring_rng_seed, shape, opts.corner_angle_threshold),
                .distance => try coloring.colorDistance(allocator, opts.coloring_rng_seed, shape, opts.corner_angle_threshold),
            };
            generate(ty, dist, px_size, shape.*, px_range, tfm);

            if (!opts.orient_contours)
                if (opts.scanline_fill_rule) |fill_rule| {
                    if (comptime n == 1)
                        try sdfSignCorrection(allocator, dist, px_size, shape.*, tfm, fill_rule)
                    else
                        try msdfSignCorrection(n, allocator, dist, px_size, shape.*, tfm, fill_rule);
                };

            if (comptime ty.requiresColoring())
                if (opts.error_correction_opts) |ec_opts|
                    try error_correction.correct(n, allocator, shape, ec_opts, opts.scanline_fill_rule != null, px_size, px_range, tfm, dist);

            const out: Image(OutPixel(ty)) = try .init(allocator, w, h, null);
            for (out.data, dist.data) |*o, d| {
                if (comptime ty == .msdf10)
                    o.* = std.mem.toBytes(Msdf10Pixel{
                        .r = @trunc(std.math.maxInt(u10) * std.math.clamp(d[0], 0.0, 1.0)),
                        .g = @trunc(std.math.maxInt(u10) * std.math.clamp(d[1], 0.0, 1.0)),
                        .b = @trunc(std.math.maxInt(u10) * std.math.clamp(d[2], 0.0, 1.0)),
                        .a = std.math.maxInt(u2),
                    })
                else for (o, d) |*v, c|
                    v.* = @trunc(std.math.maxInt(u8) * std.math.clamp(c, 0.0, 1.0));
            }
            return std.mem.sliceAsBytes(out.data);
        },
    }
}

fn sdfSignCorrection(
    allocator: Allocator,
    dist: Image([1]f64),
    scale: f64,
    shape: Shape,
    tfm: V2f64,
    fill_rule: Scanline.FillRule,
) !void {
    var scanline: Scanline = .{};
    defer scanline.intersections.deinit(allocator);
    for (0..dist.height) |y| {
        const row: u32 = @intCast(dist.height - y - 1);
        try shape.scanline(&scanline, (f64i(y) + 0.5) / scale + tfm[1], allocator);
        for (dist.row(row), 0..) |*d, x| {
            const distance = d[0];
            if ((distance > 0.5) != scanline.filled((f64i(x) + 0.5) / scale + tfm[0], fill_rule))
                d[0] = 1.0 - distance;
        }
    }
}

fn msdfSignCorrection(
    comptime n: usize,
    allocator: Allocator,
    dist: Image([n]f64),
    scale: f64,
    shape: Shape,
    tfm: V2f64,
    fill_rule: Scanline.FillRule,
) !void {
    var scanline: Scanline = .{};
    defer scanline.intersections.deinit(allocator);

    const w = dist.width;
    const h = dist.height;
    const match_map = try allocator.alloc(i32, @as(usize, w) * h);
    defer allocator.free(match_map);
    @memset(match_map, 0);

    var ambiguous = false;
    var match_idx: usize = 0;
    for (0..h) |y| {
        const row: u32 = @intCast(h - y - 1);
        try shape.scanline(&scanline, (f64i(y) + 0.5) / scale + tfm[1], allocator);
        for (dist.row(row), 0..) |*px, x| {
            const filled = scanline.filled((f64i(x) + 0.5) / scale + tfm[0], fill_rule);
            const distance = math.median(px[0], px[1], px[2]);
            if (distance == 0.5) {
                ambiguous = true;
            } else if ((distance > 0.5) != filled) {
                for (px[0..3]) |*v| v.* = 1.0 - v.*;
                match_map[match_idx] = -1;
            } else match_map[match_idx] = 1;
            if (n >= 4 and (px[3] > 0.5) != filled)
                px[3] = 1.0 - px[3];
            match_idx += 1;
        }
    }

    if (!ambiguous) return;

    match_idx = 0;
    for (0..h) |y| {
        const row: u32 = @intCast(h - y - 1);
        for (dist.row(row), 0..) |*px, x| {
            const match = match_map[match_idx];
            if (match == 0) {
                var neighbor_match: i32 = 0;
                if (x > 0) neighbor_match += match_map[match_idx - 1];
                if (x < w - 1) neighbor_match += match_map[match_idx + 1];
                if (y > 0) neighbor_match += match_map[match_idx - w];
                if (y < h - 1) neighbor_match += match_map[match_idx + w];
                if (neighbor_match < 0) {
                    for (px[0..3]) |*v| v.* = 1.0 - v.*;
                }
            }
            match_idx += 1;
        }
    }
}

fn pxRangeNorm(dist: f64, px_range: f64) f64 {
    return (dist + px_range / 2.0) / px_range;
}

pub fn findDistanceAt(
    comptime sdf_type: SdfType,
    shape: Shape,
    p: V2f64,
    px_range: f64,
) switch (sdf_type) {
    .sdf, .psdf => f64,
    inline .msdf, .msdf10, .mtsdf => |ty| [ty.numChannels()]f64,
} {
    const PsdfData = struct {
        dist: EdgeSegment.SignedDist = .init,
        edge: ?*const EdgeSegment = null,
        point_pos: EdgeSegment.PointPosition = .within_segment,
    };

    var true_ch: EdgeSegment.SignedDist = .init;
    var perp_ch: [if (sdf_type == .psdf) 1 else 3]PsdfData = @splat(.{});
    for (shape.contours.items) |contour| for (contour.edges.items) |*edge| {
        const dist, const point_pos = edge.signedDistance(p);

        switch (sdf_type) {
            .sdf, .mtsdf => {
                if (dist.lessThan(true_ch))
                    true_ch = dist;
            },
            .psdf => {
                if (dist.lessThan(perp_ch[0].dist)) perp_ch[0] = .{
                    .dist = dist,
                    .edge = edge,
                    .point_pos = point_pos,
                };
            },
            else => {},
        }

        if (sdf_type != .sdf and sdf_type != .psdf)
            for ([_]struct { channel: EdgeColor, target: *PsdfData }{
                .{ .channel = .red, .target = &perp_ch[0] },
                .{ .channel = .green, .target = &perp_ch[1] },
                .{ .channel = .blue, .target = &perp_ch[2] },
            }) |params|
                if (edge.color.hasChannel(params.channel) and dist.lessThan(params.target.dist)) {
                    params.target.* = .{
                        .dist = dist,
                        .edge = edge,
                        .point_pos = point_pos,
                    };
                };
    };

    if (sdf_type != .sdf) for (&perp_ch) |*psdf| {
        if (psdf.edge) |edge|
            edge.perpDistConvert(&psdf.dist, p, psdf.point_pos);
    };

    return switch (sdf_type) {
        .sdf => pxRangeNorm(true_ch.distance, px_range),
        .psdf => pxRangeNorm(perp_ch[0].dist.distance, px_range),
        .msdf, .msdf10 => .{
            pxRangeNorm(perp_ch[0].dist.distance, px_range),
            pxRangeNorm(perp_ch[1].dist.distance, px_range),
            pxRangeNorm(perp_ch[2].dist.distance, px_range),
        },
        .mtsdf => .{
            pxRangeNorm(perp_ch[0].dist.distance, px_range),
            pxRangeNorm(perp_ch[1].dist.distance, px_range),
            pxRangeNorm(perp_ch[2].dist.distance, px_range),
            pxRangeNorm(true_ch.distance, px_range),
        },
    };
}

fn generate(
    comptime sdf_type: SdfType,
    out: Image(DistPixel(sdf_type)),
    scale: f64,
    shape: Shape,
    px_range: f64,
    tfm: V2f64,
) void {
    for (0..out.height) |y| {
        const row: u32 = @intCast(out.height - y - 1);
        for (out.row(row), 0..) |*px, x| {
            const p = V2f64{
                (f64i(x) + 0.5),
                (f64i(y) + 0.5),
            } / math.splat(scale) + tfm;

            const dist = findDistanceAt(sdf_type, shape, p, px_range);
            px.* = if (comptime sdf_type.numChannels() == 1) .{dist} else dist;
        }
    }
}

test {
    _ = coloring;
    _ = EdgeSegment;
    _ = Scanline;
    _ = Shape;
    _ = math;
    _ = error_correction;
}
