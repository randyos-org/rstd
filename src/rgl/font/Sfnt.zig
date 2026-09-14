const std = @import("std");
const bytes = @import("../../bytes.zig");
const zigstd = @import("../../zigstd.zig");
const tables = @import("tables/__root__.zig");
const types = @import("stream_types.zig");
const aux_types = @import("aux_types.zig");
const render = @import("render/__root__.zig");
const charstring = @import("tables/cff/charstring.zig");
const glyf_shape = @import("tables/glyf/shape.zig");

const Vertex = aux_types.Vertex;
const bounds = render.bounds;
const Allocator = zigstd.Allocator;
const GlyphIndex = types.GlyphIndex;
const typeBEToNative = bytes.typeBEToNative;
const setU16BE = bytes.setU16BE;
const setU32BE = bytes.setU32BE;
const OffsetTable = tables.general.OffsetTable;
const TableRecord = tables.general.TableRecord;
const TableTag = tables.general.TableTag;
const ScalerType = tables.general.ScalerType;

const Self = @This();

/// The raw font bytes every table below borrows from. Not owned: assumed to
/// be backed by memory the caller keeps alive for at least as long as this
/// Sfnt (e.g. a memory-mapped file region), since nothing here copies out of
/// it.
font_data: []const u8 = &.{},
/// Naming table.  Required.
name: ?tables.NameTable = null,
/// Character to glyph mapping.  Required.
cmap: ?tables.CmapTable = null,
/// Font header. May be missing in bitmap-only fonts, otherwise required.
head: ?tables.HeadTable = null,
/// Horizontal header. May be missing in bitmap-only fonts, otherwise required.
hhea: ?tables.HheaTable = null,
/// Horizontal metrics. May be missing in bitmap-only fonts, otherwise required.
hmtx: ?tables.HmtxTable = null,
/// Maximum profile. Optional per TTF, required per OpenType spec.
maxp: ?tables.MaxpTable = null,
/// PostScript table. Required.
post: ?tables.PostTable = null,

// ---
// Classic glyph-based TrueType tables. Only required if not CFF.

/// Glyph data. Only required if not CFF.
glyf: ?tables.GlyfTable = null,
/// Index to location. Only required if glyf is present
loca: ?tables.LocaTable = null,

// ---
// CFF TABLES: only required if not glyf

/// Compact Font Format
cff: ?tables.CffTable = null,

// ---
// OPTIONAL TTF TABLES:

/// Kerning. Legacy pair-kerning fallback; GPOS is the modern source.
kern: ?tables.KernTable = null,

// ---
// Advanced typographic tables

/// Glyph positioning data
gpos: ?tables.GposTable = null,

// ---
// Platform-specific tables

/// OS/2 and Windows metrics. Needed for exact ascender/descender/line-height
/// parity with FreeType
os2: ?tables.Os2Table = null,

// ---
// OpenType variations tables (variable fonts)

/// Font variations: the axes and named instances a variable font supports.
fvar: ?tables.FvarTable = null,
/// Axis variations: optional remapping of each axis's default coordinate
/// normalization.
avar: ?tables.AvarTable = null,
/// Glyph variations (TrueType outlines only).
gvar: ?tables.GvarTable = null,
/// CVT variations (TrueType outlines only). Not self-describing: needs
/// fvar's axisCount, so it's resolved after the table-record loop below,
/// the same way hmtx/loca are.
cvar: ?tables.CvarTable = null,
/// Horizontal metrics variations.
hvar: ?tables.HvarTable = null,
/// Vertical metrics variations.
vvar: ?tables.VvarTable = null,
/// Metrics variations (OS/2, hhea, gasp, post, vhea deltas).
mvar: ?tables.MvarTable = null,
/// Style attributes. Required for variable fonts.
stat: ?tables.StatTable = null,

pub const SfntError = error{
    MissingRequiredTable,
    IndexMapMissing,
} || TableTag.Error || ScalerType.Error || bytes.BytesError || Allocator.Error || tables.CffError;

const InitOptions = struct {
    check: bool = false,
};

/// Every table field above is a plain view over table_data slices borrowed
/// from font_data; none of them own any heap memory (cmap subtables are the
/// one exception, but those are resolved and owned separately via
/// CmapTable.getSubtable(), not stored here). The only allocation init()
/// itself performs is a scratch table-record array, freed before it
/// returns - so there's nothing for a destroy() to release, and none exists.
pub fn init(alloc: Allocator, font_data: []const u8, options: InitOptions) SfntError!Self {
    if (options.check) try bytes.checkSize(font_data, 0, OffsetTable);
    const offset_table = typeBEToNative(OffsetTable, font_data);
    if (options.check) _ = try ScalerType.fromInt(offset_table.scaler);

    const tables_start = @sizeOf(OffsetTable);
    const num_tables = offset_table.num_tables;
    if (options.check) try bytes.checkArray(font_data, tables_start, num_tables, TableRecord);
    const table_recs = try bytes.allocArrayBEToNative(alloc, TableRecord, num_tables, font_data[tables_start..]);
    defer alloc.free(table_recs);

    var self: Self = .{ .font_data = font_data };

    // hmtx and loca aren't self-describing (their record counts come from
    // other tables), and table order in the font file isn't guaranteed, so
    // their own table_data is just captured here and parsed after the loop
    // once their dependencies are available.
    var hmtx_table_data: ?[]const u8 = null;
    var loca_table_data: ?[]const u8 = null;
    var cvar_table_data: ?[]const u8 = null;

    for (table_recs) |record| {
        const table_data = try bytes.safeSlice(font_data, record.offset, record.length, options.check);
        // Real-world fonts routinely carry tables this parser doesn't model
        // (GSUB, gasp, meta, hdmx, DSIG, vendor-private tables, ...) -
        // unrecognized tags are skipped rather than failing the whole load.
        const table_tag = TableTag.fromTagStr(record.tag) catch {
            std.debug.print("Unrecognized table tag: {s}\n", .{record.tag.strip()});
            continue;
        };
        switch (table_tag) {
            .cmap => self.cmap = try tables.CmapTable.init(table_data, options.check),
            .head => self.head = try tables.HeadTable.init(table_data, options.check),
            .hhea => self.hhea = try tables.HheaTable.init(table_data, options.check),
            .maxp => self.maxp = try tables.MaxpTable.init(table_data, options.check),
            .name => self.name = try tables.NameTable.init(table_data, options.check),
            .post => self.post = try tables.PostTable.init(table_data, options.check),
            .hmtx => hmtx_table_data = table_data,
            .loca => loca_table_data = table_data,
            .glyf => self.glyf = tables.GlyfTable.init(table_data),
            .kern => self.kern = try tables.KernTable.init(table_data, options.check),
            .GPOS => self.gpos = try tables.GposTable.init(table_data, options.check),
            .CFF => self.cff = try tables.CffTable.init(table_data, options.check),
            .@"OS/2" => self.os2 = try tables.Os2Table.init(table_data, options.check),
            .fvar => self.fvar = try tables.FvarTable.init(table_data, options.check),
            .avar => self.avar = try tables.AvarTable.init(table_data, options.check),
            .gvar => self.gvar = try tables.GvarTable.init(table_data, options.check),
            .cvar => cvar_table_data = table_data,
            .HVAR => self.hvar = try tables.HvarTable.init(table_data, options.check),
            .VVAR => self.vvar = try tables.VvarTable.init(table_data, options.check),
            .MVAR => self.mvar = try tables.MvarTable.init(table_data, options.check),
            .STAT => self.stat = try tables.StatTable.init(table_data, options.check),
            else => {},
        }
    }

    if (hmtx_table_data) |table_data| {
        if (self.hhea) |hhea| {
            if (self.maxp) |maxp| {
                self.hmtx = try tables.HmtxTable.init(
                    table_data,
                    hhea.numOfLongHorMetrics(),
                    maxp.numGlyphs(),
                    options.check,
                );
            }
        }
    }

    if (loca_table_data) |table_data| {
        if (self.head) |head| {
            if (self.maxp) |maxp| {
                self.loca = try tables.LocaTable.init(
                    table_data,
                    maxp.numGlyphs(),
                    head.indexToLocFormat() == 1,
                    options.check,
                );
            }
        }
    }

    if (cvar_table_data) |table_data| {
        if (self.fvar) |fvar| {
            self.cvar = try tables.CvarTable.init(table_data, fvar.axisCount(), options.check);
        }
    }

    return self;
}

// ---
// Higher-level API for consumers - namely a future Zig-native replacement
// for FreeType inside mist, an MSDF atlas generator (see
// ref/mist/src/Generator.zig for exactly what it calls on a FreeType Face,
// which is where this interface shape comes from). These signatures are
// locked in now, while the requirements are fresh, but several bodies are
// intentionally unimplemented placeholders: getting every table loaded
// correctly comes first, the resolution algorithms behind these come later.
// Every value here is in font design units (unscaled - divide by
// head.unitsPerEm() yourself), matching FreeType's FT_LOAD_NO_SCALE
// convention, which is the convention mist itself already loads glyphs with.

pub fn verticalMetrics(self: Self) aux_types.VMetrics {
    var ascender: types.FWord = 0;
    var descender: types.FWord = 0;
    var line_gap: types.FWord = 0;

    const use_typo = if (self.os2) |os2| os2.version() >= 4 and os2.fsSelection().use_typo_metrics else false;
    if (use_typo or self.hhea == null) {
        if (self.os2) |os2| {
            ascender = os2.sTypoAscender();
            descender = os2.sTypoDescender();
            line_gap = os2.sTypoLineGap();
        }
    } else if (self.hhea) |hhea| {
        ascender = hhea.ascent();
        descender = hhea.descent();
        line_gap = hhea.lineGap();
    }
    var vmtx: aux_types.VMetrics = .{
        .line_height = ascender - descender + line_gap,
        .ascender = ascender,
        .descender = descender,
    };

    if (self.post) |post| {
        vmtx.underline_position = post.underlinePosition();
        vmtx.underline_thickness = post.underlineThickness();
    }

    return vmtx;
}

/// Scale factor mapping font design units to a target pixel height, using
/// ascender-descender span (not unitsPerEm) as the reference height.
pub fn scaleForPixelHeight(self: Self, height: f32) f32 {
    const m = self.verticalMetrics();
    const span: f32 = @floatFromInt(m.ascender - m.descender);
    return height / span;
}

/// Format-agnostic glyph bounding box: dispatches like `glyphShape`.
pub fn getGlyphBoundingBox(self: Self, glyph: GlyphIndex) ?aux_types.GlyphBoundingBox {
    if (self.cff) |cff| return charstring.glyphBoxT2(cff, glyph);
    if (self.glyf) |glyf| {
        const loca = self.loca orelse return null;
        return glyf_shape.glyphBoxTT(loca, glyf, glyph);
    }
    return null;
}

/// Combines hmtx's advance width with the glyph outline's own bounding box.
pub fn getGlyphMetrics(self: Self, glyph: GlyphIndex) aux_types.GlyphMetrics {
    const advance_width = if (self.hmtx) |hmtx| hmtx.getMetrics(glyph).advance_width else 0;
    const bbox = self.getGlyphBoundingBox(glyph) orelse
        return .{ .advance_width = advance_width, .bearing_x = 0, .bearing_y = 0 };
    return .{ .advance_width = advance_width, .bearing_x = @intCast(bbox.x_min), .bearing_y = @intCast(bbox.y_max) };
}

/// Resolves a Unicode codepoint to a glyph index via the font's preferred
/// Unicode cmap subtable. Returns null if the font has no such subtable, or
/// the subtable has no mapping for this codepoint.
pub fn getGlyphIndex(self: Self, codepoint: u32) ?GlyphIndex {
    const cmap = self.cmap orelse return null;
    const record = cmap.findUnicodeBmpRecord() orelse return null;
    const subtable = cmap.getSubtable(record, false) catch return null;
    const glyph = subtable.getGlyphIndex(codepoint);
    return if (glyph == .notdef) null else glyph;
}

/// Format-agnostic glyph outline extraction: dispatches to the `glyf`-table
/// path or the CFF/Type2 charstring path depending on which outline format
/// this font uses.
pub fn glyphShape(self: Self, gpa: Allocator, glyph: GlyphIndex) ![]Vertex {
    if (self.cff) |cff| return charstring.glyphShapeT2(cff, gpa, glyph);
    if (self.glyf) |glyf| {
        const loca = self.loca orelse return &.{};
        return glyf_shape.glyphShapeTT(loca, glyf, gpa, glyph);
    }
    return &.{};
}

/// Whether this font has any kerning data at all (kern and/or GPOS),
/// mirroring FreeType's FT_FACE_FLAG_KERNING.
pub fn hasKerning(self: Self) bool {
    return self.kern != null or self.gpos != null;
}

/// An additional amount to advance the horizontal coordinate between the two
/// provided glyphs.
pub fn glyphKernAdvance(self: Self, a: GlyphIndex, b: GlyphIndex) i16 {
    if (self.gpos) |gpos| return gpos.glyphKernAdvance(a, b);
    if (self.kern) |kern| return kern.getKernValue(a, b) orelse 0;
    return 0;
}

/// Caller owns returned memory.
pub fn glyphBitmap(
    self: Self,
    gpa: Allocator,
    /// Appended to the list.
    /// Stored left-to-right, top-to-bottom. 8 bits per pixel. 0 is
    /// transparent, 255 is opaque.
    pixels: *std.ArrayListUnmanaged(u8),
    glyph: GlyphIndex,
    scale_x: f32,
    scale_y: f32,
) !aux_types.GlyphBitmap {
    return self.glyphBitmapSubpixel(gpa, pixels, glyph, scale_x, scale_y, 0, 0);
}

/// Caller owns returned memory.
pub fn glyphBitmapSubpixel(
    self: Self,
    gpa: Allocator,
    /// Appended to the list.
    /// Stored left-to-right, top-to-bottom. 8 bits per pixel. 0 is
    /// transparent, 255 is opaque.
    pixels: *std.ArrayListUnmanaged(u8),
    glyph: GlyphIndex,
    scale_x: f32,
    scale_y: f32,
    shift_x: f32,
    shift_y: f32,
) !aux_types.GlyphBitmap {
    const vertices = try self.glyphShape(gpa, glyph);
    defer gpa.free(vertices);

    std.debug.assert(scale_x != 0);
    std.debug.assert(scale_y != 0);

    const bbox = self.getGlyphBoundingBox(glyph) orelse return .empty;
    const box = bounds.glyphBitmapBoxSubpixel(bbox, scale_x, scale_y, shift_x, shift_y);

    const w: u32 = @intCast(box.x1 - box.x0);
    const h: u32 = @intCast(box.y1 - box.y0);

    if (w == 0 or h == 0) return .empty;

    var gbm: render.rasterize.Bitmap = .{
        .w = w,
        .h = h,
        .stride = w,
        .pixels = try pixels.addManyAsSlice(gpa, w * h),
    };
    errdefer pixels.shrinkRetainingCapacity(pixels.items.len - gbm.pixels.len);

    try render.rasterize.rasterize(gpa, &gbm, 0.35, vertices, scale_x, scale_y, shift_x, shift_y, box.x0, box.y0, true);

    return .{
        .width = @intCast(gbm.w),
        .height = @intCast(gbm.h),
        .off_x = @intCast(box.x0),
        .off_y = @intCast(box.y0),
    };
}

pub const glyphBitmapBoxSubpixel = bounds.glyphBitmapBoxSubpixel;
pub const glyphBitmapBox = bounds.glyphBitmapBox;
pub const buildTextAtlas = render.atlas.buildTextAtlas;

test "Sfnt.init rejects a font with a bad scaler tag" {
    var buf: [12]u8 = @splat(0);
    setU32BE(&buf, 0xDEADBEEF);
    try zigstd.testExpectError(error.UnsupportedScalerType, init(zigstd.testAllocator, buf[0..], .{ .check = true }));
}

test "Sfnt.init frees its scratch table-record buffer even when a later table fails validation" {
    const offset_table_size = 12;
    const rec_size = 16;
    const head_start = offset_table_size + 2 * rec_size; // 44
    const head_size = 54;
    const hhea_start = head_start + head_size; // 98
    const hhea_size = 10; // too short for hhea's 36-byte layout
    const total = hhea_start + hhea_size;

    var buf: [total]u8 = @splat(0);
    setU32BE(&buf, 0x00010000); // scaler
    setU16BE(buf[4..], 2); // num_tables

    @memcpy(buf[offset_table_size..][0..4], "head");
    setU32BE(buf[offset_table_size + 8 ..], head_start);
    setU32BE(buf[offset_table_size + 12 ..], head_size);

    @memcpy(buf[offset_table_size + rec_size ..][0..4], "hhea");
    setU32BE(buf[offset_table_size + rec_size + 8 ..], hhea_start);
    setU32BE(buf[offset_table_size + rec_size + 12 ..], hhea_size);

    // zigstd.testAllocator asserts no leaks on its own; there's nothing else
    // to release since every table is a borrowed view, so a clean error
    // return here is the whole test.
    try zigstd.testExpectError(error.EndOfStream, init(zigstd.testAllocator, &buf, .{ .check = true }));
}
