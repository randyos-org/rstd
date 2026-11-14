const std = @import("std");
const bytes = @import("../../../bytes.zig");
const zigstd = @import("../../../zigstd.zig");
const tables = @import("tables/__root__.zig");
const types = @import("types.zig");
const aux_types = @import("aux_types.zig");
const charstring = @import("tables/cff/charstring.zig");
const glyf_shape = @import("tables/glyf/shape.zig");

const Vertex = aux_types.Vertex;
const Allocator = zigstd.Allocator;
const GlyphIndex = types.GlyphIndex;
const typeBEToNative = bytes.typeBEToNative;
const OffsetTable = tables.general.OffsetTable;
const TtcHeader = tables.general.TtcHeader;
const TableRecord = tables.general.TableRecord;
const TableTag = tables.general.TableTag;
const ScalerType = tables.general.ScalerType;

const Self = @This();

/// The raw font bytes every table views (tables/__root__.zig). Not owned: the
/// caller keeps them alive as long as this Sfnt (a mapped file, say).
font_data: []const u8 = &.{},
/// Required.
name: ?tables.NameTable = null,
/// Required.
cmap: ?tables.CmapTable = null,
/// Missing only in bitmap-only fonts.
head: ?tables.HeadTable = null,
/// Missing only in bitmap-only fonts.
hhea: ?tables.HheaTable = null,
/// Missing only in bitmap-only fonts.
hmtx: ?tables.HmtxTable = null,
/// Optional per TrueType, required per OpenType.
maxp: ?tables.MaxpTable = null,
/// Required.
post: ?tables.PostTable = null,

// Outlines: glyf and loca, or CFF.

/// Only if not CFF.
glyf: ?tables.GlyfTable = null,
/// Only with glyf.
loca: ?tables.LocaTable = null,
/// Only if not glyf.
cff: ?tables.CffTable = null,

// Optional tables.

/// Legacy pair kerning, the fallback for GPOS.
kern: ?tables.KernTable = null,
gpos: ?tables.GposTable = null,
/// For FreeType-exact ascender, descender and line height.
os2: ?tables.Os2Table = null,

// Variable fonts.

fvar: ?tables.FvarTable = null,
avar: ?tables.AvarTable = null,
/// TrueType outlines only.
gvar: ?tables.GvarTable = null,
/// TrueType outlines only. Needs fvar's axisCount, so it is resolved after
/// the table-record loop, as hmtx and loca are.
cvar: ?tables.CvarTable = null,
hvar: ?tables.HvarTable = null,
vvar: ?tables.VvarTable = null,
mvar: ?tables.MvarTable = null,
/// Required for variable fonts.
stat: ?tables.StatTable = null,

pub const SfntError = error{
    MissingRequiredTable,
    IndexMapMissing,
    FontIndexOutOfRange,
} || TableTag.Error || ScalerType.Error || bytes.BytesError || Allocator.Error || tables.CffError;

const InitOptions = struct {
    check: bool = false,
    /// Which font of a collection (.ttc/.otc) to load; 0 for a single font.
    font_index: u32 = 0,
    /// Print each skipped (unrecognized) table tag to stderr, as png.zig's
    /// `DecodeOptions.debug` does.
    debug: bool = true,
};

/// How many fonts `font_data` holds: a collection's count, else 1.
pub fn numberOfFonts(font_data: []const u8) bytes.BytesError!u32 {
    try bytes.checkSize(font_data, 0, u32);
    if (bytes.getU32BE(font_data) != TtcHeader.tag) return 1;
    try bytes.checkSize(font_data, 0, TtcHeader);
    return typeBEToNative(TtcHeader, font_data).num_fonts;
}

/// Where font `index`'s OffsetTable starts: 0 for a single font, else the
/// collection's entry. Table records hold absolute file offsets either way.
pub fn fontOffset(font_data: []const u8, index: u32) SfntError!usize {
    if (index >= try numberOfFonts(font_data)) return error.FontIndexOutOfRange;
    if (bytes.getU32BE(font_data) != TtcHeader.tag) return 0;
    try bytes.checkArray(font_data, @sizeOf(TtcHeader), @as(usize, index) + 1, u32);
    return bytes.getU32BE(font_data[@sizeOf(TtcHeader) + @as(usize, index) * @sizeOf(u32) ..]);
}

/// Every table field above is a plain view over table_data slices borrowed
/// from font_data; none of them own any heap memory (cmap subtables are the
/// one exception, but those are resolved and owned separately via
/// CmapTable.getSubtable(), not stored here). The only allocation init()
/// itself performs is a scratch table-record array, freed before it
/// returns; so there's nothing for a deinit() to release, and none exists.
pub fn init(alloc: Allocator, font_data: []const u8, options: InitOptions) SfntError!Self {
    const font_start = try fontOffset(font_data, options.font_index);
    if (options.check) try bytes.checkSize(font_data, font_start, OffsetTable);
    const offset_table = typeBEToNative(OffsetTable, font_data[font_start..]);
    if (options.check) _ = try ScalerType.fromInt(offset_table.scaler);

    const tables_start = font_start + @sizeOf(OffsetTable);
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
            if (options.debug) std.debug.print("Unrecognized table tag: {s}\n", .{record.tag.strip()});
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
// Higher-level API, shaped after what mist calls on a FreeType Face
// (ref/mist/src/Generator.zig). Values are in font design units (divide by
// head.unitsPerEm()), as with FreeType's FT_LOAD_NO_SCALE.

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
pub fn getGlyphMetrics(self: Self, glyph: GlyphIndex) aux_types.GlyphMetricsDesign {
    const advance_width = if (self.hmtx) |hmtx| hmtx.getMetrics(glyph).advance_width else 0;
    const bbox = self.getGlyphBoundingBox(glyph) orelse
        return .{ .advance_width = advance_width, .bearing_x = 0, .bearing_y = 0 };
    return .{ .advance_width = advance_width, .bearing_x = @intCast(bbox.x_min), .bearing_y = @intCast(bbox.y_max) };
}

/// Resolves a Unicode codepoint to a glyph index via the font's preferred
/// Unicode cmap subtable. Prefers a full-repertoire subtable (format 12/13,
/// covering codepoints past the BMP) over a BMP-only one when both exist,
/// since the former is always a superset of the latter's mappings. Returns
/// null if the font has no Unicode subtable, or the subtable has no mapping
/// for this codepoint.
pub fn getGlyphIndex(self: Self, codepoint: u32) ?GlyphIndex {
    const cmap = self.cmap orelse return null;
    const record = cmap.findUnicodeFullRecord() orelse cmap.findUnicodeBmpRecord() orelse return null;
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

const setU16BE = bytes.setU16BE;
const setU32BE = bytes.setU32BE;

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

/// One font with only a `head` table, its OffsetTable at `base`, its head
/// at `head_start` (both absolute, as in a collection).
fn writeHeadOnlyFont(buf: []u8, base: usize, head_start: usize, units_per_em: u16) void {
    setU32BE(buf[base..], 0x00010000); // scaler
    setU16BE(buf[base + 4 ..], 1); // num_tables
    const rec = base + 12;
    @memcpy(buf[rec..][0..4], "head");
    setU32BE(buf[rec + 8 ..], @intCast(head_start));
    setU32BE(buf[rec + 12 ..], 54);
    setU16BE(buf[head_start + 18 ..], units_per_em);
}

test "Sfnt loads each font of a collection by font_index" {
    const font_size = 12 + 16 + 54;
    const fonts_start = 12 + 2 * 4;
    var buf: [fonts_start + 2 * font_size]u8 = @splat(0);
    setU32BE(&buf, TtcHeader.tag);
    setU16BE(buf[4..], 2); // majorVersion
    setU32BE(buf[8..], 2); // numFonts
    for (0..2, [_]u16{ 1000, 2048 }) |i, upem| {
        const base = fonts_start + i * font_size;
        setU32BE(buf[12 + 4 * i ..], @intCast(base)); // tableDirectoryOffsets[i]
        writeHeadOnlyFont(&buf, base, base + 12 + 16, upem);
    }

    try zigstd.testExpectEqual(@as(u32, 2), try numberOfFonts(&buf));
    for (0..2, [_]u16{ 1000, 2048 }) |i, upem| {
        const sfnt = try init(zigstd.testAllocator, &buf, .{ .check = true, .font_index = @intCast(i) });
        try zigstd.testExpectEqual(upem, sfnt.head.?.unitsPerEm());
    }
    try zigstd.testExpectError(error.FontIndexOutOfRange, init(zigstd.testAllocator, &buf, .{ .font_index = 2 }));
}

test "Sfnt treats a plain font as a one-font collection" {
    var buf: [12 + 16 + 54]u8 = @splat(0);
    writeHeadOnlyFont(&buf, 0, 12 + 16, 1000);

    try zigstd.testExpectEqual(@as(u32, 1), try numberOfFonts(&buf));
    try zigstd.testExpectEqual(@as(usize, 0), try fontOffset(&buf, 0));
    try zigstd.testExpectError(error.FontIndexOutOfRange, fontOffset(&buf, 1));
}

test "Sfnt rejects a collection header cut short" {
    var buf: [12 + 4]u8 = @splat(0);
    setU32BE(&buf, TtcHeader.tag);
    setU32BE(buf[8..], 2); // numFonts, but room for only one offset
    try zigstd.testExpectEqual(@as(usize, 0), try fontOffset(&buf, 0));
    try zigstd.testExpectError(error.EndOfStream, fontOffset(&buf, 1));
}
