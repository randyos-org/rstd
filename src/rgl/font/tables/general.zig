const std = @import("std");

const types = @import("../types.zig");
const bytes = @import("../../../bytes.zig");

const readIntBE = bytes.readIntBE;
const writeIntBE = bytes.writeIntBE;
const TagStr = types.TagStr;
const Fixed = types.Fixed;

const PlatformId = enum(u16) {
    unicode = 0,
    mac = 1,
    iso = 2,
    microsoft = 3,
};

const MicrosoftEncodingId = enum(u16) {
    symbol = 0,
    unicode_bmp = 1,
    shiftjis = 2,
    unicode_full = 10,
};

/// The list of all known table tags published in either the TTF or OTF specs.
///
/// See:
/// - https://learn.microsoft.com/en-us/typography/opentype/spec/
/// - https://developer.apple.com/fonts/TrueType-Reference-Manual/
pub const TableTag = enum {
    const Self = @This();

    /// Naming table. Required.
    name,
    /// Character to glyph mapping.  Required.
    cmap,
    /// Font header. May be missing in bitmap-only fonts, otherwise required.
    head,
    /// Horizontal header. May be missing in bitmap-only fonts, otherwise required.
    hhea,
    /// Horizontal metrics. May be missing in bitmap-only fonts, otherwise required.
    hmtx,
    /// Maximum profile. Optional per TTF, required per OpenType spec.
    maxp,
    /// PostScript table. Required.
    post,

    // ---
    // Classic glyph-based TrueType tables. Only required if not CFF.

    /// Glyph data. Only required if not CFF.
    glyf,
    /// Index to location. Only required if glyf is present
    loca,

    // ---
    // CFF TABLES: only required if not glyf

    /// Compact Font Format.
    CFF,
    /// Compact Font Format 2.
    CFF2,
    /// Vertical origin.  Optional.
    VORG,

    // ---
    // Platform-specific tables

    /// OS/2, required on OS/2 and Windows platforms only
    @"OS/2",

    // ---
    // OPTIONAL TTF TABLES:

    /// control value
    cvt,
    /// font program
    fpgm,
    /// horizontal device metrics, used on Mac only
    hdmx,
    /// kerning
    kern,
    /// control value program
    prep,
    /// Grid-fitting/scan-conversion
    gasp,

    // ---
    // SVG tables

    /// SVG table
    SVG,

    // ---
    // bitmap tables

    /// Embedded bitmap data
    EBDT,
    /// Embedded bitmap location data
    EBLC,
    /// Embedded bitmap scaling data
    EBSC,
    /// Color bitmap data
    CBDT,
    /// Color bitmap location data
    CBLC,
    /// Standard bitmap graphics
    sbix,

    // ---
    // Advanced typographic tables

    /// Baseline data
    BASE,
    /// Glyph definition data
    GDEF,
    /// Glyph positioning data
    GPOS,
    /// Glyph substitution data
    GSUB,
    /// Justification data
    JSTF,
    /// Math layout data
    MATH,

    // ---
    // OpenType variations tables

    /// Axis variations
    avar,
    /// CVT variations (TrueType outlines only)
    cvar,
    /// Font variations
    fvar,
    /// Glyph variations (TrueType outlines only)
    gvar,
    /// Horizontal metrics variations
    HVAR,
    /// Metrics variations
    MVAR,
    /// Style attributes (required for variable fonts, optional for non-variable fonts)
    STAT,
    /// Vertical metrics variations
    VVAR,

    // ---
    // Color tables

    /// Color table
    COLR,
    /// Color palette table
    CPAL,

    // ---
    // other OpenType tables

    /// Digital signature
    DSIG,
    /// Linear threshold data
    LTSH,
    /// Merge
    MERG,
    /// Metadata
    meta,
    /// PCL 5 data
    PCLT,
    /// Vertical device metrics
    VDMX,
    /// Vertical Metrics header
    vhea,
    /// Vertical Metrics
    vmtx,

    pub fn toInt(id: Self) u32 {
        var array4: TagStr = @splat(" ");
        @memcpy(&array4, @tagName(id));
        return readIntBE(u32, &array4);
    }

    pub fn fromInt(v: u32) !Self {
        var tag: TagStr = undefined;
        writeIntBE(u32, &tag, v);
        return fromTagStr(tag);
    }

    pub fn fromTagStr(tag: TagStr) !Self {
        return std.meta.stringToEnum(Self, tag.strip()) orelse error.UnsupportedTag;
    }
};

pub const ScalerType = enum {
    const Self = @This();

    ttf,
    ps,
    otf,

    pub fn fromInt(v: u32) !Self {
        const alttrue = TagStr.fromInt(0x00010000);
        const truestr = TagStr.fromStr("true");
        const typ1 = TagStr.fromStr("typ1");
        const otto = TagStr.fromStr("OTTO");
        switch (v) {
            alttrue, truestr => return .ttf,
            typ1 => return .ps,
            otto => return .otf,
            else => return error.UnsupportedScalerType,
        }
    }
};

const OffsetTable = packed struct {
    scaler: u32,
    num_tables: u16,
    search_range: u16,
    entry_selector: u16,
    range_shift: u16,
};

const TableRecord = extern struct {
    const Self = @This();

    tag: TagStr,
    checksum: u32,
    offset: u32,
    length: u32,

    pub fn calcTableChecksum(self: *Self, table_data: *[]const u8) u32 {
        var sum: u32 = 0;
        const nLongs = (self.length + 3) / 4;
        for (0..nLongs) |i| sum += @as(u32, table_data[i * 4]);
        return sum;
    }
};

const HeadTable = packed struct {
    const Self = @This();

    version: Fixed,
    font_revision: Fixed,
    checksum_adjustment: u32,
    magic_number: u32,
    flags: u16,
    units_per_em: u16,
    created: i64,
    modified: i64,
    x_min: i16,
    y_min: i16,
    x_max: i16,
    y_max: i16,
    mac_style: u16,
    lowest_rec_ppem: u16,
    font_direction_hint: i16,
    index_to_loc_format: i16,
    glyph_data_format: i16,

    pub fn calcTableChecksum(self: *Self, table_data: *[]const u8) u32 {
        var sum: u32 = self.version + self.font_revision;
        // assume checksum_adjustment is 0 for this calculation
        const nLongs = (@sizeOf(Self) + 3) / 4 - 3; // subtract 3 for the two Fixed fields and checksum_adjustment
        for (0..nLongs) |i| sum += @as(u32, table_data[12 + i * 4]);
        return sum;
    }

    /// When calculating the raw_font_checksum_adj0, the checksum_adjustment field is assumed to be 0.
    /// You may need to actually temporarily set it to zero since we don't necessarily know the packing scheme,
    /// so set it to zero, calculate the checksum for the entire font byte array with the zero checksum,
    /// and pass that value here.  Do not forget to then re-set the checksum_adjustment following the calculation.
    pub fn calcChecksumAdj(raw_font_checksum_adj0: u32) u32 {
        return 0xB1B0AFBA - raw_font_checksum_adj0;
    }
};

const MaxpTable = packed struct {
    version: Fixed,
    num_glyphs: u16,
    max_points: u16,
    max_contours: u16,
    max_component_points: u16,
    max_component_contours: u16,
    max_zones: u16,
    max_twilight_points: u16,
    max_storage: u16,
    max_function_defs: u16,
    max_instruction_defs: u16,
    maxStackElements: u16,
    maxSizeOfInstructions: u16,
    maxComponentElements: u16,
    maxComponentDepth: u16,
};

const HheaTable = packed struct {
    version: Fixed,
    ascent: i16,
    descent: i16,
    line_gap: i16,
    advance_width_max: u16,
    min_left_side_bearing: i16,
    min_right_side_bearing: i16,
    x_max_extent: i16,
    caret_slope_rise: i16,
    caret_slope_run: i16,
    caret_offset: i16,
    reserved1: i16,
    reserved2: i16,
    reserved3: i16,
    reserved4: i16,
    metric_data_format: i16,
    num_of_long_hor_metrics: u16,
};
