const types = @import("../stream_types.zig");
const bytes = @import("../../../bytes.zig");
const zigstd = @import("../../../zigstd.zig");

/// Schema for the OS/2 table across every version (0-5); each version only
/// appends fields, so one cumulative layout gives the right offset for any
/// version - see each accessor's doc comment for its minimum version.
/// ulUnicodeRange1-4/ulCodePageRange1-2 need align(2): the cumulative byte
/// count is only 2-aligned when they're reached, not the 4 a u32 wants.
const StreamLayout = extern struct {
    version: u16,
    x_avg_char_width: types.FWord,
    us_weight_class: u16,
    us_width_class: u16,
    fs_type: u16,
    y_subscript_x_size: types.FWord,
    y_subscript_y_size: types.FWord,
    y_subscript_x_offset: types.FWord,
    y_subscript_y_offset: types.FWord,
    y_superscript_x_size: types.FWord,
    y_superscript_y_size: types.FWord,
    y_superscript_x_offset: types.FWord,
    y_superscript_y_offset: types.FWord,
    y_strikeout_size: types.FWord,
    y_strikeout_position: types.FWord,
    s_family_class: i16,
    panose: [10]u8,
    ul_unicode_range_1: u32 align(2),
    ul_unicode_range_2: u32 align(2),
    ul_unicode_range_3: u32 align(2),
    ul_unicode_range_4: u32 align(2),
    ach_vend_id: types.TagStr,
    fs_selection: u16,
    us_first_char_index: u16,
    us_last_char_index: u16,
    s_typo_ascender: types.FWord,
    s_typo_descender: types.FWord,
    s_typo_line_gap: types.FWord,
    us_win_ascent: types.UFWord,
    us_win_descent: types.UFWord,
    // version 1+
    ul_code_page_range_1: u32 align(2),
    ul_code_page_range_2: u32 align(2),
    // version 2+
    sx_height: types.FWord,
    s_cap_height: types.FWord,
    us_default_char: u16,
    us_break_char: u16,
    us_max_context: u16,
    // version 5+
    us_lower_optical_point_size: u16,
    us_upper_optical_point_size: u16,
};

const v0_size = @offsetOf(StreamLayout, "ul_code_page_range_1");
const v1_size = @offsetOf(StreamLayout, "sx_height");
const v2_size = @offsetOf(StreamLayout, "us_lower_optical_point_size");
const v5_size = @sizeOf(StreamLayout);

/// fsSelection's bit flags. Bits 7-9 (use_typo_metrics/wws/oblique) are only
/// defined from version 4 onward; readers must ignore them on an earlier
/// table (per spec, fonts are required to leave them zero there anyway).
pub const FsSelection = packed struct(u16) {
    italic: bool,
    underscore: bool,
    negative: bool,
    outlined: bool,
    strikeout: bool,
    bold: bool,
    regular: bool,
    /// version 4+. If set, applications should prefer sTypoAscender/
    /// sTypoDescender/sTypoLineGap over hhea's Ascender/Descender/LineGap
    /// for default line spacing.
    use_typo_metrics: bool,
    /// version 4+.
    wws: bool,
    /// version 4+.
    oblique: bool,
    _reserved: u6 = 0,
};

/// OS/2 and Windows metrics. Some very old version-0 fonts shipped a
/// shortened table stopping at usLastCharIndex; not specially handled here
/// (init() requires the full 78-byte version-0 layout).
pub const Os2Table = struct {
    const Self = @This();

    table_data: []const u8,

    pub fn init(table_data: []const u8, check: bool) !Self {
        if (check) try bytes.checkOffset(table_data, @sizeOf(u16));
        const self: Self = .{ .table_data = table_data };

        if (check) {
            const min_size: usize = switch (self.version()) {
                0 => v0_size,
                1 => v1_size,
                2, 3, 4 => v2_size,
                5 => v5_size,
                // Unknown/future version: require at least the
                // always-present version-0 prefix.
                else => v0_size,
            };
            try bytes.checkOffset(table_data, min_size);
        }

        return self;
    }

    pub inline fn version(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "version")..]);
    }

    pub inline fn xAvgCharWidth(self: Self) types.FWord {
        return bytes.readIntBE(types.FWord, self.table_data[@offsetOf(StreamLayout, "x_avg_char_width")..]);
    }

    pub inline fn usWeightClass(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "us_weight_class")..]);
    }

    pub inline fn usWidthClass(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "us_width_class")..]);
    }

    pub inline fn fsType(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "fs_type")..]);
    }

    pub inline fn ySubscriptXSize(self: Self) types.FWord {
        return bytes.readIntBE(types.FWord, self.table_data[@offsetOf(StreamLayout, "y_subscript_x_size")..]);
    }

    pub inline fn ySubscriptYSize(self: Self) types.FWord {
        return bytes.readIntBE(types.FWord, self.table_data[@offsetOf(StreamLayout, "y_subscript_y_size")..]);
    }

    pub inline fn ySubscriptXOffset(self: Self) types.FWord {
        return bytes.readIntBE(types.FWord, self.table_data[@offsetOf(StreamLayout, "y_subscript_x_offset")..]);
    }

    pub inline fn ySubscriptYOffset(self: Self) types.FWord {
        return bytes.readIntBE(types.FWord, self.table_data[@offsetOf(StreamLayout, "y_subscript_y_offset")..]);
    }

    pub inline fn ySuperscriptXSize(self: Self) types.FWord {
        return bytes.readIntBE(types.FWord, self.table_data[@offsetOf(StreamLayout, "y_superscript_x_size")..]);
    }

    pub inline fn ySuperscriptYSize(self: Self) types.FWord {
        return bytes.readIntBE(types.FWord, self.table_data[@offsetOf(StreamLayout, "y_superscript_y_size")..]);
    }

    pub inline fn ySuperscriptXOffset(self: Self) types.FWord {
        return bytes.readIntBE(types.FWord, self.table_data[@offsetOf(StreamLayout, "y_superscript_x_offset")..]);
    }

    pub inline fn ySuperscriptYOffset(self: Self) types.FWord {
        return bytes.readIntBE(types.FWord, self.table_data[@offsetOf(StreamLayout, "y_superscript_y_offset")..]);
    }

    pub inline fn yStrikeoutSize(self: Self) types.FWord {
        return bytes.readIntBE(types.FWord, self.table_data[@offsetOf(StreamLayout, "y_strikeout_size")..]);
    }

    pub inline fn yStrikeoutPosition(self: Self) types.FWord {
        return bytes.readIntBE(types.FWord, self.table_data[@offsetOf(StreamLayout, "y_strikeout_position")..]);
    }

    pub inline fn sFamilyClass(self: Self) i16 {
        return bytes.getI16BE(self.table_data[@offsetOf(StreamLayout, "s_family_class")..]);
    }

    pub inline fn panose(self: Self) [10]u8 {
        const off = @offsetOf(StreamLayout, "panose");
        return self.table_data[off..][0..10].*;
    }

    pub inline fn ulUnicodeRange1(self: Self) u32 {
        return bytes.getU32BE(self.table_data[@offsetOf(StreamLayout, "ul_unicode_range_1")..]);
    }

    pub inline fn ulUnicodeRange2(self: Self) u32 {
        return bytes.getU32BE(self.table_data[@offsetOf(StreamLayout, "ul_unicode_range_2")..]);
    }

    pub inline fn ulUnicodeRange3(self: Self) u32 {
        return bytes.getU32BE(self.table_data[@offsetOf(StreamLayout, "ul_unicode_range_3")..]);
    }

    pub inline fn ulUnicodeRange4(self: Self) u32 {
        return bytes.getU32BE(self.table_data[@offsetOf(StreamLayout, "ul_unicode_range_4")..]);
    }

    pub inline fn achVendId(self: Self) types.TagStr {
        const off = @offsetOf(StreamLayout, "ach_vend_id");
        return .{ .data = self.table_data[off..][0..4].* };
    }

    pub inline fn fsSelection(self: Self) FsSelection {
        return @bitCast(bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "fs_selection")..]));
    }

    pub inline fn usFirstCharIndex(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "us_first_char_index")..]);
    }

    pub inline fn usLastCharIndex(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "us_last_char_index")..]);
    }

    pub inline fn sTypoAscender(self: Self) types.FWord {
        return bytes.readIntBE(types.FWord, self.table_data[@offsetOf(StreamLayout, "s_typo_ascender")..]);
    }

    pub inline fn sTypoDescender(self: Self) types.FWord {
        return bytes.readIntBE(types.FWord, self.table_data[@offsetOf(StreamLayout, "s_typo_descender")..]);
    }

    pub inline fn sTypoLineGap(self: Self) types.FWord {
        return bytes.readIntBE(types.FWord, self.table_data[@offsetOf(StreamLayout, "s_typo_line_gap")..]);
    }

    pub inline fn usWinAscent(self: Self) types.UFWord {
        return bytes.readIntBE(types.UFWord, self.table_data[@offsetOf(StreamLayout, "us_win_ascent")..]);
    }

    pub inline fn usWinDescent(self: Self) types.UFWord {
        return bytes.readIntBE(types.UFWord, self.table_data[@offsetOf(StreamLayout, "us_win_descent")..]);
    }

    /// version 1+.
    pub inline fn ulCodePageRange1(self: Self) u32 {
        return bytes.getU32BE(self.table_data[@offsetOf(StreamLayout, "ul_code_page_range_1")..]);
    }

    /// version 1+.
    pub inline fn ulCodePageRange2(self: Self) u32 {
        return bytes.getU32BE(self.table_data[@offsetOf(StreamLayout, "ul_code_page_range_2")..]);
    }

    /// version 2+.
    pub inline fn sxHeight(self: Self) types.FWord {
        return bytes.readIntBE(types.FWord, self.table_data[@offsetOf(StreamLayout, "sx_height")..]);
    }

    /// version 2+.
    pub inline fn sCapHeight(self: Self) types.FWord {
        return bytes.readIntBE(types.FWord, self.table_data[@offsetOf(StreamLayout, "s_cap_height")..]);
    }

    /// version 2+.
    pub inline fn usDefaultChar(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "us_default_char")..]);
    }

    /// version 2+.
    pub inline fn usBreakChar(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "us_break_char")..]);
    }

    /// version 2+.
    pub inline fn usMaxContext(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "us_max_context")..]);
    }

    /// version 5+.
    pub inline fn usLowerOpticalPointSize(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "us_lower_optical_point_size")..]);
    }

    /// version 5+.
    pub inline fn usUpperOpticalPointSize(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "us_upper_optical_point_size")..]);
    }
};

test "OS/2 StreamLayout matches the OpenType stream layout at every version boundary" {
    try zigstd.testExpectEqual(@as(usize, 78), v0_size);
    try zigstd.testExpectEqual(@as(usize, 86), v1_size);
    try zigstd.testExpectEqual(@as(usize, 96), v2_size);
    try zigstd.testExpectEqual(@as(usize, 100), v5_size);
    try zigstd.testExpectEqual(@as(usize, 32), @offsetOf(StreamLayout, "panose"));
    try zigstd.testExpectEqual(@as(usize, 42), @offsetOf(StreamLayout, "ul_unicode_range_1"));
    try zigstd.testExpectEqual(@as(usize, 58), @offsetOf(StreamLayout, "ach_vend_id"));
    try zigstd.testExpectEqual(@as(usize, 68), @offsetOf(StreamLayout, "s_typo_ascender"));
}

test "Os2Table.init accepts a version 0 table at exactly 78 bytes" {
    var buf: [v0_size]u8 = @splat(0);
    bytes.setU16BE(&buf, 0); // version
    bytes.setI16BE(buf[@offsetOf(StreamLayout, "s_typo_ascender")..], 1500);
    bytes.setU16BE(buf[@offsetOf(StreamLayout, "us_win_ascent")..], 2000);

    const os2 = try Os2Table.init(&buf, true);

    try zigstd.testExpectEqual(@as(u16, 0), os2.version());
    try zigstd.testExpectEqual(@as(i16, 1500), os2.sTypoAscender());
    try zigstd.testExpectEqual(@as(u16, 2000), os2.usWinAscent());
}

test "Os2Table.init rejects a version 0 table shorter than 78 bytes" {
    var buf: [40]u8 = @splat(0);
    bytes.setU16BE(&buf, 0); // version
    try zigstd.testExpectError(error.EndOfStream, Os2Table.init(&buf, true));
}

test "Os2Table.fsSelection decodes USE_TYPO_METRICS at bit 7" {
    var buf: [v0_size]u8 = @splat(0);
    bytes.setU16BE(&buf, 4); // version
    bytes.setU16BE(buf[@offsetOf(StreamLayout, "fs_selection")..], 0x00C1); // bits 0,6,7

    const os2 = try Os2Table.init(&buf, false);
    const fs = os2.fsSelection();

    try zigstd.testExpectTrue(fs.italic);
    try zigstd.testExpectTrue(fs.regular);
    try zigstd.testExpectTrue(fs.use_typo_metrics);
    try zigstd.testExpectTrue(!fs.bold);
}

test "Os2Table reads achVendId/panose as raw byte arrays" {
    var buf: [v0_size]u8 = @splat(0);
    bytes.setU16BE(&buf, 0); // version
    @memcpy(buf[@offsetOf(StreamLayout, "ach_vend_id")..][0..4], "RAND");
    @memcpy(buf[@offsetOf(StreamLayout, "panose")..][0..10], &[_]u8{ 2, 0, 6, 3, 0, 0, 0, 0, 0, 0 });

    const os2 = try Os2Table.init(&buf, true);

    try zigstd.testExpectEqualSlices(u8, "RAND", &os2.achVendId().data);
    try zigstd.testExpectEqualSlices(u8, &[_]u8{ 2, 0, 6, 3, 0, 0, 0, 0, 0, 0 }, &os2.panose());
}
