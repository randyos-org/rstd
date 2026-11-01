const bytes = @import("../../../../../bytes.zig");
const types = @import("../../stream_types.zig");
const ValueFormat = @import("PairPos.zig").ValueFormat;

const Offset16 = types.Offset16;
const Offset32 = types.Offset32;
const GlyphIndex = types.GlyphIndex;

/// Schema for the GPOS header, used only to compute field offsets/sizes at
/// comptime; never instantiated. This is the version 1.0 layout; version
/// 1.1 appends a trailing featureVariationsOffset (Offset32) we don't need
/// yet, so it's omitted here rather than modeled.
///
/// GPOS is a large table family (shared Script/Feature/Lookup list
/// machinery, 9 distinct lookup types covering mark attachment, cursive
/// positioning, contextual/chained-contextual positioning, etc.). Only
/// enough of it to describe LookupType 2 (Pair Adjustment Positioning) is
/// captured here, since that's the one used for kerning; resolving an
/// actual kerning value (walking lookups, dispatching by format, running
/// Coverage/ClassDef searches) is deferred - table-structure capture only.
const StreamLayout = extern struct {
    major_version: u16,
    minor_version: u16,
    script_list_offset: Offset16,
    feature_list_offset: Offset16,
    lookup_list_offset: Offset16,
};

const GposTable = @This();

table_data: []const u8,

pub fn init(table_data: []const u8, check: bool) !GposTable {
    if (check) try bytes.checkSize(table_data, 0, StreamLayout);
    return .{ .table_data = table_data };
}

pub inline fn majorVersion(self: GposTable) u16 {
    return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "major_version")..]);
}

pub inline fn minorVersion(self: GposTable) u16 {
    return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "minor_version")..]);
}

pub inline fn scriptListOffset(self: GposTable) Offset16 {
    return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "script_list_offset")..]);
}

pub inline fn featureListOffset(self: GposTable) Offset16 {
    return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "feature_list_offset")..]);
}

pub inline fn lookupListOffset(self: GposTable) Offset16 {
    return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "lookup_list_offset")..]);
}

pub fn numLookups(self: GposTable) u16 {
    return bytes.getU16BE(self.table_data[self.lookupListOffset()..]);
}

pub fn getLookupOffset(self: GposTable, index: usize) usize {
    const list_off = self.lookupListOffset();
    const pos = list_off + @sizeOf(u16) + index * @sizeOf(Offset16);
    return list_off + bytes.getU16BE(self.table_data[pos..]);
}

pub fn getLookup(self: GposTable, index: usize) Lookup {
    return .{ .table_data = self.table_data[self.getLookupOffset(index)..] };
}

/// Schema for a Lookup table's own header (subtableOffsets[] follows
/// immediately after); never instantiated. The optional trailing
/// markFilteringSet field (present only when lookupFlag's
/// USE_MARK_FILTERING_SET bit is set) isn't modeled since it isn't needed
/// to read pair-adjustment subtables.
const LookupHeader = extern struct {
    lookup_type: u16,
    lookup_flag: u16,
    subtable_count: u16,
};

/// One Lookup table. table_data starts at this Lookup's own beginning
/// (i.e. already resolved from the parent LookupList's offset).
pub const Lookup = struct {
    const Self = @This();

    table_data: []const u8,

    pub inline fn lookupType(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(LookupHeader, "lookup_type")..]);
    }

    pub inline fn lookupFlag(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(LookupHeader, "lookup_flag")..]);
    }

    pub inline fn subtableCount(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(LookupHeader, "subtable_count")..]);
    }

    pub fn getSubtableOffset(self: Self, index: usize) usize {
        const pos = @sizeOf(LookupHeader) + index * @sizeOf(Offset16);
        return bytes.getU16BE(self.table_data[pos..]);
    }

    pub fn getSubtableData(self: Self, index: usize) []const u8 {
        return self.table_data[self.getSubtableOffset(index)..];
    }
};

/// An additional amount to advance the horizontal coordinate between the two
/// provided glyphs, per LookupType 2 (Pair Adjustment Positioning) or
/// LookupType 9 (Extension Positioning) wrapping one. 0 if none applies.
pub fn glyphKernAdvance(self: GposTable, a: GlyphIndex, b: GlyphIndex) i16 {
    const table_data = self.table_data;
    if (self.majorVersion() != 1) return 0;
    if (self.minorVersion() > 1) return 0; // 1.1 only appends featureVariationsOffset

    const lookup_list: Offset32 = self.lookupListOffset();
    const lookup_count = self.numLookups();

    for (0..lookup_count) |i| {
        const lookup_offset: Offset32 = bytes.getU16BE(table_data[lookup_list + 2 + 2 * i ..]);
        const lookup_table = lookup_list + lookup_offset;

        const lookup_type = bytes.getU16BE(table_data[lookup_table..]);
        const sub_table_count = bytes.getU16BE(table_data[lookup_table + 4 ..]);
        const sub_table_offsets = lookup_table + 6;
        if (lookup_type != 2 and lookup_type != 9) // Pair Adjustment Positioning Subtable or Positioning Subtable Extension Subtable
            continue;

        for (0..sub_table_count) |sti| {
            const subtable_offset: Offset32 = bytes.getU16BE(table_data[sub_table_offsets + 2 * sti ..]);
            var table = lookup_table + subtable_offset;
            if (lookup_type == 9) {
                const format = bytes.getU16BE(table_data[table + 0 ..]);
                if (format != 1) continue;
                const extension_lookup_type = bytes.getU16BE(table_data[table + 2 ..]);
                if (extension_lookup_type != 2) continue;
                table += bytes.getU32BE(table_data[table + 4 ..]);
            }
            const pos_format = bytes.getU16BE(table_data[table..]);
            const coverage_offset: Offset32 = bytes.getU16BE(table_data[table + 2 ..]);
            const coverage_index = coverageIndex(table_data, table + coverage_offset, a) orelse continue;

            // Kerning is the first glyph's xAdvance; records size by both formats.
            const value_format_1: ValueFormat = @bitCast(bytes.getU16BE(table_data[table + 4 ..]));
            const value_format_2: ValueFormat = @bitCast(bytes.getU16BE(table_data[table + 6 ..]));
            const x_advance_off = value_format_1.xAdvanceOffset() orelse continue;
            const value_records_size = value_format_1.size() + value_format_2.size();

            switch (pos_format) {
                1 => {
                    const pair_set_count = bytes.getU16BE(table_data[table + 8 ..]);
                    if (coverage_index >= pair_set_count) return 0;
                    const pair_pos_offset: Offset32 = bytes.getU16BE(table_data[table + 10 + 2 * coverage_index ..]);
                    const pair_value_table = table + pair_pos_offset;
                    const pair_value_count = bytes.getU16BE(table_data[pair_value_table..]);
                    const pair_value_array = pair_value_table + 2;
                    if (pair_value_count == 0) return 0;

                    const needle = @backingInt(b);
                    var r: u32 = pair_value_count - 1;
                    var l: u32 = 0;

                    // Binary search.
                    while (l <= r) {
                        const m = (l + r) >> 1;
                        const pair_value = pair_value_array + (2 + value_records_size) * m;
                        const straw = bytes.getU16BE(table_data[pair_value..]);
                        if (needle < straw) {
                            if (m == 0) break;
                            r = m - 1;
                        } else if (needle > straw) {
                            l = m + 1;
                        } else {
                            return bytes.getI16BE(table_data[pair_value + 2 + x_advance_off ..]);
                        }
                    }
                },
                2 => {
                    const class_def10_offset: u32 = bytes.getU16BE(table_data[table + 8 ..]);
                    const class_def20_offset: u32 = bytes.getU16BE(table_data[table + 10 ..]);
                    const glyph1class = glyphClass(table_data, table + class_def10_offset, a);
                    const glyph2class = glyphClass(table_data, table + class_def20_offset, b);

                    const class1_count: u32 = bytes.getU16BE(table_data[table + 12 ..]);
                    const class2_count: u32 = bytes.getU16BE(table_data[table + 14 ..]);

                    if (glyph1class >= class1_count) return 0; // malformed
                    if (glyph2class >= class2_count) return 0; // malformed

                    const class2_record = table + 16 + value_records_size * (glyph1class * class2_count + glyph2class);
                    return bytes.getI16BE(table_data[class2_record + x_advance_off ..]);
                },
                else => continue,
            }
        }
    }

    return 0;
}

fn coverageIndex(table_data: []const u8, coverage_table: u32, glyph: GlyphIndex) ?u32 {
    const coverage_format = bytes.getU16BE(table_data[coverage_table..]);
    switch (coverage_format) {
        1 => {
            const glyph_count = bytes.getU16BE(table_data[coverage_table + 2 ..]);
            if (glyph_count == 0) return null;

            // Binary search.
            var l: u32 = 0;
            var r: u32 = glyph_count - 1;
            const needle = @backingInt(glyph);
            while (l <= r) {
                const glyph_array = coverage_table + 4;
                const m = (l + r) >> 1;
                const glyph_id = bytes.getU16BE(table_data[glyph_array + 2 * m ..]);
                const straw = glyph_id;
                if (needle < straw) {
                    if (m == 0) break;
                    r = m - 1;
                } else if (needle > straw) {
                    l = m + 1;
                } else {
                    return m;
                }
            }
        },
        2 => {
            const range_count = bytes.getU16BE(table_data[coverage_table + 2 ..]);
            if (range_count == 0) return null;
            const range_array = coverage_table + 4;

            // Binary search.
            var l: u32 = 0;
            var r: u32 = range_count - 1;
            const needle = @backingInt(glyph);
            while (l <= r) {
                const m = (l + r) >> 1;
                const range_record = range_array + 6 * m;
                const straw_start = bytes.getU16BE(table_data[range_record..]);
                const straw_end = bytes.getU16BE(table_data[range_record + 2 ..]);
                if (needle < straw_start) {
                    if (m == 0) break;
                    r = m - 1;
                } else if (needle > straw_end) {
                    l = m + 1;
                } else {
                    const start_coverage_index = bytes.getU16BE(table_data[range_record + 4 ..]);
                    return start_coverage_index + needle - straw_start;
                }
            }
        },
        else => {},
    }
    return null;
}

fn glyphClass(table_data: []const u8, class_def_table: u32, glyph: GlyphIndex) u32 {
    const glyph_int = @backingInt(glyph);
    const class_def_format = bytes.getU16BE(table_data[class_def_table..]);
    switch (class_def_format) {
        1 => {
            const start_glyph_id = bytes.getU16BE(table_data[class_def_table + 2 ..]);
            const glyph_count = bytes.getU16BE(table_data[class_def_table + 4 ..]);
            const class_def1_value_array = class_def_table + 6;

            if (glyph_int >= start_glyph_id and glyph_int < start_glyph_id + glyph_count)
                return bytes.getU16BE(table_data[class_def1_value_array + 2 * (glyph_int - start_glyph_id) ..]);
        },
        2 => {
            const class_range_count = bytes.getU16BE(table_data[class_def_table + 2 ..]);
            const class_range_records = class_def_table + 4;

            if (class_range_count == 0)
                return 0;

            // Binary search.
            var l: u32 = 0;
            var r: u32 = class_range_count - 1;
            while (l <= r) {
                const m = (l + r) >> 1;
                const class_range_record = class_range_records + 6 * m;
                const straw_start = bytes.getU16BE(table_data[class_range_record..]);
                const straw_end = bytes.getU16BE(table_data[class_range_record + 2 ..]);
                if (glyph_int < straw_start) {
                    if (m == 0) break;
                    r = m - 1;
                } else if (glyph_int > straw_end) {
                    l = m + 1;
                } else {
                    return bytes.getU16BE(table_data[class_range_record + 4 ..]);
                }
            }
        },
        else => {},
    }

    // "All glyphs not assigned to a class fall into class 0". (OpenType spec)
    return 0;
}

const zigstd = @import("../../../../../zigstd.zig");

test "GPOS StreamLayout/LookupHeader match the OpenType stream layout" {
    try zigstd.testExpectEqual(@as(usize, 10), @sizeOf(StreamLayout));
    try zigstd.testExpectEqual(@as(usize, 4), @offsetOf(StreamLayout, "script_list_offset"));
    try zigstd.testExpectEqual(@as(usize, 8), @offsetOf(StreamLayout, "lookup_list_offset"));

    try zigstd.testExpectEqual(@as(usize, 6), @sizeOf(LookupHeader));
}
