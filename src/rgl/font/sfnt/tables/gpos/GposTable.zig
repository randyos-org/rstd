const bytes = @import("../../../../../bytes.zig");
const types = @import("../../types.zig");
const pair_pos = @import("PairPos.zig");

const Offset16 = types.Offset16;
const PairPosFormat1 = pair_pos.PairPosFormat1;
const PairPosFormat2 = pair_pos.PairPosFormat2;
const GlyphIndex = types.GlyphIndex;
const zigstd = @import("../../../../../zigstd.zig");

/// The GPOS header. This is the version 1.0 layout; version 1.1 appends a
/// trailing featureVariationsOffset (Offset32) we don't need yet, so it's
/// omitted here rather than modeled.
///
/// GPOS is a large table family (shared Script/Feature/Lookup list machinery, 9
/// distinct lookup types covering mark attachment, cursive positioning,
/// contextual/chained-contextual positioning, etc.). Only LookupType 2 (Pair
/// Adjustment Positioning), alone or inside a LookupType 9 extension, is read:
/// it is what kerning uses (`glyphKernAdvance`).
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

/// A Lookup table's own header (subtableOffsets[] follows immediately after).
/// The optional trailing markFilteringSet field (present only when
/// lookupFlag's USE_MARK_FILTERING_SET bit is set) isn't modeled since it
/// isn't needed to read pair-adjustment subtables.
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
    if (self.majorVersion() != 1) return 0;
    if (self.minorVersion() > 1) return 0; // 1.1 only appends featureVariationsOffset

    for (0..self.numLookups()) |i| {
        const lookup = self.getLookup(i);
        const lookup_type = lookup.lookupType();
        if (lookup_type != 2 and lookup_type != 9) // Pair Adjustment Positioning Subtable or Positioning Subtable Extension Subtable
            continue;

        for (0..lookup.subtableCount()) |sti| {
            var subtable = lookup.getSubtableData(sti);
            if (lookup_type == 9) {
                const format = bytes.getU16BE(subtable[0..]);
                if (format != 1) continue;
                const extension_lookup_type = bytes.getU16BE(subtable[2..]);
                if (extension_lookup_type != 2) continue;
                subtable = subtable[bytes.getU32BE(subtable[4..])..];
            }
            if (pairKernAdvance(subtable, a, b)) |advance| return advance;
        }
    }

    return 0;
}

/// One PairPos subtable's kerning for (a, b); null when it doesn't apply and
/// the next subtable should be tried.
fn pairKernAdvance(subtable: []const u8, a: GlyphIndex, b: GlyphIndex) ?i16 {
    // both formats open with format, coverage and the two value formats
    const pair1: PairPosFormat1 = .{ .table_data = subtable };
    const pos_format = pair1.format();
    const coverage_index = pair1.getCoverage().index(a) orelse return null;

    // Kerning is the first glyph's xAdvance; PairSet and Class2Record sizes
    // come from both value formats.
    const value_format_1 = pair1.valueFormat1();
    const x_advance_off = value_format_1.xAdvanceOffset() orelse return null;

    switch (pos_format) {
        1 => {
            if (coverage_index >= pair1.pairSetCount()) return 0;
            const pair_set = pair1.getPairSet(coverage_index);
            const pair_value_count = pair_set.pairValueCount();
            if (pair_value_count == 0) return 0;

            const needle = @backingInt(b);
            var r: u32 = pair_value_count - 1;
            var l: u32 = 0;

            // Binary search.
            while (l <= r) {
                const m = (l + r) >> 1;
                const straw = @backingInt(pair_set.getSecondGlyph(m));
                if (needle < straw) {
                    if (m == 0) break;
                    r = m - 1;
                } else if (needle > straw) {
                    l = m + 1;
                } else {
                    return bytes.getI16BE(pair_set.table_data[pair_set.getPairValueRecordOffset(m) + 2 + x_advance_off ..]);
                }
            }
            return null;
        },
        2 => {
            const pair2: PairPosFormat2 = .{ .table_data = subtable };
            const glyph1class = pair2.getClassDef1().classOf(a);
            const glyph2class = pair2.getClassDef2().classOf(b);

            const class1_count: u32 = pair2.class1Count();
            const class2_count: u32 = pair2.class2Count();

            if (glyph1class >= class1_count) return 0; // malformed
            if (glyph2class >= class2_count) return 0; // malformed

            return bytes.getI16BE(subtable[pair2.getClass2RecordOffset(glyph1class, glyph2class) + x_advance_off ..]);
        },
        else => return null,
    }
}

test "GPOS StreamLayout/LookupHeader match the OpenType stream layout" {
    try zigstd.testExpectEqual(@as(usize, 10), @sizeOf(StreamLayout));
    try zigstd.testExpectEqual(@as(usize, 4), @offsetOf(StreamLayout, "script_list_offset"));
    try zigstd.testExpectEqual(@as(usize, 8), @offsetOf(StreamLayout, "lookup_list_offset"));

    try zigstd.testExpectEqual(@as(usize, 6), @sizeOf(LookupHeader));
}
