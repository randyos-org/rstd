const bytes = @import("../../../../../bytes.zig");
const types = @import("../../types.zig");

const GlyphIndex = types.GlyphIndex;
const zigstd = @import("../../../../../zigstd.zig");

/// A Coverage table (either format): a sorted glyph list (format 1) or sorted
/// glyph ranges (format 2); `index` maps a glyph to its coverage index.
const Coverage = @This();

table_data: []const u8,

/// Schema for a CoverageFormat2 range entry; its extern layout doubles as
/// the wire schema for @sizeOf.
pub const RangeRecord = extern struct {
    start_glyph_id: u16,
    end_glyph_id: u16,
    start_coverage_index: u16,
};

pub inline fn format(self: Coverage) u16 {
    return bytes.getU16BE(self.table_data[0..]);
}

/// format 1 only: number of entries in the glyph array.
pub inline fn glyphCount(self: Coverage) u16 {
    return bytes.getU16BE(self.table_data[@sizeOf(u16)..]);
}

/// format 1 only: glyphArray[i].
pub fn getGlyph(self: Coverage, i: usize) GlyphIndex {
    const pos = 2 * @sizeOf(u16) + i * @sizeOf(u16);
    return @fromBackingInt(@intCast(bytes.getU16BE(self.table_data[pos..])));
}

/// format 2 only: number of RangeRecords.
pub inline fn rangeCount(self: Coverage) u16 {
    return bytes.getU16BE(self.table_data[@sizeOf(u16)..]);
}

/// format 2 only: rangeRecords[i].
pub fn getRangeRecord(self: Coverage, i: usize) RangeRecord {
    const pos = 2 * @sizeOf(u16) + i * @sizeOf(RangeRecord);
    return bytes.typeBEToNative(RangeRecord, self.table_data[pos..]);
}

/// `glyph`'s coverage index, or null if it isn't covered.
pub fn index(self: Coverage, glyph: GlyphIndex) ?u32 {
    const needle = @backingInt(glyph);
    switch (self.format()) {
        1 => {
            const glyph_count = self.glyphCount();
            if (glyph_count == 0) return null;

            // Binary search.
            var l: u32 = 0;
            var r: u32 = glyph_count - 1;
            while (l <= r) {
                const m = (l + r) >> 1;
                const straw = @backingInt(self.getGlyph(m));
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
            const range_count = self.rangeCount();
            if (range_count == 0) return null;

            // Binary search.
            var l: u32 = 0;
            var r: u32 = range_count - 1;
            while (l <= r) {
                const m = (l + r) >> 1;
                const range = self.getRangeRecord(m);
                if (needle < range.start_glyph_id) {
                    if (m == 0) break;
                    r = m - 1;
                } else if (needle > range.end_glyph_id) {
                    l = m + 1;
                } else {
                    return range.start_coverage_index + needle - range.start_glyph_id;
                }
            }
        },
        else => {},
    }
    return null;
}

test "RangeRecord matches the OpenType stream layout" {
    try zigstd.testExpectEqual(@as(usize, 6), @sizeOf(RangeRecord));
}
