const bytes = @import("../../../../../bytes.zig");
const zigstd = @import("../../../../../zigstd.zig");
const types = @import("../../stream_types.zig");

const GlyphIndex = types.GlyphIndex;

/// A Coverage table (either format). Only structural array access is
/// provided here; resolving a specific glyph to its coverage index (format
/// 1: binary search a sorted glyph list; format 2: binary search sorted
/// ranges) is deferred.
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

/// format 1 only: glyphArray[index].
pub fn getGlyph(self: Coverage, index: usize) GlyphIndex {
    const pos = 2 * @sizeOf(u16) + index * @sizeOf(u16);
    return @enumFromInt(bytes.getU16BE(self.table_data[pos..]));
}

/// format 2 only: number of RangeRecords.
pub inline fn rangeCount(self: Coverage) u16 {
    return bytes.getU16BE(self.table_data[@sizeOf(u16)..]);
}

/// format 2 only: rangeRecords[index].
pub fn getRangeRecord(self: Coverage, index: usize) RangeRecord {
    const pos = 2 * @sizeOf(u16) + index * @sizeOf(RangeRecord);
    return bytes.typeBEToNative(RangeRecord, self.table_data[pos..]);
}

test "RangeRecord matches the OpenType stream layout" {
    try zigstd.testExpectEqual(@as(usize, 6), @sizeOf(RangeRecord));
}
