const bytes = @import("../../../../../bytes.zig");

/// A ClassDef table (either format). Only structural array access is
/// provided here; resolving a specific glyph to its class (format 1: direct
/// offset from startGlyphID; format 2: binary search sorted ranges) is
/// deferred.
const ClassDef = @This();

table_data: []const u8,

/// Schema for a ClassDefFormat2 range entry; its extern layout doubles as
/// the wire schema for @sizeOf.
pub const ClassRangeRecord = extern struct {
    start_glyph_id: u16,
    end_glyph_id: u16,
    class_value: u16,
};

pub inline fn format(self: ClassDef) u16 {
    return bytes.getU16BE(self.table_data[0..]);
}

/// format 1 only.
pub inline fn startGlyphId(self: ClassDef) u16 {
    return bytes.getU16BE(self.table_data[@sizeOf(u16)..]);
}

/// format 1 only: number of entries in classValueArray.
pub inline fn classCount(self: ClassDef) u16 {
    return bytes.getU16BE(self.table_data[2 * @sizeOf(u16) ..]);
}

/// format 1 only: classValueArray[index].
pub fn getClassValue(self: ClassDef, index: usize) u16 {
    const pos = 3 * @sizeOf(u16) + index * @sizeOf(u16);
    return bytes.getU16BE(self.table_data[pos..]);
}

/// format 2 only: number of ClassRangeRecords.
pub inline fn classRangeCount(self: ClassDef) u16 {
    return bytes.getU16BE(self.table_data[@sizeOf(u16)..]);
}

/// format 2 only: classRangeRecords[index].
pub fn getClassRangeRecord(self: ClassDef, index: usize) ClassRangeRecord {
    const pos = 2 * @sizeOf(u16) + index * @sizeOf(ClassRangeRecord);
    return bytes.typeBEToNative(ClassRangeRecord, self.table_data[pos..]);
}

const zigstd = @import("../../../../../zigstd.zig");

test "ClassRangeRecord matches the OpenType stream layout" {
    try zigstd.testExpectEqual(@as(usize, 6), @sizeOf(ClassRangeRecord));
}
