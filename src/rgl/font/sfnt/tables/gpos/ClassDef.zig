const bytes = @import("../../../../../bytes.zig");
const types = @import("../../types.zig");

const GlyphIndex = types.GlyphIndex;
const zigstd = @import("../../../../../zigstd.zig");

/// A ClassDef table (either format): class values from startGlyphID on
/// (format 1) or sorted glyph ranges (format 2); `classOf` maps a glyph to its
/// class.
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

/// `glyph`'s class: "All glyphs not assigned to a class fall into class 0"
/// (OpenType spec).
pub fn classOf(self: ClassDef, glyph: GlyphIndex) u32 {
    const glyph_int = @backingInt(glyph);
    switch (self.format()) {
        1 => {
            const start_glyph_id = self.startGlyphId();
            const glyph_count = self.classCount();
            if (glyph_int >= start_glyph_id and glyph_int < start_glyph_id + glyph_count)
                return bytes.getU16BE(self.table_data[@as(u32, 3 * @sizeOf(u16)) + 2 * (glyph_int - start_glyph_id) ..]);
        },
        2 => {
            const class_range_count = self.classRangeCount();
            if (class_range_count == 0)
                return 0;

            // Binary search.
            var l: u32 = 0;
            var r: u32 = class_range_count - 1;
            while (l <= r) {
                const m = (l + r) >> 1;
                const range = self.getClassRangeRecord(m);
                if (glyph_int < range.start_glyph_id) {
                    if (m == 0) break;
                    r = m - 1;
                } else if (glyph_int > range.end_glyph_id) {
                    l = m + 1;
                } else {
                    return range.class_value;
                }
            }
        },
        else => {},
    }
    return 0;
}

test "ClassRangeRecord matches the OpenType stream layout" {
    try zigstd.testExpectEqual(@as(usize, 6), @sizeOf(ClassRangeRecord));
}
