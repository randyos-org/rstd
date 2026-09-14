const bytes = @import("../../../../bytes.zig");
const zigstd = @import("../../../../zigstd.zig");
const types = @import("../../stream_types.zig");
const Coverage = @import("Coverage.zig");
const ClassDef = @import("ClassDef.zig");

const GlyphIndex = types.GlyphIndex;
const Offset16 = types.Offset16;

/// A ValueRecord's shape is entirely determined by which of these bits are
/// set: one present field occupies exactly 2 bytes (an i16, or an Offset16
/// for the four *_device fields), in the fixed field order the names are
/// declared in here (first declared field is the low bit, per Zig's packed
/// struct convention, matching the spec's 0x0001..0x0080 bit assignment).
pub const ValueFormat = packed struct(u16) {
    x_placement: bool,
    y_placement: bool,
    x_advance: bool,
    y_advance: bool,
    x_placement_device: bool,
    y_placement_device: bool,
    x_advance_device: bool,
    y_advance_device: bool,
    _reserved: u8 = 0,

    /// Total byte width of a ValueRecord using this format.
    pub fn size(self: @This()) usize {
        return @as(usize, @popCount(@as(u16, @bitCast(self)))) * @sizeOf(i16);
    }

    /// Byte offset of xAdvance within a ValueRecord using this format, or
    /// null if this format has no xAdvance field at all. xPlacement and
    /// yPlacement (if present) always precede it.
    pub fn xAdvanceOffset(self: @This()) ?usize {
        if (!self.x_advance) return null;
        var off: usize = 0;
        if (self.x_placement) off += @sizeOf(i16);
        if (self.y_placement) off += @sizeOf(i16);
        return off;
    }
};

/// Schema for PairPosFormat1's own header (pairSetOffsets[] follows
/// immediately after); never instantiated.
const Format1Header = extern struct {
    format: u16,
    coverage_offset: Offset16,
    value_format1: u16,
    value_format2: u16,
    pair_set_count: u16,
};

/// PairPosFormat1: per-glyph-pair adjustments, one PairSet per covered first
/// glyph. table_data starts at this subtable's own beginning.
pub const PairPosFormat1 = struct {
    const Self = @This();

    table_data: []const u8,

    pub inline fn format(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(Format1Header, "format")..]);
    }

    pub inline fn coverageOffset(self: Self) Offset16 {
        return bytes.getU16BE(self.table_data[@offsetOf(Format1Header, "coverage_offset")..]);
    }

    pub inline fn valueFormat1(self: Self) ValueFormat {
        return @bitCast(bytes.getU16BE(self.table_data[@offsetOf(Format1Header, "value_format1")..]));
    }

    pub inline fn valueFormat2(self: Self) ValueFormat {
        return @bitCast(bytes.getU16BE(self.table_data[@offsetOf(Format1Header, "value_format2")..]));
    }

    pub inline fn pairSetCount(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(Format1Header, "pair_set_count")..]);
    }

    pub fn getCoverage(self: Self) Coverage {
        return .{ .table_data = self.table_data[self.coverageOffset()..] };
    }

    pub fn getPairSetOffset(self: Self, index: usize) usize {
        const pos = @sizeOf(Format1Header) + index * @sizeOf(Offset16);
        return bytes.getU16BE(self.table_data[pos..]);
    }

    pub fn getPairSet(self: Self, index: usize) PairSet {
        return .{
            .table_data = self.table_data[self.getPairSetOffset(index)..],
            .value_format1 = self.valueFormat1(),
            .value_format2 = self.valueFormat2(),
        };
    }
};

/// One PairSet: every PairValueRecord for a single covered first glyph.
/// Record size is variable (depends on the parent subtable's ValueFormats),
/// so those are carried alongside table_data rather than re-read per call.
pub const PairSet = struct {
    const Self = @This();

    table_data: []const u8,
    value_format1: ValueFormat,
    value_format2: ValueFormat,

    pub inline fn pairValueCount(self: Self) u16 {
        return bytes.getU16BE(self.table_data[0..]);
    }

    inline fn recordSize(self: Self) usize {
        return @sizeOf(u16) + self.value_format1.size() + self.value_format2.size();
    }

    /// Byte offset (from this PairSet's own start) of PairValueRecord
    /// `index`'s secondGlyph field; its ValueRecords (per value_format1/2)
    /// follow immediately after.
    pub fn getPairValueRecordOffset(self: Self, index: usize) usize {
        return @sizeOf(u16) + index * self.recordSize();
    }

    pub fn getSecondGlyph(self: Self, index: usize) GlyphIndex {
        return @enumFromInt(bytes.getU16BE(self.table_data[self.getPairValueRecordOffset(index)..]));
    }
};

/// Schema for PairPosFormat2's own header (Class1Records follow immediately
/// after); never instantiated.
const Format2Header = extern struct {
    format: u16,
    coverage_offset: Offset16,
    value_format1: u16,
    value_format2: u16,
    class_def1_offset: Offset16,
    class_def2_offset: Offset16,
    class1_count: u16,
    class2_count: u16,
};

/// PairPosFormat2: adjustments keyed by (class of first glyph, class of
/// second glyph) rather than by individual glyph pairs. table_data starts
/// at this subtable's own beginning.
pub const PairPosFormat2 = struct {
    const Self = @This();

    table_data: []const u8,

    pub inline fn format(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(Format2Header, "format")..]);
    }

    pub inline fn coverageOffset(self: Self) Offset16 {
        return bytes.getU16BE(self.table_data[@offsetOf(Format2Header, "coverage_offset")..]);
    }

    pub inline fn valueFormat1(self: Self) ValueFormat {
        return @bitCast(bytes.getU16BE(self.table_data[@offsetOf(Format2Header, "value_format1")..]));
    }

    pub inline fn valueFormat2(self: Self) ValueFormat {
        return @bitCast(bytes.getU16BE(self.table_data[@offsetOf(Format2Header, "value_format2")..]));
    }

    pub inline fn classDef1Offset(self: Self) Offset16 {
        return bytes.getU16BE(self.table_data[@offsetOf(Format2Header, "class_def1_offset")..]);
    }

    pub inline fn classDef2Offset(self: Self) Offset16 {
        return bytes.getU16BE(self.table_data[@offsetOf(Format2Header, "class_def2_offset")..]);
    }

    pub inline fn class1Count(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(Format2Header, "class1_count")..]);
    }

    pub inline fn class2Count(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(Format2Header, "class2_count")..]);
    }

    pub fn getCoverage(self: Self) Coverage {
        return .{ .table_data = self.table_data[self.coverageOffset()..] };
    }

    pub fn getClassDef1(self: Self) ClassDef {
        return .{ .table_data = self.table_data[self.classDef1Offset()..] };
    }

    pub fn getClassDef2(self: Self) ClassDef {
        return .{ .table_data = self.table_data[self.classDef2Offset()..] };
    }

    /// Byte offset (from this subtable's own start) of the Class2Record for
    /// (class1, class2). A Class1Record is class2Count consecutive
    /// Class2Records; there's no separate Class1Record offset table, so
    /// this is a direct 2D-array index.
    pub fn getClass2RecordOffset(self: Self, class1: usize, class2: usize) usize {
        const class2_record_size = self.valueFormat1().size() + self.valueFormat2().size();
        const class1_record_size = self.class2Count() * class2_record_size;
        return @sizeOf(Format2Header) + class1 * class1_record_size + class2 * class2_record_size;
    }
};

test "Format1Header/Format2Header match the OpenType stream layout" {
    try zigstd.testExpectEqual(@as(usize, 10), @sizeOf(Format1Header));
    try zigstd.testExpectEqual(@as(usize, 8), @offsetOf(Format1Header, "pair_set_count"));

    try zigstd.testExpectEqual(@as(usize, 16), @sizeOf(Format2Header));
    try zigstd.testExpectEqual(@as(usize, 12), @offsetOf(Format2Header, "class1_count"));
}

test "ValueFormat.size/xAdvanceOffset match the spec's field ordering" {
    const x_only: ValueFormat = .{ .x_placement = false, .y_placement = false, .x_advance = true, .y_advance = false, .x_placement_device = false, .y_placement_device = false, .x_advance_device = false, .y_advance_device = false };
    try zigstd.testExpectEqual(@as(usize, 2), x_only.size());
    try zigstd.testExpectEqual(@as(?usize, 0), x_only.xAdvanceOffset());

    const placement_and_advance: ValueFormat = .{ .x_placement = true, .y_placement = true, .x_advance = true, .y_advance = false, .x_placement_device = false, .y_placement_device = false, .x_advance_device = false, .y_advance_device = false };
    try zigstd.testExpectEqual(@as(usize, 6), placement_and_advance.size());
    try zigstd.testExpectEqual(@as(?usize, 4), placement_and_advance.xAdvanceOffset());

    const none: ValueFormat = .{ .x_placement = true, .y_placement = false, .x_advance = false, .y_advance = false, .x_placement_device = false, .y_placement_device = false, .x_advance_device = false, .y_advance_device = false };
    try zigstd.testExpectEqual(@as(?usize, null), none.xAdvanceOffset());
}
