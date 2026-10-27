const types = @import("../stream_types.zig");
const bytes = @import("../../../../bytes.zig");

/// Schema for the STAT header, used only to compute field offsets/sizes at
/// comptime; never instantiated. v1.0 (deprecated) omits
/// elided_fallback_name_id; v1.1 adds it; v1.2 only adds a new
/// axis-value-table format. offset_to_axis_value_offsets needs align(2).
const StreamLayout = extern struct {
    major_version: u16,
    minor_version: u16,
    design_axis_size: u16,
    design_axis_count: u16,
    design_axes_offset: types.Offset32,
    axis_value_count: u16,
    offset_to_axis_value_offsets: types.Offset32 align(2),
    elided_fallback_name_id: u16,
};

/// v1.0 (deprecated) size: header without elided_fallback_name_id.
const v1_0_size = @offsetOf(StreamLayout, "elided_fallback_name_id");

/// One entry in the design axes array. getDesignAxis() strides by
/// designAxisSize() (not @sizeOf(DesignAxisRecord)) - same pattern as
/// fvar's VariationAxisRecord.
pub const DesignAxisRecord = extern struct {
    axis_tag: types.TagStr,
    axis_name_id: u16,
    axis_ordering: u16,
};

/// Flags shared by all AxisValueTable formats.
pub const AxisValueFlags = packed struct(u16) {
    older_sibling_font_attribute: bool,
    elidable_axis_value_name: bool,
    _reserved: u14 = 0,
};

pub const AxisValueFormat1 = extern struct {
    format: u16,
    axis_index: u16,
    flags: u16,
    value_name_id: u16,
    value: types.Fixed,
};

pub const AxisValueFormat2 = extern struct {
    format: u16,
    axis_index: u16,
    flags: u16,
    value_name_id: u16,
    nominal_value: types.Fixed,
    range_min_value: types.Fixed,
    range_max_value: types.Fixed,
};

pub const AxisValueFormat3 = extern struct {
    format: u16,
    axis_index: u16,
    flags: u16,
    value_name_id: u16,
    value: types.Fixed,
    linked_value: types.Fixed,
};

/// Schema for format 4's fixed header, used only to compute field
/// offsets/sizes at comptime; never instantiated. The axisValues[] array
/// that follows is genuinely variable-length (sized by axisCount), so it's
/// walked via AxisValueFormat4.getAxisValue() rather than modeled as one
/// fixed extern struct.
const AxisValueFormat4Header = extern struct {
    format: u16,
    axis_count: u16,
    flags: u16,
    value_name_id: u16,
};

/// One (axisIndex, value) pair within an AxisValueFormat4's axisValues
/// array. value's natural 4-byte alignment doesn't match its 2-mod-4 offset
/// after axis_index, hence align(2) - same situation as StreamLayout above.
pub const AxisValueRecord = extern struct {
    axis_index: u16,
    value: types.Fixed align(2),
};

/// table_data starts at this format-4 table's own beginning.
pub const AxisValueFormat4 = struct {
    const Self = @This();

    table_data: []const u8,

    pub inline fn axisCount(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(AxisValueFormat4Header, "axis_count")..]);
    }

    pub inline fn flags(self: Self) AxisValueFlags {
        return @bitCast(bytes.getU16BE(self.table_data[@offsetOf(AxisValueFormat4Header, "flags")..]));
    }

    pub inline fn valueNameId(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(AxisValueFormat4Header, "value_name_id")..]);
    }

    pub fn getAxisValue(self: Self, index: usize) AxisValueRecord {
        const off = @sizeOf(AxisValueFormat4Header) + index * @sizeOf(AxisValueRecord);
        return bytes.typeBEToNative(AxisValueRecord, self.table_data[off..]);
    }
};

/// One axis value table, in whichever of the four formats its own `format`
/// field says. table_data starts at this axis value table's own beginning.
pub const AxisValueTable = struct {
    const Self = @This();

    table_data: []const u8,

    pub inline fn format(self: Self) u16 {
        return bytes.getU16BE(self.table_data);
    }

    pub fn asFormat1(self: Self) AxisValueFormat1 {
        return bytes.typeBEToNative(AxisValueFormat1, self.table_data);
    }

    pub fn asFormat2(self: Self) AxisValueFormat2 {
        return bytes.typeBEToNative(AxisValueFormat2, self.table_data);
    }

    pub fn asFormat3(self: Self) AxisValueFormat3 {
        return bytes.typeBEToNative(AxisValueFormat3, self.table_data);
    }

    pub fn asFormat4(self: Self) AxisValueFormat4 {
        return .{ .table_data = self.table_data };
    }
};

/// Style attributes: design-axis metadata and axis-value labels used to
/// compose subfamily names. Resolving the best-match axis value table(s)
/// or composing a name string is deferred.
pub const StatTable = struct {
    const Self = @This();

    table_data: []const u8,

    pub fn init(table_data: []const u8, check: bool) !Self {
        if (check) try bytes.checkSlice(table_data, 0, v1_0_size);
        const self: Self = .{ .table_data = table_data };
        if (check and self.minorVersion() >= 1) try bytes.checkSize(table_data, 0, StreamLayout);
        return self;
    }

    pub inline fn majorVersion(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "major_version")..]);
    }

    pub inline fn minorVersion(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "minor_version")..]);
    }

    /// Size in bytes of each DesignAxisRecord in this font; use this (not
    /// @sizeOf(DesignAxisRecord)) as the stride when walking design axes,
    /// per spec.
    pub inline fn designAxisSize(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "design_axis_size")..]);
    }

    pub inline fn designAxisCount(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "design_axis_count")..]);
    }

    inline fn designAxesOffset(self: Self) types.Offset32 {
        return bytes.readIntBE(types.Offset32, self.table_data[@offsetOf(StreamLayout, "design_axes_offset")..]);
    }

    pub inline fn axisValueCount(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "axis_value_count")..]);
    }

    inline fn offsetToAxisValueOffsets(self: Self) types.Offset32 {
        return bytes.readIntBE(types.Offset32, self.table_data[@offsetOf(StreamLayout, "offset_to_axis_value_offsets")..]);
    }

    /// Null in a version 1.0 table (deprecated per spec - use of version
    /// 1.0 is discouraged and this field is simply absent from the stream).
    pub fn elidedFallbackNameId(self: Self) ?u16 {
        if (self.minorVersion() < 1) return null;
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "elided_fallback_name_id")..]);
    }

    pub fn getDesignAxis(self: Self, index: usize) DesignAxisRecord {
        const off = self.designAxesOffset() + index * self.designAxisSize();
        return bytes.typeBEToNative(DesignAxisRecord, self.table_data[off..]);
    }

    /// Axis value table offsets are relative to the start of the axis value
    /// offsets array itself (offsetToAxisValueOffsets()), not to the start
    /// of the STAT table - same indirection style as gvar's
    /// glyphVariationDataOffsets.
    pub fn getAxisValue(self: Self, index: usize) AxisValueTable {
        const array_off = self.offsetToAxisValueOffsets();
        const entry_off = array_off + index * @sizeOf(u16);
        const rel = bytes.getU16BE(self.table_data[entry_off..]);
        return .{ .table_data = self.table_data[array_off + rel ..] };
    }
};

const std = @import("std");
const zigstd = @import("../../../../zigstd.zig");

test "STAT StreamLayout/records match the OpenType stream layout" {
    try zigstd.testExpectEqual(@as(usize, 20), @sizeOf(StreamLayout));
    try zigstd.testExpectEqual(@as(usize, 18), v1_0_size);
    try zigstd.testExpectEqual(@as(usize, 14), @offsetOf(StreamLayout, "offset_to_axis_value_offsets"));
    try zigstd.testExpectEqual(@as(usize, 8), @sizeOf(DesignAxisRecord));
    try zigstd.testExpectEqual(@as(usize, 12), @sizeOf(AxisValueFormat1));
    try zigstd.testExpectEqual(@as(usize, 20), @sizeOf(AxisValueFormat2));
    try zigstd.testExpectEqual(@as(usize, 16), @sizeOf(AxisValueFormat3));
    try zigstd.testExpectEqual(@as(usize, 6), @sizeOf(AxisValueRecord));
    try zigstd.testExpectEqual(@as(usize, 2), @offsetOf(AxisValueRecord, "value"));
}

fn buildStat(alloc: zigstd.Allocator) ![]u8 {
    // 1 design axis ('wght'), 2 axis value tables: format 1 (value 400,
    // "Regular") and format 3 (value 700, "Bold", linked from nowhere).
    const header_size = @sizeOf(StreamLayout);
    const axis_size = @sizeOf(DesignAxisRecord);
    const axes_off = header_size;
    const offsets_off = axes_off + axis_size;
    const offsets_size = 2 * @sizeOf(u16);
    const values_off = offsets_off + offsets_size;
    const av1_size = @sizeOf(AxisValueFormat1);
    const av3_size = @sizeOf(AxisValueFormat3);
    const total = values_off + av1_size + av3_size;

    const buf = try alloc.alloc(u8, total);
    @memset(buf, 0);

    bytes.setU16BE(buf, 1); // majorVersion
    bytes.setU16BE(buf[2..], 2); // minorVersion
    bytes.setU16BE(buf[4..], axis_size); // designAxisSize
    bytes.setU16BE(buf[6..], 1); // designAxisCount
    bytes.setU32BE(buf[8..], axes_off); // designAxesOffset
    bytes.setU16BE(buf[12..], 2); // axisValueCount
    bytes.setU32BE(buf[14..], offsets_off); // offsetToAxisValueOffsets
    bytes.setU16BE(buf[18..], 2); // elidedFallbackNameID

    @memcpy(buf[axes_off..][0..4], "wght");
    bytes.setU16BE(buf[axes_off + 4 ..], 256); // axisNameID
    bytes.setU16BE(buf[axes_off + 6 ..], 0); // axisOrdering

    bytes.setU16BE(buf[offsets_off..], @intCast(values_off - offsets_off)); // offset to format1 table
    bytes.setU16BE(buf[offsets_off + 2 ..], @intCast(values_off + av1_size - offsets_off)); // offset to format3 table

    const av1_off = values_off;
    bytes.setU16BE(buf[av1_off..], 1); // format
    bytes.setU16BE(buf[av1_off + 2 ..], 0); // axisIndex
    bytes.setU16BE(buf[av1_off + 4 ..], 0); // flags
    bytes.setU16BE(buf[av1_off + 6 ..], 257); // valueNameID
    bytes.setI32BE(buf[av1_off + 8 ..], 400 << 16); // value = 400.0 (Fixed)

    const av3_off = av1_off + av1_size;
    bytes.setU16BE(buf[av3_off..], 3); // format
    bytes.setU16BE(buf[av3_off + 2 ..], 0); // axisIndex
    bytes.setU16BE(buf[av3_off + 4 ..], 0); // flags
    bytes.setU16BE(buf[av3_off + 6 ..], 258); // valueNameID
    bytes.setI32BE(buf[av3_off + 8 ..], 700 << 16); // value = 700.0 (Fixed)
    bytes.setI32BE(buf[av3_off + 12 ..], 700 << 16); // linkedValue = 700.0

    return buf;
}

test "StatTable reads design axes and axis value tables by format" {
    const alloc = zigstd.testAllocator;
    const buf = try buildStat(alloc);
    defer alloc.free(buf);

    const stat = try StatTable.init(buf, true);
    try zigstd.testExpectEqual(@as(u16, 1), stat.designAxisCount());
    try zigstd.testExpectEqual(@as(u16, 2), stat.axisValueCount());
    try zigstd.testExpectEqual(@as(?u16, 2), stat.elidedFallbackNameId());

    const axis0 = stat.getDesignAxis(0);
    try zigstd.testExpectEqual(@as(u16, 256), axis0.axis_name_id);
    try zigstd.testExpectEqual(true, std.mem.eql(u8, &axis0.axis_tag.data, "wght"));

    const av0 = stat.getAxisValue(0);
    try zigstd.testExpectEqual(@as(u16, 1), av0.format());
    const f1 = av0.asFormat1();
    try zigstd.testExpectEqual(@as(i32, 400 << 16), @as(i32, @bitCast(f1.value)));

    const av1 = stat.getAxisValue(1);
    try zigstd.testExpectEqual(@as(u16, 3), av1.format());
    const f3 = av1.asFormat3();
    try zigstd.testExpectEqual(@as(i32, 700 << 16), @as(i32, @bitCast(f3.value)));
    try zigstd.testExpectEqual(@as(i32, 700 << 16), @as(i32, @bitCast(f3.linked_value)));
}
