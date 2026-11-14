const types = @import("../types.zig");
const bytes = @import("../../../../bytes.zig");
const zigstd = @import("../../../../zigstd.zig");

const StreamLayout = extern struct {
    version: types.Fixed,
    ascent: types.FWord,
    descent: types.FWord,
    line_gap: types.FWord,
    advance_width_max: types.UFWord,
    min_left_side_bearing: types.FWord,
    min_right_side_bearing: types.FWord,
    x_max_extent: types.FWord,
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

pub const HheaTable = struct {
    const Self = @This();

    table_data: []const u8,

    pub fn init(table_data: []const u8, check: bool) !Self {
        if (check) try bytes.checkSize(table_data, 0, StreamLayout);
        return .{ .table_data = table_data };
    }

    pub inline fn version(self: Self) types.Fixed {
        return types.field(StreamLayout, "version", self.table_data);
    }

    pub inline fn ascent(self: Self) types.FWord {
        return types.field(StreamLayout, "ascent", self.table_data);
    }

    pub inline fn descent(self: Self) types.FWord {
        return types.field(StreamLayout, "descent", self.table_data);
    }

    pub inline fn lineGap(self: Self) types.FWord {
        return types.field(StreamLayout, "line_gap", self.table_data);
    }

    pub inline fn advanceWidthMax(self: Self) types.UFWord {
        return types.field(StreamLayout, "advance_width_max", self.table_data);
    }

    pub inline fn minLeftSideBearing(self: Self) types.FWord {
        return types.field(StreamLayout, "min_left_side_bearing", self.table_data);
    }

    pub inline fn minRightSideBearing(self: Self) types.FWord {
        return types.field(StreamLayout, "min_right_side_bearing", self.table_data);
    }

    pub inline fn xMaxExtent(self: Self) types.FWord {
        return types.field(StreamLayout, "x_max_extent", self.table_data);
    }

    pub inline fn caretSlopeRise(self: Self) i16 {
        return types.field(StreamLayout, "caret_slope_rise", self.table_data);
    }

    pub inline fn caretSlopeRun(self: Self) i16 {
        return types.field(StreamLayout, "caret_slope_run", self.table_data);
    }

    pub inline fn caretOffset(self: Self) i16 {
        return types.field(StreamLayout, "caret_offset", self.table_data);
    }

    pub inline fn reserved1(self: Self) i16 {
        return types.field(StreamLayout, "reserved1", self.table_data);
    }

    pub inline fn reserved2(self: Self) i16 {
        return types.field(StreamLayout, "reserved2", self.table_data);
    }

    pub inline fn reserved3(self: Self) i16 {
        return types.field(StreamLayout, "reserved3", self.table_data);
    }

    pub inline fn reserved4(self: Self) i16 {
        return types.field(StreamLayout, "reserved4", self.table_data);
    }

    pub inline fn metricDataFormat(self: Self) i16 {
        return types.field(StreamLayout, "metric_data_format", self.table_data);
    }

    pub inline fn numOfLongHorMetrics(self: Self) u16 {
        return types.field(StreamLayout, "num_of_long_hor_metrics", self.table_data);
    }
};

test "hhea StreamLayout matches the OpenType stream layout (36 bytes, no padding)" {
    try zigstd.testExpectEqual(@as(usize, 36), @sizeOf(StreamLayout));
    try zigstd.testExpectEqual(@as(usize, 4), @offsetOf(StreamLayout, "ascent"));
    try zigstd.testExpectEqual(@as(usize, 34), @offsetOf(StreamLayout, "num_of_long_hor_metrics"));
}

test "HheaTable getters read a real big-endian 'hhea' table correctly" {
    var buf: [36]u8 = @splat(0);
    bytes.setI16BE(buf[4..], -50); // ascent
    bytes.setU16BE(buf[34..], 7); // num_of_long_hor_metrics

    const hhea = try HheaTable.init(buf[0..], true);

    try zigstd.testExpectEqual(@as(i16, -50), hhea.ascent());
    try zigstd.testExpectEqual(@as(u16, 7), hhea.numOfLongHorMetrics());
}
