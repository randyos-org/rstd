const types = @import("../types.zig");
const bytes = @import("../../../bytes.zig");
const zigstd = @import("../../../zigstd.zig");

const Allocator = zigstd.Allocator;

pub const HheaTable = extern struct {
    const Self = @This();

    version: types.Fixed,
    ascent: i16,
    descent: i16,
    line_gap: i16,
    advance_width_max: u16,
    min_left_side_bearing: i16,
    min_right_side_bearing: i16,
    x_max_extent: i16,
    caret_slope_rise: i16,
    caret_slope_run: i16,
    caret_offset: i16,
    reserved1: i16,
    reserved2: i16,
    reserved3: i16,
    reserved4: i16,
    metric_data_format: i16,
    num_of_long_hor_metrics: u16,

    pub fn create(alloc: Allocator, table_data: []const u8, check: bool) !*Self {
        if (check) try bytes.checkSize(table_data, 0, Self);
        return bytes.allocTypeBEToNative(alloc, Self, table_data);
    }

    pub fn destroy(self: *Self, alloc: Allocator) void {
        alloc.destroy(self);
    }
};

test "HheaTable matches the OpenType wire layout (36 bytes, no padding)" {
    try zigstd.testExpectEqual(@as(usize, 36), @sizeOf(HheaTable));
    try zigstd.testExpectEqual(@as(usize, 4), @offsetOf(HheaTable, "ascent"));
    try zigstd.testExpectEqual(@as(usize, 34), @offsetOf(HheaTable, "num_of_long_hor_metrics"));
}

test "HheaTable.create reads a real big-endian 'hhea' table correctly" {
    const alloc = zigstd.testAllocator;
    var buf: [36]u8 = @splat(0);
    bytes.writeIntBE(i16, buf[4..6], -50); // ascent
    bytes.writeIntBE(u16, buf[34..36], 7); // num_of_long_hor_metrics

    const hhea = try HheaTable.create(alloc, buf[0..], true);
    defer hhea.destroy(alloc);

    try zigstd.testExpectEqual(@as(i16, -50), hhea.ascent);
    try zigstd.testExpectEqual(@as(u16, 7), hhea.num_of_long_hor_metrics);
}
