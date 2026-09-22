const types = @import("../stream_types.zig");
const bytes = @import("../../../../bytes.zig");
const zigstd = @import("../../../../zigstd.zig");

const Allocator = zigstd.Allocator;
const GlyphIndex = types.GlyphIndex;

// Horizontal metrics. Not self-describing: its record counts come from
// hhea.num_of_long_hor_metrics and maxp.num_glyphs, so init() takes both
// rather than deriving them from table_data the way other tables do. Those
// two counts are the only decoded state kept around; every metric value
// itself is read fresh from table_data on each getMetrics() call.
pub const HmtxTable = struct {
    const Self = @This();

    table_data: []const u8,
    num_h_metrics: u16,
    num_glyphs: u16,

    pub const HMetrics = extern struct {
        advance_width: types.UFWord,
        left_side_bearing: types.FWord,
    };

    pub fn init(table_data: []const u8, num_h_metrics: u16, num_glyphs: u16, check: bool) !Self {
        if (check) {
            try bytes.checkArray(table_data, 0, num_h_metrics, HMetrics);
            const lsb_start = @as(usize, num_h_metrics) * @sizeOf(HMetrics);
            const num_extra = if (num_glyphs > num_h_metrics) num_glyphs - num_h_metrics else 0;
            try bytes.checkArray(table_data, lsb_start, num_extra, types.FWord);
        }

        return .{ .table_data = table_data, .num_h_metrics = num_h_metrics, .num_glyphs = num_glyphs };
    }

    inline fn hMetricAt(self: Self, index: u16) HMetrics {
        const off = @as(usize, index) * @sizeOf(HMetrics);
        return .{
            .advance_width = bytes.readIntBE(types.UFWord, self.table_data[off..]),
            .left_side_bearing = bytes.readIntBE(types.FWord, self.table_data[off + 2 ..]),
        };
    }

    /// Glyphs past num_h_metrics share the last h_metrics entry's
    /// advance_width and carry only their own left_side_bearing.
    pub fn getMetrics(self: Self, glyph: GlyphIndex) HMetrics {
        if (self.num_h_metrics == 0) return .{ .advance_width = 0, .left_side_bearing = 0 };

        const glyph_index: u16 = @backingInt(glyph);
        if (glyph_index < self.num_h_metrics) return self.hMetricAt(glyph_index);

        const last = self.hMetricAt(self.num_h_metrics - 1);
        const lsb_index = glyph_index - self.num_h_metrics;
        const num_extra = if (self.num_glyphs > self.num_h_metrics) self.num_glyphs - self.num_h_metrics else 0;
        if (lsb_index >= num_extra) return last;

        const lsb_start = @as(usize, self.num_h_metrics) * @sizeOf(HMetrics);
        const lsb_off = lsb_start + @as(usize, lsb_index) * @sizeOf(types.FWord);
        return .{
            .advance_width = last.advance_width,
            .left_side_bearing = bytes.readIntBE(types.FWord, self.table_data[lsb_off..]),
        };
    }
};

test "HMetrics matches the stream layout (4 bytes, no padding)" {
    try zigstd.testExpectEqual(@as(usize, 4), @sizeOf(HmtxTable.HMetrics));
    try zigstd.testExpectEqual(@as(usize, 2), @offsetOf(HmtxTable.HMetrics, "left_side_bearing"));
}

fn buildHmtx(alloc: Allocator, h_metrics: []const HmtxTable.HMetrics, extra_lsbs: []const types.FWord) ![]u8 {
    const total = h_metrics.len * 4 + extra_lsbs.len * 2;
    const buf = try alloc.alloc(u8, total);
    @memset(buf, 0);

    var off: usize = 0;
    for (h_metrics) |m| {
        bytes.setU16BE(buf[off..], m.advance_width);
        bytes.setI16BE(buf[off + 2 ..], m.left_side_bearing);
        off += 4;
    }
    for (extra_lsbs) |lsb| {
        bytes.setI16BE(buf[off..], lsb);
        off += 2;
    }

    return buf;
}

test "HmtxTable.getMetrics reads full records for glyphs within num_h_metrics" {
    const alloc = zigstd.testAllocator;
    const h_metrics = [_]HmtxTable.HMetrics{
        .{ .advance_width = 500, .left_side_bearing = 10 },
        .{ .advance_width = 600, .left_side_bearing = -5 },
        .{ .advance_width = 700, .left_side_bearing = 0 },
    };
    const buf = try buildHmtx(alloc, &h_metrics, &[_]i16{ 1, 2 });
    defer alloc.free(buf);

    const hmtx = try HmtxTable.init(buf, 3, 5, true);

    const m0 = hmtx.getMetrics(@fromBackingInt(@intCast(0)));
    try zigstd.testExpectEqual(@as(u16, 500), m0.advance_width);
    try zigstd.testExpectEqual(@as(i16, 10), m0.left_side_bearing);

    const m2 = hmtx.getMetrics(@fromBackingInt(@intCast(2)));
    try zigstd.testExpectEqual(@as(u16, 700), m2.advance_width);
    try zigstd.testExpectEqual(@as(i16, 0), m2.left_side_bearing);
}

test "HmtxTable.getMetrics shares the last advance_width for glyphs past num_h_metrics" {
    const alloc = zigstd.testAllocator;
    const h_metrics = [_]HmtxTable.HMetrics{
        .{ .advance_width = 500, .left_side_bearing = 10 },
        .{ .advance_width = 700, .left_side_bearing = 0 },
    };
    const buf = try buildHmtx(alloc, &h_metrics, &[_]i16{ 3, -7 });
    defer alloc.free(buf);

    const hmtx = try HmtxTable.init(buf, 2, 4, true);

    const m2 = hmtx.getMetrics(@fromBackingInt(@intCast(2)));
    try zigstd.testExpectEqual(@as(u16, 700), m2.advance_width); // shared from the last h_metrics entry
    try zigstd.testExpectEqual(@as(i16, 3), m2.left_side_bearing); // its own trailing lsb

    const m3 = hmtx.getMetrics(@fromBackingInt(@intCast(3)));
    try zigstd.testExpectEqual(@as(u16, 700), m3.advance_width);
    try zigstd.testExpectEqual(@as(i16, -7), m3.left_side_bearing);
}

test "HmtxTable.init handles num_h_metrics == num_glyphs (no trailing lsb array)" {
    const alloc = zigstd.testAllocator;
    const h_metrics = [_]HmtxTable.HMetrics{
        .{ .advance_width = 500, .left_side_bearing = 10 },
        .{ .advance_width = 600, .left_side_bearing = -5 },
    };
    const buf = try buildHmtx(alloc, &h_metrics, &[_]i16{});
    defer alloc.free(buf);

    const hmtx = try HmtxTable.init(buf, 2, 2, true);

    const m1 = hmtx.getMetrics(@fromBackingInt(@intCast(1)));
    try zigstd.testExpectEqual(@as(u16, 600), m1.advance_width);
    try zigstd.testExpectEqual(@as(i16, -5), m1.left_side_bearing);
}

test "HmtxTable.getMetrics falls back to the last full record for an out-of-range glyph" {
    const alloc = zigstd.testAllocator;
    const h_metrics = [_]HmtxTable.HMetrics{
        .{ .advance_width = 500, .left_side_bearing = 10 },
    };
    const buf = try buildHmtx(alloc, &h_metrics, &[_]i16{});
    defer alloc.free(buf);

    const hmtx = try HmtxTable.init(buf, 1, 1, true);

    const m = hmtx.getMetrics(@fromBackingInt(@intCast(99)));
    try zigstd.testExpectEqual(@as(u16, 500), m.advance_width);
    try zigstd.testExpectEqual(@as(i16, 10), m.left_side_bearing);
}
