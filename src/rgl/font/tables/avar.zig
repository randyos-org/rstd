const std = @import("std");

const types = @import("../stream_types.zig");
const bytes = @import("../../../bytes.zig");
const zigstd = @import("../../../zigstd.zig");

/// Schema for the avar header (version 1.0 - the only version this spec
/// text defines), used only to compute field offsets/sizes at comptime;
/// never instantiated.
const StreamLayout = extern struct {
    major_version: u16,
    minor_version: u16,
    reserved: u16,
    axis_count: u16,
};

/// One AxisValueMap entry. Its extern layout doubles as the wire schema
/// for @sizeOf.
pub const AxisValueMap = extern struct {
    from_coordinate: types.F2Dot14,
    to_coordinate: types.F2Dot14,
};

/// One axis's SegmentMap: a list of (from, to) normalized-coordinate
/// correspondence pairs used to remodel that axis's default normalization
/// curve. table_data starts at this SegmentMap's own beginning.
pub const SegmentMap = struct {
    const Self = @This();

    table_data: []const u8,

    pub inline fn positionMapCount(self: Self) u16 {
        return bytes.getU16BE(self.table_data[0..]);
    }

    pub fn getAxisValueMap(self: Self, index: usize) AxisValueMap {
        const off = @sizeOf(u16) + index * @sizeOf(AxisValueMap);
        return bytes.typeBEToNative(AxisValueMap, self.table_data[off..]);
    }

    inline fn byteSize(self: Self) usize {
        return @sizeOf(u16) + @as(usize, self.positionMapCount()) * @sizeOf(AxisValueMap);
    }
};

/// Axis variations: optional piecewise-linear remapping of each axis's
/// default coordinate normalization. Applying the remapping is deferred -
/// this exposes the raw correspondence pairs.
pub const AvarTable = struct {
    const Self = @This();

    table_data: []const u8,

    pub fn init(table_data: []const u8, check: bool) !Self {
        if (check) try bytes.checkSize(table_data, 0, StreamLayout);
        const self: Self = .{ .table_data = table_data };

        if (check) {
            var off: usize = @sizeOf(StreamLayout);
            var i: u16 = 0;
            while (i < self.axisCount()) : (i += 1) {
                try bytes.checkSize(table_data, off, u16); // positionMapCount
                const map: SegmentMap = .{ .table_data = table_data[off..] };
                off += map.byteSize();
            }
            try bytes.checkOffset(table_data, off);
        }

        return self;
    }

    pub inline fn majorVersion(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "major_version")..]);
    }

    pub inline fn minorVersion(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "minor_version")..]);
    }

    pub inline fn axisCount(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "axis_count")..]);
    }

    /// Segment maps are variable-length (each has its own positionMapCount),
    /// so finding one past index 0 means walking every earlier map, the
    /// same way KernTable.getSubtableOffset() does.
    pub fn getSegmentMapOffset(self: Self, index: usize) usize {
        var off: usize = @sizeOf(StreamLayout);
        var i: usize = 0;
        while (i < index) : (i += 1) {
            const map: SegmentMap = .{ .table_data = self.table_data[off..] };
            off += map.byteSize();
        }
        return off;
    }

    pub fn getSegmentMap(self: Self, index: usize) SegmentMap {
        return .{ .table_data = self.table_data[self.getSegmentMapOffset(index)..] };
    }
};

test "avar StreamLayout/AxisValueMap match the OpenType stream layout" {
    try zigstd.testExpectEqual(@as(usize, 8), @sizeOf(StreamLayout));
    try zigstd.testExpectEqual(@as(usize, 6), @offsetOf(StreamLayout, "axis_count"));
    try zigstd.testExpectEqual(@as(usize, 4), @sizeOf(AxisValueMap));
    try zigstd.testExpectEqual(@as(usize, 2), @offsetOf(AxisValueMap, "to_coordinate"));
}

fn f2dot14Bytes(value: f64) [2]u8 {
    var buf: [2]u8 = undefined;
    const packed_val = types.F2Dot14.fromFloat(value);
    bytes.setI16BE(&buf, @bitCast(packed_val));
    return buf;
}

// Builds an avar table with two axes: the first has a modified segment map
// (three correspondence pairs), the second has none (empty segment map).
fn buildAvar(alloc: zigstd.Allocator) ![]u8 {
    const header_size = 8;
    const map0_size = @sizeOf(u16) + 3 * @sizeOf(AxisValueMap);
    const map1_size = @sizeOf(u16);
    const total = header_size + map0_size + map1_size;

    const buf = try alloc.alloc(u8, total);
    @memset(buf, 0);

    bytes.setU16BE(buf, 1); // majorVersion
    bytes.setU16BE(buf[2..], 0); // minorVersion
    bytes.setU16BE(buf[4..], 0); // reserved
    bytes.setU16BE(buf[6..], 2); // axisCount

    var off: usize = header_size;
    bytes.setU16BE(buf[off..], 3); // positionMapCount for axis 0
    off += @sizeOf(u16);
    const pairs = [_][2]f64{ .{ -1.0, -1.0 }, .{ 0.0, 0.0 }, .{ 1.0, 0.6 } };
    for (pairs) |p| {
        @memcpy(buf[off..][0..2], &f2dot14Bytes(p[0]));
        @memcpy(buf[off + 2 ..][0..2], &f2dot14Bytes(p[1]));
        off += @sizeOf(AxisValueMap);
    }

    bytes.setU16BE(buf[off..], 0); // positionMapCount for axis 1 (no remapping)
    off += @sizeOf(u16);

    return buf;
}

test "AvarTable walks variable-length segment maps and reads their AxisValueMaps" {
    const alloc = zigstd.testAllocator;
    const buf = try buildAvar(alloc);
    defer alloc.free(buf);

    const avar = try AvarTable.init(buf, true);
    try zigstd.testExpectEqual(@as(u16, 2), avar.axisCount());

    const map0 = avar.getSegmentMap(0);
    try zigstd.testExpectEqual(@as(u16, 3), map0.positionMapCount());
    const last = map0.getAxisValueMap(2);
    try zigstd.testExpectEqual(@as(f32, 1.0), last.from_coordinate.toFloat(f32));
    try std.testing.expectApproxEqAbs(@as(f32, 0.6), last.to_coordinate.toFloat(f32), 0.001);

    const map1 = avar.getSegmentMap(1);
    try zigstd.testExpectEqual(@as(u16, 0), map1.positionMapCount());
}
