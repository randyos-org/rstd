const std = @import("std");

const types = @import("../types.zig");
const bytes = @import("../../../../bytes.zig");
const ivs = @import("ivs.zig");
const zigstd = @import("../../../../zigstd.zig");

/// The MVAR header. valueRecords[] follows, strided by valueRecordSize (not
/// @sizeOf(ValueRecord)); see getValueRecord().
const StreamLayout = extern struct {
    major_version: u16,
    minor_version: u16,
    reserved: u16,
    value_record_size: u16,
    value_record_count: u16,
    item_variation_store_offset: types.Offset16,
};

/// One ValueRecord. valueTag identifies a font-wide measure from another
/// table (e.g. 'hasc' = OS/2.sTypoAscender); not modeled as an enum since
/// callers typically compare against just a few tags.
pub const ValueRecord = extern struct {
    value_tag: types.TagStr,
    delta_set_outer_index: u16,
    delta_set_inner_index: u16,
};

/// Metrics variations: deltas for font-wide metrics from OS/2, hhea, gasp,
/// post and vhea, keyed to regions in the shared Item Variation Store (see
/// ivs.zig). Resolving an actual interpolated value is deferred.
pub const MvarTable = struct {
    const Self = @This();

    table_data: []const u8,

    pub fn init(table_data: []const u8, check: bool) !Self {
        if (check) try bytes.checkSize(table_data, 0, StreamLayout);
        const self: Self = .{ .table_data = table_data };
        if (check) {
            try bytes.checkSlice(table_data, @sizeOf(StreamLayout), @as(usize, self.valueRecordCount()) * self.valueRecordSize());
        }
        return self;
    }

    pub inline fn majorVersion(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "major_version")..]);
    }

    pub inline fn minorVersion(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "minor_version")..]);
    }

    /// Size in bytes of each ValueRecord in this font; use this (not
    /// @sizeOf(ValueRecord)) as the stride when walking valueRecords, per
    /// spec.
    pub inline fn valueRecordSize(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "value_record_size")..]);
    }

    pub inline fn valueRecordCount(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "value_record_count")..]);
    }

    pub fn getItemVariationStore(self: Self) ?ivs.ItemVariationStore {
        const off = bytes.readIntBE(types.Offset16, self.table_data[@offsetOf(StreamLayout, "item_variation_store_offset")..]);
        if (off == 0) return null;
        return .{ .table_data = self.table_data[off..] };
    }

    pub fn getValueRecord(self: Self, index: usize) ValueRecord {
        const off = @sizeOf(StreamLayout) + index * self.valueRecordSize();
        return bytes.typeBEToNative(ValueRecord, self.table_data[off..]);
    }

    /// Binary search for an exact valueTag, per spec's requirement that
    /// value records be sorted in binary tag order. Null means this target
    /// item is constant across the font's variation space (no delta data
    /// at all, not just a zero delta).
    pub fn findValueRecord(self: Self, tag: types.TagStr) ?ValueRecord {
        var lo: usize = 0;
        var hi: usize = self.valueRecordCount();
        while (lo < hi) {
            const mid = lo + (hi - lo) / 2;
            const record = self.getValueRecord(mid);
            const order = std.mem.order(u8, &record.value_tag.data, &tag.data);
            switch (order) {
                .lt => lo = mid + 1,
                .gt => hi = mid,
                .eq => return record,
            }
        }
        return null;
    }
};

test "MVAR StreamLayout/ValueRecord match the OpenType stream layout" {
    try zigstd.testExpectEqual(@as(usize, 12), @sizeOf(StreamLayout));
    try zigstd.testExpectEqual(@as(usize, 6), @offsetOf(StreamLayout, "value_record_size"));
    try zigstd.testExpectEqual(@as(usize, 10), @offsetOf(StreamLayout, "item_variation_store_offset"));
    try zigstd.testExpectEqual(@as(usize, 8), @sizeOf(ValueRecord));
    try zigstd.testExpectEqual(@as(usize, 4), @offsetOf(ValueRecord, "delta_set_outer_index"));
}

fn buildMvar(alloc: zigstd.Allocator, tags: []const []const u8) ![]u8 {
    const header_size = @sizeOf(StreamLayout);
    const record_size = @sizeOf(ValueRecord);
    const total = header_size + tags.len * record_size;

    const buf = try alloc.alloc(u8, total);
    @memset(buf, 0);

    bytes.setU16BE(buf, 1); // majorVersion
    bytes.setU16BE(buf[2..], 0); // minorVersion
    bytes.setU16BE(buf[4..], 0); // reserved
    bytes.setU16BE(buf[6..], record_size);
    bytes.setU16BE(buf[8..], @intCast(tags.len));
    bytes.setU16BE(buf[10..], 0); // itemVariationStoreOffset (NULL for this test)

    for (tags, 0..) |tag, i| {
        const off = header_size + i * record_size;
        @memcpy(buf[off..][0..4], tag);
        bytes.setU16BE(buf[off + 4 ..], @intCast(i)); // outer index (arbitrary, for identification)
        bytes.setU16BE(buf[off + 6 ..], @intCast(i + 100)); // inner index
    }

    return buf;
}

test "MvarTable.findValueRecord binary-searches tags in sorted order" {
    const alloc = zigstd.testAllocator;
    // Already in binary tag order: 'cpht' < 'hasc' < 'xhgt'
    const buf = try buildMvar(alloc, &.{ "cpht", "hasc", "xhgt" });
    defer alloc.free(buf);

    const mvar = try MvarTable.init(buf, true);
    try zigstd.testExpectEqual(@as(u16, 3), mvar.valueRecordCount());
    try zigstd.testExpectEqual(@as(?ivs.ItemVariationStore, null), mvar.getItemVariationStore());

    const hasc = mvar.findValueRecord(try types.TagStr.fromStr("hasc")).?;
    try zigstd.testExpectEqual(@as(u16, 1), hasc.delta_set_outer_index);
    try zigstd.testExpectEqual(@as(u16, 101), hasc.delta_set_inner_index);

    try zigstd.testExpectEqual(@as(?ValueRecord, null), mvar.findValueRecord(try types.TagStr.fromStr("wght")));
}
