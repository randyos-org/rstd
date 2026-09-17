const types = @import("../stream_types.zig");
const bytes = @import("../../../../bytes.zig");
const zigstd = @import("../../../../zigstd.zig");
const ivs = @import("ivs.zig");

/// Schema for the HVAR header, used only to compute field offsets/sizes at
/// comptime; never instantiated.
const StreamLayout = extern struct {
    major_version: u16,
    minor_version: u16,
    item_variation_store_offset: types.Offset32 align(2),
    advance_width_mapping_offset: types.Offset32 align(2),
    lsb_mapping_offset: types.Offset32 align(2),
    rsb_mapping_offset: types.Offset32 align(2),
};

/// Horizontal metrics variations: per-glyph deltas for hmtx's advance
/// widths and side bearings, keyed to regions in the shared Item Variation
/// Store (see ivs.zig). Resolving an actual interpolated delta is
/// deferred - this exposes the raw store/mapping data.
pub const HvarTable = struct {
    const Self = @This();

    table_data: []const u8,

    pub fn init(table_data: []const u8, check: bool) !Self {
        if (check) try bytes.checkSize(table_data, 0, StreamLayout);
        return .{ .table_data = table_data };
    }

    pub inline fn majorVersion(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "major_version")..]);
    }

    pub inline fn minorVersion(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "minor_version")..]);
    }

    inline fn itemVariationStoreOffset(self: Self) types.Offset32 {
        return bytes.readIntBE(types.Offset32, self.table_data[@offsetOf(StreamLayout, "item_variation_store_offset")..]);
    }

    pub fn getItemVariationStore(self: Self) ivs.ItemVariationStore {
        return .{ .table_data = self.table_data[self.itemVariationStoreOffset()..] };
    }

    /// Maps a glyph index to a delta-set index for advance-width deltas.
    /// Null (a NULL offset per spec) means glyph indices are used directly
    /// as delta-set inner indices, with outer index 0.
    pub fn getAdvanceWidthMapping(self: Self) ?ivs.DeltaSetIndexMap {
        const off = bytes.readIntBE(types.Offset32, self.table_data[@offsetOf(StreamLayout, "advance_width_mapping_offset")..]);
        if (off == 0) return null;
        return .{ .table_data = self.table_data[off..] };
    }

    /// Maps a glyph index to a delta-set index for left-side-bearing
    /// deltas. Null if this font doesn't provide LSB variation data.
    pub fn getLsbMapping(self: Self) ?ivs.DeltaSetIndexMap {
        const off = bytes.readIntBE(types.Offset32, self.table_data[@offsetOf(StreamLayout, "lsb_mapping_offset")..]);
        if (off == 0) return null;
        return .{ .table_data = self.table_data[off..] };
    }

    /// Maps a glyph index to a delta-set index for right-side-bearing
    /// deltas. Null if this font doesn't provide RSB variation data.
    pub fn getRsbMapping(self: Self) ?ivs.DeltaSetIndexMap {
        const off = bytes.readIntBE(types.Offset32, self.table_data[@offsetOf(StreamLayout, "rsb_mapping_offset")..]);
        if (off == 0) return null;
        return .{ .table_data = self.table_data[off..] };
    }
};

test "HVAR StreamLayout matches the OpenType stream layout" {
    try zigstd.testExpectEqual(@as(usize, 20), @sizeOf(StreamLayout));
    try zigstd.testExpectEqual(@as(usize, 4), @offsetOf(StreamLayout, "item_variation_store_offset"));
    try zigstd.testExpectEqual(@as(usize, 8), @offsetOf(StreamLayout, "advance_width_mapping_offset"));
    try zigstd.testExpectEqual(@as(usize, 12), @offsetOf(StreamLayout, "lsb_mapping_offset"));
    try zigstd.testExpectEqual(@as(usize, 16), @offsetOf(StreamLayout, "rsb_mapping_offset"));
}

test "HvarTable resolves the item variation store and treats a NULL mapping offset as absent" {
    var buf: [@sizeOf(StreamLayout) + 8]u8 = @splat(0);
    bytes.setU16BE(&buf, 1); // majorVersion
    bytes.setU16BE(buf[2..], 0); // minorVersion
    bytes.setU32BE(buf[4..], @sizeOf(StreamLayout)); // itemVariationStoreOffset
    // advanceWidthMappingOffset/lsbMappingOffset/rsbMappingOffset left 0 (NULL)

    const store_off = @sizeOf(StreamLayout);
    bytes.setU16BE(buf[store_off..], 1); // ItemVariationStore.format
    bytes.setU16BE(buf[store_off + 6 ..], 0); // itemVariationDataCount

    const hvar = try HvarTable.init(&buf, true);
    try zigstd.testExpectEqual(@as(u16, 1), hvar.majorVersion());
    try zigstd.testExpectEqual(@as(u16, 1), hvar.getItemVariationStore().format());
    try zigstd.testExpectEqual(@as(?ivs.DeltaSetIndexMap, null), hvar.getAdvanceWidthMapping());
    try zigstd.testExpectEqual(@as(?ivs.DeltaSetIndexMap, null), hvar.getLsbMapping());
    try zigstd.testExpectEqual(@as(?ivs.DeltaSetIndexMap, null), hvar.getRsbMapping());
}
