const types = @import("../stream_types.zig");
const bytes = @import("../../../../bytes.zig");
const zigstd = @import("../../../../zigstd.zig");
const ivs = @import("ivs.zig");

/// Schema for the VVAR header, used only to compute field offsets/sizes at
/// comptime; never instantiated. Format is analogous to HVAR's (see
/// hvar.zig), with one addition: vOrgMappingOffset, used only by CFF2
/// outlines to vary the Y coordinate of a glyph's vertical origin.
const StreamLayout = extern struct {
    major_version: u16,
    minor_version: u16,
    item_variation_store_offset: types.Offset32 align(2),
    advance_height_mapping_offset: types.Offset32 align(2),
    tsb_mapping_offset: types.Offset32 align(2),
    bsb_mapping_offset: types.Offset32 align(2),
    v_org_mapping_offset: types.Offset32 align(2),
};

/// Vertical metrics variations: per-glyph deltas for vmtx's advance heights,
/// side bearings, and (CFF2 only) vertical origins, keyed to regions in the
/// shared Item Variation Store (see ivs.zig). Resolving an actual
/// interpolated value is deferred - this exposes the raw store/mapping data.
pub const VvarTable = struct {
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

    /// Maps a glyph index to a delta-set index for advance-height deltas.
    /// Null (a NULL offset per spec) means glyph indices are used directly
    /// as delta-set inner indices, with outer index 0.
    pub fn getAdvanceHeightMapping(self: Self) ?ivs.DeltaSetIndexMap {
        const off = bytes.readIntBE(types.Offset32, self.table_data[@offsetOf(StreamLayout, "advance_height_mapping_offset")..]);
        if (off == 0) return null;
        return .{ .table_data = self.table_data[off..] };
    }

    /// Maps a glyph index to a delta-set index for top-side-bearing deltas.
    /// Null if this font doesn't provide TSB variation data.
    pub fn getTsbMapping(self: Self) ?ivs.DeltaSetIndexMap {
        const off = bytes.readIntBE(types.Offset32, self.table_data[@offsetOf(StreamLayout, "tsb_mapping_offset")..]);
        if (off == 0) return null;
        return .{ .table_data = self.table_data[off..] };
    }

    /// Maps a glyph index to a delta-set index for bottom-side-bearing
    /// deltas. Null if this font doesn't provide BSB variation data.
    pub fn getBsbMapping(self: Self) ?ivs.DeltaSetIndexMap {
        const off = bytes.readIntBE(types.Offset32, self.table_data[@offsetOf(StreamLayout, "bsb_mapping_offset")..]);
        if (off == 0) return null;
        return .{ .table_data = self.table_data[off..] };
    }

    /// Maps a glyph index to a delta-set index for vertical-origin-Y deltas
    /// (CFF2 outlines only). Null if this font doesn't provide that data.
    pub fn getVOrgMapping(self: Self) ?ivs.DeltaSetIndexMap {
        const off = bytes.readIntBE(types.Offset32, self.table_data[@offsetOf(StreamLayout, "v_org_mapping_offset")..]);
        if (off == 0) return null;
        return .{ .table_data = self.table_data[off..] };
    }
};

test "VVAR StreamLayout matches the OpenType stream layout" {
    try zigstd.testExpectEqual(@as(usize, 24), @sizeOf(StreamLayout));
    try zigstd.testExpectEqual(@as(usize, 4), @offsetOf(StreamLayout, "item_variation_store_offset"));
    try zigstd.testExpectEqual(@as(usize, 8), @offsetOf(StreamLayout, "advance_height_mapping_offset"));
    try zigstd.testExpectEqual(@as(usize, 12), @offsetOf(StreamLayout, "tsb_mapping_offset"));
    try zigstd.testExpectEqual(@as(usize, 16), @offsetOf(StreamLayout, "bsb_mapping_offset"));
    try zigstd.testExpectEqual(@as(usize, 20), @offsetOf(StreamLayout, "v_org_mapping_offset"));
}

test "VvarTable resolves the item variation store and treats a NULL mapping offset as absent" {
    var buf: [@sizeOf(StreamLayout) + 8]u8 = @splat(0);
    bytes.setU16BE(&buf, 1); // majorVersion
    bytes.setU16BE(buf[2..], 0); // minorVersion
    bytes.setU32BE(buf[4..], @sizeOf(StreamLayout)); // itemVariationStoreOffset
    // advanceHeightMappingOffset/tsbMappingOffset/bsbMappingOffset/vOrgMappingOffset left 0 (NULL)

    const store_off = @sizeOf(StreamLayout);
    bytes.setU16BE(buf[store_off..], 1); // ItemVariationStore.format
    bytes.setU16BE(buf[store_off + 6 ..], 0); // itemVariationDataCount

    const vvar = try VvarTable.init(&buf, true);
    try zigstd.testExpectEqual(@as(u16, 1), vvar.majorVersion());
    try zigstd.testExpectEqual(@as(u16, 1), vvar.getItemVariationStore().format());
    try zigstd.testExpectEqual(@as(?ivs.DeltaSetIndexMap, null), vvar.getAdvanceHeightMapping());
    try zigstd.testExpectEqual(@as(?ivs.DeltaSetIndexMap, null), vvar.getTsbMapping());
    try zigstd.testExpectEqual(@as(?ivs.DeltaSetIndexMap, null), vvar.getBsbMapping());
    try zigstd.testExpectEqual(@as(?ivs.DeltaSetIndexMap, null), vvar.getVOrgMapping());
}
