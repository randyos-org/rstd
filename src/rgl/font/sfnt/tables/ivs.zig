const types = @import("../stream_types.zig");
const bytes = @import("../../../../bytes.zig");

// Item Variation Store: the shared variation-data format used by MVAR,
// HVAR, VVAR, BASE, GDEF and COLR. Distinct from the Tuple Variation Store
// gvar/cvar use (see tuple_variation.zig) - this organizes deltas by
// *target item* (O(1) lookups) rather than by region, and every row's size
// is directly computable, so it's fully structural (no run-length packing).

/// Schema for the ItemVariationStore header, used only to compute field
/// offsets/sizes at comptime; never instantiated. itemVariationDataOffsets[]
/// follows immediately after.
const StoreHeader = extern struct {
    format: u16,
    variation_region_list_offset: types.Offset32 align(2),
    item_variation_data_count: u16,
};

/// table_data starts at this store's own beginning (i.e. already resolved
/// from the parent HVAR/MVAR table's itemVariationStoreOffset).
pub const ItemVariationStore = struct {
    const Self = @This();

    table_data: []const u8,

    pub inline fn format(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StoreHeader, "format")..]);
    }

    pub inline fn variationRegionListOffset(self: Self) types.Offset32 {
        return bytes.readIntBE(types.Offset32, self.table_data[@offsetOf(StoreHeader, "variation_region_list_offset")..]);
    }

    pub inline fn itemVariationDataCount(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StoreHeader, "item_variation_data_count")..]);
    }

    pub fn getVariationRegionList(self: Self) VariationRegionList {
        return .{ .table_data = self.table_data[self.variationRegionListOffset()..] };
    }

    /// Null if this index has no associated subtable (a valid, explicit
    /// "no variation data" marker per spec, not a malformed-font condition).
    pub fn getItemVariationData(self: Self, index: usize) ?ItemVariationData {
        const pos = @sizeOf(StoreHeader) + index * @sizeOf(types.Offset32);
        const off = bytes.readIntBE(types.Offset32, self.table_data[pos..]);
        if (off == 0) return null;
        return .{ .table_data = self.table_data[off..] };
    }
};

/// Schema for a RegionAxisCoordinates record; its extern layout doubles as
/// the wire schema for @sizeOf. F2Dot14 (packed struct(i16)) is
/// extern-compatible since i16 is a power-of-two width, unlike the u24
/// case elsewhere in this font module.
pub const RegionAxisCoordinates = extern struct {
    start_coord: types.F2Dot14,
    peak_coord: types.F2Dot14,
    end_coord: types.F2Dot14,
};

/// Schema for the VariationRegionList header, used only to compute field
/// offsets/sizes at comptime; never instantiated.
const RegionListHeader = extern struct {
    axis_count: u16,
    region_count: u16,
};

/// table_data starts at this list's own beginning.
pub const VariationRegionList = struct {
    const Self = @This();

    table_data: []const u8,

    pub inline fn axisCount(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(RegionListHeader, "axis_count")..]);
    }

    pub inline fn regionCount(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(RegionListHeader, "region_count")..]);
    }

    /// regionAxes[axis] for the region at `region`, in the order 'fvar'
    /// defines its axes.
    pub fn getRegionAxis(self: Self, region: usize, axis: usize) RegionAxisCoordinates {
        const off = @sizeOf(RegionListHeader) +
            (region * @as(usize, self.axisCount()) + axis) * @sizeOf(RegionAxisCoordinates);
        return bytes.typeBEToNative(RegionAxisCoordinates, self.table_data[off..]);
    }
};

/// wordDeltaCount's packed flag/count field.
pub const WordDeltaCount = packed struct(u16) {
    count: u15,
    /// If set, "word" deltas are int32 (and "short" deltas are int16)
    /// rather than the usual int16/int8 split. Per spec, only ever set in
    /// the COLR table, never in MVAR/HVAR.
    long_words: bool,
};

/// Schema for the ItemVariationData header, used only to compute field
/// offsets/sizes at comptime; never instantiated. regionIndexes[] follows
/// immediately after, then the deltaSets rows.
const DataHeader = extern struct {
    item_count: u16,
    word_delta_count: u16,
    region_index_count: u16,
};

/// One ItemVariationData subtable: deltas for some set of target items
/// (rows) over some subset of the store's variation regions (columns).
/// table_data starts at this subtable's own beginning.
pub const ItemVariationData = struct {
    const Self = @This();

    table_data: []const u8,

    pub inline fn itemCount(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(DataHeader, "item_count")..]);
    }

    pub inline fn wordDeltaCount(self: Self) WordDeltaCount {
        return @bitCast(bytes.getU16BE(self.table_data[@offsetOf(DataHeader, "word_delta_count")..]));
    }

    pub inline fn regionIndexCount(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(DataHeader, "region_index_count")..]);
    }

    /// Index into the store's VariationRegionList for column `index` of
    /// this subtable's delta rows.
    pub fn getRegionIndex(self: Self, index: usize) u16 {
        const off = @sizeOf(DataHeader) + index * @sizeOf(u16);
        return bytes.getU16BE(self.table_data[off..]);
    }

    inline fn deltaSetsStart(self: Self) usize {
        return @sizeOf(DataHeader) + @as(usize, self.regionIndexCount()) * @sizeOf(u16);
    }

    inline fn rowSize(self: Self) usize {
        const wdc = self.wordDeltaCount();
        const long_size: usize = if (wdc.long_words) 4 else 2;
        const short_size: usize = if (wdc.long_words) 2 else 1;
        const short_count = @as(usize, self.regionIndexCount()) - wdc.count;
        return @as(usize, wdc.count) * long_size + short_count * short_size;
    }

    /// The delta for item (row) `item_index`, at this subtable's column
    /// `region_col` (0..regionIndexCount() - not a raw VariationRegionList
    /// index; use getRegionIndex() to map a column to one).
    pub fn getDelta(self: Self, item_index: usize, region_col: usize) i32 {
        const wdc = self.wordDeltaCount();
        const long_size: usize = if (wdc.long_words) 4 else 2;
        const short_size: usize = if (wdc.long_words) 2 else 1;
        const row_off = self.deltaSetsStart() + item_index * self.rowSize();

        if (region_col < wdc.count) {
            const off = row_off + region_col * long_size;
            return if (wdc.long_words)
                bytes.getI32BE(self.table_data[off..])
            else
                bytes.readIntBE(i16, self.table_data[off..]);
        }

        const long_bytes = @as(usize, wdc.count) * long_size;
        const off = row_off + long_bytes + (region_col - wdc.count) * short_size;
        return if (wdc.long_words)
            bytes.readIntBE(i16, self.table_data[off..])
        else
            @as(i8, @bitCast(self.table_data[off]));
    }
};

/// Maps a target item to a (outer, inner) delta-set index pair into an
/// ItemVariationStore, used by HVAR/VVAR/COLR (MVAR uses its own
/// ValueRecord array instead). Formats 0/1 differ only in mapCount's width.
pub const DeltaSetIndexMap = struct {
    const Self = @This();

    table_data: []const u8,

    pub inline fn format(self: Self) u8 {
        return self.table_data[0];
    }

    pub inline fn entryFormat(self: Self) EntryFormat {
        return @bitCast(self.table_data[1]);
    }

    pub fn mapCount(self: Self) u32 {
        return if (self.format() == 0)
            bytes.getU16BE(self.table_data[2..])
        else
            bytes.getU32BE(self.table_data[2..]);
    }

    inline fn mapDataStart(self: Self) usize {
        return if (self.format() == 0) 4 else 6;
    }

    /// Returns (outer index, inner index) for logical entry `index`. Per
    /// spec, an index >= mapCount() should be clamped to mapCount()-1 by
    /// callers rather than treated as out of range.
    pub fn getEntry(self: Self, index: usize) DeltaSetIndex {
        const ef = self.entryFormat();
        const entry_size: usize = @as(usize, ef.entrySize()) + 1;
        const off = self.mapDataStart() + index * entry_size;

        var raw: u32 = 0;
        for (0..entry_size) |i| raw = (raw << 8) | self.table_data[off + i];

        const inner_bits: u5 = @as(u5, ef.inner_index_bit_count) + 1;
        return .{
            .outer = @intCast(raw >> inner_bits),
            .inner = @intCast(raw & ((@as(u32, 1) << inner_bits) - 1)),
        };
    }
};

pub const DeltaSetIndex = struct {
    outer: u16,
    inner: u16,
};

pub const EntryFormat = packed struct(u8) {
    inner_index_bit_count: u4,
    map_entry_size_minus_one: u2,
    _reserved: u2 = 0,

    pub inline fn entrySize(self: @This()) u2 {
        return self.map_entry_size_minus_one;
    }
};

const zigstd = @import("../../../../zigstd.zig");

test "ItemVariationStore/VariationRegionList/ItemVariationData headers match the stream layout" {
    try zigstd.testExpectEqual(@as(usize, 8), @sizeOf(StoreHeader));
    try zigstd.testExpectEqual(@as(usize, 2), @offsetOf(StoreHeader, "variation_region_list_offset"));
    try zigstd.testExpectEqual(@as(usize, 4), @sizeOf(RegionListHeader));
    try zigstd.testExpectEqual(@as(usize, 6), @sizeOf(RegionAxisCoordinates));
    try zigstd.testExpectEqual(@as(usize, 6), @sizeOf(DataHeader));
}

fn f2dot14(alloc: zigstd.Allocator, values: []const f64) ![]u8 {
    const buf = try alloc.alloc(u8, values.len * 2);
    for (values, 0..) |v, i| {
        const packed_val = types.F2Dot14.fromFloat(v);
        bytes.setI16BE(buf[i * 2 ..], @bitCast(packed_val));
    }
    return buf;
}

test "VariationRegionList.getRegionAxis reads a 2-axis, 1-region list" {
    const alloc = zigstd.testAllocator;
    const header_size = 4;
    const total = header_size + 1 * 2 * 6; // 1 region, 2 axes, 6 bytes each

    const buf = try alloc.alloc(u8, total);
    defer alloc.free(buf);
    @memset(buf, 0);

    bytes.setU16BE(buf, 2); // axisCount
    bytes.setU16BE(buf[2..], 1); // regionCount

    const coords = try f2dot14(alloc, &.{ 0.0, 1.0, 1.0, -1.0, -0.5, 0.0 });
    defer alloc.free(coords);
    @memcpy(buf[header_size..][0..coords.len], coords);

    const list: VariationRegionList = .{ .table_data = buf };
    const axis0 = list.getRegionAxis(0, 0);
    try zigstd.testExpectEqual(@as(f32, 0.0), axis0.start_coord.toFloat(f32));
    try zigstd.testExpectEqual(@as(f32, 1.0), axis0.peak_coord.toFloat(f32));
    try zigstd.testExpectEqual(@as(f32, 1.0), axis0.end_coord.toFloat(f32));

    const axis1 = list.getRegionAxis(0, 1);
    try zigstd.testExpectEqual(@as(f32, -1.0), axis1.start_coord.toFloat(f32));
    try zigstd.testExpectEqual(@as(f32, -0.5), axis1.peak_coord.toFloat(f32));
    try zigstd.testExpectEqual(@as(f32, 0.0), axis1.end_coord.toFloat(f32));
}

fn buildItemVariationData(
    alloc: zigstd.Allocator,
    region_indexes: []const u16,
    word_delta_count: u15,
    rows: []const []const i32,
) ![]u8 {
    const header_size = 6;
    const regions_size = region_indexes.len * 2;
    const row_size = @as(usize, word_delta_count) * 2 + (region_indexes.len - word_delta_count) * 1;
    const total = header_size + regions_size + rows.len * row_size;

    const buf = try alloc.alloc(u8, total);
    @memset(buf, 0);

    bytes.setU16BE(buf, @intCast(rows.len)); // itemCount
    const wdc: WordDeltaCount = .{ .count = word_delta_count, .long_words = false };
    bytes.setU16BE(buf[2..], @bitCast(wdc));
    bytes.setU16BE(buf[4..], @intCast(region_indexes.len));

    var off: usize = header_size;
    for (region_indexes) |r| {
        bytes.setU16BE(buf[off..], r);
        off += 2;
    }

    for (rows) |row| {
        for (row, 0..) |v, col| {
            if (col < word_delta_count) {
                bytes.setI16BE(buf[off..], @intCast(v));
                off += 2;
            } else {
                buf[off] = @bitCast(@as(i8, @intCast(v)));
                off += 1;
            }
        }
    }

    return buf;
}

test "ItemVariationData.getDelta reads mixed word/short deltas per row" {
    const alloc = zigstd.testAllocator;
    // 3 regions, first 1 is a 16-bit ("word") delta, remaining 2 are 8-bit.
    const buf = try buildItemVariationData(alloc, &.{ 0, 2, 3 }, 1, &.{
        &.{ 1000, -5, 10 },
        &.{ -2000, 0, -128 },
    });
    defer alloc.free(buf);

    const data: ItemVariationData = .{ .table_data = buf };
    try zigstd.testExpectEqual(@as(u16, 2), data.itemCount());
    try zigstd.testExpectEqual(@as(u16, 3), data.regionIndexCount());
    try zigstd.testExpectEqual(@as(u16, 2), data.getRegionIndex(1));

    try zigstd.testExpectEqual(@as(i32, 1000), data.getDelta(0, 0));
    try zigstd.testExpectEqual(@as(i32, -5), data.getDelta(0, 1));
    try zigstd.testExpectEqual(@as(i32, 10), data.getDelta(0, 2));
    try zigstd.testExpectEqual(@as(i32, -2000), data.getDelta(1, 0));
    try zigstd.testExpectEqual(@as(i32, 0), data.getDelta(1, 1));
    try zigstd.testExpectEqual(@as(i32, -128), data.getDelta(1, 2));
}

test "DeltaSetIndexMap format 0 decodes packed outer/inner indices" {
    // entryFormat: innerIndexBitCount=8-1=7 (8 bits), entrySize=2-1=1 (2 bytes/entry)
    const ef: EntryFormat = .{ .inner_index_bit_count = 7, .map_entry_size_minus_one = 1 };
    var buf: [4 + 2 * 2]u8 = @splat(0);
    buf[0] = 0; // format
    buf[1] = @bitCast(ef);
    bytes.setU16BE(buf[2..], 2); // mapCount

    // entry 0: outer=1, inner=5 -> (1<<8)|5 = 0x0105
    bytes.setU16BE(buf[4..], 0x0105);
    // entry 1: outer=3, inner=200 -> (3<<8)|200 = 0x03C8
    bytes.setU16BE(buf[6..], 0x03C8);

    const map: DeltaSetIndexMap = .{ .table_data = &buf };
    try zigstd.testExpectEqual(@as(u32, 2), map.mapCount());

    const e0 = map.getEntry(0);
    try zigstd.testExpectEqual(@as(u16, 1), e0.outer);
    try zigstd.testExpectEqual(@as(u16, 5), e0.inner);

    const e1 = map.getEntry(1);
    try zigstd.testExpectEqual(@as(u16, 3), e1.outer);
    try zigstd.testExpectEqual(@as(u16, 200), e1.inner);
}
