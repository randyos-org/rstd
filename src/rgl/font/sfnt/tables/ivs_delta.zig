//! An Item Variation Store's delta at normalized coordinates (OpenType
//! 1.9.1, "OpenType Font Variations Common Table Formats", "Item Variation
//! Store"): each region's scalar times the item's delta for it, summed.
//! Every read is bounds-checked; `ivs.zig`'s accessors are not.

const bytes = @import("../../../../bytes.zig");
const types = @import("../types.zig");
const zigstd = @import("../../../../zigstd.zig");
const ivs = @import("ivs.zig");
const tent = @import("tent.zig").tent;
const units = @import("tent.zig").units;

const F2Dot14 = types.F2Dot14;

pub const Error = error{
    /// `coords.len` is not the store's axis count.
    WrongAxisCount,
    /// An offset, index or row runs past its data.
    MalformedVariationData,
};

/// The (outer, inner) index of `item` (a glyph ID, for HVAR) through
/// `map`: without a map, outer 0 and inner `item`; past its end, its last
/// entry, as the spec says.
pub fn mapIndex(map: ?ivs.DeltaSetIndexMap, item: u32) Error!ivs.DeltaSetIndex {
    const m = map orelse {
        if (item > 0xffff) return error.MalformedVariationData;
        return .{ .outer = 0, .inner = @intCast(item) };
    };
    const d = m.table_data;
    if (d.len < 2) return error.MalformedVariationData;
    const start: usize = if (m.format() == 0) 4 else 6;
    if (d.len < start) return error.MalformedVariationData;
    const count = m.mapCount();
    if (count == 0) return error.MalformedVariationData;
    const entry: usize = @min(item, count - 1);
    const size = @as(usize, m.entryFormat().entrySize()) + 1;
    if (start + (entry + 1) * size > d.len) return error.MalformedVariationData;
    return m.getEntry(entry);
}

/// `index`'s delta at `coords`, unrounded; 0 for an item the store has no
/// data for.
pub fn itemDelta(store: ivs.ItemVariationStore, index: ivs.DeltaSetIndex, coords: []const F2Dot14) Error!f32 {
    const s = store.table_data;
    if (s.len < 8) return error.MalformedVariationData;
    if (index.outer >= store.itemVariationDataCount()) return error.MalformedVariationData;
    const pos = 8 + @as(usize, index.outer) * 4;
    if (pos + 4 > s.len) return error.MalformedVariationData;
    const data_off = bytes.getU32BE(s[pos..]);
    if (data_off == 0) return 0;

    const rl_off = store.variationRegionListOffset();
    if (rl_off + 4 > s.len) return error.MalformedVariationData;
    const regions = store.getVariationRegionList();
    if (regions.axisCount() != coords.len) return error.WrongAxisCount;
    const region_count = regions.regionCount();
    if (rl_off + 4 + @as(usize, region_count) * coords.len * 6 > s.len) return error.MalformedVariationData;

    if (data_off + 6 > s.len) return error.MalformedVariationData;
    const data: ivs.ItemVariationData = .{ .table_data = s[data_off..] };
    if (index.inner >= data.itemCount()) return 0;
    const cols = data.regionIndexCount();
    const wdc = data.wordDeltaCount();
    if (wdc.count > cols) return error.MalformedVariationData;
    const long: usize = if (wdc.long_words) 4 else 2;
    const row = @as(usize, wdc.count) * long + (cols - @as(usize, wdc.count)) * (long / 2);
    if (6 + @as(usize, cols) * 2 + @as(usize, data.itemCount()) * row > data.table_data.len) return error.MalformedVariationData;

    var sum: f32 = 0;
    for (0..cols) |col| {
        const r = data.getRegionIndex(col);
        if (r >= region_count) return error.MalformedVariationData;
        var scalar: f32 = 1;
        for (coords, 0..) |c, axis| {
            const a = regions.getRegionAxis(r, axis);
            scalar *= tent(units(c), units(a.start_coord), units(a.peak_coord), units(a.end_coord));
            if (scalar == 0) break;
        }
        if (scalar != 0) sum += scalar * @as(f32, @floatFromInt(data.getDelta(index.inner, col)));
    }
    return sum;
}

test "itemDelta sums each region's share of an item's deltas" {
    // one axis; regions 0..1..1 and -1..-1..0; one subtable, two items,
    // one word column and one byte column
    var buf: [64]u8 = undefined;
    @memset(&buf, 0);
    bytes.setU16BE(buf[0..], 1); // format
    bytes.setU32BE(buf[2..], 12); // regionListOffset
    bytes.setU16BE(buf[6..], 1); // itemVariationDataCount
    bytes.setU32BE(buf[8..], 28); // itemVariationDataOffsets[0]
    bytes.setU16BE(buf[12..], 1); // axisCount
    bytes.setU16BE(buf[14..], 2); // regionCount
    for ([_][3]i16{ .{ 0, 16384, 16384 }, .{ -16384, -16384, 0 } }, 0..) |reg, k| {
        for (reg, 0..) |v, j| bytes.setI16BE(buf[16 + k * 6 + j * 2 ..], v);
    }
    bytes.setU16BE(buf[28..], 2); // itemCount
    bytes.setU16BE(buf[30..], 1); // wordDeltaCount
    bytes.setU16BE(buf[32..], 2); // regionIndexCount
    bytes.setU16BE(buf[34..], 0);
    bytes.setU16BE(buf[36..], 1);
    // item 0: 300, -20; item 1: -100, 7
    bytes.setI16BE(buf[38..], 300);
    buf[40] = @bitCast(@as(i8, -20));
    bytes.setI16BE(buf[41..], -100);
    buf[43] = 7;
    const store: ivs.ItemVariationStore = .{ .table_data = &buf };

    const at = struct {
        fn f(st: ivs.ItemVariationStore, inner: u16, c: f64) !f32 {
            return itemDelta(st, .{ .outer = 0, .inner = inner }, &.{F2Dot14.fromFloat(c)});
        }
    }.f;
    try zigstd.testExpectEqual(@as(f32, 150), try at(store, 0, 0.5));
    try zigstd.testExpectEqual(@as(f32, -10), try at(store, 0, -0.5));
    try zigstd.testExpectEqual(@as(f32, 7), try at(store, 1, -1.0));
    try zigstd.testExpectEqual(@as(f32, 0), try at(store, 0, 0.0));
    try zigstd.testExpectEqual(@as(f32, 0), try at(store, 5, 1.0));
    try zigstd.testExpectError(error.WrongAxisCount, itemDelta(store, .{ .outer = 0, .inner = 0 }, &.{}));
    try zigstd.testExpectError(error.MalformedVariationData, itemDelta(store, .{ .outer = 1, .inner = 0 }, &.{F2Dot14.fromFloat(0.5)}));
}

test "mapIndex: no map is the item itself; a map clamps to its last entry" {
    try zigstd.testExpectEqual(ivs.DeltaSetIndex{ .outer = 0, .inner = 9 }, try mapIndex(null, 9));
    // format 0, 1-byte entries with 4 inner bits: 0x12 is outer 1, inner 2
    const map_bytes = [_]u8{ 0, 0x03, 0, 2, 0x12, 0x34 };
    const map: ivs.DeltaSetIndexMap = .{ .table_data = &map_bytes };
    try zigstd.testExpectEqual(ivs.DeltaSetIndex{ .outer = 1, .inner = 2 }, try mapIndex(map, 0));
    try zigstd.testExpectEqual(ivs.DeltaSetIndex{ .outer = 3, .inner = 4 }, try mapIndex(map, 7));
}
