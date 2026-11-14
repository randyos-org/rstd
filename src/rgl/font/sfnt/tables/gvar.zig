const types = @import("../types.zig");
const bytes = @import("../../../../bytes.zig");
const tv = @import("tuple_variation.zig");
const zigstd = @import("../../../../zigstd.zig");

/// The gvar header. All fields land naturally aligned (the two Offset32 fields
/// each fall on a 4-byte boundary already), unlike HVAR/ivs.
const StreamLayout = extern struct {
    major_version: u16,
    minor_version: u16,
    axis_count: u16,
    shared_tuple_count: u16,
    shared_tuples_offset: types.Offset32,
    glyph_count: u16,
    flags: u16,
    glyph_variation_data_array_offset: types.Offset32,
};

/// gvar's offset-width flag: bit 0 clear means glyphVariationDataOffsets[]
/// entries are uint16 (halved, like 'loca' short format); set means uint32.
pub const GvarFlags = packed struct(u16) {
    long_offsets: bool,
    _reserved: u15 = 0,
};

/// Glyph variations: per-glyph outline-point deltas, using the Tuple
/// Variation Store format (see tuple_variation.zig). Interpolating an
/// actual instance outline is deferred; this exposes each glyph's raw
/// GlyphVariationData.
pub const GvarTable = struct {
    const Self = @This();

    table_data: []const u8,

    pub fn init(table_data: []const u8, check: bool) !Self {
        if (check) try bytes.checkSize(table_data, 0, StreamLayout);
        const self: Self = .{ .table_data = table_data };

        if (check) {
            const count = @as(usize, self.glyphCount()) + 1;
            if (self.flags().long_offsets) {
                try bytes.checkArray(table_data, @sizeOf(StreamLayout), count, types.Offset32);
            } else {
                try bytes.checkArray(table_data, @sizeOf(StreamLayout), count, u16);
            }
            try bytes.checkOffset(table_data, self.sharedTuplesOffset());
            try bytes.checkOffset(table_data, self.glyphVariationDataArrayOffset());
        }

        return self;
    }

    pub inline fn majorVersion(self: Self) u16 {
        return types.field(StreamLayout, "major_version", self.table_data);
    }

    pub inline fn minorVersion(self: Self) u16 {
        return types.field(StreamLayout, "minor_version", self.table_data);
    }

    pub inline fn axisCount(self: Self) u16 {
        return types.field(StreamLayout, "axis_count", self.table_data);
    }

    pub inline fn sharedTupleCount(self: Self) u16 {
        return types.field(StreamLayout, "shared_tuple_count", self.table_data);
    }

    pub inline fn sharedTuplesOffset(self: Self) types.Offset32 {
        return types.field(StreamLayout, "shared_tuples_offset", self.table_data);
    }

    pub inline fn glyphCount(self: Self) u16 {
        return types.field(StreamLayout, "glyph_count", self.table_data);
    }

    pub inline fn flags(self: Self) GvarFlags {
        return @bitCast(types.field(StreamLayout, "flags", self.table_data));
    }

    pub inline fn glyphVariationDataArrayOffset(self: Self) types.Offset32 {
        return types.field(StreamLayout, "glyph_variation_data_array_offset", self.table_data);
    }

    /// One of the shared tuples usable by any glyph's variation data.
    pub fn getSharedTuple(self: Self, index: usize) tv.Tuple {
        const off = self.sharedTuplesOffset() + index * @as(usize, self.axisCount()) * @sizeOf(types.F2Dot14);
        return .{ .table_data = self.table_data[off..], .axis_count = self.axisCount() };
    }

    /// Always returns an actual byte offset (from glyphVariationDataArrayOffset()),
    /// regardless of whether the stream stores it directly (long format) or
    /// halved (short format); mirrors LocaTable.getOffset().
    pub inline fn getGlyphVariationDataOffset(self: Self, glyph_index: usize) u32 {
        const array_off = @sizeOf(StreamLayout);
        if (self.flags().long_offsets) {
            return bytes.getU32BE(self.table_data[array_off + glyph_index * @sizeOf(u32) ..]);
        }
        return @as(u32, bytes.getU16BE(self.table_data[array_off + glyph_index * @sizeOf(u16) ..])) * 2;
    }

    /// Null means this glyph has no variation data (the same outline is used
    /// across the entire variation space), per spec; a valid condition, not
    /// a malformed-font one.
    pub fn getGlyphVariationData(self: Self, glyph_index: usize) ?tv.GlyphVariationData {
        const start = self.getGlyphVariationDataOffset(glyph_index);
        const end = self.getGlyphVariationDataOffset(glyph_index + 1);
        if (end <= start) return null;

        const base = self.glyphVariationDataArrayOffset();
        return .{ .table_data = self.table_data[base + start .. base + end], .axis_count = self.axisCount() };
    }
};

test "gvar StreamLayout matches the OpenType stream layout" {
    try zigstd.testExpectEqual(@as(usize, 20), @sizeOf(StreamLayout));
    try zigstd.testExpectEqual(@as(usize, 8), @offsetOf(StreamLayout, "shared_tuples_offset"));
    try zigstd.testExpectEqual(@as(usize, 16), @offsetOf(StreamLayout, "glyph_variation_data_array_offset"));
}

fn buildGvar(alloc: zigstd.Allocator) ![]u8 {
    // 1 axis, 0 shared tuples, 2 glyphs, short offsets.
    // glyph 0: one tuple variation header, no embedded tuple, 2 bytes of
    // dummy serialized data. glyph 1: no variation data (zero-size span).
    const header_size = @sizeOf(StreamLayout);
    const gvd_header = 4; // tupleVariationCount + dataOffset
    const tvh_size = 4; // fixed_size, no embedded tuple
    const glyph0_size = gvd_header + tvh_size + 2;
    // offsets array: 3 entries (glyphCount + 1), short format, must be even
    // so the halved values stay whole; pad glyph0_size up to a multiple of 2.
    const padded0 = glyph0_size + (glyph0_size % 2);
    const offsets_size = 3 * @sizeOf(u16);
    const array_off = header_size + offsets_size;
    const total = array_off + padded0;

    const buf = try alloc.alloc(u8, total);
    @memset(buf, 0);

    bytes.setU16BE(buf, 1); // majorVersion
    bytes.setU16BE(buf[2..], 0); // minorVersion
    bytes.setU16BE(buf[4..], 1); // axisCount
    bytes.setU16BE(buf[6..], 0); // sharedTupleCount
    bytes.setU32BE(buf[8..], array_off); // sharedTuplesOffset (unused, points past header)
    bytes.setU16BE(buf[12..], 2); // glyphCount
    bytes.setU16BE(buf[14..], 0); // flags: short offsets
    bytes.setU32BE(buf[16..], array_off); // glyphVariationDataArrayOffset

    bytes.setU16BE(buf[header_size..], 0); // offsets[0] halved = 0
    bytes.setU16BE(buf[header_size + 2 ..], @intCast(padded0 / 2)); // offsets[1] halved
    bytes.setU16BE(buf[header_size + 4 ..], @intCast(padded0 / 2)); // offsets[2] halved (glyph 1 empty)

    const gvd_off = array_off;
    const tvc: tv.TupleVariationCount = .{ .count = 1, .shared_point_numbers = false };
    bytes.setU16BE(buf[gvd_off..], @bitCast(tvc));
    bytes.setU16BE(buf[gvd_off + 2 ..], gvd_header + tvh_size); // dataOffset

    const tvh_off = gvd_off + gvd_header;
    bytes.setU16BE(buf[tvh_off..], 2); // variationDataSize
    const idx: tv.TupleIndex = .{ .tuple_index = 0, .private_point_numbers = false, .intermediate_region = false, .embedded_peak_tuple = false };
    bytes.setU16BE(buf[tvh_off + 2 ..], @bitCast(idx));

    return buf;
}

test "GvarTable resolves per-glyph GlyphVariationData spans via the short offset array" {
    const alloc = zigstd.testAllocator;
    const buf = try buildGvar(alloc);
    defer alloc.free(buf);

    const gvar = try GvarTable.init(buf, true);
    try zigstd.testExpectEqual(@as(u16, 1), gvar.axisCount());
    try zigstd.testExpectEqual(@as(u16, 2), gvar.glyphCount());
    try zigstd.testExpectEqual(false, gvar.flags().long_offsets);

    const g0 = gvar.getGlyphVariationData(0).?;
    try zigstd.testExpectEqual(@as(u12, 1), g0.tupleVariationCount().count);
    const h0 = g0.getTupleVariationHeader(0);
    try zigstd.testExpectEqual(@as(u16, 2), h0.variationDataSize());

    try zigstd.testExpectEqual(@as(?tv.GlyphVariationData, null), gvar.getGlyphVariationData(1));
}
