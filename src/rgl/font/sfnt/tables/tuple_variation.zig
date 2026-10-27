const types = @import("../stream_types.zig");
const bytes = @import("../../../../bytes.zig");

// Tuple Variation Store: the shared format used by 'gvar' (per-glyph) and
// 'cvar' (whole-table), grouping deltas by *region* - distinct from the
// Item Variation Store format MVAR/HVAR/VVAR use (see ivs.zig), which
// groups by target item instead. The header array here is structural; the
// packed point-number/delta data is not - see serializedData().

/// The high 4 bits of a GlyphVariationData/'cvar' header's tupleVariationCount
/// field are flags (only SHARED_POINT_NUMBERS is defined); the low 12 bits
/// are the tuple variation table count.
pub const TupleVariationCount = packed struct(u16) {
    count: u12,
    _reserved: u3 = 0,
    /// If set, the serialized data begins with a shared set of packed point
    /// numbers usable by any of this glyph's tuple variation data tables.
    shared_point_numbers: bool,
};

/// TupleVariationHeader's tupleIndex packed field.
pub const TupleIndex = packed struct(u16) {
    /// Index into the shared tuples array (only meaningful when neither
    /// embedded_peak_tuple nor a private peak tuple applies).
    tuple_index: u12,
    _reserved: u1 = 0,
    /// If set, this tuple variation table's private point numbers follow
    /// immediately in the serialized data instead of using the shared set.
    private_point_numbers: bool,
    /// If set, this header embeds an intermediateStartTuple and
    /// intermediateEndTuple after the peak tuple (or after the header's
    /// fixed part, if embedded_peak_tuple is clear).
    intermediate_region: bool,
    /// If set, this header embeds its own peak Tuple record rather than
    /// referencing the shared tuples array by tuple_index. Mandatory (always
    /// set) in 'cvar'; optional in 'gvar'.
    embedded_peak_tuple: bool,
};

/// One position in the font's variation space: axisCount F2Dot14 values, in
/// 'fvar' axis order. table_data starts at this Tuple's own beginning.
pub const Tuple = struct {
    const Self = @This();

    table_data: []const u8,
    axis_count: u16,

    pub inline fn getValue(self: Self, axis: usize) types.F2Dot14 {
        return @bitCast(bytes.getI16BE(self.table_data[axis * @sizeOf(types.F2Dot14) ..]));
    }

    pub inline fn byteSize(self: Self) usize {
        return @as(usize, self.axis_count) * @sizeOf(types.F2Dot14);
    }
};

/// One TupleVariationHeader: describes the region of applicability for one
/// tuple variation data table, plus (via headerByteSize()) how much of the
/// tupleVariationHeaders array it occupies so the next header can be found.
/// table_data starts at this header's own beginning; axis_count is needed
/// (from 'fvar', since this format never carries it itself) to size any
/// embedded tuples.
pub const TupleVariationHeader = struct {
    const Self = @This();

    table_data: []const u8,
    axis_count: u16,

    const fixed_size = @sizeOf(u16) + @sizeOf(u16);

    pub inline fn variationDataSize(self: Self) u16 {
        return bytes.getU16BE(self.table_data[0..]);
    }

    pub inline fn tupleIndex(self: Self) TupleIndex {
        return @bitCast(bytes.getU16BE(self.table_data[@sizeOf(u16)..]));
    }

    fn tupleSize(self: Self) usize {
        return @as(usize, self.axis_count) * @sizeOf(types.F2Dot14);
    }

    /// Null unless tupleIndex().embedded_peak_tuple is set.
    pub fn getPeakTuple(self: Self) ?Tuple {
        if (!self.tupleIndex().embedded_peak_tuple) return null;
        return .{ .table_data = self.table_data[fixed_size..], .axis_count = self.axis_count };
    }

    fn intermediateStart(self: Self) usize {
        return fixed_size + if (self.tupleIndex().embedded_peak_tuple) self.tupleSize() else 0;
    }

    /// Null unless tupleIndex().intermediate_region is set.
    pub fn getIntermediateStartTuple(self: Self) ?Tuple {
        if (!self.tupleIndex().intermediate_region) return null;
        return .{ .table_data = self.table_data[self.intermediateStart()..], .axis_count = self.axis_count };
    }

    /// Null unless tupleIndex().intermediate_region is set.
    pub fn getIntermediateEndTuple(self: Self) ?Tuple {
        if (!self.tupleIndex().intermediate_region) return null;
        return .{ .table_data = self.table_data[self.intermediateStart() + self.tupleSize() ..], .axis_count = self.axis_count };
    }

    /// This header's own byte size within the tupleVariationHeaders array
    /// (fixed part plus whichever tuples are embedded) - the stride to the
    /// next header.
    pub fn headerByteSize(self: Self) usize {
        const idx = self.tupleIndex();
        var size: usize = fixed_size;
        if (idx.embedded_peak_tuple) size += self.tupleSize();
        if (idx.intermediate_region) size += 2 * self.tupleSize();
        return size;
    }
};

/// Schema for a GlyphVariationData/'cvar' header's fixed part, used only to
/// compute field offsets/sizes at comptime; never instantiated.
const DataHeader = extern struct {
    tuple_variation_count: u16,
    data_offset: types.Offset16,
};

/// A tupleVariationCount/dataOffset header, a variable-length
/// tupleVariationHeaders array, then serialized point-number/delta data.
/// Used as one glyph's data in 'gvar' and as the entire 'cvar' body.
pub const GlyphVariationData = struct {
    const Self = @This();

    table_data: []const u8,
    axis_count: u16,

    pub inline fn tupleVariationCount(self: Self) TupleVariationCount {
        return @bitCast(bytes.getU16BE(self.table_data[@offsetOf(DataHeader, "tuple_variation_count")..]));
    }

    /// Offset from the start of this GlyphVariationData table to the
    /// serialized data.
    pub inline fn dataOffset(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(DataHeader, "data_offset")..]);
    }

    pub fn getTupleVariationHeaderOffset(self: Self, index: usize) usize {
        var off: usize = @sizeOf(DataHeader);
        var i: usize = 0;
        while (i < index) : (i += 1) {
            const header: TupleVariationHeader = .{ .table_data = self.table_data[off..], .axis_count = self.axis_count };
            off += header.headerByteSize();
        }
        return off;
    }

    pub fn getTupleVariationHeader(self: Self, index: usize) TupleVariationHeader {
        return .{ .table_data = self.table_data[self.getTupleVariationHeaderOffset(index)..], .axis_count = self.axis_count };
    }

    /// The raw, un-decoded serialized data: packed point numbers and packed
    /// deltas. Interpreting this run-length-encoded byte stream is deferred,
    /// matching the CFF-charstrings precedent elsewhere in this font module.
    pub fn serializedData(self: Self) []const u8 {
        return self.table_data[self.dataOffset()..];
    }
};

const zigstd = @import("../../../../zigstd.zig");

test "DataHeader matches the GlyphVariationData stream layout" {
    try zigstd.testExpectEqual(@as(usize, 4), @sizeOf(DataHeader));
    try zigstd.testExpectEqual(@as(usize, 2), @offsetOf(DataHeader, "data_offset"));
}

fn f2dot14Bytes(value: f64) [2]u8 {
    var buf: [2]u8 = undefined;
    const packed_val = types.F2Dot14.fromFloat(value);
    bytes.setI16BE(&buf, @bitCast(packed_val));
    return buf;
}

test "TupleVariationHeader sizes and reads embedded peak + intermediate tuples" {
    const idx: TupleIndex = .{
        .tuple_index = 0,
        .private_point_numbers = false,
        .intermediate_region = true,
        .embedded_peak_tuple = true,
    };

    // 2 axes: fixed(4) + peak(4) + intermediateStart(4) + intermediateEnd(4) = 16
    var buf: [16]u8 = @splat(0);
    bytes.setU16BE(&buf, 12); // variationDataSize
    bytes.setU16BE(buf[2..], @bitCast(idx));

    @memcpy(buf[4..][0..2], &f2dot14Bytes(1.0));
    @memcpy(buf[6..][0..2], &f2dot14Bytes(0.5));
    @memcpy(buf[8..][0..2], &f2dot14Bytes(0.0));
    @memcpy(buf[10..][0..2], &f2dot14Bytes(0.0));
    @memcpy(buf[12..][0..2], &f2dot14Bytes(1.0));
    @memcpy(buf[14..][0..2], &f2dot14Bytes(1.0));

    const header: TupleVariationHeader = .{ .table_data = &buf, .axis_count = 2 };
    try zigstd.testExpectEqual(@as(u16, 12), header.variationDataSize());
    try zigstd.testExpectEqual(@as(usize, 16), header.headerByteSize());

    const peak = header.getPeakTuple().?;
    try zigstd.testExpectEqual(@as(f32, 1.0), peak.getValue(0).toFloat(f32));
    try zigstd.testExpectEqual(@as(f32, 0.5), peak.getValue(1).toFloat(f32));

    const start = header.getIntermediateStartTuple().?;
    try zigstd.testExpectEqual(@as(f32, 0.0), start.getValue(0).toFloat(f32));
    try zigstd.testExpectEqual(@as(f32, 0.0), start.getValue(1).toFloat(f32));

    const end = header.getIntermediateEndTuple().?;
    try zigstd.testExpectEqual(@as(f32, 1.0), end.getValue(0).toFloat(f32));
    try zigstd.testExpectEqual(@as(f32, 1.0), end.getValue(1).toFloat(f32));
}

test "GlyphVariationData walks a variable-length tupleVariationHeaders array" {
    // Header 0: no embedded tuple, no intermediate (fixed_size=4).
    // Header 1: embedded peak tuple only, 1 axis (fixed_size=4 + 2 = 6).
    const idx0: TupleIndex = .{ .tuple_index = 0, .private_point_numbers = false, .intermediate_region = false, .embedded_peak_tuple = false };
    const idx1: TupleIndex = .{ .tuple_index = 0, .private_point_numbers = false, .intermediate_region = false, .embedded_peak_tuple = true };

    const header_fixed = @sizeOf(DataHeader);
    const h0_size = 4;
    const h1_size = 6;
    const data_off = header_fixed + h0_size + h1_size;
    const total = data_off + 3; // 3 bytes of dummy serialized data

    var buf: [header_fixed + h0_size + h1_size + 3]u8 = @splat(0);
    try zigstd.testExpectEqual(@as(usize, total), buf.len);

    const tvc: TupleVariationCount = .{ .count = 2, .shared_point_numbers = false };
    bytes.setU16BE(&buf, @bitCast(tvc));
    bytes.setU16BE(buf[2..], @intCast(data_off));

    var off: usize = header_fixed;
    bytes.setU16BE(buf[off..], 100); // header0 variationDataSize
    bytes.setU16BE(buf[off + 2 ..], @bitCast(idx0));
    off += h0_size;

    bytes.setU16BE(buf[off..], 200); // header1 variationDataSize
    bytes.setU16BE(buf[off + 2 ..], @bitCast(idx1));
    @memcpy(buf[off + 4 ..][0..2], &f2dot14Bytes(-1.0));
    off += h1_size;

    const gvd: GlyphVariationData = .{ .table_data = &buf, .axis_count = 1 };
    try zigstd.testExpectEqual(@as(u12, 2), gvd.tupleVariationCount().count);
    try zigstd.testExpectEqual(@as(u16, @intCast(data_off)), gvd.dataOffset());

    const h0 = gvd.getTupleVariationHeader(0);
    try zigstd.testExpectEqual(@as(u16, 100), h0.variationDataSize());
    try zigstd.testExpectEqual(@as(?Tuple, null), h0.getPeakTuple());

    const h1 = gvd.getTupleVariationHeader(1);
    try zigstd.testExpectEqual(@as(u16, 200), h1.variationDataSize());
    try zigstd.testExpectEqual(@as(f32, -1.0), h1.getPeakTuple().?.getValue(0).toFloat(f32));

    try zigstd.testExpectEqual(@as(usize, 3), gvd.serializedData().len);
}
