const types = @import("../../stream_types.zig");
const bytes = @import("../../../../../bytes.zig");
const zigstd = @import("../../../../../zigstd.zig");

const GlyphIndex = types.GlyphIndex;

/// Schema for the format-0 header, used only to compute field offsets/sizes
/// at comptime; never instantiated.
const StreamLayout = extern struct {
    format: u16,
    length: u16,
    language: u16,
};

const header_end = @sizeOf(StreamLayout);
const glyph_id_array_size = 256;

pub const Format0 = struct {
    const Self = @This();

    table_data: []const u8,

    pub fn init(table_data: []const u8, check: bool) !Self {
        if (check) try bytes.checkSlice(table_data, header_end, glyph_id_array_size);
        return .{ .table_data = table_data };
    }

    pub inline fn format(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "format")..]);
    }

    pub inline fn length(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "length")..]);
    }

    pub inline fn language(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "language")..]);
    }

    /// The raw glyphIdArray[codepoint] byte itself, before being widened into
    /// a GlyphIndex. getGlyphIndex() is what most callers want.
    pub inline fn getRawGlyphId(self: Self, codepoint: u8) u8 {
        return self.table_data[header_end + codepoint];
    }

    pub fn getGlyphIndex(self: Self, codepoint: u32) GlyphIndex {
        if (codepoint > 0xFF) return .notdef;
        return @fromBackingInt(self.table_data[header_end + codepoint]);
    }
};

test "Format0 stream layout matches the OpenType stream layout (262 bytes total, no padding)" {
    try zigstd.testExpectEqual(@as(usize, 6), @sizeOf(StreamLayout));
    try zigstd.testExpectEqual(@as(usize, 4), @offsetOf(StreamLayout, "language"));
    try zigstd.testExpectEqual(@as(usize, 262), header_end + glyph_id_array_size);
}

test "Format0.getGlyphIndex indexes the byte array directly" {
    var buf: [262]u8 = @splat(0);
    buf[6 + 0x41] = 10; // 'A' -> glyph 10

    const fmt0 = try Format0.init(buf[0..], true);

    try zigstd.testExpectEqual(@as(u16, 10), @backingInt(fmt0.getGlyphIndex(0x41)));
    try zigstd.testExpectEqual(GlyphIndex.notdef, fmt0.getGlyphIndex(0x42));
}

test "Format0.getGlyphIndex returns .notdef past the single-byte range" {
    var buf: [262]u8 = @splat(0);

    const fmt0 = try Format0.init(buf[0..], true);

    try zigstd.testExpectEqual(GlyphIndex.notdef, fmt0.getGlyphIndex(0x100));
}
