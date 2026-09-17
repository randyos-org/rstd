const types = @import("../../stream_types.zig");
const bytes = @import("../../../../../bytes.zig");
const zigstd = @import("../../../../../zigstd.zig");

const GlyphIndex = types.GlyphIndex;
const setU16BE = bytes.setU16BE;
const setU32BE = bytes.setU32BE;

/// Schema for the format-10 header, used only to compute field offsets/sizes
/// at comptime; never instantiated.
const StreamLayout = extern struct {
    format: u16,
    reserved: u16,
    length: u32,
    language: u32,
    start_char_code: u32,
    num_chars: u32,
};

const header_end = @sizeOf(StreamLayout);

// Trimmed array: format 6's contiguous-range trick, extended to 32-bit codes.
pub const Format10 = struct {
    const Self = @This();

    table_data: []const u8,

    pub fn init(table_data: []const u8, check: bool) !Self {
        if (check) try bytes.checkSize(table_data, 0, StreamLayout);
        const self: Self = .{ .table_data = table_data };
        if (check) try bytes.checkArray(table_data, header_end, self.numChars(), GlyphIndex);
        return self;
    }

    pub inline fn format(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "format")..]);
    }

    pub inline fn reserved(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "reserved")..]);
    }

    pub inline fn length(self: Self) u32 {
        return bytes.getU32BE(self.table_data[@offsetOf(StreamLayout, "length")..]);
    }

    pub inline fn language(self: Self) u32 {
        return bytes.getU32BE(self.table_data[@offsetOf(StreamLayout, "language")..]);
    }

    pub inline fn startCharCode(self: Self) u32 {
        return bytes.getU32BE(self.table_data[@offsetOf(StreamLayout, "start_char_code")..]);
    }

    pub inline fn numChars(self: Self) u32 {
        return bytes.getU32BE(self.table_data[@offsetOf(StreamLayout, "num_chars")..]);
    }

    /// glyphIdArray[index], with index counted from the start of the array
    /// (not from startCharCode()). getGlyphIndex() is what most callers want.
    pub fn getRawGlyphId(self: Self, index: usize) GlyphIndex {
        return @enumFromInt(bytes.getU16BE(self.table_data[header_end + index * @sizeOf(GlyphIndex) ..]));
    }

    pub fn getGlyphIndex(self: Self, codepoint: u32) GlyphIndex {
        const start_char_code = self.startCharCode();
        if (codepoint < start_char_code) return .notdef;
        const idx = codepoint - start_char_code;
        if (idx >= self.numChars()) return .notdef;

        return self.getRawGlyphId(idx);
    }
};

test "Format10 header matches the stream layout (20 bytes, no padding)" {
    try zigstd.testExpectEqual(@as(usize, 20), @sizeOf(StreamLayout));
    try zigstd.testExpectEqual(@as(usize, 12), @offsetOf(StreamLayout, "start_char_code"));
    try zigstd.testExpectEqual(@as(usize, 16), @offsetOf(StreamLayout, "num_chars"));
}

fn buildFormat10(alloc: zigstd.Allocator, start_char_code: u32, glyph_ids: []const u16) ![]u8 {
    const header_size = 20;
    const total = header_size + glyph_ids.len * 2;

    const buf = try alloc.alloc(u8, total);
    @memset(buf, 0);

    setU16BE(buf, 10); // format
    setU16BE(buf[2..], 0); // reserved
    setU32BE(buf[4..], @intCast(total)); // length
    setU32BE(buf[8..], 0); // language
    setU32BE(buf[12..], start_char_code);
    setU32BE(buf[16..], @intCast(glyph_ids.len));

    var off: usize = header_size;
    for (glyph_ids) |v| {
        setU16BE(buf[off..], v);
        off += 2;
    }

    return buf;
}

test "Format10.getGlyphIndex resolves codes within the trimmed 32-bit range" {
    const alloc = zigstd.testAllocator;
    const buf = try buildFormat10(alloc, 0x20000, &[_]u16{ 10, 11, 12 });
    defer alloc.free(buf);

    const fmt10 = try Format10.init(buf, true);

    try zigstd.testExpectEqual(@as(u16, 10), @backingInt(fmt10.getGlyphIndex(0x20000)));
    try zigstd.testExpectEqual(@as(u16, 12), @backingInt(fmt10.getGlyphIndex(0x20002)));
    try zigstd.testExpectEqual(GlyphIndex.notdef, fmt10.getGlyphIndex(0x20003));
    try zigstd.testExpectEqual(GlyphIndex.notdef, fmt10.getGlyphIndex(0x1FFFF));
}
