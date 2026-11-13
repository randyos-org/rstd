const types = @import("../../types.zig");
const bytes = @import("../../../../../bytes.zig");

const GlyphIndex = types.GlyphIndex;
const zigstd = @import("../../../../../zigstd.zig");

/// The format-6 header.
const StreamLayout = extern struct {
    format: u16,
    length: u16,
    language: u16,
    first_code: u16,
    entry_count: u16,
};

const header_end = @sizeOf(StreamLayout);

// Trimmed table mapping: a single contiguous run of character codes.
pub const Format6 = struct {
    const Self = @This();

    table_data: []const u8,

    pub fn init(table_data: []const u8, check: bool) !Self {
        if (check) try bytes.checkSize(table_data, 0, StreamLayout);
        const self: Self = .{ .table_data = table_data };
        if (check) try bytes.checkArray(table_data, header_end, self.entryCount(), GlyphIndex);
        return self;
    }

    pub inline fn format(self: Self) u16 {
        return types.field(StreamLayout, "format", self.table_data);
    }

    pub inline fn length(self: Self) u16 {
        return types.field(StreamLayout, "length", self.table_data);
    }

    pub inline fn language(self: Self) u16 {
        return types.field(StreamLayout, "language", self.table_data);
    }

    pub inline fn firstCode(self: Self) u16 {
        return types.field(StreamLayout, "first_code", self.table_data);
    }

    pub inline fn entryCount(self: Self) u16 {
        return types.field(StreamLayout, "entry_count", self.table_data);
    }

    /// glyphIdArray[index], with index counted from the start of the array
    /// (not from firstCode()). getGlyphIndex() is what most callers want.
    pub fn getRawGlyphId(self: Self, index: usize) GlyphIndex {
        return @fromBackingInt(@intCast(bytes.getU16BE(self.table_data[header_end + index * @sizeOf(GlyphIndex) ..])));
    }

    pub fn getGlyphIndex(self: Self, codepoint: u32) GlyphIndex {
        if (codepoint > 0xFFFF) return .notdef;
        const c: u16 = @intCast(codepoint);

        const first_code = self.firstCode();
        if (c < first_code) return .notdef;
        const idx = c - first_code;
        if (idx >= self.entryCount()) return .notdef;

        return self.getRawGlyphId(idx);
    }
};

const setU16BE = bytes.setU16BE;

test "Format6 header matches the stream layout (10 bytes, no padding)" {
    try zigstd.testExpectEqual(@as(usize, 10), @sizeOf(StreamLayout));
    try zigstd.testExpectEqual(@as(usize, 6), @offsetOf(StreamLayout, "first_code"));
    try zigstd.testExpectEqual(@as(usize, 8), @offsetOf(StreamLayout, "entry_count"));
}

fn buildFormat6(alloc: zigstd.Allocator, first_code: u16, glyph_ids: []const u16) ![]u8 {
    const header_size = 10;
    const total = header_size + glyph_ids.len * 2;

    const buf = try alloc.alloc(u8, total);
    @memset(buf, 0);

    setU16BE(buf, 6); // format
    setU16BE(buf[2..], @intCast(total)); // length
    setU16BE(buf[4..], 0); // language
    setU16BE(buf[6..], first_code);
    setU16BE(buf[8..], @intCast(glyph_ids.len));

    var off: usize = header_size;
    for (glyph_ids) |v| {
        setU16BE(buf[off..], v);
        off += 2;
    }

    return buf;
}

test "Format6.getGlyphIndex resolves codes within the trimmed range" {
    const alloc = zigstd.testAllocator;
    const buf = try buildFormat6(alloc, 100, &[_]u16{ 10, 11, 12 });
    defer alloc.free(buf);

    const fmt6 = try Format6.init(buf, true);

    try zigstd.testExpectEqual(@as(u16, 10), @backingInt(fmt6.getGlyphIndex(100)));
    try zigstd.testExpectEqual(@as(u16, 12), @backingInt(fmt6.getGlyphIndex(102)));
}

test "Format6.getGlyphIndex returns .notdef outside the trimmed range" {
    const alloc = zigstd.testAllocator;
    const buf = try buildFormat6(alloc, 100, &[_]u16{ 10, 11, 12 });
    defer alloc.free(buf);

    const fmt6 = try Format6.init(buf, true);

    try zigstd.testExpectEqual(GlyphIndex.notdef, fmt6.getGlyphIndex(99));
    try zigstd.testExpectEqual(GlyphIndex.notdef, fmt6.getGlyphIndex(103));
    // past the 16-bit range this format can represent at all
    try zigstd.testExpectEqual(GlyphIndex.notdef, fmt6.getGlyphIndex(0x10000));
}
