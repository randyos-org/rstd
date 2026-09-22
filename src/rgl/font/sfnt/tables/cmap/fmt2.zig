const types = @import("../../stream_types.zig");
const bytes = @import("../../../../../bytes.zig");
const zigstd = @import("../../../../../zigstd.zig");

const GlyphIndex = types.GlyphIndex;
const setU16BE = bytes.setU16BE;
const setI16BE = bytes.setI16BE;

/// Schema for the format-2 header, used only to compute field offsets/sizes
/// at comptime; never instantiated.
const StreamLayout = extern struct {
    format: u16,
    length: u16,
    language: u16,
    sub_header_keys: [256]u16,
};

const header_end = @sizeOf(StreamLayout);

/// One subheader record. Its extern layout doubles as the wire schema for
/// @offsetOf/@sizeOf.
pub const SubHeader = extern struct {
    first_code: u16,
    entry_count: u16,
    id_delta: i16,
    id_range_offset: u16,
};

/// A codepoint split into the leading byte that selects a subheader and the
/// trailing byte matched against it. First declared field gets the
/// low-order bits in a packed struct, so this lines up with @bitCast(c)
/// directly - no manual >>8/&0xFF needed.
const HiLo = packed struct(u16) {
    lo: u8,
    hi: u8,

    /// v must already fit in u16 - callers are expected to have rejected
    /// anything past whatever range their format can represent first.
    fn fromU32(v: u32) HiLo {
        return @bitCast(@as(u16, @intCast(v)));
    }
};

pub const Format2 = struct {
    const Self = @This();

    table_data: []const u8,
    /// Nothing in the spec gives an explicit subheader count; every
    /// implementation infers it as the highest sub_header_keys[] entry
    /// (divided by the 8-byte record size), plus one for the subheader that
    /// lives at that offset. Computed once at init() rather than rescanning
    /// all 256 keys on every lookup.
    num_sub_headers: usize,

    pub fn init(table_data: []const u8, check: bool) !Self {
        if (check) try bytes.checkSize(table_data, 0, StreamLayout);

        var max_key: u16 = 0;
        var hi: usize = 0;
        while (hi < 256) : (hi += 1) {
            const key = bytes.getU16BE(table_data[@offsetOf(StreamLayout, "sub_header_keys") + hi * @sizeOf(u16) ..]);
            if (key > max_key) max_key = key;
        }
        const num_sub_headers = @as(usize, max_key) / @sizeOf(SubHeader) + 1;

        if (check) try bytes.checkArray(table_data, header_end, num_sub_headers, SubHeader);

        return .{ .table_data = table_data, .num_sub_headers = num_sub_headers };
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

    pub fn getSubHeaderKey(self: Self, hi: usize) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "sub_header_keys") + hi * @sizeOf(u16) ..]);
    }

    pub fn getSubHeader(self: Self, index: usize) SubHeader {
        return bytes.typeBEToNative(SubHeader, self.table_data[header_end + index * @sizeOf(SubHeader) ..]);
    }

    pub fn getGlyphIndex(self: Self, codepoint: u32) GlyphIndex {
        if (codepoint > 0xFFFF) return .notdef;
        const c = HiLo.fromU32(codepoint);

        const sh_index = @as(usize, self.getSubHeaderKey(c.hi)) / @sizeOf(SubHeader);
        if (sh_index >= self.num_sub_headers) return .notdef;
        const sh = self.getSubHeader(sh_index);

        // The low byte is what's actually matched against the subheader's
        // range: hi == 0 (single-byte codes) routes to sub_header_keys[0],
        // which conventionally points at a subheader covering the low byte
        // directly; any other hi byte selects a subheader for that leading
        // byte's two-byte range instead.
        if (c.lo < sh.first_code) return .notdef;
        const rel = c.lo - sh.first_code;
        if (rel >= sh.entry_count) return .notdef;

        // idRangeOffset is a byte offset measured from the location of the
        // idRangeOffset field itself in the stream - the spec's own
        // pointer-following convention, usable directly now that we read
        // straight from table_data instead of separately-allocated arrays.
        const id_range_offset_field_pos = header_end + sh_index * @sizeOf(SubHeader) + @offsetOf(SubHeader, "id_range_offset");
        const glyph_array_start = header_end + self.num_sub_headers * @sizeOf(SubHeader);
        const target = id_range_offset_field_pos + @as(usize, sh.id_range_offset) + @as(usize, rel) * @sizeOf(GlyphIndex);
        if (target < glyph_array_start or target + @sizeOf(GlyphIndex) > self.table_data.len) return .notdef;

        const raw: GlyphIndex = @fromBackingInt(@intCast(bytes.getU16BE(self.table_data[target..])));
        if (raw == .notdef) return .notdef;
        const delta: u16 = @bitCast(sh.id_delta);
        return @fromBackingInt(@backingInt(raw) +% delta);
    }
};

test "Format2 StreamLayout/SubHeader match the stream layout (no padding)" {
    try zigstd.testExpectEqual(@as(usize, 518), @sizeOf(StreamLayout));
    try zigstd.testExpectEqual(@as(usize, 6), @offsetOf(StreamLayout, "sub_header_keys"));
    try zigstd.testExpectEqual(@as(usize, 8), @sizeOf(SubHeader));
    try zigstd.testExpectEqual(@as(usize, 6), @offsetOf(SubHeader, "id_range_offset"));
}

// Builds a format-2 subtable with two subheaders: subheader 0 is the
// conventional single-byte entry reached via hi byte 0, and subheader 1 is
// assigned to hi byte 0x81 for a genuine two-byte range.
fn buildFormat2(alloc: zigstd.Allocator) ![]u8 {
    const num_sub_headers = 2;
    const glyph_ids = [_]u16{ 100, 101, 0, 300, 0 };
    const glyph_array_start = header_end + num_sub_headers * @sizeOf(SubHeader);
    const total = glyph_array_start + glyph_ids.len * @sizeOf(GlyphIndex);

    const buf = try alloc.alloc(u8, total);
    @memset(buf, 0);

    setU16BE(buf, 2); // format
    setU16BE(buf[@offsetOf(StreamLayout, "length")..], @intCast(total));
    setU16BE(buf[@offsetOf(StreamLayout, "language")..], 0);

    const sub_header_keys_start = @offsetOf(StreamLayout, "sub_header_keys");
    setU16BE(buf[sub_header_keys_start + 0x81 * @sizeOf(u16) ..], @sizeOf(SubHeader)); // hi byte 0x81 -> subheader 1

    const sh0 = header_end;
    const sh1 = header_end + @sizeOf(SubHeader);

    // subheader 0: firstCode 0x41, entryCount 3, idDelta 0, idRangeOffset 10
    // (10 = distance from its own idRangeOffset field, at 524, to glyph_id_array[0] at 534)
    setU16BE(buf[sh0 + @offsetOf(SubHeader, "first_code") ..], 0x41);
    setU16BE(buf[sh0 + @offsetOf(SubHeader, "entry_count") ..], 3);
    setI16BE(buf[sh0 + @offsetOf(SubHeader, "id_delta") ..], 0);
    setU16BE(buf[sh0 + @offsetOf(SubHeader, "id_range_offset") ..], 10);

    // subheader 1: firstCode 0x30, entryCount 2, idDelta 5, idRangeOffset 8
    // (8 = distance from its own idRangeOffset field, at 532, to glyph_id_array[3] at 540)
    setU16BE(buf[sh1 + @offsetOf(SubHeader, "first_code") ..], 0x30);
    setU16BE(buf[sh1 + @offsetOf(SubHeader, "entry_count") ..], 2);
    setI16BE(buf[sh1 + @offsetOf(SubHeader, "id_delta") ..], 5);
    setU16BE(buf[sh1 + @offsetOf(SubHeader, "id_range_offset") ..], 8);

    var off: usize = glyph_array_start;
    for (glyph_ids) |v| {
        setU16BE(buf[off..], v);
        off += @sizeOf(GlyphIndex);
    }

    return buf;
}

test "Format2.getGlyphIndex resolves single-byte codes via the hi==0 subheader" {
    const alloc = zigstd.testAllocator;
    const buf = try buildFormat2(alloc);
    defer alloc.free(buf);

    const fmt2 = try Format2.init(buf, true);

    try zigstd.testExpectEqual(@as(u16, 100), @backingInt(fmt2.getGlyphIndex(0x0041)));
    try zigstd.testExpectEqual(@as(u16, 101), @backingInt(fmt2.getGlyphIndex(0x0042)));
    // raw glyph_id_array entry is 0 ("missing"), so it stays .notdef regardless of idDelta
    try zigstd.testExpectEqual(GlyphIndex.notdef, fmt2.getGlyphIndex(0x0043));
    // low byte past the subheader's entryCount
    try zigstd.testExpectEqual(GlyphIndex.notdef, fmt2.getGlyphIndex(0x0044));
}

test "Format2.getGlyphIndex resolves two-byte codes through their own subheader" {
    const alloc = zigstd.testAllocator;
    const buf = try buildFormat2(alloc);
    defer alloc.free(buf);

    const fmt2 = try Format2.init(buf, true);

    try zigstd.testExpectEqual(@as(u16, 305), @backingInt(fmt2.getGlyphIndex(0x8130))); // 300 + idDelta 5
    try zigstd.testExpectEqual(GlyphIndex.notdef, fmt2.getGlyphIndex(0x8131)); // raw 0 -> .notdef
    try zigstd.testExpectEqual(GlyphIndex.notdef, fmt2.getGlyphIndex(0x8132)); // past entryCount
}

test "Format2.getGlyphIndex returns .notdef for a hi byte without its own subheader" {
    const alloc = zigstd.testAllocator;
    const buf = try buildFormat2(alloc);
    defer alloc.free(buf);

    const fmt2 = try Format2.init(buf, true);

    // hi byte 0x82 falls back to subheader 0, whose range starts at 0x41
    try zigstd.testExpectEqual(GlyphIndex.notdef, fmt2.getGlyphIndex(0x8230));
}

test "Format2.getGlyphIndex returns .notdef for a codepoint past the 16-bit range" {
    const alloc = zigstd.testAllocator;
    const buf = try buildFormat2(alloc);
    defer alloc.free(buf);

    const fmt2 = try Format2.init(buf, true);

    try zigstd.testExpectEqual(GlyphIndex.notdef, fmt2.getGlyphIndex(0x10000));
}
