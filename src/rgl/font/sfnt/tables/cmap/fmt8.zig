const types = @import("../../stream_types.zig");
const bytes = @import("../../../../../bytes.zig");
const zigstd = @import("../../../../../zigstd.zig");

const GlyphIndex = types.GlyphIndex;
const setU16BE = bytes.setU16BE;
const setU32BE = bytes.setU32BE;

/// Schema for the format-8 header, used only to compute field offsets/sizes
/// at comptime; never instantiated.
const StreamLayout = extern struct {
    format: u16,
    reserved: u16,
    length: u32,
    language: u32,
};

/// is32 (an 8192-byte bitmap marking which 16-bit words are the leading half
/// of a 32-bit character code) sits right after the header. Per the spec,
/// this determination "can be made by looking at the 16-bit value directly,
/// with no further information required" - i.e. it's redundant with the
/// standard UTF-16 surrogate range - so it's not consulted by getGlyphIndex()
/// here, which already receives a resolved 32-bit codepoint rather than a
/// raw UTF-16 code unit. It's still exposed via isLeadingSurrogateWord() for
/// callers walking raw UTF-16 text who want the font's own bitmap.
const is32_start = @sizeOf(StreamLayout);
const is32_size = 8192;
const num_groups_pos = is32_start + is32_size;
const groups_start = num_groups_pos + @sizeOf(u32);

/// One sequential map group. Its extern layout doubles as the wire schema
/// for @offsetOf/@sizeOf.
pub const Group = extern struct {
    start_char_code: u32,
    end_char_code: u32,
    start_glyph_id: u32,
};

// Mixed 16-bit and 32-bit coverage: rare, used by a handful of CJK fonts that
// need both BMP and supplementary-plane characters in one subtable.
pub const Format8 = struct {
    const Self = @This();

    table_data: []const u8,

    pub fn init(table_data: []const u8, check: bool) !Self {
        if (check) try bytes.checkSize(table_data, num_groups_pos, u32);
        const self: Self = .{ .table_data = table_data };
        if (check) try bytes.checkArray(table_data, groups_start, self.numGroups(), Group);
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

    /// Per the spec: "To determine if a particular word (cp) is the first
    /// half of 32-bit code points, one can use an expression such as
    /// (is32[cp/8] & (1 << (7 - (cp%8))))." True means `word` is the leading
    /// half of a 32-bit character code; false means it's a standalone 16-bit
    /// character code.
    pub fn isLeadingSurrogateWord(self: Self, word: u16) bool {
        const byte = self.table_data[is32_start + word / 8];
        return (byte & (@as(u8, 1) << @intCast(7 - (word % 8)))) != 0;
    }

    pub inline fn numGroups(self: Self) u32 {
        return bytes.getU32BE(self.table_data[num_groups_pos..]);
    }

    pub fn getGroup(self: Self, index: usize) Group {
        return bytes.typeBEToNative(Group, self.table_data[groups_start + index * @sizeOf(Group) ..]);
    }

    pub fn getGlyphIndex(self: Self, codepoint: u32) GlyphIndex {
        var i: usize = 0;
        while (i < self.numGroups()) : (i += 1) {
            const group = self.getGroup(i);
            if (codepoint >= group.start_char_code and codepoint <= group.end_char_code) {
                const raw = group.start_glyph_id + (codepoint - group.start_char_code);
                return @fromBackingInt(@as(u16, @truncate(raw)));
            }
        }
        return .notdef;
    }
};

test "Format8 header and group match the stream layout (no padding)" {
    try zigstd.testExpectEqual(@as(usize, 12), @sizeOf(StreamLayout));
    try zigstd.testExpectEqual(@as(usize, 4), @offsetOf(StreamLayout, "length"));
    try zigstd.testExpectEqual(@as(usize, 8), @offsetOf(StreamLayout, "language"));
    try zigstd.testExpectEqual(@as(usize, 12), @sizeOf(Group));
    try zigstd.testExpectEqual(@as(usize, 8), @offsetOf(Group, "start_glyph_id"));
}

fn buildFormat8(alloc: zigstd.Allocator, groups: []const Group) ![]u8 {
    const header_size = 12;
    const total = groups_start + groups.len * 12;

    const buf = try alloc.alloc(u8, total);
    @memset(buf, 0);

    setU16BE(buf, 8); // format
    setU16BE(buf[2..], 0); // reserved
    setU32BE(buf[4..], @intCast(total)); // length
    setU32BE(buf[8..], 0); // language
    setU32BE(buf[header_size + is32_size ..], @intCast(groups.len)); // numGroups

    var off: usize = groups_start;
    for (groups) |g| {
        setU32BE(buf[off..], g.start_char_code);
        setU32BE(buf[off + 4 ..], g.end_char_code);
        setU32BE(buf[off + 8 ..], g.start_glyph_id);
        off += 12;
    }

    return buf;
}

test "Format8.getGlyphIndex resolves a supplementary-plane codepoint through a group" {
    const alloc = zigstd.testAllocator;
    const buf = try buildFormat8(alloc, &[_]Group{
        .{ .start_char_code = 0x10000, .end_char_code = 0x10010, .start_glyph_id = 500 },
    });
    defer alloc.free(buf);

    const fmt8 = try Format8.init(buf, true);

    try zigstd.testExpectEqual(@as(u16, 500), @backingInt(fmt8.getGlyphIndex(0x10000)));
    try zigstd.testExpectEqual(@as(u16, 505), @backingInt(fmt8.getGlyphIndex(0x10005)));
    try zigstd.testExpectEqual(GlyphIndex.notdef, fmt8.getGlyphIndex(0x10011));
}
