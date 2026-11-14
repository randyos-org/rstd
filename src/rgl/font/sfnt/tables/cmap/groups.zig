//! The cmap subtables of sorted code-range groups, formats 8, 12 and 13:
//! 8 and 12 map a group's codes to consecutive glyphs, 13 all to one. Format
//! 8 carries the `is32` bitmap between its header and group count.

const types = @import("../../types.zig");
const bytes = @import("../../../../../bytes.zig");
const zigstd = @import("../../../../../zigstd.zig");

const GlyphIndex = types.GlyphIndex;

/// The header the three formats share; 12 and 13 put the group count next.
const Header = extern struct {
    format: u16,
    reserved: u16,
    length: u32,
    language: u32,
};

/// Format 8's `is32`: 8192 bytes marking which 16-bit words lead a 32-bit
/// code. It repeats the UTF-16 surrogate ranges, so `getGlyphIndex`, which
/// takes a resolved code point, never reads it.
const is32_size = 8192;

/// One map group. Its extern layout doubles as the wire schema for
/// @offsetOf/@sizeOf.
pub const Group = extern struct {
    start_char_code: u32,
    end_char_code: u32,
    /// The first code's glyph (formats 8 and 12), or every code's (13).
    glyph_id: u32,
};

/// The format-`fmt` subtable.
pub fn Format(comptime fmt: u8) type {
    if (fmt != 8 and fmt != 12 and fmt != 13) @compileError("no group format " ++ @import("std").fmt.comptimePrint("{d}", .{fmt}));
    const count_at = @sizeOf(Header) + if (fmt == 8) is32_size else 0;
    const groups_start = count_at + @sizeOf(u32);
    return struct {
        const Self = @This();

        table_data: []const u8,

        pub fn init(table_data: []const u8, check: bool) !Self {
            if (check) try bytes.checkSize(table_data, count_at, u32);
            const self: Self = .{ .table_data = table_data };
            if (check) try bytes.checkArray(table_data, groups_start, self.numGroups(), Group);
            return self;
        }

        pub inline fn format(self: Self) u16 {
            return types.field(Header, "format", self.table_data);
        }

        pub inline fn reserved(self: Self) u16 {
            return types.field(Header, "reserved", self.table_data);
        }

        pub inline fn length(self: Self) u32 {
            return types.field(Header, "length", self.table_data);
        }

        pub inline fn language(self: Self) u32 {
            return types.field(Header, "language", self.table_data);
        }

        /// Format 8 only, per the spec: `is32[cp/8] & (1 << (7 - cp%8))`. True
        /// when `word` leads a 32-bit code, for callers walking raw UTF-16.
        pub const isLeadingSurrogateWord = if (fmt == 8) leadingWord else {};

        fn leadingWord(self: Self, word: u16) bool {
            const byte = self.table_data[@sizeOf(Header) + word / 8];
            return (byte & (@as(u8, 1) << @intCast(7 - (word % 8)))) != 0;
        }

        pub inline fn numGroups(self: Self) u32 {
            return bytes.getU32BE(self.table_data[count_at..]);
        }

        pub fn getGroup(self: Self, index: usize) Group {
            return bytes.typeBEToNative(Group, self.table_data[groups_start + index * @sizeOf(Group) ..]);
        }

        /// Groups are sorted and non-overlapping per spec.
        pub fn getGlyphIndex(self: Self, codepoint: u32) GlyphIndex {
            const Find = struct {
                t: Self,
                cp: u32,
                pub fn order(f: @This(), i: usize) types.Order {
                    const g = f.t.getGroup(i);
                    return if (f.cp < g.start_char_code) .lt else if (f.cp > g.end_char_code) .gt else .eq;
                }
            };
            const i = types.search(self.numGroups(), Find{ .t = self, .cp = codepoint }) orelse return .notdef;
            const g = self.getGroup(i);
            const raw = if (fmt == 13) g.glyph_id else g.glyph_id + (codepoint - g.start_char_code);
            return @fromBackingInt(@as(u16, @truncate(raw)));
        }
    };
}

const setU16BE = bytes.setU16BE;
const setU32BE = bytes.setU32BE;

test "group formats' header and group match the stream layout (no padding)" {
    try zigstd.testExpectEqual(@as(usize, 12), @sizeOf(Header));
    try zigstd.testExpectEqual(@as(usize, 4), @offsetOf(Header, "length"));
    try zigstd.testExpectEqual(@as(usize, 8), @offsetOf(Header, "language"));
    try zigstd.testExpectEqual(@as(usize, 12), @sizeOf(Group));
    try zigstd.testExpectEqual(@as(usize, 8), @offsetOf(Group, "glyph_id"));
}

fn build(alloc: zigstd.Allocator, comptime fmt: u8, groups: []const Group) ![]u8 {
    const count_at = @sizeOf(Header) + if (fmt == 8) is32_size else 0;
    const total = count_at + 4 + groups.len * @sizeOf(Group);
    const buf = try alloc.alloc(u8, total);
    @memset(buf, 0);
    setU16BE(buf, fmt);
    setU32BE(buf[4..], @intCast(total));
    setU32BE(buf[count_at..], @intCast(groups.len));
    var off: usize = count_at + 4;
    for (groups) |g| {
        setU32BE(buf[off..], g.start_char_code);
        setU32BE(buf[off + 4 ..], g.end_char_code);
        setU32BE(buf[off + 8 ..], g.glyph_id);
        off += @sizeOf(Group);
    }
    return buf;
}

test "Format(8).getGlyphIndex resolves a supplementary-plane codepoint through a group" {
    const alloc = zigstd.testAllocator;
    const buf = try build(alloc, 8, &.{.{ .start_char_code = 0x10000, .end_char_code = 0x10010, .glyph_id = 500 }});
    defer alloc.free(buf);
    const fmt8 = try Format(8).init(buf, true);
    try zigstd.testExpectEqual(@as(u16, 500), @backingInt(fmt8.getGlyphIndex(0x10000)));
    try zigstd.testExpectEqual(@as(u16, 505), @backingInt(fmt8.getGlyphIndex(0x10005)));
    try zigstd.testExpectEqual(GlyphIndex.notdef, fmt8.getGlyphIndex(0x10011));
}

test "Format(12).getGlyphIndex resolves through sparse, non-contiguous groups" {
    const alloc = zigstd.testAllocator;
    const buf = try build(alloc, 12, &.{
        .{ .start_char_code = 0x41, .end_char_code = 0x5A, .glyph_id = 10 },
        .{ .start_char_code = 0x10000, .end_char_code = 0x10010, .glyph_id = 500 },
    });
    defer alloc.free(buf);
    const fmt12 = try Format(12).init(buf, true);
    try zigstd.testExpectEqual(@as(u16, 10), @backingInt(fmt12.getGlyphIndex(0x41)));
    try zigstd.testExpectEqual(@as(u16, 35), @backingInt(fmt12.getGlyphIndex(0x5A)));
    try zigstd.testExpectEqual(@as(u16, 505), @backingInt(fmt12.getGlyphIndex(0x10005)));
    try zigstd.testExpectEqual(GlyphIndex.notdef, fmt12.getGlyphIndex(0x5B));
    try zigstd.testExpectEqual(GlyphIndex.notdef, fmt12.getGlyphIndex(0xFFFF));
}

test "Format(13).getGlyphIndex returns the same glyph for every code in a group" {
    const alloc = zigstd.testAllocator;
    const buf = try build(alloc, 13, &.{.{ .start_char_code = 0x3400, .end_char_code = 0x4DBF, .glyph_id = 1 }});
    defer alloc.free(buf);
    const fmt13 = try Format(13).init(buf, true);
    try zigstd.testExpectEqual(@as(u16, 1), @backingInt(fmt13.getGlyphIndex(0x3400)));
    try zigstd.testExpectEqual(@as(u16, 1), @backingInt(fmt13.getGlyphIndex(0x4000)));
    try zigstd.testExpectEqual(@as(u16, 1), @backingInt(fmt13.getGlyphIndex(0x4DBF)));
    try zigstd.testExpectEqual(GlyphIndex.notdef, fmt13.getGlyphIndex(0x4DC0));
}
