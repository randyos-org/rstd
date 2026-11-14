//! The cmap subtables of one contiguous run of codes and a glyph array:
//! format 0 (all 256 single-byte codes, byte glyphs), 6 (a 16-bit run) and
//! 10 (a 32-bit run).

const std = @import("std");
const types = @import("../../types.zig");
const bytes = @import("../../../../../bytes.zig");
const zigstd = @import("../../../../../zigstd.zig");

const GlyphIndex = types.GlyphIndex;

/// Format `fmt`'s header.
fn Layout(comptime fmt: u8) type {
    return switch (fmt) {
        0 => extern struct { format: u16, length: u16, language: u16 },
        6 => extern struct { format: u16, length: u16, language: u16, first_code: u16, entry_count: u16 },
        10 => extern struct { format: u16, reserved: u16, length: u32, language: u32, start_char_code: u32, num_chars: u32 },
        else => @compileError(std.fmt.comptimePrint("no trimmed format {d}", .{fmt})),
    };
}

/// The format-`fmt` subtable.
pub fn Format(comptime fmt: u8) type {
    const L = Layout(fmt);
    const header_end = @sizeOf(L);
    // format 0's glyphs are bytes, the others' `GlyphIndex`es
    const Entry = if (fmt == 0) u8 else u16;
    // the highest code the format can hold
    const max_code: u32 = switch (fmt) {
        0 => 0xFF,
        6 => 0xFFFF,
        else => std.math.maxInt(u32),
    };
    return struct {
        const Self = @This();

        table_data: []const u8,

        pub fn init(table_data: []const u8, check: bool) !Self {
            if (check) try bytes.checkSize(table_data, 0, L);
            const self: Self = .{ .table_data = table_data };
            if (check) try bytes.checkArray(table_data, header_end, self.count(), Entry);
            return self;
        }

        pub inline fn format(self: Self) u16 {
            return types.field(L, "format", self.table_data);
        }

        pub inline fn length(self: Self) @FieldType(L, "length") {
            return types.field(L, "length", self.table_data);
        }

        pub inline fn language(self: Self) @FieldType(L, "language") {
            return types.field(L, "language", self.table_data);
        }

        /// The first code (`firstCode`, `startCharCode`; 0 for format 0).
        pub fn first(self: Self) u32 {
            return switch (fmt) {
                0 => 0,
                6 => types.field(L, "first_code", self.table_data),
                else => types.field(L, "start_char_code", self.table_data),
            };
        }

        /// The glyph array's length (`entryCount`, `numChars`; 256 for format 0).
        pub fn count(self: Self) u32 {
            return switch (fmt) {
                0 => 256,
                6 => types.field(L, "entry_count", self.table_data),
                else => types.field(L, "num_chars", self.table_data),
            };
        }

        /// glyphIdArray[index], counted from the array's start (not from
        /// `first()`). getGlyphIndex() is what most callers want.
        pub fn getRawGlyphId(self: Self, index: usize) GlyphIndex {
            const at = self.table_data[header_end + index * @sizeOf(Entry) ..];
            return @fromBackingInt(if (fmt == 0) at[0] else bytes.getU16BE(at));
        }

        pub fn getGlyphIndex(self: Self, codepoint: u32) GlyphIndex {
            if (codepoint > max_code or codepoint < self.first()) return .notdef;
            const idx = codepoint - self.first();
            if (idx >= self.count()) return .notdef;
            return self.getRawGlyphId(idx);
        }
    };
}

const setU16BE = bytes.setU16BE;
const setU32BE = bytes.setU32BE;

test "trimmed formats' headers match the stream layout (no padding)" {
    try zigstd.testExpectEqual(@as(usize, 6), @sizeOf(Layout(0)));
    try zigstd.testExpectEqual(@as(usize, 4), @offsetOf(Layout(0), "language"));
    try zigstd.testExpectEqual(@as(usize, 10), @sizeOf(Layout(6)));
    try zigstd.testExpectEqual(@as(usize, 6), @offsetOf(Layout(6), "first_code"));
    try zigstd.testExpectEqual(@as(usize, 8), @offsetOf(Layout(6), "entry_count"));
    try zigstd.testExpectEqual(@as(usize, 20), @sizeOf(Layout(10)));
    try zigstd.testExpectEqual(@as(usize, 12), @offsetOf(Layout(10), "start_char_code"));
    try zigstd.testExpectEqual(@as(usize, 16), @offsetOf(Layout(10), "num_chars"));
}

test "Format(0).getGlyphIndex indexes the byte array directly" {
    var buf: [262]u8 = @splat(0);
    buf[6 + 0x41] = 10; // 'A' -> glyph 10
    const fmt0 = try Format(0).init(buf[0..], true);
    try zigstd.testExpectEqual(@as(u16, 10), @backingInt(fmt0.getGlyphIndex(0x41)));
    try zigstd.testExpectEqual(GlyphIndex.notdef, fmt0.getGlyphIndex(0x42));
    try zigstd.testExpectEqual(GlyphIndex.notdef, fmt0.getGlyphIndex(0x100));
}

/// A format-6 (`wide` false) or format-10 subtable of `glyph_ids` from `first`.
fn build(alloc: zigstd.Allocator, wide: bool, first: u32, glyph_ids: []const u16) ![]u8 {
    const header: usize = if (wide) 20 else 10;
    const total = header + glyph_ids.len * 2;
    const buf = try alloc.alloc(u8, total);
    @memset(buf, 0);
    if (wide) {
        setU16BE(buf, 10);
        setU32BE(buf[4..], @intCast(total));
        setU32BE(buf[12..], first);
        setU32BE(buf[16..], @intCast(glyph_ids.len));
    } else {
        setU16BE(buf, 6);
        setU16BE(buf[2..], @intCast(total));
        setU16BE(buf[6..], @intCast(first));
        setU16BE(buf[8..], @intCast(glyph_ids.len));
    }
    for (glyph_ids, 0..) |v, i| setU16BE(buf[header + 2 * i ..], v);
    return buf;
}

test "Format(6).getGlyphIndex resolves codes within the trimmed range, .notdef outside" {
    const alloc = zigstd.testAllocator;
    const buf = try build(alloc, false, 100, &.{ 10, 11, 12 });
    defer alloc.free(buf);
    const fmt6 = try Format(6).init(buf, true);
    try zigstd.testExpectEqual(@as(u16, 10), @backingInt(fmt6.getGlyphIndex(100)));
    try zigstd.testExpectEqual(@as(u16, 12), @backingInt(fmt6.getGlyphIndex(102)));
    try zigstd.testExpectEqual(GlyphIndex.notdef, fmt6.getGlyphIndex(99));
    try zigstd.testExpectEqual(GlyphIndex.notdef, fmt6.getGlyphIndex(103));
    // past the 16-bit range this format can represent at all
    try zigstd.testExpectEqual(GlyphIndex.notdef, fmt6.getGlyphIndex(0x10000));
}

test "Format(10).getGlyphIndex resolves codes within the trimmed 32-bit range" {
    const alloc = zigstd.testAllocator;
    const buf = try build(alloc, true, 0x20000, &.{ 10, 11, 12 });
    defer alloc.free(buf);
    const fmt10 = try Format(10).init(buf, true);
    try zigstd.testExpectEqual(@as(u16, 10), @backingInt(fmt10.getGlyphIndex(0x20000)));
    try zigstd.testExpectEqual(@as(u16, 12), @backingInt(fmt10.getGlyphIndex(0x20002)));
    try zigstd.testExpectEqual(GlyphIndex.notdef, fmt10.getGlyphIndex(0x20003));
    try zigstd.testExpectEqual(GlyphIndex.notdef, fmt10.getGlyphIndex(0x1FFFF));
}
