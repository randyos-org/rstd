const types = @import("../types.zig");
const bytes = @import("../../../../bytes.zig");
const zigstd = @import("../../../../zigstd.zig");

const Allocator = zigstd.Allocator;
const GlyphIndex = types.GlyphIndex;

// Index to location. Not self-describing: its entry count comes from
// maxp.num_glyphs, and whether entries are Offset16 (halved) or Offset32
// depends on head.index_to_loc_format, so init() takes both rather than
// deriving them from table_data the way other tables do. Offsets are read
// fresh from table_data on each getOffset()/getGlyphRange() call rather than
// being decoded into a stored array up front.
pub const LocaTable = struct {
    const Self = @This();

    table_data: []const u8,
    num_glyphs: u16,
    long_format: bool,

    pub const GlyphRange = struct {
        start: u32,
        end: u32,
    };

    pub fn init(table_data: []const u8, num_glyphs: u16, long_format: bool, check: bool) !Self {
        if (check) {
            const count = @as(usize, num_glyphs) + 1;
            if (long_format) {
                try bytes.checkArray(table_data, 0, count, types.Offset32);
            } else {
                try bytes.checkArray(table_data, 0, count, u16);
            }
        }

        return .{ .table_data = table_data, .num_glyphs = num_glyphs, .long_format = long_format };
    }

    /// Always returns an actual byte offset into 'glyf', regardless of
    /// whether the stream stores it directly (long format) or halved (short
    /// format).
    pub inline fn getOffset(self: Self, index: usize) u32 {
        if (self.long_format) return bytes.getU32BE(self.table_data[index * 4 ..]);
        return @as(u32, bytes.getU16BE(self.table_data[index * 2 ..])) * 2;
    }

    /// The byte range of glyph's outline data within 'glyf'. Null means the
    /// glyph has no outline (e.g. space) or the index is out of range.
    pub fn getGlyphRange(self: Self, glyph: GlyphIndex) ?GlyphRange {
        const glyph_index: usize = @backingInt(glyph);
        const count = @as(usize, self.num_glyphs) + 1;
        if (glyph_index + 1 >= count) return null;

        const start = self.getOffset(glyph_index);
        const end = self.getOffset(glyph_index + 1);
        if (end <= start) return null;

        return .{ .start = start, .end = end };
    }
};

fn buildLocaLong(alloc: Allocator, offsets: []const u32) ![]u8 {
    const buf = try alloc.alloc(u8, offsets.len * 4);
    for (offsets, 0..) |v, i| bytes.setU32BE(buf[i * 4 ..], v);
    return buf;
}

fn buildLocaShort(alloc: Allocator, halved_offsets: []const u16) ![]u8 {
    const buf = try alloc.alloc(u8, halved_offsets.len * 2);
    for (halved_offsets, 0..) |v, i| bytes.setU16BE(buf[i * 2 ..], v);
    return buf;
}

test "LocaTable.getOffset (long format) exposes offsets directly" {
    const alloc = zigstd.testAllocator;
    const buf = try buildLocaLong(alloc, &[_]u32{ 0, 40, 40, 120 });
    defer alloc.free(buf);

    const loca = try LocaTable.init(buf, 3, true, true);

    const expected = [_]u32{ 0, 40, 40, 120 };
    for (expected, 0..) |v, i| try zigstd.testExpectEqual(v, loca.getOffset(i));
}

test "LocaTable.getOffset (short format) doubles the stored halved offsets" {
    const alloc = zigstd.testAllocator;
    const buf = try buildLocaShort(alloc, &[_]u16{ 0, 20, 20, 60 });
    defer alloc.free(buf);

    const loca = try LocaTable.init(buf, 3, false, true);

    const expected = [_]u32{ 0, 40, 40, 120 };
    for (expected, 0..) |v, i| try zigstd.testExpectEqual(v, loca.getOffset(i));
}

test "LocaTable.getGlyphRange returns the byte span between consecutive offsets" {
    const alloc = zigstd.testAllocator;
    const buf = try buildLocaLong(alloc, &[_]u32{ 0, 40, 40, 120 });
    defer alloc.free(buf);

    const loca = try LocaTable.init(buf, 3, true, true);

    const r0 = loca.getGlyphRange(@fromBackingInt(@intCast(0))).?;
    try zigstd.testExpectEqual(@as(u32, 0), r0.start);
    try zigstd.testExpectEqual(@as(u32, 40), r0.end);

    // glyph 1 has no outline: loca[1] == loca[2]
    try zigstd.testExpectEqual(@as(?LocaTable.GlyphRange, null), loca.getGlyphRange(@fromBackingInt(@intCast(1))));

    const r2 = loca.getGlyphRange(@fromBackingInt(@intCast(2))).?;
    try zigstd.testExpectEqual(@as(u32, 40), r2.start);
    try zigstd.testExpectEqual(@as(u32, 120), r2.end);
}

test "LocaTable.getGlyphRange returns null for an out-of-range glyph" {
    const alloc = zigstd.testAllocator;
    const buf = try buildLocaLong(alloc, &[_]u32{ 0, 40 });
    defer alloc.free(buf);

    const loca = try LocaTable.init(buf, 1, true, true);

    try zigstd.testExpectEqual(@as(?LocaTable.GlyphRange, null), loca.getGlyphRange(@fromBackingInt(@intCast(1))));
}
