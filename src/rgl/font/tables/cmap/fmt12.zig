const types = @import("../../stream_types.zig");
const bytes = @import("../../../../bytes.zig");
const zigstd = @import("../../../../zigstd.zig");

const GlyphIndex = types.GlyphIndex;
const setU16BE = bytes.setU16BE;
const setU32BE = bytes.setU32BE;

/// Schema for the format-12 header, used only to compute field offsets/sizes
/// at comptime; never instantiated.
const StreamLayout = extern struct {
    format: u16,
    reserved: u16,
    length: u32,
    language: u32,
    num_groups: u32,
};

const groups_start = @sizeOf(StreamLayout);

/// One sequential map group. Its extern layout doubles as the wire schema
/// for @offsetOf/@sizeOf.
pub const Group = extern struct {
    start_char_code: u32,
    end_char_code: u32,
    start_glyph_id: u32,
};

// Segmented coverage: format 4's segment-plus-delta trick, extended to 32-bit
// codes and made sparse-friendly (groups need not be contiguous).
pub const Format12 = struct {
    const Self = @This();

    table_data: []const u8,

    pub fn init(table_data: []const u8, check: bool) !Self {
        if (check) try bytes.checkSize(table_data, 0, StreamLayout);
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

    pub inline fn numGroups(self: Self) u32 {
        return bytes.getU32BE(self.table_data[@offsetOf(StreamLayout, "num_groups")..]);
    }

    pub fn getGroup(self: Self, index: usize) Group {
        return bytes.typeBEToNative(Group, self.table_data[groups_start + index * @sizeOf(Group) ..]);
    }

    // Groups are sorted and non-overlapping per spec, so binary search applies.
    pub fn getGlyphIndex(self: Self, codepoint: u32) GlyphIndex {
        var lo: usize = 0;
        var hi: usize = self.numGroups();
        while (lo < hi) {
            const mid = lo + (hi - lo) / 2;
            const group = self.getGroup(mid);
            if (codepoint < group.start_char_code) {
                hi = mid;
            } else if (codepoint > group.end_char_code) {
                lo = mid + 1;
            } else {
                const raw = group.start_glyph_id + (codepoint - group.start_char_code);
                return @fromBackingInt(@as(u16, @truncate(raw)));
            }
        }
        return .notdef;
    }
};

test "Format12 header and group match the stream layout (no padding)" {
    try zigstd.testExpectEqual(@as(usize, 16), @sizeOf(StreamLayout));
    try zigstd.testExpectEqual(@as(usize, 12), @offsetOf(StreamLayout, "num_groups"));
    try zigstd.testExpectEqual(@as(usize, 12), @sizeOf(Group));
    try zigstd.testExpectEqual(@as(usize, 8), @offsetOf(Group, "start_glyph_id"));
}

fn buildFormat12(alloc: zigstd.Allocator, groups: []const Group) ![]u8 {
    const header_size = 16;
    const total = header_size + groups.len * 12;

    const buf = try alloc.alloc(u8, total);
    @memset(buf, 0);

    setU16BE(buf, 12); // format
    setU16BE(buf[2..], 0); // reserved
    setU32BE(buf[4..], @intCast(total)); // length
    setU32BE(buf[8..], 0); // language
    setU32BE(buf[12..], @intCast(groups.len));

    var off: usize = header_size;
    for (groups) |g| {
        setU32BE(buf[off..], g.start_char_code);
        setU32BE(buf[off + 4 ..], g.end_char_code);
        setU32BE(buf[off + 8 ..], g.start_glyph_id);
        off += 12;
    }

    return buf;
}

test "Format12.getGlyphIndex resolves through sparse, non-contiguous groups" {
    const alloc = zigstd.testAllocator;
    const buf = try buildFormat12(alloc, &[_]Group{
        .{ .start_char_code = 0x41, .end_char_code = 0x5A, .start_glyph_id = 10 },
        .{ .start_char_code = 0x10000, .end_char_code = 0x10010, .start_glyph_id = 500 },
    });
    defer alloc.free(buf);

    const fmt12 = try Format12.init(buf, true);

    try zigstd.testExpectEqual(@as(u16, 10), @backingInt(fmt12.getGlyphIndex(0x41)));
    try zigstd.testExpectEqual(@as(u16, 35), @backingInt(fmt12.getGlyphIndex(0x5A)));
    try zigstd.testExpectEqual(@as(u16, 505), @backingInt(fmt12.getGlyphIndex(0x10005)));
    try zigstd.testExpectEqual(GlyphIndex.notdef, fmt12.getGlyphIndex(0x5B));
    try zigstd.testExpectEqual(GlyphIndex.notdef, fmt12.getGlyphIndex(0xFFFF));
}
