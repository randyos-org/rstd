const types = @import("../types.zig");
const bytes = @import("../../../../bytes.zig");
const zigstd = @import("../../../../zigstd.zig");

/// The 'post' header.
const StreamLayout = extern struct {
    version: types.Version16Dot16,
    italic_angle: types.Fixed,
    underline_position: types.FWord,
    underline_thickness: types.FWord,
    is_fixed_pitch: u32,
    min_mem_type42: u32,
    max_mem_type42: u32,
    min_mem_type1: u32,
    max_mem_type1: u32,
};

const header_end = @sizeOf(StreamLayout);

/// PostScript information. Versions 2.0 and 2.5 additionally carry a
/// variable-length glyph-name section right after the fixed header (see
/// numGlyphs() and the accessors below it), which is why init() validates array
/// bounds based on version() rather than a single fixed size.
pub const PostTable = struct {
    const Self = @This();

    table_data: []const u8,

    pub fn init(table_data: []const u8, check: bool) !Self {
        if (check) try bytes.checkSize(table_data, 0, StreamLayout);
        const self: Self = .{ .table_data = table_data };

        if (check) {
            const version_bits: u32 = @bitCast(self.version());
            switch (version_bits) {
                0x00020000 => {
                    try bytes.checkSize(table_data, header_end, u16);
                    const num_glyphs = self.numGlyphs();
                    try bytes.checkArray(table_data, indicesStart(), num_glyphs, u16);
                    try bytes.checkOffset(table_data, stringDataStart(num_glyphs));
                },
                0x00025000 => {
                    try bytes.checkSize(table_data, header_end, u16);
                    try bytes.checkArray(table_data, indicesStart(), self.numGlyphs(), i8);
                },
                // 1.0 and 3.0 carry no glyph names at all: 1.0 implies the
                // standard Macintosh glyph set, 3.0 means none are provided.
                else => {},
            }
        }

        return self;
    }

    pub inline fn version(self: Self) types.Version16Dot16 {
        return types.Version16Dot16.fromBytes(self.table_data[@offsetOf(StreamLayout, "version")..]);
    }

    pub inline fn italicAngle(self: Self) types.Fixed {
        return types.Fixed.fromBytes(self.table_data[@offsetOf(StreamLayout, "italic_angle")..]);
    }

    pub inline fn underlinePosition(self: Self) types.FWord {
        return bytes.readIntBE(types.FWord, self.table_data[@offsetOf(StreamLayout, "underline_position")..]);
    }

    pub inline fn underlineThickness(self: Self) types.FWord {
        return bytes.readIntBE(types.FWord, self.table_data[@offsetOf(StreamLayout, "underline_thickness")..]);
    }

    pub inline fn isFixedPitch(self: Self) u32 {
        return bytes.getU32BE(self.table_data[@offsetOf(StreamLayout, "is_fixed_pitch")..]);
    }

    pub inline fn minMemType42(self: Self) u32 {
        return bytes.getU32BE(self.table_data[@offsetOf(StreamLayout, "min_mem_type42")..]);
    }

    pub inline fn maxMemType42(self: Self) u32 {
        return bytes.getU32BE(self.table_data[@offsetOf(StreamLayout, "max_mem_type42")..]);
    }

    pub inline fn minMemType1(self: Self) u32 {
        return bytes.getU32BE(self.table_data[@offsetOf(StreamLayout, "min_mem_type1")..]);
    }

    pub inline fn maxMemType1(self: Self) u32 {
        return bytes.getU32BE(self.table_data[@offsetOf(StreamLayout, "max_mem_type1")..]);
    }

    /// Valid only for version 2.0/2.5 (see version()); versions 1.0 and 3.0
    /// have no glyph-name section at all.
    pub inline fn numGlyphs(self: Self) u16 {
        return bytes.getU16BE(self.table_data[header_end..]);
    }

    inline fn indicesStart() usize {
        return header_end + @sizeOf(u16);
    }

    /// version 2.0 only: glyphNameIndex[index], see the spec for how to
    /// interpret it (either into the standard Macintosh glyph names, or
    /// offset by 258 into the Pascal strings stringData() returns).
    pub inline fn getGlyphNameIndex(self: Self, index: usize) u16 {
        return bytes.getU16BE(self.table_data[indicesStart() + index * @sizeOf(u16) ..]);
    }

    inline fn stringDataStart(num_glyphs: u16) usize {
        return indicesStart() + @as(usize, num_glyphs) * @sizeOf(u16);
    }

    /// version 2.0 only: raw Pascal-string-encoded glyph name storage.
    pub inline fn stringData(self: Self) []const u8 {
        return self.table_data[stringDataStart(self.numGlyphs())..];
    }

    /// version 2.5 (deprecated) only: per-glyph offset into the standard
    /// Macintosh glyph ordering.
    pub inline fn getMacGlyphOffset(self: Self, index: usize) i8 {
        return @bitCast(self.table_data[indicesStart() + index]);
    }
};

test "post StreamLayout matches the OpenType stream layout (32 bytes, no padding)" {
    try zigstd.testExpectEqual(@as(usize, 32), @sizeOf(StreamLayout));
    try zigstd.testExpectEqual(@as(usize, 4), @offsetOf(StreamLayout, "italic_angle"));
    try zigstd.testExpectEqual(@as(usize, 8), @offsetOf(StreamLayout, "underline_position"));
    try zigstd.testExpectEqual(@as(usize, 10), @offsetOf(StreamLayout, "underline_thickness"));
    try zigstd.testExpectEqual(@as(usize, 12), @offsetOf(StreamLayout, "is_fixed_pitch"));
    try zigstd.testExpectEqual(@as(usize, 28), @offsetOf(StreamLayout, "max_mem_type1"));
}

test "PostTable.init reads a version 3.0 header with no glyph names" {
    var buf: [32]u8 = @splat(0);
    bytes.setU32BE(buf[0..], 0x00030000); // version 3.0
    bytes.setU32BE(buf[12..], 1); // is_fixed_pitch

    const post = try PostTable.init(buf[0..], true);

    try zigstd.testExpectEqual(@as(u32, 1), post.isFixedPitch());
}

test "PostTable.init reads a version 2.0 table's glyph name indices and string data" {
    var buf: [32 + 2 + 4 + 3]u8 = @splat(0);
    bytes.setU32BE(buf[0..], 0x00020000); // version 2.0
    bytes.setU16BE(buf[32..], 2); // numGlyphs
    bytes.setU16BE(buf[34..], 65); // glyphNameIndex[0]: Macintosh standard name
    bytes.setU16BE(buf[36..], 258); // glyphNameIndex[1]: custom name 0
    buf[38] = 2; // Pascal string length
    buf[39] = 'A';
    buf[40] = 'B';

    const post = try PostTable.init(buf[0..], true);

    try zigstd.testExpectEqual(@as(u16, 65), post.getGlyphNameIndex(0));
    try zigstd.testExpectEqual(@as(u16, 258), post.getGlyphNameIndex(1));
    try zigstd.testExpectEqualSlices(u8, &[_]u8{ 2, 'A', 'B' }, post.stringData());
}

test "PostTable.init reads a version 2.5 table's per-glyph Macintosh offsets" {
    var buf: [32 + 2 + 3]u8 = @splat(0);
    bytes.setU32BE(buf[0..], 0x00025000); // version 2.5
    bytes.setU16BE(buf[32..], 3); // numGlyphs
    buf[34] = @bitCast(@as(i8, 36));
    buf[35] = @bitCast(@as(i8, 36));
    buf[36] = @bitCast(@as(i8, -1));

    const post = try PostTable.init(buf[0..], true);

    try zigstd.testExpectEqual(@as(i8, 36), post.getMacGlyphOffset(0));
    try zigstd.testExpectEqual(@as(i8, 36), post.getMacGlyphOffset(1));
    try zigstd.testExpectEqual(@as(i8, -1), post.getMacGlyphOffset(2));
}
