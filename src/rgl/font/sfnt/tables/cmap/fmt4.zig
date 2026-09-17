const types = @import("../../stream_types.zig");
const bytes = @import("../../../../../bytes.zig");
const zigstd = @import("../../../../../zigstd.zig");

const GlyphIndex = types.GlyphIndex;
const setU16BE = bytes.setU16BE;
const setI16BE = bytes.setI16BE;

/// Schema for the format-4 header, used only to compute field offsets/sizes
/// at comptime; never instantiated.
const StreamLayout = extern struct {
    format: u16,
    length: u16,
    language: u16,
    seg_count_x2: u16,
    search_range: u16,
    entry_selector: u16,
    range_shift: u16,
};

const header_end = @sizeOf(StreamLayout);

pub const Format4 = struct {
    const Self = @This();

    table_data: []const u8,
    seg_count: usize,

    pub fn init(table_data: []const u8, check: bool) !Self {
        if (check) try bytes.checkSize(table_data, 0, StreamLayout);
        const seg_count_x2 = bytes.getU16BE(table_data[@offsetOf(StreamLayout, "seg_count_x2")..]);
        const seg_count: usize = seg_count_x2 >> 1;
        const self: Self = .{ .table_data = table_data, .seg_count = seg_count };

        if (check) {
            try bytes.checkArray(table_data, header_end, seg_count, u16); // endCode
            try bytes.checkSize(table_data, self.endCodeEnd(), u16); // reservedPad
            try bytes.checkArray(table_data, self.startCodeStart(), seg_count, u16);
            try bytes.checkArray(table_data, self.idDeltaStart(), seg_count, i16);
            try bytes.checkArray(table_data, self.idRangeOffsetStart(), seg_count, u16);
            try bytes.checkOffset(table_data, self.glyphArrayStart());
        }

        return self;
    }

    /// Byte length of one of the four seg_count-sized parallel arrays
    /// (endCode/startCode/idDelta/idRangeOffset) - all 2 bytes per entry
    /// regardless of signedness, so one helper covers all four.
    inline fn segArrayBytes(self: Self) usize {
        return self.seg_count * @sizeOf(u16);
    }

    inline fn endCodeEnd(self: Self) usize {
        return header_end + self.segArrayBytes();
    }

    inline fn startCodeStart(self: Self) usize {
        return self.endCodeEnd() + @sizeOf(u16); // reservedPad
    }

    inline fn idDeltaStart(self: Self) usize {
        return self.startCodeStart() + self.segArrayBytes();
    }

    inline fn idRangeOffsetStart(self: Self) usize {
        return self.idDeltaStart() + self.segArrayBytes();
    }

    inline fn glyphArrayStart(self: Self) usize {
        return self.idRangeOffsetStart() + self.segArrayBytes();
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

    pub inline fn segCountX2(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "seg_count_x2")..]);
    }

    /// Microsoft recommends deriving the binary-search parameters (this one
    /// included) rather than trusting these header fields, since some fonts
    /// carry stale values; getGlyphIndex() always recomputes its own. Kept
    /// here only so callers can inspect what the font itself claims.
    pub inline fn searchRange(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "search_range")..]);
    }

    pub inline fn entrySelector(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "entry_selector")..]);
    }

    pub inline fn rangeShift(self: Self) u16 {
        return bytes.getU16BE(self.table_data[@offsetOf(StreamLayout, "range_shift")..]);
    }

    pub fn getEndCode(self: Self, index: usize) u16 {
        return bytes.getU16BE(self.table_data[header_end + index * @sizeOf(u16) ..]);
    }

    pub fn getStartCode(self: Self, index: usize) u16 {
        return bytes.getU16BE(self.table_data[self.startCodeStart() + index * @sizeOf(u16) ..]);
    }

    pub fn getIdDelta(self: Self, index: usize) i16 {
        return bytes.readIntBE(i16, self.table_data[self.idDeltaStart() + index * @sizeOf(i16) ..]);
    }

    pub fn getIdRangeOffset(self: Self, index: usize) u16 {
        return bytes.getU16BE(self.table_data[self.idRangeOffsetStart() + index * @sizeOf(u16) ..]);
    }

    pub fn getGlyphIndex(self: Self, codepoint: u32) GlyphIndex {
        if (codepoint > 0xffff) return .notdef;
        // Narrowed once, up front: idDelta arithmetic below is modulo 65536
        // per spec, so it must wrap at exactly 16 bits, not 32.
        const c: u16 = @intCast(codepoint);

        const seg_count = self.seg_count;

        // Microsoft recommends not using the search values in the header and
        // deriving them ourselves.  Because we do want to run on some embedded
        // systems with limited resources, we will still implement the binary search.
        const entry_selector: u16 = zigstd.log2Int(u16, @intCast(seg_count));
        // note: these are halved
        const search_range: u16 = @as(u16, 1) << @intCast(entry_selector);
        const range_shift: u16 = @as(u16, @intCast(seg_count)) - search_range;

        // If c is past the last of the "remainder" segments, the answer lives
        // in the pow2-sized block starting right after them; otherwise it's
        // somewhere in [0, pow2), which is the same size as that block and
        // still fully halvable in entry_selector steps.
        var i: usize = 0;
        if (range_shift > 0 and self.getEndCode(range_shift - 1) < c) {
            i = range_shift;
        }

        var range = search_range;
        var steps = entry_selector;
        while (steps > 0) : (steps -= 1) {
            range >>= 1;
            if (self.getEndCode(i + range - 1) < c) {
                i += range;
            }
        }

        if (self.getStartCode(i) > c) return .notdef;

        const byte_offset_from_id_offset = self.getIdRangeOffset(i);
        // idDelta arithmetic is modulo 65536 per spec: reinterpret the i16 as
        // u16 bits and use wrapping add, which is exactly mod-65536 addition.
        const delta: u16 = @bitCast(self.getIdDelta(i));
        if (byte_offset_from_id_offset == 0) {
            return @fromBackingInt(c +% delta);
        }

        // The spec's own C pointer trick applies directly here: idRangeOffset[i]
        // is a byte offset measured from the location of the idRangeOffset[i]
        // field itself, and glyphIdArray immediately follows idRangeOffset[]
        // in the actual stream bytes, since we're reading straight from
        // table_data rather than separately-allocated arrays.
        const id_range_offset_field_pos = self.idRangeOffsetStart() + i * @sizeOf(u16);
        const target = id_range_offset_field_pos + @as(usize, byte_offset_from_id_offset) + @as(usize, c - self.getStartCode(i)) * @sizeOf(GlyphIndex);
        if (target + @sizeOf(GlyphIndex) > self.table_data.len) return .notdef;

        const raw: GlyphIndex = @fromBackingInt(@intCast(bytes.getU16BE(self.table_data[target..])));
        if (raw == .notdef) return .notdef;
        return @fromBackingInt(@backingInt(raw) +% delta);
    }
};

// Builds a format-4 subtable byte buffer from parallel segment arrays plus
// an optional glyphIdArray tail, matching the stream layout exactly.
fn buildFormat4(
    alloc: zigstd.Allocator,
    end_code: []const u16,
    start_code: []const u16,
    id_delta: []const i16,
    id_range_offset: []const u16,
    glyph_id_array: []const u16,
) ![]u8 {
    const seg_count = end_code.len;
    const reserved_pad_size = @sizeOf(u16);
    // endCode + startCode + idRangeOffset (u16) + idDelta (i16) + reservedPad
    const seg_arrays_size = seg_count * @sizeOf(u16) * 3 + seg_count * @sizeOf(i16) + reserved_pad_size;
    const total = header_end + seg_arrays_size + glyph_id_array.len * @sizeOf(GlyphIndex);

    const buf = try alloc.alloc(u8, total);
    @memset(buf, 0);

    setU16BE(buf, 4); // format
    setU16BE(buf[@offsetOf(StreamLayout, "length")..], @intCast(total));
    setU16BE(buf[@offsetOf(StreamLayout, "language")..], 0);
    setU16BE(buf[@offsetOf(StreamLayout, "seg_count_x2")..], @intCast(seg_count * 2));
    setU16BE(buf[@offsetOf(StreamLayout, "search_range")..], 0); // unused by our parser, per the spec's own recommendation
    setU16BE(buf[@offsetOf(StreamLayout, "entry_selector")..], 0);
    setU16BE(buf[@offsetOf(StreamLayout, "range_shift")..], 0);

    var off: usize = header_end;
    for (end_code) |v| {
        setU16BE(buf[off..], v);
        off += @sizeOf(u16);
    }
    off += reserved_pad_size; // left 0
    for (start_code) |v| {
        setU16BE(buf[off..], v);
        off += @sizeOf(u16);
    }
    for (id_delta) |v| {
        setI16BE(buf[off..], v);
        off += @sizeOf(i16);
    }
    for (id_range_offset) |v| {
        setU16BE(buf[off..], v);
        off += @sizeOf(u16);
    }
    for (glyph_id_array) |v| {
        setU16BE(buf[off..], v);
        off += @sizeOf(GlyphIndex);
    }

    return buf;
}

// The OpenType spec's own worked example for format 4 (see "cmap Subtable
// Format 4" -> the segCountX2=8 table): segments cover 10-20, 30-90, 153-480,
// and the required trailing 0xFFFF sentinel segment. All idRangeOffset are 0,
// so every mapping resolves via idDelta alone.
test "Format4.getGlyphIndex matches the OpenType spec's own worked example" {
    const alloc = zigstd.testAllocator;
    const buf = try buildFormat4(
        alloc,
        &[_]u16{ 20, 90, 480, 0xFFFF },
        &[_]u16{ 10, 30, 153, 0xFFFF },
        &[_]i16{ -9, -18, -80, 1 },
        &[_]u16{ 0, 0, 0, 0 },
        &[_]u16{},
    );
    defer alloc.free(buf);

    const fmt4 = try Format4.init(buf, true);

    try zigstd.testExpectEqual(@as(u16, 1), @backingInt(fmt4.getGlyphIndex(10)));
    try zigstd.testExpectEqual(@as(u16, 11), @backingInt(fmt4.getGlyphIndex(20)));
    try zigstd.testExpectEqual(@as(u16, 12), @backingInt(fmt4.getGlyphIndex(30)));
    try zigstd.testExpectEqual(@as(u16, 72), @backingInt(fmt4.getGlyphIndex(90)));
    try zigstd.testExpectEqual(@as(u16, 73), @backingInt(fmt4.getGlyphIndex(153)));
    try zigstd.testExpectEqual(@as(u16, 400), @backingInt(fmt4.getGlyphIndex(480)));
    try zigstd.testExpectEqual(@as(u16, 0), @backingInt(fmt4.getGlyphIndex(0xFFFF)));
    // past the 16-bit range this format can represent at all
    try zigstd.testExpectEqual(GlyphIndex.notdef, fmt4.getGlyphIndex(0x10000));
}

test "Format4.getGlyphIndex returns .notdef for a character in a gap between segments" {
    const alloc = zigstd.testAllocator;
    const buf = try buildFormat4(
        alloc,
        &[_]u16{ 20, 90, 480, 0xFFFF },
        &[_]u16{ 10, 30, 153, 0xFFFF },
        &[_]i16{ -9, -18, -80, 1 },
        &[_]u16{ 0, 0, 0, 0 },
        &[_]u16{},
    );
    defer alloc.free(buf);

    const fmt4 = try Format4.init(buf, true);

    // 25 is between segment 0 (ends at 20) and segment 1 (starts at 30).
    try zigstd.testExpectEqual(GlyphIndex.notdef, fmt4.getGlyphIndex(25));
    // 500 is past the last real segment (480) but before the sentinel's own
    // start (0xFFFF), so it also falls in a gap.
    try zigstd.testExpectEqual(GlyphIndex.notdef, fmt4.getGlyphIndex(500));
}

// Exercises the idRangeOffset/glyphIdArray path, which the spec's own
// example (all-zero idRangeOffset) never touches. Single segment [100,102],
// idRangeOffset points straight at glyph_id_array[0..3].
test "Format4.getGlyphIndex resolves through idRangeOffset into glyph_id_array" {
    const alloc = zigstd.testAllocator;
    // idRangeOffset field sits at byte offset (header=14) + segCount*2 (endCode)
    // + 2 (reservedPad) + segCount*2 (startCode) + segCount*2 (idDelta) = for
    // segCount=1: 14+2+2+2+2 = 22. glyphIdArray starts right after it, at 24.
    // idRangeOffset[0] must equal the byte distance from ITS OWN location (22)
    // to glyphIdArray[0] (24), i.e. 2.
    const buf = try buildFormat4(
        alloc,
        &[_]u16{102},
        &[_]u16{100},
        &[_]i16{5},
        &[_]u16{2},
        &[_]u16{ 200, 0, 300 }, // codes 100,101,102 -> raw 200,0(missing),300
    );
    defer alloc.free(buf);

    const fmt4 = try Format4.init(buf, true);

    try zigstd.testExpectEqual(@as(u16, 205), @backingInt(fmt4.getGlyphIndex(100))); // 200 + idDelta 5
    try zigstd.testExpectEqual(GlyphIndex.notdef, fmt4.getGlyphIndex(101)); // raw 0 -> .notdef regardless of idDelta
    try zigstd.testExpectEqual(@as(u16, 305), @backingInt(fmt4.getGlyphIndex(102))); // 300 + idDelta 5
}

// Cross-checks the search itself against a large, deliberately non-power-of-
// two number of segments (entry_selector/pow2/remainder all non-trivial),
// comparing every possible c against an obviously-correct linear scan.
test "Format4.getGlyphIndex matches a linear scan for a non-power-of-two segment count" {
    const alloc = zigstd.testAllocator;
    const seg_count = 37;
    var end_code: [seg_count]u16 = undefined;
    var start_code: [seg_count]u16 = undefined;
    var id_delta: [seg_count]i16 = undefined;
    for (0..seg_count) |i| {
        start_code[i] = @intCast(i * 20);
        end_code[i] = @intCast(i * 20 + 10);
        id_delta[i] = 0;
    }
    end_code[seg_count - 1] = 0xFFFF;
    const id_range_offset: [seg_count]u16 = @splat(0);

    const buf = try buildFormat4(alloc, &end_code, &start_code, &id_delta, &id_range_offset, &[_]u16{});
    defer alloc.free(buf);

    const fmt4 = try Format4.init(buf, true);

    var c: u32 = 0;
    while (c <= 0xFFFF) : (c += 1) {
        var expected: GlyphIndex = .notdef;
        for (0..seg_count) |i| {
            if (end_code[i] >= c and start_code[i] <= @as(u16, @intCast(c))) {
                expected = @fromBackingInt(@as(u16, @intCast(c)));
                break;
            }
        }
        try zigstd.testExpectEqual(expected, fmt4.getGlyphIndex(@intCast(c)));
    }
}

test "Format4 stream layout matches the OpenType stream layout (14 bytes, no padding)" {
    try zigstd.testExpectEqual(@as(usize, 14), @sizeOf(StreamLayout));
    try zigstd.testExpectEqual(@as(usize, 6), @offsetOf(StreamLayout, "seg_count_x2"));
    try zigstd.testExpectEqual(@as(usize, 8), @offsetOf(StreamLayout, "search_range"));
    try zigstd.testExpectEqual(@as(usize, 10), @offsetOf(StreamLayout, "entry_selector"));
    try zigstd.testExpectEqual(@as(usize, 12), @offsetOf(StreamLayout, "range_shift"));
}
