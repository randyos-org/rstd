const types = @import("../types.zig");
const bytes = @import("../../../bytes.zig");
const zigstd = @import("../../../zigstd.zig");
const reduce = @import("../../../reduce.zig");
const Sfnt = @import("../sfnt.zig");

const Allocator = zigstd.Allocator;
const GlyphIndex = types.GlyphIndex;

const PlatformId = enum(u16) {
    unicode = 0,
    mac = 1,
    iso = 2,
    microsoft = 3,
    custom = 4,
};

const UnicodeEncodingId = enum(u16) {
    version_1_0 = 0,
    version_1_1 = 1,
    iso_10646 = 2,
    unicode_2_bmp = 3,
    unicode_2_full = 4,
    unicode_2_var = 5,
    unicode_full = 6,
};

const IsoEncodingId = enum(u16) {
    ascii = 0,
    iso10646 = 1,
    iso8859_1 = 2,
};

const MicrosoftEncodingId = enum(u16) {
    symbol = 0,
    unicode_bmp = 1,
    shiftjis = 2,
    prc = 3,
    big5 = 4,
    wansung = 5,
    johab = 6,
    reserved1 = 7,
    reserved2 = 8,
    reserved3 = 9,
    unicode_full = 10,
};

const EncodingId = union(enum) {
    unicode: UnicodeEncodingId,
    iso: IsoEncodingId,
    microsoft: MicrosoftEncodingId,
};

const CmapTableHeader = extern struct {
    version: u16,
    num_tables: u16,
};

const CmapEncRecord = struct {
    const Self = @This();

    platform_id: PlatformId,
    encoding_id: EncodingId,
    cmap_data: CmapSubtable,

    pub fn isUnicodeBmp(self: Self) bool {
        return (self.platform_id == .unicode and self.encoding_id.unicode == .unicode_2_bmp) or
            (self.platform_id == .microsoft and self.encoding_id.microsoft == .unicode_bmp);
    }
};

const CmapSubtable = union(enum) {
    const Self = @This();

    fmt0: *Format0,
    fmt2: *Format2,
    fmt4: *Format4,

    /// raw format number for formats not yet implemented
    unsupported: u16,

    pub fn create(alloc: Allocator, cmap_data: []const u8, cmap_offset: u32, check: bool) !Self {
        const table_data = cmap_data[cmap_offset..];
        const format = bytes.readIntBE(u16, table_data);
        return switch (format) {
            0 => .{ .fmt0 = try Format0.create(alloc, table_data, check) },
            2 => .{ .fmt2 = try Format2.create(alloc, table_data, check) },
            4 => .{ .fmt4 = try Format4.create(alloc, table_data, check) },
            else => .{ .unsupported = format },
        };
    }

    pub fn destroy(self: *Self, alloc: Allocator) void {
        switch (self.*) {
            .fmt0 => |fmt0| fmt0.destroy(alloc),
            .fmt2 => |fmt2| fmt2.destroy(alloc),
            .fmt4 => |fmt4| fmt4.destroy(alloc),
            .unsupported => {},
        }
    }
};

const RawCmapEncRecord = extern struct {
    const Self = @This();

    platform_id: PlatformId,
    encoding_id: u16,
    cmap_offset: u32,

    pub fn getEncodingEnum(self: Self) ?EncodingId {
        return switch (self.platform_id) {
            .unicode => EncodingId{ .unicode = @fromBackingInt(self.encoding_id) },
            .iso => EncodingId{ .iso = @fromBackingInt(self.encoding_id) },
            .microsoft => EncodingId{ .microsoft = @fromBackingInt(self.encoding_id) },
            else => null,
        };
    }
};

pub const CmapTable = struct {
    const Self = @This();

    header: CmapTableHeader,
    enc_records: []CmapEncRecord,

    pub fn destroy(self: *Self, alloc: Allocator) void {
        for (self.enc_records) |*enc_record| enc_record.cmap_data.destroy(alloc);
        alloc.free(self.enc_records);
        alloc.destroy(self);
    }

    pub fn create(alloc: Allocator, cmap_bytes: []const u8, check: bool) !*Self {
        var self = try alloc.create(Self);
        errdefer alloc.destroy(self);

        if (check) try bytes.checkSize(cmap_bytes, 0, CmapTableHeader);
        self.header = bytes.typeBEToNative(CmapTableHeader, cmap_bytes);

        const tables_start = @sizeOf(CmapTableHeader);
        const num_tables = self.header.num_tables;
        if (check) try bytes.checkArray(cmap_bytes, tables_start, num_tables, RawCmapEncRecord);
        const raw_enc_records = try bytes.allocArrayBEToNative(alloc, RawCmapEncRecord, num_tables, cmap_bytes[tables_start..]);
        defer alloc.free(raw_enc_records);

        self.enc_records = try alloc.alloc(CmapEncRecord, num_tables);
        errdefer alloc.free(self.enc_records);

        var created: usize = 0;
        errdefer for (self.enc_records[0..created]) |*enc_record| enc_record.cmap_data.destroy(alloc);
        for (raw_enc_records, 0..) |raw, i| {
            self.enc_records[i] = .{
                .platform_id = raw.platform_id,
                .encoding_id = raw.getEncodingEnum().?,
                .cmap_data = try .create(alloc, cmap_bytes, raw.cmap_offset, check),
            };
            created += 1;
        }

        return self;
    }

    // fn readIndex(self: CmapTable) Index {
    //     return fixEndianness(std.mem.bytesToValue(Index, self.cmap_bytes[0 .. @bitSizeOf(Index) / 8]));
    // }

    // fn readSubtableLookup(self: CmapTable, idx: usize) SubtableLookup {
    //     const subtable_size = @bitSizeOf(SubtableLookup) / 8;
    //     const start = @bitSizeOf(Index) / 8 + idx * subtable_size;
    //     const end = start + subtable_size;

    //     return fixEndianness(std.mem.bytesToValue(SubtableLookup, self.cmap_bytes[start..end]));
    // }

    // fn readSubtableFormat(self: CmapTable, offset: usize) u16 {
    //     return fixEndianness(std.mem.bytesToValue(u16, self.cmap_bytes[offset .. offset + 2]));
    // }

    // fn readSubtableFormat4(self: CmapTable, alloc: Allocator, offset: usize) !SubtableFormat4 {
    //     var runtime_parser = RuntimeParser{ .data = self.cmap_bytes[offset..] };
    //     const format = runtime_parser.readVal(u16);
    //     const length = runtime_parser.readVal(u16);
    //     const language = runtime_parser.readVal(u16);
    //     const seg_count_x2 = runtime_parser.readVal(u16);
    //     const search_range = runtime_parser.readVal(u16);
    //     const entry_selector = runtime_parser.readVal(u16);
    //     const range_shift = runtime_parser.readVal(u16);

    //     const end_code: []const u16 = try runtime_parser.readArray(u16, alloc, seg_count_x2 / 2);
    //     const reserved_pad = runtime_parser.readVal(u16);
    //     const start_code: []const u16 = try runtime_parser.readArray(u16, alloc, seg_count_x2 / 2);
    //     const id_delta: []const u16 = try runtime_parser.readArray(u16, alloc, seg_count_x2 / 2);
    //     const id_range_offset: []const u16 = try runtime_parser.readArray(u16, alloc, seg_count_x2 / 2);
    //     const glyph_indices: []const u16 = try runtime_parser.readArray(u16, alloc, (runtime_parser.data.len - runtime_parser.idx) / 2);

    //     return .{
    //         .format = format,
    //         .length = length,
    //         .language = language,
    //         .seg_count_x2 = seg_count_x2,
    //         .search_range = search_range,
    //         .entry_selector = entry_selector,
    //         .range_shift = range_shift,
    //         .end_code = end_code,
    //         .reserved_pad = reserved_pad,
    //         .start_code = start_code,
    //         .id_delta = id_delta,
    //         .id_range_offset = id_range_offset,
    //         .glyph_indices = glyph_indices,
    //     };
    // }
};

const Format0 = extern struct {
    const Self = @This();

    format: u16,
    length: u16,
    language: u16,
    glyph_id_array: [256]u8, // is this always 256 bytes?

    pub fn create(alloc: Allocator, table_data: []const u8, check: bool) !*Self {
        if (check) try bytes.checkSize(table_data, 0, Self);
        return bytes.allocTypeBEToNative(alloc, Self, table_data);
    }

    pub fn destroy(self: *Self, alloc: Allocator) void {
        alloc.destroy(self);
    }

    pub fn getGlyphIndex(self: *Self, codepoint: u8) GlyphIndex {
        return @fromBackingInt(self.glyph_id_array[codepoint]);
    }
};

const Format2 = struct {
    const Self = @This();

    header: Format2Header,
    sub_headers: []Format2SubHeader,
    glyph_id_array: []GlyphIndex,

    const Format2Header = extern struct {
        format: u16,
        length: u16,
        language: u16,
        sub_header_keys: [256]u16,
    };

    const Format2SubHeader = extern struct {
        first_code: u16,
        entry_count: u16,
        id_delta: i16,
        id_range_offset: u16,
    };

    pub fn create(alloc: Allocator, table_data: []const u8, check: bool) !*Self {
        if (check) try bytes.checkSize(table_data, 0, Format2Header);
        const self = try alloc.create(Self);
        errdefer alloc.destroy(self);
        self.header = bytes.typeBEToNative(Format2Header, table_data);

        // sub_header_keys[i] stores (subheader index * 8) -- the byte offset of
        // that subheader within the subHeaders array. Nothing in the spec gives
        // an explicit subheader count; every implementation infers it as the
        // highest referenced offset (divided by the 8-byte record size), plus
        // one for the subheader that lives at that offset.
        const max_key = reduce.arrayMax(self.header.sub_header_keys);
        const num_sub_headers = @as(usize, max_key) / 8 + 1;

        const header_end = @sizeOf(Format2Header);
        if (check) try bytes.checkArray(table_data, header_end, num_sub_headers, Format2SubHeader);
        self.sub_headers = try bytes.allocArrayBEToNative(alloc, Format2SubHeader, num_sub_headers, table_data[header_end..]);
        errdefer alloc.free(self.sub_headers);

        // glyphIdArray fills the rest of the subtable.
        // lookups index the remaining bytes of this table via each subheader's id_range_offset
        const glyph_array_start = header_end + num_sub_headers * @sizeOf(Format2SubHeader);
        if (check) try bytes.checkOffset(table_data, glyph_array_start);
        const glyph_array_len = (table_data.len - glyph_array_start) / @sizeOf(GlyphIndex);
        self.glyph_id_array = try bytes.allocArrayBEToNative(alloc, GlyphIndex, glyph_array_len, table_data[glyph_array_start..]);

        return self;
    }

    pub fn destroy(self: *Self, alloc: Allocator) void {
        alloc.free(self.glyph_id_array);
        alloc.free(self.sub_headers);
        alloc.destroy(self);
    }
};

// These lock in the OpenType 'cmap' wire layout: every struct here is read
// directly from file bytes via bytes.zig's memcpy+byteswap machinery, so its
// Zig-computed size and field offsets must exactly match the spec's byte
// layout with no compiler-inserted padding anywhere.
test "CmapTableHeader matches the wire layout (4 bytes, no padding)" {
    try zigstd.testExpectEqual(@as(usize, 4), @sizeOf(CmapTableHeader));
    try zigstd.testExpectEqual(@as(usize, 0), @offsetOf(CmapTableHeader, "version"));
    try zigstd.testExpectEqual(@as(usize, 2), @offsetOf(CmapTableHeader, "num_tables"));
}

test "RawCmapEncRecord matches the wire layout (8 bytes, no padding)" {
    try zigstd.testExpectEqual(@as(usize, 8), @sizeOf(RawCmapEncRecord));
    try zigstd.testExpectEqual(@as(usize, 0), @offsetOf(RawCmapEncRecord, "platform_id"));
    try zigstd.testExpectEqual(@as(usize, 2), @offsetOf(RawCmapEncRecord, "encoding_id"));
    try zigstd.testExpectEqual(@as(usize, 4), @offsetOf(RawCmapEncRecord, "cmap_offset"));
}

test "Format0 matches the wire layout (262 bytes, no padding)" {
    try zigstd.testExpectEqual(@as(usize, 262), @sizeOf(Format0));
    try zigstd.testExpectEqual(@as(usize, 6), @offsetOf(Format0, "glyph_id_array"));
}

test "Format2Header/Format2SubHeader match the wire layout (no padding)" {
    try zigstd.testExpectEqual(@as(usize, 518), @sizeOf(Format2.Format2Header));
    try zigstd.testExpectEqual(@as(usize, 6), @offsetOf(Format2.Format2Header, "sub_header_keys"));
    try zigstd.testExpectEqual(@as(usize, 8), @sizeOf(Format2.Format2SubHeader));
    try zigstd.testExpectEqual(@as(usize, 6), @offsetOf(Format2.Format2SubHeader, "id_range_offset"));
}

test "Format4Header matches the wire layout (14 bytes, no padding)" {
    try zigstd.testExpectEqual(@as(usize, 14), @sizeOf(Format4.Format4Header));
    try zigstd.testExpectEqual(@as(usize, 6), @offsetOf(Format4.Format4Header, "seg_count_x2"));
    try zigstd.testExpectEqual(@as(usize, 8), @offsetOf(Format4.Format4Header, "search_range"));
    try zigstd.testExpectEqual(@as(usize, 10), @offsetOf(Format4.Format4Header, "entry_selector"));
    try zigstd.testExpectEqual(@as(usize, 12), @offsetOf(Format4.Format4Header, "range_shift"));
}

pub const Format4 = struct {
    const Self = @This();

    header: Format4Header,
    end_code: []const u16,
    reserved_pad: u16,
    start_code: []const u16,
    id_delta: []const i16,
    id_range_offset: []const u16,
    glyph_id_array: []const GlyphIndex,

    pub fn getGlyphIndex(self: *const Self, c: u16) GlyphIndex {
        if (c > 0xffff) return .notdef;

        const seg_count_x2 = self.header.seg_count_x2;
        const seg_count = seg_count_x2 >> 1;

        // Microsoft recommends not using the search values in the header and
        // deriving them ourselves.  Because we do want to run on some embedded
        // systems with limited resources, we will still implement the binary search.
        const entry_selector: u16 = zigstd.log2Int(u16, seg_count);
        // note: these are halved
        const search_range: u16 = @as(u16, 1) << @intCast(entry_selector);
        const range_shift: u16 = seg_count - search_range;

        // If c is past the last of the "remainder" segments, the answer lives
        // in the pow2-sized block starting right after them; otherwise it's
        // somewhere in [0, pow2), which is the same size as that block and
        // still fully halvable in entry_selector steps.
        var i: usize = 0;
        if (range_shift > 0 and self.end_code[range_shift - 1] < c) {
            i = range_shift;
        }

        var range = search_range;
        var steps = entry_selector;
        while (steps > 0) : (steps -= 1) {
            range >>= 1;
            if (self.end_code[i + range - 1] < c) {
                i += range;
            }
        }

        if (self.start_code[i] > c) return .notdef;

        const byte_offset_from_id_offset = self.id_range_offset[i];
        // idDelta arithmetic is modulo 65536 per spec: reinterpret the i16 as
        // u16 bits and use wrapping add, which is exactly mod-65536 addition.
        const delta: u16 = @bitCast(self.id_delta[i]);
        if (byte_offset_from_id_offset == 0) {
            return @fromBackingInt(c +% delta);
        }

        // We apply the pointer offset a little different than the spec
        // suggests. We made individual allocations when copying/byte
        // swapping the id_range_offset and glyph_indices out of the input
        // data. This means that we can't just do pointer addition
        //
        // Instead we look at the data as follows
        //
        // [ id range ] [glyph indices ]
        //     |--offs_bytes--|
        //     ^
        //     i
        //
        // To find the index into glyph indices, we just subtract i from
        // id_range.len, and subtract that from the offset

        // The spec's C pointer trick (idRangeOffset[i]/2 + (c - startCode[i])
        // + &idRangeOffset[i]) relies on glyphIdArray immediately following
        // idRangeOffset[] in the wire bytes. We hold them as separate
        // allocations, so instead we compute the equivalent index directly
        // into our own glyph_id_array: the combined idRangeOffset+glyphIdArray
        // index is (idx + idRangeOffset[idx]/2 + (c - startCode[idx]));
        // glyphIdArray starts at combined index seg_count, so subtracting
        // seg_count lands on the right glyph_id_array index.
        const offs_from_loc: usize = @as(usize, byte_offset_from_id_offset) / 2 + (@as(usize, c) - @as(usize, self.start_code[i]));

        const dist_to_end: usize = @as(usize, seg_count) - i;
        const glyph_index_index = offs_from_loc - dist_to_end;
        const raw = self.glyph_id_array[glyph_index_index];
        if (raw == .notdef) return .notdef;
        return @fromBackingInt(@backingInt(raw) +% delta);
    }

    const Format4Header = extern struct {
        format: u16,
        length: u16,
        language: u16,
        seg_count_x2: u16,
        search_range: u16,
        entry_selector: u16,
        range_shift: u16,
    };

    pub fn create(alloc: Allocator, table_data: []const u8, check: bool) !*Self {
        if (check) try bytes.checkSize(table_data, 0, Format4Header);
        const self = try alloc.create(Self);
        errdefer alloc.destroy(self);
        self.header = bytes.typeBEToNative(Format4Header, table_data);

        const seg_count: usize = @intCast(self.header.seg_count_x2 >> 1);
        const header_end = @sizeOf(Format4Header);
        if (check) try bytes.checkArray(table_data, header_end, seg_count, u16);
        self.end_code = try bytes.allocArrayBEToNative(alloc, u16, seg_count, table_data[header_end..]);
        errdefer alloc.free(self.end_code);

        const end_code_end = header_end + seg_count * @sizeOf(u16);
        // just in case some sneaky bastard ever tries to reuse this field for something, we save it
        self.reserved_pad = bytes.readIntBE(u16, table_data[end_code_end..]);

        const start_code_start = end_code_end + @sizeOf(u16);
        if (check) try bytes.checkArray(table_data, start_code_start, seg_count, u16);
        self.start_code = try bytes.allocArrayBEToNative(alloc, u16, seg_count, table_data[start_code_start..]);
        errdefer alloc.free(self.start_code);

        const delta_start = start_code_start + seg_count * @sizeOf(u16);
        if (check) try bytes.checkArray(table_data, delta_start, seg_count, i16);
        self.id_delta = try bytes.allocArrayBEToNative(alloc, i16, seg_count, table_data[delta_start..]);
        errdefer alloc.free(self.id_delta);

        const range_offset_start = delta_start + seg_count * @sizeOf(i16);
        if (check) try bytes.checkArray(table_data, range_offset_start, seg_count, u16);
        self.id_range_offset = try bytes.allocArrayBEToNative(alloc, u16, seg_count, table_data[range_offset_start..]);
        errdefer alloc.free(self.id_range_offset);

        // glyphIdArray fills the rest of the subtable.
        // lookups index the remaining bytes of this table via each subheader's id_range_offset
        const glyph_array_start = range_offset_start + seg_count * @sizeOf(u16);
        if (check) try bytes.checkOffset(table_data, glyph_array_start);
        const glyph_array_len = (table_data.len - glyph_array_start) / @sizeOf(GlyphIndex);
        self.glyph_id_array = try bytes.allocArrayBEToNative(alloc, GlyphIndex, glyph_array_len, table_data[glyph_array_start..]);

        return self;
    }

    pub fn destroy(self: *Self, alloc: Allocator) void {
        alloc.free(self.glyph_id_array);
        alloc.free(self.id_range_offset);
        alloc.free(self.id_delta);
        alloc.free(self.start_code);
        alloc.free(self.end_code);
        alloc.destroy(self);
    }
};

fn setU16BE(dst: []u8, offset: usize, v: u16) void {
    bytes.writeIntBE(u16, dst[offset..][0..2], v);
}

fn setI16BE(dst: []u8, offset: usize, v: i16) void {
    bytes.writeIntBE(i16, dst[offset..][0..2], v);
}

// Builds a format-4 subtable byte buffer from parallel segment arrays plus
// an optional glyphIdArray tail, matching the wire layout exactly.
fn buildFormat4(
    alloc: Allocator,
    end_code: []const u16,
    start_code: []const u16,
    id_delta: []const i16,
    id_range_offset: []const u16,
    glyph_id_array: []const u16,
) ![]u8 {
    const seg_count = end_code.len;
    const header_size = 14;
    const seg_arrays_size = seg_count * 2 * 4 + 2; // endCode+startCode+idRangeOffset (u16) + idDelta (i16) + reservedPad
    const total = header_size + seg_arrays_size + glyph_id_array.len * 2;

    const buf = try alloc.alloc(u8, total);
    @memset(buf, 0);

    setU16BE(buf, 0, 4); // format
    setU16BE(buf, 2, @intCast(total)); // length
    setU16BE(buf, 4, 0); // language
    setU16BE(buf, 6, @intCast(seg_count * 2)); // segCountX2
    setU16BE(buf, 8, 0); // searchRange (unused by our parser, per the spec's own recommendation)
    setU16BE(buf, 10, 0); // entrySelector
    setU16BE(buf, 12, 0); // rangeShift

    var off: usize = header_size;
    for (end_code) |v| {
        setU16BE(buf, off, v);
        off += 2;
    }
    off += 2; // reservedPad, left 0
    for (start_code) |v| {
        setU16BE(buf, off, v);
        off += 2;
    }
    for (id_delta) |v| {
        setI16BE(buf, off, v);
        off += 2;
    }
    for (id_range_offset) |v| {
        setU16BE(buf, off, v);
        off += 2;
    }
    for (glyph_id_array) |v| {
        setU16BE(buf, off, v);
        off += 2;
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

    const fmt4 = try Format4.create(alloc, buf, true);
    defer fmt4.destroy(alloc);

    try zigstd.testExpectEqual(@as(u16, 1), @backingInt(fmt4.getGlyphIndex(10)));
    try zigstd.testExpectEqual(@as(u16, 11), @backingInt(fmt4.getGlyphIndex(20)));
    try zigstd.testExpectEqual(@as(u16, 12), @backingInt(fmt4.getGlyphIndex(30)));
    try zigstd.testExpectEqual(@as(u16, 72), @backingInt(fmt4.getGlyphIndex(90)));
    try zigstd.testExpectEqual(@as(u16, 73), @backingInt(fmt4.getGlyphIndex(153)));
    try zigstd.testExpectEqual(@as(u16, 400), @backingInt(fmt4.getGlyphIndex(480)));
    try zigstd.testExpectEqual(@as(u16, 0), @backingInt(fmt4.getGlyphIndex(0xFFFF)));
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

    const fmt4 = try Format4.create(alloc, buf, true);
    defer fmt4.destroy(alloc);

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

    const fmt4 = try Format4.create(alloc, buf, true);
    defer fmt4.destroy(alloc);

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

    const fmt4 = try Format4.create(alloc, buf, true);
    defer fmt4.destroy(alloc);

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
