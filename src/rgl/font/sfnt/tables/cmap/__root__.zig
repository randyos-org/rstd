const types = @import("../../types.zig");
const bytes = @import("../../../../../bytes.zig");

const GlyphIndex = types.GlyphIndex;
const Offset32 = types.Offset32;
const Format0 = @import("fmt0.zig").Format0;
const Format2 = @import("fmt2.zig").Format2;
const Format4 = @import("fmt4.zig").Format4;
const Format6 = @import("fmt6.zig").Format6;
const Format8 = @import("fmt8.zig").Format8;
const Format10 = @import("fmt10.zig").Format10;
const Format12 = @import("fmt12.zig").Format12;
const Format13 = @import("fmt13.zig").Format13;
const Format14 = @import("fmt14.zig").Format14;
const zigstd = @import("../../../../../zigstd.zig");

pub const PlatformId = enum(u16) {
    unicode = 0,
    mac = 1,
    iso = 2,
    microsoft = 3,
    custom = 4,
};

pub const UnicodeEncodingId = enum(u16) {
    version_1_0 = 0,
    version_1_1 = 1,
    iso_10646 = 2,
    unicode_2_bmp = 3,
    unicode_2_full = 4,
    unicode_2_var = 5,
    unicode_full = 6,
};

pub const IsoEncodingId = enum(u16) {
    ascii = 0,
    iso10646 = 1,
    iso8859_1 = 2,
};

pub const MicrosoftEncodingId = enum(u16) {
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

pub const EncodingId = union(enum) {
    unicode: UnicodeEncodingId,
    iso: IsoEncodingId,
    microsoft: MicrosoftEncodingId,
};

/// The 'cmap' header.
const StreamLayout = extern struct {
    version: u16,
    num_tables: u16,
};

/// One entry from the 'cmap' encoding record array, decoded fresh from
/// cmap_bytes on each CmapTable.getRecord() call rather than being cached -
/// records are 8 bytes each and a font rarely has more than a handful. Its
/// extern layout doubles as the wire schema for @offsetOf/@sizeOf.
pub const EncRecord = extern struct {
    const Self = @This();

    platform_id: PlatformId,
    _encoding_id: u16,
    cmap_offset: Offset32,

    pub fn encodingId(self: Self) ?EncodingId {
        return switch (self.platform_id) {
            .unicode => EncodingId{ .unicode = @fromBackingInt(@intCast(self._encoding_id)) },
            .iso => EncodingId{ .iso = @fromBackingInt(@intCast(self._encoding_id)) },
            .microsoft => EncodingId{ .microsoft = @fromBackingInt(@intCast(self._encoding_id)) },
            else => null,
        };
    }

    pub fn isUnicodeBmp(self: Self) bool {
        const encoding_id = self.encodingId() orelse return false;
        return (self.platform_id == .unicode and encoding_id.unicode == .unicode_2_bmp) or
            (self.platform_id == .microsoft and encoding_id.microsoft == .unicode_bmp);
    }

    /// A full-repertoire Unicode encoding (format 12/13-capable, i.e. not
    /// limited to the BMP); (3,10) on Windows, or (0,4)/(0,6) on the
    /// Unicode platform.
    pub fn isUnicodeFull(self: Self) bool {
        const encoding_id = self.encodingId() orelse return false;
        return (self.platform_id == .unicode and
            (encoding_id.unicode == .unicode_2_full or encoding_id.unicode == .unicode_full)) or
            (self.platform_id == .microsoft and encoding_id.microsoft == .unicode_full);
    }
};

pub const CmapSubtable = union(enum) {
    const Self = @This();

    fmt0: Format0,
    fmt2: Format2,
    fmt4: Format4,
    fmt6: Format6,
    fmt8: Format8,
    fmt10: Format10,
    fmt12: Format12,
    fmt13: Format13,
    fmt14: Format14,

    /// raw format number for formats not yet implemented
    unsupported: u16,

    /// table_data must already start at this subtable's own bytes (i.e.
    /// pre-sliced from the parent 'cmap' bytes at the record's cmap_offset).
    /// Every format is a plain view over table_data, so this never
    /// allocates.
    pub fn init(table_data: []const u8, check: bool) !Self {
        const format = bytes.readIntBE(u16, table_data);
        return switch (format) {
            0 => .{ .fmt0 = try Format0.init(table_data, check) },
            2 => .{ .fmt2 = try Format2.init(table_data, check) },
            4 => .{ .fmt4 = try Format4.init(table_data, check) },
            6 => .{ .fmt6 = try Format6.init(table_data, check) },
            8 => .{ .fmt8 = try Format8.init(table_data, check) },
            10 => .{ .fmt10 = try Format10.init(table_data, check) },
            12 => .{ .fmt12 = try Format12.init(table_data, check) },
            13 => .{ .fmt13 = try Format13.init(table_data, check) },
            14 => .{ .fmt14 = try Format14.init(table_data, check) },
            else => .{ .unsupported = format },
        };
    }

    pub fn getGlyphIndex(self: Self, codepoint: u32) GlyphIndex {
        return switch (self) {
            .fmt0 => |fmt0| fmt0.getGlyphIndex(codepoint),
            .fmt2 => |fmt2| fmt2.getGlyphIndex(codepoint),
            .fmt4 => |fmt4| fmt4.getGlyphIndex(codepoint),
            .fmt6 => |fmt6| fmt6.getGlyphIndex(codepoint),
            .fmt8 => |fmt8| fmt8.getGlyphIndex(codepoint),
            .fmt10 => |fmt10| fmt10.getGlyphIndex(codepoint),
            .fmt12 => |fmt12| fmt12.getGlyphIndex(codepoint),
            .fmt13 => |fmt13| fmt13.getGlyphIndex(codepoint),
            .fmt14, .unsupported => .notdef,
        };
    }
};

/// Character to glyph mapping.
pub const CmapTable = struct {
    const Self = @This();

    cmap_bytes: []const u8,

    pub fn init(cmap_bytes: []const u8, check: bool) !Self {
        if (check) try bytes.checkSize(cmap_bytes, 0, StreamLayout);
        const self: Self = .{ .cmap_bytes = cmap_bytes };
        if (check) try bytes.checkArray(cmap_bytes, @sizeOf(StreamLayout), self.numTables(), EncRecord);
        return self;
    }

    pub inline fn version(self: Self) u16 {
        return types.field(StreamLayout, "version", self.cmap_bytes);
    }

    pub inline fn numTables(self: Self) u16 {
        return types.field(StreamLayout, "num_tables", self.cmap_bytes);
    }

    /// Decodes one encoding record. index isn't bounds-checked against
    /// numTables() here; callers are expected to loop within 0..numTables().
    pub fn getRecord(self: Self, index: usize) EncRecord {
        const off = @sizeOf(StreamLayout) + index * @sizeOf(EncRecord);
        return bytes.typeBEToNative(EncRecord, self.cmap_bytes[off..]);
    }

    /// Scans records for a Unicode BMP-compatible encoding (Windows Unicode
    /// BMP or Unicode platform BMP); the encoding renderers almost always
    /// want. Returns null if the font has none.
    pub fn findUnicodeBmpRecord(self: Self) ?EncRecord {
        var i: usize = 0;
        while (i < self.numTables()) : (i += 1) {
            const record = self.getRecord(i);
            if (record.isUnicodeBmp()) return record;
        }
        return null;
    }

    /// Scans records for a full-repertoire Unicode encoding (format 12/13,
    /// covering codepoints past the BMP). Returns null if the font has none.
    pub fn findUnicodeFullRecord(self: Self) ?EncRecord {
        var i: usize = 0;
        while (i < self.numTables()) : (i += 1) {
            const record = self.getRecord(i);
            if (record.isUnicodeFull()) return record;
        }
        return null;
    }

    /// Parses the subtable a given record points to.
    pub fn getSubtable(self: Self, record: EncRecord, check: bool) !CmapSubtable {
        return CmapSubtable.init(self.cmap_bytes[record.cmap_offset..], check);
    }
};

const Allocator = zigstd.Allocator;

test "cmap StreamLayout/EncRecord match the OpenType stream layout (no padding)" {
    try zigstd.testExpectEqual(@as(usize, 4), @sizeOf(StreamLayout));
    try zigstd.testExpectEqual(@as(usize, 2), @offsetOf(StreamLayout, "num_tables"));
    try zigstd.testExpectEqual(@as(usize, 8), @sizeOf(EncRecord));
    try zigstd.testExpectEqual(@as(usize, 0), @offsetOf(EncRecord, "platform_id"));
    try zigstd.testExpectEqual(@as(usize, 2), @offsetOf(EncRecord, "_encoding_id"));
    try zigstd.testExpectEqual(@as(usize, 4), @offsetOf(EncRecord, "cmap_offset"));
}

test "CmapSubtable.getGlyphIndex dispatches to whichever format backend was parsed" {
    // A minimal format-6 subtable: header (10 bytes) + one glyph id.
    var buf: [12]u8 = @splat(0);
    bytes.setU16BE(&buf, 6); // format
    bytes.setU16BE(buf[2..], buf.len); // length
    bytes.setU16BE(buf[6..], 0x41); // first_code
    bytes.setU16BE(buf[8..], 1); // entry_count
    bytes.setU16BE(buf[10..], 77); // glyph_id_array[0]

    const subtable = try CmapSubtable.init(&buf, true);

    try zigstd.testExpectEqual(@as(u16, 77), @backingInt(subtable.getGlyphIndex(0x41)));
    try zigstd.testExpectEqual(GlyphIndex.notdef, subtable.getGlyphIndex(0x42));
}

test "CmapSubtable.getGlyphIndex returns .notdef for format 14 (not a plain codepoint lookup)" {
    // A minimal, empty format-14 subtable: header (10 bytes), zero records.
    var buf: [10]u8 = @splat(0);
    bytes.setU16BE(&buf, 14); // format
    bytes.setU32BE(buf[6..], 0); // numVarSelectorRecords

    const subtable = try CmapSubtable.init(&buf, true);

    try zigstd.testExpectEqual(GlyphIndex.notdef, subtable.getGlyphIndex(0x41));
}

// Builds a full 'cmap' table: header + two encoding records (a Mac record
// and a Windows Unicode BMP record) both pointing at the same format-6
// subtable appended after the record array.
fn buildCmap(alloc: Allocator) ![]u8 {
    const header_size = 4;
    const record_size = 8;
    const subtable_offset = header_size + 2 * record_size;
    const subtable_size = 12;
    const buf = try alloc.alloc(u8, subtable_offset + subtable_size);
    @memset(buf, 0);

    bytes.setU16BE(buf, 0); // version
    bytes.setU16BE(buf[2..], 2); // num_tables

    // record 0: mac platform, not a Unicode BMP encoding
    bytes.setU16BE(buf[4..], @backingInt(PlatformId.mac));
    bytes.setU16BE(buf[6..], 0);
    bytes.setU32BE(buf[8..], subtable_offset);

    // record 1: Windows Unicode BMP
    bytes.setU16BE(buf[12..], @backingInt(PlatformId.microsoft));
    bytes.setU16BE(buf[14..], @backingInt(MicrosoftEncodingId.unicode_bmp));
    bytes.setU32BE(buf[16..], subtable_offset);

    // format-6 subtable: header (10 bytes) + one glyph id
    const sub = buf[subtable_offset..];
    bytes.setU16BE(sub, 6); // format
    bytes.setU16BE(sub[2..], subtable_size); // length
    bytes.setU16BE(sub[6..], 0x41); // first_code
    bytes.setU16BE(sub[8..], 1); // entry_count
    bytes.setU16BE(sub[10..], 77); // glyph_id_array[0]

    return buf;
}

test "CmapTable.getRecord decodes platform/encoding/offset for each record" {
    const alloc = zigstd.testAllocator;
    const buf = try buildCmap(alloc);
    defer alloc.free(buf);

    const cmap = try CmapTable.init(buf, true);
    try zigstd.testExpectEqual(@as(u16, 2), cmap.numTables());

    const r0 = cmap.getRecord(0);
    try zigstd.testExpectEqual(PlatformId.mac, r0.platform_id);
    try zigstd.testExpectTrue(!r0.isUnicodeBmp());

    const r1 = cmap.getRecord(1);
    try zigstd.testExpectEqual(PlatformId.microsoft, r1.platform_id);
    try zigstd.testExpectTrue(r1.isUnicodeBmp());
}

test "CmapTable.findUnicodeBmpRecord skips non-matching records" {
    const alloc = zigstd.testAllocator;
    const buf = try buildCmap(alloc);
    defer alloc.free(buf);

    const cmap = try CmapTable.init(buf, true);
    const record = cmap.findUnicodeBmpRecord().?;
    try zigstd.testExpectEqual(PlatformId.microsoft, record.platform_id);
}

test "EncRecord.isUnicodeFull recognizes the Windows and Unicode-platform full-repertoire encodings" {
    const windows_full: EncRecord = .{
        .platform_id = .microsoft,
        ._encoding_id = @backingInt(MicrosoftEncodingId.unicode_full),
        .cmap_offset = 0,
    };
    try zigstd.testExpectTrue(windows_full.isUnicodeFull());
    try zigstd.testExpectTrue(!windows_full.isUnicodeBmp());

    const unicode_full: EncRecord = .{
        .platform_id = .unicode,
        ._encoding_id = @backingInt(UnicodeEncodingId.unicode_full),
        .cmap_offset = 0,
    };
    try zigstd.testExpectTrue(unicode_full.isUnicodeFull());

    const windows_bmp: EncRecord = .{
        .platform_id = .microsoft,
        ._encoding_id = @backingInt(MicrosoftEncodingId.unicode_bmp),
        .cmap_offset = 0,
    };
    try zigstd.testExpectTrue(!windows_bmp.isUnicodeFull());
}

test "CmapTable.findUnicodeFullRecord finds a (3,10) record and skips a BMP-only one" {
    const alloc = zigstd.testAllocator;
    const buf = try buildCmap(alloc);
    defer alloc.free(buf);

    // buildCmap() only has mac + Windows-BMP records; no full-repertoire
    // encoding present yet.
    const cmap = try CmapTable.init(buf, true);
    try zigstd.testExpectTrue(cmap.findUnicodeFullRecord() == null);

    // Same layout, but promote the Windows record to (3,10) full-repertoire.
    var full_buf = try alloc.dupe(u8, buf);
    defer alloc.free(full_buf);
    bytes.setU16BE(full_buf[14..], @backingInt(MicrosoftEncodingId.unicode_full));

    const full_cmap = try CmapTable.init(full_buf, true);
    const record = full_cmap.findUnicodeFullRecord().?;
    try zigstd.testExpectEqual(PlatformId.microsoft, record.platform_id);
    try zigstd.testExpectTrue(record.isUnicodeFull());
}

test "CmapTable.getSubtable parses the resolved record and looks up a glyph" {
    const alloc = zigstd.testAllocator;
    const buf = try buildCmap(alloc);
    defer alloc.free(buf);

    const cmap = try CmapTable.init(buf, true);
    const record = cmap.findUnicodeBmpRecord().?;

    const subtable = try cmap.getSubtable(record, true);

    try zigstd.testExpectEqual(@as(u16, 77), @backingInt(subtable.getGlyphIndex(0x41)));
}
