const types = @import("../types.zig");
const bytes = @import("../../../../bytes.zig");
const zigstd = @import("../../../../zigstd.zig");

/// The 'name' header.
const StreamLayout = extern struct {
    version: u16,
    count: u16,
    storage_offset: types.Offset16,
};

const records_start = @sizeOf(StreamLayout);

/// Naming table. getRawString() only hands back still-encoded bytes (commonly
/// UTF-16BE, but not always; see the spec): decoding text isn't this table's
/// job.
pub const NameTable = struct {
    const Self = @This();

    table_data: []const u8,

    pub const NameRecord = extern struct {
        platform_id: u16,
        encoding_id: u16,
        language_id: u16,
        name_id: u16,
        length: u16,
        string_offset: types.Offset16,
    };

    pub const LangTagRecord = extern struct {
        length: u16,
        lang_tag_offset: types.Offset16,
    };

    pub fn init(table_data: []const u8, check: bool) !Self {
        if (check) try bytes.checkSize(table_data, 0, StreamLayout);
        const self: Self = .{ .table_data = table_data };

        if (check) {
            try bytes.checkArray(table_data, records_start, self.count(), NameRecord);

            if (self.version() == 1) {
                const lang_tag_count_pos = self.langTagCountPos();
                try bytes.checkSize(table_data, lang_tag_count_pos, u16);
                try bytes.checkArray(table_data, lang_tag_count_pos + @sizeOf(u16), self.langTagCount(), LangTagRecord);
            }

            try bytes.checkOffset(table_data, self.storageOffset());
        }

        return self;
    }

    pub inline fn version(self: Self) u16 {
        return types.field(StreamLayout, "version", self.table_data);
    }

    pub inline fn count(self: Self) u16 {
        return types.field(StreamLayout, "count", self.table_data);
    }

    pub inline fn storageOffset(self: Self) types.Offset16 {
        return types.field(StreamLayout, "storage_offset", self.table_data);
    }

    /// Decodes one name record. index isn't bounds-checked against count()
    /// here; callers are expected to loop within 0..count().
    pub fn getNameRecord(self: Self, index: usize) NameRecord {
        const off = records_start + index * @sizeOf(NameRecord);
        return bytes.typeBEToNative(NameRecord, self.table_data[off..]);
    }

    inline fn langTagCountPos(self: Self) usize {
        return records_start + @as(usize, self.count()) * @sizeOf(NameRecord);
    }

    /// version 1 only.
    pub fn langTagCount(self: Self) u16 {
        return bytes.getU16BE(self.table_data[self.langTagCountPos()..]);
    }

    /// version 1 only. index isn't bounds-checked against langTagCount()
    /// here; callers are expected to loop within 0..langTagCount().
    pub fn getLangTagRecord(self: Self, index: usize) LangTagRecord {
        const off = self.langTagCountPos() + @sizeOf(u16) + index * @sizeOf(LangTagRecord);
        return bytes.typeBEToNative(LangTagRecord, self.table_data[off..]);
    }

    /// Raw string storage area; a given record's bytes are
    /// storage()[record.string_offset..][0..record.length], still in
    /// whatever encoding its platform/encoding ID specifies.
    pub inline fn storage(self: Self) []const u8 {
        return self.table_data[self.storageOffset()..];
    }

    /// The raw (still platform/encoding-specific, e.g. often UTF-16BE)
    /// bytes for one name record. Null if the record's range falls outside
    /// the storage area.
    pub fn getRawString(self: Self, record: NameRecord) ?[]const u8 {
        const store = self.storage();
        const start = record.string_offset;
        const end = @as(usize, start) + record.length;
        if (end > store.len) return null;
        return store[start..end];
    }
};

const Allocator = zigstd.Allocator;

test "name StreamLayout/NameRecord match the OpenType stream layout (no padding)" {
    try zigstd.testExpectEqual(@as(usize, 6), @sizeOf(StreamLayout));
    try zigstd.testExpectEqual(@as(usize, 2), @offsetOf(StreamLayout, "count"));
    try zigstd.testExpectEqual(@as(usize, 4), @offsetOf(StreamLayout, "storage_offset"));
    try zigstd.testExpectEqual(@as(usize, 12), @sizeOf(NameTable.NameRecord));
    try zigstd.testExpectEqual(@as(usize, 6), @offsetOf(NameTable.NameRecord, "name_id"));
    try zigstd.testExpectEqual(@as(usize, 8), @offsetOf(NameTable.NameRecord, "length"));
    try zigstd.testExpectEqual(@as(usize, 10), @offsetOf(NameTable.NameRecord, "string_offset"));
}

fn buildNameV0(alloc: Allocator, records: []const struct { name_id: u16, str: []const u8 }) ![]u8 {
    const header_size = 6;
    const record_size = 12;
    const records_size = records.len * record_size;
    const storage_offset = header_size + records_size;

    var storage_size: usize = 0;
    for (records) |r| storage_size += r.str.len;
    const total = storage_offset + storage_size;

    const buf = try alloc.alloc(u8, total);
    @memset(buf, 0);

    bytes.setU16BE(buf, 0); // version
    bytes.setU16BE(buf[2..], @intCast(records.len)); // count
    bytes.setU16BE(buf[4..], @intCast(storage_offset));

    var storage_off: u16 = 0;
    var rec_off: usize = header_size;
    for (records) |r| {
        bytes.setU16BE(buf[rec_off..], 3); // platform_id: Windows
        bytes.setU16BE(buf[rec_off + 2 ..], 1); // encoding_id: Unicode BMP
        bytes.setU16BE(buf[rec_off + 4 ..], 0x0409); // language_id: en-US
        bytes.setU16BE(buf[rec_off + 6 ..], r.name_id);
        bytes.setU16BE(buf[rec_off + 8 ..], @intCast(r.str.len));
        bytes.setU16BE(buf[rec_off + 10 ..], storage_off);
        @memcpy(buf[storage_offset + storage_off ..][0..r.str.len], r.str);
        storage_off += @intCast(r.str.len);
        rec_off += record_size;
    }

    return buf;
}

test "NameTable.init reads version 0 name records and their raw string data" {
    const alloc = zigstd.testAllocator;
    const buf = try buildNameV0(alloc, &.{
        .{ .name_id = 1, .str = "Test Font" },
        .{ .name_id = 2, .str = "Regular" },
    });
    defer alloc.free(buf);

    const name = try NameTable.init(buf, true);

    try zigstd.testExpectEqual(@as(u16, 2), name.count());
    const r0 = name.getNameRecord(0);
    try zigstd.testExpectEqual(@as(u16, 1), r0.name_id);
    try zigstd.testExpectEqual(@as(u16, 3), r0.platform_id);
    try zigstd.testExpectEqualSlices(u8, "Test Font", name.getRawString(r0).?);
    try zigstd.testExpectEqualSlices(u8, "Regular", name.getRawString(name.getNameRecord(1)).?);
}

test "NameTable.getRawString returns null when a record's range is out of bounds" {
    const alloc = zigstd.testAllocator;
    const buf = try buildNameV0(alloc, &.{.{ .name_id = 1, .str = "Hi" }});
    defer alloc.free(buf);

    const name = try NameTable.init(buf, true);

    var bad_record = name.getNameRecord(0);
    bad_record.length = 200;
    try zigstd.testExpectEqual(@as(?[]const u8, null), name.getRawString(bad_record));
}
