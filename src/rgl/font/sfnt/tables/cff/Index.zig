const bytes = @import("../../../../../bytes.zig");
const zigstd = @import("../../../../../zigstd.zig");
const std = @import("std");

/// A CFF INDEX: count + offSize header, an offset array, and concatenated
/// entry data; used for the Name/Top DICT/String/Global Subr/CharStrings/
/// Font DICT/Local Subr INDEXes alike. table_data is exactly this INDEX's
/// own extent (validated and bounded by init()), so byteSize() after that
/// is just table_data.len; callers walking several INDEXes back-to-back
/// (as CffTable.init() does) advance by that length.
const Index = @This();

table_data: []const u8,

pub fn init(table_data: []const u8, check: bool) bytes.BytesError!Index {
    if (check) try bytes.checkSize(table_data, 0, u16);
    const entry_count = bytes.getU16BE(table_data[0..]);
    if (entry_count == 0) return .{ .table_data = table_data[0..2] };

    if (check) try bytes.checkSize(table_data, 2, u8);
    const off_size = table_data[2];
    if (check) try bytes.checkSlice(table_data, offsets_start, (@as(usize, entry_count) + 1) * @as(usize, off_size));

    const self: Index = .{ .table_data = table_data };
    const total = self.byteSizeUnchecked();
    if (check) try bytes.checkSlice(table_data, 0, total);
    return .{ .table_data = table_data[0..total] };
}

pub inline fn count(self: Index) u16 {
    return bytes.getU16BE(self.table_data[0..]);
}

/// Byte width used to encode each offset in the offset array (1-4);
/// meaningless when count() == 0.
pub inline fn offSize(self: Index) u8 {
    return self.table_data[2];
}

const offsets_start = 3;

fn rawOffset(self: Index, index: usize) usize {
    const off_size: usize = self.offSize();
    const off = offsets_start + index * off_size;
    var v: usize = 0;
    for (0..off_size) |b| v = (v << 8) | self.table_data[off + b];
    return v;
}

fn dataStart(self: Index) usize {
    return offsets_start + (@as(usize, self.count()) + 1) * @as(usize, self.offSize());
}

/// Computes the total byte size from table_data directly, without
/// re-checking bounds; callers must have already validated at least
/// through the offset array (as init() does) before calling this.
fn byteSizeUnchecked(self: Index) usize {
    return self.dataStart() + self.rawOffset(self.count()) - 1;
}

/// Total byte size of this INDEX (header + offset array + data), i.e.
/// where the next structure begins; equal to table_data.len once
/// init() has bounded table_data to exactly this INDEX's extent.
pub inline fn byteSize(self: Index) usize {
    return self.table_data.len;
}

/// Entry `index`'s raw, unparsed bytes (a DICT, a charstring, or a
/// string, depending on which INDEX this is). Per spec, offsets in the
/// offset array are 1-based; callers are trusted to pass a valid index
/// (< count()), matching this module's other raw index accessors.
pub fn get(self: Index, index: usize) []const u8 {
    const data_start = self.dataStart();
    const start = self.rawOffset(index) - 1;
    const end = self.rawOffset(index + 1) - 1;
    return self.table_data[data_start + start .. data_start + end];
}

/// Test-support: builds a well-formed INDEX (offSize always 1, fine for
/// these small fixtures) containing `entries`. Reused by Dict.zig's and
/// __root__.zig's own tests to build realistic CFF byte streams, not part
/// of the table-reading API.
pub fn buildIndex(alloc: zigstd.Allocator, entries: []const []const u8) ![]u8 {
    if (entries.len == 0) {
        const buf = try alloc.alloc(u8, 2);
        bytes.setU16BE(buf, 0);
        return buf;
    }

    var total_data: usize = 0;
    for (entries) |e| total_data += e.len;
    const off_size: u8 = 1;
    const header = 3 + (entries.len + 1) * off_size;
    const buf = try alloc.alloc(u8, header + total_data);

    bytes.setU16BE(buf, @intCast(entries.len));
    buf[2] = off_size;

    var cursor: usize = 1;
    for (entries, 0..) |e, i| {
        buf[3 + i] = @intCast(cursor);
        cursor += e.len;
    }
    buf[3 + entries.len] = @intCast(cursor);

    var data_pos = header;
    for (entries) |e| {
        @memcpy(buf[data_pos..][0..e.len], e);
        data_pos += e.len;
    }

    return buf;
}

test "Index.get/count/byteSize walk a small INDEX" {
    const alloc = zigstd.testAllocator;
    const buf = try buildIndex(alloc, &.{ "ab", "cde", "" });
    defer alloc.free(buf);

    const idx = try Index.init(buf, true);
    try zigstd.testExpectEqual(@as(u16, 3), idx.count());
    try zigstd.testExpectEqual(buf.len, idx.byteSize());
    try zigstd.testExpectEqual(true, std.mem.eql(u8, "ab", idx.get(0)));
    try zigstd.testExpectEqual(true, std.mem.eql(u8, "cde", idx.get(1)));
    try zigstd.testExpectEqual(@as(usize, 0), idx.get(2).len);
}

test "Index.init treats a zero-count INDEX as an empty 2-byte structure" {
    const alloc = zigstd.testAllocator;
    const buf = try buildIndex(alloc, &.{});
    defer alloc.free(buf);

    const idx = try Index.init(buf, true);
    try zigstd.testExpectEqual(@as(u16, 0), idx.count());
    try zigstd.testExpectEqual(@as(usize, 2), idx.byteSize());
}
