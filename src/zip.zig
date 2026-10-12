//! Minimal ZIP archive writer: stored (uncompressed) entries only.
//! `std.zip` only reads/extracts archives, so the
//! writing side lives here; record signatures and the method enum still come
//! from `std.zip`. Timestamps are fixed so output is reproducible.
const std = @import("std");
const zigstd = @import("zigstd.zig");
const bytes = @import("bytes.zig");

const Allocator = zigstd.Allocator;
const ByteArrayList = bytes.ByteArrayList;
const zip = std.zip;

/// 1980-01-01 00:00:00 in MS-DOS date/time format.
const dos_date: u16 = (1 << 5) | 1;
const dos_time: u16 = 0;
/// Version 2.0: the minimum that supports stored entries in a directory tree.
const zip_version: u16 = 20;
/// General-purpose flag bit 11: file names are UTF-8.
const flag_utf8: u16 = 1 << 11;
const method_store: u16 = @backingInt(zip.CompressionMethod.store);

const Entry = struct {
    name: []const u8,
    crc32: u32,
    size: u32,
    offset: u32,
};

pub const Archive = struct {
    gpa: Allocator,
    out: ByteArrayList = .empty,
    entries: std.ArrayList(Entry) = .empty,

    pub fn init(gpa: Allocator) Archive {
        return .{ .gpa = gpa };
    }

    pub fn deinit(self: *Archive) void {
        for (self.entries.items) |e| self.gpa.free(e.name);
        self.entries.deinit(self.gpa);
        self.out.deinit(self.gpa);
    }

    fn put(self: *Archive, comptime T: type, v: T) !void {
        try bytes.appendIntLE(self.gpa, &self.out, T, v);
    }

    fn putSlice(self: *Archive, s: []const u8) !void {
        try self.out.appendSlice(self.gpa, s);
    }

    /// Appends a file. `name` is copied; `data` is written immediately.
    pub fn add(self: *Archive, name: []const u8, data: []const u8) !void {
        const offset = std.math.cast(u32, self.out.items.len) orelse return error.ZipTooLarge;
        const size = std.math.cast(u32, data.len) orelse return error.ZipTooLarge;
        const name_len = std.math.cast(u16, name.len) orelse return error.ZipNameTooLong;
        const crc32 = std.hash.Crc32.hash(data);

        try self.putSlice(&zip.local_file_header_sig);
        try self.put(u16, zip_version);
        try self.put(u16, flag_utf8);
        try self.put(u16, method_store);
        try self.put(u16, dos_time);
        try self.put(u16, dos_date);
        try self.put(u32, crc32);
        try self.put(u32, size); // compressed size
        try self.put(u32, size); // uncompressed size
        try self.put(u16, name_len);
        try self.put(u16, 0); // extra field length
        try self.putSlice(name);
        try self.putSlice(data);

        const owned_name = try self.gpa.dupe(u8, name);
        errdefer self.gpa.free(owned_name);
        try self.entries.append(self.gpa, .{
            .name = owned_name,
            .crc32 = crc32,
            .size = size,
            .offset = offset,
        });
    }

    /// Writes the central directory and returns the finished archive bytes,
    /// owned by the caller. Only `deinit` may be called afterwards.
    pub fn finish(self: *Archive) ![]u8 {
        const count = std.math.cast(u16, self.entries.items.len) orelse return error.ZipTooManyEntries;
        const cd_offset = std.math.cast(u32, self.out.items.len) orelse return error.ZipTooLarge;

        for (self.entries.items) |e| {
            try self.putSlice(&zip.central_file_header_sig);
            try self.put(u16, zip_version); // version made by (MS-DOS host)
            try self.put(u16, zip_version); // version needed
            try self.put(u16, flag_utf8);
            try self.put(u16, method_store);
            try self.put(u16, dos_time);
            try self.put(u16, dos_date);
            try self.put(u32, e.crc32);
            try self.put(u32, e.size);
            try self.put(u32, e.size);
            try self.put(u16, @intCast(e.name.len));
            try self.put(u16, 0); // extra field length
            try self.put(u16, 0); // comment length
            try self.put(u16, 0); // disk number start
            try self.put(u16, 0); // internal attributes
            try self.put(u32, 0); // external attributes
            try self.put(u32, e.offset);
            try self.putSlice(e.name);
        }

        const cd_size = std.math.cast(u32, self.out.items.len - cd_offset) orelse return error.ZipTooLarge;
        try self.putSlice(&zip.end_record_sig);
        try self.put(u16, 0); // this disk
        try self.put(u16, 0); // disk with central directory
        try self.put(u16, count);
        try self.put(u16, count);
        try self.put(u32, cd_size);
        try self.put(u32, cd_offset);
        try self.put(u16, 0); // comment length

        return self.out.toOwnedSlice(self.gpa);
    }
};

test "archive has a valid end record" {
    const gpa = zigstd.testAllocator;
    var archive: Archive = .init(gpa);
    defer archive.deinit();
    try archive.add("a.txt", "hello");
    try archive.add("dir/b.txt", "world!");
    const out = try archive.finish();
    defer gpa.free(out);

    try zigstd.testExpectEqualSlices(u8, &zip.local_file_header_sig, out[0..4]);
    // End-of-central-directory record has no comment, so it is the last 22 bytes.
    const end = out[out.len - 22 ..];
    try zigstd.testExpectEqualSlices(u8, &zip.end_record_sig, end[0..4]);
    try zigstd.testExpectEqual(@as(u16, 2), bytes.readIntLE(u16, end[10..12]));
    const cd_offset = bytes.readIntLE(u32, end[16..20]);
    try zigstd.testExpectEqualSlices(u8, &zip.central_file_header_sig, out[cd_offset..][0..4]);
    // First entry's data immediately follows its 30-byte header and name.
    try zigstd.testExpectEqualStrings("hello", out[30 + "a.txt".len ..][0..5]);
}
