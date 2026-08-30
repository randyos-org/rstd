const std = @import("std");

pub const ByteArrayList = std.ArrayList(u8);

// casts v to int_type and appends it to list in little-endian order
pub fn appendIntLE(alloc: std.mem.Allocator, list: *ByteArrayList, int_type: type, v: anytype) !void {
    return appendIntEndian(alloc, list, int_type, v, .little);
}

pub fn appendIntBE(alloc: std.mem.Allocator, list: *ByteArrayList, int_type: type, v: anytype) !void {
    return appendIntEndian(alloc, list, int_type, v, .big);
}

pub fn appendIntEndian(alloc: std.mem.Allocator, list: *ByteArrayList, int_type: type, v: anytype, endian: std.builtin.Endian) !void {
    const size = @sizeOf(int_type);
    var buf: [size]u8 = undefined;
    std.mem.writeInt(int_type, &buf, @intCast(v), endian);
    try list.appendSlice(alloc, &buf);
}

pub fn appendArrayBytes(alloc: std.mem.Allocator, list: *ByteArrayList, data: anytype) !void {
    try list.appendSlice(alloc, std.mem.sliceAsBytes(data));
}

pub inline fn readIntEndian(int_type: type, raw_bytes: []const u8, endian: std.builtin.Endian) int_type {
    return @intCast(std.mem.readInt(int_type, raw_bytes[0..@sizeOf(int_type)], endian));
}

pub inline fn readIntBE(int_type: type, raw_bytes: []const u8) int_type {
    return readIntEndian(int_type, raw_bytes, .big);
}

pub inline fn readIntLE(int_type: type, raw_bytes: []const u8) int_type {
    return readIntEndian(int_type, raw_bytes, .little);
}

// pub fn structAsBytes(comptime T: type, s: *const T) []const u8 {
//     return std.mem.sliceAsBytes(@ptrCast(s));
// }

pub const utf8ToWinString = std.unicode.utf8ToUtf16LeAlloc;
