const std = @import("std");
const zigstd = @import("zigstd.zig");
const Allocator = zigstd.Allocator;

pub const ByteArrayList = std.ArrayList(u8);

// casts v to int_type and appends it to list in little-endian order
pub fn appendIntLE(alloc: Allocator, list: *ByteArrayList, int_type: type, v: anytype) !void {
    return appendIntEndian(alloc, list, int_type, v, .little);
}

pub fn appendIntBE(alloc: Allocator, list: *ByteArrayList, int_type: type, v: anytype) !void {
    return appendIntEndian(alloc, list, int_type, v, .big);
}

pub fn appendIntEndian(alloc: Allocator, list: *ByteArrayList, int_type: type, v: anytype, endian: std.builtin.Endian) !void {
    const size = @sizeOf(int_type);
    var buf: [size]u8 = undefined;
    writeIntEndian(int_type, &buf, v, endian);
    try list.appendSlice(alloc, &buf);
}

pub fn appendArrayBytes(alloc: Allocator, list: *ByteArrayList, data: anytype) !void {
    try list.appendSlice(alloc, fromSlice(data));
}

pub inline fn readIntEndian(int_type: type, raw_bytes: *const [@sizeOf(int_type)]u8, endian: std.builtin.Endian) int_type {
    return @intCast(std.mem.readInt(int_type, raw_bytes, endian));
}

pub inline fn readIntBE(int_type: type, raw_bytes: *const [@sizeOf(int_type)]u8) int_type {
    return readIntEndian(int_type, raw_bytes, .big);
}

pub inline fn readIntLE(int_type: type, raw_bytes: *const [@sizeOf(int_type)]u8) int_type {
    return readIntEndian(int_type, raw_bytes, .little);
}

pub inline fn writeIntEndian(int_type: type, raw_bytes: *[@sizeOf(int_type)]u8, v: anytype, endian: std.builtin.Endian) void {
    std.mem.writeInt(int_type, raw_bytes, @intCast(v), endian);
}

pub inline fn writeIntBE(int_type: type, raw_bytes: *[@sizeOf(int_type)]u8, v: anytype) void {
    writeIntEndian(int_type, raw_bytes, v, .big);
}

pub inline fn writeIntLE(int_type: type, raw_bytes: *[@sizeOf(int_type)]u8, v: anytype) void {
    writeIntEndian(int_type, raw_bytes, v, .little);
}

// pub fn structAsBytes(comptime T: type, s: *const T) []const u8 {
//     return fromBytes(@ptrCast(s));
// }

pub const fromSlice = std.mem.sliceAsBytes;

pub const utf8ToWinStrConst = std.unicode.utf8ToUtf16LeStringLiteral;
pub const L = utf8ToWinStrConst;

pub const utf8ToWinStrBuf = std.unicode.utf8ToUtf16Le;
pub const utf8ToWinStrAllocZ = std.unicode.utf8ToUtf16LeAllocZ;
pub const utf8ToWinStrAlloc = std.unicode.utf8ToUtf16LeAlloc;

pub const winStrToUtf8Buf = std.unicode.utf16LeToUtf8;
pub const winStrToUtf8Alloc = std.unicode.utf16LeToUtf8Alloc;
pub const winStrToUtf8AllocZ = std.unicode.utf16LeToUtf8AllocZ;
