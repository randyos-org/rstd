const std = @import("std");

pub const ByteArrayList = std.ArrayList(u8);

// casts v to int_type and appends it to list in little-endian order
pub fn appendIntLE(alloc: std.mem.Allocator, list: *ByteArrayList, int_type: type, v: anytype) !void {
    const size = @sizeOf(int_type);
    var buf: [size]u8 = undefined;
    std.mem.writeInt(int_type, &buf, @intCast(v), .little);
    try list.appendSlice(alloc, &buf);
}

pub const utf8ToWinString = std.unicode.utf8ToUtf16LeAlloc;
