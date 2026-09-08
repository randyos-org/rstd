const sysinfo = @import("builtin");
const std = @import("std");
const zigstd = @import("zigstd.zig");
const Allocator = zigstd.Allocator;

pub const Endian = std.builtin.Endian;
pub const ByteArrayList = std.ArrayList(u8);

// casts v to int_type and appends it to list in little-endian order
pub fn appendIntLE(alloc: Allocator, list: *ByteArrayList, int_type: type, v: anytype) !void {
    return appendIntEndian(alloc, list, int_type, v, .little);
}

pub fn appendIntBE(alloc: Allocator, list: *ByteArrayList, int_type: type, v: anytype) !void {
    return appendIntEndian(alloc, list, int_type, v, .big);
}

pub fn appendIntEndian(alloc: Allocator, list: *ByteArrayList, int_type: type, v: anytype, endian: Endian) !void {
    const size = @sizeOf(int_type);
    var buf: [size]u8 = undefined;
    writeIntEndian(int_type, &buf, v, endian);
    try list.appendSlice(alloc, &buf);
}

pub fn appendArrayBytes(alloc: Allocator, list: *ByteArrayList, data: anytype) !void {
    try list.appendSlice(alloc, fromSlice(data));
}

pub inline fn readIntEndian(int_type: type, raw_bytes: []const u8, endian: Endian) int_type {
    return @intCast(std.mem.readInt(int_type, raw_bytes[0..@sizeOf(int_type)], endian));
}

pub inline fn readIntBE(int_type: type, raw_bytes: []const u8) int_type {
    return readIntEndian(int_type, raw_bytes, .big);
}

pub inline fn readIntLE(int_type: type, raw_bytes: []const u8) int_type {
    return readIntEndian(int_type, raw_bytes, .little);
}

pub inline fn writeIntEndian(int_type: type, raw_bytes: []u8, v: anytype, endian: Endian) void {
    std.mem.writeInt(int_type, raw_bytes[0..@sizeOf(int_type)], @intCast(v), endian);
}

pub inline fn writeIntBE(int_type: type, raw_bytes: []u8, v: anytype) void {
    writeIntEndian(int_type, raw_bytes, v, .big);
}

pub inline fn writeIntLE(int_type: type, raw_bytes: []u8, v: anytype) void {
    writeIntEndian(int_type, raw_bytes, v, .little);
}

pub inline fn nativeToBE(val: anytype) @TypeOf(val) {
    return swapIfNotNative(val, .big);
}

pub inline fn nativeToLE(val: anytype) @TypeOf(val) {
    return swapIfNotNative(val, .little);
}

pub fn swapIfNotNative(val: anytype, endian: Endian) @TypeOf(val) {
    if (sysinfo.cpu.arch.endian() == endian) return val;
    var ret = val;
    std.mem.byteSwap(@TypeOf(val), &ret);
    return ret;
}

pub const toNative = swapIfNotNative;

pub inline fn beToNative(val: anytype) @TypeOf(val) {
    return toNative(val, .big);
}

pub inline fn leToNative(val: anytype) @TypeOf(val) {
    return toNative(val, .little);
}

pub fn arraySwapIfNotNative(comptime Elem: type, dest: []Elem, src: []const u8, endian: Endian) void {
    const dest_bytes = fromSlice(dest);
    @memcpy(dest_bytes, src[0..dest_bytes.len]);
    if (sysinfo.cpu.arch.endian() != endian) {
        std.mem.byteSwapAllElements(Elem, dest);
    }
}

pub inline fn arrayNativeToBE(comptime Elem: type, dest: []Elem, src: []const u8) void {
    arraySwapIfNotNative(Elem, dest, src, .big);
}

pub inline fn arrayNativeToLE(comptime Elem: type, dest: []Elem, src: []const u8) void {
    arraySwapIfNotNative(Elem, dest, src, .little);
}

pub const arrayToNative = arraySwapIfNotNative;

pub inline fn arrayBEToNative(comptime Elem: type, dest: []Elem, src: []const u8) void {
    arrayToNative(Elem, dest, src, .big);
}

pub inline fn arrayLEToNative(comptime Elem: type, dest: []Elem, src: []const u8) void {
    arrayToNative(Elem, dest, src, .little);
}

pub fn typeSwapIfNotNative(comptime T: type, val: []const u8, endian: Endian) T {
    if (sysinfo.cpu.arch.endian() == endian) return viewAsValue(T, val).*;
    const size = @sizeOf(T);
    var ret: [size]u8 align(@alignOf(T)) = undefined;
    @memcpy(&ret, val[0..size]);
    std.mem.byteSwap(T, viewAsValue(T, &ret));
    return viewAsValue(T, &ret).*;
}

pub inline fn typeNativeToBE(comptime T: type, val: []const u8) T {
    return typeSwapIfNotNative(T, val, .big);
}

pub inline fn typeNativeToLE(comptime T: type, val: []const u8) T {
    return typeSwapIfNotNative(T, val, .little);
}

pub const typeToNative = typeSwapIfNotNative;

pub inline fn typeBEToNative(comptime T: type, val: []const u8) T {
    return typeToNative(T, val, .big);
}

pub inline fn typeLEToNative(comptime T: type, val: []const u8) T {
    return typeToNative(T, val, .little);
}

pub fn allocTypeSwapIfNotNative(alloc: Allocator, comptime T: type, val: []const u8, endian: Endian) !*T {
    const ret = try alloc.create(T);
    ret.* = typeSwapIfNotNative(T, val, endian);
    return ret;
}

pub inline fn allocTypeNativeToBE(alloc: Allocator, comptime T: type, val: []const u8) !*T {
    return allocTypeSwapIfNotNative(alloc, T, val, .big);
}

pub inline fn allocTypeNativeToLE(alloc: Allocator, comptime T: type, val: []const u8) !*T {
    return allocTypeSwapIfNotNative(alloc, T, val, .little);
}

pub const allocTypeToNative = allocTypeSwapIfNotNative;

pub inline fn allocTypeBEToNative(alloc: Allocator, comptime T: type, val: []const u8) !*T {
    return allocTypeToNative(alloc, T, val, .big);
}

pub inline fn allocTypeLEToNative(alloc: Allocator, comptime T: type, val: []const u8) !*T {
    return allocTypeToNative(alloc, T, val, .little);
}

pub fn allocArraySwapIfNotNative(alloc: Allocator, comptime Elem: type, count: usize, src: []const u8, endian: Endian) ![]Elem {
    const dest = try alloc.alloc(Elem, count);
    arraySwapIfNotNative(Elem, dest, src, endian);
    return dest;
}

pub inline fn allocArrayNativeToBE(alloc: Allocator, comptime Elem: type, count: usize, src: []const u8) ![]Elem {
    return allocArraySwapIfNotNative(alloc, Elem, count, src, .big);
}

pub inline fn allocArrayNativeToLE(alloc: Allocator, comptime Elem: type, count: usize, src: []const u8) ![]Elem {
    return allocArraySwapIfNotNative(alloc, Elem, count, src, .little);
}

pub const allocArrayToNative = allocArraySwapIfNotNative;

pub inline fn allocArrayBEToNative(alloc: Allocator, comptime Elem: type, count: usize, src: []const u8) ![]Elem {
    return allocArrayToNative(alloc, Elem, count, src, .big);
}

pub inline fn allocArrayLEToNative(alloc: Allocator, comptime Elem: type, count: usize, src: []const u8) ![]Elem {
    return allocArrayToNative(alloc, Elem, count, src, .little);
}

pub const viewAsValue = std.mem.bytesAsValue;
pub const viewAsSlice = std.mem.bytesAsSlice;
pub const fromSlice = std.mem.sliceAsBytes;

pub const utf8ToWinStrConst = std.unicode.utf8ToUtf16LeStringLiteral;
pub const L = utf8ToWinStrConst;

pub const utf8ToWinStrBuf = std.unicode.utf8ToUtf16Le;
pub const utf8ToWinStrAllocZ = std.unicode.utf8ToUtf16LeAllocZ;
pub const utf8ToWinStrAlloc = std.unicode.utf8ToUtf16LeAlloc;

pub const winStrToUtf8Buf = std.unicode.utf16LeToUtf8;
pub const winStrToUtf8Alloc = std.unicode.utf16LeToUtf8Alloc;
pub const winStrToUtf8AllocZ = std.unicode.utf16LeToUtf8AllocZ;

pub const BytesError = error{
    EndOfStream,
};

pub inline fn checkSize(val: []const u8, offset: usize, comptime T: type) BytesError!void {
    return checkOffset(val, offset + @sizeOf(T));
}

pub inline fn checkOffset(val: []const u8, offset: usize) BytesError!void {
    if (offset > val.len) return BytesError.EndOfStream;
}

pub inline fn checkRange(val: []const u8, start: usize, end: usize) BytesError!void {
    try checkOffset(val, start);
    try checkOffset(val, end);
}

pub inline fn checkSlice(val: []const u8, offset: usize, len: usize) BytesError!void {
    return checkRange(val, offset, offset + len);
}

pub inline fn checkArray(val: []const u8, offset: usize, count: usize, comptime T: type) BytesError!void {
    return checkSlice(val, offset, count * @sizeOf(T));
}

pub inline fn safeRange(val: []const u8, start: usize, end: usize, check: bool) BytesError![]const u8 {
    if (check) try checkRange(val, start, end);
    return val[start..end];
}

pub inline fn safeSlice(val: []const u8, offset: usize, len: usize, check: bool) BytesError![]const u8 {
    if (check) try checkSlice(val, offset, len);
    return val[offset .. offset + len];
}
