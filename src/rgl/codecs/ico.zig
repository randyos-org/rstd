const std = @import("std");

const Image = @import("../Image.zig");
const png = @import("png.zig");

fn appendU16LE(alloc: std.mem.Allocator, list: *std.ArrayList(u8), v: u16) !void {
    try list.appendSlice(alloc, &std.mem.toBytes(std.mem.nativeToLittle(u16, v)));
}
fn appendU32LE(alloc: std.mem.Allocator, list: *std.ArrayList(u8), v: u32) !void {
    try list.appendSlice(alloc, &std.mem.toBytes(std.mem.nativeToLittle(u32, v)));
}

/// writeICO packs PNG-encoded images (Vista+ supports PNG entries) into one
/// .ico.
pub fn writeICO(allocator: std.mem.Allocator, io: std.Io, name: []const u8, imgs: []const Image) !void {
    const Entry = struct { w: i32, h: i32, png: []const u8 };
    var entries: std.ArrayList(Entry) = .empty;
    defer entries.deinit(allocator);

    for (imgs) |im| {
        const png_bytes = try png.encode(allocator, &im);
        try entries.append(allocator, .{ .w = @intCast(im.width), .h = @intCast(im.height), .png = png_bytes });
    }

    var out: std.ArrayList(u8) = .empty;
    defer out.deinit(allocator);

    // ICONDIR
    try appendU16LE(allocator, &out, 0); // reserved
    try appendU16LE(allocator, &out, 1); // type: icon
    try appendU16LE(allocator, &out, @intCast(entries.items.len));

    var offset: u32 = 6 + 16 * @as(u32, @intCast(entries.items.len));
    for (entries.items) |e| {
        const dim = struct {
            fn f(v: i32) u8 {
                if (v >= 256) return 0; // 0 means 256 in the ICO spec
                return @intCast(v);
            }
        }.f;
        try out.append(allocator, dim(e.w));
        try out.append(allocator, dim(e.h));
        try out.append(allocator, 0); // color palette
        try out.append(allocator, 0); // reserved
        try appendU16LE(allocator, &out, 1); // planes
        try appendU16LE(allocator, &out, 32); // bpp
        try appendU32LE(allocator, &out, @intCast(e.png.len)); // bytes
        try appendU32LE(allocator, &out, offset); // offset
        offset += @intCast(e.png.len);
    }
    for (entries.items) |e| {
        try out.appendSlice(allocator, e.png);
        allocator.free(e.png);
    }

    try std.Io.Dir.cwd().writeFile(io, .{ .sub_path = name, .data = out.items });
}

test "writeICO produces a valid one-entry ICONDIR header" {
    const alloc = std.testing.allocator;
    var threaded: std.Io.Threaded = .init(alloc, .{});
    defer threaded.deinit();
    const io = threaded.io();

    var img = try Image.init(alloc, 4, 4, .{ .r = 0, .g = 0, .b = 0, .a = 0 });
    defer img.destroy(alloc);

    const path = "_test_ico_write_ico_output.ico";
    try writeICO(alloc, io, path, &.{img});
    defer std.Io.Dir.cwd().deleteFile(io, path) catch {};

    const written = try std.Io.Dir.cwd().readFileAlloc(io, path, alloc, .unlimited);
    defer alloc.free(written);
    // ICONDIR: reserved=0, type=1 (icon), count=1
    try std.testing.expectEqual(@as(u16, 0), std.mem.readInt(u16, written[0..2], .little));
    try std.testing.expectEqual(@as(u16, 1), std.mem.readInt(u16, written[2..4], .little));
    try std.testing.expectEqual(@as(u16, 1), std.mem.readInt(u16, written[4..6], .little));
}
