// const std = @import("std");

const color = @import("color.zig");
const Point = @import("geometry.zig").Point;
const RGBA = color.RGBA;

const bytes = @import("../bytes.zig");
const zigstd = @import("../zigstd.zig");
const Allocator = zigstd.Allocator;
const StdIo = zigstd.StdIo;

const Self = @This();
width: u32,
height: u32,
data: []RGBA,

pub fn init(alloc: Allocator, width: u32, height: u32, fill: ?RGBA) !Self {
    const data = try alloc.alloc(RGBA, width * height);
    errdefer alloc.free(data);
    if (fill) |f| @memset(data, f);
    return .{ .width = width, .height = height, .data = data };
}

pub fn destroy(self: *Self, alloc: Allocator) void {
    alloc.free(self.data);
}

pub fn setPixel(img: *Self, x: u32, y: u32, c: RGBA) void {
    img.data[(y * img.width) + x] = c;
}

pub fn setPixelPoint(img: *Self, p: Point, c: RGBA) void {
    img.data[(p.y * img.width) + p.x] = c;
}

pub fn getPixel(img: *const Self, x: u32, y: u32) RGBA {
    return img.data[(y * img.width) + x];
}

pub fn getPixelPoint(img: *const Self, p: Point) RGBA {
    return img.data[(p.y * img.width) + p.x];
}

pub fn writeRgbaFile(self: *const Self, io: StdIo, path: []const u8) !void {
    // RGBA is a packed struct of four u8 fields, so img.data is already
    // byte-for-byte identical to the output format -- no copy needed.
    try StdIo.Dir.cwd().writeFile(io, .{ .sub_path = path, .data = bytes.fromSlice(self.data) });
}

test "init/setPixel/getPixel round-trip by (x,y) and by Point" {
    const alloc = zigstd.testAllocator;
    var img = try init(alloc, 3, 2, null);
    defer img.destroy(alloc);

    img.setPixel(1, 0, .{ .r = 10, .g = 20, .b = 30, .a = 40 });
    try zigstd.testExpectTrue(img.getPixel(1, 0).equals(.{ .r = 10, .g = 20, .b = 30, .a = 40 }));

    const p: Point = .{ .x = 2, .y = 1 };
    img.setPixelPoint(p, .{ .r = 1, .g = 2, .b = 3, .a = 4 });
    try zigstd.testExpectTrue(img.getPixelPoint(p).equals(.{ .r = 1, .g = 2, .b = 3, .a = 4 }));
}

test "init fills every pixel when given a fill color" {
    const alloc = zigstd.testAllocator;
    var img = try init(alloc, 4, 4, RGBA.blank);
    defer img.destroy(alloc);
    for (img.data) |p| try zigstd.testExpectTrue(p.equals(RGBA.blank));
}

test "writeRgbaFile writes tightly-packed RGBA bytes" {
    const alloc = zigstd.testAllocator;
    var threaded: StdIo.Threaded = .init(alloc, .{});
    defer threaded.deinit();
    const io = threaded.io();

    var img = try init(alloc, 2, 1, null);
    defer img.destroy(alloc);
    img.setPixel(0, 0, .{ .r = 1, .g = 2, .b = 3, .a = 4 });
    img.setPixel(1, 0, .{ .r = 5, .g = 6, .b = 7, .a = 8 });

    const path = "_test_image_write_rgba_output.rgba";
    try img.writeRgbaFile(io, path);
    defer StdIo.Dir.cwd().deleteFile(io, path) catch {};

    const written = try StdIo.Dir.cwd().readFileAlloc(io, path, alloc, .unlimited);
    defer alloc.free(written);
    try zigstd.testExpectEqualSlices(u8, &.{ 1, 2, 3, 4, 5, 6, 7, 8 }, written);
}
