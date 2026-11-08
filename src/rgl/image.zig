// const std = @import("std");

const color = @import("color.zig");
const Point = @import("geometry.zig").Point;
const RGBA = color.RGBA;

const bytes = @import("../bytes.zig");
const zigstd = @import("../zigstd.zig");
const Allocator = zigstd.Allocator;
const StdIo = zigstd.StdIo;

/// A `width` x `height` grid of `Pixel`s, row-major with no row padding.
/// Channels belong in `Pixel` itself (`u8`, `[3]f64`, `RGBA`, ...).
pub fn Image(comptime Pixel: type) type {
    return struct {
        const Self = @This();

        width: u32,
        height: u32,
        data: []Pixel,

        pub fn init(alloc: Allocator, width: u32, height: u32, fill: ?Pixel) !Self {
            const data = try alloc.alloc(Pixel, @as(usize, width) * height);
            if (fill) |f| @memset(data, f);
            return .{ .width = width, .height = height, .data = data };
        }

        pub fn deinit(self: *Self, alloc: Allocator) void {
            alloc.free(self.data);
        }

        pub fn row(self: Self, y: u32) []Pixel {
            return self.data[@as(usize, y) * self.width ..][0..self.width];
        }

        pub fn at(self: Self, x: u32, y: u32) *Pixel {
            return &self.row(y)[x];
        }

        pub fn setPixel(self: Self, x: u32, y: u32, c: Pixel) void {
            self.at(x, y).* = c;
        }

        pub fn setPixelPoint(self: Self, p: Point, c: Pixel) void {
            self.at(p.x, p.y).* = c;
        }

        pub fn getPixel(self: Self, x: u32, y: u32) Pixel {
            return self.at(x, y).*;
        }

        pub fn getPixelPoint(self: Self, p: Point) Pixel {
            return self.at(p.x, p.y).*;
        }

        /// Copies `src`, rows of `src_width` pixels with no padding, in with
        /// its top-left corner at (x, y).
        pub fn blit(self: Self, x: u32, y: u32, src_width: u32, src: []const Pixel) void {
            for (0..src.len / src_width) |j| {
                const dst = self.row(y + @as(u32, @intCast(j)))[x..][0..src_width];
                @memcpy(dst, src[j * src_width ..][0..src_width]);
            }
        }

        /// Area-averages a square RGBA image down to `n` x `n`
        /// using premultiplied alpha so transparent edges stay clean.
        pub fn downsample(self: *const Self, alloc: Allocator, n: u32) !Self {
            var out = try Self.init(alloc, n, n, .blank);
            const scale = @as(f32, @floatFromInt(self.width)) / @as(f32, @floatFromInt(n));

            var dp: Point = .{ .x = 0, .y = 0 };
            while (dp.y < n) : (dp.y += 1) {
                dp.x = 0;
                while (dp.x < n) : (dp.x += 1) {
                    const s0 = dp.toPointF32().scaled(scale).toPoint();
                    const s1 = (Point{ .x = dp.x + 1, .y = dp.y + 1 }).toPointF32().scaled(scale).toPoint();

                    var p: color.RGBAf32 = .blank;
                    var count: f32 = 0;
                    var sp = s0;
                    while (sp.y < s1.y) : (sp.y += 1) {
                        while (sp.x < s1.x) : (sp.x += 1) {
                            const c = self.getPixelPoint(sp).toRGBAf32();
                            p.r += c.r * c.a;
                            p.g += c.g * c.a;
                            p.b += c.b * c.a;
                            p.a += c.a;
                            count += 1;
                        }
                    }
                    if (count == 0) continue;
                    if (p.a > 0) {
                        p.r /= p.a;
                        p.g /= p.a;
                        p.b /= p.a;
                        p.a /= count;
                    }
                    out.setPixelPoint(dp, p.toRGBA());
                }
            }
            return out;
        }

        /// Writes `data` as-is to `path` in `dir`: tightly packed, no header.
        pub fn writeRawFile(self: Self, io: StdIo, dir: StdIo.Dir, path: []const u8) !void {
            try dir.writeFile(io, .{ .sub_path = path, .data = bytes.fromSlice(self.data) });
        }
    };
}

pub const RgbaImage = Image(RGBA);

test "init/setPixel/getPixel round-trip by (x,y) and by Point" {
    const alloc = zigstd.testAllocator;
    var img = try RgbaImage.init(alloc, 3, 2, null);
    defer img.deinit(alloc);

    img.setPixel(1, 0, .{ .r = 10, .g = 20, .b = 30, .a = 40 });
    try zigstd.testExpectTrue(img.getPixel(1, 0).equals(.{ .r = 10, .g = 20, .b = 30, .a = 40 }));

    const p: Point = .{ .x = 2, .y = 1 };
    img.setPixelPoint(p, .{ .r = 1, .g = 2, .b = 3, .a = 4 });
    try zigstd.testExpectTrue(img.getPixelPoint(p).equals(.{ .r = 1, .g = 2, .b = 3, .a = 4 }));
}

test "init fills every pixel when given a fill color" {
    const alloc = zigstd.testAllocator;
    var img = try RgbaImage.init(alloc, 4, 4, RGBA.blank);
    defer img.deinit(alloc);
    for (img.data) |p| try zigstd.testExpectTrue(p.equals(RGBA.blank));
}

test "row and at address multi-channel pixels row-major" {
    const alloc = zigstd.testAllocator;
    var img = try Image([3]f64).init(alloc, 3, 2, @splat(0));
    defer img.deinit(alloc);

    img.at(2, 1).* = .{ 1, 2, 3 };
    try zigstd.testExpectEqual(@as(usize, 3), img.row(1).len);
    try zigstd.testExpectEqualSlices(f64, &.{ 1, 2, 3 }, &img.data[5]);
}

test "blit copies a source image into a sub-rectangle" {
    const alloc = zigstd.testAllocator;
    var dst = try Image(u8).init(alloc, 4, 3, 0);
    defer dst.deinit(alloc);
    dst.blit(1, 1, 2, &.{ 1, 2, 3, 4 });
    try zigstd.testExpectEqualSlices(u8, &.{ 0, 0, 0, 0, 0, 1, 2, 0, 0, 3, 4, 0 }, dst.data);
}

test "writeRawFile writes tightly-packed RGBA bytes" {
    const alloc = zigstd.testAllocator;
    var threaded: StdIo.Threaded = .init(alloc, .{});
    defer threaded.deinit();
    const io = threaded.io();

    var img = try RgbaImage.init(alloc, 2, 1, null);
    defer img.deinit(alloc);
    img.setPixel(0, 0, .{ .r = 1, .g = 2, .b = 3, .a = 4 });
    img.setPixel(1, 0, .{ .r = 5, .g = 6, .b = 7, .a = 8 });

    var tmp = zigstd.testTmpDir(.{});
    defer tmp.cleanup();
    const path = "out.rgba";
    try img.writeRawFile(io, tmp.dir, path);

    const written = try tmp.dir.readFileAlloc(io, path, alloc, .unlimited);
    defer alloc.free(written);
    try zigstd.testExpectEqualSlices(u8, &.{ 1, 2, 3, 4, 5, 6, 7, 8 }, written);
}

test "downsample averages a 2x2 block into one pixel" {
    const alloc = zigstd.testAllocator;
    var src = try RgbaImage.init(alloc, 2, 2, null);
    defer src.deinit(alloc);
    src.setPixel(0, 0, .{ .r = 255, .g = 0, .b = 0, .a = 255 });
    src.setPixel(1, 0, .{ .r = 0, .g = 255, .b = 0, .a = 255 });
    src.setPixel(0, 1, .{ .r = 0, .g = 0, .b = 255, .a = 255 });
    src.setPixel(1, 1, .{ .r = 255, .g = 255, .b = 255, .a = 255 });

    var out = try src.downsample(alloc, 1);
    defer out.deinit(alloc);
    const p = out.getPixel(0, 0);
    try zigstd.testExpectEqual(@as(u8, 255), p.a);
}

test "downsample of a larger image picks the right source blocks" {
    const alloc = zigstd.testAllocator;
    var src = try RgbaImage.init(alloc, 4, 4, color.RGBA.blank);
    defer src.deinit(alloc);
    // top-left quadrant red, everything else stays blank/transparent
    src.setPixel(0, 0, .{ .r = 255, .g = 0, .b = 0, .a = 255 });
    src.setPixel(1, 0, .{ .r = 255, .g = 0, .b = 0, .a = 255 });
    src.setPixel(0, 1, .{ .r = 255, .g = 0, .b = 0, .a = 255 });
    src.setPixel(1, 1, .{ .r = 255, .g = 0, .b = 0, .a = 255 });

    var out = try src.downsample(alloc, 2);
    defer out.deinit(alloc);
    try zigstd.testExpectEqual(@as(u8, 255), out.getPixel(0, 0).r);
    try zigstd.testExpectEqual(@as(u8, 0), out.getPixel(1, 1).a);
}
