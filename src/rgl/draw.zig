const std = @import("std");
const Image = @import("Image.zig");
const color = @import("color.zig");
const RGBA = color.RGBA;

pub const Point = struct {
    const Self = @This();
    x: u32,
    y: u32,

    pub const zero: Self = .{ .x = 0, .y = 0 };

    pub fn asPointF32(self: Self) PointF32 {
        return .{ .x = @floatFromInt(self.x), .y = @floatFromInt(self.y) };
    }

    pub fn asPointI32(self: Self) PointI32 {
        return .{ .x = @intCast(self.x), .y = @intCast(self.y) };
    }
};

pub const PointI32 = struct {
    const Self = @This();
    x: i32,
    y: i32,

    pub const zero: Self = .{ .x = 0, .y = 0 };

    pub fn asPointF32(self: Self) PointF32 {
        return .{ .x = @floatFromInt(self.x), .y = @floatFromInt(self.y) };
    }

    pub fn asPoint(self: Self) Point {
        return .{ .x = @intCast(self.x), .y = @intCast(self.y) };
    }

    pub fn add(self: Self, other: Self) Self {
        return .{ .x = self.x + other.x, .y = self.y + other.y };
    }
};

pub const PointF32 = struct {
    const Self = @This();
    x: f32,
    y: f32,

    pub const zero: Self = .{ .x = 0, .y = 0 };

    pub fn asPoint(self: Self) Point {
        return .{ .x = @round(self.x), .y = @round(self.y) };
    }

    pub fn asPointI32(self: Self) PointI32 {
        return .{ .x = @round(self.x), .y = @round(self.y) };
    }

    pub fn scaled(self: Self, s: f32) Self {
        return .{ .x = self.x * s, .y = self.y * s };
    }

    pub fn add(self: Self, other: Self) Self {
        return .{ .x = self.x + other.x, .y = self.y + other.y };
    }
};

/// fillRoundRect fills an axis-aligned rounded rectangle with hard edges
/// (the downsample step provides antialiasing). Alpha of col is blended
/// over dst.
pub fn fillRoundRect(dst: *Image, start: PointF32, end: PointF32, r: f32, c: RGBA) void {
    var p: PointF32 = start;
    while (p.y < end.y) : (p.y += 1) {
        p.x = start.x;
        while (p.x < end.x) : (p.x += 1) {
            const center = p.add(.{ .x = 0.5, .y = 0.5 });
            if (!insideRound(center, start, end, r)) continue;
            blendPixel(dst, p.asPoint(), c);
        }
    }
}

fn insideRound(p: PointF32, p0: PointF32, p1: PointF32, r: f32) bool {
    // Clamp point to the inner rectangle (the rounded-corner centers box).
    const cx = @max(p0.x + r, @min(p.x, p1.x - r));
    const cy = @max(p0.y + r, @min(p.y, p1.y - r));
    const dx = p.x - cx;
    const dy = p.y - cy;
    return dx * dx + dy * dy <= r * r;
}

pub fn blendPixel(img: *Image, p: Point, c: RGBA) void {
    if (c.a == 0xFF) {
        img.setPixelPoint(p, c);
        return;
    }
    const cf32 = c.asRGBAf32();
    const bgf32 = img.getPixelPoint(p).asRGBAf32();
    const blend: color.RGBAf32 = .{
        .r = cf32.r * cf32.a + bgf32.r * (1 - cf32.a),
        .g = cf32.g * cf32.a + bgf32.g * (1 - cf32.a),
        .b = cf32.b * cf32.a + bgf32.b * (1 - cf32.a),
        .a = cf32.a + bgf32.a * (1 - cf32.a),
    };
    img.setPixelPoint(p, blend.asRGBA());
}

/// area-averages a (square) source image down to an n x n image
/// using premultiplied alpha so transparent edges stay clean.
pub fn downsample(allocator: std.mem.Allocator, src: *const Image, n: u32) !Image {
    var out = try Image.init(allocator, n, n, .blank);
    const scale = @as(f32, @floatFromInt(src.width)) / @as(f32, @floatFromInt(n));

    var dp: Point = .{ .x = 0, .y = 0 };
    while (dp.y < n) : (dp.y += 1) {
        dp.x = 0;
        while (dp.x < n) : (dp.x += 1) {
            const s0 = dp.asPointF32().scaled(scale).asPoint();
            const s1 = (Point{ .x = dp.x + 1, .y = dp.y + 1 }).asPointF32().scaled(scale).asPoint();

            var p: color.RGBAf32 = .blank;
            var count: f32 = 0;
            var sp = s0;
            while (sp.y < s1.y) : (sp.y += 1) {
                while (sp.x < s1.x) : (sp.x += 1) {
                    const c = src.getPixelPoint(sp).asRGBAf32();
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
            out.setPixelPoint(dp, p.asRGBA());
        }
    }
    return out;
}

test "Point/PointF32/PointI32 conversions and arithmetic" {
    const p: Point = .{ .x = 3, .y = 5 };
    try std.testing.expectEqual(PointF32{ .x = 3, .y = 5 }, p.asPointF32());
    try std.testing.expectEqual(PointI32{ .x = 3, .y = 5 }, p.asPointI32());

    const pf: PointF32 = .{ .x = 2.6, .y = 4.4 };
    try std.testing.expectEqual(Point{ .x = 3, .y = 4 }, pf.asPoint());
    try std.testing.expectEqual(PointI32{ .x = 3, .y = 4 }, pf.asPointI32());

    const pi: PointI32 = .{ .x = 2, .y = 7 };
    try std.testing.expectEqual(PointF32{ .x = 2, .y = 7 }, pi.asPointF32());
    try std.testing.expectEqual(Point{ .x = 2, .y = 7 }, pi.asPoint());
    try std.testing.expectEqual(PointI32{ .x = 7, .y = 12 }, pi.add(.{ .x = 5, .y = 5 }));

    try std.testing.expectEqual(PointF32{ .x = 5, .y = 10 }, (PointF32{ .x = 1, .y = 2 }).scaled(5));
    try std.testing.expectEqual(PointF32{ .x = 4, .y = 6 }, (PointF32{ .x = 1, .y = 2 }).add(.{ .x = 3, .y = 4 }));
}

test "fillRoundRect fills the center and clips a fully-rounded corner" {
    const alloc = std.testing.allocator;
    var img = try Image.init(alloc, 10, 10, color.RGBA.blank);
    defer img.destroy(alloc);

    fillRoundRect(
        &img,
        .{ .x = 0, .y = 0 },
        .{ .x = 10, .y = 10 },
        3,
        .{ .r = 255, .g = 0, .b = 0, .a = 255 },
    );

    // center should be filled
    try std.testing.expect(img.getPixel(5, 5).equals(.{ .r = 255, .g = 0, .b = 0, .a = 255 }));
    // exact corner should be clipped away by the rounding radius
    try std.testing.expect(img.getPixel(0, 0).equals(color.RGBA.blank));
}

test "downsample averages a 2x2 block into one pixel" {
    const alloc = std.testing.allocator;
    var src = try Image.init(alloc, 2, 2, null);
    defer src.destroy(alloc);
    src.setPixel(0, 0, .{ .r = 255, .g = 0, .b = 0, .a = 255 });
    src.setPixel(1, 0, .{ .r = 0, .g = 255, .b = 0, .a = 255 });
    src.setPixel(0, 1, .{ .r = 0, .g = 0, .b = 255, .a = 255 });
    src.setPixel(1, 1, .{ .r = 255, .g = 255, .b = 255, .a = 255 });

    var out = try downsample(alloc, &src, 1);
    defer out.destroy(alloc);
    const p = out.getPixel(0, 0);
    try std.testing.expectEqual(@as(u8, 255), p.a);
}

test "downsample of a larger image picks the right source blocks" {
    const alloc = std.testing.allocator;
    var src = try Image.init(alloc, 4, 4, color.RGBA.blank);
    defer src.destroy(alloc);
    // top-left quadrant red, everything else stays blank/transparent
    src.setPixel(0, 0, .{ .r = 255, .g = 0, .b = 0, .a = 255 });
    src.setPixel(1, 0, .{ .r = 255, .g = 0, .b = 0, .a = 255 });
    src.setPixel(0, 1, .{ .r = 255, .g = 0, .b = 0, .a = 255 });
    src.setPixel(1, 1, .{ .r = 255, .g = 0, .b = 0, .a = 255 });

    var out = try downsample(alloc, &src, 2);
    defer out.destroy(alloc);
    try std.testing.expectEqual(@as(u8, 255), out.getPixel(0, 0).r);
    try std.testing.expectEqual(@as(u8, 0), out.getPixel(1, 1).a);
}
