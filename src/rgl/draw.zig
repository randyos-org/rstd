// const std = @import("std");
const Image = @import("Image.zig");
const color = @import("color.zig");
const RGBA = color.RGBA;
const geometry = @import("geometry.zig");
const Point = geometry.Point;
const PointF32 = geometry.PointF32;

const zigstd = @import("../zigstd.zig");
const Allocator = zigstd.Allocator;

/// fillRoundRect fills an axis-aligned rounded rectangle with hard edges
/// (the downsample step provides antialiasing). Alpha of col is blended
/// over dst. `start`/`end` may extend outside `dst` (a box can be
/// partially off-screen) -- only the iteration bounds are clamped to
/// `dst`'s size, not `start`/`end` themselves, so `insideRound`'s rounding
/// math still sees the shape's true corners rather than the clipped ones.
pub fn fillRoundRect(dst: *Image, start: PointF32, end: PointF32, r: f32, c: RGBA) void {
    const w: f32 = @floatFromInt(dst.width);
    const h: f32 = @floatFromInt(dst.height);
    const clip_start: PointF32 = .{ .x = @max(start.x, 0), .y = @max(start.y, 0) };
    const clip_end: PointF32 = .{ .x = @min(end.x, w), .y = @min(end.y, h) };

    var p: PointF32 = clip_start;
    while (p.y < clip_end.y) : (p.y += 1) {
        p.x = clip_start.x;
        while (p.x < clip_end.x) : (p.x += 1) {
            const center = p.add(.{ .x = 0.5, .y = 0.5 });
            if (!insideRound(center, start, end, r)) continue;
            blendPixel(dst, p.toPoint(), c);
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
    const cf32 = c.toRGBAf32();
    const bgf32 = img.getPixelPoint(p).toRGBAf32();
    const blend: color.RGBAf32 = .{
        .r = cf32.r * cf32.a + bgf32.r * (1 - cf32.a),
        .g = cf32.g * cf32.a + bgf32.g * (1 - cf32.a),
        .b = cf32.b * cf32.a + bgf32.b * (1 - cf32.a),
        .a = cf32.a + bgf32.a * (1 - cf32.a),
    };
    img.setPixelPoint(p, blend.toRGBA());
}

/// area-averages a (square) source image down to an n x n image
/// using premultiplied alpha so transparent edges stay clean.
pub fn downsample(alloc: Allocator, src: *const Image, n: u32) !Image {
    var out = try Image.init(alloc, n, n, .blank);
    const scale = @as(f32, @floatFromInt(src.width)) / @as(f32, @floatFromInt(n));

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
                    const c = src.getPixelPoint(sp).toRGBAf32();
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

const testExpectEqual = zigstd.testExpectEqual;

test "fillRoundRect clips a box extending past the image's edges without panicking" {
    const alloc = zigstd.testAllocator;
    var img = try Image.init(alloc, 10, 10, color.RGBA.blank);
    defer img.destroy(alloc);

    fillRoundRect(
        &img,
        .{ .x = -5, .y = -5 },
        .{ .x = 5, .y = 5 },
        0,
        .{ .r = 255, .g = 0, .b = 0, .a = 255 },
    );

    // visible corner of the box is filled...
    try zigstd.testExpectTrue(img.getPixel(0, 0).equals(.{ .r = 255, .g = 0, .b = 0, .a = 255 }));
    try zigstd.testExpectTrue(img.getPixel(4, 4).equals(.{ .r = 255, .g = 0, .b = 0, .a = 255 }));
    // ...and nothing past the box's true (partially off-screen) edge is.
    try zigstd.testExpectTrue(img.getPixel(5, 5).equals(color.RGBA.blank));
}

test "fillRoundRect fills the center and clips a fully-rounded corner" {
    const alloc = zigstd.testAllocator;
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
    try zigstd.testExpectTrue(img.getPixel(5, 5).equals(.{ .r = 255, .g = 0, .b = 0, .a = 255 }));
    // exact corner should be clipped away by the rounding radius
    try zigstd.testExpectTrue(img.getPixel(0, 0).equals(color.RGBA.blank));
}

test "downsample averages a 2x2 block into one pixel" {
    const alloc = zigstd.testAllocator;
    var src = try Image.init(alloc, 2, 2, null);
    defer src.destroy(alloc);
    src.setPixel(0, 0, .{ .r = 255, .g = 0, .b = 0, .a = 255 });
    src.setPixel(1, 0, .{ .r = 0, .g = 255, .b = 0, .a = 255 });
    src.setPixel(0, 1, .{ .r = 0, .g = 0, .b = 255, .a = 255 });
    src.setPixel(1, 1, .{ .r = 255, .g = 255, .b = 255, .a = 255 });

    var out = try downsample(alloc, &src, 1);
    defer out.destroy(alloc);
    const p = out.getPixel(0, 0);
    try testExpectEqual(@as(u8, 255), p.a);
}

test "downsample of a larger image picks the right source blocks" {
    const alloc = zigstd.testAllocator;
    var src = try Image.init(alloc, 4, 4, color.RGBA.blank);
    defer src.destroy(alloc);
    // top-left quadrant red, everything else stays blank/transparent
    src.setPixel(0, 0, .{ .r = 255, .g = 0, .b = 0, .a = 255 });
    src.setPixel(1, 0, .{ .r = 255, .g = 0, .b = 0, .a = 255 });
    src.setPixel(0, 1, .{ .r = 255, .g = 0, .b = 0, .a = 255 });
    src.setPixel(1, 1, .{ .r = 255, .g = 0, .b = 0, .a = 255 });

    var out = try downsample(alloc, &src, 2);
    defer out.destroy(alloc);
    try testExpectEqual(@as(u8, 255), out.getPixel(0, 0).r);
    try testExpectEqual(@as(u8, 0), out.getPixel(1, 1).a);
}
