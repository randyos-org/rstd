//! Pixel primitives on an `RgbaImage`, used by the software backend.

const RgbaImage = @import("image.zig").RgbaImage;
const color = @import("color.zig");
const RGBA = color.RGBA;
const geometry = @import("geometry.zig");
const Point = geometry.Point;
const PointF32 = geometry.PointF32;
const zigstd = @import("../zigstd.zig");

/// Fills an axis-aligned rounded rectangle with hard edges
/// (the downsample step provides antialiasing). Alpha of col is blended
/// over dst. `start`/`end` may extend outside `dst` (a box can be
/// partially off-screen); only the iteration bounds are clamped to
/// `dst`'s size, not `start`/`end` themselves, so `insideRound`'s rounding
/// math still sees the shape's true corners rather than the clipped ones.
pub fn fillRoundRectPixels(dst: *RgbaImage, start: PointF32, end: PointF32, r: f32, c: RGBA) void {
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

pub fn blendPixel(img: *RgbaImage, p: Point, c: RGBA) void {
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

test "fillRoundRectPixels clips a box extending past the image's edges without panicking" {
    const alloc = zigstd.testAllocator;
    var img = try RgbaImage.init(alloc, 10, 10, color.RGBA.blank);
    defer img.deinit(alloc);

    fillRoundRectPixels(
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

test "fillRoundRectPixels fills the center and clips a fully-rounded corner" {
    const alloc = zigstd.testAllocator;
    var img = try RgbaImage.init(alloc, 10, 10, color.RGBA.blank);
    defer img.deinit(alloc);

    fillRoundRectPixels(
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
