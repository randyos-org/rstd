//! Software backend: every draw call executes immediately, straight into
//! `SwSurface`'s own pixel buffer (an rgl.RgbaImage). Unlike the vulkan
//! backend, there's nothing to batch or flush; a pixel is final the
//! moment its draw call returns.

const RgbaImage = @import("image.zig").RgbaImage;
const color = @import("color.zig");
const draw = @import("draw.zig");
const RGBA = color.RGBA;
const geometry = @import("geometry.zig");
const Point = geometry.Point;
const PointF32 = geometry.PointF32;
const Rect = geometry.Rect;
const vga_font = @import("font/vga/__root__.zig");
const zigstd = @import("../zigstd.zig");
const Allocator = zigstd.Allocator;

pub const InitOptions = struct {
    fill: ?RGBA = null,
};

pub const SwSurface = struct {
    const Self = @This();
    target: RgbaImage,

    pub fn init(alloc: Allocator, width: u32, height: u32, options: InitOptions) !Self {
        return .{ .target = try RgbaImage.init(alloc, width, height, options.fill) };
    }

    pub fn deinit(self: *Self, alloc: Allocator) void {
        self.target.deinit(alloc);
    }
};

pub fn beginFrame(surface: *SwSurface, clear: ?RGBA) void {
    if (clear) |c| @memset(surface.target.data, c);
}

/// Eager backend: every pixel is already final by the time a draw call
/// returns, so there's nothing left to flush.
pub fn endFrame(surface: *SwSurface) void {
    _ = surface;
}

/// Blends `c` onto `(x, y)` if it's within `img`'s bounds, silently
/// dropping it otherwise; callers routinely pass partially- or
/// fully-off-screen coordinates (a box straddling the edge, a line to a
/// negative point), and `RgbaImage.setPixel`/`draw.blendPixel` do no bounds
/// checking of their own.
fn blendPixelClamped(img: *RgbaImage, x: i32, y: i32, c: RGBA) void {
    if (x < 0 or y < 0 or x >= img.width or y >= img.height) return;
    draw.blendPixel(img, .{ .x = @intCast(x), .y = @intCast(y) }, c);
}

pub fn drawLine(surface: *SwSurface, p0: PointF32, p1: PointF32, width: f32, c: RGBA) void {
    if (width <= 1.0) {
        drawThinLine(&surface.target, p0, p1, c);
    } else {
        fillOrientedRect(&surface.target, p0, p1, width, c);
    }
}

/// Bresenham, 1px wide.
fn drawThinLine(img: *RgbaImage, p0: PointF32, p1: PointF32, c: RGBA) void {
    var x0: i32 = @trunc(@round(p0.x));
    var y0: i32 = @trunc(@round(p0.y));
    const x1: i32 = @trunc(@round(p1.x));
    const y1: i32 = @trunc(@round(p1.y));

    const dx: i32 = @intCast(@abs(x1 - x0));
    const dy: i32 = -@as(i32, @intCast(@abs(y1 - y0)));
    const sx: i32 = if (x0 < x1) 1 else -1;
    const sy: i32 = if (y0 < y1) 1 else -1;
    var err = dx + dy;

    while (true) {
        blendPixelClamped(img, x0, y0, c);
        if (x0 == x1 and y0 == y1) break;
        const e2 = 2 * err;
        if (e2 >= dy) {
            err += dy;
            x0 += sx;
        }
        if (e2 <= dx) {
            err += dx;
            y0 += sy;
        }
    }
}

/// Fills the thin rotated rect a thick line sweeps out (a filled square
/// end-cap, for a zero-length one); the same shape shaders/shapes.zig's
/// rotation support draws on the Vulkan side, just rasterized directly
/// here instead of via a GPU pipeline.
fn fillOrientedRect(img: *RgbaImage, p0: PointF32, p1: PointF32, width: f32, c: RGBA) void {
    const d = p1.sub(p0);
    const length = @sqrt(d.x * d.x + d.y * d.y);
    const half_w = width / 2;
    if (length == 0) {
        draw.fillRoundRectPixels(img, .{ .x = p0.x - half_w, .y = p0.y - half_w }, .{ .x = p0.x + half_w, .y = p0.y + half_w }, 0, c);
        return;
    }
    const unit: PointF32 = .{ .x = d.x / length, .y = d.y / length };
    const normal: PointF32 = .{ .x = -unit.y, .y = unit.x };

    const corners = [_]PointF32{
        p0.add(normal.scaled(half_w)),
        p0.add(normal.scaled(-half_w)),
        p1.add(normal.scaled(half_w)),
        p1.add(normal.scaled(-half_w)),
    };
    var min_x = corners[0].x;
    var max_x = corners[0].x;
    var min_y = corners[0].y;
    var max_y = corners[0].y;
    for (corners[1..]) |corner| {
        min_x = @min(min_x, corner.x);
        max_x = @max(max_x, corner.x);
        min_y = @min(min_y, corner.y);
        max_y = @max(max_y, corner.y);
    }
    const img_w: f32 = @floatFromInt(img.width);
    const img_h: f32 = @floatFromInt(img.height);
    const x0: i32 = @trunc(@floor(@max(min_x, 0)));
    const y0: i32 = @trunc(@floor(@max(min_y, 0)));
    const x1: i32 = @trunc(@ceil(@min(max_x, img_w)));
    const y1: i32 = @trunc(@ceil(@min(max_y, img_h)));

    var y = y0;
    while (y < y1) : (y += 1) {
        var x = x0;
        while (x < x1) : (x += 1) {
            const p: PointF32 = .{ .x = @as(f32, @floatFromInt(x)) + 0.5, .y = @as(f32, @floatFromInt(y)) + 0.5 };
            const rel = p.sub(p0);
            const along = rel.x * unit.x + rel.y * unit.y;
            const across = rel.x * normal.x + rel.y * normal.y;
            if (along < 0 or along > length or @abs(across) > half_w) continue;
            blendPixelClamped(img, x, y, c);
        }
    }
}

pub fn fillRect(surface: *SwSurface, rect: Rect, c: RGBA) void {
    draw.fillRoundRectPixels(&surface.target, rect.topLeftF32(), rect.bottomRightF32(), 0, c);
}

pub fn strokeRect(surface: *SwSurface, rect: Rect, width: f32, c: RGBA) void {
    const tl = rect.topLeftF32();
    const br = rect.bottomRightF32();
    const tr: PointF32 = .{ .x = br.x, .y = tl.y };
    const bl: PointF32 = .{ .x = tl.x, .y = br.y };
    drawLine(surface, tl, tr, width, c);
    drawLine(surface, tr, br, width, c);
    drawLine(surface, br, bl, width, c);
    drawLine(surface, bl, tl, width, c);
}

pub fn fillRoundRect(surface: *SwSurface, rect: Rect, radius: f32, c: RGBA) void {
    draw.fillRoundRectPixels(&surface.target, rect.topLeftF32(), rect.bottomRightF32(), radius, c);
}

pub fn fillCircle(surface: *SwSurface, center: PointF32, radius: f32, c: RGBA) void {
    draw.fillRoundRectPixels(
        &surface.target,
        .{ .x = center.x - radius, .y = center.y - radius },
        .{ .x = center.x + radius, .y = center.y + radius },
        radius,
        c,
    );
}

/// Blits `text` using the embedded 8x16 VGA/CP437 font, one glyph cell per
/// byte of `text` (no UTF-8 decoding; non-CP437 bytes just show
/// whichever glyph that byte happens to map to). Background is left
/// untouched (transparent, gdi32 TRANSPARENT-mode style); only set bits
/// blend `c`.
pub fn drawText(surface: *SwSurface, pos: Point, text: []const u8, c: RGBA) void {
    const img = &surface.target;
    var pen_x: i32 = @intCast(pos.x);
    const pen_y: i32 = @intCast(pos.y);
    for (text) |ch| {
        const rows = vga_font.glyph(ch);
        for (rows.*, 0..) |row_bits, row| {
            var bit_i: i32 = 7;
            while (bit_i >= 0) : (bit_i -= 1) {
                const bit: u3 = @intCast(bit_i);
                if ((row_bits >> bit) & 1 != 0) {
                    blendPixelClamped(img, pen_x + (7 - bit_i), pen_y + @as(i32, @intCast(row)), c);
                }
            }
        }
        pen_x += vga_font.width;
    }
}

/// Straight alpha blit of `src` at `pos`, clipped to `surface`'s bounds.
/// No scaling: 1:1 source pixel to destination pixel.
pub fn drawImage(surface: *SwSurface, pos: Point, src: *const RgbaImage) void {
    const img = &surface.target;
    var sy: u32 = 0;
    while (sy < src.height) : (sy += 1) {
        var sx: u32 = 0;
        while (sx < src.width) : (sx += 1) {
            blendPixelClamped(img, @intCast(pos.x + sx), @intCast(pos.y + sy), src.getPixel(sx, sy));
        }
    }
}

const red: RGBA = .{ .r = 255, .g = 0, .b = 0, .a = 255 };

test "fillRect fills exactly the requested box" {
    const alloc = zigstd.testAllocator;
    var surface = try SwSurface.init(alloc, 10, 10, .{ .fill = color.RGBA.blank });
    defer surface.deinit(alloc);

    fillRect(&surface, .{ .x = 2, .y = 2, .w = 3, .h = 3 }, red);
    try zigstd.testExpectTrue(surface.target.getPixel(2, 2).equals(red));
    try zigstd.testExpectTrue(surface.target.getPixel(4, 4).equals(red));
    try zigstd.testExpectTrue(surface.target.getPixel(5, 5).equals(color.RGBA.blank));
}

test "drawLine (thin) draws a horizontal 1px line" {
    const alloc = zigstd.testAllocator;
    var surface = try SwSurface.init(alloc, 10, 10, .{ .fill = color.RGBA.blank });
    defer surface.deinit(alloc);

    drawLine(&surface, .{ .x = 1, .y = 5 }, .{ .x = 8, .y = 5 }, 1, red);
    try zigstd.testExpectTrue(surface.target.getPixel(1, 5).equals(red));
    try zigstd.testExpectTrue(surface.target.getPixel(8, 5).equals(red));
    try zigstd.testExpectTrue(surface.target.getPixel(5, 4).equals(color.RGBA.blank));
}

test "drawLine (thick) and off-screen coordinates don't panic" {
    const alloc = zigstd.testAllocator;
    var surface = try SwSurface.init(alloc, 10, 10, .{ .fill = color.RGBA.blank });
    defer surface.deinit(alloc);

    drawLine(&surface, .{ .x = -5, .y = 5 }, .{ .x = 5, .y = 5 }, 4, red);
    try zigstd.testExpectTrue(surface.target.getPixel(0, 5).equals(red));
}

test "fillCircle fills the center and clips the far corner" {
    const alloc = zigstd.testAllocator;
    var surface = try SwSurface.init(alloc, 10, 10, .{ .fill = color.RGBA.blank });
    defer surface.deinit(alloc);

    fillCircle(&surface, .{ .x = 5, .y = 5 }, 4, red);
    try zigstd.testExpectTrue(surface.target.getPixel(5, 5).equals(red));
    try zigstd.testExpectTrue(surface.target.getPixel(0, 0).equals(color.RGBA.blank));
}

test "drawText blits the VGA font's 'A' glyph" {
    const alloc = zigstd.testAllocator;
    var surface = try SwSurface.init(alloc, 16, 16, .{ .fill = color.RGBA.blank });
    defer surface.deinit(alloc);

    drawText(&surface, .{ .x = 0, .y = 0 }, "A", red);
    // Row 2 of 'A' (0-indexed) is 0x10 = 0b0001_0000 -> column 3 lit.
    try zigstd.testExpectTrue(surface.target.getPixel(3, 2).equals(red));
    try zigstd.testExpectTrue(surface.target.getPixel(0, 2).equals(color.RGBA.blank));
}

test "drawImage blits a source RgbaImage's pixels" {
    const alloc = zigstd.testAllocator;
    var surface = try SwSurface.init(alloc, 10, 10, .{ .fill = color.RGBA.blank });
    defer surface.deinit(alloc);

    var src = try RgbaImage.init(alloc, 2, 2, color.RGBA.blank);
    defer src.deinit(alloc);
    src.setPixel(0, 0, red);
    src.setPixel(1, 1, red);

    drawImage(&surface, .{ .x = 3, .y = 3 }, &src);
    try zigstd.testExpectTrue(surface.target.getPixel(3, 3).equals(red));
    try zigstd.testExpectTrue(surface.target.getPixel(4, 4).equals(red));
    try zigstd.testExpectTrue(surface.target.getPixel(4, 3).equals(color.RGBA.blank));
}
