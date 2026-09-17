//! Vulkan backend: batches draw calls into per-shape instance data instead
//! of executing them. Nothing here touches the Vulkan API (no `@import`
//! of the vulkan-zig bindings) -- that keeps this portable and safe to
//! compile into any rstd consumer, including headless/CLI ones and,
//! eventually, the kernel itself.
//!
//! `endFrame` is therefore a no-op: the batches just sit in `rects`/
//! `glyphs`/`images` until whoever owns the actual device, swapchain and
//! command buffers (nobody yet -- see `gpu.zig`, gated behind `-Dvulkan`)
//! uploads and records them, then calls `clearBatches`.

const std = @import("std");
const geometry = @import("../geometry.zig");
const Point = geometry.Point;
const PointF32 = geometry.PointF32;
const Rect = geometry.Rect;
const color = @import("../color.zig");
const RGBA = color.RGBA;
const Image = @import("../Image.zig");
const vga_font = @import("../font/vga/__root__.zig");
const zigstd = @import("../../zigstd.zig");
const Allocator = zigstd.Allocator;

/// Real vk-zig-backed pipeline/command-buffer code, kept in its own file so
/// it (and its `@import("vulkan")`) is only ever analyzed by a consumer
/// that actually references it -- see gpu.zig's own header.
pub const gpu = @import("gpu.zig");

/// Mirrors shaders/rgl.zig's `ShapeKind` -- keep the two in sync by hand,
/// since the shader can't be a normal Zig import here (SPIR-V-only decls).
pub const ShapeKind = enum(u32) { rect = 0, circle = 1 };

/// One instance of the rect/roundrect/circle/line pipeline (shaders/rgl.zig).
/// `rotation` (radians, about the rect's own center) is what makes an
/// arbitrary-angle `drawLine` possible as just a thin rotated rect.
pub const RectInstance = extern struct {
    rect_pos: [2]f32,
    rect_size: [2]f32,
    rotation: f32 = 0,
    corner_radius: f32 = 0,
    border_width: f32 = 0,
    fill_color: [4]f32,
    border_color: [4]f32 = .{ 0, 0, 0, 0 },
    shape_kind: u32 = @intFromEnum(ShapeKind.rect),
};

/// One instance of the glyph pipeline (shaders/glyph.zig).
pub const GlyphInstance = extern struct {
    glyph_pos: [2]f32,
    glyph_size: [2]f32,
    atlas_pos: [2]f32,
    color: [4]f32,
};

/// One instance of the image pipeline (shaders/image.zig).
pub const ImageInstance = extern struct {
    rect_pos: [2]f32,
    rect_size: [2]f32,
    color: [4]f32,
};

/// Unlike rect/glyph (which share one atlas or no texture at all), the
/// image pipeline binds one texture per draw -- so instances have to stay
/// grouped by source `Image` for batching (same image -> one draw call
/// with N instances) to mean anything.
pub const ImageDraw = struct {
    image: *const Image,
    instance: ImageInstance,
};

pub const InitOptions = struct {};

pub const VkSurface = struct {
    const Self = @This();

    alloc: Allocator,
    clear_color: ?RGBA = null,
    rects: std.ArrayList(RectInstance) = .empty,
    glyphs: std.ArrayList(GlyphInstance) = .empty,
    images: std.ArrayList(ImageDraw) = .empty,

    pub fn init(alloc: Allocator, width: u32, height: u32, options: InitOptions) !Self {
        _ = width;
        _ = height;
        _ = options;
        return .{ .alloc = alloc };
    }

    pub fn destroy(self: *Self, alloc: Allocator) void {
        self.rects.deinit(alloc);
        self.glyphs.deinit(alloc);
        self.images.deinit(alloc);
    }
};

fn rgbaToVec4(c: RGBA) [4]f32 {
    const f = c.toRGBAf32();
    return .{ f.r, f.g, f.b, f.a };
}

pub fn beginFrame(surface: *VkSurface, clear: ?RGBA) void {
    surface.clear_color = clear;
}

/// No GPU work happens here -- see this file's header.
pub fn endFrame(surface: *VkSurface) void {
    _ = surface;
}

/// Called by whoever just uploaded+recorded this frame's batches (gpu.zig),
/// once it's safe to reuse the arrays for the next frame.
pub fn clearBatches(surface: *VkSurface) void {
    surface.rects.clearRetainingCapacity();
    surface.glyphs.clearRetainingCapacity();
    surface.images.clearRetainingCapacity();
}

pub fn drawLine(surface: *VkSurface, p0: PointF32, p1: PointF32, width: f32, c: RGBA) !void {
    const d = p1.sub(p0);
    const length = @sqrt(d.x * d.x + d.y * d.y);
    const mid = p0.add(p1).scaled(0.5);
    // A zero-length line degenerates to a dot: draw it as a small square
    // rather than a zero-area rect that'd rasterize to nothing.
    const size: PointF32 = if (length == 0) .{ .x = width, .y = width } else .{ .x = length, .y = width };
    const rotation: f32 = if (length == 0) 0 else std.math.atan2(d.y, d.x);
    try surface.rects.append(surface.alloc, .{
        .rect_pos = .{ mid.x - size.x / 2, mid.y - size.y / 2 },
        .rect_size = .{ size.x, size.y },
        .rotation = rotation,
        .fill_color = rgbaToVec4(c),
    });
}

fn rectInstance(rect: Rect, c: RGBA) RectInstance {
    return .{
        .rect_pos = .{ @floatFromInt(rect.x), @floatFromInt(rect.y) },
        .rect_size = .{ @floatFromInt(rect.w), @floatFromInt(rect.h) },
        .fill_color = rgbaToVec4(c),
    };
}

pub fn fillRect(surface: *VkSurface, rect: Rect, c: RGBA) !void {
    try surface.rects.append(surface.alloc, rectInstance(rect, c));
}

pub fn strokeRect(surface: *VkSurface, rect: Rect, width: f32, c: RGBA) !void {
    var inst = rectInstance(rect, .{ .r = 0, .g = 0, .b = 0, .a = 0 });
    inst.border_width = width;
    inst.border_color = rgbaToVec4(c);
    try surface.rects.append(surface.alloc, inst);
}

pub fn fillRoundRect(surface: *VkSurface, rect: Rect, radius: f32, c: RGBA) !void {
    var inst = rectInstance(rect, c);
    inst.corner_radius = radius;
    try surface.rects.append(surface.alloc, inst);
}

pub fn fillCircle(surface: *VkSurface, center: PointF32, radius: f32, c: RGBA) !void {
    try surface.rects.append(surface.alloc, .{
        .rect_pos = .{ center.x - radius, center.y - radius },
        .rect_size = .{ radius * 2, radius * 2 },
        .fill_color = rgbaToVec4(c),
        .shape_kind = @intFromEnum(ShapeKind.circle),
    });
}

/// Lays out `text` on the embedded 8x16 VGA font's grid. No atlas texture
/// exists yet to actually sample from (see gpu.zig) -- this only records
/// where each glyph would sit in a future one-row, 256-glyph-wide strip of
/// that font, atlas_pos.x = code_point * glyph_width.
pub fn drawText(surface: *VkSurface, pos: Point, text: []const u8, c: RGBA) !void {
    const glyph_w: f32 = @floatFromInt(vga_font.width);
    const glyph_h: f32 = @floatFromInt(vga_font.height);
    var pen_x: f32 = @floatFromInt(pos.x);
    const pen_y: f32 = @floatFromInt(pos.y);
    const rgba = rgbaToVec4(c);
    for (text) |ch| {
        try surface.glyphs.append(surface.alloc, .{
            .glyph_pos = .{ pen_x, pen_y },
            .glyph_size = .{ glyph_w, glyph_h },
            .atlas_pos = .{ @as(f32, @floatFromInt(ch)) * glyph_w, 0 },
            .color = rgba,
        });
        pen_x += glyph_w;
    }
}

pub fn drawImage(surface: *VkSurface, pos: Point, img: *const Image) !void {
    try surface.images.append(surface.alloc, .{
        .image = img,
        .instance = .{
            .rect_pos = .{ @floatFromInt(pos.x), @floatFromInt(pos.y) },
            .rect_size = .{ @floatFromInt(img.width), @floatFromInt(img.height) },
            .color = .{ 1, 1, 1, 1 },
        },
    });
}

test "drawLine appends a rotated thin rect instance" {
    const alloc = zigstd.testAllocator;
    var surface = try VkSurface.init(alloc, 10, 10, .{});
    defer surface.destroy(alloc);

    try drawLine(&surface, .{ .x = 0, .y = 0 }, .{ .x = 10, .y = 0 }, 2, .{ .r = 255, .g = 0, .b = 0, .a = 255 });
    try zigstd.testExpectEqual(@as(usize, 1), surface.rects.items.len);
    const inst = surface.rects.items[0];
    try zigstd.testExpectEqual(@as(f32, 0), inst.rotation);
    try zigstd.testExpectEqual([2]f32{ 10, 2 }, inst.rect_size);
}

test "fillRect/strokeRect/fillCircle/drawText/drawImage batch without touching a GPU" {
    const alloc = zigstd.testAllocator;
    var surface = try VkSurface.init(alloc, 10, 10, .{});
    defer surface.destroy(alloc);

    try fillRect(&surface, .{ .x = 0, .y = 0, .w = 5, .h = 5 }, .{ .r = 1, .g = 2, .b = 3, .a = 255 });
    try strokeRect(&surface, .{ .x = 0, .y = 0, .w = 5, .h = 5 }, 1, .{ .r = 1, .g = 2, .b = 3, .a = 255 });
    try fillCircle(&surface, .{ .x = 5, .y = 5 }, 3, .{ .r = 1, .g = 2, .b = 3, .a = 255 });
    try zigstd.testExpectEqual(@as(usize, 3), surface.rects.items.len);

    try drawText(&surface, .{ .x = 0, .y = 0 }, "hi", .{ .r = 255, .g = 255, .b = 255, .a = 255 });
    try zigstd.testExpectEqual(@as(usize, 2), surface.glyphs.items.len);

    var img = try Image.init(alloc, 2, 2, color.RGBA.blank);
    defer img.destroy(alloc);
    try drawImage(&surface, .{ .x = 0, .y = 0 }, &img);
    try zigstd.testExpectEqual(@as(usize, 1), surface.images.items.len);

    clearBatches(&surface);
    try zigstd.testExpectEqual(@as(usize, 0), surface.rects.items.len);
    try zigstd.testExpectEqual(@as(usize, 0), surface.glyphs.items.len);
    try zigstd.testExpectEqual(@as(usize, 0), surface.images.items.len);
}
