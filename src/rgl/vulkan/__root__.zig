//! Vulkan backend: batches draw calls into per-pipeline instance data
//! instead of executing them. Nothing here touches the Vulkan API (no
//! `@import` of the vulkan-zig bindings) -- that keeps this portable and
//! safe to compile into any rstd consumer, including headless/CLI ones
//! and, eventually, the kernel itself.
//!
//! Pipelines are not limited to rgl's own three (rect/glyph/image): a
//! `Batch` never knows or cares what instance struct it's holding, only
//! its byte size (`stride`) -- an app's own custom shader (e.g. rgui's
//! `ui.zig`, a faster non-SDF path for plain axis-aligned rects) batches
//! through the exact same `drawCustom`/`Batch` machinery the built-ins
//! use, just under a `PipelineId` of its own choosing. rgl.zig's
//! SDF-based rect/roundrect/circle/line pipeline is the one guaranteed
//! fallback -- correct for any shape render.zig's box/line primitives can
//! describe, at the cost of being the most expensive per-pixel. Nothing
//! here picks a faster path automatically; that's left to whoever's
//! building on top (see rgui's eventual widget layer), which knows what
//! custom pipelines are actually registered.
//!
//! `endFrame` is therefore a no-op: the batches just sit here until
//! whoever owns the actual device, swapchain and command buffers (nobody
//! yet -- see `gpu.zig`, gated behind `-Dvulkan`) uploads and records
//! them, then calls `clearBatches`.

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
/// with N instances) to mean anything. Kept as its own field rather than
/// going through the generic `Batch` map for exactly that reason.
pub const ImageDraw = struct {
    image: *const Image,
    instance: ImageInstance,
};

/// Identifies which pipeline a batch of instances belongs to. rgl reserves
/// a handful of low values for its own built-in pipelines; anything else
/// (>= `custom_pipeline_id_start`) is free for an app or a custom shader
/// to claim via `drawCustom` -- there's no registration step here, since
/// `Batch` is created on first use and doesn't need to know anything about
/// the pipeline beyond its instance stride.
pub const PipelineId = u32;
pub const rect_pipeline_id: PipelineId = 0;
pub const glyph_pipeline_id: PipelineId = 1;
pub const image_pipeline_id: PipelineId = 2;
pub const custom_pipeline_id_start: PipelineId = 16;

/// Raw accumulated instance data for one pipeline: `stride` apart, however
/// many bytes long. Never knows the concrete instance struct -- that's
/// what lets a custom shader's instances batch exactly like a built-in
/// one's, and lets gpu.zig upload/draw any of them with one function.
pub const Batch = struct {
    stride: usize = 0,
    bytes: std.ArrayList(u8) = .empty,

    pub fn count(self: Batch) usize {
        if (self.stride == 0) return 0;
        return self.bytes.items.len / self.stride;
    }

    fn append(self: *Batch, alloc: Allocator, comptime T: type, instance: T) !void {
        zigstd.assert(self.stride == 0 or self.stride == @sizeOf(T));
        self.stride = @sizeOf(T);
        try self.bytes.appendSlice(alloc, std.mem.asBytes(&instance));
    }

    /// Reinterprets the accumulated bytes as `[]const T` -- `T` must be the
    /// same type (or at least the same size) every caller of `append` used.
    pub fn items(self: Batch, comptime T: type) []const T {
        zigstd.assert(@sizeOf(T) == self.stride);
        const ptr: [*]const T = @ptrCast(@alignCast(self.bytes.items.ptr));
        return ptr[0..self.count()];
    }

    fn clear(self: *Batch) void {
        self.bytes.clearRetainingCapacity();
    }

    fn deinit(self: *Batch, alloc: Allocator) void {
        self.bytes.deinit(alloc);
    }
};

const BatchEntry = struct { id: PipelineId, batch: Batch };

pub const InitOptions = struct {};

pub const VkSurface = struct {
    const Self = @This();

    alloc: Allocator,
    clear_color: ?RGBA = null,
    /// Small in practice (rect/glyph/image plus however many custom
    /// pipelines an app registers, not per-frame growth), so a linear
    /// scan beats pulling in a hash map for this.
    batches: std.ArrayList(BatchEntry) = .empty,
    images: std.ArrayList(ImageDraw) = .empty,

    pub fn init(alloc: Allocator, width: u32, height: u32, options: InitOptions) !Self {
        _ = width;
        _ = height;
        _ = options;
        return .{ .alloc = alloc };
    }

    pub fn destroy(self: *Self, alloc: Allocator) void {
        for (self.batches.items) |*entry| entry.batch.deinit(alloc);
        self.batches.deinit(alloc);
        self.images.deinit(alloc);
    }

    fn batchFor(self: *Self, id: PipelineId) !*Batch {
        for (self.batches.items) |*entry| {
            if (entry.id == id) return &entry.batch;
        }
        try self.batches.append(self.alloc, .{ .id = id, .batch = .{} });
        return &self.batches.items[self.batches.items.len - 1].batch;
    }

    fn appendInstance(self: *Self, id: PipelineId, comptime T: type, instance: T) !void {
        const b = try self.batchFor(id);
        try b.append(self.alloc, T, instance);
    }

    /// Read-only lookup for gpu.zig (or any other consumer) to record a
    /// pipeline's batch -- `null` if nothing's ever been drawn with `id`.
    pub fn batch(self: *const Self, id: PipelineId) ?*const Batch {
        for (self.batches.items) |*entry| {
            if (entry.id == id) return &entry.batch;
        }
        return null;
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
/// once it's safe to reuse them for the next frame.
pub fn clearBatches(surface: *VkSurface) void {
    for (surface.batches.items) |*entry| entry.batch.clear();
    surface.images.clearRetainingCapacity();
}

/// Batches one instance of an app- or custom-shader-owned pipeline, using
/// whatever instance struct that shader's own vertex attributes expect.
/// `pipeline` should be >= `custom_pipeline_id_start` to avoid colliding
/// with rgl's own built-in pipelines.
pub fn drawCustom(surface: *VkSurface, pipeline: PipelineId, comptime Instance: type, instance: Instance) !void {
    zigstd.assert(pipeline >= custom_pipeline_id_start);
    try surface.appendInstance(pipeline, Instance, instance);
}

pub fn drawLine(surface: *VkSurface, p0: PointF32, p1: PointF32, width: f32, c: RGBA) !void {
    const d = p1.sub(p0);
    const length = @sqrt(d.x * d.x + d.y * d.y);
    const mid = p0.add(p1).scaled(0.5);
    // A zero-length line degenerates to a dot: draw it as a small square
    // rather than a zero-area rect that'd rasterize to nothing.
    const size: PointF32 = if (length == 0) .{ .x = width, .y = width } else .{ .x = length, .y = width };
    const rotation: f32 = if (length == 0) 0 else std.math.atan2(d.y, d.x);
    try surface.appendInstance(rect_pipeline_id, RectInstance, .{
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
    try surface.appendInstance(rect_pipeline_id, RectInstance, rectInstance(rect, c));
}

pub fn strokeRect(surface: *VkSurface, rect: Rect, width: f32, c: RGBA) !void {
    var inst = rectInstance(rect, .{ .r = 0, .g = 0, .b = 0, .a = 0 });
    inst.border_width = width;
    inst.border_color = rgbaToVec4(c);
    try surface.appendInstance(rect_pipeline_id, RectInstance, inst);
}

pub fn fillRoundRect(surface: *VkSurface, rect: Rect, radius: f32, c: RGBA) !void {
    var inst = rectInstance(rect, c);
    inst.corner_radius = radius;
    try surface.appendInstance(rect_pipeline_id, RectInstance, inst);
}

pub fn fillCircle(surface: *VkSurface, center: PointF32, radius: f32, c: RGBA) !void {
    try surface.appendInstance(rect_pipeline_id, RectInstance, .{
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
        try surface.appendInstance(glyph_pipeline_id, GlyphInstance, .{
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
    const rects = surface.batch(rect_pipeline_id).?;
    try zigstd.testExpectEqual(@as(usize, 1), rects.count());
    const inst = rects.items(RectInstance)[0];
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
    try zigstd.testExpectEqual(@as(usize, 3), surface.batch(rect_pipeline_id).?.count());

    try drawText(&surface, .{ .x = 0, .y = 0 }, "hi", .{ .r = 255, .g = 255, .b = 255, .a = 255 });
    try zigstd.testExpectEqual(@as(usize, 2), surface.batch(glyph_pipeline_id).?.count());

    var img = try Image.init(alloc, 2, 2, color.RGBA.blank);
    defer img.destroy(alloc);
    try drawImage(&surface, .{ .x = 0, .y = 0 }, &img);
    try zigstd.testExpectEqual(@as(usize, 1), surface.images.items.len);

    clearBatches(&surface);
    try zigstd.testExpectEqual(@as(usize, 0), surface.batch(rect_pipeline_id).?.count());
    try zigstd.testExpectEqual(@as(usize, 0), surface.batch(glyph_pipeline_id).?.count());
    try zigstd.testExpectEqual(@as(usize, 0), surface.images.items.len);
}

test "drawCustom batches a caller-defined instance struct under its own pipeline id" {
    const alloc = zigstd.testAllocator;
    var surface = try VkSurface.init(alloc, 10, 10, .{});
    defer surface.destroy(alloc);

    const MyInstance = extern struct { pos: [2]f32, radius: f32 };
    const my_pipeline: PipelineId = custom_pipeline_id_start;

    try drawCustom(&surface, my_pipeline, MyInstance, .{ .pos = .{ 1, 2 }, .radius = 3 });
    try drawCustom(&surface, my_pipeline, MyInstance, .{ .pos = .{ 4, 5 }, .radius = 6 });

    const custom = surface.batch(my_pipeline).?;
    try zigstd.testExpectEqual(@as(usize, 2), custom.count());
    try zigstd.testExpectEqual(@as(f32, 6), custom.items(MyInstance)[1].radius);
    // Untouched pipelines stay absent rather than springing into existence.
    try zigstd.testExpectTrue(surface.batch(rect_pipeline_id) == null);
}
