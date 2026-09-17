//! The public, backend-agnostic drawing API -- gdi32-style primitive calls
//! (line/box/text/...) that dispatch to whichever backend `surface` was
//! created with. Callers never see `sw`/`vulkan` directly.
//!
//! sw executes every call immediately, into its own pixel buffer. vulkan
//! only records instance data (see rgl/vulkan/__root__.zig's header) --
//! that's the only reason these return `!void` rather than `void`: sw
//! itself can't fail, but vulkan's append can (allocation failure).

const Surface = @import("Surface.zig");
const geometry = @import("geometry.zig");
pub const Point = geometry.Point;
pub const PointF32 = geometry.PointF32;
pub const PointI32 = geometry.PointI32;
pub const Rect = geometry.Rect;
const color = @import("color.zig");
const RGBA = color.RGBA;
const Image = @import("Image.zig");
const sw = @import("sw/__root__.zig");
const vulkan = @import("vulkan/__root__.zig");

/// Starts a new frame. `clear`, if given, replaces every pixel (sw) or
/// discards every batched instance (vulkan) left over from the last one.
pub fn beginFrame(surface: *Surface, clear: ?RGBA) void {
    switch (surface.impl) {
        .sw => |*s| sw.beginFrame(s, clear),
        .vulkan => |*v| vulkan.beginFrame(v, clear),
    }
}

/// Ends the current frame. For sw this is a no-op (every pixel is already
/// final); for vulkan it's still a no-op today (see vulkan/__root__.zig)
/// since nothing in rstd yet owns a device to submit the batch to.
pub fn endFrame(surface: *Surface) void {
    switch (surface.impl) {
        .sw => |*s| sw.endFrame(s),
        .vulkan => |*v| vulkan.endFrame(v),
    }
}

pub fn drawLine(surface: *Surface, p0: PointF32, p1: PointF32, width: f32, c: RGBA) !void {
    switch (surface.impl) {
        .sw => |*s| sw.drawLine(s, p0, p1, width, c),
        .vulkan => |*v| try vulkan.drawLine(v, p0, p1, width, c),
    }
}

pub fn fillRect(surface: *Surface, rect: Rect, c: RGBA) !void {
    switch (surface.impl) {
        .sw => |*s| sw.fillRect(s, rect, c),
        .vulkan => |*v| try vulkan.fillRect(v, rect, c),
    }
}

pub fn strokeRect(surface: *Surface, rect: Rect, width: f32, c: RGBA) !void {
    switch (surface.impl) {
        .sw => |*s| sw.strokeRect(s, rect, width, c),
        .vulkan => |*v| try vulkan.strokeRect(v, rect, width, c),
    }
}

pub fn fillRoundRect(surface: *Surface, rect: Rect, radius: f32, c: RGBA) !void {
    switch (surface.impl) {
        .sw => |*s| sw.fillRoundRect(s, rect, radius, c),
        .vulkan => |*v| try vulkan.fillRoundRect(v, rect, radius, c),
    }
}

pub fn fillCircle(surface: *Surface, center: PointF32, radius: f32, c: RGBA) !void {
    switch (surface.impl) {
        .sw => |*s| sw.fillCircle(s, center, radius, c),
        .vulkan => |*v| try vulkan.fillCircle(v, center, radius, c),
    }
}

pub fn drawText(surface: *Surface, pos: Point, text: []const u8, c: RGBA) !void {
    switch (surface.impl) {
        .sw => |*s| sw.drawText(s, pos, text, c),
        .vulkan => |*v| try vulkan.drawText(v, pos, text, c),
    }
}

pub fn drawImage(surface: *Surface, pos: Point, img: *const Image) !void {
    switch (surface.impl) {
        .sw => |*s| sw.drawImage(s, pos, img),
        .vulkan => |*v| try vulkan.drawImage(v, pos, img),
    }
}

test "the same call sequence runs on both backends" {
    const std = @import("std");
    const zigstd = @import("../zigstd.zig");
    const alloc = zigstd.testAllocator;

    const options_list = [_]Surface.InitOptions{
        .{ .sw = .{ .fill = color.RGBA.blank } },
        .{ .vulkan = .{} },
    };
    for (options_list) |options| {
        var surface = try Surface.init(alloc, 16, 16, options);
        defer surface.destroy(alloc);

        beginFrame(&surface, color.RGBA.blank);
        try fillRect(&surface, .{ .x = 0, .y = 0, .w = 4, .h = 4 }, .{ .r = 255, .g = 0, .b = 0, .a = 255 });
        try strokeRect(&surface, .{ .x = 4, .y = 4, .w = 4, .h = 4 }, 1, .{ .r = 0, .g = 255, .b = 0, .a = 255 });
        try fillRoundRect(&surface, .{ .x = 8, .y = 0, .w = 4, .h = 4 }, 1, .{ .r = 0, .g = 0, .b = 255, .a = 255 });
        try fillCircle(&surface, .{ .x = 12, .y = 12 }, 2, .{ .r = 255, .g = 255, .b = 0, .a = 255 });
        try drawLine(&surface, .{ .x = 0, .y = 15 }, .{ .x = 15, .y = 15 }, 1, .{ .r = 255, .g = 255, .b = 255, .a = 255 });
        try drawText(&surface, .{ .x = 0, .y = 8 }, "hi", .{ .r = 255, .g = 255, .b = 255, .a = 255 });

        var img = try Image.init(alloc, 2, 2, color.RGBA.blank);
        defer img.destroy(alloc);
        try drawImage(&surface, .{ .x = 0, .y = 0 }, &img);

        endFrame(&surface);
        std.mem.doNotOptimizeAway(&surface);
    }
}
