//! A glyph a hinter loaded, its slices owned by the caller's allocator.

const std = @import("std");

const Glyph = @This();

/// 26.6, the origin at the glyph's origin (TrueType's left phantom point).
points: [][2]i32,
/// FreeType's curve tags: on the curve, a conic or a cubic control point.
tags: []u8,
/// Each contour's last point.
contours: []u16,
/// 26.6, whole pixels when hinted.
advance: i32,

pub fn deinit(g: *Glyph, allocator: std.mem.Allocator) void {
    allocator.free(g.contours);
    allocator.free(g.tags);
    allocator.free(g.points);
}
