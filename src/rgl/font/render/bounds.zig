const zigstd = @import("../../../zigstd.zig");
const aux_types = @import("../aux_types.zig");

/// Computes the pixel-space bounding box an outline's own bounding box
/// (in outline-space units) covers at the given scale/subpixel shift,
/// rounded outward to whole pixels. Y is flipped: outline-space is Y-up,
/// bitmap rows are Y-down.
pub const BitmapBox = struct {
    x0: i32,
    y0: i32,
    x1: i32,
    y1: i32,

    pub const empty: BitmapBox = .{ .x0 = 0, .y0 = 0, .x1 = 0, .y1 = 0 };
};

pub fn glyphBitmapBoxSubpixel(
    bbox: aux_types.GlyphBoundingBox,
    scale_x: f32,
    scale_y: f32,
    shift_x: f32,
    shift_y: f32,
) BitmapBox {
    return .{
        // move to integral bboxes (treating pixels as little squares, what pixels get touched)?
        .x0 = @floor(@as(f32, @floatFromInt(bbox.x_min)) * scale_x + shift_x),
        .y0 = @floor(@as(f32, @floatFromInt(-bbox.y_max)) * scale_y + shift_y),
        .x1 = @ceil(@as(f32, @floatFromInt(bbox.x_max)) * scale_x + shift_x),
        .y1 = @ceil(@as(f32, @floatFromInt(-bbox.y_min)) * scale_y + shift_y),
    };
}

pub fn glyphBitmapBox(bbox: aux_types.GlyphBoundingBox, scale_x: f32, scale_y: f32) BitmapBox {
    return glyphBitmapBoxSubpixel(bbox, scale_x, scale_y, 0, 0);
}

test "glyphBitmapBoxSubpixel scales and flips a glyph bbox into device space" {
    const bbox: aux_types.GlyphBoundingBox = .{ .x_min = 0, .y_min = -10, .x_max = 100, .y_max = 200 };
    const box = glyphBitmapBoxSubpixel(bbox, 0.5, 0.5, 0, 0);
    try zigstd.testExpectEqual(@as(i32, 0), box.x0);
    try zigstd.testExpectEqual(@as(i32, -100), box.y0);
    try zigstd.testExpectEqual(@as(i32, 50), box.x1);
    try zigstd.testExpectEqual(@as(i32, 5), box.y1);
}
