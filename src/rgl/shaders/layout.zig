//! GPU layouts shared by the shaders and the vulkan backend: the push-constant
//! globals, `ShapeKind` and the per-instance vertex data. Field order, types
//! and alignment are the wire format the pipelines read.

/// Per-draw globals, bound as a push constant by every pipeline here. `V2` is
/// `Vec2` in the shaders (shaders/lib/__root__.zig) and `[2]f32` on the host
/// (vulkan/gpu.zig).
pub fn Globals(comptime V2: type) type {
    return extern struct {
        screen_size: V2,
        /// only meaningful to `glyph.zig` (normalizes atlas texture coordinates)
        atlas_size: V2 = .{ 0, 0 },
        /// only meaningful to `glyph.zig`'s fragment shader: 0 means the bound
        /// atlas is plain single-channel coverage (sample .r directly as
        /// alpha); a nonzero value is the msdf.Generator px_range the atlas was
        /// built with, switching to median(r,g,b) signed-distance sampling.
        glyph_msdf_px_range: f32 = 0,
    };
}

/// rgl.zig's shape, per instance.
pub const ShapeKind = enum(u32) {
    /// `corner_radius` 0 is `Rectangle`, nonzero is `RoundRect`.
    rect = 0,
    /// `Ellipse`, restricted to circles; see common.zig.
    circle = 1,
};

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
    shape_kind: u32 = @backingInt(ShapeKind.rect),
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
