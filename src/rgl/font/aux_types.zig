const types = @import("stream_types.zig");

const FWord = types.FWord;
const UFWord = types.UFWord;

/// A glyph's own bounding box, in design units. Wider than FWord (glyf's
/// on-disk field width) since CFF's bounds come from actually walking a
/// Type 2 charstring's accumulated coordinates, which aren't spec-bounded
/// to i16 the way glyf's stored header fields are - the generic type should
/// be able to hold whatever either format actually produces, not force a
/// narrowing cast onto the wider one.
pub const GlyphBoundingBox = struct {
    x_min: i32,
    y_min: i32,
    x_max: i32,
    y_max: i32,
};

pub const GlyphMetrics = struct {
    advance_width: UFWord,
    /// The glyph outline's own bounding box - not necessarily identical to
    /// hmtx's left_side_bearing (though usually equal in practice).
    bearing_x: FWord,
    bearing_y: FWord,
};

pub const GlyphBitmap = struct {
    width: u16,
    height: u16,
    /// Offset in pixel space from the glyph origin to the left of the bitmap.
    off_x: i16,
    /// Offset in pixel space from the glyph origin to the top of the bitmap.
    off_y: i16,

    pub const empty: GlyphBitmap = .{
        .width = 0,
        .height = 0,
        .off_x = 0,
        .off_y = 0,
    };
};

/// The original TrueType `VerticalMetrics` (ascent/descent/line_gap),
/// augmented with OS/2-typo-metrics awareness and underline info.
pub const VMetrics = struct {
    /// hhea's Ascender+|Descender|+LineGap, or OS/2's sTypoAscender/
    /// sTypoDescender/sTypoLineGap if fsSelection.use_typo_metrics is set.
    line_height: FWord = 0,
    ascender: FWord = 0,
    descender: FWord = 0,
    underline_position: FWord = 0,
    underline_thickness: FWord = 0,
};

/// One move/line/curve command in a decomposed glyph outline (glyf's
/// vmove/vline/vcurve, or CFF's vcubic). Coordinates are in design units.
pub const Vertex = struct {
    const Self = @This();

    x: i16,
    y: i16,
    cx: i16,
    cy: i16,
    cx1: i16,
    cy1: i16,
    type: Type,

    pub const Type = enum(u8) {
        vmove = 1,
        vline = 2,
        vcurve = 3,
        vcubic = 4,
        _,
    };

    pub fn set(v: *Self, ty: Type, x: i32, y: i32, cx: i32, cy: i32) void {
        v.type = ty;
        v.x = @intCast(x);
        v.y = @intCast(y);
        v.cx = @intCast(cx);
        v.cy = @intCast(cy);
    }
};
