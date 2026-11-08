const std = @import("std");
const bytes = @import("../../../../../bytes.zig");
const zigstd = @import("../../../../../zigstd.zig");
const types = @import("../../types.zig");
const aux_types = @import("../../aux_types.zig");
const LocaTable = @import("../loca.zig").LocaTable;
const GlyfTable = @import("GlyfTable.zig");
const SimpleGlyph = @import("SimpleGlyph.zig");
const CompositeGlyph = @import("CompositeGlyph.zig");

const Allocator = zigstd.Allocator;
const GlyphIndex = types.GlyphIndex;
const Vertex = aux_types.Vertex;
const V2f32 = aux_types.V2f32;

/// Decodes the glyph occupying `range` (via GlyfTable+CompositeGlyph's
/// checked parse), then walks its points/components into vmove/vline/vcurve
/// vertices; the same state machine as the raw-cursor upstream version,
/// just fed from the already-decoded, bounds-checked types instead of a
/// raw byte cursor.
pub fn glyphShapeTT(loca: LocaTable, glyf: GlyfTable, gpa: Allocator, glyph: GlyphIndex) (Allocator.Error || bytes.BytesError || GlyfTable.Error)![]Vertex {
    const range = loca.getGlyphRange(glyph) orelse return &.{};

    var decoded = try glyf.getGlyph(gpa, range, true);
    defer decoded.deinit(gpa);

    return switch (decoded) {
        .simple => |simple| simpleGlyphVertices(gpa, simple),
        .composite => |composite| compositeGlyphVertices(loca, glyf, gpa, composite),
    };
}

/// Reads just the 10-byte GlyphHeader directly rather than decoding the
/// full point/component arrays via GlyfTable.getGlyph(), since that's all
/// this needs. Null for a glyph with no outline (e.g. space).
pub fn glyphBoxTT(loca: LocaTable, glyf: GlyfTable, glyph: GlyphIndex) ?aux_types.GlyphBoundingBox {
    const range = loca.getGlyphRange(glyph) orelse return null;

    const header_size = @sizeOf(GlyfTable.GlyphHeader);
    if (range.end - range.start < header_size) return null;
    const header_bytes = glyf.table_data[range.start..][0..header_size];
    const header = bytes.typeBEToNative(GlyfTable.GlyphHeader, header_bytes);
    return .{ .x_min = header.x_min, .y_min = header.y_min, .x_max = header.x_max, .y_max = header.y_max };
}

fn simpleGlyphVertices(gpa: Allocator, simple: SimpleGlyph) ![]Vertex {
    var vertices: std.ArrayList(Vertex) = .empty;
    defer vertices.deinit(gpa);

    const n = simple.numPoints();
    if (n == 0) return &.{};

    var next_move: usize = 0;
    var j: usize = 0;
    var sx: i32 = 0;
    var sy: i32 = 0;
    var cx: i32 = 0;
    var cy: i32 = 0;
    var scx: i32 = 0;
    var scy: i32 = 0;
    var start_off: bool = false;
    var was_off: bool = false;

    var i: usize = 0;
    while (i < n) : (i += 1) {
        const on_curve = simple.flags[i].on_curve_point;
        const x: i32 = simple.x_coordinates[i];
        const y: i32 = simple.y_coordinates[i];

        if (next_move == i) {
            if (i != 0)
                try closeShape(&vertices, gpa, was_off, start_off, sx, sy, scx, scy, cx, cy);

            // now start the new one
            start_off = !on_curve;
            if (start_off) {
                // if we start off with an off-curve point, then when we need to find a point on the curve
                // where we can start, and we need to save some state for when we wraparound.
                scx = x;
                scy = y;
                if (i + 1 < n and simple.flags[i + 1].on_curve_point) {
                    // next point is also a curve point, so interpolate an on-point curve
                    sx = (x + simple.x_coordinates[i + 1]) >> 1;
                    sy = (y + simple.y_coordinates[i + 1]) >> 1;
                } else {
                    // otherwise just use the next point as our start point
                    sx = if (i + 1 < n) simple.x_coordinates[i + 1] else x;
                    sy = if (i + 1 < n) simple.y_coordinates[i + 1] else y;
                    i += 1; // we're using point i+1 as the starting point, so skip it
                }
            } else {
                sx = x;
                sy = y;
            }
            try appendVertex(&vertices, gpa, .vmove, sx, sy, 0, 0);
            was_off = false;
            next_move = @as(usize, simple.end_pts_of_contours[j]) + 1;
            j += 1;
        } else {
            if (!on_curve) { // if it's a curve
                if (was_off) {
                    // two off-curve control points in a row means interpolate an on-curve midpoint
                    try appendVertex(&vertices, gpa, .vcurve, (cx + x) >> 1, (cy + y) >> 1, cx, cy);
                }
                cx = x;
                cy = y;
                was_off = true;
            } else {
                if (was_off)
                    try appendVertex(&vertices, gpa, .vcurve, x, y, cx, cy)
                else
                    try appendVertex(&vertices, gpa, .vline, x, y, 0, 0);
                was_off = false;
            }
        }
    }
    try closeShape(&vertices, gpa, was_off, start_off, sx, sy, scx, scy, cx, cy);

    return vertices.toOwnedSlice(gpa);
}

fn compositeGlyphVertices(loca: LocaTable, glyf: GlyfTable, gpa: Allocator, composite: CompositeGlyph) ![]Vertex {
    var vertices: std.ArrayList(Vertex) = .empty;
    defer vertices.deinit(gpa);

    for (composite.components) |c| {
        // x' = xscale*x + scale10*y + dx; y' = scale01*x + yscale*y + dy
        //
        // Upstream (stb_truetype) additionally multiplies this whole sum by
        // an extra magnitude scale m=sqrt(xscale^2+scale01^2) (and n for y),
        // citing Apple's composite-glyph spec. For a plain axis-aligned
        // scale (scale01=scale10=0) that squares the scale actually applied
        // to x/y (xscale=0.5 ends up scaling by 0.25) while dx/dy
        // only pick up one factor of it. That's wrong: the spec's own
        // formula (used here) is just the plain affine transform, no extra
        // magnitude multiply.
        const row_x: V2f32 = .{ c.xscale, c.scale10 };
        const row_y: V2f32 = .{ c.scale01, c.yscale };
        const translate: V2f32 = .{ c.dx, c.dy };

        // Get indexed glyph.
        const comp_verts = try glyphShapeTT(loca, glyf, gpa, c.glyph_index);
        defer gpa.free(comp_verts);
        if (comp_verts.len > 0) {
            // Transform vertices.
            for (comp_verts) |*v| {
                const p0 = transformPoint(row_x, row_y, translate, v.x, v.y);
                v.x = p0[0];
                v.y = p0[1];
                const p1 = transformPoint(row_x, row_y, translate, v.cx, v.cy);
                v.cx = p1[0];
                v.cy = p1[1];
            }
            try vertices.appendSlice(gpa, comp_verts);
        }
    }
    return vertices.toOwnedSlice(gpa);
}

/// 2D affine transform of one (x, y) point via two matrix rows + a
/// translation, as `V2f32` dot products; mirrors how mist's own
/// EdgeSegment/Shape geometry represents points, so composite-glyph
/// transforms and the MSDF math read the same way.
fn transformPoint(row_x: V2f32, row_y: V2f32, translate: V2f32, x: i16, y: i16) [2]i16 {
    const p: V2f32 = .{ @floatFromInt(x), @floatFromInt(y) };
    const tx = @reduce(.Add, row_x * p) + translate[0];
    const ty = @reduce(.Add, row_y * p) + translate[1];
    return .{ @trunc(tx), @trunc(ty) };
}

fn closeShape(
    vertices: *std.ArrayList(Vertex),
    gpa: Allocator,
    was_off: bool,
    start_off: bool,
    sx: i32,
    sy: i32,
    scx: i32,
    scy: i32,
    cx: i32,
    cy: i32,
) !void {
    if (start_off) {
        if (was_off) {
            try appendVertex(vertices, gpa, .vcurve, (cx + scx) >> 1, (cy + scy) >> 1, cx, cy);
        }
        try appendVertex(vertices, gpa, .vcurve, sx, sy, scx, scy);
    } else {
        if (was_off)
            try appendVertex(vertices, gpa, .vcurve, sx, sy, cx, cy)
        else
            try appendVertex(vertices, gpa, .vline, sx, sy, 0, 0);
    }
}

fn appendVertex(vertices: *std.ArrayList(Vertex), gpa: Allocator, ty: Vertex.Type, x: i32, y: i32, cx: i32, cy: i32) Allocator.Error!void {
    try vertices.append(gpa, undefined);
    vertices.items[vertices.items.len - 1].set(ty, x, y, cx, cy);
}
