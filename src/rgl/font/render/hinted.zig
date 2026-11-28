//! Hinted glyph bitmaps: a hinted outline (`hint.Font`, 26.6 pixels) drawn
//! by the coverage rasterizer, in a box sized as FreeType 2.14.3 sizes it
//! for a gray target (`src/base/ftobjs.c` `ft_glyphslot_preset_bitmap`:
//! the control box's floor and ceiling), the outline walked as
//! `FT_Outline_Decompose` (`src/base/ftoutln.c`) walks it.

const std = @import("std");
const aux_types = @import("../sfnt/aux_types.zig");
const hint = @import("../hint/__root__.zig");
const rasterize = @import("rasterize.zig");

const Allocator = std.mem.Allocator;
const Vertex = aux_types.Vertex;
const Coverage = @import("../../image.zig").Image(u8);

pub const Error = error{ GlyphTooLarge, InvalidOutline } || hint.Font.Error;

pub const HintedBitmap = struct {
    bitmap: aux_types.GlyphBitmap,
    /// 26.6 pixels: whole pixels unless unhinted.
    advance: i32,
    lsb_delta: i32,
    rsb_delta: i32,
};

const tag_on = 1;
const tag_cubic = 2;

fn vertex(kind: Vertex.Type, p: [2]i32, c0: [2]i32, c1: [2]i32) Error!Vertex {
    const fits = struct {
        fn f(v: i32) Error!i16 {
            return std.math.cast(i16, v) orelse error.GlyphTooLarge;
        }
    }.f;
    return .{ .type = kind, .x = try fits(p[0]), .y = try fits(p[1]), .cx = try fits(c0[0]), .cy = try fits(c0[1]), .cx1 = try fits(c1[0]), .cy1 = try fits(c1[1]) };
}

fn mid(a: [2]i32, b: [2]i32) [2]i32 {
    return .{ @divFloor(a[0] + b[0], 2), @divFloor(a[1] + b[1], 2) };
}

/// The outline as moves, lines, conics and cubics, in 26.6 units; the
/// caller owns it. `FT_Outline_Decompose`: a contour starting off the
/// curve starts at its last point, or between its first and last.
pub fn toVertices(allocator: Allocator, g: hint.Glyph) Error![]Vertex {
    var out: std.ArrayList(Vertex) = .empty;
    errdefer out.deinit(allocator);
    var first: usize = 0;
    for (g.contours) |end_u16| {
        const end: usize = end_u16;
        defer first = end + 1;
        if (end < first) continue;
        const pts = g.points[first .. end + 1];
        const tags = g.tags[first .. end + 1];
        var limit = pts.len;
        var at: usize = 0;
        var start = pts[0];
        if (tags[0] & 3 != tag_on) {
            // FreeType rejects a contour starting on a cubic point
            if (tags[0] & 3 == tag_cubic) continue;
            if (tags[limit - 1] & 3 == tag_on) {
                start = pts[limit - 1];
                limit -= 1;
            } else start = mid(pts[0], pts[limit - 1]);
            // the first point is a control point: keep it
            at = 0;
        } else at = 1;
        try out.append(allocator, try vertex(.vmove, start, .{ 0, 0 }, .{ 0, 0 }));
        var i = at;
        while (i < limit) {
            const p = pts[i];
            switch (tags[i] & 3) {
                tag_on => {
                    try out.append(allocator, try vertex(.vline, p, .{ 0, 0 }, .{ 0, 0 }));
                    i += 1;
                },
                tag_cubic => {
                    // two cubic controls, then an on point or the start
                    if (i + 1 >= limit or tags[i + 1] & 3 != tag_cubic) return error.InvalidOutline;
                    const c1 = pts[i + 1];
                    const to = if (i + 2 < limit) pts[i + 2] else start;
                    try out.append(allocator, try vertex(.vcubic, to, p, c1));
                    i += 3;
                },
                else => {
                    // conics: an implied on point between two of them
                    var control = p;
                    i += 1;
                    while (true) {
                        if (i >= limit) {
                            try out.append(allocator, try vertex(.vcurve, start, control, .{ 0, 0 }));
                            break;
                        }
                        const q = pts[i];
                        if (tags[i] & 3 == tag_on) {
                            try out.append(allocator, try vertex(.vcurve, q, control, .{ 0, 0 }));
                            i += 1;
                            break;
                        }
                        try out.append(allocator, try vertex(.vcurve, mid(control, q), control, .{ 0, 0 }));
                        control = q;
                        i += 1;
                    }
                    continue;
                },
            }
        }
        // back to the start
        try out.append(allocator, try vertex(.vline, start, .{ 0, 0 }, .{ 0, 0 }));
    }
    return out.toOwnedSlice(allocator);
}

/// Glyph `gid` hinted at `ppem` per `opts`, its coverage appended to
/// `pixels` (left to right, top to bottom, 0 transparent).
pub fn glyphBitmapHinted(font: *hint.Font, gpa: Allocator, pixels: *std.ArrayList(u8), gid: u16, ppem: u16, opts: hint.Font.Options) Error!HintedBitmap {
    var g = try font.load(gpa, gid, ppem, opts);
    defer g.deinit(gpa);
    var r: HintedBitmap = .{ .bitmap = .empty, .advance = g.advance, .lsb_delta = g.lsb_delta, .rsb_delta = g.rsb_delta };
    if (g.points.len == 0) return r;
    var lo: [2]i32 = g.points[0];
    var hi: [2]i32 = g.points[0];
    for (g.points) |p| for (0..2) |k| {
        lo[k] = @min(lo[k], p[k]);
        hi[k] = @max(hi[k], p[k]);
    };
    const x0 = lo[0] >> 6;
    const y0 = lo[1] >> 6;
    const x1 = (hi[0] + 63) >> 6;
    const y1 = (hi[1] + 63) >> 6;
    const w: u32 = @intCast(x1 - x0);
    const h: u32 = @intCast(y1 - y0);
    if (w == 0 or h == 0) return r;
    const vertices = try toVertices(gpa, g);
    defer gpa.free(vertices);
    const bmp: Coverage = .{ .width = w, .height = h, .data = try pixels.addManyAsSlice(gpa, w * h) };
    errdefer pixels.shrinkRetainingCapacity(pixels.items.len - bmp.data.len);
    try rasterize.rasterize(gpa, bmp, vertices, .{
        .flatness_in_pixels = 0.35,
        .scale_x = 1.0 / 64.0,
        .scale_y = 1.0 / 64.0,
        .shift_x = 0,
        .shift_y = 0,
        .off_x = x0,
        .off_y = -y1,
        .invert = true,
    });
    r.bitmap = .{ .width = @intCast(w), .height = @intCast(h), .off_x = @intCast(x0), .off_y = @intCast(-y1) };
    return r;
}
