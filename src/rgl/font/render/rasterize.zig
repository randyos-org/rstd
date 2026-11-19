const std = @import("std");
const zigstd = @import("../../../zigstd.zig");
const Allocator = zigstd.Allocator;

const Sfnt = @import("../sfnt/Sfnt.zig");
const types = @import("../sfnt/types.zig");
const aux_types = @import("../sfnt/aux_types.zig");
const tessellate = @import("tessellate.zig");
const sweep = @import("sweep.zig");
const bounds = @import("bounds.zig");

const Vertex = aux_types.Vertex;
const GlyphIndex = types.GlyphIndex;
const Point = tessellate.Point;
const flattenCurves = tessellate.flattenCurves;
const glyphBitmapBoxSubpixel = bounds.glyphBitmapBoxSubpixel;
const Edge = sweep.Edge;
const rasterizeSortedEdges = sweep.rasterizeSortedEdges;

/// 8bpp coverage: 0 transparent, 255 opaque.
const Coverage = @import("../../image.zig").Image(u8);

/// How `rasterize` maps outline space onto its `Coverage`.
pub const RasterizeOptions = struct {
    flatness_in_pixels: f32,
    scale_x: f32,
    scale_y: f32,
    shift_x: f32,
    shift_y: f32,
    /// The bitmap's top-left corner in device space.
    off_x: i32,
    off_y: i32,
    /// Flip Y (outline space is Y-up, bitmap rows Y-down).
    invert: bool,
};

pub fn rasterize(gpa: Allocator, result: Coverage, vertices: []Vertex, o: RasterizeOptions) Allocator.Error!void {
    const scale = @min(o.scale_x, o.scale_y);
    var windings = try flattenCurves(gpa, vertices, o.flatness_in_pixels / scale);
    defer windings.deinit(gpa);
    try rasterizeInner(gpa, result, windings.points, windings.contour_lengths, o.scale_x, o.scale_y, o.shift_x, o.shift_y, o.off_x, o.off_y, o.invert);
}

/// Caller owns returned memory.
pub fn glyphBitmap(
    sfnt: Sfnt,
    gpa: Allocator,
    /// Appended to the list.
    /// Stored left-to-right, top-to-bottom. 8 bits per pixel. 0 is
    /// transparent, 255 is opaque.
    pixels: *std.ArrayList(u8),
    glyph: GlyphIndex,
    scale_x: f32,
    scale_y: f32,
) !aux_types.GlyphBitmap {
    return glyphBitmapSubpixel(sfnt, gpa, pixels, glyph, scale_x, scale_y, 0, 0);
}

/// Caller owns returned memory.
pub fn glyphBitmapSubpixel(
    sfnt: Sfnt,
    gpa: Allocator,
    /// Appended to the list.
    /// Stored left-to-right, top-to-bottom. 8 bits per pixel. 0 is
    /// transparent, 255 is opaque.
    pixels: *std.ArrayList(u8),
    glyph: GlyphIndex,
    scale_x: f32,
    scale_y: f32,
    shift_x: f32,
    shift_y: f32,
) !aux_types.GlyphBitmap {
    const vertices = try sfnt.glyphShape(gpa, glyph);
    defer gpa.free(vertices);

    std.debug.assert(scale_x != 0);
    std.debug.assert(scale_y != 0);

    const bbox = sfnt.getGlyphBoundingBox(glyph) orelse return .empty;
    const box = glyphBitmapBoxSubpixel(bbox, scale_x, scale_y, shift_x, shift_y);

    const w: u32 = @intCast(box.x1 - box.x0);
    const h: u32 = @intCast(box.y1 - box.y0);

    if (w == 0 or h == 0) return .empty;

    const gbm: Coverage = .{
        .width = w,
        .height = h,
        .data = try pixels.addManyAsSlice(gpa, w * h),
    };
    errdefer pixels.shrinkRetainingCapacity(pixels.items.len - gbm.data.len);

    try rasterize(gpa, gbm, vertices, .{
        .flatness_in_pixels = 0.35,
        .scale_x = scale_x,
        .scale_y = scale_y,
        .shift_x = shift_x,
        .shift_y = shift_y,
        .off_x = box.x0,
        .off_y = box.y0,
        .invert = true,
    });

    return .{
        .width = @intCast(gbm.width),
        .height = @intCast(gbm.height),
        .off_x = @intCast(box.x0),
        .off_y = @intCast(box.y0),
    };
}

fn rasterizeInner(
    gpa: Allocator,
    result: Coverage,
    pts: []Point,
    wcount: []u32,
    scale_x: f32,
    scale_y: f32,
    shift_x: f32,
    shift_y: f32,
    off_x: i32,
    off_y: i32,
    invert: bool,
) Allocator.Error!void {
    const y_scale_inv: f32 = if (invert) -scale_y else scale_y;

    // now we have to blow out the windings into explicit edge lists
    const edge_alloc_n = n: {
        var n: u32 = 1; // Add an extra one as a sentinel.
        for (wcount) |elem| n += elem;
        break :n n;
    };

    const e = try gpa.alloc(Edge, edge_alloc_n);
    defer gpa.free(e);

    var n: u32 = 0;
    var m: u32 = 0;
    for (wcount) |wcount_elem| {
        const p: []Point = pts[m..];
        m += wcount_elem;
        var j: u32 = wcount_elem - 1;
        var k: u32 = 0;
        while (k < wcount_elem) : ({
            j = k;
            k += 1;
        }) {
            var a = k;
            var b = j;
            // skip the edge if horizontal
            if (p[j].y == p[k].y)
                continue;
            // add edge from j to k to the list
            e[n].invert = false;
            if (if (invert) p[j].y > p[k].y else p[j].y < p[k].y) {
                e[n].invert = true;
                a = j;
                b = k;
            }
            e[n].x0 = p[a].x * scale_x + shift_x;
            e[n].y0 = (p[a].y * y_scale_inv + shift_y);
            e[n].x1 = p[b].x * scale_x + shift_x;
            e[n].y1 = (p[b].y * y_scale_inv + shift_y);
            n += 1;
        }
    }
    // now sort the edges by their highest point (should snap to integer, and then by x)
    std.mem.sortUnstable(Edge, e[0..n], Edge.Sort{}, Edge.Sort.lessThan);

    // now, traverse the scanlines and find the intersections on each scanline, use xor winding rule
    try rasterizeSortedEdges(gpa, result, e[0 .. n + 1], off_x, off_y);
}

test "rasterize fills a fully-covered square with full coverage and misses outside it" {
    const alloc = zigstd.testAllocator;
    // A 4x4 square from (2,2) to (6,6) in outline space (Y-up).
    var outline: [5]Vertex = undefined;
    outline[0].set(.vmove, 2, 2, 0, 0);
    outline[1].set(.vline, 6, 2, 0, 0);
    outline[2].set(.vline, 6, 6, 0, 0);
    outline[3].set(.vline, 2, 6, 0, 0);
    outline[4].set(.vline, 2, 2, 0, 0);

    // With invert=true, device_y = -outline_y, so the square's y in [2,6]
    // lands at device_y in [-6,-2]. off_x/off_y (this bitmap's own top-left
    // corner, in that same device space) must actually cover that range -
    // computeBitmapBox() is exactly how a real caller derives this from a
    // glyph's bounding box; here it's given directly so the test is legible.
    const off_x: i32 = 2;
    const off_y: i32 = -6;

    const pixels = try alloc.alloc(u8, 8 * 8);
    defer alloc.free(pixels);
    @memset(pixels, 0);
    const bitmap: Coverage = .{ .width = 8, .height = 8, .data = pixels };

    try rasterize(alloc, bitmap, &outline, .{ .flatness_in_pixels = 0.3, .scale_x = 1, .scale_y = 1, .shift_x = 0, .shift_y = 0, .off_x = off_x, .off_y = off_y, .invert = true });

    // Local (1,1) is device (3,-5); inside the square (outline x in [2,6],
    // device_y in [-6,-2]).
    try zigstd.testExpectEqual(@as(u8, 255), pixels[1 * 8 + 1]);
    // Local (6,6) is device (8,0); outside the square (x=8 > 6).
    try zigstd.testExpectEqual(@as(u8, 0), pixels[6 * 8 + 6]);
}
