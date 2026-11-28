//! A variable glyph's advance-width delta at normalized coordinates
//! (OpenType 1.9.1, "HVAR", and "gvar", "Phantom points"): HVAR's when the
//! font has one, else how far `gvar` moves the second phantom point
//! against the first.

const bytes = @import("../../../bytes.zig");
const zigstd = @import("../../../zigstd.zig");
const types = @import("types.zig");
const tables = @import("tables/__root__.zig");
const Sfnt = @import("Sfnt.zig");
const gvar_example = @import("tables/_test/gvar_example.zig");
const CompositeGlyph = @import("tables/glyf/CompositeGlyph.zig");

const Allocator = zigstd.Allocator;
const F2Dot14 = types.F2Dot14;
const GlyphIndex = types.GlyphIndex;
const gvar_deltas = tables.gvar_deltas;
const ivs_delta = tables.ivs_delta;

pub const Error = ivs_delta.Error || gvar_deltas.Error || bytes.BytesError || tables.GlyfTable.Error;

/// The delta to add to `glyph`'s `hmtx` advance at `coords`, unrounded; 0
/// with no coordinates, or a font with neither HVAR nor `gvar`.
pub fn advanceDelta(sfnt: Sfnt, gpa: Allocator, glyph: GlyphIndex, coords: []const F2Dot14) Error!f32 {
    if (coords.len == 0) return 0;
    const id = @backingInt(glyph);
    if (sfnt.hvar) |hvar| {
        const index = try ivs_delta.mapIndex(hvar.getAdvanceWidthMapping(), id);
        return ivs_delta.itemDelta(hvar.getItemVariationStore(), index, coords);
    }
    const gvar = sfnt.gvar orelse return 0;
    const glyf = sfnt.glyf orelse return 0;
    const loca = sfnt.loca orelse return 0;
    var decoded: ?tables.GlyfTable.Glyph = if (loca.getGlyphRange(glyph)) |r| try glyf.getGlyph(gpa, r, true) else null;
    defer if (decoded) |*d| d.deinit(gpa);
    return phantomDelta(gpa, gvar, id, coords, decoded);
}

/// The second phantom point's x delta less the first's: the advance's change.
fn phantomDelta(gpa: Allocator, gvar: tables.GvarTable, id: u16, coords: []const F2Dot14, glyph: ?tables.GlyfTable.Glyph) Error!f32 {
    var outline: gvar_deltas.Outline = .{};
    var n: usize = 0;
    if (glyph) |g| switch (g) {
        .simple => |s| {
            n = s.numPoints();
            if (n > s.x_coordinates.len or n > s.y_coordinates.len) return error.EndOfStream;
            outline = .{ .xs = s.x_coordinates[0..n], .ys = s.y_coordinates[0..n], .ends = s.end_pts_of_contours };
        },
        .composite => |c| n = c.components.len,
    };
    const d = try gpa.alloc(f32, 2 * (n + 4));
    defer gpa.free(d);
    try gvar_deltas.glyphDeltas(gpa, gvar, id, coords, outline, d[0 .. n + 4], d[n + 4 ..]);
    return d[n + 1] - d[n];
}

test "phantomDelta is the advance phantom's move against the origin's" {
    // the test gvar's glyph as four components: tuple 2 moves every point,
    // phantoms too, by x 1, so the advance keeps its width
    var buf: [128]u8 = undefined;
    const gvar = try tables.GvarTable.init(gvar_example.twoTuples(&buf), true);
    const gpa = zigstd.testAllocator;
    const none: CompositeGlyph.Component = .{ .glyph_index = .notdef, .dx = 0, .dy = 0, .xscale = 1, .scale01 = 0, .scale10 = 0, .yscale = 1, .round_xy_to_grid = false, .use_my_metrics = false, .overlap_compound = false, .scaled_component_offset = false, .unscaled_component_offset = false };
    var comps: [4]@TypeOf(none) = @splat(none);
    const header: tables.GlyfTable.GlyphHeader = .{ .number_of_contours = -1, .x_min = 0, .y_min = 0, .x_max = 0, .y_max = 0 };
    const glyph: tables.GlyfTable.Glyph = .{ .composite = .{ .header = header, .components = &comps, .instructions = &.{} } };
    try zigstd.testExpectEqual(@as(f32, 0), try phantomDelta(gpa, gvar, 0, &.{F2Dot14.fromFloat(0.5)}, glyph));
    // with no glyph data the point list is the phantoms alone, which tuple 2's 8 deltas overrun
    try zigstd.testExpectError(error.MalformedVariationData, phantomDelta(gpa, gvar, 0, &.{F2Dot14.fromFloat(0.5)}, null));
}
