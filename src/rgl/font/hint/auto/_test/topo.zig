//! Tests: the autohinter's outline analysis, segments and edges on small
//! hand-made outlines.

const std = @import("std");
const GlyphHints = @import("../GlyphHints.zig");
const edges = @import("../edges.zig");

const testing = std.testing;
const flag = GlyphHints.flag;

/// 12 ppem at 1000 units per em, in 16.16 (`FT_DivFix(12 * 64, 1000)`).
const scale: i32 = 50332;

fn analyze(h: *GlyphHints, points: []const [2]i32, contours: []const u16) !void {
    const tags = try testing.allocator.alloc(u8, points.len);
    defer testing.allocator.free(tags);
    @memset(tags, 1);
    h.x_scale = scale;
    h.y_scale = scale;
    try h.reload(.{ .points = points, .tags = tags, .contours = contours }, 1000);
    for ([_]GlyphHints.Dim{ .horz, .vert }) |dim| try edges.detectFeatures(h, dim, 0, 10, false, 1000);
}

test "autohint directions: a dominant arm or none" {
    try testing.expectEqual(GlyphHints.dir_left, GlyphHints.direction(-100, 0));
    try testing.expectEqual(GlyphHints.dir_right, GlyphHints.direction(100, 0));
    try testing.expectEqual(GlyphHints.dir_down, GlyphHints.direction(0, -100));
    try testing.expectEqual(GlyphHints.dir_up, GlyphHints.direction(7, 100));
    // within about 4.1 degrees only
    try testing.expectEqual(GlyphHints.dir_none, GlyphHints.direction(8, 100));
}

test "autohint orientation: TrueType clockwise, PostScript counterclockwise" {
    const cw = [_][2]i32{ .{ 0, 0 }, .{ 0, 700 }, .{ 500, 700 }, .{ 500, 0 } };
    const ccw = [_][2]i32{ .{ 0, 0 }, .{ 500, 0 }, .{ 500, 700 }, .{ 0, 700 } };
    const tags = [_]u8{ 1, 1, 1, 1 };
    try testing.expectEqual(GlyphHints.Orientation.truetype, GlyphHints.orientation(.{ .points = &cw, .tags = &tags, .contours = &.{3} }));
    try testing.expectEqual(GlyphHints.Orientation.postscript, GlyphHints.orientation(.{ .points = &ccw, .tags = &tags, .contours = &.{3} }));
}

test "autohint topology: a rectangle's stems pair into linked edges" {
    var h: GlyphHints = .{ .allocator = testing.allocator };
    defer h.deinit();
    // a near point on the left side is weak, the corners strong
    try analyze(&h, &.{ .{ 0, 0 }, .{ 0, 5 }, .{ 0, 700 }, .{ 500, 700 }, .{ 500, 0 } }, &.{4});
    const pts = h.points.items;
    try testing.expect(pts[1].flags & flag.weak != 0);
    for ([_]usize{ 0, 2, 3, 4 }) |i| try testing.expect(pts[i].flags & flag.weak == 0);
    try testing.expectEqual(GlyphHints.dir_up, pts[0].out_dir);
    try testing.expectEqual(@as(i32, 384), pts[3].ox); // FT_MulFix(500, scale)

    // vertical stems: the left side up, the right side down, linked both ways
    const hs = h.axis[0].segments.items;
    try testing.expectEqual(@as(usize, 2), hs.len);
    try testing.expectEqual(GlyphHints.dir_up, hs[0].dir);
    try testing.expectEqual(@as(i32, 0), hs[0].pos);
    try testing.expectEqual(@as(i32, 700), hs[0].height);
    try testing.expectEqual(@as(i32, 500), hs[1].pos);
    try testing.expectEqual(@as(u32, 1), hs[0].link);
    try testing.expectEqual(@as(u32, 0), hs[1].link);
    // 500 for the distance, 6000 * 1000 / 2048 over 700 for the overlap
    try testing.expectEqual(@as(i32, 504), hs[0].score);
    const he = h.axis[0].edges.items;
    try testing.expectEqual(@as(usize, 2), he.len);
    try testing.expectEqual(@as(i32, 0), he[0].fpos);
    try testing.expectEqual(@as(i32, 500), he[1].fpos);
    try testing.expectEqual(@as(i32, 384), he[1].opos);
    try testing.expectEqual(@as(u32, 1), he[0].link);
    try testing.expectEqual(@as(u32, 0), he[1].link);
    try testing.expectEqual(@as(u8, 0), he[0].flags);

    // horizontal stems: bottom (leftward, the major direction) and top
    const vs = h.axis[1].segments.items;
    try testing.expectEqual(@as(usize, 2), vs.len);
    try testing.expectEqual(GlyphHints.dir_right, vs[0].dir);
    try testing.expectEqual(@as(i32, 700), vs[0].pos);
    try testing.expectEqual(GlyphHints.dir_left, vs[1].dir);
    try testing.expectEqual(@as(i32, 0), vs[1].pos);
    const ve = h.axis[1].edges.items;
    try testing.expectEqual(@as(usize, 2), ve.len);
    try testing.expectEqual(@as(i32, 0), ve[0].fpos);
    try testing.expectEqual(@as(u32, 1), ve[0].link);
}

test "autohint topology: a round bowl's segments are round, a one-point contour has no direction" {
    var h: GlyphHints = .{ .allocator = testing.allocator };
    defer h.deinit();
    const tags = [_]u8{ 1, 0, 1, 0, 1, 0, 1, 0, 1 };
    // an octagon-ish bowl of conics, then a lone point
    const pts = [_][2]i32{ .{ 250, 0 }, .{ 0, 0 }, .{ 0, 350 }, .{ 0, 700 }, .{ 250, 700 }, .{ 500, 700 }, .{ 500, 350 }, .{ 500, 0 }, .{ 900, 900 } };
    h.x_scale = scale;
    h.y_scale = scale;
    try h.reload(.{ .points = &pts, .tags = &tags, .contours = &.{ 7, 8 } }, 1000);
    try edges.detectFeatures(&h, .horz, 0, 10, false, 1000);
    const hs = h.axis[0].segments.items;
    // the two sides, each a control point between on points, and the lone point
    try testing.expectEqual(@as(usize, 3), hs.len);
    try testing.expect(hs[0].flags & GlyphHints.edge_flag.round != 0);
    try testing.expectEqual(GlyphHints.dir_none, hs[2].dir);
    try testing.expectEqual(@as(i32, 900), hs[2].pos);
}
