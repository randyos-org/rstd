//! Tests: walking a hinted outline into vertices as `FT_Outline_Decompose`.

const std = @import("std");
const hinted = @import("../hinted.zig");
const Glyph = @import("../../hint/Glyph.zig");

const testing = std.testing;

test "hinted outline: a contour starting off the curve, conics with implied on points" {
    // off, off, on, off: starts between the last and first (both off)
    var points = [_][2]i32{ .{ 0, 0 }, .{ 128, 0 }, .{ 128, 128 }, .{ 0, 128 } };
    var tags = [_]u8{ 0, 0, 1, 0 };
    var contours = [_]u16{3};
    const g: Glyph = .{ .points = &points, .tags = &tags, .contours = &contours, .advance = 0 };
    const v = try hinted.toVertices(testing.allocator, g);
    defer testing.allocator.free(v);
    try testing.expectEqual(@as(usize, 5), v.len);
    // the start: between the last point and the first
    try testing.expectEqual(@as(i16, 0), v[0].x);
    try testing.expectEqual(@as(i16, 64), v[0].y);
    // a conic through the first control to the implied point, then to the on point
    try testing.expectEqual(@as(i16, 64), v[1].x);
    try testing.expectEqual(@as(i16, 0), v[1].cx);
    try testing.expectEqual(@as(i16, 128), v[2].y);
    // the last control closes the contour back to the start
    try testing.expectEqual(@as(i16, 0), v[3].cx);
    try testing.expectEqual(@as(i16, 128), v[3].cy);
    try testing.expectEqual(@as(i16, 64), v[3].y);
}
