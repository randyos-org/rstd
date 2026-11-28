//! Tests: the charstring evaluator, the glyph path and the outline builder
//! on small hand-made charstrings.

const std = @import("std");
const charstring = @import("../charstring.zig");
const private = @import("../private.zig");
const Blues = @import("../blues.zig").Blues;
const Builder = @import("../Builder.zig");
const GlyphPath = @import("../GlyphPath.zig");
const Index = @import("../../../sfnt/tables/cff/Index.zig");

const testing = std.testing;

/// Draws `code` unhinted at one font unit per 26.6 unit (scale 1/1024).
fn draw(out: *Builder, code: []const u8) !?i32 {
    const priv: private.Private = .{};
    const blues: Blues = .init(&priv, 0x400);
    var path: GlyphPath = .init(&blues, out, 0x400 * 64, 0x400 * 64, false);
    const empty: Index = try .init(&.{ 0, 0 }, true);
    const r = try charstring.evaluate(.{ .global = empty, .local = null }, code, &path, .{ 0, 0 });
    out.closeContour();
    return r.width;
}

/// `n` as a one-byte charstring integer.
fn n(comptime v: i32) u8 {
    return @intCast(v + 139);
}

test "cff path: a square closes without repeating its first point" {
    var out: Builder = .{ .allocator = testing.allocator };
    defer out.deinit();
    // 50 10 20 rmoveto 100 hlineto 100 vlineto -100 hlineto endchar
    const width = try draw(&out, &.{ n(50), n(10), n(20), 21, n(100), 6, n(100), 7, n(-100), 6, 14 });
    try testing.expectEqual(@as(?i32, 50 << 16), width);
    // a unit is 64 here: points are in 26.6
    try testing.expectEqualSlices([2]i32, &.{ .{ 640, 1280 }, .{ 7040, 1280 }, .{ 7040, 7680 }, .{ 640, 7680 } }, out.points.items);
    try testing.expectEqualSlices(u8, &.{ 1, 1, 1, 1 }, out.tags.items);
    try testing.expectEqualSlices(u16, &.{3}, out.contours.items);
}

test "cff path: a curve's control points are cubic, a zero-length line is dropped" {
    var out: Builder = .{ .allocator = testing.allocator };
    defer out.deinit();
    // 0 0 rmoveto 0 0 rlineto 10 0 10 10 0 10 rrcurveto endchar
    const width = try draw(&out, &.{ n(0), n(0), 21, n(0), n(0), 5, n(10), n(0), n(10), n(10), n(0), n(10), 8, 14 });
    try testing.expectEqual(@as(?i32, null), width);
    try testing.expectEqualSlices(u8, &.{ 1, 2, 2, 1 }, out.tags.items);
    try testing.expectEqualSlices([2]i32, &.{ .{ 0, 0 }, .{ 640, 0 }, .{ 1280, 640 }, .{ 1280, 1280 } }, out.points.items);
}

test "cff charstring: arithmetic and subroutine limits" {
    var out: Builder = .{ .allocator = testing.allocator };
    defer out.deinit();
    // 7 2 div 3 add: 6.5, then its width; then a callsubr with no subrs
    try testing.expectEqual(@as(?i32, 0x68000), try draw(&out, &.{ n(7), n(2), 12, 12, n(3), 12, 10, n(0), 22, 14 }));
    try testing.expectError(error.InvalidCharstring, draw(&out, &.{ n(0), 10 }));
}
