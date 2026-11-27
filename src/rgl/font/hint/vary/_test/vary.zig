//! Tests: packed point numbers and deltas, tuple scalars, inferred deltas,
//! and (under `-Dtest-fonts`) Cascadia Code's normalized coordinates.

const std = @import("std");
const packed_data = @import("../packed.zig");
const iup = @import("../iup.zig");
const Blend = @import("../Blend.zig");
const Sfnt = @import("../../../sfnt/Sfnt.zig");
const test_fonts = @import("rstd_test_fonts");

const testing = std.testing;
const Cursor = packed_data.Cursor;

test "vary packed points: all, byte and word runs, and a short run fails" {
    const al = testing.allocator;
    var c: Cursor = .{ .data = &.{0} };
    try testing.expect((try packed_data.readPoints(al, &c)).? == .all);
    // three byte increments, then two word ones
    c = .{ .data = &.{ 5, 0x02, 1, 2, 3, 0x81, 0x01, 0x00, 0x00, 0x05 } };
    const p = (try packed_data.readPoints(al, &c)).?;
    defer p.deinit(al);
    try testing.expectEqualSlices(u16, &.{ 1, 3, 6, 262, 267 }, p.list);
    try testing.expectEqual(@as(usize, 10), c.at);
    c = .{ .data = &.{ 4, 0x03, 1, 2 } };
    try testing.expect(try packed_data.readPoints(al, &c) == null);
}

test "vary packed deltas: zero, word and byte runs in 16.16, and a short run fails" {
    const al = testing.allocator;
    var c: Cursor = .{ .data = &.{ 0x81, 0x41, 0xFF, 0xFE, 0x00, 0x03, 0x01, 0xFF, 0x05 } };
    const d = (try packed_data.readDeltas(al, &c, 6)).?;
    defer al.free(d);
    try testing.expectEqualSlices(i32, &.{ 0, 0, -2 << 16, 3 << 16, -1 << 16, 5 << 16 }, d);
    c = .{ .data = &.{ 0x41, 0x00 } };
    try testing.expect(try packed_data.readDeltas(al, &c, 2) == null);
}

test "vary tuple scalars: peaks, axes a tuple ignores, and intermediate regions" {
    const half = [_]i32{0x8000};
    try testing.expectEqual(@as(i32, 0x8000), Blend.tupleScalar(&half, &.{0x10000}, null));
    try testing.expectEqual(@as(i32, 0x10000), Blend.tupleScalar(&half, &.{0}, null));
    try testing.expectEqual(@as(i32, 0x10000), Blend.tupleScalar(&half, &.{0x8000}, null));
    try testing.expectEqual(@as(i32, 0), Blend.tupleScalar(&half, &.{-0x10000}, null));
    try testing.expectEqual(@as(i32, 0), Blend.tupleScalar(&.{0}, &.{0x10000}, null));
    // past the peak, falling to the end
    try testing.expectEqual(@as(i32, 0x8000), Blend.tupleScalar(&.{0xC000}, &.{0x8000}, .{ &.{0}, &.{0x10000} }));
    try testing.expectEqual(@as(i32, 0), Blend.tupleScalar(&.{0x10000}, &.{0x8000}, .{ &.{0}, &.{0x10000} }));
}

test "vary inferred deltas: a lone delta shifts its contour, two interpolate" {
    const org = [_]iup.Point{ .{ 0, 0 }, .{ 100, 0 }, .{ 100, 100 }, .{ 0, 100 }, .{ 50, 50 } };
    var out = org;
    out[1][0] += 10;
    var has = [_]bool{ false, true, false, false, true };
    // contour 0..3 with one delta; point 4 its own contour
    iup.interpolate(&.{ 3, 4 }, &out, &org, &has);
    for (out[0..4], org[0..4]) |o, r| try testing.expectEqual(r[0] + 10, o[0]);
    // deltas at x = 0 (+0) and x = 100 (+20): halfway scales by 1.2
    const org2 = [_]iup.Point{ .{ 0, 0 }, .{ 50, 0 }, .{ 100, 0 } };
    var out2 = org2;
    out2[2][0] += 20;
    var has2 = [_]bool{ true, false, true };
    iup.interpolate(&.{2}, &out2, &org2, &has2);
    try testing.expectEqual(@as(i32, 60), out2[1][0]);
}

test "vary instances: Cascadia Code's weights normalized through its avar (-Dtest-fonts)" {
    const bytes = test_fonts.cascadia_code orelse return error.SkipZigTest;
    const al = testing.allocator;
    var arena: std.heap.ArenaAllocator = .init(al);
    defer arena.deinit();
    const sfnt = try Sfnt.init(arena.allocator(), bytes, .{});
    // the default weight is the static font
    try testing.expect(try Blend.init(al, &sfnt, &.{400 << 16}) == null);
    try testing.expect(try Blend.init(al, &sfnt, &.{}) == null);
    // -0.5 maps to avar's -0.666687; 0.5 lands between 0 and 0.666687's
    // 0.599976; past the maximum clamps to 1
    for ([_][2]i32{ .{ 300, -10923 * 4 }, .{ 550, 29489 }, .{ 700, 0x10000 }, .{ 900, 0x10000 } }) |c| {
        var b = (try Blend.init(al, &sfnt, &.{c[0] << 16})).?;
        defer b.deinit();
        try testing.expectEqualSlices(i32, &.{c[1]}, b.coords);
        try testing.expect(b.gvar != null and b.hvar != null);
    }
}
