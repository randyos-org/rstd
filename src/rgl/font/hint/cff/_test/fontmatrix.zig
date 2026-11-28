//! Tests: FontMatrix parsing and normalization (`FontMatrix.zig`), with
//! FreeType's rounding worked by hand, on a DICT and (under `-Dtest-fonts`)
//! on `tools/cfffonts.zig`'s fonts.

const std = @import("std");
const FontMatrix = @import("../FontMatrix.zig");
const Sfnt = @import("../../../sfnt/Sfnt.zig");
const test_fonts = @import("rstd_test_fonts");

const testing = std.testing;

/// `[0.001 0 0.00025 0.001 0.01 0] FontMatrix` as DICT reals.
const oblique = [_]u8{ 30, 0x0A, 0x00, 0x1F, 139, 30, 0x0A, 0x00, 0x02, 0x5F, 30, 0x0A, 0x00, 0x1F, 30, 0x0A, 0x01, 0xFF, 139, 12, 7 };

test "cff FontMatrix: the six reals scaled by the largest one's power of ten" {
    const m = try FontMatrix.parse(&oblique);
    // 0.001 and 0.00025 at 1/100: (65536 + 5) / 10 and (25 * 65536 + 500) / 1000
    try testing.expect(m.has);
    try testing.expectEqual([4]i32{ 6554, 0, 1638, 6554 }, m.matrix);
    try testing.expectEqual([2]i32{ 65536, 0 }, m.offset);
    try testing.expectEqual(@as(i32, 100), m.upm);
    // no FontMatrix: the default
    const d = try FontMatrix.parse(&.{ 139, 15 });
    try testing.expect(!d.has);
    try testing.expectEqual(FontMatrix.identity, d.matrix);
}

test "cff FontMatrix: normalized to yy = 1, a CID subfont's on the top font's (-Dtest-fonts)" {
    var arena: std.heap.ArenaAllocator = .init(testing.allocator);
    defer arena.deinit();
    const seac = test_fonts.cff_seac orelse return error.SkipZigTest;
    const s = try Sfnt.init(arena.allocator(), seac, .{});
    // 1638 / 6554 of 65536 rounds to 16379, and the offset 10 truncates to 9
    const m = try FontMatrix.forGlyph(s.cff.?, 1000, @fromBackingInt(1));
    try testing.expectEqual([4]i32{ 65536, 0, 16379, 65536 }, m.matrix);
    try testing.expectEqual([2]i32{ 9, 0 }, m.offset);
    try testing.expectEqual(@as(i32, 1000), m.top_upm);
    try testing.expectEqual(@as(i32, 1000), m.sub_upm);
    const cid = test_fonts.cff_cid orelse return error.SkipZigTest;
    const c = try Sfnt.init(arena.allocator(), cid, .{});
    // FD 1's half scale: 2000 units per em, its 0.1 shear 0.2, its offset lost
    const fd1 = try FontMatrix.forGlyph(c.cff.?, 1000, @fromBackingInt(2));
    try testing.expectEqual([4]i32{ 65536, 0, 13108, 65536 }, fd1.matrix);
    try testing.expectEqual([2]i32{ 0, 0 }, fd1.offset);
    try testing.expectEqual(@as(i32, 2000), fd1.sub_upm);
    // FD 0's identity, concatenated: the top font's units
    const fd0 = try FontMatrix.forGlyph(c.cff.?, 1000, @fromBackingInt(1));
    try testing.expectEqual(FontMatrix.identity, fd0.matrix);
    try testing.expectEqual(@as(i32, 1000), fd0.sub_upm);
}
