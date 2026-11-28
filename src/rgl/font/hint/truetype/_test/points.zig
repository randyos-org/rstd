//! Tests: the point-moving instructions on the rig's glyph zone, point `i`
//! at x = 50 + 5 i (26.6) from 100 + 10 i font units, scale 0.5.

const std = @import("std");
const Rig = @import("rig.zig");
const zone = @import("../zone.zig");

const testing = std.testing;

fn xs(r: *const Rig) [8]i32 {
    var out: [8]i32 = undefined;
    for (r.cur, 0..) |p, i| out[i] = p[0];
    return out;
}

test "truetype points: MDAP rounds, MIAP to the CVT with its cut-in, the reference points" {
    var r: Rig = .{};
    const e = r.init();
    try Rig.run(e, .glyph, &.{ 0x01, 0xB0, 1, 0x2F });
    try testing.expectEqual(@as(i32, 64), r.cur[1][0]);
    try testing.expect(r.tags[1] & zone.touched_x != 0);
    try testing.expectEqual(@as(u16, 1), e.gs.rp0);
    try testing.expectEqual(@as(u16, 1), e.gs.rp1);
    try Rig.run(e, .glyph, &.{ 0x01, 0xB1, 3, 2, 0x3E });
    try testing.expectEqual(@as(i32, 128), r.cur[3][0]);
    // |448 - 70| is past the cut-in: the outline's distance, rounded
    try Rig.run(e, .glyph, &.{ 0x01, 0xB1, 4, 7, 0x3F });
    try testing.expectEqual(@as(i32, 64), r.cur[4][0]);
}

test "truetype points: MDRP and MIRP by flags; MSIRP; SHP by a reference's move" {
    var r: Rig = .{};
    const e = r.init();
    // 10 units of 0.5 rounded to 0, the minimum distance making it 64
    try Rig.run(e, .glyph, &.{ 0x01, 0xB0, 0, 0x10, 0xB0, 2, 0xCC });
    try testing.expectEqual(@as(i32, 114), r.cur[2][0]);
    try testing.expectEqual(@as(u16, 0), e.gs.rp1);
    try testing.expectEqual(@as(u16, 2), e.gs.rp2);
    try testing.expectEqual(@as(u16, 0), e.gs.rp0);
    // MIRP unrounded takes CVT 3, 192; rp0 moves to the point
    try Rig.run(e, .glyph, &.{ 0x01, 0xB0, 0, 0x10, 0xB1, 4, 3, 0xF0 });
    try testing.expectEqual(@as(i32, 242), r.cur[4][0]);
    try testing.expectEqual(@as(u16, 4), e.gs.rp0);
    // rounded, the cut-in falls back to the outline's 20, the minimum to 64
    try Rig.run(e, .glyph, &.{ 0x01, 0xB0, 0, 0x10, 0xB1, 5, 3, 0xEC });
    try testing.expectEqual(@as(i32, 114), r.cur[5][0]);
    try Rig.run(e, .glyph, &.{ 0x01, 0xB0, 0, 0x10, 0xB1, 6, 64, 0x3A });
    try testing.expectEqual(@as(i32, 114), r.cur[6][0]);
    // MDAP moves point 1 by 9; SHP[1] shifts point 7 by rp1's move
    try Rig.run(e, .glyph, &.{ 0x01, 0xB0, 1, 0x2F, 0xB0, 7, 0x33 });
    try testing.expectEqual(@as(i32, 94), r.cur[7][0]);
}

test "truetype points: IP keeps the original fraction; IUP interpolates and shifts" {
    var r: Rig = .{};
    const e = r.init();
    r.cur[4][0] += 40;
    try Rig.run(e, .glyph, &.{ 0x01, 0xB1, 4, 0, 0x11, 0x12, 0xB0, 2, 0x39 });
    // half of the way in font units (20 of 40), so half of 60
    try testing.expectEqual(@as(i32, 80), r.cur[2][0]);
    var s: Rig = .{};
    const f = s.init();
    s.cur[0][0] += 10;
    s.cur[4][0] += 20;
    s.tags[0] |= zone.touched_x;
    s.tags[4] |= zone.touched_x;
    s.tags[2] = zone.on_curve;
    try Rig.run(f, .glyph, &.{0x31});
    try testing.expectEqual([8]i32{ 60, 68, 75, 83, 90, 95, 100, 105 }, xs(&s));
    // one touched point shifts the rest of its contour
    var t: Rig = .{};
    const g = t.init();
    t.cur[3][1] += 64;
    t.tags[3] |= zone.touched_y;
    try Rig.run(g, .glyph, &.{0x30});
    for (t.cur, t.org) |c, o| try testing.expectEqual(o[1] + 64, c[1]);
}

test "truetype points: SHPIX, SCFS, ALIGNRP, ALIGNPTS, ISECT, DELTAP" {
    var r: Rig = .{};
    const e = r.init();
    try Rig.run(e, .glyph, &.{ 0x01, 0xB1, 1, 64, 0x38 });
    try testing.expectEqual(@as(i32, 119), r.cur[1][0]);
    try Rig.run(e, .glyph, &.{ 0x01, 0xB1, 2, 200, 0x48 });
    try testing.expectEqual(@as(i32, 200), r.cur[2][0]);
    try Rig.run(e, .glyph, &.{ 0x01, 0xB0, 0, 0x10, 0xB0, 3, 0x3C });
    try testing.expectEqual(@as(i32, 50), r.cur[3][0]);
    try Rig.run(e, .glyph, &.{ 0x01, 0xB1, 4, 6, 0x27 });
    try testing.expectEqual(@as(i32, 75), r.cur[4][0]);
    try testing.expectEqual(@as(i32, 75), r.cur[6][0]);
    // a horizontal line through (0, 0) and (256, 0), a vertical one at x = 128
    r.cur[0] = .{ 0, 0 };
    r.cur[1] = .{ 256, 0 };
    r.cur[2] = .{ 128, -128 };
    r.cur[3] = .{ 128, 128 };
    try Rig.run(e, .glyph, &.{ 0xB4, 5, 0, 1, 2, 3, 0x0F });
    try testing.expectEqual(zone.Point{ 128, 0 }, r.cur[5]);
    try testing.expectEqual(zone.touched_x | zone.touched_y, r.tags[5] & (zone.touched_x | zone.touched_y));
    // at 12 ppem, delta base 9: point 7 one step (1/8 px) right
    try Rig.run(e, .glyph, &.{ 0x01, 0xB2, 0x38, 7, 1, 0x5D });
    try testing.expectEqual(@as(i32, 93), r.cur[7][0]);
}

test "truetype points: v40 backward compatibility freezes x, then y after IUP both ways" {
    var r: Rig = .{};
    const e = r.init();
    e.backward_compatibility = 4;
    try Rig.run(e, .glyph, &.{ 0x01, 0xB1, 1, 64, 0x38, 0x00, 0xB1, 1, 64, 0x38 });
    // SHPIX moves no point untouched in y, in x or y
    try testing.expectEqual(@as(i32, 55), r.cur[1][0]);
    try testing.expectEqual(@as(i32, -25), r.cur[1][1]);
    // MDAP moves point 1 in y; IUP in y shifts point 2 with it; once IUP
    // ran both ways, SCFS cannot move point 2
    try Rig.run(e, .glyph, &.{ 0x00, 0xB0, 1, 0x2F, 0x30, 0x31, 0xB1, 2, 0x80, 0x48 });
    try testing.expectEqual(@as(i32, 0), r.cur[1][1]);
    try testing.expectEqual(@as(i32, -25), r.cur[2][1]);
    try testing.expectEqual(@as(u8, 7), e.backward_compatibility);
}
