//! Tests: the interpreter on synthetic programs, instruction by instruction,
//! with FreeType 2.14.3's answers for bad input (`ttinterp.c`).

const std = @import("std");
const fixed = @import("../../fixed.zig");
const Rig = @import("rig.zig");
const zone = @import("../zone.zig");

const testing = std.testing;
const expectStack = Rig.expectStack;

test "truetype arithmetic: ADD to MIN in 26.6, DIV by 0" {
    var r: Rig = .{};
    const e = r.init();
    try expectStack(e, &.{ 0xB1, 5, 3, 0x60 }, &.{8});
    try expectStack(e, &.{ 0xB1, 5, 3, 0x61 }, &.{2});
    // 2.0 * 3.0, 6.0 / 2.0
    try expectStack(e, &.{ 0xB9, 0, 128, 0, 192, 0x63 }, &.{384});
    try expectStack(e, &.{ 0xB9, 1, 0x80, 0, 128, 0x62 }, &.{192});
    try expectStack(e, &.{ 0xB8, 0xFF, 0xFB, 0x64 }, &.{5});
    try expectStack(e, &.{ 0xB0, 5, 0x65 }, &.{-5});
    try expectStack(e, &.{ 0xB0, 100, 0x66 }, &.{64});
    try expectStack(e, &.{ 0xB0, 100, 0x67 }, &.{128});
    try expectStack(e, &.{ 0xB1, 3, 9, 0x8B, 0xB1, 3, 9, 0x8C }, &.{ 9, 3 });
    try testing.expectError(error.DivideByZero, Rig.run(e, .font, &.{ 0xB1, 5, 0, 0x62 }));
}

test "truetype logic: comparisons, AND, OR, NOT, ODD and EVEN by the round state" {
    var r: Rig = .{};
    const e = r.init();
    try expectStack(e, &.{ 0xB1, 3, 5, 0x50, 0xB1, 5, 3, 0x50, 0xB1, 5, 5, 0x51, 0xB1, 5, 3, 0x52 }, &.{ 1, 0, 1, 1 });
    try expectStack(e, &.{ 0xB1, 5, 5, 0x53, 0xB1, 4, 5, 0x54, 0xB1, 4, 5, 0x55 }, &.{ 1, 0, 1 });
    try expectStack(e, &.{ 0xB1, 1, 0, 0x5A, 0xB1, 1, 0, 0x5B, 0xB0, 0, 0x5C }, &.{ 0, 1, 1 });
    try expectStack(e, &.{ 0xB0, 64, 0x56, 0xB0, 128, 0x56, 0xB0, 128, 0x57 }, &.{ 1, 0, 1 });
}

test "truetype stack: DUP to ROLL, the pushes, CINDEX's 0 for a bad index" {
    var r: Rig = .{};
    const e = r.init();
    try expectStack(e, &.{ 0xB1, 1, 2, 0x23 }, &.{ 2, 1 });
    try expectStack(e, &.{ 0xB0, 7, 0x20 }, &.{ 7, 7 });
    try expectStack(e, &.{ 0xB2, 1, 2, 3, 0x24 }, &.{ 1, 2, 3, 3 });
    try expectStack(e, &.{ 0xB3, 10, 20, 30, 2, 0x25 }, &.{ 10, 20, 30, 20 });
    try expectStack(e, &.{ 0xB3, 10, 20, 30, 3, 0x26 }, &.{ 20, 30, 10 });
    try expectStack(e, &.{ 0xB2, 1, 2, 3, 0x8A }, &.{ 2, 3, 1 });
    try expectStack(e, &.{ 0xB1, 1, 2, 0x22 }, &.{});
    try expectStack(e, &.{ 0x41, 2, 0xFF, 0xFF, 0x01, 0x00 }, &.{ -1, 256 });
    try expectStack(e, &.{ 0x40, 3, 200, 0, 9, 0x21 }, &.{ 200, 0 });
    try expectStack(e, &.{ 0xB1, 10, 5, 0x25 }, &.{ 10, 0 });
    e.pedantic = true;
    try testing.expectError(error.InvalidReference, Rig.run(e, .font, &.{ 0xB1, 10, 5, 0x25 }));
}

test "truetype stack: too few values zero-fill as FreeType does, unless pedantic; overflow; truncated pushes" {
    var r: Rig = .{};
    const e = r.init();
    // ADD with one value sees two zeros
    try expectStack(e, &.{ 0xB0, 9, 0x60 }, &.{0});
    try testing.expectError(error.CodeOverflow, Rig.run(e, .font, &.{ 0xB1, 5 }));
    var big: [2 + 70]u8 = @splat(1);
    big[0] = 0x40;
    big[1] = 70;
    try testing.expectError(error.StackOverflow, Rig.run(e, .font, &big));
    e.pedantic = true;
    try testing.expectError(error.TooFewArguments, Rig.run(e, .font, &.{ 0xB0, 9, 0x60 }));
}

test "truetype flow: IF, ELSE and nested EIF; JMPR, JROT, JROF; bad and endless jumps" {
    var r: Rig = .{};
    const e = r.init();
    try expectStack(e, &.{ 0xB0, 1, 0x58, 0xB0, 7, 0x1B, 0xB0, 8, 0x59 }, &.{7});
    try expectStack(e, &.{ 0xB0, 0, 0x58, 0xB0, 7, 0x1B, 0xB0, 8, 0x59 }, &.{8});
    try expectStack(e, &.{ 0xB0, 0, 0x58, 0xB0, 1, 0x58, 0xB0, 5, 0x59, 0xB0, 6, 0x1B, 0xB0, 9, 0x59 }, &.{9});
    try expectStack(e, &.{ 0xB0, 3, 0x1C, 0xB0, 1, 0xB0, 2 }, &.{2});
    try expectStack(e, &.{ 0xB1, 3, 1, 0x78, 0xB0, 1, 0xB0, 2 }, &.{2});
    try expectStack(e, &.{ 0xB1, 3, 0, 0x78, 0xB0, 1, 0xB0, 2 }, &.{ 1, 2 });
    try expectStack(e, &.{ 0xB1, 3, 0, 0x79, 0xB0, 1, 0xB0, 2 }, &.{2});
    try testing.expectError(error.BadArgument, Rig.run(e, .font, &.{ 0xB0, 0, 0x1C }));
    try testing.expectError(error.BadArgument, Rig.run(e, .font, &.{ 0xB8, 0xFF, 0x00, 0x1C }));
    // back to the push, forever, until the backward-jump limit
    try testing.expectError(error.ExecutionTooLong, Rig.run(e, .font, &.{ 0xB8, 0xFF, 0xFD, 0x1C }));
    try testing.expectError(error.CodeOverflow, Rig.run(e, .font, &.{ 0xB0, 0, 0x58, 0xB0, 1 }));
}

test "truetype definitions: FDEF and CALL, LOOPCALL, IDEF, and their errors" {
    var r: Rig = .{};
    const e = r.init();
    try expectStack(e, &.{ 0xB0, 0, 0x2C, 0xB0, 42, 0x2D, 0xB0, 0, 0x2B }, &.{42});
    // defined in fpgm, called from prep
    try Rig.run(e, .control_value, &.{ 0xB0, 0, 0x2B, 0xB0, 0, 0x2B });
    try testing.expectEqualSlices(i32, &.{ 42, 42 }, Rig.values(e));
    try expectStack(e, &.{ 0xB0, 1, 0x2C, 0xB0, 1, 0x2D, 0xB1, 3, 1, 0x2A }, &.{ 1, 1, 1 });
    try expectStack(e, &.{ 0xB0, 1, 0x2C, 0xB0, 1, 0x2D, 0xB1, 0, 1, 0x2A }, &.{});
    try expectStack(e, &.{ 0xB0, 0x93, 0x89, 0xB0, 77, 0x2D, 0x93 }, &.{77});
    try testing.expectError(error.UnhandledOpcode, Rig.run(e, .font, &.{0x93}));
    try testing.expectError(error.InvalidDefinition, Rig.run(e, .font, &.{ 0xB0, 5, 0x2B }));
    try testing.expectError(error.DefinitionInGlyphProgram, Rig.run(e, .glyph, &.{ 0xB0, 0, 0x2C, 0x2D }));
    try testing.expectError(error.NestedDefinition, Rig.run(e, .font, &.{ 0xB0, 0, 0x2C, 0xB0, 1, 0x2C, 0x2D, 0x2D }));
    try testing.expectError(error.EndfWithoutCall, Rig.run(e, .font, &.{0x2D}));
    try testing.expectError(error.CodeOverflow, Rig.run(e, .font, &.{ 0xB0, 0, 0x2C, 0xB0, 1 }));
    // nine functions where maxp allows eight
    var many: [9 * 4]u8 = undefined;
    for (0..9) |i| many[4 * i ..][0..4].* = .{ 0xB0, @intCast(i), 0x2C, 0x2D };
    try testing.expectError(error.TooManyDefinitions, Rig.run(e, .font, &many));
}

test "truetype storage and CVT: reads, writes, a glyph's copy, WCVTF in place, DELTAC" {
    var r: Rig = .{};
    const e = r.init();
    try expectStack(e, &.{ 0xB1, 3, 99, 0x42, 0xB0, 3, 0x43 }, &.{99});
    try testing.expectEqual(@as(i32, 99), r.store[3]);
    // a glyph's writes go to its copy
    try Rig.run(e, .glyph, &.{ 0xB1, 2, 55, 0x42, 0xB1, 4, 7, 0x44 });
    try testing.expectEqual(@as(i32, 0), r.store[2]);
    try testing.expectEqual(@as(i32, 55), r.store_copy[2]);
    try testing.expectEqual(@as(i32, 256), r.cvt[4]);
    try testing.expectEqual(@as(i32, 7), r.cvt_copy[4]);
    // WCVTF before any other write reaches the size's CVT
    try Rig.run(e, .glyph, &.{ 0xB1, 1, 100, 0x70 });
    try testing.expectEqual(@as(i32, 50), r.cvt[1]);
    try expectStack(e, &.{ 0xB0, 2, 0x45, 0xB0, 20, 0x45, 0xB0, 20, 0x43 }, &.{ 128, 0, 0 });
    // at 12 ppem, delta base 9: one step of 1/8 px on CVT 2
    try Rig.run(e, .control_value, &.{ 0xB2, 0x38, 2, 1, 0x73 });
    try testing.expectEqual(@as(i32, 136), r.cvt[2]);
    // another ppem's pair is consumed but not applied
    try Rig.run(e, .control_value, &.{ 0xB2, 0x48, 2, 1, 0x73 });
    try testing.expectEqual(@as(i32, 136), r.cvt[2]);
    try testing.expectEqual(@as(usize, 0), e.top);
    e.pedantic = true;
    try testing.expectError(error.InvalidReference, Rig.run(e, .font, &.{ 0xB0, 20, 0x45 }));
    try testing.expectError(error.InvalidReference, Rig.run(e, .font, &.{ 0xB1, 20, 1, 0x42 }));
}

test "truetype graphics state: vectors, super rounding, INSTCTRL, SCANCTRL and bad settings" {
    var r: Rig = .{};
    const e = r.init();
    try expectStack(e, &.{ 0x00, 0x0C, 0x0D }, &.{ 0, 0x4000, 0, 0x4000 });
    const n = fixed.normalize14(1, 1);
    try expectStack(e, &.{ 0xB1, 1, 1, 0x0A, 0x0C }, &.{ n[0], n[1] });
    // a zero vector leaves the old one
    try expectStack(e, &.{ 0x01, 0xB1, 0, 0, 0x0A, 0x0C }, &.{ 0x4000, 0 });
    try Rig.run(e, .font, &.{ 0xB0, 0x48, 0x76 });
    try testing.expectEqual(@as(i32, 64), e.super_round.period);
    try testing.expectEqual(@as(i32, 0), e.super_round.phase);
    try testing.expectEqual(@as(i32, 32), e.super_round.threshold);
    try expectStack(e, &.{ 0xB0, 0x48, 0x76, 0xB0, 90, 0x68 }, &.{64});
    try Rig.run(e, .control_value, &.{ 0xB1, 4, 3, 0x8E });
    try testing.expectEqual(@as(u8, 4), e.gs.instruct_control);
    e.backward_compatibility = 4;
    try Rig.run(e, .glyph, &.{ 0xB1, 4, 3, 0x8E });
    try testing.expectEqual(@as(u8, 0), e.backward_compatibility);
    try Rig.run(e, .font, &.{ 0xB8, 0x01, 0x0C, 0x85 });
    try testing.expect(e.gs.scan_control);
    try testing.expectError(error.BadArgument, Rig.run(e, .font, &.{ 0xB8, 0xFF, 0xFF, 0x17 }));
    try testing.expectError(error.BadArgument, Rig.run(e, .font, &.{ 0xB0, 7, 0x5F }));
    try testing.expectError(error.DebugOpcode, Rig.run(e, .font, &.{ 0xB0, 0, 0x4F }));
    try Rig.run(e, .font, &.{ 0xB0, 2, 0x13 });
    try testing.expectEqual(zone.Pointer.glyph, e.zp[0]);
    e.pedantic = true;
    try testing.expectError(error.InvalidReference, Rig.run(e, .font, &.{ 0xB0, 2, 0x13 }));
    try testing.expectError(error.InvalidReference, Rig.run(e, .font, &.{ 0xB1, 4, 3, 0x8E }));
}

test "truetype queries: GETINFO by version and target, MPPEM, MPS, GETVARIATION, GETDATA" {
    var r: Rig = .{};
    const e = r.init();
    const v40_bits: i32 = 40 | (1 << 13) | (1 << 17) | (1 << 18) | (1 << 19);
    try expectStack(e, &.{ 0xB8, 0x1C, 0x41, 0x88 }, &.{v40_bits});
    e.target = .mono;
    try expectStack(e, &.{ 0xB8, 0x1C, 0x41, 0x88 }, &.{40});
    e.target = .normal;
    e.version = .v35;
    try expectStack(e, &.{ 0xB0, 33, 0x88, 0x4B, 0x4C }, &.{ 35 | (1 << 12), 12, 12 });
    e.version = .v40;
    try expectStack(e, &.{0x4C}, &.{768});
    try testing.expectError(error.UnhandledOpcode, Rig.run(e, .font, &.{0x91}));
    e.axis_count = 2;
    e.coords = &.{0x10000};
    try expectStack(e, &.{ 0x91, 0x92 }, &.{ 0x4000, 0, 17 });
}

test "truetype measures: GC and MD on the glyph and twilight zones" {
    var r: Rig = .{};
    const e = r.init();
    // x: point 2 at 60; MD from point 1 to 3, current and from font units
    try expectStack(e, &.{ 0x01, 0xB0, 2, 0x46, 0xB0, 2, 0x47, 0xB1, 3, 1, 0x4A, 0xB1, 3, 1, 0x49 }, &.{ 60, 60, 10, 10 });
    // y; with zp0 on the twilight zone MD measures original points: (0, 0) to point 1
    try expectStack(e, &.{ 0x00, 0xB1, 3, 1, 0x4A, 0xB0, 0, 0x13, 0xB1, 3, 1, 0x49 }, &.{ -50, 25 });
    try expectStack(e, &.{ 0xB0, 50, 0x46 }, &.{0});
}

test "truetype tags: FLIPPT by the loop, FLIPRGON, FLIPRGOFF, UTP" {
    var r: Rig = .{};
    const e = r.init();
    try Rig.run(e, .glyph, &.{ 0xB2, 1, 2, 2, 0x17, 0x80 });
    try testing.expectEqualSlices(u8, &.{ 1, 0, 0, 1, 1, 1, 1, 1 }, &r.tags);
    try testing.expectEqual(@as(i32, 1), e.gs.loop);
    try Rig.run(e, .glyph, &.{ 0xB1, 4, 6, 0x82, 0xB1, 1, 1, 0x81 });
    try testing.expectEqualSlices(u8, &.{ 1, 1, 0, 1, 0, 0, 0, 1 }, &r.tags);
    r.tags[5] |= zone.touched_x | zone.touched_y;
    try Rig.run(e, .glyph, &.{ 0x01, 0xB0, 5, 0x29 });
    try testing.expectEqual(zone.touched_y, r.tags[5]);
}
