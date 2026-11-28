//! Tests: `fixed.zig` and `round.zig` against FreeType's own results, the
//! vectors fontations computed with `FT_MulDiv_No_Round`, `TT_MulFix14`
//! and `FT_Vector_NormLen` (`skrifa/src/outline/glyf/hint/math.rs` and
//! `round.rs`, commit 4f115482, MIT).

const std = @import("std");
const fixed = @import("../../fixed.zig");
const round = @import("../round.zig");

const testing = std.testing;

test "fixed mulDivNoRound: FreeType's results" {
    const cases = [_][4]i32{
        .{ -326, -11474, 9942, 376 },      .{ -6781, 13948, 11973, -7899 },     .{ 3517, 15622, 8075, 6804 },
        .{ -6127, 15026, 2276, -40450 },   .{ 11257, 14828, 2542, 65664 },      .{ -12797, -16280, -9086, -22929 },
        .{ -7994, -3340, 9583, 2786 },     .{ -16101, -13780, -1427, -155481 }, .{ 10304, -16331, 15480, -10870 },
        .{ -15879, 11912, -4650, 40677 },  .{ -5015, 6382, -15977, 2003 },      .{ 2080, -11930, -15457, 1605 },
        .{ -11071, 13350, 16138, -9158 },  .{ 16084, -13564, -770, 283329 },    .{ 14304, -10377, -21, 7068219 },
        .{ -14056, -8853, -5488, -22674 }, .{ -10319, 14797, 8554, -17850 },    .{ -7820, 6826, 10555, -5057 },
        .{ 7257, 15928, 8159, 14167 },     .{ 14929, 11579, -13204, -13091 },   .{ 2808, 12070, -14697, -2306 },
        .{ -13818, 8544, -1649, 71595 },   .{ 3265, 7325, -1373, -17418 },      .{ 14832, 10586, -6440, -24380 },
        .{ 4123, 8274, -2022, -16871 },    .{ 4645, -4149, -7242, 2661 },       .{ -3891, 8366, 5771, -5640 },
        .{ -15447, -3428, -9335, -5672 },  .{ 13670, -14311, -11122, 17589 },   .{ 12590, -6592, 13159, -6306 },
        .{ -8369, -10193, 5051, 16888 },   .{ -9539, 5167, 2595, -18993 },
    };
    for (cases) |c| try testing.expectEqual(c[3], fixed.mulDivNoRound(c[0], c[1], c[2]));
    try testing.expectEqual(std.math.minInt(i32), fixed.mulDivNoRound(std.math.minInt(i32), 1, 1));
    try testing.expectEqual(std.math.minInt(i32), fixed.mulDivNoRound(std.math.minInt(i32), -1, 1));
}

test "fixed mul14: FreeType's results" {
    const cases = [_][3]i32{
        .{ 6236, -10078, -3836 },  .{ -6803, -5405, 2244 },    .{ -10006, -12852, 7849 },  .{ -15434, -4102, 3864 },
        .{ -8681, 9269, -4911 },   .{ 9449, -9130, -5265 },    .{ 12643, 2161, 1668 },     .{ -6115, 9284, -3465 },
        .{ 316, 3390, 65 },        .{ 15077, -12901, -11872 }, .{ -12182, 11613, -8635 },  .{ -7213, 8246, -3630 },
        .{ 13482, 8096, 6662 },    .{ 5690, 15016, 5215 },     .{ -5991, 12613, -4612 },   .{ 13112, -8404, -6726 },
        .{ 13524, 6786, 5601 },    .{ 7156, 3291, 1437 },      .{ -2978, 353, -64 },       .{ -1755, 14626, -1567 },
        .{ 14402, 7886, 6932 },    .{ 7124, 15730, 6840 },     .{ -12679, 14830, -11476 }, .{ -9374, -12999, 7437 },
        .{ 12301, -4685, -3517 },  .{ 5324, 2066, 671 },       .{ 6783, -4946, -2048 },    .{ 12078, -968, -714 },
        .{ -10137, 14116, -8734 }, .{ -13946, 11585, -9861 },  .{ -678, -2205, 91 },       .{ -2629, -3319, 533 },
    };
    for (cases) |c| try testing.expectEqual(c[2], fixed.mul14(c[0], c[1]));
}

test "fixed normalize14: FreeType's unit vectors" {
    const cases = [_][4]i32{
        .{ -13660, 11807, -12395, 10713 }, .{ -10763, 9293, -12401, 10707 },   .{ -3673, 673, -16115, 2952 },
        .{ 15886, -2964, 16106, -3005 },   .{ 15442, -2871, 16108, -2994 },    .{ -6308, 5744, -12114, 11031 },
        .{ 9410, -10415, 10983, -12156 },  .{ -10620, -14856, -9528, -13328 }, .{ -9372, 12029, -10069, 12924 },
        .{ -1272, -1261, -11635, -11534 }, .{ -7076, -5517, -12920, -10074 },  .{ -10297, 179, -16381, 284 },
        .{ 9256, -13235, 9389, -13426 },   .{ 5315, -12449, 6433, -15068 },    .{ 8064, 15213, 7673, 14476 },
        .{ -8665, 41, -16383, 77 },        .{ -3455, -4720, -9677, -13220 },   .{ 13449, -5152, 15299, -5861 },
        .{ -15605, 8230, -14492, 7643 },   .{ 4716, -13690, 5336, -15490 },    .{ 12904, -11422, 12268, -10859 },
        .{ 2825, -6396, 6619, -14987 },    .{ 4654, 15245, 4783, 15670 },      .{ -14769, 15133, -11443, 11725 },
        .{ -8090, -9057, -10914, -12219 }, .{ -472, 1953, -3848, 15925 },      .{ -12563, 1040, -16328, 1351 },
        .{ -7938, 15587, -7435, 14599 },   .{ -9701, 5356, -14343, 7919 },     .{ -642, -14484, -725, -16367 },
        .{ 12963, -9690, 13123, -9809 },   .{ 7067, 5361, 13053, 9902 },       .{ 0x4000, 0, 0x4000, 0 },
        .{ 0, 0x4000, 0, 0x4000 },         .{ -0x4000, 0, -0x4000, 0 },        .{ 0, -0x4000, 0, -0x4000 },
    };
    for (cases) |c| try testing.expectEqual([2]i32{ c[2], c[3] }, fixed.normalize14(c[0], c[1]));
}

test "fixed rounding wraps at the range's end" {
    try testing.expectEqual(std.math.minInt(i32), fixed.round(std.math.maxInt(i32)));
    try testing.expectEqual(std.math.minInt(i32), fixed.ceil(std.math.maxInt(i32)));
    try testing.expectEqual(std.math.minInt(i32), fixed.roundPad(std.math.maxInt(i32), 32));
}

fn expectRound(mode: round.Mode, cases: []const [2]i32) !void {
    for (cases) |c| try testing.expectEqual(c[1], round.round(mode, .{}, c[0]));
}

test "round: each mode, FreeType 2.14.3's half-grid fallback to 32" {
    try expectRound(.grid, &.{ .{ 0, 0 }, .{ 32, 64 }, .{ -32, -64 }, .{ 64, 64 }, .{ 50, 64 } });
    try expectRound(.half_grid, &.{ .{ 0, 32 }, .{ 32, 32 }, .{ -32, -32 }, .{ 64, 96 }, .{ 50, 32 } });
    try expectRound(.double_grid, &.{ .{ 0, 0 }, .{ 32, 32 }, .{ -32, -32 }, .{ 64, 64 }, .{ 50, 64 } });
    try expectRound(.down_to_grid, &.{ .{ 0, 0 }, .{ 32, 0 }, .{ -32, 0 }, .{ 64, 64 }, .{ 50, 0 } });
    try expectRound(.up_to_grid, &.{ .{ 0, 0 }, .{ 32, 64 }, .{ -32, -64 }, .{ 64, 64 }, .{ 50, 64 } });
    try expectRound(.off, &.{ .{ 0, 0 }, .{ 32, 32 }, .{ -32, -32 }, .{ 64, 64 }, .{ 50, 50 } });
    // at the range's end: half grid stays positive, grid wraps and falls back to 0
    try testing.expectEqual(@as(i32, 2147483616), round.round(.half_grid, .{}, std.math.maxInt(i32)));
    try testing.expectEqual(@as(i32, 0), round.round(.grid, .{}, std.math.maxInt(i32)));
}

test "round: super modes at FreeType's extremes do not trap" {
    for ([_]i32{ std.math.maxInt(i32), std.math.minInt(i32), 0, 100, -100 }) |d| {
        _ = round.round(.super, .{ .threshold = std.math.maxInt(i32), .phase = std.math.minInt(i32), .period = std.math.minInt(i32) }, d);
        _ = round.round(.super45, .{ .threshold = std.math.maxInt(i32), .phase = std.math.minInt(i32), .period = -1 }, d);
        _ = round.round(.grid, .{}, d);
    }
    // SROUND's half grid, phase 0, threshold half: 1.5 px
    try testing.expectEqual(@as(i32, 96), round.round(.super, .{ .period = 32, .phase = 0, .threshold = 16 }, 90));
}
