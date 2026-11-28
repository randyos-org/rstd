//! Tests: walking a charmap's mapped codes in order (`nextChar`).

const std = @import("std");
const bytes = @import("../../../../../../bytes.zig");
const CmapSubtable = @import("../__root__.zig").CmapSubtable;

test "cmap nextChar: codes in order, unmapped ones skipped, across groups" {
    // format 12: 0x41-0x43 to glyphs 5-7, 0x1F600-0x1F601 to glyphs 0 and 1
    var buf: [16 + 2 * 12]u8 = @splat(0);
    bytes.setU16BE(&buf, 12);
    bytes.setU32BE(buf[4..], buf.len);
    bytes.setU32BE(buf[12..], 2);
    for ([_][3]u32{ .{ 0x41, 0x43, 5 }, .{ 0x1F600, 0x1F601, 0 } }, 0..) |g, i| {
        for (g, 0..) |v, k| bytes.setU32BE(buf[16 + i * 12 + k * 4 ..], v);
    }
    const sub = try CmapSubtable.init(&buf, true);
    var cp: u32 = 0;
    var seen: [4]u32 = undefined;
    var n: usize = 0;
    while (sub.nextChar(cp)) |m| : (cp = m.codepoint) {
        seen[n] = m.codepoint;
        n += 1;
    }
    // 0x1F600 maps to glyph 0
    try std.testing.expectEqualSlices(u32, &.{ 0x41, 0x42, 0x43, 0x1F601 }, seen[0..n]);
    try std.testing.expectEqual(@as(u16, 7), @backingInt(sub.nextChar(0x42).?.glyph));
}
