//! Tests: FreeType's table checksum for its tricky-font list.

const std = @import("std");
const tricky = @import("../tricky.zig");

test "tricky fonts: the sfnt checksum sums big-endian words, a tail left-aligned" {
    try std.testing.expectEqual(@as(u32, 0x01020304 +% 0x05060708), tricky.checksum(&.{ 1, 2, 3, 4, 5, 6, 7, 8 }));
    try std.testing.expectEqual(@as(u32, 0x01020304 +% 0x05060000), tricky.checksum(&.{ 1, 2, 3, 4, 5, 6 }));
    try std.testing.expectEqual(@as(u32, 0), tricky.checksum(&.{}));
}
