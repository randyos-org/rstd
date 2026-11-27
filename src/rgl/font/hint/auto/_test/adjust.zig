//! Tests: the adjustment database's lookup, as FreeType's.

const std = @import("std");
const adjust = @import("../adjust.zig");

const testing = std.testing;

test "autohint adjustments: FreeType's lookup, misses and all" {
    try testing.expectEqual(adjust.up | adjust.no_height_check, adjust.lookup('!'));
    try testing.expectEqual(adjust.up | adjust.tilde_top, adjust.lookup(0xC3));
    try testing.expectEqual(adjust.up, adjust.lookup('i'));
    try testing.expectEqual(adjust.up | adjust.no_height_check, adjust.lookup('?'));
    // `Q` is out of order in FreeType's table: its binary search misses it
    try testing.expectEqual(@as(u32, 0), adjust.lookup('Q'));
    // `~` is in an `#if 0` block
    try testing.expectEqual(@as(u32, 0), adjust.lookup('~'));
    try testing.expectEqual(@as(u32, 0), adjust.lookup(' '));
}
