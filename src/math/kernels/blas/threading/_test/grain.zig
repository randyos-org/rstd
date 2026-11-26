//! `threading.crew`'s tests.

const std = @import("std");
const types = @import("../../types.zig");
const threading = @import("../../threading.zig");

const crew = threading.crew;
const product = threading.product;
const testing = std.testing;

test "crew: one worker per grain, a complex product four, never more than asked" {
    const grain: u64 = 1_000_000;
    try testing.expectEqual(@as(u32, 1), crew(f64, 12, grain, product(&.{ 64, 64, 64 })));
    try testing.expectEqual(@as(u32, 1), crew(f64, 12, grain, product(&.{ 100, 100, 199 })));
    try testing.expectEqual(@as(u32, 2), crew(f64, 12, grain, product(&.{ 100, 100, 200 })));
    try testing.expectEqual(@as(u32, 4), crew(types.Complex64, 12, grain, product(&.{ 100, 100, 100 })));
    try testing.expectEqual(@as(u32, 12), crew(f32, 12, grain, product(&.{ 512, 512, 512 })));
    try testing.expectEqual(@as(u32, 1), crew(f32, 1, grain, product(&.{ 512, 512, 512 })));
    try testing.expectEqual(@as(u32, 3), crew(f32, 3, 0, 8));
    try testing.expectEqual(@as(u32, 5), crew(f64, 5, grain, std.math.maxInt(u64)));
    try testing.expectEqual(@as(u32, 5), crew(types.Complex64, 5, grain, std.math.maxInt(u64)));
}

test "share: one worker per grain of work, at least one, never more than asked" {
    const share = threading.share;
    const grain = threading.matrix_grain;
    try testing.expectEqual(@as(u32, 1), share(12, grain, 64 * 64 * 16));
    try testing.expectEqual(@as(u32, 6), share(12, grain, 512 * 512 * 4));
    try testing.expectEqual(@as(u32, 12), share(12, grain, 4096 * 4096 * 4));
    try testing.expectEqual(@as(u32, 2), share(2, grain, 4096 * 4096 * 4));
    try testing.expectEqual(@as(u32, 3), share(3, 0, 8));
}
