//! `histogramInto` over strided views, and its errors.

const std = @import("std");
const strided = @import("../../blas/strided.zig");
const stats = @import("../__root__.zig");

const testing = std.testing;
const ConstStrided = strided.ConstStrided;

test "histogramInto: a strided view, NumPy's bins, the last closed" {
    // every other element, backward: NaN, 3, 2, 1, 0.5, 0
    const data = [_]f64{ 0, 9, 0.5, 9, 1, 9, 2, 9, 3, 9, std.math.nan(f64) };
    const x = try ConstStrided(f64).init(&data, 6, -2);
    var counts: [2]u64 = .{ 1, 0 };
    try stats.histogramInto(f64, &.{ 0, 1, 3 }, x, &counts);
    // 0 and 0.5 in [0, 1), 1, 2 and 3 in [1, 3]; the NaN nowhere; added to 1
    try testing.expectEqualSlices(u64, &.{ 3, 3 }, &counts);
}

test "histogramInto: too few edges, unsorted edges, too few counts" {
    const x: ConstStrided(f32) = .of(&.{ 1, 2 });
    var counts: [2]u64 = @splat(0);
    try testing.expectError(error.ZeroBins, stats.histogramInto(f32, &.{0}, x, &counts));
    try testing.expectError(error.EdgesNotSorted, stats.histogramInto(f32, &.{ 0, 2, 1 }, x, &counts));
    try testing.expectError(error.EdgesNotSorted, stats.histogramInto(f32, &.{ 0, std.math.nan(f32) }, x, &counts));
    try testing.expectError(error.SliceTooShort, stats.histogramInto(f32, &.{ 0, 1, 2, 3 }, x, &counts));
}
