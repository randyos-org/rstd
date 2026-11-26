//! One axis's share of a variation region (OpenType 1.9.1, "OpenType Font
//! Variations Common Table Formats", "Algorithm for interpolation of
//! instance values"), shared by `gvar` tuples and Item Variation Stores.
//! Values are raw F2Dot14 units (`units`).

const types = @import("../types.zig");
const zigstd = @import("../../../../zigstd.zig");

/// 1 at `peak`, falling linearly to 0 at `start` and `end`, 0 outside. A
/// region whose start, peak and end are out of order, or that crosses 0,
/// or peaks at 0, does not depend on this axis: 1.
pub fn tent(v: f32, start: f32, peak: f32, end: f32) f32 {
    if (peak == 0 or start > peak or peak > end or (start < 0 and end > 0)) return 1;
    if (v < start or v > end) return 0;
    if (v == peak) return 1;
    return if (v < peak) (v - start) / (peak - start) else (end - v) / (end - peak);
}

/// An F2Dot14 in raw units, as the region math takes it.
pub fn units(v: types.F2Dot14) f32 {
    return @floatFromInt(@as(i16, @bitCast(v)));
}

test "tent rises to its peak and falls to its end; bad regions ignore the axis" {
    try zigstd.testExpectEqual(@as(f32, 0.5), tent(0.5, 0, 1, 1));
    try zigstd.testExpectEqual(@as(f32, 1), tent(1, 0, 1, 1));
    try zigstd.testExpectEqual(@as(f32, 0), tent(-0.5, 0, 1, 1));
    try zigstd.testExpectEqual(@as(f32, 0.5), tent(0.75, 0.5, 0.5, 1));
    try zigstd.testExpectEqual(@as(f32, 1), tent(0.3, -1, 1, 1));
    try zigstd.testExpectEqual(@as(f32, 1), tent(0.3, 0.5, 0.25, 1));
}
