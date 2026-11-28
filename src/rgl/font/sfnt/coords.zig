//! A variable font's normalized axis coordinates from user axis values
//! (OpenType 1.9.1, "OpenType Font Variations Overview", "Coordinate
//! scales and normalization"): clamp to the axis's range, normalize
//! against its default (`fvar`), remap through `avar`, round to F2Dot14.

const std = @import("std");
const tables = @import("tables/__root__.zig");
const types = @import("types.zig");
const VariationAxisRecord = @import("tables/fvar.zig").VariationAxisRecord;
const SegmentMap = @import("tables/avar.zig").SegmentMap;
const zigstd = @import("../../../zigstd.zig");
const bytes = @import("../../../bytes.zig");
const fvar_example = @import("tables/_test/fvar_example.zig");

const F2Dot14 = types.F2Dot14;
const FvarTable = tables.FvarTable;
const AvarTable = tables.AvarTable;

/// A user-scale value for the axis tagged `tag` ("wght", "wdth", ...).
pub const AxisValue = struct {
    tag: []const u8,
    value: f64,
};

pub const Error = error{
    /// `values` names a tag the font has no axis for.
    UnknownAxis,
    /// `out` is shorter than the font's axis count.
    SliceTooShort,
};

/// Fills `out[0..fvar.axisCount()]` with each axis's normalized
/// coordinate: `values`' entry for it (the last, if several), else its
/// default, 0.
pub fn normalize(fvar: FvarTable, avar: ?AvarTable, values: []const AxisValue, out: []F2Dot14) Error!void {
    const count = fvar.axisCount();
    if (out.len < count) return error.SliceTooShort;
    // with a different axis count, avar does not describe these axes
    const maps = if (avar) |a| (if (a.axisCount() == count) a else null) else null;
    for (values) |v| {
        for (0..count) |i| {
            if (std.mem.eql(u8, fvar.getAxis(i).axis_tag.strip(), v.tag)) break;
        } else return error.UnknownAxis;
    }
    for (0..count, out[0..count]) |i, *c| {
        const axis = fvar.getAxis(i);
        var user = axis.default_value.toFloat(f64);
        for (values) |v| {
            if (std.mem.eql(u8, axis.axis_tag.strip(), v.tag)) user = v.value;
        }
        var n = defaultNormalized(axis, user);
        if (maps) |a| n = remap(a.getSegmentMap(i), n);
        c.* = F2Dot14.fromFloat(@round(n * 16384) / 16384);
    }
}

/// The default normalization: -1 at min, 0 at default, 1 at max, linear between.
fn defaultNormalized(axis: VariationAxisRecord, user: f64) f64 {
    const min = axis.min_value.toFloat(f64);
    const default = axis.default_value.toFloat(f64);
    const max = axis.max_value.toFloat(f64);
    const v = std.math.clamp(user, @min(min, default), @max(max, default));
    if (v < default) return -(default - v) / (default - min);
    if (v > default) return (v - default) / (max - default);
    return 0;
}

/// `n` through one axis's segment map, piecewise linear between its
/// pairs. A map without the required -1, 0 and 1 pairs leaves `n` as is.
fn remap(map: SegmentMap, n: f64) f64 {
    const count = map.positionMapCount();
    if (count < 3) return n;
    var prev = map.getAxisValueMap(0);
    if (n <= prev.from_coordinate.toFloat(f64)) return prev.to_coordinate.toFloat(f64);
    for (1..count) |k| {
        const cur = map.getAxisValueMap(k);
        const from0 = prev.from_coordinate.toFloat(f64);
        const from1 = cur.from_coordinate.toFloat(f64);
        if (n < from1) {
            const to0 = prev.to_coordinate.toFloat(f64);
            const to1 = cur.to_coordinate.toFloat(f64);
            return to0 + (to1 - to0) * (n - from0) / (from1 - from0);
        }
        if (n == from1) return cur.to_coordinate.toFloat(f64);
        prev = cur;
    }
    return prev.to_coordinate.toFloat(f64);
}

fn expectCoords(want: []const f64, got: []const F2Dot14) !void {
    for (want, got) |w, g| try zigstd.testExpectApproxEqAbs(w, g.toFloat(f64), 1.0 / 16384.0);
}

test "normalize maps fvar's example axes: min to -1, default to 0, max to 1, clamped" {
    const gpa = zigstd.testAllocator;
    const buf = try fvar_example.build(gpa);
    defer gpa.free(buf);
    const fvar = try FvarTable.init(buf, true);
    var out: [2]F2Dot14 = undefined;
    // wght 300..400..700, wdth 62.5..100..150
    try normalize(fvar, null, &.{}, &out);
    try expectCoords(&.{ 0, 0 }, &out);
    try normalize(fvar, null, &.{ .{ .tag = "wght", .value = 550 }, .{ .tag = "wdth", .value = 81.25 } }, &out);
    try expectCoords(&.{ 0.5, -0.5 }, &out);
    try normalize(fvar, null, &.{ .{ .tag = "wght", .value = 100 }, .{ .tag = "wdth", .value = 999 } }, &out);
    try expectCoords(&.{ -1, 1 }, &out);
    try zigstd.testExpectError(error.UnknownAxis, normalize(fvar, null, &.{.{ .tag = "slnt", .value = 0 }}, &out));
    try zigstd.testExpectError(error.SliceTooShort, normalize(fvar, null, &.{}, out[0..1]));
}

test "normalize remaps through avar between its pairs" {
    const gpa = zigstd.testAllocator;
    const buf = try fvar_example.build(gpa);
    defer gpa.free(buf);
    const fvar = try FvarTable.init(buf, true);
    // axis 0: -1 to -1, 0 to 0, 0.5 to 0.25, 1 to 1; axis 1: no map
    var avar_buf: [8 + 2 + 4 * 4 + 2]u8 = undefined;
    bytes.setU16BE(avar_buf[0..], 1);
    bytes.setU16BE(avar_buf[2..], 0);
    bytes.setU16BE(avar_buf[4..], 0);
    bytes.setU16BE(avar_buf[6..], 2);
    bytes.setU16BE(avar_buf[8..], 4);
    const pairs = [_][2]i16{ .{ -16384, -16384 }, .{ 0, 0 }, .{ 8192, 4096 }, .{ 16384, 16384 } };
    for (pairs, 0..) |p, k| {
        bytes.setI16BE(avar_buf[10 + k * 4 ..], p[0]);
        bytes.setI16BE(avar_buf[12 + k * 4 ..], p[1]);
    }
    bytes.setU16BE(avar_buf[26..], 0);
    const avar = try AvarTable.init(&avar_buf, true);
    var out: [2]F2Dot14 = undefined;
    // wght 550 is 0.5, mapped to 0.25; 625 is 0.75, halfway from 0.25 to 1
    try normalize(fvar, avar, &.{ .{ .tag = "wght", .value = 550 }, .{ .tag = "wdth", .value = 125 } }, &out);
    try expectCoords(&.{ 0.25, 0.5 }, &out);
    try normalize(fvar, avar, &.{.{ .tag = "wght", .value = 625 }}, &out);
    try expectCoords(&.{ 0.625, 0 }, &out);
}
