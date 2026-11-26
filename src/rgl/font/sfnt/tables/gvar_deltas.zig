//! A glyph's `gvar` deltas at normalized coordinates (OpenType 1.9.1, "gvar"
//! and "OpenType Font Variations Common Table Formats", "Tuple Variation
//! Store"): each tuple's scalar, its packed point numbers and deltas, the
//! points it leaves out inferred (`iup.zig`), summed per point.

const std = @import("std");
const bytes = @import("../../../../bytes.zig");
const types = @import("../types.zig");
const zigstd = @import("../../../../zigstd.zig");
const GvarTable = @import("gvar.zig").GvarTable;
const tv = @import("tuple_variation.zig");
const iup = @import("iup.zig");
const gvar_example = @import("_test/gvar_example.zig");

const Allocator = zigstd.Allocator;
const F2Dot14 = types.F2Dot14;

pub const Error = Allocator.Error || error{
    /// `coords.len` is not the font's axis count.
    WrongAxisCount,
    GlyphOutOfRange,
    /// A header, point list or delta run runs past its data.
    MalformedVariationData,
};

/// A glyph's original points, for inferring the deltas of the points a
/// tuple leaves out. A composite has none: its components get no inferred deltas.
pub const Outline = struct {
    xs: []const i16 = &.{},
    ys: []const i16 = &.{},
    /// each contour's last point index
    ends: []const u16 = &.{},
};

/// Writes `glyph`'s deltas at `coords` to `dx` and `dy`, one per point: a
/// simple glyph's points or a composite's components, then the four
/// phantom points. Zero where the glyph has no variation data.
pub fn glyphDeltas(gpa: Allocator, gvar: GvarTable, glyph: usize, coords: []const F2Dot14, outline: Outline, dx: []f32, dy: []f32) Error!void {
    std.debug.assert(dx.len == dy.len);
    @memset(dx, 0);
    @memset(dy, 0);
    if (coords.len != gvar.axisCount()) return error.WrongAxisCount;
    if (glyph >= gvar.glyphCount()) return error.GlyphOutOfRange;
    const data = gvar.getGlyphVariationData(glyph) orelse return;
    const all = data.table_data;
    if (all.len < 4 or data.dataOffset() > all.len) return error.MalformedVariationData;
    const count = data.tupleVariationCount();

    var serialized: Reader = .{ .data = data.serializedData() };
    const shared = if (count.shared_point_numbers) try readPoints(gpa, &serialized) else null;
    defer if (shared) |s| gpa.free(s);

    const n = dx.len;
    const scratch = try gpa.alloc(f32, 2 * n);
    defer gpa.free(scratch);
    const touched = try gpa.alloc(bool, n);
    defer gpa.free(touched);

    var header_at: usize = 4;
    var tuple_at = serialized.at;
    for (0..count.count) |_| {
        if (header_at + 4 > all.len) return error.MalformedVariationData;
        const header: tv.TupleVariationHeader = .{ .table_data = all[header_at..], .axis_count = data.axis_count };
        const header_len = header.headerByteSize();
        if (header_at + header_len > all.len) return error.MalformedVariationData;
        header_at += header_len;
        const size = header.variationDataSize();
        if (tuple_at + size > serialized.data.len) return error.MalformedVariationData;
        defer tuple_at += size;

        const index = header.tupleIndex();
        const peak = header.getPeakTuple() orelse blk: {
            if (index.tuple_index >= gvar.sharedTupleCount()) return error.MalformedVariationData;
            const tuple_end = gvar.sharedTuplesOffset() + (@as(usize, index.tuple_index) + 1) * @as(usize, data.axis_count) * @sizeOf(F2Dot14);
            if (tuple_end > gvar.table_data.len) return error.MalformedVariationData;
            break :blk gvar.getSharedTuple(index.tuple_index);
        };
        const s = scalar(coords, peak, header.getIntermediateStartTuple(), header.getIntermediateEndTuple());
        if (s == 0) continue;

        var r: Reader = .{ .data = serialized.data[tuple_at..][0..size] };
        const private = if (index.private_point_numbers) try readPoints(gpa, &r) else null;
        defer if (private) |p| gpa.free(p);
        const points: ?[]const u16 = if (index.private_point_numbers) private else shared;

        const tx = scratch[0..n];
        const ty = scratch[n..];
        if (points) |pts| {
            @memset(tx, 0);
            @memset(ty, 0);
            @memset(touched, false);
            const deltas = try gpa.alloc(f32, 2 * pts.len);
            defer gpa.free(deltas);
            try readDeltas(&r, deltas[0..pts.len]);
            try readDeltas(&r, deltas[pts.len..]);
            for (pts, deltas[0..pts.len], deltas[pts.len..]) |p, x, y| {
                if (p >= n) continue; // past the phantom points: ignored, as FreeType does
                tx[p] = x;
                ty[p] = y;
                touched[p] = true;
            }
            const real = outline.xs.len;
            if (real > 0 and real <= n) iup.infer(outline.ends, outline.xs, outline.ys, tx[0..real], ty[0..real], touched[0..real]);
        } else {
            try readDeltas(&r, tx);
            try readDeltas(&r, ty);
        }
        for (dx, dy, tx, ty) |*ox, *oy, x, y| {
            ox.* += s * x;
            oy.* += s * y;
        }
    }
}

/// How much one tuple applies at `coords`: the product over its axes of
/// each axis's tent, 0 outside its region. An intermediate region whose
/// start, peak and end are out of order, or that crosses 0, is ignored on
/// that axis, as FreeType does.
fn scalar(coords: []const F2Dot14, peak: tv.Tuple, start: ?tv.Tuple, end: ?tv.Tuple) f32 {
    var s: f32 = 1;
    for (coords, 0..) |c, i| {
        const p = raw(peak.getValue(i));
        if (p == 0) continue;
        const v = raw(c);
        if (start) |st_tuple| {
            const st = raw(st_tuple.getValue(i));
            const en = raw(end.?.getValue(i));
            if (st > p or p > en or (st < 0 and en > 0)) continue;
            if (v < st or v > en) return 0;
            if (v == p) continue;
            s *= if (v < p) (v - st) / (p - st) else (en - v) / (en - p);
        } else {
            if (v == 0 or v < @min(0, p) or v > @max(0, p)) return 0;
            if (v == p) continue;
            s *= v / p;
        }
    }
    return s;
}

fn raw(v: F2Dot14) f32 {
    return @floatFromInt(@as(i16, @bitCast(v)));
}

const Reader = struct {
    data: []const u8,
    at: usize = 0,

    fn byte(r: *Reader) Error!u8 {
        if (r.at >= r.data.len) return error.MalformedVariationData;
        defer r.at += 1;
        return r.data[r.at];
    }

    fn int(r: *Reader, comptime T: type) Error!T {
        if (r.at + @sizeOf(T) > r.data.len) return error.MalformedVariationData;
        defer r.at += @sizeOf(T);
        return bytes.readIntBE(T, r.data[r.at..]);
    }
};

/// Packed point numbers; null for "every point" (a count of 0).
fn readPoints(gpa: Allocator, r: *Reader) Error!?[]u16 {
    const first = try r.byte();
    if (first == 0) return null;
    const count: usize = if (first & 0x80 != 0) (@as(usize, first & 0x7f) << 8) | try r.byte() else first;
    const points = try gpa.alloc(u16, count);
    errdefer gpa.free(points);
    var last: u16 = 0;
    var i: usize = 0;
    while (i < count) {
        const control = try r.byte();
        const run = @as(usize, control & 0x7f) + 1;
        if (i + run > count) return error.MalformedVariationData;
        for (points[i..][0..run]) |*p| {
            last +%= if (control & 0x80 != 0) try r.int(u16) else try r.byte();
            p.* = last;
        }
        i += run;
    }
    return points;
}

/// Packed deltas, exactly `out.len` of them.
fn readDeltas(r: *Reader, out: []f32) Error!void {
    var i: usize = 0;
    while (i < out.len) {
        const control = try r.byte();
        const run = @as(usize, control & 0x3f) + 1;
        if (i + run > out.len) return error.MalformedVariationData;
        for (out[i..][0..run]) |*d| d.* = switch (control & 0xc0) {
            0x80 => 0,
            0x40 => @floatFromInt(try r.int(i16)),
            0xc0 => @floatFromInt(try r.int(i32)),
            else => @floatFromInt(@as(i8, @bitCast(try r.byte()))),
        };
        i += run;
    }
}

test "glyphDeltas sums two tuples at half their peak, inferring untouched points" {
    var buf: [128]u8 = undefined;
    const gvar = try GvarTable.init(gvar_example.twoTuples(&buf), true);
    const outline: Outline = .{
        .xs = &.{ 0, 100, 100, 0 },
        .ys = &.{ 0, 0, 100, 100 },
        .ends = &.{3},
    };
    var dx: [8]f32 = undefined;
    var dy: [8]f32 = undefined;
    const gpa = zigstd.testAllocator;
    try glyphDeltas(gpa, gvar, 0, &.{F2Dot14.fromFloat(0.5)}, outline, &dx, &dy);
    // tuple 1 at 0.5: (5, 2) (10, 2) (10, 4) (5, 4); tuple 2: x + 0.5 everywhere
    try zigstd.testExpectEqualSlices(f32, &.{ 5.5, 10.5, 10.5, 5.5, 0.5, 0.5, 0.5, 0.5 }, &dx);
    try zigstd.testExpectEqualSlices(f32, &.{ 2, 2, 4, 4, 0, 0, 0, 0 }, &dy);
    // the default, and the other side of it, see none
    for ([_]f64{ 0, -0.5 }) |c| {
        try glyphDeltas(gpa, gvar, 0, &.{F2Dot14.fromFloat(c)}, outline, &dx, &dy);
        const zero: [8]f32 = @splat(0);
        try zigstd.testExpectEqualSlices(f32, &zero, &dx);
    }
    try zigstd.testExpectError(error.WrongAxisCount, glyphDeltas(gpa, gvar, 0, &.{}, outline, &dx, &dy));
    try zigstd.testExpectError(error.GlyphOutOfRange, glyphDeltas(gpa, gvar, 1, &.{F2Dot14.fromFloat(1.0)}, outline, &dx, &dy));
}

test "scalar: a tent from 0 to its peak, and an intermediate region" {
    var peak_bytes: [2]u8 = undefined;
    bytes.setI16BE(&peak_bytes, 0x4000);
    const peak: tv.Tuple = .{ .table_data = &peak_bytes, .axis_count = 1 };
    const half = F2Dot14.fromFloat(0.5);
    try zigstd.testExpectEqual(@as(f32, 0.5), scalar(&.{half}, peak, null, null));
    try zigstd.testExpectEqual(@as(f32, 1), scalar(&.{F2Dot14.fromFloat(1.0)}, peak, null, null));
    try zigstd.testExpectEqual(@as(f32, 0), scalar(&.{F2Dot14.fromFloat(-0.25)}, peak, null, null));
    // region 0.25 .. 1 .. 1: at 0.625, halfway up from 0.25
    var start_bytes: [2]u8 = undefined;
    bytes.setI16BE(&start_bytes, 0x1000);
    const start: tv.Tuple = .{ .table_data = &start_bytes, .axis_count = 1 };
    try zigstd.testExpectEqual(@as(f32, 0.5), scalar(&.{F2Dot14.fromFloat(0.625)}, peak, start, peak));
    try zigstd.testExpectEqual(@as(f32, 0), scalar(&.{F2Dot14.fromFloat(0.125)}, peak, start, peak));
}
