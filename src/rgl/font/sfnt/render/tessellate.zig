const std = @import("std");
const zigstd = @import("../../../../zigstd.zig");
const Allocator = zigstd.Allocator;

const Vertex = @import("../aux_types.zig").Vertex;

pub const Point = struct { x: f32, y: f32 };

pub const FlattenedCurves = struct {
    points: []Point,
    /// Point count per contour (a moveTo starts a new one); sums to
    /// points.len.
    contour_lengths: []u32,

    pub const empty: FlattenedCurves = .{
        .points = &.{},
        .contour_lengths = &.{},
    };

    pub fn deinit(fc: *FlattenedCurves, gpa: Allocator) void {
        gpa.free(fc.points);
        gpa.free(fc.contour_lengths);
        fc.* = undefined;
    }
};

pub fn flattenCurves(
    gpa: Allocator,
    vertices: []const Vertex,
    objspace_flatness: f32,
) error{OutOfMemory}!FlattenedCurves {
    var points: std.ArrayList(Point) = .empty;
    defer points.deinit(gpa);
    var contour_lengths: std.ArrayList(u32) = .empty;
    defer contour_lengths.deinit(gpa);

    const objspace_flatness_squared = objspace_flatness * objspace_flatness;

    var start: u32 = 0;
    var x: f32 = 0;
    var y: f32 = 0;
    for (vertices) |v| {
        sw: switch (v.type) {
            .vmove => {
                if (points.items.len > 0) {
                    try contour_lengths.append(gpa, @intCast(points.items.len - start));
                    start = @intCast(points.items.len);
                }

                continue :sw .vline;
            },
            .vline => {
                x = @floatFromInt(v.x);
                y = @floatFromInt(v.y);
                try points.append(gpa, .{ .x = x, .y = y });
            },
            .vcurve => {
                try tesselateCurve(
                    gpa,
                    &points,
                    x,
                    y,
                    @floatFromInt(v.cx),
                    @floatFromInt(v.cy),
                    @floatFromInt(v.x),
                    @floatFromInt(v.y),
                    objspace_flatness_squared,
                    0,
                );
                x = @floatFromInt(v.x);
                y = @floatFromInt(v.y);
            },
            .vcubic => {
                try tesselateCubic(
                    gpa,
                    &points,
                    x,
                    y,
                    @floatFromInt(v.cx),
                    @floatFromInt(v.cy),
                    @floatFromInt(v.cx1),
                    @floatFromInt(v.cy1),
                    @floatFromInt(v.x),
                    @floatFromInt(v.y),
                    objspace_flatness_squared,
                    0,
                );
                x = @floatFromInt(v.x);
                y = @floatFromInt(v.y);
            },
            _ => continue,
        }
    }
    try contour_lengths.append(gpa, @intCast(points.items.len - start));

    return .{
        .points = try points.toOwnedSlice(gpa),
        .contour_lengths = try contour_lengths.toOwnedSlice(gpa),
    };
}

/// tessellate until threshold p is happy... @TODO warped to compensate for non-linear stretching
fn tesselateCurve(
    gpa: Allocator,
    points: *std.ArrayList(Point),
    x0: f32,
    y0: f32,
    x1: f32,
    y1: f32,
    x2: f32,
    y2: f32,
    objspace_flatness_squared: f32,
    n: u32,
) Allocator.Error!void {
    // midpoint
    const mx: f32 = (x0 + 2 * x1 + x2) / 4;
    const my: f32 = (y0 + 2 * y1 + y2) / 4;
    // versus directly drawn line
    const dx: f32 = (x0 + x2) / 2 - mx;
    const dy: f32 = (y0 + y2) / 2 - my;
    if (n > 16) // 65536 segments on one curve better be enough!
        return;
    if (dx * dx + dy * dy > objspace_flatness_squared) { // half-pixel error allowed... need to be smaller if AA
        try tesselateCurve(gpa, points, x0, y0, (x0 + x1) / 2.0, (y0 + y1) / 2.0, mx, my, objspace_flatness_squared, n + 1);
        try tesselateCurve(gpa, points, mx, my, (x1 + x2) / 2.0, (y1 + y2) / 2.0, x2, y2, objspace_flatness_squared, n + 1);
    } else {
        try points.append(gpa, .{ .x = x2, .y = y2 });
    }
}

fn tesselateCubic(
    gpa: Allocator,
    points: *std.ArrayList(Point),
    x0: f32,
    y0: f32,
    x1: f32,
    y1: f32,
    x2: f32,
    y2: f32,
    x3: f32,
    y3: f32,
    objspace_flatness_squared: f32,
    n: u32,
) Allocator.Error!void {
    // According to Dougall Johnson, this "flatness" calculation is just
    // made-up nonsense that seems to work well enough.
    const dx0 = x1 - x0;
    const dy0 = y1 - y0;
    const dx1 = x2 - x1;
    const dy1 = y2 - y1;
    const dx2 = x3 - x2;
    const dy2 = y3 - y2;
    const dx = x3 - x0;
    const dy = y3 - y0;
    const longlen = @sqrt(dx0 * dx0 + dy0 * dy0) + @sqrt(dx1 * dx1 + dy1 * dy1) + @sqrt(dx2 * dx2 + dy2 * dy2);
    const shortlen = @sqrt(dx * dx + dy * dy);
    const flatness_squared = longlen * longlen - shortlen * shortlen;

    if (n > 16) // 65536 segments on one curve better be enough!
        return;

    if (flatness_squared > objspace_flatness_squared) {
        const x01 = (x0 + x1) / 2;
        const y01 = (y0 + y1) / 2;
        const x12 = (x1 + x2) / 2;
        const y12 = (y1 + y2) / 2;
        const x23 = (x2 + x3) / 2;
        const y23 = (y2 + y3) / 2;

        const xa = (x01 + x12) / 2;
        const ya = (y01 + y12) / 2;
        const xb = (x12 + x23) / 2;
        const yb = (y12 + y23) / 2;

        const mx = (xa + xb) / 2;
        const my = (ya + yb) / 2;

        try tesselateCubic(gpa, points, x0, y0, x01, y01, xa, ya, mx, my, objspace_flatness_squared, n + 1);
        try tesselateCubic(gpa, points, mx, my, xb, yb, x23, y23, x3, y3, objspace_flatness_squared, n + 1);
    } else {
        try points.append(gpa, .{ .x = x3, .y = y3 });
    }
}

test "flattenCurves traces a simple triangle as 3 points in one contour" {
    const alloc = zigstd.testAllocator;
    var segments: [4]Vertex = undefined;
    segments[0].set(.vmove, 0, 0, 0, 0);
    segments[1].set(.vline, 10, 0, 0, 0);
    segments[2].set(.vline, 5, 10, 0, 0);
    segments[3].set(.vline, 0, 0, 0, 0);

    var flattened = try flattenCurves(alloc, &segments, 0.1);
    defer flattened.deinit(alloc);

    try zigstd.testExpectEqual(@as(usize, 1), flattened.contour_lengths.len);
    try zigstd.testExpectEqual(@as(u32, 4), flattened.contour_lengths[0]);
    try zigstd.testExpectEqual(@as(f32, 5), flattened.points[2].x);
    try zigstd.testExpectEqual(@as(f32, 10), flattened.points[2].y);
}

test "flattenCurves subdivides a curved vcurve into more than its two endpoints" {
    const alloc = zigstd.testAllocator;
    var segments: [2]Vertex = undefined;
    segments[0].set(.vmove, 0, 0, 0, 0);
    segments[1].set(.vcurve, 100, 0, 50, 100); // end (100,0), control (50,100)

    var flattened = try flattenCurves(alloc, &segments, 0.1);
    defer flattened.deinit(alloc);

    try zigstd.testExpectEqual(@as(usize, 1), flattened.contour_lengths.len);
    try zigstd.testExpectTrue(flattened.contour_lengths[0] > 2);
    try zigstd.testExpectEqual(@as(f32, 0), flattened.points[0].x);
    const last = flattened.points[flattened.points.len - 1];
    try zigstd.testExpectEqual(@as(f32, 100), last.x);
    try zigstd.testExpectEqual(@as(f32, 0), last.y);
}
