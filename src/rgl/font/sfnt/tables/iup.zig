//! Inferred deltas for the points a `gvar` tuple leaves out (OpenType 1.9.1,
//! "gvar", "Inferred deltas for un-referenced point numbers"), contour by
//! contour, each axis on its own, as FreeType's `tt_interpolate_deltas`.

const zigstd = @import("../../../../zigstd.zig");

/// Fills `dx` and `dy` for the points `touched` leaves false, from the
/// touched points around them on their contour. `xs` and `ys` are the
/// original coordinates; `ends` each contour's last point index. A contour
/// with no touched point keeps its zero deltas; with one, all take its delta.
pub fn infer(ends: []const u16, xs: []const i16, ys: []const i16, dx: []f32, dy: []f32, touched: []const bool) void {
    var start: usize = 0;
    for (ends) |e| {
        const end: usize = e;
        if (end < start or end >= xs.len) return; // malformed: stop at the last good contour
        defer start = end + 1;
        const first = for (start..end + 1) |p| {
            if (touched[p]) break p;
        } else continue;
        // walk each run between consecutive touched points, wrapping at the contour's end
        var a = first;
        while (true) {
            var b = next(a, start, end);
            while (!touched[b]) b = next(b, start, end);
            var p = next(a, start, end);
            while (p != b) : (p = next(p, start, end)) {
                dx[p] = interp(xs, dx, a, b, p);
                dy[p] = interp(ys, dy, a, b, p);
            }
            if (b == first) break;
            a = b;
        }
    }
}

fn next(p: usize, start: usize, end: usize) usize {
    return if (p == end) start else p + 1;
}

/// Point `p`'s delta along one axis from touched points `a` and `b`: the
/// nearer one's beyond their range, linear between, none if they coincide
/// with different deltas. With `a == b` (one touched point) it is `a`'s.
fn interp(orig: []const i16, d: []const f32, a: usize, b: usize, p: usize) f32 {
    var in1: f32 = @floatFromInt(orig[a]);
    var d1 = d[a];
    var in2: f32 = @floatFromInt(orig[b]);
    var d2 = d[b];
    if (in1 > in2) {
        const t = in1;
        in1 = in2;
        in2 = t;
        const u = d1;
        d1 = d2;
        d2 = u;
    }
    if (in1 == in2) return if (d1 == d2) d1 else 0;
    const c: f32 = @floatFromInt(orig[p]);
    if (c <= in1) return d1;
    if (c >= in2) return d2;
    return d1 + (c - in1) * (d2 - d1) / (in2 - in1);
}

test "infer interpolates between touched points, takes the nearer delta past them" {
    // a square (0,0) (100,0) (100,100) (0,100), then a lone point
    const ends = [_]u16{ 3, 4 };
    const xs = [_]i16{ 0, 100, 100, 0, 50 };
    const ys = [_]i16{ 0, 0, 100, 100, 50 };
    var dx = [_]f32{ 10, 0, 20, 0, 0 };
    var dy = [_]f32{ 4, 0, 8, 0, 0 };
    const touched = [_]bool{ true, false, true, false, false };
    infer(&ends, &xs, &ys, &dx, &dy, &touched);
    try zigstd.testExpectEqualSlices(f32, &.{ 10, 20, 20, 10, 0 }, &dx);
    try zigstd.testExpectEqualSlices(f32, &.{ 4, 4, 8, 8, 0 }, &dy);
}

test "infer gives a contour's one touched delta to all its points, and splits a range" {
    const ends = [_]u16{2};
    const xs = [_]i16{ 0, 50, 100 };
    const ys = [_]i16{ 0, 0, 0 };
    var dx = [_]f32{ 0, 6, 0 };
    var dy = [_]f32{ 0, 2, 0 };
    infer(&ends, &xs, &ys, &dx, &dy, &.{ false, true, false });
    try zigstd.testExpectEqualSlices(f32, &.{ 6, 6, 6 }, &dx);
    try zigstd.testExpectEqualSlices(f32, &.{ 2, 2, 2 }, &dy);
    // touched ends of a line: the middle point lies halfway
    var ex = [_]f32{ 0, 0, 10 };
    var ey = [_]f32{ 0, 0, 0 };
    infer(&ends, &xs, &ys, &ex, &ey, &.{ true, false, true });
    try zigstd.testExpectEqualSlices(f32, &.{ 0, 5, 10 }, &ex);
}
