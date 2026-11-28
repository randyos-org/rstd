//! The deltas a `gvar` tuple leaves out, inferred from its neighbors as
//! FreeType 2.14.3 infers them (`src/truetype/ttgxvar.c`
//! `tt_interpolate_deltas`, `tt_delta_interpolate`, `tt_delta_shift`), in
//! 16.16 font units.

const fixed = @import("../fixed.zig");

pub const Point = [2]i32;

/// Moves the points of each contour (`ends`, their last points) that have
/// no delta, between `org` and `out`: shifted with a contour's only point
/// that has one, else interpolated between its neighbors that do.
pub fn interpolate(ends: []const u16, out: []Point, org: []const Point, has_delta: []const bool) void {
    var point: usize = 0;
    for (ends) |end_u16| {
        const end: usize = end_u16;
        const first_point = point;
        while (point <= end and !has_delta[point]) point += 1;
        if (point <= end) {
            const first_delta = point;
            var cur_delta = point;
            point += 1;
            while (point <= end) : (point += 1) {
                if (has_delta[point]) {
                    interp(cur_delta + 1, point, cur_delta, point, org, out);
                    cur_delta = point;
                }
            }
            if (cur_delta == first_delta) {
                shift(first_point, end, cur_delta, org, out);
            } else {
                interp(cur_delta + 1, end + 1, cur_delta, first_delta, org, out);
                if (first_delta > 0) interp(first_point, first_delta, cur_delta, first_delta, org, out);
            }
        }
        // FreeType walks on from wherever this contour's search stopped
        point = @max(point, end + 1);
    }
}

/// `tt_delta_shift`: points `p1..=p2` but `ref` moved by `ref`'s delta.
fn shift(p1: usize, p2: usize, ref: usize, org: []const Point, out: []Point) void {
    const d: Point = .{ out[ref][0] -% org[ref][0], out[ref][1] -% org[ref][1] };
    if (d[0] == 0 and d[1] == 0) return;
    for (p1..p2 + 1) |p| {
        if (p == ref) continue;
        out[p] = .{ out[p][0] +% d[0], out[p][1] +% d[1] };
    }
}

/// `tt_delta_interpolate`: points `from..to` (exclusive) placed between
/// `ref1` and `ref2`, each axis on its own.
fn interp(from: usize, to: usize, ref1_in: usize, ref2_in: usize, org: []const Point, out: []Point) void {
    if (from >= to) return;
    for (0..2) |i| {
        var ref1 = ref1_in;
        var ref2 = ref2_in;
        if (org[ref1][i] > org[ref2][i]) {
            ref1 = ref2_in;
            ref2 = ref1_in;
        }
        const in1 = org[ref1][i];
        const in2 = org[ref2][i];
        const out1 = out[ref1][i];
        const out2 = out[ref2][i];
        // same coordinate but different deltas: no inferred delta
        if (in1 == in2 and out1 != out2) continue;
        const scale = if (in1 != in2) fixed.div(out2 -% out1, in2 -% in1) else 0;
        for (from..to) |p| {
            const v = org[p][i];
            out[p][i] = if (v <= in1) v +% (out1 -% in1) else if (v >= in2) v +% (out2 -% in2) else out1 +% fixed.mul(v -% in1, scale);
        }
    }
}
