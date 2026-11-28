//! FreeType 2.14's diacritic handling in the Latin hinter: tildes
//! stretched to stay legible, accents kept a pixel from their base, blue
//! zones ignored above or below the base. FreeType 2.14.3
//! `src/autofit/aflatin.c` (`af_move_contour_vertically` through
//! `af_glyph_hints_apply_vertical_separation_adjustments`), its quirks kept.

const std = @import("std");
const GlyphHints = @import("../GlyphHints.zig");
const Metrics = @import("Metrics.zig");
const adjust = @import("../adjust.zig");

const flag = GlyphHints.flag;
const blue_flag = Metrics.blue_flag;
const long_max = std.math.maxInt(i64);
const long_min = std.math.minInt(i64);

/// Each contour's vertical extent, 64-bit as `FT_Pos`
/// (`contour_y_minima`, `contour_y_maxima`).
pub const Extrema = struct {
    min: []i64,
    max: []i64,

    /// `af_compute_vertical_extrema`: contours of under three points span
    /// nothing.
    pub fn compute(e: Extrema, h: *const GlyphHints) void {
        const pts = h.points.items;
        for (h.contours.items, 0..) |first, c| {
            var lo: i64 = long_max;
            var hi: i64 = long_min;
            if (pts[pts[first].next].next != first) {
                var p = first;
                while (true) {
                    lo = @min(lo, pts[p].y);
                    hi = @max(hi, pts[p].y);
                    p = pts[p].next;
                    if (p == first) break;
                }
            }
            e.min[c] = lo;
            e.max[c] = hi;
        }
    }
};

/// Calls `f` on each point of contour `c`, its first point last.
fn eachPoint(h: *GlyphHints, c: usize, ctx: anytype, comptime f: fn (@TypeOf(ctx), *GlyphHints.Point) void) void {
    const pts = h.points.items;
    const first = h.contours.items[c];
    var p = first;
    while (true) {
        p = pts[p].next;
        f(ctx, &pts[p]);
        if (p == first) break;
    }
}

fn moveContour(h: *GlyphHints, c: usize, movement: i64) void {
    eachPoint(h, c, movement, struct {
        fn f(m: i64, p: *GlyphHints.Point) void {
            p.y +%= @truncate(m);
        }
    }.f);
}

/// `af_move_contours_up`: every contour above `limit`.
fn moveUp(h: *GlyphHints, e: Extrema, limit: i64, delta: i64) void {
    for (e.min, e.max, 0..) |lo, hi, c| {
        if (lo < hi and lo > limit) moveContour(h, c, delta);
    }
}

fn moveDown(h: *GlyphHints, e: Extrema, limit: i64, delta: i64) void {
    for (e.min, e.max, 0..) |lo, hi, c| {
        if (lo < hi and hi < limit) moveContour(h, c, -delta);
    }
}

/// `af_find_highest_contour`: the highest top, the shorter on a tie.
fn highest(e: Extrema) usize {
    var best: usize = 0;
    var best_min: i64 = long_max;
    var best_max: i64 = long_min;
    for (e.min, e.max, 0..) |lo, hi, c| {
        if (hi > best_max or (hi == best_max and lo > best_min)) {
            best_min = lo;
            best_max = hi;
            best = c;
        }
    }
    return best;
}

/// `af_find_second_highest_contour`: the highest one reaching below the
/// highest contour's bottom.
fn secondHighest(e: Extrema) usize {
    if (e.min.len < 3) return 0;
    const top = highest(e);
    const top_min = e.min[top];
    var best: usize = 0;
    var best_max: i64 = long_min;
    for (e.min, e.max, 0..) |lo, hi, c| {
        if (c == top) continue;
        if (hi > best_max and lo < top_min) {
            best_max = hi;
            best = c;
        }
    }
    return best;
}

fn lowest(e: Extrema) usize {
    var best: usize = 0;
    var best_min: i64 = long_max;
    var best_max: i64 = long_min;
    for (e.min, e.max, 0..) |lo, hi, c| {
        if (lo < best_min or (lo == best_min and hi < best_max)) {
            best_min = lo;
            best_max = hi;
            best = c;
        }
    }
    return best;
}

fn secondLowest(e: Extrema) usize {
    if (e.min.len < 3) return 0;
    const bottom = lowest(e);
    const bottom_max = e.max[bottom];
    var best: usize = 0;
    var best_min: i64 = long_max;
    for (e.min, e.max, 0..) |lo, hi, c| {
        if (c == bottom) continue;
        if (lo < best_min and hi > bottom_max) {
            best_min = lo;
            best = c;
        }
    }
    return best;
}

/// `af_latin_get_base_glyph_blues`: the active capital (or small) top and
/// bottom zones, if any.
pub fn baseBlues(m: *const Metrics, capital: bool) [2]?*const Metrics.Blue {
    var r: [2]?*const Metrics.Blue = .{ null, null };
    const zones = m.axis[1].blueSlice();
    const want = [2]u32{
        (if (capital) blue_flag.top else blue_flag.adjustment) | blue_flag.active,
        (if (capital) blue_flag.bottom else blue_flag.bottom_small) | blue_flag.active,
    };
    for (want, 0..) |w, k| {
        for (zones) |*z| if (z.flags & w == w) {
            r[k] = z;
            break;
        };
    }
    return r;
}

/// `af_latin_ignore_top`: edges above the base, past a seventh of its
/// height and a quarter pixel, align to no blue zone.
pub fn ignoreTop(h: *GlyphHints, top: *const Metrics.Blue, bottom: *const Metrics.Blue) void {
    const height = top.shoot.cur - bottom.shoot.cur;
    const limit = top.shoot.cur + @divTrunc(height, 7) + 16;
    for (h.axis[1].edges.items) |*e| if (e.pos > limit) {
        e.flags |= GlyphHints.edge_flag.no_blue;
    };
}

pub fn ignoreBottom(h: *GlyphHints, top: *const Metrics.Blue, bottom: *const Metrics.Blue) void {
    const height = top.shoot.cur - bottom.shoot.cur;
    const limit = bottom.shoot.cur - @divTrunc(height, 7) - 16;
    for (h.axis[1].edges.items) |*e| if (e.pos < limit) {
        e.flags |= GlyphHints.edge_flag.no_blue;
    };
}

/// `af_touch_contour`: kept out of segments and vertical hinting.
fn touchContour(h: *GlyphHints, c: usize) void {
    eachPoint(h, c, {}, struct {
        fn f(_: void, p: *GlyphHints.Point) void {
            p.flags |= flag.ignore;
            if (p.flags & flag.control == 0) p.flags |= flag.touch_y;
        }
    }.f);
}

/// The lowest point of a tilde's wiggles off its extremes, as a distance
/// from the extreme it faces; 0 if none.
fn tildeMeasurement(h: *GlyphHints, c: usize, lo: i64, hi: i64) i64 {
    const pts = h.points.items;
    const first = h.contours.items[c];
    const threshold = @divTrunc(hi - lo, 8);
    var min_measurement: i64 = long_max;
    var taken = false;
    var p = first;
    while (true) {
        p = pts[p].next;
        const pt = pts[p];
        const prev = pts[pt.prev];
        const next = pts[pt.next];
        if (pt.flags & flag.control == 0 and prev.y == pt.y and next.y == pt.y and pt.y != lo and pt.y != hi and
            prev.flags & flag.control != 0 and next.flags & flag.control != 0)
        candidate: {
            var prev_on = pt.prev;
            while (pts[prev_on].flags & flag.control != 0) prev_on = pts[prev_on].prev;
            var next_on = pt.next;
            while (pts[next_on].flags & flag.control != 0) next_on = pts[next_on].next;
            const ny = pts[next_on].y;
            const py = pts[prev_on].y;
            const measurement: i64 = if (ny > pt.y and py > pt.y) pt.y - lo else if (ny < pt.y and py < pt.y) hi - pt.y else break :candidate;
            if (measurement < threshold) break :candidate;
            if (!taken or measurement < min_measurement) {
                taken = true;
                min_measurement = measurement;
            }
        }
        if (p == first) break;
    }
    return if (taken) min_measurement else 0;
}

/// `af_latin_stretch_top_tilde` (`top`) or `_bottom_tilde`: the tilde made
/// a pixel taller than its wiggle, it and the contours beyond it touched;
/// the height it gained.
fn stretchTilde(h: *GlyphHints, e: Extrema, c: usize, top: bool) i64 {
    const lo = e.min[c];
    const hi = e.max[c];
    if (lo == hi) return 0;
    const height = hi - lo;
    const min_measurement = tildeMeasurement(h, c, lo, hi);
    // the tilde and what lies beyond it keep this shape
    const limit = e.min[c];
    for (e.min, e.max, 0..) |clo, chi, k| {
        if (clo < chi and (if (top) clo >= limit else chi <= limit)) touchContour(h, k);
    }
    const target = min_measurement + 64;
    if (height >= target) return 0;
    const pts = h.points.items;
    const first = h.contours.items[c];
    const base = if (top) lo else hi;
    var p = first;
    while (true) {
        p = pts[p].next;
        pts[p].y = @intCast(@divTrunc((pts[p].y - base) * target, height) + base);
        if (p == first) break;
    }
    return target - height;
}

/// `af_latin_align_top_tilde` (`top`) or `_bottom_tilde`: the tilde's base
/// (or, under three pixels, its center) on the grid; the shift.
fn alignTilde(h: *GlyphHints, c: usize, top: bool) i64 {
    const pts = h.points.items;
    const first = h.contours.items[c];
    var lo: i64 = pts[first].y;
    var hi: i64 = lo;
    var p = first;
    while (true) {
        p = pts[p].next;
        lo = @min(lo, pts[p].y);
        hi = @max(hi, pts[p].y);
        if (p == first) break;
    }
    const height = hi - lo;
    var delta: i64 = if (top) pixRound(lo) - lo else pixRound(hi) - hi;
    if (height < 64 * 3) {
        const centering = @divTrunc(pixRound(height) - height, 2);
        delta = if (top) delta + centering else delta - centering;
    }
    moveContour(h, c, delta);
    return delta;
}

fn pixRound(x: i64) i64 {
    return (x + 32) & ~@as(i64, 63);
}

/// The tildes, before the vertical features are found: inner ones first.
pub fn stretchTildes(h: *GlyphHints, e: Extrema, flags: u32) void {
    if (flags & adjust.tilde_top2 != 0) {
        const c = secondHighest(e);
        var y_offset = stretchTilde(h, e, c, true);
        y_offset += alignTilde(h, c, true);
        moveUp(h, e, e.min[c], y_offset);
    }
    if (flags & adjust.tilde_bottom2 != 0) {
        const c = secondLowest(e);
        var y_offset = stretchTilde(h, e, c, false);
        y_offset -= alignTilde(h, c, false);
        moveDown(h, e, e.max[c], y_offset);
    }
    if (flags & adjust.tilde_top != 0) {
        const c = highest(e);
        _ = stretchTilde(h, e, c, true);
        _ = alignTilde(h, c, true);
    }
    if (flags & adjust.tilde_bottom != 0) {
        const c = lowest(e);
        _ = stretchTilde(h, e, c, false);
        _ = alignTilde(h, c, false);
    }
}

/// `af_check_contour_horizontal_overlap`: whether contour `ci` overlaps the
/// others' horizontal extent.
fn overlapsOthers(h: *const GlyphHints, ci: usize) bool {
    const pts = h.points.items;
    var c_min: i64 = long_max;
    var c_max: i64 = long_min;
    var o_min: i64 = long_max;
    var o_max: i64 = long_min;
    for (h.contours.items, 0..) |first, c| {
        if (pts[pts[first].next].next == first) continue;
        var p = first;
        while (true) {
            p = pts[p].next;
            const x: i64 = pts[p].x;
            if (c == ci) {
                c_min = @min(c_min, x);
                c_max = @max(c_max, x);
            } else {
                o_min = @min(o_min, x);
                o_max = @max(o_max, x);
            }
            if (p == first) break;
        }
    }
    return (o_min <= c_max and c_max <= o_max) or (o_min <= c_min and c_min <= o_max) or (c_max >= o_max and c_min <= o_min);
}

/// `af_glyph_hints_apply_vertical_separation_adjustments`: an accent above
/// (or below) its base pushed to a pixel away from it.
pub fn separate(h: *GlyphHints, e: Extrema, flags: u32, accent_height_limit: i64) void {
    const top = flags & adjust.up != 0;
    const below_top = flags & adjust.up2 != 0;
    const bottom = flags & adjust.down != 0;
    const above_bottom = flags & adjust.down2 != 0;
    const n = h.contours.items.len;
    if (((top or bottom) and n >= 2) or ((below_top or above_bottom) and n >= 3)) e.compute(h);
    if ((top and n >= 2) or (below_top and n >= 3)) {
        const high = if (below_top) secondHighest(e) else highest(e);
        if (!overlapsOthers(h, high)) return;
        const high_min = e.min[high];
        const high_height = e.max[high] - high_min;
        if (high_height > accent_height_limit) return;
        var min_distance: i64 = 64;
        for (e.min, e.max, 0..) |lo, hi, c| {
            if (c == high) continue;
            // the other contour must reach below the accent
            const distance = high_min - hi;
            if (distance < 64 and distance < min_distance and lo < high_min) min_distance = distance;
        }
        const amount = 64 - min_distance;
        const is_tilde = flags & adjust.tilde_top != 0;
        const is_below_tilde = flags & adjust.tilde_top2 != 0;
        var centering: i64 = 0;
        if (is_tilde or is_below_tilde) {
            const t = if (top or is_below_tilde) high else highest(e);
            const t_height = e.max[t] - e.min[t];
            var pos = high_min + amount;
            if (below_top and is_tilde) pos += high_height;
            if (@rem(pos, 64) == 0 and t_height < 3 * 64) centering = @divTrunc(pixRound(t_height) - t_height, 2);
        }
        const calculated = if ((top and is_tilde) or (below_top and is_below_tilde)) amount + centering else amount;
        // 2/64 pixel of slack for rounding
        if (calculated != 0 and calculated >= -2 and (calculated <= 66 or amount <= 66)) {
            const min_y_limit = high_min - @divTrunc(high_height, 8);
            moveUp(h, e, min_y_limit, amount);
            if (below_top and is_tilde) moveUp(h, e, min_y_limit + high_height, centering);
        }
    }
    if ((bottom and n >= 2) or (above_bottom and n >= 3)) {
        const low = if (above_bottom) secondLowest(e) else lowest(e);
        if (!overlapsOthers(h, low)) return;
        const low_max = e.max[low];
        const low_height = low_max - e.min[low];
        if (low_height > accent_height_limit) return;
        var min_distance: i64 = 64;
        for (e.min, e.max, 0..) |lo, hi, c| {
            if (c == low) continue;
            const distance = lo - low_max;
            if (distance < 64 and distance < min_distance and hi > low_max) min_distance = distance;
        }
        const amount = 64 - min_distance;
        const is_tilde = flags & adjust.tilde_bottom != 0;
        const is_above_tilde = flags & adjust.tilde_bottom2 != 0;
        var centering: i64 = 0;
        if (is_tilde or is_above_tilde) {
            const t = if (bottom or is_above_tilde) low else lowest(e);
            const t_height = e.max[t] - e.min[t];
            var pos = low_max - amount;
            if (above_bottom and is_tilde) pos -= low_height;
            if (@rem(pos, 64) == 0 and t_height < 3 * 64) centering = @divTrunc(pixRound(t_height) - t_height, 2);
        }
        const calculated = if ((bottom and is_tilde) or (above_bottom and is_above_tilde)) amount + centering else amount;
        if (calculated != 0 and calculated >= -2 and (calculated <= 66 or amount <= 66)) {
            const max_y_limit = low_max + @divTrunc(low_height, 8);
            moveDown(h, e, max_y_limit, amount);
            if (above_bottom and is_tilde) moveDown(h, e, max_y_limit - low_height, centering);
        }
    }
}
