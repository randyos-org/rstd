//! A Latin-system style's blue zones: for each blue string, the extremes
//! of its characters' outlines, flat ones the reference and round ones the
//! overshoot. FreeType 2.14.3 `src/autofit/aflatin.c`
//! (`af_latin_metrics_init_blues`, `af_latin_sort_blue`), without HarfBuzz
//! (no y offsets), its quirks kept.

const std = @import("std");
const Metrics = @import("Metrics.zig");
const Source = @import("../Source.zig");
const Glyph = @import("../../Glyph.zig");
const shaper = @import("../shaper.zig");
const sortPos = @import("../widths.zig").sortPos;
const blue_strings = @import("../data/blues.zig").blue_strings;
const styles = @import("../data/styles.zig").styles;

const Allocator = std.mem.Allocator;
const blue_flag = Metrics.blue_flag;

/// `AF_BLUE_PROPERTY_LATIN_*`.
const prop = struct {
    const top: u16 = 1 << 0;
    const sub_top: u16 = 1 << 1;
    const neutral: u16 = 1 << 2;
    const x_height: u16 = 1 << 3;
    const long: u16 = 1 << 4;
    const capital_bottom: u16 = 1 << 5;
    const small_bottom: u16 = 1 << 6;
};

/// `AF_BLUE_STRING_MAX_LEN`.
const max_len = 51;
const tag_on = 1;

/// Whether any zone was found.
pub fn init(allocator: Allocator, m: *Metrics, src: *Source) Source.Error!bool {
    const axis = &m.axis[1];
    const flat_threshold = @divTrunc(@as(i32, m.units_per_em), 14);
    var bs: usize = styles[m.style].blue_stringset;
    while (blue_strings[bs].chars.len != 0) : (bs += 1) {
        const props = blue_strings[bs].properties;
        const is_top = props & (prop.top | prop.sub_top) != 0;
        var flats: [max_len]i32 = undefined;
        var rounds: [max_len]i32 = undefined;
        var num_flats: usize = 0;
        var num_rounds: usize = 0;
        var ascender: i32 = 0;
        var descender: i32 = 0;
        const chars = blue_strings[bs].chars;
        var at: usize = 0;
        while (at < chars.len) {
            const c = shaper.cluster(src.sfnt, chars, &at);
            if (c.count == 0 or c.glyph == 0) continue;
            var g = src.load(allocator, c.glyph) catch |e| switch (e) {
                error.OutOfMemory => return e,
                else => continue,
            };
            defer g.deinit(allocator);
            // glyphs drawing nothing are skipped
            if (g.points.len <= 2) continue;
            const r = extremum(m, &g, props, &ascender, &descender, flat_threshold) orelse continue;
            if (r.round) {
                rounds[num_rounds] = r.y;
                num_rounds += 1;
            } else {
                flats[num_flats] = r.y;
                num_flats += 1;
            }
        }
        // no glyph for the zone: none
        if (num_flats == 0 and num_rounds == 0) continue;
        // the medians
        sortPos(rounds[0..num_rounds]);
        sortPos(flats[0..num_flats]);
        const blue = &axis.blues[axis.blue_count];
        axis.blue_count += 1;
        if (num_flats == 0) {
            blue.ref.org = rounds[num_rounds / 2];
            blue.shoot.org = blue.ref.org;
        } else if (num_rounds == 0) {
            blue.ref.org = flats[num_flats / 2];
            blue.shoot.org = blue.ref.org;
        } else {
            blue.ref.org = flats[num_flats / 2];
            blue.shoot.org = rounds[num_rounds / 2];
        }
        // an overshoot on the wrong side of its reference: the mean of both
        if (blue.shoot.org != blue.ref.org) {
            const ref = blue.ref.org;
            const shoot = blue.shoot.org;
            if (is_top != (shoot > ref)) {
                blue.ref.org = @divTrunc(shoot + ref, 2);
                blue.shoot.org = blue.ref.org;
            }
        }
        blue.ascender = ascender;
        blue.descender = descender;
        blue.flags = 0;
        if (props & prop.top != 0) blue.flags |= blue_flag.top;
        if (props & prop.sub_top != 0) blue.flags |= blue_flag.sub_top;
        if (props & prop.neutral != 0) blue.flags |= blue_flag.neutral;
        if (props & prop.capital_bottom != 0) blue.flags |= blue_flag.bottom;
        if (props & prop.small_bottom != 0) blue.flags |= blue_flag.bottom_small;
        if (props & prop.x_height != 0) blue.flags |= blue_flag.adjustment;
    }
    if (axis.blue_count == 0) return false;
    // zones must not overlap: sorted by bottom, each top lowered to the next
    var sorted: [Metrics.max_blues]*Metrics.Blue = undefined;
    for (axis.blues[0..axis.blue_count], 0..) |*b, i| sorted[i] = b;
    sortBlue(sorted[0..axis.blue_count]);
    for (0..axis.blue_count - 1) |i| {
        const a = edge(sorted[i], true);
        const b = edge(sorted[i + 1], true);
        if (a.* > b.*) a.* = b.*;
    }
    return true;
}

/// A top zone's overshoot or a bottom zone's reference, for the overlap
/// fix; or the other one, for sorting.
fn edge(b: *Metrics.Blue, overlap: bool) *i32 {
    const top = b.flags & (blue_flag.top | blue_flag.sub_top) != 0;
    return if (top == overlap) &b.shoot.org else &b.ref.org;
}

/// `af_latin_sort_blue`: bottom to top.
fn sortBlue(table: []*Metrics.Blue) void {
    for (1..@max(table.len, 1)) |i| {
        var j = i;
        while (j > 0) : (j -= 1) {
            const a = keyOf(table[j - 1]);
            const b = keyOf(table[j]);
            if (b >= a) break;
            std.mem.swap(*Metrics.Blue, &table[j], &table[j - 1]);
        }
    }
}

fn keyOf(b: *Metrics.Blue) i32 {
    return if (b.flags & (blue_flag.top | blue_flag.sub_top) != 0) b.ref.org else b.shoot.org;
}

const Extremum = struct { y: i32, round: bool };

/// The glyph's highest (or lowest) point and whether its segment is round;
/// null to skip the glyph. Updates the zone's `ascender` and `descender`.
fn extremum(m: *const Metrics, g: *const Glyph, props: u16, ascender: *i32, descender: *i32, flat_threshold: i32) ?Extremum {
    const is_top = props & (prop.top | prop.sub_top) != 0;
    const pts = g.points;
    const tags = g.tags;
    var best_point: i32 = -1;
    var best_first: i32 = -1;
    var best_last: i32 = -1;
    var best_y: i32 = 0;
    var last: i32 = -1;
    for (g.contours) |end| {
        const first = last + 1;
        last = end;
        // single points are never drawn: mark anchors far off the outline
        if (last <= first) continue;
        var pp = first;
        while (pp <= last) : (pp += 1) {
            const y = pts[@intCast(pp)][1];
            if (is_top) {
                if (best_point < 0 or y > best_y) {
                    best_point = pp;
                    best_y = y;
                    ascender.* = @max(ascender.*, best_y);
                } else descender.* = @min(descender.*, y);
            } else {
                if (best_point < 0 or y < best_y) {
                    best_point = pp;
                    best_y = y;
                    descender.* = @min(descender.*, best_y);
                } else ascender.* = @max(ascender.*, y);
            }
        }
        if (best_point > best_last) {
            best_first = first;
            best_last = last;
        }
    }
    if (best_point < 0) return .{ .y = best_y, .round = false };

    const P = struct {
        p: []const [2]i32,
        first: i32,
        last: i32,
        fn x(s: @This(), i: i32) i32 {
            return s.p[@intCast(i)][0];
        }
        fn y(s: @This(), i: i32) i32 {
            return s.p[@intCast(i)][1];
        }
        fn prev(s: @This(), i: i32) i32 {
            return if (i > s.first) i - 1 else s.last;
        }
        fn next(s: @This(), i: i32) i32 {
            return if (i < s.last) i + 1 else s.first;
        }
    };
    const c: P = .{ .p = pts, .first = best_first, .last = best_last };
    const on = struct {
        fn f(t: []const u8, i: i32) bool {
            return t[@intCast(i)] & 3 == tag_on;
        }
    }.f;
    const best_x = c.x(best_point);
    var seg_first = best_point;
    var seg_last = best_point;
    var on_first: i32 = if (on(tags, best_point)) best_point else -1;
    var on_last: i32 = on_first;

    // the neighbors on about the same height (within 5 units, or 2.9 degrees)
    var prev = best_point;
    var next = best_point;
    while (true) {
        prev = c.prev(prev);
        const dist: i32 = @intCast(@abs(c.y(prev) - best_y));
        if (dist > 5 and @abs(c.x(prev) - best_x) <= 20 * dist) break;
        seg_first = prev;
        if (on(tags, prev)) {
            on_first = prev;
            if (on_last < 0) on_last = prev;
        }
        if (prev == best_point) break;
    }
    var dist: i32 = 0;
    while (true) {
        next = c.next(next);
        dist = @intCast(@abs(c.y(next) - best_y));
        if (dist > 5 and @abs(c.x(next) - best_x) <= 20 * dist) break;
        seg_last = next;
        if (on(tags, next)) {
            on_last = next;
            if (on_first < 0) on_first = next;
        }
        if (next == best_point) break;
    }

    if (props & prop.long != 0) {
        // a long zone needs a segment longer than a 25th of the em, on the
        // extreme contour, near the extremum and of the same direction
        const length_threshold = @divTrunc(@as(i32, m.units_per_em), 25);
        const seg_len: i32 = @intCast(@abs(c.x(seg_last) - c.x(seg_first)));
        if (seg_len < length_threshold and seg_last - seg_first + 2 <= best_last - best_first) {
            const height_threshold = @divTrunc(@as(i32, m.units_per_em), 4);
            prev = best_point;
            while (true) {
                prev = c.prev(prev);
                if (c.x(prev) != best_x) break;
                if (prev == best_point) break;
            }
            // degenerate: the glyph is skipped
            if (prev == best_point) return null;
            const left2right = c.x(prev) < c.x(best_point);
            var first = seg_last;
            var lst = first;
            var hit = false;
            var p_first: i32 = 0;
            var p_last: i32 = 0;
            search: while (true) {
                body: {
                    if (!hit) {
                        first = lst;
                        if (on(tags, first)) {
                            p_first = first;
                            p_last = first;
                        } else {
                            p_first = -1;
                            p_last = -1;
                        }
                        hit = true;
                    }
                    lst = c.next(lst);
                    if (@abs(best_y - c.y(first)) > height_threshold) {
                        hit = false;
                        break :body;
                    }
                    dist = @intCast(@abs(c.y(lst) - c.y(first)));
                    if (dist > 5 and @abs(c.x(lst) - c.x(first)) <= 20 * dist) {
                        hit = false;
                        break :body;
                    }
                    if (on(tags, lst)) {
                        p_last = lst;
                        if (p_first < 0) p_first = lst;
                    }
                    const l2r = c.x(first) < c.x(lst);
                    const d: i32 = @intCast(@abs(c.x(lst) - c.x(first)));
                    if (l2r == left2right and d >= length_threshold) {
                        // found: run on to the segment's end
                        while (true) {
                            lst = c.next(lst);
                            const dy: i32 = @intCast(@abs(c.y(lst) - c.y(first)));
                            // FreeType tests `next` and `dist` here, not `last` and `d`
                            if (dy > 5 and @abs(c.x(next) - c.x(first)) <= 20 * dist) {
                                lst = c.prev(lst);
                                break;
                            }
                            p_last = lst;
                            if (on(tags, lst)) {
                                p_last = lst;
                                if (p_first < 0) p_first = lst;
                            }
                            if (lst == seg_first) break;
                        }
                        best_y = c.y(first);
                        seg_first = first;
                        seg_last = lst;
                        on_first = p_first;
                        on_last = p_last;
                        break :search;
                    }
                }
                if (lst == seg_first) break;
            }
        }
    }

    // flat if its on points span more than the threshold, else round if
    // either end is off the curve
    const round = if (on_first >= 0 and on_last >= 0 and @abs(c.x(on_last) - c.x(on_first)) > flat_threshold)
        false
    else
        !on(tags, seg_first) or !on(tags, seg_last);
    // a neutral zone takes flat segments only
    if (round and props & prop.neutral != 0) return null;
    return .{ .y = best_y, .round = round };
}
