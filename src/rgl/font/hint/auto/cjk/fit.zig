//! Grid fitting a dimension's edges the CJK way: blue zones, then stems
//! kept apart, then serifs and the rest interpolated. FreeType 2.14.3
//! `src/autofit/afcjk.c` (`af_cjk_compute_stem_width`,
//! `af_cjk_align_linked_edge`, `af_hint_normal_stem`, `af_cjk_hint_edges`,
//! `af_cjk_align_edge_points`).

const fixed = @import("../../fixed.zig");
const GlyphHints = @import("../GlyphHints.zig");
const Metrics = @import("Metrics.zig");
const Width = @import("../widths.zig").Width;
const latin_fit = @import("../latin/fit.zig");

const Edge = GlyphHints.Edge;
const none = GlyphHints.none;
const ef = GlyphHints.edge_flag;
const Mode = latin_fit.Mode;
const pixRound = latin_fit.pixRound;

/// The light mode's limits (`AF_LIGHT_MODE_MAX_*`).
const max_horz_gap = 9;
const max_vert_gap = 15;
const max_delta_abs = 14;

/// `af_cjk_compute_stem_width`.
fn stemWidth(m: *const Metrics, mode: Mode, dim: usize, width: i32) i32 {
    const axis = &m.axis[dim];
    if (!mode.stem_adjust) return width;
    var dist = width;
    var sign = false;
    if (dist < 0) {
        dist = -width;
        sign = true;
    }
    const vertical = dim == 1;
    if ((vertical and !mode.vert_snap) or (!vertical and !mode.horz_snap)) {
        done: {
            // the widths' `cur` is never set: a stem near 0 takes 48
            if (axis.width_count > 0 and @abs(dist - axis.widths[0].cur) < 40) {
                dist = @max(axis.widths[0].cur, 48);
                break :done;
            }
            if (dist < 54) {
                dist += @divTrunc(54 - dist, 2);
            } else if (dist < 3 * 64) {
                const delta = dist & 63;
                dist &= -64;
                if (delta < 10) {
                    dist += delta;
                } else if (delta < 22) {
                    dist += 10;
                } else if (delta < 42) {
                    dist += delta;
                } else if (delta < 54) {
                    dist += 54;
                } else dist += delta;
            }
        }
    } else {
        dist = latin_fit.snapWidth(axis.widths[0..axis.width_count], dist);
        if (vertical) {
            dist = if (dist >= 64) (dist + 16) & ~@as(i32, 63) else 64;
        } else if (mode.mono) {
            dist = if (dist < 64) 64 else (dist + 32) & ~@as(i32, 63);
        } else if (dist < 48) {
            dist = (dist + 64) >> 1;
        } else if (dist < 128) {
            dist = (dist + 22) & ~@as(i32, 63);
        } else dist = (dist + 32) & ~@as(i32, 63);
    }
    return if (sign) -dist else dist;
}

fn alignLinked(m: *const Metrics, mode: Mode, dim: usize, edges: []Edge, base: u32, stem: u32) void {
    const dist = edges[stem].opos - edges[base].opos;
    edges[stem].pos = edges[base].pos + stemWidth(m, mode, dim, dist);
}

/// `af_hint_normal_stem`: the stem centered, then nudged so an edge lands
/// on the grid without shrinking the stem past a threshold; the nudge.
fn hintNormalStem(m: *const Metrics, mode: Mode, edges: []Edge, e1: u32, e2: u32, anchor: i32, dim: usize) i32 {
    const edge = &edges[e1];
    const edge2 = &edges[e2];
    var threshold: i32 = 64;
    if (!mode.stem_adjust) {
        const both_round = edge.flags & ef.round != 0 and edge2.flags & ef.round != 0;
        const gap: i32 = if (dim == 1) max_horz_gap else max_vert_gap;
        threshold = if (both_round) 64 - gap else 64 - @divTrunc(gap, 3);
    }
    const org_len = edge2.opos - edge.opos;
    const cur_len = stemWidth(m, mode, dim, org_len);
    const org_center = @divTrunc(edge.opos + edge2.opos, 2) + anchor;
    var cur_pos1 = org_center - @divTrunc(cur_len, 2);
    const cur_pos2 = cur_pos1 + cur_len;
    var d_off1 = cur_pos1 - (cur_pos1 & ~@as(i32, 63));
    var d_off2 = cur_pos2 - (cur_pos2 & ~@as(i32, 63));
    var u_off1 = 64 - d_off1;
    var u_off2 = 64 - d_off2;
    var delta: i32 = 0;
    exit: {
        if (d_off1 == 0 or d_off2 == 0) break :exit;
        if (cur_len <= threshold) {
            if (d_off2 < cur_len) delta = if (u_off1 <= d_off2) u_off1 else -d_off2;
            break :exit;
        }
        if (threshold < 64 and (d_off1 >= threshold or u_off1 >= threshold or d_off2 >= threshold or u_off2 >= threshold)) break :exit;
        var offset = cur_len & 63;
        if (offset < 32) {
            if (u_off1 <= offset or d_off2 <= offset) break :exit;
        } else offset = 64 - threshold;
        d_off1 = threshold - u_off1;
        u_off1 = u_off1 - offset;
        u_off2 = threshold - d_off2;
        d_off2 = d_off2 - offset;
        if (d_off1 <= u_off1) u_off1 = -d_off1;
        if (d_off2 <= u_off2) u_off2 = -d_off2;
        delta = if (@abs(u_off1) <= @abs(u_off2)) u_off1 else u_off2;
    }
    if (!mode.stem_adjust) delta = @max(-max_delta_abs, @min(delta, max_delta_abs));
    cur_pos1 += delta;
    if (edge.opos < edge2.opos) {
        edge.pos = cur_pos1;
        edge2.pos = cur_pos1 + cur_len;
    } else {
        edge.pos = cur_pos1 + cur_len;
        edge2.pos = cur_pos1;
    }
    return delta;
}

/// `af_cjk_hint_edges`. `blues` holds each zone's reference and undershoot.
pub fn hintEdges(h: *GlyphHints, m: *const Metrics, mode: Mode, dim: usize, blues: []const Width) void {
    const edges = h.axis[dim].edges.items;
    var anchor: u32 = none;
    var delta: i32 = 0;
    var skipped: i32 = 0;
    var has_last_stem = false;
    var last_stem_pos: i32 = 0;

    for (edges, 0..) |*edge, ei| {
        if (edge.flags & ef.done != 0) continue;
        var blue = edge.blue_edge;
        var edge1: u32 = none;
        var edge2 = edge.link;
        if (blue != none) {
            edge1 = @intCast(ei);
        } else if (edge2 != none and edges[edge2].blue_edge != none) {
            blue = edges[edge2].blue_edge;
            edge1 = edge2;
            edge2 = @intCast(ei);
        }
        if (edge1 == none) continue;
        edges[edge1].pos = blues[blue].fit;
        edges[edge1].flags |= ef.done;
        if (edge2 != none and edges[edge2].blue_edge == none) {
            alignLinked(m, mode, dim, edges, edge1, edge2);
            edges[edge2].flags |= ef.done;
        }
        if (anchor == none) anchor = @intCast(ei);
    }

    for (edges, 0..) |*edge, ei_usize| {
        const ei: u32 = @intCast(ei_usize);
        if (edge.flags & ef.done != 0) continue;
        const e2 = edge.link;
        if (e2 == none) {
            skipped += 1;
            continue;
        }
        // stems too near the last one are interpolated at the end instead
        if (has_last_stem and (edge.pos < last_stem_pos + 64 or edges[e2].pos < last_stem_pos + 64)) {
            skipped += 1;
            continue;
        }
        if (edges[e2].blue_edge != none) {
            alignLinked(m, mode, dim, edges, e2, ei);
            edge.flags |= ef.done;
            continue;
        }
        if (e2 < ei) {
            alignLinked(m, mode, dim, edges, e2, ei);
            edge.flags |= ef.done;
            has_last_stem = true;
            last_stem_pos = edge.pos;
            continue;
        }
        if (dim != 1 and anchor == none) {
            delta = hintNormalStem(m, mode, edges, ei, e2, 0, 0);
        } else {
            _ = hintNormalStem(m, mode, edges, ei, e2, delta, dim);
        }
        anchor = ei;
        edge.flags |= ef.done;
        edges[e2].flags |= ef.done;
        has_last_stem = true;
        last_stem_pos = edges[e2].pos;
    }

    // a lowercase m keeps its symmetry
    const n_edges = edges.len;
    if (dim == 0 and (n_edges == 6 or n_edges == 12)) {
        const idx: [3]u32 = if (n_edges == 6) .{ 0, 2, 4 } else .{ 1, 5, 9 };
        const span = @abs((edges[idx[1]].opos - edges[idx[0]].opos) - (edges[idx[2]].opos - edges[idx[1]].opos));
        if (edges[idx[0]].link == idx[0] + 1 and edges[idx[1]].link == idx[1] + 1 and edges[idx[2]].link == idx[2] + 1 and span < 8) {
            const e3 = &edges[idx[2]];
            const d = e3.pos - (2 * edges[idx[1]].pos - edges[idx[0]].pos);
            e3.pos -= d;
            if (e3.link != none) edges[e3.link].pos -= d;
            if (n_edges == 12) {
                edges[8].pos -= d;
                edges[11].pos -= d;
            }
            e3.flags |= ef.done;
            if (e3.link != none) edges[e3.link].flags |= ef.done;
        }
    }

    if (skipped == 0) return;
    for (edges) |*edge| {
        if (edge.flags & ef.done != 0) continue;
        if (edge.serif != none) {
            const base = edges[edge.serif];
            edge.pos = base.pos + (edge.opos - base.opos);
            edge.flags |= ef.done;
            skipped -= 1;
        }
    }
    if (skipped == 0) return;
    for (edges, 0..) |*edge, ei| {
        if (edge.flags & ef.done != 0) continue;
        var before: isize = @as(isize, @intCast(ei)) - 1;
        while (before >= 0 and edges[@intCast(before)].flags & ef.done == 0) before -= 1;
        var after: usize = ei + 1;
        while (after < edges.len and edges[after].flags & ef.done == 0) after += 1;
        if (before < 0 and after >= edges.len) continue;
        if (before < 0) {
            const a = edges[after];
            edge.pos = a.pos + (edge.opos - a.opos);
        } else if (after >= edges.len) {
            const b = edges[@intCast(before)];
            edge.pos = b.pos + (edge.opos - b.opos);
        } else {
            const b = edges[@intCast(before)];
            const a = edges[after];
            edge.pos = if (a.fpos == b.fpos) b.pos else b.pos + fixed.mulDiv(edge.fpos - b.fpos, a.pos - b.pos, a.fpos - b.fpos);
        }
    }
}

/// `af_cjk_align_edge_points`: each edge's points moved with it, or onto
/// it when snapping.
pub fn alignEdgePoints(h: *GlyphHints, mode: Mode, dim: usize) void {
    const pts = h.points.items;
    const segs = h.axis[dim].segments.items;
    const snapping = (dim == 0 and mode.horz_snap) or (dim == 1 and mode.vert_snap);
    const touch = if (dim == 0) GlyphHints.flag.touch_x else GlyphHints.flag.touch_y;
    for (h.axis[dim].edges.items) |edge| {
        var s = edge.first;
        const delta = edge.pos - edge.opos;
        while (true) {
            var p = segs[s].first;
            while (true) {
                const c = if (dim == 0) &pts[p].x else &pts[p].y;
                c.* = if (snapping) edge.pos else c.* + delta;
                pts[p].flags |= touch;
                if (p == segs[s].last) break;
                p = pts[p].next;
            }
            s = segs[s].edge_next;
            if (s == edge.first) break;
        }
    }
}
