//! The CJK system's segments, links and edges: FreeType 2.14.3
//! `src/autofit/afcjk.c` (`af_cjk_hints_compute_segments`,
//! `af_cjk_hints_link_segments`, `af_cjk_hints_compute_edges`,
//! `af_cjk_hints_detect_features`, `af_cjk_hints_compute_blue_edges`),
//! on the Latin segments.

const std = @import("std");
const fixed = @import("../../fixed.zig");
const GlyphHints = @import("../GlyphHints.zig");
const segments = @import("../segments.zig");
const Metrics = @import("Metrics.zig");

const Allocator = std.mem.Allocator;
const Segment = GlyphHints.Segment;
const Edge = GlyphHints.Edge;
const none = GlyphHints.none;
const ef = GlyphHints.edge_flag;

/// `AF_SEGMENT_DIST`.
fn segDist(a: Segment, b: Segment) i32 {
    return if (a.pos > b.pos) a.pos - b.pos else b.pos - a.pos;
}

/// The Latin segments. FreeType's CJK version then redoes their roundness,
/// but over the segment count read before they were made, which the
/// outline's reload had set to 0: the loop never runs, the Latin roundness
/// stays, and so it does here.
fn computeSegments(h: *GlyphHints, dim: GlyphHints.Dim, units_per_em: u16) Allocator.Error!void {
    try segments.compute(h, dim, units_per_em);
}

/// `af_cjk_hints_link_segments`: the nearest opposite segment of enough
/// overlap, then serifs at the ends of strokes wider there.
fn linkSegments(h: *GlyphHints, dim: GlyphHints.Dim, units_per_em: u16) void {
    const axis = &h.axis[@backingInt(dim)];
    const segs = axis.segments.items;
    const len_threshold = segments.latinConstant(units_per_em, 8);
    const dist_threshold = fixed.div(64 * 3, if (dim == .horz) h.x_scale else h.y_scale);
    for (segs, 0..) |*seg1, ix1| {
        if (seg1.dir != axis.major_dir) continue;
        for (segs, 0..) |*seg2, ix2| {
            if (ix2 == ix1 or seg1.dir + seg2.dir != 0) continue;
            const dist = seg2.pos - seg1.pos;
            if (dist < 0) continue;
            const len = @min(seg1.max_coord, seg2.max_coord) - @max(seg1.min_coord, seg2.min_coord);
            if (len < len_threshold) continue;
            if (dist * 8 < seg1.score * 9 and (dist * 8 < seg1.score * 7 or seg1.len < len)) {
                seg1.score = dist;
                seg1.len = len;
                seg1.link = @intCast(ix2);
            }
            if (dist * 8 < seg2.score * 9 and (dist * 8 < seg2.score * 7 or seg2.len < len)) {
                seg2.score = dist;
                seg2.len = len;
                seg2.link = @intCast(ix1);
            }
        }
    }
    // a stroke wider at an end: the ends are serifs, or the link goes
    for (segs, 0..) |*seg1, ix1| {
        const l1 = seg1.link;
        if (l1 == none or segs[l1].link != ix1 or segs[l1].pos <= seg1.pos) continue;
        if (seg1.score >= dist_threshold) continue;
        for (segs, 0..) |*seg2, ix2| {
            if (seg2.pos > seg1.pos or ix1 == ix2) continue;
            const l2 = seg2.link;
            if (l2 == none or segs[l2].link != ix2 or segs[l2].pos < segs[l1].pos) continue;
            if (seg1.pos == seg2.pos and segs[l1].pos == segs[l2].pos) continue;
            if (seg2.score <= seg1.score or seg1.score * 4 <= seg2.score) continue;
            // seg2 < seg1 < link1 < link2
            if (seg1.len >= seg2.len * 3) {
                for (segs) |*seg| {
                    if (seg.link == ix2) {
                        seg.link = none;
                        seg.serif = l1;
                    } else if (seg.link == l2) {
                        seg.link = none;
                        seg.serif = @intCast(ix1);
                    }
                }
            } else {
                segs[l1].link = none;
                seg1.link = none;
                break;
            }
        }
    }
    for (segs, 0..) |*seg1, ix1| {
        const s2 = seg1.link;
        if (s2 == none) continue;
        if (segs[s2].link != ix1) {
            seg1.link = none;
            if (segs[s2].score < dist_threshold or seg1.score < segs[s2].score * 4) seg1.serif = segs[s2].link;
        }
    }
}

/// `af_cjk_hints_compute_edges`.
fn computeEdges(h: *GlyphHints, m: *const Metrics, dim: GlyphHints.Dim) Allocator.Error!void {
    const d = @backingInt(dim);
    const axis = &h.axis[d];
    axis.edges.clearRetainingCapacity();
    const segs = axis.segments.items;
    const scale = if (dim == .horz) h.x_scale else h.y_scale;
    const laxis_edt = m.axis[d].edge_distance_threshold;
    // at most a quarter pixel apart
    const edt = if (fixed.mul(laxis_edt, scale) > 64 / 4) fixed.div(64 / 4, scale) else laxis_edt;
    for (segs, 0..) |seg, si| {
        var found: ?usize = null;
        var best: i32 = 0xFFFF;
        for (axis.edges.items, 0..) |edge, ee| {
            if (edge.dir != seg.dir) continue;
            const dist: i32 = @intCast(@abs(seg.pos - edge.fpos));
            if (dist >= edt or dist >= best) continue;
            // the edge's linked segments must all stay near this one's link
            if (seg.link != none) {
                var s1 = edge.first;
                var dist2: i32 = 0;
                while (true) {
                    const link1 = segs[s1].link;
                    if (link1 != none) {
                        dist2 = segDist(segs[seg.link], segs[link1]);
                        if (dist2 >= edt) break;
                    }
                    s1 = segs[s1].edge_next;
                    if (s1 == edge.first) break;
                }
                if (dist2 >= edt) continue;
            }
            best = dist;
            found = ee;
        }
        if (found) |f| {
            const e = &axis.edges.items[f];
            segs[si].edge_next = e.first;
            segs[e.last].edge_next = @intCast(si);
            e.last = @intCast(si);
        } else {
            const at = try newEdge(h, axis, seg.pos, seg.dir);
            const opos = fixed.mul(seg.pos, scale);
            axis.edges.items[at] = .{ .first = @intCast(si), .last = @intCast(si), .dir = seg.dir, .fpos = seg.pos, .opos = opos, .pos = opos };
            segs[si].edge_next = @intCast(si);
        }
    }
    const edges = axis.edges.items;
    for (edges, 0..) |e, ei| {
        var s = e.first;
        if (s == none) continue;
        while (true) {
            segs[s].edge = @intCast(ei);
            s = segs[s].edge_next;
            if (s == e.first) break;
        }
    }
    for (edges, 0..) |*e, ei| {
        var round: u32 = 0;
        var straight: u32 = 0;
        var s = e.first;
        if (s != none) while (true) {
            const seg = segs[s];
            if (seg.flags & ef.round != 0) round += 1 else straight += 1;
            const is_serif = seg.serif != none and segs[seg.serif].edge != ei;
            if (seg.link != none or is_serif) {
                const seg2 = if (is_serif) seg.serif else seg.link;
                var edge2 = if (is_serif) e.serif else e.link;
                if (edge2 != none) {
                    if (segDist(seg, segs[seg2]) < @abs(e.fpos - edges[edge2].fpos)) edge2 = segs[seg2].edge;
                } else edge2 = segs[seg2].edge;
                if (is_serif) {
                    e.serif = edge2;
                    edges[edge2].flags |= ef.serif;
                } else e.link = edge2;
            }
            s = seg.edge_next;
            if (s == e.first) break;
        };
        e.flags = 0;
        if (round > 0 and round >= straight) e.flags |= ef.round;
        if (e.serif != none and e.link != none) e.serif = none;
    }
}

/// `af_axis_hints_new_edge`, bottom to top.
fn newEdge(h: *GlyphHints, axis: *GlyphHints.Axis, fpos: i32, dir: GlyphHints.Dir) Allocator.Error!usize {
    try axis.edges.append(h.allocator, .{});
    const edges = axis.edges.items;
    var at = edges.len - 1;
    while (at > 0) : (at -= 1) {
        const before = edges[at - 1].fpos;
        if (before < fpos) break;
        if (before == fpos and dir == axis.major_dir) break;
        edges[at] = edges[at - 1];
    }
    edges[at] = .{};
    return at;
}

pub fn detectFeatures(h: *GlyphHints, m: *const Metrics, dim: GlyphHints.Dim) Allocator.Error!void {
    try computeSegments(h, dim, m.units_per_em);
    linkSegments(h, dim, m.units_per_em);
    try computeEdges(h, m, dim);
}

/// `af_cjk_hints_compute_blue_edges`: each edge's nearest zone edge, the
/// zone's reference or undershoot, whichever is nearer. A blue edge is an
/// index: twice the zone's, plus 1 for its undershoot.
pub fn computeBlueEdges(h: *GlyphHints, m: *const Metrics, dim: GlyphHints.Dim) void {
    const d = @backingInt(dim);
    const axis = &h.axis[d];
    const cjk = &m.axis[d];
    const scale = cjk.scale;
    const best_dist0 = @min(fixed.mul(@divTrunc(@as(i32, m.units_per_em), 40), scale), 64 / 2);
    for (axis.edges.items) |*edge| {
        var best: u32 = none;
        var best_dist = best_dist0;
        for (cjk.blueSlice(), 0..) |blue, bb| {
            if (blue.flags & Metrics.blue_flag.active == 0) continue;
            const is_top_right = blue.flags & Metrics.blue_flag.top != 0;
            const is_major = edge.dir == axis.major_dir;
            if (is_top_right == is_major) continue;
            const shoot_nearer = @abs(edge.fpos - blue.ref.org) > @abs(edge.fpos - blue.shoot.org);
            const cmp = if (shoot_nearer) blue.shoot.org else blue.ref.org;
            const dist = fixed.mul(@intCast(@abs(edge.fpos - cmp)), scale);
            if (dist < best_dist) {
                best_dist = dist;
                best = @intCast(2 * bb + @intFromBool(shoot_nearer));
            }
        }
        if (best != none) edge.blue_edge = best;
    }
}
