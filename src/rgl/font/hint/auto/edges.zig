//! Edges: segments at one position gathered, with their roundness, stem
//! links and serifs; FreeType 2.14.3 `src/autofit/aflatin.c`
//! (`af_latin_hints_compute_edges`, `af_latin_hints_detect_features`) and
//! `afhints.c` (`af_axis_hints_new_edge`); after fontations
//! `skrifa/src/outline/autohint/topo/edges.rs` (commit 4f115482, MIT).

const std = @import("std");
const fixed = @import("../fixed.zig");
const GlyphHints = @import("GlyphHints.zig");
const segments = @import("segments.zig");

const Allocator = std.mem.Allocator;
const Edge = GlyphHints.Edge;
const edge_flag = GlyphHints.edge_flag;
const none = GlyphHints.none;

/// `af_axis_hints_new_edge`: a zeroed edge at `fpos`, its index; edges stay
/// sorted by position, a minor-direction edge before a major one at the
/// same place.
fn newEdge(h: *GlyphHints, axis: *GlyphHints.Axis, fpos: i32, dir: GlyphHints.Dir, top_to_bottom: bool) Allocator.Error!usize {
    try axis.edges.append(h.allocator, .{});
    const edges = axis.edges.items;
    var at = edges.len - 1;
    while (at > 0) : (at -= 1) {
        const before = edges[at - 1].fpos;
        if (if (top_to_bottom) before > fpos else before < fpos) break;
        if (before == fpos and dir == axis.major_dir) break;
        edges[at] = edges[at - 1];
    }
    edges[at] = .{};
    return at;
}

/// The edge within `threshold` of `pos` (and of direction `dir`, unless
/// `dir` is null).
fn findEdge(edges: []const Edge, pos: i32, threshold: i32, dir: ?GlyphHints.Dir) ?usize {
    for (edges, 0..) |e, i| {
        if (@abs(pos - e.fpos) < threshold and (dir == null or e.dir == dir.?)) return i;
    }
    return null;
}

/// Adds segment `s` to edge `e`'s ring.
fn addToEdge(segs: []GlyphHints.Segment, e: *Edge, s: u32) void {
    segs[s].edge_next = e.first;
    segs[e.last].edge_next = s;
    e.last = s;
}

/// `af_latin_hints_compute_edges`. `distance_threshold` is the metrics'
/// edge distance threshold in font units; `top_to_bottom`, the script's
/// vertical hinting order.
pub fn compute(h: *GlyphHints, dim: GlyphHints.Dim, distance_threshold: i32, top_to_bottom: bool) Allocator.Error!void {
    const axis = &h.axis[@backingInt(dim)];
    axis.edges.clearRetainingCapacity();
    const segs = axis.segments.items;
    const scale = if (dim == .horz) h.x_scale else h.y_scale;
    const ttb = dim == .vert and top_to_bottom;
    // segments under a pixel long, or over a pixel wide, are ignored
    const length_threshold: i32 = if (dim == .horz) fixed.div(64, h.y_scale) else 0;
    const width_threshold = fixed.div(32, scale);
    // at most a quarter pixel apart
    const edge_distance = fixed.div(@min(fixed.mul(distance_threshold, scale), 64 / 4), scale);

    for (segs, 0..) |s, i| {
        if (s.height < length_threshold or s.delta > width_threshold or s.dir == GlyphHints.dir_none) continue;
        // serifs under 1.5 pixels too
        if (s.serif != none and 2 * s.height < 3 * length_threshold) continue;
        if (findEdge(axis.edges.items, s.pos, edge_distance, s.dir)) |found| {
            addToEdge(segs, &axis.edges.items[found], @intCast(i));
        } else {
            const at = try newEdge(h, axis, s.pos, s.dir, ttb);
            const opos = fixed.mul(s.pos, scale);
            axis.edges.items[at] = .{ .first = @intCast(i), .last = @intCast(i), .dir = s.dir, .fpos = s.pos, .opos = opos, .pos = opos };
            segs[i].edge_next = @intCast(i);
        }
    }
    // one-point segments of no direction join an edge near them, if any
    for (segs, 0..) |s, i| {
        if (s.dir != GlyphHints.dir_none) continue;
        if (findEdge(axis.edges.items, s.pos, edge_distance, null)) |found| addToEdge(segs, &axis.edges.items[found], @intCast(i));
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
        while (true) {
            const seg = segs[s];
            if (seg.flags & edge_flag.round != 0) round += 1 else straight += 1;
            // a serif link overrides the stem link
            const is_serif = seg.serif != none and segs[seg.serif].edge != none and segs[seg.serif].edge != ei;
            if ((seg.link != none and segs[seg.link].edge != none) or is_serif) {
                const seg2 = if (is_serif) seg.serif else seg.link;
                var edge2 = if (is_serif) e.serif else e.link;
                if (edge2 != none) {
                    // the nearer of the two
                    if (@abs(seg.pos - segs[seg2].pos) < @abs(e.fpos - edges[edge2].fpos)) edge2 = segs[seg2].edge;
                } else {
                    edge2 = segs[seg2].edge;
                }
                if (is_serif) {
                    e.serif = edge2;
                    edges[edge2].flags |= edge_flag.serif;
                } else {
                    e.link = edge2;
                }
            }
            s = seg.edge_next;
            if (s == e.first) break;
        }
        e.flags = 0;
        if (round > 0 and round >= straight) e.flags |= edge_flag.round;
        // a linked edge drops its serif
        if (e.serif != none and e.link != none) e.serif = none;
    }
}

/// `af_latin_hints_detect_features`: segments, their links and edges for
/// `dim`.
pub fn detectFeatures(h: *GlyphHints, dim: GlyphHints.Dim, max_width: i32, distance_threshold: i32, top_to_bottom: bool, units_per_em: u16) Allocator.Error!void {
    try segments.compute(h, dim, units_per_em);
    segments.link(h, dim, max_width, units_per_em);
    try compute(h, dim, distance_threshold, top_to_bottom);
}
