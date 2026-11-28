//! Grid fitting a dimension's edges: blue zones first, then stems in order,
//! then serifs and lone edges. FreeType 2.14.3 `src/autofit/aflatin.c`
//! (`af_latin_snap_width`, `af_latin_compute_stem_width`,
//! `af_latin_align_linked_edge`, `af_latin_align_serif_edge`,
//! `af_latin_hint_edges`).

const fixed = @import("../../fixed.zig");
const GlyphHints = @import("../GlyphHints.zig");
const Metrics = @import("Metrics.zig");
const Width = @import("../widths.zig").Width;

const Edge = GlyphHints.Edge;
const none = GlyphHints.none;
const ef = GlyphHints.edge_flag;

/// What the render mode asks of the Latin hinter (`AF_LATIN_HINTS_*`).
pub const Mode = struct {
    horz_snap: bool = false,
    vert_snap: bool = false,
    stem_adjust: bool = false,
    mono: bool = false,
    /// The size, for the double-rounding correction.
    ppem: u16,
};

/// `FT_PIX_ROUND`.
pub fn pixRound(x: i32) i32 {
    return (x +% 32) & ~@as(i32, 63);
}

/// `af_latin_snap_width`: to the nearest standard width within its
/// rounding's reach.
pub fn snapWidth(widths: []const Width, width: i32) i32 {
    var best: i32 = 64 + 32 + 2;
    var reference = width;
    for (widths) |w| {
        const dist: i32 = @intCast(@abs(width - w.cur));
        if (dist < best) {
            best = dist;
            reference = w.cur;
        }
    }
    const scaled = pixRound(reference);
    if (width >= reference) {
        if (width < scaled + 48) return reference;
    } else if (width > scaled - 48) return reference;
    return width;
}

/// `af_latin_compute_stem_width`: the width a stem is fitted to.
pub fn stemWidth(m: *const Metrics, mode: Mode, dim: usize, width: i32, base_delta: i32, base_flags: u8, stem_flags: u8) i32 {
    const axis = &m.axis[dim];
    if (!mode.stem_adjust or axis.extra_light) return width;
    var dist = width;
    var sign = false;
    if (dist < 0) {
        dist = -width;
        sign = true;
    }
    const vertical = dim == 1;
    if ((vertical and !mode.vert_snap) or (!vertical and !mode.horz_snap)) {
        // smooth hinting: the width lightly quantized
        done: {
            // serifs keep their width
            if (stem_flags & ef.serif != 0 and vertical and dist < 3 * 64) break :done;
            if (base_flags & ef.round != 0) {
                if (dist < 80) dist = 64;
            } else if (dist < 56) dist = 56;
            if (axis.width_count == 0) break :done;
            // near the standard width: it
            if (@abs(dist - axis.widths[0].cur) < 40) {
                dist = @max(axis.widths[0].cur, 48);
                break :done;
            }
            if (dist < 3 * 64) {
                const delta = dist & 63;
                dist &= -64;
                if (delta < 10) {
                    dist += delta;
                } else if (delta < 32) {
                    dist += 10;
                } else if (delta < 54) {
                    dist += 54;
                } else dist += delta;
            } else {
                // a stem end is rounded twice, as its start and its
                // length: at small sizes, the length gives way
                var bdelta: i32 = 0;
                if ((width > 0 and base_delta > 0) or (width < 0 and base_delta < 0)) {
                    const ppem: i32 = mode.ppem;
                    if (ppem < 10) {
                        bdelta = base_delta;
                    } else if (ppem < 30) bdelta = @divTrunc(base_delta * (30 - ppem), 20);
                    bdelta = @intCast(@abs(bdelta));
                }
                dist = (dist - bdelta + 32) & ~@as(i32, 63);
            }
        }
    } else {
        // strong hinting: whole pixels
        const org_dist = dist;
        dist = snapWidth(axis.widthSlice(), dist);
        if (vertical) {
            dist = if (dist >= 64) (dist + 16) & ~@as(i32, 63) else 64;
        } else if (mode.mono) {
            dist = if (dist < 64) 64 else (dist + 32) & ~@as(i32, 63);
        } else if (dist < 48) {
            dist = (dist + 64) >> 1;
        } else if (dist < 128) {
            // whole pixels only for a distortion under a quarter pixel
            dist = (dist + 22) & ~@as(i32, 63);
            if (@abs(dist - org_dist) >= 16) {
                dist = org_dist;
                if (dist < 48) dist = (dist + 64) >> 1;
            }
        } else dist = (dist + 32) & ~@as(i32, 63);
    }
    return if (sign) -dist else dist;
}

/// `af_latin_align_linked_edge`.
fn alignLinked(m: *const Metrics, mode: Mode, dim: usize, edges: []Edge, base: u32, stem: u32) void {
    const dist = edges[stem].opos - edges[base].opos;
    const base_delta = edges[base].pos - edges[base].opos;
    edges[stem].pos = edges[base].pos + stemWidth(m, mode, dim, dist, base_delta, edges[base].flags, edges[stem].flags);
}

/// The stem width's centering offsets: up and down.
fn offsets(cur_len: i32) [2]i32 {
    return if (cur_len <= 64) .{ 32, 32 } else .{ 38, 26 };
}

/// `af_latin_hint_edges`.
pub fn hintEdges(h: *GlyphHints, m: *const Metrics, mode: Mode, dim: usize, top_to_bottom: bool, blues: []Width) void {
    const edges = h.axis[dim].edges.items;
    const ttb = dim == 1 and top_to_bottom;
    var anchor: u32 = none;
    var has_non_stem_edges = false;

    // stems on blue zones first
    if (dim == 1) for (edges, 0..) |*edge, ei| {
        if (edge.flags & ef.done != 0) continue;
        var edge1: u32 = none;
        var edge2 = edge.link;
        // a stem's neutral zone gives way to the other one
        if (edge.blue_edge != none and edge2 != none and edges[edge2].blue_edge != none) {
            if (edges[edge2].flags & ef.neutral != 0) {
                edges[edge2].blue_edge = none;
                edges[edge2].flags &= ~ef.neutral;
            } else if (edge.flags & ef.neutral != 0) {
                edge.blue_edge = none;
                edge.flags &= ~ef.neutral;
            }
        }
        var blue = edge.blue_edge;
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
    };

    // the other stems, in order
    for (edges, 0..) |*edge, ei_usize| {
        const ei: u32 = @intCast(ei_usize);
        if (edge.flags & ef.done != 0) continue;
        const e2 = edge.link;
        if (e2 == none) {
            has_non_stem_edges = true;
            continue;
        }
        const edge2 = &edges[e2];
        if (edge2.blue_edge != none) {
            // FreeType's "should not happen"
            alignLinked(m, mode, dim, edges, e2, ei);
            edge.flags |= ef.done;
            continue;
        }
        if (anchor == none) {
            // the first stem: its center, then its edges
            const org_len = edge2.opos - edge.opos;
            const cur_len = stemWidth(m, mode, dim, org_len, 0, edge.flags, edge2.flags);
            const off = offsets(cur_len);
            if (cur_len < 96) {
                const org_center = edge.opos + (org_len >> 1);
                var cur_pos1 = pixRound(org_center);
                const error1 = @abs(org_center - (cur_pos1 - off[0]));
                const error2 = @abs(org_center - (cur_pos1 + off[1]));
                if (error1 < error2) cur_pos1 -= off[0] else cur_pos1 += off[1];
                edge.pos = cur_pos1 - @divTrunc(cur_len, 2);
                edge2.pos = edge.pos + cur_len;
            } else edge.pos = pixRound(edge.opos);
            anchor = ei;
            edge.flags |= ef.done;
            alignLinked(m, mode, dim, edges, ei, e2);
        } else {
            const a = edges[anchor];
            const org_pos = a.pos + (edge.opos - a.opos);
            const org_len = edge2.opos - edge.opos;
            const org_center = org_pos + (org_len >> 1);
            const cur_len = stemWidth(m, mode, dim, org_len, 0, edge.flags, edge2.flags);
            if (edge2.flags & ef.done != 0) {
                edge.pos = edge2.pos - cur_len;
            } else if (cur_len < 96) {
                const off = offsets(cur_len);
                var cur_pos1 = pixRound(org_center);
                const delta1 = @abs(org_center - (cur_pos1 - off[0]));
                const delta2 = @abs(org_center - (cur_pos1 + off[1]));
                if (delta1 < delta2) cur_pos1 -= off[0] else cur_pos1 += off[1];
                edge.pos = cur_pos1 - @divTrunc(cur_len, 2);
                edge2.pos = cur_pos1 + @divTrunc(cur_len, 2);
            } else {
                const cur_pos1 = pixRound(org_pos);
                const delta1 = @abs(cur_pos1 + (cur_len >> 1) - org_center);
                const cur_pos2 = pixRound(org_pos + org_len) - cur_len;
                const delta2 = @abs(cur_pos2 + (cur_len >> 1) - org_center);
                edge.pos = if (delta1 < delta2) cur_pos1 else cur_pos2;
                edge2.pos = edge.pos + cur_len;
            }
            edge.flags |= ef.done;
            edge2.flags |= ef.done;
            bound(edges, ei, ttb);
        }
    }

    // a lowercase m keeps its symmetry: three equally spaced stems
    const n_edges = edges.len;
    if (dim == 0 and (n_edges == 6 or n_edges == 12)) {
        const idx: [3]usize = if (n_edges == 6) .{ 0, 2, 4 } else .{ 1, 5, 9 };
        const e1 = &edges[idx[0]];
        const e2 = &edges[idx[1]];
        const e3 = &edges[idx[2]];
        const span = @abs((e2.opos - e1.opos) - (e3.opos - e2.opos));
        if (span < 8) {
            const delta = e3.pos - (2 * e2.pos - e1.pos);
            e3.pos -= delta;
            if (e3.link != none) edges[e3.link].pos -= delta;
            // the serifs move with the stem
            if (n_edges == 12) {
                edges[8].pos -= delta;
                edges[11].pos -= delta;
            }
            e3.flags |= ef.done;
            if (e3.link != none) edges[e3.link].flags |= ef.done;
        }
    }

    if (!has_non_stem_edges and anchor != none) return;
    // serifs and single edges
    for (edges, 0..) |*edge, ei_usize| {
        const ei: u32 = @intCast(ei_usize);
        if (edge.flags & ef.done != 0) continue;
        var delta: i32 = 1000;
        if (edge.serif != none) {
            const serif = edges[edge.serif];
            delta = serif.opos - edge.opos;
            var top: usize = edge.serif;
            var bottom: usize = ei;
            if (delta < 0) {
                delta = -delta;
                top = ei;
                bottom = edge.serif;
            }
            if (delta < 64 + 32) {
                // a real serif has no edge between it and its stem that overlaps them
                const lo = @min(extent(h, dim, edges[top], .min), extent(h, dim, edges[bottom], .min));
                const hi = @max(extent(h, dim, edges[top], .max), extent(h, dim, edges[bottom], .max));
                var e = bottom + 1;
                while (e < top) : (e += 1) {
                    const e_min = extent(h, dim, edges[e], .min);
                    const e_max = extent(h, dim, edges[e], .max);
                    if (!((e_min < lo and e_max < lo) or (e_min > hi and e_max > hi))) {
                        delta = 1000;
                        break;
                    }
                }
                if (delta == 1000) continue;
            }
        }
        if (delta < 64 + 16) {
            const base = edges[edge.serif];
            edge.pos = base.pos + (edge.opos - base.opos);
        } else if (anchor == none) {
            edge.pos = pixRound(edge.opos);
            anchor = ei;
        } else {
            var before: isize = @as(isize, ei) - 1;
            while (before >= 0 and edges[@intCast(before)].flags & ef.done == 0) before -= 1;
            var after: usize = ei + 1;
            while (after < edges.len and edges[after].flags & ef.done == 0) after += 1;
            if (before >= 0 and after < edges.len) {
                const b = edges[@intCast(before)];
                const af = edges[after];
                edge.pos = if (af.opos == b.opos) b.pos else b.pos + fixed.mulDiv(edge.opos - b.opos, af.pos - b.pos, af.opos - b.opos);
            } else {
                const a = edges[anchor];
                edge.pos = a.pos + ((edge.opos - a.opos + 16) & ~@as(i32, 31));
            }
        }
        edge.flags |= ef.done;
        bound(edges, ei, ttb);
        // FreeType compares with the previous edge here too
        if (ei + 1 < edges.len and edges[ei + 1].flags & ef.done != 0 and
            (if (ttb) edge.pos < edges[ei + 1].pos else edge.pos > edges[ei + 1].pos))
        {
            if (edge.link != none and ei > 0 and @abs(edges[edge.link].pos - edges[ei - 1].pos) > 16) edge.pos = edges[ei + 1].pos;
        }
    }
}

/// An edge past its predecessor is held back, unless its stem would
/// (almost) vanish.
fn bound(edges: []Edge, ei: u32, ttb: bool) void {
    if (ei == 0) return;
    const edge = &edges[ei];
    const prev = edges[ei - 1];
    if (!(if (ttb) edge.pos > prev.pos else edge.pos < prev.pos)) return;
    if (edge.link != none and @abs(edges[edge.link].pos - prev.pos) > 16) edge.pos = prev.pos;
}

/// The least or greatest `v` of an edge's first and last segments' ends.
fn extent(h: *const GlyphHints, dim: usize, e: Edge, which: enum { min, max }) i32 {
    const pts = h.points.items;
    const s = h.axis[dim].segments.items;
    const vs = [4]i32{ pts[s[e.first].first].v, pts[s[e.first].last].v, pts[s[e.last].first].v, pts[s[e.last].last].v };
    var r = vs[0];
    for (vs[1..]) |v| r = if (which == .min) @min(r, v) else @max(r, v);
    return r;
}
