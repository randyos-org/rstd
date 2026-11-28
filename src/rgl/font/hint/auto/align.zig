//! Moving points after their edges: FreeType 2.14.3
//! `src/autofit/afhints.c` (`af_glyph_hints_align_edge_points`,
//! `af_glyph_hints_align_strong_points`, `af_iup_shift`, `af_iup_interp`,
//! `af_glyph_hints_align_weak_points`).

const fixed = @import("../fixed.zig");
const GlyphHints = @import("GlyphHints.zig");

const flag = GlyphHints.flag;
const none = GlyphHints.none;

fn touchFlag(dim: usize) u16 {
    return if (dim == 0) flag.touch_x else flag.touch_y;
}

/// Each segment's points onto its edge.
pub fn edgePoints(h: *GlyphHints, dim: usize) void {
    const pts = h.points.items;
    const edges = h.axis[dim].edges.items;
    for (h.axis[dim].segments.items) |seg| {
        if (seg.edge == none) continue;
        const pos = edges[seg.edge].pos;
        var p = seg.first;
        while (true) {
            if (dim == 0) pts[p].x = pos else pts[p].y = pos;
            pts[p].flags |= touchFlag(dim);
            if (p == seg.last) break;
            p = pts[p].next;
        }
    }
}

/// Strong points between the edges they lie between, as TrueType's IP.
pub fn strongPoints(h: *GlyphHints, dim: usize) void {
    const edges = h.axis[dim].edges.items;
    if (edges.len == 0) return;
    const touch = touchFlag(dim);
    for (h.points.items) |*p| {
        if (p.flags & touch != 0 or p.flags & flag.weak != 0) continue;
        var u = if (dim == 1) p.fy else p.fx;
        const ou = if (dim == 1) p.oy else p.ox;
        const fu = u;
        store: {
            // before the first edge or after the last: shifted with it
            if (edges[0].fpos - u >= 0) {
                u = edges[0].pos - (edges[0].opos - ou);
                break :store;
            }
            const last = edges[edges.len - 1];
            if (u - last.fpos >= 0) {
                u = last.pos + (ou - last.opos);
                break :store;
            }
            var min: usize = 0;
            var max: usize = edges.len;
            if (max <= 8) {
                var nn: usize = 0;
                while (nn < max and edges[nn].fpos < u) nn += 1;
                if (edges[nn].fpos == u) {
                    u = edges[nn].pos;
                    break :store;
                }
                min = nn;
            } else while (min < max) {
                const mid = (max + min) >> 1;
                const fpos = edges[mid].fpos;
                if (u < fpos) {
                    max = mid;
                } else if (u > fpos) {
                    min = mid + 1;
                } else {
                    u = edges[mid].pos;
                    break :store;
                }
            }
            const before = &edges[min - 1];
            const after = edges[min];
            if (before.scale == 0) before.scale = fixed.div(after.pos - before.pos, after.fpos - before.fpos);
            u = before.pos + fixed.mul(fu - before.fpos, before.scale);
        }
        if (dim == 0) p.x = u else p.y = u;
        p.flags |= touch;
    }
}

/// `af_iup_shift`.
fn iupShift(pts: []GlyphHints.Point, p1: usize, p2: usize, ref: usize) void {
    const delta = pts[ref].u - pts[ref].v;
    if (delta == 0) return;
    var p = p1;
    while (p < ref) : (p += 1) pts[p].u = pts[p].v + delta;
    p = ref + 1;
    while (p <= p2) : (p += 1) pts[p].u = pts[p].v + delta;
}

/// `af_iup_interp`.
fn iupInterp(pts: []GlyphHints.Point, p1: usize, p2: usize, ref1_in: usize, ref2_in: usize) void {
    if (p1 > p2) return;
    var ref1 = ref1_in;
    var ref2 = ref2_in;
    if (pts[ref1].v > pts[ref2].v) {
        ref1 = ref2_in;
        ref2 = ref1_in;
    }
    const v1 = pts[ref1].v;
    const v2 = pts[ref2].v;
    const ua = pts[ref1].u;
    const ub = pts[ref2].u;
    const d1 = ua - v1;
    const d2 = ub - v2;
    const flat = ua == ub or v1 == v2;
    const scale = if (flat) 0 else fixed.div(ub - ua, v2 - v1);
    var p = p1;
    while (p <= p2) : (p += 1) {
        const u = pts[p].v;
        pts[p].u = if (u <= v1) u + d1 else if (u >= v2) u + d2 else if (flat) ua else ua + fixed.mul(u - v1, scale);
    }
}

/// Weak points between their touched neighbors, as TrueType's IUP.
pub fn weakPoints(h: *GlyphHints, dim: usize) void {
    const pts = h.points.items;
    const touch = touchFlag(dim);
    for (pts) |*p| {
        p.u = if (dim == 0) p.x else p.y;
        p.v = if (dim == 0) p.ox else p.oy;
    }
    for (h.contours.items) |c| {
        const first_point: usize = c;
        const end_point: usize = pts[c].prev;
        var point = first_point;
        while (point <= end_point and pts[point].flags & touch == 0) point += 1;
        // no touched point in the contour
        if (point > end_point) continue;
        const first_touched = point;
        var last_touched = point;
        while (true) {
            while (point < end_point and pts[point + 1].flags & touch != 0) point += 1;
            last_touched = point;
            point += 1;
            while (point <= end_point and pts[point].flags & touch == 0) point += 1;
            if (point > end_point) break;
            iupInterp(pts, last_touched + 1, point - 1, last_touched, point);
        }
        if (last_touched == first_touched) {
            iupShift(pts, first_point, end_point, first_touched);
        } else {
            if (last_touched < end_point) iupInterp(pts, last_touched + 1, end_point, last_touched, first_touched);
            // FreeType tests the outline's first point here, not the contour's
            if (first_touched > 0) iupInterp(pts, first_point, first_touched -| 1, last_touched, first_touched);
        }
    }
    for (pts) |*p| {
        if (dim == 0) p.x = p.u else p.y = p.u;
    }
}
