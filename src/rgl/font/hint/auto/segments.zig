//! Segments: runs of a contour along one dimension's major direction, and
//! the stems and serifs they pair into; FreeType 2.14.3
//! `src/autofit/aflatin.c` (`af_latin_hints_compute_segments`,
//! `af_latin_hints_link_segments`) and `afhints.c`
//! (`af_axis_hints_new_segment`); after fontations
//! `skrifa/src/outline/autohint/topo/segments.rs` (commit 4f115482, MIT).

const std = @import("std");
const GlyphHints = @import("GlyphHints.zig");

const Allocator = std.mem.Allocator;
const Segment = GlyphHints.Segment;
const flag = GlyphHints.flag;
const none = GlyphHints.none;

/// FreeType's cap; a glyph with more gets no segments.
const max_segments = 1000;

/// A Latin constant, given for 2048 units per em (`AF_LATIN_CONSTANT`).
pub fn latinConstant(units_per_em: u16, c: i32) i32 {
    return @intCast(@divTrunc(@as(i64, c) * units_per_em, 2048));
}

/// The running extent of a segment being built.
const Extent = struct {
    min_pos: i32 = 32000,
    max_pos: i32 = -32000,
    min_coord: i32 = 32000,
    max_coord: i32 = -32000,
    min_flags: u16 = 0,
    max_flags: u16 = 0,
    min_on: i32 = 32000,
    max_on: i32 = -32000,

    fn round(e: Extent, flat_threshold: i32) bool {
        return (e.min_flags | e.max_flags) & flag.control != 0 and e.max_on - e.min_on < flat_threshold;
    }

    fn setPos(e: Extent, s: *Segment) void {
        s.pos = short((e.min_pos + e.max_pos) >> 1);
        s.delta = short((e.max_pos - e.min_pos) >> 1);
    }

    fn setCoords(e: Extent, s: *Segment) void {
        s.min_coord = short(e.min_coord);
        s.max_coord = short(e.max_coord);
        s.height = short(s.max_coord - s.min_coord);
    }
};

/// `(FT_Short)`.
fn short(v: i32) i32 {
    return @as(i16, @truncate(v));
}

/// `af_latin_hints_compute_segments`.
pub fn compute(h: *GlyphHints, dim: GlyphHints.Dim, units_per_em: u16) Allocator.Error!void {
    const axis = &h.axis[@backingInt(dim)];
    const pts = h.points.items;
    const segs = &axis.segments;
    segs.clearRetainingCapacity();
    const flat_threshold = @divTrunc(@as(i32, units_per_em), 14);
    const major_dir: GlyphHints.Dir = @intCast(@abs(axis.major_dir));
    var segment_dir = major_dir;
    for (pts) |*p| {
        p.u = if (dim == .horz) p.fx else p.fy;
        p.v = if (dim == .horz) p.fy else p.fx;
    }
    for (h.contours.items) |c| {
        var point: usize = c;
        var last: usize = pts[point].prev;
        var on_edge = false;
        var e: Extent = .{};
        var prev: Extent = e;
        var prev_segment: u32 = none;
        var segment: u32 = none;
        if (pts[point].flags & flag.ignore == 0 and @abs(pts[last].out_dir) == major_dir and @abs(pts[point].out_dir) == major_dir) {
            // already on an edge: find its start
            last = point;
            while (true) {
                point = pts[point].prev;
                if (@abs(pts[point].out_dir) != major_dir) {
                    point = pts[point].next;
                    break;
                }
                if (point == last) break;
            }
        }
        last = point;
        var passed = false;
        while (true) {
            const p = pts[point];
            if (on_edge) {
                e.min_pos = @min(e.min_pos, p.u);
                e.max_pos = @max(e.max_pos, p.u);
                if (p.v < e.min_coord) {
                    e.min_coord = p.v;
                    e.min_flags = p.flags;
                }
                if (p.v > e.max_coord) {
                    e.max_coord = p.v;
                    e.max_flags = p.flags;
                }
                if (p.flags & flag.control == 0) {
                    e.min_on = @min(e.min_on, p.v);
                    e.max_on = @max(e.max_on, p.v);
                }
                if (p.flags & flag.ignore != 0 or p.out_dir != segment_dir or point == last) {
                    const s = &segs.items[segment];
                    if (p.flags & flag.ignore != 0 or prev_segment == none or s.first != segs.items[prev_segment].last) {
                        // leaving an edge: a new segment
                        s.last = @intCast(point);
                        e.setPos(s);
                        if (e.round(flat_threshold)) s.flags |= GlyphHints.edge_flag.round;
                        e.setCoords(s);
                        prev_segment = segment;
                        prev = e;
                    } else {
                        // it starts where the previous one ends: merge them
                        const ps = &segs.items[prev_segment];
                        if (pts[ps.last].in_dir == p.in_dir) {
                            // the same direction, as a zig-zag degenerate outline has
                            e.min_pos = @min(e.min_pos, prev.min_pos);
                            e.max_pos = @max(e.max_pos, prev.max_pos);
                            if (prev.min_coord < e.min_coord) {
                                e.min_coord = prev.min_coord;
                                e.min_flags = prev.min_flags;
                            }
                            if (prev.max_coord > e.max_coord) {
                                e.max_coord = prev.max_coord;
                                e.max_flags = prev.max_flags;
                            }
                            e.min_on = @min(e.min_on, prev.min_on);
                            e.max_on = @max(e.max_on, prev.max_on);
                            ps.last = @intCast(point);
                            e.setPos(ps);
                            if (e.round(flat_threshold)) ps.flags |= GlyphHints.edge_flag.round else ps.flags &= ~GlyphHints.edge_flag.round;
                            e.setCoords(ps);
                        } else if (@abs(prev.max_coord - prev.min_coord) > @abs(e.max_coord - e.min_coord)) {
                            // the longer one's properties: the current one is dropped
                            prev.min_pos = @min(prev.min_pos, e.min_pos);
                            prev.max_pos = @max(prev.max_pos, e.max_pos);
                            ps.last = @intCast(point);
                            prev.setPos(ps);
                        } else {
                            e.min_pos = @min(e.min_pos, prev.min_pos);
                            e.max_pos = @max(e.max_pos, prev.max_pos);
                            s.last = @intCast(point);
                            e.setPos(s);
                            if (e.round(flat_threshold)) s.flags |= GlyphHints.edge_flag.round;
                            e.setCoords(s);
                            ps.* = s.*;
                            prev = e;
                        }
                        segs.items.len -= 1;
                    }
                    on_edge = false;
                    segment = none;
                }
            }
            if (point == last) {
                if (passed) break;
                passed = true;
            }
            if (p.flags & flag.ignore == 0 and !on_edge and (@abs(p.out_dir) == major_dir or point == p.prev)) {
                if (segs.items.len > max_segments) {
                    segs.clearRetainingCapacity();
                    return;
                }
                // a new segment starts
                segment_dir = p.out_dir;
                segment = @intCast(segs.items.len);
                try segs.append(h.allocator, .{ .dir = segment_dir, .first = @intCast(point), .last = @intCast(point) });
                e = .{ .min_pos = p.u, .max_pos = p.u, .min_coord = p.v, .max_coord = p.v, .min_flags = p.flags, .max_flags = p.flags };
                if (p.flags & flag.control == 0) {
                    e.min_on = p.v;
                    e.max_on = p.v;
                }
                on_edge = true;
                if (point == p.prev) {
                    // a contour of one point, of no direction
                    const s = &segs.items[segment];
                    s.pos = short(e.min_pos);
                    if (p.flags & flag.control != 0) s.flags |= GlyphHints.edge_flag.round;
                    s.min_coord = short(p.v);
                    s.max_coord = short(p.v);
                    s.height = 0;
                    on_edge = false;
                    segment = none;
                }
            }
            point = p.next;
        }
    }

    // a little more height where the contour goes on beyond the ends, to tell serifs
    for (segs.items) |*s| {
        const first_v = pts[s.first].v;
        const last_v = pts[s.last].v;
        const before = pts[pts[s.first].prev].v;
        const after = pts[pts[s.last].next].v;
        if (first_v < last_v) {
            if (before < first_v) s.height = short(s.height + ((first_v - before) >> 1));
            if (after > last_v) s.height = short(s.height + ((after - last_v) >> 1));
        } else {
            if (before > first_v) s.height = short(s.height + ((before - first_v) >> 1));
            if (after < last_v) s.height = short(s.height + ((last_v - after) >> 1));
        }
    }
}

/// `af_latin_hints_link_segments`: each major-direction segment pairs with
/// the best opposite one to its right; `max_width`, the widest standard
/// width, 0 if none, weighs the distance.
pub fn link(h: *GlyphHints, dim: GlyphHints.Dim, max_width: i32, units_per_em: u16) void {
    const axis = &h.axis[@backingInt(dim)];
    const segs = axis.segments.items;
    const len_threshold = @max(latinConstant(units_per_em, 8), 1);
    const len_score = latinConstant(units_per_em, 6000);
    const dist_score = 3000;
    for (segs, 0..) |*seg1, ix1| {
        if (seg1.dir != axis.major_dir) continue;
        for (segs, 0..) |*seg2, ix2| {
            const pos1 = seg1.pos;
            const pos2 = seg2.pos;
            if (seg1.dir + seg2.dir != 0 or pos2 <= pos1) continue;
            // how much they overlap
            const len = @min(seg1.max_coord, seg2.max_coord) - @max(seg1.min_coord, seg2.min_coord);
            if (len < len_threshold) continue;
            const dist = pos2 - pos1;
            var dist_demerit: i32 = dist;
            if (max_width != 0) {
                // in multiples of the widest stem, times 1024
                const delta = @divTrunc(dist << 10, max_width) - (1 << 10);
                dist_demerit = if (delta > 10000) 32000 else if (delta > 0) @divTrunc(delta * delta, dist_score) else 0;
            }
            const score = dist_demerit + @divTrunc(len_score, len);
            if (score < seg1.score) {
                seg1.score = score;
                seg1.link = @intCast(ix2);
            }
            if (score < seg2.score) {
                seg2.score = score;
                seg2.link = @intCast(ix1);
            }
        }
    }
    // a one-sided link is a serif of the segment's partner's partner
    for (segs, 0..) |*seg1, ix1| {
        if (seg1.link == none) continue;
        const seg2 = segs[seg1.link];
        if (seg2.link != ix1) {
            seg1.link = none;
            seg1.serif = seg2.link;
        }
    }
}
