//! A glyph's outline as the autohinter analyzes it: FreeType 2.14.3
//! `src/autofit/afhints.c` (`af_direction_compute`, `af_glyph_hints_reload`)
//! and `afhints.h`, with `FT_Outline_Get_Orientation` (`src/base/ftoutln.c`)
//! and `ft_corner_is_flat` (`ftcalc.c`); laid out after fontations
//! `skrifa/src/outline/autohint/outline.rs` (commit 4f115482, MIT). Points,
//! segments and edges refer to each other by index.

const std = @import("std");
const fixed = @import("../fixed.zig");

const Allocator = std.mem.Allocator;
const GlyphHints = @This();

/// A direction: opposite ones sum to 0 (`AF_Direction`).
pub const Dir = i8;
pub const dir_none: Dir = 4;
pub const dir_right: Dir = 1;
pub const dir_left: Dir = -1;
pub const dir_up: Dir = 2;
pub const dir_down: Dir = -2;

/// Point flags (`AF_FLAG_*`).
pub const flag = struct {
    pub const conic: u16 = 1 << 0;
    pub const cubic: u16 = 1 << 1;
    pub const control: u16 = conic | cubic;
    pub const touch_x: u16 = 1 << 2;
    pub const touch_y: u16 = 1 << 3;
    pub const weak: u16 = 1 << 4;
    pub const near: u16 = 1 << 5;
    /// Kept out of segments.
    pub const ignore: u16 = 1 << 6;
};

/// Segment and edge flags (`AF_EDGE_*`).
pub const edge_flag = struct {
    pub const round: u8 = 1 << 0;
    pub const serif: u8 = 1 << 1;
    pub const done: u8 = 1 << 2;
    pub const neutral: u8 = 1 << 3;
    pub const no_blue: u8 = 1 << 4;
};

pub const Dim = enum(u1) { horz = 0, vert = 1 };

/// No point, segment or edge.
pub const none = std.math.maxInt(u32);

pub const Point = struct {
    flags: u16 = 0,
    in_dir: Dir = dir_none,
    out_dir: Dir = dir_none,
    /// Scaled.
    ox: i32 = 0,
    oy: i32 = 0,
    /// Font units, as 16 bits.
    fx: i32 = 0,
    fy: i32 = 0,
    /// Hinted.
    x: i32 = 0,
    y: i32 = 0,
    /// (x, y) or (y, x) by dimension; indices while analyzing the outline.
    u: i32 = 0,
    v: i32 = 0,
    next: u32 = 0,
    prev: u32 = 0,
};

pub const Segment = struct {
    flags: u8 = 0,
    dir: Dir = 0,
    /// Font units, as 16 bits like the rest.
    pos: i32 = 0,
    delta: i32 = 0,
    min_coord: i32 = 0,
    max_coord: i32 = 0,
    height: i32 = 0,
    edge: u32 = none,
    edge_next: u32 = none,
    link: u32 = none,
    serif: u32 = none,
    score: i32 = 32000,
    len: i32 = 0,
    first: u32 = 0,
    last: u32 = 0,
};

pub const Edge = struct {
    fpos: i32 = 0,
    opos: i32 = 0,
    pos: i32 = 0,
    flags: u8 = 0,
    dir: Dir = 0,
    scale: i32 = 0,
    /// The blue zone width it aligns to (`blue_edge`), as an index.
    blue_edge: u32 = none,
    link: u32 = none,
    serif: u32 = none,
    score: i32 = 0,
    first: u32 = none,
    last: u32 = none,
};

pub const Axis = struct {
    segments: std.ArrayList(Segment) = .empty,
    edges: std.ArrayList(Edge) = .empty,
    major_dir: Dir = dir_none,
};

/// An outline in font units, as FreeType's `FT_Outline` holds it.
pub const Outline = struct {
    points: []const [2]i32,
    /// `FT_CURVE_TAG_*`: on 1, conic 0, cubic 2.
    tags: []const u8,
    contours: []const u16,
};

allocator: Allocator,
points: std.ArrayList(Point) = .empty,
/// Each contour's first point.
contours: std.ArrayList(u32) = .empty,
axis: [2]Axis = .{ .{}, .{} },
x_scale: i32 = 1 << 16,
y_scale: i32 = 1 << 16,
x_delta: i32 = 0,
y_delta: i32 = 0,

pub fn deinit(h: *GlyphHints) void {
    h.points.deinit(h.allocator);
    h.contours.deinit(h.allocator);
    for (&h.axis) |*a| {
        a.segments.deinit(h.allocator);
        a.edges.deinit(h.allocator);
    }
}

/// `af_direction_compute`: none unless one arm is 14 times the other.
pub fn direction(dx: i32, dy: i32) Dir {
    var dir: Dir = undefined;
    var ll: i32 = undefined;
    var ss: i32 = undefined;
    if (dy >= dx) {
        if (dy >= -dx) {
            dir = dir_up;
            ll = dy;
            ss = dx;
        } else {
            dir = dir_left;
            ll = -dx;
            ss = dy;
        }
    } else if (dy >= -dx) {
        dir = dir_right;
        ll = dx;
        ss = dy;
    } else {
        dir = dir_down;
        ll = -dy;
        ss = dx;
    }
    if (@as(i64, ll) <= 14 * @as(i64, @abs(ss))) return dir_none;
    return dir;
}

/// `FT_HYPOT`'s approximation.
fn hypot(x: i32, y: i32) i32 {
    const ax: i32 = @intCast(@abs(x));
    const ay: i32 = @intCast(@abs(y));
    return if (ax > ay) ax + ((3 * ay) >> 3) else ay + ((3 * ax) >> 3);
}

/// `ft_corner_is_flat`: one of the vectors much the longer.
fn cornerIsFlat(in_x: i32, in_y: i32, out_x: i32, out_y: i32) bool {
    const d_in = hypot(in_x, in_y);
    const d_out = hypot(out_x, out_y);
    const d_hypot = hypot(in_x + out_x, in_y + out_y);
    return d_in + d_out - d_hypot < (d_hypot >> 4);
}

pub const Orientation = enum { truetype, postscript, none };

/// `FT_Outline_Get_Orientation`: the polygon's signed area, by nonzero winding.
pub fn orientation(o: Outline) Orientation {
    if (o.points.len == 0) return .truetype;
    var min: [2]i32 = o.points[0];
    var max: [2]i32 = o.points[0];
    for (o.points) |p| for (0..2) |k| {
        min[k] = @min(min[k], p[k]);
        max[k] = @max(max[k], p[k]);
    };
    if (min[0] == max[0] or min[1] == max[1]) return .none;
    const lim = 0x1000000;
    if (min[0] < -lim or min[1] < -lim or max[0] > lim or max[1] > lim) return .none;
    const xshift: u5 = @intCast(@max(msb(@as(u32, @abs(max[0])) | @abs(min[0])) - 14, 0));
    const yshift: u5 = @intCast(@max(msb(@intCast(max[1] - min[1])) - 14, 0));
    var area: i64 = 0;
    var first: usize = 0;
    for (o.contours) |end| {
        const last: usize = end;
        var prev: [2]i64 = .{ o.points[last][0] >> xshift, o.points[last][1] >> yshift };
        for (o.points[first .. last + 1]) |p| {
            const cur: [2]i64 = .{ p[0] >> xshift, p[1] >> yshift };
            area +%= (cur[1] - prev[1]) * (cur[0] + prev[0]);
            prev = cur;
        }
        first = last + 1;
    }
    return if (area > 0) .postscript else if (area < 0) .truetype else .none;
}

/// `FT_MSB`.
fn msb(x: u32) i32 {
    return 31 - @as(i32, @clz(x));
}

/// `af_glyph_hints_reload`: the points, their directions, and which are
/// weak, from `o`, scaled by the hints' scales and deltas.
pub fn reload(h: *GlyphHints, o: Outline, units_per_em: u16) Allocator.Error!void {
    h.points.clearRetainingCapacity();
    h.contours.clearRetainingCapacity();
    for (&h.axis) |*a| {
        a.segments.clearRetainingCapacity();
        a.edges.clearRetainingCapacity();
    }
    h.axis[0].major_dir = dir_up;
    h.axis[1].major_dir = dir_left;
    if (orientation(o) == .postscript) {
        h.axis[0].major_dir = dir_down;
        h.axis[1].major_dir = dir_right;
    }
    if (o.points.len == 0) return;
    try h.points.resize(h.allocator, o.points.len);
    try h.contours.resize(h.allocator, o.contours.len);
    const pts = h.points.items;
    @memset(pts, .{});
    const near_limit: i32 = @divTrunc(20 * @as(i32, units_per_em), 2048);

    // coordinates, flags, and the links between points
    var contour: usize = 0;
    var end: usize = o.contours[0];
    var prev: usize = end;
    for (pts, o.points, o.tags, 0..) |*p, v, tag, i| {
        // field by field, as FreeType: the last point's flags are reset
        // after its near flag was set, its link kept
        p.in_dir = dir_none;
        p.out_dir = dir_none;
        p.fx = @as(i16, @truncate(v[0]));
        p.fy = @as(i16, @truncate(v[1]));
        p.flags = switch (tag & 3) {
            0 => flag.conic,
            2 => flag.cubic,
            else => 0,
        };
        p.ox = fixed.mul(v[0], h.x_scale) +% h.x_delta;
        p.x = p.ox;
        p.oy = fixed.mul(v[1], h.y_scale) +% h.y_delta;
        p.y = p.oy;
        // the contour's last point, needed before it is reached
        pts[end].fx = @as(i16, @truncate(o.points[end][0]));
        pts[end].fy = @as(i16, @truncate(o.points[end][1]));
        if (@abs(p.fx - pts[prev].fx) + @abs(p.fy - pts[prev].fy) < near_limit) pts[prev].flags |= flag.near;
        p.prev = @intCast(prev);
        pts[prev].next = @intCast(i);
        prev = i;
        if (i == end) {
            contour += 1;
            if (contour < o.contours.len) {
                end = o.contours[contour];
                prev = end;
            }
        }
    }
    var start: u32 = 0;
    for (h.contours.items, o.contours) |*c, e| {
        c.* = start;
        start = @as(u32, e) + 1;
    }

    // in and out directions; runs of near points take their sum's direction
    const near_limit2 = 2 * near_limit - 1;
    for (h.contours.items) |c| {
        var point: usize = c;
        var back: usize = pts[c].prev;
        while (back != c) {
            if (@abs(pts[point].fx - pts[back].fx) + @abs(pts[point].fy - pts[back].fy) >= near_limit2) break;
            point = back;
            back = pts[back].prev;
        }
        const first = point;
        var curr = first;
        // u and v: the next and previous non-near points
        pts[curr].u = @intCast(first);
        pts[first].v = @intCast(curr);
        var out_x: i32 = 0;
        var out_y: i32 = 0;
        var next = first;
        while (true) {
            point = next;
            next = pts[point].next;
            out_x += pts[next].fx - pts[point].fx;
            out_y += pts[next].fy - pts[point].fy;
            if (@abs(out_x) + @abs(out_y) < near_limit) {
                pts[next].flags |= flag.weak;
            } else {
                pts[curr].u = @intCast(next);
                pts[next].v = @intCast(curr);
                const out_dir = direction(out_x, out_y);
                pts[curr].out_dir = out_dir;
                curr = pts[curr].next;
                while (curr != next) : (curr = pts[curr].next) {
                    pts[curr].in_dir = out_dir;
                    pts[curr].out_dir = out_dir;
                }
                pts[next].in_dir = out_dir;
                pts[curr].u = @intCast(first);
                pts[first].v = @intCast(curr);
                out_x = 0;
                out_y = 0;
            }
            if (next == first) break;
        }
    }

    // vectors of no major direction into one quadrant act as one: the point between is weak
    for (pts) |*p| {
        if (p.flags & flag.weak != 0) continue;
        if (p.in_dir == dir_none and p.out_dir == dir_none) {
            const nu: usize = @intCast(p.u);
            const pv: usize = @intCast(p.v);
            const in_x = p.fx - pts[pv].fx;
            const in_y = p.fy - pts[pv].fy;
            const out_x = pts[nu].fx - p.fx;
            const out_y = pts[nu].fy - p.fy;
            if ((in_x ^ out_x) >= 0 and (in_y ^ out_y) >= 0) {
                p.flags |= flag.weak;
                pts[pv].u = @intCast(nu);
                pts[nu].v = @intCast(pv);
            }
        }
    }

    // control points, points inside a straight run, flat corners and spikes are weak
    for (pts) |*p| {
        if (p.flags & flag.weak != 0) continue;
        var weak = false;
        if (p.flags & flag.control != 0) {
            weak = true;
        } else if (p.out_dir == p.in_dir) {
            if (p.out_dir != dir_none) {
                weak = true;
            } else {
                const nu: usize = @intCast(p.u);
                const pv: usize = @intCast(p.v);
                if (cornerIsFlat(p.fx - pts[pv].fx, p.fy - pts[pv].fy, pts[nu].fx - p.fx, pts[nu].fy - p.fy)) {
                    pts[pv].u = @intCast(nu);
                    pts[nu].v = @intCast(pv);
                    weak = true;
                }
            }
        } else if (p.in_dir == -p.out_dir) {
            weak = true;
        }
        if (weak) p.flags |= flag.weak;
    }
}
