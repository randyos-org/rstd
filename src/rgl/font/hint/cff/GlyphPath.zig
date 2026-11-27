//! The charstring's path, hinted and drawn as FreeType 2.14.3's CFF engine
//! draws it (`src/psaux/pshints.c` `cf2_glyphpath_*`, the stem and mask
//! handling of `psintrp.c` `cf2_interpT2CharString`), with stem darkening
//! off: each element is hinted one element late, with the hint map in force
//! when the next one arrives, and a contour closes with the map of its
//! start. The evaluator's sink (`charstring.zig`).

const std = @import("std");
const fixed = @import("../fixed.zig");
const hintmap = @import("hintmap.zig");
const Blues = @import("blues.zig").Blues;
const Builder = @import("Builder.zig");

const Allocator = std.mem.Allocator;
const HintMap = hintmap.HintMap;
const HintMask = hintmap.HintMask;
const StemHint = hintmap.StemHint;
const GlyphPath = @This();
const Point = [2]i32;

blues: *const Blues,
/// The x scale (`innerTransform.a`); y's is the hint maps'.
scale_x: i32,
out: *Builder,
err: ?Allocator.Error = null,
too_many_hints: bool = false,

stems: [hintmap.max_hints]StemHint = undefined,
n_hstems: usize = 0,
n_vstems: usize = 0,
mask: HintMask = .{},
mask_valid: bool = false,
mask_new: bool = false,
initial: HintMap,
first: HintMap,
map: HintMap,

move_pending: bool = true,
path_open: bool = false,
path_closing: bool = false,
queued: bool = false,
prev_cubic: bool = false,
prev: [4]Point = undefined,
current_cs: Point = .{ 0, 0 },
start: Point = .{ 0, 0 },
current_ds: Point = .{ 0, 0 },
offset_start0: Point = .{ 0, 0 },
offset_start1: Point = .{ 0, 0 },

pub fn init(blues: *const Blues, out: *Builder, scale_x: i32, scale_y: i32, hinted: bool) GlyphPath {
    const empty: HintMap = .{ .scale = scale_y, .hinted = hinted };
    return .{ .blues = blues, .scale_x = scale_x, .out = out, .initial = empty, .first = empty, .map = empty };
}

fn hstems(p: *GlyphPath) []StemHint {
    return p.stems[0..p.n_hstems];
}

/// A mask's bit count, false past FreeType's limit (`cf2_hintmask_setCounts`).
fn countOk(p: *GlyphPath) bool {
    const n = p.n_hstems + p.n_vstems;
    if (n > hintmap.max_hints) p.too_many_hints = true;
    return n > 0 and n <= hintmap.max_hints;
}

/// `cf2_hintmap_build` for the glyph's own maps.
fn build(p: *GlyphPath, m: *HintMap) void {
    if (!p.mask_valid) {
        // no mask yet: every hint, from now on fixed
        if (!p.countOk()) return;
        p.mask = .{};
        p.mask_valid = true;
    }
    if (!p.initial.valid and !p.countOk()) return;
    m.build(p.blues, p.hstems(), p.mask, &p.initial, 0);
    p.mask_new = false;
}

pub fn hstem(p: *GlyphPath, min: i32, max: i32) void {
    // hints after the mask is set are ignored
    if (p.mask_valid) return;
    if (p.n_hstems + p.n_vstems >= hintmap.max_hints) {
        p.too_many_hints = true;
        return;
    }
    p.stems[p.n_hstems] = .{ .min = min, .max = max };
    p.n_hstems += 1;
}

pub fn vstem(p: *GlyphPath, min: i32, max: i32) void {
    _ = min;
    _ = max;
    if (!p.mask_valid) p.n_vstems += 1;
}

pub fn hintMask(p: *GlyphPath, bytes: []const u8) void {
    if (!p.countOk()) return;
    p.mask = HintMask.fromBytes(bytes[0..@min(bytes.len, hintmap.max_hints / 8)]);
    p.mask_valid = true;
    p.mask_new = true;
}

/// Places and locks the stems a counter mask names, in a map of its own.
pub fn counterMask(p: *GlyphPath, bytes: []const u8) void {
    if (!p.countOk()) return;
    var m: HintMap = .{ .scale = p.map.scale, .hinted = p.map.hinted };
    m.build(p.blues, p.hstems(), HintMask.fromBytes(bytes[0..@min(bytes.len, hintmap.max_hints / 8)]), &p.initial, 0);
}

fn hintPoint(p: *GlyphPath, m: *const HintMap, pt: Point) Point {
    return .{ fixed.mul(p.scale_x, pt[0]), m.map(pt[1]) };
}

pub fn moveTo(p: *GlyphPath, x: i32, y: i32) void {
    p.closeOpenPath();
    p.current_cs = .{ x, y };
    p.start = .{ x, y };
    p.move_pending = true;
    if (!p.map.valid or p.mask_new) p.build(&p.map);
    p.first = p.map;
}

pub fn lineTo(p: *GlyphPath, x: i32, y: i32) void {
    const new_map = p.mask_new and !p.path_closing;
    // a line of no length in character space, under the same map, draws nothing
    if (p.current_cs[0] == x and p.current_cs[1] == y and !new_map) return;
    var p0 = p.current_cs;
    const p1: Point = .{ x, y };
    p.beginElement(p0, p1);
    if (p.queued) p.pushPrev(&p.map, &p0, p1, false);
    p.queued = true;
    p.prev_cubic = false;
    p.prev[0] = p0;
    p.prev[1] = p1;
    if (new_map) p.build(&p.map);
    p.current_cs = p1;
}

pub fn curveTo(p: *GlyphPath, x1: i32, y1: i32, x2: i32, y2: i32, x3: i32, y3: i32) void {
    var p0 = p.current_cs;
    const p1: Point = .{ x1, y1 };
    p.beginElement(p0, p1);
    if (p.queued) p.pushPrev(&p.map, &p0, p1, false);
    p.queued = true;
    p.prev_cubic = true;
    p.prev = .{ p0, p1, .{ x2, y2 }, .{ x3, y3 } };
    if (p.mask_new) p.build(&p.map);
    p.current_cs = .{ x3, y3 };
}

/// The contour's move, drawn with its first element.
fn beginElement(p: *GlyphPath, p0: Point, p1: Point) void {
    if (!p.move_pending) return;
    // a first contour with no move starts at the origin (`cf2_glyphpath_pushMove`)
    if (!p.map.valid) p.moveTo(p.start[0], p.start[1]);
    const pt = p.hintPoint(&p.map, p0);
    p.out.moveTo();
    p.current_ds = pt;
    p.offset_start0 = p0;
    p.move_pending = false;
    p.path_open = true;
    p.offset_start1 = p1;
}

/// The queued element, drawn (`cf2_glyphpath_pushPrevElem`); with no
/// darkening offset, elements always join.
fn pushPrev(p: *GlyphPath, m: *const HintMap, next_p0: *Point, next_p1: Point, closing: bool) void {
    _ = next_p1;
    if (p.err != null) return;
    const close_map: *const HintMap = if (closing) &p.first else m;
    if (!p.prev_cubic) {
        const pt1 = p.hintPoint(close_map, p.prev[1]);
        if (pt1[0] != p.current_ds[0] or pt1[1] != p.current_ds[1]) {
            p.out.lineTo(p.current_ds, pt1) catch |e| {
                p.err = e;
            };
            p.current_ds = pt1;
        }
    } else {
        const pt1 = p.hintPoint(m, p.prev[1]);
        const pt2 = p.hintPoint(m, p.prev[2]);
        const pt3 = p.hintPoint(m, p.prev[3]);
        p.out.cubeTo(p.current_ds, pt1, pt2, pt3) catch |e| {
            p.err = e;
        };
        p.current_ds = pt3;
    }
    // the join to the next element, of no length unless the maps differ
    const pt = p.hintPoint(close_map, next_p0.*);
    if (pt[0] != p.current_ds[0] or pt[1] != p.current_ds[1]) {
        p.out.lineTo(p.current_ds, pt) catch |e| {
            p.err = e;
        };
        p.current_ds = pt;
    }
}

pub fn closeOpenPath(p: *GlyphPath) void {
    if (!p.path_open) return;
    p.path_closing = true;
    p.lineTo(p.start[0], p.start[1]);
    if (p.queued) p.pushPrev(&p.map, &p.offset_start0, p.offset_start1, true);
    p.move_pending = true;
    p.path_open = false;
    p.path_closing = false;
    p.queued = false;
}

/// The evaluator closes contours at each move; the path does it itself.
pub fn close(p: *GlyphPath) void {
    _ = p;
}

pub fn finish(p: *GlyphPath) void {
    p.closeOpenPath();
    p.out.closeContour();
}
