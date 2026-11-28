//! Loading a TrueType glyph scaled and hinted as FreeType 2.14.3 loads it
//! (`src/truetype/ttgload.c`: `load_truetype_glyph`, `tt_get_metrics`,
//! `TT_Get_VMetrics`, `tt_loader_set_pp`, `TT_Process_Simple_Glyph`,
//! `TT_Process_Composite_Component`, `TT_Process_Composite_Glyph`,
//! `TT_Hint_Glyph`, `compute_glyph_metrics`; `src/base/ftobjs.c`'s grid
//! fitting of the advance): the outline in 26.6 pixels with its origin at
//! the left phantom point, and the advance; a variable font's at the
//! instance's blend (`../../vary/`).

const std = @import("std");
const fixed = @import("../../fixed.zig");
const trig = @import("../../trig.zig");
const glyf = @import("glyf.zig");
const zone = @import("../zone.zig");
const Instance = @import("../Instance.zig");
const Engine = @import("../Engine.zig");
const vary = @import("../../vary/__root__.zig");
const GlyphIndex = @import("../../../sfnt/types.zig").GlyphIndex;

const Allocator = std.mem.Allocator;
const Point = zone.Point;
const Loader = @This();

pub const Error = error{ InvalidGlyphIndex, InvalidComposite, InvalidOutline, TooManyHints } || vary.Gvar.Error || Engine.Error || Allocator.Error;

pub const Glyph = @import("../../Glyph.zig");

allocator: Allocator,
inst: *Instance,
/// The load hints, as `FT_LOAD_NO_HINTING` is off.
hinted: bool,
/// The programs run: `hinted`, unless `prep` turned hinting off.
run: bool,
points: std.ArrayList(Point) = .empty,
tags: std.ArrayList(u8) = .empty,
contours: std.ArrayList(u16) = .empty,
/// Phantom points: left side, advance, top, bottom.
pp: [4]Point = @splat(.{ 0, 0 }),
x_max: i16 = 0,
y_max: i16 = 0,
x_min: i16 = 0,
left_bearing: i16 = 0,
advance: u16 = 0,
top_bearing: i16 = 0,
vadvance: u16 = 0,

/// Glyph `gid` at `inst`'s size, hinted or only scaled. A failed `fpgm`
/// or `prep` fails a hinted load; a failed glyph program keeps what it did
/// unless the instance is pedantic, as in FreeType.
pub fn load(allocator: Allocator, inst: *Instance, gid: u16, hinted: bool) Error!Glyph {
    if (hinted) if (inst.failed) |err| return err;
    var l: Loader = .{ .allocator = allocator, .inst = inst, .hinted = hinted, .run = hinted and !inst.hinting_off };
    defer l.points.deinit(allocator);
    defer l.tags.deinit(allocator);
    defer l.contours.deinit(allocator);
    inst.engine.backward_compatibility = if (hinted) inst.backward_compatibility else 0;
    try l.glyph(gid, 0);
    const shift = l.pp[0][0];
    if (shift != 0) for (l.points.items) |*p| {
        p[0] -%= shift;
    };
    var advance = l.pp[1][0] -% l.pp[0][0];
    if (hinted) {
        if (l.run) if (inst.widths) |w| {
            advance = @as(i32, w[gid]) * 64;
        };
        advance = fixed.round(advance);
    }
    const tags = try l.tags.toOwnedSlice(allocator);
    errdefer allocator.free(tags);
    for (tags) |*t| t.* &= glyf.on_curve;
    const contours = try l.contours.toOwnedSlice(allocator);
    errdefer allocator.free(contours);
    return .{ .points = try l.points.toOwnedSlice(allocator), .tags = tags, .contours = contours, .advance = advance };
}

fn sfntData(l: *const Loader, gid: u16) []const u8 {
    const s = l.inst.sfnt;
    const range = s.loca.?.getGlyphRange(@fromBackingInt(gid)) orelse return &.{};
    const data = s.glyf.?.table_data;
    if (range.end > data.len) return &.{};
    return data[range.start..range.end];
}

/// `load_truetype_glyph`.
fn glyph(l: *Loader, gid: u16, depth: u32) Error!void {
    if (depth > 100) return error.InvalidComposite;
    const s = l.inst.sfnt;
    if (gid >= s.maxp.?.numGlyphs()) return error.InvalidGlyphIndex;
    const data = l.sfntData(gid);
    var n_contours: i16 = 0;
    l.x_min = 0;
    l.x_max = 0;
    l.y_max = 0;
    if (data.len > 0) {
        const h = try glyf.header(data);
        n_contours = h.n_contours;
        if (n_contours != 0) {
            l.x_min = h.x_min;
            l.x_max = h.x_max;
            l.y_max = h.y_max;
        }
    }
    l.metrics(gid);
    l.setPhantoms();
    if (data.len == 0 or n_contours == 0) {
        if (l.inst.options.blend != null) {
            var pts = l.pp;
            var unrounded: [4]Point = undefined;
            try l.applyDeltas(gid, &pts, &.{}, &unrounded);
        }
        return l.scalePhantoms();
    }
    if (n_contours > 0) return l.simple(gid, data, @intCast(n_contours));
    return l.composite(gid, data, depth);
}

/// The blend's `HVAR`, if FreeType uses one.
fn hvar(l: *const Loader) ?*const vary.Hvar {
    const b = l.inst.options.blend orelse return null;
    return if (b.hvar) |*h| h else null;
}

/// `TT_Vary_Apply_Glyph_Deltas` on `points` (font units, the last four
/// the phantoms), which then give the phantoms but with `HVAR`.
fn applyDeltas(l: *Loader, gid: u16, points: []Point, ends: []const u16, unrounded: []Point) Error!void {
    const b = l.inst.options.blend.?;
    const g = if (b.gvar) |*g| g else {
        for (unrounded, points) |*u, p| u.* = .{ p[0] * 64, p[1] * 64 };
        return;
    };
    const has_hvar = l.hvar() != null;
    try g.apply(l.allocator, gid, b.coords, points, ends, has_hvar, unrounded);
    const n = points.len;
    if (!has_hvar) l.pp[0..2].* = points[n - 4 ..][0..2].*;
    l.pp[2..4].* = points[n - 2 ..][0..2].*;
}

/// `tt_get_metrics`: `hmtx`, then the vertical metrics from `OS/2`'s typo
/// values or `hhea` (no font here has `vmtx`).
fn metrics(l: *Loader, gid: u16) void {
    const s = l.inst.sfnt;
    if (s.hmtx) |h| {
        const m = h.getMetrics(@fromBackingInt(gid));
        l.left_bearing = m.left_side_bearing;
        l.advance = m.advance_width;
    }
    // `HVAR` varies the advance, not the bearing (`tt_face_get_metrics`)
    if (l.hvar()) |h| l.advance = @truncate(@as(u32, @bitCast(@as(i32, l.advance) +% h.advanceDelta(gid, l.inst.options.blend.?.coords))));
    const asc: i32, const desc: i32 = if (s.os2) |o| .{ o.sTypoAscender(), o.sTypoDescender() } else if (s.hhea) |h| .{ h.ascent(), h.descent() } else .{ 0, 0 };
    l.top_bearing = @truncate(asc - l.y_max);
    l.vadvance = @truncate(@abs(asc - desc));
}

/// `tt_loader_set_pp`; v40 centers the vertical phantoms but for mono
/// and LCD targets.
fn setPhantoms(l: *Loader) void {
    const x0: i32 = @as(i32, l.x_min) - l.left_bearing;
    l.pp[0] = .{ x0, 0 };
    l.pp[1] = .{ x0 + l.advance, 0 };
    const top: i32 = @as(i32, l.y_max) + l.top_bearing;
    l.pp[2] = .{ 0, top };
    l.pp[3] = .{ 0, top - l.vadvance };
    const o = l.inst.options;
    if (l.hinted and o.version == .v40 and o.target != .mono and o.target != .lcd and o.target != .lcd_v) {
        const half: i32 = l.advance / 2;
        l.pp[2][0] = half;
        l.pp[3][0] = half;
    }
}

fn scalePhantoms(l: *Loader) void {
    const sc = l.inst.scale;
    l.pp[0][0] = fixed.mul(l.pp[0][0], sc);
    l.pp[1][0] = fixed.mul(l.pp[1][0], sc);
    inline for (2..4) |i| {
        l.pp[i][0] = fixed.mul(l.pp[i][0], sc);
        l.pp[i][1] = fixed.mul(l.pp[i][1], sc);
    }
}

/// `TT_Process_Simple_Glyph`: the points and phantoms varied and scaled,
/// hinted on a zone of their own, then added to the outline.
fn simple(l: *Loader, gid: u16, data: []const u8, n_contours: usize) Error!void {
    var g = try glyf.simple(l.allocator, data, n_contours);
    defer g.deinit(l.allocator);
    const n = g.points.len;
    const z = try l.zoneBuffers(n + 4);
    defer l.freeZone(z);
    @memcpy(z.cur[0..n], g.points);
    @memcpy(z.cur[n..], &l.pp);
    @memcpy(z.tags[0..n], g.tags);
    @memset(z.tags[n..], 0);
    const unrounded: ?[]Point = if (l.inst.options.blend != null) try l.allocator.alloc(Point, n + 4) else null;
    defer if (unrounded) |u| l.allocator.free(u);
    if (unrounded) |u| try l.applyDeltas(gid, z.cur, g.contours, u);
    // the hinter's font units are the varied points, rounded
    @memcpy(z.orus, z.cur);
    const sc = l.inst.scale;
    if (unrounded) |u| {
        // scaled from the unrounded 26.6 points; FT_LOAD_NO_SCALE keeps the rounded ones
        if (!l.inst.unscaled) for (z.cur, u) |*p, q| {
            p.* = .{ (fixed.mul(q[0], sc) +% 32) >> 6, (fixed.mul(q[1], sc) +% 32) >> 6 };
        };
    } else for (z.cur) |*p| {
        p.* = .{ fixed.mul(p[0], sc), fixed.mul(p[1], sc) };
    }
    if (l.hvar() != null and l.hinted) {
        // the `HVAR` advance, scaled as is
        l.pp[0][0] = fixed.mul(l.pp[0][0], sc);
        l.pp[1][0] = fixed.mul(l.pp[1][0], sc);
        l.pp[2..4].* = z.cur[n + 2 ..][0..2].*;
    } else l.pp = z.cur[n..][0..4].*;
    if (l.run) try l.hint(z, g.contours, g.instructions, false);
    const base: u16 = @intCast(l.points.items.len);
    try l.points.appendSlice(l.allocator, z.cur[0..n]);
    try l.tags.appendSlice(l.allocator, z.tags[0..n]);
    for (g.contours) |c| try l.contours.append(l.allocator, base + c);
}

const ZoneBuffers = struct { cur: []Point, org: []Point, orus: []Point, tags: []u8 };

fn zoneBuffers(l: *Loader, n: usize) Allocator.Error!ZoneBuffers {
    const pts = try l.allocator.alloc(Point, 3 * n);
    errdefer l.allocator.free(pts);
    return .{ .cur = pts[0..n], .org = pts[n .. 2 * n], .orus = pts[2 * n ..], .tags = try l.allocator.alloc(u8, n) };
}

fn freeZone(l: *Loader, z: ZoneBuffers) void {
    l.allocator.free(z.tags);
    l.allocator.free(z.cur.ptr[0 .. 3 * z.cur.len]);
}

/// `TT_Hint_Glyph`: the phantoms rounded, the program run, the phantoms
/// kept unless backward compatibility keeps the old ones.
fn hint(l: *Loader, z: ZoneBuffers, contours: []const u16, code: []const u8, composite_glyph: bool) Error!void {
    const n = z.cur.len;
    if (code.len > 0) @memcpy(z.org, z.cur);
    if (composite_glyph) @memcpy(z.orus, z.cur);
    z.cur[n - 4][0] = fixed.round(z.cur[n - 4][0]);
    z.cur[n - 3][0] = fixed.round(z.cur[n - 3][0]);
    z.cur[n - 2][1] = fixed.round(z.cur[n - 2][1]);
    z.cur[n - 1][1] = fixed.round(z.cur[n - 1][1]);
    const e = &l.inst.engine;
    if (code.len > 0) {
        e.zones[1] = .{ .orus = z.orus, .org = z.org, .cur = z.cur, .tags = z.tags, .contours = contours };
        e.code[2] = code;
        e.composite = composite_glyph;
        e.gs = l.inst.gs;
        e.run(.glyph) catch |err| if (e.pedantic) return err;
        e.zones[1] = .{};
    }
    if (e.backward_compatibility != 0) return;
    l.pp = z.cur[n - 4 ..][0..4].*;
}

/// The composite's components, each loaded and placed in turn, then its
/// own program over all their points.
fn composite(l: *Loader, gid: u16, data: []const u8, depth: u32) Error!void {
    const start_point = l.points.items.len;
    const start_contour = l.contours.items.len;
    const c = try glyf.composite(l.allocator, data, l.inst.sfnt.maxp.?.numGlyphs());
    defer l.allocator.free(c.components);
    if (l.inst.options.blend != null) try l.varyComponents(gid, c.components);
    l.scalePhantoms();
    var num_points = start_point;
    for (c.components) |comp| {
        const pp = l.pp;
        const vadvance = l.vadvance;
        const num_base = l.points.items.len;
        try l.glyph(comp.glyph, depth + 1);
        if (comp.flags & glyf.use_my_metrics == 0) {
            l.pp = pp;
            l.vadvance = vadvance;
        }
        num_points = l.points.items.len;
        if (num_points == num_base) continue;
        try l.place(comp, start_point, num_base);
    }
    const last = c.components[c.components.len - 1];
    if (l.run and last.flags & glyf.we_have_instructions != 0 and num_points > start_point)
        try l.compositeProgram(data, c.instructions_at, start_point, start_contour);
}

/// `load_truetype_glyph`'s composite deltas: one point per component (its
/// offset; a contour of its own), then the phantoms; offsets that are
/// points to match stay as they are.
fn varyComponents(l: *Loader, gid: u16, components: []glyf.Component) Error!void {
    const k = components.len;
    const pts = try l.allocator.alloc(Point, 2 * (k + 4));
    defer l.allocator.free(pts);
    const ends = try l.allocator.alloc(u16, k);
    defer l.allocator.free(ends);
    for (components, pts[0..k], ends, 0..) |comp, *p, *e, i| {
        p.* = .{ comp.arg1, comp.arg2 };
        e.* = @intCast(i);
    }
    pts[k..][0..4].* = l.pp;
    try l.applyDeltas(gid, pts[0 .. k + 4], ends, pts[k + 4 ..]);
    for (components, pts[0..k]) |*comp, p| if (comp.flags & glyf.args_are_xy_values != 0) {
        comp.arg1 = @as(i16, @truncate(p[0]));
        comp.arg2 = @as(i16, @truncate(p[1]));
    };
}

/// `TT_Process_Composite_Component`: the transform, then the offset, by
/// points to match or scaled (and with ROUND_XY_TO_GRID rounded) values.
fn place(l: *Loader, c: glyf.Component, start_point: usize, num_base: usize) Error!void {
    const pts = l.points.items;
    const cur = pts[num_base..];
    const have_scale = c.flags & (glyf.we_have_a_scale | glyf.we_have_an_xy_scale | glyf.we_have_a_2x2) != 0;
    if (have_scale) for (cur) |*p| {
        p.* = .{ fixed.mul(p[0], c.xx) +% fixed.mul(p[1], c.xy), fixed.mul(p[0], c.yx) +% fixed.mul(p[1], c.yy) };
    };
    var x: i32 = undefined;
    var y: i32 = undefined;
    if (c.flags & glyf.args_are_xy_values == 0) {
        const k = @as(usize, @intCast(c.arg1)) + start_point;
        const m = @as(usize, @intCast(c.arg2)) + num_base;
        if (k >= num_base or m >= pts.len) return error.InvalidComposite;
        x = pts[k][0] -% pts[m][0];
        y = pts[k][1] -% pts[m][1];
    } else {
        x = c.arg1;
        y = c.arg2;
        if (x == 0 and y == 0) return;
        // FreeType's default: offsets unscaled unless the flag asks
        if (have_scale and c.flags & glyf.scaled_component_offset != 0) {
            x = fixed.mul(x, trig.vectorLength(c.xx, c.xy));
            y = fixed.mul(y, trig.vectorLength(c.yy, c.yx));
        }
        const sc = l.inst.scale;
        x = fixed.mul(x, sc);
        y = fixed.mul(y, sc);
        if (c.flags & glyf.round_xy_to_grid != 0 and l.run) {
            if (l.inst.engine.backward_compatibility == 0) x = fixed.round(x);
            y = fixed.round(y);
        }
    }
    if (x != 0 or y != 0) for (cur) |*p| {
        p.* = .{ p[0] +% x, p[1] +% y };
    };
}

/// `TT_Process_Composite_Glyph`: the composite's program over its
/// components' points, untouched, and the phantoms.
fn compositeProgram(l: *Loader, data: []const u8, at: usize, start_point: usize, start_contour: usize) Error!void {
    if (at + 2 > data.len) return error.InvalidOutline;
    const n_ins = std.mem.readInt(u16, data[at..][0..2], .big);
    if (n_ins == 0) return;
    if (n_ins > data.len or at + 2 + n_ins > data.len) return error.TooManyHints;
    const code = data[at + 2 ..][0..n_ins];
    const n = l.points.items.len - start_point;
    const z = try l.zoneBuffers(n + 4);
    defer l.freeZone(z);
    @memcpy(z.cur[0..n], l.points.items[start_point..]);
    @memcpy(z.cur[n..], &l.pp);
    for (z.tags[0..n], l.tags.items[start_point..]) |*t, s| t.* = s & ~(zone.touched_x | zone.touched_y);
    @memset(z.tags[n..], 0);
    const contours = try l.allocator.alloc(u16, l.contours.items.len - start_contour);
    defer l.allocator.free(contours);
    for (contours, l.contours.items[start_contour..]) |*d, s| d.* = s - @as(u16, @intCast(start_point));
    try l.hint(z, contours, code, true);
    @memcpy(l.points.items[start_point..], z.cur[0..n]);
    @memcpy(l.tags.items[start_point..], z.tags[0..n]);
}
