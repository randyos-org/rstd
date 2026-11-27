//! A variable TrueType font at one instance, as FreeType 2.14.3 sets one
//! from design coordinates (`FT_Set_Var_Design_Coordinates`;
//! `src/truetype/ttgxvar.c` `TT_Get_MM_Var`, `TT_Set_Var_Design`,
//! `ft_var_to_normalized`, `ft_var_load_avar`, `tt_set_mm_blend`): each
//! axis normalized to -1..1 in 16.16, then mapped through `avar`, with the
//! `gvar` and `HVAR` data its glyphs vary by.

const std = @import("std");
const fixed = @import("../fixed.zig");
const Sfnt = @import("../../sfnt/Sfnt.zig");
const packed_data = @import("packed.zig");
const ItemStore = @import("ItemStore.zig");
const Gvar = @import("Gvar.zig");
const Hvar = @import("Hvar.zig");

const Allocator = std.mem.Allocator;
const Cursor = packed_data.Cursor;
const Blend = @This();

pub const Error = error{
    /// The font has no valid `fvar`.
    NotVariable,
    /// `gvar` is present but malformed.
    InvalidTable,
} || Allocator.Error;

allocator: Allocator,
/// Normalized, 16.16, one per axis.
coords: []i32,
gvar: ?Gvar,
/// Null without a valid `HVAR`; the phantom points' deltas give the advance then.
hvar: ?Hvar,

/// The instance at `design` (16.16 axis values in `fvar` order; missing
/// ones at their defaults), or null at the default instance, which
/// FreeType loads as a static font.
pub fn init(allocator: Allocator, sfnt: *const Sfnt, design: []const i32) Error!?Blend {
    const fvar = (sfnt.fvar orelse return error.NotVariable).table_data;
    const axes = try readAxes(allocator, fvar);
    defer allocator.free(axes);
    const coords = try allocator.alloc(i32, axes.len);
    errdefer allocator.free(coords);
    normalize(axes, design, coords);
    if (sfnt.avar) |a| try mapAvar(allocator, a.table_data, coords);
    if (std.mem.allEqual(i32, coords, 0)) {
        allocator.free(coords);
        return null;
    }
    var gvar: ?Gvar = if (sfnt.gvar) |g| try Gvar.init(allocator, g.table_data, @intCast(axes.len)) else null;
    errdefer if (gvar) |*g| g.deinit(allocator);
    var hvar: ?Hvar = if (sfnt.hvar) |h| try Hvar.init(allocator, h.table_data, @intCast(axes.len)) else null;
    errdefer if (hvar) |*h| h.deinit(allocator);
    return .{ .allocator = allocator, .coords = coords, .gvar = gvar, .hvar = hvar };
}

pub fn deinit(b: *Blend) void {
    if (b.hvar) |*h| h.deinit(b.allocator);
    if (b.gvar) |*g| g.deinit(b.allocator);
    b.allocator.free(b.coords);
}

/// An axis's range, 16.16.
const Axis = struct { min: i32, def: i32, max: i32 };

/// `fvar`'s axes, as `sfnt_init_face` validates and `TT_Get_MM_Var` reads
/// them: a range out of order collapses to its default.
fn readAxes(allocator: Allocator, fvar: []const u8) Error![]Axis {
    var c: Cursor = .{ .data = fvar };
    const version = c.ulong();
    const offset = c.ushort();
    _ = c.ushort();
    const n = c.ushort();
    const axis_size = c.ushort();
    const n_inst = c.ushort();
    const inst_size = c.ushort();
    if (fvar.len < 20 or version != 0x10000 or axis_size != 20 or n == 0 or n > 0x3FFE or
        !(inst_size == 4 + 4 * @as(u32, n) or inst_size == 6 + 4 * @as(u32, n)) or n_inst > 0x7EFF or
        @as(usize, offset) + 20 * @as(usize, n) + @as(usize, inst_size) * n_inst > fvar.len) return error.NotVariable;
    const axes = try allocator.alloc(Axis, n);
    c.seek(offset);
    for (axes) |*a| {
        _ = c.ulong();
        const min: i32 = @bitCast(c.ulong());
        const def: i32 = @bitCast(c.ulong());
        const max: i32 = @bitCast(c.ulong());
        _ = c.ulong();
        a.* = if (min > def or def > max) .{ .min = def, .def = def, .max = def } else .{ .min = min, .def = def, .max = max };
    }
    return axes;
}

/// `ft_var_to_normalized`'s first stage: -1 at the minimum, 0 at the
/// default, 1 at the maximum, clamped.
fn normalize(axes: []const Axis, design: []const i32, out: []i32) void {
    for (axes, out, 0..) |a, *n, i| {
        const v = if (i < design.len) design[i] else a.def;
        n.* = if (v > a.def)
            (if (v >= a.max) 0x10000 else fixed.div(v -% a.def, a.max -% a.def))
        else if (v < a.def)
            (if (v <= a.min) -0x10000 else fixed.div(v -% a.def, a.def -% a.min))
        else
            0;
    }
}

/// `ft_var_load_avar` and the second stage: each axis through its
/// segment map, then (`avar` 2) the item store's deltas, clamped to -1..1.
/// A table that does not match `fvar` is ignored, as FreeType ignores it.
fn mapAvar(allocator: Allocator, avar: []const u8, coords: []i32) Allocator.Error!void {
    var c: Cursor = .{ .data = avar };
    const version = c.ulong();
    if (version != 0x10000 and version != 0x20000) return;
    if (c.ulong() != coords.len) return;
    // segment maps, each pair F2Dot14, all checked before any is used
    const maps_at = c.at;
    for (0..coords.len) |_| {
        const pairs = c.ushort();
        if (@as(usize, pairs) * 4 > avar.len) return;
        c.seek(c.at + 4 * @as(usize, pairs));
    }
    const after_maps = c.at;
    c.seek(maps_at);
    for (coords) |*v| {
        const pairs = c.ushort();
        const map = c;
        c.seek(c.at + 4 * @as(usize, pairs));
        v.* = mapSegment(map, pairs, v.*);
    }
    if (version < 0x20000) return;
    c.seek(after_maps);
    const map_offset = c.ulong();
    const store_offset = c.ulong();
    if (store_offset == 0) return;
    var store = ItemStore.init(allocator, avar, store_offset, @intCast(coords.len)) catch |e| switch (e) {
        error.OutOfMemory => return error.OutOfMemory,
        error.InvalidTable => return,
    };
    defer store.deinit(allocator);
    var map: ?ItemStore.IndexMap = null;
    if (map_offset != 0) map = ItemStore.IndexMap.init(allocator, avar, map_offset, &store) catch |e| switch (e) {
        error.OutOfMemory => return error.OutOfMemory,
        error.InvalidTable => return,
    };
    defer if (map) |*m| m.deinit(allocator);
    if (store.data.len == 0) return;
    // the deltas see the half-mapped coordinates
    const half = try allocator.dupe(i32, coords);
    defer allocator.free(half);
    for (coords, 0..) |*v, i| {
        const direct: ItemStore.Index = .{ .outer = 0, .inner = @intCast(i) };
        const idx = if (map) |m| m.get(i) orelse direct else direct;
        const d = store.delta(idx, half);
        v.* = std.math.clamp(half[i] +% d *% 4, -0x10000, 0x10000);
    }
}

/// One axis's segment map: the first segment ending past `v`, linear.
fn mapSegment(map: Cursor, pairs: u16, v: i32) i32 {
    var c = map;
    var prev: [2]i32 = .{ @as(i32, c.short()) * 4, @as(i32, c.short()) * 4 };
    for (1..pairs) |_| {
        const cur: [2]i32 = .{ @as(i32, c.short()) * 4, @as(i32, c.short()) * 4 };
        if (v < cur[0]) return fixed.mulDiv(v -% prev[0], cur[1] -% prev[1], cur[0] -% prev[0]) +% prev[1];
        prev = cur;
    }
    return v;
}

/// `ft_var_apply_tuple`: how much a tuple (its peak, and for an
/// intermediate one its start and end, 16.16) applies at `coords`.
pub fn tupleScalar(coords: []const i32, peak: []const i32, intermediate: ?[2][]const i32) i32 {
    var apply: i32 = 0x10000;
    for (coords, peak, 0..) |ncv, p, i| {
        if (p == 0) continue;
        if (ncv == 0) return 0;
        if (p == ncv) continue;
        if (intermediate) |im| {
            const start = im[0][i];
            const end = im[1][i];
            if (ncv <= start or ncv >= end) return 0;
            apply = if (ncv < p) fixed.mulDiv(apply, ncv -% start, p -% start) else fixed.mulDiv(apply, end -% ncv, end -% p);
        } else {
            if (!((p > ncv and ncv > 0) or (p < ncv and ncv < 0))) return 0;
            apply = fixed.mulDiv(apply, ncv, p);
        }
    }
    return apply;
}
