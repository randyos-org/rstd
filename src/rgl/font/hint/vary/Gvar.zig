//! `gvar`, read and applied as FreeType 2.14.3 does
//! (`src/truetype/ttgxvar.c` `ft_var_load_gvar`,
//! `TT_Vary_Apply_Glyph_Deltas`): each applying tuple's deltas scaled in
//! 16.16 and summed, the points a tuple leaves out inferred (`iup.zig`).

const std = @import("std");
const fixed = @import("../fixed.zig");
const packed_data = @import("packed.zig");
const tuples = @import("tuples.zig");
const iup = @import("iup.zig");

const Allocator = std.mem.Allocator;
const Cursor = packed_data.Cursor;
const Gvar = @This();

pub const Error = tuples.Error;
pub const Point = iup.Point;

table: []const u8,
/// Each glyph's data, `offsets[g]..offsets[g + 1]` in `table`.
offsets: []u32,
/// The shared tuples' peaks, 16.16, one per axis each.
shared: []i32,

pub fn init(allocator: Allocator, table: []const u8, axis_count: u16) Error!Gvar {
    var c: Cursor = .{ .data = table };
    if (c.left() < 20 or c.ulong() != 0x10000) return error.InvalidTable;
    if (c.ushort() != axis_count) return error.InvalidTable;
    const shared_count = c.ushort();
    const shared_at = c.ulong();
    const glyph_count = c.ushort();
    const long = c.ushort() & 1 != 0;
    const data_at = c.ulong();
    if (@as(usize, shared_count) * axis_count > table.len / 2) return error.InvalidTable;
    const offsets_len = (@as(usize, glyph_count) + 1) * @as(usize, if (long) 4 else 2);
    if (offsets_len > c.left()) return error.InvalidTable;
    const offsets = try allocator.alloc(u32, @as(usize, glyph_count) + 1);
    errdefer allocator.free(offsets);
    // monotonic, and clamped to the table
    var max: u64 = 0;
    for (offsets) |*o| {
        const raw: u64 = if (long) c.ulong() else @as(u64, c.ushort()) * 2;
        var v = @as(u64, data_at) + raw;
        if (v >= max) max = v else v = max;
        o.* = @intCast(@min(v, table.len));
    }
    c.seek(shared_at);
    const n = @as(usize, shared_count) * axis_count;
    if (shared_at > table.len or c.left() < 2 * n) return error.InvalidTable;
    const shared = try allocator.alloc(i32, n);
    for (shared) |*v| v.* = @as(i32, c.short()) * 4;
    return .{ .table = table, .offsets = offsets, .shared = shared };
}

pub fn deinit(g: *Gvar, allocator: Allocator) void {
    allocator.free(g.shared);
    allocator.free(g.offsets);
}

/// `TT_Vary_Apply_Glyph_Deltas`: glyph `gid`'s deltas at `coords` added to
/// `points` (font units: its points, then the four phantom points) rounded
/// to whole units, and to `unrounded` (26.6 font units, set from `points`
/// first). `ends` are the contours' last points; with `hvar` the phantom
/// points stay put.
pub fn apply(g: *const Gvar, allocator: Allocator, gid: u16, coords: []const i32, points: []Point, ends: []const u16, hvar: bool, unrounded: []Point) Error!void {
    for (unrounded, points) |*u, p| u.* = .{ p[0] * 64, p[1] * 64 };
    if (@as(usize, gid) + 1 >= g.offsets.len) return;
    const frame = g.table[g.offsets[gid]..g.offsets[gid + 1]];
    if (frame.len == 0) return;
    var head: Cursor = .{ .data = frame };
    const count = head.ushort() & 0x0FFF;
    const data_offset = head.ushort();
    if (data_offset > frame.len or @as(usize, count) * 4 > frame.len) return error.InvalidTable;

    const n = points.len;
    const buf = try allocator.alloc(Point, 3 * n);
    defer allocator.free(buf);
    const acc = buf[0..n];
    const org = buf[n .. 2 * n];
    const out = buf[2 * n ..];
    const has_delta = try allocator.alloc(bool, n);
    defer allocator.free(has_delta);
    @memset(acc, .{ 0, 0 });
    for (org, points) |*o, p| o.* = .{ p[0] << 16, p[1] << 16 };

    var it = try tuples.Iterator.init(allocator, frame, 0, coords, g.shared);
    defer it.deinit();
    while (try it.next()) |t| {
        defer it.release(t);
        const pts = t.points orelse continue;
        const len = switch (pts) {
            .all => n,
            .list => |l| l.len,
        };
        var d = t.deltas;
        const dx = try packed_data.readDeltas(allocator, &d, len);
        defer if (dx) |v| allocator.free(v);
        const dy = try packed_data.readDeltas(allocator, &d, len);
        defer if (dy) |v| allocator.free(v);
        if (dx == null or dy == null) continue;
        switch (pts) {
            .all => for (acc, dx.?, dy.?) |*a, x, y| {
                a.* = .{ a[0] +% fixed.mul(x, t.scalar), a[1] +% fixed.mul(y, t.scalar) };
            },
            .list => |l| {
                @memset(has_delta, false);
                @memcpy(out, org);
                for (l, dx.?, dy.?) |idx, x, y| {
                    if (idx >= n) continue;
                    has_delta[idx] = true;
                    out[idx] = .{ out[idx][0] +% fixed.mul(x, t.scalar), out[idx][1] +% fixed.mul(y, t.scalar) };
                }
                // phantom points are in no contour, so never inferred
                iup.interpolate(ends, out, org, has_delta);
                for (acc, out, org) |*a, o, r| a.* = .{ a[0] +% (o[0] -% r[0]), a[1] +% (o[1] -% r[1]) };
            },
        }
    }
    if (hvar) @memset(acc[n - 4 .. n - 2], .{ 0, 0 });
    for (acc, points, unrounded) |a, *p, *u| for (0..2) |k| {
        // FT_fixedToFdot6 and FT_fixedToInt
        u[k] +%= (a[k] +% 0x200) >> 10;
        p[k] +%= (a[k] +% 0x8000) >> 16;
    };
}
