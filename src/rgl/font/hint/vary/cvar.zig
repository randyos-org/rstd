//! `cvar`'s deltas on the control values, as FreeType 2.14.3 applies them
//! (`src/truetype/ttgxvar.c` `tt_face_vary_cvt`): summed in 16.16 per
//! value, then rounded to 26.6 font units. A missing or malformed table
//! leaves the values as they are.

const std = @import("std");
const fixed = @import("../fixed.zig");
const packed_data = @import("packed.zig");
const tuples = @import("tuples.zig");
const Cursor = packed_data.Cursor;

const Allocator = std.mem.Allocator;

/// Adds the deltas at `coords` to `cvt` (26.6 font units); `shared` are
/// `gvar`'s shared tuples, which FreeType lets `cvar` use.
pub fn apply(allocator: Allocator, cvar: []const u8, coords: []const i32, shared: []const i32, cvt: []i32) Allocator.Error!void {
    var c: Cursor = .{ .data = cvar };
    if (c.left() < 8 or c.ulong() != 0x10000) return;
    const count = c.ushort() & 0x0FFF;
    const data_offset = c.ushort();
    if (@as(usize, data_offset) + @as(usize, count) * 4 > cvar.len) return;
    const deltas = try allocator.alloc(i32, cvt.len);
    defer allocator.free(deltas);
    @memset(deltas, 0);
    var it = tuples.Iterator.init(allocator, cvar, 4, coords, shared) catch |e| switch (e) {
        error.OutOfMemory => return error.OutOfMemory,
        error.InvalidTable => return,
    };
    defer it.deinit();
    while (it.next() catch |e| switch (e) {
        error.OutOfMemory => return error.OutOfMemory,
        // FreeType keeps no deltas from a table it fails on
        error.InvalidTable => return,
    }) |t| {
        defer it.release(t);
        const pts = t.points orelse continue;
        var d = t.deltas;
        const len = switch (pts) {
            .all => cvt.len,
            .list => |l| l.len,
        };
        const dv = try packed_data.readDeltas(allocator, &d, len) orelse continue;
        defer allocator.free(dv);
        switch (pts) {
            .all => for (deltas, dv) |*a, v| {
                a.* +%= fixed.mul(v, t.scalar);
            },
            .list => |l| for (l, dv) |idx, v| {
                if (idx >= cvt.len) continue;
                deltas[idx] +%= fixed.mul(v, t.scalar);
            },
        }
    }
    for (cvt, deltas) |*v, d| v.* +%= (d +% 0x200) >> 10;
}
