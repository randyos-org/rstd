//! The tuples of a `gvar` glyph's or `cvar`'s variation data that apply at
//! an instance, walked as FreeType 2.14.3 walks them
//! (`src/truetype/ttgxvar.c` `TT_Vary_Apply_Glyph_Deltas`,
//! `tt_face_vary_cvt`): each one's scalar, its points and its deltas.

const std = @import("std");
const packed_data = @import("packed.zig");
const Blend = @import("Blend.zig");

const Allocator = std.mem.Allocator;
const Cursor = packed_data.Cursor;
const Points = packed_data.Points;

pub const Error = error{InvalidTable} || Allocator.Error;

const embedded_peak = 0x8000;
const intermediate = 0x4000;
const private_points = 0x2000;
const index_mask = 0x0FFF;
const shared_points = 0x8000;
const count_mask = 0x0FFF;

/// One applying tuple; `deltas` is at its serialized data.
pub const Tuple = struct {
    scalar: i32,
    /// Null if they could not be read: the tuple is then skipped.
    points: ?Points,
    deltas: Cursor,
    owned: bool,
};

pub const Iterator = struct {
    allocator: Allocator,
    /// The whole frame: the glyph's data, or `cvar`.
    frame: Cursor,
    /// At the next tuple header.
    headers: Cursor,
    /// Where the next tuple's serialized data starts.
    data_at: usize,
    left: usize,
    coords: []const i32,
    shared_tuples: []const i32,
    shared: ?Points = null,
    scratch: []i32,

    /// The variation data in `frame` whose count and data offset are at
    /// `at`; `shared_tuples` are `gvar`'s, `coords.len` per tuple.
    pub fn init(allocator: Allocator, frame: []const u8, at: usize, coords: []const i32, shared_tuples: []const i32) Error!Iterator {
        var c: Cursor = .{ .data = frame };
        c.seek(at);
        const count = c.ushort();
        // the data offset counts from the frame's start
        var data_at: usize = c.ushort();
        var it: Iterator = .{
            .allocator = allocator,
            .frame = .{ .data = frame },
            .headers = c,
            .data_at = data_at,
            .left = count & count_mask,
            .coords = coords,
            .shared_tuples = shared_tuples,
            .scratch = try allocator.alloc(i32, 3 * coords.len),
        };
        errdefer allocator.free(it.scratch);
        if (count & shared_points != 0) {
            var d: Cursor = .{ .data = frame };
            d.seek(data_at);
            it.shared = try packed_data.readPoints(allocator, &d);
            data_at = d.at;
            it.data_at = data_at;
        }
        return it;
    }

    pub fn deinit(it: *Iterator) void {
        if (it.shared) |s| s.deinit(it.allocator);
        it.allocator.free(it.scratch);
    }

    /// The next tuple that applies at `coords`, or null after the last.
    pub fn next(it: *Iterator) Error!?Tuple {
        const n = it.coords.len;
        while (it.left > 0) {
            it.left -= 1;
            const h = &it.headers;
            if (h.left() < 4) return error.InvalidTable;
            const size = h.ushort();
            const index = h.ushort();
            const at = it.data_at;
            it.data_at += size;
            var peak: []const i32 = undefined;
            if (index & embedded_peak != 0) {
                if (h.left() < 2 * n) return error.InvalidTable;
                for (it.scratch[0..n]) |*v| v.* = @as(i32, h.short()) * 4;
                peak = it.scratch[0..n];
            } else if (index & index_mask < it.shared_tuples.len / @max(n, 1)) {
                peak = it.shared_tuples[(index & index_mask) * n ..][0..n];
            } else return error.InvalidTable;
            var range: ?[2][]const i32 = null;
            if (index & intermediate != 0) {
                if (h.left() < 4 * n) return error.InvalidTable;
                for (it.scratch[n..]) |*v| v.* = @as(i32, h.short()) * 4;
                range = .{ it.scratch[n .. 2 * n], it.scratch[2 * n ..] };
            }
            const scalar = Blend.tupleScalar(it.coords, peak, range);
            if (scalar == 0) continue;
            var d = it.frame;
            d.seek(at);
            if (index & private_points != 0) {
                const local = try packed_data.readPoints(it.allocator, &d);
                return .{ .scalar = scalar, .points = local, .deltas = d, .owned = true };
            }
            return .{ .scalar = scalar, .points = it.shared, .deltas = d, .owned = false };
        }
        return null;
    }

    /// Frees a tuple's own points.
    pub fn release(it: *Iterator, t: Tuple) void {
        if (t.owned) if (t.points) |p| p.deinit(it.allocator);
    }
};
