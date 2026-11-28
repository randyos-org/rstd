//! A CFF2 font's variation store, its regions and each item variation
//! data's region list (CFF2 has no delta sets), and the blend vector a
//! `blend` weighs its deltas by, as FreeType 2.14.3 reads and builds them
//! (`src/cff/cffload.c` `cff_vstore_load`, `cff_blend_build_vector`).

const std = @import("std");
const fixed = @import("../fixed.zig");
const Cursor = @import("../vary/packed.zig").Cursor;

const Allocator = std.mem.Allocator;
const VarStore = @This();

pub const Error = error{
    /// The store does not parse, or a vector asks for what it lacks.
    InvalidVariationStore,
} || Allocator.Error;

axis_count: u16,
/// Start, peak and end per axis per region, 16.16; an invalid range has
/// peak 0, so it ignores its axis.
regions: []i32,
/// Each item variation data's region indexes, in one array.
indexes: []u16,
/// `indexes[starts[i]..starts[i + 1]]` are data `i`'s.
starts: []usize,

/// The store at `offset` in the CFF2 table: a 16-bit length, then an item
/// variation store.
pub fn init(allocator: Allocator, table: []const u8, offset: usize) Error!VarStore {
    if (offset + 2 > table.len) return error.InvalidVariationStore;
    const base = offset + 2;
    var c: Cursor = .{ .data = table[base..] };
    if (c.left() < 8 or c.ushort() != 1) return error.InvalidVariationStore;
    const region_list = c.ulong();
    const data_count = c.ushort();
    if (c.left() < 4 * @as(usize, data_count)) return error.InvalidVariationStore;
    const data_at = c.at;
    c.seek(region_list);
    if (c.left() < 4) return error.InvalidVariationStore;
    const axis_count = c.ushort();
    const region_count = c.ushort();
    const n = @as(usize, region_count) * axis_count;
    if (c.left() < 6 * n) return error.InvalidVariationStore;
    const regions = try allocator.alloc(i32, 3 * n);
    errdefer allocator.free(regions);
    for (0..n) |i| {
        const start: i32 = c.short();
        var peak: i32 = c.short();
        const end: i32 = c.short();
        if ((start < 0 and end > 0) or start > peak or peak > end) peak = 0;
        regions[3 * i ..][0..3].* = .{ start * 4, peak * 4, end * 4 };
    }
    const starts = try allocator.alloc(usize, @as(usize, data_count) + 1);
    errdefer allocator.free(starts);
    var list: std.ArrayList(u16) = .empty;
    errdefer list.deinit(allocator);
    for (0..data_count) |i| {
        starts[i] = list.items.len;
        c.seek(data_at + 4 * i);
        c.seek(c.ulong());
        // itemCount and wordDeltaCount: CFF2 has no delta sets
        if (c.left() < 6) return error.InvalidVariationStore;
        c.seek(c.at + 4);
        const k = c.ushort();
        if (c.left() < 2 * @as(usize, k)) return error.InvalidVariationStore;
        for (0..k) |_| try list.append(allocator, c.ushort());
    }
    starts[data_count] = list.items.len;
    return .{ .axis_count = axis_count, .regions = regions, .indexes = try list.toOwnedSlice(allocator), .starts = starts };
}

pub fn deinit(v: *VarStore, allocator: Allocator) void {
    allocator.free(v.starts);
    allocator.free(v.indexes);
    allocator.free(v.regions);
}

/// `cff_blend_build_vector`: 1, then each of data `vsindex`'s regions'
/// scalars at `coords` (16.16, normalized; none for the default instance,
/// whose vector is 1, 0, 0, ...). The caller owns it.
pub fn blendVector(v: *const VarStore, allocator: Allocator, vsindex: usize, coords: []const i32) Error![]i32 {
    if (coords.len != 0 and coords.len != v.axis_count) return error.InvalidVariationStore;
    if (vsindex + 1 >= v.starts.len) return error.InvalidVariationStore;
    const idx = v.indexes[v.starts[vsindex]..v.starts[vsindex + 1]];
    const bv = try allocator.alloc(i32, idx.len + 1);
    errdefer allocator.free(bv);
    bv[0] = 0x10000;
    const region_count = v.regions.len / 3 / @max(v.axis_count, 1);
    for (idx, bv[1..]) |r, *w| {
        if (r >= region_count) return error.InvalidVariationStore;
        if (coords.len == 0) {
            w.* = 0;
            continue;
        }
        w.* = 0x10000;
        const axes = v.regions[@as(usize, r) * v.axis_count * 3 ..][0 .. v.axis_count * 3];
        for (coords, 0..) |ncv, j| {
            const start, const peak, const end = axes[3 * j ..][0..3].*;
            if (peak == ncv or peak == 0) continue;
            if (ncv <= start or ncv >= end) {
                w.* = 0;
                break;
            }
            w.* = if (ncv < peak) fixed.mulDiv(w.*, ncv -% start, peak -% start) else fixed.mulDiv(w.*, end -% ncv, end -% peak);
        }
    }
    return bv;
}

/// The blend vector of the last `vsindex` asked for, rebuilt when it
/// changes (`cff_blend_check_vector`).
pub const Vectors = struct {
    allocator: Allocator,
    store: *const VarStore,
    coords: []const i32,
    vector: ?[]i32 = null,
    vsindex: usize = 0,

    pub fn get(s: *Vectors, vsindex: usize) Error![]const i32 {
        if (s.vector) |v| {
            if (s.vsindex == vsindex) return v;
            s.allocator.free(v);
            s.vector = null;
        }
        s.vector = try s.store.blendVector(s.allocator, vsindex, s.coords);
        s.vsindex = vsindex;
        return s.vector.?;
    }

    pub fn deinit(s: *Vectors) void {
        if (s.vector) |v| s.allocator.free(v);
    }
};
