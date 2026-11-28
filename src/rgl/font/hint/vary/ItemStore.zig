//! An item variation store (`HVAR`, `avar` 2) and its delta-set index map,
//! read and evaluated as FreeType 2.14.3 does (`src/truetype/ttgxvar.c`
//! `tt_var_load_item_variation_store`, `tt_var_load_delta_set_index_mapping`,
//! `tt_calculate_scalar`, `tt_var_get_item_delta`): regions in 16.16, an
//! item's deltas summed in 64 bits and rounded once.

const std = @import("std");
const fixed = @import("../fixed.zig");
const Cursor = @import("packed.zig").Cursor;

const Allocator = std.mem.Allocator;
const ItemStore = @This();

pub const Error = error{InvalidTable} || Allocator.Error;

pub const Index = struct { outer: u16, inner: u16 };

const VarData = struct {
    item_count: u16,
    word_count: u16,
    long_words: bool,
    regions: []u16,
    /// `item_count` rows of deltas, borrowed from the table.
    deltas: []const u8,
};

axis_count: u16,
/// Each region's start, peak and end per axis, 16.16; an invalid range
/// has peak 0, so it ignores its axis.
regions: []i32,
data: []VarData,

/// The store at `offset` in `table`, for a font of `axis_count` axes.
pub fn init(allocator: Allocator, table: []const u8, offset: usize, axis_count: u16) Error!ItemStore {
    if (offset >= table.len) return error.InvalidTable;
    var c: Cursor = .{ .data = table[offset..] };
    if (c.left() < 8 or c.ushort() != 1) return error.InvalidTable;
    const region_offset = c.ulong();
    const data_count = c.ushort();
    if (data_count == 0 or c.left() < 4 * @as(usize, data_count)) return error.InvalidTable;
    const data_offsets = c.at;

    c.seek(region_offset);
    if (c.left() < 4) return error.InvalidTable;
    if (c.ushort() != axis_count) return error.InvalidTable;
    const region_count = c.ushort();
    if (region_count >= 32768 or c.left() < @as(usize, region_count) * axis_count * 6) return error.InvalidTable;
    const regions = try allocator.alloc(i32, @as(usize, region_count) * axis_count * 3);
    errdefer allocator.free(regions);
    for (0..@as(usize, region_count) * axis_count) |i| {
        const start: i32 = c.short();
        var peak: i32 = c.short();
        const end: i32 = c.short();
        if ((start < 0 and end > 0) or start > peak or peak > end) peak = 0;
        regions[3 * i ..][0..3].* = .{ start * 4, peak * 4, end * 4 };
    }

    const data = try allocator.alloc(VarData, data_count);
    var done: usize = 0;
    errdefer {
        for (data[0..done]) |d| allocator.free(d.regions);
        allocator.free(data);
    }
    for (data, 0..) |*d, i| {
        c.seek(data_offsets + 4 * i);
        c.seek(c.ulong());
        if (c.left() < 6) return error.InvalidTable;
        const item_count = c.ushort();
        const word_raw = c.ushort();
        const region_idx_count = c.ushort();
        const word_count = word_raw & 0x7FFF;
        const long_words = word_raw & 0x8000 != 0;
        if (word_count > region_idx_count or region_idx_count > region_count) return error.InvalidTable;
        if (c.left() < 2 * @as(usize, region_idx_count)) return error.InvalidTable;
        const idx = try allocator.alloc(u16, region_idx_count);
        for (idx) |*r| {
            r.* = c.ushort();
            if (r.* >= region_count) {
                allocator.free(idx);
                return error.InvalidTable;
            }
        }
        var row: usize = @as(usize, word_count) + region_idx_count;
        if (long_words) row *= 2;
        if (c.left() < row * item_count) {
            allocator.free(idx);
            return error.InvalidTable;
        }
        d.* = .{ .item_count = item_count, .word_count = word_count, .long_words = long_words, .regions = idx, .deltas = c.data[c.at..][0 .. row * item_count] };
        done += 1;
    }
    return .{ .axis_count = axis_count, .regions = regions, .data = data };
}

pub fn deinit(s: *ItemStore, allocator: Allocator) void {
    for (s.data) |d| allocator.free(d.regions);
    allocator.free(s.data);
    allocator.free(s.regions);
}

/// `tt_calculate_scalar`: how much a region applies at `coords`.
fn regionScalar(s: *const ItemStore, region: usize, coords: []const i32) i32 {
    var scalar: i32 = 0x10000;
    const axes = s.regions[region * s.axis_count * 3 ..][0 .. s.axis_count * 3];
    for (0..s.axis_count) |j| {
        const start, const peak, const end = axes[3 * j ..][0..3].*;
        const ncv = coords[j];
        if (peak == ncv or peak == 0) continue;
        if (ncv <= start or ncv >= end) return 0;
        scalar = if (ncv < peak) fixed.mulDiv(scalar, ncv -% start, peak -% start) else fixed.mulDiv(scalar, end -% ncv, end -% peak);
    }
    return scalar;
}

/// `tt_var_get_item_delta`: an item's delta at `coords`, 0 out of range.
pub fn delta(s: *const ItemStore, idx: Index, coords: []const i32) i32 {
    if (idx.outer == 0xFFFF and idx.inner == 0xFFFF) return 0;
    if (idx.outer >= s.data.len) return 0;
    const d = s.data[idx.outer];
    if (idx.inner >= d.item_count or d.regions.len == 0) return 0;
    var row: usize = @as(usize, d.word_count) + d.regions.len;
    if (d.long_words) row *= 2;
    var c: Cursor = .{ .data = d.deltas[row * idx.inner ..][0..row] };
    var sum: i64 = 0;
    for (d.regions, 0..) |r, m| {
        const wide = m < d.word_count;
        const v: i32 = if (d.long_words)
            (if (wide) @bitCast(c.ulong()) else c.short())
        else
            (if (wide) c.short() else @as(i8, @bitCast(c.byte())));
        const scalar = s.regionScalar(r, coords);
        if (scalar != 0) sum += @as(i64, v) * scalar;
    }
    // FreeType truncates to 32 bits before the shift
    return @as(i32, @truncate(sum + 0x8000)) >> 16;
}

/// A delta-set index map: an item's outer and inner index by position.
pub const IndexMap = struct {
    indices: []Index,

    pub fn init(allocator: Allocator, table: []const u8, offset: usize, store: *const ItemStore) Error!IndexMap {
        if (offset >= table.len) return error.InvalidTable;
        var c: Cursor = .{ .data = table[offset..] };
        if (c.left() < 2) return error.InvalidTable;
        const format = c.byte();
        const entry_format = c.byte();
        const count: usize = switch (format) {
            0 => if (c.left() < 2) return error.InvalidTable else c.ushort(),
            1 => if (c.left() < 4) return error.InvalidTable else c.ulong(),
            else => return error.InvalidTable,
        };
        if (entry_format & 0xC0 != 0) return error.InvalidTable;
        const entry_size: usize = ((entry_format & 0x30) >> 4) + 1;
        const inner_bits: u5 = @intCast((entry_format & 0x0F) + 1);
        if (count * entry_size > table.len or c.left() < count * entry_size) return error.InvalidTable;
        const indices = try allocator.alloc(Index, count);
        errdefer allocator.free(indices);
        for (indices) |*e| {
            var v: u32 = 0;
            for (0..entry_size) |_| v = (v << 8) | c.byte();
            if (v == 0xFFFFFFFF) {
                e.* = .{ .outer = 0xFFFF, .inner = 0xFFFF };
                continue;
            }
            const outer = v >> inner_bits;
            const inner = v & ((@as(u32, 1) << inner_bits) - 1);
            if (outer >= store.data.len or inner >= store.data[outer].item_count) return error.InvalidTable;
            e.* = .{ .outer = @intCast(outer), .inner = @intCast(inner) };
        }
        return .{ .indices = indices };
    }

    pub fn deinit(m: *IndexMap, allocator: Allocator) void {
        allocator.free(m.indices);
    }

    /// Positions past the end take the last entry; an empty map, none
    /// (FreeType allocates nothing for it).
    pub fn get(m: IndexMap, i: usize) ?Index {
        if (m.indices.len == 0) return null;
        return m.indices[@min(i, m.indices.len - 1)];
    }
};
