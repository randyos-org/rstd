//! `HVAR`'s advance width deltas, as FreeType 2.14.3 reads and applies
//! them (`src/truetype/ttgxvar.c` `ft_var_load_hvvar`,
//! `tt_hvadvance_adjust`): an invalid table is not used, and the phantom
//! points' `gvar` deltas give the advance instead.

const std = @import("std");
const ItemStore = @import("ItemStore.zig");
const Cursor = @import("packed.zig").Cursor;

const Allocator = std.mem.Allocator;
const Hvar = @This();

store: ItemStore,
widths: ?ItemStore.IndexMap,

/// The table, or null if FreeType would reject it.
pub fn init(allocator: Allocator, table: []const u8, axis_count: u16) Allocator.Error!?Hvar {
    var c: Cursor = .{ .data = table };
    if (c.left() < 12 or c.ushort() != 1) return null;
    _ = c.ushort();
    const store_offset = c.ulong();
    const map_offset = c.ulong();
    var store = ItemStore.init(allocator, table, store_offset, axis_count) catch |e| switch (e) {
        error.OutOfMemory => return error.OutOfMemory,
        error.InvalidTable => return null,
    };
    errdefer store.deinit(allocator);
    var widths: ?ItemStore.IndexMap = null;
    if (map_offset != 0) widths = ItemStore.IndexMap.init(allocator, table, map_offset, &store) catch |e| switch (e) {
        error.OutOfMemory => return error.OutOfMemory,
        error.InvalidTable => {
            store.deinit(allocator);
            return null;
        },
    };
    return .{ .store = store, .widths = widths };
}

pub fn deinit(h: *Hvar, allocator: Allocator) void {
    if (h.widths) |*w| w.deinit(allocator);
    h.store.deinit(allocator);
}

/// Glyph `gid`'s advance delta at `coords`, in font units.
pub fn advanceDelta(h: *const Hvar, gid: u16, coords: []const i32) i32 {
    const direct: ItemStore.Index = .{ .outer = 0, .inner = gid };
    const idx = if (h.widths) |w| w.get(gid) orelse direct else direct;
    return h.store.delta(idx, coords);
}
