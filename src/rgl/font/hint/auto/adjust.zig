//! The diacritic adjustments' lookups: FreeType 2.14.3
//! `src/autofit/afadjust.c` (`af_adjustment_database_lookup`, and
//! `af_reverse_character_map_new` without HarfBuzz) and the flags of
//! `afadjust.h`.

const std = @import("std");
const Sfnt = @import("../../sfnt/Sfnt.zig");
const entries = @import("data/adjust.zig").entries;

const Allocator = std.mem.Allocator;

pub const up: u32 = 0x01;
pub const down: u32 = 0x02;
pub const up2: u32 = 0x04;
pub const down2: u32 = 0x08;
pub const tilde_top: u32 = 0x10;
pub const tilde_bottom: u32 = 0x20;
pub const tilde_top2: u32 = 0x40;
pub const tilde_bottom2: u32 = 0x80;
pub const ignore_capital_top: u32 = 0x100;
pub const ignore_capital_bottom: u32 = 0x200;
pub const ignore_small_top: u32 = 0x400;
pub const ignore_small_bottom: u32 = 0x800;
pub const no_height_check: u32 = 0x1000;

/// `af_adjustment_database_lookup`: FreeType's binary search over its
/// table, which misses the few entries out of order.
pub fn lookup(codepoint: u32) u32 {
    var low: usize = 0;
    var high: usize = entries.len - 1;
    while (high >= low) {
        const mid = (low + high) / 2;
        const c = entries[mid][0];
        if (c < codepoint) {
            low = mid + 1;
        } else if (c > codepoint) {
            // FreeType's unsigned index would wrap; no code it looks up gets here
            if (mid == 0) return 0;
            high = mid - 1;
        } else return entries[mid][1];
    }
    return 0;
}

/// Each glyph the database's code points map to, and its code point: the
/// last entry in the table's order wins.
pub const ReverseMap = struct {
    map: std.AutoHashMapUnmanaged(u16, u32) = .empty,

    pub fn init(allocator: Allocator, sfnt: *const Sfnt) Allocator.Error!ReverseMap {
        var r: ReverseMap = .{};
        errdefer r.map.deinit(allocator);
        const cmap = sfnt.unicodeCmap() orelse return r;
        for (entries) |e| {
            const g = @backingInt(cmap.getGlyphIndex(e[0]));
            if (g == 0) continue;
            try r.map.put(allocator, g, e[0]);
        }
        return r;
    }

    pub fn deinit(r: *ReverseMap, allocator: Allocator) void {
        r.map.deinit(allocator);
    }

    /// The glyph's adjustment flags, 0 for none.
    pub fn flags(r: *const ReverseMap, glyph: u16) u32 {
        const cp = r.map.get(glyph) orelse return 0;
        return lookup(cp);
    }
};
