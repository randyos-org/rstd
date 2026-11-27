//! Each glyph's style, from the Unicode charmap: FreeType 2.14.3
//! `src/autofit/afglobal.c` (`af_face_globals_compute_style_coverage`)
//! without HarfBuzz, where only the default coverages count.

const std = @import("std");
const Sfnt = @import("../../sfnt/Sfnt.zig");
const scripts = @import("data/scripts.zig").scripts;
const styles = @import("data/styles.zig").styles;
const Cmap = @import("../../sfnt/tables/cmap/__root__.zig").CmapSubtable;

const Allocator = std.mem.Allocator;

/// The style's index in a glyph's entry (`AF_STYLE_MASK`).
pub const style_mask: u16 = 0x1FFF;
pub const unassigned: u16 = style_mask;
pub const digit: u16 = 0x8000;
pub const nonbase: u16 = 0x4000;
pub const has_cmap_entry: u16 = 0x2000;

/// The style of glyphs no script covers (`AF_STYLE_FALLBACK`, CJK's, as
/// FreeType is built with `AF_CONFIG_OPTION_CJK`).
pub const fallback_style: u16 = styleIndex("hani_dflt");
/// The style of glyphs left unhinted (`AF_STYLE_NONE_DFLT`).
pub const none_style: u16 = styleIndex("none_dflt");

fn styleIndex(comptime name: []const u8) u16 {
    for (styles, 0..) |s, i| if (std.mem.eql(u8, s.name, name)) return i;
    @compileError("no style " ++ name);
}

/// A glyph style entry per glyph, owned by the caller.
pub fn compute(allocator: Allocator, sfnt: *const Sfnt) Allocator.Error![]u16 {
    const count: u32 = if (sfnt.maxp) |m| m.numGlyphs() else 0;
    const gstyles = try allocator.alloc(u16, count);
    @memset(gstyles, unassigned);
    if (sfnt.unicodeCmap()) |cmap| {
        for (styles, 0..) |style, ss_usize| {
            if (!style.default_coverage) continue;
            const ss: u16 = @intCast(ss_usize);
            const script = scripts[style.script];
            for (script.ranges) |r| {
                var it = rangeIter(cmap, r);
                while (it.next()) |g| {
                    if (g < count and gstyles[g] & style_mask == unassigned) gstyles[g] = ss | has_cmap_entry;
                }
            }
            for (script.nonbase) |r| {
                var it = rangeIter(cmap, r);
                while (it.next()) |g| {
                    if (g < count and gstyles[g] & style_mask == ss) gstyles[g] |= nonbase;
                }
            }
        }
        var c: u32 = 0x30;
        while (c <= 0x39) : (c += 1) {
            const g = @backingInt(cmap.getGlyphIndex(c));
            if (g != 0 and g < count) gstyles[g] |= digit;
        }
    }
    for (gstyles) |*s| {
        if (s.* & style_mask == unassigned) s.* = (s.* & ~style_mask) | fallback_style;
    }
    return gstyles;
}

/// The glyphs of a range's mapped codes, as `FT_Get_Char_Index` then
/// `FT_Get_Next_Char` walk it.
const RangeIter = struct {
    cmap: Cmap,
    code: u32,
    last: u32,
    started: bool = false,

    fn next(it: *RangeIter) ?u32 {
        if (!it.started) {
            it.started = true;
            const g = @backingInt(it.cmap.getGlyphIndex(it.code));
            if (g != 0) return g;
        }
        const m = it.cmap.nextChar(it.code) orelse return null;
        if (m.codepoint > it.last) return null;
        it.code = m.codepoint;
        return @backingInt(m.glyph);
    }
};

fn rangeIter(cmap: Cmap, r: [2]u32) RangeIter {
    return .{ .cmap = cmap, .code = r[0], .last = r[1] };
}
