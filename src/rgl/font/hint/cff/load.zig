//! Loading a CFF glyph scaled and hinted as FreeType 2.14.3 loads it
//! (`src/cff/cffgload.c` `cff_slot_load`, `src/psaux/psft.c`
//! `cf2_decoder_parse_charstrings` and `cf2_getScaleAndHintFlag`,
//! `src/base/ftobjs.c`'s grid fitting of the advance): hinted, the engine
//! draws in 26.6 pixels; unhinted, in font units scaled afterwards. A
//! FontMatrix other than the default one and SEAC accents are not
//! supported.

const std = @import("std");
const fixed = @import("../fixed.zig");
const charstring = @import("charstring.zig");
const private = @import("private.zig");
const Blues = @import("blues.zig").Blues;
const Builder = @import("Builder.zig");
const GlyphPath = @import("GlyphPath.zig");
const Instance = @import("Instance.zig");
const Glyph = @import("../Glyph.zig");
const GlyphIndex = @import("../../sfnt/types.zig").GlyphIndex;

const Allocator = std.mem.Allocator;

pub const Error = error{ InvalidGlyphIndex, TooManyHints } || charstring.Error || private.Error || Allocator.Error;

/// `glyph` at `inst`'s size, hinted or not; the caller owns its slices.
pub fn load(allocator: Allocator, inst: *const Instance, glyph: GlyphIndex, hinted: bool) Error!Glyph {
    const cff = inst.cff;
    if (@backingInt(glyph) >= cff.charstrings_index.count()) return error.InvalidGlyphIndex;
    const priv: private.Private = if (cff.privateDictForGlyph(glyph)) |d| try private.parse(d.data) else .{};
    // the engine's scale has no factor of 64; unhinted, it draws in font units
    const scale: i32 = if (hinted) @divTrunc(inst.scale +% 32, 64) else 0x400;
    const blues: Blues = .init(&priv, scale);
    var out: Builder = .{ .allocator = allocator };
    defer out.deinit();
    var path: GlyphPath = .init(&blues, &out, scale, scale, hinted);
    const ctx: charstring.Context = .{ .global = cff.global_subrs_index, .local = cff.getLocalSubrIndexForGlyph(glyph) };
    const width = try charstring.evaluate(ctx, cff.getCharString(glyph), &path);
    if (path.err) |e| return e;
    if (path.too_many_hints) return error.TooManyHints;
    if (!hinted) for (out.points.items) |*p| {
        p[0] = fixed.mul(p[0], inst.scale);
        p[1] = fixed.mul(p[1], inst.scale);
    };
    // hmtx's advance, else the charstring's
    var advance: i32 = if (inst.sfnt.hmtx) |hmtx|
        hmtx.getMetrics(glyph).advance_width
    else if (width) |w|
        (w +% priv.nominal_width +% 0x8000) >> 16
    else
        (priv.default_width orelse 0) >> 16;
    advance = fixed.mul(advance, inst.scale);
    if (hinted) advance = (advance + 32) & ~@as(i32, 63);
    const points = try allocator.dupe([2]i32, out.points.items);
    errdefer allocator.free(points);
    const tags = try allocator.dupe(u8, out.tags.items);
    errdefer allocator.free(tags);
    return .{ .points = points, .tags = tags, .contours = try allocator.dupe(u16, out.contours.items), .advance = advance };
}
