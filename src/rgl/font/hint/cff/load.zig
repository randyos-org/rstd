//! Loading a CFF glyph scaled and hinted as FreeType 2.14.3 loads it
//! (`src/cff/cffgload.c` `cff_slot_load`, `src/psaux/psft.c`
//! `cf2_decoder_parse_charstrings` and `cf2_getScaleAndHintFlag`,
//! `src/base/ftobjs.c`'s grid fitting of the advance): hinted, the engine
//! draws in 26.6 pixels; unhinted, in font units scaled afterwards. A
//! SEAC accent is drawn, then its base, each with stem hints of its own
//! (`psintrp.c` `cf2_cmdENDCHAR`, `psft.c` `cf2_getSeacComponent`). The
//! FontMatrix (`FontMatrix.zig`) then transforms the drawn outline and
//! offsets it in font units, and a CID subfont's units per em rescale it.

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
const charset = @import("../../sfnt/tables/cff/charset.zig");
const FontMatrix = @import("FontMatrix.zig");

const Allocator = std.mem.Allocator;

pub const Error = error{
    InvalidGlyphIndex,
    TooManyHints,
    /// A SEAC names a code no glyph has, or a component is a SEAC itself.
    InvalidSeac,
} || charstring.Error || private.Error || Allocator.Error;

/// `glyph` at `inst`'s size, hinted or not; the caller owns its slices.
pub fn load(allocator: Allocator, inst: *const Instance, glyph: GlyphIndex, hinted: bool) Error!Glyph {
    const cff = inst.cff;
    if (@backingInt(glyph) >= cff.charstrings_index.count()) return error.InvalidGlyphIndex;
    const priv: private.Private = if (cff.privateDictForGlyph(glyph)) |d| try private.parse(d.data) else .{};
    const face_upem: u16 = if (inst.sfnt.head) |h| h.unitsPerEm() else 1000;
    const fm = try FontMatrix.forGlyph(cff, face_upem, glyph);
    // a subfont of other units per em scales by their ratio
    const glyph_scale = if (fm.top_upm != fm.sub_upm) fixed.mulDiv(inst.scale, fm.top_upm, fm.sub_upm) else inst.scale;
    // the engine's scale has no factor of 64; unhinted, it draws in font units
    const scale: i32 = if (hinted) @divTrunc(glyph_scale +% 32, 64) else 0x400;
    const blues: Blues = .init(&priv, scale);
    var out: Builder = .{ .allocator = allocator };
    defer out.deinit();
    const ctx: charstring.Context = .{ .global = cff.global_subrs_index, .local = cff.getLocalSubrIndexForGlyph(glyph) };
    const result = try draw(ctx, cff.getCharString(glyph), &blues, &out, scale, hinted, .{ 0, 0 });
    if (result.seac) |seac| {
        // the accent at its offset, its stems counted from it, then the base
        const parts = [_]struct { code: i32, at: [2]i32 }{ .{ .code = seac.accent, .at = .{ seac.adx, seac.ady } }, .{ .code = seac.base, .at = .{ 0, 0 } } };
        for (parts) |part| {
            const gid = charset.glyphForStdCode(cff, part.code) orelse return error.InvalidSeac;
            const r = try draw(ctx, cff.getCharString(@fromBackingInt(gid)), &blues, &out, scale, hinted, part.at);
            if (r.seac != null) return error.InvalidSeac;
        }
    }
    out.closeContour();
    const width = result.width;
    const transformed = !std.mem.eql(i32, &fm.matrix, &FontMatrix.identity);
    // on the drawn outline, hinted or not, as `cff_slot_load` applies them
    for (out.points.items) |*p| {
        if (transformed) p.* = FontMatrix.transform(p.*, fm.matrix);
        p.* = .{ p[0] +% fm.offset[0], p[1] +% fm.offset[1] };
        if (!hinted) p.* = .{ fixed.mul(p[0], glyph_scale), fixed.mul(p[1], glyph_scale) };
    }
    // hmtx's advance, else the charstring's
    var advance: i32 = if (inst.sfnt.hmtx) |hmtx|
        hmtx.getMetrics(glyph).advance_width
    else if (width) |w|
        (w +% priv.nominal_width +% 0x8000) >> 16
    else
        (priv.default_width orelse 0) >> 16;
    if (transformed) advance = fixed.mul(advance, fm.matrix[0]);
    advance = fixed.mul(advance +% fm.offset[0], glyph_scale);
    if (hinted) advance = (advance + 32) & ~@as(i32, 63);
    const points = try allocator.dupe([2]i32, out.points.items);
    errdefer allocator.free(points);
    const tags = try allocator.dupe(u8, out.tags.items);
    errdefer allocator.free(tags);
    return .{ .points = points, .tags = tags, .contours = try allocator.dupe(u16, out.contours.items), .advance = advance };
}

/// One charstring on a path of its own (`cf2_interpT2CharString`), from
/// `start`, its stems counted from `start`'s y.
fn draw(ctx: charstring.Context, code: []const u8, blues: *const Blues, out: *Builder, scale: i32, hinted: bool, start: [2]i32) Error!charstring.Result {
    var path: GlyphPath = .init(blues, out, scale, scale, hinted);
    path.origin_y = start[1];
    const result = try charstring.evaluate(ctx, code, &path, start);
    if (path.err) |e| return e;
    if (path.too_many_hints) return error.TooManyHints;
    return result;
}
