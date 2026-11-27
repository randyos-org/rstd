//! The autohinter for one font: its glyph styles and each style's metrics,
//! and loading a glyph hinted as `FT_LOAD_FORCE_AUTOHINT` loads it.
//! FreeType 2.14.3 `src/autofit/afloader.c` (`af_loader_load_glyph`),
//! `afglobal.c` (`af_face_globals_get_metrics`), `afdummy.c` and
//! `aflatin.c` (`af_latin_hints_init`); stem darkening is off, as FreeType
//! leaves it. The Latin and CJK writing systems; Indic is issue #815.

const std = @import("std");
const fixed = @import("../fixed.zig");
const Sfnt = @import("../../sfnt/Sfnt.zig");
const Glyph = @import("../Glyph.zig");
const GlyphHints = @import("GlyphHints.zig");
const Source = @import("Source.zig");
const coverage = @import("coverage.zig");
const adjust = @import("adjust.zig");
const LatinMetrics = @import("latin/Metrics.zig");
const latin = @import("latin/apply.zig");
const CjkMetrics = @import("cjk/Metrics.zig");
const cjk = @import("cjk/apply.zig");
const styles = @import("data/styles.zig").styles;

const Allocator = std.mem.Allocator;
const Hinter = @This();

pub const RenderMode = latin.RenderMode;
pub const Error = error{ InvalidGlyphIndex, UnsupportedWritingSystem } || Source.Error;

allocator: Allocator,
src: Source,
gstyles: []u16,
reverse: adjust.ReverseMap,
/// Each style's metrics, made on first use.
metrics: [styles.len]?StyleMetrics = @splat(null),
fixed_width: bool,
italic: bool,

pub fn init(allocator: Allocator, sfnt: *const Sfnt) (Allocator.Error || Source.InitError)!Hinter {
    var src = try Source.init(allocator, sfnt);
    errdefer src.deinit();
    const gstyles = try coverage.compute(allocator, sfnt);
    errdefer allocator.free(gstyles);
    const reverse = try adjust.ReverseMap.init(allocator, sfnt);
    const italic = if (sfnt.os2) |os2| os2.fsSelection().oblique or os2.fsSelection().italic else if (sfnt.head) |head| head.macStyle() & 2 != 0 else false;
    return .{
        .allocator = allocator,
        .src = src,
        .gstyles = gstyles,
        .reverse = reverse,
        .fixed_width = if (sfnt.post) |post| post.isFixedPitch() != 0 else false,
        .italic = italic,
    };
}

pub fn deinit(h: *Hinter) void {
    h.reverse.deinit(h.allocator);
    h.allocator.free(h.gstyles);
    h.src.deinit();
}

const StyleMetrics = union(enum) {
    latin: LatinMetrics,
    cjk: CjkMetrics,

    fn scales(s: *const StyleMetrics) [2]i32 {
        return switch (s.*) {
            inline else => |*m| .{ m.axis[0].scale, m.axis[1].scale },
        };
    }

    fn digitsSameWidth(s: *const StyleMetrics) bool {
        return switch (s.*) {
            inline else => |*m| m.digits_have_same_width,
        };
    }
};

/// `af_face_globals_get_metrics`: the glyph's style's metrics, null for
/// the dummy system; a Latin style with no blue zones falls to it.
fn metricsFor(h: *Hinter, gid: u16) Error!?*StyleMetrics {
    while (true) {
        const style = h.gstyles[gid] & coverage.unassigned;
        if (h.metrics[style]) |*m| return m;
        h.metrics[style] = switch (styles[style].writing_system) {
            .dummy => return null,
            .latin => .{ .latin = LatinMetrics.init(h.allocator, &h.src, style, h.gstyles) catch |e| switch (e) {
                error.NoBlueZones => continue,
                else => |x| return x,
            } },
            .cjk => .{ .cjk = try CjkMetrics.init(h.allocator, &h.src, style) },
            .indic => return error.UnsupportedWritingSystem,
        };
    }
}

fn pixRound(x: i32) i32 {
    return (x +% 32) & ~@as(i32, 63);
}

/// Glyph `gid` at `ppem`, autohinted for `mode`; the caller owns it.
pub fn load(h: *Hinter, allocator: Allocator, gid: u16, ppem: u16, mode: RenderMode) Error!Glyph {
    if (gid >= h.gstyles.len) return error.InvalidGlyphIndex;
    const scale = fixed.div(@as(i32, ppem) << 6, h.src.unitsPerEm());
    const metrics = try h.metricsFor(gid);
    const opts: latin.Options = .{ .mode = mode, .ppem = ppem, .italic = h.italic };
    var hints: GlyphHints = .{ .allocator = allocator };
    defer hints.deinit();
    if (metrics) |m| {
        switch (m.*) {
            inline else => |*x| x.scale(scale, scale, 0, 0),
        }
        hints.x_scale, hints.y_scale = m.scales();
    } else {
        hints.x_scale = scale;
        hints.y_scale = scale;
    }
    var g = try h.src.load(allocator, gid);
    errdefer g.deinit(allocator);
    const advance_units = h.src.advance(gid);
    var pp1x: i32 = 0;
    var pp2x = fixed.mul(advance_units, hints.x_scale);
    var lsb_delta: i32 = 0;
    var rsb_delta: i32 = 0;
    if (g.points.len != 0) {
        try hints.reload(.{ .points = g.points, .tags = g.tags, .contours = g.contours }, h.src.unitsPerEm());
        var fit_advance = false;
        if (metrics) |m| switch (m.*) {
            .latin => |*lm| {
                try latin.apply(allocator, &hints, lm, opts, h.gstyles[gid], h.reverse.flags(gid));
                fit_advance = true;
            },
            // CJK sets `AF_SCALER_FLAG_NO_ADVANCE`
            .cjk => |*cm| try cjk.apply(&hints, cm, mode, ppem),
        };
        latin.save(&hints, g.points, g.tags);
        // the side bearings after the hinting moved the outline
        const edges = hints.axis[0].edges.items;
        if (mode != .light and edges.len > 1 and fit_advance) {
            const e1 = edges[0];
            const e2 = edges[edges.len - 1];
            const old_rsb = pp2x - e2.opos;
            const old_lsb = e1.opos;
            const new_lsb = e1.pos;
            var pp1x_uh = new_lsb - old_lsb;
            var pp2x_uh = e2.pos + old_rsb;
            // too much room rather than too little at small sizes
            if (old_lsb < 24) pp1x_uh -= 8;
            if (old_rsb < 24) pp2x_uh += 8;
            pp1x = pixRound(pp1x_uh);
            pp2x = pixRound(pp2x_uh);
            if (pp1x >= new_lsb and old_lsb > 0) pp1x -= 64;
            if (pp2x <= e2.pos and old_rsb > 0) pp2x += 64;
            lsb_delta = pp1x - pp1x_uh;
            rsb_delta = pp2x - pp2x_uh;
        } else {
            const a = pp1x;
            const b = pp2x;
            pp1x = pixRound(a);
            pp2x = pixRound(b);
            lsb_delta = pp1x - a;
            rsb_delta = pp2x - b;
        }
    }
    if (pp1x != 0) for (g.points) |*p| {
        p[0] -= pp1x;
    };
    // a fixed-width face, or a digit of a set of one width, keeps its advance
    const keep = mode != .light and (h.fixed_width or (h.gstyles[gid] & coverage.digit != 0 and (if (metrics) |m| m.digitsSameWidth() else false)));
    var advance: i32 = 0;
    if (keep) {
        advance = fixed.mul(advance_units, hints.x_scale);
        lsb_delta = 0;
        rsb_delta = 0;
    } else if (advance_units != 0) {
        advance = pp2x - pp1x;
    }
    g.advance = pixRound(advance);
    g.lsb_delta = lsb_delta;
    g.rsb_delta = rsb_delta;
    return g;
}
