//! Hinting one glyph with the Latin writing system: FreeType 2.14.3
//! `src/autofit/aflatin.c` (`af_latin_hints_compute_blue_edges`,
//! `af_latin_hints_init`, `af_latin_hints_apply`) and `afhints.c`
//! (`af_glyph_hints_save`).

const std = @import("std");
const fixed = @import("../../fixed.zig");
const GlyphHints = @import("../GlyphHints.zig");
const Metrics = @import("Metrics.zig");
const edges = @import("../edges.zig");
const align_ = @import("../align.zig");
const adjust = @import("../adjust.zig");
const coverage = @import("../coverage.zig");
const fit = @import("fit.zig");
const diacritics = @import("diacritics.zig");
const Width = @import("../widths.zig").Width;
const scripts = @import("../data/scripts.zig").scripts;
const styles = @import("../data/styles.zig").styles;

const Allocator = std.mem.Allocator;
const blue_flag = Metrics.blue_flag;
const none = GlyphHints.none;

pub const RenderMode = enum { normal, light, mono, lcd, lcd_v };

/// What `af_latin_hints_init` derives from the render mode.
pub const Options = struct {
    mode: RenderMode,
    ppem: u16,
    italic: bool = false,

    fn fitMode(o: Options) fit.Mode {
        return .{
            .horz_snap = o.mode == .mono or o.mode == .lcd,
            .vert_snap = o.mode == .mono or o.mode == .lcd_v,
            .stem_adjust = o.mode != .light and o.mode != .lcd,
            .mono = o.mode == .mono,
            .ppem = o.ppem,
        };
    }

    /// Light and LCD modes, and italic faces, are hinted vertically only.
    pub fn horizontal(o: Options) bool {
        return !(o.mode == .light or o.mode == .lcd or o.italic);
    }
};

/// `af_latin_hints_compute_blue_edges`: each horizontal edge's nearest blue
/// zone edge within half a pixel. A blue edge is an index, twice the zone's
/// plus 1 for its overshoot.
fn computeBlueEdges(h: *GlyphHints, m: *const Metrics) void {
    const axis = &h.axis[1];
    const latin = &m.axis[1];
    const scale = latin.scale;
    for (axis.edges.items) |*edge| {
        if (edge.flags & GlyphHints.edge_flag.no_blue != 0) continue;
        var best_dist = @min(fixed.mul(@divTrunc(@as(i32, m.units_per_em), 40), scale), 64 / 2);
        var best: u32 = none;
        var best_neutral = false;
        for (latin.blueSlice(), 0..) |blue, bb| {
            if (blue.flags & blue_flag.active == 0) continue;
            const is_top = blue.flags & (blue_flag.top | blue_flag.sub_top) != 0;
            const is_neutral = blue.flags & blue_flag.neutral != 0;
            const is_major = edge.dir == axis.major_dir;
            // top zones take edges against the major direction, bottom ones with it
            if (is_top == is_major and !is_neutral) continue;
            var dist = fixed.mul(@intCast(@abs(edge.fpos - blue.ref.org)), scale);
            if (dist < best_dist) {
                best_dist = dist;
                best = @intCast(2 * bb);
                best_neutral = is_neutral;
            }
            // a round edge beyond the reference may take the overshoot
            if (edge.flags & GlyphHints.edge_flag.round != 0 and dist != 0 and !is_neutral) {
                const is_under_ref = edge.fpos < blue.ref.org;
                if (is_top != is_under_ref) {
                    dist = fixed.mul(@intCast(@abs(edge.fpos - blue.shoot.org)), scale);
                    if (dist < best_dist) {
                        best_dist = dist;
                        best = @intCast(2 * bb + 1);
                        best_neutral = is_neutral;
                    }
                }
            }
        }
        if (best != none) {
            edge.blue_edge = best;
            if (best_neutral) edge.flags |= GlyphHints.edge_flag.neutral;
        }
    }
}

/// `af_latin_hints_apply`: hints `h`, already loaded from the outline
/// (`af_glyph_hints_reload`) at the metrics' scales.
pub fn apply(allocator: Allocator, h: *GlyphHints, m: *const Metrics, opts: Options, glyph_style: u16, adjust_flags: u32) Allocator.Error!void {
    const style = styles[m.style];
    const top_to_bottom = scripts[style.script].top_to_bottom;
    const horizontal = opts.horizontal();
    if (horizontal) {
        const a = &m.axis[0];
        try edges.detectFeatures(h, .horz, maxWidth(a), a.edge_distance_threshold, top_to_bottom, m.units_per_em);
    }
    const nc = h.contours.items.len;
    const ext_buf = try allocator.alloc(i64, 2 * nc);
    defer allocator.free(ext_buf);
    const ext: diacritics.Extrema = .{ .min = ext_buf[0..nc], .max = ext_buf[nc..] };
    var accent_height_limit: i64 = 0;

    // vertical
    {
        const any_tilde = adjust_flags & (adjust.tilde_top | adjust.tilde_bottom | adjust.tilde_top2 | adjust.tilde_bottom2) != 0;
        if (any_tilde) ext.compute(h);
        diacritics.stretchTildes(h, ext, adjust_flags);
        const a = &m.axis[1];
        try edges.detectFeatures(h, .vert, maxWidth(a), a.edge_distance_threshold, top_to_bottom, m.units_per_em);
        const capital = diacritics.baseBlues(m, true);
        const small = diacritics.baseBlues(m, false);
        if (adjust_flags != 0 and adjust_flags & adjust.no_height_check == 0) {
            // the accent at most 2/3 of the small letters, or half the
            // capitals, or half the PostScript ascender (8/10 em)
            accent_height_limit = if (small[0] != null and small[1] != null)
                @divTrunc(2 * (small[0].?.shoot.cur - small[1].?.shoot.cur), 3)
            else if (capital[0] != null and capital[1] != null)
                @divTrunc(capital[0].?.shoot.cur - capital[1].?.shoot.cur, 2)
            else
                fixed.mul(@divTrunc(@as(i32, m.units_per_em) * 4, 10), m.axis[1].scale);
        }
        if (capital[0] != null and capital[1] != null) {
            if (adjust_flags & adjust.ignore_capital_top != 0) diacritics.ignoreTop(h, capital[0].?, capital[1].?);
            if (adjust_flags & adjust.ignore_capital_bottom != 0) diacritics.ignoreBottom(h, capital[0].?, capital[1].?);
        }
        if (small[0] != null and small[1] != null) {
            if (adjust_flags & adjust.ignore_small_top != 0) diacritics.ignoreTop(h, small[0].?, small[1].?);
            if (adjust_flags & adjust.ignore_small_bottom != 0) diacritics.ignoreBottom(h, small[0].?, small[1].?);
        }
        // blue zones apply to base characters only
        if (glyph_style & coverage.nonbase == 0) computeBlueEdges(h, m);
    }

    // the zones' edges for the grid fitting: ref then overshoot
    var blue_widths: [2 * Metrics.max_blues]Width = undefined;
    for (m.axis[1].blueSlice(), 0..) |b, i| {
        blue_widths[2 * i] = b.ref;
        blue_widths[2 * i + 1] = b.shoot;
    }
    const fm = opts.fitMode();
    for (0..2) |dim| {
        if (dim == 0 and !horizontal) continue;
        fit.hintEdges(h, m, fm, dim, top_to_bottom, &blue_widths);
        align_.edgePoints(h, dim);
        align_.strongPoints(h, dim);
        align_.weakPoints(h, dim);
        if (dim == 1) diacritics.separate(h, ext, adjust_flags, accent_height_limit);
    }
}

/// The widest standard width, for scoring stem links.
fn maxWidth(a: *const Metrics.Axis) i32 {
    return if (a.width_count > 0) a.widths[a.width_count - 1].org else 0;
}

/// `af_glyph_hints_save`: the hinted points and their tags.
pub fn save(h: *const GlyphHints, points: [][2]i32, tags: []u8) void {
    for (h.points.items, points, tags) |p, *v, *t| {
        v.* = .{ p.x, p.y };
        t.* = if (p.flags & GlyphHints.flag.conic != 0) 0 else if (p.flags & GlyphHints.flag.cubic != 0) 2 else 1;
    }
}
