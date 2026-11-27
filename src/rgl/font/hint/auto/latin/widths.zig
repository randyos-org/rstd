//! A Latin-system style's standard stem widths, from the first of its
//! standard characters the font has: FreeType 2.14.3
//! `src/autofit/aflatin.c` (`af_latin_metrics_init_widths`).

const std = @import("std");
const Metrics = @import("Metrics.zig");
const Source = @import("../Source.zig");
const GlyphHints = @import("../GlyphHints.zig");
const segments = @import("../segments.zig");
const shaper = @import("../shaper.zig");
const sortAndQuantize = @import("../widths.zig").sortAndQuantize;
const scripts = @import("../data/scripts.zig").scripts;
const styles = @import("../data/styles.zig").styles;

const Allocator = std.mem.Allocator;

pub fn init(allocator: Allocator, m: *Metrics, src: *Source) Source.Error!void {
    m.axis[0].width_count = 0;
    m.axis[1].width_count = 0;
    defer finish(m);
    const chars = scripts[styles[m.style].script].standard_chars;
    var at: usize = 0;
    var glyph: u16 = 0;
    while (at < chars.len) {
        const c = shaper.cluster(src.sfnt, chars, &at);
        if (c.count > 1) continue;
        glyph = c.glyph;
        if (glyph != 0) break;
    }
    // no standard character: the fallback widths
    if (glyph == 0) return;
    var g = src.load(allocator, glyph) catch |e| switch (e) {
        error.OutOfMemory => return e,
        else => return,
    };
    defer g.deinit(allocator);
    if (g.points.len == 0) return;
    var h: GlyphHints = .{ .allocator = allocator };
    defer h.deinit();
    try h.reload(.{ .points = g.points, .tags = g.tags, .contours = g.contours }, m.units_per_em);
    for ([_]GlyphHints.Dim{ .horz, .vert }) |d| {
        const dim = @backingInt(d);
        try segments.compute(&h, d, m.units_per_em);
        segments.link(&h, d, 0, m.units_per_em);
        const axis = &m.axis[dim];
        var n: usize = 0;
        const segs = h.axis[dim].segments.items;
        for (segs, 0..) |seg, i| {
            // stems only: mutual links, each counted once
            if (seg.link == GlyphHints.none or segs[seg.link].link != i or seg.link <= i) continue;
            if (n < Metrics.max_widths) {
                axis.widths[n].org = @intCast(@abs(seg.pos - segs[seg.link].pos));
                n += 1;
            }
        }
        // near-equal widths merge (100 is FreeType's heuristic)
        axis.width_count = sortAndQuantize(&axis.widths, n, @divTrunc(@as(i32, m.units_per_em), 100));
    }
}

/// The standard width, the first one or the fallback, and the edge distance
/// threshold, a fifth of it.
fn finish(m: *Metrics) void {
    for (&m.axis) |*axis| {
        const stdw: i32 = if (axis.width_count > 0) axis.widths[0].org else segments.latinConstant(m.units_per_em, 50);
        axis.edge_distance_threshold = @divTrunc(stdw, 5);
        axis.standard_width = stdw;
        axis.extra_light = false;
    }
}
