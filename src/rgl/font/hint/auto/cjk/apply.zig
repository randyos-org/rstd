//! Hinting one glyph with the CJK writing system: FreeType 2.14.3
//! `src/autofit/afcjk.c` (`af_cjk_hints_init`, `af_cjk_hints_apply`).
//! Unlike the Latin system, light mode still hints horizontally, and the
//! advance is never fitted to the edges (`AF_SCALER_FLAG_NO_ADVANCE`).

const std = @import("std");
const GlyphHints = @import("../GlyphHints.zig");
const Metrics = @import("Metrics.zig");
const topo = @import("topo.zig");
const fit = @import("fit.zig");
const align_ = @import("../align.zig");
const latin_fit = @import("../latin/fit.zig");
const RenderMode = @import("../latin/apply.zig").RenderMode;
const Width = @import("../widths.zig").Width;

const Allocator = std.mem.Allocator;

fn fitMode(mode: RenderMode, ppem: u16) latin_fit.Mode {
    return .{
        .horz_snap = mode == .mono or mode == .lcd,
        .vert_snap = mode == .mono or mode == .lcd_v,
        .stem_adjust = mode != .light and mode != .lcd,
        .mono = mode == .mono,
        .ppem = ppem,
    };
}

/// `af_cjk_hints_apply`: hints `h`, loaded at the metrics' scales.
pub fn apply(h: *GlyphHints, m: *const Metrics, mode: RenderMode, ppem: u16) Allocator.Error!void {
    for ([_]GlyphHints.Dim{ .horz, .vert }) |d| {
        try topo.detectFeatures(h, m, d);
        topo.computeBlueEdges(h, m, d);
    }
    const fm = fitMode(mode, ppem);
    for (0..2) |dim| {
        // each zone's reference then undershoot
        var blue_widths: [16]Width = undefined;
        for (m.axis[dim].blueSlice(), 0..) |b, i| {
            blue_widths[2 * i] = b.ref;
            blue_widths[2 * i + 1] = b.shoot;
        }
        fit.hintEdges(h, m, fm, dim, &blue_widths);
        fit.alignEdgePoints(h, fm, dim);
        align_.strongPoints(h, dim);
        align_.weakPoints(h, dim);
    }
}
