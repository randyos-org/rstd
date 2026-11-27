//! A CJK-system style's global metrics: standard stem widths, and blue
//! zones on both axes (tops and bottoms, lefts and rights). FreeType 2.14.3
//! `src/autofit/afcjk.h` and `afcjk.c` (`af_cjk_metrics_init`,
//! `af_cjk_metrics_init_blues`, `af_cjk_metrics_scale_dim`,
//! `af_cjk_metrics_scale`), without HarfBuzz. As FreeType's, its widths
//! are never scaled: their `cur` stays 0.

const std = @import("std");
const fixed = @import("../../fixed.zig");
const Width = @import("../widths.zig").Width;
const sortPos = @import("../widths.zig").sortPos;
const Source = @import("../Source.zig");
const shaper = @import("../shaper.zig");
const stdwidths = @import("../stdwidths.zig");
const digits = @import("../digits.zig");
const blue_strings = @import("../data/blues.zig").blue_strings;
const styles = @import("../data/styles.zig").styles;

const Allocator = std.mem.Allocator;
const Metrics = @This();

pub const max_widths = 16;
const max_blues = 8;
/// `AF_BLUE_STRING_MAX_LEN`.
const max_len = 51;

/// `AF_CJK_BLUE_*`.
pub const blue_flag = struct {
    pub const active: u32 = 1 << 0;
    /// Top, or right on the horizontal axis.
    pub const top: u32 = 1 << 1;
};

/// `AF_BLUE_PROPERTY_CJK_*`.
const prop = struct {
    const top: u16 = 1 << 0;
    const horiz: u16 = 1 << 1;
};

pub const Blue = struct {
    ref: Width = .{},
    /// The undershoot.
    shoot: Width = .{},
    flags: u32 = 0,
};

pub const Axis = struct {
    scale: i32 = 0,
    delta: i32 = 0,
    width_count: usize = 0,
    widths: [max_widths]Width = @splat(.{}),
    edge_distance_threshold: i32 = 0,
    standard_width: i32 = 0,
    extra_light: bool = false,
    blue_count: usize = 0,
    blues: [max_blues]Blue = @splat(.{}),
    org_scale: i32 = 0,
    org_delta: i32 = 0,

    pub fn blueSlice(a: *const Axis) []const Blue {
        return a.blues[0..a.blue_count];
    }
};

style: u16,
units_per_em: u16,
axis: [2]Axis = .{ .{}, .{} },
digits_have_same_width: bool = true,

/// `af_cjk_metrics_init`; unlike Latin's, zones or none, it succeeds.
pub fn init(allocator: Allocator, src: *Source, style: u16) Source.Error!Metrics {
    var m: Metrics = .{ .style = style, .units_per_em = src.unitsPerEm() };
    if (src.sfnt.unicodeCmap() == null) return m;
    try stdwidths.init(allocator, &m, src);
    try m.initBlues(allocator, src);
    m.digits_have_same_width = digits.sameWidth(src);
    return m;
}

/// `af_cjk_metrics_init_blues`: per string, the medians of the fill
/// extremes (before `|`) and the flat ones (after).
fn initBlues(m: *Metrics, allocator: Allocator, src: *Source) Source.Error!void {
    var bs: usize = styles[m.style].blue_stringset;
    while (blue_strings[bs].chars.len != 0) : (bs += 1) {
        const props = blue_strings[bs].properties;
        const horiz = props & prop.horiz != 0;
        const top = props & prop.top != 0;
        const axis = &m.axis[if (horiz) 0 else 1];
        var fills: [max_len]i32 = undefined;
        var flats: [max_len]i32 = undefined;
        var num_fills: usize = 0;
        var num_flats: usize = 0;
        var fill = true;
        const chars = blue_strings[bs].chars;
        var at: usize = 0;
        while (at < chars.len) {
            while (at < chars.len and chars[at] == ' ') at += 1;
            if (at < chars.len and chars[at] == '|') {
                fill = false;
                at += 1;
                continue;
            }
            const c = shaper.cluster(src.sfnt, chars, &at);
            if (c.count > 1 or c.glyph == 0) continue;
            var g = src.load(allocator, c.glyph) catch |e| switch (e) {
                error.OutOfMemory => return e,
                else => continue,
            };
            defer g.deinit(allocator);
            if (g.points.len <= 2) continue;
            // the extreme point, single points left out
            var best_point: i32 = -1;
            var best_pos: i32 = 0;
            var last: i32 = -1;
            for (g.contours) |end| {
                const first = last + 1;
                last = end;
                if (last <= first) continue;
                var pp = first;
                while (pp <= last) : (pp += 1) {
                    const v = g.points[@intCast(pp)][if (horiz) 0 else 1];
                    if (best_point < 0 or (if (top) v > best_pos else v < best_pos)) {
                        best_point = pp;
                        best_pos = v;
                    }
                }
            }
            if (fill) {
                fills[num_fills] = best_pos;
                num_fills += 1;
            } else {
                flats[num_flats] = best_pos;
                num_flats += 1;
            }
        }
        if (num_flats == 0 and num_fills == 0) continue;
        sortPos(fills[0..num_fills]);
        sortPos(flats[0..num_flats]);
        const blue = &axis.blues[axis.blue_count];
        axis.blue_count += 1;
        if (num_flats == 0) {
            blue.ref.org = fills[num_fills / 2];
            blue.shoot.org = blue.ref.org;
        } else if (num_fills == 0) {
            blue.ref.org = flats[num_flats / 2];
            blue.shoot.org = blue.ref.org;
        } else {
            blue.ref.org = fills[num_fills / 2];
            blue.shoot.org = flats[num_flats / 2];
        }
        // the reference must lie beyond the undershoot: else their mean
        if (blue.shoot.org != blue.ref.org) {
            const ref = blue.ref.org;
            const shoot = blue.shoot.org;
            if (top != (shoot < ref)) {
                blue.ref.org = @divTrunc(shoot + ref, 2);
                blue.shoot.org = blue.ref.org;
            }
        }
        blue.flags = if (top) blue_flag.top else 0;
    }
}

/// `af_cjk_metrics_scale_dim`: the blue zones scaled; a zone under 3/4
/// pixel tall is active, its undershoot a whole pixel or none away.
fn scaleDim(m: *Metrics, dim: usize, sc: i32, delta: i32) void {
    const axis = &m.axis[dim];
    if (axis.org_scale == sc and axis.org_delta == delta) return;
    axis.org_scale = sc;
    axis.org_delta = delta;
    axis.scale = sc;
    axis.delta = delta;
    for (axis.blues[0..axis.blue_count]) |*blue| {
        blue.ref.cur = fixed.mul(blue.ref.org, sc) + delta;
        blue.ref.fit = blue.ref.cur;
        blue.shoot.cur = fixed.mul(blue.shoot.org, sc) + delta;
        blue.shoot.fit = blue.shoot.cur;
        blue.flags &= ~blue_flag.active;
        const dist = fixed.mul(blue.ref.org - blue.shoot.org, sc);
        if (dist <= 48 and dist >= -48) {
            blue.ref.fit = (blue.ref.cur + 32) & ~@as(i32, 63);
            const delta1 = fixed.div(blue.ref.fit, sc) - blue.shoot.org;
            var delta2 = fixed.mul(@intCast(@abs(delta1)), sc);
            delta2 = if (delta2 < 32) 0 else (delta2 + 32) & ~@as(i32, 63);
            if (delta1 < 0) delta2 = -delta2;
            blue.shoot.fit = blue.ref.fit - delta2;
            blue.flags |= blue_flag.active;
        }
    }
}

/// `af_cjk_metrics_scale`: the scales as given, no x-height tuning.
pub fn scale(m: *Metrics, x_scale: i32, y_scale: i32, x_delta: i32, y_delta: i32) void {
    m.scaleDim(0, x_scale, x_delta);
    m.scaleDim(1, y_scale, y_delta);
}
