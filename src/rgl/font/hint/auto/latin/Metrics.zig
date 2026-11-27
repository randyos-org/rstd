//! A Latin-system style's global metrics: its standard stem widths and
//! blue zones in font units, scaled per size. FreeType 2.14.3
//! `src/autofit/aflatin.h` and `aflatin.c` (`af_latin_metrics_init`,
//! `af_latin_metrics_check_digits`, `af_latin_metrics_scale_dim`,
//! `af_latin_metrics_scale`), without HarfBuzz; after fontations
//! `skrifa/src/outline/autohint/metrics/` (commit 4f115482, MIT).

const std = @import("std");
const fixed = @import("../../fixed.zig");
const Width = @import("../widths.zig").Width;
const Source = @import("../Source.zig");
const shaper = @import("../shaper.zig");
const coverage = @import("../coverage.zig");
const widths = @import("widths.zig");
const blues = @import("blues.zig");

const Allocator = std.mem.Allocator;
const Metrics = @This();

pub const max_widths = 16;
/// `AF_BLUE_STRINGSET_MAX_LEN`.
pub const max_blues = 8;

/// `AF_LATIN_BLUE_*`.
pub const blue_flag = struct {
    /// At most 3/4 pixel tall at this size.
    pub const active: u32 = 1 << 0;
    pub const top: u32 = 1 << 1;
    pub const sub_top: u32 = 1 << 2;
    pub const neutral: u32 = 1 << 3;
    /// The x height, which tunes the vertical scale.
    pub const adjustment: u32 = 1 << 4;
    pub const bottom: u32 = 1 << 5;
    pub const bottom_small: u32 = 1 << 6;
};

pub const Blue = struct {
    ref: Width = .{},
    shoot: Width = .{},
    ascender: i32 = 0,
    descender: i32 = 0,
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
    /// Vertical only.
    blue_count: usize = 0,
    blues: [max_blues]Blue = @splat(.{}),
    org_scale: i32 = 0,
    org_delta: i32 = 0,

    pub fn widthSlice(a: *const Axis) []const Width {
        return a.widths[0..a.width_count];
    }

    pub fn blueSlice(a: *const Axis) []const Blue {
        return a.blues[0..a.blue_count];
    }
};

pub const InitError = error{NoBlueZones} || Source.Error;

style: u16,
units_per_em: u16,
axis: [2]Axis = .{ .{}, .{} },
digits_have_same_width: bool = true,

/// `af_latin_metrics_init`. With no blue zones, the style's glyphs in
/// `gstyles` go unhinted and `error.NoBlueZones` returns, as FreeType's -1.
pub fn init(allocator: Allocator, src: *Source, style: u16, gstyles: []u16) InitError!Metrics {
    var m: Metrics = .{ .style = style, .units_per_em = src.unitsPerEm() };
    if (src.sfnt.unicodeCmap() == null) return m;
    try widths.init(allocator, &m, src);
    if (!try blues.init(allocator, &m, src)) {
        for (gstyles) |*s| {
            if (s.* & coverage.style_mask == style) s.* = coverage.none_style;
        }
        return error.NoBlueZones;
    }
    m.checkDigits(src);
    return m;
}

/// `af_latin_metrics_check_digits`: whether 0 to 9 share one advance.
fn checkDigits(m: *Metrics, src: *const Source) void {
    const digits = "0 1 2 3 4 5 6 7 8 9";
    var at: usize = 0;
    var started = false;
    var old_advance: i32 = 0;
    m.digits_have_same_width = true;
    while (at < digits.len) {
        const c = shaper.cluster(src.sfnt, digits, &at);
        if (c.count > 1 or c.glyph == 0) continue;
        const advance = src.advance(c.glyph);
        if (started) {
            if (advance != old_advance) {
                m.digits_have_same_width = false;
                break;
            }
        } else {
            old_advance = advance;
            started = true;
        }
    }
}

/// `af_latin_metrics_scale_dim`: the scale (the vertical one nudged so the
/// x height fits the grid), the scaled widths and blue zones.
fn scaleDim(m: *Metrics, dim: usize, scale_in: i32, delta: i32) void {
    var sc = scale_in;
    const axis = &m.axis[dim];
    if (axis.org_scale == sc and axis.org_delta == delta) return;
    axis.org_scale = sc;
    axis.org_delta = delta;
    const vert = &m.axis[1];
    for (vert.blueSlice()) |blue| {
        if (blue.flags & blue_flag.adjustment == 0) continue;
        const scaled = fixed.mul(blue.shoot.org, sc);
        // `increase-x-height` is off by default: the threshold is 40
        const fitted = (scaled + 40) & ~@as(i32, 63);
        if (scaled != fitted and dim == 1) {
            const new_scale = fixed.mulDiv(sc, fitted, scaled);
            // at most two pixels' change
            var max_height: i32 = m.units_per_em;
            for (vert.blueSlice()) |b| {
                max_height = @max(max_height, b.ascender);
                max_height = @max(max_height, -b.descender);
            }
            const dist = fixed.mul(max_height, new_scale - sc);
            if (-128 < dist and dist < 128) sc = new_scale;
        }
        break;
    }
    axis.scale = sc;
    axis.delta = delta;
    for (axis.widths[0..axis.width_count]) |*w| {
        w.cur = fixed.mul(w.org, sc);
        w.fit = w.cur;
    }
    // under 5/8 pixel
    axis.extra_light = fixed.mul(axis.standard_width, sc) < 32 + 8;
    if (dim != 1) return;
    for (axis.blues[0..axis.blue_count]) |*blue| {
        blue.ref.cur = fixed.mul(blue.ref.org, sc) + delta;
        blue.ref.fit = blue.ref.cur;
        blue.shoot.cur = fixed.mul(blue.shoot.org, sc) + delta;
        blue.shoot.fit = blue.shoot.cur;
        blue.flags &= ~blue_flag.active;
        const dist = fixed.mul(blue.ref.org - blue.shoot.org, sc);
        if (dist <= 48 and dist >= -48) {
            // discrete zone heights: 0, half a pixel or one
            const mag: i32 = @intCast(@abs(dist));
            var delta2: i32 = if (mag < 32) 0 else if (mag < 48) 32 else 64;
            if (dist < 0) delta2 = -delta2;
            blue.ref.fit = fixed.roundPad(blue.ref.cur, 64);
            blue.shoot.fit = blue.ref.fit - delta2;
            blue.flags |= blue_flag.active;
        }
    }
    // a sub-top zone overlapping another active one is dropped
    for (axis.blues[0..axis.blue_count]) |*blue| {
        if (blue.flags & blue_flag.sub_top == 0 or blue.flags & blue_flag.active == 0) continue;
        for (axis.blues[0..axis.blue_count]) |b| {
            if (b.flags & blue_flag.sub_top != 0 or b.flags & blue_flag.active == 0) continue;
            if (b.ref.fit <= blue.shoot.fit and b.shoot.fit >= blue.ref.fit) {
                blue.flags &= ~blue_flag.active;
                break;
            }
        }
    }
}

/// `af_latin_metrics_scale`.
pub fn scale(m: *Metrics, x_scale: i32, y_scale: i32, x_delta: i32, y_delta: i32) void {
    m.scaleDim(0, x_scale, x_delta);
    m.scaleDim(1, y_scale, y_delta);
}
