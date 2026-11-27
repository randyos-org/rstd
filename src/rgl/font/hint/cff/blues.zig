//! A subfont's blue zones at one scale, and the capture of hint edges by
//! them: FreeType 2.14.3 `src/psaux/psblues.c` (`cf2_blues_init`,
//! `cf2_blues_capture`), with stem darkening off as FreeType's CFF driver
//! leaves it; after fontations `skrifa/src/outline/cff/hint.rs`
//! (`HintState`, commit 4f115482, MIT).

const std = @import("std");
const fixed = @import("../fixed.zig");
const private = @import("private.zig");
const Hint = @import("hintmap.zig").Hint;

/// The em box FreeType assumes for ideographic fonts (`CF2_ICF_Top`).
pub const icf_top: i32 = 880 << 16;
pub const icf_bottom: i32 = -120 << 16;
/// Room kept beside a hint edge (`CF2_MIN_COUNTER`).
pub const min_counter: i32 = 0x8000;
const max_zones = 12;

const Zone = struct {
    bottom: bool,
    cs_bottom: i32,
    cs_top: i32,
    cs_flat: i32,
    ds_flat: i32,
};

pub const Blues = struct {
    scale: i32,
    blue_shift: i32,
    blue_fuzz: i32,
    suppress_overshoot: bool = false,
    em_box: bool = false,
    zones: [max_zones]Zone = undefined,
    count: usize = 0,

    pub fn init(p: *const private.Private, scale: i32) Blues {
        var b: Blues = .{ .scale = scale, .blue_shift = p.blue_shift, .blue_fuzz = p.blue_fuzz };
        const blues = p.blues.slice();
        if (p.language_group == 1 and (blues.len == 0 or (blues.len == 2 and
            blues[0][0] < icf_bottom and blues[0][1] < icf_bottom and blues[1][0] > icf_top and blues[1][1] > icf_top)))
        {
            b.em_box = true;
            return b;
        }
        var max_height: i32 = 0;
        for (blues, 0..) |z, i| {
            const height = z[1] -% z[0];
            if (height < 0) continue;
            max_height = @max(max_height, height);
            // the first BlueValues pair is the bottom zone, the rest top ones
            b.zones[b.count] = .{ .bottom = i == 0, .cs_bottom = z[0], .cs_top = z[1], .cs_flat = if (i == 0) z[1] else z[0], .ds_flat = 0 };
            b.count += 1;
        }
        for (p.other_blues.slice()) |z| {
            const height = z[1] -% z[0];
            if (height < 0) continue;
            max_height = @max(max_height, height);
            b.zones[b.count] = .{ .bottom = true, .cs_bottom = z[0], .cs_top = z[1], .cs_flat = z[1], .ds_flat = 0 };
            b.count += 1;
        }
        // a FamilyBlues edge within a pixel replaces the flat edge
        const units_per_pixel = fixed.div(1 << 16, scale);
        const family = p.family_blues.slice();
        for (b.zones[0..b.count]) |*z| {
            const flat = z.cs_flat;
            var min_diff: i32 = std.math.maxInt(i32);
            if (z.bottom) {
                for (p.family_other_blues.slice()) |f| {
                    const diff: i32 = @intCast(@abs(flat -% f[1]));
                    if (diff < min_diff and diff < units_per_pixel) {
                        z.cs_flat = f[1];
                        min_diff = diff;
                        if (diff == 0) break;
                    }
                }
                if (family.len >= 1) {
                    const diff: i32 = @intCast(@abs(flat -% family[0][1]));
                    if (diff < min_diff and diff < units_per_pixel) z.cs_flat = family[0][1];
                }
            } else if (family.len > 1) for (family[1..]) |f| {
                const diff: i32 = @intCast(@abs(flat -% f[0]));
                if (diff < min_diff and diff < units_per_pixel) {
                    z.cs_flat = f[0];
                    min_diff = diff;
                    if (diff == 0) break;
                }
            };
        }
        var blue_scale = p.blue_scale;
        if (max_height > 0 and blue_scale > fixed.div(1 << 16, max_height)) blue_scale = fixed.div(1 << 16, max_height);
        // below BlueScale, overshoot is suppressed and the zones boosted
        var boost: i32 = 0;
        if (scale < blue_scale) {
            b.suppress_overshoot = true;
            const six_tenths: i32 = 39322; // cf2_doubleToFixed(0.6)
            boost = @min(six_tenths - fixed.mulDiv(six_tenths, scale, blue_scale), 0x7FFF);
        }
        for (b.zones[0..b.count]) |*z| {
            const v = fixed.mul(z.cs_flat, scale) +% (if (z.bottom) -boost else boost);
            z.ds_flat = round(v);
        }
        return b;
    }

    /// The ghost edges at the em box's bottom and top (`emBoxBottomEdge`).
    pub fn emBoxEdges(b: *const Blues) [2]Hint {
        const bottom_cs = icf_bottom - 1;
        const top_cs = icf_top + 1;
        return .{
            .{ .flags = Hint.ghost_bottom | Hint.locked | Hint.synthetic, .cs = bottom_cs, .ds = round(fixed.mul(bottom_cs, b.scale)) - min_counter, .scale = b.scale },
            .{ .flags = Hint.ghost_top | Hint.locked | Hint.synthetic, .cs = top_cs, .ds = round(fixed.mul(top_cs, b.scale)) + min_counter, .scale = b.scale },
        };
    }

    /// Locks `bottom` and `top` to the first zone that captures one; true
    /// if one did.
    pub fn capture(b: *const Blues, bottom: *Hint, top: *Hint) bool {
        var adjustment: i32 = 0;
        var captured = false;
        for (b.zones[0..b.count]) |z| {
            if (z.bottom and bottom.isBottom() and z.cs_bottom -% b.blue_fuzz <= bottom.cs and bottom.cs <= z.cs_top +% b.blue_fuzz) {
                adjustment = if (b.suppress_overshoot)
                    z.ds_flat
                else if (z.cs_top -% bottom.cs >= b.blue_shift)
                    // at least a pixel of overshoot
                    @min(round(bottom.ds), z.ds_flat - (1 << 16))
                else
                    round(bottom.ds);
                adjustment -%= bottom.ds;
                captured = true;
                break;
            }
            if (!z.bottom and top.isTop() and z.cs_bottom -% b.blue_fuzz <= top.cs and top.cs <= z.cs_top +% b.blue_fuzz) {
                adjustment = if (b.suppress_overshoot)
                    z.ds_flat
                else if (top.cs -% z.cs_bottom >= b.blue_shift)
                    @max(round(top.ds), z.ds_flat + (1 << 16))
                else
                    round(top.ds);
                adjustment -%= top.ds;
                captured = true;
                break;
            }
        }
        if (captured) {
            if (bottom.flags != 0) {
                bottom.ds +%= adjustment;
                bottom.flags |= Hint.locked;
            }
            if (top.flags != 0) {
                top.ds +%= adjustment;
                top.flags |= Hint.locked;
            }
        }
        return captured;
    }
};

/// `cf2_fixedRound`.
pub fn round(x: i32) i32 {
    return (x +% 0x8000) & ~@as(i32, 0xFFFF);
}
