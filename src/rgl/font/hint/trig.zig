//! FreeType's CORDIC trigonometry for the hinters: `FT_Vector_Length` (the
//! loader's `FT_Hypot`) and `FT_Atan2` (the autohinter's), from FreeType
//! 2.14.3 `src/base/fttrigon.c` (`ft_trig_downscale`'s 64-bit form,
//! `ft_trig_prenorm`, `ft_trig_pseudo_polarize`). Angles are 16.16 degrees.

/// 180 degrees (`FT_ANGLE_PI`).
pub const pi: i32 = 180 << 16;
pub const pi2: i32 = 90 << 16;
pub const pi4: i32 = 45 << 16;

/// The CORDIC gain's reciprocal, 0.858785336480436 in 0.32.
const trig_scale: u32 = 0xDBD95B16;
/// The highest bit a normalized component may have without overflow.
const safe_msb = 29;
const max_iters = 23;

/// `atan(2^-i)` in 16.16 degrees, `i` from 0.
const arctan_table = [_]i32{ 1740967, 919879, 466945, 234379, 117304, 58666, 29335, 14668, 7334, 3667, 1833, 917, 458, 229, 115, 57, 29, 14, 7, 4, 2, 1 };

/// `val` times the CORDIC shrink factor, rounded.
fn downscale(val: i32) i32 {
    const neg = val < 0;
    const v: u32 = @abs(val);
    const r: u64 = (@as(u64, v) * trig_scale + 0x40000000) >> 32;
    const out: i32 = @intCast(r);
    return if (neg) -out else out;
}

/// Scales `v` so its larger component's top bit is bit 29; the shift taken.
fn prenorm(v: *[2]i32) i32 {
    const m: u32 = @as(u32, @abs(v[0])) | @as(u32, @abs(v[1]));
    const msb: i32 = 31 - @as(i32, @clz(m));
    if (msb <= safe_msb) {
        const s: u5 = @intCast(safe_msb - msb);
        v[0] = @bitCast(@as(u32, @bitCast(v[0])) << s);
        v[1] = @bitCast(@as(u32, @bitCast(v[1])) << s);
        return safe_msb - msb;
    }
    const s: u5 = @intCast(msb - safe_msb);
    v[0] >>= s;
    v[1] >>= s;
    return -(msb - safe_msb);
}

/// Rotates `v` onto the x axis: `v[0]` its length times the CORDIC gain,
/// `v[1]` its angle.
fn pseudoPolarize(v: *[2]i32) void {
    var x = v[0];
    var y = v[1];
    var theta: i32 = 0;
    if (y > x) {
        if (y > -x) {
            theta = pi2;
            const t = y;
            y = -x;
            x = t;
        } else {
            theta = if (y > 0) pi else -pi;
            x = -x;
            y = -y;
        }
    } else if (y < -x) {
        theta = -pi2;
        const t = -y;
        y = x;
        x = t;
    }
    var b: i32 = 1;
    for (1..max_iters) |i| {
        const s: u5 = @intCast(i);
        const xt = if (y > 0) x + ((y + b) >> s) else x - ((y + b) >> s);
        if (y > 0) {
            y = y - ((x + b) >> s);
            theta += arctan_table[i - 1];
        } else {
            y = y + ((x + b) >> s);
            theta -= arctan_table[i - 1];
        }
        x = xt;
        b <<= 1;
    }
    theta = if (theta >= 0) (theta + 8) & ~@as(i32, 15) else -((-theta + 8) & ~@as(i32, 15));
    v[0] = x;
    v[1] = theta;
}

/// The length of `(x, y)` (`FT_Vector_Length`, `FT_Hypot`).
pub fn vectorLength(x: i32, y: i32) i32 {
    if (x == 0) return @intCast(@abs(y));
    if (y == 0) return @intCast(@abs(x));
    var v = [2]i32{ x, y };
    const shift = prenorm(&v);
    pseudoPolarize(&v);
    const len = downscale(v[0]);
    if (shift > 0) return (len + (@as(i32, 1) << @intCast(shift - 1))) >> @intCast(shift);
    return @bitCast(@as(u32, @bitCast(len)) << @intCast(-shift));
}

/// The angle of `(dx, dy)` in 16.16 degrees, 0 for the zero vector
/// (`FT_Atan2`).
pub fn atan2(dx: i32, dy: i32) i32 {
    if (dx == 0 and dy == 0) return 0;
    var v = [2]i32{ dx, dy };
    _ = prenorm(&v);
    pseudoPolarize(&v);
    return v[1];
}

test "trig: lengths and angles near their true values" {
    const std = @import("std");
    try std.testing.expect(@abs(vectorLength(3 << 16, 4 << 16) - (5 << 16)) <= 1);
    try std.testing.expectEqual(@as(i32, 7), vectorLength(0, -7));
    try std.testing.expect(@abs(vectorLength(0x4000, 0x4000) - 23170) <= 1);
    try std.testing.expect(@abs(atan2(1 << 16, 1 << 16) - (45 << 16)) <= 16);
    // FreeType takes -180 degrees for the negative x axis
    try std.testing.expect(@abs(atan2(-1 << 16, 0) + pi) <= 16);
    try std.testing.expect(@abs(atan2(0, -1 << 16) + pi2) <= 16);
    try std.testing.expectEqual(@as(i32, 0), atan2(0, 0));
}
