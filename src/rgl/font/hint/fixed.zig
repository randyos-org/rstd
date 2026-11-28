//! FreeType's fixed-point arithmetic for the hinters: 26.6 pixels, 16.16
//! scales and 2.14 unit vectors, all in `i32` with FreeType's wrapping and
//! rounding, so every target gets the same bits (issue #373).
//!
//! Ported from fontations `font-types/src/fixed.rs` (`Fixed`'s `*`, `/` and
//! `mul_div`) and `skrifa/src/outline/glyf/hint/math.rs`, commit
//! 4f115482ebf24d1989e2c4446bf877d4812a268b (MIT): FreeType's `FT_MulFix`,
//! `FT_DivFix`, `FT_MulDiv`, `FT_MulDiv_No_Round` and `FT_Vector_NormLen`.

/// 26.6 floor.
pub fn floor(x: i32) i32 {
    return x & ~@as(i32, 63);
}

/// 26.6 round half up.
pub fn round(x: i32) i32 {
    return floor(x +% 32);
}

/// 26.6 ceiling.
pub fn ceil(x: i32) i32 {
    return floor(x +% 63);
}

/// `x` rounded to a multiple of `n`, a power of two.
pub fn roundPad(x: i32, n: i32) i32 {
    return (x +% @divTrunc(n, 2)) & ~(n - 1);
}

/// `a * b` in 16.16 (`FT_MulFix`).
pub fn mul(a: i32, b: i32) i32 {
    const ab = @as(i64, a) * b;
    return @truncate((ab + 0x8000 - @intFromBool(ab < 0)) >> 16);
}

/// `a / b` in 16.16 (`FT_DivFix`); 0x7FFFFFFF, signed, for `b == 0`.
pub fn div(a: i32, b: i32) i32 {
    const au: u64 = @abs(a);
    const bu: u64 = @abs(b);
    const q: u32 = if (bu == 0) 0x7FFFFFFF else @truncate(((au << 16) + (bu >> 1)) / bu);
    return signed(q, (a < 0) != (b < 0));
}

/// `a * b / c`, rounded (`FT_MulDiv`); 0x7FFFFFFF, signed, for `c == 0`.
pub fn mulDiv(a: i32, b: i32, c: i32) i32 {
    const p = unsignedParts(a, b, c);
    const r: u64 = if (p.c > 0) (p.a *% p.b +% (p.c >> 1)) / p.c else 0x7FFFFFFF;
    return signed(@truncate(r), p.neg);
}

/// `a * b / c`, truncated (`FT_MulDiv_No_Round`); 0x7FFFFFFF for `c == 0`.
pub fn mulDivNoRound(a: i32, b: i32, c: i32) i32 {
    const p = unsignedParts(a, b, c);
    const r: u64 = if (p.c > 0) (p.a *% p.b) / p.c else 0x7FFFFFFF;
    return signed(@truncate(r), p.neg);
}

/// `a * b` in 2.14, rounded half away from zero.
pub fn mul14(a: i32, b: i32) i32 {
    var v = @as(i64, a) * b;
    v += 0x2000 + (v >> 63);
    return @truncate(v >> 14);
}

/// The dot product of `(ax, ay)` and the 2.14 vector `(bx, by)`.
pub fn dot14(ax: i32, ay: i32, bx: i32, by: i32) i32 {
    var v = @as(i64, ax) * bx +% @as(i64, ay) * by;
    v +%= 0x2000 + (v >> 63);
    return @truncate(v >> 14);
}

/// `(x, y)` scaled to unit length in 2.14 (`FT_Vector_NormLen`, then `/ 4`).
pub fn normalize14(x: i32, y: i32) [2]i32 {
    var ux: u32 = @bitCast(x);
    var uy: u32 = @bitCast(y);
    var sx: i32 = 1;
    var sy: i32 = 1;
    if (x < 0) {
        ux = 0 -% ux;
        sx = -1;
    }
    if (y < 0) {
        uy = 0 -% uy;
        sy = -1;
    }
    if (ux == 0) return .{ @divTrunc(x, 4), if (uy > 0) @divTrunc(sy *% 0x10000, 4) else 0 };
    if (uy == 0) return .{ if (ux > 0) @divTrunc(sx *% 0x10000, 4) else 0, @divTrunc(y, 4) };
    var len = if (ux > uy) ux +% (uy >> 1) else uy +% (ux >> 1);
    var shift: i32 = @clz(len);
    shift -= 15 + @as(i32, @intFromBool(len >= (@as(u32, 0xAAAAAAAA) >> @intCast(shift))));
    if (shift > 0) {
        ux <<= @intCast(shift);
        uy <<= @intCast(shift);
        len = if (ux > uy) ux +% (uy >> 1) else uy +% (ux >> 1);
    } else {
        const s: u5 = @intCast(-shift);
        ux >>= s;
        uy >>= s;
        len >>= s;
    }
    var b: i32 = 0x10000 -% @as(i32, @bitCast(len));
    const xs: i32 = @bitCast(ux);
    const ys: i32 = @bitCast(uy);
    var u: u32 = undefined;
    var v: u32 = undefined;
    while (true) {
        u = @bitCast(xs +% ((xs *% b) >> 16));
        v = @bitCast(ys +% ((ys *% b) >> 16));
        var z: i32 = @divTrunc(-%@as(i32, @bitCast(u *% u +% v *% v)), 0x200);
        z = @divTrunc(z *% ((0x10000 +% b) >> 8), 0x10000);
        b +%= z;
        if (z <= 0) break;
    }
    return .{ @divTrunc(@as(i32, @bitCast(u)) *% sx, 4), @divTrunc(@as(i32, @bitCast(v)) *% sy, 4) };
}

const Parts = struct { a: u64, b: u64, c: u64, neg: bool };

/// The magnitudes of `a`, `b` and `c` as FreeType takes them (`0 - x`,
/// wrapping), and whether the result is negative.
fn unsignedParts(a: i32, b: i32, c: i32) Parts {
    return .{
        .a = magnitude(a),
        .b = magnitude(b),
        .c = magnitude(c),
        .neg = ((a < 0) != (b < 0)) != (c < 0),
    };
}

fn magnitude(x: i32) u64 {
    const w: u64 = @bitCast(@as(i64, x));
    return if (x < 0) 0 -% w else w;
}

fn signed(q: u32, neg: bool) i32 {
    const s: i32 = @bitCast(q);
    return if (neg) -%s else s;
}
