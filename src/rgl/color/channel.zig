//! Byte and float channel conversions, and the packed colors' byte order.

const std = @import("std");
const sysinfo = @import("builtin");
const clamp = std.math.clamp;

/// Rounds an already byte-scaled value (e.g. `128 + cb * 255`) to the nearest `u8`,
/// clamping first so an out-of-range result can't overflow the cast.
pub inline fn roundToByte(value: f32) u8 {
    return @trunc(clamp(@round(value), 0, 255));
}

/// Converts a `u8` channel to its normalized `[0, 1]` float value.
pub inline fn byteToFloat(value: u8) f32 {
    return @as(f32, @floatFromInt(value)) / 255.0;
}

/// Converts a normalized (and clamped) `[0, 1]` float channel to a `u8`.
pub inline fn floatToByte(value: f32) u8 {
    return roundToByte(clamp(value, 0, 1) * 255.0);
}

/// Converts a `u8` channel encoded as `offset + value * scale` back to `value`
/// (e.g. luma or chroma, whatever domain `value` is naturally in; there's no
/// clamp here since the arithmetic is just the inverse affine map).
pub inline fn byteToScaledOffset(value: u8, scale: f32, offset: f32) f32 {
    return (@as(f32, @floatFromInt(value)) - offset) / scale;
}

/// Encodes `value` as `offset + value * scale`, rounded to the nearest `u8`. There's
/// no domain-specific clamp on `value` itself (luma's is `[0, 1]`, chroma's is
/// `[-0.5, 0.5]`); `roundToByte`'s clamp to `[0, 255]` already saturates correctly
/// in the right direction for any of them, since the map is a monotonic affine one.
pub inline fn scaledOffsetToByte(value: f32, scale: f32, offset: f32) u8 {
    return roundToByte(offset + value * scale);
}

/// A big-endian target stores a packed color's fields last to first, so its
/// bytes in memory are not in field order; `toNative` fixes that.
/// `u32` values keep one meaning everywhere: `r` is `RGBA`'s low byte.
pub const swapped = sysinfo.target.cpu.arch.endian() == .big;
