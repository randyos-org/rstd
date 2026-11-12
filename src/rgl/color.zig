const std = @import("std");
const math = std.math;
const clamp = math.clamp;
const pow = math.pow;
const zigstd = @import("../zigstd.zig");
const sysinfo = @import("builtin");

/// The linear-segment plus power-law transfer curve every gamma standard uses.
/// Encode: `linear_slope * L` for `L <= threshold`, else
/// `power_scale * L^gamma - (power_scale - 1)`; that bias makes `V(1) = 1`.
/// Decode is the inverse, with threshold `threshold * linear_slope`.
pub const GammaLawConsts = struct {
    threshold: f32,
    linear_slope: f32,
    power_scale: f32,
    encode_gamma: f32,
};

pub const YCbCrConsts = struct {
    kr: f32,
    kb: f32,
};

/// Byte-quantization convention for Y/Cb/Cr. Full range spans the entire 0-255 byte
/// (the JPEG/JFIF convention); safe range reserves footroom/headroom the way
/// broadcast video does (luma 16-235, chroma 16-240).
/// Source: ITU-R BT.601-7 (2011) §2.5.3.
pub const YCbCrRange = enum {
    full,
    safe,
};

/// The `offset`/`scale` pair for luma and the shared chroma `scale` (chroma is
/// always centered on byte 128) that a `YCbCrRange` maps onto.
pub const YCbCrRangeConsts = struct {
    luma_scale: f32,
    luma_offset: f32,
    chroma_scale: f32,
};

pub fn rangeConsts(range: YCbCrRange) YCbCrRangeConsts {
    return switch (range) {
        .full => .{ .luma_scale = 255.0, .luma_offset = 0.0, .chroma_scale = 255.0 },
        .safe => .{ .luma_scale = 219.0, .luma_offset = 16.0, .chroma_scale = 224.0 },
    };
}

/// BT.601 SD video (1982).
/// Source: ITU-R BT.601-7 (2011) §2.5.
pub const Rec601Std = struct {
    pub const gamma_law: GammaLawConsts = .{
        .threshold = 0.018,
        .linear_slope = 4.5,
        .power_scale = 1.099,
        .encode_gamma = 0.45,
    };
    pub const ycbcr: YCbCrConsts = .{ .kr = 0.299, .kb = 0.114 };
};

/// BT.709 HD video (1990).
/// Only the YCbCr matrix actually differs, reflecting BT.709's different color primaries.
/// Source: ITU-R BT.709-6 (2015) §3 (kr/kb).
pub const Rec709Std = struct {
    pub const gamma_law: GammaLawConsts = Rec601Std.gamma_law;
    pub const ycbcr: YCbCrConsts = .{ .kr = 0.2126, .kb = 0.0722 };
};

/// BT.2020 UHD video. `.gamma_law` uses the recommendation's full, unrounded
/// precision (rounded 10-/12-bit variants are also permitted, but not used here).
/// Source: ITU-R BT.2020-2 (2015) table 4 (gamma); §4 (kr/kb).
pub const Rec2020Std = struct {
    pub const gamma_law: GammaLawConsts = .{
        .threshold = 0.018053968510807,
        .linear_slope = 4.5,
        .power_scale = 1.09929682680944,
        .encode_gamma = 0.45,
    };
    pub const ycbcr: YCbCrConsts = .{ .kr = 0.2627, .kb = 0.0593 };
};

/// sRGB display transfer function. sRGB itself defines no YCbCr matrix; `.ycbcr`
/// here just borrows `Rec601Std`'s, per the JPEG/JFIF (ITU-T T.871) convention of
/// pairing sRGB-gamma stills with the BT.601 matrix. That's a convention, not an
/// inherent property of sRGB.
/// Source: IEC 61966-2-1:1999.
pub const SrgbStd = struct {
    pub const gamma_law: GammaLawConsts = .{
        .threshold = 0.0031308,
        .linear_slope = 12.92,
        .power_scale = 1.055,
        // The commonly-cited "sRGB gamma" of 2.4 is the *decode* exponent
        .encode_gamma = 1.0 / 2.4,
    };
    pub const ycbcr: YCbCrConsts = Rec601Std.ycbcr;
};

pub fn gammaEncodeChannel(law: GammaLawConsts, linear: f32) f32 {
    const value = clamp(linear, 0, 1);
    const power_bias = law.power_scale - 1;
    return if (value <= law.threshold)
        value * law.linear_slope
    else
        law.power_scale * pow(f32, value, law.encode_gamma) - power_bias;
}

pub fn gammaDecodeChannel(law: GammaLawConsts, encoded: f32) f32 {
    const value = clamp(encoded, 0, 1);
    const power_bias = law.power_scale - 1;
    const encoded_threshold = law.threshold * law.linear_slope;
    return if (value <= encoded_threshold)
        value / law.linear_slope
    else
        pow(f32, (value + power_bias) / law.power_scale, 1.0 / law.encode_gamma);
}

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

pub const RGBA = packed struct {
    const Self = @This();
    r: u8,
    g: u8,
    b: u8,
    a: u8,

    pub const blank: Self = .{ .r = 0, .g = 0, .b = 0, .a = 0 };

    pub inline fn equals(self: Self, other: Self) bool {
        return self.toU32() == other.toU32();
    }

    pub fn print(self: Self) void {
        std.debug.print("<RGBA {d} {d} {d} {d}>", .{ self.r, self.g, self.b, self.a });
    }

    pub fn printHex(self: Self) void {
        std.debug.print("#{X}{X}{X}{X}", .{ self.r, self.g, self.b, self.a });
    }

    pub fn printRgbHex(self: Self) void {
        std.debug.print("#{X}{X}{X}", .{ self.r, self.g, self.b });
    }

    pub inline fn toU32(self: Self) u32 {
        return @bitCast(self);
    }

    /// `self` with its bytes in field order in memory, for files, codecs and
    /// GPU uploads: byte-swapped where `swapped`, its fields then reading
    /// wrong. Its own inverse, so it also reads such bytes back.
    pub inline fn toNative(self: Self) Self {
        return if (swapped) @bitCast(@byteSwap(self.toU32())) else self;
    }

    pub inline fn toArray(self: Self) [4]u8 {
        return .{ self.r, self.g, self.b, self.a };
    }

    pub inline fn toRGBAf32(self: Self) RGBAf32 {
        return .{
            .r = byteToFloat(self.r),
            .g = byteToFloat(self.g),
            .b = byteToFloat(self.b),
            .a = byteToFloat(self.a),
        };
    }

    pub inline fn toBGRA(self: Self) BGRA {
        return .{ .b = self.b, .g = self.g, .r = self.r, .a = self.a };
    }
};

/// used by GLSL
pub const RGBAf32 = packed struct {
    const Self = @This();
    r: f32,
    g: f32,
    b: f32,
    a: f32,

    pub const blank: Self = .{ .r = 0, .g = 0, .b = 0, .a = 0 };

    pub inline fn equals(self: Self, other: Self) bool {
        return self.r == other.r and self.g == other.g and self.b == other.b and self.a == other.a;
    }

    pub fn print(self: Self) void {
        std.debug.print("<RGBAf32 {f} {f} {f} {f}>", .{ self.r, self.g, self.b, self.a });
    }

    pub inline fn toArray(self: Self) [4]f32 {
        return .{ self.r, self.g, self.b, self.a };
    }

    pub inline fn toRGBA(self: Self) RGBA {
        return .{
            .r = floatToByte(self.r),
            .g = floatToByte(self.g),
            .b = floatToByte(self.b),
            .a = floatToByte(self.a),
        };
    }

    pub fn toHSLAf32(self: Self) HSLAf32 {
        const red = clamp(self.r, 0, 1);
        const green = clamp(self.g, 0, 1);
        const blue = clamp(self.b, 0, 1);
        const max_channel = @max(red, green, blue);
        const min_channel = @min(red, green, blue);
        const chroma = max_channel - min_channel;
        var hue_prime: f32 = 0;
        if (chroma != 0) {
            if (max_channel == red) {
                hue_prime = @mod((green - blue) / chroma, 6.0);
            } else if (max_channel == green) {
                hue_prime = (blue - red) / chroma + 2;
            } else if (max_channel == blue) {
                hue_prime = (red - green) / chroma + 4;
            }
        }
        const lightness = (min_channel + max_channel) / 2;
        const saturation: f32 = if (chroma == 0) 0 else chroma / (1 - @abs(2 * lightness - 1));
        return .{ .h = hue_prime / 6, .s = saturation, .l = lightness, .a = clamp(self.a, 0, 1) };
    }

    /// Converts gamma-encoded (display-ready) R'G'B' to linear-light RGB. Alpha is
    /// always linear and passes through unchanged.
    pub fn toLinear(self: Self, law: GammaLawConsts) Self {
        return .{
            .r = gammaDecodeChannel(law, self.r),
            .g = gammaDecodeChannel(law, self.g),
            .b = gammaDecodeChannel(law, self.b),
            .a = self.a,
        };
    }

    /// Converts linear-light RGB to a gamma-encoded (display-ready) R'G'B' signal.
    /// Alpha is always linear and passes through unchanged.
    pub fn toGammaEncoded(self: Self, law: GammaLawConsts) Self {
        return .{
            .r = gammaEncodeChannel(law, self.r),
            .g = gammaEncodeChannel(law, self.g),
            .b = gammaEncodeChannel(law, self.b),
            .a = self.a,
        };
    }

    /// Derives luma/chroma from gamma-encoded R'G'B' per `consts`' matrix.
    pub fn toYCbCrAf32(self: Self, consts: YCbCrConsts) YCbCrAf32 {
        const kr = consts.kr;
        const kb = consts.kb;
        const kg = 1 - kr - kb;
        const r = clamp(self.r, 0, 1);
        const g = clamp(self.g, 0, 1);
        const b = clamp(self.b, 0, 1);
        const y = kr * r + kg * g + kb * b;
        const cb = (b - y) / (2 * (1 - kb));
        const cr = (r - y) / (2 * (1 - kr));
        return .{ .y = y, .cb = cb, .cr = cr, .a = clamp(self.a, 0, 1) };
    }
};

pub const BGRA = packed struct {
    const Self = @This();
    b: u8,
    g: u8,
    r: u8,
    a: u8,

    pub inline fn equals(self: Self, other: Self) bool {
        return self.toU32() == other.toU32();
    }

    pub fn print(self: Self) void {
        std.debug.print("<BGRA {d} {d} {d} {d}>", .{ self.b, self.g, self.r, self.a });
    }

    pub inline fn toU32(self: Self) u32 {
        return @bitCast(self);
    }

    /// `self` with its bytes in field order in memory, for files, codecs and
    /// GPU uploads: byte-swapped where `swapped`, its fields then reading
    /// wrong. Its own inverse, so it also reads such bytes back.
    pub inline fn toNative(self: Self) Self {
        return if (swapped) @bitCast(@byteSwap(self.toU32())) else self;
    }

    pub inline fn toRGBA(self: Self) RGBA {
        return .{ .r = self.r, .g = self.g, .b = self.b, .a = self.a };
    }
};

pub const HSLAf32 = packed struct {
    const Self = @This();
    h: f32,
    s: f32,
    l: f32,
    a: f32,

    pub inline fn equals(self: Self, other: Self) bool {
        return self.h == other.h and self.s == other.s and self.l == other.l and self.a == other.a;
    }

    pub fn print(self: Self) void {
        std.debug.print("<HSLAf32 {f} {f} {f} {f}>", .{ self.h, self.s, self.l, self.a });
    }

    pub inline fn toArray(self: Self) [4]f32 {
        return .{ self.h, self.s, self.l, self.a };
    }

    pub fn toRGBAf32(self: Self) RGBAf32 {
        const s = clamp(self.s, 0, 1);
        const l = clamp(self.l, 0, 1);
        const chroma = (1 - @abs(2 * l - 1)) * s;
        const hprime = @mod(self.h, 1.0) * 6.0;
        const second_largest = chroma * (1 - @abs(@mod(hprime, 2.0) - 1));
        const rgbprime: [3]f32 = if (math.isNan(self.h))
            .{ 0, 0, 0 }
        else if (hprime <= 1)
            .{ chroma, second_largest, 0 }
        else if (hprime <= 2)
            .{ second_largest, chroma, 0 }
        else if (hprime <= 3)
            .{ 0, chroma, second_largest }
        else if (hprime <= 4)
            .{ 0, second_largest, chroma }
        else if (hprime <= 5)
            .{ second_largest, 0, chroma }
        else
            .{ chroma, 0, second_largest };
        const lightness_match = l - chroma * 0.5;
        return .{
            .r = clamp(rgbprime[0] + lightness_match, 0, 1),
            .g = clamp(rgbprime[1] + lightness_match, 0, 1),
            .b = clamp(rgbprime[2] + lightness_match, 0, 1),
            .a = clamp(self.a, 0, 1),
        };
    }

    pub fn toRGBA(self: Self) RGBA {
        return self.toRGBAf32().toRGBA();
    }

    pub fn toHSLA(self: Self) HSLA {
        return .{
            .h = floatToByte(self.h),
            .s = floatToByte(self.s),
            .l = floatToByte(self.l),
            .a = floatToByte(self.a),
        };
    }
};

pub const HSLA = packed struct {
    const Self = @This();
    h: u8,
    s: u8,
    l: u8,
    a: u8,

    pub inline fn equals(self: Self, other: Self) bool {
        return self.toU32() == other.toU32();
    }

    pub fn print(self: Self) void {
        std.debug.print("<HSLA {d} {d} {d} {d}>", .{ self.h, self.s, self.l, self.a });
    }

    pub inline fn toU32(self: Self) u32 {
        return @bitCast(self);
    }

    /// `self` with its bytes in field order in memory, for files, codecs and
    /// GPU uploads: byte-swapped where `swapped`, its fields then reading
    /// wrong. Its own inverse, so it also reads such bytes back.
    pub inline fn toNative(self: Self) Self {
        return if (swapped) @bitCast(@byteSwap(self.toU32())) else self;
    }

    pub inline fn toHSLAf32(self: Self) HSLAf32 {
        return .{
            .h = byteToFloat(self.h),
            .s = byteToFloat(self.s),
            .l = byteToFloat(self.l),
            .a = byteToFloat(self.a),
        };
    }

    pub inline fn toRGBAf32(self: Self) RGBAf32 {
        return self.toHSLAf32().toRGBAf32();
    }

    pub inline fn toRGBA(self: Self) RGBA {
        return self.toRGBAf32().toRGBA();
    }
};

/// Luma (`y`, in `[0, 1]`) and chroma differences (`cb`/`cr`, in `[-0.5, 0.5]`)
/// derived from gamma-encoded R'G'B'.
pub const YCbCrAf32 = packed struct {
    const Self = @This();
    y: f32,
    cb: f32,
    cr: f32,
    a: f32,

    pub inline fn equals(self: Self, other: Self) bool {
        return self.y == other.y and self.cb == other.cb and self.cr == other.cr and self.a == other.a;
    }

    pub fn print(self: Self) void {
        std.debug.print("<YCbCrAf32 {f} {f} {f} {f}>", .{ self.y, self.cb, self.cr, self.a });
    }

    pub inline fn toArray(self: Self) [4]f32 {
        return .{ self.y, self.cb, self.cr, self.a };
    }

    pub fn toRGBAf32(self: Self, consts: YCbCrConsts) RGBAf32 {
        const kr = consts.kr;
        const kb = consts.kb;
        const kg = 1 - kr - kb;
        const y = self.y;
        const r = y + 2 * (1 - kr) * self.cr;
        const b = y + 2 * (1 - kb) * self.cb;
        const g = (y - kr * r - kb * b) / kg;
        return .{
            .r = clamp(r, 0, 1),
            .g = clamp(g, 0, 1),
            .b = clamp(b, 0, 1),
            .a = clamp(self.a, 0, 1),
        };
    }

    pub fn toYCbCrA(self: Self, range: YCbCrRange) YCbCrA {
        const r = rangeConsts(range);
        return .{
            .y = scaledOffsetToByte(self.y, r.luma_scale, r.luma_offset),
            .cb = scaledOffsetToByte(self.cb, r.chroma_scale, 128.0),
            .cr = scaledOffsetToByte(self.cr, r.chroma_scale, 128.0),
            .a = floatToByte(self.a),
        };
    }

    pub fn toRGBA(self: Self, consts: YCbCrConsts) RGBA {
        return self.toRGBAf32(consts).toRGBA();
    }
};

/// Byte-quantized Y'CbCr, per `YCbCrRange` (full or safe). Which `YCbCrConsts`
/// matrix it was derived from is not stored; pass the same one back in when
/// converting back to RGB.
pub const YCbCrA = packed struct {
    const Self = @This();
    y: u8,
    cb: u8,
    cr: u8,
    a: u8,

    pub inline fn equals(self: Self, other: Self) bool {
        return self.toU32() == other.toU32();
    }

    pub fn print(self: Self) void {
        std.debug.print("<YCbCrA {d} {d} {d} {d}>", .{ self.y, self.cb, self.cr, self.a });
    }

    pub inline fn toU32(self: Self) u32 {
        return @bitCast(self);
    }

    /// `self` with its bytes in field order in memory, for files, codecs and
    /// GPU uploads: byte-swapped where `swapped`, its fields then reading
    /// wrong. Its own inverse, so it also reads such bytes back.
    pub inline fn toNative(self: Self) Self {
        return if (swapped) @bitCast(@byteSwap(self.toU32())) else self;
    }

    pub inline fn toArray(self: Self) [4]u8 {
        return .{ self.y, self.cb, self.cr, self.a };
    }

    pub fn toYCbCrAf32(self: Self, range: YCbCrRange) YCbCrAf32 {
        const r = rangeConsts(range);
        return .{
            .y = byteToScaledOffset(self.y, r.luma_scale, r.luma_offset),
            .cb = byteToScaledOffset(self.cb, r.chroma_scale, 128.0),
            .cr = byteToScaledOffset(self.cr, r.chroma_scale, 128.0),
            .a = byteToFloat(self.a),
        };
    }

    pub fn toRGBAf32(self: Self, consts: YCbCrConsts, range: YCbCrRange) RGBAf32 {
        return self.toYCbCrAf32(range).toRGBAf32(consts);
    }

    pub fn toRGBA(self: Self, consts: YCbCrConsts, range: YCbCrRange) RGBA {
        return self.toRGBAf32(consts, range).toRGBA();
    }
};

// Reference test vectors for sRGB <-> HSL conversion.
//
// Source: Wikipedia, "HSL and HSV" article, Examples section table (the full 19-row
// table, R, G, B, H, S_HSL, L columns): https://en.wikipedia.org/wiki/HSL_and_HSV#Examples
// (retrieved 2026-08-23). `h_deg` is `0` for the achromatic rows (white/gray/black),
// where the table lists hue as undefined ("N/A").
const _HslRgbSample = struct {
    r: f32,
    g: f32,
    b: f32,
    h_deg: f32,
    s: f32,
    l: f32,
};

const _hsl_rgb_samples = [_]_HslRgbSample{
    .{ .r = 1.000, .g = 1.000, .b = 1.000, .h_deg = 0.0, .s = 0.000, .l = 1.000 }, // #FFFFFF
    .{ .r = 0.500, .g = 0.500, .b = 0.500, .h_deg = 0.0, .s = 0.000, .l = 0.500 }, // #808080
    .{ .r = 0.000, .g = 0.000, .b = 0.000, .h_deg = 0.0, .s = 0.000, .l = 0.000 }, // #000000
    .{ .r = 1.000, .g = 0.000, .b = 0.000, .h_deg = 0.0, .s = 1.000, .l = 0.500 }, // #FF0000
    .{ .r = 0.750, .g = 0.750, .b = 0.000, .h_deg = 60.0, .s = 1.000, .l = 0.375 }, // #BFBF00
    .{ .r = 0.000, .g = 0.500, .b = 0.000, .h_deg = 120.0, .s = 1.000, .l = 0.250 }, // #008000
    .{ .r = 0.500, .g = 1.000, .b = 1.000, .h_deg = 180.0, .s = 1.000, .l = 0.750 }, // #80FFFF
    .{ .r = 0.500, .g = 0.500, .b = 1.000, .h_deg = 240.0, .s = 1.000, .l = 0.750 }, // #8080FF
    .{ .r = 0.750, .g = 0.250, .b = 0.750, .h_deg = 300.0, .s = 0.500, .l = 0.500 }, // #BF40BF
    .{ .r = 0.628, .g = 0.643, .b = 0.142, .h_deg = 61.8, .s = 0.638, .l = 0.393 }, // #A0A424
    .{ .r = 0.255, .g = 0.104, .b = 0.918, .h_deg = 251.1, .s = 0.832, .l = 0.511 }, // #411BEA
    .{ .r = 0.931, .g = 0.463, .b = 0.316, .h_deg = 14.3, .s = 0.817, .l = 0.624 }, // #ED7651
    .{ .r = 0.116, .g = 0.675, .b = 0.255, .h_deg = 134.9, .s = 0.707, .l = 0.396 }, // #1EAC41
    .{ .r = 0.941, .g = 0.785, .b = 0.053, .h_deg = 49.5, .s = 0.893, .l = 0.497 }, // #F0C80E
    .{ .r = 0.704, .g = 0.187, .b = 0.897, .h_deg = 283.7, .s = 0.775, .l = 0.542 }, // #B430E5
    .{ .r = 0.998, .g = 0.974, .b = 0.532, .h_deg = 56.9, .s = 0.991, .l = 0.765 }, // #FEF888
    .{ .r = 0.099, .g = 0.795, .b = 0.591, .h_deg = 162.4, .s = 0.779, .l = 0.447 }, // #19CB97
    .{ .r = 0.211, .g = 0.149, .b = 0.597, .h_deg = 248.3, .s = 0.601, .l = 0.373 }, // #362698
    .{ .r = 0.495, .g = 0.493, .b = 0.721, .h_deg = 240.5, .s = 0.290, .l = 0.607 }, // #7E7EB8
};

test "RGBAf32.toHSLAf32 matches the Wikipedia HSL reference table" {
    for (_hsl_rgb_samples) |sample| {
        const rgba = RGBAf32{ .r = sample.r, .g = sample.g, .b = sample.b, .a = 1 };
        const hsl = rgba.toHSLAf32();
        if (sample.s != 0) {
            // hue is meaningless (and left at 0) for achromatic colors
            try zigstd.testExpectApproxEqAbs(sample.h_deg / 360.0, hsl.h, 0.003);
        }
        try zigstd.testExpectApproxEqAbs(sample.s, hsl.s, 0.003);
        try zigstd.testExpectApproxEqAbs(sample.l, hsl.l, 0.003);
    }
}

test "HSLAf32.toRGBAf32 matches the Wikipedia HSL reference table" {
    for (_hsl_rgb_samples) |sample| {
        const hsl = HSLAf32{ .h = sample.h_deg / 360.0, .s = sample.s, .l = sample.l, .a = 1 };
        const rgba = hsl.toRGBAf32();
        try zigstd.testExpectApproxEqAbs(sample.r, rgba.r, 0.003);
        try zigstd.testExpectApproxEqAbs(sample.g, rgba.g, 0.003);
        try zigstd.testExpectApproxEqAbs(sample.b, rgba.b, 0.003);
    }
}

test "a u32 holds r low on every target, and toNative puts the bytes in field order" {
    const c: RGBA = .{ .r = 1, .g = 2, .b = 3, .a = 4 };
    try zigstd.testExpectEqual(@as(u32, 0x04030201), c.toU32());
    const from_u32: RGBA = @bitCast(@as(u32, 0x04030201));
    try zigstd.testExpectTrue(from_u32.equals(c));
    try zigstd.testExpectEqualSlices(u8, &.{ 1, 2, 3, 4 }, std.mem.asBytes(&c.toNative()));
    try zigstd.testExpectTrue(c.toNative().toNative().equals(c));
}

test "RGBA <-> RGBAf32 round-trips every byte value exactly" {
    var channel: u16 = 0;
    while (channel <= 255) : (channel += 1) {
        const byte: u8 = @intCast(channel);
        const original = RGBA{ .r = byte, .g = byte, .b = byte, .a = byte };
        try zigstd.testExpectEqual(original, original.toRGBAf32().toRGBA());
    }
}

test "HSLA <-> HSLAf32 round-trips every byte value exactly" {
    var channel: u16 = 0;
    while (channel <= 255) : (channel += 1) {
        const byte: u8 = @intCast(channel);
        const original = HSLA{ .h = byte, .s = byte, .l = byte, .a = byte };
        try zigstd.testExpectEqual(original, original.toHSLAf32().toHSLA());
    }
}

// Gamma (transfer function) reference points and round-trip checks.
//
// The exact 0/1 endpoints and the sRGB midpoint are well-known fixed points of these
// curves (e.g. sRGB-encoded 0.5 corresponds to ~21.4% linear light; see e.g.
// https://en.wikipedia.org/wiki/SRGB#Transfer_function_(%22gamma%22)); the rest is
// a round-trip property check against the formulas cited above, not an external table.
const _gamma_laws = [_]GammaLawConsts{ SrgbStd.gamma_law, Rec601Std.gamma_law, Rec709Std.gamma_law, Rec2020Std.gamma_law };
const _gamma_round_trip_samples = [_]f32{ 0.0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0 };

test "gamma encode/decode fixed points (0, 1, and sRGB midpoint)" {
    for (_gamma_laws) |law| {
        try zigstd.testExpectApproxEqAbs(@as(f32, 0.0), gammaEncodeChannel(law, 0.0), 1e-6);
        try zigstd.testExpectApproxEqAbs(@as(f32, 1.0), gammaEncodeChannel(law, 1.0), 1e-6);
        try zigstd.testExpectApproxEqAbs(@as(f32, 0.0), gammaDecodeChannel(law, 0.0), 1e-6);
        try zigstd.testExpectApproxEqAbs(@as(f32, 1.0), gammaDecodeChannel(law, 1.0), 1e-6);
    }
    try zigstd.testExpectApproxEqAbs(@as(f32, 0.214041), gammaDecodeChannel(SrgbStd.gamma_law, 0.5), 0.0005);
}

test "RGBAf32.toLinear/toGammaEncoded round-trip for every standard" {
    for (_gamma_laws) |law| {
        for (_gamma_round_trip_samples) |value| {
            const rgba = RGBAf32{ .r = value, .g = value, .b = value, .a = 1 };
            const round_tripped = rgba.toGammaEncoded(law).toLinear(law);
            try zigstd.testExpectApproxEqAbs(value, round_tripped.r, 1e-4);
            try zigstd.testExpectApproxEqAbs(value, round_tripped.g, 1e-4);
            try zigstd.testExpectApproxEqAbs(value, round_tripped.b, 1e-4);
        }
    }
}

test "rec601 and rec709 share the same gamma curve" {
    for (_gamma_round_trip_samples) |value| {
        try zigstd.testExpectEqual(gammaEncodeChannel(Rec601Std.gamma_law, value), gammaEncodeChannel(Rec709Std.gamma_law, value));
        try zigstd.testExpectEqual(gammaDecodeChannel(Rec601Std.gamma_law, value), gammaDecodeChannel(Rec709Std.gamma_law, value));
    }
}

// YCbCr reference vectors: BT.601 full range from Wikipedia, "YCbCr", "JPEG
// conversion" (https://en.wikipedia.org/wiki/YCbCr, retrieved 2026-08-23),
// applied here to the primaries (matching the usual JPEG values). `.srgb`
// shares the BT.601 matrix (`SrgbStd`), so they apply to it too.
const _YCbCr601FullRangeSample = struct {
    r: f32,
    g: f32,
    b: f32,
    y: u8,
    cb: u8,
    cr: u8,
};

const _ycbcr_601_full_range_samples = [_]_YCbCr601FullRangeSample{
    .{ .r = 1.000, .g = 1.000, .b = 1.000, .y = 255, .cb = 128, .cr = 128 }, // white
    .{ .r = 0.000, .g = 0.000, .b = 0.000, .y = 0, .cb = 128, .cr = 128 }, // black
    .{ .r = 1.000, .g = 0.000, .b = 0.000, .y = 76, .cb = 85, .cr = 255 }, // red
    .{ .r = 0.000, .g = 1.000, .b = 0.000, .y = 150, .cb = 44, .cr = 21 }, // green
    .{ .r = 0.000, .g = 0.000, .b = 1.000, .y = 29, .cb = 255, .cr = 107 }, // blue
};

test "RGBAf32.toYCbCrAf32 (rec601, full range) matches known JPEG conversion values" {
    for (_ycbcr_601_full_range_samples) |sample| {
        const rgba = RGBAf32{ .r = sample.r, .g = sample.g, .b = sample.b, .a = 1 };
        const ycbcr = rgba.toYCbCrAf32(Rec601Std.ycbcr).toYCbCrA(.full);
        try zigstd.testExpectEqual(sample.y, ycbcr.y);
        try zigstd.testExpectEqual(sample.cb, ycbcr.cb);
        try zigstd.testExpectEqual(sample.cr, ycbcr.cr);
    }
}

test "YCbCrA.toRGBA (rec601, full range) matches known JPEG conversion values" {
    for (_ycbcr_601_full_range_samples) |sample| {
        const ycbcr = YCbCrA{ .y = sample.y, .cb = sample.cb, .cr = sample.cr, .a = 255 };
        const rgba = ycbcr.toRGBA(Rec601Std.ycbcr, .full);
        try zigstd.testExpectApproxEqAbs(sample.r, rgba.toRGBAf32().r, 0.004);
        try zigstd.testExpectApproxEqAbs(sample.g, rgba.toRGBAf32().g, 0.004);
        try zigstd.testExpectApproxEqAbs(sample.b, rgba.toRGBAf32().b, 0.004);
    }
}

test "srgb's YCbCr matrix matches rec601's" {
    for (_ycbcr_601_full_range_samples) |sample| {
        const rgba = RGBAf32{ .r = sample.r, .g = sample.g, .b = sample.b, .a = 1 };
        try zigstd.testExpectEqual(rgba.toYCbCrAf32(Rec601Std.ycbcr), rgba.toYCbCrAf32(SrgbStd.ycbcr));
    }
}

const ycbcr_consts = [_]YCbCrConsts{ SrgbStd.ycbcr, Rec601Std.ycbcr, Rec709Std.ycbcr, Rec2020Std.ycbcr };

test "safe-range Y hits the standard 16 (black) / 235 (white) video levels" {
    for (ycbcr_consts) |consts| {
        const black = (RGBAf32{ .r = 0, .g = 0, .b = 0, .a = 1 }).toYCbCrAf32(consts).toYCbCrA(.safe);
        const white = (RGBAf32{ .r = 1, .g = 1, .b = 1, .a = 1 }).toYCbCrAf32(consts).toYCbCrA(.safe);
        try zigstd.testExpectEqual(@as(u8, 16), black.y);
        try zigstd.testExpectEqual(@as(u8, 128), black.cb);
        try zigstd.testExpectEqual(@as(u8, 128), black.cr);
        try zigstd.testExpectEqual(@as(u8, 235), white.y);
        try zigstd.testExpectEqual(@as(u8, 128), white.cb);
        try zigstd.testExpectEqual(@as(u8, 128), white.cr);
    }
}

test "RGBAf32 <-> YCbCrAf32 round-trips for every standard and range" {
    const ranges = [_]YCbCrRange{ .full, .safe };
    const colors = [_][3]f32{
        .{ 0.9, 0.2, 0.1 },
        .{ 0.1, 0.8, 0.3 },
        .{ 0.2, 0.3, 0.9 },
        .{ 0.5, 0.5, 0.5 },
    };
    for (ycbcr_consts) |consts| {
        for (ranges) |range| {
            for (colors) |color| {
                const rgba = RGBAf32{ .r = color[0], .g = color[1], .b = color[2], .a = 1 };
                const round_tripped = rgba.toYCbCrAf32(consts).toYCbCrA(range).toYCbCrAf32(range).toRGBAf32(consts);
                // Byte quantization (safe range in particular, with its narrower
                // 219/224-step scale) introduces rounding noise that the matrix
                // inversion can amplify, so this needs a looser tolerance than the
                // exact-formula tests above.
                try zigstd.testExpectApproxEqAbs(color[0], round_tripped.r, 0.01);
                try zigstd.testExpectApproxEqAbs(color[1], round_tripped.g, 0.01);
                try zigstd.testExpectApproxEqAbs(color[2], round_tripped.b, 0.01);
            }
        }
    }
}

test "f32 colors' equals compares every component, alpha too" {
    const rgba: RGBAf32 = .{ .r = 0.25, .g = 0.5, .b = 0.75, .a = 1 };
    var rgba2 = rgba;
    try zigstd.testExpectTrue(rgba.equals(rgba2));
    rgba2.a = 0.5;
    try zigstd.testExpectTrue(!rgba.equals(rgba2));
    const hsla: HSLAf32 = .{ .h = 0.25, .s = 0.5, .l = 0.75, .a = 1 };
    try zigstd.testExpectTrue(hsla.equals(hsla));
    try zigstd.testExpectTrue(!hsla.equals(.{ .h = 0.25, .s = 0.5, .l = 0.75, .a = 0.5 }));
    const ycbcr: YCbCrAf32 = .{ .y = 0.25, .cb = 0.5, .cr = 0.75, .a = 1 };
    try zigstd.testExpectTrue(ycbcr.equals(ycbcr));
    try zigstd.testExpectTrue(!ycbcr.equals(.{ .y = 0.25, .cb = 0.5, .cr = 0.75, .a = 0.5 }));
}
