//! The gamma transfer curves and YCbCr constants of sRGB and BT.601, 709 and 2020.

const std = @import("std");
const clamp = std.math.clamp;
const pow = std.math.pow;

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
