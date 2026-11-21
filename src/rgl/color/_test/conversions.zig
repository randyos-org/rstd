//! The color conversions against reference tables, and their round trips.

const std = @import("std");
const zigstd = @import("../../../zigstd.zig");
const rgl_color = @import("../__root__.zig");
const GammaLawConsts = rgl_color.GammaLawConsts;
const YCbCrConsts = rgl_color.YCbCrConsts;
const YCbCrRange = rgl_color.YCbCrRange;
const gammaEncodeChannel = rgl_color.gammaEncodeChannel;
const gammaDecodeChannel = rgl_color.gammaDecodeChannel;
const SrgbStd = rgl_color.SrgbStd;
const Rec601Std = rgl_color.Rec601Std;
const Rec709Std = rgl_color.Rec709Std;
const Rec2020Std = rgl_color.Rec2020Std;
const RGBA = rgl_color.RGBA;
const RGBAf32 = rgl_color.RGBAf32;
const HSLA = rgl_color.HSLA;
const HSLAf32 = rgl_color.HSLAf32;
const YCbCrA = rgl_color.YCbCrA;
const YCbCrAf32 = rgl_color.YCbCrAf32;

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
