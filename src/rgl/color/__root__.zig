//! Colors: the standards' constants, channel conversions, and the RGB, HSL
//! and Y'CbCr types.

const standards = @import("standards.zig");
const channel = @import("channel.zig");
const rgb = @import("rgb.zig");
const hsl = @import("hsl.zig");
const ycbcr = @import("ycbcr.zig");

pub const GammaLawConsts = standards.GammaLawConsts;
pub const YCbCrConsts = standards.YCbCrConsts;
pub const YCbCrRange = standards.YCbCrRange;
pub const YCbCrRangeConsts = standards.YCbCrRangeConsts;
pub const rangeConsts = standards.rangeConsts;
pub const Rec601Std = standards.Rec601Std;
pub const Rec709Std = standards.Rec709Std;
pub const Rec2020Std = standards.Rec2020Std;
pub const SrgbStd = standards.SrgbStd;
pub const gammaEncodeChannel = standards.gammaEncodeChannel;
pub const gammaDecodeChannel = standards.gammaDecodeChannel;

pub const roundToByte = channel.roundToByte;
pub const byteToFloat = channel.byteToFloat;
pub const floatToByte = channel.floatToByte;
pub const byteToScaledOffset = channel.byteToScaledOffset;
pub const scaledOffsetToByte = channel.scaledOffsetToByte;
pub const swapped = channel.swapped;

pub const RGBA = rgb.RGBA;
pub const RGBAf32 = rgb.RGBAf32;
pub const BGRA = rgb.BGRA;
pub const HSLA = hsl.HSLA;
pub const HSLAf32 = hsl.HSLAf32;
pub const YCbCrA = ycbcr.YCbCrA;
pub const YCbCrAf32 = ycbcr.YCbCrAf32;

test {
    _ = @import("_test/conversions.zig");
}
