//! RGB colors: `RGBA` and `BGRA` bytes, `RGBAf32` floats.

const std = @import("std");
const standards = @import("standards.zig");
const channel = @import("channel.zig");
const hsl = @import("hsl.zig");
const ycbcr = @import("ycbcr.zig");
const clamp = std.math.clamp;
const GammaLawConsts = standards.GammaLawConsts;
const YCbCrConsts = standards.YCbCrConsts;
const gammaEncodeChannel = standards.gammaEncodeChannel;
const gammaDecodeChannel = standards.gammaDecodeChannel;
const byteToFloat = channel.byteToFloat;
const floatToByte = channel.floatToByte;
const swapped = channel.swapped;
const HSLAf32 = hsl.HSLAf32;
const YCbCrAf32 = ycbcr.YCbCrAf32;

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
