//! Y'CbCr colors, in bytes and floats.

const std = @import("std");
const standards = @import("standards.zig");
const channel = @import("channel.zig");
const rgb = @import("rgb.zig");
const clamp = std.math.clamp;
const YCbCrConsts = standards.YCbCrConsts;
const YCbCrRange = standards.YCbCrRange;
const rangeConsts = standards.rangeConsts;
const byteToFloat = channel.byteToFloat;
const floatToByte = channel.floatToByte;
const byteToScaledOffset = channel.byteToScaledOffset;
const scaledOffsetToByte = channel.scaledOffsetToByte;
const swapped = channel.swapped;
const RGBA = rgb.RGBA;
const RGBAf32 = rgb.RGBAf32;

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
