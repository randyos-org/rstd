//! HSL colors, in bytes and floats.

const std = @import("std");
const channel = @import("channel.zig");
const rgb = @import("rgb.zig");
const math = std.math;
const clamp = std.math.clamp;
const byteToFloat = channel.byteToFloat;
const floatToByte = channel.floatToByte;
const swapped = channel.swapped;
const RGBA = rgb.RGBA;
const RGBAf32 = rgb.RGBAf32;

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
