const std = @import("std");

const bytes = @import("../../bytes.zig");
const zigstd = @import("../../zigstd.zig");

const getI16BE = bytes.getI16BE;
const setI16BE = bytes.setI16BE;
const getI32BE = bytes.getI32BE;
const setI32BE = bytes.setI32BE;
const getU32BE = bytes.getU32BE;
const setU32BE = bytes.setU32BE;

const testExpectEqual = zigstd.testExpectEqual;
const testExpectEqualSlices = zigstd.testExpectEqualSlices;

pub const ShortFrac = packed struct(i16) {
    const Self = @This();
    const scalar = 1 << 14;

    fraction: u14,
    integer: i2,

    pub fn toFloat(self: Self, comptime T: type) T {
        const raw: i16 = @bitCast(self);
        const denom: T = @floatFromInt(scalar);
        return @as(T, @floatFromInt(raw)) / denom;
    }

    pub fn fromFloat(value: anytype) Self {
        const scaled = value * scalar;
        const clamped = std.math.clamp(scaled, std.math.minInt(i16), std.math.maxInt(i16));
        const raw: i16 = @trunc(clamped);
        return @bitCast(raw);
    }

    pub fn toBytes(self: Self) [2]u8 {
        var buf: [2]u8 = undefined;
        setI16BE(&buf, @as(i16, @bitCast(self)));
        return buf;
    }

    pub fn fromBytes(raw_bytes: []const u8) Self {
        const raw = getI16BE(raw_bytes);
        return @bitCast(raw);
    }
};
pub const F2Dot14 = ShortFrac;

test "ShortFrac: Boundary and conversion testing" {
    // 1. Test 1.0 (0x4000)
    const one_bytes = [_]u8{ 0x40, 0x00 };
    const sf_one = ShortFrac.fromBytes(&one_bytes);
    try std.testing.expectApproxEqAbs(@as(f32, 1.0), sf_one.toFloat(f32), 0.00001);
    try testExpectEqual(@as(i2, 1), sf_one.integer);
    try testExpectEqual(@as(u14, 0), sf_one.fraction);

    // 2. Test Max Value: 1.99993896... (0x7FFF)
    const max_bytes = [_]u8{ 0x7F, 0xFF };
    const sf_max = ShortFrac.fromBytes(&max_bytes);
    try std.testing.expectApproxEqAbs(@as(f32, 1.99993896), sf_max.toFloat(f32), 0.00001);

    // 3. Test 0.0 (0x0000)
    const zero_bytes = [_]u8{ 0x00, 0x00 };
    const sf_zero = ShortFrac.fromBytes(&zero_bytes);
    try std.testing.expectApproxEqAbs(@as(f32, 0.0), sf_zero.toFloat(f32), 0.00001);

    // 4. Test -1.0 (0xC000)
    const neg_one_bytes = [_]u8{ 0xC0, 0x00 };
    const sf_neg_one = ShortFrac.fromBytes(&neg_one_bytes);
    try std.testing.expectApproxEqAbs(@as(f32, -1.0), sf_neg_one.toFloat(f32), 0.00001);
    try testExpectEqual(@as(i2, -1), sf_neg_one.integer);

    // 5. Test Min Value: -2.0 (0x8000)
    const min_bytes = [_]u8{ 0x80, 0x00 };
    const sf_min = ShortFrac.fromBytes(&min_bytes);
    try std.testing.expectApproxEqAbs(@as(f32, -2.0), sf_min.toFloat(f32), 0.00001);
    try testExpectEqual(@as(i2, -2), sf_min.integer);

    // 6. Test float -> shortFrac conversion & serialization cycle
    // 0.25 is exactly 1/4 of 1.0, which means 16384 / 4 = 4096 (0x1000)
    const sf_quarter = ShortFrac.fromFloat(0.25);
    const quarter_bytes = sf_quarter.toBytes();
    const expected_quarter = [_]u8{ 0x10, 0x00 };
    try testExpectEqualSlices(u8, &expected_quarter, &quarter_bytes);
}

pub const Fixed = packed struct(i32) {
    const Self = @This();
    const scalar = 1 << 16;

    fraction: u16,
    integer: i16,

    pub fn toFloat(self: Self, comptime T: type) T {
        const raw: i32 = @bitCast(self);
        const denom: T = @floatFromInt(scalar);
        return @as(T, @floatFromInt(raw)) / denom;
    }

    pub fn fromFloat(value: anytype) Self {
        const scaled = value * scalar;
        const clamped = std.math.clamp(scaled, std.math.minInt(i32), std.math.maxInt(i32));
        const raw: i32 = @trunc(clamped);
        return @bitCast(raw);
    }

    pub fn toBytes(self: Self) [4]u8 {
        var buf: [4]u8 = undefined;
        setI32BE(&buf, @as(i32, @bitCast(self)));
        return buf;
    }

    pub fn fromBytes(raw_bytes: []const u8) Self {
        const raw = getI32BE(raw_bytes);
        return @bitCast(raw);
    }
};

test "TTF Fixed type serialization and math" {
    // 1. Test reading 1.0 (0x00010000)
    const one_bytes = [4]u8{ 0x00, 0x01, 0x00, 0x00 };
    const fixed_one = Fixed.fromBytes(&one_bytes);
    try std.testing.expectApproxEqAbs(1.0, fixed_one.toFloat(f32), 0.00001);

    // 2. Test reading a typical font version like 1.000 (often 0x00010000 in TTF headers)
    try testExpectEqual(@as(i16, 1), fixed_one.integer);
    try testExpectEqual(@as(u16, 0), fixed_one.fraction);

    // 3. Test conversion from float (e.g., -1.5)
    const fixed_neg = Fixed.fromFloat(-1.5);
    const neg_bytes = fixed_neg.toBytes();
    const expected_neg_bytes = [4]u8{ 0xFF, 0xFE, 0x80, 0x00 }; // -1 integer part, 0x8000 fraction part
    try testExpectEqualSlices(u8, &expected_neg_bytes, &neg_bytes);
}

pub const Version16Dot16 = packed struct(u32) {
    minor: u16,
    major: u16,

    pub fn fromBytes(raw_bytes: []const u8) Version16Dot16 {
        const raw = getU32BE(raw_bytes);
        return @bitCast(raw);
    }

    /// Formats the version correctly as a string without float rounding issues.
    /// This handles the Microsoft quirk by treating the minor version as a literal integer.
    pub fn format(
        self: Version16Dot16,
        comptime fmt: []const u8,
        options: std.fmt.FormatOptions,
        writer: anytype,
    ) !void {
        _ = fmt;
        _ = options;

        // If minor is 0x0050, it prints "1.50" or "1.5", not "1.00122"
        if (self.minor == 0) {
            try writer.print("{d}.0", .{self.major});
        } else {
            // Standard OpenType approach: print major and the raw minor integer
            try writer.print("{d}.{d}", .{ self.major, self.minor });
        }
    }
};

test "Parsing Version16Dot16" {
    // Version 1.0 (0x00010000)
    const v1_bytes = [_]u8{ 0x00, 0x01, 0x00, 0x00 };
    const version = Version16Dot16.fromBytes(&v1_bytes);

    try testExpectEqual(@as(u16, 1), version.major);
    try testExpectEqual(@as(u16, 0), version.minor);
}

pub const FWord = i16;
pub const UFWord = u16;
pub const LongDateTime = i64;
pub const Offset16 = u16;
pub const Offset32 = u32;

pub const TagStr = extern struct {
    data: [4]u8,

    pub fn strip(self: *const TagStr) []const u8 {
        var len: usize = 4;
        while (len > 0 and self.data[len - 1] == ' ') : (len -= 1) {}
        return self.data[0..len];
    }

    pub fn fromInt(v: u32) TagStr {
        var tag: TagStr = undefined;
        setU32BE(&tag.data, v);
        return tag;
    }

    pub fn fromStr(s: []const u8) !TagStr {
        if (s.len > 4) return error.StringTooLong;
        var tag: TagStr = .{ .data = @splat(' ') };
        @memcpy(tag.data[0..s.len], s);
        return tag;
    }
};

test "TagStr.fromStr pads short tags with spaces, strip() removes them" {
    const short = try TagStr.fromStr("cvt");
    try testExpectEqualSlices(u8, "cvt ", &short.data);
    try testExpectEqualSlices(u8, "cvt", short.strip());

    const full = try TagStr.fromStr("cmap");
    try testExpectEqualSlices(u8, "cmap", &full.data);
    try testExpectEqualSlices(u8, "cmap", full.strip());
}

test "TagStr.fromInt/fromStr agree on the same 4 bytes" {
    const from_str = try TagStr.fromStr("cmap");
    const from_int = TagStr.fromInt(getU32BE(&from_str.data));
    try testExpectEqualSlices(u8, &from_str.data, &from_int.data);
}

pub const GlyphIndex = enum(u16) {
    notdef = 0,
    _,
};
