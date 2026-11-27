//! The Private DICT values the CFF hinter and loader use, parsed as FreeType
//! 2.14.3 parses them (`src/cff/cffparse.c` `cff_parse_real` and the
//! `cfftoken.h` table, `BlueScale` scaled by 1000), after fontations
//! `read-fonts/src/ps/cff/dict.rs`, `num.rs` (`BcdComponents`) and
//! `hinting.rs` (commit 4f115482, MIT). Values are 16.16.

const std = @import("std");
const fixed = @import("../fixed.zig");

pub const Error = error{InvalidDict};

/// Up to seven zones, each a bottom and a top edge (`Blues`).
pub const Blues = struct {
    values: [7][2]i32 = undefined,
    len: usize = 0,

    pub fn slice(b: *const Blues) []const [2]i32 {
        return b.values[0..b.len];
    }

    /// Pairs from `v`, which is already summed from deltas.
    fn init(v: []const i32) Blues {
        var b: Blues = .{};
        var i: usize = 0;
        while (i + 1 < v.len and b.len < 7) : (i += 2) {
            b.values[b.len] = .{ v[i], v[i + 1] };
            b.len += 1;
        }
        return b;
    }
};

pub const Private = struct {
    blues: Blues = .{},
    family_blues: Blues = .{},
    other_blues: Blues = .{},
    family_other_blues: Blues = .{},
    /// 0.039625, the spec's default.
    blue_scale: i32 = fromF64(0.039625),
    blue_shift: i32 = 7 << 16,
    blue_fuzz: i32 = 1 << 16,
    language_group: i32 = 0,
    /// The local Subrs' offset from the Private DICT's start.
    subrs: ?usize = null,
    default_width: ?i32 = null,
    nominal_width: i32 = 0,
};

fn fromF64(comptime v: f64) i32 {
    return @intFromFloat(@round(v * 65536.0));
}

const max_stack = 48;

/// A number on the DICT stack: an integer, or 16.16 from a real.
const Operand = struct {
    v: i32,
    real: bool,
    /// The real's digits, for BlueScale's rescale.
    bcd: ?Bcd = null,

    fn fixed(o: Operand) i32 {
        return if (o.real) o.v else o.v *% 65536;
    }
};

pub fn parse(data: []const u8) Error!Private {
    var p: Private = .{};
    var stack: [max_stack]Operand = undefined;
    var n: usize = 0;
    var at: usize = 0;
    while (at < data.len) {
        const b0 = data[at];
        at += 1;
        if (b0 == 30) {
            const bcd = try Bcd.parse(data, &at);
            if (n < max_stack) {
                stack[n] = .{ .v = bcd.value(false), .real = true, .bcd = bcd };
                n += 1;
            }
            continue;
        }
        if (b0 == 28 or b0 == 29 or b0 >= 32) {
            const v = try integer(data, &at, b0);
            if (n < max_stack) {
                stack[n] = .{ .v = v, .real = false };
                n += 1;
            }
            continue;
        }
        var op: u16 = b0;
        if (b0 == 12) {
            if (at >= data.len) return error.InvalidDict;
            op = 0x0C00 | @as(u16, data[at]);
            at += 1;
        }
        const args = stack[0..n];
        switch (op) {
            6, 7, 8, 9 => {
                var vals: [14]i32 = undefined;
                var sum: i32 = 0;
                // FreeType keeps 14 numbers of the blues, 10 of the others
                const k = @min(args.len, @as(usize, if (op == 7 or op == 9) 10 else 14));
                for (args[0..k], 0..) |a, i| {
                    sum +%= a.fixed();
                    vals[i] = sum;
                }
                const blues = Blues.init(vals[0..k]);
                switch (op) {
                    6 => p.blues = blues,
                    7 => p.other_blues = blues,
                    8 => p.family_blues = blues,
                    else => p.family_other_blues = blues,
                }
            },
            0x0C09 => if (args.len > 0) {
                const a = args[args.len - 1];
                p.blue_scale = if (a.bcd) |b| b.value(true) else a.fixed();
            },
            0x0C0A => if (args.len > 0) {
                p.blue_shift = args[args.len - 1].fixed() & ~@as(i32, 0xFFFF);
            },
            0x0C0B => if (args.len > 0) {
                p.blue_fuzz = args[args.len - 1].fixed() & ~@as(i32, 0xFFFF);
            },
            0x0C11 => if (args.len > 0) {
                p.language_group = args[args.len - 1].fixed() >> 16;
            },
            19 => if (args.len > 0) {
                const v = args[args.len - 1].fixed() >> 16;
                if (v >= 0) p.subrs = @intCast(v);
            },
            20 => if (args.len > 0) {
                p.default_width = args[args.len - 1].fixed() & ~@as(i32, 0xFFFF);
            },
            21 => if (args.len > 0) {
                p.nominal_width = args[args.len - 1].fixed() & ~@as(i32, 0xFFFF);
            },
            else => {},
        }
        n = 0;
    }
    return p;
}

/// An integer operand: `28` (16 bits), `29` (32 bits) or the short forms.
fn integer(data: []const u8, at: *usize, b0: u8) Error!i32 {
    switch (b0) {
        32...246 => return @as(i32, b0) - 139,
        247...250 => {
            if (at.* >= data.len) return error.InvalidDict;
            defer at.* += 1;
            return (@as(i32, b0) - 247) * 256 + data[at.*] + 108;
        },
        251...254 => {
            if (at.* >= data.len) return error.InvalidDict;
            defer at.* += 1;
            return -(@as(i32, b0) - 251) * 256 - @as(i32, data[at.*]) - 108;
        },
        28 => {
            if (at.* + 2 > data.len) return error.InvalidDict;
            defer at.* += 2;
            return std.mem.readInt(i16, data[at.*..][0..2], .big);
        },
        29 => {
            if (at.* + 4 > data.len) return error.InvalidDict;
            defer at.* += 4;
            return std.mem.readInt(i32, data[at.*..][0..4], .big);
        },
        else => return error.InvalidDict,
    }
}

const power_tens = [_]i32{ 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000 };
const number_limit: i32 = 0xCCCCCCC;
const integer_limit: i32 = 0x7FFF;
const overflow: i32 = 0x7FFFFFFF;

/// A real's binary-coded decimal digits (`BcdComponents`).
const Bcd = struct {
    err: ?i32 = null,
    number: i32 = 0,
    sign: i32 = 1,
    exponent: i32 = 0,
    exponent_add: i32 = 0,
    integer_len: i32 = 0,
    fraction_len: i32 = 0,

    fn parse(data: []const u8, at: *usize) Error!Bcd {
        const Phase = enum { integer, fraction, exponent };
        var phase: Phase = .integer;
        var b: Bcd = .{};
        var exp_sign: i32 = 1;
        var too_big = false;
        outer: while (true) {
            if (at.* >= data.len) return error.InvalidDict;
            const byte = data[at.*];
            at.* += 1;
            for ([2]u8{ byte >> 4, byte & 0xF }) |nib| {
                switch (phase) {
                    .integer => switch (nib) {
                        0...9 => if (b.number >= number_limit) {
                            b.exponent_add += 1;
                        } else if (nib != 0 or b.number != 0) {
                            b.number = b.number * 10 + nib;
                            b.integer_len += 1;
                        },
                        0xE => b.sign = -1,
                        0xA => phase = .fraction,
                        0xB => phase = .exponent,
                        0xC => {
                            phase = .exponent;
                            exp_sign = -1;
                        },
                        else => break :outer,
                    },
                    .fraction => switch (nib) {
                        0...9 => if (nib == 0 and b.number == 0) {
                            b.exponent_add -= 1;
                        } else if (b.number < number_limit and b.fraction_len < 9) {
                            b.number = b.number * 10 + nib;
                            b.fraction_len += 1;
                        },
                        0xB => phase = .exponent,
                        0xC => {
                            phase = .exponent;
                            exp_sign = -1;
                        },
                        else => break :outer,
                    },
                    .exponent => switch (nib) {
                        0...9 => if (b.exponent > 1000) {
                            too_big = true;
                        } else {
                            b.exponent = b.exponent * 10 + nib;
                        },
                        else => break :outer,
                    },
                }
            }
        }
        if (too_big and b.number != 0) return .{ .err = if (exp_sign == -1) 0 else overflow *% b.sign };
        b.exponent *= exp_sign;
        return b;
    }

    /// The value in 16.16, scaled down by 1000 for `BlueScale`.
    fn value(b: Bcd, by_1000: bool) i32 {
        if (b.err) |e| return e;
        var number = b.number;
        if (number == 0) return 0;
        const exponent = b.exponent + b.exponent_add + @as(i32, if (by_1000) 3 else 0);
        const integer_len = b.integer_len + exponent;
        var fraction_len = b.fraction_len - exponent;
        if (integer_len > 5) return overflow;
        if (integer_len < -5) return 0;
        if (integer_len < 0) {
            number = @divTrunc(number, power_tens[@intCast(-integer_len)]);
            fraction_len += integer_len;
        }
        if (fraction_len == 10) {
            number = @divTrunc(number, 10);
            fraction_len -= 1;
        }
        var result: i32 = undefined;
        if (fraction_len > 0) {
            const d = power_tens[@intCast(fraction_len)];
            result = if (@divTrunc(number, d) > integer_limit) 0 else fixed.div(number, d);
        } else {
            number *%= power_tens[@intCast(-fraction_len)];
            if (number > integer_limit) return overflow;
            result = number << 16;
        }
        if (by_1000) result = fixed.div(result, 1000 << 16);
        return result *% b.sign;
    }
};

test "cff private: a BlueScale real, deltas into zones, widths floored" {
    // BlueValues -20 20 480 20 (deltas), BlueScale 0.039625 as BCD, nominalWidthX 500.5
    const data = [_]u8{ 139 - 20, 139 + 20, 248, 116, 139 + 20, 6, 30, 0x0A, 0x03, 0x96, 0x25, 0xFF, 12, 9, 30, 0x50, 0x0A, 0x5F, 21 };
    const p = try parse(&data);
    try std.testing.expectEqual(@as(usize, 2), p.blues.len);
    try std.testing.expectEqual([2]i32{ -20 << 16, 0 }, p.blues.values[0]);
    try std.testing.expectEqual([2]i32{ 480 << 16, 500 << 16 }, p.blues.values[1]);
    try std.testing.expect(@abs(p.blue_scale - fromF64(0.039625)) <= 1);
    try std.testing.expectEqual(@as(i32, 500 << 16), p.nominal_width);
}
