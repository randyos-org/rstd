//! A CFF DICT's operands and operators, read as FreeType 2.14.3 reads them
//! (`src/cff/cffparse.c`: `cff_parse_integer`, `cff_parse_real` with and
//! without its scaling, `cff_parse_fixed_dynamic`), after fontations
//! `read-fonts/src/ps/cff/dict.rs` and `num.rs` (`BcdComponents`, commit
//! 4f115482, MIT). The Private DICT (`private.zig`) and FontMatrix
//! (`FontMatrix.zig`) read through it.

const std = @import("std");
const fixed_point = @import("../fixed.zig");

pub const Error = error{InvalidDict};

const max_stack = 48;

/// A number on the DICT stack: an integer, or 16.16 from a real.
pub const Operand = struct {
    v: i32,
    real: bool,
    /// The real's digits, for rescaling.
    bcd: ?Bcd = null,

    pub fn fixed(o: Operand) i32 {
        return if (o.real) o.v else o.v *% 65536;
    }

    /// `cff_parse_fixed_dynamic`: the number as precisely as 16.16 holds
    /// it, and the power of ten it was divided by.
    pub fn dynamic(o: Operand) struct { i32, i32 } {
        if (o.bcd) |b| return b.scaled();
        const number = o.v;
        if (number <= 0x7FFF) return .{ number << 16, 0 };
        var len: usize = 5;
        while (len < 10) : (len += 1) {
            if (number < power_tens[len]) break;
        }
        if (@divTrunc(number, power_tens[len - 5]) > 0x7FFF) return .{ fixed_point.div(number, power_tens[len - 4]), @intCast(len - 4) };
        return .{ fixed_point.div(number, power_tens[len - 5]), @intCast(len - 5) };
    }
};

/// The DICT's operators in order, each with its operands.
pub const Parser = struct {
    data: []const u8,
    at: usize = 0,
    stack: [max_stack]Operand = undefined,

    pub const Op = struct { op: u16, args: []const Operand };

    /// The next operator; null at the end. Operands past FreeType's stack
    /// depth are dropped.
    pub fn next(p: *Parser) Error!?Op {
        var n: usize = 0;
        while (p.at < p.data.len) {
            const b0 = p.data[p.at];
            p.at += 1;
            if (b0 == 30) {
                const bcd = try Bcd.parse(p.data, &p.at);
                if (n < max_stack) {
                    p.stack[n] = .{ .v = bcd.value(false), .real = true, .bcd = bcd };
                    n += 1;
                }
                continue;
            }
            if (b0 == 28 or b0 == 29 or b0 >= 32) {
                const v = try integer(p.data, &p.at, b0);
                if (n < max_stack) {
                    p.stack[n] = .{ .v = v, .real = false };
                    n += 1;
                }
                continue;
            }
            var op: u16 = b0;
            if (b0 == 12) {
                if (p.at >= p.data.len) return error.InvalidDict;
                op = 0x0C00 | @as(u16, p.data[p.at]);
                p.at += 1;
            }
            return .{ .op = op, .args = p.stack[0..n] };
        }
        return null;
    }
};

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
pub const Bcd = struct {
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
    pub fn value(b: Bcd, by_1000: bool) i32 {
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
            result = if (@divTrunc(number, d) > integer_limit) 0 else fixed_point.div(number, d);
        } else {
            number *%= power_tens[@intCast(-fraction_len)];
            if (number > integer_limit) return overflow;
            result = number << 16;
        }
        if (by_1000) result = fixed_point.div(result, 1000 << 16);
        return result *% b.sign;
    }

    /// `cff_parse_real` with its scaling: as many digits as 16.16 holds,
    /// and the power of ten taken out.
    pub fn scaled(b: Bcd) struct { i32, i32 } {
        if (b.err) |e| return .{ e, 0 };
        var number = b.number;
        if (number == 0) return .{ 0, 0 };
        var exponent = b.exponent + b.exponent_add + b.integer_len;
        const fraction_len = b.fraction_len + b.integer_len;
        var result: i32 = undefined;
        var scaling: i32 = undefined;
        if (fraction_len <= 5) {
            if (number > integer_limit) {
                result = fixed_point.div(number, 10);
                scaling = exponent - fraction_len + 1;
            } else {
                // as small a scaling as the digits allow
                const shift = @min(exponent, 5) - fraction_len;
                if (exponent > 0 and shift > 0) {
                    exponent -= @min(exponent, 5);
                    number *= power_tens[@intCast(shift)];
                    if (number > integer_limit) {
                        number = @divTrunc(number, 10);
                        exponent += 1;
                    }
                } else exponent -= fraction_len;
                result = number << 16;
                scaling = exponent;
            }
        } else if (@divTrunc(number, power_tens[@intCast(fraction_len - 5)]) > integer_limit) {
            result = fixed_point.div(number, power_tens[@intCast(fraction_len - 4)]);
            scaling = exponent - 4;
        } else {
            result = fixed_point.div(number, power_tens[@intCast(fraction_len - 5)]);
            scaling = exponent - 5;
        }
        return .{ result *% b.sign, scaling };
    }
};
