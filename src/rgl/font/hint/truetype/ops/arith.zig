//! Arithmetic, logic and rounding: FreeType 2.14.3 `ttinterp.c` `Ins_ADD`
//! to `Ins_MIN`, `Ins_LT` to `Ins_NOT`, `Ins_ODD`, `Ins_EVEN`, `Ins_ROUND`.

const fixed = @import("../../fixed.zig");
const Engine = @import("../Engine.zig");

const Error = Engine.Error;

pub const Binary = enum { add, sub, div, mul, max, min, lt, lteq, gt, gteq, eq, neq, @"and", @"or" };

/// `args[0] op args[1]` into `args[0]`.
pub fn binary(e: *Engine, comptime op: Binary) Error!void {
    const a = e.arg(0).*;
    const b = e.arg(1).*;
    e.arg(0).* = switch (op) {
        .add => a +% b,
        .sub => a -% b,
        .div => if (b == 0) return error.DivideByZero else fixed.mulDivNoRound(a, 64, b),
        .mul => fixed.mulDiv(a, b, 64),
        .max => @max(a, b),
        .min => @min(a, b),
        .lt => @intFromBool(a < b),
        .lteq => @intFromBool(a <= b),
        .gt => @intFromBool(a > b),
        .gteq => @intFromBool(a >= b),
        .eq => @intFromBool(a == b),
        .neq => @intFromBool(a != b),
        .@"and" => @intFromBool(a != 0 and b != 0),
        .@"or" => @intFromBool(a != 0 or b != 0),
    };
}

pub const Unary = enum { abs, neg, floor, ceiling, not, odd, even, round };

pub fn unary(e: *Engine, comptime op: Unary) void {
    const a = e.arg(0).*;
    e.arg(0).* = switch (op) {
        .abs => if (a < 0) 0 -% a else a,
        .neg => 0 -% a,
        .floor => fixed.floor(a),
        .ceiling => fixed.ceil(a),
        .not => @intFromBool(a == 0),
        .odd => @intFromBool(e.roundDistance(a) & 64 == 64),
        .even => @intFromBool(e.roundDistance(a) & 64 == 0),
        .round => e.roundDistance(a),
    };
}
