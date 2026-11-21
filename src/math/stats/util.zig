//! Type helpers shared by `rstd.stats`' routines.

const std = @import("std");

pub fn Wide(comptime F: type) type {
    return if (@bitSizeOf(F) < 64) f64 else F;
}

/// The element type of a slice, array or pointer to one.
pub fn Elem(comptime S: type) type {
    return std.meta.Elem(S);
}

/// `std.math.Complex` and other `re`, `im` structs.
pub fn isComplex(comptime T: type) bool {
    return @typeInfo(T) == .@"struct" and @hasField(T, "re") and @hasField(T, "im");
}

pub fn Real(comptime R: type) type {
    return if (isComplex(R)) @FieldType(R, "re") else R;
}

/// A real element as the float `W`.
pub inline fn to(comptime W: type, x: anytype) W {
    return switch (@typeInfo(@TypeOf(x))) {
        .int, .comptime_int => @floatFromInt(x),
        .float, .comptime_float => @floatCast(x),
        else => @compileError("rstd.stats: " ++ @typeName(@TypeOf(x)) ++ " is not a real number type"),
    };
}
/// `R` a float for real `E` (an integer or a float), a `Complex` of a float
/// `R` a float for real `E` (an integer or a float), a `Complex` of a float
/// for complex `E`.
pub fn checkPair(comptime R: type, comptime E: type) void {
    if (isComplex(E)) {
        if (!isComplex(R)) @compileError("rstd.stats: complex " ++ @typeName(E) ++ " needs a Complex result, not " ++ @typeName(R));
        checkFloat(Real(R));
    } else {
        checkReal(E);
        checkFloat(R);
    }
}

pub fn checkReal(comptime E: type) void {
    switch (@typeInfo(E)) {
        .int, .float => {},
        else => @compileError("rstd.stats: " ++ @typeName(E) ++ " is not an integer or float type"),
    }
}

pub fn checkFloat(comptime F: type) void {
    if (@typeInfo(F) != .float) @compileError("rstd.stats: " ++ @typeName(F) ++ " is not a float type");
}
