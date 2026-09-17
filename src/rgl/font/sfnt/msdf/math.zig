const std = @import("std");
const aux_types = @import("../aux_types.zig");
const linalg = @import("../../../../linalg/__root__.zig");
const rsimd = @import("../../../../simd/__root__.zig");

const vec = linalg.vector;
const zigvec = rsimd.zigvec;
const V2f64 = zigvec.V2f64;

pub const EquationParams = struct {
    a: f64,
    b: f64,
    c: f64,
    /// Keeping this as NaN implies a quadratic equation
    d: f64 = std.math.nan(f64),
};

pub const dot2 = vec.dot2;
pub const length = vec.length;
pub const dot = vec.dot;
pub const cross = vec.cross2;
pub const mix = rsimd.mix;
pub const splatV2f64 = zigvec.splatV2f64;

pub fn ortho(v: V2f64) V2f64 {
    return .{ -v[1], v[0] };
}

fn boolToF64(b: bool) f64 {
    return @floatFromInt(@intFromBool(b));
}

pub fn normal(v: V2f64, disallow_zero: bool) V2f64 {
    return vec.norm(v, .{ 0.0, boolToF64(disallow_zero) });
}

pub fn orthonormal(v: V2f64, disallow_zero: bool) V2f64 {
    const len = length(v);
    if (len == 0.0) return .{ 0.0, boolToF64(disallow_zero) };
    return ortho(v) / splatV2f64(len);
}

pub fn median(a: anytype, b: anytype, c: anytype) @TypeOf(a) {
    return @max(@min(a, b), @min(@max(a, b), c));
}

pub fn boolSign(b: bool) i2 {
    return @intCast(@as(i3, @intFromBool(b)) * 2 - 1);
}

pub fn nonZeroSign(n: anytype) @TypeOf(n) {
    const sign = boolSign(n >= 0);
    switch (@typeInfo(@TypeOf(n))) {
        .float, .comptime_float => return @floatFromInt(sign),
        .int, .comptime_int => return sign,
        else => @compileError("Invalid type, only floats and ints are supported"),
    }
}

/// This populates the input `buf` buffer and returns a slice from it,
/// so you should be mindful of its lifetime.
pub fn solveEquation(buf: []f64, params: *const EquationParams) []f64 {
    const quadratic_specified = std.math.isNan(params.d);
    if (quadratic_specified or params.a == 0.0 or @abs(params.b / params.a) > 1e6) {
        const a, const b, const c = if (quadratic_specified)
            .{ params.a, params.b, params.c }
        else
            .{ params.b, params.c, params.d };

        if (a == 0 or @abs(b) > 1e12 * @abs(a)) {
            if (b == 0) return &.{};
            buf[0] = -c / b;
            return buf[0..1];
        }

        const discriminant = b * b - 4.0 * a * c;
        if (discriminant > 0) {
            const discriminant_sqrt = @sqrt(discriminant);
            buf[0] = (-b + discriminant_sqrt) / (2 * a);
            buf[1] = (-b - discriminant_sqrt) / (2 * a);
            return buf[0..2];
        }

        if (discriminant == 0) {
            buf[0] = -b / (2 * a);
            return buf[0..1];
        }

        return &.{};
    }

    const a = params.b / params.a;
    const b = params.c / params.a;
    const c = params.d / params.a;
    const a_sqr = a * a;
    const q = (a_sqr - 3.0 * b) / 9.0;
    const r = (a * (2.0 * a_sqr - 9.0 * b) + 27.0 * c) / 54.0;
    const r_sqr = r * r;
    const q_cube = q * q * q;
    const a_third = a / 3.0;
    if (r_sqr < q_cube) {
        const t = std.math.acos(std.math.clamp(r / @sqrt(q_cube), -1.0, 1.0));
        const sq = -2 * @sqrt(q);
        buf[0] = sq * @cos(t / 3.0) - a_third;
        buf[1] = sq * @cos((t + std.math.tau) / 3.0) - a_third;
        buf[2] = sq * @cos((t + 2 * std.math.tau) / 3.0) - a_third;
        return buf[0..3];
    }

    const inv: f64 = if (r < 0.0) 1.0 else -1.0;
    const u = inv * std.math.pow(f64, @abs(r) + @sqrt(r_sqr - q_cube), 1.0 / 3.0);
    if (u == 0.0) {
        buf[0] = -a_third;
        buf[1] = -a_third;
        return buf[0..2];
    }

    const v = q / u;
    buf[0] = (u + v) - a_third;
    if (u == v or @abs(u - v) < 1e-12 * @abs(u + v)) {
        buf[1] = -0.5 * (u + v) - a_third;
        return buf[0..2];
    }

    return buf[0..1];
}
