//! Vector norms by kind (MATLAB's `norm(v, p)`), as `matrixNorm` takes one
//! for matrices: 1, 2 (`norm`'s, the BLAS's `nrm2`), inf, and a general
//! `p`, scaled by the largest entry so it neither overflows nor underflows
//! early (issue #327).

const std = @import("std");
const vector = @import("../vector.zig");
const reductions = @import("reductions.zig");
const util = @import("util.zig");

const Io = std.Io;
const isComplex = vector.isComplex;
const Real = vector.Real;
const ElementOf = util.ElementOf;
const constOf = util.constOf;
const strided = util.strided;

pub const VectorNorm = union(enum) { one, two, inf, p: f64 };

pub const VectorNormError = reductions.VectorError || error{
    /// `p` is not positive, or NaN.
    OrderNotPositive,
};

/// `x`'s norm of the given kind: `Σ|xᵢ|`, `‖x‖₂`, `max|xᵢ|` or
/// `(Σ|xᵢ|ᵖ)^(1/p)`, a complex `|xᵢ|` its modulus; 0 when empty, NaN when
/// `x` holds a NaN.
pub fn vectorNorm(io: Io, x: anytype, which: VectorNorm) VectorNormError!Real(ElementOf(@TypeOf(x))) {
    const T = ElementOf(@TypeOf(x));
    const F = Real(T);
    const kind: VectorNorm = switch (which) {
        .p => |p| if (!(p > 0)) return error.OrderNotPositive else if (p == 1) .one else if (p == 2) .two else if (std.math.isInf(p)) .inf else which,
        else => which,
    };
    if (kind == .two) return reductions.norm(io, x);
    const s = try strided(constOf(x));
    var big: F = 0;
    for (0..s.n) |i| {
        const a = abs(T, s.ptr[i * s.inc]);
        if (a != a) return a;
        big = @max(big, a);
    }
    switch (kind) {
        .inf => return big,
        .one => {
            var sum: F = 0;
            for (0..s.n) |i| sum += abs(T, s.ptr[i * s.inc]);
            return sum;
        },
        .p => |p| {
            if (big == 0 or std.math.isInf(big)) return big;
            const pf: F = @floatCast(p);
            var sum: F = 0;
            for (0..s.n) |i| sum += std.math.pow(F, abs(T, s.ptr[i * s.inc]) / big, pf);
            return big * std.math.pow(F, sum, 1 / pf);
        },
        .two => unreachable,
    }
}

fn abs(comptime T: type, x: T) Real(T) {
    return if (comptime isComplex(T)) std.math.hypot(x.re, x.im) else @abs(x);
}

// ---- tests ----

const testing = std.testing;

test "vectorNorm: each kind, scaled p, special values" {
    const io = testing.io;
    const v = [_]f64{ 3, -4 };
    try testing.expectEqual(7, try vectorNorm(io, @as([]const f64, &v), .one));
    try testing.expectEqual(5, try vectorNorm(io, @as([]const f64, &v), .two));
    try testing.expectEqual(4, try vectorNorm(io, @as([]const f64, &v), .inf));
    try testing.expectApproxEqRel(@as(f64, std.math.pow(f64, 27 + 64, 1.0 / 3.0)), try vectorNorm(io, @as([]const f64, &v), .{ .p = 3 }), 1e-15);
    // no overflow at 1e300 with p = 3
    const huge = [_]f64{ 1e300, 1e300 };
    try testing.expectApproxEqRel(1e300 * std.math.pow(f64, 2, 1.0 / 3.0), try vectorNorm(io, @as([]const f64, &huge), .{ .p = 3 }), 1e-14);
    try testing.expectEqual(5, try vectorNorm(io, @as([]const f64, &v), .{ .p = 2 }));
    const C = std.math.Complex(f64);
    const z = [_]C{ .init(3, 4), .init(0, 1) };
    try testing.expectEqual(6, try vectorNorm(io, @as([]const C, &z), .one));
    try testing.expectError(error.OrderNotPositive, vectorNorm(io, @as([]const f64, &v), .{ .p = 0 }));
    const nan = [_]f64{ 1, std.math.nan(f64) };
    try testing.expect(std.math.isNan(try vectorNorm(io, @as([]const f64, &nan), .inf)));
}
