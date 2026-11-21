//! The 1-norm reciprocal condition estimate (MATLAB's `rcond`): `getrf`,
//! `lange` and `gecon`, the cheap companion to `cond`'s SVD (issue #324).

const std = @import("std");
const lapack = @import("../../math/kernels/lapack/__root__.zig");
const matrix = @import("../matrix.zig");
const vector = @import("../vector.zig");
const ops = @import("../../math/kernels/lapack/ops.zig");
const util = @import("util.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Real = vector.Real;
const ElementOf = util.ElementOf;
const constOf = util.constOf;
const copyOf = util.copyOf;

pub const RcondError = error{NotSquare} || Allocator.Error || Io.Cancelable;

/// An estimate of `1 / (‖A‖₁ ‖A⁻¹‖₁)`: near 1 for a well-conditioned `A`,
/// near `eps` for a badly conditioned one, 0 when an LU pivot is exactly 0.
/// As MATLAB's: infinite for an empty `A`, NaN when `A` holds a NaN, 0 when
/// it holds an infinity.
pub fn rcond(allocator: Allocator, io: Io, a: anytype) RcondError!Real(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const F = Real(T);
    const ca = constOf(a);
    if (ca.rows != ca.cols) return error.NotSquare;
    const n = ca.rows;
    if (n == 0) return std.math.inf(F);
    const f = try copyOf(T, allocator, ca);
    defer f.deinit(allocator);
    const anorm = lapack.lange(allocator, io, T, .o, n, n, f.data, n) catch |e| return narrow(e);
    if (std.math.isNan(anorm)) return std.math.nan(F);
    if (std.math.isInf(anorm)) return 0;
    const ipiv = try allocator.alloc(usize, n);
    defer allocator.free(ipiv);
    if (lapack.getrf(allocator, io, T, n, n, f.data, n, ipiv) catch |e| return narrow(e)) |_| return 0;
    return lapack.gecon(allocator, io, T, .o, n, f.data, n, anorm) catch |e| switch (e) {
        // an LU too ill-scaled to estimate: as singular
        error.EstimateFailed => 0,
        else => return narrow(e),
    };
}

fn narrow(e: anyerror) RcondError {
    return ops.narrow(RcondError, e);
}

// ---- tests ----

const testing = std.testing;

test "rcond: identity, singular, ill-conditioned, special values" {
    const al = testing.allocator;
    const io = testing.io;
    const V = matrix.ConstView(f64);
    const id = [_]f64{ 1, 0, 0, 1 };
    try testing.expectEqual(1, try rcond(al, io, try V.init(&id, 2, 2, 2)));
    const sing = [_]f64{ 1, 2, 2, 4 };
    try testing.expectEqual(0, try rcond(al, io, try V.init(&sing, 2, 2, 2)));
    // [1 1; 1 1 + 1e-10]: 1-norm condition about 4e10
    const ill = [_]f64{ 1, 1, 1, 1 + 1e-10 };
    const r = try rcond(al, io, try V.init(&ill, 2, 2, 2));
    try testing.expect(r > 1e-12 and r < 1e-9);
    try testing.expect(std.math.isInf(try rcond(al, io, try V.init(&.{}, 0, 0, 1))));
    const nan = [_]f64{ std.math.nan(f64), 0, 0, 1 };
    try testing.expect(std.math.isNan(try rcond(al, io, try V.init(&nan, 2, 2, 2))));
    try testing.expectError(error.NotSquare, rcond(al, io, try V.init(&id, 1, 2, 1)));
}
