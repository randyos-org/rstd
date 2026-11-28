//! `minres`: SciPy's case, indefinite real and complex Hermitian systems,
//! preconditioned or not, and the errors.

const std = @import("std");
const check = @import("../../blas/_test/check.zig");
const sparse = @import("../__root__.zig");
const operators = @import("operators.zig");
const types = @import("../../blas/types.zig");

const testing = std.testing;
const Float = types.Float;
const C64 = std.math.Complex(f64);
const Dense = operators.Dense;
const Jacobi = operators.Jacobi;
const tridiag = operators.tridiag;

/// Solves at `tol`, then checks `‖b - A x‖ <= bound ‖b‖`: MINRES's test is
/// against `‖A‖ ‖x‖`, not `‖b‖`.
fn expectMinres(comptime T: type, d: *const Dense(T), m: ?sparse.LinearOperator(T), tol: Float(T), bound: Float(T), b: []const T, x: []T) !sparse.SolveReport {
    const rep = try sparse.minres(testing.allocator, check.io, T, d.op(), m, tol, 1000, b, x);
    try testing.expect(try d.residual(b, x) <= bound * operators.norm(T, b));
    try testing.expect(rep.iterations > 0);
    return rep;
}

test "minres: SciPy's docstring example, indefinite" {
    const a = [_]f64{ 6, 3, 0, 3, -2, 5, 0, 5, 2 };
    const d: Dense(f64) = .{ .n = 3, .a = &a };
    const b = [_]f64{ 2, 4, -1 };
    var x: [3]f64 = @splat(0);
    const rep = try expectMinres(f64, &d, null, 1e-12, 1e-9, &b, &x);
    try testing.expect(rep.iterations <= 4);
}

test "minres: an indefinite Laplacian, f64 and f32, and complex Hermitian" {
    const n = 50;
    var a: [n * n]f64 = undefined;
    var diag: [n]f64 = undefined;
    // eigenvalues `0.9 - 2 cos(k pi / 51)`: both sides of 0, none at it
    tridiag(f64, &a, &diag, n, -1, 0.9, -1);
    const d: Dense(f64) = .{ .n = n, .a = &a };
    var b: [n]f64 = undefined;
    for (&b, 0..) |*bi, i| bi.* = @sin(@as(f64, @floatFromInt(i)));
    var x: [n]f64 = @splat(0);
    _ = try expectMinres(f64, &d, null, 1e-12, 1e-8, &b, &x);
    // from the solution: the estimate is met at once
    const again = try sparse.minres(testing.allocator, check.io, f64, d.op(), null, 1e-6, 1000, &b, &x);
    try testing.expect(again.iterations <= 1);

    var a32: [n * n]f32 = undefined;
    var d32: [n]f32 = undefined;
    tridiag(f32, &a32, &d32, n, -1, 0.9, -1);
    const dense32: Dense(f32) = .{ .n = n, .a = &a32 };
    var b32: [n]f32 = undefined;
    for (&b32, b) |*o, v| o.* = @floatCast(v);
    var x32: [n]f32 = @splat(0);
    _ = try expectMinres(f32, &dense32, null, 1e-5, 1e-3, &b32, &x32);

    const m = 30;
    var c: [m * m]C64 = undefined;
    var cd: [m]C64 = undefined;
    tridiag(C64, &c, &cd, m, .init(1, -1), .init(4, 0), .init(1, 1));
    for (0..m) |i| {
        const s: f64 = if (i % 2 == 0) 1 else -1;
        c[i * m + i] = .init(4 * s, 0);
        cd[i] = .init(4, 0);
    }
    const dc: Dense(C64) = .{ .n = m, .a = &c };
    const jacobi: Jacobi(C64) = .{ .diag = &cd };
    var bc: [m]C64 = undefined;
    for (&bc, 0..) |*bi, i| bi.* = .init(@cos(@as(f64, @floatFromInt(i))), 1);
    var xc: [m]C64 = @splat(.init(0, 0));
    _ = try expectMinres(C64, &dc, null, 1e-12, 1e-8, &bc, &xc);
    @memset(&xc, .init(0, 0));
    _ = try expectMinres(C64, &dc, jacobi.op(), 1e-12, 1e-8, &bc, &xc);
}

test "minres: a zero b, DimensionMismatch, NoConvergence, NotPositiveDefinite" {
    const al = testing.allocator;
    const n = 40;
    var a: [n * n]f64 = undefined;
    var diag: [n]f64 = undefined;
    tridiag(f64, &a, &diag, n, -1, 2, -1);
    const lap: Dense(f64) = .{ .n = n, .a = &a };
    const ones: [n]f64 = @splat(1);
    var y: [n]f64 = @splat(0);
    try testing.expectError(error.DimensionMismatch, sparse.minres(al, check.io, f64, lap.op(), null, 1e-10, 100, &ones, y[0..3]));
    try testing.expectError(error.NoConvergence, sparse.minres(al, check.io, f64, lap.op(), null, 1e-12, 2, &ones, &y));
    @memset(&y, 0);
    try testing.expectError(error.NoConvergence, sparse.minres(al, check.io, f64, lap.op(), null, 1e-12, 0, &ones, &y));
    const negative: [n]f64 = @splat(-1);
    const bad: Jacobi(f64) = .{ .diag = &negative };
    @memset(&y, 0);
    try testing.expectError(error.NotPositiveDefinite, sparse.minres(al, check.io, f64, lap.op(), bad.op(), 1e-12, 100, &ones, &y));
    const zero: [n]f64 = @splat(0);
    @memset(&y, 5);
    const rep = try sparse.minres(al, check.io, f64, lap.op(), null, 1e-10, 100, &zero, &y);
    try testing.expectEqual(0, rep.iterations);
    try testing.expectEqualSlices(f64, &zero, &y);
}
