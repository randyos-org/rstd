//! `bicgstab`: faer's case, complex and real nonsymmetric systems,
//! preconditioned or not, and the errors.

const std = @import("std");
const check = @import("../../blas/_test/check.zig");
const sparse = @import("../__root__.zig");
const operators = @import("operators.zig");

const testing = std.testing;
const C64 = std.math.Complex(f64);
const Dense = operators.Dense;
const Jacobi = operators.Jacobi;
const tridiag = operators.tridiag;
const expectSolves = operators.expectSolves;

test "bicgstab: faer's test_bicgstab" {
    const a = [_]f64{ 2.5, 1.0, -1.0, 3.1 };
    const d: Dense(f64) = .{ .n = 2, .a = &a };
    const sol = [_]f64{ 2.1, 4.1 };
    var b: [2]f64 = undefined;
    try Dense(f64).apply(&d, check.io, &sol, &b);
    const diag = [_]f64{ 1.7, 2.3 };
    const jacobi: Jacobi(f64) = .{ .diag = &diag };
    var x: [2]f64 = @splat(0);
    const tol = std.math.floatEps(f64) * 128;
    const rep = try sparse.bicgstab(testing.allocator, check.io, f64, d.op(), jacobi.op(), tol, 10, &b, &x);
    try testing.expect(rep.iterations <= 2);
    try testing.expect(try d.residual(&b, &x) <= tol * operators.norm(f64, &b));
}

test "bicgstab: complex Hermitian and nonsymmetric, real nonsymmetric" {
    const n = 30;
    var a: [n * n]C64 = undefined;
    var diag: [n]C64 = undefined;
    tridiag(C64, &a, &diag, n, .init(1, -1), .init(4, 0), .init(1, 1));
    const d: Dense(C64) = .{ .n = n, .a = &a };
    var b: [n]C64 = undefined;
    for (&b, 0..) |*bi, i| bi.* = .init(@cos(@as(f64, @floatFromInt(i))), 1);
    var x: [n]C64 = @splat(.init(0, 0));
    _ = try expectSolves(C64, .bicgstab, &d, null, 1e-10, &b, &x);

    tridiag(C64, &a, &diag, n, .init(-1.3, 0.2), .init(2, 0.5), .init(-0.7, 0));
    const jacobi: Jacobi(C64) = .{ .diag = &diag };
    @memset(&x, .init(0, 0));
    _ = try expectSolves(C64, .bicgstab, &d, null, 1e-10, &b, &x);
    @memset(&x, .init(0, 0));
    _ = try expectSolves(C64, .bicgstab, &d, jacobi.op(), 1e-10, &b, &x);

    var ar: [n * n]f64 = undefined;
    var dr: [n]f64 = undefined;
    tridiag(f64, &ar, &dr, n, -1.3, 2.5, -0.7);
    const dense: Dense(f64) = .{ .n = n, .a = &ar };
    var br: [n]f64 = undefined;
    for (&br, 0..) |*bi, i| bi.* = @floatFromInt(i % 5);
    var xr: [n]f64 = @splat(0);
    _ = try expectSolves(f64, .bicgstab, &dense, null, 1e-10, &br, &xr);
    // from the solution: no iteration
    const again = try sparse.bicgstab(testing.allocator, check.io, f64, dense.op(), null, 1e-8, 100, &br, &xr);
    try testing.expectEqual(0, again.iterations);
}

test "bicgstab: a zero b, DimensionMismatch, NoConvergence" {
    const al = testing.allocator;
    const n = 40;
    var a: [n * n]f64 = undefined;
    var diag: [n]f64 = undefined;
    tridiag(f64, &a, &diag, n, -1, 2, -1);
    const lap: Dense(f64) = .{ .n = n, .a = &a };
    const ones: [n]f64 = @splat(1);
    var y: [n]f64 = @splat(0);
    try testing.expectError(error.DimensionMismatch, sparse.bicgstab(al, check.io, f64, lap.op(), null, 1e-10, 100, &ones, y[0..3]));
    try testing.expectError(error.NoConvergence, sparse.bicgstab(al, check.io, f64, lap.op(), null, 1e-12, 2, &ones, &y));
    @memset(&y, 0);
    try testing.expectError(error.NoConvergence, sparse.bicgstab(al, check.io, f64, lap.op(), null, 1e-12, 0, &ones, &y));
    const zero: [n]f64 = @splat(0);
    @memset(&y, 5);
    const rep = try sparse.bicgstab(al, check.io, f64, lap.op(), null, 1e-10, 100, &zero, &y);
    try testing.expectEqual(0, rep.iterations);
    try testing.expectEqualSlices(f64, &zero, &y);
}
