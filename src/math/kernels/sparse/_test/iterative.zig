//! `cg` and `gmres`: faer's and SciPy's cases, then Laplacians (`bicgstab`
//! too), complex Hermitian and nonsymmetric systems, preconditioned or not,
//! and the errors.

const std = @import("std");
const check = @import("../../blas/_test/check.zig");
const ops = @import("../../lapack/ops.zig");
const types = @import("../../blas/types.zig");
const sparse = @import("../__root__.zig");
const operators = @import("operators.zig");

const testing = std.testing;
const Io = std.Io;
const Float = types.Float;
const C64 = std.math.Complex(f64);
const Dense = operators.Dense;
const Jacobi = operators.Jacobi;
const Solver = operators.Solver;
const norm = operators.norm;
const tridiag = operators.tridiag;
const value = operators.value;
const expectSolves = operators.expectSolves;

test "cg: faer's test_cg and codeberg_326" {
    const al = testing.allocator;
    const a = [_]f64{ 2.5, -1.0, -1.0, 3.1 };
    const d: Dense(f64) = .{ .n = 2, .a = &a };
    const sol = [_]f64{ 2.1, 4.1 };
    var b: [2]f64 = undefined;
    try Dense(f64).apply(&d, check.io, &sol, &b);
    var x: [2]f64 = @splat(0);
    const tol = std.math.floatEps(f64) * 128;
    const rep = try sparse.cg(al, check.io, f64, d.op(), null, tol, 10, &b, &x);
    try testing.expect(rep.iterations <= 2);
    try testing.expect(try d.residual(&b, &x) <= tol * norm(f64, &b));

    const one = [_]f64{2};
    const d1: Dense(f64) = .{ .n = 1, .a = &one };
    var x1 = [_]f64{0};
    _ = try sparse.cg(al, check.io, f64, d1.op(), null, tol, 10, &[_]f64{2}, &x1);
    try testing.expectApproxEqAbs(1.0, x1[0], 1e-10);
}

test "gmres: SciPy's docstring example" {
    const a = [_]f64{ 3, 1, 0, 2, -1, 5, 0, 0, 1 };
    const d: Dense(f64) = .{ .n = 3, .a = &a };
    const b = [_]f64{ 2, 4, -1 };
    var x: [3]f64 = @splat(0);
    _ = try expectSolves(f64, .{ .gmres = 10 }, &d, null, 1e-10, &b, &x);
    try testing.expectApproxEqAbs(2.0, x[0], 1e-8);
    try testing.expectApproxEqAbs(-2.0, x[1], 1e-8);
    try testing.expectApproxEqAbs(9.0, x[2], 1e-8);
}

fn laplacian(comptime T: type, solver: Solver, tol: Float(T)) !void {
    const n = 50;
    var a: [n * n]T = undefined;
    var diag: [n]T = undefined;
    // the 1D Laplacian, scaled down the diagonal so Jacobi has work to do
    tridiag(T, &a, &diag, n, value(T, -1, 0), value(T, 2, 0), value(T, -1, 0));
    for (0..n) |i| {
        const s = value(T, 1 + @as(f64, @floatFromInt(i)) / 10, 0);
        for (0..n) |j| {
            a[j * n + i] = ops.mul(T, a[j * n + i], s);
            a[i * n + j] = ops.mul(T, a[i * n + j], s);
        }
        diag[i] = a[i * n + i];
    }
    const d: Dense(T) = .{ .n = n, .a = &a };
    const jacobi: Jacobi(T) = .{ .diag = &diag };
    var b: [n]T = undefined;
    for (&b, 0..) |*bi, i| bi.* = value(T, @sin(@as(f64, @floatFromInt(i))), 0);
    var x: [n]T = @splat(ops.zero(T));
    const plain = try expectSolves(T, solver, &d, null, tol, &b, &x);
    @memset(&x, ops.zero(T));
    _ = try expectSolves(T, solver, &d, jacobi.op(), tol, &b, &x);
    // from a nonzero guess, and from the solution itself
    for (&x, 0..) |*xi, i| xi.* = value(T, @floatFromInt(i % 3), 0);
    _ = try expectSolves(T, solver, &d, jacobi.op(), tol, &b, &x);
    const again = switch (solver) {
        .cg => try sparse.cg(testing.allocator, check.io, T, d.op(), null, tol, 1000, &b, &x),
        .gmres => try sparse.gmres(testing.allocator, check.io, T, d.op(), null, tol, 1000, 0, &b, &x),
        .bicgstab => try sparse.bicgstab(testing.allocator, check.io, T, d.op(), null, tol, 1000, &b, &x),
    };
    try testing.expect(again.iterations <= 1);
    try testing.expect(plain.iterations > 1);
    // unrestarted, GMRES ends within `n` iterations
    if (solver == .gmres) try testing.expect(plain.iterations <= n);
}

test "cg, gmres and bicgstab: a Laplacian, f64 and f32" {
    try laplacian(f64, .cg, 1e-10);
    try laplacian(f32, .cg, 1e-4);
    try laplacian(f64, .{ .gmres = 50 }, 1e-10);
    try laplacian(f32, .{ .gmres = 50 }, 1e-4);
    try laplacian(f64, .bicgstab, 1e-10);
    try laplacian(f32, .bicgstab, 1e-4);
}

test "cg and gmres: complex Hermitian, gmres nonsymmetric" {
    const n = 30;
    var a: [n * n]C64 = undefined;
    var diag: [n]C64 = undefined;
    tridiag(C64, &a, &diag, n, .init(1, -1), .init(4, 0), .init(1, 1));
    const d: Dense(C64) = .{ .n = n, .a = &a };
    var b: [n]C64 = undefined;
    for (&b, 0..) |*bi, i| bi.* = .init(@cos(@as(f64, @floatFromInt(i))), 1);
    var x: [n]C64 = @splat(.init(0, 0));
    _ = try expectSolves(C64, .cg, &d, null, 1e-10, &b, &x);
    @memset(&x, .init(0, 0));
    _ = try expectSolves(C64, .{ .gmres = 10 }, &d, null, 1e-10, &b, &x);

    // convection-diffusion: nonsymmetric, complex and real
    tridiag(C64, &a, &diag, n, .init(-1.3, 0.2), .init(2, 0.5), .init(-0.7, 0));
    const jacobi: Jacobi(C64) = .{ .diag = &diag };
    @memset(&x, .init(0, 0));
    _ = try expectSolves(C64, .{ .gmres = 10 }, &d, jacobi.op(), 1e-10, &b, &x);
    var ar: [n * n]f64 = undefined;
    var dr: [n]f64 = undefined;
    tridiag(f64, &ar, &dr, n, -1.3, 2, -0.7);
    const dense: Dense(f64) = .{ .n = n, .a = &ar };
    var br: [n]f64 = undefined;
    for (&br, 0..) |*bi, i| bi.* = @floatFromInt(i % 5);
    var xr: [n]f64 = @splat(0);
    _ = try expectSolves(f64, .{ .gmres = 10 }, &dense, null, 1e-10, &br, &xr);
}

test "cg and gmres: a zero b, DimensionMismatch, NoConvergence, NotPositiveDefinite" {
    const al = testing.allocator;
    const neg = [_]f64{ -1, 0, 0, -2 };
    const d: Dense(f64) = .{ .n = 2, .a = &neg };
    const b = [_]f64{ 1, 2 };
    var x = [_]f64{ 5, 5 };
    try testing.expectError(error.DimensionMismatch, sparse.cg(al, check.io, f64, d.op(), null, 1e-10, 100, &b, x[0..1]));
    try testing.expectError(error.DimensionMismatch, sparse.gmres(al, check.io, f64, d.op(), null, 1e-10, 100, 20, &b, x[0..1]));
    try testing.expectError(error.NotPositiveDefinite, sparse.cg(al, check.io, f64, d.op(), null, 1e-10, 100, &b, &x));
    // gmres needs no definiteness
    @memset(&x, 0);
    _ = try expectSolves(f64, .{ .gmres = 10 }, &d, null, 1e-12, &b, &x);

    const zero = [_]f64{ 0, 0 };
    x = .{ 5, 5 };
    const rep = try sparse.cg(al, check.io, f64, d.op(), null, 1e-10, 100, &zero, &x);
    try testing.expectEqual(0, rep.iterations);
    try testing.expectEqualSlices(f64, &zero, &x);
    x = .{ 5, 5 };
    _ = try sparse.gmres(al, check.io, f64, d.op(), null, 1e-10, 100, 20, &zero, &x);
    try testing.expectEqualSlices(f64, &zero, &x);

    const n = 40;
    var a: [n * n]f64 = undefined;
    var diag: [n]f64 = undefined;
    tridiag(f64, &a, &diag, n, -1, 2, -1);
    const lap: Dense(f64) = .{ .n = n, .a = &a };
    const ones: [n]f64 = @splat(1);
    var y: [n]f64 = @splat(0);
    try testing.expectError(error.NoConvergence, sparse.cg(al, check.io, f64, lap.op(), null, 1e-12, 3, &ones, &y));
    @memset(&y, 0);
    try testing.expectError(error.NoConvergence, sparse.gmres(al, check.io, f64, lap.op(), null, 1e-12, 3, 20, &ones, &y));
    @memset(&y, 0);
    try testing.expectError(error.NoConvergence, sparse.gmres(al, check.io, f64, lap.op(), null, 1e-12, 0, 20, &ones, &y));
}
