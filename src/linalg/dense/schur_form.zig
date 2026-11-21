//! The Hessenberg and Schur forms (MATLAB's `hess` and `schur`): `A = P H
//! Pᴴ` by `gehrd` and `orghr`, `A = U T Uᴴ` by `gees`. A real `A` gets the
//! real Schur form, its 2 x 2 blocks standardized, or with `schurComplex`
//! the complex one. Replaces `matlab/schur.zig` (issue #322).

const std = @import("std");
const lapack = @import("../../math/kernels/lapack/__root__.zig");
const matrix = @import("../matrix.zig");
const vector = @import("../vector.zig");
const ops = @import("../../math/kernels/lapack/ops.zig");
const util = @import("util.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = matrix.Mat;
const ConstView = matrix.ConstView;
const isComplex = vector.isComplex;
const Real = vector.Real;
const ElementOf = util.ElementOf;
const constOf = util.constOf;
const copyOf = util.copyOf;

pub const HessError = error{NotSquare} || Allocator.Error || Io.Cancelable;

pub fn Hess(comptime T: type) type {
    return struct {
        /// Unitary.
        p: Mat(T),
        /// Upper Hessenberg: zero below the subdiagonal.
        h: Mat(T),

        pub fn deinit(f: @This(), allocator: Allocator) void {
            f.p.deinit(allocator);
            f.h.deinit(allocator);
        }
    };
}

/// `A = P H Pᴴ`, `H` upper Hessenberg (MATLAB's `[P, H] = hess(A)`).
pub fn hess(allocator: Allocator, io: Io, a: anytype) HessError!Hess(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const ca = constOf(a);
    if (ca.rows != ca.cols) return error.NotSquare;
    const n = ca.rows;
    const h = try copyOf(T, allocator, ca);
    errdefer h.deinit(allocator);
    const p = try Mat(T).init(allocator, n, n);
    errdefer p.deinit(allocator);
    if (n == 0) return .{ .p = p, .h = h };
    const tau = try allocator.alloc(T, n - 1);
    defer allocator.free(tau);
    lapack.gehrd(allocator, io, T, n, 0, n, h.data, n, tau) catch |e| return ops.narrow(HessError, e);
    @memcpy(p.data, h.data);
    lapack.orghr(allocator, io, T, n, 0, n, p.data, n, tau) catch |e| return ops.narrow(HessError, e);
    for (0..n -| 2) |j| for (j + 2..n) |i| {
        h.data[i + j * n] = ops.zero(T);
    };
    return .{ .p = p, .h = h };
}

/// `schur`'s and `schurComplex`'s: `NoConvergence` also for a NaN in `A`.
pub const SchurError = error{ NotSquare, NoConvergence } || Allocator.Error || Io.Cancelable;

pub fn Schur(comptime T: type) type {
    return struct {
        /// Unitary (orthogonal for real `T`).
        u: Mat(T),
        /// Upper triangular, or for real `T` quasi-triangular.
        t: Mat(T),

        pub fn deinit(f: @This(), allocator: Allocator) void {
            f.u.deinit(allocator);
            f.t.deinit(allocator);
        }
    };
}

/// `A = U T Uᴴ` (MATLAB's `[U, T] = schur(A)`): for a real `A` the real
/// Schur form, a 2 x 2 block on the diagonal per complex pair.
pub fn schur(allocator: Allocator, io: Io, a: anytype) SchurError!Schur(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const ca = constOf(a);
    if (ca.rows != ca.cols) return error.NotSquare;
    return gees(T, allocator, io, try copyOf(T, allocator, ca));
}

/// `A = U T Uᴴ` with `T` triangular even for a real `A` (MATLAB's
/// `schur(A, 'complex')`).
pub fn schurComplex(allocator: Allocator, io: Io, a: anytype) SchurError!Schur(std.math.Complex(Real(ElementOf(@TypeOf(a))))) {
    const T = ElementOf(@TypeOf(a));
    const C = std.math.Complex(Real(T));
    const ca = constOf(a);
    if (ca.rows != ca.cols) return error.NotSquare;
    const n = ca.rows;
    const t = try Mat(C).init(allocator, n, n);
    for (0..n) |j| for (t.data[j * n ..][0..n], ca.col(j)) |*d, x| {
        d.* = if (comptime isComplex(T)) x else C.init(x, 0);
    };
    return gees(C, allocator, io, t);
}

/// `gees` on `t` (taken), with `U` from its Schur vectors.
fn gees(comptime T: type, allocator: Allocator, io: Io, t: Mat(T)) SchurError!Schur(T) {
    errdefer t.deinit(allocator);
    const n = t.rows;
    const u = try Mat(T).init(allocator, n, n);
    errdefer u.deinit(allocator);
    if (n == 0) return .{ .u = u, .t = t };
    const w = try allocator.alloc(T, if (comptime isComplex(T)) n else 2 * n);
    defer allocator.free(w);
    const ev: lapack.Eigenvalues(T) = if (comptime isComplex(T)) w else .{ .re = w[0..n], .im = w[n..] };
    _ = lapack.gees(allocator, io, T, .compute, null, n, t.data, n, ev, u.data, n, null) catch |e| return switch (e) {
        error.NanInA => error.NoConvergence,
        else => ops.narrow(SchurError, e),
    };
    return .{ .u = u, .t = t };
}

// ---- tests ----

const testing = std.testing;

/// `‖Q X Qᴴ - A‖` over the entries, real `T`.
fn residual(a: []const f64, q: Mat(f64), x: Mat(f64)) f64 {
    const n = q.rows;
    var worst: f64 = 0;
    for (0..n) |i| for (0..n) |j| {
        var s: f64 = 0;
        for (0..n) |k| for (0..n) |l| {
            s += q.at(i, k) * x.at(k, l) * q.at(j, l);
        };
        worst = @max(worst, @abs(s - a[i + j * n]));
    };
    return worst;
}

test "hess and schur: P H Pᵀ and U T Uᵀ give A back, forms as promised" {
    const al = testing.allocator;
    const io = testing.io;
    // magic(4)
    const a = [_]f64{ 16, 5, 9, 4, 2, 11, 7, 14, 3, 10, 6, 15, 13, 8, 12, 1 };
    const v: ConstView(f64) = try .init(&a, 4, 4, 4);
    const h = try hess(al, io, v);
    defer h.deinit(al);
    try testing.expect(residual(&a, h.p, h.h) < 1e-12 * 34);
    for (0..2) |j| for (j + 2..4) |i| try testing.expectEqual(0, h.h.at(i, j));

    // a rotation block: the real form keeps a 2 x 2 block
    const r = [_]f64{ 0, -1, 0, 1, 0, 0, 0, 0, 2 };
    const s = try schur(al, io, try ConstView(f64).init(&r, 3, 3, 3));
    defer s.deinit(al);
    try testing.expect(residual(&r, s.u, s.t) < 1e-13);
    const sc = try schurComplex(al, io, try ConstView(f64).init(&r, 3, 3, 3));
    defer sc.deinit(al);
    for (0..3) |j| for (j + 1..3) |i| try testing.expectEqual(0, sc.t.at(i, j).re);
    var trace: f64 = 0;
    for (0..3) |i| trace += sc.t.at(i, i).re;
    try testing.expectApproxEqAbs(2, trace, 1e-13);
}
