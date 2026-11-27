//! The Hessenberg and Schur forms (MATLAB's `hess` and `schur`): `A = P H
//! Pᴴ` by `gehrd` and `orghr`, `A = U T Uᴴ` by `gees`. A real `A` gets the
//! real Schur form, its 2 x 2 blocks standardized, or with `schurComplex`
//! the complex one. Replaces `matlab/schur.zig` (issue #322).

const std = @import("std");
const lapack = @import("../../kernels/lapack/__root__.zig");
const matrix = @import("../matrix.zig");
const vector = @import("../vector.zig");
const ops = @import("../../kernels/lapack/ops.zig");
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

/// `qz`'s and `qzComplex`'s: `NoConvergence` also for a NaN or infinity.
pub const QzError = error{ NotSquare, DimensionMismatch, NoConvergence } || Allocator.Error || Io.Cancelable;

pub fn Qz(comptime T: type) type {
    return struct {
        /// Upper triangular, or for real `T` quasi-triangular.
        aa: Mat(T),
        /// Upper triangular, its diagonal real and nonnegative.
        bb: Mat(T),
        /// Unitary, `Q A Z = AA` and `Q B Z = BB`.
        q: Mat(T),
        z: Mat(T),

        pub fn deinit(f: @This(), allocator: Allocator) void {
            f.aa.deinit(allocator);
            f.bb.deinit(allocator);
            f.q.deinit(allocator);
            f.z.deinit(allocator);
        }
    };
}

/// The generalized Schur form `Q A Z = AA`, `Q B Z = BB` (MATLAB's `[AA, BB,
/// Q, Z] = qz(A, B)`), by `gges`: for real `A` and `B` the real form, a 2 x 2
/// block of `AA` per complex pair.
pub fn qz(allocator: Allocator, io: Io, a: anytype, b: anytype) QzError!Qz(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const ca = constOf(a);
    const cb = constOf(b);
    if (ca.rows != ca.cols or cb.rows != cb.cols) return error.NotSquare;
    if (cb.rows != ca.rows) return error.DimensionMismatch;
    return gges(T, allocator, io, try copyOf(T, allocator, ca), try copyOf(T, allocator, cb));
}

/// `qz` with `AA` triangular even for real `A` and `B` (MATLAB's `qz(A, B,
/// "complex")`).
pub fn qzComplex(allocator: Allocator, io: Io, a: anytype, b: anytype) QzError!Qz(std.math.Complex(Real(ElementOf(@TypeOf(a))))) {
    const T = ElementOf(@TypeOf(a));
    const C = std.math.Complex(Real(T));
    const ca = constOf(a);
    const cb = constOf(b);
    if (ca.rows != ca.cols or cb.rows != cb.cols) return error.NotSquare;
    if (cb.rows != ca.rows) return error.DimensionMismatch;
    const n = ca.rows;
    const f = try Mat(C).init(allocator, n, n);
    errdefer f.deinit(allocator);
    const g = try Mat(C).init(allocator, n, n);
    errdefer g.deinit(allocator);
    for (0..n) |j| for (f.data[j * n ..][0..n], g.data[j * n ..][0..n], ca.col(j), cb.col(j)) |*d, *e, x, y| {
        d.* = if (comptime isComplex(T)) x else C.init(x, 0);
        e.* = if (comptime isComplex(T)) y else C.init(y, 0);
    };
    return gges(C, allocator, io, f, g);
}

/// `gges` on `aa` and `bb` (taken), `Q` returned as LAPACK's `Qᴴ`.
fn gges(comptime T: type, allocator: Allocator, io: Io, aa: Mat(T), bb: Mat(T)) QzError!Qz(T) {
    errdefer aa.deinit(allocator);
    errdefer bb.deinit(allocator);
    const n = aa.rows;
    const ld = @max(n, 1);
    const q = try Mat(T).init(allocator, n, n);
    errdefer q.deinit(allocator);
    const z = try Mat(T).init(allocator, n, n);
    errdefer z.deinit(allocator);
    const w = try allocator.alloc(T, 3 * n);
    defer allocator.free(w);
    if (comptime isComplex(T)) {
        lapack.ggesComplex(allocator, io, T, .compute, .compute, n, aa.data, ld, bb.data, ld, w[0..n], w[n .. 2 * n], q.data, ld, z.data, ld) catch |e| return ops.narrow(QzError, e);
    } else lapack.ggesReal(allocator, io, T, .compute, .compute, n, aa.data, ld, bb.data, ld, w[0..n], w[n .. 2 * n], w[2 * n ..], q.data, ld, z.data, ld) catch |e| return ops.narrow(QzError, e);
    // MATLAB's `Q` is LAPACK's `Qᴴ`
    for (0..n) |j| for (0..j + 1) |i| {
        const x = ops.conj(T, q.data[i + j * n]);
        q.data[i + j * n] = ops.conj(T, q.data[j + i * n]);
        q.data[j + i * n] = x;
    };
    return .{ .aa = aa, .bb = bb, .q = q, .z = z };
}

/// `ordschur`'s: `DimensionMismatch` unless `u`, `t` and `select` agree in
/// size; `EigenvaluesTooClose` when a swap would lose accuracy (real `T`:
/// `trsen`'s `SwapRejected`, MATLAB's error there).
pub const OrdschurError = error{ NotSquare, DimensionMismatch, EigenvaluesTooClose } || Allocator.Error || Io.Cancelable;

/// The Schur form `A = U T Uᴴ` reordered so the eigenvalues `select` picks
/// lead `T` (MATLAB's `[US, TS] = ordschur(U, T, select)`), by `trsen`. For
/// a real `T`, either row of a 2 x 2 block picks its pair.
pub fn ordschur(allocator: Allocator, io: Io, u: anytype, t: anytype, select: []const bool) OrdschurError!Schur(ElementOf(@TypeOf(t))) {
    const T = ElementOf(@TypeOf(t));
    const cu = constOf(u);
    const ct = constOf(t);
    if (ct.rows != ct.cols) return error.NotSquare;
    const n = ct.rows;
    if (cu.rows != n or cu.cols != n or select.len != n) return error.DimensionMismatch;
    const rt = try copyOf(T, allocator, ct);
    errdefer rt.deinit(allocator);
    const ru = try copyOf(T, allocator, cu);
    errdefer ru.deinit(allocator);
    if (n == 0) return .{ .u = ru, .t = rt };
    const w = try allocator.alloc(T, if (comptime isComplex(T)) n else 2 * n);
    defer allocator.free(w);
    const ev: lapack.Eigenvalues(T) = if (comptime isComplex(T)) w else .{ .re = w[0..n], .im = w[n..] };
    _ = lapack.trsen(allocator, io, T, .none, .update, select, n, rt.data, n, ru.data, n, ev) catch |e| return switch (e) {
        error.SwapRejected => error.EigenvaluesTooClose,
        else => ops.narrow(OrdschurError, e),
    };
    return .{ .u = ru, .t = rt };
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

test "ordschur: the selected eigenvalues lead T, U T Uᴴ unchanged" {
    const al = testing.allocator;
    const io = testing.io;
    // triangular with U = I: 3 and 4 move to the front, in their order
    const a = [_]f64{ 1, 0, 0, 0, 0.5, 2, 0, 0, -1, 0.25, 3, 0, 2, 1, -0.5, 4 };
    const v: ConstView(f64) = try .init(&a, 4, 4, 4);
    const eye = [_]f64{ 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 };
    const id: ConstView(f64) = try .init(&eye, 4, 4, 4);
    const s = try ordschur(al, io, id, v, &.{ false, false, true, true });
    defer s.deinit(al);
    try testing.expect(residual(&a, s.u, s.t) < 1e-13 * 8);
    try testing.expectApproxEqAbs(3, s.t.at(0, 0), 1e-13);
    try testing.expectApproxEqAbs(4, s.t.at(1, 1), 1e-13);
    for (0..4) |j| for (j + 1..4) |i| try testing.expectEqual(0, s.t.at(i, j));

    // a real pair ahead of 2: selecting 2 moves it past the 2 x 2 block
    const r = [_]f64{ 0, -1, 0, 1, 0, 0, 0, 0, 2 };
    const rv = try ConstView(f64).init(&r, 3, 3, 3);
    const sr = try schur(al, io, rv);
    defer sr.deinit(al);
    var at2: usize = 0;
    for (0..3) |i| {
        if (@abs(sr.t.at(i, i) - 2) < 1e-12) at2 = i;
    }
    var pick = [_]bool{ false, false, false };
    pick[at2] = true;
    const o = try ordschur(al, io, sr.u, sr.t, &pick);
    defer o.deinit(al);
    try testing.expect(residual(&r, o.u, o.t) < 1e-13 * 8);
    try testing.expectApproxEqAbs(2, o.t.at(0, 0), 1e-13);

    // complex: the same with schurComplex's triangular form
    const sc = try schurComplex(al, io, rv);
    defer sc.deinit(al);
    var cpick = [_]bool{ false, false, false };
    for (0..3) |i| cpick[i] = @abs(sc.t.at(i, i).re - 2) < 1e-12;
    const oc = try ordschur(al, io, sc.u, sc.t, &cpick);
    defer oc.deinit(al);
    try testing.expectApproxEqAbs(2, oc.t.at(0, 0).re, 1e-13);
    try testing.expectApproxEqAbs(0, oc.t.at(0, 0).im, 1e-13);

    try testing.expectError(error.DimensionMismatch, ordschur(al, io, id, v, &.{ true, false }));
}
