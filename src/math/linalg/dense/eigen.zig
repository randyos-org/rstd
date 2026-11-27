//! Eigenvalues and eigenvectors.

const sysinfo = @import("builtin");
const std = @import("std");
const Allocator = std.mem.Allocator;
const Io = std.Io;
const blas = @import("../../kernels/blas/__root__.zig");
const ext = @import("../../kernels/blas/ext/interface.zig");
const l3 = @import("../../kernels/blas/blas3/core.zig");
const l2 = @import("../../kernels/blas/blas2/core.zig");
const l1 = @import("../../kernels/blas/blas1/interface.zig");
const lapack = @import("../../kernels/lapack/__root__.zig");
const vector = @import("../vector.zig");
const matrix = @import("../matrix.zig");
const ops = @import("../../kernels/lapack/ops.zig");

const View = matrix.View;
const ConstView = matrix.ConstView;
const Mat = matrix.Mat;
const isComplex = vector.isComplex;
const Real = vector.Real;
const eMul = matrix.eMul;
const eMulAdd = matrix.eMulAdd;

const singular = @import("singular.zig");
const svdError = singular.svdError;
const util = @import("util.zig");
const ElementOf = util.ElementOf;
const constOf = util.constOf;
const copyOf = util.copyOf;

/// `eig`'s result: Hermitian `A` (decided by its values alone, exactly)
/// has real eigenvalues, nondecreasing, and orthonormal eigenvectors;
/// otherwise they are complex, a real `A`'s conjugate pairs together with
/// the positive imaginary part first, and each eigenvector (`A v = λ v`) of
/// unit 2-norm with its largest entry real (design.md, decision 6).
pub fn Eig(comptime T: type) type {
    const C = std.math.Complex(Real(T));
    return union(enum) {
        hermitian: struct { values: []Real(T), vectors: ?Mat(T) },
        general: struct { values: []C, vectors: ?Mat(C) },

        pub fn deinit(e: @This(), allocator: Allocator) void {
            switch (e) {
                inline else => |r| {
                    allocator.free(r.values);
                    if (r.vectors) |v| v.deinit(allocator);
                },
            }
        }
    };
}

pub const EigError = error{ NotSquare, NoConvergence } || Allocator.Error || Io.Cancelable;

/// Eigenvalues, and with `want_vectors` eigenvectors (MATLAB's `eig`):
/// `syevd` (`heevd`) for Hermitian `A`, else `geev`. A NaN in `A` gives
/// `error.NoConvergence`.
pub fn eig(allocator: Allocator, io: Io, a: anytype, want_vectors: bool) EigError!Eig(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const C = std.math.Complex(Real(T));
    const ca = constOf(a);
    if (ca.rows != ca.cols) return error.NotSquare;
    const n = ca.rows;
    const ld = @max(n, 1);
    const jobz: lapack.VectorJob = if (want_vectors) .compute else .none;
    const f = try copyOf(T, allocator, ca);
    if (isHermitian(ca)) {
        errdefer f.deinit(allocator);
        const w = try allocator.alloc(Real(T), n);
        errdefer allocator.free(w);
        lapack.syevd(allocator, io, T, jobz, .u, n, f.data, ld, w) catch |e| return svdError(e);
        if (!want_vectors) f.deinit(allocator);
        return .{ .hermitian = .{ .values = w, .vectors = if (want_vectors) f else null } };
    }
    defer f.deinit(allocator);
    const w = try allocator.alloc(C, n);
    errdefer allocator.free(w);
    const vr: ?Mat(T) = if (want_vectors) try Mat(T).init(allocator, n, n) else null;
    var none: [1]T = undefined;
    const pr: []T = if (vr) |x| x.data else &none;
    if (comptime isComplex(T)) {
        errdefer if (vr) |x| x.deinit(allocator);
        lapack.geev(allocator, io, T, .none, jobz, n, f.data, ld, w, &none, 1, pr, ld, null) catch |e| return svdError(e);
        return .{ .general = .{ .values = w, .vectors = vr } };
    }
    defer if (vr) |x| x.deinit(allocator);
    const wr = try allocator.alloc(T, 2 * n);
    defer allocator.free(wr);
    lapack.geev(allocator, io, T, .none, jobz, n, f.data, ld, .{ .re = wr[0..n], .im = wr[n..] }, &none, 1, pr, ld, null) catch |e| return svdError(e);
    for (w, wr[0..n], wr[n..]) |*d, re, im| d.* = .init(re, im);
    const x = vr orelse return .{ .general = .{ .values = w, .vectors = null } };
    // `v(:, j) ± i v(:, j + 1)` for a pair, `wi(j) > 0` first
    const v = try Mat(C).init(allocator, n, n);
    var j: usize = 0;
    while (j < n) : (j += 1) {
        const c = x.data[j * n ..][0..n];
        if (wr[n + j] == 0) {
            for (v.data[j * n ..][0..n], c) |*d, re| d.* = .init(re, 0);
            continue;
        }
        const ci = x.data[(j + 1) * n ..][0..n];
        for (v.data[j * n ..][0..n], v.data[(j + 1) * n ..][0..n], c, ci) |*p, *q, re, im| {
            p.* = .init(re, im);
            q.* = .init(re, -im);
        }
        j += 1;
    }
    return .{ .general = .{ .values = w, .vectors = v } };
}

pub const EigGeneralizedError = error{ NotSquare, DimensionMismatch, NoConvergence } || Allocator.Error || Io.Cancelable;

/// The generalized eigenvalues of `A x = λ B x` (MATLAB's `eig(A, B)`), with
/// `want_vectors` the right eigenvectors. Hermitian `A` and Hermitian
/// positive definite `B`: `sygvd` (`hegvd`), `Eig`'s `.hermitian`, real
/// values nondecreasing and vectors of unit `B`-norm (`Vᴴ B V = I`),
/// MATLAB's `chol` algorithm. Any other pencil, a Hermitian `B` that turns
/// out indefinite too: `ggev` (MATLAB's `qz`), `.general`, `λ = alpha /
/// beta` (Inf where `beta` is 0, NaN for `0 / 0`), vectors as `eig`'s. A NaN
/// gives `error.NoConvergence`.
pub fn eigGeneralized(allocator: Allocator, io: Io, a: anytype, b: anytype, want_vectors: bool) EigGeneralizedError!Eig(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const ca = constOf(a);
    const cb = constOf(b);
    if (ca.rows != ca.cols or cb.rows != cb.cols) return error.NotSquare;
    if (cb.rows != ca.rows) return error.DimensionMismatch;
    const n = ca.rows;
    const ld = @max(n, 1);
    if (isHermitian(ca) and isHermitian(cb)) chol: {
        const f = try copyOf(T, allocator, ca);
        errdefer f.deinit(allocator);
        const g = try copyOf(T, allocator, cb);
        defer g.deinit(allocator);
        const w = try allocator.alloc(Real(T), n);
        errdefer allocator.free(w);
        const jobz: lapack.VectorJob = if (want_vectors) .compute else .none;
        const solve = if (comptime isComplex(T)) lapack.hegvd else lapack.sygvd;
        solve(allocator, io, T, .ax_equals_lambda_bx, jobz, .u, n, f.data, ld, g.data, ld, w) catch |e| switch (e) {
            error.NotPositiveDefinite => {
                f.deinit(allocator);
                allocator.free(w);
                break :chol;
            },
            else => return ops.narrow(EigGeneralizedError, e),
        };
        if (!want_vectors) f.deinit(allocator);
        return .{ .hermitian = .{ .values = w, .vectors = if (want_vectors) f else null } };
    }
    return generalPencil(T, allocator, io, ca, cb, want_vectors);
}

/// `eigGeneralized` by `ggev`.
fn generalPencil(comptime T: type, allocator: Allocator, io: Io, ca: ConstView(T), cb: ConstView(T), want_vectors: bool) EigGeneralizedError!Eig(T) {
    const C = std.math.Complex(Real(T));
    const n = ca.rows;
    const ld = @max(n, 1);
    const f = try copyOf(T, allocator, ca);
    defer f.deinit(allocator);
    const g = try copyOf(T, allocator, cb);
    defer g.deinit(allocator);
    const vr = try Mat(T).init(allocator, n, n);
    defer vr.deinit(allocator);
    const jobvr: lapack.VectorJob = if (want_vectors) .compute else .none;
    var none: [1]T = undefined;
    const values = try allocator.alloc(C, n);
    errdefer allocator.free(values);
    const w = try allocator.alloc(T, 3 * n);
    defer allocator.free(w);
    if (comptime isComplex(T)) {
        lapack.ggevComplex(allocator, io, T, .none, jobvr, n, f.data, ld, g.data, ld, w[0..n], w[n .. 2 * n], &none, 1, vr.data, ld) catch |e| return ops.narrow(EigGeneralizedError, e);
        for (values, w[0..n], w[n .. 2 * n]) |*v, al, be| v.* = ratio(Real(T), al, be.re);
    } else {
        lapack.ggevReal(allocator, io, T, .none, jobvr, n, f.data, ld, g.data, ld, w[0..n], w[n .. 2 * n], w[2 * n ..], &none, 1, vr.data, ld) catch |e| return ops.narrow(EigGeneralizedError, e);
        for (values, w[0..n], w[n .. 2 * n], w[2 * n ..]) |*v, re, im, be| v.* = ratio(Real(T), C.init(re, im), be);
    }
    if (!want_vectors) return .{ .general = .{ .values = values, .vectors = null } };
    const v = try Mat(C).init(allocator, n, n);
    errdefer v.deinit(allocator);
    var j: usize = 0;
    while (j < n) : (j += 1) {
        const c = vr.data[j * n ..][0..n];
        if (comptime isComplex(T)) {
            @memcpy(v.data[j * n ..][0..n], c);
        } else if (w[n + j] == 0) {
            for (v.data[j * n ..][0..n], c) |*d, re| d.* = .init(re, 0);
        } else {
            // a pair: `v(:, j) ± i v(:, j + 1)`
            const ci = vr.data[(j + 1) * n ..][0..n];
            for (v.data[j * n ..][0..n], v.data[(j + 1) * n ..][0..n], c, ci) |*p, *q, re, im| {
                p.* = .init(re, im);
                q.* = .init(re, -im);
            }
            j += 1;
        }
    }
    for (0..n) |k| largestReal(C, v.data[k * n ..][0..n]);
    return .{ .general = .{ .values = values, .vectors = v } };
}

/// `alpha / beta`, `beta` real: Inf where `beta` is 0, NaN for `0 / 0`.
fn ratio(comptime F: type, alpha: std.math.Complex(F), beta: F) std.math.Complex(F) {
    if (beta != 0) return .init(alpha.re / beta, alpha.im / beta);
    if (alpha.re == 0 and alpha.im == 0) return .init(std.math.nan(F), 0);
    return .init(std.math.inf(F), 0);
}

/// A unit vector turned so its largest entry is real and positive, as
/// `eig`'s `geev` vectors (design.md, decision 6).
fn largestReal(comptime C: type, x: []C) void {
    var best: usize = 0;
    var m: @TypeOf(x[0].re) = -1;
    for (x, 0..) |e, i| {
        const a = std.math.hypot(e.re, e.im);
        if (a > m) {
            m = a;
            best = i;
        }
    }
    if (m <= 0) return;
    const phase = C.init(x[best].re / m, -x[best].im / m);
    for (x) |*e| e.* = e.mul(phase);
    x[best].im = 0;
}

/// `a == aᴴ` exactly, the diagonal real.
pub fn isHermitian(a: anytype) bool {
    const T = @TypeOf(a).Element;
    for (0..a.rows) |j| {
        if (comptime isComplex(T)) if (a.at(j, j).im != 0) return false;
        for (j + 1..a.rows) |i| {
            const x = a.at(i, j);
            const y = ops.conj(T, a.at(j, i));
            const same = if (comptime isComplex(T)) x.re == y.re and x.im == y.im else x == y;
            if (!same) return false;
        }
    }
    return true;
}
