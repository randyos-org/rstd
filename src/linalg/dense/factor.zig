//! LU, Cholesky and QR.

const sysinfo = @import("builtin");
const std = @import("std");
const Allocator = std.mem.Allocator;
const Io = std.Io;
const blas = @import("../blas/__root__.zig");
const ext = @import("../blas/ext/interface.zig");
const l3 = @import("../blas/blas3/core.zig");
const l2 = @import("../blas/blas2/core.zig");
const l1 = @import("../blas/blas1/interface.zig");
const lapack = @import("../lapack/__root__.zig");
const vector = @import("../vector.zig");
const matrix = @import("../matrix.zig");
const ops = @import("../lapack/ops.zig");

const View = matrix.View;
const ConstView = matrix.ConstView;
const Mat = matrix.Mat;
const isComplex = vector.isComplex;
const Real = vector.Real;
const eMul = matrix.eMul;
const eMulAdd = matrix.eMulAdd;

const util = @import("util.zig");
const ElementOf = util.ElementOf;
const Error = util.Error;
const constOf = util.constOf;
const copyOf = util.copyOf;
const narrow = util.narrow;

/// `lu`'s result, `P A = L U` for the `m` by `n` `A`, `getrf`'s factors.
pub fn Lu(comptime T: type) type {
    return struct {
        /// `L` below the diagonal (its unit diagonal implied), `U` on and
        /// above it.
        factors: Mat(T),
        /// Row `i` was interchanged with row `ipiv[i]`, for `i` from 0 up.
        ipiv: []usize,
        /// The first `i` with `U(i, i)` exactly 0; a singular matrix still
        /// factors completely (design.md, decision 5).
        zero_pivot: ?usize,

        const Self = @This();

        pub fn deinit(f: Self, allocator: Allocator) void {
            f.factors.deinit(allocator);
            allocator.free(f.ipiv);
        }

        /// `L`, `m` by `min(m, n)`, unit lower trapezoidal.
        pub fn l(f: Self, allocator: Allocator) Allocator.Error!Mat(T) {
            const k = @min(f.factors.rows, f.factors.cols);
            const r = try Mat(T).zeros(allocator, f.factors.rows, k);
            for (0..k) |j| {
                r.data[j + j * r.rows] = ops.one(T);
                for (j + 1..r.rows) |i| r.data[i + j * r.rows] = f.factors.at(i, j);
            }
            return r;
        }

        /// `U`, `min(m, n)` by `n`, upper trapezoidal.
        pub fn u(f: Self, allocator: Allocator) Allocator.Error!Mat(T) {
            const k = @min(f.factors.rows, f.factors.cols);
            const r = try Mat(T).zeros(allocator, k, f.factors.cols);
            for (0..r.cols) |j| for (0..@min(j + 1, k)) |i| {
                r.data[i + j * k] = f.factors.at(i, j);
            };
            return r;
        }
    };
}

/// `P A = L U` with partial pivoting (MATLAB's `lu`), by `getrf`.
pub fn lu(allocator: Allocator, io: Io, a: anytype) (Allocator.Error || Io.Cancelable)!Lu(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const ca = constOf(a);
    const k = @min(ca.rows, ca.cols);
    const f = try copyOf(T, allocator, ca);
    errdefer f.deinit(allocator);
    const ipiv = try allocator.alloc(usize, k);
    errdefer allocator.free(ipiv);
    const zero_pivot = lapack.getrf(allocator, io, T, ca.rows, ca.cols, f.data, @max(ca.rows, 1), ipiv) catch |e| return narrow(Error, e);
    return .{ .factors = f, .ipiv = ipiv, .zero_pivot = zero_pivot };
}

pub const CholError = error{ NotSquare, NotPositiveDefinite } || Allocator.Error || Io.Cancelable;

/// The upper `R` with `A = Rᴴ R` (MATLAB's `chol`), by `potrf` on `a`'s
/// upper triangle (the lower is not read); zeros below the diagonal.
pub fn chol(allocator: Allocator, io: Io, a: anytype) CholError!Mat(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const ca = constOf(a);
    if (ca.rows != ca.cols) return error.NotSquare;
    const n = ca.rows;
    const r = try copyOf(T, allocator, ca);
    errdefer r.deinit(allocator);
    lapack.potrf(allocator, io, T, .u, n, r.data, @max(n, 1)) catch |e| return narrow(error{NotPositiveDefinite} || Error, e);
    for (0..n) |j| @memset(r.data[j * n + j + 1 .. (j + 1) * n], ops.zero(T));
    return r;
}

/// `qr`'s shapes for the `m` by `n` `A`, `k = min(m, n)`: `.economy` `Q` `m`
/// by `k` and `R` `k` by `n`; `.full` `Q` `m` by `m` and `R` `m` by `n`.
pub const QrMode = enum { economy, full };

/// `qr`'s result.
pub fn Qr(comptime T: type) type {
    return struct {
        /// Orthonormal columns.
        q: Mat(T),
        /// Upper trapezoidal.
        r: Mat(T),

        pub fn deinit(f: @This(), allocator: Allocator) void {
            f.q.deinit(allocator);
            f.r.deinit(allocator);
        }
    };
}

/// `A = Q R` (MATLAB's `qr`, unpivoted), by `geqrf` and `orgqr` (`ungqr`);
/// `R`'s diagonal keeps LAPACK's signs.
pub fn qr(allocator: Allocator, io: Io, a: anytype, mode: QrMode) (Allocator.Error || Io.Cancelable)!Qr(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const ca = constOf(a);
    const m = ca.rows;
    const n = ca.cols;
    const k = @min(m, n);
    const qc = if (mode == .full) m else k;
    const f = try copyOf(T, allocator, ca);
    defer f.deinit(allocator);
    const tau = try allocator.alloc(T, k);
    defer allocator.free(tau);
    const ld = @max(m, 1);
    lapack.geqrf(allocator, io, T, m, n, f.data, ld, tau) catch |e| return narrow(Error, e);
    const r = try Mat(T).zeros(allocator, qc, n);
    errdefer r.deinit(allocator);
    for (0..n) |j| for (0..@min(j + 1, k)) |i| {
        r.data[i + j * qc] = f.at(i, j);
    };
    const q = try Mat(T).zeros(allocator, m, qc);
    errdefer q.deinit(allocator);
    for (0..k) |j| @memcpy(q.data[j * m ..][0..m], f.data[j * m ..][0..m]);
    if (m > 0) lapack.orgqr(allocator, io, T, m, qc, k, q.data, ld, tau) catch |e| return narrow(Error, e);
    return .{ .q = q, .r = r };
}
