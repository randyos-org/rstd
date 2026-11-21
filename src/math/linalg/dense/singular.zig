//! Matrix norms, the SVD and what it gives: `pinv`, `rank`, `cond`.

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

const util = @import("util.zig");
const ElementOf = util.ElementOf;
const Error = util.Error;
const constOf = util.constOf;
const copyOf = util.copyOf;
const narrow = util.narrow;
const toIsize = util.toIsize;

/// `matrixNorm`'s norms: the largest column sum of `|a(i, j)|`, the largest
/// row sum, Frobenius, and the largest singular value.
pub const MatrixNorm = enum { one, inf, fro, two };

pub const MatrixNormError = error{NoConvergence} || Allocator.Error || Io.Cancelable;

/// A matrix norm (MATLAB's `norm(A, p)`); 0 when empty. `lange`, or `.two`
/// the first singular value from `gesvd`.
pub fn matrixNorm(allocator: Allocator, io: Io, a: anytype, which: MatrixNorm) MatrixNormError!Real(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const ca = constOf(a);
    const m = ca.rows;
    const n = ca.cols;
    if (m == 0 or n == 0) return 0;
    const kind: lapack.NormType = switch (which) {
        .one => .one,
        .inf => .infinity,
        .fro => .frobenius,
        .two => {
            const s = try singularValues(allocator, io, ca);
            defer allocator.free(s);
            return s[0];
        },
    };
    return lapack.lange(allocator, io, T, kind, m, n, ca.data, ca.ld) catch |e| return narrow(Allocator.Error, e);
}

/// `svd`'s outputs for the `m` by `n` `A`, `k = min(m, n)`: `.values` `s`
/// alone; `.economy` `U` `m` by `k` and `V` `n` by `k`; `.full` `U` `m` by
/// `m` and `V` `n` by `n`.
pub const SvdMode = enum { values, economy, full };

/// `svd`'s result, `A = U Σ Vᴴ`.
pub fn Svd(comptime T: type) type {
    return struct {
        /// The `k` singular values, nonincreasing.
        s: []Real(T),
        /// Orthonormal columns; null for `.values`.
        u: ?Mat(T),
        v: ?Mat(T),

        pub fn deinit(f: @This(), allocator: Allocator) void {
            allocator.free(f.s);
            if (f.u) |u| u.deinit(allocator);
            if (f.v) |v| v.deinit(allocator);
        }
    };
}

/// `svd`'s, `pinv`'s, `rank`'s and `cond`'s: `NoConvergence` also for a NaN
/// in `A`.
pub const SvdError = error{NoConvergence} || Allocator.Error || Io.Cancelable;

/// `NanInA` is `NoConvergence`.
pub fn svdError(e: anyerror) SvdError {
    return if (e == error.NanInA) error.NoConvergence else narrow(SvdError, e);
}

/// `A = U Σ Vᴴ` (MATLAB's `svd`): the values alone by `gesvd` (QR
/// iteration, which keeps a graded matrix's small values), with vectors by
/// `gesdd` (divide and conquer), so the two may differ in the last bits.
/// Empty `A` with `.full`: `U` and `V` are identities.
pub fn svd(allocator: Allocator, io: Io, a: anytype, mode: SvdMode) SvdError!Svd(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const ca = constOf(a);
    if (mode == .values) return .{ .s = try singularValues(allocator, io, ca), .u = null, .v = null };
    const m = ca.rows;
    const n = ca.cols;
    const k = @min(m, n);
    const uc = if (mode == .full) m else k;
    const vc = if (mode == .full) n else k;
    const s = try allocator.alloc(Real(T), k);
    errdefer allocator.free(s);
    const u = try Mat(T).zeros(allocator, m, uc);
    errdefer u.deinit(allocator);
    const v = try Mat(T).zeros(allocator, n, vc);
    errdefer v.deinit(allocator);
    if (k == 0) {
        for (0..uc) |i| u.data[i * (m + 1)] = ops.one(T);
        for (0..vc) |i| v.data[i * (n + 1)] = ops.one(T);
        return .{ .s = s, .u = u, .v = v };
    }
    const f = try copyOf(T, allocator, ca);
    defer f.deinit(allocator);
    const vt = try allocator.alloc(T, vc * n);
    defer allocator.free(vt);
    lapack.gesdd(allocator, io, T, if (mode == .full) .all else .leading, m, n, f.data, m, s, u.data, m, vt, vc) catch |e| return svdError(e);
    for (0..n) |j| for (0..vc) |i| {
        v.data[j + i * n] = ops.conj(T, vt[i + j * vc]);
    };
    return .{ .s = s, .u = u, .v = v };
}

/// `a`'s `min(m, n)` singular values by `gesvd`, owned by the caller.
pub fn singularValues(allocator: Allocator, io: Io, a: anytype) SvdError![]Real(@TypeOf(a).Element) {
    const T = @TypeOf(a).Element;
    const m = a.rows;
    const n = a.cols;
    const s = try allocator.alloc(Real(T), @min(m, n));
    errdefer allocator.free(s);
    if (s.len == 0) return s;
    const f = try copyOf(T, allocator, a);
    defer f.deinit(allocator);
    lapack.gesvd(allocator, io, T, .none, .none, m, n, f.data, m, s, f.data[0..0], 1, f.data[0..0], 1, null) catch |e| return svdError(e);
    return s;
}

/// The spacing of floats at `x` (MATLAB's `eps(x)`).
pub fn epsOf(x: anytype) @TypeOf(x) {
    const F = @TypeOf(x);
    if (x == 0) return std.math.floatTrueMin(F);
    const e = std.math.frexp(x).exponent;
    return @max(std.math.ldexp(std.math.floatEps(F), e - 1), std.math.floatTrueMin(F));
}

/// How many of the nonincreasing `s` exceed `tol`, by default
/// `max(m, n) eps(σ₁)` (MATLAB's `rank` and `pinv`).
pub fn rankOf(comptime F: type, s: []const F, m: usize, n: usize, tol: ?F) usize {
    if (s.len == 0) return 0;
    const t = tol orelse (@as(F, @floatFromInt(@max(m, n))) * epsOf(s[0]));
    var r: usize = 0;
    while (r < s.len and s[r] > t) r += 1;
    return r;
}

/// The Moore-Penrose pseudoinverse (MATLAB's `pinv`), `n` by `m`, owned by
/// the caller: `V Σ⁺ Uᴴ` over the singular values above `tol`, by default
/// `max(m, n) eps(σ₁)`.
pub fn pinv(allocator: Allocator, io: Io, a: anytype, tol: ?Real(ElementOf(@TypeOf(a)))) SvdError!Mat(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const ca = constOf(a);
    const m = ca.rows;
    const n = ca.cols;
    const x = try Mat(T).zeros(allocator, n, m);
    errdefer x.deinit(allocator);
    const f = try svd(allocator, io, ca, .economy);
    defer f.deinit(allocator);
    const r = rankOf(Real(T), f.s, m, n, tol);
    if (r == 0) return x;
    const u = f.u.?;
    const v = f.v.?;
    for (f.s[0..r], 0..) |sj, j| for (v.data[j * n ..][0..n]) |*e| {
        e.* = if (comptime isComplex(T)) T.init(e.re / sj, e.im / sj) else e.* / sj;
    };
    const t: blas.Trans = if (comptime isComplex(T)) .c else .t;
    l3.gemm(allocator, io, T, .n, t, toIsize(n), toIsize(m), toIsize(r), ops.one(T), v.data.ptr, toIsize(n), u.data.ptr, toIsize(m), ops.zero(T), x.data.ptr, toIsize(n)) catch |e| return narrow(Error, e);
    return x;
}

/// The number of singular values above `tol`, by default `max(m, n)
/// eps(σ₁)` (MATLAB's `rank`); 0 when empty.
pub fn rank(allocator: Allocator, io: Io, a: anytype, tol: ?Real(ElementOf(@TypeOf(a)))) SvdError!usize {
    const ca = constOf(a);
    const s = try singularValues(allocator, io, ca);
    defer allocator.free(s);
    return rankOf(Real(ElementOf(@TypeOf(a))), s, ca.rows, ca.cols, tol);
}

/// The 2-norm condition number `σ₁ / σₖ` (MATLAB's `cond`); infinite when
/// `σₖ` is 0, 0 when empty.
pub fn cond(allocator: Allocator, io: Io, a: anytype) SvdError!Real(ElementOf(@TypeOf(a))) {
    const s = try singularValues(allocator, io, constOf(a));
    defer allocator.free(s);
    if (s.len == 0) return 0;
    const last = s[s.len - 1];
    return if (last == 0) std.math.inf(@TypeOf(last)) else s[0] / last;
}
