//! Eigenvalues and eigenvectors.

const sysinfo = @import("builtin");
const std = @import("std");
const Allocator = std.mem.Allocator;
const Io = std.Io;
const blas = @import("../../math/kernels/blas/__root__.zig");
const ext = @import("../../math/kernels/blas/ext/interface.zig");
const l3 = @import("../../math/kernels/blas/blas3/core.zig");
const l2 = @import("../../math/kernels/blas/blas2/core.zig");
const l1 = @import("../../math/kernels/blas/blas1/interface.zig");
const lapack = @import("../../math/kernels/lapack/__root__.zig");
const vector = @import("../vector.zig");
const matrix = @import("../matrix.zig");
const ops = @import("../../math/kernels/lapack/ops.zig");

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
    const jobz: lapack.Jobz = if (want_vectors) .v else .n;
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
        lapack.geev(allocator, io, T, .n, jobz, n, f.data, ld, w, &none, 1, pr, ld, null) catch |e| return svdError(e);
        return .{ .general = .{ .values = w, .vectors = vr } };
    }
    defer if (vr) |x| x.deinit(allocator);
    const wr = try allocator.alloc(T, 2 * n);
    defer allocator.free(wr);
    lapack.geev(allocator, io, T, .n, jobz, n, f.data, ld, .{ .re = wr[0..n], .im = wr[n..] }, &none, 1, pr, ld, null) catch |e| return svdError(e);
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
