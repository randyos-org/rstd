//! EIG (`syevd`/`heevd` or `geev`) over LAPACK (M2), in place of EISPACK's
//! HTRIDI, IMTQL2, HTRIBK, CORTH and COMQR3 (`MATFN2`).

const std = @import("std");
const zigstd = @import("../../zigstd.zig");
const scalar = @import("scalar.zig");
const CxMat = @import("CxMat.zig");
const lapack = @import("../lapack/__root__.zig");
const ops = @import("../lapack/ops.zig");

const Cx = scalar.Cx;
const Allocator = std.mem.Allocator;
const Io = std.Io;

/// Eigenvalues, and optionally eigenvectors, of an `n`x`n` complex matrix.
/// All slices are owned by the caller and must be freed with the same
/// allocator used to create them (via `deinit`).
pub const EigResult = struct {
    values_r: []f64,
    values_i: []f64,
    /// Column-major `n`x`n`, present iff `eig` was asked to compute
    /// vectors: orthonormal for a Hermitian `A`, else each of unit 2-norm
    /// with its largest entry real (`geev`'s).
    vectors_r: ?[]f64,
    vectors_i: ?[]f64,

    pub fn deinit(self: EigResult, allocator: Allocator) void {
        allocator.free(self.values_r);
        allocator.free(self.values_i);
        if (self.vectors_r) |v| allocator.free(v);
        if (self.vectors_i) |v| allocator.free(v);
    }
};

/// `NoConvergence` for a QR iteration that did not converge, or an
/// infinity or NaN in the matrix.
pub const EigError = error{NoConvergence} || Allocator.Error || Io.Cancelable;

fn isHermitian(a: CxMat) bool {
    var j: usize = 1;
    while (j <= a.n) : (j += 1) {
        var i: usize = 1;
        while (i <= j) : (i += 1) {
            if (a.r[a.idx(i, j)] != a.r[a.idx(j, i)] or a.i[a.idx(i, j)] != -a.i[a.idx(j, i)]) return false;
        }
    }
    return true;
}

/// An eigen routine's result, its arguments taken from a square `CxMat`:
/// `InvalidA` (a NaN) is `NoConvergence`.
pub fn converged(r: anytype) EigError!void {
    r catch |e| {
        const x: anyerror = e;
        if (x == error.InvalidA) return error.NoConvergence;
        return ops.narrow(EigError, x);
    };
}

/// EIG: eigenvalues, and with `want_vectors` eigenvectors, of `a` (not
/// mutated): `syevd` (nondecreasing) for Hermitian `a`, else `geev`.
pub fn eig(allocator: Allocator, io: Io, a: CxMat, want_vectors: bool) EigError!EigResult {
    const n = a.n;
    const ni: isize = @intCast(n);
    const ld: isize = @intCast(@max(1, n));
    const values_r = try allocator.alloc(f64, n);
    errdefer allocator.free(values_r);
    const values_i = try allocator.alloc(f64, n);
    errdefer allocator.free(values_i);
    @memset(values_i, 0);
    const vectors: ?CxMat = if (want_vectors) blk: {
        const vr = try allocator.alloc(f64, n * n);
        errdefer allocator.free(vr);
        break :blk CxMat.initPacked(vr, try allocator.alloc(f64, n * n), n, n);
    } else null;
    errdefer if (vectors) |v| v.deinit(allocator);
    const jobz: lapack.Jobz = if (want_vectors) .v else .n;

    switch (a.isReal()) {
        inline else => |real| {
            const T = if (real) f64 else Cx;
            const x = try a.toPacked(T, allocator);
            defer allocator.free(x);
            if (isHermitian(a)) {
                try converged(lapack.syevd(allocator, io, T, jobz, .u, ni, x.ptr, ld, values_r.ptr));
                if (vectors) |v| v.fromPacked(T, x);
            } else {
                var none: [1]T = undefined;
                const vr = if (want_vectors) try allocator.alloc(T, n * n) else &none;
                defer if (want_vectors) allocator.free(vr);
                const ldvr: isize = if (want_vectors) ld else 1;
                if (real) {
                    try converged(lapack.geev(allocator, io, T, .n, jobz, ni, x.ptr, ld, .{ .re = values_r.ptr, .im = values_i.ptr }, &none, 1, vr.ptr, ldvr, null));
                    if (vectors) |v| unpackPairs(v, values_i, vr);
                } else {
                    const w = try allocator.alloc(Cx, n);
                    defer allocator.free(w);
                    try converged(lapack.geev(allocator, io, T, .n, jobz, ni, x.ptr, ld, w.ptr, &none, 1, vr.ptr, ldvr, null));
                    for (w, values_r, values_i) |e, *re, *im| {
                        re.* = e.re;
                        im.* = e.im;
                    }
                    if (vectors) |v| v.fromPacked(T, vr);
                }
            }
        },
    }
    return .{
        .values_r = values_r,
        .values_i = values_i,
        .vectors_r = if (vectors) |v| v.r else null,
        .vectors_i = if (vectors) |v| v.i else null,
    };
}

/// `dgeev`'s real `vr` into `v`: a pair `j`, `j + 1` (`wi[j] > 0`) is
/// `vr(:, j) ± i vr(:, j + 1)`.
fn unpackPairs(v: CxMat, wi: []const f64, vr: []const f64) void {
    const n = v.n;
    var j: usize = 0;
    while (j < n) : (j += 1) {
        const re = vr[j * n ..][0..n];
        const c = v.col(j + 1);
        if (wi[j] == 0) {
            @memcpy(c.r, re);
            @memset(c.i, 0);
            continue;
        }
        const im = vr[(j + 1) * n ..][0..n];
        const d = v.col(j + 2);
        for (c.r, c.i, d.r, d.i, re, im) |*cr, *ci, *dr, *di, x, y| {
            cr.* = x;
            ci.* = y;
            dr.* = x;
            di.* = -y;
        }
        j += 1;
    }
}

// ---- tests ------------------------------------------------------------------

const testing = std.testing;

test "eig dispatches to the Hermitian path for a symmetric matrix" {
    // A = [[2,1],[1,2]]: eigenvalues 1 and 3.
    var ar = [_]f64{ 2.0, 1.0, 1.0, 2.0 };
    var ai: [4]f64 = @splat(0.0);
    const result = try eig(testing.allocator, testing.io, CxMat.initPacked(&ar, &ai, 2, 2), true);
    defer result.deinit(testing.allocator);

    try zigstd.testExpectApproxEqAbs(1.0, result.values_r[0], 1e-15);
    try zigstd.testExpectApproxEqAbs(3.0, result.values_r[1], 1e-15);
    try zigstd.testExpectTrue(result.vectors_r != null);
}

test "eig dispatches to the general path for a non-symmetric matrix" {
    var ar = [_]f64{ 1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 10.0 };
    var ai: [9]f64 = @splat(0.0);
    const result = try eig(testing.allocator, testing.io, CxMat.initPacked(&ar, &ai, 3, 3), false);
    defer result.deinit(testing.allocator);

    try zigstd.testExpectTrue(result.vectors_r == null);
    var trace: f64 = 0.0;
    for (result.values_r) |v| trace += v;
    try zigstd.testExpectApproxEqAbs(16.0, trace, 1e-13);
}

test "eig's vectors satisfy A v = lambda v: real pairs, complex general, Hermitian" {
    // a rotation-like real matrix (a conjugate pair and 2), a complex
    // general one, and a complex Hermitian one
    const cases = [_]struct { r: [9]f64, i: [9]f64 }{
        .{ .r = .{ 0, -1, 0, 1, 0, 0, 0, 0, 2 }, .i = @splat(0) },
        .{ .r = .{ 1, 0, 1, 2, 3, 0, 0, 1, 1 }, .i = .{ 1, 0, 0, 0, -1, 0, 2, 0, 0 } },
        .{ .r = .{ 2, 1, 0, 1, 3, 1, 0, 1, 4 }, .i = .{ 0, -1, 0, 1, 0, 2, 0, -2, 0 } },
    };
    for (cases) |c| {
        var ar = c.r;
        var ai = c.i;
        const a = CxMat.initPacked(&ar, &ai, 3, 3);
        const result = try eig(testing.allocator, testing.io, a, true);
        defer result.deinit(testing.allocator);
        const v = CxMat.initPacked(result.vectors_r.?, result.vectors_i.?, 3, 3);
        for (1..4) |j| {
            const lambda = Cx.init(result.values_r[j - 1], result.values_i[j - 1]);
            for (1..4) |i| {
                var s = Cx.init(0, 0);
                for (1..4) |k| s = s.add(a.get(i, k).mul(v.get(k, j)));
                const d = s.sub(lambda.mul(v.get(i, j)));
                try zigstd.testExpectTrue(d.magnitude() < 1e-14);
            }
        }
    }
    // the real case's pair is exactly conjugate
    var ar = cases[0].r;
    var ai = cases[0].i;
    const result = try eig(testing.allocator, testing.io, CxMat.initPacked(&ar, &ai, 3, 3), false);
    defer result.deinit(testing.allocator);
    const j = std.mem.indexOfScalar(f64, result.values_i, 0).?;
    const p: usize = if (j == 0) 1 else 0;
    try zigstd.testExpectApproxEqAbs(2.0, result.values_r[j], 1e-15);
    try zigstd.testExpectEqual(result.values_r[p], result.values_r[p + 1]);
    try zigstd.testExpectEqual(result.values_i[p], -result.values_i[p + 1]);
}
