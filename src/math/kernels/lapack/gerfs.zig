//! `xGERFS`: iterative refinement of a `getrs` solution with error bounds, Reference-LAPACK
//! 3.12.0 (`lapack-netlib/SRC/dgerfs.f`, `zgerfs.f`, OpenBLAS v0.3.34-317; modified BSD).

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const interface = @import("interface.zig");
const ops = @import("ops.zig");
const blas_calls = @import("blas_calls.zig");
const lacn2_estimate = @import("lacn2.zig");
const lu = @import("lu.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const isComplex = types.isComplex;
const Trans = types.Trans;
const Error = ops.Error;
const abs1 = ops.abs1;
const lacn2 = lacn2_estimate.lacn2;

/// `ITMAX` (`dgerfs.f`): the most refinement steps per right-hand side.
const itmax = 5;

/// What `gerfs` returns: a `trans` other than `n`, `t` or `c`, a pivot that
/// is not a row of `A`, a slice too short or a leading dimension too small,
/// and the allocator's and cancellation.
pub const GerfsError = interface.TransError || interface.PivotError || interface.ShapeError || Allocator.Error || Io.Cancelable;

/// `?GERFS` (`dgerfs.f`, `zgerfs.f`): improves each column of `X`
/// (`n x nrhs`), a solution of `op(A) X = B`, by iterative refinement with
/// `af` and `ipiv`, `getrf`'s factor of `A`, and bounds its error.
/// `trans` is `n`, `t` or `c` (real `c` is `t`). For each column `j`,
/// `berr[j]` is the componentwise backward error, the smallest relative
/// change to the entries of `A` and `B` that makes `X(:, j)` exact, and
/// `ferr[j]` an estimated bound on `‖X(:, j) - Xtrue‖∞ / ‖X(:, j)‖∞`, by
/// `lacn2`. A column is refined while its backward error exceeds `eps`,
/// halves at least at each step and has had at most `itmax` steps, the
/// reference's stopping test.
pub fn gerfs(allocator: Allocator, io: Io, comptime T: type, trans: Trans, n: usize, nrhs: usize, a: []const T, lda: usize, af: []const T, ldaf: usize, ipiv: []const usize, b: []const T, ldb: usize, x: []T, ldx: usize, ferr: []Float(T), berr: []Float(T)) GerfsError!void {
    comptime ops.requireFloat(T, "gerfs");
    const F = Float(T);
    const notran = switch (trans) {
        .n => true,
        .t, .c => false,
        .r, .o, .u, .s, .d => return error.UnsupportedTranspose,
    };
    try util.holds(a.len, n, n, lda);
    try util.holds(af.len, n, n, ldaf);
    try util.holds(b.len, n, nrhs, ldb);
    try util.holds(x.len, n, nrhs, ldx);
    if (ferr.len < nrhs or berr.len < nrhs) return error.SliceTooShort;
    try interface.pivotsIn(ipiv, n);
    if (n == 0 or nrhs == 0) {
        @memset(ferr[0..nrhs], 0);
        @memset(berr[0..nrhs], 0);
        return;
    }
    // `zgerfs.f`'s `TRANSN` and `TRANST`: `op(A)` and its adjoint for
    // `lacn2`; `dgerfs.f` has `TRANS` and `T` or `N`
    const transn: Trans = if (notran) .n else if (comptime isComplex(T)) .c else .t;
    const transt: Trans = if (!notran) .n else if (comptime isComplex(T)) .c else .t;
    const nz: F = @floatFromInt(n + 1);
    const eps = std.math.floatEps(F);
    const safmin = std.math.floatMin(F);
    const safe1 = nz * safmin;
    const safe2 = safe1 / eps;
    // the residual and `lacn2`'s `x`, its `v`, and `|op(A)| |X| + |B|`
    const work = try allocator.alloc(T, 2 * n);
    defer allocator.free(work);
    const r = work[0..n];
    const v = work[n..];
    const bound = try allocator.alloc(F, n);
    defer allocator.free(bound);
    const am: ops.Mat(T) = .{ .ptr = @constCast(a.ptr), .rows = n, .cols = n, .ld = lda };
    const ctx: Ferr(T) = .{ .al = allocator, .io = io, .n = n, .af = af, .ldaf = ldaf, .ipiv = ipiv, .bound = bound, .transn = transn, .transt = transt };
    for (0..nrhs) |j| {
        const bj = b[j * ldb ..][0..n];
        const xj = x[j * ldx ..][0..n];
        var count: usize = 1;
        var lstres: F = 3;
        while (true) {
            // `r = B - op(A) X`
            @memcpy(r, bj);
            try gemv(allocator, io, T, trans, am, xj, r);
            for (bound, bj) |*w, e| w.* = abs1(T, e);
            if (notran) {
                for (0..n) |k| {
                    const xk = abs1(T, xj[k]);
                    for (0..n) |i| bound[i] += abs1(T, a[i + k * lda]) * xk;
                }
            } else {
                for (0..n) |k| {
                    var t: F = 0;
                    for (0..n) |i| t += abs1(T, a[i + k * lda]) * abs1(T, xj[i]);
                    bound[k] += t;
                }
            }
            var s: F = 0;
            for (bound, r) |w, e| {
                s = @max(s, if (w > safe2) abs1(T, e) / w else (abs1(T, e) + safe1) / (w + safe1));
            }
            berr[j] = s;
            if (!(s > eps and 2 * s <= lstres and count <= itmax)) break;
            try solve(allocator, io, T, trans, n, af, ldaf, ipiv, r);
            try blas_calls.axpy(io, T, n, ops.one(T), r.ptr, 1, xj.ptr, 1);
            lstres = s;
            count += 1;
        }
        // `|inv(op(A))| (|r| + nz eps (|op(A)| |X| + |B|))`, its infinity
        // norm by `lacn2` over `diag(bound) inv(op(A))ᴴ`
        for (bound, r) |*w, e| {
            w.* = abs1(T, e) + nz * eps * w.* + (if (w.* > safe2) 0 else safe1);
        }
        var c = ctx;
        ferr[j] = lacn2(allocator, T, GerfsError, n, v, r, &c, Ferr(T).apply) catch |e| switch (e) {
            // sized above
            error.SliceTooShort => unreachable,
            else => |err| return err,
        };
        var xmax: F = 0;
        for (xj) |e| xmax = @max(xmax, abs1(T, e));
        if (xmax != 0) ferr[j] /= xmax;
    }
}

/// `r := r - op(A) x` (`?GEMV` with `alpha = -1`, `beta = 1`).
fn gemv(al: Allocator, io: Io, comptime T: type, trans: Trans, a: ops.Mat(T), x: []T, r: []T) Error!void {
    try blas_calls.gemv(al, io, T, trans, ops.scalar(T, -1), a, x.ptr, 1, ops.one(T), r.ptr, 1);
}

/// `x := op(A)⁻¹ x` by `getrs` on one right-hand side.
fn solve(al: Allocator, io: Io, comptime T: type, trans: Trans, n: usize, af: []const T, ldaf: usize, ipiv: []const usize, x: []T) GerfsError!void {
    try lu.getrs(al, io, T, trans, n, 1, af, ldaf, ipiv, x, n);
}

/// `B = diag(bound) inv(op(A))ᴴ` for `lacn2` (`zgerfs.f`'s `KASE` loop):
/// `B x` solves with the adjoint, then scales; `Bᴴ x` scales, then solves.
fn Ferr(comptime T: type) type {
    return struct {
        al: Allocator,
        io: Io,
        n: usize,
        af: []const T,
        ldaf: usize,
        ipiv: []const usize,
        bound: []const Float(T),
        transn: Trans,
        transt: Trans,

        fn apply(s: *@This(), adjoint: bool, x: [*]T) GerfsError!void {
            const v = x[0..s.n];
            if (!adjoint) try solve(s.al, s.io, T, s.transt, s.n, s.af, s.ldaf, s.ipiv, v);
            for (v, s.bound) |*e, w| e.* = ops.mulReal(T, e.*, w);
            if (adjoint) try solve(s.al, s.io, T, s.transn, s.n, s.af, s.ldaf, s.ipiv, v);
        }
    };
}
