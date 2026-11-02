//! SVD-based conveniences (`COND`, `RANK`, `PINV`) and `ORTH`, over rstd's
//! LAPACK `gesvd`, `geqrf` and `orgqr` (M2) in place of WSVDC, WQRDC and
//! WQRSL; a real matrix goes to the real routines, so `pinv` and `orth`
//! stay real for it.

const std = @import("std");
const zigstd = @import("../../zigstd.zig");
const scalar = @import("scalar.zig");
const CxMat = @import("CxMat.zig");
const lapack = @import("../lapack/__root__.zig");

const Cx = scalar.Cx;
const Complex = CxMat.Complex;
const Allocator = std.mem.Allocator;
const Io = std.Io;
const ix = scalar.ix;

/// `NoConvergence` from `gesvd` (an infinity or NaN in the matrix).
pub const Error = error{NoConvergence} || Allocator.Error || Io.Cancelable;

/// `gesvd` of `x` (`m x n`, destroyed): `s` (`min(m, n)`), and with `u`
/// and `vt` non-null the thin `U` (`m x k`) and `Vᴴ` (`k x n`).
fn gesvd(comptime T: type, allocator: Allocator, io: Io, m: usize, n: usize, x: []T, s: []f64, u: ?[]T, vt: ?[]T) Error!void {
    const k = @min(m, n);
    var none: [1]T = undefined;
    const job: lapack.Job = if (u != null) .s else .n;
    const up: [*]T = if (u) |y| y.ptr else &none;
    const vp: [*]T = if (vt) |y| y.ptr else &none;
    lapack.gesvd(allocator, io, T, job, job, @intCast(m), @intCast(n), x.ptr, @intCast(@max(1, m)), s.ptr, up, @intCast(if (u != null) @max(1, m) else 1), vp, @intCast(if (vt != null) @max(1, k) else 1), null) catch |e| switch (e) {
        error.NoConvergence, error.OutOfMemory, error.Canceled => |y| return y,
        // the arguments come from the shape
        else => unreachable,
    };
}

/// The singular values of `a` (not mutated), nonincreasing. Caller frees.
pub fn singularValues(a: CxMat, allocator: Allocator, io: Io) Error![]f64 {
    const s = try allocator.alloc(f64, @min(a.m, a.n));
    errdefer allocator.free(s);
    switch (a.isReal()) {
        inline else => |real| {
            const T = if (real) f64 else Complex;
            const x = try a.toPacked(T, allocator);
            defer allocator.free(x);
            try gesvd(T, allocator, io, a.m, a.n, x, s, null, null);
        },
    }
    return s;
}

/// COND: ratio of the largest to smallest singular value of `a` (2-norm
/// condition number). `a` is not mutated.
pub fn cond(allocator: Allocator, io: Io, a: CxMat) Error!f64 {
    const s = try singularValues(a, allocator, io);
    defer allocator.free(s);
    const smallest = s[s.len - 1];
    if (smallest == 0.0) return std.math.inf(f64);
    return s[0] / smallest;
}

/// The singular values above `tol`, `null` for the original's default
/// `max(m, n) * eps * s[0]`.
fn rankOf(s: []const f64, m: usize, n: usize, tol: ?f64) usize {
    const threshold = tol orelse (@as(f64, @floatFromInt(@max(m, n))) * std.math.floatEps(f64) * s[0]);
    for (s, 0..) |sv, k| if (sv <= threshold) return k;
    return s.len;
}

/// RANK: number of singular values of `a` that exceed `tol` (`null`: the
/// original's `max(m, n) * eps * largest_singular_value`). `a` is not
/// mutated.
pub fn rank(allocator: Allocator, io: Io, a: CxMat, tol: ?f64) Error!usize {
    const s = try singularValues(a, allocator, io);
    defer allocator.free(s);
    return rankOf(s, a.m, a.n, tol);
}

/// PINV: Moore-Penrose pseudoinverse of `a`, via SVD
/// (`pinv(A) = V * pinv(S) * ctrans(U)`, treating singular values at or
/// below `tol` as zero - `null` uses the same default as `rank`). Writes
/// the `n`x`m` result into caller-supplied `result_r`/`result_i`
/// (column-major, leading dimension `n`). `a` is not mutated.
pub fn pinv(allocator: Allocator, io: Io, a: CxMat, tol: ?f64, result_r: []f64, result_i: []f64) Error!void {
    const m = a.m;
    const n = a.n;
    const k = @min(m, n);
    const s = try allocator.alloc(f64, k);
    defer allocator.free(s);
    switch (a.isReal()) {
        inline else => |real| {
            const T = if (real) f64 else Complex;
            const x = try a.toPacked(T, allocator);
            defer allocator.free(x);
            const u = try allocator.alloc(T, m * k);
            defer allocator.free(u);
            const vt = try allocator.alloc(T, k * n);
            defer allocator.free(vt);
            try gesvd(T, allocator, io, m, n, x, s, u, vt);
            const r = rankOf(s, m, n, tol);
            const out = CxMat.initPacked(result_r, result_i, n, m);
            // pinv(i, j) = sum over l of conj(Vᴴ(l, i)) conj(U(j, l)) / s(l)
            for (1..n + 1) |i| {
                for (1..m + 1) |j| {
                    var acc = Cx.init(0, 0);
                    for (0..r) |l| {
                        const v = cx(vt[(i - 1) * k + l]).conjugate();
                        const w = cx(u[l * m + j - 1]).conjugate();
                        const t = v.mul(w);
                        acc = acc.add(Cx.init(t.re / s[l], t.im / s[l]));
                    }
                    out.set(i, j, acc);
                }
            }
        },
    }
}

fn cx(x: anytype) Cx {
    return if (@TypeOf(x) == f64) Cx.init(x, 0) else x;
}

/// ORTH: `Q` from a QR decomposition of `a` (`m`x`n`, `n <= m`, not
/// mutated), by `geqrf` and `orgqr` - an orthonormal basis (as columns) for
/// its column space. Caller owns the returned `m`x`n` slices.
pub fn orth(allocator: Allocator, io: Io, a: CxMat) (Allocator.Error || Io.Cancelable)!struct { r: []f64, i: []f64 } {
    const m = a.m;
    const n = a.n;
    const result_r = try allocator.alloc(f64, m * n);
    errdefer allocator.free(result_r);
    const result_i = try allocator.alloc(f64, m * n);
    errdefer allocator.free(result_i);
    switch (a.isReal()) {
        inline else => |real| {
            const T = if (real) f64 else Complex;
            const x = try a.toPacked(T, allocator);
            defer allocator.free(x);
            const tau = try allocator.alloc(T, n);
            defer allocator.free(tau);
            const mi: isize = @intCast(m);
            const ni: isize = @intCast(n);
            const lda: isize = @intCast(@max(1, m));
            lapack.geqrf(allocator, io, T, mi, ni, x.ptr, lda, tau.ptr) catch |e| switch (e) {
                error.OutOfMemory, error.Canceled => |y| return y,
                else => unreachable,
            };
            lapack.orgqr(allocator, io, T, mi, ni, ni, x.ptr, lda, tau.ptr) catch |e| switch (e) {
                error.OutOfMemory, error.Canceled => |y| return y,
                // `n <= m` is the caller's
                else => unreachable,
            };
            CxMat.initPacked(result_r, result_i, m, n).fromPacked(T, x);
        },
    }
    return .{ .r = result_r, .i = result_i };
}

const testing = std.testing;

test "cond, rank and pinv agree for a known rank-deficient matrix" {
    // A = [[1,2],[2,4]] (rank 1, columns parallel).
    var ar = [_]f64{ 1.0, 2.0, 2.0, 4.0 };
    var ai: [4]f64 = @splat(0.0);
    const a = CxMat.initPacked(&ar, &ai, 2, 2);

    const c = try cond(testing.allocator, testing.io, a);
    try zigstd.testExpectTrue(c > 1e10 or std.math.isInf(c));

    const r = try rank(testing.allocator, testing.io, a, null);
    try zigstd.testExpectEqual(@as(usize, 1), r);

    var pr: [4]f64 = @splat(0.0);
    var pi: [4]f64 = @splat(1.0);
    try pinv(testing.allocator, testing.io, a, null, &pr, &pi);
    // a real `A` gives an exactly real pinv
    try zigstd.testExpectEqualSlices(f64, &@as([4]f64, @splat(0.0)), &pi);
    // A * pinv(A) * A should reconstruct A (pseudoinverse identity).
    var apa_full: [4]f64 = @splat(0.0);
    for (0..2) |row| {
        for (0..2) |col| {
            var s: f64 = 0.0;
            for (0..2) |m| {
                var ap: f64 = 0.0;
                for (0..2) |l| ap += ar[ix(2, row + 1, l + 1)] * pr[ix(2, l + 1, m + 1)];
                s += ap * ar[ix(2, m + 1, col + 1)];
            }
            apa_full[ix(2, row + 1, col + 1)] = s;
        }
    }
    for (0..4) |idx| try zigstd.testExpectApproxEqAbs(ar[idx], apa_full[idx], 1e-12);
}

test "pinv of a complex rectangular matrix satisfies A pinv(A) A = A" {
    // A = [[1, i], [0, 1], [1, 0]] (3x2)
    var ar = [_]f64{ 1, 0, 1, 0, 1, 0 };
    var ai = [_]f64{ 0, 0, 0, 1, 0, 0 };
    const a = CxMat.initPacked(&ar, &ai, 3, 2);
    var pr: [6]f64 = undefined;
    var pi: [6]f64 = undefined;
    try pinv(testing.allocator, testing.io, a, null, &pr, &pi);
    const p = CxMat.initPacked(&pr, &pi, 2, 3);
    for (1..4) |i| {
        for (1..3) |j| {
            var s = Cx.init(0, 0);
            for (1..3) |l| {
                for (1..4) |q| s = s.add(a.get(i, l).mul(p.get(l, q)).mul(a.get(q, j)));
            }
            try zigstd.testExpectApproxEqAbs(a.get(i, j).re, s.re, 1e-14);
            try zigstd.testExpectApproxEqAbs(a.get(i, j).im, s.im, 1e-14);
        }
    }
}

test "orth returns an orthonormal basis spanning the column space" {
    // A = [[1,1],[0,1],[1,0]] (3x2, full column rank), then with an
    // imaginary entry
    for ([_]f64{ 0, 1 }) |im| {
        var ar = [_]f64{ 1.0, 0.0, 1.0, 1.0, 1.0, 0.0 };
        var ai = [_]f64{ 0, 0, 0, 0, im, 0 };
        const q = try orth(testing.allocator, testing.io, CxMat.initPacked(&ar, &ai, 3, 2));
        defer testing.allocator.free(q.r);
        defer testing.allocator.free(q.i);
        const qm = CxMat.initPacked(q.r, q.i, 3, 2);

        for (1..3) |a_col| {
            for (1..3) |b_col| {
                var dot = Cx.init(0, 0);
                for (1..4) |row| dot = dot.add(qm.get(row, a_col).conjugate().mul(qm.get(row, b_col)));
                const expected: f64 = if (a_col == b_col) 1.0 else 0.0;
                try zigstd.testExpectApproxEqAbs(expected, dot.re, 1e-14);
                try zigstd.testExpectApproxEqAbs(0.0, dot.im, 1e-14);
            }
        }
    }
}
