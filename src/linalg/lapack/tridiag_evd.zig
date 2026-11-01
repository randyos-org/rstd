//! The eigendecomposition of a real symmetric tridiagonal matrix, from faer
//! v0.24.4 (<https://codeberg.org/sarah-quinones/faer>, tag `faer-v0.24.4`,
//! MIT): `faer/src/linalg/evd/tridiag_evd.rs`. `diag` holds the diagonal,
//! `offdiag` the subdiagonal (`offdiag[n-1]`, when present, is 0).
//!
//! - `qrAlgorithm`: implicit Wilkinson-shift QR sweeps with Givens
//!   rotations; the reference's `xSTEQR` method.
//! - `divideAndConquer`: Cuppen's divide and conquer in Gu and Eisenstat's
//!   stable form ("A divide-and-conquer algorithm for the symmetric
//!   tridiagonal eigenproblem", SIAM J. Matrix Anal. Appl. 16(1), 1995), the
//!   reference's `xSTEDC`: split at `n / 2` by a rank-one tear, solve the
//!   halves, deflate small `z` entries and (by a Householder reflection per
//!   run) equal eigenvalues, solve the secular equation with the SVD's root
//!   finder (`bidiag_svd.rootFinder`), rebuild `z` from the roots (the
//!   reference's `xLAED3`, `sqrt(-W)`), and merge the vectors with `gemm`.
//!   Below `threshold` it is `qrAlgorithm`.
//!
//! Faer runs the halves on two threads (`join_raw`); here they run in turn,
//! the products threading inside the BLAS, as `bidiag_svd.zig` does.

const std = @import("std");
const ops = @import("ops.zig");
const hh = @import("householder.zig");
const bsvd = @import("bidiag_svd.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = ops.Mat;

/// What the solvers return: `NoConvergence` once `qrAlgorithm` has run its
/// iterations out.
pub const Error = ops.Error || error{NoConvergence};

fn normMax(comptime F: type, x: []const F) F {
    var m: F = 0;
    for (x) |v| m = @max(m, @abs(v));
    return m;
}

/// `x *= 1 / max`, by division where that reciprocal overflows (faer
/// multiplies, which a denormal `max` makes infinite; issues.md).
fn scaleDown(comptime F: type, x: []F, max: F) void {
    if (max >= std.math.floatMin(F)) {
        const inv = 1 / max;
        for (x) |*v| v.* *= inv;
    } else for (x) |*v| {
        v.* /= max;
    }
}

fn swapCols(comptime F: type, m: Mat(F), i: usize, j: usize) void {
    for (0..m.rows) |r| std.mem.swap(F, m.at(r, i), m.at(r, j));
}

/// `qr_algorithm`: the eigenvalues into `diag`, nondecreasing; `u`
/// (`n x n`), when given, set to the eigenvectors.
pub fn qrAlgorithm(comptime F: type, diag: []F, offdiag: []F, u: ?Mat(F)) error{NoConvergence}!void {
    const R = hh.Rotation(F);
    const n = diag.len;
    if (u) |um| um.identity();
    if (n <= 1) return;
    const eps = std.math.floatEps(F);
    const sml = std.math.floatMin(F);
    if (n == 2) {
        const a = diag[0];
        const d = diag[1];
        const b = offdiag[0];
        const t0 = std.math.hypot(a - d, b * 2) * 0.5;
        const t1 = (a + d) * 0.5;
        const r0 = t1 - t0;
        const r1 = t1 + t0;
        const tol = @max(@abs(r0), @abs(r1)) * eps;
        if (u) |um| {
            if (r1 - r0 <= tol) {
                um.identity();
            } else if (@abs(b) <= tol) {
                if (!(diag[0] < diag[1])) {
                    um.fill(0);
                    um.at(1, 0).* = 1;
                    um.at(0, 1).* = 1;
                }
            } else {
                const tau = ((d - a) / b) * 0.5;
                var t = 1 / (@abs(tau) + std.math.hypot(tau, @as(F, 1)));
                if (tau < 0) t = -t;
                var c = std.math.hypot(t, @as(F, 1));
                var s = c * t;
                const r = std.math.hypot(c, s);
                c /= r;
                s /= r;
                const r0_r = (c * a - s * b) / c;
                if (@abs(r0 - r0_r) < r1 - r0_r) {
                    um.at(0, 0).* = c;
                    um.at(1, 0).* = -s;
                    um.at(0, 1).* = s;
                    um.at(1, 1).* = c;
                } else {
                    um.at(0, 1).* = c;
                    um.at(1, 1).* = -s;
                    um.at(0, 0).* = s;
                    um.at(1, 0).* = c;
                }
            }
        }
        diag[0] = r0;
        diag[1] = r1;
        return;
    }
    const max = @max(normMax(F, diag), normMax(F, offdiag));
    if (max == 0) return;
    scaleDown(F, diag, max);
    scaleDown(F, offdiag, max);
    var end = n - 1;
    var start: usize = 0;
    const max_iters = @max(30, (std.math.floatMantissaBits(F) + 1) / 2) *| n *| n;
    var iter: usize = 0;
    while (iter < max_iters) : (iter += 1) {
        for (start..end) |i| {
            if (@abs(offdiag[i]) < sml or @abs(offdiag[i]) < eps * std.math.hypot(diag[i], diag[i + 1])) offdiag[i] = 0;
        }
        while (end > 0 and offdiag[end - 1] == 0) end -= 1;
        if (end == 0) break;
        if (iter + 1 == max_iters) {
            for (diag) |*x| x.* *= max;
            for (offdiag) |*x| x.* *= max;
            return error.NoConvergence;
        }
        start = end - 1;
        while (start > 0 and offdiag[start - 1] != 0) start -= 1;
        const td = (diag[end - 1] - diag[end]) * 0.5;
        const e = offdiag[end - 1];
        var mu = diag[end];
        if (td == 0) {
            mu -= @abs(e);
        } else if (e != 0) {
            const e2 = e * e;
            var hd = std.math.hypot(td, e);
            if (!(td > 0)) hd = -hd;
            mu = if (e2 == 0) mu - e / ((td + hd) / e) else mu - e2 / (td + hd);
        }
        var x = diag[start] - mu;
        var z = offdiag[start];
        var k = start;
        while (k < end and z != 0) : (k += 1) {
            const rot = R.make(x, z);
            const sdk = rot.s * diag[k] + rot.c * offdiag[k];
            const dkp1 = rot.s * offdiag[k] + rot.c * diag[k + 1];
            diag[k] = rot.c * (rot.c * diag[k] - rot.s * offdiag[k]) - rot.s * (rot.c * offdiag[k] - rot.s * diag[k + 1]);
            diag[k + 1] = rot.s * sdk + rot.c * dkp1;
            offdiag[k] = rot.c * sdk - rot.s * dkp1;
            if (k > start) offdiag[k - 1] = rot.c * offdiag[k - 1] - rot.s * z;
            x = offdiag[k];
            if (k < end - 1) {
                z = -rot.s * offdiag[k + 1];
                offdiag[k + 1] = rot.c * offdiag[k + 1];
            }
            if (u) |um| rot.right(um.rows, um.col(k + 1), um.col(k));
        }
    }
    for (0..n - 1) |i| {
        var idx = i;
        var min = diag[i];
        for (i + 1..n) |kk| if (diag[kk] < min) {
            idx = kk;
            min = diag[kk];
        };
        if (idx != i) {
            std.mem.swap(F, &diag[i], &diag[idx]);
            if (u) |um| swapCols(F, um, i, idx);
        }
    }
    for (diag) |*x| x.* *= max;
}

/// The secular function `1 / rho + Σ z² / ((d - shift) - mu)` over the
/// undeflated entries.
fn Secular(comptime F: type) type {
    return struct {
        d: []const F,
        z: []const F,
        rho_recip: F,

        const Self = @This();

        /// `secular_eq`.
        pub fn eval(e: Self, shift: F, mu: F) F {
            var res = e.rho_recip;
            for (e.d, e.z) |d, z| res += z * (z / ((d - shift) - mu));
            return res;
        }

        /// `batch_secular_eq`.
        pub fn batch(e: Self, shift: [4]F, mu: [4]F) [4]F {
            var res: [4]F = @splat(e.rho_recip);
            for (e.d, e.z) |d, z| {
                for (&res, mu, shift) |*r, m, s| r.* += z * (z / ((d - s) - m));
            }
            return res;
        }
    };
}

/// `compute_eigenvalues`: each root as `shift + mu`, `shift` its nearer
/// pole; the deflated past `nd` are their own.
fn eigenvalues(comptime F: type, mus: []F, shifts: []F, d: []const F, z: []const F, rho: F, nd: usize) void {
    var zz: F = 0;
    for (z) |x| zz += x * x;
    const sec: Secular(F) = .{ .d = d[0..nd], .z = z[0..nd], .rho_recip = 1 / rho };
    for (0..nd) |i| {
        const last = i == nd - 1;
        const right = if (last) d[i] + rho * zz else d[i + 1];
        const r = bsvd.rootFinder(F, sec, d[i], right, last);
        shifts[i] = r[0];
        mus[i] = r[1];
    }
    for (nd..d.len) |i| {
        shifts[i] = 0;
        mus[i] = d[i];
    }
}

/// `divide_and_conquer`'s workspace, split in halves for the halves.
fn Work(comptime F: type) type {
    return struct {
        pl_before: []usize,
        pl_after: []usize,
        pr: []usize,
        run_info: []usize,
        z: []F,
        permuted_diag: []F,
        permuted_z: []F,
        householder: []F,
        mus: []F,
        shifts: []F,
        repaired_u: Mat(F),
        tmp: Mat(F),

        const Self = @This();

        fn half(w: Self, lo: usize, hi: usize) Self {
            var r: Self = undefined;
            inline for (.{ "pl_before", "pl_after", "pr", "run_info", "z", "permuted_diag", "permuted_z", "householder", "mus", "shifts" }) |f| {
                @field(r, f) = @field(w, f)[lo..hi];
            }
            r.repaired_u = w.repaired_u.sub(lo, lo, hi - lo, hi - lo);
            r.tmp = w.tmp.sub(lo, lo, hi - lo, hi - lo);
            return r;
        }
    };
}

/// `divide_and_conquer`: the eigenvalues into `diag`, nondecreasing, and
/// the eigenvectors into `u` (`n x n`).
pub fn divideAndConquer(al: Allocator, io: Io, comptime F: type, diag: []F, offdiag: []F, u: Mat(F), threshold: usize) Error!void {
    const n = diag.len;
    const idx = try al.alloc(usize, 4 * n);
    defer al.free(idx);
    const vec = try al.alloc(F, 6 * n);
    defer al.free(vec);
    const ru = try ops.alloc(F, al, n, n);
    defer ops.free(F, al, ru);
    const tmp = try ops.alloc(F, al, n, n);
    defer ops.free(F, al, tmp);
    const w: Work(F) = .{
        .pl_before = idx[0..n],
        .pl_after = idx[n .. 2 * n],
        .pr = idx[2 * n .. 3 * n],
        .run_info = idx[3 * n ..],
        .z = vec[0..n],
        .permuted_diag = vec[n .. 2 * n],
        .permuted_z = vec[2 * n .. 3 * n],
        .householder = vec[3 * n .. 4 * n],
        .mus = vec[4 * n .. 5 * n],
        .shifts = vec[5 * n ..],
        .repaired_u = ru,
        .tmp = tmp,
    };
    try recurse(al, io, F, diag, offdiag, u, w, threshold);
}

/// `divide_and_conquer_recurse`.
fn recurse(al: Allocator, io: Io, comptime F: type, diag: []F, offdiag: []F, u: Mat(F), w: Work(F), threshold: usize) Error!void {
    const n = diag.len;
    if (n <= @max(threshold, 4)) return qrAlgorithm(F, diag, offdiag, u);
    const n1 = n / 2;
    var rho = offdiag[n1 - 1];
    diag[n1 - 1] -= @abs(rho);
    diag[n1] -= @abs(rho);
    const ua = u.sub(0, 0, n1, n1);
    const ub = u.sub(n1, n1, n - n1, n - n1);
    try recurse(al, io, F, diag[0..n1], offdiag[0 .. n1 - 1], ua, w.half(0, n1), threshold);
    try recurse(al, io, F, diag[n1..], offdiag[n1..], ub, w.half(n1, n), threshold);
    var repaired_u = w.repaired_u;
    var tmp = w.tmp;
    const z = w.z;
    const pd = w.permuted_diag;
    const pz = w.permuted_z;
    for (0..n1) |i| z[i] = ua.get(n1 - 1, i);
    for (0..n - n1) |i| z[n1 + i] = if (rho < 0) -ub.get(0, i) else ub.get(0, i);
    const inv_sqrt2 = @sqrt(@as(F, 0.5));
    for (z) |*x| x.* *= inv_sqrt2;
    rho = @abs(rho) * 2;
    // the halves' eigenvalues merged, nondecreasing
    {
        var i: usize = 0;
        var j: usize = 0;
        for (w.pl_before) |*p| {
            if (i == n1) {
                p.* = n1 + j;
                j += 1;
            } else if (j == n - n1 or diag[i] < diag[n1 + j]) {
                p.* = i;
                i += 1;
            } else {
                p.* = n1 + j;
                j += 1;
            }
        }
    }
    for (w.pl_before, 0..) |p, i| {
        pd[i] = diag[p];
        pz[i] = z[p];
    }
    const dmax = normMax(F, pd);
    const zmax = normMax(F, pz);
    const tol = 8 * std.math.floatEps(F) * @max(dmax, zmax);
    if (rho * zmax <= tol) {
        // all deflated: the halves' vectors, merged
        tmp.fill(0);
        tmp.sub(0, 0, n1, n1).copyFrom(ua);
        tmp.sub(n1, n1, n - n1, n - n1).copyFrom(ub);
        for (w.pl_before, 0..) |p, j| @memcpy(u.col(j)[0..n], tmp.col(p)[0..n]);
        @memcpy(diag, pd);
        return;
    }
    for (pz) |*x| if (@abs(rho * x.*) <= tol) {
        x.* = 0;
    };
    // runs of equal eigenvalues: one reflection each folds their z into one
    var applied_householder = false;
    var idx: usize = 0;
    while (idx < n) {
        var run_len: usize = 1;
        const d_prev = pd[idx];
        while (idx + run_len < n) {
            if (pd[idx + run_len] - d_prev <= tol) {
                pd[idx + run_len] = d_prev;
                run_len += 1;
            } else break;
        }
        @memset(w.run_info[idx..][0..run_len], run_len);
        if (run_len > 1) {
            applied_householder = true;
            const hv = w.householder[idx..][0..run_len];
            const zr = pz[idx..][0..run_len];
            @memcpy(hv, zr);
            // faer reverses the tail; the reflection's vector is the same
            const info = hh.make(io, F, &hv[run_len - 1], run_len - 1, hv.ptr, 1, null, 0);
            @memset(zr, 0);
            zr[run_len - 1] = hv[run_len - 1];
            hv[run_len - 1] = info.tau;
        }
        idx += run_len;
    }
    // deflated entries to the back
    var nd: usize = 0;
    var deflated: usize = 0;
    for (0..n) |reader| {
        const zr = pz[reader];
        const d = pd[reader];
        if (zr == 0) {
            diag[deflated] = d;
            w.pr[deflated] = reader;
            deflated += 1;
        } else {
            pz[nd] = zr;
            pd[nd] = d;
            w.pl_after[nd] = reader;
            nd += 1;
        }
    }
    for (0..deflated) |i| {
        pd[nd + i] = diag[i];
        w.pl_after[nd + i] = w.pr[i];
    }
    eigenvalues(F, w.mus, w.shifts, pd, pz, rho, nd);
    // z from the roots (Gu and Eisenstat; the reference's `xLAED3`)
    for (0..nd) |i| {
        const di = pd[i];
        var prod = w.mus[i] + (w.shifts[i] - di);
        for (0..nd) |kk| {
            if (kk == i) continue;
            prod *= (w.mus[kk] + (w.shifts[kk] - di)) / (pd[kk] - di);
        }
        // faer takes |prod|, where `xLAED3` takes sqrt(-W) (NaN for a
        // spuriously positive W; OpenBLAS's laed3 clamps it to 0)
        const p = @sqrt(@abs(prod));
        pz[i] = if (pz[i] < 0) -p else p;
    }
    const eigenvals = z;
    for (0..n) |i| eigenvals[i] = w.mus[i] + w.shifts[i];
    for (w.pr, 0..) |*p, i| p.* = i;
    std.mem.sort(usize, w.pr, eigenvals, lessNanLast(F));
    if (!applied_householder) {
        for (w.pl_after) |*p| p.* = w.pl_before[p.*];
    }
    for (w.pr, 0..) |pj, j| {
        const col = repaired_u.col(j);
        if (pj >= nd) {
            @memset(col[0..n], 0);
            col[w.pl_after[pj]] = 1;
        } else {
            const mu_j = w.mus[pj];
            const shift_j = w.shifts[pj];
            for (w.pl_after[0..nd], 0..) |pa, i| col[pa] = pz[i] / ((pd[i] - shift_j) - mu_j);
            for (w.pl_after[nd .. nd + deflated]) |pa| col[pa] = 0;
            const inv_norm = 1 / ops.nrm2(io, F, n, col, 1);
            for (col[0..n]) |*x| x.* *= inv_norm;
        }
    }
    if (applied_householder) {
        idx = 0;
        while (idx < n) {
            const run_len = w.run_info[idx];
            if (run_len > 1) {
                const hv = w.householder[idx..][0..run_len];
                const tau_inv = 1 / hv[run_len - 1];
                hv[run_len - 1] = 1;
                for (0..n) |j| {
                    const col = repaired_u.col(j) + idx;
                    var dot: F = 0;
                    for (0..run_len) |i| dot += hv[i] * col[i];
                    dot *= tau_inv;
                    for (0..run_len) |i| col[i] -= dot * hv[i];
                }
            }
            idx += run_len;
        }
        for (0..n) |j| for (w.pl_before, 0..) |p, i| {
            tmp.at(p, j).* = repaired_u.get(i, j);
        };
        std.mem.swap(Mat(F), &repaired_u, &tmp);
    }
    try ops.gemm(al, io, F, .n, .n, 1, ua, repaired_u.sub(0, 0, n1, n), 0, tmp.sub(0, 0, n1, n));
    try ops.gemm(al, io, F, .n, .n, 1, ub, repaired_u.sub(n1, 0, n - n1, n), 0, tmp.sub(n1, 0, n - n1, n));
    u.copyFrom(tmp);
    for (0..n) |i| diag[i] = w.mus[w.pr[i]] + w.shifts[w.pr[i]];
}

/// `pr`'s order: by value, NaN last.
fn lessNanLast(comptime F: type) fn ([]const F, usize, usize) bool {
    return struct {
        fn f(v: []const F, i: usize, j: usize) bool {
            const a = v[i];
            const b = v[j];
            if (std.math.isNan(a)) return false;
            if (std.math.isNan(b)) return true;
            return a < b;
        }
    }.f;
}

// ---- tests -----------------------------------------------------------------------

const testing = std.testing;
const check = @import("../blas/_test/blas1/check.zig");

const Solver = enum { qr, dc4, dc8 };

/// `U diag(s) Uᵀ` against the tridiagonal per element, faer's
/// `ApproxEq::eps()` (`128 eps`) times `max(1, ‖T‖)`; `U` orthogonal
/// alike, `s` nondecreasing.
fn expectTridiag(comptime F: type, al: Allocator, diag: []const F, off: []const F, solver: Solver) !void {
    const n = diag.len;
    const d = try al.dupe(F, diag);
    defer al.free(d);
    const e = try al.alloc(F, n);
    defer al.free(e);
    @memcpy(e[0 .. n - 1], off);
    e[n - 1] = 0;
    const u = try ops.alloc(F, al, n, n);
    defer ops.free(F, al, u);
    u.fill(std.math.nan(F));
    switch (solver) {
        .qr => try qrAlgorithm(F, d, e, u),
        .dc4 => try divideAndConquer(al, check.io, F, d, e, u, 4),
        .dc8 => try divideAndConquer(al, check.io, F, d, e, u, 8),
    }
    var tmax: F = 0;
    for (diag) |x| tmax = @max(tmax, @abs(x));
    for (off) |x| tmax = @max(tmax, @abs(x));
    const tol = 128 * std.math.floatEps(F) * @max(tmax, 1);
    for (1..n) |i| try testing.expect(d[i - 1] <= d[i]);
    for (0..n) |j| for (0..n) |i| {
        var s: F = 0;
        var g: F = 0;
        for (0..n) |k| {
            s += u.get(i, k) * d[k] * u.get(j, k);
            g += u.get(k, i) * u.get(k, j);
        }
        const want: F = if (i == j) diag[j] else if (i == j + 1) off[j] else if (j == i + 1) off[i] else 0;
        const id: F = if (i == j) 1 else 0;
        if (!(@abs(s - want) <= tol) or !(@abs(g - id) <= tol * @sqrt(@as(F, @floatFromInt(n))))) {
            check.say("tridiag evd {s} {t} n={d}: ({d},{d}) U S Ut {e} want {e}, Ut U {e}\n", .{ @typeName(F), solver, n, i, j, s, want, g });
            return error.TestUnexpectedResult;
        }
    };
}

test "faer tridiag evd: faer's evd_qr_tests through QR iteration and divide and conquer" {
    const al = testing.allocator;
    const x = 1e-6;
    const cases = [_]struct { d: []const f64, e: []const f64 }{
        .{ .d = &.{ 1, 1 }, .e = &.{0} },
        .{ .d = &.{ 1, 1 }, .e = &.{0.5213289} },
        .{ .d = &.{ 1.79069356, 1.20930644, 1 }, .e = &.{ -4.06813537e-01, 0 } },
        .{ .d = &.{ 1.95069537, 2.44845332, 2.56957029, 3.03128102, 1 }, .e = &.{ -7.02200909e-01, -1.11661820e+00, -6.81418803e-01, 0 } },
        .{ .d = &.{ 3, 2, 1, 0, 1, 2, 3 }, .e = &.{ 1, 1, 1, 1, 1, 1 } },
        .{ .d = &.{ 3, 2, 1, 0, 1, 2, 3, 3, 2, 1, 0, 1, 2, 3 }, .e = &.{ 1, 1, 1, 1, 1, 1, x, 1, 1, 1, 1, 1, 1 } },
        .{ .d = &.{ 0, 0, 1.0769230769230773, -0.4290761869709236, -0.8278050499098524, 0.07994922044020283, -0.35579623371016944, 0.6487378508167678, -0.9347442346214521, -0.08624745233962683, -0.4999243909534632, 1.3708277457481026, -0.2592167303689501, -0.5929351972647323, -0.5863220906879729, 0.15069873027683844, 0.2449309426221532, 0.5599151389028441, 0.440084861097156, 9.811634162559901e-17 }, .e = &.{ 1.7320508075688772, 2.081665999466133, 2.0303418353670932, 1.2463948607107287, 1.5895840148470526, 1.3810057029812097, 1.265168346300635, 0.8941431038915991, 1.007512301091709, 0.5877505835309086, 1.0370970338888965, 0.8628932798233644, 1.1935059937001073, 1.1614143449715744, 0.41040224297074174, 0.561318309959268, 3.1807090401145072e-15, 0.4963971959331084, 1.942890293094024e-16 } },
        .{ .d = &.{ 1.0000000000000002, 1.0000000000000002 }, .e = &.{7.216449660063518e-16} },
    };
    for (cases) |c| for ([_]Solver{ .qr, .dc4, .dc8 }) |s| try expectTridiag(f64, al, c.d, c.e, s);
}

test "faer tridiag evd: random, clustered, graded and zero tridiagonals, f32 and f64" {
    const al = testing.allocator;
    inline for (.{ f32, f64 }) |F| {
        for ([_]usize{ 3, 9, 33, 100 }) |n| for (0..5) |kind| {
            const d = try al.alloc(F, n);
            defer al.free(d);
            const e = try al.alloc(F, n - 1);
            defer al.free(e);
            check.fill(F, d, n + kind);
            check.fill(F, e, 2 * n + kind);
            switch (kind) {
                0 => {},
                // clustered: equal diagonal, tiny coupling
                1 => {
                    @memset(d, 1);
                    for (e) |*x| x.* *= std.math.floatEps(F);
                },
                // graded over the type's range of significant digits
                2 => for (d, 0..) |*x, i| {
                    x.* = std.math.pow(F, 10, -@as(F, @floatFromInt(i % 12)));
                },
                // a zero block in the middle
                3 => {
                    @memset(d[n / 3 .. n / 2], 0);
                    @memset(e[n / 3 .. n / 2], 0);
                },
                // all zero
                else => {
                    @memset(d, 0);
                    @memset(e, 0);
                },
            }
            for ([_]Solver{ .qr, .dc4, .dc8 }) |s| try expectTridiag(F, al, d, e, s);
        };
    }
}
