//! The eigendecomposition of a real symmetric tridiagonal matrix, from faer
//! v0.24.4 (<https://codeberg.org/sarah-quinones/faer>, tag `faer-v0.24.4`,
//! MIT): `faer/src/linalg/evd/tridiag_evd.rs`, QR iteration (`xSTEQR`'s
//! method) and Gu and Eisenstat's divide and conquer (`xSTEDC`'s).

const std = @import("std");
const mem = @import("../../../mem.zig");
const ops = @import("ops.zig");
const blas_calls = @import("blas_calls.zig");
const hh = @import("householder.zig");
const bsvd = @import("bidiag_svd.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = ops.Mat;

/// What the solvers return: `NoConvergence` once `qrAlgorithm` has run its
/// iterations out.
pub const Error = ops.Error || error{NoConvergence};

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
    const max = @max(ops.absMax(F, diag), ops.absMax(F, offdiag));
    if (max == 0) return;
    ops.scaleDown(F, diag, max);
    ops.scaleDown(F, offdiag, max);
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
            mem.swap(F, &diag[i], &diag[idx]);
            if (u) |um| ops.swapCols(F, um, i, idx);
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
    const dmax = ops.absMax(F, pd);
    const zmax = ops.absMax(F, pz);
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
            const inv_norm = 1 / blas_calls.nrm2(io, F, n, col, 1);
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
        mem.swap(Mat(F), &repaired_u, &tmp);
    }
    try blas_calls.gemm(al, io, F, .n, .n, 1, ua, repaired_u.sub(0, 0, n1, n), 0, tmp.sub(0, 0, n1, n));
    try blas_calls.gemm(al, io, F, .n, .n, 1, ub, repaired_u.sub(n1, 0, n - n1, n), 0, tmp.sub(n1, 0, n - n1, n));
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
