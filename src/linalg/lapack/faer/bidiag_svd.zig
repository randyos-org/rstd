//! The SVD of a real bidiagonal matrix, from faer v0.24.4
//! (<https://codeberg.org/sarah-quinones/faer>, tag `faer-v0.24.4`, MIT):
//! `faer/src/linalg/svd/bidiag_svd.rs`. The matrix is lower bidiagonal,
//! `n + 1` by `n`: `diag` on its diagonal, `subdiag` below it (`subdiag[n-1]`
//! in row `n`), so `U` is `(n + 1) x (n + 1)` and `V` `n x n`.
//!
//! - `qrAlgorithm`: implicit-shift QR sweeps with Givens rotations, after
//!   Golub and Kahan; the sweeps take the Wilkinson shift of the trailing
//!   `2 x 2` of `Bᵀ B`.
//! - `divideAndConquer`: Gu and Eisenstat's divide and conquer ("A
//!   divide-and-conquer algorithm for the bidiagonal SVD", SIAM J. Matrix
//!   Anal. Appl. 16(1), 1995), as Eigen's `BDCSVD`: split at `n / 2`, solve
//!   the halves, deflate, solve the secular equation, and rebuild the
//!   vectors with `gemm`. Below `threshold` it is `qrAlgorithm`.
//!
//! Faer runs the two halves on two threads (`join_raw`); here they run in
//! turn, and the products that merge them thread inside the BLAS. Rstd's
//! threading splits flat ranges into chunks (`threading.md`), not a tree of
//! tasks; the halves' own products already fill the workers.

const std = @import("std");
const ops = @import("ops.zig");
const hh = @import("householder.zig");

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

/// `nbits`: the significand's bits.
fn nbits(comptime F: type) usize {
    return std.math.floatMantissaBits(F) + 1;
}

/// `qr_algorithm`: `diag` and `subdiag` (`n` each; `subdiag[n-1]` must be 0)
/// to the singular values, nonnegative and nonincreasing, rotating `u`'s and
/// `v`'s columns (`u.rows`, `v.rows` rows; `n` columns used).
pub fn qrAlgorithm(comptime F: type, diag: []F, subdiag: []F, u: ?Mat(F), v: ?Mat(F)) error{NoConvergence}!void {
    const R = hh.Rotation(F);
    const n = diag.len;
    const eps = std.math.floatEps(F);
    const sml = std.math.floatMin(F);
    if (n == 0) return;
    const max_iters = @max(30, nbits(F) / 2) *| n *| n;
    const last = n - 1;
    const max = @max(normMax(F, diag), normMax(F, subdiag));
    if (max == 0) return;
    scaleDown(F, diag, max);
    scaleDown(F, subdiag, max);
    const eps2 = eps * eps;
    var iter: usize = 0;
    while (iter < max_iters) : (iter += 1) {
        for (0..last) |ia| {
            if (subdiag[ia] * subdiag[ia] <= eps2 * @abs(diag[ia] * diag[ia + 1]) + sml) subdiag[ia] = 0;
        }
        var end = n;
        while (end >= 2 and subdiag[end - 2] * subdiag[end - 2] <= sml) end -= 1;
        if (end == 1) break;
        var start = end - 1;
        while (start >= 1 and subdiag[start - 1] != 0) start -= 1;
        var found_zero_diag = false;
        for (start..end - 1) |i| {
            if (@abs(diag[i]) <= eps) {
                found_zero_diag = true;
                var val = subdiag[i];
                subdiag[i] = 0;
                for (i + 1..end) |j| {
                    const rot = R.make(diag[j], val);
                    diag[j] = rot.c * diag[j] - rot.s * val;
                    if (j + 1 < end) {
                        val = rot.s * subdiag[j];
                        subdiag[j] = rot.c * subdiag[j];
                    }
                    if (v) |vm| rot.adjoint().right(vm.rows, vm.col(j), vm.col(i));
                }
            }
        }
        if (found_zero_diag) {
            if (iter + 1 == max_iters) return unscaled(F, diag, subdiag, max);
            continue;
        }
        const end2 = end - 2;
        const end1 = end - 1;
        const t00 = if (end - start == 2) diag[end2] * diag[end2] else diag[end2] * diag[end2] + subdiag[end - 3] * subdiag[end - 3];
        const t11 = diag[end1] * diag[end1] + subdiag[end2] * subdiag[end2];
        const t01 = diag[end2] * subdiag[end2];
        const t01_2 = t01 * t01;
        var mu: F = undefined;
        if (t01_2 > sml) {
            const d = (t00 - t11) * 0.5;
            var delta = @sqrt(d * d + t01_2);
            if (d < 0) delta = -delta;
            mu = t11 - t01_2 / (d + delta);
        } else {
            mu = t11;
        }
        var y = diag[start] * diag[start] - mu;
        var z = diag[start] * subdiag[start];
        for (start..end1) |k| {
            var rot = R.make(y, z);
            if (k > start) subdiag[k - 1] = @abs(rot.c * y - rot.s * z);
            var diag_k = diag[k];
            const t0 = rot.c * diag_k - rot.s * subdiag[k];
            const t1 = rot.s * diag_k + rot.c * subdiag[k];
            diag_k = t0;
            subdiag[k] = t1;
            const k1 = k + 1;
            y = diag_k;
            z = -rot.s * diag[k1];
            diag[k1] = rot.c * diag[k1];
            if (u) |um| rot.right(um.rows, um.col(k1), um.col(k));
            rot = R.make(y, z);
            diag_k = rot.c * y - rot.s * z;
            diag[k] = diag_k;
            const s0 = rot.c * subdiag[k] - rot.s * diag[k1];
            const s1 = rot.s * subdiag[k] + rot.c * diag[k1];
            subdiag[k] = s0;
            diag[k1] = s1;
            if (k < end - 2) {
                y = subdiag[k];
                z = -rot.s * subdiag[k1];
                subdiag[k1] = rot.c * subdiag[k1];
            }
            if (v) |vm| rot.right(vm.rows, vm.col(k1), vm.col(k));
        }
        if (iter + 1 == max_iters) return unscaled(F, diag, subdiag, max);
    }
    for (0..n) |j| if (diag[j] < 0) {
        diag[j] = -diag[j];
        if (v) |vm| for (vm.col(j)[0..vm.rows]) |*x| {
            x.* = -x.*;
        };
    };
    for (0..n) |k| {
        var m: F = 0;
        var idx = k;
        for (k..n) |kk| if (diag[kk] > m) {
            m = diag[kk];
            idx = kk;
        };
        if (k != idx) {
            std.mem.swap(F, &diag[k], &diag[idx]);
            if (u) |um| swapCols(F, um, k, idx);
            if (v) |vm| swapCols(F, vm, k, idx);
        }
    }
    for (diag) |*x| x.* *= max;
    for (subdiag) |*x| x.* *= max;
}

/// `x /= max`, as faer's `x *= 1 / max`, but by division where that
/// reciprocal overflows (`max` below the safe minimum).
fn scaleDown(comptime F: type, x: []F, max: F) void {
    if (max >= std.math.floatMin(F)) {
        const inv = 1 / max;
        for (x) |*v| v.* *= inv;
    } else for (x) |*v| {
        v.* /= max;
    }
}

/// The scaling undone, and `NoConvergence`.
fn unscaled(comptime F: type, diag: []F, subdiag: []F, max: F) error{NoConvergence} {
    for (diag) |*x| x.* *= max;
    for (subdiag) |*x| x.* *= max;
    return error.NoConvergence;
}

fn swapCols(comptime F: type, m: Mat(F), i: usize, j: usize) void {
    for (0..m.rows) |r| std.mem.swap(F, m.at(r, i), m.at(r, j));
}

/// `U` as `divide_and_conquer` takes it: whole, or only its first and last
/// rows (`2 x (n + 1)`), or none but scratch for those two rows.
pub fn MatU(comptime F: type) type {
    return union(enum) {
        full: Mat(F),
        two_rows: Mat(F),
        two_rows_storage: Mat(F),
    };
}

/// `divide_and_conquer`: the SVD of the `(n + 1) x n` lower bidiagonal
/// matrix, singular values into `diag`, nonincreasing; `subdiag` is
/// scratch. `v`, when given, is `n x n`.
pub fn divideAndConquer(al: Allocator, io: Io, comptime F: type, diag: []F, subdiag: []F, u: MatU(F), v: ?Mat(F), threshold: usize) Error!void {
    const R = hh.Rotation(F);
    const thr = @max(threshold, 4);
    const n = diag.len;
    if (n < thr) {
        const ua = try ops.alloc(F, al, n + 1, n + 1);
        defer ops.free(F, al, ua);
        ua.identity();
        if (v) |vm| vm.identity();
        var val = subdiag[n - 1];
        subdiag[n - 1] = 0;
        var i = n;
        while (i > 0) {
            i -= 1;
            const rot = R.make(diag[i], val);
            diag[i] = rot.c * diag[i] - rot.s * val;
            if (i > 0) {
                val = rot.s * subdiag[i - 1];
                subdiag[i - 1] = rot.c * subdiag[i - 1];
            }
            rot.adjoint().right(n + 1, ua.col(i), ua.col(n));
        }
        try qrAlgorithm(F, diag, subdiag, ua, v);
        switch (u) {
            .full => |um| um.sub(0, 0, n + 1, n + 1).copyFrom(ua),
            .two_rows => |um| for (0..n + 1) |j| {
                um.at(0, j).* = ua.get(0, j);
                um.at(1, j).* = ua.get(n, j);
            },
            .two_rows_storage => {},
        }
        return;
    }
    const max = @max(normMax(F, diag), normMax(F, subdiag));
    if (max == 0) {
        switch (u) {
            .full => |um| um.identity(),
            .two_rows => |um| {
                um.fill(0);
                um.at(0, 0).* = 1;
                // faer stores 0 here: the last row of the identity is eₙ
                um.at(1, n).* = 1;
            },
            .two_rows_storage => {},
        }
        if (v) |vm| vm.identity();
        return;
    }
    scaleDown(F, diag, max);
    scaleDown(F, subdiag, max);
    const k = n / 2;
    const rem = n - k - 1;
    const alpha = diag[k];
    const beta = subdiag[k];
    {
        const u_top: MatU(F), const u_bot: MatU(F) = switch (u) {
            .full => |um| .{ .{ .full = um.sub(0, 1, k + 1, k + 1) }, .{ .full = um.sub(k + 1, k + 1, rem + 1, rem + 1) } },
            .two_rows, .two_rows_storage => |um| .{ .{ .two_rows = um.sub(0, 0, 2, k + 1) }, .{ .two_rows = um.sub(0, k + 1, 2, rem + 1) } },
        };
        const v_top: ?Mat(F) = if (v) |vm| vm.sub(0, 1, k, k) else null;
        const v_bot: ?Mat(F) = if (v) |vm| vm.sub(k + 1, k + 1, rem, rem) else null;
        try divideAndConquer(al, io, F, diag[0..k], subdiag[0..k], u_top, v_top, threshold);
        try divideAndConquer(al, io, F, diag[k + 1 ..], subdiag[k + 1 ..], u_bot, v_bot, threshold);
        switch (u_top) {
            .two_rows => |um| {
                var i = k;
                while (i > 0) {
                    i -= 1;
                    swapCols(F, um, i, i + 1);
                }
            },
            .full, .two_rows_storage => {},
        }
    }
    if (v) |vm| vm.at(k, 0).* = 1;
    {
        var i = k;
        while (i > 0) {
            i -= 1;
            diag[i + 1] = diag[i];
        }
    }
    const lambda = switch (u) {
        .full => |um| um.get(k, k + 1),
        .two_rows, .two_rows_storage => |um| um.get(1, 0),
    };
    const phi = switch (u) {
        .full => |um| um.get(k + 1, n),
        .two_rows, .two_rows_storage => |um| um.get(0, n),
    };
    const al_ = alpha * lambda;
    const bp = beta * phi;
    const r0 = @sqrt(al_ * al_ + bp * bp);
    const c0: F, const s0: F = if (r0 == 0) .{ 1, 0 } else .{ al_ / r0, bp / r0 };
    const col0 = subdiag;
    diag[0] = r0;
    col0[0] = r0;
    switch (u) {
        .full => |um| {
            for (1..k + 1) |j| col0[j] = alpha * um.get(k, j);
            for (k + 1..n) |j| col0[j] = beta * um.get(k + 1, j);
            for (0..k + 1) |i| {
                const a = um.get(i, k + 1);
                um.at(i, 0).* = c0 * a;
                um.at(i, n).* = -s0 * a;
            }
            for (k + 1..n + 1) |i| {
                const xn = um.get(i, n);
                um.at(i, 0).* = s0 * xn;
                um.at(i, n).* = c0 * xn;
            }
        },
        .two_rows => |um| {
            for (1..k + 1) |j| {
                col0[j] = alpha * um.get(1, j);
                um.at(1, j).* = 0;
            }
            for (k + 1..n) |j| {
                col0[j] = beta * um.get(0, j);
                um.at(0, j).* = 0;
            }
            const q10 = um.get(0, 0);
            const q21 = um.get(1, n);
            um.at(0, 0).* = c0 * q10;
            um.at(0, n).* = -s0 * q10;
            um.at(1, 0).* = s0 * q21;
            um.at(1, n).* = c0 * q21;
        },
        .two_rows_storage => |um| {
            for (1..k + 1) |j| col0[j] = alpha * um.get(1, j);
            for (k + 1..n) |j| col0[j] = beta * um.get(0, j);
        },
    }
    const perm = try al.alloc(usize, n);
    defer al.free(perm);
    for (perm, 0..) |*p, i| p.* = i;
    const jacobi_coeff = try al.alloc(R, n);
    defer al.free(jacobi_coeff);
    @memset(jacobi_coeff, .{ .c = 0, .s = 0 });
    const jacobi_idx = try al.alloc(usize, n);
    defer al.free(jacobi_idx);
    for (jacobi_idx, 0..) |*p, i| p.* = i;
    const jacobi_0i: usize, const jacobi_ij: usize = blk: {
        const transpositions = try al.alloc(usize, n);
        defer al.free(transpositions);
        for (transpositions, 0..) |*p, i| p.* = i;
        break :blk try deflate(al, F, diag, col0, jacobi_coeff, jacobi_idx, transpositions, perm, k);
    };
    const um = try ops.alloc(F, al, n + 1, n + 1);
    defer ops.free(F, al, um);
    um.fill(0);
    const vm: ?Mat(F) = if (v != null) try ops.alloc(F, al, n, n) else null;
    defer if (vm) |x| ops.free(F, al, x);
    if (vm) |x| x.fill(0);
    try svdOfM(al, io, F, um, vm, diag, col0, perm);
    @memset(col0, 0);
    {
        var t = jacobi_ij;
        while (t > 0) {
            t -= 1;
            const rot = jacobi_coeff[jacobi_0i + t];
            const ii = jacobi_idx[jacobi_0i + t];
            const ai = perm[ii - 1];
            const aj = perm[ii];
            rot.left(n, um.ptr + aj, um.ld, um.ptr + ai, um.ld);
            if (vm) |x| rot.left(n, x.ptr + aj, x.ld, x.ptr + ai, x.ld);
        }
        t = jacobi_0i;
        while (t > 0) {
            t -= 1;
            const rot = jacobi_coeff[t];
            const ii = jacobi_idx[t];
            rot.left(n, um.ptr + ii, um.ld, um.ptr, um.ld);
        }
    }
    if (v) |vl| try updateV(al, io, F, vl, vm.?, k, rem);
    switch (u) {
        .two_rows_storage => {},
        .two_rows => |ul| {
            const cu = try ops.alloc(F, al, 2, n + 1);
            defer ops.free(F, al, cu);
            try ops.gemm(al, io, F, .n, .n, 1, ul, um, 0, cu);
            ul.copyFrom(cu);
        },
        .full => |ul| try updateU(al, io, F, ul, um, k, rem),
    }
    for (diag) |*x| x.* *= max;
}

/// `update_v`: `V := V Vm`, by its blocks.
fn updateV(al: Allocator, io: Io, comptime F: type, v: Mat(F), vm: Mat(F), k: usize, rem: usize) ops.Error!void {
    const n = v.rows;
    const cv = try ops.alloc(F, al, n, n);
    defer ops.free(F, al, cv);
    try ops.gemm(al, io, F, .n, .n, 1, v.sub(0, 1, k, k), vm.sub(1, 0, k, n), 0, cv.sub(0, 0, k, n));
    try ops.gemm(al, io, F, .n, .n, 1, v.sub(k + 1, k + 1, rem, rem), vm.sub(k + 1, 0, rem, n), 0, cv.sub(k + 1, 0, rem, n));
    try ops.gemm(al, io, F, .n, .n, 1, v.sub(k, 0, 1, 1), vm.sub(0, 0, 1, n), 0, cv.sub(k, 0, 1, n));
    v.copyFrom(cv);
}

/// `update_u`: `U := U Um`, by its blocks.
fn updateU(al: Allocator, io: Io, comptime F: type, u: Mat(F), um: Mat(F), k: usize, rem: usize) ops.Error!void {
    const n = u.rows - 1;
    const cu = try ops.alloc(F, al, n + 1, n + 1);
    defer ops.free(F, al, cu);
    const c1 = cu.sub(0, 0, k + 1, n + 1);
    const c2 = cu.sub(k + 1, 0, rem + 1, n + 1);
    try ops.gemm(al, io, F, .n, .n, 1, u.sub(0, 0, k + 1, k + 1), um.sub(0, 0, k + 1, n + 1), 0, c1);
    try ops.gemm(al, io, F, .n, .n, 1, u.sub(0, n, k + 1, 1), um.sub(n, 0, 1, n + 1), 1, c1);
    try ops.gemm(al, io, F, .n, .n, 1, u.sub(k + 1, k + 1, rem + 1, rem + 1), um.sub(k + 1, 0, rem + 1, n + 1), 0, c2);
    try ops.gemm(al, io, F, .n, .n, 1, u.sub(k + 1, 0, rem + 1, 1), um.sub(0, 0, 1, n + 1), 1, c2);
    u.copyFrom(cu);
}

/// `compute_svd_of_m`: the SVD of the deflated `M` (`col0` in its first
/// column, `diag` on its diagonal), vectors into `um` and `vm` (zeroed).
fn svdOfM(al: Allocator, io: Io, comptime F: type, um: Mat(F), vm: ?Mat(F), diag: []F, col0: []const F, outer_perm: []const usize) ops.Error!void {
    const n = diag.len;
    diag[0] = 0;
    var actual_n = n;
    while (actual_n > 1 and diag[actual_n - 1] == 0) actual_n -= 1;
    const perm_buf = try al.alloc(usize, actual_n);
    defer al.free(perm_buf);
    var o: usize = 0;
    for (0..actual_n) |i| if (@abs(col0[i]) != 0) {
        perm_buf[o] = i;
        o += 1;
    };
    const perm = perm_buf[0..o];
    const scratch = try al.alloc(F, 2 * o + 4 * n);
    defer al.free(scratch);
    const col0_perm = scratch[0..o];
    const diag_perm = scratch[o .. 2 * o];
    const shifts = scratch[2 * o ..][0..n];
    const mus = scratch[2 * o + n ..][0..n];
    const s = scratch[2 * o + 2 * n ..][0..n];
    const zhat = scratch[2 * o + 3 * n ..][0..n];
    for (perm, 0..) |p, kk| {
        col0_perm[kk] = col0[p];
        diag_perm[kk] = diag[p];
    }
    const sec: Secular(F) = .{ .col0_perm = col0_perm, .diag_perm = diag_perm };
    singularValues(io, F, shifts, mus, s, diag, col0, sec);
    perturbCol0(F, zhat, col0, diag, perm, s, shifts, mus);
    const col_perm = try al.alloc(usize, 2 * actual_n);
    defer al.free(col_perm);
    const col_perm_inv = col_perm[actual_n..];
    for (col_perm[0..actual_n], 0..) |*p, i| p.* = i;
    for (0..actual_n - 1) |ia| {
        if (s[ia] > s[ia + 1]) {
            std.mem.swap(F, &s[ia], &s[ia + 1]);
            std.mem.swap(usize, &col_perm[ia], &col_perm[ia + 1]);
        }
    }
    for (col_perm[0..actual_n], 0..) |p, i| col_perm_inv[p] = i;
    singularVectors(io, F, um, vm, zhat, diag, perm, outer_perm, col_perm_inv, actual_n, shifts, mus);
    for (0..actual_n) |idx| diag[idx] = s[actual_n - idx - 1];
    for (actual_n..n) |idx| diag[idx] = s[idx];
}

/// `compute_singular_vectors`.
fn singularVectors(io: Io, comptime F: type, um: Mat(F), vm: ?Mat(F), zhat: []const F, diag: []const F, perm: []const usize, outer_perm: []const usize, col_perm_inv: []const usize, actual_n: usize, shifts: []const F, mus: []const F) void {
    const n = diag.len;
    for (0..n) |k| {
        const actual_k = if (k >= actual_n) k else actual_n - col_perm_inv[k] - 1;
        um.at(n, actual_k).* = 0;
        const u = um.col(actual_k);
        const v: ?[*]F = if (vm) |x| x.col(actual_k) else null;
        if (zhat[k] == 0) {
            u[outer_perm[k]] = 1;
            if (v) |vc| vc[outer_perm[k]] = 1;
            continue;
        }
        const mu = mus[k];
        const shift = shifts[k];
        for (perm) |i| u[outer_perm[i]] = (zhat[i] / ((diag[i] - shift) - mu)) / (diag[i] + (shift + mu));
        const nu = 1 / ops.nrm2(io, F, n, u, 1);
        for (u[0..n]) |*x| x.* *= nu;
        if (v) |vc| {
            vc[outer_perm[0]] = -1;
            for (perm[1..]) |i| vc[outer_perm[i]] = ((diag[i] * zhat[i]) / ((diag[i] - shift) - mu)) / (diag[i] + (shift + mu));
            const nv = 1 / ops.nrm2(io, F, n, vc, 1);
            for (vc[0..n]) |*x| x.* *= nv;
        }
    }
    um.at(n, n).* = 1;
}

/// `perturb_col0`: `col0` recomputed from the computed singular values,
/// so the vectors come out orthogonal (Gu and Eisenstat).
fn perturbCol0(comptime F: type, zhat: []F, col0: []const F, diag: []const F, perm: []const usize, s: []const F, shifts: []const F, mus: []const F) void {
    const n = diag.len;
    const m = perm.len;
    if (m == 0) {
        @memset(zhat, 0);
        return;
    }
    const last_idx = perm[m - 1];
    for (0..n) |k| {
        if (col0[k] == 0) {
            zhat[k] = 0;
            continue;
        }
        const dk = diag[k];
        var prod = (s[last_idx] + dk) * (mus[last_idx] + (shifts[last_idx] - dk));
        for (perm, 0..) |i, l| {
            if (i == k) continue;
            if (i >= k and l == 0) {
                prod = 0;
                break;
            }
            const j = if (i < k) i else if (l > 0) perm[l - 1] else i;
            const term = ((s[j] + dk) / (diag[i] + dk)) * ((mus[j] + (shifts[j] - dk)) / (diag[i] - dk));
            prod = prod * term;
        }
        const tmp = @sqrt(prod);
        zhat[k] = if (col0[k] > 0) tmp else -tmp;
    }
}

/// The secular function's data: the deflated `col0` and `diag` entries.
fn Secular(comptime F: type) type {
    return struct {
        col0_perm: []const F,
        diag_perm: []const F,

        const Self = @This();

        /// `secular_eq`: `1 + Σ c² / ((d - shift - mu) (d + shift + mu))`.
        fn eval(e: Self, shift: F, mu: F) F {
            var res: F = 1;
            for (e.col0_perm, e.diag_perm) |c, d| res += (c / ((d - shift) - mu)) * (c / ((d + shift) + mu));
            return res;
        }

        /// `batch_secular_eq`: four at once, its own way.
        fn batch(e: Self, shift: [4]F, mu: [4]F) [4]F {
            var res: [4]F = @splat(1);
            for (e.col0_perm, e.diag_perm) |c, d| {
                for (&res, mu, shift) |*r, m, s| r.* += (c * c) / (((d - s) - m) * ((d + s) + m));
            }
            return res;
        }
    };
}

/// `compute_singular_values`: each root as `shift + mu`, `shift` the nearer
/// pole, for accuracy.
fn singularValues(io: Io, comptime F: type, shifts: []F, mus: []F, s: []F, diag: []const F, col0: []const F, sec: Secular(F)) void {
    const n = diag.len;
    var actual_n = n;
    while (actual_n > 1 and col0[actual_n - 1] == 0) actual_n -= 1;
    for (0..n) |k| {
        s[k] = 0;
        shifts[k] = 0;
        mus[k] = 0;
        if (col0[k] == 0 or actual_n == 1) {
            s[k] = if (k == 0) col0[0] else diag[k];
            shifts[k] = s[k];
            mus[k] = 0;
            continue;
        }
        const last_k = k == actual_n - 1;
        const left = diag[k];
        const right = if (last_k) diag[actual_n - 1] + ops.nrm2(io, F, n, @constCast(col0.ptr), 1) else blk: {
            var l = k + 1;
            while (col0[l] == 0) l += 1;
            break :blk diag[l];
        };
        const r = rootFinder(F, sec, left, right, last_k);
        s[k] = r[0] + r[1];
        shifts[k] = r[0];
        mus[k] = r[1];
    }
}

/// `secular_eq_root_finder`: the root in `(left, right)` of the increasing
/// secular function, as `.{ shift, mu }`: bracketing and bisection at
/// geometric midpoints, then the secant method, then bisection if it leaves
/// the bracket.
fn rootFinder(comptime F: type, sec: Secular(F), left: F, right: F, last: bool) [2]F {
    const two: F = 2;
    const eight: F = 8;
    const half: F = 0.5;
    const epsilon = std.math.floatEps(F);
    const inf = std.math.inf(F);
    const mid = left + (right - left) * half;
    const f0 = sec.batch(.{ 0, left, left, right }, .{ mid, if (last) right - left else (right - left) * half, (right - left) * half, (left - right) * half });
    var f_mid = f0[0];
    const f_max = f0[1];
    const f_mid_left_shift = f0[2];
    const f_mid_right_shift = f0[3];
    var shift: F, const mu0: F = if (last or f_mid > 0) .{ left, (right - left) * half } else .{ right, (left - right) * half };
    if (f_mid_left_shift <= 0 and f_mid_right_shift > 0) return .{ shift, mu0 };
    if (!last) {
        if (shift == left) {
            if (f_mid_left_shift < 0) {
                shift = right;
                f_mid = f_mid_right_shift;
            }
        } else if (f_mid_right_shift > 0) {
            shift = left;
            f_mid = f_mid_left_shift;
        }
    }
    var left_shifted: F = undefined;
    var f_left: F = undefined;
    var right_shifted: F = undefined;
    var f_right: F = undefined;
    if (shift == left) {
        left_shifted = 0;
        f_left = -inf;
        right_shifted = if (last) right - left else (right - left) * half;
        f_right = if (last) f_max else f_mid;
    } else {
        left_shifted = (left - right) * half;
        f_left = f_mid;
        right_shifted = 0;
        f_right = inf;
    }
    var iteration_count: usize = 0;
    var f_prev = f_mid;
    const half0 = half;
    const half1 = half0 * half0;
    const half2 = half1 * half1;
    const half3 = half2 * half2;
    const base = if (shift == left) right_shifted else left_shifted;
    const mu_values: [4]F = .{ base * half3, base * half2, base * half1, base * half0 };
    const f_values = sec.batch(@splat(shift), mu_values);
    if (shift == left) {
        var i: usize = 0;
        for (mu_values, f_values, 0..) |m, f, idx| if (f < 0) {
            left_shifted = m;
            f_left = f;
            i = idx + 1;
        };
        if (i < 4) {
            right_shifted = mu_values[i];
            f_right = f_values[i];
        }
    } else {
        var i: usize = 0;
        for (mu_values, f_values, 0..) |m, f, idx| if (f > 0) {
            right_shifted = m;
            f_right = f;
            i = idx + 1;
        };
        if (i < 4) {
            left_shifted = mu_values[i];
            f_left = f_values[i];
        }
    }
    while (right_shifted - left_shifted > two * epsilon * @max(@abs(left_shifted), @abs(right_shifted))) {
        const mid_arith = (left_shifted + right_shifted) * half;
        var mid_geo = @sqrt(@abs(left_shifted)) * @sqrt(@abs(right_shifted));
        if (left_shifted < 0) mid_geo = -mid_geo;
        const mid_shifted = if (mid_geo == 0) mid_arith else mid_geo;
        const fm = sec.eval(shift, mid_shifted);
        if (fm == 0) {
            return .{ shift, mid_shifted };
        } else if (fm > 0) {
            right_shifted = mid_shifted;
            f_prev = f_right;
            f_right = fm;
        } else {
            left_shifted = mid_shifted;
            f_prev = f_left;
            f_left = fm;
        }
        if (iteration_count == 4) break;
        iteration_count += 1;
    }
    const args: [4]F = if (left_shifted == 0)
        .{ right_shifted * two, right_shifted, f_prev, f_right }
    else if (right_shifted == 0)
        .{ left_shifted * two, left_shifted, f_prev, f_left }
    else
        .{ left_shifted, right_shifted, f_left, f_right };
    const sc = secant(F, sec, shift, left, right, args[0], args[1], args[2], args[3]);
    var mu_cur = sc.mu_cur;
    if (sc.left) |lc| if (sc.right) |rc| if (lc < rc) {
        if (lc > left_shifted) left_shifted = lc;
        if (rc < right_shifted) right_shifted = rc;
    };
    if (sc.use_bisection) {
        while (right_shifted - left_shifted > two * epsilon * @max(@abs(left_shifted), @abs(right_shifted))) {
            const m = (left_shifted + right_shifted) * half;
            const fm = sec.eval(shift, m);
            if (fm == 0) {
                break;
            } else if (fm > 0) {
                right_shifted = m;
            } else {
                left_shifted = m;
            }
        }
        mu_cur = (left_shifted + right_shifted) * half;
    }
    _ = eight;
    return .{ shift, mu_cur };
}

fn SecantResult(comptime F: type) type {
    return struct { use_bisection: bool, mu_cur: F, left: ?F, right: ?F };
}

/// `rootFinder`'s `secant` closure.
fn secant(comptime F: type, sec: Secular(F), shift: F, left: F, right: F, mu_cur_in: F, mu_prev_in: F, f_cur_in: F, f_prev_in: F) SecantResult(F) {
    const epsilon = std.math.floatEps(F);
    var mu_cur = mu_cur_in;
    var mu_prev = mu_prev_in;
    var f_cur = f_cur_in;
    var f_prev = f_prev_in;
    if (@abs(f_prev) < @abs(f_cur)) {
        std.mem.swap(F, &f_prev, &f_cur);
        std.mem.swap(F, &mu_prev, &mu_cur);
    }
    var left_candidate: ?F = null;
    var right_candidate: ?F = null;
    var use_bisection = false;
    const same_sign = (f_prev > 0) == (f_cur > 0);
    if (!same_sign) {
        left_candidate = @min(mu_cur, mu_prev);
        right_candidate = @max(mu_cur, mu_prev);
        if (mu_cur < mu_prev) {
            left_candidate = mu_cur;
            right_candidate = mu_prev;
        } else {
            left_candidate = mu_prev;
            right_candidate = mu_cur;
        }
    }
    while (f_cur != 0 and @abs(mu_cur - mu_prev) > 8 * epsilon * @max(@abs(mu_cur), @abs(mu_prev)) and @abs(f_cur - f_prev) > epsilon and !use_bisection) {
        const a = (f_cur - f_prev) * (mu_prev * mu_cur) / (mu_prev - mu_cur);
        const b = f_cur - a / mu_cur;
        const mu_zero = -a / b;
        const f_zero = sec.eval(shift, mu_zero);
        if (f_zero < 0) {
            left_candidate = mu_zero;
        } else {
            right_candidate = mu_zero;
        }
        mu_prev = mu_cur;
        f_prev = f_cur;
        mu_cur = mu_zero;
        f_cur = f_zero;
        if (shift == left and (mu_cur < 0 or mu_cur > right - left)) use_bisection = true;
        if (shift == right and (mu_cur > 0 or mu_cur < left - right)) use_bisection = true;
        if (@abs(f_cur) > @abs(f_prev)) {
            var k: F = 1;
            for (0..4) |_| {
                const mu_opposite = -a / (k * f_zero + b);
                const f_opposite = sec.eval(shift, mu_opposite);
                if (f_zero < 0 and f_opposite >= 0) {
                    right_candidate = mu_opposite;
                    break;
                }
                if (f_zero > 0 and f_opposite <= 0) {
                    left_candidate = mu_opposite;
                    break;
                }
                k *= 2;
            }
            use_bisection = true;
        }
    }
    // Faer stops here with whatever the secant reached when `f` stops changing
    // by more than `eps`, which happens far from a root that lies within
    // rounding of a pole (a tiny `col0` entry): then `mu` is wrong by orders
    // of magnitude, and `perturbCol0` with it. Bisect instead (issues.md).
    if (!use_bisection and f_cur != 0 and !(@abs(f_cur - f_prev) > epsilon) and
        @abs(mu_cur - mu_prev) > 8 * epsilon * @max(@abs(mu_cur), @abs(mu_prev))) use_bisection = true;
    return .{ .use_bisection = use_bisection, .mu_cur = mu_cur, .left = left_candidate, .right = right_candidate };
}

/// `deflate`: Gu and Eisenstat's deflations; returns how many rotations
/// each kind recorded (`jacobi_0i` of `deflation_43`, then `jacobi_ij` of
/// `deflation_44`), and `perm` the order `svdOfM` reads.
fn deflate(al: Allocator, comptime F: type, diag: []F, col0: []F, jacobi_coeff: []hh.Rotation(F), jacobi_idx: []usize, transpositions: []usize, perm: []usize, k: usize) Allocator.Error![2]usize {
    const n = diag.len;
    var jacobi_0i: usize = 0;
    var jacobi_ij: usize = 0;
    const max_diag = normMax(F, diag[1..]);
    const max_col0 = normMax(F, col0);
    const mx = @max(max_diag, max_col0);
    const eps = std.math.floatEps(F);
    const sml = std.math.floatMin(F);
    const eps_strict = @max(eps * max_diag, sml);
    const eps_coarse = 8 * eps * mx;
    if (diag[0] < eps_coarse) {
        diag[0] = eps_coarse;
        col0[0] = eps_coarse;
    }
    for (1..n) |i| if (@abs(col0[i]) < eps_strict) {
        col0[i] = 0;
    };
    for (1..n) |i| if (diag[i] < eps_coarse) {
        if (deflation43(F, diag, col0, i)) |rot| {
            jacobi_coeff[jacobi_0i] = rot;
            jacobi_idx[jacobi_0i] = i;
            jacobi_0i += 1;
        }
    };
    var total_deflation = true;
    for (1..n) |i| if (!(@abs(col0[i]) < sml)) {
        total_deflation = false;
        break;
    };
    perm[0] = 0;
    var p: usize = 1;
    for (1..n) |i| if (@abs(diag[i]) < sml) {
        perm[p] = i;
        p += 1;
    };
    {
        var i: usize = 1;
        var j: usize = k + 1;
        for (perm[p..]) |*slot| {
            if (i >= k + 1) {
                slot.* = j;
                j += 1;
            } else if (j >= n) {
                slot.* = i;
                i += 1;
            } else if (diag[i] < diag[j]) {
                slot.* = j;
                j += 1;
            } else {
                slot.* = i;
                i += 1;
            }
        }
    }
    if (total_deflation) {
        for (1..n) |i| {
            const pi = perm[i];
            if (@abs(diag[pi]) < sml or diag[pi] > diag[0]) {
                perm[i - 1] = perm[i];
            } else {
                perm[i - 1] = 0;
                break;
            }
        }
    }
    const real = try al.alloc(usize, 2 * n);
    defer al.free(real);
    const real_ind = real[0..n];
    const real_col = real[n..];
    for (real_ind, real_col, 0..) |*a, *b, i| {
        a.* = i;
        b.* = i;
    }
    for ((if (total_deflation) @as(usize, 0) else 1)..n) |i| {
        const pi = perm[n - if (total_deflation) i + 1 else i];
        const j = real_col[pi];
        std.mem.swap(F, &diag[i], &diag[j]);
        if (i != 0 and j != 0) std.mem.swap(F, &col0[i], &col0[j]);
        transpositions[i] = j;
        const real_i = real_ind[i];
        real_col[real_i] = j;
        real_col[pi] = i;
        real_ind[j] = real_i;
        real_ind[i] = pi;
    }
    col0[0] = diag[0];
    for (perm, 0..) |*x, i| x.* = i;
    for (transpositions, 0..) |j, i| std.mem.swap(usize, &perm[i], &perm[j]);
    var i = n - 1;
    while (i > 0 and (@abs(diag[i]) < sml or @abs(col0[i]) < sml)) i -= 1;
    while (i > 1) {
        const ib = i - 1;
        if (diag[i] - diag[ib] < eps_strict) {
            if (deflation44(F, diag, col0, ib, i)) |rot| {
                jacobi_coeff[jacobi_0i + jacobi_ij] = rot;
                jacobi_idx[jacobi_0i + jacobi_ij] = i;
                jacobi_ij += 1;
            }
        }
        i = ib;
    }
    return .{ jacobi_0i, jacobi_ij };
}

/// `deflation_43`: a tiny `diag[i]` folds `col0[i]` into `col0[0]`.
fn deflation43(comptime F: type, diag: []F, col0: []F, i: usize) ?hh.Rotation(F) {
    const p = col0[0];
    const q = col0[i];
    if (p == 0 and q == 0) return null;
    const rot = hh.Rotation(F).make(p, q);
    const r = rot.c * p - rot.s * q;
    col0[0] = r;
    diag[0] = r;
    col0[i] = 0;
    return rot;
}

/// `deflation_44`: equal `diag[i]` and `diag[j]` fold `col0[j]` into
/// `col0[i]`.
fn deflation44(comptime F: type, diag: []F, col0: []F, i: usize, j: usize) ?hh.Rotation(F) {
    const p = col0[i];
    const q = col0[j];
    if (p == 0 and q == 0) {
        diag[i] = diag[j];
        return null;
    }
    const rot = hh.Rotation(F).make(p, q);
    const r = rot.c * p - rot.s * q;
    col0[i] = r;
    col0[j] = 0;
    diag[i] = diag[j];
    return rot;
}

// ---- tests -----------------------------------------------------------------------

const testing = std.testing;
const check = @import("../../blas/_test/blas1/check.zig");

/// `U [S; 0] Vᵀ` against the `(n + 1) x n` (or with `square`, `n x n`)
/// lower bidiagonal matrix, within faer's `test_divide_and_conquer` (and
/// `test_qr_algorithm`) tolerance: `128 eps` times `max(|d|, |e|) sqrt n`,
/// times 10 for divide and conquer; and `Uᵀ U = I`, `Vᵀ V = I` within it.
fn expectBidiagSvd(comptime F: type, diag: []const F, sub: []const F, d: []const F, u: Mat(F), v: Mat(F), square: bool, factor: F) !void {
    const n = diag.len;
    const rows = if (square) n else n + 1;
    const scale = @max(normMax(F, diag), normMax(F, sub));
    const tol = 128 * std.math.floatEps(F) * @max(scale, 1) * @sqrt(@as(F, @floatFromInt(n))) * factor;
    for (0..n) |j| for (0..rows) |i| {
        var x: F = 0;
        for (0..n) |kk| x += u.get(i, kk) * d[kk] * v.get(j, kk);
        const want: F = if (i == j) diag[j] else if (i == j + 1) sub[j] else 0;
        if (!(@abs(x - want) <= tol)) {
            check.say("bidiag svd n={d}: ({d},{d}) = {e}, want {e}\n", .{ n, i, j, x, want });
            return error.TestUnexpectedResult;
        }
    };
    for (0..n) |j| for (0..n) |i| {
        var x: F = 0;
        for (0..n) |kk| x += v.get(kk, i) * v.get(kk, j);
        try testing.expect(@abs(x - @as(F, if (i == j) 1 else 0)) <= tol);
    };
    for (0..rows) |j| for (0..rows) |i| {
        var x: F = 0;
        for (0..rows) |kk| x += u.get(kk, i) * u.get(kk, j);
        try testing.expect(@abs(x - @as(F, if (i == j) 1 else 0)) <= tol);
    };
    for (1..n) |i| try testing.expect(d[i - 1] >= d[i]);
    for (d) |x| try testing.expect(x >= 0);
}

/// A bidiagonal with entries spread over several magnitudes, and clusters.
fn testBidiag(comptime F: type, diag: []F, sub: []F, seed: u64) void {
    var prng = std.Random.DefaultPrng.init(seed);
    const r = prng.random();
    for (diag, sub, 0..) |*d, *e, i| {
        d.* = (r.float(F) * 2 - 1) * std.math.pow(F, 10, @floatFromInt(@as(i32, @intCast(i % 5)) - 2));
        e.* = (r.float(F) * 2 - 1);
        if (i % 7 == 3) d.* = diag[i - 1];
    }
}

test "faer bidiag svd: qr_algorithm and divide_and_conquer reconstruct B, both thresholds" {
    const al = testing.allocator;
    inline for (.{ f64, f32 }) |F| {
        for ([_]usize{ 1, 2, 3, 7, 20, 64, 129 }) |n| {
            const diag = try al.alloc(F, n);
            defer al.free(diag);
            const sub = try al.alloc(F, n);
            defer al.free(sub);
            testBidiag(F, diag, sub, n);
            const d = try al.alloc(F, n);
            defer al.free(d);
            const e = try al.alloc(F, n);
            defer al.free(e);
            // qr_algorithm on the square part
            {
                @memcpy(d, diag);
                @memcpy(e, sub);
                e[n - 1] = 0;
                const sq = try al.dupe(F, e);
                defer al.free(sq);
                const u = try ops.alloc(F, al, n, n);
                defer ops.free(F, al, u);
                const v = try ops.alloc(F, al, n, n);
                defer ops.free(F, al, v);
                u.identity();
                v.identity();
                try qrAlgorithm(F, d, e, u, v);
                for (e) |x| try testing.expectEqual(@as(F, 0), x);
                try expectBidiagSvd(F, diag, sq, d, u, v, true, 1);
            }
            // divide and conquer, faer's threshold of 4 and 8, and 128
            for ([_]usize{ 4, 8, 128 }) |thr| {
                @memcpy(d, diag);
                @memcpy(e, sub);
                const u = try ops.alloc(F, al, n + 1, n + 1);
                defer ops.free(F, al, u);
                const v = try ops.alloc(F, al, n, n);
                defer ops.free(F, al, v);
                try divideAndConquer(al, check.io, F, d, e, .{ .full = u }, v, thr);
                expectBidiagSvd(F, diag, sub, d, u, v, false, 10) catch |err| {
                    check.say("{s} n={d} threshold {d}\n", .{ @typeName(F), n, thr });
                    return err;
                };
                // the values alone, and without U, agree
                const d2 = try al.dupe(F, diag);
                defer al.free(d2);
                const e2 = try al.dupe(F, sub);
                defer al.free(e2);
                const u_bot = try ops.alloc(F, al, 2, n + 1);
                defer ops.free(F, al, u_bot);
                try divideAndConquer(al, check.io, F, d2, e2, .{ .two_rows_storage = u_bot }, null, thr);
                for (d, d2) |x, y| try testing.expect(@abs(x - y) <= 1024 * std.math.floatEps(F) * @max(d[0], 1));
            }
        }
    }
}

test "faer bidiag svd: a zero block deflates, with or without U; zink's bidiagonal keeps its smallest value" {
    const al = testing.allocator;
    // the zero block is the second half of the first half, so its `TwoRows`
    // last row (faer: 0, the identity: eₙ) becomes the coupling `lambda`
    const n = 24;
    var diag: [n]f64 = undefined;
    var sub: [n]f64 = undefined;
    for (0..n) |i| {
        const zero = i >= 7 and i < 12;
        diag[i] = if (zero) 0 else 1 + @as(f64, @floatFromInt(i));
        sub[i] = if (zero) 0 else 0.5;
    }
    var d = diag;
    var e = sub;
    const u = try ops.alloc(f64, al, n + 1, n + 1);
    defer ops.free(f64, al, u);
    const v = try ops.alloc(f64, al, n, n);
    defer ops.free(f64, al, v);
    try divideAndConquer(al, check.io, f64, &d, &e, .{ .full = u }, v, 4);
    try expectBidiagSvd(f64, &diag, &sub, &d, u, v, false, 10);
    var d2 = diag;
    var e2 = sub;
    const u_bot = try ops.alloc(f64, al, 2, n + 1);
    defer ops.free(f64, al, u_bot);
    try divideAndConquer(al, check.io, f64, &d2, &e2, .{ .two_rows_storage = u_bot }, null, 4);
    // within a few ulps; faer's 0 is off by 1.5e-13 at d[5]
    for (d, d2) |x, y| try testing.expectApproxEqAbs(x, y, 8 * std.math.floatEps(f64) * d[0]);
    // faer's `test_zink` (svd/mod.rs)
    var zd = [_]f64{ -9.931833701529301, -10.920807536026027, -52.33647796311243, 2.3685025127736967, 2.421701994236093, -0.5051763005624579, -0.04808263896606017, -0.003875251886338955, -0.0006413264967716465, -0.003381944152463707, 2.981152313236375e-5, 5.4290648208388795e-6, -6.329275972084404e-7, -6.879142344209158e-7, -5.265228263479126e-9, -2.941999902335516e-9, -1.3060984997930294e-10, 7.07516117218088e-12, 1.8657003929029376e-12, -6.216080089659131e-14 };
    var zs = [_]f64{ -57.8029649868477, 17.67263066467847, 8.884153814270894, -9.01998231080713, -1.028638150814966, 0.22247719217200435, 0.016389886745811315, -0.004090989452162578, 0.00036818904090536926, -0.0031394146217732367, -7.571300829706796e-6, 3.0045718718618155e-6, 2.1329796886727743e-6, 9.259701025627789e-8, 2.2291214755992877e-9, -2.3017207713252894e-9, 6.807967994979358e-11, 2.1677299575405587e-12, -3.07282771050034e-13, 0.0 };
    const zu = try ops.alloc(f64, al, 2, zd.len + 1);
    defer ops.free(f64, al, zu);
    try divideAndConquer(al, check.io, f64, &zd, &zs, .{ .two_rows_storage = zu }, null, 8);
    try testing.expect(zd[zd.len - 1] != 0);
}

test "faer bidiag svd: deflation 4.3 and 4.4, faer's cases" {
    // `test_deflation43`, `test_deflation44`, `test_both_deflation`
    const al = testing.allocator;
    const cases = [_]struct { diag: [7]f32, want: [2]usize }{
        .{ .diag = .{ 1, 5, 3, 1e-7, 4, 2, 2e-7 }, .want = .{ 2, 0 } },
        .{ .diag = .{ 1, 5, 3, 1, 4, 2, 1 }, .want = .{ 0, 1 } },
        .{ .diag = .{ 1, 5, 3, 2, 4, 2, 0 }, .want = .{ 1, 1 } },
    };
    for (cases) |c| {
        var diag = c.diag;
        var col: [7]f32 = @splat(1);
        var coeff: [7]hh.Rotation(f32) = undefined;
        var idx: [7]usize = undefined;
        var perm: [7]usize = undefined;
        var tr: [7]usize = .{ 0, 1, 2, 3, 4, 5, 6 };
        const got = try deflate(al, f32, &diag, &col, &coeff, &idx, &tr, &perm, 3);
        try testing.expectEqual(c.want, got);
        // M' in the permuted places, then the rotations undone, is M: the
        // diagonal and a first column of ones
        var m: [7][7]f32 = @splat(@splat(0));
        for (0..7) |i| {
            m[perm[i]][perm[i]] = diag[i];
            m[perm[i]][perm[0]] = col[i];
        }
        var t = got[1];
        while (t > 0) {
            t -= 1;
            const rot = coeff[got[0] + t];
            const pi = perm[idx[got[0] + t] - 1];
            const pj = perm[idx[got[0] + t]];
            rot.left(7, &m[pj], 1, &m[pi], 1);
            rot.left(7, @as([*]f32, @ptrCast(&m)) + pj, 7, @as([*]f32, @ptrCast(&m)) + pi, 7);
        }
        t = got[0];
        while (t > 0) {
            t -= 1;
            coeff[t].left(7, &m[idx[t]], 1, &m[0], 1);
        }
        for (0..7) |i| for (0..7) |j| {
            const want: f32 = if (j == 0) 1 else if (i == j) c.diag[i] else 0;
            try testing.expect(@abs(m[i][j] - want) <= 128 * std.math.floatEps(f32) * @max(1, @abs(want)));
        };
    }
}
