//! The SVD of a real lower bidiagonal matrix, from faer v0.24.4
//! (<https://codeberg.org/sarah-quinones/faer>, tag `faer-v0.24.4`, MIT):
//! `faer/src/linalg/svd/bidiag_svd.rs`'s divide and conquer, its secular
//! solver (`bidiag_secular.zig`) and deflation (`bidiag_deflate.zig`), with
//! Reference-LAPACK's `xBDSQR` as `qrAlgorithm` (`bdsqr.zig`, T126).

const std = @import("std");
const mem = @import("../../mem.zig");
const ops = @import("ops.zig");
const blas_calls = @import("blas_calls.zig");
const hh = @import("householder.zig");
const bdsqr = @import("bdsqr.zig");
const bidiag_secular = @import("bidiag_secular.zig");
const bidiag_deflate = @import("bidiag_deflate.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = ops.Mat;

pub const qrAlgorithm = bdsqr.qrAlgorithm;
const Secular = bidiag_secular.Secular;
const singularValues = bidiag_secular.singularValues;
pub const rootFinder = bidiag_secular.rootFinder;
const deflate = bidiag_deflate.deflate;

/// What the solvers return: `NoConvergence` once `qrAlgorithm` has run its
/// iterations out.
pub const Error = ops.Error || error{NoConvergence};

/// `nbits`: the significand's bits.
fn nbits(comptime F: type) usize {
    return std.math.floatMantissaBits(F) + 1;
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
    const max = @max(ops.absMax(F, diag), ops.absMax(F, subdiag));
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
    ops.scaleDown(F, diag, max);
    ops.scaleDown(F, subdiag, max);
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
                    ops.swapCols(F, um, i, i + 1);
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
            try blas_calls.gemm(al, io, F, .n, .n, 1, ul, um, 0, cu);
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
    try blas_calls.gemm(al, io, F, .n, .n, 1, v.sub(0, 1, k, k), vm.sub(1, 0, k, n), 0, cv.sub(0, 0, k, n));
    try blas_calls.gemm(al, io, F, .n, .n, 1, v.sub(k + 1, k + 1, rem, rem), vm.sub(k + 1, 0, rem, n), 0, cv.sub(k + 1, 0, rem, n));
    try blas_calls.gemm(al, io, F, .n, .n, 1, v.sub(k, 0, 1, 1), vm.sub(0, 0, 1, n), 0, cv.sub(k, 0, 1, n));
    v.copyFrom(cv);
}

/// `update_u`: `U := U Um`, by its blocks.
fn updateU(al: Allocator, io: Io, comptime F: type, u: Mat(F), um: Mat(F), k: usize, rem: usize) ops.Error!void {
    const n = u.rows - 1;
    const cu = try ops.alloc(F, al, n + 1, n + 1);
    defer ops.free(F, al, cu);
    const c1 = cu.sub(0, 0, k + 1, n + 1);
    const c2 = cu.sub(k + 1, 0, rem + 1, n + 1);
    try blas_calls.gemm(al, io, F, .n, .n, 1, u.sub(0, 0, k + 1, k + 1), um.sub(0, 0, k + 1, n + 1), 0, c1);
    try blas_calls.gemm(al, io, F, .n, .n, 1, u.sub(0, n, k + 1, 1), um.sub(n, 0, 1, n + 1), 1, c1);
    try blas_calls.gemm(al, io, F, .n, .n, 1, u.sub(k + 1, k + 1, rem + 1, rem + 1), um.sub(k + 1, 0, rem + 1, n + 1), 0, c2);
    try blas_calls.gemm(al, io, F, .n, .n, 1, u.sub(k + 1, 0, rem + 1, 1), um.sub(0, 0, 1, n + 1), 1, c2);
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
    sortAscending(F, s[0..actual_n], col_perm[0..actual_n]);
    for (col_perm[0..actual_n], 0..) |p, i| col_perm_inv[p] = i;
    singularVectors(io, F, um, vm, zhat, diag, perm, outer_perm, col_perm_inv, actual_n, shifts, mus);
    for (0..actual_n) |idx| diag[idx] = s[actual_n - idx - 1];
    for (actual_n..n) |idx| diag[idx] = s[idx];
}

/// `s` ascending, `perm` moved with it. Faer makes one bubble pass,
/// assuming deflation leaves at most the last value out of place; heavy
/// deflation can leave more (T130, issue #254), so a stable insertion sort,
/// the same order wherever one pass was enough.
pub fn sortAscending(comptime F: type, s: []F, perm: []usize) void {
    for (1..s.len) |i| {
        var j = i;
        while (j > 0 and s[j - 1] > s[j]) : (j -= 1) {
            mem.swap(F, &s[j - 1], &s[j]);
            mem.swap(usize, &perm[j - 1], &perm[j]);
        }
    }
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
        const nu = 1 / blas_calls.nrm2(io, F, n, u, 1);
        for (u[0..n]) |*x| x.* *= nu;
        if (v) |vc| {
            vc[outer_perm[0]] = -1;
            for (perm[1..]) |i| vc[outer_perm[i]] = ((diag[i] * zhat[i]) / ((diag[i] - shift) - mu)) / (diag[i] + (shift + mu));
            const nv = 1 / blas_calls.nrm2(io, F, n, vc, 1);
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
