//! Deflation in `bidiag_svd.zig`'s divide and conquer, from faer v0.24.4
//! (`faer/src/linalg/svd/bidiag_svd.rs`, MIT).

const std = @import("std");
const mem = @import("../../../mem.zig");
const ops = @import("ops.zig");
const hh = @import("householder.zig");

const Allocator = std.mem.Allocator;

/// `deflate`: Gu and Eisenstat's deflations; returns how many rotations
/// each kind recorded (`jacobi_0i` of `deflation_43`, then `jacobi_ij` of
/// `deflation_44`), and `perm` the order `svdOfM` reads.
pub fn deflate(al: Allocator, comptime F: type, diag: []F, col0: []F, jacobi_coeff: []hh.Rotation(F), jacobi_idx: []usize, transpositions: []usize, perm: []usize, k: usize) Allocator.Error![2]usize {
    const n = diag.len;
    var jacobi_0i: usize = 0;
    var jacobi_ij: usize = 0;
    const max_diag = ops.absMax(F, diag[1..]);
    const max_col0 = ops.absMax(F, col0);
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
        mem.swap(F, &diag[i], &diag[j]);
        if (i != 0 and j != 0) mem.swap(F, &col0[i], &col0[j]);
        transpositions[i] = j;
        const real_i = real_ind[i];
        real_col[real_i] = j;
        real_col[pi] = i;
        real_ind[j] = real_i;
        real_ind[i] = pi;
    }
    col0[0] = diag[0];
    for (perm, 0..) |*x, i| x.* = i;
    for (transpositions, 0..) |j, i| mem.swap(usize, &perm[i], &perm[j]);
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
