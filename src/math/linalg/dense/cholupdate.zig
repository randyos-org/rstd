//! A Cholesky factor's rank-1 update and downdate (MATLAB's `cholupdate`),
//! by faer v0.24.4's rank-r update with `r = 1`
//! (`faer/src/linalg/cholesky/llt/update.rs`, `RankRUpdate::run`, MIT), on
//! `L = Rᴴ`.

const std = @import("std");
const matrix = @import("../matrix.zig");
const vector = @import("../vector.zig");
const ops = @import("../../kernels/lapack/ops.zig");
const util = @import("util.zig");

const Allocator = std.mem.Allocator;
const Mat = matrix.Mat;
const Real = vector.Real;
const ElementOf = util.ElementOf;
const constOf = util.constOf;

/// `cholupdate`'s: `DimensionMismatch` unless `x` is a vector of `R`'s
/// order; `NotPositiveDefinite` for a diagonal of `R` that is not positive
/// or a downdate that leaves no Cholesky factor (MATLAB's `p > 0`).
pub const CholupdateError = error{ NotSquare, DimensionMismatch, NotPositiveDefinite } || Allocator.Error;

/// MATLAB's `'+'` and `'-'`.
pub const CholupdateOp = enum { update, downdate };

/// The upper `R1` with `R1ᴴ R1 = Rᴴ R + x xᴴ` (`.update`) or `Rᴴ R - x xᴴ`
/// (`.downdate`), from the upper `R` with a positive diagonal (MATLAB's
/// `cholupdate(R, x, '+' | '-')`); `R`'s lower triangle is not read. O(n²).
pub fn cholupdate(allocator: Allocator, r: anytype, x: anytype, op: CholupdateOp) CholupdateError!Mat(ElementOf(@TypeOf(r))) {
    const T = ElementOf(@TypeOf(r));
    const F = Real(T);
    const cr = constOf(r);
    if (cr.rows != cr.cols) return error.NotSquare;
    const n = cr.rows;
    const xs = util.strided(constOf(x)) catch return error.DimensionMismatch;
    if (xs.n != n) return error.DimensionMismatch;
    // `L = Rᴴ`, packed column-major
    const l = try allocator.alloc(T, n * n);
    defer allocator.free(l);
    for (0..n) |j| for (j..n) |i| {
        l[i + j * n] = ops.conj(T, cr.at(j, i));
    };
    const w = try allocator.alloc(T, n);
    defer allocator.free(w);
    for (w, 0..) |*v, i| v.* = xs.ptr[i * xs.inc];
    var alpha: F = if (op == .update) 1 else -1;
    for (0..n) |j| {
        const d = ops.re(T, l[j + j * n]);
        if (!(d > 0)) return error.NotPositiveDefinite;
        // the new pivot, and the column's rotation `p`, `beta`, `gamma`
        const alpha_conj_p = ops.mulReal(T, ops.conj(T, w[j]), alpha);
        const new_d2 = d * d + ops.re(T, ops.mul(T, alpha_conj_p, w[j]));
        if (!(new_d2 > 0)) return error.NotPositiveDefinite;
        const new_d = @sqrt(new_d2);
        const gamma = new_d / d;
        const beta = ops.mulReal(T, alpha_conj_p, 1 / new_d);
        const p = ops.neg(T, ops.mulReal(T, w[j], 1 / d));
        alpha -= ops.re(T, ops.mul(T, beta, ops.conj(T, beta)));
        l[j + j * n] = ops.scalar(T, new_d);
        for (j + 1..n) |i| {
            w[i] = ops.add(T, ops.mul(T, p, l[i + j * n]), w[i]);
            l[i + j * n] = ops.add(T, ops.mul(T, beta, w[i]), ops.mulReal(T, l[i + j * n], gamma));
        }
    }
    // `R1 = Lᴴ`, zeros below the diagonal
    const r1 = try Mat(T).zeros(allocator, n, n);
    for (0..n) |j| for (j..n) |i| {
        r1.data[j + i * n] = ops.conj(T, l[i + j * n]);
    };
    return r1;
}
