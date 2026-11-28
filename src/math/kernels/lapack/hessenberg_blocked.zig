//! The blocked Hessenberg reduction, from faer v0.24.4
//! (<https://codeberg.org/sarah-quinones/faer>, tag `faer-v0.24.4`, MIT):
//! `faer/src/linalg/evd/hessenberg.rs`'s `hessenberg_gqvdg_blocked` and
//! `hessenberg_gqvdg_unblocked` (G. Quintana-Ortí and R. van de Geijn's
//! algorithm, "Improving the performance of reduction to Hessenberg form",
//! ACM TOMS 32(2), 2006): each block's reflections and their block factor
//! by the unblocked pass, then the rest of `A` by `gemm` and `trsm`. The
//! factors are faer's (`H = I - V T⁻¹ Vᴴ`, `T` upper, `tau` on its
//! diagonal), so `hessenberg.zig`'s callers read them as the unblocked
//! form's. A strictly triangular operand is copied with a zero diagonal.

const std = @import("std");
const ops = @import("ops.zig");
const blas_calls = @import("blas_calls.zig");
const hh = @import("householder.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = ops.Mat;
const Error = ops.Error;
const one = ops.one;
const zero = ops.zero;

/// `hessenberg_gqvdg_blocked`: `a` is `n x n`, `h` `b x (n - 1)`, its
/// rows the block size.
pub fn blocked(al: Allocator, io: Io, comptime T: type, a: Mat(T), h: Mat(T)) Error!void {
    const n = a.rows;
    const b = h.rows;
    const z = try ops.alloc(T, al, n, b);
    defer ops.free(T, al, z);
    const beta = try al.alloc(T, b);
    defer al.free(beta);
    // `X0` (`j x bs_u`), `X2`'s transpose (`bs_u x rest`), strict triangles of `U1`
    const x0 = try ops.alloc(T, al, n, b);
    defer ops.free(T, al, x0);
    const xt = try ops.alloc(T, al, b, n);
    defer ops.free(T, al, xt);
    const l1 = try ops.alloc(T, al, b, b);
    defer ops.free(T, al, l1);
    var j: usize = 0;
    while (j < n) {
        const bs = @min(b, n - j);
        const bu = @min(bs, n - j - 1);
        const t11 = h.sub(0, j, bu, bu);
        try unblocked(al, io, T, a.sub(j, j, n - j, n - j), z.sub(j, 0, n - j, bs), t11, beta);
        const rest = n - j - bu;
        const a01 = a.sub(0, j, j, bu);
        const a02 = a.sub(0, j + bu, j, rest);
        const a12 = a.sub(j, j + bu, bu, rest);
        const a22 = a.sub(j + bu, j + bu, rest, rest);
        const v1 = a.sub(j, j, bu, bu);
        const v2 = a.sub(j + bu, j, rest, bu);
        const z1 = z.sub(j, 0, bu, bu);
        const z2 = z.sub(j + bu, 0, rest, bu);
        // `U1`'s strictly lower triangle
        const l = l1.sub(0, 0, bu, bu);
        for (0..bu) |c| for (0..bu) |r| {
            l.at(r, c).* = if (r > c) v1.get(r, c) else zero(T);
        };
        // the rows above: `A0 := A0 (I - U T⁻¹ Uᴴ)`
        const x = x0.sub(0, 0, j, bu);
        try blas_calls.gemm(al, io, T, .n, .n, one(T), a01, l, zero(T), x);
        try blas_calls.gemm(al, io, T, .n, .n, one(T), a02, v2, one(T), x);
        try blas_calls.triRight(al, io, T, true, true, .n, false, one(T), t11, x);
        try blas_calls.gemm(al, io, T, .n, .c, ops.neg(T, one(T)), x, l, one(T), a01);
        try blas_calls.gemm(al, io, T, .n, .c, ops.neg(T, one(T)), x, v2, one(T), a02);
        // the right: `A := A - Z T⁻¹ Uᴴ`, then from the left `(I - U T⁻ᴴ Uᴴ) A`
        try blas_calls.triRight(al, io, T, true, true, .n, false, one(T), t11, z1);
        try blas_calls.triRight(al, io, T, true, true, .n, false, one(T), t11, z2);
        try blas_calls.gemm(al, io, T, .n, .c, ops.neg(T, one(T)), z1, v2, one(T), a12);
        try blas_calls.gemm(al, io, T, .n, .c, ops.neg(T, one(T)), z2, v2, one(T), a22);
        const y = xt.sub(0, 0, bu, rest);
        try blas_calls.gemm(al, io, T, .c, .n, one(T), l, a12, zero(T), y);
        try blas_calls.gemm(al, io, T, .c, .n, one(T), v2, a22, one(T), y);
        try blas_calls.tri(al, io, T, true, true, .c, false, one(T), t11, y);
        try blas_calls.gemm(al, io, T, .n, .n, ops.neg(T, one(T)), l, y, one(T), a12);
        try blas_calls.gemm(al, io, T, .n, .n, ops.neg(T, one(T)), v2, y, one(T), a22);
        for (0..bs) |k| if (j + k + 1 < n) {
            a.at(j + k + 1, j + k).* = beta[k];
        };
        j += bs;
    }
}

/// `hessenberg_gqvdg_unblocked`: the first `h.rows` reflections of the
/// `n x n` `a` (`n` being what is left), `z` (`n x h.rows`) their `A V`
/// columns, `h` their block factor, `beta` the subdiagonal they leave.
fn unblocked(al: Allocator, io: Io, comptime T: type, a: Mat(T), z: Mat(T), h: Mat(T), beta: []T) Error!void {
    const n = a.rows;
    const b = h.rows;
    const x = try al.alloc(T, n);
    defer al.free(x);
    for (0..b) |k| {
        const x0 = x[0..k];
        const t00 = h.sub(0, 0, k, k);
        const a1 = a.col(k);
        const v0 = a.sub(0, 0, n, k);
        const v20 = a.sub(k + 1, 0, n - k - 1, k);
        // `x0 := T00⁻¹ U10ᴴ`, then `A1 -= Z0 x0`
        for (x0, 0..) |*e, i| e.* = ops.conj(T, v0.get(k, i));
        if (k > 0) {
            try blas_calls.tri(al, io, T, true, true, .n, false, one(T), t00, column(T, x0));
            try blas_calls.gemv(al, io, T, .n, ops.neg(T, one(T)), z.sub(0, 0, n, k), x0.ptr, 1, one(T), a1, 1);
        }
        const a01 = a1[0..k];
        const a11 = &a1[k];
        const a21 = a1[k + 1 .. n];
        if (k > 0) {
            // `x0 := strict(U00ᴴ) A01 + U10ᴴ A11 + U20ᴴ A21`
            for (x0, 0..) |*e, i| {
                var s = zero(T);
                for (i + 1..k) |r| s = ops.add(T, s, ops.mul(T, ops.conj(T, v0.get(r, i)), a01[r]));
                e.* = ops.add(T, s, ops.mul(T, a11.*, ops.conj(T, v0.get(k, i))));
            }
            if (a21.len > 0) try blas_calls.gemv(al, io, T, .c, one(T), v20, a21.ptr, 1, one(T), x0.ptr, 1);
            try blas_calls.tri(al, io, T, true, true, .c, false, one(T), t00, column(T, x0));
            // `A01 -= strict(U00) x0`, `A11 -= U10 x0`, `A21 -= U20 x0`
            var r = k;
            while (r > 0) {
                r -= 1;
                var s = zero(T);
                for (0..r) |i| s = ops.add(T, s, ops.mul(T, v0.get(r, i), x0[i]));
                a01[r] = ops.sub(T, a01[r], s);
            }
            var s = zero(T);
            for (x0, 0..) |e, i| s = ops.add(T, s, ops.mul(T, v0.get(k, i), e));
            a11.* = ops.sub(T, a11.*, s);
            if (a21.len > 0) try blas_calls.gemv(al, io, T, .n, ops.neg(T, one(T)), v20, x0.ptr, 1, one(T), a21.ptr, 1);
        }
        if (k + 1 < n) {
            const info = hh.make(io, T, &a21[0], a21.len - 1, a21.ptr + 1, 1, null, 0);
            beta[k] = a21[0];
            a21[0] = one(T);
            h.at(k, k).* = ops.scalar(T, info.tau);
        } else {
            h.at(k, k).* = ops.scalar(T, std.math.inf(ops.Float(T)));
        }
        // `Z1 := A2 A21`, `T01 := U20ᴴ A21`
        const zk = z.col(k);
        try blas_calls.gemv(al, io, T, .n, one(T), a.sub(0, k + 1, n, n - k - 1), a21.ptr, 1, zero(T), zk, 1);
        if (k > 0) try blas_calls.gemv(al, io, T, .c, one(T), v20, a21.ptr, 1, zero(T), h.col(k), 1);
    }
}

/// `v` as a one-column matrix.
fn column(comptime T: type, v: []T) Mat(T) {
    return .{ .ptr = v.ptr, .rows = v.len, .cols = 1, .ld = @max(v.len, 1) };
}
