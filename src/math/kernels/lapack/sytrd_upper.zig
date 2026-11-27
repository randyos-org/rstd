//! `?SYTRD` and `?HETRD` on the upper triangle, from the bottom right:
//! Reference-LAPACK as OpenBLAS ships it (`lapack-netlib/SRC/dsytrd.f`,
//! `zhetrd.f`, `dlatrd.f`, `zlatrd.f`, `dsytd2.f`, `zhetd2.f`, their `UPLO
//! = 'U'` branches; OpenBLAS v0.3.34-317; modified BSD), as `sytrd.zig`
//! takes the lower.

const std = @import("std");
const ops = @import("ops.zig");
const auxiliary = @import("auxiliary.zig");
const blas_calls = @import("blas_calls.zig");
const strided = @import("../blas/strided.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = ops.Mat;
const Float = ops.Float;
const Error = ops.Error;
const isComplex = ops.isComplex;
const Strided = strided.Strided;
const zero = ops.zero;
const one = ops.one;
const scalar = ops.scalar;
const re = ops.re;
const add = ops.add;
const sub = ops.sub;
const mul = ops.mul;
const conj = ops.conj;
const mulReal = ops.mulReal;

/// The upper triangle's reduction (`dsytrd.f`'s `UPPER` branch): blocks of
/// `nb` columns from the right through `latrd`, each followed by one
/// rank-2k update of the rest, the first `kk` columns through `sytd2`.
/// `al` is `sytrd.zig`'s kept allocator.
pub fn reduce(al: Allocator, io: Io, comptime T: type, a: Mat(T), d: []Float(T), e: []Float(T), tau: []T, nb: usize, nx: usize) Error!void {
    const n = a.rows;
    var kk = n;
    if (n > nx) {
        kk = n - ((n - nx + nb - 1) / nb) * nb;
        const w = try ops.alloc(T, al, n, nb);
        defer ops.free(T, al, w);
        const tmp = try al.alloc(T, nb);
        defer al.free(tmp);
        // `DO I = N - NB + 1, KK + 1, -NB`, `i` 0-based
        var i = n - nb;
        while (true) {
            try latrd(al, io, T, a.sub(0, 0, i + nb, i + nb), e, tau, w.sub(0, 0, i + nb, nb), nb, tmp);
            try blas_calls.rank2k(al, io, T, true, a.sub(0, i, i, nb), w.sub(0, 0, i, nb), a.sub(0, 0, i, i));
            for (i..i + nb) |j| {
                a.at(j - 1, j).* = scalar(T, e[j - 1]);
                d[j] = re(T, a.get(j, j));
            }
            if (i < kk + nb) break;
            i -= nb;
        }
    }
    try sytd2(al, io, T, a.sub(0, 0, kk, kk), d, e, tau);
}

/// `?LATRD` upper: reduces the last `nb` columns of the `n x n` `a`, leaving
/// in `w` (`n x nb`) what the rank-2k update of the rest takes: `A11 -= V
/// Wᴴ + W Vᴴ`. `tmp` holds `nb` elements.
fn latrd(al: Allocator, io: Io, comptime T: type, a: Mat(T), e: []Float(T), tau: []T, w: Mat(T), nb: usize, tmp: []T) Error!void {
    const n = a.rows;
    const cx = comptime isComplex(T);
    const minus = scalar(T, -1);
    // `DO I = N, N - NB + 1, -1`, `i` 0-based, `iw` its column of `w`
    var i = n;
    while (i > n - nb) {
        i -= 1;
        const iw = i - (n - nb);
        const rest = n - 1 - i;
        if (rest > 0) {
            if (cx) a.at(i, i).* = scalar(T, re(T, a.get(i, i)));
            // A(0:i+1, i) -= A(0:i+1, i+1:n) conj(W(i, iw+1:nb))ᵀ + W(0:i+1, iw+1:nb) conj(A(i, i+1:n))ᵀ
            for (tmp[0..rest], 0..) |*t, j| t.* = conj(T, w.get(i, iw + 1 + j));
            try blas_calls.gemv(al, io, T, .n, minus, a.sub(0, i + 1, i + 1, rest), tmp.ptr, 1, one(T), a.col(i), 1);
            for (tmp[0..rest], 0..) |*t, j| t.* = conj(T, a.get(i, i + 1 + j));
            try blas_calls.gemv(al, io, T, .n, minus, w.sub(0, iw + 1, i + 1, rest), tmp.ptr, 1, one(T), a.col(i), 1);
            if (cx) a.at(i, i).* = scalar(T, re(T, a.get(i, i)));
        }
        if (i == 0) continue;
        // the reflector taking A(0:i-1, i) to zero against A(i-1, i)
        const v = a.col(i);
        var alpha = v[i - 1];
        tau[i - 1] = try auxiliary.larfg(io, T, &alpha, Strided(T).of(v[0 .. i - 1]));
        e[i - 1] = re(T, alpha);
        v[i - 1] = one(T);
        // W(0:i, iw) := A11 v - V W(i+1:n, iw)... as `dlatrd.f`
        const wi = w.col(iw);
        try blas_calls.symv(al, io, T, true, one(T), a.sub(0, 0, i, i), v, wi);
        if (rest > 0) {
            const lower = wi + i + 1;
            try blas_calls.gemv(al, io, T, .c, one(T), w.sub(0, iw + 1, i, rest), v, 1, zero(T), lower, 1);
            try blas_calls.gemv(al, io, T, .n, minus, a.sub(0, i + 1, i, rest), lower, 1, one(T), wi, 1);
            try blas_calls.gemv(al, io, T, .c, one(T), a.sub(0, i + 1, i, rest), v, 1, zero(T), lower, 1);
            try blas_calls.gemv(al, io, T, .n, minus, w.sub(0, iw + 1, i, rest), lower, 1, one(T), wi, 1);
        }
        if (cx) try blas_calls.scalComplex(io, T, i, tau[i - 1], wi, 1) else try blas_calls.scal(io, T, i, tau[i - 1], wi, 1);
        // W(0:i, iw) += -tau / 2 (wᴴ v) v
        const wv = try blas_calls.dot(io, T, true, i, wi, 1, v, 1);
        const half = mul(T, mulReal(T, tau[i - 1], -0.5), wv);
        try blas_calls.axpy(io, T, i, half, v, 1, wi, 1);
    }
}

/// `?SYTD2`, `?HETD2` upper: the unblocked reduction of the `n x n` `a`,
/// from its last column.
fn sytd2(al: Allocator, io: Io, comptime T: type, a: Mat(T), d: []Float(T), e: []Float(T), tau: []T) Error!void {
    const n = a.rows;
    const cx = comptime isComplex(T);
    if (n == 0) return;
    if (cx) a.at(n - 1, n - 1).* = scalar(T, re(T, a.get(n - 1, n - 1)));
    if (n > 1) {
        const x = try al.alloc(T, n - 1);
        defer al.free(x);
        // `DO I = N - 1, 1, -1`, `i` 0-based: the column `i + 1` above row `i`
        var i = n - 1;
        while (i > 0) {
            i -= 1;
            const v = a.col(i + 1);
            var alpha = v[i];
            const taui = try auxiliary.larfg(io, T, &alpha, Strided(T).of(v[0..i]));
            e[i] = re(T, alpha);
            if (!ops.isZero(T, taui)) {
                v[i] = one(T);
                // x := taui A11 v; x -= taui / 2 (xᴴ v) v; A11 -= v xᴴ + x vᴴ
                const a11 = a.sub(0, 0, i + 1, i + 1);
                try blas_calls.symv(al, io, T, true, taui, a11, v, x.ptr);
                const xv = try blas_calls.dot(io, T, true, i + 1, x.ptr, 1, v, 1);
                const s = mul(T, mulReal(T, taui, -0.5), xv);
                try blas_calls.axpy(io, T, i + 1, s, v, 1, x.ptr, 1);
                rank2Upper(T, a11, v, x.ptr);
            } else if (cx) {
                a.at(i, i).* = scalar(T, re(T, a.get(i, i)));
            }
            v[i] = scalar(T, e[i]);
            d[i + 1] = re(T, a.get(i + 1, i + 1));
            tau[i] = taui;
        }
    }
    d[0] = re(T, a.get(0, 0));
}

/// `?SYR2`, `?HER2` upper with `alpha` -1: `A -= v xᴴ + x vᴴ` on the upper
/// triangle (`dsyr2.f`, `zher2.f`: by columns, `x(j)` and `v(j)` scaled
/// once each; a complex diagonal comes out real).
fn rank2Upper(comptime T: type, a: Mat(T), v: [*]const T, x: [*]const T) void {
    const m = a.rows;
    for (0..m) |j| {
        const t1 = conj(T, x[j]);
        const t2 = conj(T, v[j]);
        const col = a.col(j);
        for (0..j + 1) |r| col[r] = sub(T, col[r], add(T, mul(T, v[r], t1), mul(T, x[r], t2)));
        if (comptime isComplex(T)) col[j] = scalar(T, re(T, col[j]));
    }
}
