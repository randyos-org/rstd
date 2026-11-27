//! `?SYTRD` and `?HETRD`: the blocked reduction of a symmetric (complex:
//! Hermitian) matrix to tridiagonal form, `Qᴴ A Q = T`, on the lower
//! triangle here and the upper in `sytrd_upper.zig`,
//! from Reference-LAPACK as OpenBLAS ships it (`lapack-netlib/SRC/dsytrd.f`,
//! `zhetrd.f`, `dlatrd.f`, `zlatrd.f`, `dsytd2.f`, `zhetd2.f`; OpenBLAS
//! v0.3.34-794-gdbf2c24c9). Blocks of `nb` columns go through `latrd`, the
//! rest of the matrix through one rank-2k update each, the last `nx` columns
//! through `sytd2`. The BLAS calls give the same bits for any worker count.

const std = @import("std");
const ops = @import("ops.zig");
const auxiliary = @import("auxiliary.zig");
const blas_calls = @import("blas_calls.zig");
const strided = @import("../blas/strided.zig");
const Kept = @import("Kept.zig");
const util = @import("../blas/util.zig");
const types = @import("../blas/types.zig");
const interface = @import("interface.zig");
const upper_half = @import("sytrd_upper.zig");

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

/// `ILAENV(1, 'DSYTRD')`, the block (`ilaenv.f`: 32 for `SYTRD`, `HETRD`).
pub const nb = 32;
/// `ILAENV(3, 'DSYTRD')`, the order below which the rest is unblocked
/// (`ilaenv.f`: 32).
pub const nx = 32;

pub const SytrdError = interface.ShapeError || Allocator.Error || Io.Cancelable;

/// `?SYTRD`, `?HETRD` (`dsytrd.f`, `zhetrd.f`): reduces the `n x n`
/// symmetric (Hermitian) `A`, its `uplo` triangle in `a` (the other not
/// read), to the real tridiagonal `T = Qᴴ A Q`: `d` (`n`) its diagonal, `e`
/// (`n - 1`) its off-diagonal, `tau` (`n - 1`) the reflectors' factors, the
/// reflectors in `a` as the reference leaves them (above the superdiagonal
/// for `.u`, `Q = H(n-1) … H(1)`; below the subdiagonal for `.l`, `Q =
/// H(1) … H(n-1)`).
pub fn sytrd(allocator: Allocator, io: Io, comptime T: type, uplo: types.Uplo, n: usize, a: []T, lda: usize, d: []Float(T), e: []Float(T), tau: []T) SytrdError!void {
    comptime ops.requireFloat(T, "sytrd");
    try util.holds(a.len, n, n, lda);
    if (d.len < n or e.len < n -| 1 or tau.len < n -| 1) return error.SliceTooShort;
    const m: Mat(T) = .{ .ptr = a.ptr, .rows = n, .cols = n, .ld = @max(lda, 1) };
    try tridiagonalize(allocator, io, T, uplo == .u, m, d, e, tau);
}

/// `sytrd` after its checks, on `a`'s upper or lower triangle.
pub fn tridiagonalize(child: Allocator, io: Io, comptime T: type, upper: bool, a: Mat(T), d: []Float(T), e: []Float(T), tau: []T) Error!void {
    const n = a.rows;
    if (n == 0) return;
    // one set of packing space and partial vectors for every call (`Kept.zig`)
    var kept: Kept = .init(child);
    defer kept.deinit();
    const al = kept.allocator();
    if (upper) return upper_half.reduce(al, io, T, a, d, e, tau, nb, nx);
    var i: usize = 0;
    if (n > nx) {
        const w = try ops.alloc(T, al, n, nb);
        defer ops.free(T, al, w);
        const tmp = try al.alloc(T, nb);
        defer al.free(tmp);
        // `DO I = 1, N - NX, NB`
        while (i < n - nx) : (i += nb) {
            const m = n - i;
            try latrd(al, io, T, a.sub(i, i, m, m), e[i..], tau[i..], w.sub(0, 0, m, nb), tmp);
            const rest = m - nb;
            try blas_calls.rank2k(al, io, T, false, a.sub(i + nb, i, rest, nb), w.sub(nb, 0, rest, nb), a.sub(i + nb, i + nb, rest, rest));
            for (i..i + nb) |j| {
                a.at(j + 1, j).* = scalar(T, e[j]);
                d[j] = re(T, a.get(j, j));
            }
        }
    }
    try sytd2(al, io, T, a.sub(i, i, n - i, n - i), d[i..], e[i..], tau[i..]);
}

/// `?LATRD` lower: reduces the first `nb` columns of the `n x n` `a`,
/// leaving in `w` (`n x nb`) what the rank-2k update of the rest takes:
/// `A22 -= V Wᴴ + W Vᴴ`. `tmp` holds `nb` elements.
fn latrd(al: Allocator, io: Io, comptime T: type, a: Mat(T), e: []Float(T), tau: []T, w: Mat(T), tmp: []T) Error!void {
    const n = a.rows;
    const cx = comptime isComplex(T);
    const minus = scalar(T, -1);
    for (0..nb) |i| {
        if (cx) a.at(i, i).* = scalar(T, re(T, a.get(i, i)));
        if (i > 0) {
            // A(i:n, i) -= A(i:n, 0:i) conj(W(i, 0:i))ᵀ + W(i:n, 0:i) conj(A(i, 0:i))ᵀ
            const col = a.col(i) + i;
            rowConj(T, w, i, tmp[0..i]);
            try blas_calls.gemv(al, io, T, .n, minus, a.sub(i, 0, n - i, i), tmp.ptr, 1, one(T), col, 1);
            rowConj(T, a, i, tmp[0..i]);
            try blas_calls.gemv(al, io, T, .n, minus, w.sub(i, 0, n - i, i), tmp.ptr, 1, one(T), col, 1);
            if (cx) a.at(i, i).* = scalar(T, re(T, a.get(i, i)));
        }
        if (i + 1 == n) continue;
        // the reflector taking A(i+1:n, i) to e(i) e0
        const m = n - i - 1;
        const v = a.col(i) + i + 1;
        var alpha = v[0];
        tau[i] = try auxiliary.larfg(io, T, &alpha, Strided(T).of(v[1..m]));
        e[i] = re(T, alpha);
        v[0] = one(T);
        // W(i+1:n, i) := A22 v - V W(i+1:n, 0:i)ᴴ v - W V(i+1:n, 0:i)ᴴ v
        const wi = w.col(i) + i + 1;
        try blas_calls.symv(al, io, T, false, one(T), a.sub(i + 1, i + 1, m, m), v, wi);
        if (i > 0) {
            const top = w.col(i);
            try blas_calls.gemv(al, io, T, .c, one(T), w.sub(i + 1, 0, m, i), v, 1, zero(T), top, 1);
            try blas_calls.gemv(al, io, T, .n, minus, a.sub(i + 1, 0, m, i), top, 1, one(T), wi, 1);
            try blas_calls.gemv(al, io, T, .c, one(T), a.sub(i + 1, 0, m, i), v, 1, zero(T), top, 1);
            try blas_calls.gemv(al, io, T, .n, minus, w.sub(i + 1, 0, m, i), top, 1, one(T), wi, 1);
        }
        if (cx) try blas_calls.scalComplex(io, T, m, tau[i], wi, 1) else try blas_calls.scal(io, T, m, tau[i], wi, 1);
        // W(i+1:n, i) += -tau / 2 (wᴴ v) v
        const wv = try blas_calls.dot(io, T, true, m, wi, 1, v, 1);
        const half = mul(T, mulReal(T, tau[i], -0.5), wv);
        try blas_calls.axpy(io, T, m, half, v, 1, wi, 1);
    }
}

/// `conj(M(i, 0:k))` into `out` (`k` long), as `?LACGV` leaves the row for
/// `?GEMV`; plain for a real `T`.
fn rowConj(comptime T: type, m: Mat(T), i: usize, out: []T) void {
    for (out, 0..) |*o, j| o.* = conj(T, m.get(i, j));
}

/// `?SYTD2`, `?HETD2` lower: the unblocked reduction of the `n x n` `a`.
fn sytd2(al: Allocator, io: Io, comptime T: type, a: Mat(T), d: []Float(T), e: []Float(T), tau: []T) Error!void {
    const n = a.rows;
    const cx = comptime isComplex(T);
    if (n == 0) return;
    if (cx) a.at(0, 0).* = scalar(T, re(T, a.get(0, 0)));
    if (n > 1) {
        const x = try al.alloc(T, n - 1);
        defer al.free(x);
        for (0..n - 1) |i| {
            const m = n - i - 1;
            const v = a.col(i) + i + 1;
            var alpha = v[0];
            const taui = try auxiliary.larfg(io, T, &alpha, Strided(T).of(v[1..m]));
            e[i] = re(T, alpha);
            if (!ops.isZero(T, taui)) {
                v[0] = one(T);
                // x := taui A22 v; x -= taui / 2 (xᴴ v) v; A22 -= v xᴴ + x vᴴ
                const a22 = a.sub(i + 1, i + 1, m, m);
                try blas_calls.symv(al, io, T, false, taui, a22, v, x.ptr);
                const xv = try blas_calls.dot(io, T, true, m, x.ptr, 1, v, 1);
                const s = mul(T, mulReal(T, taui, -0.5), xv);
                try blas_calls.axpy(io, T, m, s, v, 1, x.ptr, 1);
                rank2Lower(T, a22, v, x.ptr);
            } else if (cx) {
                a.at(i + 1, i + 1).* = scalar(T, re(T, a.get(i + 1, i + 1)));
            }
            v[0] = scalar(T, e[i]);
            d[i] = re(T, a.get(i, i));
            tau[i] = taui;
        }
    }
    d[n - 1] = re(T, a.get(n - 1, n - 1));
}

/// `?SYR2`, `?HER2` lower with `alpha` -1: `A -= v xᴴ + x vᴴ` on the lower
/// triangle (`dsyr2.f`, `zher2.f`: by columns, `x(j)` and `v(j)` scaled
/// once each; a complex diagonal comes out real).
fn rank2Lower(comptime T: type, a: Mat(T), v: [*]const T, x: [*]const T) void {
    const m = a.rows;
    for (0..m) |j| {
        const t1 = conj(T, x[j]);
        const t2 = conj(T, v[j]);
        const col = a.col(j);
        for (j..m) |r| col[r] = sub(T, col[r], add(T, mul(T, v[r], t1), mul(T, x[r], t2)));
        if (comptime isComplex(T)) col[j] = scalar(T, re(T, col[j]));
    }
}
