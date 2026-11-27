//! The unblocked real QR and RQ of `dtgex2`'s blocks (at most 4 x 4),
//! Reference-LAPACK 3.12.0 (`lapack-netlib/SRC/dgeqr2.f`, `dorg2r.f`,
//! `dorm2r.f`, `dgerq2.f`, `dorgr2.f`, `dormr2.f`, `dlarf.f`, `dlarf1f.f`,
//! `dlarf1l.f`, OpenBLAS v0.3.34-317; modified BSD): reflectors one at a
//! time, in the reference's order, the unit element written as `dlarf.f`'s
//! callers do.

const std = @import("std");
const strided = @import("../blas/strided.zig");
const ops = @import("ops.zig");
const auxiliary = @import("auxiliary.zig");

const Io = std.Io;
const Error = ops.Error;
const Strided = strided.Strided;

/// `n` elements of `p`, `inc` apart, as `larfg` reads them.
fn vec(comptime F: type, p: [*]F, n: usize, inc: usize) Strided(F) {
    return .{ .data = p[0..if (n == 0) 0 else (n - 1) * inc + 1], .n = n, .stride = @intCast(inc) };
}

/// `DLARF` (`dlarf.f`), all of `v`: `C := H C` (`left`, `C` `m x n`, `v`
/// of `m`) or `C := C H` (`v` of `n`), `H = I - tau v vᵀ`, `v` `incv` apart.
pub fn larf(comptime F: type, left: bool, m: usize, n: usize, v: [*]const F, incv: usize, tau: F, c: [*]F, ldc: usize) void {
    if (tau == 0) return;
    if (left) {
        for (0..n) |j| {
            var w: F = 0;
            for (0..m) |i| w += c[i + j * ldc] * v[i * incv];
            for (0..m) |i| c[i + j * ldc] -= tau * v[i * incv] * w;
        }
    } else {
        for (0..m) |i| {
            var w: F = 0;
            for (0..n) |j| w += c[i + j * ldc] * v[j * incv];
            for (0..n) |j| c[i + j * ldc] -= tau * w * v[j * incv];
        }
    }
}

/// `DGEQR2` (`dgeqr2.f`): `A = Q R`, the reflectors below `R`.
pub fn geqr2(io: Io, comptime F: type, m: usize, n: usize, a: [*]F, lda: usize, tau: []F) Error!void {
    for (0..@min(m, n)) |i| {
        const d = a + i + i * lda;
        tau[i] = try auxiliary.larfg(io, F, &d[0], vec(F, a + @min(i + 1, m - 1) + i * lda, m - i - 1, 1));
        if (i + 1 < n) {
            const aii = d[0];
            d[0] = 1;
            larf(F, true, m - i, n - i - 1, d, 1, tau[i], d + lda, lda);
            d[0] = aii;
        }
    }
}

/// `DORG2R` (`dorg2r.f`): the first `n` columns of `Q = H(1) … H(k)`.
pub fn org2r(comptime F: type, m: usize, n: usize, k: usize, a: [*]F, lda: usize, tau: []const F) void {
    for (k..n) |j| {
        for (0..m) |l| a[l + j * lda] = 0;
        a[j + j * lda] = 1;
    }
    var i = k;
    while (i > 0) {
        i -= 1;
        const d = a + i + i * lda;
        if (i + 1 < n) {
            d[0] = 1;
            larf(F, true, m - i, n - i - 1, d, 1, tau[i], d + lda, lda);
        }
        for (1..m - i) |l| d[l] *= -tau[i];
        d[0] = 1 - tau[i];
        for (0..i) |l| a[l + i * lda] = 0;
    }
}

/// `DORM2R` (`dorm2r.f`), `nq x nq` square `Q = H(1) … H(k)` from `geqr2`:
/// `C := op(Q) C` (`left`) or `C op(Q)`, `op` the transpose with `trans`.
pub fn orm2r(comptime F: type, left: bool, trans: bool, m: usize, n: usize, k: usize, a: [*]F, lda: usize, tau: []const F, c: [*]F, ldc: usize) void {
    const forward = left == trans;
    for (0..k) |s| {
        const i = if (forward) s else k - 1 - s;
        const d = a + i + i * lda;
        const aii = d[0];
        d[0] = 1;
        if (left) larf(F, true, m - i, n, d, 1, tau[i], c + i, ldc) else larf(F, false, m, n - i, d, 1, tau[i], c + i * ldc, ldc);
        d[0] = aii;
    }
}

/// `DGERQ2` (`dgerq2.f`): `A = R Q`, the reflectors left of `R` in `A`'s
/// last rows.
pub fn gerq2(io: Io, comptime F: type, m: usize, n: usize, a: [*]F, lda: usize, tau: []F) Error!void {
    const k = @min(m, n);
    var i = k;
    while (i > 0) {
        i -= 1;
        const row = m - k + i;
        const col = n - k + i;
        const d = a + row + col * lda;
        tau[i] = try auxiliary.larfg(io, F, &d[0], vec(F, a + row, col, lda));
        const aii = d[0];
        d[0] = 1;
        larf(F, false, row, col + 1, a + row, lda, tau[i], a, lda);
        d[0] = aii;
    }
}

/// `DORGR2` (`dorgr2.f`): the last `m` rows of `Q = H(1) … H(k)` from
/// `gerq2`, `m <= n`.
pub fn orgr2(comptime F: type, m: usize, n: usize, k: usize, a: [*]F, lda: usize, tau: []const F) void {
    if (k < m) for (0..n) |j| {
        for (0..m - k) |l| a[l + j * lda] = 0;
        if (j + m >= n and j + k < n) a[j + m - n + j * lda] = 1;
    };
    for (0..k) |i| {
        const ii = m - k + i;
        const col = n - m + ii;
        a[ii + col * lda] = 1;
        larf(F, false, ii, col + 1, a + ii, lda, tau[i], a, lda);
        for (0..col) |l| a[ii + l * lda] *= -tau[i];
        a[ii + col * lda] = 1 - tau[i];
        for (col + 1..n) |l| a[ii + l * lda] = 0;
    }
}

/// `DORMR2` (`dormr2.f`), `nq x nq` square `Q = H(1) … H(k)` from `gerq2`
/// (rows `0..k` of `a`): `C := op(Q) C` (`left`) or `C op(Q)`.
pub fn ormr2(comptime F: type, left: bool, trans: bool, m: usize, n: usize, k: usize, a: [*]F, lda: usize, tau: []const F, c: [*]F, ldc: usize) void {
    const nq = if (left) m else n;
    const forward = left == trans;
    for (0..k) |s| {
        const i = if (forward) s else k - 1 - s;
        const d = a + i + (nq - k + i) * lda;
        const aii = d[0];
        d[0] = 1;
        if (left) larf(F, true, m - k + i + 1, n, a + i, lda, tau[i], c, ldc) else larf(F, false, m, n - k + i + 1, a + i, lda, tau[i], c, ldc);
        d[0] = aii;
    }
}
