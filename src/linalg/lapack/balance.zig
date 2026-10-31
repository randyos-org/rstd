//! Balancing, from Reference-LAPACK 3.12.0 as OpenBLAS v0.3.34-317 ships it
//! (`lapack-netlib/SRC/`, modified BSD): `xGEBAL` permutes `A` to isolate
//! eigenvalues and scales rows and columns by powers of 2 toward equal
//! norms; `xGEBAK` takes eigenvectors of the balanced matrix back. Neither
//! faer nor libflame has its own (libflame's are the reference through f2c);
//! `dgeev.f` balances before the Hessenberg reduction, and so does `geev`
//! (`evd.zig`, through `balance` and `back`, job `B`).
//!
//! `scale` and `ilo`, `ihi` are the reference's: 1-based, a permuted row's
//! partner stored as a float index. Row and column norms are rstd's `nrm2`;
//! swaps and scalings go through its BLAS.

const std = @import("std");
const blas = @import("../blas/__root__.zig");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const ops = @import("ops.zig");
const norm = @import("norm.zig");
const chol = @import("cholesky.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const isComplex = types.isComplex;
const Side = types.Side;
const Error = norm.Error;
const Mat = ops.Mat;

/// `JOB`: nothing (`n`), permute only (`p`), scale only (`s`), or both
/// (`b`).
pub const BalanceJob = enum { n, p, s, b };

/// `SCLFAC`: scalings are powers of 2.
const sclfac = 2;
/// `FACTOR`: a scaling must cut `c + r` below this fraction of it.
const factor = 0.95;

/// What `gebal` returns: `InvalidA` for a NaN met while scaling (the
/// reference's `INFO = -3`, `A` then part balanced).
pub const GebalError = error{ InvalidN, InvalidLda, InvalidA } || Io.Cancelable;

/// `?GEBAL` (`dgebal.f`, `zgebal.f`): balances the `n x n` `A` in place,
/// `D⁻¹ Pᵀ A P D`. Rows and columns outside `ilo..ihi` (1-based) are then
/// upper triangular. `scale` (`n`) gets, outside `ilo..ihi`, the index
/// each row was swapped with (in the reference's order: from `n` down to
/// `ihi + 1`, then up from 1), and inside it `D`'s entries.
pub fn gebal(allocator: Allocator, io: Io, comptime T: type, job: BalanceJob, n: isize, a: [*]T, lda: isize, ilo: *isize, ihi: *isize, scale: [*]Float(T)) GebalError!void {
    comptime chol.floatType(T, "gebal");
    _ = allocator;
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, n)) return error.InvalidLda;
    if (n == 0) {
        ilo.* = 1;
        ihi.* = 0;
        return;
    }
    if (!util.fitsMatrix(n, n, lda, @sizeOf(T))) return error.InvalidLda;
    const lo, const hi = balanceJob(io, T, job, @intCast(n), a, @intCast(lda), scale) catch |e| return switch (e) {
        error.InvalidA => error.InvalidA,
        error.OutOfMemory => unreachable,
        error.Canceled => error.Canceled,
    };
    ilo.* = @intCast(lo);
    ihi.* = @intCast(hi);
}

/// `gebal` after its checks: returns `ilo`, `ihi` (1-based).
pub fn balanceJob(io: Io, comptime T: type, job: BalanceJob, n: usize, a: [*]T, lda: usize, scale: [*]Float(T)) (Error || error{InvalidA})![2]usize {
    const F = Float(T);
    if (job == .n) {
        for (0..n) |i| scale[i] = 1;
        return .{ 1, n };
    }
    // `k`, `l`: the 1-based active window
    var k: usize = 1;
    var l: usize = n;
    if (job != .s) {
        // rows isolating an eigenvalue go to the bottom
        var noconv = true;
        while (noconv) {
            noconv = false;
            var i = l;
            while (i >= 1) : (i -= 1) {
                const can = for (1..l + 1) |j| {
                    if (i != j and !ops.isZero(T, at(a, lda, i, j))) break false;
                } else true;
                if (can) {
                    scale[l - 1] = @floatFromInt(i);
                    if (i != l) {
                        try swap(io, T, l, a + (i - 1) * lda, 1, a + (l - 1) * lda, 1);
                        try swap(io, T, n - k + 1, a + (i - 1) + (k - 1) * lda, lda, a + (l - 1) + (k - 1) * lda, lda);
                    }
                    noconv = true;
                    if (l == 1) return .{ 1, 1 };
                    l -= 1;
                }
                if (i == 1) break;
            }
        }
        // columns isolating one go to the left
        noconv = true;
        while (noconv) {
            noconv = false;
            const last = l;
            var j = k;
            while (j <= last) : (j += 1) {
                const can = for (k..l + 1) |i| {
                    if (i != j and !ops.isZero(T, at(a, lda, i, j))) break false;
                } else true;
                if (can) {
                    scale[k - 1] = @floatFromInt(j);
                    if (j != k) {
                        try swap(io, T, l, a + (j - 1) * lda, 1, a + (k - 1) * lda, 1);
                        try swap(io, T, n - k + 1, a + (j - 1) + (k - 1) * lda, lda, a + (k - 1) + (k - 1) * lda, lda);
                    }
                    noconv = true;
                    k += 1;
                }
            }
        }
    }
    for (k - 1..l) |i| scale[i] = 1;
    if (job == .p) return .{ k, l };
    // scale toward equal row and column norms
    const sfmin1 = std.math.floatMin(F) / std.math.floatEps(F);
    const sfmax1 = 1 / sfmin1;
    const sfmin2 = sfmin1 * sclfac;
    const sfmax2 = 1 / sfmin2;
    var noconv = true;
    while (noconv) {
        noconv = false;
        for (k..l + 1) |i| {
            const w = l - k + 1;
            var c = ops.nrm2(io, T, w, a + (k - 1) + (i - 1) * lda, 1);
            var r = ops.nrm2(io, T, w, a + (i - 1) + (k - 1) * lda, lda);
            const ica = norm.iamax(T, l, a + (i - 1) * lda, 1);
            var ca = ops.abs(T, a[ica + (i - 1) * lda]);
            const ira = norm.iamax(T, n - k + 1, a + (i - 1) + (k - 1) * lda, lda);
            var ra = ops.abs(T, a[(i - 1) + (ira + k - 1) * lda]);
            if (c == 0 or r == 0) continue;
            if (std.math.isNan(c + ca + r + ra)) return error.InvalidA;
            var g = r / sclfac;
            var f: F = 1;
            const s = c + r;
            while (c < g and @max(f, c, ca) < sfmax2 and @min(r, g, ra) > sfmin2) {
                f *= sclfac;
                c *= sclfac;
                ca *= sclfac;
                r /= sclfac;
                g /= sclfac;
                ra /= sclfac;
            }
            g = c / sclfac;
            while (g >= r and @max(r, ra) < sfmax2 and @min(f, c, g, ca) > sfmin2) {
                f /= sclfac;
                c /= sclfac;
                g /= sclfac;
                ca /= sclfac;
                r *= sclfac;
                ra *= sclfac;
            }
            if (c + r >= factor * s) continue;
            if (f < 1 and scale[i - 1] < 1 and f * scale[i - 1] <= sfmin1) continue;
            if (f > 1 and scale[i - 1] > 1 and scale[i - 1] >= sfmax1 / f) continue;
            g = 1 / f;
            scale[i - 1] *= f;
            noconv = true;
            try norm.scal(io, T, n - k + 1, g, a + (i - 1) + (k - 1) * lda, lda);
            try norm.scal(io, T, l, f, a + (i - 1) * lda, 1);
        }
    }
    return .{ k, l };
}

/// `A(i, j)`, 1-based.
fn at(a: anytype, lda: usize, i: usize, j: usize) @typeInfo(@TypeOf(a)).pointer.child {
    return a[(i - 1) + (j - 1) * lda];
}

/// `?SWAP` through rstd's BLAS.
fn swap(io: Io, comptime T: type, n: usize, x: [*]T, incx: usize, y: [*]T, incy: usize) Error!void {
    if (n == 0) return;
    blas.swap(io, T, @intCast(n), x, @intCast(incx), y, @intCast(incy)) catch |e| return norm.narrow(e);
}

/// What `gebak` returns.
pub const GebakError = error{ InvalidN, InvalidIlo, InvalidIhi, InvalidM, InvalidLdv } || Io.Cancelable;

/// `?GEBAK` (`dgebak.f`, `zgebak.f`): `V := P D V` for right eigenvectors
/// (`side` `r`) or `V := P D⁻¹ V` for left ones (`l`), the `n x m` `V`'s
/// rows, from `gebal`'s `job`, `ilo`, `ihi` and `scale`.
pub fn gebak(allocator: Allocator, io: Io, comptime T: type, job: BalanceJob, side: Side, n: isize, ilo: isize, ihi: isize, scale: [*]const Float(T), m: isize, v: [*]T, ldv: isize) GebakError!void {
    comptime chol.floatType(T, "gebak");
    _ = allocator;
    if (n < 0) return error.InvalidN;
    if (ilo < 1 or ilo > @max(1, n)) return error.InvalidIlo;
    if (ihi < @min(ilo, n) or ihi > n) return error.InvalidIhi;
    if (m < 0) return error.InvalidM;
    if (ldv < @max(1, n)) return error.InvalidLdv;
    if (n == 0 or m == 0 or job == .n) return;
    if (!util.fitsMatrix(n, m, ldv, @sizeOf(T))) return error.InvalidLdv;
    backJob(io, T, job, side, @intCast(n), @intCast(ilo), @intCast(ihi), scale, @intCast(m), v, @intCast(ldv)) catch |e| return switch (e) {
        error.OutOfMemory => unreachable,
        error.Canceled => error.Canceled,
    };
}

/// `gebak` after its checks (`ilo`, `ihi` 1-based).
pub fn backJob(io: Io, comptime T: type, job: BalanceJob, side: Side, n: usize, ilo: usize, ihi: usize, scale: [*]const Float(T), m: usize, v: [*]T, ldv: usize) Error!void {
    if (ilo != ihi and (job == .s or job == .b)) {
        for (ilo..ihi + 1) |i| {
            const s = if (side == .r) scale[i - 1] else 1 / scale[i - 1];
            try norm.scal(io, T, m, s, v + (i - 1), ldv);
        }
    }
    if (job == .p or job == .b) {
        for (1..n + 1) |ii| {
            var i = ii;
            if (i >= ilo and i <= ihi) continue;
            if (i < ilo) i = ilo - ii;
            const k: usize = @intFromFloat(scale[i - 1]);
            if (k == i) continue;
            try swap(io, T, m, v + (i - 1), ldv, v + (k - 1), ldv);
        }
    }
}

// ---- geev's job B -------------------------------------------------------------------

/// What `balance` found: `a`'s rows and columns outside `ilo..ihi`
/// (0-based, `ihi` past the end) are already triangular.
pub const Info = struct { ilo: usize, ihi: usize };

/// `xGEBAL` with job `B` on the `n x n` `a`, which holds no NaN; `scale`
/// (`n`) as `gebal`'s.
pub fn balance(io: Io, comptime T: type, a: Mat(T), scale: []Float(T)) Error!Info {
    const r = balanceJob(io, T, .b, a.rows, a.ptr, a.ld, scale.ptr) catch |e| switch (e) {
        error.InvalidA => unreachable,
        error.OutOfMemory, error.Canceled => |x| return x,
    };
    return .{ .ilo = r[0] - 1, .ihi = r[1] };
}

/// `xGEBAK` with job `B`: `v`'s rows (the balanced matrix's right vectors,
/// or with `left` its left ones) back to `a`'s.
pub fn back(io: Io, comptime T: type, info: Info, scale: []const Float(T), left: bool, v: Mat(T)) Error!void {
    if (v.rows == 0 or v.cols == 0) return;
    try backJob(io, T, .b, if (left) .l else .r, v.rows, info.ilo + 1, info.ihi, scale.ptr, v.cols, v.ptr, v.ld);
}

// ---- tests -----------------------------------------------------------------------

const testing = std.testing;
const check = @import("../blas/_test/blas1/check.zig");

const test_types = chol.test_types;

/// `TESTING/dbal.in` of Reference-LAPACK 3.12.0 (OpenBLAS v0.3.34-317's
/// `lapack-netlib/`), verbatim: `dchkbl.f`'s cases for `xGEBAL` with job
/// `B`, each `n`, `A` by rows, the expected `ilo`, `ihi`, balanced `A` and
/// `scale`. Copy it again rather than editing it if it goes stale.
const dbal_in =
    \\DBL:  Tests DGEBAL
    \\  5
    \\  0.1000D+01  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00
    \\  0.0000D+00  0.2000D+01  0.0000D+00  0.0000D+00  0.0000D+00
    \\  0.0000D+00  0.0000D+00  0.3000D+01  0.0000D+00  0.0000D+00
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.4000D+01  0.0000D+00
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00  0.5000D+01
    \\
    \\   1   1
    \\  0.1000D+01  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00
    \\  0.0000D+00  0.2000D+01  0.0000D+00  0.0000D+00  0.0000D+00
    \\  0.0000D+00  0.0000D+00  0.3000D+01  0.0000D+00  0.0000D+00
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.4000D+01  0.0000D+00
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00  0.5000D+01
    \\
    \\  0.1000D+01  0.2000D+01  0.3000D+01  0.4000D+01  0.5000D+01
    \\
    \\  5
    \\  0.1000D+01  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00
    \\  0.1000D+01  0.2000D+01  0.0000D+00  0.0000D+00  0.0000D+00
    \\  0.1000D+01  0.2000D+01  0.3000D+01  0.0000D+00  0.0000D+00
    \\  0.1000D+01  0.2000D+01  0.3000D+01  0.4000D+01  0.0000D+00
    \\  0.1000D+01  0.2000D+01  0.3000D+01  0.4000D+01  0.5000D+01
    \\
    \\   1   1
    \\  0.5000D+01  0.4000D+01  0.3000D+01  0.2000D+01  0.1000D+01
    \\  0.0000D+00  0.4000D+01  0.3000D+01  0.2000D+01  0.1000D+01
    \\  0.0000D+00  0.0000D+00  0.3000D+01  0.2000D+01  0.1000D+01
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.2000D+01  0.1000D+01
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00  0.1000D+01
    \\
    \\  0.1000D+01  0.2000D+01  0.3000D+01  0.2000D+01  0.1000D+01
    \\
    \\  5
    \\  0.1000D+01  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00
    \\  0.1000D+01  0.1000D+01  0.0000D+00  0.0000D+00  0.0000D+00
    \\  0.0000D+00  0.1000D+01  0.1000D+01  0.0000D+00  0.0000D+00
    \\  0.0000D+00  0.0000D+00  0.1000D+01  0.1000D+01  0.0000D+00
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.1000D+01  0.1000D+01
    \\
    \\   1   1
    \\  0.1000D+01  0.1000D+01  0.0000D+00  0.0000D+00  0.0000D+00
    \\  0.0000D+00  0.1000D+01  0.1000D+01  0.0000D+00  0.0000D+00
    \\  0.0000D+00  0.0000D+00  0.1000D+01  0.1000D+01  0.0000D+00
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.1000D+01  0.1000D+01
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00  0.1000D+01
    \\
    \\  0.1000D+01  0.2000D+01  0.3000D+01  0.2000D+01  0.1000D+01
    \\
    \\  4
    \\  0.0000D+00  0.2000D+01  0.1000D+00  0.0000D+00
    \\  0.2000D+01  0.0000D+00  0.0000D+00  0.1000D+00
    \\  0.1000D+03  0.0000D+00  0.0000D+00  0.2000D+01
    \\  0.0000D+00  0.1000D+03  0.2000D+01  0.0000D+00
    \\
    \\   1   4
    \\  0.0000D-03  2.0000D+00  3.2000D+00  0.0000D-03
    \\  2.0000D+00  0.0000D-03  0.0000D-03  3.2000D+00
    \\  3.1250D+00  0.0000D-03  0.0000D-03  2.0000D+00
    \\  0.0000D-03  3.1250D+00  2.0000D+00  0.0000D-03
    \\
    \\  62.5000D-03 62.5000D-03  2.0000D+00  2.0000D+00
    \\
    \\  6
    \\  0.2000D+01  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00  0.1024D+04
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00  0.1280D+03
    \\  0.0000D+00  0.2000D+01  0.3000D+04  0.0000D+00  0.0000D+00  0.2000D+01
    \\  0.1280D+03  0.4000D+01  0.4000D-02  0.5000D+01  0.6000D+03  0.8000D+01
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00  0.2000D-02  0.2000D+01
    \\  0.8000D+01  0.8192D+04  0.0000D+00  0.0000D+00  0.0000D+00  0.2000D+01
    \\
    \\   4   6
    \\  0.5000D+01  0.4000D-02  0.6000D+03  0.1024D+04  0.5000D+00  0.8000D+01
    \\  0.0000D+00  0.3000D+04  0.0000D+00  0.0000D+00  0.2500D+00  0.2000D+01
    \\  0.0000D+00  0.0000D+00  0.2000D-02  0.0000D+00  0.0000D+00  0.2000D+01
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.2000D+01  0.0000D+00  0.1280D+03
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00  0.1024D+04
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.6400D+02  0.1024D+04  0.2000D+01
    \\
    \\  0.4000D+01  0.3000D+01  0.5000D+01  0.8000D+01  0.1250D+00  0.1000D+01
    \\
    \\  5
    \\  0.1000D+01  0.0000D+00  0.0000D+00  0.0000D+00  0.8000D+01
    \\  0.0000D+00  0.2000D+01  0.8192D+04  0.2000D+01  0.4000D+01
    \\  0.2500D-03  0.1250D-03  0.4000D+01  0.0000D+00  0.6400D+02
    \\  0.0000D+00  0.2000D+01  0.1024D+04  0.4000D+01  0.8000D+01
    \\  0.0000D+00  0.8192D+04  0.0000D+00  0.0000D+00  0.8000D+01
    \\
    \\   1   5
    \\  1.0000D+00     0.0000D-03     0.0000D-03     0.0000D-03     2.0000D+00
    \\  0.0000D-03     2.0000D+00     1.0240D+03    16.0000D+00    16.0000D+00
    \\  3.2000D-02     1.0000D-03     4.0000D+00     0.0000D-03     2.0480D+03
    \\  0.0000D-03   250.0000D-03    16.0000D+00     4.0000D+00     4.0000D+00
    \\  0.0000D-03     2.0480D+03     0.0000D-03     0.0000D-03     8.0000D+00
    \\
    \\  8.0000D+00  500.0000D-03  62.5000D-03  4.0000D+00  2.0000D+00
    \\
    \\  4
    \\  0.1000D+01  0.1000D+07  0.1000D+07  0.1000D+07
    \\ -0.2000D+07  0.3000D+01  0.2000D-05  0.3000D-05
    \\ -0.3000D+07  0.0000D+00  0.1000D-05  0.2000D+01
    \\  0.1000D+07  0.0000D+00  0.3000D-05  0.4000D+07
    \\
    \\   1   4
    \\  1.0000D+00     1.0000D+06     2.0000D+06     1.0000D+06
    \\ -2.0000D+06     3.0000D+00     4.0000D-06     3.0000D-06
    \\ -1.5000D+06     0.0000D-03     1.0000D-06     1.0000D+00
    \\  1.0000D+06     0.0000D-03     6.0000D-06     4.0000D+06
    \\
    \\  1.0000D+00  1.0000D+00 2.0000D+00  1.0000D+00
    \\
    \\   4
    \\  0.1000D+01  0.1000D+05  0.1000D+05  0.1000D+05
    \\ -0.2000D+05  0.3000D+01  0.2000D-02  0.3000D-02
    \\  0.0000D+00  0.2000D+01  0.0000D+00 -0.3000D+05
    \\  0.0000D+00  0.0000D+00  0.1000D+05  0.0000D+00
    \\
    \\   1   4
    \\  1.0000D+00    10.0000D+03    10.0000D+03     5.0000D+03
    \\ -20.0000D+03     3.0000D+00     2.0000D-03     1.5000D-03
    \\  0.0000D-03     2.0000D+00     0.0000D-03   -15.0000D+03
    \\  0.0000D-03     0.0000D-03    20.0000D+03     0.0000D-03
    \\
    \\   1.0000D+00     1.0000D+00     1.0000D+00   500.0000D-03
    \\
    \\  5
    \\  0.1000D+01  0.5120D+03  0.4096D+04  3.2768D+04  2.62144D+05
    \\  0.8000D+01  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00
    \\  0.0000D+00  0.8000D+01  0.0000D+00  0.0000D+00  0.0000D+00
    \\  0.0000D+00  0.0000D+00  0.8000D+01  0.0000D+00  0.0000D+00
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.8000D+01  0.0000D+00
    \\
    \\   1   5
    \\  1.0000D+00    32.0000D+00   32.0000D+00  32.0000D+000 32.0000D+00
    \\  128.0000D+00  0.0000D-03    0.0000D-03   0.0000D-003  0.0000D-03
    \\  0.0000D-03    64.0000D+00   0.0000D-03   0.0000D-003  0.0000D-03
    \\  0.0000D-03    0.0000D-03    64.0000D+00  0.0000D-003  0.0000D-03
    \\  0.0000D-03    0.0000D-03    0.0000D-03   64.0000D+000 0.0000D-03
    \\
    \\  256.0000D+00  16.0000D+00  2.0000D+00  250.0000D-03  31.2500D-03
    \\
    \\  6
    \\  0.1000D+01  0.1000D+01  0.0000D+00  0.1000D+01  0.1000D+01  0.1000D+01
    \\  0.1000D+01  0.1000D+01  0.0000D+00  0.1000D+01  0.1000D+01  0.1000D+01
    \\  0.1000D+01  0.1000D+01  0.1000D+01  0.1000D+01  0.1000D+01  0.1000D+01
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.1000D+01  0.0000D+00  0.0000D+00
    \\  0.1000D+01  0.1000D+01  0.0000D+00  0.1000D+01  0.1000D+01  0.1000D+01
    \\  0.1000D+01  0.1000D+01  0.0000D+00  0.1000D+01  0.1000D+01  0.1000D+01
    \\
    \\   2   5
    \\  0.1000D+01  0.1000D+01  0.1000D+01  0.1000D+01  0.1000D+01  0.1000D+01
    \\  0.0000D+00  0.1000D+01  0.1000D+01  0.1000D+01  0.1000D+01  0.1000D+01
    \\  0.0000D+00  0.1000D+01  0.1000D+01  0.1000D+01  0.1000D+01  0.1000D+01
    \\  0.0000D+00  0.1000D+01  0.1000D+01  0.1000D+01  0.1000D+01  0.1000D+01
    \\  0.0000D+00  0.1000D+01  0.1000D+01  0.1000D+01  0.1000D+01  0.1000D+01
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00  0.1000D+01
    \\
    \\  0.3000D+01  0.1000D+01  0.1000D+01  0.1000D+01  0.1000D+01  0.4000D+01
    \\
    \\  7
    \\  0.6000D+01  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00  0.1000D+01 0.0000D+00
    \\  0.0000D+00  0.4000D+01  0.0000D+00  0.2500D-03  0.1250D-01  0.2000D-01 0.1250D+00
    \\  0.1000D+01  0.1280D+03  0.6400D+02  0.0000D+00  0.0000D+00 -0.2000D+01 0.1600D+02
    \\  0.0000D+00  1.6384D+04  0.0000D+00  0.1000D+01 -0.4000D+03  0.2560D+03 -0.4000D+04
    \\ -0.2000D+01 -0.2560D+03  0.0000D+00  0.1250D-01  0.2000D+01  0.2000D+01 0.3200D+02
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00 0.0000D+00
    \\  0.0000D+00  0.8000D+01  0.0000D+00  0.4000D-02  0.1250D+00 -0.2000D+00 0.3000D+01
    \\
    \\  2   5
    \\  6.4000D+01   1.0000D+00   5.00000D-01   0.0000D+00   0.0000D+00   1.0000D+00  -2.0000D+00
    \\  0.0000D+00   4.0000D+00   5.00000D-01   1.0240D+00   8.0000D-01   0.0000D+00   2.5600D+00
    \\  0.0000D+00   2.0000D+00   3.00000D+00   4.0960D+00   2.0000D+00   0.0000D+00  -6.4000D+00
    \\  0.0000D+00   4.0000D+00  -3.90625D+00   1.0000D+00  -6.2500D+00   0.0000D+00   8.0000D+00
    \\  0.0000D+00  -4.0000D+00   2.00000D+00   8.0000D-01   2.0000D+00  -4.0000D+00   4.0000D+00
    \\  0.0000D+00   0.0000D+00   0.00000D+00   0.0000D+00   0.0000D+00   6.0000D+00   1.0000D+00
    \\  0.0000D+00   0.0000D+00   0.00000D+00   0.0000D+00   0.0000D+00   0.0000D+00   0.0000D+00
    \\
    \\  3.0000D+00  7.812500D-03  3.1250D-02  3.2000D+01  5.0000D-01  1.0000D+00 6.0000D+00
    \\
    \\  5
    \\  0.1000D+04  0.2000D+01  0.3000D+01  0.4000D+01  0.5000D+06
    \\  0.9000D+01  0.0000D+00  0.2000D-03  0.1000D+01  0.3000D+01
    \\  0.0000D+00 -0.3000D+03  0.2000D+01  0.1000D+01  0.1000D+01
    \\  0.9000D+01  0.2000D-02  0.1000D+01  0.1000D+01 -0.1000D+04
    \\  0.6000D+01  0.2000D+03  0.1000D+01  0.6000D+03  0.3000D+01
    \\
    \\  1   5
    \\  1.0000D+03     3.1250D-02     3.7500D-01     3.1250D-02     1.95312500D+03
    \\  5.7600D+02     0.0000D+00     1.6000D-03     5.0000D-01     7.50000000D-01
    \\  0.0000D+00    -3.7500D+01     2.0000D+00     6.2500D-02     3.12500000D-02
    \\  1.1520D+03     4.0000D-03     1.6000D+01     1.0000D+00    -5.00000000D+02
    \\  1.5360D+03     8.0000D+02     3.2000D+01     1.2000D+03     3.00000000D+00
    \\
    \\  3.2000D+01  5.0000D-01   4.0000D+00   2.5000D-01  1.2500D-01
    \\
    \\  6
    \\  1.0000D+00  1.0000D+120 0.0000D+00  0.0000D+00  0.0000D+00 0.0000D+00
    \\  1.0000D-120 1.0000D+00  1.0000D+120 0.0000D+00  0.0000D+00 0.0000D+00
    \\  0.0000D+00  1.0000D-120 1.0000D+00  1.0000D+120 0.0000D+00 0.0000D+00
    \\  0.0000D+00  0.0000D+00  1.0000D-120 1.0000D+00  1.0000D+120 0.0000D+00
    \\  0.0000D+00  0.0000D+00  0.0000D+00  1.0000D-120 1.0000D+00 1.0000D+120
    \\  0.0000D+00  0.0000D+00  0.0000D+00  0.0000D+00  1.0000D-120 1.0000D+00
    \\
    \\   1   6
    \\0.10000000000000000000D+01 0.63448545932891229313D+04 0.00000000000000000000D+00 0.00000000000000000000D+00 0.00000000000000000000D+00 0.00000000000000000000D+00
    \\0.15760802478557791348D-03 0.10000000000000000000D+01 0.63448545932891229313D+04 0.00000000000000000000D+00 0.00000000000000000000D+00 0.00000000000000000000D+00
    \\0.00000000000000000000D+00 0.15760802478557791348D-03 0.10000000000000000000D+01 0.31724272966445614657D+04 0.00000000000000000000D+00 0.00000000000000000000D+00
    \\0.00000000000000000000D+00 0.00000000000000000000D+00 0.31521604957115582695D-03 0.10000000000000000000D+01 0.15862136483222807328D+04 0.00000000000000000000D+00
    \\0.00000000000000000000D+00 0.00000000000000000000D+00 0.00000000000000000000D+00 0.63043209914231165391D-03 0.10000000000000000000D+01 0.79310682416114036641D+03
    \\0.00000000000000000000D+00 0.00000000000000000000D+00 0.00000000000000000000D+00 0.00000000000000000000D+00 0.12608641982846233078D-02 0.10000000000000000000D+01
    \\
    \\  2.494800386918399765D+291  1.582914569427869018D+175  1.004336277661868922D+59  3.186183822264904554D-58  5.053968264940243633D-175  0.40083367200179455560D-291
    \\
    \\0
;

/// The next number of `dbal_in`, Fortran's `D` exponent read as `e`.
fn nextNumber(it: *std.mem.TokenIterator(u8, .any)) f64 {
    const tok = it.next().?;
    var buf: [64]u8 = undefined;
    const s = buf[0..tok.len];
    for (tok, s) |c, *d| d.* = if (c == 'D' or c == 'd') 'e' else c;
    return std.fmt.parseFloat(f64, s) catch unreachable;
}

/// `gebak` against `gebal`'s `A` and `B = D⁻¹ Pᵀ A P D`: `A X = X B` with
/// `X = P D` (right), `Yᵀ A = B Yᵀ` with `Y = P D⁻¹` (left), within
/// `4 n eps ‖A‖ ‖X‖` per entry.
fn expectBack(comptime T: type, al: Allocator, job: BalanceJob, n: usize, a: []const T, b: []const T, ilo: isize, ihi: isize, scale: []const Float(T)) !void {
    const F = Float(T);
    const x = try al.alloc(T, n * n);
    defer al.free(x);
    const ni: isize = @intCast(n);
    var amax: F = 0;
    for (a, b) |p, q| amax = @max(amax, ops.abs(T, p), ops.abs(T, q));
    for ([_]Side{ .r, .l }) |side| {
        @memset(x, ops.zero(T));
        for (0..n) |i| x[i * (n + 1)] = ops.one(T);
        try gebak(al, check.io, T, job, side, ni, ilo, ihi, scale.ptr, ni, x.ptr, ni);
        var xmax: F = 0;
        for (x) |v| xmax = @max(xmax, ops.abs(T, v));
        const tol = 4 * @as(F, @floatFromInt(n)) * std.math.floatEps(F) * amax * xmax;
        for (0..n) |i| for (0..n) |j| {
            var l = ops.zero(T);
            var r = ops.zero(T);
            for (0..n) |k| {
                if (side == .r) {
                    l = ops.add(T, l, ops.mul(T, a[i + k * n], x[k + j * n]));
                    r = ops.add(T, r, ops.mul(T, x[i + k * n], b[k + j * n]));
                } else {
                    l = ops.add(T, l, ops.mul(T, x[k + i * n], a[k + j * n]));
                    r = ops.add(T, r, ops.mul(T, b[i + k * n], x[j + k * n]));
                }
            }
            const d = ops.abs(T, ops.sub(T, l, r));
            if (!(d <= tol)) {
                check.say("gebak {s} job {s} side {s} n {d}: ({d}, {d}) off by {e}, tolerance {e}\n", .{ @typeName(T), @tagName(job), @tagName(side), n, i, j, @as(f64, @floatCast(d)), @as(f64, @floatCast(tol)) });
                return error.TestUnexpectedResult;
            }
        };
    }
}

test "lapack gebal: the reference's dbal.in cases, and gebak takes vectors back" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        var it = std.mem.tokenizeAny(u8, dbal_in, " \n\r");
        // the title line
        for (0..4) |_| _ = it.next();
        var case: usize = 1;
        while (true) : (case += 1) {
            const n: usize = @intFromFloat(nextNumber(&it));
            if (n == 0) break;
            var a: [7 * 7]f64 = undefined;
            var want: [7 * 7]f64 = undefined;
            var wscale: [7]f64 = undefined;
            for (0..n) |i| for (0..n) |j| {
                a[i + j * n] = nextNumber(&it);
            };
            const wlo: isize = @intFromFloat(nextNumber(&it));
            const whi: isize = @intFromFloat(nextNumber(&it));
            for (0..n) |i| for (0..n) |j| {
                want[i + j * n] = nextNumber(&it);
            };
            for (0..n) |i| wscale[i] = nextNumber(&it);
            // case 13's 1e120 and its scalings are past f32's range
            var fits = true;
            for (a[0 .. n * n]) |v| if (v != 0 and (@abs(v) > std.math.floatMax(F) or @abs(v) < std.math.floatMin(F))) {
                fits = false;
            };
            for (wscale[0..n]) |v| if (@abs(v) > std.math.floatMax(F) or @abs(v) < std.math.floatMin(F)) {
                fits = false;
            };
            if (!fits) continue;
            var got: [7 * 7]T = undefined;
            for (0..n * n) |i| got[i] = chol.fromF64(T, a[i]);
            const a0 = got;
            var scale: [7]F = undefined;
            var ilo: isize = 0;
            var ihi: isize = 0;
            try gebal(al, check.io, T, .b, @intCast(n), &got, @intCast(n), &ilo, &ihi, &scale);
            try testing.expectEqual(wlo, ilo);
            try testing.expectEqual(whi, ihi);
            // dchkbl.f's relative difference; the file's digits are exact
            // but for case 13's twenty
            const tol: F = @max(8 * std.math.floatEps(F), 1e-15);
            for (0..n * n) |i| {
                const g = got[i];
                const w: F = @floatCast(want[i]);
                const d = if (comptime isComplex(T)) @abs(g.re - w) + @abs(g.im) else @abs(g - w);
                if (!(d <= tol * @max(ops.abs(T, g), @abs(w), std.math.floatMin(F)))) {
                    check.say("gebal {s} dbal.in case {d}: A[{d}] {e} against {e}\n", .{ @typeName(T), case, i, @as(f64, @floatCast(ops.re(T, g))), want[i] });
                    return error.TestUnexpectedResult;
                }
            }
            for (0..n) |i| {
                const w: F = @floatCast(wscale[i]);
                if (!(@abs(scale[i] - w) <= tol * @max(scale[i], w))) {
                    check.say("gebal {s} dbal.in case {d}: scale[{d}] {e} against {e}\n", .{ @typeName(T), case, i, @as(f64, @floatCast(scale[i])), wscale[i] });
                    return error.TestUnexpectedResult;
                }
            }
            try expectBack(T, al, .b, n, a0[0 .. n * n], got[0 .. n * n], ilo, ihi, scale[0..n]);
        }
    }
}

test "lapack gebal: every job on random and structured matrices, with gebak" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        for ([_]usize{ 1, 2, 5, 12 }) |n| for ([_]BalanceJob{ .n, .p, .s, .b }) |job| for (0..3) |kind| {
            const a = try al.alloc(T, n * n);
            defer al.free(a);
            check.fill(T, a, n * 31 + kind);
            // badly scaled rows and columns
            for (0..n) |i| for (0..n) |j| {
                const e: F = @floatFromInt(@as(isize, @intCast(i)) - @as(isize, @intCast(j)));
                a[i + j * n] = ops.mulReal(T, a[i + j * n], @as(F, @floatCast(std.math.pow(f64, 4, @floatCast(e)))));
            };
            // isolated eigenvalues: a row, then a column, zero off the diagonal
            if (kind >= 1 and n > 2) for (0..n) |j| if (j != 1) {
                a[1 + j * n] = ops.zero(T);
            };
            if (kind == 2 and n > 3) for (0..n) |i| if (i != n - 2) {
                a[i + (n - 2) * n] = ops.zero(T);
            };
            const b = try al.alloc(T, n * n);
            defer al.free(b);
            @memcpy(b, a);
            const scale = try al.alloc(F, n);
            defer al.free(scale);
            var ilo: isize = 0;
            var ihi: isize = 0;
            const ni: isize = @intCast(n);
            try gebal(al, check.io, T, job, ni, b.ptr, ni, &ilo, &ihi, scale.ptr);
            try testing.expect(1 <= ilo and ilo <= ihi + 1 and ihi <= ni);
            if (job == .n) {
                try testing.expect(ilo == 1 and ihi == ni);
                try testing.expectEqualSlices(u8, std.mem.sliceAsBytes(a), std.mem.sliceAsBytes(b));
            }
            if (job == .s) try testing.expect(ilo == 1 and ihi == ni);
            // below the diagonal outside ilo..ihi: zero
            for (0..n) |j| for (j + 1..n) |i| {
                if (j + 1 < ilo or i + 1 > ihi) try testing.expect(ops.isZero(T, b[i + j * n]));
            };
            // the scalings are powers of 2
            if (job == .s or job == .b) for (@intCast(ilo - 1)..@intCast(ihi)) |i| {
                try testing.expectEqual(@as(F, 0.5), std.math.frexp(scale[i]).significand);
            };
            try expectBack(T, al, job, n, a, b, ilo, ihi, scale);
        };
    }
}

test "lapack gebal: n = 0 and 1, NaN, bad arguments; gebak's" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        var ilo: isize = 9;
        var ihi: isize = 9;
        var scale = [_]F{ 7, 7 };
        var a = [_]T{ ops.scalar(T, 2), ops.scalar(T, 1), ops.scalar(T, 1), ops.scalar(T, 3) };
        try gebal(al, check.io, T, .b, 0, &a, 1, &ilo, &ihi, &scale);
        try testing.expect(ilo == 1 and ihi == 0 and scale[0] == 7);
        try gebal(al, check.io, T, .b, 1, &a, 1, &ilo, &ihi, &scale);
        try testing.expect(ilo == 1 and ihi == 1 and scale[0] == 1);
        a[1] = chol.poisoned(T);
        try testing.expectError(error.InvalidA, gebal(al, check.io, T, .s, 2, &a, 2, &ilo, &ihi, &scale));
    }
    var a = [_]f64{1};
    var s = [_]f64{1};
    var lo: isize = 0;
    var hi: isize = 0;
    try testing.expectError(error.InvalidN, gebal(al, check.io, f64, .b, -1, &a, 1, &lo, &hi, &s));
    try testing.expectError(error.InvalidLda, gebal(al, check.io, f64, .b, 2, &a, 1, &lo, &hi, &s));
    try testing.expectError(error.InvalidN, gebak(al, check.io, f64, .b, .r, -1, 1, 0, &s, 1, &a, 1));
    try testing.expectError(error.InvalidIlo, gebak(al, check.io, f64, .b, .r, 1, 2, 1, &s, 1, &a, 1));
    try testing.expectError(error.InvalidIhi, gebak(al, check.io, f64, .b, .l, 2, 2, 1, &s, 1, &a, 2));
    try testing.expectError(error.InvalidM, gebak(al, check.io, f64, .b, .l, 1, 1, 1, &s, -1, &a, 1));
    try testing.expectError(error.InvalidLdv, gebak(al, check.io, f64, .b, .l, 2, 1, 2, &s, 1, &a, 1));
    try gebak(al, check.io, f64, .b, .l, 0, 1, 0, &s, 1, &a, 1);
}

test "faer balance: D⁻¹ Pᵀ A P D is the balanced matrix, isolated eigenvalues outside ilo..ihi" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex32 }) |T| {
        const F = Float(T);
        const n = 9;
        const a0 = try ops.alloc(T, al, n, n);
        defer ops.free(T, al, a0);
        check.fill(T, a0.ptr[0 .. n * n], 4);
        // graded rows and columns, and two isolated eigenvalues: row 2 and
        // column 6 zero off the diagonal
        for (0..n) |j| for (0..n) |i| {
            const d = @as(F, @floatCast(std.math.pow(f64, 2, @as(f64, @floatFromInt(i)) * 3 - @as(f64, @floatFromInt(j)) * 3)));
            a0.at(i, j).* = ops.mulReal(T, a0.get(i, j), d);
            if ((i == 2 and j != 2) or (j == 6 and i != 6)) a0.at(i, j).* = ops.zero(T);
        };
        const b = try ops.alloc(T, al, n, n);
        defer ops.free(T, al, b);
        b.copyFrom(a0);
        var scale: [n]F = undefined;
                const info = try balance(check.io, T, b, &scale);
        try testing.expect(info.ilo >= 1 and info.ihi <= n - 1);
        // X = P D, the identity backed: then A X = X B
        const x = try ops.alloc(T, al, n, n);
        defer ops.free(T, al, x);
        x.identity();
        try back(check.io, T, info, &scale, false, x);
        const ax = try ops.alloc(T, al, n, n);
        defer ops.free(T, al, ax);
        const xb = try ops.alloc(T, al, n, n);
        defer ops.free(T, al, xb);
        try ops.gemm(al, check.io, T, .n, .n, ops.one(T), a0, x, ops.zero(T), ax);
        try ops.gemm(al, check.io, T, .n, .n, ops.one(T), x, b, ops.zero(T), xb);
        for (0..n) |j| for (0..n) |i| {
            const err = ops.abs(T, ops.sub(T, ax.get(i, j), xb.get(i, j)));
            try testing.expect(err <= 16 * std.math.floatEps(F) * @max(ops.abs(T, ax.get(i, j)), std.math.floatMin(F)));
        };
        // outside the window: triangular
        for (0..n) |j| for (0..n) |i| {
            if (i > j and (j < info.ilo or i >= info.ihi)) try testing.expect(ops.isZero(T, b.get(i, j)));
        };
    }
}
