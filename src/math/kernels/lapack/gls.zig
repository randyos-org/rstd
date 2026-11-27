//! Generalized least squares, `?GGLSE` and `?GGGLM`: Reference-LAPACK
//! 3.12.0 (`lapack-netlib/SRC/dgglse.f`, `zgglse.f`, `dggglm.f`,
//! `zggglm.f`, with their power-of-two scaling, OpenBLAS v0.3.34-794;
//! modified BSD), over `ggrqf` and `ggqrf`.

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const interface = @import("interface.zig");
const ops = @import("ops.zig");
const blas_calls = @import("blas_calls.zig");
const auxiliary = @import("auxiliary.zig");
const norm = @import("norm.zig");
const qr = @import("qr.zig");
const rq = @import("rq.zig");
const gqr = @import("gqr.zig");
const triangular = @import("triangular.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = ops.Mat;
const Float = types.Float;

const Common = interface.ShapeError || Allocator.Error || Io.Cancelable;

/// What `gglse` returns: `TooManyConstraints` for `p > n`, `TooFewRows`
/// for `m + p < n`; `BRankDeficient` when `rank(B) < p` and
/// `ABRankDeficient` when `rank([A; B]) < n` (the reference's `INFO` 1
/// and 2); a slice too short (`a`, `b`, `c` under `m`, `d` under `p`, `x`
/// under `n`) or a leading dimension too small; and the allocator's and
/// cancellation.
pub const GglseError = error{ TooManyConstraints, TooFewRows, BRankDeficient, ABRankDeficient } || Common;

/// What `ggglm` returns: `TooManyColumns` for `m > n`, `TooFewColumns` for
/// `m + p < n`; `ABRankDeficient` when `rank([A B]) < n` and
/// `ARankDeficient` when `rank(A) < m` (the reference's `INFO` 1 and 2); a
/// slice too short (`a`, `b`, `d` under `n`, `x` under `m`, `y` under `p`)
/// or a leading dimension too small; and the allocator's and cancellation.
pub const GgglmError = error{ TooManyColumns, TooFewColumns, ABRankDeficient, ARankDeficient } || Common;

/// `?GGLSE`: `x` (`n`) minimizing `‖c - A x‖` subject to `B x = d`, for
/// the `m x n` `A` and `p x n` `B`, `p <= n <= m + p`, by `ggrqf`. `A`,
/// `B` and `d` are destroyed; `c[n - p..m]` gets a residual whose norm is
/// the minimum `‖c - A x‖`.
pub fn gglse(allocator: Allocator, io: Io, comptime T: type, m: usize, n: usize, p: usize, a: []T, lda: usize, b: []T, ldb: usize, c: []T, d: []T, x: []T) GglseError!void {
    comptime ops.requireFloat(T, "gglse");
    const F = Float(T);
    if (p > n) return error.TooManyConstraints;
    if (m + p < n) return error.TooFewRows;
    try util.holds(a.len, m, n, lda);
    try util.holds(b.len, p, n, ldb);
    if (c.len < m or d.len < p or x.len < n) return error.SliceTooShort;
    if (n == 0) return;
    const am: Mat(T) = .{ .ptr = a.ptr, .rows = m, .cols = n, .ld = lda };
    const bm: Mat(T) = .{ .ptr = b.ptr, .rows = p, .cols = n, .ld = ldb };
    const cv = c[0..m];
    const dv = d[0..p];
    const xv = x[0..n];
    const r: Range(F) = .{};
    const ka = try r.scaleMat(allocator, T, am);
    const kb = try r.scaleMat(allocator, T, bm);
    // c carries 2^(ka + kt), d 2^(kb + kt), x 2^kt; kt brings the larger
    // right-hand side into range
    const ic = r.exponent(maxAbs(T, cv), ka);
    const id = r.exponent(maxAbs(T, dv), kb);
    var kt: i32 = 0;
    if (ic != null or id != null) {
        const iu = @max(ic orelse std.math.minInt(i32), id orelse std.math.minInt(i32));
        if (iu > r.ibig) kt = r.ibig - iu else if (iu < r.isml) kt = r.isml - iu;
    }
    scaleVec(T, cv, ka + kt);
    scaleVec(T, dv, kb + kt);
    // B Qᴴ = (0 T12), Zᴴ A Qᴴ = (R11 R12; 0 R22)
    const mn = @min(m, n);
    const taub = try allocator.alloc(T, p);
    defer allocator.free(taub);
    const taua = try allocator.alloc(T, mn);
    defer allocator.free(taua);
    try gqr.ggrqf(allocator, io, T, p, m, n, b, ldb, taub, a, lda, taua);
    // c := Zᴴ c
    if (mn > 0) qr.ormqr(allocator, io, T, .l, .c, m, 1, mn, a, lda, taua, cv, @max(m, 1)) catch |e| return ops.narrow(GglseError, e);
    if (p > 0) {
        // T12 x2 = d
        solveUpper(allocator, io, T, bm.sub(0, n - p, p, p), dv) catch |e| return mapSingular(GglseError, error.BRankDeficient, e);
        @memcpy(xv[n - p ..], dv);
        // c1 -= R12 x2
        try blas_calls.gemv(allocator, io, T, .n, ops.scalar(T, -1), am.sub(0, n - p, n - p, p), dv.ptr, 1, ops.one(T), cv.ptr, 1);
    }
    if (n > p) {
        // R11 x1 = c1
        solveUpper(allocator, io, T, am.sub(0, 0, n - p, n - p), cv[0 .. n - p]) catch |e| return mapSingular(GglseError, error.ABRankDeficient, e);
        @memcpy(xv[0 .. n - p], cv[0 .. n - p]);
    }
    // the residual: c2 -= R22 x2
    var nr = p;
    if (m < n) {
        nr = m + p - n;
        if (nr > 0) try blas_calls.gemv(allocator, io, T, .n, ops.scalar(T, -1), am.sub(n - p, m, nr, n - m), dv[nr..].ptr, 1, ops.one(T), cv[n - p ..].ptr, 1);
    }
    if (nr > 0) {
        upperTimes(T, am.sub(n - p, n - p, nr, nr), dv[0..nr]);
        for (cv[n - p ..][0..nr], dv[0..nr]) |*ci, di| ci.* = ops.sub(T, ci.*, di);
    }
    // x := Qᴴ x
    if (p > 0) rq.ormrq(allocator, io, T, .l, .c, n, 1, p, b, ldb, taub, xv, n) catch |e| return ops.narrow(GglseError, e);
    scaleVec(T, xv, -kt);
    if (m + p > n) scaleVec(T, cv[n - p ..], -(ka + kt));
}

/// `?GGGLM`: `x` (`m`) and `y` (`p`) minimizing `‖y‖` subject to `d = A x
/// + B y`, for the `n x m` `A` and `n x p` `B`, `m <= n <= m + p`, by
/// `ggqrf`. `A`, `B` and `d` are destroyed.
pub fn ggglm(allocator: Allocator, io: Io, comptime T: type, n: usize, m: usize, p: usize, a: []T, lda: usize, b: []T, ldb: usize, d: []T, x: []T, y: []T) GgglmError!void {
    comptime ops.requireFloat(T, "ggglm");
    const F = Float(T);
    if (m > n) return error.TooManyColumns;
    if (m + p < n) return error.TooFewColumns;
    try util.holds(a.len, n, m, lda);
    try util.holds(b.len, n, p, ldb);
    if (d.len < n or x.len < m or y.len < p) return error.SliceTooShort;
    const xv = x[0..m];
    const yv = y[0..p];
    if (n == 0) {
        @memset(xv, ops.zero(T));
        @memset(yv, ops.zero(T));
        return;
    }
    const am: Mat(T) = .{ .ptr = a.ptr, .rows = n, .cols = m, .ld = lda };
    const bm: Mat(T) = .{ .ptr = b.ptr, .rows = n, .cols = p, .ld = ldb };
    const dv = d[0..n];
    // x then carries 2^(kd - ka), y 2^(kd - kb)
    const r: Range(F) = .{};
    const ka = try r.scaleMat(allocator, T, am);
    const kb = try r.scaleMat(allocator, T, bm);
    const kd = r.shift(maxAbs(T, dv));
    scaleVec(T, dv, kd);
    // Qᴴ A = (R11; 0), Qᴴ B Zᴴ = (T11 T12; 0 T22)
    const np = @min(n, p);
    const taua = try allocator.alloc(T, m);
    defer allocator.free(taua);
    const taub = try allocator.alloc(T, np);
    defer allocator.free(taub);
    try gqr.ggqrf(allocator, io, T, n, m, p, a, lda, taua, b, ldb, taub);
    // d := Qᴴ d
    if (m > 0) qr.ormqr(allocator, io, T, .l, .c, n, 1, m, a, lda, taua, dv, n) catch |e| return ops.narrow(GgglmError, e);
    const y1 = m + p - n;
    if (n > m) {
        // T22 y2 = d2
        solveUpper(allocator, io, T, bm.sub(m, y1, n - m, n - m), dv[m..]) catch |e| return mapSingular(GgglmError, error.ABRankDeficient, e);
        @memcpy(yv[y1..], dv[m..]);
    }
    @memset(yv[0..y1], ops.zero(T));
    // d1 -= T12 y2
    if (n > m and m > 0) try blas_calls.gemv(allocator, io, T, .n, ops.scalar(T, -1), bm.sub(0, y1, m, n - m), yv[y1..].ptr, 1, ops.one(T), dv.ptr, 1);
    if (m > 0) {
        // R11 x = d1
        solveUpper(allocator, io, T, am.sub(0, 0, m, m), dv[0..m]) catch |e| return mapSingular(GgglmError, error.ARankDeficient, e);
        @memcpy(xv, dv[0..m]);
    }
    // y := Zᴴ y, Z's reflections in B's last np rows
    if (np > 0) rq.ormrq(allocator, io, T, .l, .c, p, 1, np, b[n - np ..], ldb, taub, yv, p) catch |e| return ops.narrow(GgglmError, e);
    scaleVec(T, xv, ka - kd);
    scaleVec(T, yv, kb - kd);
}

/// The drivers' range `[smlnum, bignum)`, `smlnum = DLAMCH('S') /
/// DLAMCH('P')`, as the exponents `[isml, ibig]` (Fortran's `EXPONENT`,
/// `frexp`'s).
fn Range(comptime F: type) type {
    return struct {
        const Self = @This();

        isml: i32 = std.math.frexp(std.math.floatMin(F) / std.math.floatEps(F)).exponent,
        ibig: i32 = std.math.frexp(std.math.floatEps(F) / std.math.floatMin(F)).exponent - 1,

        /// `nrm`'s exponent plus `k`; null for zero, infinity or NaN.
        fn exponent(r: Self, nrm: F, k: i32) ?i32 {
            _ = r;
            if (!(nrm > 0 and nrm <= std.math.floatMax(F))) return null;
            return std.math.frexp(nrm).exponent + k;
        }

        /// The power of two bringing `nrm` into range.
        fn shift(r: Self, nrm: F) i32 {
            const e = r.exponent(nrm, 0) orelse return 0;
            if (e < r.isml) return r.isml - e;
            if (e > r.ibig) return r.ibig - e;
            return 0;
        }

        /// `a` scaled by the shift of its largest entry, returned.
        fn scaleMat(r: Self, al: Allocator, comptime T: type, a: Mat(T)) Allocator.Error!i32 {
            if (a.rows == 0 or a.cols == 0) return 0;
            const k = r.shift(try norm.matNorm(al, T, .max_abs, a.rows, a.cols, a.ptr, a.ld));
            // `zlascl` with both ends at or above one, as the reference
            if (k != 0) auxiliary.lascl(T, false, std.math.ldexp(@as(F, 1), @max(-k, 0)), std.math.ldexp(@as(F, 1), @max(k, 0)), a);
            return k;
        }
    };
}

/// `zlange('M')` of a vector: NaN if an entry is.
fn maxAbs(comptime T: type, v: []const T) Float(T) {
    var big: Float(T) = 0;
    for (v) |e| {
        const s = ops.abs(T, e);
        if (std.math.isNan(s)) return s;
        big = @max(big, s);
    }
    return big;
}

/// `v *= 2^k`, exact.
fn scaleVec(comptime T: type, v: []T, k: i32) void {
    if (k == 0) return;
    const s = std.math.ldexp(@as(Float(T), 1), k);
    for (v) |*e| e.* = ops.mulReal(T, e.*, s);
}

/// `U z = v` for the upper triangular, non-unit `U`, in place.
fn solveUpper(al: Allocator, io: Io, comptime T: type, u: Mat(T), v: []T) triangular.TrtrsError!void {
    const k = u.rows;
    return triangular.trtrs(al, io, T, .u, .n, .n, k, 1, u.ptr[0 .. (k - 1) * u.ld + k], u.ld, v[0..k], k);
}

/// `trtrs`'s `Singular` as `as`, its others narrowed to `E`.
fn mapSingular(comptime E: type, as: E, e: triangular.TrtrsError) E {
    return if (e == error.Singular) as else ops.narrow(E, e);
}

/// `v := U v` for the upper triangular, non-unit `U` (`ztrmv.f`'s
/// `Upper`, `No transpose`).
fn upperTimes(comptime T: type, u: Mat(T), v: []T) void {
    for (0..v.len) |i| {
        var s = ops.zero(T);
        for (i..v.len) |j| s = ops.add(T, s, ops.mul(T, u.get(i, j), v[j]));
        v[i] = s;
    }
}
