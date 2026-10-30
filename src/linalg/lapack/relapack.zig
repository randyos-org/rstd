//! ReLAPACK's recursive algorithms (Elmar Peise, MIT), from OpenBLAS's copy,
//! `relapack/src/` (OpenBLAS v0.3.34-317), which its `BUILD_RELAPACK` option
//! builds. Each routine splits its problem in two by `split` (`?REC_SPLIT`),
//! recurses on the halves and joins them with level-3 BLAS, down to
//! `crossover` (`CROSSOVER`, `relapack_config.h`), where LAPACK's unblocked
//! kernel runs.
//!
//! `potrf`, `getrf`, `trtri` and `lauum` replace OpenBLAS's only under
//! `-Drelapack` (`__root__.zig`), as OpenBLAS's `RELAPACK_REPLACE`; their
//! base cases are rstd's `potf2`, `getf2`, `trti2` and `lauu2`. The rest
//! (`sytrf`, `hetrf` and their `_rook` forms, `sygst`, `hegst`, `pbtrf`,
//! `gbtrf`, `trsyl`, `tgsyl`) OpenBLAS takes from LAPACK, so ReLAPACK's are
//! the default. Their unblocked kernels are Reference-LAPACK's, as OpenBLAS
//! ships them in `lapack-netlib/SRC/` (faer has none with LAPACK's storage):
//! `indefinite.zig`, `sylvester.zig` and the `?2` loops here.
//!
//! The BLAS calls are the level-3 drivers with the interface's worker rules
//! (`blas3/interface.zig`), so they thread as a ReLAPACK build's public BLAS
//! calls do, with the same bits at any worker count; Accelerate never runs
//! them. Workspace ReLAPACK takes from `malloc` comes from the allocator.
//!
//! Where ReLAPACK and the reference differ, the reference wins, and the
//! routine says so.

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const threading = @import("../blas/threading.zig");
const driver = @import("../blas/blas3/driver.zig");
const threaded = @import("../blas/blas3/threaded.zig");
const laswp = @import("../blas/ext/laswp/kernel.zig");
const chol = @import("cholesky.zig");
const lu = @import("lu.zig");
const tri = @import("triangular.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Uplo = types.Uplo;
const Diag = types.Diag;
const Trans = types.Trans;
const Side = types.Side;
const Float = types.Float;
const isComplex = types.isComplex;
const scalar = chol.scalar;

const Error = driver.Error;

/// `CROSSOVER` (`relapack_config.h`): every routine's unblocked size.
pub const crossover = 24;

/// `SREC_SPLIT` and the rest (`relapack.h`): half of `n`, from `2 * u` on
/// rounded so a half is a multiple of `u`, a 64-byte cache line of `T`
/// (16 for `s`, 8 for `d` and `c`, 4 for `z`).
pub fn split(comptime T: type, n: usize) usize {
    const u = @max(1, 64 / @sizeOf(T));
    return if (n >= 2 * u) (n + u) / (2 * u) * u else n / 2;
}

// ---- the BLAS, as ReLAPACK calls it ------------------------------------------------

/// `?GEMM`, already checked; nothing for an empty `C`.
pub fn gemm(al: Allocator, io: Io, comptime T: type, ta: Trans, tb: Trans, m: usize, n: usize, k: usize, alpha: T, a: [*]const T, lda: usize, b: [*]const T, ldb: usize, beta: T, c: [*]T, ldc: usize) Error!void {
    if (m == 0 or n == 0) return;
    const conj: driver.Conj = .of(isComplex(T) and ta.conjA(), isComplex(T) and tb.conjA());
    const w = threaded.workers(.gemm, T, threading.product(&.{ m, n, k }));
    try driver.gemm(al, io, T, chol.table(T), ta.transposed(), tb.transposed(), conj, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, w);
}

/// `?TRSM` (`solve`) or `?TRMM`, already checked.
pub fn trsm(al: Allocator, io: Io, comptime T: type, comptime solve: bool, side: Side, uplo: Uplo, ta: Trans, diag: Diag, m: usize, n: usize, alpha: T, a: [*]const T, lda: usize, b: [*]T, ldb: usize) Error!void {
    if (m == 0 or n == 0) return;
    const f = if (solve) driver.trsm else driver.trmm;
    const w = threaded.workers(.trsm, T, threading.product(&.{ m, n }));
    try f(al, io, T, chol.table(T), side == .l, uplo == .u, ta.transposed(), isComplex(T) and ta.conjA(), diag == .u, m, n, alpha, a, lda, b, ldb, w);
}

/// `?SYRK`, or with `herm` `?HERK` (`alpha`, `beta` real); `?SYR2K` and
/// `?HER2K` with `b`. `trans` is `t` or `c` for `op(A) = A'`.
pub fn syrk(al: Allocator, io: Io, comptime T: type, comptime herm: bool, uplo: Uplo, trans: Trans, n: usize, k: usize, alpha: T, a: [*]const T, lda: usize, b: ?[*]const T, ldb: usize, beta: T, c: [*]T, ldc: usize) Error!void {
    if (n == 0) return;
    const w = if (b == null)
        threaded.workers(.syrk, T, threading.product(&.{ n + 1, n, k }))
    else
        threaded.workers(.syr2k, T, threading.product(&.{ n, k }));
    try driver.syrk(al, io, T, chol.table(T), herm and isComplex(T), uplo == .u, trans != .n, n, k, alpha, a, lda, b, ldb, beta, c, ldc, w);
}

/// `?SYMM`, or with `herm` `?HEMM`.
pub fn symm(al: Allocator, io: Io, comptime T: type, comptime herm: bool, side: Side, uplo: Uplo, m: usize, n: usize, alpha: T, a: [*]const T, lda: usize, b: [*]const T, ldb: usize, beta: T, c: [*]T, ldc: usize) Error!void {
    if (m == 0 or n == 0) return;
    const w = threaded.workers(.symm, T, threading.product(&.{ 2, m, m, n }));
    try driver.symm(al, io, T, chol.table(T), herm and isComplex(T), side == .l, uplo == .u, m, n, alpha, a, lda, b, ldb, beta, c, ldc, w);
}

/// The Hermitian forms (`herk`, `C`) for complex `T`, the symmetric ones
/// (`syrk`, `T`) for real.
fn ctrans(comptime T: type) Trans {
    return if (comptime isComplex(T)) .c else .t;
}

// ---- potrf (`-Drelapack`) ------------------------------------------------------------

/// `?POTRF` by ReLAPACK's recursion (`relapack/src/zpotrf.c`, `dpotrf.c`
/// and the rest; OpenBLAS v0.3.34-317), as `cholesky.potrf` takes it, with
/// its errors and `minor`. The base case is `potf2`.
pub fn potrf(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, a: [*]T, lda: isize, minor: ?*usize) chol.PotrfError!void {
    comptime chol.floatType(T, "potrf");
    if (!try chol.checkSquare(T, uplo, n, lda)) return;
    const k = try potrfRec(allocator, io, T, uplo, @intCast(n), a, @intCast(lda));
    if (k == 0) return;
    if (minor) |m| m.* = k;
    return error.NotPositiveDefinite;
}

/// `RELAPACK_?potrf_rec`: the first minor that is not positive definite, or 0.
pub fn potrfRec(al: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, a: [*]T, lda: usize) Error!usize {
    if (n <= crossover) return chol.unblocked(T, uplo == .u, n, a, lda);
    const one = scalar(T, 1);
    const minus = scalar(T, -1);
    const n1 = split(T, n);
    const n2 = n - n1;
    const br = a + n1 * lda + n1;
    const k = try potrfRec(al, io, T, uplo, n1, a, lda);
    if (k != 0) return k;
    if (uplo == .l) {
        // A_BL = A_BL / A_TL'; A_BR -= A_BL A_BL'
        try trsm(al, io, T, true, .r, .l, ctrans(T), .n, n2, n1, one, a, lda, a + n1, lda);
        try syrk(al, io, T, true, .l, .n, n2, n1, minus, a + n1, lda, null, 0, one, br, lda);
    } else {
        // A_TR = A_TL' \ A_TR; A_BR -= A_TR' A_TR
        try trsm(al, io, T, true, .l, .u, ctrans(T), .n, n1, n2, one, a, lda, a + n1 * lda, lda);
        try syrk(al, io, T, true, .u, ctrans(T), n2, n1, minus, a + n1 * lda, lda, null, 0, one, br, lda);
    }
    const k2 = try potrfRec(al, io, T, uplo, n2, br, lda);
    return if (k2 != 0) k2 + n1 else 0;
}

// ---- getrf (`-Drelapack`) -----------------------------------------------------------

/// `?GETRF` by ReLAPACK's recursion (`relapack/src/zgetrf.c`, `dgetrf.c`
/// and the rest; OpenBLAS v0.3.34-317), as `lu.getrf` takes it, with its
/// errors and `pivot`. The base case is `lu.getf2` (ReLAPACK calls
/// `xGETRF2`, whose pivots and semantics it shares).
///
/// The reference wins: a zero pivot in the left half does not stop the
/// recursion, which ReLAPACK's `if (*info) return` does, so the
/// factorization completes, as `xGETRF`'s.
pub fn getrf(allocator: Allocator, io: Io, comptime T: type, m: isize, n: isize, a: [*]T, lda: isize, ipiv: [*]isize, pivot: ?*usize) lu.GetrfError!void {
    comptime chol.floatType(T, "getrf");
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, m)) return error.InvalidLda;
    if (m == 0 or n == 0) return;
    if (!util.fitsMatrix(m, n, lda, @sizeOf(T))) return error.InvalidLda;
    const um: usize = @intCast(m);
    const un: usize = @intCast(n);
    const la: usize = @intCast(lda);
    const info = try getrfRec(allocator, io, T, um, @min(um, un), a, la, ipiv);
    // the right remainder: `A_R = A_L \ P A_R`
    if (um < un) {
        const ar = a + la * um;
        laswp.swapRows(T, un - um, ar, la, 1, um, ipiv, 1);
        try trsm(allocator, io, T, true, .l, .l, .n, .u, um, un - um, scalar(T, 1), a, la, ar, la);
    }
    return lu.singular(info, pivot);
}

/// `RELAPACK_?getrf_rec` on `m x n`, `m >= n`: the first zero pivot, or 0.
fn getrfRec(al: Allocator, io: Io, comptime T: type, m: usize, n: usize, a: [*]T, lda: usize, ipiv: [*]isize) Error!usize {
    if (m == 0 or n == 0) return 0;
    if (n <= crossover) return lu.getf2(T, chol.table(T), m, n, a, lda, ipiv);
    const n1 = split(T, n);
    const n2 = n - n1;
    const m2 = m - n1;
    const ar = a + lda * n1;
    const bl = a + n1;
    const br = ar + n1;
    const info = try getrfRec(al, io, T, m, n1, a, lda, ipiv);
    laswp.swapRows(T, n2, ar, lda, 1, n1, ipiv, 1);
    // A_TR = A_TL \ A_TR; A_BR -= A_BL A_TR
    try trsm(al, io, T, true, .l, .l, .n, .u, n1, n2, scalar(T, 1), a, lda, ar, lda);
    try gemm(al, io, T, .n, .n, m2, n2, n1, scalar(T, -1), bl, lda, ar, lda, scalar(T, 1), br, lda);
    const info2 = try getrfRec(al, io, T, m2, n2, br, lda, ipiv + n1);
    laswp.swapRows(T, n1, bl, lda, 1, n2, ipiv + n1, 1);
    for (ipiv[n1..n]) |*p| p.* += @intCast(n1);
    return if (info != 0) info else if (info2 != 0) info2 + n1 else 0;
}

// ---- trtri, lauum (`-Drelapack`) -------------------------------------------------------

/// `?TRTRI` by ReLAPACK's recursion (`relapack/src/ztrtri.c`, `dtrtri.c`
/// and the rest; OpenBLAS v0.3.34-317), as `triangular.trtri` takes it,
/// with its errors and `pivot`. The base case is `trti2`.
///
/// The reference wins: a zero `A(i, i)` reports `i` 1-based, where
/// ReLAPACK's singularity check sets `INFO = i - 1`, so a zero first
/// diagonal element reports success.
pub fn trtri(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, diag: Diag, n: isize, a: [*]T, lda: isize, pivot: ?*usize) tri.TrtriError!void {
    comptime chol.floatType(T, "trtri");
    if (!try chol.checkSquare(T, uplo, n, lda)) return;
    const un: usize = @intCast(n);
    const la: usize = @intCast(lda);
    if (diag == .n) if (tri.firstZero(T, un, a, la)) |k| {
        if (pivot) |p| p.* = k;
        return error.Singular;
    };
    try trtriRec(allocator, io, T, uplo, diag, un, a, la);
}

fn trtriRec(al: Allocator, io: Io, comptime T: type, uplo: Uplo, diag: Diag, n: usize, a: [*]T, lda: usize) Error!void {
    if (n <= crossover) return tri.unblockedInverse(T, uplo == .u, diag == .u, n, a, lda);
    const one = scalar(T, 1);
    const minus = scalar(T, -1);
    const n1 = split(T, n);
    const n2 = n - n1;
    const br = a + n1 * lda + n1;
    try trtriRec(al, io, T, uplo, diag, n1, a, lda);
    if (uplo == .l) {
        // A_BL = -A_BL A_TL; A_BL = A_BR \ A_BL
        try trsm(al, io, T, false, .r, .l, .n, diag, n2, n1, minus, a, lda, a + n1, lda);
        try trsm(al, io, T, true, .l, .l, .n, diag, n2, n1, one, br, lda, a + n1, lda);
    } else {
        // A_TR = -A_TL A_TR; A_TR = A_TR / A_BR
        try trsm(al, io, T, false, .l, .u, .n, diag, n1, n2, minus, a, lda, a + n1 * lda, lda);
        try trsm(al, io, T, true, .r, .u, .n, diag, n1, n2, one, br, lda, a + n1 * lda, lda);
    }
    try trtriRec(al, io, T, uplo, diag, n2, br, lda);
}

/// `?LAUUM` by ReLAPACK's recursion (`relapack/src/zlauum.c`, `dlauum.c`
/// and the rest; OpenBLAS v0.3.34-317), as `triangular.lauum` takes it. The
/// base case is `lauu2`.
pub fn lauum(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, a: [*]T, lda: isize) tri.LauumError!void {
    comptime chol.floatType(T, "lauum");
    if (!try chol.checkSquare(T, uplo, n, lda)) return;
    try lauumRec(allocator, io, T, uplo, @intCast(n), a, @intCast(lda));
}

fn lauumRec(al: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, a: [*]T, lda: usize) Error!void {
    if (n <= crossover) return tri.unblockedProduct(T, uplo == .u, n, a, lda);
    const one = scalar(T, 1);
    const n1 = split(T, n);
    const n2 = n - n1;
    const br = a + n1 * lda + n1;
    try lauumRec(al, io, T, uplo, n1, a, lda);
    if (uplo == .l) {
        // A_TL += A_BL' A_BL; A_BL = A_BR' A_BL
        try syrk(al, io, T, true, .l, ctrans(T), n1, n2, one, a + n1, lda, null, 0, one, a, lda);
        try trsm(al, io, T, false, .l, .l, ctrans(T), .n, n2, n1, one, br, lda, a + n1, lda);
    } else {
        // A_TL += A_TR A_TR'; A_TR = A_TR A_BR'
        try syrk(al, io, T, true, .u, .n, n1, n2, one, a + n1 * lda, lda, null, 0, one, a, lda);
        try trsm(al, io, T, false, .r, .u, ctrans(T), .n, n1, n2, one, br, lda, a + n1 * lda, lda);
    }
    try lauumRec(al, io, T, uplo, n2, br, lda);
}

// ---- tests ---------------------------------------------------------------------

const testing = std.testing;
const blas = @import("../blas/__root__.zig");
const check = @import("../blas/_test/blas1/check.zig");
const Vtable = @import("../blas/Vtable.zig");
const options = @import("lapack_options");

pub const test_types = chol.test_types;
const isNanT = chol.isNanT;

/// Around the crossover and `split`'s rounding for every element size.
pub const test_sizes = [_]usize{ 1, 2, 24, 25, 33, 49, 70, 100 };

pub fn uploOf(upper: bool) Uplo {
    return if (upper) .u else .l;
}

fn potrfCases(vt: *const Vtable) anyerror!void {
    const al = testing.allocator;
    const saved = blas.testUse(vt);
    defer blas.testRestore(saved);
    inline for (test_types) |T| for ([_]bool{ true, false }) |upper| for (test_sizes) |n| {
        if (check.rig() and n > 70) continue;
        var p = try chol.Problem(T).init(al, n, upper, n * 37 + @intFromBool(upper));
        defer p.deinit(al);
        try potrf(al, check.io, T, uploOf(upper), @intCast(n), p.a.ptr, @intCast(p.lda), null);
        p.expectFactor() catch |e| {
            check.say("relapack potrf {s} n={d}\n", .{ @typeName(T), n });
            return e;
        };
    };
}

test "relapack potrf: A = L Lᴴ from a known L, both triangles, on every linked core" {
    try check.eachCore(potrfCases);
}

test "relapack potrf: the first minor that is not positive definite, either half" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex32 }) |T| for ([_]bool{ true, false }) |upper| {
        const n = 49;
        var p = try chol.Problem(T).init(al, n, upper, 5);
        defer p.deinit(al);
        for ([_]usize{ 0, 23, 24, 25, 40, 48 }) |k| {
            p.store();
            p.a[k * (p.lda + 1)] = scalar(T, -100);
            var minor: usize = 0;
            try testing.expectError(error.NotPositiveDefinite, potrf(al, check.io, T, uploOf(upper), n, p.a.ptr, @intCast(p.lda), &minor));
            try testing.expectEqual(k + 1, minor);
        }
    };
}

const getrf_shapes = [_][2]usize{ .{ 1, 1 }, .{ 24, 24 }, .{ 25, 25 }, .{ 50, 50 }, .{ 70, 33 }, .{ 33, 70 }, .{ 100, 60 }, .{ 26, 100 } };

fn getrfCases(vt: *const Vtable) anyerror!void {
    const al = testing.allocator;
    const saved = blas.testUse(vt);
    defer blas.testRestore(saved);
    inline for (test_types) |T| for (getrf_shapes) |s| {
        if (check.rig() and s[0] * s[1] > 3000) continue;
        var p = try lu.Problem(T).init(al, s[0], s[1], s[0] * 101 + s[1]);
        defer p.deinit(al);
        try getrf(al, check.io, T, @intCast(s[0]), @intCast(s[1]), p.a.ptr, @intCast(p.lda), p.ipiv.ptr, null);
        p.expectFactor(al) catch |e| {
            check.say("relapack getrf {s}\n", .{@typeName(T)});
            return e;
        };
    };
}

test "relapack getrf: P L U = A, square, tall and wide, on every linked core" {
    try check.eachCore(getrfCases);
}

test "relapack getrf: a zero column gives that pivot, the factorization completed, either half" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64 }) |T| for ([_]usize{ 0, 10, 30, 45 }) |k| {
        var p = try lu.Problem(T).init(al, 50, 50, 9);
        defer p.deinit(al);
        for (0..50) |i| p.aw[i + k * 50] = check.wide(T, scalar(T, 0));
        p.store();
        var piv: usize = 0;
        try testing.expectError(error.Singular, getrf(al, check.io, T, 50, 50, p.a.ptr, @intCast(p.lda), p.ipiv.ptr, &piv));
        try testing.expectEqual(k + 1, piv);
        try p.expectFactor(al);
    };
}

fn triCases(vt: *const Vtable) anyerror!void {
    const al = testing.allocator;
    const saved = blas.testUse(vt);
    defer blas.testRestore(saved);
    inline for (test_types) |T| for ([_]bool{ true, false }) |upper| for ([_]bool{ false, true }) |unit| for (test_sizes) |n| {
        if (check.rig() and n > 70) continue;
        var t = try tri.Tri(T).init(al, n, upper, unit, n * 13 + @intFromBool(upper));
        defer t.deinit(al);
        try trtri(al, check.io, T, uploOf(upper), if (unit) .u else .n, @intCast(n), t.a.ptr, @intCast(t.lda), null);
        try t.expectInverse();
        if (unit) continue;
        t.store();
        try lauum(al, check.io, T, uploOf(upper), @intCast(n), t.a.ptr, @intCast(t.lda));
        try t.expectProduct();
    };
}

test "relapack trtri and lauum: A X = I and U Uᴴ, Lᴴ L, both triangles, on every linked core" {
    try check.eachCore(triCases);
}

test "relapack trtri: a zero first diagonal element is singular at 1, nothing written" {
    const al = testing.allocator;
    var t = try tri.Tri(f64).init(al, 30, true, false, 3);
    defer t.deinit(al);
    t.a[0] = 0;
    const before = try al.dupe(f64, t.a);
    defer al.free(before);
    var piv: usize = 0;
    try testing.expectError(error.Singular, trtri(al, check.io, f64, .u, .n, 30, t.a.ptr, @intCast(t.lda), &piv));
    try testing.expectEqual(1, piv);
    for (before, t.a) |x, y| try testing.expect(chol.sameBits(f64, x, y));
}

test "relapack potrf, getrf, trtri, lauum: the first bad argument, nothing read when empty" {
    const al = testing.allocator;
    var a: [4]f64 = @splat(chol.poisoned(f64));
    var ipiv: [2]isize = @splat(-1);
    try testing.expectError(error.InvalidUplo, potrf(al, check.io, f64, .v, 2, &a, 2, null));
    try testing.expectError(error.InvalidN, potrf(al, check.io, f64, .u, -1, &a, 2, null));
    try testing.expectError(error.InvalidLda, potrf(al, check.io, f64, .u, 2, &a, 1, null));
    try testing.expectError(error.InvalidM, getrf(al, check.io, f64, -1, 2, &a, 2, &ipiv, null));
    try testing.expectError(error.InvalidN, getrf(al, check.io, f64, 2, -1, &a, 2, &ipiv, null));
    try testing.expectError(error.InvalidLda, getrf(al, check.io, f64, 2, 2, &a, 1, &ipiv, null));
    try testing.expectError(error.InvalidUplo, trtri(al, check.io, f64, .m, .n, 2, &a, 2, null));
    try testing.expectError(error.InvalidLda, lauum(al, check.io, f64, .l, 2, &a, 1));
    try potrf(al, check.io, f64, .u, 0, &a, 1, null);
    try getrf(al, check.io, f64, 0, 2, &a, 1, &ipiv, null);
    try trtri(al, check.io, f64, .l, .n, 0, &a, 1, null);
    try lauum(al, check.io, f64, .l, 0, &a, 1);
    for (a) |x| try testing.expect(std.math.isNan(x));
    try testing.expectEqual(-1, ipiv[0]);
}

/// `x` and `y` within `64 n` ulps of `y`'s largest element, NaN where both
/// are.
pub fn expectAgree(comptime T: type, what: []const u8, n: usize, x: []const T, y: []const T) !void {
    var big: f128 = 0;
    for (y) |v| if (!isNanT(T, v)) {
        big = @max(big, check.abs1(check.wide(T, v)));
    };
    const tol = 64 * @as(f128, @floatFromInt(n)) * std.math.floatEps(Float(T)) * big;
    for (x, y, 0..) |u, v, i| {
        if (isNanT(T, u) and isNanT(T, v)) continue;
        const d = check.abs1(check.subW(check.wide(T, u), check.wide(T, v)));
        if (!(d <= tol)) {
            check.say("{s} {s} n={d}: element {d} differs by {e}, allowed {e}\n", .{ what, @typeName(T), n, i, @as(f64, @floatCast(d)), @as(f64, @floatCast(tol)) });
            return error.TestUnexpectedResult;
        }
    }
}

test "with -Drelapack, ReLAPACK's potrf, getrf, trtri and lauum agree with OpenBLAS's" {
    if (!options.relapack) return error.SkipZigTest;
    const root = @import("__root__.zig");
    const al = testing.allocator;
    const io = check.io;
    inline for (.{ f32, f64, types.Complex32, types.Complex64 }) |T| for ([_]bool{ true, false }) |upper| {
        const n = 100;
        const uplo = uploOf(upper);
        var p = try chol.Problem(T).init(al, n, upper, 77);
        defer p.deinit(al);
        const x = try al.dupe(T, p.a);
        defer al.free(x);
        try root.potrf(al, io, T, uplo, n, x.ptr, @intCast(p.lda), null);
        try chol.potrf(al, io, T, uplo, n, p.a.ptr, @intCast(p.lda), null);
        try expectAgree(T, "potrf", n, x, p.a);

        var t = try tri.Tri(T).init(al, n, upper, false, 78);
        defer t.deinit(al);
        const z = try al.dupe(T, t.a);
        defer al.free(z);
        try root.trtri(al, io, T, uplo, .n, n, z.ptr, @intCast(t.lda), null);
        try tri.trtri(al, io, T, uplo, .n, n, t.a.ptr, @intCast(t.lda), null);
        try expectAgree(T, "trtri", n, z, t.a);
        t.store();
        @memcpy(z, t.a);
        try root.lauum(al, io, T, uplo, n, z.ptr, @intCast(t.lda));
        try tri.lauum(al, io, T, uplo, n, t.a.ptr, @intCast(t.lda));
        try expectAgree(T, "lauum", n, z, t.a);

        var g = try lu.Problem(T).init(al, n, n, 79 + @as(u64, @intFromBool(upper)));
        defer g.deinit(al);
        const y = try al.dupe(T, g.a);
        defer al.free(y);
        const ip = try al.alloc(isize, n);
        defer al.free(ip);
        try root.getrf(al, io, T, n, n, y.ptr, @intCast(g.lda), ip.ptr, null);
        try lu.getrf(al, io, T, n, n, g.a.ptr, @intCast(g.lda), g.ipiv.ptr, null);
        try testing.expectEqualSlices(isize, g.ipiv, ip);
        try expectAgree(T, "getrf", n, y, g.a);
    };
}

/// `f(ctx, io, out)` on `Io.failing` (tasks inline) and on `Io.Threaded`,
/// bit for bit; the BLAS calls take their workers by size.
pub fn expectIoBits(comptime T: type, ctx: anytype, comptime f: fn (@TypeOf(ctx), Io, []T) anyerror!void, len: usize, what: []const u8) !void {
    if (comptime !threading.hosted) return;
    const al = testing.allocator;
    const want = try al.alloc(T, len);
    defer al.free(want);
    const got = try al.alloc(T, len);
    defer al.free(got);
    try f(ctx, Io.failing, want);
    var th: Io.Threaded = .init(al, .{});
    defer th.deinit();
    try f(ctx, th.io(), got);
    for (want, got, 0..) |x, y, i| if (!chol.sameBits(T, x, y)) {
        check.say("{s} {s}: threaded differs from inline at element {d}\n", .{ what, @typeName(T), i });
        return error.TestUnexpectedResult;
    };
}

/// A diagonally dominant symmetric `n x n`, positive definite.
pub fn dominant(al: Allocator, n: usize, seed: u64) ![]f64 {
    const a = try al.alloc(f64, n * n);
    check.fill(f64, a, seed);
    for (0..n) |j| for (0..n) |i| {
        if (i > j) a[i + j * n] = a[j + i * n];
        if (i == j) a[i + j * n] += @floatFromInt(2 * n);
    };
    return a;
}

const big_n = 300;

const Big = struct {
    a: []const f64,

    fn potrfRun(self: Big, io: Io, out: []f64) anyerror!void {
        @memcpy(out, self.a);
        try potrf(testing.allocator, io, f64, .l, big_n, out.ptr, big_n, null);
    }

    fn getrfRun(self: Big, io: Io, out: []f64) anyerror!void {
        @memcpy(out, self.a);
        var ipiv: [big_n]isize = undefined;
        try getrf(testing.allocator, io, f64, big_n, big_n, out.ptr, big_n, &ipiv, null);
    }
};

test "relapack potrf and getrf: the same bits inline and threaded" {
    const al = testing.allocator;
    const a = try dominant(al, big_n, 1);
    defer al.free(a);
    try expectIoBits(f64, Big{ .a = a }, Big.potrfRun, big_n * big_n, "potrf");
    try expectIoBits(f64, Big{ .a = a }, Big.getrfRun, big_n * big_n, "getrf");
}
