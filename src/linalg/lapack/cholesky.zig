//! Cholesky factorization, solve and inverse. `potrf` is OpenBLAS's
//! `interface/lapack/potrf.c` and `zpotrf.c` over `lapack/potrf/potrf_U_single.c`
//! and `potrf_L_single.c`, or past upstream's size threshold
//! `potrf_U_parallel.c` and `potrf_L_parallel.c`; the base case is
//! `lapack/potf2/potf2_U.c`, `potf2_L.c`, `zpotf2_U.c` and `zpotf2_L.c`
//! (OpenBLAS v0.3.34-317), public as `potf2` (`interface/lapack/potf2.c`,
//! `zpotf2.c`). The semantics are Reference-LAPACK's `xPOTRF2`
//! (`lapack-netlib/SRC/zpotrf2.f` in that checkout).
//!
//! Where they differ, the reference wins:
//! - a NaN pivot stops the factorization (`DISNAN` in `xPOTRF2`), where
//!   upstream's `ajj <= 0` lets it through and reports success;
//! - the base case divides by the pivot, as rstd's `trsm` does, where upstream
//!   multiplies by its reciprocal (infinite for a denormal pivot).
//!
//! The trailing updates are the level-3 drivers' `trsm` and `herk` (`syrk`
//! for real `T`), called directly as upstream calls its kernels, so
//! Accelerate never runs them. Upstream takes the parallel algorithm (other
//! blocks) when it has more than one thread; here the size alone picks it
//! (`threads`), and the drivers' workers keep the bits, so no result depends
//! on the worker count. On macOS `potrf` goes to Accelerate's `?potrf_` when
//! routing is on (`apple_accel.potrf`).
//!
//! `potrs` is Reference-LAPACK's `xPOTRS` (`lapack-netlib/SRC/zpotrs.f`),
//! which OpenBLAS uses too: two triangular solves. `potri` is
//! `interface/lapack/potri.c` and `zpotri.c`: `trtri`, then `lauum`
//! (`triangular.zig`), as `xPOTRI`.

const std = @import("std");
const sysinfo = @import("builtin");
const blas = @import("../blas/__root__.zig");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const portable = @import("../blas/portable.zig");
const threading = @import("../blas/threading.zig");
const driver = @import("../blas/blas3/driver.zig");
const accel = @import("../blas/apple_accel.zig");
const Vtable = @import("../blas/Vtable.zig");
const tri = @import("triangular.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Uplo = types.Uplo;
const Float = types.Float;
const isComplex = types.isComplex;

/// The reference's `xerbla` argument, by name.
pub const ArgError = error{
    InvalidUplo,
    InvalidN,
    InvalidNrhs,
    InvalidLda,
    InvalidLdb,
};

/// The arguments of a routine on one square matrix (`checkSquare`).
pub const SquareError = error{ InvalidUplo, InvalidN, InvalidLda };

/// What `potrs` returns.
pub const PotrsError = ArgError || Allocator.Error || Io.Cancelable;

/// What `potrf` and `potf2` return: `NotPositiveDefinite` when a leading
/// minor is not positive definite.
pub const PotrfError = ArgError || Allocator.Error || Io.Cancelable || error{NotPositiveDefinite};

/// What `potri` returns: `Singular` when the factor has a zero on its
/// diagonal.
pub const PotriError = ArgError || Allocator.Error || Io.Cancelable || error{Singular};

/// `dtb_entries` of a hand-built vtable that has none (`blas2/triangular.zig`).
pub const default_dtb = 64;

pub fn floatType(comptime T: type, comptime name: []const u8) void {
    if (@typeInfo(Float(T)) != .float) @compileError(name ++ ": " ++ @typeName(T) ++ " is not a float type or its Complex");
}

/// The vtable for `T`; a type without slots runs the portable kernels and
/// never reads it (as `blas3/interface.zig`'s `table`).
pub fn table(comptime T: type) *const Vtable {
    return if (comptime portable.letter(T) != null) blas.vtable() else &no_table;
}

const no_table: Vtable = .{};

/// `dtb_entries` of `vt`, or `default_dtb`.
pub fn dtbOf(vt: *const Vtable) usize {
    return if (vt.dtb_entries > 0) @intCast(vt.dtb_entries) else default_dtb;
}

// ---- threads -------------------------------------------------------------------

/// The routines whose upstream interface picks a thread count by size.
pub const Routine = enum { getrf, gesv, getrs, potrf, potri, lauum, trtri, trtrs };

/// Upstream's interface rule for `r` over `T` at `size` (`m * n` for
/// `getrf`, `n * n` for `gesv`, else `n`): null below its threshold, where
/// upstream runs one thread and its single algorithm, else the workers it
/// allows (`num_cpu_avail`, capped as it caps them; `threading.workers`). The
/// size alone picks the algorithm, so one worker past the threshold still
/// runs the parallel one, and the bits never depend on the worker count.
pub fn threads(comptime r: Routine, comptime T: type, size: u64) ?u32 {
    const double = Float(T) == f64;
    const complex = comptime isComplex(T);
    const any = std.math.maxInt(u32);
    const cap: u64 = switch (r) {
        // `getrf.c`: below `nmax` one thread, else at most `m * n / nmax`;
        // `zgetrf.c`: below 10000
        .getrf => if (complex) (if (size < 10000) return null else any) else blk: {
            const nmax: u64 = if (double) 10000 else 40000;
            if (size < nmax) return null;
            break :blk size / nmax;
        },
        // `gesv.c`: `m * m` below 40000 (`DOUBLE` 10000); MSVC arm64 builds
        // up to 500 (complex 300), then 4 up to 1000
        .gesv => if (sysinfo.os.tag == .windows and sysinfo.cpu.arch == .aarch64)
            (if (size <= @as(u64, if (complex) 300 else 500)) return null else if (size <= 1000) 4 else any)
        else if (size < @as(u64, if (double) 10000 else 40000)) return null else any,
        // `getrs.c`, `trtrs.c`, `lauum.c` and their `z` forms: `n < 64`
        .getrs, .trtrs, .lauum => if (size < 64) return null else any,
        // `potrf.c`: below `nmax` one thread, else at most `n / nmax`;
        // `zpotrf.c`: below 64
        .potrf => if (complex) (if (size < 64) return null else any) else blk: {
            const nmax: u64 = if (double) 64 else 128;
            if (size < nmax) return null;
            break :blk size / nmax;
        },
        // `potri.c`: below 180; `zpotri.c`: below 200 (`DOUBLE` 150)
        .potri => if (size < @as(u64, if (!complex) 180 else if (double) 150 else 200)) return null else any,
        // `trtri.c`: up to 150; `ztrtri.c`: below 64
        .trtri => if (if (complex) size < 64 else size <= 150) return null else any,
    };
    return workers(cap);
}

/// `threading.workers` for a size cap, at most `threading.chunks_max`.
pub fn workers(size_max: u64) u32 {
    // selects the core first, which fills `blas.config`
    _ = blas.coreName();
    const w = threading.workers(&blas.config, std.math.cast(u32, size_max) orelse std.math.maxInt(u32));
    return @min(w, threading.chunks_max);
}

// ---- potrf, potf2 ---------------------------------------------------------------

/// `?POTRF`: the Cholesky factor of the Hermitian (real: symmetric)
/// positive-definite `n x n` `A`, from and into its `uplo` triangle:
/// `A = Uᴴ U` for `u`, `A = L Lᴴ` for `l`. The other triangle and the
/// imaginary parts of the diagonal are never read; the factor's diagonal
/// comes out real and positive.
///
/// When the leading minor of order `k` is not positive definite (its pivot
/// is not positive, or is NaN), returns `error.NotPositiveDefinite` and sets
/// `minor.*` to `k` (the reference's `INFO`). `A` then holds the factor's
/// first `k - 1` columns (`u`: rows), the failed pivot at `A(k, k)`, and
/// partial updates elsewhere. `minor` is written only then.
///
/// Routed to Accelerate (macOS), a failure leaves its partial factor in `A`.
pub fn potrf(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, a: [*]T, lda: isize, minor: ?*usize) PotrfError!void {
    comptime floatType(T, "potrf");
    if (!try checkSquare(T, uplo, n, lda)) return;
    if (accel.potrf(T, uplo, n, a, lda)) |info| return failed(@intCast(info), minor);
    const vt = table(T);
    const blk = driver.Blocking.of(T, vt);
    const upper = uplo == .u;
    const un: usize = @intCast(n);
    const la: usize = @intCast(lda);
    const k = if (threads(.potrf, T, un)) |w|
        try factorParallel(allocator, io, T, vt, upper, blk, un, a, la, w)
    else
        try factor(allocator, io, T, vt, upper, dtbOf(vt), blk.q, un, a, la);
    return failed(k, minor);
}

/// `?POTF2`: `potrf` unblocked (`potf2_U.c` and the rest, as loops), with
/// its errors and its `minor`.
pub fn potf2(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, a: [*]T, lda: isize, minor: ?*usize) PotrfError!void {
    comptime floatType(T, "potf2");
    _ = allocator;
    _ = io;
    if (!try checkSquare(T, uplo, n, lda)) return;
    return failed(unblocked(T, uplo == .u, @intCast(n), a, @intCast(lda)), minor);
}

/// The reference's checks of `uplo`, `n` and `lda` in its order; false for
/// `n = 0`, where nothing is read.
pub fn checkSquare(comptime T: type, uplo: Uplo, n: isize, lda: isize) SquareError!bool {
    if (uplo.isRev()) return error.InvalidUplo;
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, n)) return error.InvalidLda;
    if (n == 0) return false;
    if (!util.fitsMatrix(n, n, lda, @sizeOf(T))) return error.InvalidLda;
    return true;
}

/// `error.NotPositiveDefinite` for a nonzero `k`, written to `minor`.
fn failed(k: usize, minor: ?*usize) error{NotPositiveDefinite}!void {
    if (k == 0) return;
    if (minor) |m| m.* = k;
    return error.NotPositiveDefinite;
}

/// `potrf_U_single.c` and `potrf_L_single.c`: the order of the first minor
/// that is not positive definite, or 0. Upstream fuses the panel's `trsm`
/// with the update's `herk` column block by column block; here each is one
/// driver call, the same arithmetic.
fn factor(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, upper: bool, dtb: usize, q: usize, n: usize, a: [*]T, lda: usize) driver.Error!usize {
    if (n <= dtb / 2) return unblocked(T, upper, n, a, lda);
    // `(n + 3) / 4` and `n / 4` as upstream has them
    const blocking = if (n > 4 * q) q else if (upper) (n + 3) / 4 else @max(1, n / 4);
    var j: usize = 0;
    while (j < n) : (j += blocking) {
        const bk = @min(n - j, blocking);
        const d = a + j * (lda + 1);
        const k = try factor(allocator, io, T, vt, upper, dtb, q, bk, d, lda);
        if (k != 0) return k + j;
        try trailing(allocator, io, T, vt, upper, n - j - bk, bk, d, lda, 1);
    }
    return 0;
}

/// `potrf_U_parallel.c` and `potrf_L_parallel.c`: blocks of half the order
/// rounded up to `unroll_n`, at most `Q`, each factored recursively, until
/// `4 * unroll_n` is left for `factor`; the drivers take up to `want`
/// workers. Upstream drops to `factor` for one thread too; here only the
/// size does (`threads`).
fn factorParallel(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, upper: bool, blk: driver.Blocking, n: usize, a: [*]T, lda: usize, want: u32) driver.Error!usize {
    if (n <= 4 * blk.un) return factor(allocator, io, T, vt, upper, dtbOf(vt), blk.q, n, a, lda);
    const blocking = @min((n / 2 + blk.un - 1) / blk.un * blk.un, blk.q);
    var i: usize = 0;
    while (i < n) : (i += blocking) {
        const bk = @min(n - i, blocking);
        const d = a + i * (lda + 1);
        const k = try factorParallel(allocator, io, T, vt, upper, blk, bk, d, lda, want);
        if (k != 0) return k + i;
        try trailing(allocator, io, T, vt, upper, n - i - bk, bk, d, lda, want);
    }
    return 0;
}

/// The `rest` columns (`u`: rows) past the `bk` just factored at `d`:
/// `A12 := U11^-H A12; A22 -= A12^H A12`, or `A21 := A21 L11^-H;
/// A22 -= A21 A21^H` (`TRSM_LCUN`, `HERK_UC`; `TRSM_RCLN`, `HERK_LN`).
fn trailing(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, upper: bool, rest: usize, bk: usize, d: [*]T, lda: usize, want: u32) driver.Error!void {
    if (rest == 0) return;
    const complex = comptime isComplex(T);
    const one = scalar(T, 1);
    const minus = scalar(T, -1);
    const p = if (upper) d + bk * lda else d + bk;
    const c = d + bk * (lda + 1);
    if (upper) {
        try driver.trsm(allocator, io, T, vt, true, true, true, complex, false, bk, rest, one, d, lda, p, lda, want);
    } else {
        try driver.trsm(allocator, io, T, vt, false, false, true, complex, false, rest, bk, one, d, lda, p, lda, want);
    }
    try driver.syrk(allocator, io, T, vt, complex, upper, upper, rest, bk, minus, p, lda, null, 0, one, c, lda, want);
}

/// `?POTRS`: solves `A X = B` for the `n x nrhs` `X`, in place of `B`,
/// given `potrf`'s factor of `A` in its `uplo` triangle: `Uᴴ U X = B` or
/// `L Lᴴ X = B`, as two `trsm`s. The other triangle is never read.
pub fn potrs(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, nrhs: isize, a: [*]const T, lda: isize, b: [*]T, ldb: isize) PotrsError!void {
    comptime floatType(T, "potrs");
    if (uplo.isRev()) return error.InvalidUplo;
    if (n < 0) return error.InvalidN;
    if (nrhs < 0) return error.InvalidNrhs;
    if (lda < @max(1, n)) return error.InvalidLda;
    if (ldb < @max(1, n)) return error.InvalidLdb;
    if (n == 0 or nrhs == 0) return;
    if (!util.fitsMatrix(n, n, lda, @sizeOf(T))) return error.InvalidLda;
    if (!util.fitsMatrix(n, nrhs, ldb, @sizeOf(T))) return error.InvalidLdb;
    const vt = table(T);
    const complex = comptime isComplex(T);
    const m: usize = @intCast(n);
    const r: usize = @intCast(nrhs);
    const la: usize = @intCast(lda);
    const lb: usize = @intCast(ldb);
    const one = scalar(T, 1);
    const upper = uplo == .u;
    // `u`: `Uᴴ Y = B`, then `U X = Y`; `l`: `L Y = B`, then `Lᴴ X = Y`
    try driver.trsm(allocator, io, T, vt, true, upper, upper, upper and complex, false, m, r, one, a, la, b, lb, 1);
    try driver.trsm(allocator, io, T, vt, true, upper, !upper, !upper and complex, false, m, r, one, a, la, b, lb, 1);
}

/// `potf2_U.c`, `potf2_L.c` and their `z` forms: the unblocked factor, with
/// upstream's `dot`, `gemv` and `scal` as loops. A pivot that is not
/// positive is stored, and its order returned.
fn unblocked(comptime T: type, upper: bool, n: usize, a: [*]T, lda: usize) usize {
    const F = Float(T);
    for (0..n) |j| {
        const djj = &a[j * (lda + 1)];
        // `upper`: column `j` above the diagonal; else row `j` left of it
        const s: usize = if (upper) 1 else lda;
        const v = if (upper) a + j * lda else a + j;
        var sum: F = 0;
        for (0..j) |i| sum += abs2(T, v[i * s]);
        const ajj = re(T, djj.*) - sum;
        // `!(ajj > 0)`: NaN fails too, as `xPOTRF2`'s `DISNAN`
        if (!(ajj > 0)) {
            djj.* = scalar(T, ajj);
            return j + 1;
        }
        const d = @sqrt(ajj);
        djj.* = scalar(T, d);
        if (upper) {
            // row `j` right of the diagonal: `(A(j, k) - U(:j, j)^H U(:j, k)) / d`
            for (j + 1..n) |k| {
                const c = a + k * lda;
                var t: T = scalar(T, 0);
                for (0..j) |i| t = add(T, t, mulConj(T, v[i], c[i]));
                c[j] = divReal(T, sub(T, c[j], t), d);
            }
        } else {
            // column `j` below the diagonal: `(A(k, j) - L(k, :j) L(j, :j)^H) / d`
            const c = a + j * lda;
            for (0..j) |i| {
                const l = a + i * lda;
                const x = v[i * lda];
                for (j + 1..n) |k| c[k] = sub(T, c[k], mulConj(T, x, l[k]));
            }
            for (j + 1..n) |k| c[k] = divReal(T, c[k], d);
        }
    }
    return 0;
}

// ---- potri ---------------------------------------------------------------------

/// `?POTRI`: `A^-1` of the Hermitian (real: symmetric) positive-definite
/// `A` from `potrf`'s factor in its `uplo` triangle, into that triangle:
/// `U^-1 U^-H` or `L^-H L^-1`, as `trtri` (non-unit) then `lauum`. The other
/// triangle is never read.
///
/// A zero on the factor's diagonal returns `error.Singular` before `A` is
/// touched, its index (1-based) in `pivot`, as the reference's `xTRTRI`
/// reports it; upstream's `potri.c` calls `trtri`'s algorithm without that
/// check and divides by the zero.
pub fn potri(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, a: [*]T, lda: isize, pivot: ?*usize) PotriError!void {
    comptime floatType(T, "potri");
    if (!try checkSquare(T, uplo, n, lda)) return;
    const un: usize = @intCast(n);
    const la: usize = @intCast(lda);
    if (tri.firstZero(T, un, a, la)) |k| {
        if (pivot) |p| p.* = k;
        return error.Singular;
    }
    // upstream runs both parallel algorithms, or both single ones
    const w = threads(.potri, T, un);
    try tri.inverse(allocator, io, T, uplo == .u, false, un, a, la, w);
    try tri.product(allocator, io, T, uplo == .u, un, a, la, w);
}

// ---- element arithmetic ---------------------------------------------------------

pub fn scalar(comptime T: type, v: Float(T)) T {
    return if (comptime isComplex(T)) .init(v, 0) else v;
}

pub fn re(comptime T: type, v: T) Float(T) {
    return if (comptime isComplex(T)) v.re else v;
}

pub fn abs2(comptime T: type, v: T) Float(T) {
    return if (comptime isComplex(T)) v.re * v.re + v.im * v.im else v * v;
}

pub fn add(comptime T: type, x: T, y: T) T {
    return if (comptime isComplex(T)) .init(x.re + y.re, x.im + y.im) else x + y;
}

pub fn sub(comptime T: type, x: T, y: T) T {
    return if (comptime isComplex(T)) .init(x.re - y.re, x.im - y.im) else x - y;
}

pub fn mul(comptime T: type, x: T, y: T) T {
    if (comptime !isComplex(T)) return x * y;
    return .init(x.re * y.re - x.im * y.im, x.re * y.im + x.im * y.re);
}

pub fn conj(comptime T: type, v: T) T {
    return if (comptime isComplex(T)) .init(v.re, -v.im) else v;
}

pub fn neg(comptime T: type, v: T) T {
    return if (comptime isComplex(T)) .init(-v.re, -v.im) else -v;
}

/// `conj(x) * y`.
pub fn mulConj(comptime T: type, x: T, y: T) T {
    if (comptime !isComplex(T)) return x * y;
    return .init(x.re * y.re + x.im * y.im, x.re * y.im - x.im * y.re);
}

/// `x / d`, part by part: a real `d`, so each part's quotient is exact.
pub fn divReal(comptime T: type, x: T, d: Float(T)) T {
    return if (comptime isComplex(T)) .init(x.re / d, x.im / d) else x / d;
}

/// `x * d`, part by part, `d` real (`?dscal`).
pub fn mulReal(comptime T: type, x: T, d: Float(T)) T {
    return if (comptime isComplex(T)) .init(x.re * d, x.im * d) else x * d;
}

pub fn isZero(comptime T: type, v: T) bool {
    return if (comptime isComplex(T)) v.re == 0 and v.im == 0 else v == 0;
}

// ---- tests ---------------------------------------------------------------------

const testing = std.testing;
const check = @import("../blas/_test/blas1/check.zig");
const guard = @import("../blas/_test/blas2/guard.zig");

const nan = std.math.nan(f64);

/// Every element type the tests take: the four precisions, `XDouble` where
/// it is its own, and `f80` where it has no slots.
pub const test_types = .{ f32, f64, types.Complex32, types.Complex64 } ++ check.x_real ++ check.x_complex ++ check.f80_unslotted;

pub fn W(comptime T: type) type {
    return check.Wide(T);
}

pub fn fromF64(comptime T: type, v: f64) T {
    return if (comptime isComplex(T)) .init(@floatCast(v), 0) else @floatCast(v);
}

pub fn toT(comptime T: type, v: W(T)) T {
    return if (comptime isComplex(T)) .init(@floatCast(v.re), @floatCast(v.im)) else @floatCast(v);
}

/// NaN in every part.
pub fn poisoned(comptime T: type) T {
    const x = std.math.nan(Float(T));
    return if (comptime isComplex(T)) .init(x, x) else x;
}

/// Equal bits, part by part (`f80`'s padding aside).
pub fn sameBits(comptime T: type, x: T, y: T) bool {
    const F = Float(T);
    const U = @Int(.unsigned, @bitSizeOf(F));
    if (comptime isComplex(T)) return sameBits(F, x.re, y.re) and sameBits(F, x.im, y.im);
    return @as(U, @bitCast(x)) == @as(U, @bitCast(y));
}

pub fn isNanT(comptime T: type, v: T) bool {
    return if (comptime isComplex(T)) std.math.isNan(v.re) and std.math.isNan(v.im) else std.math.isNan(v);
}

/// The worker counts the threaded tests take past one.
pub const test_workers = [_]u32{ 2, 3, 8 };

/// `run(ctx, io, want, out)` at one worker on `Io.failing` against
/// `test_workers`, on `Io.failing` (tasks inline) and on a threaded `Io`,
/// bit for bit over the `len` elements it writes.
pub fn expectSameBits(comptime T: type, ctx: anytype, comptime run: fn (@TypeOf(ctx), Io, u32, []T) anyerror!void, len: usize, what: []const u8) !void {
    const al = testing.allocator;
    const want = try al.alloc(T, len);
    defer al.free(want);
    const got = try al.alloc(T, len);
    defer al.free(got);
    try run(ctx, Io.failing, 1, want);
    var th: if (threading.hosted) Io.Threaded else void = if (comptime threading.hosted) .init(al, .{}) else {};
    defer if (comptime threading.hosted) th.deinit();
    for (test_workers) |w| for (0..if (threading.hosted) 2 else 1) |k| {
        const io = if (comptime threading.hosted) (if (k == 0) Io.failing else th.io()) else Io.failing;
        try run(ctx, io, w, got);
        for (want, got, 0..) |x, y, i| if (!sameBits(T, x, y)) {
            check.say("{s} {s}: {d} workers ({s}) differ from one at element {d}: {any} against {any}\n", .{ what, @typeName(T), w, if (k == 0) "inline" else "threaded", i, y, x });
            return error.TestUnexpectedResult;
        };
    };
}

/// A test matrix: `L`, lower with a diagonal in [2, 4) and the rest within
/// `2 / n` (so `L` and `A = L Lᴴ` are well conditioned), and `A` built from
/// it in `W(T)`, stored in `uplo`'s triangle of `n x n` with `lda = n + 2`.
/// The other triangle, the padding rows and a complex diagonal's imaginary
/// parts hold NaN.
fn Problem(comptime T: type) type {
    return struct {
        n: usize,
        lda: usize,
        upper: bool,
        l: []T,
        a: []T,
        aw: []W(T),

        const Self = @This();

        fn init(al: Allocator, n: usize, upper: bool, seed: u64) !Self {
            const lda = n + 2;
            const l = try al.alloc(T, n * n);
            errdefer al.free(l);
            const a = try al.alloc(T, @max(1, lda * n));
            errdefer al.free(a);
            const aw = try al.alloc(W(T), n * n);
            check.fill(T, l, seed);
            const scale: Float(T) = @floatCast(1 / @as(f64, @floatFromInt(@max(n, 1))));
            for (0..n) |j| for (0..n) |i| {
                const v = &l[i + j * n];
                if (i < j) {
                    v.* = scalar(T, 0);
                } else if (i == j) {
                    v.* = scalar(T, 3 + re(T, v.*) / 2);
                } else if (comptime isComplex(T)) {
                    v.* = .init(v.re * scale, v.im * scale);
                } else {
                    v.* = v.* * scale;
                }
            };
            for (0..n) |j| for (0..n) |i| {
                var s: W(T) = check.wide(T, scalar(T, 0));
                for (0..@min(i, j) + 1) |k| s = check.addW(s, check.mulW(check.wide(T, l[i + k * n]), check.conjW(check.wide(T, l[j + k * n]))));
                aw[i + j * n] = s;
            };
            var self: Self = .{ .n = n, .lda = lda, .upper = upper, .l = l, .a = a, .aw = aw };
            self.store();
            return self;
        }

        /// `A` into `a`, poisoned where `potrf` must not read.
        fn store(self: *Self) void {
            @memset(self.a, poisoned(T));
            for (0..self.n) |j| for (0..self.n) |i| {
                if (if (self.upper) i > j else i < j) continue;
                const v = toT(T, self.aw[i + j * self.n]);
                self.a[i + j * self.lda] = if (comptime isComplex(T)) (if (i == j) .init(v.re, std.math.nan(Float(T))) else v) else v;
            };
        }

        fn deinit(self: *Self, al: Allocator) void {
            al.free(self.l);
            al.free(self.a);
            al.free(self.aw);
        }

        /// The computed factor's `(i, j)` as lower, `i >= j`: `L(i, j)`, or
        /// `conj(U(j, i))`.
        fn got(self: *const Self, i: usize, j: usize) T {
            if (!self.upper) return self.a[i + j * self.lda];
            const v = self.a[j + i * self.lda];
            return if (comptime isComplex(T)) v.conjugate() else v;
        }

        fn factorIt(self: *Self, al: Allocator, minor: ?*usize) PotrfError!void {
            return potrf(al, check.io, T, if (self.upper) .u else .l, @intCast(self.n), self.a.ptr, @intCast(self.lda), minor);
        }

        /// The factor against `A` (backward, within `2 (n + 2)` ulps of
        /// `|L| |L|ᴴ`) and against the known `L` (forward, loosely), with
        /// the poisoned elements untouched.
        fn expectFactor(self: *const Self) !void {
            const n = self.n;
            const F = Float(T);
            const eps: f128 = std.math.floatEps(F);
            const ulps: f128 = @floatFromInt(2 * (n + 2));
            for (0..n) |j| for (j..n) |i| {
                var s: W(T) = check.wide(T, scalar(T, 0));
                var mag: f128 = 0;
                for (0..j + 1) |k| {
                    const x = check.wide(T, self.got(i, k));
                    const y = check.wide(T, self.got(j, k));
                    s = check.addW(s, check.mulW(x, check.conjW(y)));
                    mag += check.abs1(x) * check.abs1(y);
                }
                const d = check.abs1(check.subW(s, self.aw[i + j * n]));
                if (!(d <= ulps * eps * mag)) {
                    check.say("potrf {s} n={d}: (L Lᴴ)({d},{d}) off by {e}, allowed {e}\n", .{ if (self.upper) "u" else "l", n, i, j, @as(f64, @floatCast(d)), @as(f64, @floatCast(ulps * eps * mag)) });
                    return error.TestUnexpectedResult;
                }
                const f = check.abs1(check.subW(check.wide(T, self.got(i, j)), check.wide(T, self.l[i + j * n])));
                if (!(f <= 64 * @as(f128, @floatFromInt(n)) * eps * 4)) {
                    check.say("potrf {s} n={d}: L({d},{d}) off by {e}\n", .{ if (self.upper) "u" else "l", n, i, j, @as(f64, @floatCast(f)) });
                    return error.TestUnexpectedResult;
                }
                if (comptime isComplex(T)) if (i == j and self.got(i, i).im != 0) {
                    check.say("potrf n={d}: L({d},{d}) has imaginary part {e}\n", .{ n, i, i, self.got(i, i).im });
                    return error.TestUnexpectedResult;
                };
            };
            for (0..n) |j| for (0..self.lda) |i| {
                const stored = i < n and (if (self.upper) i <= j else i >= j);
                if (!stored and !isNanT(T, self.a[i + j * self.lda])) {
                    check.say("potrf {s} n={d}: wrote ({d},{d}), outside its triangle\n", .{ if (self.upper) "u" else "l", n, i, j });
                    return error.TestUnexpectedResult;
                }
            };
        }
    };
}

/// `vt` with every precision's `q` and `dtb_entries` shrunk, so `potrf`'s
/// base case ends at 2 columns and its blocks are 3 wide past 12; static, as
/// a vtable is large for a rig's stack.
var shrunk: Vtable = undefined;

pub fn shrink(vt: *const Vtable) *const Vtable {
    shrunk = vt.*;
    shrunk.dtb_entries = 4;
    inline for (.{ "s", "d", "q", "c", "z", "x" }) |p| @field(shrunk, p ++ "gemm_q") = 3;
    return &shrunk;
}

/// Sizes around every block edge: on the shrunk blocking, and on a core's
/// own (`dtb_entries / 2` is 16 to 128 by core, so 130 is blocked on every
/// core, twice over where it is 32).
const shrunk_sizes = [_]usize{ 1, 2, 3, 4, 5, 7, 12, 13, 14, 17, 26 };
const own_sizes = [_]usize{ 1, 31, 33, 70, 130 };

fn factorCases(vt: *const Vtable) anyerror!void {
    const al = testing.allocator;
    inline for (test_types) |T| for ([_]bool{ true, false }) |upper| {
        for ([_]*const Vtable{ shrink(vt), vt }, 0..) |v, vi| {
            const saved = blas.testUse(v);
            defer blas.testRestore(saved);
            const sizes: []const usize = if (vi == 0) &shrunk_sizes else &own_sizes;
            for (sizes) |n| {
                if (check.rig() and n > 70) continue;
                var p = try Problem(T).init(al, n, upper, n * 31 + @as(u64, @intFromBool(upper)));
                defer p.deinit(al);
                var poison: guard.Poison = .{ .child = al };
                p.factorIt(poison.allocator(), null) catch |e| {
                    check.say("potrf {s} {s} n={d}\n", .{ @typeName(T), if (upper) "u" else "l", n });
                    return e;
                };
                p.expectFactor() catch |e| {
                    check.say("{s}, blocking {s}\n", .{ @typeName(T), if (vi == 0) "shrunk" else "own" });
                    return e;
                };
            }
        }
    };
}

test "potrf: A = L Lᴴ from a known L, both triangles, on every linked core" {
    try check.eachCore(factorCases);
}

/// `potrf` on `p` must fail at minor `want`.
fn expectMinor(comptime T: type, p: *Problem(T), want: usize, what: []const u8) !void {
    var minor: usize = 0;
    var poison: guard.Poison = .{ .child = testing.allocator };
    if (p.factorIt(poison.allocator(), &minor)) |_| {
        check.say("{s} {s} {s} n={d}: no error, want minor {d}\n", .{ what, @typeName(T), if (p.upper) "u" else "l", p.n, want });
        return error.TestExpectedError;
    } else |e| if (e != error.NotPositiveDefinite or minor != want) {
        check.say("{s} {s} {s} n={d}: {s} at minor {d}, want {d}\n", .{ what, @typeName(T), if (p.upper) "u" else "l", p.n, @errorName(e), minor, want });
        return error.TestUnexpectedResult;
    }
}

/// Positions of a bad element: the first, last and either side of block
/// edges on both blockings.
const bad_at = [_]usize{ 0, 1, 2, 3, 6, 12, 13, 25, 32, 33, 69 };

fn failCases(vt: *const Vtable) anyerror!void {
    const al = testing.allocator;
    const inf = std.math.inf(f64);
    inline for (test_types) |T| for ([_]bool{ true, false }) |upper| {
        for ([_]*const Vtable{ shrink(vt), vt }, 0..) |v, vi| {
            const saved = blas.testUse(v);
            defer blas.testRestore(saved);
            const n: usize = if (vi == 0) 26 else 70;
            var p = try Problem(T).init(al, n, upper, 7 + @as(u64, @intFromBool(upper)));
            defer p.deinit(al);
            for (bad_at) |k| {
                if (k >= n) continue;
                // a pivot made negative, NaN, or a NaN or infinity left of it
                for ([_]f64{ -100, nan }) |x| {
                    p.store();
                    p.a[k * (p.lda + 1)] = fromF64(T, x);
                    try expectMinor(T, &p, k + 1, "diagonal");
                }
                if (k == 0) continue;
                for ([_]f64{ nan, inf, -inf }) |x| for ([_]usize{ 0, k / 2, k - 1 }) |c| {
                    p.store();
                    // `(k, c)` in the lower triangle, `(c, k)` in the upper
                    p.a[if (upper) c + k * p.lda else k + c * p.lda] = fromF64(T, x);
                    try expectMinor(T, &p, k + 1, "off the diagonal");
                };
            }
        }
    };
}

test "potrf: the first minor that is not positive definite, by a pivot, NaN or infinity, on every linked core" {
    try check.eachCore(failCases);
}

test "potrf: a failed pivot is stored, real; minor is written only on failure" {
    // rstd's own (Accelerate's keeps a complex pivot's imaginary part)
    const was = accel.routing();
    _ = accel.route(false);
    defer _ = accel.route(was);
    inline for (.{ f64, types.Complex64 }) |T| for ([_]Uplo{ .u, .l }) |uplo| {
        // the identity with -3 at (2, 2): exact pivots
        var a: [16]T = @splat(scalar(T, 0));
        for (0..4) |i| a[i * 5] = scalar(T, 1);
        a[10] = if (comptime isComplex(T)) .init(-3, 5) else -3;
        var minor: usize = 99;
        try testing.expectError(error.NotPositiveDefinite, potrf(testing.allocator, check.io, T, uplo, 4, &a, 4, &minor));
        try testing.expectEqual(@as(usize, 3), minor);
        try testing.expectEqual(scalar(T, -3), a[10]);
        try testing.expectEqual(scalar(T, 1), a[5]);
        // a zero pivot fails too
        a[10] = scalar(T, 0);
        try testing.expectError(error.NotPositiveDefinite, potrf(testing.allocator, check.io, T, uplo, 4, &a, 4, &minor));
        try testing.expectEqual(@as(usize, 3), minor);
        a[10] = scalar(T, 4);
        minor = 99;
        try potrf(testing.allocator, check.io, T, uplo, 4, &a, 4, &minor);
        try testing.expectEqual(@as(usize, 99), minor);
        try testing.expectEqual(scalar(T, 2), a[10]);
    };
}

test "potrf: reports the first bad argument and reads nothing for n = 0" {
    const fa = testing.failing_allocator;
    var a: [4]f64 = @splat(nan);
    try testing.expectError(error.InvalidUplo, potrf(fa, check.io, f64, .v, -1, &a, 0, null));
    try testing.expectError(error.InvalidUplo, potrf(fa, check.io, f64, .m, -1, &a, 0, null));
    try testing.expectError(error.InvalidN, potrf(fa, check.io, f64, .u, -1, &a, 0, null));
    try testing.expectError(error.InvalidLda, potrf(fa, check.io, f64, .l, 2, &a, 1, null));
    try testing.expectError(error.InvalidLda, potrf(fa, check.io, f64, .u, 0, &a, 0, null));
    try testing.expectError(error.InvalidLda, potrf(fa, check.io, f64, .u, 2, &a, std.math.maxInt(isize) / 2, null));
    try potrf(fa, check.io, f64, .u, 0, &a, 1, null);
    for (a) |v| try testing.expect(std.math.isNan(v));
}

/// `potrs` on `potrf`'s factor of `p.a`, `nrhs` right-hand sides with
/// `ldb = n + 1` (the padding row NaN), then the residual `A X - B` in `W(T)`
/// within `8 (n + 2)` ulps of `|L| |L|ᴴ |X|`, and the padding untouched. With
/// `nan_at`, `B(i, j)` is NaN: column `j` must come out all NaN, the others
/// as without it, bit for bit.
fn solveCase(comptime T: type, al: Allocator, p: *Problem(T), nrhs: usize, seed: u64, nan_at: ?[2]usize) !void {
    const n = p.n;
    const ldb = n + 1;
    const b = try al.alloc(T, ldb * nrhs);
    defer al.free(b);
    const b0 = try al.alloc(T, n * nrhs);
    defer al.free(b0);
    check.fill(T, b0, seed);
    for (0..nrhs) |j| {
        @memcpy(b[j * ldb ..][0..n], b0[j * n ..][0..n]);
        b[j * ldb + n] = poisoned(T);
    }
    var poison: guard.Poison = .{ .child = al };
    const uplo: Uplo = if (p.upper) .u else .l;
    try potrs(poison.allocator(), check.io, T, uplo, @intCast(n), @intCast(nrhs), p.a.ptr, @intCast(p.lda), b.ptr, @intCast(ldb));
    const F = Float(T);
    const eps: f128 = std.math.floatEps(F);
    const ulps: f128 = @floatFromInt(8 * (n + 2));
    const t = try al.alloc(f128, n);
    defer al.free(t);
    for (0..nrhs) |j| {
        const x = b[j * ldb ..][0..n];
        try testing.expect(isNanT(T, b[j * ldb + n]));
        // `|L|ᴴ |x|`, then `|L|` times it
        for (0..n) |k| {
            t[k] = 0;
            for (k..n) |i| t[k] += check.abs1(check.wide(T, p.got(i, k))) * check.abs1(check.wide(T, x[i]));
        }
        for (0..n) |i| {
            var s: W(T) = check.wide(T, b0[i + j * n]);
            var mag: f128 = 0;
            for (0..n) |k| s = check.subW(s, check.mulW(p.aw[i + k * n], check.wide(T, x[k])));
            for (0..i + 1) |k| mag += check.abs1(check.wide(T, p.got(i, k))) * t[k];
            const d = check.abs1(s);
            if (!(d <= ulps * eps * mag)) {
                check.say("potrs {s} {s} n={d} nrhs={d}: residual ({d},{d}) {e}, allowed {e}\n", .{ @typeName(T), @tagName(uplo), n, nrhs, i, j, @as(f64, @floatCast(d)), @as(f64, @floatCast(ulps * eps * mag)) });
                return error.TestUnexpectedResult;
            }
        }
    }
    const at = nan_at orelse return;
    const bn = try al.alloc(T, ldb * nrhs);
    defer al.free(bn);
    for (0..nrhs) |j| @memcpy(bn[j * ldb ..][0..n], b0[j * n ..][0..n]);
    bn[at[0] + at[1] * ldb] = poisoned(T);
    try potrs(poison.allocator(), check.io, T, uplo, @intCast(n), @intCast(nrhs), p.a.ptr, @intCast(p.lda), bn.ptr, @intCast(ldb));
    for (0..nrhs) |j| for (0..n) |i| {
        const got = bn[i + j * ldb];
        const ok = if (j == at[1]) isNanT(T, got) else sameBits(T, got, b[i + j * ldb]);
        if (!ok) {
            check.say("potrs {s} {s} n={d}: NaN at B({d},{d}), X({d},{d}) is {any}\n", .{ @typeName(T), @tagName(uplo), n, at[0], at[1], i, j, got });
            return error.TestUnexpectedResult;
        }
    };
}

fn solveCases(vt: *const Vtable) anyerror!void {
    const al = testing.allocator;
    inline for (test_types) |T| for ([_]bool{ true, false }) |upper| {
        for ([_]*const Vtable{ shrink(vt), vt }, 0..) |v, vi| {
            const saved = blas.testUse(v);
            defer blas.testRestore(saved);
            const sizes: []const usize = if (vi == 0) &shrunk_sizes else &own_sizes;
            for (sizes) |n| {
                if (check.rig() and n > 70) continue;
                var p = try Problem(T).init(al, n, upper, n * 17 + @as(u64, @intFromBool(upper)));
                defer p.deinit(al);
                try p.factorIt(al, null);
                for ([_]usize{ 1, 3 }) |nrhs| {
                    const nan_at: [2]usize = .{ (n * 7 + nrhs) % n, nrhs / 2 };
                    solveCase(T, al, &p, nrhs, n + nrhs, nan_at) catch |e| {
                        check.say("{s}, blocking {s}\n", .{ @typeName(T), if (vi == 0) "shrunk" else "own" });
                        return e;
                    };
                }
            }
        }
    };
}

test "potrs: A X = B on potrf's factor, both triangles, NaN in one column of B, on every linked core" {
    try check.eachCore(solveCases);
}

test "potrs: reports the first bad argument and reads nothing when empty" {
    const fa = testing.failing_allocator;
    var a: [4]f64 = @splat(nan);
    var b: [4]f64 = @splat(nan);
    try testing.expectError(error.InvalidUplo, potrs(fa, check.io, f64, .v, -1, -1, &a, 0, &b, 0));
    try testing.expectError(error.InvalidN, potrs(fa, check.io, f64, .u, -1, -1, &a, 0, &b, 0));
    try testing.expectError(error.InvalidNrhs, potrs(fa, check.io, f64, .l, 0, -1, &a, 0, &b, 0));
    try testing.expectError(error.InvalidLda, potrs(fa, check.io, f64, .l, 2, 1, &a, 1, &b, 0));
    try testing.expectError(error.InvalidLdb, potrs(fa, check.io, f64, .l, 2, 1, &a, 2, &b, 1));
    try testing.expectError(error.InvalidLdb, potrs(fa, check.io, f64, .u, 0, 1, &a, 1, &b, 0));
    try testing.expectError(error.InvalidLdb, potrs(fa, check.io, f64, .u, 2, 2, &a, 2, &b, std.math.maxInt(isize) / 2));
    try potrs(fa, check.io, f64, .u, 0, 2, &a, 1, &b, 1);
    try potrs(fa, check.io, f64, .l, 2, 0, &a, 2, &b, 2);
    for (b) |v| try testing.expect(std.math.isNan(v));
}

fn parallelCases(vt: *const Vtable) anyerror!void {
    const al = testing.allocator;
    inline for (test_types) |T| for ([_]bool{ true, false }) |upper| {
        for ([_]*const Vtable{ shrink(vt), vt }, 0..) |v, vi| {
            const saved = blas.testUse(v);
            defer blas.testRestore(saved);
            const blk = driver.Blocking.of(T, v);
            const sizes: []const usize = if (vi == 0) &shrunk_sizes else &own_sizes;
            for (sizes) |n| {
                if (check.rig() and n > 70) continue;
                var p = try Problem(T).init(al, n, upper, n * 29 + @as(u64, @intFromBool(upper)));
                defer p.deinit(al);
                var poison: guard.Poison = .{ .child = al };
                const k = try factorParallel(poison.allocator(), check.io, T, table(T), upper, blk, n, p.a.ptr, p.lda, 1);
                try testing.expectEqual(@as(usize, 0), k);
                p.expectFactor() catch |e| {
                    check.say("{s}, blocking {s}, parallel algorithm\n", .{ @typeName(T), if (vi == 0) "shrunk" else "own" });
                    return e;
                };
                if (n < 13) continue;
                // a failure past the first block
                p.store();
                p.a[(n - 2) * (p.lda + 1)] = fromF64(T, -100);
                try testing.expectEqual(n - 1, try factorParallel(al, check.io, T, table(T), upper, blk, n, p.a.ptr, p.lda, 1));
            }
        }
    };
}

test "potrf's parallel algorithm: A = L Lᴴ, and the failing minor, on every linked core" {
    try check.eachCore(parallelCases);
}

/// Inputs of a threaded-bits case: `A`, copied into each run's `out`.
fn Threaded(comptime T: type) type {
    return struct {
        p: *const Problem(T),

        fn run(c: @This(), io: Io, want: u32, out: []T) anyerror!void {
            const vt = table(T);
            @memcpy(out, c.p.a);
            _ = try factorParallel(testing.allocator, io, T, vt, c.p.upper, driver.Blocking.of(T, vt), c.p.n, out.ptr, c.p.lda, want);
        }

        /// The public `potri`, threading by `blas.config`'s workers.
        fn inverse(c: @This(), io: Io, want: u32, out: []T) anyerror!void {
            const saved = blas.config.cpu_count;
            defer blas.config.cpu_count = saved;
            blas.config.cpu_count = want;
            @memcpy(out, c.p.a);
            try potri(testing.allocator, io, T, if (c.p.upper) .u else .l, @intCast(c.p.n), out.ptr, @intCast(c.p.lda), null);
        }
    };
}

test "potrf's parallel algorithm and potri give one worker's bits at any worker count" {
    _ = blas.coreName();
    const was = accel.routing();
    _ = accel.route(false);
    defer _ = accel.route(was);
    const al = testing.allocator;
    inline for (.{ f32, f64, types.Complex32, types.Complex64 }) |T| for ([_]bool{ true, false }) |upper| {
        // past every type's `potri` threshold
        const n: usize = if (check.rig()) 70 else 210;
        var p = try Problem(T).init(al, n, upper, 41);
        defer p.deinit(al);
        const c: Threaded(T) = .{ .p = &p };
        try expectSameBits(T, c, Threaded(T).run, p.a.len, "potrf");
        try p.factorIt(al, null);
        try expectSameBits(T, c, Threaded(T).inverse, p.a.len, "potri");
    };
}

test "potf2: the unblocked factor, its failing minor, its first bad argument" {
    const al = testing.allocator;
    inline for (test_types) |T| for ([_]bool{ true, false }) |upper| {
        const n = 26;
        var p = try Problem(T).init(al, n, upper, 5);
        defer p.deinit(al);
        try potf2(al, check.io, T, if (upper) .u else .l, n, p.a.ptr, @intCast(p.lda), null);
        try p.expectFactor();
        p.store();
        p.a[7 * (p.lda + 1)] = fromF64(T, nan);
        try expectMinorWith(T, &p, 8, potf2);
    };
    const fa = testing.failing_allocator;
    var a: [4]f64 = @splat(nan);
    try testing.expectError(error.InvalidUplo, potf2(fa, check.io, f64, .v, -1, &a, 0, null));
    try testing.expectError(error.InvalidN, potf2(fa, check.io, f64, .u, -1, &a, 0, null));
    try testing.expectError(error.InvalidLda, potf2(fa, check.io, f64, .l, 2, &a, 1, null));
    try potf2(fa, check.io, f64, .u, 0, &a, 1, null);
    for (a) |v| try testing.expect(std.math.isNan(v));
}

fn expectMinorWith(comptime T: type, p: *Problem(T), want: usize, comptime f: anytype) !void {
    var minor: usize = 0;
    const r = f(testing.allocator, check.io, T, if (p.upper) .u else .l, @intCast(p.n), p.a.ptr, @intCast(p.lda), &minor);
    try testing.expectError(error.NotPositiveDefinite, r);
    try testing.expectEqual(want, minor);
}

/// `potri` on `potrf`'s factor of `p`: `A X = I` within `32 (n + 2)` ulps
/// of `|A| |X|`, `X` Hermitian from its stored triangle, nothing outside
/// it written.
fn inverseCase(comptime T: type, al: Allocator, p: *Problem(T)) !void {
    const n = p.n;
    try p.factorIt(al, null);
    var poison: guard.Poison = .{ .child = al };
    try potri(poison.allocator(), check.io, T, if (p.upper) .u else .l, @intCast(n), p.a.ptr, @intCast(p.lda), null);
    const eps: f128 = std.math.floatEps(Float(T));
    const ulps: f128 = @floatFromInt(32 * (n + 2));
    for (0..n) |j| for (0..n) |i| {
        var s: W(T) = check.wide(T, scalar(T, if (i == j) -1 else 0));
        var mag: f128 = 0;
        for (0..n) |k| {
            const inside = if (p.upper) k <= j else k >= j;
            const x = check.wide(T, if (inside) p.a[k + j * p.lda] else conj(T, p.a[j + k * p.lda]));
            s = check.addW(s, check.mulW(p.aw[i + k * n], x));
            mag += check.abs1(p.aw[i + k * n]) * check.abs1(x);
        }
        const d = check.abs1(s);
        if (!(d <= ulps * eps * mag)) {
            check.say("potri {s} {s} n={d}: (A X - I)({d},{d}) = {e}, allowed {e}\n", .{ @typeName(T), if (p.upper) "u" else "l", n, i, j, @as(f64, @floatCast(d)), @as(f64, @floatCast(ulps * eps * mag)) });
            return error.TestUnexpectedResult;
        }
    };
    for (0..n) |j| for (0..p.lda) |i| {
        const stored = i < n and (if (p.upper) i <= j else i >= j);
        if (!stored and !isNanT(T, p.a[i + j * p.lda])) {
            check.say("potri {s} n={d}: wrote ({d},{d}), outside its triangle\n", .{ @typeName(T), n, i, j });
            return error.TestUnexpectedResult;
        }
    };
}

fn inverseCases(vt: *const Vtable) anyerror!void {
    const al = testing.allocator;
    inline for (test_types) |T| for ([_]bool{ true, false }) |upper| {
        for ([_]*const Vtable{ shrink(vt), vt }, 0..) |v, vi| {
            const saved = blas.testUse(v);
            defer blas.testRestore(saved);
            const sizes: []const usize = if (vi == 0) &shrunk_sizes else &own_sizes;
            for (sizes) |n| {
                if (check.rig() and n > 70) continue;
                var p = try Problem(T).init(al, n, upper, n * 11 + @as(u64, @intFromBool(upper)));
                defer p.deinit(al);
                inverseCase(T, al, &p) catch |e| {
                    check.say("{s}, blocking {s}\n", .{ @typeName(T), if (vi == 0) "shrunk" else "own" });
                    return e;
                };
            }
        }
    };
}

test "potri: A X = I on potrf's factor, both triangles, on every linked core" {
    try check.eachCore(inverseCases);
}

test "potri: a zero on the factor's diagonal before A is touched; the first bad argument" {
    const fa = testing.failing_allocator;
    inline for (.{ f64, types.Complex64 }) |T| for ([_]Uplo{ .u, .l }) |uplo| {
        var a: [9]T = @splat(scalar(T, 1));
        a[8] = scalar(T, 0);
        const before = a;
        var pivot: usize = 99;
        try testing.expectError(error.Singular, potri(fa, check.io, T, uplo, 3, &a, 3, &pivot));
        try testing.expectEqual(@as(usize, 3), pivot);
        try testing.expectEqualSlices(T, &before, &a);
    };
    var a: [4]f64 = @splat(nan);
    try testing.expectError(error.InvalidUplo, potri(fa, check.io, f64, .m, -1, &a, 0, null));
    try testing.expectError(error.InvalidN, potri(fa, check.io, f64, .u, -1, &a, 0, null));
    try testing.expectError(error.InvalidLda, potri(fa, check.io, f64, .l, 2, &a, 1, null));
    try potri(fa, check.io, f64, .u, 0, &a, 1, null);
    for (a) |v| try testing.expect(std.math.isNan(v));
}

test "potrf goes to Accelerate where it routes, and agrees with rstd's" {
    _ = blas.coreName();
    if (!accel.routing()) return error.SkipZigTest;
    const al = testing.allocator;
    inline for (.{ f32, f64, types.Complex32, types.Complex64 }) |T| for ([_]bool{ true, false }) |upper| {
        for ([_]usize{ 1, 7, 70, 130 }) |n| {
            var p = try Problem(T).init(al, n, upper, n * 3 + @as(u64, @intFromBool(upper)));
            defer p.deinit(al);
            const before = accel.calls.load(.monotonic);
            try p.factorIt(al, null);
            try testing.expect(accel.calls.load(.monotonic) > before);
            p.expectFactor() catch |e| {
                check.say("Accelerate's {s}potrf_\n", .{@typeName(T)});
                return e;
            };
            if (n < 7) continue;
            p.store();
            p.a[(n / 2) * (p.lda + 1)] = fromF64(T, nan);
            try expectMinor(T, &p, n / 2 + 1, "Accelerate's NaN pivot");
        }
    };
}
