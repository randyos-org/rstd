//! ReLAPACK's recursive algorithms (Elmar Peise, MIT), from OpenBLAS
//! v0.3.34-317's `relapack/src/`: `sygst`, `hegst` (`relapack_sygst.zig`)
//! and `pbtrf` (`relapack_pbtrf.zig`), and under
//! `-Drelapack` `potrf`, `getrf`, `trtri` and `lauum`, each split by `split`
//! (`?REC_SPLIT`) down to `crossover`, where the unblocked kernel runs.

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const threading = @import("../blas/threading.zig");
const driver = @import("../blas/blas3/driver.zig");
const threaded = @import("../blas/blas3/threaded.zig");
const laswp = @import("../blas/ext/laswp/kernel.zig");
const ops = @import("ops.zig");
const interface = @import("interface.zig");
const chol = @import("cholesky.zig");
const lu = @import("lu.zig");
const tri = @import("triangular.zig");
const relapack_sygst = @import("relapack_sygst.zig");
const relapack_pbtrf = @import("relapack_pbtrf.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Uplo = types.Uplo;
const Diag = types.Diag;
const Trans = types.Trans;
const Side = types.Side;
const isComplex = types.isComplex;
const scalar = ops.scalar;

const Error = driver.Error;

pub const Problem = relapack_sygst.Problem;
pub const SygstError = relapack_sygst.SygstError;
pub const sygst = relapack_sygst.sygst;
pub const hegst = relapack_sygst.hegst;
pub const PbtrfError = relapack_pbtrf.PbtrfError;
pub const pbtrf = relapack_pbtrf.pbtrf;

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
    try driver.gemm(al, io, T, interface.table(T), ta.transposed(), tb.transposed(), conj, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, w);
}

/// `?TRSM` (`solve`) or `?TRMM`, already checked.
pub fn trsm(al: Allocator, io: Io, comptime T: type, comptime solve: bool, side: Side, uplo: Uplo, ta: Trans, diag: Diag, m: usize, n: usize, alpha: T, a: [*]const T, lda: usize, b: [*]T, ldb: usize) Error!void {
    if (m == 0 or n == 0) return;
    const f = if (solve) driver.trsm else driver.trmm;
    const w = threaded.workers(.trsm, T, threading.product(&.{ m, n }));
    try f(al, io, T, interface.table(T), side == .l, uplo == .u, ta.transposed(), isComplex(T) and ta.conjA(), diag == .u, m, n, alpha, a, lda, b, ldb, w);
}

/// `?SYRK`, or with `herm` `?HERK` (`alpha`, `beta` real); `?SYR2K` and
/// `?HER2K` with `b`. `trans` is `t` or `c` for `op(A) = A'`.
pub fn syrk(al: Allocator, io: Io, comptime T: type, comptime herm: bool, uplo: Uplo, trans: Trans, n: usize, k: usize, alpha: T, a: [*]const T, lda: usize, b: ?[*]const T, ldb: usize, beta: T, c: [*]T, ldc: usize) Error!void {
    if (n == 0) return;
    const w = if (b == null)
        threaded.workers(.syrk, T, threading.product(&.{ n + 1, n, k }))
    else
        threaded.workers(.syr2k, T, threading.product(&.{ n, k }));
    try driver.syrk(al, io, T, interface.table(T), herm and isComplex(T), uplo == .u, trans != .n, n, k, alpha, a, lda, b, ldb, beta, c, ldc, w);
}

/// `?SYMM`, or with `herm` `?HEMM`.
pub fn symm(al: Allocator, io: Io, comptime T: type, comptime herm: bool, side: Side, uplo: Uplo, m: usize, n: usize, alpha: T, a: [*]const T, lda: usize, b: [*]const T, ldb: usize, beta: T, c: [*]T, ldc: usize) Error!void {
    if (m == 0 or n == 0) return;
    const w = threaded.workers(.symm, T, threading.product(&.{ 2, m, m, n }));
    try driver.symm(al, io, T, interface.table(T), herm and isComplex(T), side == .l, uplo == .u, m, n, alpha, a, lda, b, ldb, beta, c, ldc, w);
}

/// The Hermitian forms (`herk`, `C`) for complex `T`, the symmetric ones
/// (`syrk`, `T`) for real.
pub fn ctrans(comptime T: type) Trans {
    return if (comptime isComplex(T)) .c else .t;
}

// ---- potrf (`-Drelapack`) ------------------------------------------------------------

/// `?POTRF` by ReLAPACK's recursion (`relapack/src/zpotrf.c`, `dpotrf.c`
/// and the rest; OpenBLAS v0.3.34-317), as `cholesky.potrf` takes it, with
/// its errors. The base case is `potf2`.
pub fn potrf(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, a: []T, lda: usize) chol.PotrfError!void {
    comptime ops.requireFloat(T, "potrf");
    try interface.holdsTriangle(uplo, n, a.len, lda);
    if (n == 0) return;
    if (try potrfRec(allocator, io, T, uplo, n, a.ptr, lda) != null) return error.NotPositiveDefinite;
}

/// `RELAPACK_?potrf_rec`: the first pivot that is not positive (0-based),
/// or null.
pub fn potrfRec(al: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, a: [*]T, lda: usize) Error!?usize {
    if (n <= crossover) return chol.unblocked(T, uplo == .u, n, a, lda);
    const one = scalar(T, 1);
    const minus = scalar(T, -1);
    const n1 = split(T, n);
    const n2 = n - n1;
    const br = a + n1 * lda + n1;
    if (try potrfRec(al, io, T, uplo, n1, a, lda)) |k| return k;
    if (uplo == .l) {
        // A_BL = A_BL / A_TL'; A_BR -= A_BL A_BL'
        try trsm(al, io, T, true, .r, .l, ctrans(T), .n, n2, n1, one, a, lda, a + n1, lda);
        try syrk(al, io, T, true, .l, .n, n2, n1, minus, a + n1, lda, null, 0, one, br, lda);
    } else {
        // A_TR = A_TL' \ A_TR; A_BR -= A_TR' A_TR
        try trsm(al, io, T, true, .l, .u, ctrans(T), .n, n1, n2, one, a, lda, a + n1 * lda, lda);
        try syrk(al, io, T, true, .u, ctrans(T), n2, n1, minus, a + n1 * lda, lda, null, 0, one, br, lda);
    }
    return if (try potrfRec(al, io, T, uplo, n2, br, lda)) |k| k + n1 else null;
}

// ---- getrf (`-Drelapack`) -----------------------------------------------------------

/// `?GETRF` by ReLAPACK's recursion (`relapack/src/zgetrf.c`, `dgetrf.c`
/// and the rest; OpenBLAS v0.3.34-317), as `lu.getrf` takes it, with its
/// errors and result. The base case is `lu.getf2` (ReLAPACK calls
/// `xGETRF2`, whose pivots and semantics it shares).
///
/// The reference wins: a zero pivot in the left half does not stop the
/// recursion, which ReLAPACK's `if (*info) return` does, so the
/// factorization completes, as `xGETRF`'s.
pub fn getrf(allocator: Allocator, io: Io, comptime T: type, m: usize, n: usize, a: []T, lda: usize, ipiv: []usize) lu.GetrfError!?usize {
    comptime ops.requireFloat(T, "getrf");
    try util.holds(a.len, m, n, lda);
    if (ipiv.len < @min(m, n)) return error.SliceTooShort;
    if (m == 0 or n == 0) return null;
    const zero = try getrfRec(allocator, io, T, m, @min(m, n), a.ptr, lda, ipiv.ptr);
    // the right remainder: `A_R = A_L \ P A_R`
    if (m < n) {
        const ar = a.ptr + lda * m;
        laswp.swapRows(T, n - m, ar, lda, 0, m, ipiv.ptr, false);
        try trsm(allocator, io, T, true, .l, .l, .n, .u, m, n - m, scalar(T, 1), a.ptr, lda, ar, lda);
    }
    return zero;
}

/// `RELAPACK_?getrf_rec` on `m x n`, `m >= n`: the first zero pivot, or null.
fn getrfRec(al: Allocator, io: Io, comptime T: type, m: usize, n: usize, a: [*]T, lda: usize, ipiv: [*]usize) Error!?usize {
    if (m == 0 or n == 0) return null;
    if (n <= crossover) return lu.getf2(T, interface.table(T), m, n, a, lda, ipiv);
    const n1 = split(T, n);
    const n2 = n - n1;
    const m2 = m - n1;
    const ar = a + lda * n1;
    const bl = a + n1;
    const br = ar + n1;
    const zero = try getrfRec(al, io, T, m, n1, a, lda, ipiv);
    laswp.swapRows(T, n2, ar, lda, 0, n1, ipiv, false);
    // A_TR = A_TL \ A_TR; A_BR -= A_BL A_TR
    try trsm(al, io, T, true, .l, .l, .n, .u, n1, n2, scalar(T, 1), a, lda, ar, lda);
    try gemm(al, io, T, .n, .n, m2, n2, n1, scalar(T, -1), bl, lda, ar, lda, scalar(T, 1), br, lda);
    const zero2 = try getrfRec(al, io, T, m2, n2, br, lda, ipiv + n1);
    laswp.swapRows(T, n1, bl, lda, 0, n2, ipiv + n1, false);
    for (ipiv[n1..n]) |*p| p.* += n1;
    return interface.firstOf(zero, zero2, n1);
}

// ---- trtri, lauum (`-Drelapack`) -------------------------------------------------------

/// `?TRTRI` by ReLAPACK's recursion (`relapack/src/ztrtri.c`, `dtrtri.c`
/// and the rest; OpenBLAS v0.3.34-317), as `triangular.trtri` takes it,
/// with its errors. The base case is `trti2`.
///
/// The reference wins: any zero `A(i, i)` is `error.Singular`, where
/// ReLAPACK's singularity check sets `INFO = i - 1`, so a zero first
/// diagonal element reports success.
pub fn trtri(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, diag: Diag, n: usize, a: []T, lda: usize) tri.TrtriError!void {
    comptime ops.requireFloat(T, "trtri");
    try interface.holdsTriangle(uplo, n, a.len, lda);
    if (n == 0) return;
    if (diag == .n and tri.firstZero(T, n, a.ptr, lda) != null) return error.Singular;
    try trtriRec(allocator, io, T, uplo, diag, n, a.ptr, lda);
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
pub fn lauum(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, a: []T, lda: usize) tri.LauumError!void {
    comptime ops.requireFloat(T, "lauum");
    try interface.holdsTriangle(uplo, n, a.len, lda);
    if (n == 0) return;
    try lauumRec(allocator, io, T, uplo, n, a.ptr, lda);
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

/// A column-major matrix read as it is or transposed.
pub fn Tview(comptime T: type) type {
    return struct {
        p: [*]T,
        ld: usize,
        t: bool,

        pub inline fn at(v: @This(), i: usize, j: usize) *T {
            return &v.p[if (v.t) j + i * v.ld else i + j * v.ld];
        }
    };
}
