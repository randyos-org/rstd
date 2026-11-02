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
const indefinite = @import("indefinite.zig");
const sylvester = @import("sylvester.zig");

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

// ---- sytrf, hetrf and their rook forms -----------------------------------------------

/// What `sytrf`, `hetrf` and their `_rook` forms return: `Singular` when a
/// diagonal block of `D` is exactly zero.
pub const SytrfError = chol.SquareError || error{Singular} || Allocator.Error || Io.Cancelable;

/// `?SYTRF`: `A = U D Uᵀ` (`u`) or `L D Lᵀ` (`l`) of the symmetric
/// (complex: complex symmetric) `n x n` `A` in its `uplo` triangle, by Bunch
/// and Kaufman's diagonal pivoting; `D` and the multipliers in that triangle
/// and `ipiv` as LAPACK stores them: `ipiv[k]` positive for a 1 x 1 block
/// whose rows `k + 1` and `ipiv[k]` (1-based) were interchanged, equal and
/// negative on both rows of a 2 x 2 block.
///
/// When a block of `D` is exactly zero (or a NaN pivot is met in `tf2`), the
/// factorization still completes, then `error.Singular` sets `pivot.*` to
/// the reference's `INFO`, that block's first index (1-based), the first met:
/// the lowest for `l`, the highest for `u`.
///
/// ReLAPACK's recursion (`relapack/src/dsytrf.c` and the rest; OpenBLAS
/// v0.3.34-317) does not factor `A` past one level: it updates the rest
/// with panels already put back in LAPACK's form, defers updates that its
/// pivot search then reads stale, and its `_rook` forms call the non-rook
/// `xSYTF2` (`issues.md`). So this is Reference-LAPACK's blocked `xSYTRF`
/// (`lapack-netlib/SRC/dsytrf.f`, `zsytrf.f`, `zhetrf.f` and the `_rook`
/// forms in that checkout): panels of `sytrf_nb` columns by ReLAPACK's
/// `rec2` (`indefinite.panel`), each followed by `xLASYF`'s update of the
/// rest, here ReLAPACK's `gemmt`, then its interchanges undone; the last
/// columns by `indefinite.tf2`. Workspace (`n x sytrf_nb`) comes from
/// `allocator`.
pub fn sytrf(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, a: [*]T, lda: isize, ipiv: [*]isize, pivot: ?*usize) SytrfError!void {
    comptime chol.floatType(T, "sytrf");
    return indefiniteImpl(allocator, io, T, false, false, uplo, n, a, lda, ipiv, pivot);
}

/// `?HETRF` (`zhetrf.f`, `chetrf.f`; ReLAPACK's `zhetrf_rec2.c` panels):
/// `sytrf` for the Hermitian `A`, `A = U D Uᴴ` or `L D Lᴴ`; the diagonal's
/// imaginary parts are never read and `D`'s diagonal comes out real.
/// Complex only.
pub fn hetrf(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, a: [*]T, lda: isize, ipiv: [*]isize, pivot: ?*usize) SytrfError!void {
    comptime complexType(T, "hetrf");
    return indefiniteImpl(allocator, io, T, true, false, uplo, n, a, lda, ipiv, pivot);
}

/// `?SYTRF_ROOK` (`dsytrf_rook.f` and the rest; ReLAPACK's
/// `dsytrf_rook_rec2.c` panels): `sytrf` by the bounded Bunch-Kaufman
/// ("rook") pivoting, whose `ipiv` names both interchanges of a 2 x 2 block
/// (`-ipiv[k]` and `-ipiv[k + 1]` for `l`). No NaN test: only an exactly
/// zero block is singular, as the reference's.
pub fn sytrfRook(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, a: [*]T, lda: isize, ipiv: [*]isize, pivot: ?*usize) SytrfError!void {
    comptime chol.floatType(T, "sytrfRook");
    return indefiniteImpl(allocator, io, T, false, true, uplo, n, a, lda, ipiv, pivot);
}

/// `?HETRF_ROOK` (`zhetrf_rook.f`, `chetrf_rook.f`; ReLAPACK's
/// `zhetrf_rook_rec2.c` panels): `sytrfRook` for the Hermitian `A`.
/// Complex only.
pub fn hetrfRook(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, a: [*]T, lda: isize, ipiv: [*]isize, pivot: ?*usize) SytrfError!void {
    comptime complexType(T, "hetrfRook");
    return indefiniteImpl(allocator, io, T, true, true, uplo, n, a, lda, ipiv, pivot);
}

fn complexType(comptime T: type, comptime name: []const u8) void {
    if (!isComplex(T)) @compileError(name ++ ": " ++ @typeName(T) ++ " is not complex");
}

/// `ILAENV`'s block for `xSYTRF` and `xHETRF` and their `_ROOK` forms
/// (`lapack-netlib/SRC/ilaenv.f`); at or past `n`, `tf2` takes all of `A`.
pub const sytrf_nb = 64;

fn indefiniteImpl(al: Allocator, io: Io, comptime T: type, comptime herm: bool, comptime rook: bool, uplo: Uplo, n: isize, a: [*]T, lda: isize, ipiv: [*]isize, pivot: ?*usize) SytrfError!void {
    if (!try chol.checkSquare(T, uplo, n, lda)) return;
    const un: usize = @intCast(n);
    const blocked = sytrf_nb < un;
    const work = try al.alloc(T, if (blocked) un * sytrf_nb else 0);
    defer al.free(work);
    const info = try sytrfBlocked(al, io, T, herm, rook, uplo == .u, un, a, @intCast(lda), ipiv, work.ptr, un);
    return lu.singular(info, pivot);
}

/// `xSYTRF`'s loop: panels from the first column (upper: the last), each
/// `panel`, the update of the rest, `undoPanel`; the reference's `INFO`.
fn sytrfBlocked(al: Allocator, io: Io, comptime T: type, comptime herm: bool, comptime rook: bool, upper: bool, n: usize, a: [*]T, lda: usize, ipiv: [*]isize, w: [*]T, ldw: usize) Error!usize {
    const one = scalar(T, 1);
    const minus = scalar(T, -1);
    const nb = sytrf_nb;
    var info: usize = 0;
    if (!upper) {
        var k: usize = 0;
        while (k < n) {
            const m = n - k;
            const d = a + k * (lda + 1);
            const p = ipiv + k;
            var kb = m;
            var iinfo: usize = 0;
            if (k + nb < n) {
                const r = indefinite.panel(T, herm, rook, false, m, nb, d, lda, p, w, ldw);
                kb = r.kb;
                iinfo = r.info;
                // A22 -= L21 W21ᵀ (`xLASYF`), then L21 in LAPACK's form
                try gemmtRec(al, io, T, .l, m - kb, kb, minus, d + kb, lda, w + kb, ldw, one, d + kb * (lda + 1), lda);
                indefinite.undoPanel(T, rook, false, m, d, lda, p, kb);
            } else iinfo = indefinite.tf2(T, herm, rook, false, m, d, lda, p);
            if (info == 0 and iinfo != 0) info = iinfo + k;
            const off: isize = @intCast(k);
            for (p[0..kb]) |*v| v.* += if (v.* > 0) off else -off;
            k += kb;
        }
        return info;
    }
    var k = n;
    while (k > 0) {
        var kb = k;
        var iinfo: usize = 0;
        if (k > nb) {
            const r = indefinite.panel(T, herm, rook, true, k, nb, a, lda, ipiv, w, ldw);
            kb = r.kb;
            iinfo = r.info;
            // A11 -= U12 W12ᵀ (`xLASYF`), W's last `kb` columns
            try gemmtRec(al, io, T, .u, k - kb, kb, minus, a + (k - kb) * lda, lda, w + (nb - kb) * ldw, ldw, one, a, lda);
            indefinite.undoPanel(T, rook, true, k, a, lda, ipiv, kb);
        } else iinfo = indefinite.tf2(T, herm, rook, true, k, a, lda, ipiv);
        if (info == 0 and iinfo != 0) info = iinfo;
        k -= kb;
    }
    return info;
}

/// `RELAPACK_?gemmt` as its `sytrf` calls it (`relapack/src/dgemmt.c` and
/// the rest, `HAVE_XGEMMT` 0): `C := alpha A Bᵀ + beta C` on the `uplo`
/// triangle of `n x n` `C`, `A` and `B` `n x k`; halves by `gemm` down to the
/// crossover (`CROSSOVER_XGEMMT`, `sytrf`'s), then `gemv` per column.
fn gemmtRec(al: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, k: usize, alpha: T, a: [*]const T, lda: usize, b: [*]const T, ldb: usize, beta: T, c: [*]T, ldc: usize) Error!void {
    if (n <= crossover) {
        // `RELAPACK_?gemmt_rec2`: `xGEMV` `N` per column, `beta` first
        for (0..n) |i| {
            const lo = if (uplo == .l) i else 0;
            const hi = if (uplo == .l) n else i + 1;
            const y = c + i * ldc;
            if (!chol.isZero(T, chol.sub(T, beta, scalar(T, 1)))) for (lo..hi) |r| {
                y[r] = if (chol.isZero(T, beta)) scalar(T, 0) else chol.mul(T, beta, y[r]);
            };
            for (0..k) |j| {
                const t = chol.mul(T, alpha, b[i + j * ldb]);
                for (lo..hi) |r| y[r] = chol.add(T, y[r], chol.mul(T, t, a[r + j * lda]));
            }
        }
        return;
    }
    const n1 = split(T, n);
    const n2 = n - n1;
    try gemmtRec(al, io, T, uplo, n1, k, alpha, a, lda, b, ldb, beta, c, ldc);
    if (uplo == .l) {
        try gemm(al, io, T, .n, .t, n2, n1, k, alpha, a + n1, lda, b, ldb, beta, c + n1, ldc);
    } else {
        try gemm(al, io, T, .n, .t, n1, n2, k, alpha, a, lda, b + n1, ldb, beta, c + n1 * ldc, ldc);
    }
    try gemmtRec(al, io, T, uplo, n2, k, alpha, a + n1, lda, b + n1, ldb, beta, c + n1 * ldc + n1, ldc);
}

// ---- sygst, hegst ---------------------------------------------------------------------

/// What `sygst` and `hegst` return: bad arguments in the reference's order,
/// `error.OutOfMemory` or `error.Canceled`.
pub const SygstError = error{ InvalidItype, InvalidUplo, InvalidN, InvalidLda, InvalidLdb } || Allocator.Error || Io.Cancelable;

/// `?SYGST` by ReLAPACK's recursion (`relapack/src/dsygst.c`, `ssygst.c`;
/// OpenBLAS v0.3.34-317): the symmetric-definite `A x = λ B x` (`itype` 1),
/// `A B x = λ x` (2) or `B A x = λ x` (3) to standard form, `B` holding
/// `potrf`'s factor in the same `uplo` triangle: `A := inv(Uᵀ) A inv(U)` or
/// `inv(L) A inv(Lᵀ)` for 1, `U A Uᵀ` or `Lᵀ A L` for 2 and 3, into `A`'s
/// `uplo` triangle. Real `T`; the other triangles are never read.
///
/// The base case is Reference-LAPACK's `xSYGS2` (`gs2`). The workspace
/// ReLAPACK takes from `malloc` (`n1 (n - n1)`, `XSYGST_ALLOW_MALLOC`) comes
/// from `allocator`; it holds `-1/2 B21 A11` (or its like) once, where
/// otherwise `symm` runs twice. Every level uses it, as `zhegst.c`'s `>=`
/// has it; `dsygst.c`, `ssygst.c` and `chegst.c` test `>`, which leaves the
/// top level's exactly sized buffer unused.
pub fn sygst(allocator: Allocator, io: Io, comptime T: type, itype: isize, uplo: Uplo, n: isize, a: [*]T, lda: isize, b: [*]const T, ldb: isize) SygstError!void {
    comptime chol.floatType(T, "sygst");
    if (comptime isComplex(T)) @compileError("sygst: " ++ @typeName(T) ++ " is complex; hegst takes it");
    return gstImpl(allocator, io, T, itype, uplo, n, a, lda, b, ldb);
}

/// `?HEGST` (`relapack/src/zhegst.c`, `chegst.c`): `sygst` for the
/// Hermitian-definite problem, `ᴴ` for `ᵀ`; `A`'s diagonal comes out real.
/// Complex only.
pub fn hegst(allocator: Allocator, io: Io, comptime T: type, itype: isize, uplo: Uplo, n: isize, a: [*]T, lda: isize, b: [*]const T, ldb: isize) SygstError!void {
    comptime complexType(T, "hegst");
    return gstImpl(allocator, io, T, itype, uplo, n, a, lda, b, ldb);
}

fn gstImpl(al: Allocator, io: Io, comptime T: type, itype: isize, uplo: Uplo, n: isize, a: [*]T, lda: isize, b: [*]const T, ldb: isize) SygstError!void {
    if (itype < 1 or itype > 3) return error.InvalidItype;
    if (uplo.isRev()) return error.InvalidUplo;
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, n)) return error.InvalidLda;
    if (ldb < @max(1, n)) return error.InvalidLdb;
    if (n == 0) return;
    if (!util.fitsMatrix(n, n, lda, @sizeOf(T))) return error.InvalidLda;
    if (!util.fitsMatrix(n, n, ldb, @sizeOf(T))) return error.InvalidLdb;
    const un: usize = @intCast(n);
    const n1 = split(T, un);
    const work = try al.alloc(T, n1 * (un - n1));
    defer al.free(work);
    try gstRec(al, io, T, itype == 1, uplo, un, a, @intCast(lda), b, @intCast(ldb), work);
}

/// `RELAPACK_?sygst_rec`; `inv` for `itype` 1.
fn gstRec(al: Allocator, io: Io, comptime T: type, inv: bool, uplo: Uplo, n: usize, a: [*]T, lda: usize, b: [*]const T, ldb: usize, work: []T) Error!void {
    if (n <= crossover) return gs2(T, inv, uplo == .u, n, a, lda, b, ldb);
    const zero = scalar(T, 0);
    const one = scalar(T, 1);
    const minus = scalar(T, -1);
    const half = scalar(T, 0.5);
    const mhalf = scalar(T, -0.5);
    const ct = ctrans(T);
    const n1 = split(T, n);
    const n2 = n - n1;
    const a_tr = a + lda * n1;
    const a_bl = a + n1;
    const a_br = a + lda * n1 + n1;
    const b_tr = b + ldb * n1;
    const b_bl = b + n1;
    const b_br = b + ldb * n1 + n1;
    try gstRec(al, io, T, inv, uplo, n1, a, lda, b, ldb, work);
    // `T`, computed once into the workspace, or `symm` twice
    const use = work.len >= n2 * n1;
    const lower = uplo == .l;
    // the off-diagonal block, its rows and columns
    const x = if (lower) a_bl else a_tr;
    const xm = if (lower) n2 else n1;
    const xn = if (lower) n1 else n2;
    if (inv) {
        if (lower) {
            // A_BL = A_BL / B_TL'
            try trsm(al, io, T, true, .r, .l, ct, .n, n2, n1, one, b, ldb, a_bl, lda);
        } else {
            // A_TR = B_TL' \ A_TR
            try trsm(al, io, T, true, .l, .u, ct, .n, n1, n2, one, b, ldb, a_tr, lda);
        }
        // X -= 1/2 B_BL A_TL (`u`: A_TL B_TR), twice around the rank-2 update
        const side: Side = if (lower) .r else .l;
        const s = if (lower) b_bl else b_tr;
        if (use) {
            try symm(al, io, T, true, side, uplo, xm, xn, mhalf, a, lda, s, ldb, zero, work.ptr, xm);
            addBlock(T, xm, xn, work.ptr, xm, x, lda);
        } else try symm(al, io, T, true, side, uplo, xm, xn, mhalf, a, lda, s, ldb, one, x, lda);
        // A_BR -= A_BL B_BL' + B_BL A_BL' (`u`: A_TR' B_TR + B_TR' A_TR)
        try syrk(al, io, T, true, uplo, if (lower) .n else ct, n2, n1, minus, x, lda, s, ldb, one, a_br, lda);
        if (use) {
            addBlock(T, xm, xn, work.ptr, xm, x, lda);
        } else try symm(al, io, T, true, side, uplo, xm, xn, mhalf, a, lda, s, ldb, one, x, lda);
        if (lower) {
            // A_BL = B_BR \ A_BL
            try trsm(al, io, T, true, .l, .l, .n, .n, n2, n1, one, b_br, ldb, a_bl, lda);
        } else {
            // A_TR = A_TR / B_BR
            try trsm(al, io, T, true, .r, .u, .n, .n, n1, n2, one, b_br, ldb, a_tr, lda);
        }
    } else {
        if (lower) {
            // A_BL = A_BL B_TL
            try trsm(al, io, T, false, .r, .l, .n, .n, n2, n1, one, b, ldb, a_bl, lda);
        } else {
            // A_TR = B_TL A_TR
            try trsm(al, io, T, false, .l, .u, .n, .n, n1, n2, one, b, ldb, a_tr, lda);
        }
        // X += 1/2 A_BR B_BL (`u`: B_TR A_BR), twice around the rank-2 update
        const side: Side = if (lower) .l else .r;
        const s = if (lower) b_bl else b_tr;
        if (use) {
            try symm(al, io, T, true, side, uplo, xm, xn, half, a_br, lda, s, ldb, zero, work.ptr, xm);
            addBlock(T, xm, xn, work.ptr, xm, x, lda);
        } else try symm(al, io, T, true, side, uplo, xm, xn, half, a_br, lda, s, ldb, one, x, lda);
        // A_TL += A_BL' B_BL + B_BL' A_BL (`u`: A_TR B_TR' + B_TR A_TR')
        try syrk(al, io, T, true, uplo, if (lower) ct else .n, n1, n2, one, x, lda, s, ldb, one, a, lda);
        if (use) {
            addBlock(T, xm, xn, work.ptr, xm, x, lda);
        } else try symm(al, io, T, true, side, uplo, xm, xn, half, a_br, lda, s, ldb, one, x, lda);
        if (lower) {
            // A_BL = B_BR' A_BL
            try trsm(al, io, T, false, .l, .l, ct, .n, n2, n1, one, b_br, ldb, a_bl, lda);
        } else {
            // A_TR = A_TR B_BR'
            try trsm(al, io, T, false, .r, .u, ct, .n, n1, n2, one, b_br, ldb, a_tr, lda);
        }
    }
    try gstRec(al, io, T, inv, uplo, n2, a_br, lda, b_br, ldb, work);
}

/// `Y += X` on `m x n` (`?AXPY` with 1, a column at a time).
fn addBlock(comptime T: type, m: usize, n: usize, x: [*]const T, ldx: usize, y: [*]T, ldy: usize) void {
    for (0..n) |j| for (0..m) |i| {
        y[i + j * ldy] = chol.add(T, y[i + j * ldy], x[i + j * ldx]);
    };
}

/// A column-major matrix read as it is or transposed.
fn Tview(comptime T: type) type {
    return struct {
        p: [*]T,
        ld: usize,
        t: bool,

        inline fn at(v: @This(), i: usize, j: usize) *T {
            return &v.p[if (v.t) j + i * v.ld else i + j * v.ld];
        }
    };
}

/// Reference-LAPACK's `xSYGS2` and `xHEGS2` (`lapack-netlib/SRC/dsygs2.f`,
/// `zhegs2.f`; OpenBLAS v0.3.34-317), the unblocked `sygst`, with its
/// `syr2`/`her2`, `trsv`, `trmv`, `axpy` and `scal` as loops in the
/// reference BLAS's order. `upper` runs the lower algorithm on the
/// transposes, the conjugates of the reference's upper vectors (it conjugates
/// a row to get a column), so the same arithmetic up to their loop orders.
fn gs2(comptime T: type, inv: bool, upper: bool, n: usize, a_: [*]T, lda: usize, b_: [*]const T, ldb: usize) void {
    const herm = comptime isComplex(T);
    const a: Tview(T) = .{ .p = a_, .ld = lda, .t = upper };
    const b: Tview(T) = .{ .p = @constCast(b_), .ld = ldb, .t = upper };
    const cj = struct {
        fn f(v: T) T {
            return if (comptime herm) v.conjugate() else v;
        }
    }.f;
    const mul = chol.mul;
    const add = chol.add;
    for (0..n) |k| {
        const akk0 = chol.re(T, a.at(k, k).*);
        const bkk = chol.re(T, b.at(k, k).*);
        if (inv) {
            // inv(L) A inv(Lᴴ): the column below `k`, then the rest
            const akk = akk0 / (bkk * bkk);
            a.at(k, k).* = scalar(T, akk);
            if (k + 1 == n) continue;
            const r = 1 / bkk;
            for (k + 1..n) |i| a.at(i, k).* = chol.mulReal(T, a.at(i, k).*, r);
            const c = -0.5 * akk;
            for (k + 1..n) |i| a.at(i, k).* = add(T, a.at(i, k).*, chol.mulReal(T, b.at(i, k).*, c));
            rank2(T, a, k + 1, n, k, b, k, scalar(T, -1));
            for (k + 1..n) |i| a.at(i, k).* = add(T, a.at(i, k).*, chol.mulReal(T, b.at(i, k).*, c));
            // x := inv(B22) x, `xTRSV` `N`
            for (k + 1..n) |j| {
                const xj = a.at(j, k);
                if (chol.isZero(T, xj.*)) continue;
                xj.* = divide(T, xj.*, b.at(j, j).*);
                const t = xj.*;
                for (j + 1..n) |i| a.at(i, k).* = chol.sub(T, a.at(i, k).*, mul(T, t, b.at(i, j).*));
            }
        } else {
            // Lᴴ A L: row `k` left of the diagonal, as a conjugated column
            // x := B11ᴴ x, `xTRMV` `C`
            for (0..k) |j| {
                var t = cj(a.at(k, j).*);
                t = mul(T, t, cj(b.at(j, j).*));
                for (j + 1..k) |i| t = add(T, t, mul(T, cj(b.at(i, j).*), cj(a.at(k, i).*)));
                a.at(k, j).* = cj(t);
            }
            const c = 0.5 * akk0;
            for (0..k) |j| a.at(k, j).* = cj(add(T, cj(a.at(k, j).*), chol.mulReal(T, cj(b.at(k, j).*), c)));
            rank2Row(T, a, k, b, scalar(T, 1));
            for (0..k) |j| a.at(k, j).* = cj(chol.mulReal(T, add(T, cj(a.at(k, j).*), chol.mulReal(T, cj(b.at(k, j).*), c)), bkk));
            a.at(k, k).* = scalar(T, akk0 * (bkk * bkk));
        }
    }
}

/// `xTRSV`'s division by the diagonal.
fn divide(comptime T: type, x: T, d: T) T {
    if (comptime !isComplex(T)) return x / d;
    const z = util.cdiv(Float(T), x.re, x.im, d.re, d.im);
    return .init(z[0], z[1]);
}

/// `xSYR2` or `xHER2` (lower) with `alpha` on `A(lo:hi, lo:hi)`, `x` column
/// `xc` of `A` and `y` column `yc` of `B`, both from row `lo`.
fn rank2(comptime T: type, a: Tview(T), lo: usize, hi: usize, xc: usize, b: Tview(T), yc: usize, alpha: T) void {
    const herm = comptime isComplex(T);
    for (lo..hi) |j| {
        const xj = a.at(j, xc).*;
        const yj = b.at(j, yc).*;
        if (chol.isZero(T, xj) and chol.isZero(T, yj)) {
            if (comptime herm) a.at(j, j).* = scalar(T, a.at(j, j).re);
            continue;
        }
        const t1 = chol.mul(T, alpha, chol.conj(T, yj));
        const t2 = chol.conj(T, chol.mul(T, alpha, xj));
        if (comptime herm) {
            const d = chol.add(T, chol.mul(T, xj, t1), chol.mul(T, yj, t2));
            a.at(j, j).* = scalar(T, a.at(j, j).re + d.re);
        }
        for (if (herm) j + 1 else j..hi) |i| {
            a.at(i, j).* = chol.add(T, a.at(i, j).*, chol.add(T, chol.mul(T, a.at(i, xc).*, t1), chol.mul(T, b.at(i, yc).*, t2)));
        }
    }
}

/// `rank2` with `x` and `y` the conjugated rows `k` of `A` and `B` left of
/// the diagonal, on `A(:k, :k)`.
fn rank2Row(comptime T: type, a: Tview(T), k: usize, b: Tview(T), alpha: T) void {
    const herm = comptime isComplex(T);
    const cj = chol.conj;
    for (0..k) |j| {
        const xj = cj(T, a.at(k, j).*);
        const yj = cj(T, b.at(k, j).*);
        if (chol.isZero(T, xj) and chol.isZero(T, yj)) {
            if (comptime herm) a.at(j, j).* = scalar(T, a.at(j, j).re);
            continue;
        }
        const t1 = chol.mul(T, alpha, cj(T, yj));
        const t2 = cj(T, chol.mul(T, alpha, xj));
        if (comptime herm) {
            const d = chol.add(T, chol.mul(T, xj, t1), chol.mul(T, yj, t2));
            a.at(j, j).* = scalar(T, a.at(j, j).re + d.re);
        }
        for (if (herm) j + 1 else j..k) |i| {
            a.at(i, j).* = chol.add(T, a.at(i, j).*, chol.add(T, chol.mul(T, cj(T, a.at(k, i).*), t1), chol.mul(T, cj(T, b.at(k, i).*), t2)));
        }
    }
}

// ---- pbtrf -----------------------------------------------------------------------------

/// What `pbtrf` returns: bad arguments in the reference's order,
/// `NotPositiveDefinite` (with the minor, as `potrf`'s), `error.OutOfMemory`
/// or `error.Canceled`.
pub const PbtrfError = error{ InvalidUplo, InvalidN, InvalidKd, InvalidLdab, NotPositiveDefinite } || Allocator.Error || Io.Cancelable;

/// `?PBTRF` by ReLAPACK's recursion (`relapack/src/dpbtrf.c`, `zpbtrf.c`
/// and the rest; OpenBLAS v0.3.34-317): the Cholesky factor of the Hermitian
/// (real: symmetric) positive-definite band matrix with `kd` off-diagonals
/// in LAPACK's band storage `ab` (`ldab >= kd + 1`; `u`: `A(i, j)` at row
/// `kd + i - j` of column `j`, `l`: at row `i - j`), into the same storage.
/// Blocks of up to `kd` columns are `potrf`'s (this file's), the part below
/// them a `trsm` and `herk`, the triangle past the band's edge through a
/// workspace, as ReLAPACK; the base case is Reference-LAPACK's `xPBTF2`
/// (`pbtf2`). A failure is `potrf`'s: `minor.*` gets the order of the first
/// leading minor that is not positive definite (NaN included, as `potrf`'s).
///
/// The reference wins where ReLAPACK factors the rest as a dense matrix
/// whenever `kd` exceeds the block (`if (*kd > n1)`), which reads and
/// writes outside the band once the rest is wider than `kd + 1`; there it
/// recurses on the band. The workspace is `n1 x n1` for the top level's
/// block `n1`, which every level's `n22 x n1` fits (`dpbtrf.c`'s `n22`;
/// `zpbtrf.c` bounds it by `kd`, past the block).
pub fn pbtrf(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, kd: isize, ab: [*]T, ldab: isize, minor: ?*usize) PbtrfError!void {
    comptime chol.floatType(T, "pbtrf");
    if (uplo.isRev()) return error.InvalidUplo;
    if (n < 0) return error.InvalidN;
    if (kd < 0) return error.InvalidKd;
    if (ldab < kd + 1) return error.InvalidLdab;
    if (n == 0) return;
    if (!util.fitsMatrix(kd + 1, n, ldab, @sizeOf(T))) return error.InvalidLdab;
    const un: usize = @intCast(n);
    const ukd: usize = @intCast(kd);
    const n1 = @min(split(T, un), ukd);
    // zero, as `?laset` leaves the triangle the copies never write
    const work = try allocator.alloc(T, n1 * n1);
    defer allocator.free(work);
    @memset(work, scalar(T, 0));
    const k = try pbtrfRec(allocator, io, T, uplo, un, ukd, ab, @intCast(ldab), work.ptr, @max(1, n1));
    if (k == 0) return;
    if (minor) |m| m.* = k;
    return error.NotPositiveDefinite;
}

/// `RELAPACK_?pbtrf_rec`: the first minor that is not positive definite, or 0.
fn pbtrfRec(al: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, kd: usize, ab: [*]T, ldab: usize, work: [*]T, ldw: usize) Error!usize {
    if (n <= crossover or ldab == 1) return pbtf2(T, uplo == .u, n, kd, ab, ldab);
    const one = scalar(T, 1);
    const minus = scalar(T, -1);
    const ct = ctrans(T);
    // the band as a dense matrix of leading dimension `ldab - 1`
    const ld = ldab - 1;
    const a = ab + (if (uplo == .l) 0 else kd);
    const n1 = @min(split(T, n), kd);
    const n2 = n - n1;
    const k = try potrfRec(al, io, T, uplo, n1, a, ld);
    if (k != 0) return k;
    const n21 = @min(n2, kd - n1);
    const n22 = @min(n2 - n21, n1);
    const a_br = a + ld * n1 + n1;
    if (uplo == .l) {
        const blt = a + n1;
        const blb = blt + n21;
        // A_BLt = A_BLt / A_TL'; A_BRtl -= A_BLt A_BLt'
        try trsm(al, io, T, true, .r, .l, ct, .n, n21, n1, one, a, ld, blt, ld);
        try syrk(al, io, T, true, .l, .n, n21, n1, minus, blt, ld, null, 0, one, a_br, ld);
        // the triangle of A_BLb in the band, through the workspace
        copyTri(T, true, n22, n1, blb, ld, work, ldw);
        try trsm(al, io, T, true, .r, .l, ct, .n, n22, n1, one, a, ld, work, ldw);
        try gemm(al, io, T, .n, ct, n22, n21, n1, minus, work, ldw, blt, ld, one, a_br + n21, ld);
        try syrk(al, io, T, true, .l, .n, n22, n1, minus, work, ldw, null, 0, one, a_br + ld * n21 + n21, ld);
        copyTri(T, true, n22, n1, work, ldw, blb, ld);
    } else {
        const trl = a + ld * n1;
        const trr = trl + ld * n21;
        // A_TRl = A_TL' \ A_TRl; A_BRtl -= A_TRl' A_TRl
        try trsm(al, io, T, true, .l, .u, ct, .n, n1, n21, one, a, ld, trl, ld);
        try syrk(al, io, T, true, .u, ct, n21, n1, minus, trl, ld, null, 0, one, a_br, ld);
        copyTri(T, false, n1, n22, trr, ld, work, ldw);
        try trsm(al, io, T, true, .l, .u, ct, .n, n1, n22, one, a, ld, work, ldw);
        try gemm(al, io, T, ct, .n, n21, n22, n1, minus, trl, ld, work, ldw, one, a_br + ld * n21, ld);
        try syrk(al, io, T, true, .u, ct, n22, n1, minus, work, ldw, null, 0, one, a_br + ld * n21 + n21, ld);
        copyTri(T, false, n1, n22, work, ldw, trr, ld);
    }
    // the rest, dense where it lies inside the band, else banded
    const k2 = if (kd > n1 and n2 <= kd + 1)
        try potrfRec(al, io, T, uplo, n2, a_br, ld)
    else
        try pbtrfRec(al, io, T, uplo, n2, kd, ab + ldab * n1, ldab, work, ldw);
    return if (k2 != 0) k2 + n1 else 0;
}

/// `?LACPY` of the upper (`upper`) or lower triangle of `m x n`.
fn copyTri(comptime T: type, upper: bool, m: usize, n: usize, src: [*]const T, lds: usize, dst: [*]T, ldd: usize) void {
    for (0..n) |j| {
        const lo = if (upper) 0 else @min(j, m);
        const hi = if (upper) @min(j + 1, m) else m;
        for (lo..hi) |i| dst[i + j * ldd] = src[i + j * lds];
    }
}

/// Reference-LAPACK's `xPBTF2` (`lapack-netlib/SRC/dpbtf2.f`, `zpbtf2.f`;
/// OpenBLAS v0.3.34-317): column by column, `xSYR` or `xHER` on the
/// `kd`-wide window below (`u`: right of) the pivot, in the reference
/// BLAS's order. `upper` runs the lower algorithm on the transpose (the
/// reference conjugates the row it takes). A pivot that is not positive
/// fails, NaN too, as `potf2`'s (the reference's `AJJ.LE.ZERO` passes NaN).
fn pbtf2(comptime T: type, upper: bool, n: usize, kd: usize, ab: [*]T, ldab: usize) usize {
    const herm = comptime isComplex(T);
    const a: Tview(T) = .{ .p = ab + (if (upper) kd else 0), .ld = ldab - 1, .t = upper };
    for (0..n) |j| {
        const ajj0 = chol.re(T, a.at(j, j).*);
        if (!(ajj0 > 0)) {
            a.at(j, j).* = scalar(T, ajj0);
            return j + 1;
        }
        const ajj = @sqrt(ajj0);
        a.at(j, j).* = scalar(T, ajj);
        const kn = @min(kd, n - j - 1);
        if (kn == 0) continue;
        const r = 1 / ajj;
        for (j + 1..j + 1 + kn) |i| a.at(i, j).* = chol.mulReal(T, a.at(i, j).*, r);
        // `xSYR`, `xHER` with -1 on the window
        for (j + 1..j + 1 + kn) |c| {
            const x = a.at(c, j).*;
            if (chol.isZero(T, x)) {
                if (comptime herm) a.at(c, c).* = scalar(T, a.at(c, c).re);
                continue;
            }
            const t = chol.neg(T, chol.conj(T, x));
            if (comptime herm) a.at(c, c).* = scalar(T, a.at(c, c).re + chol.mul(T, x, t).re);
            for (if (herm) c + 1 else c..j + 1 + kn) |i| a.at(i, c).* = chol.add(T, a.at(i, c).*, chol.mul(T, a.at(i, j).*, t));
        }
    }
    return 0;
}

// ---- gbtrf -----------------------------------------------------------------------------

/// What `gbtrf` returns: bad arguments in the reference's order, `Singular`
/// (an exactly zero `U(i, i)`, the factorization completed), or
/// `error.OutOfMemory` and `error.Canceled`.
pub const GbtrfError = error{ InvalidM, InvalidN, InvalidKl, InvalidKu, InvalidLdab, Singular } || Allocator.Error || Io.Cancelable;

/// `?GBTRF`: `A = P L U` of the `m x n` band matrix with `kl` sub- and `ku`
/// superdiagonals in LAPACK's band storage `ab` (`ldab >= 2 kl + ku + 1`,
/// `A(i, j)` at row `kl + ku + i - j` of column `j`; the first `kl` rows
/// take `U`'s fill-in), by partial pivoting; `ipiv` and the multipliers as
/// LAPACK stores them, each interchange applied to the columns right of its
/// step only. An exactly zero `U(i, i)` returns `error.Singular` with the
/// first `i` (1-based) in `pivot`, the factorization completed.
///
/// ReLAPACK's (`relapack/src/dgbtrf.c` and the rest; OpenBLAS v0.3.34-317)
/// runs the reference's `xGBTF2` (`gbtf2`) unless `n <= kl`, where its
/// recursion gives a wrong factor, as upstream's full recursion does on most
/// shapes (`issues.md`); so `gbtrf` is `gbtf2` throughout. `allocator` and
/// `io` are unused, kept for the blocked form to come.
pub fn gbtrf(allocator: Allocator, io: Io, comptime T: type, m: isize, n: isize, kl: isize, ku: isize, ab: [*]T, ldab: isize, ipiv: [*]isize, pivot: ?*usize) GbtrfError!void {
    comptime chol.floatType(T, "gbtrf");
    _ = allocator;
    _ = io;
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    if (kl < 0) return error.InvalidKl;
    if (ku < 0) return error.InvalidKu;
    if (ldab < 2 * kl + ku + 1) return error.InvalidLdab;
    if (m == 0 or n == 0) return;
    if (!util.fitsMatrix(2 * kl + ku + 1, n, ldab, @sizeOf(T))) return error.InvalidLdab;
    const info = gbtf2(T, @intCast(m), @intCast(n), @intCast(kl), @intCast(ku), ab, @intCast(ldab), ipiv);
    return lu.singular(info, pivot);
}

/// Reference-LAPACK's `xGBTF2` (`lapack-netlib/SRC/dgbtf2.f`, `zgbtf2.f`;
/// OpenBLAS v0.3.34-317): column by column, the pivot the first NaN or else
/// the first largest `|re| + |im|` (rstd's `iamax`), the interchange over the
/// columns the step reaches (`JU`), the multipliers times the pivot's
/// reciprocal, and `xGER` (`xGERU`) in the reference's order, skipping a
/// zero `y(j)`. The first zero pivot (1-based), or 0.
fn gbtf2(comptime T: type, m: usize, n: usize, kl: usize, ku: usize, ab: [*]T, ldab: usize, ipiv: [*]isize) usize {
    const F = Float(T);
    const kv = ku + kl;
    const zero = scalar(T, 0);
    // `AB(r, c)`, a place of the storage
    const at = struct {
        inline fn f(p: [*]T, ld: usize, r: usize, c: usize) *T {
            return &p[r + c * ld];
        }
    }.f;
    // `A(r, c)` of the matrix, `c - kv <= r`
    const el = struct {
        inline fn f(p: [*]T, ld: usize, v: usize, r: usize, c: usize) *T {
            return &p[v + r - c + c * ld];
        }
    }.f;
    // fill-in of columns `ku + 1` to `kv - 1` (0-based), zeroed
    for (ku + 1..@max(ku + 1, @min(kv, n))) |j| for (kv - j..kl) |i| {
        at(ab, ldab, i, j).* = zero;
    };
    var info: usize = 0;
    var ju: usize = 0;
    for (0..@min(m, n)) |j| {
        if (j + kv < n) for (0..kl) |i| {
            at(ab, ldab, i, j + kv).* = zero;
        };
        const km = @min(kl, m - j - 1);
        // the pivot among the `km + 1` rows from the diagonal
        var jp: usize = 0;
        var big: F = -1;
        for (0..km + 1) |i| {
            const v = at(ab, ldab, kv + i, j).*;
            const x = if (comptime isComplex(T)) @abs(v.re) + @abs(v.im) else @abs(v);
            if (std.math.isNan(x)) {
                jp = i;
                break;
            }
            if (x > big) {
                big = x;
                jp = i;
            }
        }
        ipiv[j] = @intCast(jp + j + 1);
        const piv = at(ab, ldab, kv + jp, j).*;
        if (chol.isZero(T, piv)) {
            if (info == 0) info = j + 1;
            continue;
        }
        ju = @max(ju, @min(j + ku + jp, n - 1));
        // rows `j` and `j + jp` over columns `j..ju`, along the band's diagonals
        if (jp != 0) for (j..ju + 1) |c| {
            std.mem.swap(T, el(ab, ldab, kv, j + jp, c), el(ab, ldab, kv, j, c));
        };
        if (km == 0) continue;
        const r = recipOf(T, at(ab, ldab, kv, j).*);
        for (1..km + 1) |i| at(ab, ldab, kv + i, j).* = chol.mul(T, at(ab, ldab, kv + i, j).*, r);
        // A(j+1:, j+1:ju) -= l u, `xGER` with -1
        for (j + 1..ju + 1) |c| {
            const y = el(ab, ldab, kv, j, c).*;
            if (chol.isZero(T, y)) continue;
            const t = chol.neg(T, y);
            for (1..km + 1) |i| {
                const dst = el(ab, ldab, kv, j + i, c);
                dst.* = chol.add(T, dst.*, chol.mul(T, at(ab, ldab, kv + i, j).*, t));
            }
        }
    }
    return info;
}

/// `ONE / p` as Fortran divides (complex: Smith's, `util.cdiv`).
fn recipOf(comptime T: type, p: T) T {
    if (comptime !isComplex(T)) return 1 / p;
    const z = util.cdiv(Float(T), 1, 0, p.re, p.im);
    return .init(z[0], z[1]);
}

// ---- trsyl -----------------------------------------------------------------------------

/// What `trsyl` returns: bad arguments in the reference's order;
/// `Perturbed` (the reference's `INFO = 1`: `A` and `-isgn B` have close or
/// common eigenvalues, so perturbed values were used), with `C` and `scale`
/// written; `error.OutOfMemory` or `error.Canceled`.
pub const TrsylError = error{ InvalidTrana, InvalidTranb, InvalidIsgn, InvalidM, InvalidN, InvalidLda, InvalidLdb, InvalidLdc, Perturbed } || Allocator.Error || Io.Cancelable;

/// `?TRSYL` by ReLAPACK's recursion (`relapack/src/dtrsyl.c`, `ztrsyl.c`
/// and the rest; OpenBLAS v0.3.34-317): `op(A) X + isgn X op(B) = scale C`
/// for `X`, into `C` (`m x n`), `A` (`m x m`) and `B` (`n x n`) upper
/// triangular, real ones in Schur form (2 x 2 diagonal blocks); `op` is `n`
/// or the (conjugate) transpose, `t` or `c` for real `T` and `c` only for
/// complex, as the reference. `scale <= 1` keeps `X` from overflowing.
///
/// The larger side splits in two (never inside a 2 x 2 block), the halves
/// joined by `gemm`; the base case is `sylvester.trsy2`, the reference's
/// `xTRSYL`, which reads the whole of `A` and `B` (`xLANGE`).
///
/// The reference wins where ReLAPACK, solving the top half first (`op(A)`
/// transposed), rescales the bottom half by the second solve's `scale`; the
/// top half, solved first, takes it here.
pub fn trsyl(allocator: Allocator, io: Io, comptime T: type, trana: Trans, tranb: Trans, isgn: isize, m: isize, n: isize, a: [*]const T, lda: isize, b: [*]const T, ldb: isize, c: [*]T, ldc: isize, scale: *Float(T)) TrsylError!void {
    comptime chol.floatType(T, "trsyl");
    const ta = sylTrans(T, trana) orelse return error.InvalidTrana;
    const tb = sylTrans(T, tranb) orelse return error.InvalidTranb;
    if (isgn != 1 and isgn != -1) return error.InvalidIsgn;
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, m)) return error.InvalidLda;
    if (ldb < @max(1, n)) return error.InvalidLdb;
    if (ldc < @max(1, m)) return error.InvalidLdc;
    scale.* = 1;
    if (m == 0 or n == 0) return;
    if (!util.fitsMatrix(m, m, lda, @sizeOf(T))) return error.InvalidLda;
    if (!util.fitsMatrix(n, n, ldb, @sizeOf(T))) return error.InvalidLdb;
    if (!util.fitsMatrix(m, n, ldc, @sizeOf(T))) return error.InvalidLdc;
    const r = try trsylRec(allocator, io, T, ta, tb, isgn, @intCast(m), @intCast(n), a, @intCast(lda), b, @intCast(ldb), c, @intCast(ldc));
    scale.* = r.scale;
    if (r.perturbed) return error.Perturbed;
}

/// Whether `t` transposes `A` or `B` for `trsyl`, or null where the
/// reference rejects it (complex `t`).
fn sylTrans(comptime T: type, t: Trans) ?bool {
    return switch (t) {
        .n => false,
        .c => true,
        .t => if (comptime isComplex(T)) null else true,
        .r, .o, .u, .s, .d => null,
    };
}

/// `RELAPACK_?trsyl_rec`.
fn trsylRec(al: Allocator, io: Io, comptime T: type, ta: bool, tb: bool, isgn: isize, m: usize, n: usize, a: [*]const T, lda: usize, b: [*]const T, ldb: usize, c: [*]T, ldc: usize) Error!sylvester.Solved(Float(T)) {
    if (m <= crossover and n <= crossover) return sylvester.trsy2(T, ta, tb, isgn, m, n, a, lda, b, ldb, c, ldc);
    const minus = scalar(T, -1);
    const msgn = scalar(T, @floatFromInt(-isgn));
    const opa: Trans = if (ta) ctrans(T) else .n;
    const opb: Trans = if (tb) ctrans(T) else .n;
    var r1: sylvester.Solved(Float(T)) = undefined;
    var r2: sylvester.Solved(Float(T)) = undefined;
    if (m > n) {
        var m1 = split(T, m);
        // never inside a 2 x 2 block
        if (comptime !isComplex(T)) if (a[m1 + lda * (m1 - 1)] != 0) {
            m1 += 1;
        };
        const m2 = m - m1;
        const a_tr = a + lda * m1;
        const a_br = a_tr + m1;
        const c_b = c + m1;
        if (!ta) {
            // C_T = scale1 C_T - A_TR C_B
            r1 = try trsylRec(al, io, T, ta, tb, isgn, m2, n, a_br, lda, b, ldb, c_b, ldc);
            try gemm(al, io, T, .n, .n, m1, n, m2, minus, a_tr, lda, c_b, ldc, scalar(T, r1.scale), c, ldc);
            r2 = try trsylRec(al, io, T, ta, tb, isgn, m1, n, a, lda, b, ldb, c, ldc);
            if (r2.scale != 1) scaleBlock(T, m2, n, r2.scale, c_b, ldc);
        } else {
            // C_B = scale1 C_B - A_TR' C_T
            r1 = try trsylRec(al, io, T, ta, tb, isgn, m1, n, a, lda, b, ldb, c, ldc);
            try gemm(al, io, T, opa, .n, m2, n, m1, minus, a_tr, lda, c, ldc, scalar(T, r1.scale), c_b, ldc);
            r2 = try trsylRec(al, io, T, ta, tb, isgn, m2, n, a_br, lda, b, ldb, c_b, ldc);
            if (r2.scale != 1) scaleBlock(T, m1, n, r2.scale, c, ldc);
        }
    } else {
        var n1 = split(T, n);
        if (comptime !isComplex(T)) if (b[n1 + ldb * (n1 - 1)] != 0) {
            n1 += 1;
        };
        const n2 = n - n1;
        const b_tr = b + ldb * n1;
        const b_br = b_tr + n1;
        const c_r = c + ldc * n1;
        if (!tb) {
            // C_R = scale1 C_R -/+ C_L B_TR
            r1 = try trsylRec(al, io, T, ta, tb, isgn, m, n1, a, lda, b, ldb, c, ldc);
            try gemm(al, io, T, .n, .n, m, n2, n1, msgn, c, ldc, b_tr, ldb, scalar(T, r1.scale), c_r, ldc);
            r2 = try trsylRec(al, io, T, ta, tb, isgn, m, n2, a, lda, b_br, ldb, c_r, ldc);
            if (r2.scale != 1) scaleBlock(T, m, n1, r2.scale, c, ldc);
        } else {
            // C_L = scale1 C_L -/+ C_R B_TR'
            r1 = try trsylRec(al, io, T, ta, tb, isgn, m, n2, a, lda, b_br, ldb, c_r, ldc);
            try gemm(al, io, T, .n, opb, m, n1, n2, msgn, c_r, ldc, b_tr, ldb, scalar(T, r1.scale), c, ldc);
            r2 = try trsylRec(al, io, T, ta, tb, isgn, m, n1, a, lda, b, ldb, c, ldc);
            if (r2.scale != 1) scaleBlock(T, m, n2, r2.scale, c_r, ldc);
        }
    }
    return .{ .scale = r1.scale * r2.scale, .perturbed = r1.perturbed or r2.perturbed };
}

/// `?LASCL('G', 1, s)`: `C *= s`, part by part.
fn scaleBlock(comptime T: type, m: usize, n: usize, s: Float(T), c: [*]T, ldc: usize) void {
    for (0..n) |j| for (0..m) |i| {
        c[i + j * ldc] = chol.mulReal(T, c[i + j * ldc], s);
    };
}

// ---- tgsyl -----------------------------------------------------------------------------

/// What `tgsyl` returns: bad arguments in the reference's order (`m` or `n`
/// below 1 is invalid, as there); `Perturbed` (the reference's `INFO > 0`:
/// `(A, D)` and `(B, E)` have common or close eigenvalues, so perturbed
/// values were used), with `C`, `F`, `scale` and `dif` written;
/// `error.OutOfMemory` or `error.Canceled`.
pub const TgsylError = error{ InvalidTrans, InvalidIjob, InvalidM, InvalidN, InvalidLda, InvalidLdb, InvalidLdc, InvalidLdd, InvalidLde, InvalidLdf, Perturbed } || Allocator.Error || Io.Cancelable;

/// `?TGSYL` by ReLAPACK's recursion (`relapack/src/dtgsyl.c`, `ztgsyl.c`
/// and the rest; OpenBLAS v0.3.34-317): the generalized Sylvester equation
/// `A R - L B = scale C`, `D R - L E = scale F` for `R` and `L` (`m x n`),
/// into `C` and `F`; `(A, D)` (`m x m`) and `(B, E)` (`n x n`) in
/// generalized Schur form, `A`, `B` upper quasi-triangular (complex:
/// triangular), `D`, `E` upper triangular. `trans` `t` (real) or `c`
/// (complex) solves the transposed system, `Aᴴ R + Dᴴ L = scale C`,
/// `R Bᴴ + L Eᴴ = -scale F`. `scale <= 1` keeps `R` and `L` from
/// overflowing.
///
/// `ijob` (read only for `trans` `n`): 0 solves; 3 and 4 only estimate
/// `Dif[(A, D), (B, E)]`, the smallest singular value of the system's
/// matrix, into `dif` (an upper bound on it), by the look-ahead of
/// `xLATDF` (3) or its `xGECON` null vectors (4), leaving that estimate's
/// `R` and `L` in `C` and `F`; 1 and 2 do both, the solution in `C` and
/// `F`. `dif` is not written for `ijob` 0 or `trans`.
///
/// The larger side splits in two (never inside a 2 x 2 block), the halves
/// joined by `gemm`; the base case is `sylvester.tgsy2`, the reference's
/// `xTGSY2`. Where ReLAPACK differs from the reference, the reference wins
/// ([issues](../../../issues.md#relapack-in-openblas)): the
/// transposed updates read `F` and `E` with their own leading dimensions,
/// and the real `dif` from `ijob` 2 and 4 counts every base case's block
/// pairs, not the last one's.
pub fn tgsyl(allocator: Allocator, io: Io, comptime T: type, trans: Trans, ijob: isize, m: isize, n: isize, a: [*]const T, lda: isize, b: [*]const T, ldb: isize, c: [*]T, ldc: isize, d: [*]const T, ldd: isize, e: [*]const T, lde: isize, f: [*]T, ldf: isize, scale: *Float(T), dif: *Float(T)) TgsylError!void {
    comptime chol.floatType(T, "tgsyl");
    const F = Float(T);
    // `N`, or `T` for real and `C` for complex, as `dtgsyl.f` and `ztgsyl.f`
    const tr = switch (trans) {
        .n => false,
        .t => if (comptime isComplex(T)) return error.InvalidTrans else true,
        .c => if (comptime isComplex(T)) true else return error.InvalidTrans,
        .r, .o, .u, .s, .d => return error.InvalidTrans,
    };
    if (!tr and (ijob < 0 or ijob > 4)) return error.InvalidIjob;
    if (m <= 0) return error.InvalidM;
    if (n <= 0) return error.InvalidN;
    if (lda < m) return error.InvalidLda;
    if (ldb < n) return error.InvalidLdb;
    if (ldc < m) return error.InvalidLdc;
    if (ldd < m) return error.InvalidLdd;
    if (lde < n) return error.InvalidLde;
    if (ldf < m) return error.InvalidLdf;
    inline for (.{ .{ m, lda, error.InvalidLda }, .{ n, ldb, error.InvalidLdb }, .{ m, ldd, error.InvalidLdd }, .{ n, lde, error.InvalidLde } }) |x| {
        if (!util.fitsMatrix(x[0], x[0], x[1], @sizeOf(T))) return x[2];
    }
    if (!util.fitsMatrix(m, n, ldc, @sizeOf(T))) return error.InvalidLdc;
    if (!util.fitsMatrix(m, n, ldf, @sizeOf(T))) return error.InvalidLdf;
    const um: usize = @intCast(m);
    const un: usize = @intCast(n);
    const p: Gsyl(T) = .{ .tr = tr, .a = a, .lda = @intCast(lda), .b = b, .ldb = @intCast(ldb), .c = c, .ldc = @intCast(ldc), .d = d, .ldd = @intCast(ldd), .e = e, .lde = @intCast(lde), .f = f, .ldf = @intCast(ldf) };
    const job: usize = if (tr) 0 else @intCast(ijob);
    // `ijob` 1 and 2 solve, then estimate on zeros; 3 and 4 only estimate
    var ifunc: u2 = if (job >= 3) @intCast(job - 2) else 0;
    const isolve: usize = if (job == 1 or job == 2) 2 else 1;
    if (job >= 3) zeroBoth(T, p, um, un);
    const saved = try allocator.alloc(T, if (isolve == 2) 2 * um * un else 0);
    defer allocator.free(saved);
    var perturbed = false;
    var scale2: F = 1;
    for (0..isolve) |round| {
        var acc: sylvester.Dif(F) = .{};
        const r = try tgsylRec(allocator, io, T, ifunc, um, un, p, &acc);
        scale.* = r.scale;
        if (r.perturbed) perturbed = true;
        if (acc.scale != 0) {
            const num: F = @floatFromInt(if (job == 1 or job == 3) 2 * um * un else acc.pq);
            dif.* = @sqrt(num) / (acc.scale * @sqrt(acc.sum));
        }
        if (isolve == 2 and round == 0) {
            ifunc = @intCast(job);
            scale2 = scale.*;
            copyBlock(T, um, un, c, p.ldc, saved.ptr, um);
            copyBlock(T, um, un, f, p.ldf, saved.ptr + um * un, um);
            zeroBoth(T, p, um, un);
        } else if (isolve == 2) {
            copyBlock(T, um, un, saved.ptr, um, c, p.ldc);
            copyBlock(T, um, un, saved.ptr + um * un, um, f, p.ldf);
            scale.* = scale2;
        }
    }
    if (perturbed) return error.Perturbed;
}

/// `tgsyl`'s six matrices, as one recursion level sees them.
fn Gsyl(comptime T: type) type {
    return struct {
        tr: bool,
        a: [*]const T,
        lda: usize,
        b: [*]const T,
        ldb: usize,
        c: [*]T,
        ldc: usize,
        d: [*]const T,
        ldd: usize,
        e: [*]const T,
        lde: usize,
        f: [*]T,
        ldf: usize,

        /// Rows `i..` of `C`, `F` and the `A`, `D` blocks from `(i, i)`.
        fn rows(s: @This(), i: usize) @This() {
            var x = s;
            x.a += i * s.lda + i;
            x.d += i * s.ldd + i;
            x.c += i;
            x.f += i;
            return x;
        }

        /// Columns `j..` of `C`, `F` and the `B`, `E` blocks from `(j, j)`.
        fn cols(s: @This(), j: usize) @This() {
            var x = s;
            x.b += j * s.ldb + j;
            x.e += j * s.lde + j;
            x.c += j * s.ldc;
            x.f += j * s.ldf;
            return x;
        }
    };
}

fn zeroBoth(comptime T: type, p: Gsyl(T), m: usize, n: usize) void {
    for (0..n) |j| for (0..m) |i| {
        p.c[i + j * p.ldc] = scalar(T, 0);
        p.f[i + j * p.ldf] = scalar(T, 0);
    };
}

/// `?LACPY('F')`.
fn copyBlock(comptime T: type, m: usize, n: usize, src: [*]const T, lds: usize, dst: [*]T, ldd: usize) void {
    for (0..n) |j| @memcpy(dst[j * ldd ..][0..m], src[j * lds ..][0..m]);
}

/// `RELAPACK_?tgsyl_rec`: the halves of the larger side, `dif`'s sums
/// carried through.
fn tgsylRec(al: Allocator, io: Io, comptime T: type, ifunc: u2, m: usize, n: usize, p: Gsyl(T), acc: *sylvester.Dif(Float(T))) Error!sylvester.Solved(Float(T)) {
    if (m <= crossover and n <= crossover) return sylvester.tgsy2(al, io, T, p.tr, ifunc, m, n, p.a, p.lda, p.b, p.ldb, p.c, p.ldc, p.d, p.ldd, p.e, p.lde, p.f, p.ldf, acc);
    const one = scalar(T, 1);
    const minus = scalar(T, -1);
    const ct = ctrans(T);
    var r1: sylvester.Solved(Float(T)) = undefined;
    var r2: sylvester.Solved(Float(T)) = undefined;
    if (m > n) {
        var m1 = split(T, m);
        // never inside a 2 x 2 block
        if (comptime !isComplex(T)) if (p.a[m1 + p.lda * (m1 - 1)] != 0) {
            m1 += 1;
        };
        const m2 = m - m1;
        const bot = p.rows(m1);
        const a_tr = p.a + p.lda * m1;
        const d_tr = p.d + p.ldd * m1;
        if (!p.tr) {
            // C_T = scale1 C_T - A_TR C_B, F_T = scale1 F_T - D_TR C_B
            r1 = try tgsylRec(al, io, T, ifunc, m2, n, bot, acc);
            const s1 = scalar(T, r1.scale);
            try gemm(al, io, T, .n, .n, m1, n, m2, minus, a_tr, p.lda, bot.c, p.ldc, s1, p.c, p.ldc);
            try gemm(al, io, T, .n, .n, m1, n, m2, minus, d_tr, p.ldd, bot.c, p.ldc, s1, p.f, p.ldf);
            r2 = try tgsylRec(al, io, T, ifunc, m1, n, p, acc);
            if (r2.scale != 1) {
                scaleBlock(T, m2, n, r2.scale, bot.c, p.ldc);
                scaleBlock(T, m2, n, r2.scale, bot.f, p.ldf);
            }
        } else {
            // C_B = scale1 C_B - A_TRᴴ C_T - D_TRᴴ F_T, F_B = scale1 F_B
            r1 = try tgsylRec(al, io, T, ifunc, m1, n, p, acc);
            if (r1.scale != 1) scaleBlock(T, m2, n, r1.scale, bot.f, p.ldf);
            try gemm(al, io, T, ct, .n, m2, n, m1, minus, a_tr, p.lda, p.c, p.ldc, scalar(T, r1.scale), bot.c, p.ldc);
            try gemm(al, io, T, ct, .n, m2, n, m1, minus, d_tr, p.ldd, p.f, p.ldf, one, bot.c, p.ldc);
            r2 = try tgsylRec(al, io, T, ifunc, m2, n, bot, acc);
            if (r2.scale != 1) {
                scaleBlock(T, m1, n, r2.scale, p.c, p.ldc);
                scaleBlock(T, m1, n, r2.scale, p.f, p.ldf);
            }
        }
    } else {
        var n1 = split(T, n);
        if (comptime !isComplex(T)) if (p.b[n1 + p.ldb * (n1 - 1)] != 0) {
            n1 += 1;
        };
        const n2 = n - n1;
        const right = p.cols(n1);
        const b_tr = p.b + p.ldb * n1;
        const e_tr = p.e + p.lde * n1;
        if (!p.tr) {
            // C_R = scale1 C_R + F_L B_TR, F_R = scale1 F_R + F_L E_TR
            r1 = try tgsylRec(al, io, T, ifunc, m, n1, p, acc);
            const s1 = scalar(T, r1.scale);
            try gemm(al, io, T, .n, .n, m, n2, n1, one, p.f, p.ldf, b_tr, p.ldb, s1, right.c, p.ldc);
            try gemm(al, io, T, .n, .n, m, n2, n1, one, p.f, p.ldf, e_tr, p.lde, s1, right.f, p.ldf);
            r2 = try tgsylRec(al, io, T, ifunc, m, n2, right, acc);
            if (r2.scale != 1) {
                scaleBlock(T, m, n1, r2.scale, p.c, p.ldc);
                scaleBlock(T, m, n1, r2.scale, p.f, p.ldf);
            }
        } else {
            // F_L = scale1 F_L + C_R B_TRᴴ + F_R E_TRᴴ, C_L = scale1 C_L
            r1 = try tgsylRec(al, io, T, ifunc, m, n2, right, acc);
            if (r1.scale != 1) scaleBlock(T, m, n1, r1.scale, p.c, p.ldc);
            try gemm(al, io, T, .n, ct, m, n1, n2, one, right.c, p.ldc, b_tr, p.ldb, scalar(T, r1.scale), p.f, p.ldf);
            try gemm(al, io, T, .n, ct, m, n1, n2, one, right.f, p.ldf, e_tr, p.lde, one, p.f, p.ldf);
            r2 = try tgsylRec(al, io, T, ifunc, m, n1, p, acc);
            if (r2.scale != 1) {
                scaleBlock(T, m, n2, r2.scale, right.c, p.ldc);
                scaleBlock(T, m, n2, r2.scale, right.f, p.ldf);
            }
        }
    }
    return .{ .scale = r1.scale * r2.scale, .perturbed = r1.perturbed or r2.perturbed };
}

// ---- tests ---------------------------------------------------------------------

const testing = std.testing;
const blas = @import("../blas/__root__.zig");
const check = @import("../blas/_test/blas1/check.zig");
const Vtable = @import("../blas/Vtable.zig");
const options = @import("lapack_options");
const faer_svd = @import("svd.zig");

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
        potrf(al, check.io, T, uploOf(upper), @intCast(n), p.a.ptr, @intCast(p.lda), null) catch |e| {
            check.say("relapack potrf {s} n={d} upper={}\n", .{ @typeName(T), n, upper });
            return e;
        };
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

// ---- sytrf tests ------------------------------------------------------------------

const W = chol.W;

/// `sytrf`'s factor in `a` (`n x n`, `lda`) and `ipiv` multiplied out in
/// `W(T)`: `L D Lᵀ` or `U D Uᵀ` (`ᴴ` with `herm`), as LAPACK's `xLAVSY` and
/// `xLAVHE` read the format, applied to the identity; `mag` gets the same
/// products in magnitude (`|L| |D| |L|ᵀ`, permuted).
pub fn ldlProduct(comptime T: type, comptime herm: bool, comptime rook: bool, al: Allocator, upper: bool, n: usize, a: []const T, lda: usize, ipiv: []const isize, out: []W(T), mag: []f128) !void {
    // `P = S(p1) S(p2)` per step, in the order the factorization took them
    const Step = struct { b: usize, s: usize, p1: [2]usize, p2: ?[2]usize };
    var steps: std.ArrayList(Step) = .empty;
    defer steps.deinit(al);
    const at = struct {
        fn neg(v: isize) usize {
            return @intCast(-v - 1);
        }
    };
    if (!upper) {
        var k: usize = 0;
        while (k < n) {
            if (ipiv[k] > 0) {
                try steps.append(al, .{ .b = k, .s = 1, .p1 = .{ k, @intCast(ipiv[k] - 1) }, .p2 = null });
                k += 1;
            } else {
                try steps.append(al, if (rook)
                    .{ .b = k, .s = 2, .p1 = .{ k, at.neg(ipiv[k]) }, .p2 = .{ k + 1, at.neg(ipiv[k + 1]) } }
                else
                    .{ .b = k, .s = 2, .p1 = .{ k + 1, at.neg(ipiv[k]) }, .p2 = null });
                k += 2;
            }
        }
    } else {
        var k: usize = n;
        while (k > 0) {
            const j = k - 1;
            if (ipiv[j] > 0) {
                try steps.append(al, .{ .b = j, .s = 1, .p1 = .{ j, @intCast(ipiv[j] - 1) }, .p2 = null });
                k -= 1;
            } else {
                try steps.append(al, if (rook)
                    .{ .b = j - 1, .s = 2, .p1 = .{ j, at.neg(ipiv[j]) }, .p2 = .{ j - 1, at.neg(ipiv[j - 1]) } }
                else
                    .{ .b = j - 1, .s = 2, .p1 = .{ j - 1, at.neg(ipiv[j]) }, .p2 = null });
                k -= 2;
            }
        }
    }
    const F = check.Wide(T);
    const zero = check.wide(T, scalar(T, 0));
    const Ctx = struct {
        n: usize,
        x: []F,
        m: []f128,
        fn swap(c: @This(), r: usize, s: usize) void {
            if (r == s) return;
            for (0..c.n) |j| {
                std.mem.swap(F, &c.x[r + j * c.n], &c.x[s + j * c.n]);
                std.mem.swap(f128, &c.m[r + j * c.n], &c.m[s + j * c.n]);
            }
        }
    };
    for (out, mag, 0..) |*o, *g, i| {
        const one = i % (n + 1) == 0;
        o.* = if (one) check.wide(T, scalar(T, 1)) else zero;
        g.* = if (one) 1 else 0;
    }
    const c: Ctx = .{ .n = n, .x = out, .m = mag };
    const el = struct {
        fn of(v: T) F {
            return check.wide(T, v);
        }
        fn cjw(v: F) F {
            return if (comptime herm) check.conjW(v) else v;
        }
    };
    // the multiplier rows of step `st`
    const Rows = struct {
        fn lo(st: Step, up: bool) usize {
            return if (up) 0 else st.b + st.s;
        }
        fn hi(st: Step, up: bool, nn: usize) usize {
            return if (up) st.b else nn;
        }
    };
    // X := Lᵀ X (`ᴴ`): each step's `Pᵀ` (p1, then p2), then its multipliers
    for (steps.items) |st| {
        c.swap(st.p1[0], st.p1[1]);
        if (st.p2) |q| c.swap(q[0], q[1]);
        for (st.b..st.b + st.s) |col| for (Rows.lo(st, upper)..Rows.hi(st, upper, n)) |r| {
            const l = el.cjw(el.of(a[r + col * lda]));
            for (0..n) |j| {
                c.x[col + j * n] = check.addW(c.x[col + j * n], check.mulW(l, c.x[r + j * n]));
                c.m[col + j * n] += check.abs1(l) * c.m[r + j * n];
            }
        };
    }
    // X := D X
    for (steps.items) |st| {
        const b = st.b;
        if (st.s == 1) {
            const d = el.of(a[b + b * lda]);
            const dd = if (comptime herm) check.wide(T, scalar(T, a[b + b * lda].re)) else d;
            for (0..n) |j| {
                c.x[b + j * n] = check.mulW(dd, c.x[b + j * n]);
                c.m[b + j * n] *= check.abs1(dd);
            }
            continue;
        }
        const d11 = if (comptime herm) check.wide(T, scalar(T, a[b + b * lda].re)) else el.of(a[b + b * lda]);
        const d22 = if (comptime herm) check.wide(T, scalar(T, a[b + 1 + (b + 1) * lda].re)) else el.of(a[b + 1 + (b + 1) * lda]);
        const d21 = if (upper) el.cjw(el.of(a[b + (b + 1) * lda])) else el.of(a[b + 1 + b * lda]);
        const d12 = el.cjw(d21);
        for (0..n) |j| {
            const x0 = c.x[b + j * n];
            const x1 = c.x[b + 1 + j * n];
            const m0 = c.m[b + j * n];
            const m1 = c.m[b + 1 + j * n];
            c.x[b + j * n] = check.addW(check.mulW(d11, x0), check.mulW(d12, x1));
            c.x[b + 1 + j * n] = check.addW(check.mulW(d21, x0), check.mulW(d22, x1));
            c.m[b + j * n] = check.abs1(d11) * m0 + check.abs1(d12) * m1;
            c.m[b + 1 + j * n] = check.abs1(d21) * m0 + check.abs1(d22) * m1;
        }
    }
    // X := L X: steps in reverse, multipliers then `P` (p2, then p1)
    var i = steps.items.len;
    while (i > 0) {
        i -= 1;
        const st = steps.items[i];
        for (st.b..st.b + st.s) |col| for (Rows.lo(st, upper)..Rows.hi(st, upper, n)) |r| {
            const l = el.of(a[r + col * lda]);
            for (0..n) |j| {
                c.x[r + j * n] = check.addW(c.x[r + j * n], check.mulW(l, c.x[col + j * n]));
                c.m[r + j * n] += check.abs1(l) * c.m[col + j * n];
            }
        };
        if (st.p2) |q| c.swap(q[0], q[1]);
        c.swap(st.p1[0], st.p1[1]);
    }
}

/// A random symmetric (`herm`: Hermitian) `n x n` in its `upper` or lower
/// triangle at `lda = n + 2`, `zero_diag` with a zero diagonal (2 x 2
/// pivots); the other triangle, the padding and a Hermitian diagonal's
/// imaginary parts NaN. `w` keeps the whole matrix.
pub fn Indefinite(comptime T: type, comptime herm: bool) type {
    return struct {
        n: usize,
        lda: usize,
        upper: bool,
        a: []T,
        w: []W(T),
        ipiv: []isize,

        const Self = @This();

        pub fn init(al: Allocator, n: usize, upper: bool, zero_diag: bool, seed: u64) !Self {
            const lda = n + 2;
            const v = try al.alloc(T, n * n);
            defer al.free(v);
            check.fill(T, v, seed);
            const w = try al.alloc(W(T), n * n);
            errdefer al.free(w);
            for (0..n) |j| for (0..n) |i| {
                var x = check.wide(T, v[@max(i, j) + @min(i, j) * n]);
                if (i < j and herm) x = check.conjW(x);
                if (i == j) x = if (zero_diag) check.wide(T, scalar(T, 0)) else if (comptime herm) check.wide(T, scalar(T, v[i + i * n].re)) else x;
                w[i + j * n] = x;
            };
            const a = try al.alloc(T, @max(1, lda * n));
            errdefer al.free(a);
            const ipiv = try al.alloc(isize, @max(1, n));
            var self: Self = .{ .n = n, .lda = lda, .upper = upper, .a = a, .w = w, .ipiv = ipiv };
            self.store();
            return self;
        }

        pub fn store(self: *Self) void {
            @memset(self.a, chol.poisoned(T));
            for (0..self.n) |j| for (0..self.n) |i| {
                if (if (self.upper) i > j else i < j) continue;
                var x = chol.toT(T, self.w[i + j * self.n]);
                if (comptime herm) if (i == j) {
                    x.im = std.math.nan(Float(T));
                };
                self.a[i + j * self.lda] = x;
            };
            @memset(self.ipiv, 0);
        }

        pub fn deinit(self: *Self, al: Allocator) void {
            al.free(self.a);
            al.free(self.w);
            al.free(self.ipiv);
        }

        /// The factor against `A` within `32 (n + 2)` ulps of `|A| + |L| |D|
        /// |L|ᵀ`, a Hermitian `D`'s diagonal real, nothing outside the
        /// triangle written.
        pub fn expectFactor(self: *const Self, comptime rook: bool, al: Allocator, what: []const u8) !void {
            const n = self.n;
            const got = try al.alloc(W(T), n * n);
            defer al.free(got);
            const mag = try al.alloc(f128, n * n);
            defer al.free(mag);
            try ldlProduct(T, herm, rook, al, self.upper, n, self.a, self.lda, self.ipiv[0..n], got, mag);
            const eps: f128 = std.math.floatEps(Float(T));
            const ulps: f128 = @floatFromInt(32 * (n + 2));
            for (0..n) |j| for (0..n) |i| {
                const want = self.w[i + j * n];
                const d = check.abs1(check.subW(got[i + j * n], want));
                const allowed = ulps * eps * (check.abs1(want) + mag[i + j * n]) + std.math.floatMin(Float(T));
                if (!(d <= allowed)) {
                    check.say("{s} {s} {s} n={d}: (L D Lᵀ)({d},{d}) off by {e}, allowed {e}\n", .{ what, @typeName(T), if (self.upper) "u" else "l", n, i, j, @as(f64, @floatCast(d)), @as(f64, @floatCast(allowed)) });
                    return error.TestUnexpectedResult;
                }
            };
            for (0..n) |j| for (0..self.lda) |i| {
                const x = self.a[i + j * self.lda];
                const inside = i < n and (if (self.upper) i <= j else i >= j);
                if (!inside and !isNanT(T, x)) {
                    check.say("{s} {s} n={d}: wrote ({d},{d}), outside its triangle\n", .{ what, @typeName(T), n, i, j });
                    return error.TestUnexpectedResult;
                }
                if (comptime herm) if (i == j and i < n and x.im != 0) {
                    check.say("{s} {s} n={d}: D({d},{d}) has imaginary part {e}\n", .{ what, @typeName(T), n, i, i, x.im });
                    return error.TestUnexpectedResult;
                };
            };
        }
    };
}

/// `sytrf`, `hetrf` (`herm`) or a `_rook` form, by its flags.
pub fn indefiniteRoutine(comptime T: type, comptime herm: bool, comptime rook: bool) fn (Allocator, Io, comptime type, Uplo, isize, [*]T, isize, [*]isize, ?*usize) SytrfError!void {
    return if (herm) (if (rook) hetrfRook else hetrf) else (if (rook) sytrfRook else sytrf);
}

/// Around the panel width (`sytrf_nb`), past it the blocked loop.
const sytrf_sizes = [_]usize{ 1, 2, 3, 24, 57, 65, 66, 100 };

fn sytrfCases(comptime T: type, comptime herm: bool, comptime rook: bool) !void {
    const al = testing.allocator;
    const f = comptime indefiniteRoutine(T, herm, rook);
    const what = if (herm) (if (rook) "hetrfRook" else "hetrf") else (if (rook) "sytrfRook" else "sytrf");
    for ([_]bool{ true, false }) |upper| for ([_]bool{ false, true }) |zd| for (sytrf_sizes) |n| {
        // a zero 1 x 1 is singular
        if ((check.rig() and n > 66) or (zd and n == 1)) continue;
        var p = try Indefinite(T, herm).init(al, n, upper, zd, n * 7 + @intFromBool(upper) + 2 * @as(u64, @intFromBool(zd)));
        defer p.deinit(al);
        f(al, check.io, T, uploOf(upper), @intCast(n), p.a.ptr, @intCast(p.lda), p.ipiv.ptr, null) catch |e| {
            check.say("{s} {s} n={d}\n", .{ what, @typeName(T), n });
            return e;
        };
        try p.expectFactor(rook, al, what);
    };
}

test "relapack sytrf, sytrfRook: A = L D Lᵀ and U D Uᵀ, random and zero-diagonal, every type" {
    inline for (test_types) |T| inline for (.{ false, true }) |rook| try sytrfCases(T, false, rook);
}

test "relapack hetrf, hetrfRook: A = L D Lᴴ and U D Uᴴ, random and zero-diagonal, complex types" {
    inline for (test_types) |T| if (comptime isComplex(T)) inline for (.{ false, true }) |rook| try sytrfCases(T, true, rook);
}

test "relapack sytrf and the rest: three panels, f64 and c64" {
    if (check.rig()) return error.SkipZigTest;
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64 }) |T| inline for (.{ false, true }) |herm| inline for (.{ false, true }) |rook| {
        if (comptime !herm or isComplex(T)) {
            const f = comptime indefiniteRoutine(T, herm, rook);
            for ([_]bool{ true, false }) |upper| {
                var p = try Indefinite(T, herm).init(al, 150, upper, false, 99);
                defer p.deinit(al);
                try f(al, check.io, T, uploOf(upper), 150, p.a.ptr, @intCast(p.lda), p.ipiv.ptr, null);
                try p.expectFactor(rook, al, "sytrf 150");
            }
        }
    };
}

test "relapack sytrf and the rest: a zero row and column is singular there, the factorization completed" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64 }) |T| inline for (.{ false, true }) |herm| inline for (.{ false, true }) |rook| {
        if (comptime !herm or isComplex(T)) {
            const f = comptime indefiniteRoutine(T, herm, rook);
            for ([_]bool{ true, false }) |upper| for ([_]usize{ 0, 5, 70, 99 }) |k| {
                const n = 100;
                var p = try Indefinite(T, herm).init(al, n, upper, false, 3 + k);
                defer p.deinit(al);
                for (0..n) |i| {
                    p.w[i + k * n] = check.wide(T, scalar(T, 0));
                    p.w[k + i * n] = check.wide(T, scalar(T, 0));
                }
                p.store();
                var piv: usize = 0;
                try testing.expectError(error.Singular, f(al, check.io, T, uploOf(upper), n, p.a.ptr, @intCast(p.lda), p.ipiv.ptr, &piv));
                // the zero row may move with the interchanges (the reference's rook
                // search moves it too): `D` is zero where `pivot` says
                try testing.expect(chol.isZero(T, p.a[(piv - 1) * (p.lda + 1)]));
                try p.expectFactor(rook, al, "singular");
            };
            // all zero: the first column met, the last for `u`
            for ([_]bool{ true, false }) |upper| for ([_]usize{ 3, 70 }) |n| {
                var p = try Indefinite(T, herm).init(al, n, upper, false, 1);
                defer p.deinit(al);
                @memset(p.w, check.wide(T, scalar(T, 0)));
                p.store();
                var piv: usize = 0;
                try testing.expectError(error.Singular, f(al, check.io, T, uploOf(upper), @intCast(n), p.a.ptr, @intCast(p.lda), p.ipiv.ptr, &piv));
                try testing.expectEqual(if (upper) n else 1, piv);
            };
        }
    };
}

test "relapack sytrf, hetrf: the first bad argument; nothing read for n = 0" {
    const al = testing.allocator;
    const C = types.Complex64;
    var a: [4]C = @splat(chol.poisoned(C));
    var ipiv: [2]isize = @splat(-7);
    try testing.expectError(error.InvalidUplo, sytrf(al, check.io, C, .v, 2, &a, 2, &ipiv, null));
    try testing.expectError(error.InvalidN, hetrf(al, check.io, C, .u, -1, &a, 2, &ipiv, null));
    try testing.expectError(error.InvalidLda, hetrfRook(al, check.io, C, .l, 2, &a, 1, &ipiv, null));
    try sytrfRook(al, check.io, C, .l, 0, &a, 1, &ipiv, null);
    for (a) |x| try testing.expect(isNanT(C, x));
    try testing.expectEqual(-7, ipiv[0]);
}

const BigSym = struct {
    a: []const f64,

    fn run(self: BigSym, io: Io, out: []f64) anyerror!void {
        @memcpy(out, self.a);
        var ipiv: [big_n]isize = undefined;
        try sytrf(testing.allocator, io, f64, .l, big_n, out.ptr, big_n, &ipiv, null);
        @memcpy(out, self.a);
        try sytrfRook(testing.allocator, io, f64, .u, big_n, out.ptr, big_n, &ipiv, null);
    }
};

test "relapack sytrf: the same bits inline and threaded" {
    const al = testing.allocator;
    const a = try al.alloc(f64, big_n * big_n);
    defer al.free(a);
    check.fill(f64, a, 4);
    try expectIoBits(f64, BigSym{ .a = a }, BigSym.run, big_n * big_n, "sytrf");
}

// ---- sygst tests ------------------------------------------------------------------

/// `M X Mᴴ` of `n x n` wide matrices, and its magnitude `|M| |X| |M|ᴴ`.
pub fn sandwich(comptime T: type, al: Allocator, n: usize, m: []const W(T), x: []const W(T), out: []W(T), mag: []f128) !void {
    const t = try al.alloc(W(T), n * n);
    defer al.free(t);
    const tm = try al.alloc(f128, n * n);
    defer al.free(tm);
    const zero = check.wide(T, scalar(T, 0));
    for (0..n) |j| for (0..n) |i| {
        var s = zero;
        var g: f128 = 0;
        for (0..n) |k| {
            s = check.addW(s, check.mulW(m[i + k * n], x[k + j * n]));
            g += check.abs1(m[i + k * n]) * check.abs1(x[k + j * n]);
        }
        t[i + j * n] = s;
        tm[i + j * n] = g;
    };
    for (0..n) |j| for (0..n) |i| {
        var s = zero;
        var g: f128 = 0;
        for (0..n) |k| {
            const c = check.conjW(m[j + k * n]);
            s = check.addW(s, check.mulW(t[i + k * n], c));
            g += tm[i + k * n] * check.abs1(c);
        }
        out[i + j * n] = s;
        mag[i + j * n] = g;
    };
}

fn gstCase(comptime T: type, al: Allocator, itype: isize, upper: bool, n: usize, seed: u64) !void {
    // `B` from a known Cholesky factor, `A` Hermitian, both in `uplo`'s triangle
    var p = try chol.Problem(T).init(al, n, upper, seed);
    defer p.deinit(al);
    var h = try Indefinite(T, isComplex(T)).init(al, n, upper, false, seed + 1);
    defer h.deinit(al);
    const ldb = n + 1;
    const b = try al.alloc(T, @max(1, ldb * n));
    defer al.free(b);
    @memset(b, chol.poisoned(T));
    // `L`, or `U = Lᴴ`, whole in `W(T)`
    const f = try al.alloc(W(T), n * n);
    defer al.free(f);
    for (0..n) |j| for (0..n) |i| {
        const l = check.wide(T, p.l[i + j * n]);
        f[i + j * n] = if (upper) check.conjW(check.wide(T, p.l[j + i * n])) else l;
        if (if (upper) i <= j else i >= j) b[i + j * ldb] = chol.toT(T, f[i + j * n]);
    };
    const routine = if (comptime isComplex(T)) hegst else sygst;
    try routine(al, check.io, T, itype, uploOf(upper), @intCast(n), h.a.ptr, @intCast(h.lda), b.ptr, @intCast(ldb));
    // the result, whole
    const c = try al.alloc(W(T), n * n);
    defer al.free(c);
    for (0..n) |j| for (0..n) |i| {
        const s = if (upper) i <= j else i >= j;
        const v = check.wide(T, if (s) h.a[i + j * h.lda] else h.a[j + i * h.lda]);
        c[i + j * n] = if (s) v else check.conjW(v);
    };
    // `M X Mᴴ` against the other side: `L C Lᴴ` or `Uᴴ C U` against `A`
    // (itype 1); `C` against `Lᴴ A L` or `U A Uᴴ`
    const m = try al.alloc(W(T), n * n);
    defer al.free(m);
    for (0..n) |j| for (0..n) |i| {
        const flip = (itype == 1) == upper;
        m[i + j * n] = if (flip) check.conjW(f[j + i * n]) else f[i + j * n];
    };
    const r = try al.alloc(W(T), n * n);
    defer al.free(r);
    const mag = try al.alloc(f128, n * n);
    defer al.free(mag);
    try sandwich(T, al, n, m, if (itype == 1) c else h.w, r, mag);
    const other = if (itype == 1) h.w else c;
    const eps: f128 = std.math.floatEps(Float(T));
    const ulps: f128 = @floatFromInt(32 * (n + 2));
    for (0..n) |j| for (0..n) |i| {
        const d = check.abs1(check.subW(r[i + j * n], other[i + j * n]));
        const allowed = ulps * eps * (mag[i + j * n] + check.abs1(other[i + j * n])) + std.math.floatMin(Float(T));
        if (!(d <= allowed)) {
            check.say("sygst {s} itype {d} {s} n={d}: ({d},{d}) off by {e}, allowed {e}\n", .{ @typeName(T), itype, if (upper) "u" else "l", n, i, j, @as(f64, @floatCast(d)), @as(f64, @floatCast(allowed)) });
            return error.TestUnexpectedResult;
        }
    };
    for (0..n) |j| for (0..h.lda) |i| {
        const x = h.a[i + j * h.lda];
        const inside = i < n and (if (upper) i <= j else i >= j);
        if (!inside and !isNanT(T, x)) {
            check.say("sygst {s} n={d}: wrote ({d},{d}), outside its triangle\n", .{ @typeName(T), n, i, j });
            return error.TestUnexpectedResult;
        }
        if (comptime isComplex(T)) if (i == j and i < n and x.im != 0) {
            check.say("hegst {s} n={d}: ({d},{d}) has imaginary part {e}\n", .{ @typeName(T), n, i, i, x.im });
            return error.TestUnexpectedResult;
        };
    };
}

test "relapack sygst, hegst: every itype, both triangles, against the products, every type" {
    const al = testing.allocator;
    inline for (test_types) |T| for ([_]isize{ 1, 2, 3 }) |itype| for ([_]bool{ true, false }) |upper| {
        for ([_]usize{ 1, 2, 24, 25, 40, 57 }) |n| {
            if (check.rig() and n > 40) continue;
            gstCase(T, al, itype, upper, n, n * 5 + @as(u64, @intCast(itype))) catch |e| {
                check.say("sygst {s} itype {d} n={d}\n", .{ @typeName(T), itype, n });
                return e;
            };
        }
    };
}

test "relapack sygst, hegst: past the crossover twice over, f64 and c64" {
    if (check.rig()) return error.SkipZigTest;
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64 }) |T| for ([_]isize{ 1, 2 }) |itype| for ([_]bool{ true, false }) |upper| {
        try gstCase(T, al, itype, upper, 100, 7);
    };
}

test "relapack sygst, hegst: the first bad argument, nothing read for n = 0" {
    const al = testing.allocator;
    var a: [4]f64 = @splat(chol.poisoned(f64));
    const b: [4]f64 = @splat(chol.poisoned(f64));
    try testing.expectError(error.InvalidItype, sygst(al, check.io, f64, 0, .v, 2, &a, 2, &b, 2));
    try testing.expectError(error.InvalidItype, sygst(al, check.io, f64, 4, .u, 2, &a, 2, &b, 2));
    try testing.expectError(error.InvalidUplo, sygst(al, check.io, f64, 1, .m, 2, &a, 2, &b, 2));
    try testing.expectError(error.InvalidN, sygst(al, check.io, f64, 1, .u, -1, &a, 2, &b, 2));
    try testing.expectError(error.InvalidLda, sygst(al, check.io, f64, 2, .u, 2, &a, 1, &b, 2));
    try testing.expectError(error.InvalidLdb, sygst(al, check.io, f64, 3, .l, 2, &a, 2, &b, 1));
    try sygst(al, check.io, f64, 1, .l, 0, &a, 1, &b, 1);
    for (a) |x| try testing.expect(std.math.isNan(x));
    var z: [1]types.Complex32 = @splat(chol.poisoned(types.Complex32));
    try testing.expectError(error.InvalidLdb, hegst(al, check.io, types.Complex32, 1, .u, 1, &z, 1, &z, 0));
}

const BigGst = struct {
    a: []const f64,
    b: []const f64,

    fn run(self: BigGst, io: Io, out: []f64) anyerror!void {
        @memcpy(out, self.a);
        try sygst(testing.allocator, io, f64, 1, .l, big_n, out.ptr, big_n, self.b.ptr, big_n);
    }
};

test "relapack sygst: the same bits inline and threaded" {
    const al = testing.allocator;
    const a = try dominant(al, big_n, 5);
    defer al.free(a);
    const b = try dominant(al, big_n, 6);
    defer al.free(b);
    try potrf(al, check.io, f64, .l, big_n, b.ptr, big_n, null);
    try expectIoBits(f64, BigGst{ .a = a, .b = b }, BigGst.run, big_n * big_n, "sygst");
}

// ---- pbtrf tests ------------------------------------------------------------------

/// A Hermitian (real: symmetric) positive-definite band matrix `A = L Lᴴ`,
/// `L` lower with `kd` subdiagonals, its diagonal in [2, 4) and the rest
/// within `1 / (kd + 1)`; stored in `upper`'s band form at `ldab = kd + 2`,
/// every unused place NaN.
fn Band(comptime T: type) type {
    return struct {
        n: usize,
        kd: usize,
        ldab: usize,
        upper: bool,
        ab: []T,
        aw: []W(T),

        const Self = @This();

        fn init(al: Allocator, n: usize, kd: usize, upper: bool, seed: u64) !Self {
            const v = try al.alloc(T, n * n);
            defer al.free(v);
            check.fill(T, v, seed);
            const scale: Float(T) = @floatCast(1 / @as(f64, @floatFromInt(kd + 1)));
            const l = try al.alloc(W(T), n * n);
            defer al.free(l);
            for (0..n) |j| for (0..n) |i| {
                const x = v[i + j * n];
                l[i + j * n] = check.wide(T, if (i < j or i > j + kd) scalar(T, 0) else if (i == j) scalar(T, 3 + chol.re(T, x) / 2) else chol.mulReal(T, x, scale));
            };
            const aw = try al.alloc(W(T), n * n);
            errdefer al.free(aw);
            for (0..n) |j| for (0..n) |i| {
                var s = check.wide(T, scalar(T, 0));
                for (0..@min(i, j) + 1) |k| s = check.addW(s, check.mulW(l[i + k * n], check.conjW(l[j + k * n])));
                aw[i + j * n] = s;
            };
            const ldab = kd + 2;
            const ab = try al.alloc(T, @max(1, ldab * n));
            var self: Self = .{ .n = n, .kd = kd, .ldab = ldab, .upper = upper, .ab = ab, .aw = aw };
            self.store();
            return self;
        }

        fn deinit(self: *Self, al: Allocator) void {
            al.free(self.ab);
            al.free(self.aw);
        }

        /// The storage place of `(i, j)`, `i >= j` (`l`) or `i <= j` (`u`),
        /// inside the band, else null.
        fn place(self: *const Self, i: usize, j: usize) ?usize {
            if (self.upper) {
                if (i > j or j - i > self.kd) return null;
                return self.kd + i - j + j * self.ldab;
            }
            if (i < j or i - j > self.kd) return null;
            return i - j + j * self.ldab;
        }

        fn store(self: *Self) void {
            @memset(self.ab, chol.poisoned(T));
            for (0..self.n) |j| for (0..self.n) |i| if (self.place(i, j)) |at| {
                self.ab[at] = chol.toT(T, self.aw[i + j * self.n]);
            };
        }

        /// The factor's `(i, j)` as lower, `i >= j`, zero outside the band.
        fn got(self: *const Self, i: usize, j: usize) W(T) {
            const at = (if (self.upper) self.place(j, i) else self.place(i, j)) orelse return check.wide(T, scalar(T, 0));
            const v = check.wide(T, self.ab[at]);
            return if (self.upper) check.conjW(v) else v;
        }

        /// `F Fᴴ` against `A` within `2 (n + 2)` ulps of `|F| |F|ᴴ`, the
        /// unused places untouched.
        fn expectFactor(self: *const Self) !void {
            const n = self.n;
            const eps: f128 = std.math.floatEps(Float(T));
            const ulps: f128 = @floatFromInt(2 * (n + 2));
            for (0..n) |j| for (j..n) |i| {
                var s = check.wide(T, scalar(T, 0));
                var mag: f128 = 0;
                const lo = @min(if (i > self.kd) i - self.kd else 0, j + 1);
                for (lo..j + 1) |k| {
                    const x = self.got(i, k);
                    const y = check.conjW(self.got(j, k));
                    s = check.addW(s, check.mulW(x, y));
                    mag += check.abs1(x) * check.abs1(y);
                }
                const d = check.abs1(check.subW(s, self.aw[i + j * n]));
                if (!(d <= ulps * eps * mag + std.math.floatMin(Float(T)))) {
                    check.say("pbtrf {s} {s} n={d} kd={d}: (F Fᴴ)({d},{d}) off by {e}\n", .{ @typeName(T), if (self.upper) "u" else "l", n, self.kd, i, j, @as(f64, @floatCast(d)) });
                    return error.TestUnexpectedResult;
                }
            };
            for (0..n) |j| for (0..self.ldab) |r| {
                const used = r <= self.kd and (if (self.upper) r + j >= self.kd else r + j < n);
                if (!used and !isNanT(T, self.ab[r + j * self.ldab])) {
                    check.say("pbtrf {s} n={d} kd={d}: wrote unused ({d},{d})\n", .{ @typeName(T), n, self.kd, r, j });
                    return error.TestUnexpectedResult;
                }
            };
        }
    };
}

const band_cases = [_][2]usize{ .{ 1, 0 }, .{ 5, 0 }, .{ 30, 1 }, .{ 30, 5 }, .{ 26, 25 }, .{ 50, 10 }, .{ 60, 30 }, .{ 40, 60 } };

test "relapack pbtrf: A = F Fᴴ in band storage, both triangles, every type" {
    const al = testing.allocator;
    inline for (test_types) |T| for ([_]bool{ true, false }) |upper| for (band_cases) |c| {
        var p = try Band(T).init(al, c[0], c[1], upper, c[0] * 3 + c[1]);
        defer p.deinit(al);
        pbtrf(al, check.io, T, uploOf(upper), @intCast(c[0]), @intCast(c[1]), p.ab.ptr, @intCast(p.ldab), null) catch |e| {
            check.say("pbtrf {s} n={d} kd={d}\n", .{ @typeName(T), c[0], c[1] });
            return e;
        };
        try p.expectFactor();
    };
}

test "relapack pbtrf: a rest wider than the band past the block (ReLAPACK's dense case), f64 and c64" {
    if (check.rig()) return error.SkipZigTest;
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64 }) |T| for ([_]bool{ true, false }) |upper| for ([_][2]usize{ .{ 100, 50 }, .{ 100, 60 }, .{ 120, 30 } }) |c| {
        var p = try Band(T).init(al, c[0], c[1], upper, 11);
        defer p.deinit(al);
        try pbtrf(al, check.io, T, uploOf(upper), @intCast(c[0]), @intCast(c[1]), p.ab.ptr, @intCast(p.ldab), null);
        try p.expectFactor();
    };
}

test "relapack pbtrf: the first minor that is not positive definite, NaN too" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex32 }) |T| for ([_]bool{ true, false }) |upper| {
        var p = try Band(T).init(al, 60, 10, upper, 2);
        defer p.deinit(al);
        for ([_]usize{ 0, 23, 30, 59 }) |k| for ([_]f64{ -100, std.math.nan(f64) }) |x| {
            p.store();
            p.ab[p.place(k, k).?] = chol.fromF64(T, x);
            var minor: usize = 0;
            try testing.expectError(error.NotPositiveDefinite, pbtrf(al, check.io, T, uploOf(upper), 60, 10, p.ab.ptr, @intCast(p.ldab), &minor));
            try testing.expectEqual(k + 1, minor);
        };
    };
}

test "relapack pbtrf: the first bad argument, nothing read for n = 0" {
    const al = testing.allocator;
    var ab: [6]f64 = @splat(chol.poisoned(f64));
    try testing.expectError(error.InvalidUplo, pbtrf(al, check.io, f64, .v, 2, 1, &ab, 2, null));
    try testing.expectError(error.InvalidN, pbtrf(al, check.io, f64, .u, -1, 1, &ab, 2, null));
    try testing.expectError(error.InvalidKd, pbtrf(al, check.io, f64, .u, 2, -1, &ab, 2, null));
    try testing.expectError(error.InvalidLdab, pbtrf(al, check.io, f64, .l, 2, 1, &ab, 1, null));
    try pbtrf(al, check.io, f64, .l, 0, 1, &ab, 2, null);
    for (ab) |x| try testing.expect(std.math.isNan(x));
}

// ---- gbtrf tests ------------------------------------------------------------------

/// `gbtrf` on a random `m x n` band matrix (`kl`, `ku`) at `ldab = 2 kl + ku
/// + 2`, the unused places NaN: `P1 L1 P2 L2 ... U` rebuilt in `W(T)` within
/// `2 (min(m, n) + 2)` ulps of the same product in magnitude, the multipliers
/// at most 1 in `|re| + |im|`. `zero_col` zeroes that column first (singular
/// there).
fn gbCase(comptime T: type, al: Allocator, m: usize, n: usize, kl: usize, ku: usize, seed: u64, zero_col: ?usize) !void {
    const kv = kl + ku;
    const ldab = 2 * kl + ku + 2;
    const ab = try al.alloc(T, ldab * n);
    defer al.free(ab);
    @memset(ab, chol.poisoned(T));
    const aw = try al.alloc(W(T), m * n);
    defer al.free(aw);
    const v = try al.alloc(T, m * n);
    defer al.free(v);
    check.fill(T, v, seed);
    for (0..n) |j| for (0..m) |i| {
        const inside = i <= j + kl and j <= i + ku and zero_col != j;
        const x = if (inside) v[i + j * m] else scalar(T, 0);
        aw[i + j * m] = check.wide(T, x);
        if (i + ku >= j and i <= j + kl) ab[kv + i - j + j * ldab] = x;
    };
    // the fill-in rows, which `gbtrf` zeroes where it uses them
    for (0..n) |j| for (0..kl) |r| {
        ab[r + j * ldab] = chol.poisoned(T);
    };
    const ipiv = try al.alloc(isize, @max(1, @min(m, n)));
    defer al.free(ipiv);
    var piv: usize = 0;
    gbtrf(al, check.io, T, @intCast(m), @intCast(n), @intCast(kl), @intCast(ku), ab.ptr, @intCast(ldab), ipiv.ptr, &piv) catch |e| {
        if (e != error.Singular or zero_col == null) return e;
    };
    const mn = @min(m, n);
    if (zero_col) |z| try testing.expectEqual(z + 1, piv);
    // X := U, then X := P_j L_j X from the last step back
    const x = try al.alloc(W(T), m * n);
    defer al.free(x);
    const g = try al.alloc(f128, m * n);
    defer al.free(g);
    for (0..n) |j| for (0..m) |i| {
        const u = if (i <= j and j <= i + kv and i < mn) check.wide(T, ab[kv + i - j + j * ldab]) else check.wide(T, scalar(T, 0));
        x[i + j * m] = u;
        g[i + j * m] = check.abs1(u);
    };
    var s = mn;
    while (s > 0) {
        s -= 1;
        const km = @min(kl, m - s - 1);
        for (1..km + 1) |r| {
            const l = check.wide(T, ab[kv + r + s * ldab]);
            // the pivot largest in `|re| + |im|`: `|l|` at most `sqrt 2` complex
            const l2 = if (comptime isComplex(T)) l.re * l.re + l.im * l.im else l * l;
            if (!(l2 <= (if (comptime isComplex(T)) 2 else 1) * (1 + 8 * std.math.floatEps(Float(T))))) {
                check.say("gbtrf {s}: |L({d},{d})|^2 = {e}\n", .{ @typeName(T), s + r, s, @as(f64, @floatCast(l2)) });
                return error.TestUnexpectedResult;
            }
            for (0..n) |j| {
                x[s + r + j * m] = check.addW(x[s + r + j * m], check.mulW(l, x[s + j * m]));
                g[s + r + j * m] += check.abs1(l) * g[s + j * m];
            }
        }
        const p: usize = @intCast(ipiv[s] - 1);
        if (p != s) for (0..n) |j| {
            std.mem.swap(W(T), &x[s + j * m], &x[p + j * m]);
            std.mem.swap(f128, &g[s + j * m], &g[p + j * m]);
        };
    }
    const eps: f128 = std.math.floatEps(Float(T));
    const ulps: f128 = @floatFromInt(2 * (mn + 2));
    for (0..n) |j| for (0..m) |i| {
        const d = check.abs1(check.subW(x[i + j * m], aw[i + j * m]));
        if (!(d <= ulps * eps * g[i + j * m] + std.math.floatMin(Float(T)))) {
            check.say("gbtrf {s} {d}x{d} kl={d} ku={d}: (P L U)({d},{d}) off by {e}\n", .{ @typeName(T), m, n, kl, ku, i, j, @as(f64, @floatCast(d)) });
            return error.TestUnexpectedResult;
        }
    };
    // the padding row, never used
    for (0..n) |j| if (!isNanT(T, ab[ldab - 1 + j * ldab])) {
        check.say("gbtrf {s}: wrote the padding of column {d}\n", .{ @typeName(T), j });
        return error.TestUnexpectedResult;
    };
}

const gb_cases = [_][4]usize{ .{ 1, 1, 0, 0 }, .{ 5, 5, 0, 0 }, .{ 30, 30, 2, 3 }, .{ 60, 60, 10, 5 }, .{ 50, 40, 45, 3 }, .{ 40, 50, 3, 45 }, .{ 100, 100, 30, 20 }, .{ 80, 100, 90, 7 }, .{ 3, 7, 1, 1 } };

test "relapack gbtrf: P L U = A in band storage, tall, wide and square, every type" {
    const al = testing.allocator;
    inline for (test_types) |T| for (gb_cases) |c| {
        if (check.rig() and c[0] * c[1] > 4000) continue;
        gbCase(T, al, c[0], c[1], c[2], c[3], c[0] + c[2] * 7, null) catch |e| {
            check.say("gbtrf {s} {d}x{d} kl={d} ku={d}\n", .{ @typeName(T), c[0], c[1], c[2], c[3] });
            return e;
        };
    };
}

test "relapack gbtrf: a zero column is singular there, the factorization completed" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64 }) |T| for ([_]usize{ 0, 7, 39 }) |z| {
        try gbCase(T, al, 40, 40, 6, 4, 3, z);
    };
}

test "relapack gbtrf: the first bad argument, nothing read when empty" {
    const al = testing.allocator;
    var ab: [8]f64 = @splat(chol.poisoned(f64));
    var ipiv: [2]isize = @splat(-1);
    try testing.expectError(error.InvalidM, gbtrf(al, check.io, f64, -1, 2, 0, 0, &ab, 1, &ipiv, null));
    try testing.expectError(error.InvalidN, gbtrf(al, check.io, f64, 2, -1, 0, 0, &ab, 1, &ipiv, null));
    try testing.expectError(error.InvalidKl, gbtrf(al, check.io, f64, 2, 2, -1, 0, &ab, 1, &ipiv, null));
    try testing.expectError(error.InvalidKu, gbtrf(al, check.io, f64, 2, 2, 0, -1, &ab, 1, &ipiv, null));
    try testing.expectError(error.InvalidLdab, gbtrf(al, check.io, f64, 2, 2, 1, 1, &ab, 3, &ipiv, null));
    try gbtrf(al, check.io, f64, 0, 2, 1, 1, &ab, 4, &ipiv, null);
    for (ab) |x| try testing.expect(std.math.isNan(x));
    try testing.expectEqual(-1, ipiv[0]);
}

// ---- trsyl tests ------------------------------------------------------------------

/// An `n x n` upper triangular (real: quasi-triangular, a 2 x 2 block every
/// `every` rows) matrix whose eigenvalues' real parts lie in `[lo, lo + 2)`,
/// the rest within `1 / n`; zero below.
fn schur(comptime T: type, al: Allocator, n: usize, lo: f64, every: usize, seed: u64) ![]T {
    const a = try al.alloc(T, n * n);
    check.fill(T, a, seed);
    const s: Float(T) = @floatCast(1 / @as(f64, @floatFromInt(@max(n, 1))));
    for (0..n) |j| for (0..n) |i| {
        const v = &a[i + j * n];
        if (i > j) {
            v.* = scalar(T, 0);
        } else if (i == j) {
            v.* = chol.add(T, scalar(T, @floatCast(lo + 1)), chol.mulReal(T, v.*, 0.5));
        } else v.* = chol.mulReal(T, v.*, s);
    };
    if (comptime !isComplex(T)) {
        var i: usize = 0;
        while (i + 1 < n) : (i += every) {
            // `[d, 0.7; -x, d]`, eigenvalues `d +- i sqrt(0.7 x)`
            a[i + 1 + i * n] = -0.5;
            a[i + (i + 1) * n] = 0.7;
            a[i + 1 + (i + 1) * n] = a[i + i * n];
        }
    }
    return a;
}

fn trsylCase(comptime T: type, al: Allocator, ta: bool, tb: bool, isgn: isize, m: usize, n: usize, seed: u64) !void {
    const a = try schur(T, al, m, 2, 3, seed);
    defer al.free(a);
    // `-isgn B`'s spectrum away from `A`'s
    const b = try schur(T, al, n, if (isgn == 1) 2 else -4, 4, seed + 1);
    defer al.free(b);
    const c0 = try al.alloc(T, m * n);
    defer al.free(c0);
    check.fill(T, c0, seed + 2);
    const x = try al.dupe(T, c0);
    defer al.free(x);
    const op = struct {
        fn of(comptime U: type, t: bool) Trans {
            return if (!t) .n else if (comptime isComplex(U)) .c else .t;
        }
    };
    var scale: Float(T) = 0;
    try trsyl(al, check.io, T, op.of(T, ta), op.of(T, tb), isgn, @intCast(m), @intCast(n), a.ptr, @intCast(m), b.ptr, @intCast(n), x.ptr, @intCast(m), &scale);
    try testing.expect(scale > 0 and scale <= 1);
    const eps: f128 = std.math.floatEps(Float(T));
    const ulps: f128 = @floatFromInt(32 * (m + n + 2));
    const sg: f128 = @floatFromInt(isgn);
    for (0..n) |j| for (0..m) |i| {
        var r = check.scaleW(-@as(f128, scale), check.wide(T, c0[i + j * m]));
        var mag: f128 = scale * check.abs1(check.wide(T, c0[i + j * m]));
        for (0..m) |k| {
            const av = check.wide(T, if (ta) a[k + i * m] else a[i + k * m]);
            const aw = if (ta) check.conjW(av) else av;
            const xv = check.wide(T, x[k + j * m]);
            r = check.addW(r, check.mulW(aw, xv));
            mag += check.abs1(aw) * check.abs1(xv);
        }
        for (0..n) |k| {
            const bv = check.wide(T, if (tb) b[j + k * n] else b[k + j * n]);
            const bw = if (tb) check.conjW(bv) else bv;
            const xv = check.wide(T, x[i + k * m]);
            r = check.addW(r, check.scaleW(sg, check.mulW(xv, bw)));
            mag += check.abs1(bw) * check.abs1(xv);
        }
        const d = check.abs1(r);
        if (!(d <= ulps * eps * mag + std.math.floatMin(Float(T)))) {
            check.say("trsyl {s} {d}x{d} {any}{any} isgn {d}: residual ({d},{d}) {e}, allowed {e}\n", .{ @typeName(T), m, n, ta, tb, isgn, i, j, @as(f64, @floatCast(d)), @as(f64, @floatCast(ulps * eps * mag)) });
            return error.TestUnexpectedResult;
        }
    };
}

const trsyl_shapes = [_][2]usize{ .{ 1, 1 }, .{ 5, 3 }, .{ 24, 24 }, .{ 30, 10 }, .{ 10, 40 }, .{ 60, 50 }, .{ 100, 70 } };

test "relapack trsyl: op(A) X + isgn X op(B) = scale C, every op and sign, quasi-triangular real, every type" {
    const al = testing.allocator;
    inline for (test_types) |T| for (trsyl_shapes) |s| for ([_]bool{ false, true }) |ta| for ([_]bool{ false, true }) |tb| for ([_]isize{ 1, -1 }) |sg| {
        if (check.rig() and s[0] * s[1] > 1000) continue;
        trsylCase(T, al, ta, tb, sg, s[0], s[1], s[0] * 13 + s[1]) catch |e| {
            check.say("trsyl {s} {d}x{d}\n", .{ @typeName(T), s[0], s[1] });
            return e;
        };
    };
}

test "relapack trsyl: common eigenvalues are perturbed, C and scale written; the first bad argument" {
    const al = testing.allocator;
    // A X - X A = C, A = 2 I: singular
    var a = [_]f64{ 2, 0, 0, 2 };
    var c = [_]f64{ 1, 2, 3, 4 };
    var scale: f64 = 0;
    try testing.expectError(error.Perturbed, trsyl(al, check.io, f64, .n, .n, -1, 2, 2, &a, 2, &a, 2, &c, 2, &scale));
    try testing.expect(scale > 0 and scale <= 1 and std.math.isFinite(c[0]));
    const z: [4]types.Complex64 = @splat(chol.poisoned(types.Complex64));
    var zc = z;
    try testing.expectError(error.InvalidTrana, trsyl(al, check.io, types.Complex64, .t, .n, 1, 2, 2, &z, 2, &z, 2, &zc, 2, &scale));
    try testing.expectError(error.InvalidTranb, trsyl(al, check.io, f64, .n, .r, 1, 2, 2, &a, 2, &a, 2, &c, 2, &scale));
    try testing.expectError(error.InvalidIsgn, trsyl(al, check.io, f64, .n, .n, 0, 2, 2, &a, 2, &a, 2, &c, 2, &scale));
    try testing.expectError(error.InvalidM, trsyl(al, check.io, f64, .n, .n, 1, -1, 2, &a, 2, &a, 2, &c, 2, &scale));
    try testing.expectError(error.InvalidN, trsyl(al, check.io, f64, .n, .n, 1, 2, -1, &a, 2, &a, 2, &c, 2, &scale));
    try testing.expectError(error.InvalidLda, trsyl(al, check.io, f64, .n, .n, 1, 2, 2, &a, 1, &a, 2, &c, 2, &scale));
    try testing.expectError(error.InvalidLdb, trsyl(al, check.io, f64, .n, .n, 1, 2, 2, &a, 2, &a, 1, &c, 2, &scale));
    try testing.expectError(error.InvalidLdc, trsyl(al, check.io, f64, .n, .n, 1, 2, 2, &a, 2, &a, 2, &c, 1, &scale));
    scale = 0;
    try trsyl(al, check.io, types.Complex64, .c, .c, 1, 0, 2, &z, 1, &z, 2, &zc, 1, &scale);
    try testing.expectEqual(1, scale);
    for (zc) |v| try testing.expect(isNanT(types.Complex64, v));
}

const BigSyl = struct {
    a: []const f64,
    b: []const f64,
    c: []const f64,

    fn run(self: BigSyl, io: Io, out: []f64) anyerror!void {
        @memcpy(out, self.c);
        var scale: f64 = 0;
        try trsyl(testing.allocator, io, f64, .n, .t, 1, big_n, big_n, self.a.ptr, big_n, self.b.ptr, big_n, out.ptr, big_n, &scale);
    }
};

test "relapack trsyl: the same bits inline and threaded" {
    const al = testing.allocator;
    const a = try schur(f64, al, big_n, 2, 3, 1);
    defer al.free(a);
    const b = try schur(f64, al, big_n, 2, 5, 2);
    defer al.free(b);
    const c = try al.alloc(f64, big_n * big_n);
    defer al.free(c);
    check.fill(f64, c, 3);
    try expectIoBits(f64, BigSyl{ .a = a, .b = b, .c = c }, BigSyl.run, big_n * big_n, "trsyl");
}

// ---- tgsyl tests ------------------------------------------------------------------

/// An `n x n` upper triangular matrix like `schur`'s, without its 2 x 2
/// blocks: `D` and `E` of a generalized Schur form.
fn upperTri(comptime T: type, al: Allocator, n: usize, lo: f64, seed: u64) ![]T {
    const a = try schur(T, al, n, lo, n + 1, seed);
    if (comptime !isComplex(T)) for (0..n -| 1) |i| {
        a[i + 1 + i * n] = 0;
    };
    return a;
}

/// `tgsyl`'s operands, each at its own padded leading dimension so a mixed
/// up one shows: `(A, D)` with eigenvalues in `[0.8, 8]`, `(B, E)` in
/// `[-8, -0.8]`.
fn Gen(comptime T: type) type {
    return struct {
        m: usize,
        n: usize,
        a: []T,
        b: []T,
        d: []T,
        e: []T,
        c0: []T,
        f0: []T,
        c: []T,
        f: []T,
        const ld = .{ .a = 1, .b = 2, .c = 3, .d = 2, .e = 1, .f = 5 };

        const Self = @This();

        fn init(al: Allocator, m: usize, n: usize, seed: u64) !Self {
            var g: Self = undefined;
            g.m = m;
            g.n = n;
            g.a = try padded(al, m, ld.a, try schur(T, al, m, 2, 3, seed));
            g.d = try padded(al, m, ld.d, try upperTri(T, al, m, 0.5, seed + 1));
            g.b = try padded(al, n, ld.b, try schur(T, al, n, -4, 4, seed + 2));
            g.e = try padded(al, n, ld.e, try upperTri(T, al, n, 0.5, seed + 3));
            g.c0 = try al.alloc(T, (m + ld.c) * n);
            check.fill(T, g.c0, seed + 4);
            g.f0 = try al.alloc(T, (m + ld.f) * n);
            check.fill(T, g.f0, seed + 5);
            g.c = try al.dupe(T, g.c0);
            g.f = try al.dupe(T, g.f0);
            return g;
        }

        /// `x` (`k x k`, packed) at leading dimension `k + pad`, freed.
        fn padded(al: Allocator, k: usize, pad: usize, x: []T) ![]T {
            defer al.free(x);
            const y = try al.alloc(T, (k + pad) * k);
            @memset(y, chol.poisoned(T));
            for (0..k) |j| @memcpy(y[j * (k + pad) ..][0..k], x[j * k ..][0..k]);
            return y;
        }

        fn deinit(g: Self, al: Allocator) void {
            for ([_][]T{ g.a, g.b, g.d, g.e, g.c0, g.f0, g.c, g.f }) |x| al.free(x);
        }

        fn reset(g: Self) void {
            @memcpy(g.c, g.c0);
            @memcpy(g.f, g.f0);
        }

        fn run(g: Self, al: Allocator, io: Io, trans: Trans, ijob: isize, scale: *Float(T), dif: *Float(T)) TgsylError!void {
            const m = g.m;
            const n = g.n;
            return tgsyl(al, io, T, trans, ijob, @intCast(m), @intCast(n), g.a.ptr, @intCast(m + ld.a), g.b.ptr, @intCast(n + ld.b), g.c.ptr, @intCast(m + ld.c), g.d.ptr, @intCast(m + ld.d), g.e.ptr, @intCast(n + ld.e), g.f.ptr, @intCast(m + ld.f), scale, dif);
        }

        fn at(x: []const T, ldx: usize, i: usize, j: usize) Wt {
            return check.wide(T, x[i + j * ldx]);
        }
        const Wt = @TypeOf(check.wide(T, scalar(T, 0)));

        /// Both equations' residuals within `32 (m + n + 2)` ulps of their
        /// terms' magnitudes; padding untouched.
        fn expectSolved(g: Self, tr: bool, scale: Float(T)) !void {
            const m = g.m;
            const n = g.n;
            const lda = m + ld.a;
            const ldb = n + ld.b;
            const ldc = m + ld.c;
            const ldd = m + ld.d;
            const lde = n + ld.e;
            const ldf = m + ld.f;
            const eps: f128 = std.math.floatEps(Float(T));
            const ulps: f128 = @floatFromInt(32 * (m + n + 2));
            const cjw = struct {
                fn f(x: Wt, c: bool) Wt {
                    return if (c) check.conjW(x) else x;
                }
            }.f;
            for (0..n) |j| for (0..m) |i| for (0..2) |eq| {
                const rhs = if (eq == 0) at(g.c0, ldc, i, j) else at(g.f0, ldf, i, j);
                // notrans: A R - L B, D R - L E; trans: Aᴴ R + Dᴴ L, R Bᴴ + L Eᴴ
                var r = check.scaleW(if (tr and eq == 1) @as(f128, scale) else -@as(f128, scale), rhs);
                var mag: f128 = scale * check.abs1(rhs);
                for (0..@max(m, n)) |k| {
                    var terms: [2][2]Wt = undefined;
                    var used: [2]bool = .{ false, false };
                    if (!tr) {
                        if (k < m) {
                            terms[0] = .{ at(if (eq == 0) g.a else g.d, if (eq == 0) lda else ldd, i, k), at(g.c, ldc, k, j) };
                            used[0] = true;
                        }
                        if (k < n) {
                            terms[1] = .{ check.scaleW(-1, at(g.f, ldf, i, k)), at(if (eq == 0) g.b else g.e, if (eq == 0) ldb else lde, k, j) };
                            used[1] = true;
                        }
                    } else if (eq == 0) {
                        if (k < m) {
                            terms[0] = .{ cjw(at(g.a, lda, k, i), true), at(g.c, ldc, k, j) };
                            terms[1] = .{ cjw(at(g.d, ldd, k, i), true), at(g.f, ldf, k, j) };
                            used = .{ true, true };
                        }
                    } else if (k < n) {
                        terms[0] = .{ at(g.c, ldc, i, k), cjw(at(g.b, ldb, j, k), true) };
                        terms[1] = .{ at(g.f, ldf, i, k), cjw(at(g.e, lde, j, k), true) };
                        used = .{ true, true };
                    }
                    for (terms, used) |t, u| if (u) {
                        r = check.addW(r, check.mulW(t[0], t[1]));
                        mag += check.abs1(t[0]) * check.abs1(t[1]);
                    };
                }
                const dv = check.abs1(r);
                if (!(dv <= ulps * eps * mag + std.math.floatMin(Float(T)))) {
                    check.say("tgsyl {s} {d}x{d} trans {any}: equation {d} ({d},{d}) residual {e}, allowed {e}\n", .{ @typeName(T), m, n, tr, eq + 1, i, j, @as(f64, @floatCast(dv)), @as(f64, @floatCast(ulps * eps * mag)) });
                    return error.TestUnexpectedResult;
                }
            };
        }
    };
}

fn transOf(comptime T: type, tr: bool) Trans {
    return if (!tr) .n else if (comptime isComplex(T)) .c else .t;
}

const tgsyl_shapes = [_][2]usize{ .{ 1, 1 }, .{ 2, 3 }, .{ 5, 4 }, .{ 24, 24 }, .{ 25, 10 }, .{ 10, 40 }, .{ 60, 50 }, .{ 33, 70 } };

test "relapack tgsyl: A R - L B = scale C, D R - L E = scale F and the transposed system, every type" {
    const al = testing.allocator;
    inline for (test_types) |T| for (tgsyl_shapes) |s| {
        if (check.rig() and s[0] * s[1] > 1000) continue;
        const g = try Gen(T).init(al, s[0], s[1], s[0] * 7 + s[1]);
        defer g.deinit(al);
        for ([_]bool{ false, true }) |tr| {
            g.reset();
            var scale: Float(T) = 0;
            var dif: Float(T) = -1;
            try g.run(al, check.io, transOf(T, tr), 0, &scale, &dif);
            try testing.expect(scale > 0 and scale <= 1 and dif == -1);
            g.expectSolved(tr, scale) catch |e| {
                check.say("tgsyl {s} {d}x{d}\n", .{ @typeName(T), s[0], s[1] });
                return e;
            };
        }
        // ijob 1 and 2 solve as 0 does, bit for bit, and estimate
        const x = try al.dupe(T, g.c);
        defer al.free(x);
        const y = try al.dupe(T, g.f);
        defer al.free(y);
        g.reset();
        var scale: Float(T) = 0;
        var dif: Float(T) = -1;
        try g.run(al, check.io, .n, 0, &scale, &dif);
        @memcpy(x, g.c);
        @memcpy(y, g.f);
        for ([_]isize{ 1, 2 }) |job| {
            g.reset();
            try g.run(al, check.io, .n, job, &scale, &dif);
            try testing.expect(dif > 0 and std.math.isFinite(dif));
            // by value: f80's padding bytes are not its bits
            for (x, g.c) |u, v| try testing.expect(chol.sameBits(T, u, v));
            for (y, g.f) |u, v| try testing.expect(chol.sameBits(T, u, v));
        }
    };
}

/// `σ_min(Z)` of `dtgsyl.f`'s (2), from `gesvd` in `f64` or `c64`.
fn sigmaMin(comptime T: type, al: Allocator, g: Gen(T)) !f64 {
    const C = if (comptime isComplex(T)) types.Complex64 else f64;
    const m = g.m;
    const n = g.n;
    const k = m * n;
    const z = try al.alloc(C, 4 * k * k);
    defer al.free(z);
    @memset(z, scalar(C, 0));
    const lda = m + Gen(T).ld.a;
    const ldb = n + Gen(T).ld.b;
    const ldd = m + Gen(T).ld.d;
    const lde = n + Gen(T).ld.e;
    const to = struct {
        fn f(v: T) C {
            return if (comptime isComplex(T)) .init(@floatCast(v.re), @floatCast(v.im)) else @floatCast(v);
        }
    }.f;
    for (0..n) |j| for (0..m) |i| {
        const r = i + m * j;
        for (0..m) |kk| {
            z[r + 2 * k * (kk + m * j)] = to(g.a[i + kk * lda]);
            z[k + r + 2 * k * (kk + m * j)] = to(g.d[i + kk * ldd]);
        }
        for (0..n) |kk| {
            z[r + 2 * k * (k + i + m * kk)] = chol.neg(C, to(g.b[kk + j * ldb]));
            z[k + r + 2 * k * (k + i + m * kk)] = chol.neg(C, to(g.e[kk + j * lde]));
        }
    };
    const s = try al.alloc(f64, 2 * k);
    defer al.free(s);
    var dummy: [1]C = undefined;
    try faer_svd.gesvd(al, check.io, C, .n, .n, @intCast(2 * k), @intCast(2 * k), z.ptr, @intCast(2 * k), s.ptr, &dummy, 1, &dummy, 1, null);
    return s[2 * k - 1];
}

test "relapack tgsyl: dif estimates σ_min of the Kronecker form; ijob 1 and 3, 2 and 4 agree" {
    const al = testing.allocator;
    inline for (test_types) |T| for ([_][2]usize{ .{ 1, 1 }, .{ 2, 2 }, .{ 3, 5 }, .{ 6, 4 } }) |s| {
        const g = try Gen(T).init(al, s[0], s[1], s[0] * 11 + s[1]);
        defer g.deinit(al);
        const sig = try sigmaMin(T, al, g);
        var difs: [5]Float(T) = undefined;
        for (1..5) |job| {
            g.reset();
            var scale: Float(T) = 0;
            difs[job] = -1;
            try g.run(al, check.io, .n, @intCast(job), &scale, &difs[job]);
            // an upper bound on σ_min, within sqrt(2 m n) of it as the
            // estimate's right-hand sides are about that long
            const d: f64 = @floatCast(difs[job]);
            const lim = 4 * @sqrt(2 * @as(f64, @floatFromInt(s[0] * s[1])));
            if (!(d >= sig * (1 - 1e-3) / lim and d <= sig * lim)) {
                check.say("tgsyl {s} {d}x{d} ijob {d}: dif {e}, sigma_min {e}\n", .{ @typeName(T), s[0], s[1], job, d, sig });
                return error.TestUnexpectedResult;
            }
        }
        try testing.expectEqual(difs[1], difs[3]);
        try testing.expectEqual(difs[2], difs[4]);
    };
}

test "relapack tgsyl: the recursion's dif is the unblocked xTGSY2's, every base case's blocks counted" {
    // ReLAPACK keeps only the last base case's `PQ`: its `ijob` 2 and 4
    // estimate would come out about sqrt(base cases) too small
    const al = testing.allocator;
    inline for (.{ f64, types.Complex32 }) |T| {
        const F = Float(T);
        const g = try Gen(T).init(al, 40, 30, 3);
        defer g.deinit(al);
        for ([_]isize{ 3, 4 }) |job| {
            g.reset();
            var scale: F = 0;
            var dif: F = 0;
            try g.run(al, check.io, .n, job, &scale, &dif);
            // the whole problem as one `xTGSY2`
            g.reset();
            for (0..g.n) |j| for (0..g.m) |i| {
                g.c[i + j * (g.m + 3)] = scalar(T, 0);
                g.f[i + j * (g.m + 5)] = scalar(T, 0);
            };
            var acc: sylvester.Dif(F) = .{};
            _ = try sylvester.tgsy2(al, check.io, T, false, @intCast(job - 2), g.m, g.n, g.a.ptr, g.m + 1, g.b.ptr, g.n + 2, g.c.ptr, g.m + 3, g.d.ptr, g.m + 2, g.e.ptr, g.n + 1, g.f.ptr, g.m + 5, &acc);
            const num: F = @floatFromInt(if (job == 3) 2 * g.m * g.n else acc.pq);
            const want = @sqrt(num) / (acc.scale * @sqrt(acc.sum));
            try testing.expect(@abs(dif - want) <= 1e3 * std.math.floatEps(F) * want);
        }
    }
}

test "relapack tgsyl: common eigenvalues are perturbed, written and finite; scale below 1; the first bad argument" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64 }) |T| {
        const F = Float(T);
        // (A, D) = (B, E): Z is singular; a huge right-hand side scales
        const a = [_]T{ scalar(T, 2), scalar(T, 0), scalar(T, 1), scalar(T, 2) };
        const d = [_]T{ scalar(T, 1), scalar(T, 0), scalar(T, 0.5), scalar(T, 1) };
        for ([_]F{ 1, std.math.floatMax(F) / 4 }) |big| {
            var c: [4]T = @splat(scalar(T, big));
            var f: [4]T = @splat(scalar(T, big));
            var scale: F = 0;
            var dif: F = -1;
            try testing.expectError(error.Perturbed, tgsyl(al, check.io, T, .n, 0, 2, 2, &a, 2, &a, 2, &c, 2, &d, 2, &d, 2, &f, 2, &scale, &dif));
            try testing.expect(scale > 0 and scale <= 1);
            if (big > 1) try testing.expect(scale < 1);
            for (c ++ f) |x| try testing.expect(if (comptime isComplex(T)) std.math.isFinite(x.re) and std.math.isFinite(x.im) else std.math.isFinite(x));
        }
        var c: [4]T = @splat(chol.poisoned(T));
        var scale: F = 0;
        var dif: F = 0;
        const t: Trans = if (comptime isComplex(T)) .t else .c;
        const z = &c;
        try testing.expectError(error.InvalidTrans, tgsyl(al, check.io, T, t, 0, 2, 2, z, 2, z, 2, z, 2, z, 2, z, 2, z, 2, &scale, &dif));
        try testing.expectError(error.InvalidIjob, tgsyl(al, check.io, T, .n, 5, 2, 2, z, 2, z, 2, z, 2, z, 2, z, 2, z, 2, &scale, &dif));
        try testing.expectError(error.InvalidM, tgsyl(al, check.io, T, transOf(T, true), 9, 0, 2, z, 2, z, 2, z, 2, z, 2, z, 2, z, 2, &scale, &dif));
        try testing.expectError(error.InvalidN, tgsyl(al, check.io, T, .n, 0, 2, 0, z, 2, z, 2, z, 2, z, 2, z, 2, z, 2, &scale, &dif));
        try testing.expectError(error.InvalidLda, tgsyl(al, check.io, T, .n, 0, 2, 2, z, 1, z, 2, z, 2, z, 2, z, 2, z, 2, &scale, &dif));
        try testing.expectError(error.InvalidLdb, tgsyl(al, check.io, T, .n, 0, 2, 2, z, 2, z, 1, z, 2, z, 2, z, 2, z, 2, &scale, &dif));
        try testing.expectError(error.InvalidLdc, tgsyl(al, check.io, T, .n, 0, 2, 2, z, 2, z, 2, z, 1, z, 2, z, 2, z, 2, &scale, &dif));
        try testing.expectError(error.InvalidLdd, tgsyl(al, check.io, T, .n, 0, 2, 2, z, 2, z, 2, z, 2, z, 1, z, 2, z, 2, &scale, &dif));
        try testing.expectError(error.InvalidLde, tgsyl(al, check.io, T, .n, 0, 2, 2, z, 2, z, 2, z, 2, z, 2, z, 1, z, 2, &scale, &dif));
        try testing.expectError(error.InvalidLdf, tgsyl(al, check.io, T, .n, 0, 2, 2, z, 2, z, 2, z, 2, z, 2, z, 2, z, 1, &scale, &dif));
        for (c) |x| try testing.expect(isNanT(T, x));
    }
}

const BigGsyl = struct {
    g: Gen(f64),

    fn run(self: BigGsyl, io: Io, out: []f64) anyerror!void {
        self.g.reset();
        var scale: f64 = 0;
        var dif: f64 = 0;
        try self.g.run(testing.allocator, io, .n, 1, &scale, &dif);
        @memcpy(out[0..self.g.c.len], self.g.c);
        @memcpy(out[self.g.c.len..][0..self.g.f.len], self.g.f);
        out[out.len - 1] = dif;
    }
};

test "relapack tgsyl: the same bits inline and threaded" {
    const al = testing.allocator;
    const g = try Gen(f64).init(al, big_n, big_n, 9);
    defer g.deinit(al);
    try expectIoBits(f64, BigGsyl{ .g = g }, BigGsyl.run, g.c.len + g.f.len + 1, "tgsyl");
}
