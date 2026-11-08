//! ReLAPACK's recursive algorithms (Elmar Peise, MIT), from OpenBLAS
//! v0.3.34-317's `relapack/src/`: `sygst`, `hegst` and `pbtrf`, and under
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

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Uplo = types.Uplo;
const Diag = types.Diag;
const Trans = types.Trans;
const Side = types.Side;
const isComplex = types.isComplex;
const scalar = ops.scalar;

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
/// its errors and `minor`. The base case is `potf2`.
pub fn potrf(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, a: [*]T, lda: isize, minor: ?*usize) chol.PotrfError!void {
    comptime ops.requireFloat(T, "potrf");
    if (!try interface.checkSquare(T, uplo, n, lda)) return;
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
    comptime ops.requireFloat(T, "getrf");
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
    if (n <= crossover) return lu.getf2(T, interface.table(T), m, n, a, lda, ipiv);
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
    comptime ops.requireFloat(T, "trtri");
    if (!try interface.checkSquare(T, uplo, n, lda)) return;
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
    comptime ops.requireFloat(T, "lauum");
    if (!try interface.checkSquare(T, uplo, n, lda)) return;
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
    comptime ops.requireFloat(T, "sygst");
    if (comptime isComplex(T)) @compileError("sygst: " ++ @typeName(T) ++ " is complex; hegst takes it");
    return gstImpl(allocator, io, T, itype, uplo, n, a, lda, b, ldb);
}

/// `?HEGST` (`relapack/src/zhegst.c`, `chegst.c`): `sygst` for the
/// Hermitian-definite problem, `ᴴ` for `ᵀ`; `A`'s diagonal comes out real.
/// Complex only.
pub fn hegst(allocator: Allocator, io: Io, comptime T: type, itype: isize, uplo: Uplo, n: isize, a: [*]T, lda: isize, b: [*]const T, ldb: isize) SygstError!void {
    comptime ops.requireComplex(T, "hegst");
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
        y[i + j * ldy] = ops.add(T, y[i + j * ldy], x[i + j * ldx]);
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
    const mul = ops.mul;
    const add = ops.add;
    for (0..n) |k| {
        const akk0 = ops.re(T, a.at(k, k).*);
        const bkk = ops.re(T, b.at(k, k).*);
        if (inv) {
            // inv(L) A inv(Lᴴ): the column below `k`, then the rest
            const akk = akk0 / (bkk * bkk);
            a.at(k, k).* = scalar(T, akk);
            if (k + 1 == n) continue;
            const r = 1 / bkk;
            for (k + 1..n) |i| a.at(i, k).* = ops.mulReal(T, a.at(i, k).*, r);
            const c = -0.5 * akk;
            for (k + 1..n) |i| a.at(i, k).* = add(T, a.at(i, k).*, ops.mulReal(T, b.at(i, k).*, c));
            rank2(T, a, k + 1, n, k, b, k, scalar(T, -1));
            for (k + 1..n) |i| a.at(i, k).* = add(T, a.at(i, k).*, ops.mulReal(T, b.at(i, k).*, c));
            // x := inv(B22) x, `xTRSV` `N`
            for (k + 1..n) |j| {
                const xj = a.at(j, k);
                if (ops.isZero(T, xj.*)) continue;
                xj.* = ops.div(T, xj.*, b.at(j, j).*);
                const t = xj.*;
                for (j + 1..n) |i| a.at(i, k).* = ops.sub(T, a.at(i, k).*, mul(T, t, b.at(i, j).*));
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
            for (0..k) |j| a.at(k, j).* = cj(add(T, cj(a.at(k, j).*), ops.mulReal(T, cj(b.at(k, j).*), c)));
            rank2Row(T, a, k, b, scalar(T, 1));
            for (0..k) |j| a.at(k, j).* = cj(ops.mulReal(T, add(T, cj(a.at(k, j).*), ops.mulReal(T, cj(b.at(k, j).*), c)), bkk));
            a.at(k, k).* = scalar(T, akk0 * (bkk * bkk));
        }
    }
}

/// `xSYR2` or `xHER2` (lower) with `alpha` on `A(lo:hi, lo:hi)`, `x` column
/// `xc` of `A` and `y` column `yc` of `B`, both from row `lo`.
fn rank2(comptime T: type, a: Tview(T), lo: usize, hi: usize, xc: usize, b: Tview(T), yc: usize, alpha: T) void {
    const herm = comptime isComplex(T);
    for (lo..hi) |j| {
        const xj = a.at(j, xc).*;
        const yj = b.at(j, yc).*;
        if (ops.isZero(T, xj) and ops.isZero(T, yj)) {
            if (comptime herm) a.at(j, j).* = scalar(T, a.at(j, j).re);
            continue;
        }
        const t1 = ops.mul(T, alpha, ops.conj(T, yj));
        const t2 = ops.conj(T, ops.mul(T, alpha, xj));
        if (comptime herm) {
            const d = ops.add(T, ops.mul(T, xj, t1), ops.mul(T, yj, t2));
            a.at(j, j).* = scalar(T, a.at(j, j).re + d.re);
        }
        for (if (herm) j + 1 else j..hi) |i| {
            a.at(i, j).* = ops.add(T, a.at(i, j).*, ops.add(T, ops.mul(T, a.at(i, xc).*, t1), ops.mul(T, b.at(i, yc).*, t2)));
        }
    }
}

/// `rank2` with `x` and `y` the conjugated rows `k` of `A` and `B` left of
/// the diagonal, on `A(:k, :k)`.
fn rank2Row(comptime T: type, a: Tview(T), k: usize, b: Tview(T), alpha: T) void {
    const herm = comptime isComplex(T);
    const cj = ops.conj;
    for (0..k) |j| {
        const xj = cj(T, a.at(k, j).*);
        const yj = cj(T, b.at(k, j).*);
        if (ops.isZero(T, xj) and ops.isZero(T, yj)) {
            if (comptime herm) a.at(j, j).* = scalar(T, a.at(j, j).re);
            continue;
        }
        const t1 = ops.mul(T, alpha, cj(T, yj));
        const t2 = cj(T, ops.mul(T, alpha, xj));
        if (comptime herm) {
            const d = ops.add(T, ops.mul(T, xj, t1), ops.mul(T, yj, t2));
            a.at(j, j).* = scalar(T, a.at(j, j).re + d.re);
        }
        for (if (herm) j + 1 else j..k) |i| {
            a.at(i, j).* = ops.add(T, a.at(i, j).*, ops.add(T, ops.mul(T, cj(T, a.at(k, i).*), t1), ops.mul(T, cj(T, b.at(k, i).*), t2)));
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
    comptime ops.requireFloat(T, "pbtrf");
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
        const ajj0 = ops.re(T, a.at(j, j).*);
        if (!(ajj0 > 0)) {
            a.at(j, j).* = scalar(T, ajj0);
            return j + 1;
        }
        const ajj = @sqrt(ajj0);
        a.at(j, j).* = scalar(T, ajj);
        const kn = @min(kd, n - j - 1);
        if (kn == 0) continue;
        const r = 1 / ajj;
        for (j + 1..j + 1 + kn) |i| a.at(i, j).* = ops.mulReal(T, a.at(i, j).*, r);
        // `xSYR`, `xHER` with -1 on the window
        for (j + 1..j + 1 + kn) |c| {
            const x = a.at(c, j).*;
            if (ops.isZero(T, x)) {
                if (comptime herm) a.at(c, c).* = scalar(T, a.at(c, c).re);
                continue;
            }
            const t = ops.neg(T, ops.conj(T, x));
            if (comptime herm) a.at(c, c).* = scalar(T, a.at(c, c).re + ops.mul(T, x, t).re);
            for (if (herm) c + 1 else c..j + 1 + kn) |i| a.at(i, c).* = ops.add(T, a.at(i, c).*, ops.mul(T, a.at(i, j).*, t));
        }
    }
    return 0;
}
