//! The LAPACK routines' interface rules: the argument checks (the
//! reference's, and the Zig-idiomatic ones of issue #313 and after), the
//! vtable, and upstream's thread counts by size, kept so the bits never
//! depend on the worker count.

const std = @import("std");
const sysinfo = @import("builtin");
const blas = @import("../blas/__root__.zig");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const strided = @import("../blas/strided.zig");
const threading = @import("../blas/threading.zig");
const Vtable = @import("../blas/Vtable.zig");

const Float = types.Float;
const isComplex = types.isComplex;
const Uplo = types.Uplo;

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

/// A slice too short for its matrix or vector, or a leading dimension below
/// its rows.
pub const ShapeError = strided.SliceError || strided.LeadingDimensionError;
/// A reversed `uplo` (`v`, `m`), which no LAPACK routine takes.
pub const UploError = error{UnsupportedUplo};
/// A `trans` letter `o` to `d`, which no LAPACK routine takes.
pub const TransError = error{UnsupportedTranspose};
/// A pivot that is not a row of the matrix.
pub const PivotError = error{PivotOutOfRange};
/// An exact zero on a triangular factor's diagonal, which a solve or an
/// inverse cannot divide by.
pub const SingularError = error{Singular};
/// The checks of a routine on one triangle of a square matrix.
pub const TriangleError = UploError || ShapeError;

/// `uplo` is one the routine takes and `a` holds the `n` by `n` matrix.
pub fn holdsTriangle(uplo: Uplo, n: usize, a_len: usize, lda: usize) TriangleError!void {
    if (uplo.isRev()) return error.UnsupportedUplo;
    try util.holds(a_len, n, n, lda);
}

/// `trans` is one the routine takes.
pub fn takes(trans: types.Trans) TransError!void {
    if (trans.conjX()) return error.UnsupportedTranspose;
}

/// More reflections than the matrix they act on holds.
pub const ReflectionsError = error{TooManyReflections};

/// `trans` is `n`, `t` or `c` (complex: `n` or `c`), as the orthogonal and
/// unitary routines take it.
pub fn takesAdjoint(comptime T: type, trans: types.Trans) TransError!void {
    switch (trans) {
        .n, .c => {},
        .t => if (comptime types.isComplex(T)) return error.UnsupportedTranspose,
        .r, .o, .u, .s, .d => return error.UnsupportedTranspose,
    }
}

/// A NaN in `A`.
pub const NanError = error{NanInA};

/// A block `lo..hi` of rows and columns that is reversed or past `n`.
pub const BlockError = error{BlockOutOfRange};

/// `lo..hi` (0-based, `hi` past its end) lies within the order `n` matrix.
pub fn blockIn(lo: usize, hi: usize, n: usize) BlockError!void {
    if (lo > hi or hi > n) return error.BlockOutOfRange;
}

/// `ipiv` holds `n` pivots, each a row of the order `n` matrix.
pub fn pivotsIn(ipiv: []const usize, n: usize) (strided.SliceError || PivotError)!void {
    if (ipiv.len < n) return error.SliceTooShort;
    for (ipiv[0..n]) |p| if (p >= n) return error.PivotOutOfRange;
}

/// The first of two zero pivots: `first`, else `later`'s at `offset`.
pub fn firstOf(first: ?usize, later: ?usize, offset: usize) ?usize {
    return first orelse if (later) |k| k + offset else null;
}

/// `dtb_entries` of a hand-built vtable that has none (`blas2/triangular.zig`).
pub const default_dtb = 64;

pub const table = blas.tableOf;

/// `dtb_entries` of `vt`, or `default_dtb`.
pub fn dtbOf(vt: *const Vtable) usize {
    return if (vt.dtb_entries > 0) @intCast(vt.dtb_entries) else default_dtb;
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

/// `threading.team_grain` where the choice of an algorithm reads it, so the
/// bits are one on every system.
pub const block_grain: u64 = 1_000_000;

/// `threading.workers` for a size cap, at most `threading.chunks_max`.
pub fn workers(size_max: u64) u32 {
    // selects the core first, which fills `blas.config`
    _ = blas.coreName();
    const w = threading.workers(&blas.config, std.math.cast(u32, size_max) orelse std.math.maxInt(u32));
    return @min(w, threading.chunks_max);
}

/// LAPACK's `JOBZ`, `JOBVL`, `JOBVR`, `JOBVS`: vectors (`V`) or not (`N`).
pub const VectorJob = enum { none, compute };

/// LAPACK's `EQUED`: what `laqge` and `gesvx` scaled, `diag(R) A diag(C)`.
pub const Equilibration = enum { none, rows, columns, both };

/// LAPACK's `FACT`: `gesvx` factors `A` (`N`), takes its given factor
/// (`F`), or equilibrates `A` first (`E`).
pub const FactorizationInput = enum { factor, factored, equilibrate_and_factor };

/// LAPACK's `RANGE`: every eigenvalue (`A`), those in `(lo, hi]` (`V`), or
/// the `first`-th to `last`-th smallest, 0-based and inclusive (`I`).
pub fn SpectrumRange(comptime F: type) type {
    return union(enum) { all, values: struct { lo: F, hi: F }, indices: struct { first: usize, last: usize } };
}

/// `stebz`'s `ORDER`: every eigenvalue ascending (`E`), or by split block,
/// each block's ascending (`B`).
pub const SpectrumOrder = enum { ascending, by_block };
