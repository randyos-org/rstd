//! The LAPACK routines' interface rules: the reference's argument checks,
//! the vtable, and upstream's thread counts by size, kept so the bits never
//! depend on the worker count.

const std = @import("std");
const sysinfo = @import("builtin");
const blas = @import("../blas/__root__.zig");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const threading = @import("../blas/threading.zig");
const portable = @import("../blas/portable.zig");
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

/// `dtb_entries` of a hand-built vtable that has none (`blas2/triangular.zig`).
pub const default_dtb = 64;

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

/// `threading.workers` for a size cap, at most `threading.chunks_max`.
pub fn workers(size_max: u64) u32 {
    // selects the core first, which fills `blas.config`
    _ = blas.coreName();
    const w = threading.workers(&blas.config, std.math.cast(u32, size_max) orelse std.math.maxInt(u32));
    return @min(w, threading.chunks_max);
}
