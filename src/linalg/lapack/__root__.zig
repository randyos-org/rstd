//! LAPACK, over rstd's BLAS: OpenBLAS's own LAPACK routines (`lapack/`,
//! `interface/lapack/`) where it has them, else Reference-LAPACK's (as
//! OpenBLAS ships it in `lapack-netlib/`). The conventions are the BLAS
//! routines' (`../blas/readme.md`): Zig functions over column-major data in
//! the Fortran argument order, the element type explicit, sizes `isize`,
//! `(allocator, io, T, ...)`, the BLAS letter enums, and named errors for bad
//! arguments in the reference's `xerbla` order.

const cholesky = @import("cholesky.zig");

pub const Uplo = @import("../blas/types.zig").Uplo;

/// `potrf`: bad arguments, `error.NotPositiveDefinite`, `error.OutOfMemory`
/// or `error.Canceled`.
pub const PotrfError = cholesky.PotrfError;
/// `potrs`: bad arguments, `error.OutOfMemory` or `error.Canceled`.
pub const PotrsError = cholesky.PotrsError;

pub const potrf = cholesky.potrf;
pub const potrs = cholesky.potrs;

test {
    _ = cholesky;
}
