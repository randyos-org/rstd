//! LAPACK, over rstd's BLAS: OpenBLAS's own LAPACK routines (`lapack/`,
//! `interface/lapack/`) where it has them, else Reference-LAPACK's (as
//! OpenBLAS ships it in `lapack-netlib/`). The conventions are the BLAS
//! routines' (`../blas/readme.md`): Zig functions over column-major data in
//! the Fortran argument order, the element type explicit, sizes `isize`,
//! `(allocator, io, T, ...)`, the BLAS letter enums, and named errors for bad
//! arguments in the reference's `xerbla` order.

const cholesky = @import("cholesky.zig");
const lu = @import("lu.zig");

pub const Uplo = @import("../blas/types.zig").Uplo;
pub const Trans = @import("../blas/types.zig").Trans;

/// `potrf`: bad arguments, `error.NotPositiveDefinite`, `error.OutOfMemory`
/// or `error.Canceled`.
pub const PotrfError = cholesky.PotrfError;
/// `potrs`: bad arguments, `error.OutOfMemory` or `error.Canceled`.
pub const PotrsError = cholesky.PotrsError;

/// `getrf`: bad arguments, `error.Singular`, `error.OutOfMemory` or
/// `error.Canceled`.
pub const GetrfError = lu.GetrfError;
/// `getrs`: bad arguments, a pivot outside `1..n` (`error.InvalidIpiv`),
/// `error.OutOfMemory` or `error.Canceled`.
pub const GetrsError = lu.GetrsError;

pub const potrf = cholesky.potrf;
pub const potrs = cholesky.potrs;
pub const getrf = lu.getrf;
pub const getrs = lu.getrs;

test {
    _ = cholesky;
    _ = lu;
}
