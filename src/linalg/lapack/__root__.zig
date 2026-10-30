//! LAPACK, over rstd's BLAS: OpenBLAS's own LAPACK routines (`lapack/`,
//! `interface/lapack/`) where it has them, else Reference-LAPACK's (as
//! OpenBLAS ships it in `lapack-netlib/`). The conventions are the BLAS
//! routines' (`../blas/readme.md`): Zig functions over column-major data in
//! the Fortran argument order, the element type explicit, sizes `isize`,
//! `(allocator, io, T, ...)`, the BLAS letter enums, and named errors for bad
//! arguments in the reference's `xerbla` order.

const cholesky = @import("cholesky.zig");
const lu = @import("lu.zig");
const triangular = @import("triangular.zig");

pub const Uplo = @import("../blas/types.zig").Uplo;
pub const Trans = @import("../blas/types.zig").Trans;
pub const Diag = @import("../blas/types.zig").Diag;

/// `potrf`, `potf2`: bad arguments, `error.NotPositiveDefinite`,
/// `error.OutOfMemory` or `error.Canceled`.
pub const PotrfError = cholesky.PotrfError;
/// `potrs`: bad arguments, `error.OutOfMemory` or `error.Canceled`.
pub const PotrsError = cholesky.PotrsError;
/// `potri`: bad arguments, `error.Singular`, `error.OutOfMemory` or
/// `error.Canceled`.
pub const PotriError = cholesky.PotriError;

/// `getrf`: bad arguments, `error.Singular`, `error.OutOfMemory` or
/// `error.Canceled`.
pub const GetrfError = lu.GetrfError;
/// `getrs`: bad arguments, a pivot outside `1..n` (`error.InvalidIpiv`),
/// `error.OutOfMemory` or `error.Canceled`.
pub const GetrsError = lu.GetrsError;
/// `gesv`: bad arguments, `error.Singular`, `error.OutOfMemory` or
/// `error.Canceled`.
pub const GesvError = lu.GesvError;

/// `trtri`: bad arguments, `error.Singular`, `error.OutOfMemory` or
/// `error.Canceled`.
pub const TrtriError = triangular.TrtriError;
/// `trti2`, `lauum`, `lauu2`: bad arguments, `error.OutOfMemory` or
/// `error.Canceled`.
pub const LauumError = triangular.LauumError;
/// `trtrs`: bad arguments, `error.Singular`, `error.OutOfMemory` or
/// `error.Canceled`.
pub const TrtrsError = triangular.TrtrsError;

pub const potrf = cholesky.potrf;
pub const potf2 = cholesky.potf2;
pub const potrs = cholesky.potrs;
pub const potri = cholesky.potri;
pub const getrf = lu.getrf;
pub const getrs = lu.getrs;
pub const gesv = lu.gesv;
pub const trtri = triangular.trtri;
pub const trti2 = triangular.trti2;
pub const lauum = triangular.lauum;
pub const lauu2 = triangular.lauu2;
pub const trtrs = triangular.trtrs;

test {
    _ = cholesky;
    _ = lu;
    _ = triangular;
}
