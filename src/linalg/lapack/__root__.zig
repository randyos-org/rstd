//! LAPACK, over rstd's BLAS: OpenBLAS's own LAPACK routines (`lapack/`,
//! `interface/lapack/`) where it has them; past those, faer's methods
//! (`faer/`) behind Reference-LAPACK's interfaces. The conventions are the BLAS
//! routines' (`../blas/readme.md`): Zig functions over column-major data in
//! the Fortran argument order, the element type explicit, sizes `isize`,
//! `(allocator, io, T, ...)`, the BLAS letter enums, and named errors for bad
//! arguments in the reference's `xerbla` order.
//!
//! ReLAPACK's recursive algorithms (`relapack.zig`) serve the routines
//! OpenBLAS takes from the reference, and with `-Drelapack` also `potrf`,
//! `getrf`, `trtri` and `lauum`, as OpenBLAS's `BUILD_RELAPACK` with
//! `RELAPACK_REPLACE`; OpenBLAS's own stay the default.

const cholesky = @import("cholesky.zig");
const lu = @import("lu.zig");
const triangular = @import("triangular.zig");
const svd = @import("faer/svd.zig");
const faer_ops = @import("faer/ops.zig");
const faer_householder = @import("faer/householder.zig");
const faer_qr = @import("faer/qr.zig");
const faer_qr_solve = @import("faer/qr_solve.zig");
const faer_tridiag = @import("faer/tridiag.zig");
const faer_tridiag_evd = @import("faer/tridiag_evd.zig");
const evd = @import("faer/evd.zig");
const faer_hessenberg = @import("faer/hessenberg.zig");
const faer_real_schur = @import("faer/real_schur.zig");
const faer_complex_schur = @import("faer/complex_schur.zig");
const faer_balance = @import("faer/balance.zig");
const faer_bidiag = @import("faer/bidiag.zig");
const faer_bidiag_svd = @import("faer/bidiag_svd.zig");
const relapack = @import("relapack.zig");
const options = @import("lapack_options");

pub const Uplo = @import("../blas/types.zig").Uplo;
pub const Trans = @import("../blas/types.zig").Trans;
pub const Diag = @import("../blas/types.zig").Diag;
pub const Side = @import("../blas/types.zig").Side;

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

/// `gesvd`: bad arguments, `error.NoConvergence`, `error.OutOfMemory` or
/// `error.Canceled`.
pub const GesvdError = svd.GesvdError;
/// `gesdd`: bad arguments (a NaN in `A` is `error.InvalidA`),
/// `error.NoConvergence`, `error.OutOfMemory` or `error.Canceled`.
pub const GesddError = svd.GesddError;
/// `gesvd`'s and `gesdd`'s jobs.
pub const Job = svd.Job;

/// `geqrf`: bad arguments, `error.OutOfMemory` or `error.Canceled`.
pub const GeqrfError = faer_qr.GeqrfError;
/// `orgqr`, `ungqr`: bad arguments, `error.OutOfMemory` or `error.Canceled`.
pub const OrgqrError = faer_qr.OrgqrError;
/// `ormqr`, `unmqr`: bad arguments, `error.OutOfMemory` or `error.Canceled`.
pub const OrmqrError = faer_qr.OrmqrError;
/// `gels`: bad arguments, `error.Singular` (rank deficient),
/// `error.OutOfMemory` or `error.Canceled`.
pub const GelsError = faer_qr_solve.GelsError;

/// `syevd`, `heevd`, `syev`, `heev`: bad arguments, `error.NoConvergence`
/// (also for a NaN or infinity in `A`), `error.OutOfMemory` or
/// `error.Canceled`.
pub const SyevError = evd.SyevError;
/// The eigensolvers' `jobz`, `jobvl`, `jobvr`.
pub const Jobz = evd.Jobz;
/// `geev`: bad arguments (a NaN in `A` is `error.InvalidA`),
/// `error.NoConvergence` (also for an infinity in `A`), `error.OutOfMemory`
/// or `error.Canceled`.
pub const GeevError = evd.GeevError;
/// `geev`'s eigenvalues: `.{ .re = wr, .im = wi }` real, `w` complex.
pub const Eigenvalues = evd.Eigenvalues;

pub const potrf = if (options.relapack) relapack.potrf else cholesky.potrf;
pub const potf2 = cholesky.potf2;
pub const potrs = cholesky.potrs;
pub const potri = cholesky.potri;
pub const getrf = if (options.relapack) relapack.getrf else lu.getrf;
pub const getrs = lu.getrs;
pub const gesv = lu.gesv;
pub const trtri = if (options.relapack) relapack.trtri else triangular.trtri;
pub const trti2 = triangular.trti2;
pub const lauum = if (options.relapack) relapack.lauum else triangular.lauum;
pub const lauu2 = triangular.lauu2;
pub const trtrs = triangular.trtrs;
pub const gesvd = svd.gesvd;
pub const gesdd = svd.gesdd;
pub const geqrf = faer_qr.geqrf;
pub const orgqr = faer_qr.orgqr;
pub const ungqr = faer_qr.orgqr;
pub const ormqr = faer_qr.ormqr;
pub const unmqr = faer_qr.ormqr;
pub const gels = faer_qr_solve.gels;
pub const syevd = evd.syevd;
pub const heevd = evd.syevd;
pub const syev = evd.syev;
pub const heev = evd.syev;
pub const geev = evd.geev;

test {
    _ = cholesky;
    _ = lu;
    _ = triangular;
    _ = faer_ops;
    _ = faer_householder;
    _ = faer_qr;
    _ = faer_qr_solve;
    _ = faer_tridiag;
    _ = faer_tridiag_evd;
    _ = evd;
    _ = faer_hessenberg;
    _ = faer_real_schur;
    _ = faer_complex_schur;
    _ = faer_balance;
    _ = faer_bidiag;
    _ = faer_bidiag_svd;
    _ = svd;
    _ = relapack;
}
