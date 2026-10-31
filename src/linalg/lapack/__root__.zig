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
const svd = @import("svd.zig");
const faer_ops = @import("ops.zig");
const faer_householder = @import("householder.zig");
const faer_qr = @import("qr.zig");
const faer_qr_solve = @import("qr_solve.zig");
const faer_tridiag = @import("tridiag.zig");
const faer_tridiag_evd = @import("tridiag_evd.zig");
const evd = @import("evd.zig");
const faer_hessenberg = @import("hessenberg.zig");
const faer_real_schur = @import("real_schur.zig");
const faer_complex_schur = @import("complex_schur.zig");
const faer_balance = @import("balance.zig");
const faer_bidiag = @import("bidiag.zig");
const faer_bidiag_svd = @import("bidiag_svd.zig");
const relapack = @import("relapack.zig");
const norm = @import("norm.zig");
const schur = @import("schur.zig");
const lq = @import("lq.zig");
const qr_pivot = @import("qr_pivot.zig");
const lsq = @import("lsq.zig");
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
/// `getri`: bad arguments, `error.Singular`, `error.OutOfMemory` or
/// `error.Canceled`.
pub const GetriError = lu.GetriError;

/// `trtri`: bad arguments, `error.Singular`, `error.OutOfMemory` or
/// `error.Canceled`.
pub const TrtriError = triangular.TrtriError;
/// `trti2`, `lauum`, `lauu2`: bad arguments, `error.OutOfMemory` or
/// `error.Canceled`.
pub const LauumError = triangular.LauumError;
/// `trtrs`: bad arguments, `error.Singular`, `error.OutOfMemory` or
/// `error.Canceled`.
pub const TrtrsError = triangular.TrtrsError;

/// `gesvd`: bad arguments, `error.NoConvergence` (its count of unconverged
/// superdiagonals in the last argument), `error.OutOfMemory` or
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
/// `error.NoConvergence` (also for an infinity in `A`; the count of
/// eigenvalues not found in the last argument), `error.OutOfMemory` or
/// `error.Canceled`.
pub const GeevError = evd.GeevError;
/// `geev`'s eigenvalues: `.{ .re = wr, .im = wi }` real, `w` complex.
pub const Eigenvalues = evd.Eigenvalues;
/// `lange`'s and `gecon`'s `NORM`.
pub const Norm = norm.Norm;
/// `lange`: a negative size or short `lda`, or `error.OutOfMemory`.
pub const LangeError = norm.LangeError;
/// `gecon`: bad arguments (a NaN or infinite `anorm` is
/// `error.InvalidAnorm`), `error.EstimateFailed`, `error.OutOfMemory` or
/// `error.Canceled`.
pub const GeconError = norm.GeconError;
/// `laset`: bad arguments.
pub const LasetError = norm.LasetError;
/// `gebal`'s and `gebak`'s `JOB`.
/// `gehrd`, `orghr`: bad arguments, `error.OutOfMemory` or `error.Canceled`.
pub const GehrdError = faer_hessenberg.GehrdError;
pub const OrghrError = faer_hessenberg.OrghrError;
pub const BalanceJob = faer_balance.BalanceJob;
/// `gebal`: bad arguments, a NaN in `A` (`error.InvalidA`) or
/// `error.Canceled`.
pub const GebalError = faer_balance.GebalError;
/// `gebak`: bad arguments or `error.Canceled`.
pub const GebakError = faer_balance.GebakError;
/// `trexc`'s and `trsen`'s `COMPQ`.
pub const Compq = schur.Compq;
/// `trsen`'s `JOB`.
pub const TrsenJob = schur.TrsenJob;
/// `trexc`: bad arguments, `error.SwapRejected` (real: a swap too
/// ill-conditioned to make) or `error.Canceled`.
pub const TrexcError = schur.TrexcError;
/// `trsen`: bad arguments, `error.SwapRejected`, `error.OutOfMemory` or
/// `error.Canceled`.
pub const TrsenError = schur.TrsenError;
/// `gelqf`: bad arguments, `error.OutOfMemory` or `error.Canceled`.
pub const GelqfError = lq.GelqfError;
/// `orglq`: bad arguments, `error.OutOfMemory` or `error.Canceled`.
pub const OrglqError = lq.OrglqError;
/// `ormlq`: bad arguments, `error.OutOfMemory` or `error.Canceled`.
pub const OrmlqError = lq.OrmlqError;
/// `geqp3`: bad arguments, `error.OutOfMemory` or `error.Canceled`.
pub const Geqp3Error = qr_pivot.Geqp3Error;
/// `gelsy`: bad arguments, `error.OutOfMemory` or `error.Canceled`.
pub const GelsyError = lsq.GelsyError;

/// `sytrf`, `hetrf`, `sytrfRook`, `hetrfRook`: bad arguments,
/// `error.Singular` (a zero block of `D`, its index in the last argument),
/// `error.OutOfMemory` or `error.Canceled`.
pub const SytrfError = relapack.SytrfError;
/// `sygst`, `hegst`: bad arguments, `error.OutOfMemory` or `error.Canceled`.
pub const SygstError = relapack.SygstError;
/// `pbtrf`: bad arguments, `error.NotPositiveDefinite` (the minor's order in
/// the last argument), `error.OutOfMemory` or `error.Canceled`.
pub const PbtrfError = relapack.PbtrfError;
/// `gbtrf`: bad arguments, `error.Singular` (the first zero `U(i, i)` in the
/// last argument), `error.OutOfMemory` or `error.Canceled`.
pub const GbtrfError = relapack.GbtrfError;
/// `trsyl`: bad arguments, `error.Perturbed` (close eigenvalues; `C` and
/// `scale` written), `error.OutOfMemory` or `error.Canceled`.
pub const TrsylError = relapack.TrsylError;
/// `tgsyl`: bad arguments, `error.Perturbed` (close eigenvalues; `C`, `F`,
/// `scale` and `dif` written), `error.OutOfMemory` or `error.Canceled`.
pub const TgsylError = relapack.TgsylError;

pub const potrf = if (options.relapack) relapack.potrf else cholesky.potrf;
pub const potf2 = cholesky.potf2;
pub const potrs = cholesky.potrs;
pub const potri = cholesky.potri;
pub const getrf = if (options.relapack) relapack.getrf else lu.getrf;
pub const getrs = lu.getrs;
pub const gesv = lu.gesv;
pub const getri = lu.getri;
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
pub const lange = norm.lange;
pub const gecon = norm.gecon;
pub const laset = norm.laset;
pub const lapy2 = schur.lapy2;
pub const gehrd = faer_hessenberg.gehrd;
pub const orghr = faer_hessenberg.orghr;
pub const unghr = faer_hessenberg.orghr;
pub const gebal = faer_balance.gebal;
pub const gebak = faer_balance.gebak;
pub const trexc = schur.trexc;
pub const trsen = schur.trsen;
pub const gelqf = lq.gelqf;
pub const orglq = lq.orglq;
pub const unglq = lq.orglq;
pub const ormlq = lq.ormlq;
pub const unmlq = lq.ormlq;
pub const geqp3 = qr_pivot.geqp3;
pub const gelsy = lsq.gelsy;

pub const sytrf = relapack.sytrf;
pub const hetrf = relapack.hetrf;
pub const sytrfRook = relapack.sytrfRook;
pub const hetrfRook = relapack.hetrfRook;
pub const sygst = relapack.sygst;
pub const hegst = relapack.hegst;
pub const pbtrf = relapack.pbtrf;
pub const gbtrf = relapack.gbtrf;
pub const trsyl = relapack.trsyl;
pub const tgsyl = relapack.tgsyl;

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
    _ = norm;
    _ = schur;
    _ = lq;
    _ = qr_pivot;
    _ = lsq;
}

const std = @import("std");
const testing = std.testing;
const types = @import("../blas/types.zig");
const check = @import("../blas/_test/blas1/check.zig");

test "lapack: ReLAPACK's routines through their public names, 2 x 2 smoke tests" {
    const al = testing.allocator;
    const io = check.io;
    const near = struct {
        fn f(want: []const f64, got: []const f64) !void {
            for (want, got) |w, g| try testing.expectApproxEqAbs(w, g, 1e-14);
        }
    }.f;
    // sytrf: A = [4 1; 1 3] is positive definite, so 1 x 1 pivots, no swap
    var s = [_]f64{ 4, 1, 1, 3 };
    var ipiv: [2]isize = undefined;
    try sytrf(al, io, f64, .l, 2, &s, 2, &ipiv, null);
    try near(&.{ 4, 0.25, 1, 3 - 0.25 }, &s);
    try testing.expectEqualSlices(isize, &.{ 1, 2 }, &ipiv);
    var h = [_]types.Complex64{ .init(4, 0), .init(1, 1), .init(1, -1), .init(3, 0) };
    try hetrf(al, io, types.Complex64, .l, 2, &h, 2, &ipiv, null);
    try testing.expectApproxEqAbs(@as(f64, 2.5), h[3].re, 1e-14);
    s = .{ 4, 1, 1, 3 };
    try sytrfRook(al, io, f64, .u, 2, &s, 2, &ipiv, null);
    try testing.expectApproxEqAbs(@as(f64, 3), s[3], 1e-14);
    h = .{ .init(4, 0), .init(1, 1), .init(1, -1), .init(3, 0) };
    try hetrfRook(al, io, types.Complex64, .u, 2, &h, 2, &ipiv, null);
    try testing.expectApproxEqAbs(@as(f64, 3), h[3].re, 1e-14);
    // sygst, itype 1: B = 2 I factored is L = sqrt(2) I, so C = A / 2
    var a = [_]f64{ 4, 2, 2, 6 };
    const l = [_]f64{ @sqrt(2.0), 0, 0, @sqrt(2.0) };
    try sygst(al, io, f64, 1, .l, 2, &a, 2, &l, 2);
    try near(&.{ 2, 1 }, &.{ a[0], a[1] });
    try testing.expectApproxEqAbs(@as(f64, 3), a[3], 1e-14);
    var ha = [_]types.Complex64{ .init(4, 0), .init(2, 0), .init(2, 0), .init(6, 0) };
    const hl = [_]types.Complex64{ .init(@sqrt(2.0), 0), .init(0, 0), .init(0, 0), .init(@sqrt(2.0), 0) };
    try hegst(al, io, types.Complex64, 1, .u, 2, &ha, 2, &hl, 2);
    try testing.expectApproxEqAbs(@as(f64, 3), ha[3].re, 1e-14);
    // pbtrf: [4 2; 2 5] in upper band storage (kd 1) is U = [2 1; 0 2]
    var ab = [_]f64{ 0, 4, 2, 5 };
    try pbtrf(al, io, f64, .u, 2, 1, &ab, 2, null);
    try near(&.{ 2, 1, 2 }, &.{ ab[1], ab[2], ab[3] });
    // gbtrf: [1 2; 3 4] (kl = ku = 1, ldab 4) swaps its rows: U = [3 4; 0 2/3]
    var gb = [_]f64{ 0, 0, 1, 3, 0, 2, 4, 0 };
    try gbtrf(al, io, f64, 2, 2, 1, 1, &gb, 4, &ipiv, null);
    try testing.expectEqualSlices(isize, &.{ 2, 2 }, &ipiv);
    try near(&.{ 3, 1.0 / 3.0, 4, 2.0 / 3.0 }, &.{ gb[2], gb[3], gb[5], gb[6] });
    // trsyl: 2 x + x 3 = 10
    var scale: f64 = 0;
    var x = [_]f64{10};
    try trsyl(al, io, f64, .n, .n, 1, 1, 1, &[_]f64{2}, 1, &[_]f64{3}, 1, &x, 1, &scale);
    try near(&.{ 1, 2 }, &.{ scale, x[0] });
    // tgsyl: 2 r + l = 3, r - l = 0 (B = -1, E = 1)
    var c = [_]f64{3};
    var f = [_]f64{0};
    var dif: f64 = 0;
    try tgsyl(al, io, f64, .n, 0, 1, 1, &[_]f64{2}, 1, &[_]f64{-1}, 1, &c, 1, &[_]f64{1}, 1, &[_]f64{1}, 1, &f, 1, &scale, &dif);
    try near(&.{ 1, 1, 1 }, &.{ scale, c[0], f[0] });
}
