//! LAPACK over rstd's BLAS: OpenBLAS's own routines where it has them,
//! ReLAPACK's, faer's methods and Reference-LAPACK's past those, behind
//! Reference-LAPACK's interfaces in the BLAS's conventions
//! (`../blas/readme.md`). `-Drelapack` takes ReLAPACK's `potrf`, `getrf`,
//! `trtri` and `lauum`, as OpenBLAS's `RELAPACK_REPLACE`.

const options = @import("lapack_options");
const types = @import("../blas/types.zig");
const ops = @import("ops.zig");
const auxiliary = @import("auxiliary.zig");
const cholesky = @import("cholesky.zig");
const lu = @import("lu.zig");
const gbtrf_factor = @import("gbtrf.zig");
const triangular = @import("triangular.zig");
const relapack = @import("relapack.zig");
const indefinite = @import("indefinite.zig");
const sytrf_factor = @import("sytrf.zig");
const householder = @import("householder.zig");
const qr = @import("qr.zig");
const qr_solve = @import("qr_solve.zig");
const lq = @import("lq.zig");
const qr_pivot = @import("qr_pivot.zig");
const lsq = @import("lsq.zig");
const tridiag = @import("tridiag.zig");
const tridiag_evd = @import("tridiag_evd.zig");
const evd = @import("evd.zig");
const balance = @import("balance.zig");
const hessenberg = @import("hessenberg.zig");
const real_schur = @import("real_schur.zig");
const complex_schur = @import("complex_schur.zig");
const schur = @import("schur.zig");
const trevc_vectors = @import("trevc.zig");
const sylvester = @import("sylvester.zig");
const trsyl_solve = @import("trsyl.zig");
const tgsyl_solve = @import("tgsyl.zig");
const bidiag = @import("bidiag.zig");
const bidiag_svd = @import("bidiag_svd.zig");
const svd = @import("svd.zig");
const norm = @import("norm.zig");

pub const Uplo = types.Uplo;
pub const Trans = types.Trans;
pub const Diag = types.Diag;
pub const Side = types.Side;

// Each `...Error` is its routine's bad arguments in the reference's `xerbla`
// order, `error.OutOfMemory` and `error.Canceled` where it can meet them, and
// the failures its doc names.
/// `potrf`, `potf2`: `error.NotPositiveDefinite`.
pub const PotrfError = cholesky.PotrfError;
pub const PotrsError = cholesky.PotrsError;
/// `potri`: `error.Singular`.
pub const PotriError = cholesky.PotriError;
/// `pbtrf`: `error.NotPositiveDefinite`, the minor's order in the last
/// argument.
pub const PbtrfError = relapack.PbtrfError;
pub const potrf = if (options.relapack) relapack.potrf else cholesky.potrf;
pub const potf2 = cholesky.potf2;
pub const potrs = cholesky.potrs;
pub const potri = cholesky.potri;
pub const pbtrf = relapack.pbtrf;

/// `getrf`: `error.Singular`.
pub const GetrfError = lu.GetrfError;
/// `getrs`: a pivot outside `1..n` is `error.InvalidIpiv`.
pub const GetrsError = lu.GetrsError;
/// `gesv`: `error.Singular`.
pub const GesvError = lu.GesvError;
/// `getri`: `error.Singular`.
pub const GetriError = lu.GetriError;
/// `gbtrf`: `error.Singular`, the first zero `U(i, i)` in the last argument.
pub const GbtrfError = gbtrf_factor.GbtrfError;
pub const getrf = if (options.relapack) relapack.getrf else lu.getrf;
pub const getrs = lu.getrs;
pub const gesv = lu.gesv;
pub const getri = lu.getri;
pub const gbtrf = gbtrf_factor.gbtrf;

/// `trtri`: `error.Singular`.
pub const TrtriError = triangular.TrtriError;
/// `trti2`, `lauum`, `lauu2`.
pub const LauumError = triangular.LauumError;
/// `trtrs`: `error.Singular`.
pub const TrtrsError = triangular.TrtrsError;
pub const trtri = if (options.relapack) relapack.trtri else triangular.trtri;
pub const trti2 = triangular.trti2;
pub const lauum = if (options.relapack) relapack.lauum else triangular.lauum;
pub const lauu2 = triangular.lauu2;
pub const trtrs = triangular.trtrs;

/// `sytrf`, `hetrf`, `sytrfRook`, `hetrfRook`: `error.Singular`, a zero block
/// of `D`, its index in the last argument.
pub const SytrfError = sytrf_factor.SytrfError;
/// `sygst`, `hegst`.
pub const SygstError = relapack.SygstError;
pub const sytrf = sytrf_factor.sytrf;
pub const hetrf = sytrf_factor.hetrf;
pub const sytrfRook = sytrf_factor.sytrfRook;
pub const hetrfRook = sytrf_factor.hetrfRook;
pub const sygst = relapack.sygst;
pub const hegst = relapack.hegst;

pub const GeqrfError = qr.GeqrfError;
/// `orgqr`, `ungqr`.
pub const OrgqrError = qr.OrgqrError;
/// `ormqr`, `unmqr`.
pub const OrmqrError = qr.OrmqrError;
/// `gels`: `error.Singular` when rank deficient.
pub const GelsError = qr_solve.GelsError;
pub const GelqfError = lq.GelqfError;
pub const OrglqError = lq.OrglqError;
pub const OrmlqError = lq.OrmlqError;
pub const Geqp3Error = qr_pivot.Geqp3Error;
pub const GelsyError = lsq.GelsyError;
pub const geqrf = qr.geqrf;
pub const orgqr = qr.orgqr;
pub const ungqr = qr.orgqr;
pub const ormqr = qr.ormqr;
pub const unmqr = qr.ormqr;
pub const gels = qr_solve.gels;
pub const gelqf = lq.gelqf;
pub const orglq = lq.orglq;
pub const unglq = lq.orglq;
pub const ormlq = lq.ormlq;
pub const unmlq = lq.ormlq;
pub const geqp3 = qr_pivot.geqp3;
pub const gelsy = lsq.gelsy;

/// `syevd`, `heevd`, `syev`, `heev`: `error.NoConvergence`, also for a NaN or
/// infinity in `A`.
pub const SyevError = evd.SyevError;
/// The eigensolvers' `jobz`, `jobvl`, `jobvr`.
pub const Jobz = evd.Jobz;
pub const syevd = evd.syevd;
pub const heevd = evd.syevd;
pub const syev = evd.syev;
pub const heev = evd.syev;

/// `geev`: a NaN in `A` is `error.InvalidA`; `error.NoConvergence` (also for
/// an infinity in `A`), the count of eigenvalues not found in the last
/// argument.
pub const GeevError = evd.GeevError;
/// `geev`'s eigenvalues: `.{ .re = wr, .im = wi }` real, `w` complex.
pub const Eigenvalues = evd.Eigenvalues;
/// `hseqr`: `error.NoConvergence`.
pub const HseqrError = evd.HseqrError;
pub const HseqrJob = evd.HseqrJob;
pub const Compz = evd.Compz;
/// `gees`: a NaN in `A` is `error.InvalidA`; `error.NoConvergence`,
/// `error.SelectChanged`, `error.SwapRejected`.
pub const GeesError = evd.GeesError;
pub const GeesSelect = evd.GeesSelect;
pub const Trevc3Error = trevc_vectors.Trevc3Error;
pub const TrevcSide = trevc_vectors.Side;
pub const Howmny = trevc_vectors.Howmny;
/// `gebal`'s and `gebak`'s `JOB`.
pub const BalanceJob = balance.BalanceJob;
/// `gebal`: a NaN in `A` is `error.InvalidA`.
pub const GebalError = balance.GebalError;
pub const GebakError = balance.GebakError;
/// `gehrd`, `orghr`.
pub const GehrdError = hessenberg.GehrdError;
pub const OrghrError = hessenberg.OrghrError;
/// `trexc`'s and `trsen`'s `COMPQ`.
pub const Compq = schur.Compq;
/// `trsen`'s `JOB`.
pub const TrsenJob = schur.TrsenJob;
/// `trexc`: `error.SwapRejected`, real only: a swap too ill-conditioned to
/// make.
pub const TrexcError = schur.TrexcError;
/// `trsen`: `error.SwapRejected`.
pub const TrsenError = schur.TrsenError;
/// `trsyl`: `error.Perturbed` for close eigenvalues, `C` and `scale` written.
pub const TrsylError = trsyl_solve.TrsylError;
/// `tgsyl`: `error.Perturbed` for close eigenvalues, `C`, `F`, `scale` and
/// `dif` written.
pub const TgsylError = tgsyl_solve.TgsylError;
pub const geev = evd.geev;
pub const hseqr = evd.hseqr;
pub const gees = evd.gees;
pub const trevc3 = trevc_vectors.trevc3;
pub const trevc = trevc_vectors.trevc3;
pub const gebal = balance.gebal;
pub const gebak = balance.gebak;
pub const gehrd = hessenberg.gehrd;
pub const orghr = hessenberg.orghr;
pub const unghr = hessenberg.orghr;
pub const trexc = schur.trexc;
pub const trsen = schur.trsen;
pub const trsyl = trsyl_solve.trsyl;
pub const tgsyl = tgsyl_solve.tgsyl;

/// `gesvd`: `error.NoConvergence`, its count of unconverged superdiagonals in
/// the last argument.
pub const GesvdError = svd.GesvdError;
/// `gesdd`: a NaN in `A` is `error.InvalidA`; `error.NoConvergence`.
pub const GesddError = svd.GesddError;
/// `gesvd`'s and `gesdd`'s jobs.
pub const Job = svd.Job;
pub const gesvd = svd.gesvd;
pub const gesdd = svd.gesdd;

/// `lange`'s and `gecon`'s `NORM`.
pub const Norm = norm.Norm;
/// `lange`: a negative size or short `lda`, or `error.OutOfMemory`.
pub const LangeError = norm.LangeError;
/// `gecon`: a NaN or infinite `anorm` is `error.InvalidAnorm`;
/// `error.EstimateFailed`.
pub const GeconError = norm.GeconError;
pub const LasetError = auxiliary.LasetError;
pub const lange = norm.lange;
pub const gecon = norm.gecon;
pub const laset = auxiliary.laset;
pub const lapy2 = auxiliary.lapy2;
test {
    _ = cholesky;
    _ = lu;
    _ = gbtrf_factor;
    _ = triangular;
    _ = ops;
    _ = householder;
    _ = qr;
    _ = qr_solve;
    _ = tridiag;
    _ = tridiag_evd;
    _ = evd;
    _ = hessenberg;
    _ = real_schur;
    _ = complex_schur;
    _ = balance;
    _ = trevc_vectors;
    _ = bidiag;
    _ = bidiag_svd;
    _ = svd;
    _ = relapack;
    _ = norm;
    _ = schur;
    _ = lq;
    _ = qr_pivot;
    _ = lsq;
    _ = sylvester;
    _ = trsyl_solve;
    _ = tgsyl_solve;
    _ = indefinite;
    _ = sytrf_factor;
    _ = relapack_tests;
    _ = evd_tests;
}

const std = @import("std");
const testing = std.testing;
const check = @import("../blas/_test/check.zig");
const relapack_tests = @import("_test/relapack.zig");
const evd_tests = @import("_test/evd.zig");

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
