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

// Each `...Error` is its routine's bad arguments by plain names (issues
// #313 to #319), `error.OutOfMemory` and `error.Canceled` where it can
// meet them, and the failures its doc names.
/// `potrf`: `error.NotPositiveDefinite`.
pub const PotrfError = cholesky.PotrfError;
/// `potf2`: `error.NotPositiveDefinite`; it allocates nothing.
pub const Potf2Error = cholesky.Potf2Error;
pub const PotrsError = cholesky.PotrsError;
/// `potri`: `error.Singular`.
pub const PotriError = cholesky.PotriError;
/// `pbtrf`: `error.NotPositiveDefinite`.
pub const PbtrfError = relapack.PbtrfError;
pub const potrf = if (options.relapack) relapack.potrf else cholesky.potrf;
pub const potf2 = cholesky.potf2;
pub const potrs = cholesky.potrs;
pub const potri = cholesky.potri;
pub const pbtrf = relapack.pbtrf;

/// `getrf`: returns the first zero `U(i, i)`'s `i` (0-based), or null.
pub const GetrfError = lu.GetrfError;
/// `getrs`: a pivot outside `0..n - 1` (0-based) is `error.PivotOutOfRange`.
pub const GetrsError = lu.GetrsError;
/// `gesv`: `error.Singular`.
pub const GesvError = lu.GesvError;
/// `getri`: `error.Singular`.
pub const GetriError = lu.GetriError;
/// `gbtrf`: returns the first zero `U(i, i)`'s `i` (0-based), or null.
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

/// `sytrf`, `hetrf`, `sytrfRook`, `hetrfRook`: they return the first row of
/// the first zero block of `D` (0-based), or null.
pub const SytrfError = sytrf_factor.SytrfError;
/// One entry of `sytrf`'s pivots: a row (0-based) and its block of `D`.
pub const Pivot = sytrf_factor.Pivot;
/// `sygst`'s and `hegst`'s eigenproblem, LAPACK's `itype`.
pub const Problem = relapack.Problem;
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

/// `geev`: a NaN in `A` is `error.NanInA`; `error.NoConvergence` (also for
/// an infinity in `A`), the count of eigenvalues not found in the last
/// argument.
pub const GeevError = evd.GeevError;
/// `geev`'s eigenvalues: `.{ .re = wr, .im = wi }` real, `w` complex.
pub const Eigenvalues = schur.Eigenvalues;
/// `hseqr`: its block `lo..hi` 0-based, `hi` past the end, as `gebal`
/// returns it; `error.NoConvergence`.
pub const HseqrError = schur.HseqrError;
pub const HseqrJob = schur.HseqrJob;
pub const Compz = schur.Compz;
/// `gees`: returns `sdim`; a NaN in `A` is `error.NanInA`;
/// `error.NoConvergence`, `error.SelectChanged`, `error.SwapRejected`.
pub const GeesError = schur.GeesError;
pub const GeesSelect = schur.GeesSelect;
/// `trevc3`: returns the columns used (the reference's `m`).
pub const Trevc3Error = trevc_vectors.Trevc3Error;
pub const TrevcSide = trevc_vectors.Side;
pub const Howmny = trevc_vectors.Howmny;
/// `gebal`'s and `gebak`'s `JOB`.
pub const BalanceJob = balance.BalanceJob;
/// `gebal`'s result: the block `lo..hi` (0-based, `hi` past the end)
/// outside which `A` is triangular.
pub const BalanceBlock = balance.Block;
/// `gebal`: returns its block; a NaN in `A` is `error.NanInA`. `scale`'s
/// swapped rows are 0-based.
pub const GebalError = balance.GebalError;
/// `gebak`: a swapped row in `scale` past `V`'s is `error.PivotOutOfRange`.
pub const GebakError = balance.GebakError;
/// `gehrd`, `orghr`.
pub const GehrdError = hessenberg.GehrdError;
pub const OrghrError = hessenberg.OrghrError;
/// `trexc`'s and `trsen`'s `COMPQ`.
pub const Compq = schur.Compq;
/// `trsen`'s `JOB`.
pub const TrsenJob = schur.TrsenJob;
/// `trexc`: `ifst` and `ilst` 0-based; `error.SwapRejected`, real only: a
/// swap too ill-conditioned to make.
pub const TrexcError = schur.TrexcError;
/// `trsen`: returns `Reordered`; `error.SwapRejected`.
pub const TrsenError = schur.TrsenError;
/// `trsen`'s result: the cluster's size, and `s` and `sep` as `job` asks.
pub const Reordered = schur.Reordered;
/// `trsyl`: returns `Solved`, `scale` and `perturbed` for close eigenvalues.
pub const TrsylError = trsyl_solve.TrsylError;
/// `trsyl`'s result.
pub const Solved = trsyl_solve.Solved;
/// `tgsyl`: returns `SolvedWithDif`, `scale`, `dif` as `ijob` asks and
/// `perturbed` for close eigenvalues.
pub const TgsylError = tgsyl_solve.TgsylError;
/// `tgsyl`'s result.
pub const SolvedWithDif = tgsyl_solve.SolvedWithDif;
pub const geev = evd.geev;
pub const hseqr = schur.hseqr;
pub const gees = schur.gees;
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

/// `gesvd`: both jobs `o` is `error.ConflictingJobs`; `error.NoConvergence`,
/// its count of unconverged superdiagonals in the last argument.
pub const GesvdError = svd.GesvdError;
/// `gesdd`: a NaN in `A` is `error.NanInA`; `error.NoConvergence`.
pub const GesddError = svd.GesddError;
/// `gesvd`'s and `gesdd`'s jobs.
pub const Job = svd.Job;
pub const gesvd = svd.gesvd;
pub const gesdd = svd.gesdd;

/// `lange`'s and `gecon`'s `NORM`.
pub const Norm = norm.Norm;
/// `lange`: a slice too short or short `lda`, or `error.OutOfMemory`.
pub const LangeError = norm.LangeError;
/// `gecon`: returns `rcond`; a negative, NaN or infinite `anorm` is
/// `error.NormOutOfRange`; `error.EstimateFailed`.
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
    _ = norm_tests;
    _ = trexc_tests;
    _ = bidiag_svd_tests;
    _ = lu_tests;
    _ = lu_solve_tests;
    _ = svd_tests;
    _ = balance_tests;
    _ = cholesky_tests;
    _ = trevc_tests;
    _ = triangular_tests;
    _ = tgsyl_tests;
    _ = geev_tests;
    _ = hseqr_tests;
    _ = stebr_tests;
}

const std = @import("std");
const testing = std.testing;
const check = @import("../blas/_test/check.zig");
const relapack_tests = @import("_test/relapack.zig");
const evd_tests = @import("_test/evd.zig");
const norm_tests = @import("_test/norm.zig");
const trexc_tests = @import("_test/trexc.zig");
const bidiag_svd_tests = @import("_test/bidiag_svd.zig");
const lu_tests = @import("_test/lu.zig");
const lu_solve_tests = @import("_test/lu_solve.zig");
const svd_tests = @import("_test/svd.zig");
const balance_tests = @import("_test/balance.zig");
const cholesky_tests = @import("_test/cholesky.zig");
const trevc_tests = @import("_test/trevc.zig");
const triangular_tests = @import("_test/triangular.zig");
const tgsyl_tests = @import("_test/tgsyl.zig");
const geev_tests = @import("_test/geev.zig");
const hseqr_tests = @import("_test/hseqr.zig");
const stebr_tests = @import("_test/stebr.zig");

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
    var ipiv: [2]Pivot = undefined;
    try testing.expectEqual(@as(?usize, null), try sytrf(al, io, f64, .l, 2, &s, 2, &ipiv));
    try near(&.{ 4, 0.25, 1, 3 - 0.25 }, &s);
    try testing.expectEqualSlices(Pivot, &.{ .of(0, .one), .of(1, .one) }, &ipiv);
    var h = [_]types.Complex64{ .init(4, 0), .init(1, 1), .init(1, -1), .init(3, 0) };
    _ = try hetrf(al, io, types.Complex64, .l, 2, &h, 2, &ipiv);
    try testing.expectApproxEqAbs(@as(f64, 2.5), h[3].re, 1e-14);
    s = .{ 4, 1, 1, 3 };
    _ = try sytrfRook(al, io, f64, .u, 2, &s, 2, &ipiv);
    try testing.expectApproxEqAbs(@as(f64, 3), s[3], 1e-14);
    h = .{ .init(4, 0), .init(1, 1), .init(1, -1), .init(3, 0) };
    _ = try hetrfRook(al, io, types.Complex64, .u, 2, &h, 2, &ipiv);
    try testing.expectApproxEqAbs(@as(f64, 3), h[3].re, 1e-14);
    // sygst, `A x = λ B x`: B = 2 I factored is L = sqrt(2) I, so C = A / 2
    var a = [_]f64{ 4, 2, 2, 6 };
    const l = [_]f64{ @sqrt(2.0), 0, 0, @sqrt(2.0) };
    try sygst(al, io, f64, .ax_lbx, .l, 2, &a, 2, &l, 2);
    try near(&.{ 2, 1 }, &.{ a[0], a[1] });
    try testing.expectApproxEqAbs(@as(f64, 3), a[3], 1e-14);
    var ha = [_]types.Complex64{ .init(4, 0), .init(2, 0), .init(2, 0), .init(6, 0) };
    const hl = [_]types.Complex64{ .init(@sqrt(2.0), 0), .init(0, 0), .init(0, 0), .init(@sqrt(2.0), 0) };
    try hegst(al, io, types.Complex64, .ax_lbx, .u, 2, &ha, 2, &hl, 2);
    try testing.expectApproxEqAbs(@as(f64, 3), ha[3].re, 1e-14);
    // pbtrf: [4 2; 2 5] in upper band storage (kd 1) is U = [2 1; 0 2]
    var ab = [_]f64{ 0, 4, 2, 5 };
    try pbtrf(al, io, f64, .u, 2, 1, &ab, 2);
    try near(&.{ 2, 1, 2 }, &.{ ab[1], ab[2], ab[3] });
    // gbtrf: [1 2; 3 4] (kl = ku = 1, ldab 4) swaps its rows: U = [3 4; 0 2/3]
    var gb = [_]f64{ 0, 0, 1, 3, 0, 2, 4, 0 };
    var gpiv: [2]usize = undefined;
    try testing.expectEqual(@as(?usize, null), try gbtrf(f64, 2, 2, 1, 1, &gb, 4, &gpiv));
    try testing.expectEqualSlices(usize, &.{ 1, 1 }, &gpiv);
    try near(&.{ 3, 1.0 / 3.0, 4, 2.0 / 3.0 }, &.{ gb[2], gb[3], gb[5], gb[6] });
    // trsyl: 2 x + x 3 = 10
    var x = [_]f64{10};
    const ts = try trsyl(al, io, f64, .n, .n, 1, 1, 1, &[_]f64{2}, 1, &[_]f64{3}, 1, &x, 1);
    try near(&.{ 1, 2 }, &.{ ts.scale, x[0] });
    // tgsyl: 2 r + l = 3, r - l = 0 (B = -1, E = 1)
    var c = [_]f64{3};
    var f = [_]f64{0};
    const gs = try tgsyl(al, io, f64, .n, 0, 1, 1, &[_]f64{2}, 1, &[_]f64{-1}, 1, &c, 1, &[_]f64{1}, 1, &[_]f64{1}, 1, &f, 1);
    try near(&.{ 1, 1, 1 }, &.{ gs.scale, c[0], f[0] });
}
