// Eigen behind the Fortran BLAS and LAPACK names blas_bench calls, for the
// bench only (tools/blasbench/readme.md). Eigen (MPL 2.0) comes from
// upstream through nixpkgs (eigen.nix); rstd never links this, holds no
// Eigen source, and ports nothing from it. This file is rstd's own. Each
// routine takes the arguments the bench passes and aborts on any other.
// One thread: built without OpenMP.
#include <Eigen/Dense>
#include <algorithm>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace {

template <class T> using Mat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template <class T> using Vec = Eigen::Matrix<T, Eigen::Dynamic, 1>;
using Stride = Eigen::OuterStride<>;
template <class T> using Map = Eigen::Map<Mat<T>, 0, Stride>;
template <class T> using CMap = Eigen::Map<const Mat<T>, 0, Stride>;
template <class T> using Ref = Eigen::Ref<Mat<T>, 0, Stride>;
template <class T> using Real = typename Eigen::NumTraits<T>::Real;

void need(bool ok, const char *what) {
  if (!ok) {
    std::fprintf(stderr, "eigen shim: only %s\n", what);
    std::abort();
  }
}

bool is(const char *c, char want) { return (*c & ~0x20) == want; }

template <class T> Map<T> map(T *a, int r, int c, int ld) { return Map<T>(a, r, c, Stride(ld)); }
template <class T> CMap<T> cmap(const T *a, int r, int c, int ld) { return CMap<T>(a, r, c, Stride(ld)); }

// LAPACK's info for a factor with a zero on its diagonal: the first, from 1
template <class T> int zeroDiag(const Map<T> &a) {
  for (int k = 0; k < a.rows(); ++k)
    if (a(k, k) == T(0)) return k + 1;
  return 0;
}

// LAPACK's row exchanges, from 1, for Eigen's permutation P (P A = L U)
template <class P> void exchanges(const P &perm, int n, int *ipiv) {
  // want[r]: the row of A at row r of P A; cur and pos track the exchanges
  std::vector<int> want(n), cur(n), pos(n);
  for (int i = 0; i < n; ++i) want[perm.indices()[i]] = i;
  for (int i = 0; i < n; ++i) cur[i] = pos[i] = i;
  for (int k = 0; k < n; ++k) {
    const int p = pos[want[k]];
    ipiv[k] = p + 1;
    std::swap(cur[k], cur[p]);
    pos[cur[k]] = k;
    pos[cur[p]] = p;
  }
}

template <class T>
void gemm(const char *ta, const char *tb, const int *m, const int *n, const int *k, const T *alpha, const T *a,
          const int *lda, const T *b, const int *ldb, const T *beta, T *c, const int *ldc) {
  need(is(ta, 'N') && is(tb, 'N'), "gemm NN");
  auto C = map(c, *m, *n, *ldc);
  if (*beta != T(1)) C *= *beta;
  C.noalias() += *alpha * cmap(a, *m, *k, *lda) * cmap(b, *k, *n, *ldb);
}

// real only: Eigen's self-adjoint view is Hermitian, not complex symmetric
template <class T>
void symm(const char *side, const char *uplo, const int *m, const int *n, const T *alpha, const T *a, const int *lda,
          const T *b, const int *ldb, const T *beta, T *c, const int *ldc) {
  need(is(side, 'L') && is(uplo, 'U'), "symm left upper");
  auto C = map(c, *m, *n, *ldc);
  if (*beta != T(1)) C *= *beta;
  C.noalias() += *alpha * cmap(a, *m, *m, *lda).template selfadjointView<Eigen::Upper>() * cmap(b, *m, *n, *ldb);
}

template <class T>
void syrk(const char *uplo, const char *trans, const int *n, const int *k, const T *alpha, const T *a, const int *lda,
          const T *beta, T *c, const int *ldc) {
  need(is(uplo, 'U') && is(trans, 'N'), "syrk upper N");
  auto C = map(c, *n, *n, *ldc);
  if (*beta != T(1)) C.template triangularView<Eigen::Upper>() *= *beta;
  C.template selfadjointView<Eigen::Upper>().rankUpdate(cmap(a, *n, *k, *lda), *alpha);
}

void triArgs(const char *side, const char *uplo, const char *transa, const char *diag) {
  need(is(side, 'L') && is(uplo, 'U') && is(transa, 'N') && is(diag, 'N'), "left upper N non-unit");
}

template <class T>
void trmm(const char *side, const char *uplo, const char *transa, const char *diag, const int *m, const int *n,
          const T *alpha, const T *a, const int *lda, T *b, const int *ldb) {
  triArgs(side, uplo, transa, diag);
  auto B = map(b, *m, *n, *ldb);
  const Mat<T> x = cmap(a, *m, *m, *lda).template triangularView<Eigen::Upper>() * B;
  B = x;
  if (*alpha != T(1)) B *= *alpha;
}

template <class T>
void trsm(const char *side, const char *uplo, const char *transa, const char *diag, const int *m, const int *n,
          const T *alpha, const T *a, const int *lda, T *b, const int *ldb) {
  triArgs(side, uplo, transa, diag);
  auto B = map(b, *m, *n, *ldb);
  if (*alpha != T(1)) B *= *alpha;
  cmap(a, *m, *m, *lda).template triangularView<Eigen::Upper>().solveInPlace(B);
}

template <class T> void potrf(const char *uplo, const int *n, T *a, const int *lda, int *info) {
  need(is(uplo, 'U'), "potrf upper");
  auto M = map(a, *n, *n, *lda);
  Ref<T> A(M);
  Eigen::LLT<Ref<T>, Eigen::Upper> llt(A);
  *info = llt.info() == Eigen::Success ? 0 : 1;
}

template <class T> void getrf(const int *m, const int *n, T *a, const int *lda, int *ipiv, int *info) {
  need(*m == *n, "getrf square");
  auto M = map(a, *n, *n, *lda);
  Ref<T> A(M);
  Eigen::PartialPivLU<Ref<T>> lu(A);
  exchanges(lu.permutationP(), *n, ipiv);
  *info = zeroDiag<T>(M);
}

template <class T>
void getrs(const char *trans, const int *n, const int *nrhs, const T *a, const int *lda, const int *ipiv, T *b,
           const int *ldb, int *info) {
  need(is(trans, 'N'), "getrs N");
  auto A = cmap(a, *n, *n, *lda);
  auto B = map(b, *n, *nrhs, *ldb);
  for (int k = 0; k < *n; ++k)
    if (ipiv[k] - 1 != k) B.row(k).swap(B.row(ipiv[k] - 1));
  A.template triangularView<Eigen::UnitLower>().solveInPlace(B);
  A.template triangularView<Eigen::Upper>().solveInPlace(B);
  *info = 0;
}

template <class T>
void gesv(const int *n, const int *nrhs, T *a, const int *lda, int *ipiv, T *b, const int *ldb, int *info) {
  auto M = map(a, *n, *n, *lda);
  Ref<T> A(M);
  Eigen::PartialPivLU<Ref<T>> lu(A);
  auto B = map(b, *n, *nrhs, *ldb);
  const Mat<T> x = lu.solve(B);
  B = x;
  exchanges(lu.permutationP(), *n, ipiv);
  *info = zeroDiag<T>(M);
}

template <class T> void trtri(const char *uplo, const char *diag, const int *n, T *a, const int *lda, int *info) {
  need(is(uplo, 'U') && is(diag, 'N'), "trtri upper non-unit");
  auto A = map(a, *n, *n, *lda);
  if ((*info = zeroDiag<T>(A)) != 0) return;
  Mat<T> x = Mat<T>::Identity(*n, *n);
  A.template triangularView<Eigen::Upper>().solveInPlace(x);
  A.template triangularView<Eigen::Upper>() = x;
}

// inv(A) = inv(U) inv(U)ᴴ from potrf's U
template <class T> void potri(const char *uplo, const int *n, T *a, const int *lda, int *info) {
  need(is(uplo, 'U'), "potri upper");
  auto A = map(a, *n, *n, *lda);
  if ((*info = zeroDiag<T>(A)) != 0) return;
  Mat<T> ui = Mat<T>::Identity(*n, *n);
  A.template triangularView<Eigen::Upper>().solveInPlace(ui);
  const Mat<T> x = ui.template triangularView<Eigen::Upper>() * ui.adjoint();
  A.template triangularView<Eigen::Upper>() = x;
}

// Eigen's divide-and-conquer SVD, as LAPACK's gesdd: `jobz` S
template <class T>
void gesdd(const char *jobz, const int *m, const int *n, T *a, const int *lda, Real<T> *s, T *u, const int *ldu, T *vt,
           const int *ldvt, T *work, const int *lwork, int *info) {
  *info = 0;
  if (*lwork == -1) {
    work[0] = T(1);
    return;
  }
  need(is(jobz, 'S'), "gesdd S");
  const int mn = std::min(*m, *n);
  Eigen::BDCSVD<Mat<T>> svd(cmap(a, *m, *n, *lda), Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::Map<Vec<Real<T>>>(s, mn) = svd.singularValues();
  map(u, *m, mn, *ldu) = svd.matrixU();
  map(vt, mn, *n, *ldvt) = svd.matrixV().adjoint();
}

template <class T>
void geqrf(const int *m, const int *n, T *a, const int *lda, T *tau, T *work, const int *lwork, int *info) {
  *info = 0;
  if (*lwork == -1) {
    work[0] = T(1);
    return;
  }
  auto M = map(a, *m, *n, *lda);
  Ref<T> A(M);
  Eigen::HouseholderQR<Ref<T>> qr(A);
  Eigen::Map<Vec<T>>(tau, std::min(*m, *n)) = qr.hCoeffs();
}

// Eigen reads the lower triangle: the upper one is made whole first
template <class T>
void eigh(const char *jobz, const char *uplo, const int *n, T *a, const int *lda, Real<T> *w, int *info) {
  need(is(jobz, 'V') && is(uplo, 'U'), "syevd/heevd V upper");
  auto A = map(a, *n, *n, *lda);
  const Mat<T> full = A.template selfadjointView<Eigen::Upper>();
  Eigen::SelfAdjointEigenSolver<Mat<T>> es(full, Eigen::ComputeEigenvectors);
  *info = es.info() == Eigen::Success ? 0 : 1;
  Eigen::Map<Vec<Real<T>>>(w, *n) = es.eigenvalues();
  A = es.eigenvectors();
}

} // namespace

#define STR2(x) #x
#define STR(x) STR2(x)

// the only symbols the library exports (built with -fvisibility=hidden)
#if defined(_WIN32)
#define EXPORT extern "C" __declspec(dllexport)
#else
#define EXPORT extern "C" __attribute__((visibility("default")))
#endif

EXPORT const char *rstd_bench_eigen_version() {
  return STR(EIGEN_WORLD_VERSION) "." STR(EIGEN_MAJOR_VERSION) "." STR(EIGEN_MINOR_VERSION);
}

// the routines every type has
#define ALL(p, T)                                                                                                     \
  EXPORT void p##gemm_(const char *ta, const char *tb, const int *m, const int *n, const int *k, const T *alpha,  \
                           const T *a, const int *lda, const T *b, const int *ldb, const T *beta, T *c,              \
                           const int *ldc) {                                                                          \
    gemm(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);                                                       \
  }                                                                                                                   \
  EXPORT void p##trmm_(const char *side, const char *uplo, const char *ta, const char *diag, const int *m,        \
                           const int *n, const T *alpha, const T *a, const int *lda, T *b, const int *ldb) {         \
    trmm(side, uplo, ta, diag, m, n, alpha, a, lda, b, ldb);                                                          \
  }                                                                                                                   \
  EXPORT void p##trsm_(const char *side, const char *uplo, const char *ta, const char *diag, const int *m,        \
                           const int *n, const T *alpha, const T *a, const int *lda, T *b, const int *ldb) {         \
    trsm(side, uplo, ta, diag, m, n, alpha, a, lda, b, ldb);                                                          \
  }                                                                                                                   \
  EXPORT void p##potrf_(const char *uplo, const int *n, T *a, const int *lda, int *info) {                       \
    potrf(uplo, n, a, lda, info);                                                                                     \
  }                                                                                                                   \
  EXPORT void p##getrf_(const int *m, const int *n, T *a, const int *lda, int *ipiv, int *info) {                \
    getrf(m, n, a, lda, ipiv, info);                                                                                  \
  }                                                                                                                   \
  EXPORT void p##getrs_(const char *trans, const int *n, const int *nrhs, const T *a, const int *lda,             \
                            const int *ipiv, T *b, const int *ldb, int *info) {                                       \
    getrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info);                                                                \
  }                                                                                                                   \
  EXPORT void p##gesv_(const int *n, const int *nrhs, T *a, const int *lda, int *ipiv, T *b, const int *ldb,     \
                           int *info) {                                                                               \
    gesv(n, nrhs, a, lda, ipiv, b, ldb, info);                                                                        \
  }                                                                                                                   \
  EXPORT void p##trtri_(const char *uplo, const char *diag, const int *n, T *a, const int *lda, int *info) {     \
    trtri(uplo, diag, n, a, lda, info);                                                                               \
  }                                                                                                                   \
  EXPORT void p##potri_(const char *uplo, const int *n, T *a, const int *lda, int *info) {                       \
    potri(uplo, n, a, lda, info);                                                                                     \
  }                                                                                                                   \
  EXPORT void p##geqrf_(const int *m, const int *n, T *a, const int *lda, T *tau, T *work, const int *lwork,     \
                            int *info) {                                                                              \
    geqrf(m, n, a, lda, tau, work, lwork, info);                                                                      \
  }

// real: symm, syrk (complex symmetric has no Eigen view), gesdd, syevd
#define REAL(p, T)                                                                                                    \
  ALL(p, T)                                                                                                           \
  EXPORT void p##symm_(const char *side, const char *uplo, const int *m, const int *n, const T *alpha,           \
                           const T *a, const int *lda, const T *b, const int *ldb, const T *beta, T *c,              \
                           const int *ldc) {                                                                          \
    symm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc);                                                      \
  }                                                                                                                   \
  EXPORT void p##syrk_(const char *uplo, const char *trans, const int *n, const int *k, const T *alpha,          \
                           const T *a, const int *lda, const T *beta, T *c, const int *ldc) {                        \
    syrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc);                                                             \
  }                                                                                                                   \
  EXPORT void p##gesdd_(const char *jobz, const int *m, const int *n, T *a, const int *lda, T *s, T *u,          \
                            const int *ldu, T *vt, const int *ldvt, T *work, const int *lwork, int *iwork,           \
                            int *info) {                                                                              \
    (void)iwork;                                                                                                      \
    gesdd(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);                                                \
  }                                                                                                                   \
  EXPORT void p##syevd_(const char *jobz, const char *uplo, const int *n, T *a, const int *lda, T *w, T *work,  \
                            const int *lwork, int *iwork, const int *liwork, int *info) {                            \
    if (*lwork == -1 || *liwork == -1) {                                                                              \
      work[0] = 1;                                                                                                    \
      iwork[0] = 1;                                                                                                   \
      *info = 0;                                                                                                      \
      return;                                                                                                         \
    }                                                                                                                 \
    eigh(jobz, uplo, n, a, lda, w, info);                                                                             \
  }

// complex: gesdd and heevd, with their real workspaces
#define CPLX(p, T, R)                                                                                                 \
  ALL(p, T)                                                                                                           \
  EXPORT void p##gesdd_(const char *jobz, const int *m, const int *n, T *a, const int *lda, R *s, T *u,          \
                            const int *ldu, T *vt, const int *ldvt, T *work, const int *lwork, R *rwork,             \
                            int *iwork, int *info) {                                                                  \
    (void)rwork;                                                                                                      \
    (void)iwork;                                                                                                      \
    gesdd(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);                                                \
  }                                                                                                                   \
  EXPORT void p##heevd_(const char *jobz, const char *uplo, const int *n, T *a, const int *lda, R *w, T *work,  \
                            const int *lwork, R *rwork, const int *lrwork, int *iwork, const int *liwork,            \
                            int *info) {                                                                              \
    if (*lwork == -1 || *lrwork == -1 || *liwork == -1) {                                                             \
      work[0] = 1;                                                                                                    \
      rwork[0] = 1;                                                                                                   \
      iwork[0] = 1;                                                                                                   \
      *info = 0;                                                                                                      \
      return;                                                                                                         \
    }                                                                                                                 \
    eigh(jobz, uplo, n, a, lda, w, info);                                                                             \
  }

REAL(s, float)
REAL(d, double)
CPLX(c, std::complex<float>, float)
CPLX(z, std::complex<double>, double)
