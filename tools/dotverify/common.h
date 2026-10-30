/* The few names OpenBLAS's x86_64 dot drivers take from its common.h, so they
   compile alone for tools/dotverify/run.sh. */
#ifndef DOTVERIFY_COMMON_H
#define DOTVERIFY_COMMON_H
#include <stddef.h>
typedef long long BLASLONG;
#ifdef DOUBLE
#define FLOAT double
#else
#define FLOAT float
#endif
typedef struct { FLOAT real, imag; } OPENBLAS_COMPLEX_FLOAT_T;
#define OPENBLAS_COMPLEX_FLOAT OPENBLAS_COMPLEX_FLOAT_T
#define OPENBLAS_MAKE_COMPLEX_FLOAT(r, i) ((OPENBLAS_COMPLEX_FLOAT){(r), (i)})
#define CREAL(z) ((z).real)
#define CIMAG(z) ((z).imag)
#endif
