# matlab88 linalg port

Linear algebra kernels ported to Zig from `matlab88.f90`, the 1988 Argonne National
Laboratory "MATLAB" Fortran source (see `matlab88_doc.txt`, `matlab88_help.txt`, and
`readme.txt` in this directory for the original documentation). Only the
LINPACK/EISPACK-derived numerical core was ported - the original REPL (parser,
tokenizer, `<>` stack machine, `PLOT`/`SAVLOD` file I/O) was left behind as
interpreter-specific, not general linalg functionality.

## Importing

From elsewhere in `rstd`, the whole port is reachable as `rstd.linalg.matlab`:

```zig
const rstd = @import("rstd");
const matlab = rstd.linalg.matlab;
```

Within `linalg/` itself (or when depending on this directory directly), import
`matlab/__root__.zig`:

```zig
const matlab = @import("matlab/__root__.zig");
```

`matlab.Cx` is `std.math.Complex(f64)` - the type used for individual complex scalars
throughout this port (see "Data layout" below for how matrices and vectors are
represented).

## Two layers

- **`matlab.api`** - the friendly layer, named after the original REPL's commands
  (`det`, `inv`, `eig`, `chol`, ...). Use this unless you specifically need a
  lower-level primitive. A few of these (`eig`, `cond`, `rank`, `pinv`) take a
  `std.mem.Allocator` because they need scratch space sized at runtime; the rest work
  directly on caller-supplied buffers, same as the low-level kernels.
- **Everything else** (`matlab.lu`, `matlab.svd`, `matlab.qr`, `matlab.blas1`, ...) -
  the low-level kernels, named after their original Fortran routines (lowercased;
  Fortran is case-insensitive, so this renames nothing about the source). Use these
  directly when you want full control over workspace and job flags, or when `api`
  doesn't cover what you need.

`api` is a thin dispatch layer on top of the kernels, not a separate implementation -
e.g. `api.eig` is exactly the Hermitian-vs-general dispatch the original REPL's `EIG`
command did, calling straight into `tridiagonal`/`eigen_symmetric` or
`hessenberg`/`eigen_general` beneath it.

## Data layout

Complex matrices and vectors are **split real/imaginary arrays** (`ar: []f64, ai:
[]f64`), not an interleaved `Complex(f64)` slice - this matches the original Fortran's
`AR`/`AI` convention, and it's what lets the BLAS1 kernels vectorize the common
contiguous-stride case with `@Vector` (real and imaginary lanes are contiguous
runs, not interleaved). A single complex *scalar* (a multiplier, a temporary, a
return value) is `Cx = std.math.Complex(f64)`.

Matrices are **column-major**, with an explicit leading dimension (`lda`), matching
Fortran array storage - `A(i, j)` in the original is `ar[blas1.ix(lda, i, j)]` here.
Row/column indices passed to these functions are **1-based**, again matching the
Fortran source; `blas1.ix(lda, i, j)` converts a 1-based `(i, j)` to the 0-based flat
offset.

Vector arguments to the BLAS1-style kernels (`waxpy`, `wdotcr`, `wnrm2`, ...) take an
explicit stride (`incx`/`incy`, may be negative) and expect the slice to already start
at the first accessed element - exactly like Fortran handing a subroutine the address
of `AR(k, j)`. Most call sites in this port use a stride of 1 (a matrix column) or
`lda` (a matrix row); non-unit strides fall back to a scalar loop instead of the
`@Vector` fast path.

None of the low-level kernels or most of `api` allocate - they read/write
caller-supplied buffers, sized as documented on each function. `api.eig`, `api.cond`,
`api.rank`, and `api.pinv` take a `std.mem.Allocator` for scratch space (and, for
`eig`, the returned eigenvalues/eigenvectors) since that sizing depends on which
runtime code path (Hermitian vs. general) is taken.

## Function reference

### `matlab.api` - friendly (REPL-named) layer

| Name | Signature (elided) | What it does |
| --- | --- | --- |
| `det` | `(ar, ai, lda, n, ipvt) Determinant` | Determinant via LU. `A` destroyed. |
| `inv` | `(ar, ai, lda, n, ipvt, zr, zi) f64` | In-place inverse via LU; returns `rcond`. |
| `chol` | `(ar, ai, lda, n) usize` | Cholesky factor into a clean upper-triangular `R`. |
| `rref` | `(ar, ai, lda, m, n, eps) void` | Reduced row-echelon form in place. Alias of `generators.rref`. |
| `hilb` | `(a, lda, n) void` | Exact inverse Hilbert matrix (real only). Alias of `generators.hilber`. |
| `lu` | `(ar, ai, lda, n, ipvt) usize` | Packed LU factorization. Alias of `lu.wgefa`. |
| `svd` | see `svd.wsvdc` | Alias of `svd.wsvdc`. |
| `qr` | see `qr.wqrdc` | Alias of `qr.wqrdc`. |
| `qrsl` | see `qr.wqrsl` | Alias of `qr.wqrsl`. |
| `eig` | `(allocator, ar, ai, n, want_vectors) !EigResult` | Eigenvalues (+ optional eigenvectors); dispatches Hermitian vs. general. |
| `cond` | `(allocator, ar, ai, n, p) !f64` | 2-norm condition number via SVD. |
| `rank` | `(allocator, ar, ai, n, p, tol) !usize` | Numerical rank via SVD. |
| `pinv` | `(allocator, ar, ai, n, p, tol, result_r, result_i) !void` | Moore-Penrose pseudoinverse via SVD. |

`EigResult` has `values_r`/`values_i` (length `n`) and optional `vectors_r`/`vectors_i`
(`n`x`n`, column-major); call `.deinit(allocator)` to free it.

### Low-level kernels, by module

| Module | Function | Ported from | What it does |
| --- | --- | --- | --- |
| `blas1` | `pythag` | `PYTHAG` | Overflow-safe `sqrt(a^2+b^2)` (aliases `std.math.hypot`). |
| `blas1` | `wmul`, `wdiv` | `WMUL`, `WDIV` | Complex scalar multiply/divide (`wdiv` is overflow-scaled). |
| `blas1` | `wsign`, `wsqrt`, `wlog`, `watan` | `WSIGN`, `WSQRT`, `WLOG`, `WATAN` | Complex scalar sign-copy, sqrt, log, arctangent. |
| `blas1` | `wnrm2`, `wasum`, `iwamax` | `WNRM2`, `WASUM`, `IWAMAX` | Complex vector 2-norm, 1-norm, max-magnitude index. |
| `blas1` | `wrscal`, `wscal`, `waxpy` | `WRSCAL`, `WSCAL`, `WAXPY` | Complex vector scale (real/complex scalar) and `y += s*x`. |
| `blas1` | `wdotur`, `wdotui`, `wdotcr`, `wdotci` | `WDOTUR`/`UI`/`CR`/`CI` | Complex dot products (unconjugated/conjugated, real/imag part). |
| `blas1` | `wcopy`, `wset`, `wswap` | `WCOPY`, `WSET`, `WSWAP` | Complex vector copy, fill, swap. |
| `blas1` | `rset`, `rswap`, `rrot`, `rrotg` | `RSET`, `RSWAP`, `RROT`, `RROTG` | Real vector fill/swap and Givens rotation apply/construct. |
| `blas1` | `ix` | *(new)* | 1-based `(i, j)` → 0-based flat column-major offset. |
| `lu` | `wgefa` | `WGEFA` | Complex LU factorization with partial pivoting. |
| `lu` | `wgesl` | `WGESL` | Solve `A*x=b` or `ctrans(A)*x=b` from `wgefa`'s factors. |
| `lu` | `wgeco` | `WGECO` | `wgefa` + reciprocal condition number estimate. |
| `lu` | `wgedi` | `WGEDI` | Determinant and/or inverse from `wgefa`/`wgeco`'s factors. |
| `cholesky` | `wpofa` | `WPOFA` | Complex Cholesky factorization. |
| `generators` | `rref` | `RREF` | Reduced row-echelon form. |
| `generators` | `hilber` | `HILBER` | Exact inverse Hilbert matrix (closed form). |
| `tridiagonal` | `htridi` | `HTRIDI` | Reduce a Hermitian matrix to real symmetric tridiagonal form. |
| `tridiagonal` | `htribk` | `HTRIBK` | Back-transform tridiagonal eigenvectors to the original Hermitian matrix's. |
| `eigen_symmetric` | `imtql2` | `IMTQL2` | Eigenvalues/eigenvectors of a real symmetric tridiagonal matrix. |
| `hessenberg` | `corth` | `CORTH` | Reduce a general complex matrix to upper Hessenberg form. |
| `eigen_general` | `comqr3` | `COMQR3` | Eigenvalues/Schur form/eigenvectors of a complex Hessenberg matrix (shifted QR). |
| `svd` | `wsvdc` | `WSVDC` | Singular value decomposition of a general complex matrix. |
| `qr` | `wqrdc` | `WQRDC` | QR decomposition with optional column pivoting. |
| `qr` | `wqrsl` | `WQRSL` | Apply `wqrdc`'s output: `Q*y`, `ctrans(Q)*y`, least squares, residual, projection. |

See each function's doc comment for the full parameter contract (buffer sizes, `job`
flag encodings, what's destroyed vs. read-only) - the table above is a map of what
exists, not a substitute for those.

## Example: eigenvalues and eigenvectors

```zig
const std = @import("std");
const matlab = @import("matlab/__root__.zig");

pub fn eigExample(allocator: std.mem.Allocator) !void {
    // A = [[2, 1], [1, 2]] (Hermitian; api.eig will pick the symmetric path).
    const ar = [_]f64{ 2.0, 1.0, 1.0, 2.0 }; // column-major
    const ai = [_]f64{ 0.0, 0.0, 0.0, 0.0 };

    const result = try matlab.api.eig(allocator, &ar, &ai, 2, true);
    defer result.deinit(allocator);

    // result.values_r/values_i: eigenvalues (1 and 3, here).
    // result.vectors_r/vectors_i: matching eigenvectors as columns, if requested.
}
```

## Example: solve a linear system with the low-level kernels

```zig
const matlab = @import("matlab/__root__.zig");
const Cx = matlab.Cx;

pub fn solveExample() void {
    // A = [[1+i, 2], [3, 4-i]], solve A*x = b in place.
    var ar = [_]f64{ 1.0, 3.0, 2.0, 4.0 }; // column-major
    var ai = [_]f64{ 1.0, 0.0, 0.0, -1.0 };
    var ipvt: [2]usize = undefined;

    _ = matlab.lu.wgefa(&ar, &ai, 2, 2, &ipvt);

    var br = [_]f64{ 5.0, 6.0 };
    var bi = [_]f64{ 0.0, 0.0 };
    matlab.lu.wgesl(&ar, &ai, 2, 2, &ipvt, &br, &bi, 0);
    // br/bi now hold the solution x.
}
```
