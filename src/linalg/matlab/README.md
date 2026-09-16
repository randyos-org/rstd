# matlab88 linalg port

A Zig port of `matlab88.f90`'s function library - the 1988 Argonne National Laboratory
"MATLAB" Fortran source (see `ref/matlab88_doc.txt`, `ref/matlab88_help.txt`, and
`ref/readme.txt` for the original documentation). The goal is every function a MATLAB88
user could call, minus the REPL itself: the parser, tokenizer, `<>` stack machine,
`FOR`/`WHILE`/`IF`, `SAVE`/`LOAD`/`DIARY`/`PLOT`/`PRINT` file and terminal I/O, and
`HELP` text are all left out, since none of that has meaning outside an interactive
session.

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
  (`det`, `inv`, `eig`, `chol`, `kron`, `roots`, ...). Use this unless you specifically
  need a lower-level primitive. A handful of these (`eig`, `cond`, `rank`, `pinv`,
  `orth`, `roots`, and the `expm`/`logm`/... matrix-function family) take a
  `std.mem.Allocator` because they need scratch space sized at runtime; the rest work
  directly on caller-supplied buffers, same as the low-level kernels.
- **Everything else** (`matlab.lu`, `matlab.svd`, `matlab.qr`, `matlab.blas1`,
  `matlab.generators`, ...) - the low-level kernels, named after their original Fortran
  routines (lowercased; Fortran is case-insensitive, so this renames nothing about the
  source). Use these directly when you want full control over workspace and job flags,
  or when `api` doesn't cover what you need.

`api` is a thin dispatch layer on top of the kernels, not a separate implementation -
e.g. `api.eig` is exactly the Hermitian-vs-general dispatch the original REPL's `EIG`
command did, calling straight into `tridiagonal`/`eigen_symmetric` or
`hessenberg`/`eigen_general` beneath it. `api`'s own implementation is spread across a
few files kept under this project's line-count guideline - `eig.zig`, `divide.zig`,
`svdapi.zig`, `matrixfun.zig` - all re-exported from `api.zig` so that every name still
resolves as `api.<name>` regardless of which file it actually lives in.

One deliberate naming departure from the original: MATLAB88 overloaded `SIN`, `COS`,
`EXP`, `SQRT`, `LOG`, and `ATAN` to mean *elementwise* application for a vector but the
*matrix function* (via eigendecomposition) for a square matrix, picking one or the
other based on the argument's shape. This port keeps the elementwise versions under
those names (`api.sin`, `api.cos`, ...) and gives the matrix-valued versions their own
names (`api.sinm`, `api.cosm`, `api.expm`, `api.logm`, `api.sqrtm`, `api.atanm`) -
matching the distinction real MATLAB itself later adopted (`expm` vs. `exp`), and
avoiding a shape-dependent meaning change that's easy to trip over in a typed language.

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

Most of this port doesn't allocate - functions read/write caller-supplied buffers,
sized as documented on each one. The exceptions, all in `api`, are `eig`, `cond`,
`rank`, `pinv`, `orth`, `roots`, and the `expm`/`logm`/`sqrtm`/`sinm`/`cosm`/`atanm`
family, which take a `std.mem.Allocator` because their scratch requirements depend on
a runtime code path (e.g. `eig`'s Hermitian-vs-general dispatch) or on an
eigendecomposition computed internally.

## What's *not* here

`matlab88.f90` is half REPL, half numerics. Left out entirely: the parser/tokenizer,
the `<>` expression stack machine, `FOR`/`WHILE`/`IF`/macros, `SAVE`/`LOAD`/`DIARY`
(the REPL's own binary variable-persistence format), `PRINT`/`DISP` (terminal
rendering - `rstd` has its own I/O for that), `PLOT` (line-printer/terminal plotting),
`EXEC` (script execution), `USER` (linking custom Fortran into the interpreter), and
the `HELP` text database. `CHAR` and `CHOP` are REPL input/display settings with no
meaning outside a session. `SIZE` is REPL introspection a typed caller already has
(the dimensions it passed in).

## Function reference

### `matlab.api` - friendly (REPL-named) layer

| Name | Signature (elided) | What it does |
| --- | --- | --- |
| `det` | `(ar, ai, lda, n, ipvt) Determinant` | Determinant via LU. `A` destroyed. |
| `inv` | `(ar, ai, lda, n, ipvt, zr, zi) f64` | In-place inverse via LU; returns `rcond`. |
| `chol` | `(ar, ai, lda, n) usize` | Cholesky factor into a clean upper-triangular `R`. |
| `rcond` | see `lu.wgeco` | LU-based reciprocal condition number. Alias of `lu.wgeco`. |
| `rref` | `(ar, ai, lda, m, n, eps) void` | Reduced row-echelon form in place. Alias of `generators.rref`. |
| `hilb` | `(a, lda, n) void` | Exact inverse Hilbert matrix (real only). Alias of `generators.hilber`. |
| `magic` | `(a, lda, n) void` | Magic square generator. Alias of `generators.magic`. |
| `eye` | `(ar, ai, lda, m, n) void` | Rectangular identity. Alias of `generators.eye`. |
| `ones` | `(ar, ai, lda, m, n) void` | All-ones matrix. Alias of `generators.ones`. |
| `tril` / `triu` | `(ar, ai, lda, m, n, k) void` | Zero the opposite triangle in place. Alias of `generators.tril`/`triu`. |
| `kron` | see `generators.kron` | Kronecker (tensor) product. Alias of `generators.kron`. |
| `diagFromVector` / `diagToVector` | see `generators.diag*` | `DIAG`'s two overloads: build a diagonal matrix from a vector, or extract one. |
| `base` | `(allocator, x, b, eps) ![]u8` | Format `x` in base `b` as an ASCII string. Alias of `base.base`. |
| `rat` | `(x, max_term, d) Fraction` | Continued-fraction rational approximation. Alias of `rat.rat`. |
| `lu` | `(ar, ai, lda, n, ipvt) usize` | Packed LU factorization. Alias of `lu.wgefa`. |
| `svd` | see `svd.wsvdc` | Alias of `svd.wsvdc`. |
| `qr` / `qrsl` | see `qr.wqrdc`/`wqrsl` | Aliases of `qr.wqrdc`/`qr.wqrsl`. |
| `mrdivide` | `(ar, ai, lda, m, n, br, bi, ldb, ipvt, zr, zi, rowr, rowi) f64` | Right matrix division `A/B` (square `B`) via LU. |
| `eig` | `(allocator, ar, ai, n, want_vectors) !EigResult` | Eigenvalues (+ optional eigenvectors); dispatches Hermitian vs. general. |
| `cond` | `(allocator, ar, ai, n, p) !f64` | 2-norm condition number via SVD. |
| `rank` | `(allocator, ar, ai, n, p, tol) !usize` | Numerical rank via SVD. |
| `pinv` | `(allocator, ar, ai, n, p, tol, result_r, result_i) !void` | Moore-Penrose pseudoinverse via SVD. |
| `orth` | `(allocator, ar, ai, lda, m, n) !struct{r,i}` | Orthonormal basis for the column space, via QR. |
| `poly` | `(rootsr, rootsi, n, cr, ci) void` | Characteristic polynomial from `n` roots. |
| `roots` | `(allocator, coeffr, coeffi, wr, wi) !usize` | Roots of a polynomial, via its companion matrix's eigenvalues. |
| `sum` / `prod` | `(ar, ai) Cx` | Sum/product of all entries (any shape). |
| `sin`/`cos`/`exp`/`sqrt`/`log`/`atan` | `(ar, ai) void` | Elementwise, in place (any shape). |
| `sinm`/`cosm`/`expm`/`logm`/`sqrtm`/`atanm` | `(allocator, ar, ai, n, out_r, out_i) !void` | The *matrix* function, via eigendecomposition (`n`x`n` only). |

`EigResult` has `values_r`/`values_i` (length `n`) and optional `vectors_r`/`vectors_i`
(`n`x`n`, column-major); call `.deinit(allocator)` to free it. `ABS`/`ROUND`/`REAL`/
`IMAG`/`CONJG` weren't given wrapper functions - they're one-liners you can inline
directly: `blas1.pythag(re, im)` (ABS), `elementary.round(re)` (ROUND, real part only -
the original ignored the imaginary part), `re`/`im` themselves (REAL/IMAG), and
`Cx.init(re, -im)` (CONJG).

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
| `elementary` | `wsin`, `wcos`, `wexp`, `wpow` | `MATFN2`'s `SIN`/`COS`/`EXP`/power bodies | Elementary complex scalar functions. |
| `elementary` | `round` | `ROUND` | Round to the nearest integer via exact binary decomposition. |
| `lu` | `wgefa` | `WGEFA` | Complex LU factorization with partial pivoting. |
| `lu` | `wgesl` | `WGESL` | Solve `A*x=b` or `ctrans(A)*x=b` from `wgefa`'s factors. |
| `lu` | `wgeco` | `WGECO` | `wgefa` + reciprocal condition number estimate. |
| `lu` | `wgedi` | `WGEDI` | Determinant and/or inverse from `wgefa`/`wgeco`'s factors. |
| `cholesky` | `wpofa` | `WPOFA` | Complex Cholesky factorization. |
| `generators` | `rref` | `RREF` | Reduced row-echelon form. |
| `generators` | `hilber` | `HILBER` | Exact inverse Hilbert matrix (closed form). |
| `generators` | `eye`, `ones` | `MATFN6`'s `EYE`/`ONES` bodies | Identity / all-ones matrix generators. |
| `generators` | `magic` | `MAGIC` | Magic square (odd/singly-even/doubly-even). |
| `generators` | `diagFromVector`, `diagToVector` | `MATFN6`'s `DIAG` body | `DIAG`'s vector-to-matrix and matrix-to-vector forms. |
| `generators` | `tril`, `triu` | `MATFN6`'s `TRIL`/`TRIU` bodies | Zero the opposite triangle. |
| `generators` | `kron` | `MATFN6`'s `KRON` body | Kronecker product (and its `./`/`.\` elementwise-block variants). |
| `tridiagonal` | `htridi` | `HTRIDI` | Reduce a Hermitian matrix to real symmetric tridiagonal form. |
| `tridiagonal` | `htribk` | `HTRIBK` | Back-transform tridiagonal eigenvectors to the original Hermitian matrix's. |
| `eigen_symmetric` | `imtql2` | `IMTQL2` | Eigenvalues/eigenvectors of a real symmetric tridiagonal matrix. |
| `hessenberg` | `corth` | `CORTH` | Reduce a general complex matrix to upper Hessenberg form. |
| `eigen_general` | `comqr3` | `COMQR3` | Eigenvalues/Schur form/eigenvectors of a complex Hessenberg matrix (shifted QR). |
| `svd` | `wsvdc` | `WSVDC` | Singular value decomposition of a general complex matrix. |
| `qr` | `wqrdc` | `WQRDC` | QR decomposition with optional column pivoting. |
| `qr` | `wqrsl` | `WQRSL` | Apply `wqrdc`'s output: `Q*y`, `ctrans(Q)*y`, least squares, residual, projection. |
| `rat` | `rat` | `RAT` | Continued-fraction rational approximation. |
| `base` | `base` | `BASE` | Format a number in an arbitrary base (ASCII output; see the module doc comment). |
| `random` | `Urand` | `URAND` + `MATFN6`'s `RAND` body | Legacy 1970s LCG (`.next()`) plus a Marsaglia-polar normal draw (`.nextNormal()`) and matrix fills. Kept for REPL-command parity and reproducible sequences, not statistical quality - reach for `std.Random` elsewhere. |

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

## Example: a random matrix, and the matrix exponential

```zig
const std = @import("std");
const matlab = @import("matlab/__root__.zig");

pub fn randomAndExpm(allocator: std.mem.Allocator) !void {
    var rng = matlab.random.Urand.init(42);
    var ar: [9]f64 = undefined;
    var ai: [9]f64 = @splat(0.0);
    rng.fillUniform(&ar, 3, 3); // ar now holds a 3x3 uniform-random real matrix

    // expm needs a diagonalizable matrix; a random one usually is.
    var out_r: [9]f64 = undefined;
    var out_i: [9]f64 = undefined;
    try matlab.api.expm(allocator, &ar, &ai, 3, &out_r, &out_i);
}
```
