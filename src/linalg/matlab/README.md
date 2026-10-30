# matlab88 linalg port

A Zig port of `matlab88.f90`'s function library - Cleve Moler's Fortran MATLAB, written mostly at
the University of New Mexico on LINPACK and EISPACK routines (`ref/matlab88_doc.txt`: the
MATLAB Users' Guide of May 1981, its acknowledgement and references [1]-[3]; the source
prints "Version of 05/25/82", `ref/matlab88.f90:2902`). `ref/readme.txt` is a later bulletin
board note; its "developed by Argonne National Laboratories" is the distributor's, not
Moler's. See those files and `ref/matlab88_help.txt` for the original documentation. The goal is every function a MATLAB88
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

## Two surfaces

- **`matlab.CxMat`** - bundles a complex matrix's storage (`r`/`i`), leading dimension
  (`lda`), and shape (`m`/`n`) into one value, and carries every REPL command that
  operates on one matrix as its primary subject as a method: `mat.det(ipvt)`,
  `mat.eye()`, `mat.hess(allocator)`, `mat.eig(allocator, true)`, and so on. Use this
  for anything self-shaped. See its own doc comment for `lda`/striding details and
  `init`/`initPacked`/`clone`.
- **`matlab.<name>`** (top-level) - the REPL commands that *aren't* self-shaped: a
  real-only generator with no `CxMat` shape to bundle (`hilb`, `magic`), a
  flat-array/vector operation (`sum`, `prod`, `poly`, `roots`, vector `norm`, the
  elementwise `sin`/`cos`/... family), or a standalone command (`base`, `rat`, `Rand`).
  `roots` takes a `std.mem.Allocator` because it needs scratch space sized at runtime;
  the rest work directly on caller-supplied buffers, same as `CxMat`'s methods.

Everything else in this directory (`lu`, `svd`, `qr`, `blas1`, `generators`, ...) - the
low-level kernels, named after their original Fortran routines (lowercased; Fortran is
case-insensitive, so this renames nothing about the source) - is an internal
implementation detail behind those two surfaces, not part of the public API. Every
matrix-level (not vector-level BLAS1) kernel among these takes its primary matrix as a
`CxMat` first parameter, so `CxMat`'s methods are exactly these calls with dot-call
sugar, not a second implementation - `mat.wgefa(ipvt)` *is* `lu.wgefa(mat, ipvt)`. See
"Internal kernels" below if you need to read the source for one of these.

One deliberate naming departure from the original: MATLAB88 overloaded `SIN`, `COS`,
`EXP`, `SQRT`, `LOG`, and `ATAN` to mean *elementwise* application for a vector but the
*matrix function* (via eigendecomposition) for a square matrix, picking one or the
other based on the argument's shape. This port keeps the elementwise versions under
those top-level names (`matlab.sin`, `matlab.cos`, ...) and gives the matrix-valued
versions their own names on `CxMat` (`mat.sinm(allocator, out)`, `mat.cosm(...)`,
`mat.expm(...)`, `mat.logm(...)`, `mat.sqrtm(...)`, `mat.atanm(...)`) - matching the
distinction real MATLAB itself later adopted (`expm` vs. `exp`,
<https://www.mathworks.com/help/matlab/ref/expm.html>), and avoiding a
shape-dependent meaning change that's easy to trip over in a typed language.

## Operators (no operator overloading in Zig)

Every REPL operator with matrix/vector meaning has a named `CxMat` method or
top-level function standing in for it, since Zig has no `+`/`*`/`'` overloading:

| REPL operator | Meaning | Stands in as |
| --- | --- | --- |
| `X'` | Conjugate transpose | `mat.ctranspose(out)` (plain `mat.transpose(out)` also available) |
| `X + Y` | Addition | `mat.add(other, out)` |
| `X - Y` | Subtraction | `mat.sub(other, out)` |
| `X*Y` | Matrix multiplication | `mat.matmul(other, out, scratch_r, scratch_i)` |
| `s*X` | Scalar multiply | `mat.scale(s)` (in place) |
| `A\B` | Left division | `mat.mldivide(b, ipvt, zr, zi)` |
| `B/A` | Right division | `mat.mrdivide(other, ipvt, zr, zi, rowr, rowi)` |
| `X .* Y` | Entrywise multiply | `mat.dotMul(other, out)` |
| `X ./ Y` | Entrywise divide | `mat.dotDiv(other, out)` |
| `X .\ Y` | Entrywise left divide | `mat.dotRdiv(other, out)` |
| `X .*. Y` | Kronecker product | `mat.kron(other, out, .mul)` |
| `X ./. Y` | Kronecker quotient | `mat.kron(other, out, .div)` |
| `X .\. Y` | Kronecker left quotient | `mat.kron(other, out, .rdiv)` |
| `X**p` | Power (matrix `X`, scalar `p`) | `mat.powm(allocator, p, out)` |
| `x**y` | Power (scalar `x`) | `matlab.pow(x, y)` |

`ABS`/`ROUND`/`REAL`/`IMAG`/`CONJG` are trivial enough to stay as inline one-liners
(see the function reference below) rather than named operator stand-ins.

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
offset. `lda` is the storage stride between columns - normally equal to the matrix's
row count `m` (`CxMat.initPacked` assumes this), but it may exceed `m`, letting a
`CxMat` describe a submatrix view into a larger allocated block without copying.

`CxMat` (`r`/`i`/`lda`/`m`/`n`) bundles exactly these five values so matrix-level
functions take one parameter instead of five; `CxMat.get`/`.set`/`.col`/`.rowView` give
indexed/strided access, and `.idx` is `blas1.ix(self.lda, ...)`. Vector-level BLAS1
kernels are unaffected - they take a 1D strided view (a row, column, or arbitrary
vector with an explicit `incx`), not a whole matrix.

Vector arguments to the BLAS1-style kernels (`waxpy`, `wdotcr`, `wnrm2`, ...) take an
explicit stride (`incx`/`incy`, may be negative) and expect the slice to already start
at the first accessed element - exactly like Fortran handing a subroutine the address
of `AR(k, j)`. Most call sites in this port use a stride of 1 (a matrix column) or
`lda` (a matrix row); non-unit strides fall back to a scalar loop instead of the
`@Vector` fast path.

Most of this port doesn't allocate - functions read/write caller-supplied buffers,
sized as documented on each one. The exceptions are `CxMat`'s `eig`, `hess`, `schur`,
`cond`, `rank`, `pinv`, `orth`, and the `expm`/`logm`/`sqrtm`/`sinm`/`cosm`/`atanm`
family, plus the top-level `roots` - all take a `std.mem.Allocator` because their
scratch requirements depend on a runtime code path (e.g. `eig`'s Hermitian-vs-general
dispatch) or on an eigendecomposition/QR factorization computed internally.

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

### `matlab.CxMat` - self-shaped matrix methods

| Name | Signature (elided) | What it does |
| --- | --- | --- |
| `det` | `(self, ipvt) Determinant` | Determinant via LU. `self` destroyed. |
| `inv` | `(self, ipvt, zr, zi) f64` | In-place inverse via LU; returns `rcond`. |
| `chol` | `(self) usize` | Cholesky factor into a clean upper-triangular `R`. |
| `wgeco` | `(self, ipvt, zr, zi) f64` | `wgefa` + reciprocal condition number estimate (`RCOND`). |
| `rref` | `(self, eps) void` | Reduced row-echelon form in place. |
| `eye` | `(self) void` | Rectangular identity. |
| `ones` | `(self) void` | All-ones matrix. |
| `tril` / `triu` | `(self, k) void` | Zero the opposite triangle in place. |
| `kron` | `(a, b, out, op) void` | Kronecker (tensor) product. |
| `diagFromVector` | `(self, vr, vi, n, k) void` | `DIAG`: build a diagonal matrix from a vector. |
| `diagToVector` | `(self, k, dr, di) usize` | `DIAG`: extract a diagonal into a vector. |
| `wgefa` | `(self, ipvt) usize` | Packed LU factorization. |
| `wgesl` | `(self, ipvt, br, bi, job) void` | Solve `A*x=b` or `ctrans(A)*x=b` from `wgefa`'s factors. |
| `wgedi` | `(self, ipvt, workr, worki, job, det_out) void` | Determinant and/or inverse from `wgefa`/`wgeco`'s factors. |
| `wpofa` | `(self) usize` | Complex Cholesky factorization (raw, lower triangle left untouched). |
| `mrdivide` | `(a, b, ipvt, zr, zi, rowr, rowi) f64` | Right matrix division `A/B` (square `B`) via LU. |
| `mldivide` | `(self, b, ipvt, zr, zi) f64` | Left matrix division `A\B` (square `A`) via LU. |
| `matmul` | `(a, b, out, scratch_r, scratch_i) void` | `C = A*B`. |
| `ctranspose` | `(self, out) void` | `X'`: conjugate transpose. |
| `transpose` | `(self, out) void` | Plain (non-conjugating) transpose. |
| `scale` | `(self, s) void` | Multiply every entry by a scalar in place (`s*X`/`X*s`). |
| `add` / `sub` | `(a, b, out) void` | `X + Y` / `X - Y`, entrywise. |
| `dotMul` / `dotDiv` / `dotRdiv` | `(a, b, out) void` | `.*` / `./` / `.\`, entrywise (not Kronecker - see `kron` above). |
| `powm` | `(self, allocator, p, out) !void` | `X**p` (matrix case): via eigendecomposition, like `expm`/... |
| `eig` | `(self, allocator, want_vectors) !EigResult` | Eigenvalues (+ optional eigenvectors); dispatches Hermitian vs. general. |
| `hess` | `(self, allocator) !HessResult` | Unitary similarity `P` and Hessenberg form `H` (`A = P*H*P'`). |
| `schur` | `(self, allocator) !SchurResult` | Unitary similarity `U` and Schur form `T` (`A = U*T*U'`). |
| `corth` | `(self, low, igh, ortr, orti) void` | Reduce a submatrix to upper Hessenberg form (raw kernel behind `hess`/`schur`/`eig`). |
| `htridi` / `htribk` | see `tridiagonal.zig` | Reduce a Hermitian matrix to tridiagonal form / back-transform its eigenvectors. |
| `comqr3` | `(self, low, igh, ortr, orti, wr, wi, z, job) usize` | Eigenvalues/Schur form/eigenvectors of a Hessenberg matrix (shifted QR). |
| `sinm`/`cosm`/`expm`/`logm`/`sqrtm`/`atanm` | `(self, allocator, out) !void` | The *matrix* function, via eigendecomposition. |
| `wsvdc` | `(self, sr, si, er, ei, u, v, workr, worki, job) usize` | Singular value decomposition. `self` destroyed. |
| `wqrdc` | `(self, qrauxr, qrauxi, jpvt, workr, worki, job) void` | QR decomposition with optional column pivoting. |
| `wqrsl` | `(self, k, qrauxr, qrauxi, yr, yi, qyr, qyi, qtyr, qtyi, br, bi, rsdr, rsdi, xbr, xbi, job) usize` | Apply `wqrdc`'s output: `Q*y`, `ctrans(Q)*y`, least squares, residual, projection. |
| `cond` | `(self, allocator) !f64` | 2-norm condition number via SVD. |
| `rank` | `(self, allocator, tol) !usize` | Numerical rank via SVD. |
| `pinv` | `(self, allocator, tol, result_r, result_i) !void` | Moore-Penrose pseudoinverse via SVD. |
| `orth` | `(self, allocator) !struct{r,i}` | Orthonormal basis for the column space, via QR. |
| `matNorm2` | `(self, sr, si, er, ei, workr, worki) usize` | Largest singular value. |
| `matNorm1` / `matNormInf` / `matNormFro` | `(self) f64` | 1-norm / infinity-norm / Frobenius norm. |

`EigResult`/`HessResult`/`SchurResult` own caller-freed slices - call `.deinit(allocator)`.

### `matlab.<name>` (top-level) - REPL commands that aren't `CxMat`-shaped

| Name | Signature (elided) | What it does |
| --- | --- | --- |
| `hilb` | `(a, lda, n) void` | Exact inverse Hilbert matrix (real only - no `CxMat`, see its own doc comment). Alias of `generators.hilber`. |
| `magic` | `(a, lda, n) void` | Magic square generator (real only). Alias of `generators.magic`. |
| `base` | `(allocator, x, b, eps) ![]u8` | Format `x` in base `b` as an ASCII string. Alias of `base.base`. |
| `rat` | `(x, max_term, d) Fraction` | Continued-fraction rational approximation. Alias of `rat.rat`. |
| `poly` | `(rootsr, rootsi, n, cr, ci) void` | Characteristic polynomial from `n` roots. |
| `roots` | `(allocator, coeffr, coeffi, wr, wi) !usize` | Roots of a polynomial, via its companion matrix's eigenvalues. |
| `sum` / `prod` | `(ar, ai) Cx` | Sum/product of all entries (any shape). |
| `sin`/`cos`/`exp`/`sqrt`/`log`/`atan` | `(ar, ai) void` | Elementwise, in place (any shape). |
| `pow` | see `elementary.wpow` | `x^y` (principal branch). Alias of `elementary.wpow`. |
| `Rand` | see `random.Urand` | Legacy LCG RNG. Alias of `random.Urand`. |
| `vecNorm2`/`vecNormP`/`vecNormInf` | see `norm.zig` | Vector norms (the matrix-norm forms are `CxMat` methods - see above). |

`ABS`/`ROUND`/`REAL`/`IMAG`/`CONJG` weren't given wrapper functions - they're
one-liners you can inline directly: `blas1.pythag(re, im)` (ABS), `@round(re)` (ROUND,
real part only - the original ignored the imaginary part), `re`/`im` themselves
(REAL/IMAG), and `Cx.init(re, -im)` (CONJG). Likewise, `SIN`/`COS`/`EXP`/`SQRT` route
straight through `std.math.complex` (via `matrixfun.zig`'s
`applyElementwise`/`applyMatrixFunction`) rather than through a named wrapper in this
port - std's versions carry overflow/edge-case handling the original Fortran formulas
didn't have.

### Internal kernels, by module (not part of the public API)

These modules aren't reachable as `matlab.<module>` - they're implementation details
behind `CxMat`'s methods and the top-level functions above (see "Two surfaces"). Listed
here for anyone reading the source or wanting the exact original Fortran routine name.

| Module | Function | Ported from | What it does |
| --- | --- | --- | --- |
| `blas1` | `pythag` | `PYTHAG` | Overflow-safe `sqrt(a^2+b^2)` (aliases `std.math.hypot`). |
| `blas1` | `wmul`, `wdiv` | `WMUL`, `WDIV` | Complex scalar multiply/divide (`wdiv` is overflow-scaled). |
| `blas1` | `wsign`, `wlog`, `watan` | `WSIGN`, `WLOG`, `WATAN` | Complex scalar sign-copy, log, arctangent. `WSQRT` is `std.math.complex.sqrt` directly - no wrapper. |
| `blas1` | `wnrm2`, `wasum`, `iwamax` | `WNRM2`, `WASUM`, `IWAMAX` | Complex vector 2-norm, 1-norm, max-magnitude index. |
| `blas1` | `wrscal`, `wscal`, `waxpy` | `WRSCAL`, `WSCAL`, `WAXPY` | Complex vector scale (real/complex scalar) and `y += s*x`. |
| `blas1` | `wdotur`, `wdotui`, `wdotcr`, `wdotci` | `WDOTUR`/`UI`/`CR`/`CI` | Complex dot products (unconjugated/conjugated, real/imag part). |
| `blas1` | `wcopy`, `wset`, `wswap` | `WCOPY`, `WSET`, `WSWAP` | Complex vector copy, fill, swap. |
| `blas1` | `rset`, `rswap`, `rrot`, `rrotg` | `RSET`, `RSWAP`, `RROT`, `RROTG` | Real vector fill/swap and Givens rotation apply/construct. |
| `blas1` | `ix` | *(new)* | 1-based `(i, j)` → 0-based flat column-major offset. |
| `elementary` | `wpow` | `MATFN2`'s power body | `x^y` (principal branch). `SIN`/`COS`/`EXP`/`SQRT` are `std.math.complex` calls directly - no wrapper. `ROUND` is `@round` directly. |
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
| `random` | `Urand` | `URAND` + `MATFN6`'s `RAND` body | Legacy LCG after Knuth (1969), as `URAND`'s comment says (`.next()`), plus `RAND`'s polar (rejection on the unit disk) normal draw (`.nextNormal()`) and matrix fills. Kept for REPL-command parity and reproducible sequences, not statistical quality - reach for `std.Random` elsewhere. |

See each function's doc comment for the full parameter contract (buffer sizes, `job`
flag encodings, what's destroyed vs. read-only) - the table above is a map of what
exists, not a substitute for those.

## Example: eigenvalues and eigenvectors

```zig
const std = @import("std");
const matlab = @import("matlab/__root__.zig");
const CxMat = matlab.CxMat;

pub fn eigExample(allocator: std.mem.Allocator) !void {
    // A = [[2, 1], [1, 2]] (Hermitian; .eig will pick the symmetric path).
    var ar = [_]f64{ 2.0, 1.0, 1.0, 2.0 }; // column-major
    var ai = [_]f64{ 0.0, 0.0, 0.0, 0.0 };
    const a = CxMat.initPacked(&ar, &ai, 2, 2);

    const result = try a.eig(allocator, true);
    defer result.deinit(allocator);

    // result.values_r/values_i: eigenvalues (1 and 3, here).
    // result.vectors_r/vectors_i: matching eigenvectors as columns, if requested.
}
```

## Example: solve a linear system via LU factorization

```zig
const matlab = @import("matlab/__root__.zig");
const Cx = matlab.Cx;
const CxMat = matlab.CxMat;

pub fn solveExample() void {
    // A = [[1+i, 2], [3, 4-i]], solve A*x = b in place.
    var ar = [_]f64{ 1.0, 3.0, 2.0, 4.0 }; // column-major
    var ai = [_]f64{ 1.0, 0.0, 0.0, -1.0 };
    var ipvt: [2]usize = undefined;
    const a = CxMat.initPacked(&ar, &ai, 2, 2);

    _ = a.wgefa(&ipvt);

    var br = [_]f64{ 5.0, 6.0 };
    var bi = [_]f64{ 0.0, 0.0 };
    a.wgesl(&ipvt, &br, &bi, 0);
    // br/bi now hold the solution x.
}
```

## Example: a random matrix, and the matrix exponential

```zig
const std = @import("std");
const matlab = @import("matlab/__root__.zig");
const CxMat = matlab.CxMat;

pub fn randomAndExpm(allocator: std.mem.Allocator) !void {
    var rng = matlab.Rand.init(42);
    var ar: [9]f64 = undefined;
    var ai: [9]f64 = @splat(0.0);
    rng.fillUniform(&ar, 3, 3); // ar now holds a 3x3 uniform-random real matrix
    const a = CxMat.initPacked(&ar, &ai, 3, 3);

    // expm needs a diagonalizable matrix; a random one usually is.
    var out_r: [9]f64 = undefined;
    var out_i: [9]f64 = undefined;
    try a.expm(allocator, CxMat.initPacked(&out_r, &out_i, 3, 3));
}
```
