# matlab88 linalg port

A Zig port of `matlab88.f90`'s function library - Cleve Moler's Fortran MATLAB, written mostly at the University of New Mexico on LINPACK and EISPACK routines (`ref/matlab88_doc.txt`: the MATLAB Users' Guide of May 1981, its acknowledgement and references [1]-[3]; the source prints "Version of 05/25/82", `ref/matlab88.f90:2902`). `ref/readme.txt` is a later bulletin board note; its "developed by Argonne National Laboratories" is the distributor's, not Moler's. See those files and `ref/matlab88_help.txt` for the original documentation. The goal is every function a MATLAB88 user could call, minus the REPL itself: the parser, tokenizer, `<>` stack machine, `FOR`/`WHILE`/`IF`, `SAVE`/`LOAD`/`DIARY`/`PLOT`/`PRINT` file and terminal I/O, and `HELP` text are all left out, since none of that has meaning outside an interactive session.

## Importing

From elsewhere in `rstd`, the whole port is reachable as `rstd.linalg.matlab`:

```zig
const rstd = @import("rstd");
const matlab = rstd.linalg.matlab;
```

Within `linalg/` itself (or when depending on this directory directly), import `matlab/__root__.zig`:

```zig
const matlab = @import("matlab/__root__.zig");
```

`matlab.Cx` is `std.math.Complex(f64)` - the type used for individual complex scalars throughout this port (see "Data layout" below for how matrices and vectors are represented).

## Two surfaces

- **`matlab.CxMat`** - bundles a complex matrix's storage (`r`/`i`), leading dimension (`lda`), and shape (`m`/`n`) into one value, and carries every REPL command that operates on one matrix as its primary subject as a method: `mat.det(allocator, io)`, `mat.eye()`, `mat.chol(allocator, io, null)`, `CxMat.eig(allocator, io, mat, true)`, and so on. Use this for anything self-shaped. See its own doc comment for `lda`/striding details and `init`/`initPacked`/`clone`.
- **`matlab.<name>`** (top-level) - the REPL commands that *aren't* self-shaped: a real-only generator with no `CxMat` shape to bundle (`hilb`, `magic`), a flat-array/vector operation (`sum`, `prod`, `poly`, `roots`, vector `norm`, the elementwise `sin`/`cos`/... family), or a standalone command (`base`, `rat`, `Rand`). `roots` takes a `std.mem.Allocator` because it needs scratch space sized at runtime; the rest work directly on caller-supplied buffers, same as `CxMat`'s methods.

Everything else in this directory (`lu`, `cholesky`, `svd`, `scalar`, `generators`, ...) - the low-level kernels, named after their original Fortran routines (lowercased; Fortran is case-insensitive, so this renames nothing about the source) - is an internal implementation detail behind those two surfaces, not part of the public API. Every matrix-level (not vector-level BLAS1) kernel among these takes its primary matrix as a `CxMat` first parameter, so `CxMat`'s methods are exactly these calls with dot-call sugar, not a second implementation - `mat.det(allocator, io)` *is* `lu.det(mat, allocator, io)`. See "Internal kernels" below if you need to read the source for one of these.

One deliberate naming departure from the original: MATLAB88 overloaded `SIN`, `COS`, `EXP`, `SQRT`, `LOG`, and `ATAN` to mean *elementwise* application for a vector but the *matrix function* (via eigendecomposition) for a square matrix, picking one or the other based on the argument's shape. This port keeps the elementwise versions under those top-level names (`matlab.sin`, `matlab.cos`, ...) and gives the matrix-valued versions their own names on `CxMat` (`CxMat.sinm(allocator, io, a, out)`, `cosm`, `expm`, `logm`, `sqrtm`, `atanm`) - matching the distinction real MATLAB itself later adopted (`expm` vs. `exp`, <https://www.mathworks.com/help/matlab/ref/expm.html>), and avoiding a shape-dependent meaning change that's easy to trip over in a typed language.

## Operators (no operator overloading in Zig)

Every REPL operator with matrix/vector meaning has a named `CxMat` method or top-level function standing in for it, since Zig has no `+`/`*`/`'` overloading:

| REPL operator | Meaning | Stands in as |
| --- | --- | --- |
| `X'` | Conjugate transpose | `mat.ctranspose(out)` (plain `mat.transpose(out)` also available) |
| `X + Y` | Addition | `mat.add(other, out)` |
| `X - Y` | Subtraction | `mat.sub(other, out)` |
| `X*Y` | Matrix multiplication | `mat.matmul(other, out, allocator, io)` |
| `s*X` | Scalar multiply | `mat.scale(io, s)` (in place) |
| `A\B` | Left division | `mat.mldivide(b, allocator, io)` |
| `B/A` | Right division | `mat.mrdivide(other, allocator, io)` |
| `X .* Y` | Entrywise multiply | `mat.dotMul(other, out)` |
| `X ./ Y` | Entrywise divide | `mat.dotDiv(other, out)` |
| `X .\ Y` | Entrywise left divide | `mat.dotRdiv(other, out)` |
| `X .*. Y` | Kronecker product | `mat.kron(other, out, .mul)` |
| `X ./. Y` | Kronecker quotient | `mat.kron(other, out, .div)` |
| `X .\. Y` | Kronecker left quotient | `mat.kron(other, out, .rdiv)` |
| `X**p` | Power (matrix `X`, scalar `p`) | `CxMat.powm(allocator, io, a, p, out)` |
| `x**y` | Power (scalar `x`) | `matlab.pow(x, y)` |

`ABS`/`ROUND`/`REAL`/`IMAG`/`CONJG` are trivial enough to stay as inline one-liners (see the function reference below) rather than named operator stand-ins.

## Data layout

Complex matrices and vectors are **split real/imaginary arrays** (`ar: []f64, ai: []f64`), the original Fortran's `AR`/`AI`, not an interleaved `Complex(f64)` slice. The numerics are rstd's BLAS and LAPACK (M2), which convert at the boundary (`CxMat.isReal`, `toPacked`, `fromPacked`): a real matrix (`ai` all zero) goes to the real routines, in place where it can; a complex one through one packed interleaved copy. A single complex *scalar* (a multiplier, a temporary, a return value) is `Cx = std.math.Complex(f64)`.

Matrices are **column-major**, with an explicit leading dimension (`lda`), matching Fortran array storage - `A(i, j)` in the original is `ar[scalar.ix(lda, i, j)]` here. Row/column indices passed to these functions are **1-based**, again matching the Fortran source; `scalar.ix(lda, i, j)` converts a 1-based `(i, j)` to the 0-based flat offset. `lda` is the storage stride between columns - normally equal to the matrix's row count `m` (`CxMat.initPacked` assumes this), but it may exceed `m`, letting a `CxMat` describe a submatrix view into a larger allocated block without copying.

`CxMat` (`r`/`i`/`lda`/`m`/`n`) bundles exactly these five values so matrix-level functions take one parameter instead of five; `CxMat.get`/`.set`/`.col`/`.rowView` give indexed/strided access, and `.idx` is `scalar.ix(self.lda, ...)`.

Vector arguments (`vecNorm2`, ...) take an explicit stride (`incx`, may be negative) as BLAS does: the slice starts at the storage's first element, and a negative stride walks it backward.

Most of this port doesn't allocate - functions read/write caller-supplied buffers, sized as documented on each one. The exceptions are `CxMat`'s `det`, `inv`, `mldivide`, `mrdivide`, `chol`, `rref`, `matmul`, `matNorm2`, `eig`, `hess`, `schur`, `cond`, `rank`, `pinv`, `orth`, and the `expm`/`logm`/`sqrtm`/`sinm`/`cosm`/`atanm` family, plus the top-level `roots` - all take a `std.mem.Allocator` because their scratch requirements depend on a runtime code path (e.g. `eig`'s Hermitian-vs-general dispatch), on an eigendecomposition/QR factorization computed internally, or on LAPACK (M2: a complex matrix goes through an interleaved copy); those over LAPACK also take a `std.Io`.

## What's *not* here

`matlab88.f90` is half REPL, half numerics. Left out entirely: the parser/tokenizer, the `<>` expression stack machine, `FOR`/`WHILE`/`IF`/macros, `SAVE`/`LOAD`/`DIARY` (the REPL's own binary variable-persistence format), `PRINT`/`DISP` (terminal rendering - `rstd` has its own I/O for that), `PLOT` (line-printer/terminal plotting), `EXEC` (script execution), `USER` (linking custom Fortran into the interpreter), and the `HELP` text database. `CHAR` and `CHOP` are REPL input/display settings with no meaning outside a session. `SIZE` is REPL introspection a typed caller already has (the dimensions it passed in).

## Function reference

### `matlab.CxMat` - self-shaped matrix methods

| Name | Signature (elided) | What it does |
| --- | --- | --- |
| `det` | `(self, allocator, io) !Determinant` | Determinant via `getrf`, as a mantissa times `10^exponent10`. `self` destroyed. |
| `inv` | `(self, allocator, io) !f64` | In-place inverse via `getri`; returns `gecon`'s `rcond`; `error.Singular` for a zero pivot. |
| `chol` | `(self, allocator, io, minor) !void` | Cholesky factor into a clean upper-triangular `R`, by `potrf(.u)`. |
| `rref` | `(self, allocator, io, eps) !void` | Reduced row-echelon form in place, by BLAS 1. |
| `eye` | `(self) void` | Rectangular identity. |
| `ones` | `(self) void` | All-ones matrix. |
| `tril` / `triu` | `(self, k) void` | Zero the opposite triangle in place. |
| `kron` | `(a, b, out, op) void` | Kronecker (tensor) product. |
| `diagFromVector` | `(self, vr, vi, n, k) void` | `DIAG`: build a diagonal matrix from a vector. |
| `diagToVector` | `(self, k, dr, di) usize` | `DIAG`: extract a diagonal into a vector. |
| `mrdivide` | `(a, b, allocator, io) !f64` | Right matrix division `A/B` (square `B`) via `getrs`; returns `B`'s `rcond`. |
| `mldivide` | `(self, b, allocator, io) !f64` | Left matrix division `A\B` (square `A`) via `getrs`; returns `A`'s `rcond`. |
| `matmul` | `(a, b, out, allocator, io) !void` | `C = A*B`, by `gemm`. |
| `ctranspose` | `(self, out) void` | `X'`: conjugate transpose. |
| `transpose` | `(self, out) void` | Plain (non-conjugating) transpose. |
| `scale` | `(self, io, s) !void` | Multiply every entry by a scalar in place (`s*X`/`X*s`), by `drot`. |
| `add` / `sub` | `(a, b, out) void` | `X + Y` / `X - Y`, entrywise. |
| `dotMul` / `dotDiv` / `dotRdiv` | `(a, b, out) void` | `.*` / `./` / `.\`, entrywise (not Kronecker - see `kron` above). |
| `powm` | `(allocator, io, a, p, out) !void` | `X**p` (matrix case): via eigendecomposition, like `expm`/... |
| `eig` | `(allocator, io, a, want_vectors) !EigResult` | Eigenvalues (+ optional eigenvectors); dispatches Hermitian (`syevd`) vs. general (`geev`). |
| `hess` | `(allocator, io, a) !HessResult` | Unitary similarity `P` and Hessenberg form `H` (`A = P*H*P'`), by `gehrd` and `orghr`. |
| `schur` | `(allocator, io, a) !SchurResult` | Unitary similarity `U` and triangular Schur form `T` (`A = U*T*U'`), by `zgees`. |
| `sinm`/`cosm`/`expm`/`logm`/`sqrtm`/`atanm` | `(allocator, io, a, out) !void` | The *matrix* function, via eigendecomposition. |
| `cond` | `(allocator, io, a) !f64` | 2-norm condition number via `gesvd`. |
| `rank` | `(allocator, io, a, tol) !usize` | Numerical rank via `gesvd`. |
| `pinv` | `(allocator, io, a, tol, result_r, result_i) !void` | Moore-Penrose pseudoinverse via `gesvd`. |
| `orth` | `(allocator, io, a) !struct{r,i}` | Orthonormal basis for the column space, via `geqrf` and `orgqr`. |
| `matNorm2` | `(self, allocator, io) !f64` | Largest singular value, via `gesvd`. |
| `matNorm1` / `matNormInf` / `matNormFro` | `(self, io) !f64` | 1-norm / infinity-norm / Frobenius norm, by `dasum` and `dnrm2`. |

`EigResult`/`HessResult`/`SchurResult` own caller-freed slices - call `.deinit(allocator)`.

### `matlab.<name>` (top-level) - REPL commands that aren't `CxMat`-shaped

| Name | Signature (elided) | What it does |
| --- | --- | --- |
| `hilb` | `(a, lda, n) void` | Exact inverse Hilbert matrix (real only - no `CxMat`, see its own doc comment). Alias of `generators.hilber`. |
| `magic` | `(io, a, lda, n) !void` | Magic square generator (real only). Alias of `generators.magic`. |
| `base` | `(allocator, x, b, eps) ![]u8` | Format `x` in base `b` as an ASCII string. Alias of `base.base`. |
| `rat` | `(x, max_term, d) Fraction` | Continued-fraction rational approximation. Alias of `rat.rat`. |
| `poly` | `(rootsr, rootsi, n, cr, ci) void` | Characteristic polynomial from `n` roots. |
| `roots` | `(allocator, io, coeffr, coeffi, wr, wi) !usize` | Roots of a polynomial, via its companion matrix's eigenvalues (`geev`, balanced). |
| `sum` / `prod` | `(ar, ai) Cx` | Sum/product of all entries (any shape). |
| `sin`/`cos`/`exp`/`sqrt`/`log`/`atan` | `(ar, ai) void` | Elementwise, in place (any shape). |
| `pow` | see `elementary.wpow` | `x^y` (principal branch). Alias of `elementary.wpow`. |
| `Rand` | see `random.Urand` | Legacy LCG RNG. Alias of `random.Urand`. |
| `vecNorm2`/`vecNormP`/`vecNormInf` | see `norm.zig` | Vector norms, `vecNorm2` by `dnrm2` (the matrix-norm forms are `CxMat` methods - see above). |

`ABS`/`ROUND`/`REAL`/`IMAG`/`CONJG` weren't given wrapper functions - they're one-liners you can inline directly: `std.math.hypot(re, im)` (ABS), `@round(re)` (ROUND, real part only - the original ignored the imaginary part), `re`/`im` themselves (REAL/IMAG), and `Cx.init(re, -im)` (CONJG). Likewise, `SIN`/`COS`/`EXP`/`SQRT` route straight through `std.math.complex` (via `matrixfun.zig`'s `applyElementwise`/`applyMatrixFunction`) rather than through a named wrapper in this port - std's versions carry overflow/edge-case handling the original Fortran formulas didn't have.

### Internal kernels, by module (not part of the public API)

These modules aren't reachable as `matlab.<module>` - they're implementation details behind `CxMat`'s methods and the top-level functions above (see "Two surfaces"). Listed here for anyone reading the source or wanting the exact original Fortran routine name.

| Module | Function | Ported from | What it does |
| --- | --- | --- | --- |
| `scalar` | `wdiv` | `WDIV` | Complex scalar divide, overflow-scaled. `PYTHAG` is `std.math.hypot` and `WMUL` is `Cx.mul`, no wrappers. |
| `scalar` | `wsign`, `wlog`, `watan` | `WSIGN`, `WLOG`, `WATAN` | Complex scalar sign-copy, log, arctangent. `WSQRT` is `std.math.complex.sqrt` directly - no wrapper. |
| `scalar` | `level1` | `WNRM2`, `WASUM`, `IWAMAX`, `WSCAL`, `WAXPY`, `WSWAP`, `RSWAP` | rstd's BLAS level 1 since M2; `level1` drops its stride errors. |
| `scalar` | `ix` | *(new)* | 1-based `(i, j)` → 0-based flat column-major offset. |
| `elementary` | `wpow` | `MATFN2`'s power body | `x^y` (principal branch). `SIN`/`COS`/`EXP`/`SQRT` are `std.math.complex` calls directly - no wrapper. `ROUND` is `@round` directly. |
| `rat` | `rat` | `RAT` | Continued-fraction rational approximation. |
| `base` | `base` | `BASE` | Format a number in an arbitrary base (ASCII output; see the module doc comment). |
| `lu` | `Lu`, `lange1`, `rcond` | `WGEFA`, `WGESL`, `WGECO` | LAPACK's `getrf`, `getrs`, `lange` and `gecon` at the split-layout boundary (M2). |
| `lu` | `det`, `inv` | `WGEDI` | `WGEDI`'s scaled determinant over `getrf`'s factors; the inverse by `getri`. |
| `divide` | `mldivide`, `mrdivide` | `WGECO`, `WGESL` | `A\B` and `B/A` by `getrf` and `getrs`. |
| `cholesky` | `chol` | `WPOFA` | LAPACK's `potrf(.u)` at the split-layout boundary (M2). |
| `generators` | `hilber` | `HILBER` | Exact inverse Hilbert matrix (closed form). |
| `generators` | `eye`, `ones` | `MATFN6`'s `EYE`/`ONES` bodies | Identity / all-ones matrix generators. |
| `generators` | `magic` | `MAGIC` | Magic square (odd/singly-even/doubly-even). |
| `generators` | `diagFromVector`, `diagToVector` | `MATFN6`'s `DIAG` body | `DIAG`'s vector-to-matrix and matrix-to-vector forms. |
| `generators` | `tril`, `triu` | `MATFN6`'s `TRIL`/`TRIU` bodies | Zero the opposite triangle. |
| `generators` | `rref` | `RREF` | Reduced row-echelon form, by BLAS 1. |
| `generators` | `kron` | `MATFN6`'s `KRON` body | Kronecker product (and its `./`/`.\` elementwise-block variants). |
| `matmul` | `matmul` | `STACK2`'s `*` | `C = A*B` by `gemm` (M2). |
| `elementwise` | `add`, `sub`, `dotMul`, `dotDiv`, `dotRdiv` | `STACK2`'s entrywise operators | `+`, `-`, `.*`, `./`, `.\`. |
| `__root__` | `sum`, `prod` | `MATFN6`'s `SUM`/`PROD` bodies | Sum/product of all entries. |
| `eig` | `eig` | `HTRIDI`, `IMTQL2`, `HTRIBK`, `CORTH`, `COMQR3` | LAPACK's `syevd` and `geev` (M2). |
| `schur` | `hess`, `schur` | `MATFN4`/`MATFN5`'s `HESS`/`SCHUR` bodies, `CORTH`, `COMQR3` | LAPACK's `gehrd`, `orghr` and `gees` (M2). |
| `svd` | `singularValues` | `WSVDC`, `WQRDC`, `WQRSL` | LAPACK's `gesvd`, `geqrf` and `orgqr` behind `cond`, `rank`, `pinv`, `orth` and `matNorm2` (M2). |
| `matrixfun` | `poly`, `roots`, `sin` ... `atan`, `expm` ... `powm` | `MATFN2` | Polynomials, and the elementwise and matrix elementary functions. |
| `norm` | `vecNorm2` ... `matNormFro` | `MATFN3`'s `NORM` body | Vector and matrix norms. |
| `random` | `Urand` | `URAND` + `MATFN6`'s `RAND` body | Legacy LCG after Knuth (1969), as `URAND`'s comment says (`.next()`), plus `RAND`'s polar (rejection on the unit disk) normal draw (`.nextNormal()`) and matrix fills. Kept for REPL-command parity and reproducible sequences, not statistical quality - reach for `std.Random` elsewhere. |

See each function's doc comment for the full parameter contract (buffer sizes, `job` flag encodings, what's destroyed vs. read-only) - the table above is a map of what exists, not a substitute for those.

## Example: eigenvalues and eigenvectors

```zig
const std = @import("std");
const matlab = @import("matlab/__root__.zig");
const CxMat = matlab.CxMat;

pub fn eigExample(allocator: std.mem.Allocator, io: std.Io) !void {
    // A = [[2, 1], [1, 2]] (Hermitian; eig will pick the symmetric path).
    var ar = [_]f64{ 2.0, 1.0, 1.0, 2.0 }; // column-major
    var ai = [_]f64{ 0.0, 0.0, 0.0, 0.0 };
    const a = CxMat.initPacked(&ar, &ai, 2, 2);

    const result = try CxMat.eig(allocator, io, a, true);
    defer result.deinit(allocator);

    // result.values_r/values_i: eigenvalues (1 and 3, here).
    // result.vectors_r/vectors_i: matching eigenvectors as columns, if requested.
}
```

## Example: solve a linear system via LU factorization

```zig
const std = @import("std");
const matlab = @import("matlab/__root__.zig");
const CxMat = matlab.CxMat;

pub fn solveExample(allocator: std.mem.Allocator, io: std.Io) !void {
    // A = [[1+i, 2], [3, 4-i]], solve A*x = b in place.
    var ar = [_]f64{ 1.0, 3.0, 2.0, 4.0 }; // column-major
    var ai = [_]f64{ 1.0, 0.0, 0.0, -1.0 };
    const a = CxMat.initPacked(&ar, &ai, 2, 2);

    var br = [_]f64{ 5.0, 6.0 };
    var bi = [_]f64{ 0.0, 0.0 };
    _ = try a.mldivide(CxMat.initPacked(&br, &bi, 2, 1), allocator, io);
    // br/bi now hold the solution x.
}
```

## Example: a random matrix, and the matrix exponential

```zig
const std = @import("std");
const matlab = @import("matlab/__root__.zig");
const CxMat = matlab.CxMat;

pub fn randomAndExpm(allocator: std.mem.Allocator, io: std.Io) !void {
    var rng = matlab.Rand.init(42);
    var ar: [9]f64 = undefined;
    var ai: [9]f64 = @splat(0.0);
    rng.fillUniform(&ar, 3, 3); // ar now holds a 3x3 uniform-random real matrix
    const a = CxMat.initPacked(&ar, &ai, 3, 3);

    // expm needs a diagonalizable matrix; a random one usually is.
    var out_r: [9]f64 = undefined;
    var out_i: [9]f64 = undefined;
    try CxMat.expm(allocator, io, a, CxMat.initPacked(&out_r, &out_i, 3, 3));
}
```
