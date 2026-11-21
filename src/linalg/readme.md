# Linear algebra

`rstd.linalg` is dense linear algebra for Zig: vectors and matrices of any float or complex element type, the operations MATLAB users expect under the names they expect, and results that are bit-identical for any thread count. Below it sit [`linalg.blas`](../math/kernels/blas/readme.md) and `linalg.lapack`, for tools, tests and experts; most code needs only this layer. How it is built, and why, is in [`design.md`](design.md).

## Matrices and views

| Type | What it is |
| --- | --- |
| `Mat(T)` | a matrix on the heap, column-major; freed with `deinit(allocator)` |
| `View(T)`, `ConstView(T)` | a window on any column-major storage: a slice, its rows and columns, and the distance between columns |
| `StaticMat(r, c, T)` | a fixed-size matrix by value, for small shapes |
| `Vec(n, T)` | a fixed-size vector |

`T` is `f16`, `f32`, `f64`, `f80`, `f128` or `std.math.Complex` of one. Every function takes its operands as any of a `Mat`, a view, a pointer to a `StaticMat`, or a slice (a column vector), in any mix; a vector is a view of one column or one row. `view.sub(i, j, rows, cols)` takes a block without copying.

## Conventions

- A function that allocates takes the allocator first, then an `std.Io`: `(allocator, io, ...)`. A result matrix comes back as an owned `Mat(T)`; most have an `...Into` form that writes into storage you pass.
- Inputs are never changed.
- Errors are exact and plainly named: `DimensionMismatch`, `NotSquare`, `Singular`, `NotPositiveDefinite`, `NoConvergence`, plus `OutOfMemory` where a function allocates and `Canceled` from `io`. Nothing panics on bad input.
- Up to `inline_limit` (8 by default, `Options`) a shape takes the fixed-size code path, so a `Mat` and a `StaticMat` holding the same numbers give the same bits; past it, the BLAS and LAPACK.
- Results are the same bits on every run and for any number of threads.

```zig
const linalg = rstd.linalg;

// x := A \ b, then A x again
const x = try linalg.solve(allocator, io, a, b);
defer x.deinit(allocator);
const r = try linalg.mul(allocator, io, a, x);
defer r.deinit(allocator);
```

## Functions

| Family | Functions |
| --- | --- |
| Vectors | `dot`, `norm`, `vectorNorm` (1, 2, inf or `p`), `sum`, `max`, `min`, `axpy`, `scale` |
| Products | `mul`, `transpose`, `ctranspose`, and their `...Into` forms |
| Element-wise | `add`, `sub`, `times` (`.*`), `rdivide` (`./`), `prod`, `map`, `realPart`, `imagPart`, `conj` |
| Building | `identity`, `ones`, `tril`, `triu`, `diagMatrix`, `diagonal`, `kron` |
| Solving | `solve` (square, triangular, Hermitian or least squares, chosen from the values), `inv`, `det`, `rcond`, `rref` |
| Factorizations | `lu`, `chol`, `qr`, `svd`, `eig`, `hess`, `schur`, `schurComplex` |
| From the SVD | `pinv`, `rank`, `cond`, `orth`, `nullSpace`, `matrixNorm` (1, inf, Frobenius, 2) |
| Matrix functions | `expm`, `logm`, `sqrtm`, `powm`, `funm`, `sinm`, `cosm`, `atanm` |
| Polynomials | `poly.roots`, `poly.poly`, `poly.charPoly` |
| Test matrices | `gallery.magic`, `gallery.hilbert`, `gallery.invHilbert`, `gallery.uniform`, `gallery.normal` |

A real matrix stays real wherever the result is real: `schur` gives the real Schur form, and `sqrtm`, `logm` and `powm` work on it, returning `NoRealSquareRoot` or `NoRealLogarithm` when the answer is complex (pass the matrix as complex then).

## Below this layer

- [`linalg.blas`](../math/kernels/blas/readme.md): the BLAS, OpenBLAS ported to Zig with its per-core kernels, on slices and strided views.
- `linalg.lapack`: the LAPACK routines this layer uses, ported from faer, Reference-LAPACK and ReLAPACK, on slices.
- `linalg.matrix` and `linalg.vector`: the fixed-size types and their operations.
