# `rstd.linalg` public layer: design proposal (M3)

Status: proposal for review (2026-10-30); no code yet. It builds on [`blas/design.md`](blas/design.md),
Settled decisions ("`rstd.linalg` is the public API", "an anthology", "Bit guarantees are the
CPU's", "Shaders keep their own API") and [`blas/todo.md`](blas/todo.md) items M3, M4, S1 and S2.
Open questions for the user are at the end, each with a recommendation.

## Today

| Where | What |
| --- | --- |
| `linalg.vector` | `dot` (a fused `mulAdd` chain), `dot2`, `length` (`@sqrt(dot2)`, `hypot` fallback off SPIR-V), `distance`, `norm(v, fallback)` (a safe normalize), `cross2`, `cross3`; generic over any `@Vector` |
| `rstd.simd` | `mix`, `mulAdd` (`blas/fma.zig`'s policy), `vectorLen`; `zigvec` (`V2f32` to `V4f64`, `SimdVecF32`/`F64`, 14 splat and zero helpers); `reduce` (`reduce`, `array{Max,Min,Sum,Product,And,Or,Xor}`, `mapReduce`) |
| callers | bare `@Vector` types everywhere; no matrix type in use (shaderlib's `Mat2` to `Mat4` are unused `[N]VecN`; glyf's `transformPoint` writes its 2x2 rows by hand) |
| SPIR-V (0.18.0-dev.751) | `@sqrt`, `@mulAdd`, `@abs`, `@floor`, `@min`/`@max` lower to GLSL.std.450 (`Sqrt`, `Fma`, `FAbs`, `Floor`, `FMin`/`FMax`); `@Vector` is `OpTypeVector` for lengths 2 to 4 (8 and 16 with `Vector16`), otherwise unrolled; `OpDot` and `OpExtInst` are reachable only through inline assembly; there is no matrix type: `@SpirvType` has none and the assembler's `OpTypeMatrix` is a TODO (`src/codegen/spirv/Assembler.zig`, `processTypeInstruction`) |

## 1. Types

| Type | Definition | Notes |
| --- | --- | --- |
| `Vec(N, T)` | `@Vector(N, T)` | the builtin itself, no wrapper: `+ - * /`, `@splat`, `@reduce` and indexing keep working for every caller |
| `Mat(R, C, T)` | `struct { cols: [C]Vec(R, T) }` | `R` rows, `C` columns, column-major: `m.cols[j]` is column `j` |
| aliases | `V2f32` to `V4f64` (from `zigvec`); `M2f32` to `M4f64` (square) | |

`T` is any int or float type `@Vector` takes; float-only operations say so. Complex is not in the
first set, since `@Vector` cannot hold a struct (question 8).

- **`@Vector`.** `Vec` is `@Vector`, so `*` is element-wise (MATLAB's `.*`). Zig has no operator
  overloading, so matrix products are `mul`.
- **Shaderlib.** `Vec2` to `Vec4` stay aliases of `V2f32` to `V4f32`; `Mat2` to `Mat4` become
  `Mat(N, N, f32)` (the same memory as today's `[N]VecN`, plus a field name).
- **SPIR-V.** `Vec` of 2 to 4 is `OpTypeVector`. `Mat` maps to `OpTypeMatrix` (`C` columns of
  `OpTypeVector(R)`) once Zig can declare one; until then it is a struct of an array of vectors,
  and the matrix instructions wait (section 3).
- **Layout: column-major.** The BLAS (Settled: column-major only), GLSL (`m[j]` is a column) and
  `OpTypeMatrix` all are; there is no row-major form. `fromRows` and `transpose` take row data.
- **Column stride** is `@sizeOf(Vec(R, T))`: for `R = 3` a padded column (16 bytes for `f32` on
  the tier 1 and 2 targets; checked by a comptime assert), which is also std140's and std430's
  `mat3` column stride. The BLAS path passes `ld = @sizeOf(Vec(R, T)) / @sizeOf(T)`.

## 2. Names

One file, `linalg/vector.zig`, free functions over `anytype` (as today), inferring `N`, `R`, `C`
and `T`. MATLAB names where natural; shaderlib re-exports them under GLSL names.

### Vectors

| Name | Meaning | MATLAB / GLSL | SPIR-V |
| --- | --- | --- | --- |
| `splat(V, x)`, `zero(V)` | every element `x`, or 0 | `repmat` / `vec3(x)` | `OpCompositeConstruct` |
| `dot(a, b)` | `Σ a[i] b[i]` | `dot` / `dot` | `OpDot` |
| `dot2(v)` | `dot(v, v)` | | `OpDot` |
| `length(v)`, `norm(v)` | Euclidean length, overflow-safe | `norm` / `length` | `Length` |
| `distance(a, b)` | `length(b - a)` | `norm(b - a)` / `distance` | `Distance` |
| `normalizeOr(v, fallback)` | `v / length(v)`, or `fallback` for 0 | / `normalize` plus select | `Normalize`, select |
| `cross(a, b)` | 3D cross product | `cross` / `cross` | `Cross` |
| `cross2(a, b)` | 2D perp-dot, a scalar | | none |
| `outer(a, b)` | `a * bᵀ`, a `Mat` | `a * b'` / `outerProduct` | `OpOuterProduct` (†) |
| `mix(a, b, t)` | `a (1 - t) + b t` | / `mix` | `FMix` |
| `mulAdd(T, a, b, c)` | `a b + c`, fused per policy | / `fma` | `Fma` |
| `sum(v)`, `prod(v)` | reductions | `sum`, `prod` | adds, multiplies |
| `max(v)`, `min(v)` | reductions; first NaN wins | `max`, `min` | `FMax`/`FMin` chain |

### Matrices

| Name | Meaning | MATLAB / GLSL | SPIR-V |
| --- | --- | --- | --- |
| `identity(M)`, `fromCols`, `fromRows`, `row(m, i)` | construction, access | `eye` / `mat3(1)` | |
| `transpose(m)` | `mᵀ` | `m.'` / `transpose` | `OpTranspose` (†) |
| `mul(a, b)` | `Mat x Mat`, `Mat x Vec`, `Vec x Mat` (`a` as a row) | `a * b` / `a * b` | `OpMatrixTimesMatrix`, `OpMatrixTimesVector`, `OpVectorTimesMatrix` (†) |
| `scale(m, s)`, `add(a, b)`, `sub(a, b)` | per element | `s * m`, `a + b` | `OpMatrixTimesScalar` (†), `OpFAdd` per column |
| `trace(m)` | `Σ m[i][i]` | `trace` | none |
| `det(m)` | determinant | `det` / `determinant` | `Determinant` (†) |
| `inv(m)` | `error{Singular}!Mat` | `inv` / `inverse` | `MatrixInverse` (†) |
| `solve(a, b)` | `a \ b`, `error{Singular}!Vec` or `!Mat` | `mldivide` | none |

(†) Needs `OpTypeMatrix`; column code until Zig can declare one.

- **Renames.** `cross3` becomes `cross` (shaderlib's alias already is); `norm(v, fallback)` becomes
  `normalizeOr`, and `norm(v)` is MATLAB's 2-norm. The arity changes, so an old call fails to
  compile rather than changing meaning.
- **`length`** keeps GLSL's meaning (the Settled decision's name), not MATLAB's element count.
- **Plain `normalize`** (NaN for 0) stays in shaderlib, per the Settled decision (question 4).
- **Shaderlib's GLSL names:** `inverse = inv`, `determinant = det`, `outerProduct = outer`,
  `matrixCompMult` (element-wise per column), plus the ones it already has.

## 3. Dispatch

Comptime, per operation, element type, size and target, so a call takes the same path every time:

1. SPIR-V: the native instruction where one exists (table below), else inline code.
2. CPU, size at or below the crossover: inline `@Vector` code.
3. CPU, above it: the BLAS kernel, through a serial entry (below).

Element-wise operations (`+`, `mix`, `mulAdd`, `splat`, `cross`, `outer`, `transpose`, `scale`,
`add`) are always inline: the data is by value, in registers, and a BLAS call would store and
reload it. `det`, `inv` and `solve` on fixed sizes are always inline too (unblocked LU with
partial pivoting, `getf2`'s order); LAPACK serves the dynamic layer (section 4).

### Provisional crossovers

| Operation | Inline while | Above: BLAS |
| --- | --- | --- |
| `dot`, `dot2`, `sum` | `N <= 16` | `dot`, `sum` |
| `length`, `norm`, `distance` | `N <= 16` | `nrm2` |
| `max`, `min` | `N <= 16` | `extreme(.max)`, `extreme(.min)` |
| `mul(Mat, Vec)`, `mul(Vec, Mat)` | `R * C <= 16` (4x4) | `gemv` (`n`, `t`) |
| `mul(Mat, Mat)` | `R * C * K <= 64` (4x4x4) | `gemm` (small-matrix kernels) |

M4 replaces these from the M6 and Zen 3 timings, per operation and type; one table for every CPU
target unless M4 finds an architecture that needs its own row.

- **Serial entry.** The BLAS routines take an `Io` and return `Error!T`; `linalg.vector` takes
  neither (value in, value out, as GLSL). Fixed sizes sit far below every level-1 and level-2
  serial threshold (`threading.zig`: 10000 and more for level 1), so the call never threads. An
  internal `blas.fixed` calls the kernel slot directly (`portable.fetch`) with comptime-valid
  sizes: no argument errors, no `Io`, and no allocator (gemv's buffer and the small-matrix gemm
  fit a stack array). It reads `blas.vtable()`, which selects the core on first use, as every
  routine does. It skips the Accelerate routing (which lives in the public interfaces), so macOS
  runs rstd's kernels here.
- **Configuration.** The crossovers are a comptime `Options` in `linalg/vector.zig`; a root
  declaration overrides any of them, as `blas_fuse_fma` is read from root today:

  ```zig
  pub const rstd_linalg_options: rstd.linalg.Options = .{ .crossover = .{ .dot = 32 } };
  ```

  `0` means always the BLAS, `maxInt` always inline. Comptime only: the path never depends on
  run-time state.
- **Tests** call each form directly (an internal `dotWith(.inline_, ...)`, `.blas`) at every size 1
  to 32 against the wide reference, so both paths are covered whatever the table says.

### Dispatch under SPIR-V

No BLAS and no crossover; each operation is its native instruction, used as is.

| Operation | Instruction | How |
| --- | --- | --- |
| `dot`, `dot2` | `OpDot` | inline assembly |
| `length`, `norm`, `distance`, `normalizeOr`, `cross`, `mix` | `Length`, `Distance`, `Normalize`, `Cross`, `FMix` | inline assembly, `OpExtInst` on `GLSL.std.450` |
| `mulAdd` | `Fma` | `@mulAdd` (the backend's lowering) |
| `@sqrt`, `@abs`, `@min`/`@max` and the rest | GLSL.std.450 | already, through the backend |
| matrix operations (†) | `OpMatrixTimesVector` and the rest | once Zig has `OpTypeMatrix`; column code meanwhile (`Σ cols[j] * splat(v[j])`) |

- `N` outside 2 to 4 is not an `OpTypeVector`, so those sizes take the inline code.
- `mulAdd` changes: today `blas/fma.zig` says unfused on SPIR-V; the Settled decision lists `Fma`,
  so `linalg.vector.mulAdd` takes `@mulAdd` there, and the CPU keeps `fma.zig`'s policy.
- The missing matrix type goes to `blas/issues.md` (not filed upstream).
- S1 checks each form by disassembly; S2's build test keeps them native.

## 4. Larger objects

### Types

| Type | Fields | Use |
| --- | --- | --- |
| `Matrix(T)` | `data: []T`, `rows`, `cols`; column-major, `ld = rows` | owned: `init(allocator, rows, cols)`, `deinit(allocator)`, `zeros`, `eye`, `clone`, `view()` |
| `View(T)`, `ConstView(T)` | `data`, `rows`, `cols`, `ld` | borrowed: `init(slice, rows, cols, ld)` (checks the slice covers it), `sub`, `col` |
| vectors | `[]T`, `[]const T`; `Strided(T)` (`data`, `len`, `inc`) for a row | |

`T` is any float type or its `std.math.Complex` (interleaved), as the BLAS takes. A fixed `Mat`
gives a `ConstView` (`ld` from its column stride), so the dynamic functions take it too.

### Conventions

From `agent_rules.md` and design.md's "Planned: a Zig-idiomatic interface":

- `(allocator, io, ...)` where a function allocates, `(io, ...)` otherwise; allocator first.
- A new matrix comes back owned (`Matrix(T)`), freed with `deinit(allocator)` on the same
  allocator; each such function has an `...Into` form that writes caller storage.
- Inputs are const; nothing is overwritten unless the name says `InPlace`.
- Exact, named error sets per function, plain names; `OutOfMemory` only where it allocates,
  `Canceled` from `io`. Shared names: `DimensionMismatch`, `NotSquare`, `SliceTooShort`,
  `Singular`, `NotPositiveDefinite`, `NoConvergence` (LAPACK's `INFO > 0` in SVD and eigen).
- LAPACK's argument errors (`InvalidLda`) cannot arise, since `linalg` forms the arguments; one
  would be a bug, so it is `unreachable`.

### Functions

Re-exported flat at `rstd.linalg` from `linalg/dense.zig` (`Matrix.zig` holds the types). Errors
listed are beyond `OutOfMemory` and `Canceled`.

| Function | Result | Backed by | Errors |
| --- | --- | --- | --- |
| `dot(io, x, y)` | `T` (`xᴴ y` for complex, as MATLAB) | `dot`, `dotc` | `DimensionMismatch` |
| `norm(io, x)`, `sum`, `max`, `min` | scalar | `nrm2`, `sum`, `extreme` | |
| `axpy(io, a, x, y)`, `scale(io, a, x)` | in place | `axpy`, `scal` | `DimensionMismatch` |
| `mul(allocator, io, a, b)` | `Matrix` | `gemm`; `gemv` for a vector | `DimensionMismatch` |
| `transpose`, `ctranspose` | `Matrix` | `omatcopy` | |
| `solve(allocator, io, a, b)` | `Matrix`: `A \ B` | triangular: `trtrs`; Hermitian with a positive real diagonal: `potrf`, `potrs`, LU if not positive definite; square: `gesv`; rectangular: `gelsy` (QR with column pivoting, as MATLAB) | `DimensionMismatch`, `Singular` |
| `lu` | `Lu` (factors, `ipiv`, `l()`, `u()`, `zero_pivot: ?usize`) | `getrf` | none; singular is not an error (question 5) |
| `chol` | `Matrix` (upper `R`, `A = Rᴴ R`) | `potrf` | `NotSquare`, `NotPositiveDefinite` |
| `qr(..., .economy or .full)` | `Qr` (`q`, `r`) | `geqrf`, `orgqr` | |
| `svd(..., .values, .economy or .full)` | `Svd` (`s`, `u`, `v`) | `gesdd` or `gesvd`, per `lapack/sources.md` | `NoConvergence` |
| `eig(..., want_vectors)` | `Eig` (`values: []Complex`, `vectors`) | `geev`; Hermitian input: `syevd`/`heevd` | `NotSquare`, `NoConvergence` |
| `eigh` | `EigH` (real `values`, `vectors`) | `syevd`, `heevd` | `NotSquare`, `NoConvergence` |
| `inv` | `Matrix` | `getrf`, `getri` (T117; `getrs` on `I` until then) | `NotSquare`, `Singular` |
| `pinv(..., tol: ?F)` | `Matrix` | SVD; default `tol` `max(m, n) eps(σ₁)`, as MATLAB | `NoConvergence` |
| `rank(..., tol: ?F)` | `usize` | SVD values | `NoConvergence` |
| `det` | `T`; 0 when singular | `getrf` | `NotSquare` |
| `cond` | `σ₁ / σₙ` | SVD values | `NoConvergence` |
| `matrixNorm(..., .one, .inf, .fro or .two)` | `F` | `lange`; `.two` by SVD | |

`solve`'s structure tests read only the values, so the same input always takes the same route.
The 1980s port, `linalg.matlab`, keeps its own API (split complex) until M2; this layer does not
depend on it.

## 5. The `rstd.simd` split

| Function | Goes to | Note |
| --- | --- | --- |
| `mix` | `linalg.vector.mix` | |
| `mulAdd` | `linalg.vector.mulAdd` | the BLAS's own use imports `blas/fma.zig`, so the BLAS never imports `linalg.vector` |
| `vectorLen` | stays | SIMD width, not math |
| `zigvec.V2f32` to `V4f64` | `linalg.vector` | |
| `zigvec.SimdVecF32`, `SimdVecF64`, `f32_simd_lanes`, `f64_simd_lanes` | stay | width types |
| `zigvec`'s 14 splat and zero helpers | `linalg.vector.splat`, `zero`, or plain `@splat` | deleted |
| `reduce.reduce` | stays | |
| `arraySum` | float slice: `linalg.sum` (BLAS `sum`); vector: `linalg.vector.sum` | ints and arrays stay |
| `arrayMax`, `arrayMin` | float slice: `linalg.max`, `min` (BLAS `extreme`); vector: `linalg.vector.max`, `min` | NaN rule changes (below); ints stay |
| `arrayProduct` | stays | no BLAS counterpart; `linalg.vector.prod` is inline |
| `arrayAnd`, `arrayOr`, `arrayXor` | stay | |
| `mapReduce` | stays | its matlab callers become BLAS calls in M2 |

**NaN rule.** `@reduce(.Max)` follows `@max`, which returns "the biggest non-NaN value included"
(langref, `@max`). The BLAS's extremes and `linalg`'s `max`/`min` on both paths return NaN when
any element is NaN, the first one (`blas/comparisons.md`, (3)). A caller that relied on NaN being
skipped must filter first. No caller outside `src/simd/` uses them today.

### Callers that move

| File | Uses | Becomes |
| --- | --- | --- |
| `linalg/vector.zig` | `zigvec` types, `splatAs`, `zeroAs`, `mulAdd` | defines them |
| `linalg/matlab/blas1.zig` | `mulAdd`; `f64_simd_lanes`, `SimdVecF64`; `mapReduce` (`wnrm2`, `wasum`, `wdot*`) | `linalg.vector.mulAdd`; width types stay; M2 makes the reductions BLAS `nrm2`, `asum`, `dot`, `dotc` (after its split-complex decision) |
| `linalg/matlab/hessenberg.zig` | `mapReduce(cabs1)` | M2: BLAS `asum` (complex `asum` is `Σ abs(re) + abs(im)`) |
| `linalg/blas/blas1/axpy/__root__.zig` | `vectorLen`, `mulAdd` | `vectorLen` stays; `fma.mulAdd` |
| `rgl/font/sfnt/aux_types.zig` | `zigvec.V2f32`, `V2f64` | `linalg.vector` |
| `rgl/font/sfnt/msdf/math.zig` | `zigvec.V2f64`, `splatV2f64`, `rsimd.mix` | `linalg.vector.V2f64`, `splat`, `mix` |
| `rgl/font/sfnt/msdf/` `EdgeSegment`, `Shape`, `ErrorCorrection`, `coloring`, `__root__` | `math.splatV2f64`, `math.mix` | through `msdf/math.zig`; `splatV2f64` becomes `splat` |
| `rgl/shaders/types.zig` | `zigvec.V2f32` to `V4f32`, `splatVecF32`; `Mat2` to `Mat4` | `linalg.vector`; `Mat(N, N, f32)` |
| `rgl/shaders/lib/math.zig` | `rsimd.mix`, `reduce.arraySum`; `vec.norm`, `vec.cross3` | `linalg.vector.mix`, `sum`; `normalizeOr`, `cross` |
| `rgl/shaders/lib/coords.zig`, `sdf.zig`, `shaders/rgl.zig` | `types.splatVec` | `types`' re-export of `splat` |
| `_test_qemu.zig`, `__root__.zig` | `simd` root | unchanged |

Not moved: glyf's `transformPoint` (`tables/glyf/shape.zig`): hand-written 2x2 rows with
`@reduce`, not a `simd` caller, and `mul` could fuse differently and move truncated glyph
coordinates. Nothing in RandyOS or rgui outside rstd imports `rstd.simd` or `linalg.vector`;
zfls is not checked here.

## 6. Bit guarantees

### Bits on the CPU

| Path | Promise |
| --- | --- |
| inline | a function of the inputs, operation, type, size and build target (`fma.zig`'s policy): the same bits every call, thread and run-time core |
| BLAS | the BLAS's: the same bits for any thread count and schedule, for a given build and core; one contract on every core; bits may differ between cores |
| both | one contract: first NaN wins in `max`/`min`; `length`/`norm` overflow- and underflow-safe; element-wise operations IEEE per element; the caller's floating-point environment untouched; no address dependence; for `N <= 4` the same operation order as today, so font and SDF bits do not change |

Crossing a crossover may change rounding (summation order), never meaning. The dynamic layer
carries the BLAS's and LAPACK's promises.

### Bits under SPIR-V

No bit promise. Vulkan's precision rules apply
(<https://registry.khronos.org/vulkan/specs/latest/html/vkspec.html#spirvenv-precision-operation>):
the driver sets `OpDot`'s, `Length`'s and the rest's precision, may fuse `a*b + c` and may flush
denormals; `FMin`/`FMax` with a NaN operand are undefined, so the first-NaN rule does not hold.
The promise there is the same API and meaning, compiled to the native instruction.

## 7. Migration plan

Each step is its own commit for review; tests per step as `agent_rules.md` says (native
`x86_64-windows`, compile checks of the tier-2 rigs, an rgl shader build for SPIR-V).

1. **`linalg.vector`, inline only** (crossovers off): types, names, `splat`/`zero`, `mix`,
   `mulAdd`; every caller in section 5 moves in the same commit, with no shims left in
   `rstd.simd`. Tests: today's, plus bit-for-bit for `N` 2 to 4 against today's formulas.
2. **Fixed matrices:** `Mat`, `mul`, `transpose`, `outer`, `det`, `inv`, `solve`; shaderlib's
   `Mat2` to `Mat4`.
3. **The BLAS path:** `blas.fixed`, the crossover table and `Options`; both paths tested. M4
   then times them and sets the table.
4. **SPIR-V forms**, after S1; shaderlib's aliases (S2).
5. **Dynamic layer:** `Matrix`, `View`, the level 1 to 3 wrappers, then `solve`, `lu`, `chol`,
   `qr`, `svd`, `eig`, `eigh`, `inv`, `pinv`, `rank`, `det`, `cond`, `matrixNorm`.
6. **Float-slice reductions** to `linalg` (the NaN change), and `rstd.simd` trimmed to section 5.
7. **M2** moves the matlab port onto this layer.

## Open questions

1. **Matrix storage:** `[C]Vec(R, T)` (padded `mat3` columns; matches std140/std430 and
   `OpTypeMatrix`) or packed `[C][R]T`? **Recommend** `[C]Vec(R, T)`.
2. **Core-dependent bits above a crossover:** a fixed-size call routed to the BLAS gives the
   selected core's bits. Accept that, or keep fixed sizes always inline and the BLAS for slices
   only? **Recommend** accept (the Settled decision routes by size); the inline path stays
   core-independent.
3. **Crossover override:** a root declaration or a build option? **Recommend** the root
   declaration: no build plumbing, comptime, and tests force each path anyway.
4. **Plain `normalize`:** the Settled decision keeps it in shaderlib. Keep that, with
   `linalg.vector.normalizeOr` for CPU callers? **Recommend** yes; no reason to reopen it.
5. **`lu` of a singular matrix:** return the factors with `zero_pivot` set (MATLAB's `lu` does
   not fail) while `solve`, `inv` return `error.Singular`? **Recommend** yes; the factorization is
   complete and valid, unlike `potrf`'s.
6. **`eig`'s result type:** always complex values (`geev`, or `syevd` with zero imaginary parts),
   plus `eigh` for real Hermitian results? **Recommend** both.
7. **`linalg.blas` and `linalg.lapack`:** keep them public? **Recommend** keep them reachable for
   tools, tests and experts, with the docs pointing at `rstd.linalg`'s own layer.
8. **Complex fixed-size vectors:** none in the first set? **Recommend** defer; the dynamic layer
   takes complex.
9. **`inv` under SPIR-V:** `MatrixInverse` reports nothing for a singular matrix. Keep one
   signature by testing `Determinant == 0` for `error.Singular` there? **Recommend** yes.
10. **Order against todo.md:** M3 says "Needs: M1, S1". Only step 4 needs S1, and only step 7
    needs M1. Start steps 1 to 3 and 5 now? **Recommend** yes, after the tier-1 BLAS work, as the
    MVP note says.
11. **No matrix type in Zig's SPIR-V backend:** record it in `blas/issues.md` and lower matrix
    products as column code until then? **Recommend** yes.
