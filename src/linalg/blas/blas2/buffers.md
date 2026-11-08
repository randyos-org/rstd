# Level-2 kernel buffers

What each ported level-2 kernel does with its `buffer` argument, read from its OpenBLAS source (v0.3.34-317), and what `interface.zig` gives it. The kernel tests (`_test/guard.zig`) hand every kernel exactly this buffer, with canaries on both sides, for shapes up to 129 (and `gemv`'s past 4096). The interfaces' own stack buffers hold NaN in tests, as the allocators do. The level-2 design is in [`design.md`](../../design.md#level-2).

## What the interface allocates

| routine | reference (`interface/`) | here |
| --- | --- | --- |
| `gemv`, `zgemv` | `STACK_ALLOC` of `cs*(m+n) + 128/sizeof + 3 & ~3` elements, plus `160/sizeof` on Windows; 2048 bytes on the stack (`aligned(0x20)`), else a page-aligned block | the same size (`symmetric.gemvLen`); stack up to 2048 bytes (32-aligned), else the allocator, 4096-aligned |
| `ger`, `zger` | unit strides and `m*n <= 2048*4`: NULL; else `STACK_ALLOC` of `cs*m` | unit `x`: the unused stack array, no capacity promised; strided: `cs*m` elements, stack or heap as `gemv`, where the interface copies `x` once and passes the copy, at unit stride, to every kernel call |
| `symv`, `zsymv`, `zhemv` | `blas_memory_alloc`, page-aligned, `BUFFER_SIZE` (16 to 128 MB by arch; 128 MB at `common_x86_64.h:260`) | always the allocator, 4096-aligned, the larger of the two contracts below |
| `trmv` and the other drivers | own copies of `x` and `y`, gemv's buffer after them | same, 2048 bytes on the stack, else the allocator |

`interface.stack_bytes` (`MAX_STACK_ALLOC`) is the one stack limit for `interface.zig`, `triangular.zig` and `banded.zig`.

## What each kernel needs

`cs` is 2 for complex, `es` an element's bytes (`cs * sizeof(F)`).

| Kernels (cores) | use of `buffer` | needs |
| --- | --- | --- |
| x86_64 `dgemv_n.S`, `dgemv_n_atom.S`, `dgemv_n_bulldozer.S` | `y` accumulated there when `incy != 1` or `y` is misaligned: `16*SIZE` head, a `SIZE` shift, then `m` | `m + 17` elements; 16-byte aligned (`movapd`) |
| x86_64 `dgemv_t_atom.S`, `dgemv_t_bulldozer.S` | `x` copied there, `16*SIZE` head | `m + 17` elements; 16-byte aligned |
| x86_64 `sgemv_n_4.c`, `dgemv_n_4.c` (`gemv/driver.zig`) | `y` summed there, zeroed first, when `incy != 1`: one row block (`NB`, at most `m`) | `m` elements |
| x86_64 `cgemv_n_4.c`, `zgemv_n_4.c` (`gemv/driver.zig`) | `y` summed there for each row block (`NB`, at most `m`), zeroed first | `2m` floats |
| x86_64 `cgemv_t_4.c`, `zgemv_t_4.c` (`gemv/driver.zig`) | `x` copied there when `incx != 1` (`NB`) | `2m` floats |
| x86_64 `sgemv_t_4.c`, `dgemv_t_4.c` (`gemv/driver.zig`) | `x` copied there when `incx != 1` (`NB`), and from `min(m, NBMAX)` on four sums per four columns, up to `NBMAX` of them | `min(m, NBMAX) + min(n, NBMAX)` elements, at most `m + n` |
| x86_64 `zgemv_n_atom.S`, `zgemv_n_dup.S`, `zgemv_t_atom.S`, `zgemv_t.S` (bobcat, opteron) | `y` (n) or `x` (t) of `2m`, `16*SIZE` head | `2m + 16` elements; 16-byte aligned; any contents (below) |
| x86_64 `symv_[LU]_sse.S`, `symv_[LU]_sse2.S`, `zsymv_[LU]_sse2.S` (also `zhemv`) | `x` copied (times `alpha`), the pointer rounded up to the next 512-byte boundary, `y` copied, blocks of 4 read past each | `2*m*es + 1280` bytes (`sseBytes`), any alignment of 16; no NaN in it (they read past the copies; the interface zeroes it) |
| x86 `gemv_*.S`, `zgemv_*.S` (x87 and SSE) | `x` or `y` copy of `m` | `cs*m` elements |
| arm `*gemv_*_vfp*.S`, arm64 `gemv_*.S`, `zgemv_*.S`, power `zgemv_n.S` | none | nothing |
| loongarch64 `gemv_[nt].S`, `zgemv_[nt].S` | `y` (n) or `x` (t) copy of `m` when its stride is not 1 | `cs*m` elements |
| loongarch64 `*gemv_*_lsx.S`, `*_lasx.S`, `[ds]symv_*` | defined, never used | nothing |
| power `gemv_[nt].S`, `zgemv_t.S`, `*_ppc440.S` | `y` or `x` copy of `m`, addressed as `CO = BUFFER - SIZE` and stored from `1*SIZE(CO)`, the buffer's first element | `cs*m` elements |
| sparc `gemv_[nt].S`, `zgemv_[nt].S` | `y` or `x` copy of `m` | `cs*m` elements |
| generic `ger.c`, `zger.c` (every core's `ger`) | `x` copied there when `incx != 1` | `cs*m` elements for a strided `x`, nothing otherwise |
| x86_64 `sgemv_n_4.c` on Skylake-X, its few-column path (`m <= 16384`, `n <= 48`, `n != 4`; `gemv/driver.zig`) | `x` copied to the start when `incx != 1`, `y` after `n` elements when `incy != 1` | `m + n` elements |
| x86_64 `sgemv_t_4.c` on Skylake-X, its path for `lda == m <= 8`, `n <= 16384` (`gemv/driver.zig`) | `x` copied to the start when `incx != 1`, `y` after `m` elements when `incy != 1` | `m + n` elements |
| x86_64 `?symv_[LU].c` (`symv/driver.zig`) | none | nothing |
| arm64 `symv_[LU]_asimd_4x4.c` (`symv/driver.zig`) | `y` copied there when `incy != 1`, then `x` when `incx != 1` (after `y`'s copy, else at the start) | `2m` elements at most |
| x86_64 `sger.c`, `dger.c` (Sandy Bridge; `ger/x86_64/kernel.zig`) | `x` copied there when `incx != 1` | `m` elements for a strided `x`, nothing otherwise |
| generic `symv_k.c`, `zsymv_k.c`, `zhemv_k.c` (`symmetric.Compose`) | a `SYMV_P` square block, then `y` and `x` copies (strided only) and gemv's buffer, each from a page boundary | from a page-aligned buffer, `composeBytes` |

Every gemv need is inside the interface's `cs*(m+n) + 128/sizeof` elements, so the interface size covers every core. No ported kernel needs more than 16-byte alignment except the generic symv layout, which rounds to pages itself and is sized for a page-aligned start.

The x86_64 complex `t` forms (`zgemv_t.S`, `zgemv_t_atom.S`) copy `x` to the buffer's start and load it ahead of use, so past the last row they load what the buffer holds; their row tails (`.L15`, `.L17` and kin) multiply only what was loaded for rows that exist, so those values never reach `y`. The reference clears its x86_64 stack buffer for `cgemv_t.S` (`interface/zgemv.c:250`, <https://github.com/OpenMathLib/OpenBLAS/issues/746>), which no core here takes ([`coverage.md`](coverage.md)); the interface does not, and the kernel tests hand every form a NaN buffer.

The `ger` size is the reference's own. No ported `ger` kernel touches the buffer for a unit stride, which the interface relies on: every call it makes takes a unit-stride `x`, the caller's or its copy in the buffer. The per-core ones (x86_64 `sger.c` and `dger.c` for Sandy Bridge) copy only a strided `x`, so none needs a buffer flag. arm64's `ger_sve_v1x3.c` is not ported yet.
