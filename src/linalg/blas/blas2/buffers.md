# Level-2 kernel buffers

What each ported level-2 kernel does with its `buffer` argument, read from
its OpenBLAS source (v0.3.34-317), and what `interface.zig` gives it. The
kernel tests (`guard.zig`) hand every kernel exactly this buffer, with
canaries on both sides, for shapes up to 129.

## What the interface allocates

| routine | reference (`interface/`) | here |
| --- | --- | --- |
| `gemv`, `zgemv` | `STACK_ALLOC` of `cs*(m+n) + 128/sizeof + 3 & ~3` elements, plus `160/sizeof` on Windows; 2048 bytes on the stack (`aligned(0x20)`), else a page-aligned block | the same size (`symmetric.gemvLen`); stack up to 2048 bytes (32-aligned), else the allocator, 4096-aligned |
| `ger`, `zger` | unit strides and `m*n <= 2048*4`: NULL; else `STACK_ALLOC` of `cs*m` | unit and small: the unused stack array, no capacity promised; else `cs*m` elements, stack or heap as `gemv` |
| `symv`, `zsymv`, `zhemv` | `blas_memory_alloc`, page-aligned, `BUFFER_SIZE` (32 MB) | always the allocator, 4096-aligned, the larger of the two contracts below |
| `trmv` and the other drivers | own copies of `x` and `y`, gemv's buffer after them | same, 2048 bytes on the stack, else the allocator |

`interface.stack_bytes` (`MAX_STACK_ALLOC`) is the one stack limit for
`interface.zig`, `triangular.zig` and `banded.zig`.

## What each kernel needs

`cs` is 2 for complex, `es` an element's bytes (`cs * sizeof(F)`).

| Kernels (cores) | use of `buffer` | needs |
| --- | --- | --- |
| x86_64 `dgemv_n.S`, `dgemv_n_atom.S`, `dgemv_n_bulldozer.S` | `y` accumulated there when `incy != 1` or `y` is misaligned: `16*SIZE` head, a `SIZE` shift, then `m` | `m + 17` elements; 16-byte aligned (`movapd`) |
| x86_64 `dgemv_t_atom.S`, `dgemv_t_bulldozer.S` | `x` copied there, `16*SIZE` head | `m + 17` elements; 16-byte aligned |
| x86_64 `zgemv_n_atom.S`, `zgemv_n_dup.S`, `zgemv_t_atom.S`, `zgemv_t.S` (bobcat, opteron) | `y` (n) or `x` (t) of `2m`, `16*SIZE` head | `2m + 16` elements; 16-byte aligned; no NaN or infinity in it for `t` forms (the interface zeroes it) |
| x86_64 `symv_[LU]_sse.S`, `symv_[LU]_sse2.S`, `zsymv_[LU]_sse2.S` (also `zhemv`) | `x` copied (times `alpha`), the pointer rounded up to the next 512-byte boundary, `y` copied, blocks of 4 read past each | `2*m*es + 1280` bytes (`sseBytes`), any alignment of 16; no NaN in it (they read past the copies; the interface zeroes it) |
| x86 `gemv_*.S`, `zgemv_*.S` (x87 and SSE) | `x` or `y` copy of `m` | `cs*m` elements |
| arm `*gemv_*_vfp*.S`, arm64 `gemv_*.S`, `zgemv_*.S`, power `zgemv_n.S` | none | nothing |
| loongarch64 `gemv_[nt].S`, `zgemv_[nt].S` | `y` (n) or `x` (t) copy of `m` when its stride is not 1 | `cs*m` elements |
| loongarch64 `*gemv_*_lsx.S`, `*_lasx.S`, `[ds]symv_*` | defined, never used | nothing |
| power `gemv_[nt].S`, `zgemv_t.S`, `*_ppc440.S` | `y` or `x` copy of `m`, addressed as `CO = BUFFER - SIZE` and stored from `1*SIZE(CO)`, the buffer's first element | `cs*m` elements |
| sparc `gemv_[nt].S`, `zgemv_[nt].S` | `y` or `x` copy of `m` | `cs*m` elements |
| generic `ger.c`, `zger.c` (every core's `ger`) | `x` copied there when `incx != 1` | `cs*m` elements for a strided `x`, nothing otherwise |
| generic `symv_k.c`, `zsymv_k.c`, `zhemv_k.c` (`symmetric.Compose`) | a `SYMV_P` square block, then `y` and `x` copies (strided only) and gemv's buffer, each from a page boundary | from a page-aligned buffer, `composeBytes` |

Every gemv need is inside the interface's `cs*(m+n) + 128/sizeof` elements,
so the interface size covers every core. No ported kernel needs more than
16-byte alignment except the generic symv layout, which rounds to pages
itself and is sized for a page-aligned start.

The `ger` size is the reference's own. No ported `ger` kernel touches the
buffer for a unit stride, which is what the small unit-stride case relies on.
No core has its own `ger` kernel yet (x86_64 `sger.c` and `dger.c` for
Sandy Bridge, arm64 `ger_sve_v1x3.c`), so none needs a buffer flag.
