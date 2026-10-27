# BLAS threading design

This document specifies how rstd's BLAS routines use multiple threads: the public API, how
work is divided, how results stay reproducible, and how the design degrades on targets
without threads. Plan step 1 is implemented (`threading.zig`); no routine threads yet. Upstream references are to the
OpenBLAS checkout that `tools/blas_asm` converts (v0.3.34-317); std references are to Zig's
`lib/std` at 0.18.0-dev.751. The library's overall structure is in [`design.md`](design.md);
how this rule compares with MKL, BLIS and OpenBLAS is in [`comparisons.md`](comparisons.md) (8).

## Requirements

- **Reproducible results.** For a given build and core selection, every routine returns the
  same bits for any worker count, any `Io`, any scheduling, and in threaded and single-threaded
  builds alike. Run time and the order in which work executes may vary; results may not.
  Reproducibility takes precedence over speed.
- **Adaptive parallelism.** The degree of parallelism follows the machine. No caller has to
  state a thread count for good performance.
- **std conventions.** Concurrency flows through a caller-supplied `std.Io`, passed like an
  `Allocator`.
- **Every target.** The same API serves Windows, Linux, macOS, wasm, UEFI and freestanding.
  Targets without threads run serially with no special cases in callers.
- **No hidden failure modes.** Allocation happens before any task starts; failures surface
  as errors, never panics, and are never silently converted into slower paths.
- **Memory safety.** No task outlives the call that spawned it, and nothing a task touches is
  freed or unloaded while it runs ([Guardrails](#guardrails)).

Across core types, results may differ in the last bits: different cores run different
kernels, with different fusion, accumulation order and blocking. OpenBLAS behaves the same way.
Beyond rounding they agree ([`readme.md`](readme.md#one-contract-on-every-core)), and on one
core no result depends on where its operands lie in memory. Pinning the core type
(`OPENBLAS_CORETYPE`) makes machines that can run that core agree.

## Background

### OpenBLAS

- **Server.** `driver/others/blas_server.c` keeps a global pthread pool. Workers spin
  (`YIELDING`) for `thread_timeout` ticks, then sleep on a condvar. `exec_blas` queues
  `queue[1..]` and runs `queue[0]` on the caller. Each worker has private `sa`/`sb` buffers
  from `blas_memory_alloc`. `blas_server_win32.c` mirrors this; `blas_server_omp.c` uses
  OpenMP.
- **Level 1.** `blas_l1_thread.c` splits `n` into near-equal ranges, one per thread.
  Reductions are threaded inside per-arch kernels (`kernel/x86_64/ddot.c`, `dasum.c`,
  `drot.c`, `kernel/arm64/dot.c`) as per-thread partials summed in thread order.
- **Level 2.** `driver/level2/*_thread.c`: `gemv` splits rows (N) or columns (T), or, when
  `m` is short, splits `x` into per-thread `y` buffers that are then summed. `ger` splits
  columns. `symv`/`trmv` use sqrt-balanced triangle ranges and per-thread `y` buffers combined
  with `AXPYU_K`.
- **Level 3.** `level3_thread.c` builds a `t_m x t_n` grid (`CNAME`: halve `t_m` until
  `m >= t_m * switch_ratio`, set `t_n = ceil(n / (switch_ratio * t_m))`, then rebalance
  toward square tiles). Each thread packs its slice of `B` into `DIVIDE_RATE` (2) buffers and
  publishes them through `job[].working[][]` flags; peers in its N group consume them with
  their own packed `A`. Producers spin until consumers clear the flags, so the protocol
  requires every thread to run concurrently. `trsm` uses the synchronization-free
  `gemm_thread_n`/`gemm_thread_m`.
- **Nesting.** `num_cpu_avail` returns 1 inside a parallel region.

Because levels 1 and 2 split by thread count, OpenBLAS's reductions and per-thread-buffer
routines change their low bits when the thread count changes.

### std.Io

- `Io.async`, `Io.concurrent`, `Io.Future(T)`, and `Io.Group` with `async`, `concurrent`,
  `await` and `cancel` (`Io.zig`). A `Group` task returns something coercible to
  `Cancelable!void`; other errors are discarded, so tasks cannot report `OutOfMemory`.
- `Io.Threaded.groupAsync` (`Io/Threaded.zig`) runs a task inline (`groupAsyncEager`) when
  the build is `single_threaded`, when task allocation fails, when `busy_count >=
  async_limit`, or when a thread spawn fails. `async` degrades the same way. Nested or
  oversubscribed use therefore neither deadlocks nor errors.
- `groupConcurrent` never runs inline: it returns `error.ConcurrencyUnavailable` in those
  cases. Its `concurrent_limit` defaults to `.unlimited`.
- `Threaded.init` sets `async_limit` to `getCpuCount() - 1`; `setAsyncLimit` changes it.
  `start.zig` passes such an instance to `main` as `std.process.Init.io`. `groupAwait` blocks
  on a futex and does not run queued work.
- `Io.failing` (`noAsync`/`noGroupAsync`) and rstd's UEFI `Io` (`src/_os/uefi/io/async.zig`)
  run tasks inline and fail `groupConcurrent`.
- The `Io` interface does not expose its parallelism; the limits are fields of the concrete
  `Io.Threaded`.
- `futexWait` is cancelable. `std.Thread.Pool` no longer exists, and `WaitGroup` is private
  to `Io/Threaded.zig`.
- wasm32/wasm64 default to single-threaded (`defaultSingleThreaded` in the compiler's
  `src/target.zig`).

## Design

### API

Every public BLAS routine takes `io: Io` immediately after its allocator, or first when it
has none: `gemm(allocator, io, T, ...)`, `dot(io, T, ...)`. This holds for routines that
never thread as well, so signatures change once and threading can be added to any routine
later without breaking callers. Callers pass the `Io` they already have. A single-threaded
`Io` (`Io.failing`, `Io.Threaded.global_single_threaded`) gives serial execution.

There is no per-call thread count. Thread policy belongs to the owner of the `Io`: a caller
wanting less parallelism configures it (`Io.Threaded.setAsyncLimit`) or passes a smaller
dedicated instance.

### Worker count

`blas.init` caches `std.Thread.getCpuCount() catch 1` in `blas.config`. A call uses

```text
w = min(cpu_count, ceiling, size_max)
```

workers, where `ceiling` is `OPENBLAS_NUM_THREADS`, else `GOTO_NUM_THREADS`, else unbounded,
and `size_max` is the most workers the routine's threshold justifies for the problem size.
`OMP_NUM_THREADS`, `OMP_ADAPTIVE` and `thread_timeout` remain parsed in `Config` but have no
effect, since there is no spinning server. The calling thread is worker 0.

The worker count affects only speed. Targets that report one CPU (freestanding builds, the
qemu rigs through their `Thread` override) run serially without special cases. When the `Io`
admits fewer than `w` workers, the surplus run inline on the caller.

### Levels 1 and 2: fixed chunks, dynamic claiming, ordered combine

Above a routine's threshold the problem is cut into chunks whose boundaries depend only on
the problem size and the routine's chunk size, never on `w`:

- **Chunks.** A chunk is about the most work OpenBLAS's threshold runs serially: a call of
  `work` units above `serial_max` makes `ceil(work / serial_max)` chunks, at most 64. Bounds
  split the range evenly, or by square roots where row costs grow or shrink (`Shape`), and
  round to a quantum the routine picks.
- **Claiming.** `w` workers start through `Group.async`; each claims the next chunk from an
  atomic counter until none remain. The work spreads over whatever actually runs: all
  workers when the `Io` runs them all, fewer when it admits fewer (late starters find nothing
  left), and the remaining workers absorb the load when a core is busy.
- **Combining.** Each chunk's partial result is determined by its boundaries, not by the
  worker that computed it. Reductions keep partials in a stack `[64]T` (for `nrm2`, a
  `(scale, ssq)` pair) and fold them in chunk order on the caller; `i?amax` combines exactly
  (first index wins). Level 2 routines with partial vectors (`gemv` with short `m`, `symv`,
  `trmv`) keep one per chunk and combine them in chunk order.

Element-wise routines (`axpy`, `scal`, `swap`, `rot`, `gemv` over rows, `ger`) compute each
output in exactly one chunk, so the split cannot change their results. Below the threshold a
routine makes one kernel call, as in a single-threaded build, which gives the same bits.

### Level 3: OpenBLAS's shared packed B

`gemm` and the other level 3 routines follow `level3_thread.c`: the `CNAME` grid, with
threads in an N group sharing each other's packed `B`. Because producers wait for consumers,
every worker must run concurrently, and `Group.async`'s inline fallback would deadlock. The
protocol is therefore built on `Group.concurrent`:

- **All-or-nothing start.** The caller spawns `w - 1` workers, each of which first waits on a
  start barrier. The caller counts the workers that actually started, builds the grid for
  that count, then releases them. With none started it runs serially. No worker ever waits
  for a peer that does not exist.
- **Cap.** `concurrent_limit` defaults to unlimited, so `w` is the cap.
- **Waiting.** Workers wait with `io.futexWait`/`futexWake` rather than spinning, so a
  fiber-based `Io` (`Io.Evented`) is not blocked.
- **Reproducibility.** The grid depends on the worker count, but every element of `C` is
  computed by one worker and the `k` blocking (`Blocking.split`) does not depend on the
  split, so results match the serial path bit for bit.

`trsm`, `syrk` and `symm` use the same mechanisms, splitting as `gemm_thread_n`/
`gemm_thread_m` and the triangle-balanced variants do.

### Thresholds

| Routine | OpenBLAS threshold | Division |
| --- | --- | --- |
| real `gemm` | `MNK > 65536 * 4`; `t = min(ncpu, MNK / 262144)` | `CNAME` grid, shared packed `B` |
| complex `gemm` | `MNK > 8192 * 4` | same |
| `syrk`, `symm` | `(n+1)nk` over `109944 * 4`; `2m^2 n` over `65536 * 4` | triangle-balanced grid |
| `trsm` | `m * n >= SMP_FACTOR (256 real, 128 cplx) * 4` | left: N panels; right: M panels |
| `gemv` | `m * n >= 115200 * 4` | N: row chunks; T: column chunks |
| `ger` | `m*n > 2048*4`: up to 3; `> 20480*4`: all | column chunks |
| `symv`, `trmv` | `symv n >= 200`; `trmv n < 50`: 1, `< 500`: 2 | sqrt-balanced row chunks |
| `syr`, `spr`, `syr2`, `spr2` (real and complex symmetric) | `n < 50`: 1, `< 500`: 2 | sqrt-balanced column chunks |
| `her`, `hpr`, `her2`, `hpr2` | `n < 100`: 1 | same |
| `axpy` / `scal` / `swap` | `n > 10000` / `n > 1048576` / `n >= 8 MiB / sizeof(T)` | range chunks |
| `dot`, `asum`, `sum`, `rot` | x86_64 kernels: `n > 10000` (dot), `> 100000` | range chunks, ordered combine |

The generic thresholds apply with `GEMM_MULTITHREAD_THRESHOLD = 4` as a comptime constant,
together with OpenBLAS's per-core overrides (`neoversev1`/`v2` in `interface/gemm.c` and
`gemv.c`, `a64fx` and `neoversev1` in `arm64/dot.c`). Complex routines have their own
constants; `threading.zig`'s `rule` holds them all. These are defaults; per-core tables derived from rstd's own benchmarks will
replace them, and the chunk sizes, later. `switch_ratio` already lives in the vtable
(`blas3/gemm/params.zig`). M ranges round to `unroll_m` and a cache line to avoid false
sharing in `C`. The rank-update cut-offs only pick the thread count; their serial small-`n`
paths are already ported ([`design.md`](design.md), Level 2).

### Memory

All scratch is allocated by the caller before any task starts, in one `alignedAlloc` sized
for `w` workers: `gemm` packing buffers (`driver.Work(T).lens`, about `(P + min(n_t, R)) * Q`
elements per worker, aligned to `Work(T).alignment`) and level 2 per-chunk vectors. Tasks are
therefore infallible, and the allocator need not be thread safe. Level 1 needs none. Task
arguments are a pointer to the shared job plus a worker index; `Group.async` copies them with
the `Io`'s own allocator and runs the task inline if that allocation fails.

### Errors and cancellation

- **Out of memory.** `error.OutOfMemory` propagates to the caller with nothing started. The
  routine never retries with fewer workers or serially: the caller decides, per Zig's
  allocator convention, and memory pressure never shows up as an unexplained slowdown.
- **Cancellation.** Every routine that takes `io` has `error.Canceled` in its error set (all
  but the scalar `rotg`, `zrotg` and `rotmg`). Workers check for cancellation at chunk
  boundaries and while waiting; on cancellation the routine stops its workers, waits for them
  to finish, and returns `error.Canceled`. The output is then
  undefined (for `gemm`, `C` may be partly scaled and partly computed), and the routine docs
  say so. Level 3 wakes peers blocked on the shared `B` flags so they can exit.
- **Inline tasks.** A task the `Io` runs inline shares the caller's cancellation, and `Group`
  drops a task's `error.Canceled`, so a worker that takes one calls `io.recancel()`. An `Io`
  whose `checkCancel` is `Io.unreachableCheckCancel` (`Io.failing`, UEFI's) has no
  cancellation; workers skip the check there.

### Single-threaded builds and freestanding targets

Threading compiles out only under `builtin.single_threaded`: such builds take a comptime
serial path that never references `Io.Group`, over the same chunks. Freestanding and UEFI are not forced serial, so
a multithreaded freestanding build (for example a RandyOS kernel with its own `Io`) can
thread. Today's freestanding rigs build `single_threaded` and report one CPU; UEFI's `Io` runs
tasks inline.

### Vendor backends

A backend that delegates to a vendor library (Apple's Accelerate on macOS) must meet the
same reproducibility requirement. If the library's own threading changes its results, rstd
runs it single-threaded. Accelerate's does (its threaded `gemv` rounds otherwise), so
`apple_accel.load` makes it single-threaded for the process and the loading thread; every other
thread sets its own model (`singleThreadHere`) on its first routed call
([`design.md`](design.md#apple-accelerate)). rstd's own workers are to parallelize around it.

### Guardrails

Reproducibility:

- **Bit-identical on a core.** On a given build and core, a routine's bits depend on its
  arguments only. A change that breaks this is a bug, whatever it gains in speed.
- **One contract across cores.** Each routine has one numerical contract on every core; cores
  differ in rounding only. A threaded path keeps its serial path's contract.
- **Split inputs.** Chunk sizes, thresholds and the level-3 `Q` read the problem size and the
  core only, never the CPU count, `OPENBLAS_NUM_THREADS` or the `Io`'s limits. Thread-count
  blocking sets (LA464's) may change `P` and `R`, never `Q`.
- **Same chunks everywhere.** Above a threshold, a `single_threaded` build and an `Io` that runs
  tasks inline combine the same chunks in the same order as a threaded run.
- **Partition edges.** Level-3 partitions round `M` ranges to `unroll_m` and `N` ranges to
  `unroll_n`, so no element moves between a full tile and a tail path.
- **No address dependence.** No kernel's bits follow an operand's alignment (`vector.peels`), so
  a chunk may start anywhere.

Memory safety for concurrent work:

- **Await before free.** A routine waits for every task it spawned, on every exit path
  (success, error, cancellation), before it frees or returns anything a task touches. `await`
  may return `error.Canceled` early; `cancel` returns only once every task has run.
- **`defer` order.** Scratch is allocated and its free deferred before the group exists, and
  `group.cancel` is deferred after, so it runs first:

  ```zig
  const work = try allocator.alloc(T, len); // before any task starts
  defer allocator.free(work); // runs last
  var group: Io.Group = .init;
  defer group.cancel(io); // runs first; every task has run
  ```

- **Caller-owned buffers.** The routine that opens the parallel region owns every buffer its
  tasks use: operands, scratch and per-chunk partials. Tasks borrow slices of them, allocate
  nothing and free nothing; a task's arguments point only at memory that lives until the region
  closes.
- **Disjoint writes.** Each output element and each partial has one writer; shared packed
  buffers are written by their producer before the flag that publishes them and read only after
  it.
- **Never unload a library.** A dynamically loaded library is never closed: Accelerate's table
  points into it, and any thread may be calling through it (`src/dll.zig`).

## Alternatives considered

| Alternative | Reason not chosen |
| --- | --- |
| A `blas.Threads { io, max }` parameter | No std precedent; the cap it carried belongs to the `Io` or `OPENBLAS_NUM_THREADS` |
| A process-global `Io` (`blas.setThreads`) | The owner may deinit it; parallel tests contend for it; hides a dependency |
| Adding `io` only as routines gain threading | Breaks signatures repeatedly; one change instead |
| Splitting levels 1 and 2 by thread count, as OpenBLAS | Results then depend on the thread count and machine |
| Sizing work from the `Io`'s limits | The `Io` interface does not expose them; reading `Io.Threaded` fields breaks other `Io`s |
| Independent level 3 tiles via `Group.async` | Redundant packing; OpenBLAS's shared-`B` protocol is followed instead |
| Retrying serially on `OutOfMemory` | Hides memory pressure and makes speed unpredictable; the serial path can fail too |
| Blocking cancellation for the call | Long calls could not be abandoned |
| Forcing serial on freestanding and UEFI | Rules out multithreaded freestanding builds |

## Verification

- **Everywhere, rigs included:** unit tests for chunking, the worker count and the ordered
  combine, and a test that the rigs resolve to one worker.
- **Hosted** (native Windows; `qemu-test-linux-*` at `-smp 4`; `ssh-test-aarch64-macos`):
  threaded results equal serial results bit-for-bit for worker counts 1 to 9 (a test hook
  overrides the cached CPU count), all transposes and odd sizes, on the shrunk-blocking
  vtable `_test/blas3/check.zig` uses and its NaN-poisoning allocator. Each case also runs on
  `Io.Threaded` with `async_limit` 0, 3 and 8 (oversubscribed) and on `Io.failing` as an
  inline baseline. Reductions additionally meet `check.zig`'s error bound.
- **Races:** the hosted tests loop about 1000 times on the weak-memory rigs (aarch64,
  ppc64le, riscv64) and the M6, with guard bands asserting nothing outside the output is
  written. TSan (`sanitize_thread`) on native Linux follows once its support in this
  toolchain is confirmed.
- **Cancellation:** canceled calls return `error.Canceled`, leave no worker running, and leave
  memory outside the output untouched; the allocator sees every free after the last task ends.
- **Split inputs:** a test that `Q` and every chunk size are the same for every thread count
  and CPU count.
- **Freestanding:** the kernel rigs build `single_threaded` (`build/qemutest.zig`,
  `addKernelRun`), and `rstd-check-uefi` keeps UEFI compiling.

## Level 3 threading prototype

`blas3/thread_proto.zig` (internal, `wip/l3proto`) threads real `gemm` both ways over the serial
driver's packing, kernels and `Blocking.split`. Both start through `threading.team`: workers
start with `Group.concurrent`, wait on a futex until the caller has counted them, and the plan
is sized for those that started; with none, the caller runs alone. Neither waits on a worker
that does not exist, so neither can deadlock when the `Io` admits fewer.

- **(A) OpenBLAS**, a port of `level3_thread.c`: the `CNAME` grid; each worker packs two halves
  of its share of its N group's `B`, publishes each with a step number, and its readers
  release it through a count. Per slice a worker waits up to `2 * (tm - 1)` times for halves
  and twice for releases. The halves need fixed offsets, since a slice can be deeper than the
  one before (`split` rounds up to `unroll_m`). 172 lines.
- **(B) BLIS**, from `docs/Multithreading.md`: `jc x ic` ways; the `ic` workers of a `jc` group
  pack each slice of `B` together into one of two buffers and meet at one futex barrier
  (`threading.Barrier`) per slice, where BLIS meets twice around one buffer. 93 lines; twice
  the `B` buffer memory of (A).

Both match `driver.gemm` bit for bit at 1 to 9 workers (`f32`, `f64`, every transpose, shapes
around every edge of the shrunk blocking, poisoned scratch) on `Io.failing` and on
`Io.Threaded` with `async_limit` and `concurrent_limit` 0, 3 and 8. A canceled call returns
`error.Canceled` with no worker left running. Row ranges must start on `unroll_m`: unrounded
ones change the last bits (row tails narrower than a vector go through `thin`).

GFLOPS on this Zen 3 (8 cores, 16 threads, core `zen`, portable Zig tiles), `C = A * B`
square, `NN`, by workers; mean of two quiet runs of medians of 9 to 25 interleaved calls:

| Type | n | serial | A 2 | B 2 | A 4 | B 4 | A 8 | B 8 | A 16 | B 16 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `f64` | 256 | 36 | 63 | 62 | 90 | 88 | 94 | 91 | 79 | 76 |
| `f64` | 512 | 36 | 70 | 71 | 124 | 113 | 143 | 141 | 147 | 156 |
| `f64` | 1000 | 39 | 81 | 76 | 136 | 120 | 134 | 134 | 186 | 155 |
| `f64` | 2000 | 50 | 82 | 82 | 129 | 126 | 154 | 147 | 162 | 155 |
| `f32` | 256 | 74 | 127 | 126 | 181 | 172 | 157 | 149 | 112 | 104 |
| `f32` | 512 | 82 | 147 | 154 | 253 | 264 | 295 | 294 | 359 | 343 |
| `f32` | 1000 | 78 | 175 | 172 | 312 | 274 | 324 | 310 | 422 | 403 |
| `f32` | 2000 | 101 | 191 | 186 | 319 | 294 | 384 | 349 | 414 | 373 |

(A) is ahead by 0 to 5% at most points and by 10 to 17% at a few (4 or 16 workers, `n` of
1000 and 2000); (B) leads at some 512s. Both scale alike past 4 workers, so the limit there is
shared (blocking, kernels, memory), not the protocol.

**Recommendation: (B)** for step 2, with (A)'s gap as its tuning target. Half the code, one
synchronization point per slice instead of `2 * tm` flag waits, one primitive (a barrier) with
a one-line deadlock argument, and room for BLIS's `jr` ways where threads share an L2 (SMT
pairs, Apple clusters). (A)'s lead likely comes from multiplying each half as it is packed;
(B) can do the same with each worker's own chunks before the barrier (not yet measured).

## Implementation plan

1. **Signatures and plumbing** (done): `io: Io` on every public routine but the scalar
   `rotg`, `zrotg` and `rotmg`, with `error.Canceled` in their error sets; in `threading.zig`,
   the cached CPU count (`blas.config.cpu_count`) and ceiling (`Config.threadCeiling`),
   thresholds (`rule`), chunking (`Chunks`, `plan`), the claim loop (`forEach`) and ordered
   combine (`reduce`), with their unit tests.
2. **Real `gemm`:** the `CNAME` grid and shared packed `B` through `Group.concurrent`, with
   the start barrier, futex waiting, per-worker packing buffers and the serial fallback.
3. **Level 1:** `axpy`, `scal`, `swap`, `rot` and the chunked reductions.
4. **Level 2:** `gemv`, `ger`, then `symv`/`trmv`.
5. **Remaining level 3:** `trsm`, `syrk`, `symm` and complex, as OpenBLAS divides them.

## Open issues

- Behavior of the futex-based level 3 protocol on `Io.Evented` is unverified.
- TSan support with this toolchain is unverified.
- Level 3 protocol for step 2: OpenBLAS's flags or BLIS's barriers
  ([Level 3 threading prototype](#level-3-threading-prototype)).
- Thresholds and chunk sizes are OpenBLAS's defaults until rstd's own per-core benchmarks
  exist.
