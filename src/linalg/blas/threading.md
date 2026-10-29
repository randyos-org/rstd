# BLAS threading design

This document specifies how rstd's BLAS routines use multiple threads: the public API, how
work is divided, how results stay reproducible, and how the design degrades on targets
without threads. The plan is implemented (`threading.zig`, `blas3/threaded.zig`, the level 1
and 2 interfaces): every level 3 routine but `gemmt` threads, and levels 1 and 2 thread where
the split keeps the serial bits ([Levels 1 and 2](#levels-1-and-2-fixed-chunks-dynamic-claiming-ordered-combine)).
Upstream references are to the OpenBLAS checkout that `tools/blas_asm` converts
(v0.3.34-317), BLIS references to its tag 2.1, std references to Zig's `lib/std` at
0.18.0-dev.751. The library's overall structure is in [`design.md`](design.md); how this rule
compares with MKL, BLIS and OpenBLAS is in [`comparisons.md`](comparisons.md) (8).

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
- **Nesting.** Under OpenMP, `num_cpu_avail` (`common_thread.h:151`) returns 1 inside a parallel
  region unless `openblas_set_num_threads_local` raised `blas_omp_threads_local` (default 1).

Because levels 1 and 2 split by thread count, OpenBLAS's reductions and per-thread-buffer
routines change their low bits when the thread count changes.

### BLIS

- **Loops.** Five loops around the microkernel: `jc` (`NC` columns), `pc` (`KC` depth), `ic`
  (`MC` rows), `jr` and `ir` (register tiles). Each `pc` step packs a `KC x NC` panel of `B`,
  each `ic` step an `MC x KC` block of `A`.
- **Threads.** Ways of parallelism per loop; `pc` is not split for `gemm`. Given a thread
  count, `bli_rntm_factorize` splits it `ic x jc` by `bli_thread_partition_2x2` (`m / ic`
  near `n / jc`), then moves up to 4 of `jc` to `jr`. Threads that share a packed buffer pack
  it together and meet at a barrier before and after (`bli_packm_int.c`).
- **Waiting.** `bli_thrcomm_barrier_atomic` is a sense-reversing barrier; waiters spin on an
  empty loop without bound.

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
  worker that computed it. Reductions keep partials in a stack `[64]T` and fold them in chunk
  order on the caller; the extremes fold exactly (first index wins). No level 2 routine keeps
  partial vectors (below).

Element-wise routines compute each output in exactly one chunk, so the split cannot change
their results; with one worker they make one kernel call, as below the threshold. A reduction
above its threshold folds its chunks whatever the worker count, in a single-threaded build too.

What threads, as landed:

- **Element-wise level 1:** `copy`, `swap`, `scal`, `rot`, `axpy` and `axpyc`. Chunks start on
  multiples of 1024 elements, so every chunk but the last is whole kernel blocks, the last ends
  where the vector does, and each keeps the vector's alignment: an element takes the body or
  tail path one call gives it. A zero stride runs serially, as OpenBLAS's do.
- **Reductions:** the dots (`dsdot` too), `asum` and `sum` at `f32` and `f64`, real and complex.
  Partials of `f32` fold in `f64` and round once; `sdot`'s partials are `dsdot`'s, so it keeps
  its contract (products and sum in `f64`, rounded once). On the x87 each step of the fold is
  rounded to its float, so no partial keeps a register's width. Wider types do not thread their sums:
  they give the bits of OpenBLAS's x87 kernels (`readme.md`), which no fold keeps.
- **Extremes:** values and indices, every type, folded exactly: the first NaN, else the first
  extreme, as one kernel call finds them. Where a complex key overflows, each chunk's kernel has
  picked by `|re|/4 + |im|/4`, and the fold does too between two such chunks, so the index is
  the serial one.
- **`gemv`:** rows without a transpose, columns with one, each chunk a serial kernel call with a
  buffer of its own. Rows split only on 4096-row blocks: x86_64's `gemv_n` kernels take rows in
  `NBMAX` blocks and run the last block's remainder first, so a finer split moves rows between
  paths (measured: `sgemv_n` on Sandy Bridge, Haswell and Zen). Columns start on multiples of 64.
- **`ger`:** columns. Its runs of nonzero `y` are cut on multiples of 64 past each run's start,
  so a run's parts make its kernel call's column blocks.
- **Rank updates** (`syr`, `spr`, `syr2`, `spr2`, `her`, `hpr`, `her2`, `hpr2`): columns, each one
  `axpy` call as in the serial driver, balanced by length (`Shape.rising`, `falling`).
- **`gbmv` with a transpose:** columns, each element of `y` one column's dot.

What stays serial, and why:

- **`nrm2`:** chunks' norms would fold with a rounding each, where `s` and `c` round the root
  once and `d` and `z` keep Blue's three sums inside the kernel.
- **`symv`, `hemv`, `sbmv`, `hbmv`, `spmv`, `hpmv`, `trmv`, `tbmv`, `tpmv` and `gbmv` without a
  transpose:** each element of `y` adds up over several kernel calls (`axpy` columns, `gemv`
  blocks, or a core's own `symv` kernel). A split of `y` cuts those calls elsewhere, and kernels
  round a call's body, tail and remainder apart, so the bits would follow the split. OpenBLAS's
  per-thread partial vectors fold a sum per thread, so theirs follow the thread count.
- **Not threaded upstream:** `axpby`, `rotm`, the triangular solves, `sbdot` and the bfloat16
  conversions (OpenBLAS threads those two on x86_64 only, `sbdot.c`, `tobf16.c`; left for T05).

### Level 3: BLIS's shared packed B

`gemm` follows BLIS's model (Van Zee and van de Geijn, "BLIS: A Framework for Rapidly
Instantiating BLAS Functionality", ACM TOMS 41(3), 2015, <https://doi.org/10.1145/2764454>; Smith
et al., "Anatomy of High-Performance Many-Threaded Matrix Multiplication", IPDPS 2014,
<https://doi.org/10.1109/IPDPS.2014.110>; `docs/Multithreading.md`), written from that design rather than ported
(`blas3/threaded.zig`). Its loops, packing, kernels and depth slices are the serial driver's.

- **Ways.** `w` workers split `C` `ic x jc` ways: the most workers that leave each part a row
  panel and a column panel, then the least packing, `n * ic + m * jc` (as BLIS's
  `bli_thread_partition_2x2` balances `m / ic` against `n / jc`), more `ic` on a tie. Row
  ranges start on `unroll_m`, column ranges on `unroll_n`. BLIS's `jr` and `ir` ways pay where
  threads share an L2, which `Io.Threaded` cannot arrange without affinity; not used.
- **Packing.** Each `jc` group owns one buffer for a depth slice of its columns of `op(B)`, at
  most `R` wide; its `ic` workers pack equal shares of the slice's column chunks into it. Each
  worker packs its own `P`-row blocks of `op(A)` into a buffer of its own.
- **Two barriers per slice**, as BLIS's `bli_packm_int.c`: one before packing (every peer is
  done with the last slice) and one after (the slice is whole). The first slice skips the first.
- **Fused packing.** A worker multiplies each of its own chunks with its first row block of
  `op(A)` as it packs it, before the second barrier, while the chunk is in cache; OpenBLAS's
  `level3_thread.c` does the same with its own share of `B`.
- **Start.** `threading.team` spawns `w - 1` workers through `Group.concurrent`; each waits
  until the caller has counted those that started, and the ways are sized for that count.
  None started, the caller runs alone. No worker waits for a peer that does not exist, so the
  call cannot deadlock when the `Io` admits fewer; `Group.async`'s inline fallback would.
- **Waiting.** BLIS and OpenBLAS both spin without bound inside a call: BLIS's
  `bli_thrcomm_barrier_atomic` on an empty loop, OpenBLAS's `level3_thread.c` with `YIELDING`
  (eight `nop`s on x86_64, arm, POWER and riscv64; `YieldProcessor` under MSVC;
  `sched_yield` or `SwitchToThread` elsewhere). Rstd spins too, checking the word around each
  `std.atomic.spinLoopHint` (`pause`, `isb`), but at most `spin_limit` (4096) times, about
  90 us on Zen 3, before it sleeps on `io.futexWait`; a waker calls `futexWake` only when
  someone sleeps (`threading.Signal`). Unbounded, a wait on a fiber `Io` (`Io.Evented` with
  peers on the waiter's thread) would never yield to the peer it waits for, and a canceled
  waiter would never reach a cancelable call. The bound outlasts the usual imbalance at a
  barrier, so a running team rarely sleeps, and stays well under a slice's work, so a fiber
  `Io` loses little per barrier.
- **Reproducibility.** The ways depend on the worker count, but every element of `C` is
  computed by one worker with the serial kernel calls in the serial order, and the depth
  slices (`Blocking.split`) do not depend on the split, so results match the serial path bit
  for bit.

`symm` and `hemm` run the same loops, their symmetric operand packed whole as in the serial
driver. The rank-k updates (`syrk`, `herk`, `syr2k`, `her2k`) keep one triangle of `C`, so
their `jc` groups are balanced by cost (`Shape.rising`, `falling`), as OpenBLAS's
`syrk_thread` balances its ranges. Tiles across the diagonal go through a scratch tile and
are added, which rounds differently from the kernel's own update, so a triangle's column cuts
fall where the serial driver's chunks start, and its rows follow the serial column block
(`drive`'s `r_lo`, `r_hi`); a full region's columns may split anywhere. `syr2k` runs its
second pass over the same split. `gemmt` stays serial, as upstream (`interface/gemmt.c`).

`trmm` and `trsm` split `B` into independent parts, as OpenBLAS's `gemm_thread_n` (left side,
column panels) and `gemm_thread_m` (right side, row panels) do, and as BLIS splits only its
`jc` and `ic` loops there: one part per worker, starting on `unroll_n` columns or `unroll_m`
rows, each a serial solve or product on packing space of its own. Nothing waits, so the parts
run through `forEach` (`Group.async`, inline when the `Io` admits fewer).

### Thresholds

| Routine | OpenBLAS threshold | Division |
| --- | --- | --- |
| real `gemm` | `MNK > 65536 * 4`; `t = min(ncpu, MNK / 262144)` | `ic x jc` ways, shared packed `B` |
| complex `gemm` | `MNK > 8192 * 4` | same |
| `syrk`, `symm` | `(n+1)nk` over `109944 * 4`; `2m^2 n` over `65536 * 4` | triangle-balanced ways; `gemm`'s |
| `syr2k` | `n * k >= 1000` | triangle-balanced ways |
| `trsm`, `trmm` | `m * n >= SMP_FACTOR (256 real, 128 cplx) * 4` | left: N panels; right: M panels |
| `gemv` | `m * n >= 115200 * 4` | N: 4096-row chunks; T: column chunks |
| `gbmv` | `m * n >= 250000` (125000 complex), `kl + ku >= 15` | T: column chunks; N serial |
| `ger` | `m*n > 2048*4`: up to 3; `> 20480*4`: all | column chunks |
| `symv`, `trmv` | `symv n >= 200`; `trmv n < 50`: 1, `< 500`: 2 | serial here (above) |
| `syr`, `spr`, `syr2`, `spr2` (real and complex symmetric) | `n < 50`: 1, `< 500`: 2 | sqrt-balanced column chunks |
| `her`, `hpr`, `her2`, `hpr2` | `n < 100`: 1 | same |
| `axpy` / `scal` / `swap` | `n > 10000` / `n > 1048576` / `n >= 8 MiB / sizeof(T)` | range chunks |
| `copy` | `copy_thunderx2t99.c`: `n > 10000` | range chunks |
| `dot`, `asum`, `sum`, `rot` | x86_64 kernels: `n > 10000` (dot), `> 100000` | range chunks, ordered combine |
| extremes, indices | `iamax_thunderx2t99.c`, `izamax_thunderx2t99.c`: `n > 10000` | range chunks, exact fold |

The `copy` and extremes rows are OpenBLAS's ThunderX2 kernels' own, taken as defaults for every
core; the others are its generic interfaces' or x86_64 kernels'. The generic thresholds apply
with `GEMM_MULTITHREAD_THRESHOLD = 4` as a comptime constant,
together with OpenBLAS's per-core overrides (`neoversev1`/`v2` in `interface/gemm.c` and
`gemv.c`, `a64fx` and `neoversev1` in `arm64/dot.c`). Complex routines have their own
constants; `threading.zig`'s `rule` holds them all. These are defaults; per-core tables derived from rstd's own benchmarks will
replace them, and the chunk sizes, later. `switch_ratio` already lives in the vtable
(`blas3/gemm/params.zig`). M ranges round to `unroll_m` and a cache line to avoid false
sharing in `C`. The rank-update cut-offs only pick the thread count; their serial small-`n`
paths are already ported ([`design.md`](design.md), Level 2).

### Memory

All scratch is allocated by the caller before any task starts, in one `alignedAlloc` sized
for `w` workers: level 3 packing buffers (`P x Q` elements per worker, and per `jc` group a
slice of `B` as wide as the widest group, at most `R`, each on a cache line; for `trmm` and
`trsm`, the serial driver's space per part) and `gemv`'s buffer per chunk. Tasks are therefore
infallible, and the allocator need not be thread safe. Level 1, `ger`, `gbmv` and the rank
updates need none beyond their serial copies of strided vectors, made before any task starts. Task
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
  say so. Level 3 workers waiting at a barrier are canceled with the group; the bounded spin
brings each to its cancelable futex wait.
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
([`design.md`](design.md#apple-accelerate)). Rstd's own workers are to parallelize around it.

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
- **The caller's control state.** Threads do not inherit floating-point control state (x86
  Windows starts them at 53-bit x87 precision, 0x27f, where a MinGW main thread has 64-bit,
  0x37f). Every task enters through `team`'s worker or `Region.work`; `team` and `forEachIn`
  capture `util.FpState` on the caller, and each task's first step adopts it, restoring its
  thread's own when it ends. Only control bits move: exception flags stay per thread, and the
  caller's state is never changed ([`design.md`](design.md), Settled decisions).

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
- **Disjoint writes.** Each output element and each partial has one writer; a shared packed
  buffer is written between two barriers and read only after the second.
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
| Independent level 3 tiles via `Group.async` | Redundant packing; BLIS's shared packed `B` is followed instead |
| OpenBLAS's `level3_thread.c` flags for level 3 | Prototyped: within a few percent of BLIS's barriers, at twice the code and `2 * tm` waits per slice |
| Two `B` buffers and one barrier per slice | BLIS meets twice around one buffer: half the `B` memory |
| Unbounded spinning, as BLIS and OpenBLAS | Stalls a fiber `Io` and cancellation |
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

## Level 3 measurements

Two prototypes (`wip/l3proto`) threaded real `gemm` before step 2, both over the serial
driver's packing, kernels and `Blocking.split` and both started through `threading.team`:

- **(A) OpenBLAS**, a port of `level3_thread.c`: the `CNAME` grid; each worker packs two halves
  of its share of its N group's `B`, publishes each with a step number, and its readers
  release it through a count. Per slice a worker waits up to `2 * (tm - 1)` times for halves
  and twice for releases. 172 lines.
- **(B) BLIS**, from `docs/Multithreading.md`: `jc x ic` ways; the `ic` workers of a `jc` group
  pack each slice of `B` together into one of two buffers and meet at one barrier per slice.
  93 lines.

Both matched `driver.gemm` bit for bit. Row ranges must start on `unroll_m`: unrounded ones
change the last bits (row tails narrower than a vector go through `thin`); columns may split
anywhere. (B), reworked to BLIS's one buffer and two barriers, fused and with bounded spinning,
is `blas3/threaded.zig`; (A) was dropped.

GFLOPS of the prototypes on this Zen 3 (8 cores, 16 threads, core `zen`, portable Zig tiles),
`C = A * B` square, `NN`, by workers; mean of two quiet runs of medians of 9 to 25 interleaved
calls:

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

(A) led by 0 to 5% at most points and by 10 to 17% at a few (4 or 16 workers, `n` of 1000 and
2000). Both scaled alike past 4 workers, so the limit there is shared (blocking, kernels,
memory), not the protocol.

`blas3/threaded.zig` as shipped, the same machine and bench (`l3thread bench`, built with
`zig test -O fast`, since `rstd-test` builds the tests `debug` whatever `-Doptimize` says),
2026-10-29; medians of four quiet runs, each a median of 9 to 25 interleaved calls:

| Type | n | serial | 2 | 4 | 8 | 16 |
| --- | --- | --- | --- | --- | --- | --- |
| `f64` | 256 | 36 | 65 | 91 | 96 | 74 |
| `f64` | 512 | 36 | 69 | 107 | 136 | 165 |
| `f64` | 1000 | 39 | 75 | 106 | 133 | 160 |
| `f64` | 2000 | 44 | 53 | 86 | 126 | 145 |
| `f32` | 256 | 68 | 91 | 113 | 113 | 89 |
| `f32` | 512 | 79 | 92 | 150 | 210 | 268 |
| `f32` | 1000 | 77 | 109 | 179 | 252 | 413 |
| `f32` | 2000 | 94 | 109 | 174 | 259 | 357 |

At 8 and 16 workers this is within 5% of (B) at most points (`f64` 1000 at 16: 160 against
155; `f32` 2000 at 16: 357 against 373), above it at 512, and below at `f64` 2000 at 8 (126
against 147). At 2 and 4 workers it is well below (B) here, but not in every run: two earlier
quiet runs that day had 2 and 4 workers 40 to 80% faster (`f64` 2000: 77 and 121) and 8 and 16
alike. The four runs agree with each other, so the difference is the machine's state, most
likely where Windows puts two to four threads (one core's SMT pair or two cores), which
`Io.Threaded` does not control. Fused packing against packing first, over those six runs, is a
wash to 4% ahead (geometric mean of the 32 cells per run).

## Implementation plan

1. **Signatures and plumbing** (done): `io: Io` on every public routine but the scalar
   `rotg`, `zrotg` and `rotmg`, with `error.Canceled` in their error sets; in `threading.zig`,
   the cached CPU count (`blas.config.cpu_count`) and ceiling (`Config.threadCeiling`),
   thresholds (`rule`), chunking (`Chunks`, `plan`), the claim loop (`forEach`) and ordered
   combine (`reduce`), with their unit tests.
2. **`gemm`** (done), every precision: BLIS's ways and shared packed `B` through
   `Group.concurrent`, with the start barrier, bounded spinning, the serial fallback and the
   public interface's threshold.
3. **Level 1** (done): `copy`, `swap`, `scal`, `rot`, `axpy`, the dots, `asum`, `sum` and the
   extremes; `nrm2` stays serial.
4. **Level 2** (done): `gemv`, `ger`, the rank updates and transposed `gbmv`; `symv`, `trmv`
   and the other banded and packed products stay serial.
5. **Remaining level 3** (done): `symm`, `hemm` and the rank-k updates over the same
   driver; `trmm` and `trsm` by panels of `B`.

## Open issues

- Level 3 on `Io.Evented`: `gemm` keeps its bits on `Uring` with every fiber on the caller's
  thread (`thread_limit = 0`, `qemu-test-linux-x86_64`), the case an unbounded spin would
  hang. With worker threads `Uring` hangs even a plain futex barrier ([`issues.md`](issues.md));
  `Dispatch` (macOS) has no thread limit and is untested. Windows has no `Evented`.
- TSan support with this toolchain is unverified.
- At 64-bit x87 precision (0x37f), the scalar `gemm` kernel's full tile keeps its sums in
  registers and its partial tile in memory, so they round `f64` apart; where a partition puts
  a column in the other one (`symm`, `pentium3`, the tests' shrunk `R`), threaded bits differ.
  At 53-bit they agree. T04's x87 audit.
- Thresholds and chunk sizes are OpenBLAS's defaults until rstd's own per-core benchmarks
  exist.
- `symv`, `trmv` and the banded and packed products could split `y` keeping their bits if every
  core's `axpy` (and `gemv`) gave an element the same bits wherever it sits in a call; untested
  across cores, so they stay serial. Folding fixed per-size partial vectors instead (OpenBLAS's
  way, with chunks in place of threads) would keep the rule too, but change their serial bits
  and cost `n` per chunk.
