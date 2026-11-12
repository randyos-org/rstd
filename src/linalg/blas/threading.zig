//! The machinery of `threading.md`'s levels 1 and 2: thresholds, the worker
//! count, fixed chunks, the claim loop and the ordered combine. Chunks follow
//! the problem size and the core only, never the worker count, and partials
//! fold in chunk order, so the bits are the same for any worker count, `Io` or
//! scheduling. Level 3 uses the thresholds and the worker count.

const std = @import("std");
const sysinfo = @import("builtin");
const blas = @import("__root__.zig");
const types = @import("types.zig");
const Config = @import("Config.zig");
const FpState = @import("FpState.zig");
const util = @import("util.zig");
const test_io = @import("_test/io.zig");
const zigstd = @import("../../zigstd.zig");
const coretypes = @import("cpu/coretypes.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;

/// The most chunks per call; a reduction keeps one partial each on the stack.
pub const chunks_max = 64;

/// OpenBLAS's `GEMM_MULTITHREAD_THRESHOLD` default.
const mt = 4;

// ---- thresholds ------------------------------------------------------------------

/// A routine family, as OpenBLAS thresholds it, with the `work` its rule reads.
/// OpenBLAS runs a zero-stride level-1 call serially; so must the routine.
pub const Family = enum {
    /// `axpy`, `axpby`: `n`.
    axpy,
    /// `scal`: `n`.
    scal,
    /// `swap`: `n`.
    swap,
    /// `copy`: `n`.
    copy,
    /// `rot`, `rotm`: `n`.
    rot,
    /// `dot`, `dotc` and the mixed dots: `n`.
    dot,
    /// `asum`, `sum`: `n`.
    asum,
    /// `index`, `extreme`: `n`.
    index,
    /// `gemv`: `m * n`.
    gemv,
    /// `gbmv`: `m * n`.
    gbmv,
    /// `ger`: `m * n`.
    ger,
    /// `symv`, `spmv`, `hpmv`: `n`.
    symv,
    /// `hemv`: `n`.
    hemv,
    /// `trmv`: `n`.
    trmv,
    /// `tbmv`, `tpmv`: `n`.
    tpmv,
    /// `syr`, `spr`, `syr2`, `spr2`, real and complex symmetric: `n`.
    syr,
    /// `her`, `hpr`, `her2`, `hpr2`: `n`.
    her,
    /// `gemm`: `m * n * k`.
    gemm,
    /// `syrk`, `herk`: `(n + 1) * n * k`.
    syrk,
    /// `syr2k`, `her2k`: `n * k`.
    syr2k,
    /// `symm`: `2 * m * m * n`.
    symm,
    /// `trsm`: `m * n`.
    trsm,
};

/// Work below `below` takes at most `max` workers.
pub const Tier = struct { below: u64, max: u32 };

/// Work up to `x` takes at most `max` workers.
fn upTo(x: u64, max: u32) Tier {
    return .{ .below = x + 1, .max = max };
}

/// A family's threshold for one element type and core.
pub const Rule = struct {
    /// The most work run serially, as one chunk.
    serial_max: u64,
    /// Past the tiers, at most `work / serial_max` workers rather than any.
    per_unit: bool = false,
    /// Caps above `serial_max`, by ascending `below`.
    tiers: []const Tier = &.{},

    /// The most workers `work` justifies.
    pub fn sizeMax(r: Rule, work: u64) u32 {
        if (work <= r.serial_max) return 1;
        for (r.tiers) |t| if (work < t.below) return t.max;
        if (!r.per_unit) return std.math.maxInt(u32);
        return std.math.cast(u32, work / @max(r.serial_max, 1)) orelse std.math.maxInt(u32);
    }

    /// One chunk up to `serial_max`; above, chunks of about `serial_max`, at
    /// most `chunks_max`.
    pub fn chunkCount(r: Rule, work: u64) usize {
        if (work <= r.serial_max) return 1;
        const per = @max(r.serial_max, 1);
        return @intCast(@min(work / per + @intFromBool(work % per != 0), chunks_max));
    }
};

/// `family`'s rule for `T` on `core` (`blas.coreType()`; any architecture's
/// core enum, or null for none): OpenBLAS's
/// interface or kernel thresholds, with its per-core overrides as a
/// `DYNAMIC_ARCH` build applies them.
pub fn rule(comptime family: Family, comptime T: type, core: anytype) Rule {
    const complex = comptime types.isComplex(T);
    const F = types.Float(T);
    const s: u64 = @sizeOf(F);
    return switch (family) {
        .axpy => .{ .serial_max = 10000 },
        .scal => .{ .serial_max = 1048576 },
        .swap => .{ .serial_max = @as(u64, if (complex) 1048576 else 2097152) * mt / s - 1 },
        // `kernel/arm64/copy_thunderx2t99.c`: `n <= 10000`
        .copy => .{ .serial_max = 10000 },
        .rot => .{ .serial_max = 100000 },
        .dot => (if (complex) null else dotCore(F, core)) orelse .{ .serial_max = 10000 },
        .asum => .{ .serial_max = if (complex) 10000 else 100000, .per_unit = true },
        // `kernel/arm64/iamax_thunderx2t99.c`, `izamax_thunderx2t99.c`: `n <= 10000`
        .index => .{ .serial_max = 10000 },
        .gemv => if (complex) .{ .serial_max = 1024 * mt - 1 } else gemvCore(F, core) orelse .{ .serial_max = 115200 * mt - 1 },
        // `gbmv.c`, `zgbmv.c`: `m * n < 250000` (125000 complex); their `kl + ku < 15` is the interface's
        .gbmv => .{ .serial_max = if (complex) 124999 else 249999 },
        .ger => if (complex) .{ .serial_max = 36 * s * s * mt } else .{ .serial_max = 2048 * mt, .tiers = comptime &.{upTo(20480 * mt, 3)} },
        // `ztrmv`'s `n*n > 36*s*s*4` and `n*n < 64*s*s*4`, as bounds on `n`
        .trmv => if (complex) .{ .serial_max = 12 * s, .tiers = comptime &.{.{ .below = 16 * s, .max = 2 }} } else .{ .serial_max = 49, .tiers = comptime &.{.{ .below = 500, .max = 2 }} },
        // `syr.c`, `spr.c`, `syr2.c`, `spr2.c`, `zsyr.c`, `zspr.c`: `n < 50`, `n < 500`
        .syr => .{ .serial_max = 49, .tiers = comptime &.{.{ .below = 500, .max = 2 }} },
        // `zher.c`, `zhpr.c`, `zher2.c`, `zhpr2.c`: `n < 100`
        .her => .{ .serial_max = 99 },
        // `symv.c`: `n < 200`; `zsymv.c`, `zhpmv.c`: `n < 50`, `n < 500`
        .symv => if (complex) .{ .serial_max = 49, .tiers = comptime &.{.{ .below = 500, .max = 2 }} } else .{ .serial_max = 199 },
        // `zhemv.c`: `n < MULTI_THREAD_MINIMAL`, 362
        .hemv => .{ .serial_max = 361 },
        // `tbmv.c`, `tpmv.c`, `ztbmv.c`, `ztpmv.c`: `n < 50`, `n < 500`
        .tpmv => .{ .serial_max = 49, .tiers = comptime &.{.{ .below = 500, .max = 2 }} },
        .gemm => (if (complex) null else gemmCore(F, core)) orelse .{ .serial_max = @as(u64, if (complex) 8192 else 65536) * mt, .per_unit = true },
        .syrk => .{ .serial_max = @as(u64, if (complex) 14824 else 109944) * mt },
        // `syr2k.c`: `n * k < 1000`
        .syr2k => .{ .serial_max = 999 },
        .symm => .{ .serial_max = @as(u64, if (complex) 8192 else 65536) * mt },
        .trsm => .{ .serial_max = @as(u64, if (complex) 128 else 256) * mt - 1 },
    };
}

/// Whether `core` (`rule`'s) is arm64's `c`.
fn is(core: anytype, comptime c: coretypes.Arm64CoreType) bool {
    const C = @TypeOf(core);
    if (C == coretypes.Arm64CoreType) return core == c;
    if (C == ?coretypes.Arm64CoreType) return if (core) |k| k == c else false;
    return false;
}

/// `kernel/arm64/dot.c`'s tables.
fn dotCore(comptime F: type, core: anytype) ?Rule {
    if (is(core, .neoversev1)) return switch (F) {
        f64 => .{ .serial_max = 64500, .tiers = comptime &.{ upTo(100000, 2), upTo(150000, 4), upTo(260000, 8), upTo(360000, 16), upTo(520000, 24), upTo(1010000, 56) } },
        f32 => .{ .serial_max = 110000, .tiers = comptime &.{ upTo(200000, 2), upTo(280000, 4), upTo(520000, 8), upTo(830000, 16), upTo(1010000, 24) } },
        else => null,
    };
    if (is(core, .a64fx)) return switch (F) {
        f64 => .{ .serial_max = 11000, .tiers = comptime &.{ upTo(20000, 2), upTo(35000, 4), upTo(50000, 6), upTo(440000, 8), upTo(880000, 16), upTo(1020000, 24) } },
        f32 => .{ .serial_max = 22000, .tiers = comptime &.{ upTo(39000, 2), upTo(79000, 4), upTo(120000, 6), upTo(1020000, 8) } },
        else => null,
    };
    return null;
}

/// `interface/gemv.c`'s tables.
fn gemvCore(comptime F: type, core: anytype) ?Rule {
    if (F != f32 and F != f64) return null;
    if (is(core, .neoversev1)) return switch (F) {
        f64 => .{ .serial_max = 8099, .tiers = comptime &.{ .{ .below = 12100, .max = 2 }, .{ .below = 36100, .max = 4 }, .{ .below = 84100, .max = 8 }, .{ .below = 348100, .max = 16 }, .{ .below = 435600, .max = 24 }, .{ .below = 810000, .max = 32 }, .{ .below = 1050625, .max = 40 } } },
        else => .{ .serial_max = 25599, .tiers = comptime &.{ .{ .below = 63001, .max = 4 }, .{ .below = 202500, .max = 8 }, .{ .below = 806404, .max = 16 }, .{ .below = 1638400, .max = 32 } } },
    };
    if (is(core, .neoversev2)) return .{ .serial_max = 24963, .tiers = comptime &.{ .{ .below = 145924, .max = 8 }, .{ .below = 692224, .max = 16 }, .{ .below = 1638400, .max = 32 } } };
    return null;
}

/// `interface/gemm.c`'s tables, single precision only.
fn gemmCore(comptime F: type, core: anytype) ?Rule {
    if (F != f32) return null;
    if (is(core, .neoversev1)) return .{ .serial_max = 262143, .tiers = comptime &.{ .{ .below = 1124864, .max = 6 }, .{ .below = 7880599, .max = 12 }, .{ .below = 17173512, .max = 16 }, .{ .below = 33386248, .max = 20 }, .{ .below = 57066625, .max = 24 }, .{ .below = 91733851, .max = 32 }, .{ .below = 265847707, .max = 40 }, .{ .below = 458314011, .max = 48 }, .{ .below = 729000000, .max = 56 } } };
    if (is(core, .neoversev2)) return .{ .serial_max = 124999, .tiers = comptime &.{ .{ .below = 1092727, .max = 6 }, .{ .below = 2628072, .max = 8 }, .{ .below = 8000000, .max = 12 }, .{ .below = 20346417, .max = 16 }, .{ .below = 57066625, .max = 24 }, .{ .below = 91125000, .max = 28 }, .{ .below = 238328000, .max = 40 }, .{ .below = 454756609, .max = 48 }, .{ .below = 857375000, .max = 56 }, .{ .below = 1073741824, .max = 64 } } };
    return null;
}

/// `dims` multiplied, saturating: a rule's `work`.
pub fn product(dims: []const u64) u64 {
    var p: u64 = 1;
    for (dims) |d| p *|= d;
    return p;
}

// ---- worker count ----------------------------------------------------------------

/// Whether `std.Thread.getCpuCount` exists here: std's own threads, or the
/// root's override where std has none (Zig 0.17 has no overrides).
const has_cpu_count = !sysinfo.single_threaded and (thread_override or switch (sysinfo.target.os.tag) {
    .windows, .linux, .wasi => true,
    else => std.Thread.use_pthreads,
});

const thread_override = if (@hasDecl(std, "os_options")) std.os_options.Thread != null else false;

/// `std.Thread.getCpuCount() catch 1`; 1 without threads. `blas.init` caches it.
pub fn cpuCount() u32 {
    if (comptime !has_cpu_count) return 1;
    const n = std.Thread.getCpuCount() catch return 1;
    return std.math.cast(u32, n) orelse std.math.maxInt(u32);
}

/// `min(cpu_count, ceiling, size_max)`, at least 1; 1 in `single_threaded` builds.
pub fn workers(cfg: *const Config, size_max: u32) u32 {
    if (comptime sysinfo.single_threaded) return 1;
    var w = @min(cfg.cpu_count, size_max);
    if (cfg.threadCeiling()) |c| w = @min(w, c);
    return @max(w, 1);
}

// ---- chunks ------------------------------------------------------------------------

/// How cost spreads over a range's units.
pub const Shape = enum {
    /// Every unit costs the same.
    even,
    /// Unit `i` costs about `i`, as a triangle's growing rows.
    rising,
    /// Unit `i` costs about `len - i`.
    falling,
};

/// Units `start..end` of a range, chunk `index` of its split.
pub const Chunk = struct { index: usize, start: usize, end: usize };

/// A range split into balanced, non-empty chunks by its length, the chunk
/// count asked for, shape and quantum only.
pub const Chunks = struct {
    count: usize,
    bounds: [chunks_max + 1]usize,

    /// `want` chunks (clamped to 1 and `chunks_max`) of `len` units, inner
    /// bounds rounded down to multiples of `quantum`; empty ones dropped. An
    /// empty range is one empty chunk.
    pub fn init(len: usize, want: usize, shape: Shape, quantum: usize) Chunks {
        const q = @max(quantum, 1);
        const c = std.math.clamp(want, 1, chunks_max);
        var r: Chunks = .{ .count = 0, .bounds = undefined };
        r.bounds[0] = 0;
        for (1..c + 1) |k| {
            const b = if (k == c) len else blk: {
                const x = cut(len, c, k, shape);
                break :blk x - x % q;
            };
            if (b > r.bounds[r.count]) {
                r.count += 1;
                r.bounds[r.count] = b;
            }
        }
        if (r.count == 0) {
            r.count = 1;
            r.bounds[1] = 0;
        }
        return r;
    }

    pub fn get(r: *const Chunks, i: usize) Chunk {
        return .{ .index = i, .start = r.bounds[i], .end = r.bounds[i + 1] };
    }

    /// Bound `k` of `c` before rounding: `k/c` of the cost.
    fn cut(len: usize, c: usize, k: usize, shape: Shape) usize {
        return switch (shape) {
            .even => @intCast(fraction(len, c, k)),
            .rising => @intCast(std.math.sqrt(fraction(@as(u128, len) * len, c, k))),
            .falling => len - @as(usize, @intCast(std.math.sqrt(fraction(@as(u128, len) * len, c, c - k)))),
        };
    }

    /// `floor(x * k / c)` for `k <= c`, without overflow.
    fn fraction(x: u128, c: usize, k: usize) u128 {
        return (x / c) * k + (x % c) * k / c;
    }
};

/// A call's chunks and the most workers that run them.
pub const Plan = struct {
    chunks: Chunks,
    workers: u32,
};

/// The least `serial_max` of `family`'s rule for `T` on any core.
fn serialFloor(comptime family: Family, comptime T: type) u64 {
    comptime {
        @setEvalBranchQuota(1_000_000);
        var m = rule(family, T, null).serial_max;
        for (coretypes.all) |C| {
            for (std.meta.tags(C)) |c| m = @min(m, rule(family, T, c).serial_max);
        }
        return m;
    }
}

/// Whether `work` is one chunk on one worker on every core: `plan` then needs
/// neither the core nor the config, which small calls would pay for (T100).
pub fn serial(comptime family: Family, comptime T: type, work: u64) bool {
    if (sysinfo.is_test and test_serial_max != null) return false;
    return work <= comptime serialFloor(family, T);
}

/// `family`'s plan for `work` over `len` units, from `blas.config` and the
/// selected core.
pub fn plan(comptime family: Family, comptime T: type, work: u64, len: usize, shape: Shape, quantum: usize) Plan {
    if (serial(family, T, work)) return .{ .chunks = .init(len, 1, shape, quantum), .workers = 1 };
    // selects the core first, which fills `blas.config`
    var r = rule(family, T, blas.coreType());
    if (sysinfo.is_test) if (test_serial_max) |m| {
        r = .{ .serial_max = m };
    };
    return planWith(&blas.config, r, work, len, shape, quantum);
}

/// Tests set this to replace every rule of `plan` by one with this
/// `serial_max`, so small problems make many chunks.
pub var test_serial_max: ?u64 = null;

/// `plan` from `cfg` and `r`; only `workers` reads `cfg`.
pub fn planWith(cfg: *const Config, r: Rule, work: u64, len: usize, shape: Shape, quantum: usize) Plan {
    const chunks: Chunks = .init(len, r.chunkCount(work), shape, quantum);
    return .{ .chunks = chunks, .workers = workers(cfg, @min(r.sizeMax(work), @as(u32, @intCast(chunks.count)))) };
}

// ---- the claim loop and the ordered combine ---------------------------------------

/// `f(ctx, chunk)` for every chunk of `p`, on up to `p.workers` workers that
/// claim chunks in turn; the caller is worker 0, and the rest start through
/// `Group.async` (`Tasks`), so an `Io` that admits fewer runs them inline,
/// and adopt the caller's floating-point control state (`FpState.zig`).
/// Returns once every task has finished. `error.Canceled` leaves chunks undone.
pub fn forEach(io: Io, p: *const Plan, ctx: anytype, comptime f: fn (@TypeOf(ctx), Chunk) void) Io.Cancelable!void {
    return forEachIn(io, p, ctx, f, null);
}

/// `forEach` for a reduction: `partial(ctx, chunk)` per chunk, folded in chunk
/// order on the caller, so the result does not depend on who computed what.
pub fn reduce(
    io: Io,
    p: *const Plan,
    comptime R: type,
    ctx: anytype,
    comptime partial: fn (@TypeOf(ctx), Chunk) R,
    comptime fold: fn (R, R) R,
) Io.Cancelable!R {
    return reduceIn(io, p, R, ctx, partial, fold, null);
}

/// `forEach` with scratch: one allocation, cut into `len(ctx, chunk)`
/// elements per chunk in chunk order; `f` gets its chunk's part.
pub fn forEachPart(
    allocator: Allocator,
    io: Io,
    p: *const Plan,
    comptime E: type,
    comptime alignment: usize,
    ctx: anytype,
    comptime len: fn (@TypeOf(ctx), Chunk) usize,
    comptime f: fn (@TypeOf(ctx), Chunk, []E) void,
) (Allocator.Error || Io.Cancelable)!void {
    var at: [chunks_max + 1]usize = undefined;
    at[0] = 0;
    for (0..p.chunks.count) |i| at[i + 1] = at[i] + len(ctx, p.chunks.get(i));
    const buf = try allocator.alignedAlloc(E, .fromByteUnits(alignment), at[p.chunks.count]);
    defer allocator.free(buf);
    const Parts = struct {
        ctx: @TypeOf(ctx),
        buf: [*]E,
        at: *const [chunks_max + 1]usize,

        fn one(s: @This(), c: Chunk) void {
            f(s.ctx, c, s.buf[s.at[c.index]..s.at[c.index + 1]]);
        }
    };
    return forEach(io, p, Parts{ .ctx = ctx, .buf = buf.ptr, .at = &at }, Parts.one);
}

/// `forEach`, claiming chunks in `order` if given (tests).
fn forEachIn(io: Io, p: *const Plan, ctx: anytype, comptime f: fn (@TypeOf(ctx), Chunk) void, order: ?[]const u8) Io.Cancelable!void {
    var region: Region(@TypeOf(ctx), f) = .{ .ctx = ctx, .chunks = &p.chunks, .order = order };
    if (comptime !sysinfo.single_threaded) {
        const w = @min(p.workers, p.chunks.count);
        if (w > 1) return region.parallel(io, w);
    }
    return region.claim(io);
}

fn reduceIn(
    io: Io,
    p: *const Plan,
    comptime R: type,
    ctx: anytype,
    comptime partial: fn (@TypeOf(ctx), Chunk) R,
    comptime fold: fn (R, R) R,
    order: ?[]const u8,
) Io.Cancelable!R {
    const Into = struct {
        ctx: @TypeOf(ctx),
        partials: *[chunks_max]R,

        fn one(into: @This(), c: Chunk) void {
            into.partials[c.index] = partial(into.ctx, c);
        }
    };
    var partials: [chunks_max]R = undefined;
    try forEachIn(io, p, Into{ .ctx = ctx, .partials = &partials }, Into.one, order);
    var acc = partials[0];
    for (partials[1..p.chunks.count]) |x| acc = fold(acc, x);
    return acc;
}

/// One parallel region's shared state; tasks get a pointer to it.
fn Region(comptime Ctx: type, comptime f: fn (Ctx, Chunk) void) type {
    return struct {
        const Self = @This();

        ctx: Ctx,
        chunks: *const Chunks,
        order: ?[]const u8,
        next: std.atomic.Value(usize) = .init(0),
        /// The caller's, which `work` adopts.
        fp: FpState = undefined,

        fn claim(r: *Self, io: Io) Io.Cancelable!void {
            while (true) {
                try checkCancel(io);
                const k = r.next.fetchAdd(1, .monotonic);
                if (k >= r.chunks.count) return;
                f(r.ctx, r.chunks.get(if (r.order) |o| o[k] else k));
            }
        }

        /// A worker task. It re-arms a cancellation it took: a task run inline
        /// shares the caller's, and `Group` drops the error.
        fn work(r: *Self, io: Io) void {
            const own = r.fp.adopt();
            defer if (own) |s| s.apply();
            r.claim(io) catch |err| switch (err) {
                error.Canceled => io.recancel(),
            };
        }

        fn parallel(r: *Self, io: Io, w: usize) Io.Cancelable!void {
            r.fp = .capture();
            var tasks: Tasks = .init(io);
            // runs before the caller's frame goes; every task has then run
            defer {
                tasks.cancel(io);
                // the caller's state where a fiber `Io` resumed it (`Lent`)
                _ = r.fp.adopt();
            }
            for (1..w) |_| tasks.async(io, Self.work, .{ r, io });
            try r.claim(io);
            try tasks.await(io);
        }
    };
}

// ---- worker tasks ------------------------------------------------------------------

/// Plain futures instead of `Io.Group` on `Io.Dispatch`, whose `Group` crashes
/// with worker threads (issue #208). Once `ssh-repro-group-aarch64-macos`
/// passes, set this false and delete `Futures`.
pub const dispatch_futures = true;

/// `Io.Evented` for the fiber tests; void where it does not build: Zig 0.17's
/// `Uring`, and its `Dispatch`, whose vtable names a field 0.17's `Io.VTable`
/// lacks (version-gated, std exposes nothing to detect); SVE targets, whose
/// LLVM rejects the fiber switch's `ffr` clobber (issue #212).
pub const Evented = if ((Io.Evented == Io.Uring or Io.Evented == Io.Dispatch) and sysinfo.zig_version.major == 0 and sysinfo.zig_version.minor == 17) void else if (sysinfo.cpu.arch.isAarch64() and std.Target.aarch64.featureSetHas(sysinfo.cpu.features, .sve)) void else Io.Evented;

/// Whether `io`'s worker tasks run on `Futures`.
pub fn usesFutures(io: Io) bool {
    // `Evented`, not `Io.Evented`: Zig 0.17.0 cannot instantiate `Dispatch`
    if (comptime !dispatch_futures or Evented != Io.Dispatch) return false;
    const S = struct {
        var ev: Io.Dispatch = undefined;
    };
    // `Dispatch.io` returns one constant vtable
    return io.vtable == S.ev.io().vtable;
}

/// A region's worker tasks: an `Io.Group`, or `Futures` (`usesFutures`).
const Tasks = struct {
    group: Io.Group = .init,
    futures: Futures = .{},
    on_futures: bool,

    fn init(io: Io) Tasks {
        return .{ .on_futures = usesFutures(io) };
    }

    fn async(t: *Tasks, io: Io, function: anytype, args: std.meta.ArgsTuple(@TypeOf(function))) void {
        if (t.on_futures) return t.futures.async(io, function, args);
        t.group.async(io, function, args);
    }

    fn concurrent(t: *Tasks, io: Io, function: anytype, args: std.meta.ArgsTuple(@TypeOf(function))) Io.ConcurrentError!void {
        if (t.on_futures) return t.futures.concurrent(io, function, args);
        return t.group.concurrent(io, function, args);
    }

    fn await(t: *Tasks, io: Io) Io.Cancelable!void {
        if (t.on_futures) return t.futures.await(io);
        return t.group.await(io);
    }

    fn cancel(t: *Tasks, io: Io) void {
        if (t.on_futures) return t.futures.cancel(io);
        t.group.cancel(io);
    }
};

/// Tasks as plain futures, awaited and canceled as `Io.Group`'s are: `await`
/// waits cancelably for them all, then `cancel` (the caller's `defer`)
/// cancels and awaits each.
const Futures = struct {
    list: [cap]Io.Future(void) = undefined,
    len: usize = 0,
    /// Tasks not yet finished, plus one the caller holds until `await`.
    left: std.atomic.Value(u32) = .init(1),
    done: Signal = .{},

    const cap = if (dispatch_futures and Evented == Io.Dispatch) chunks_max else 0;

    fn Task(comptime function: anytype) type {
        return struct {
            fn run(fs: *Futures, io: Io, args: std.meta.ArgsTuple(@TypeOf(function))) void {
                @call(.auto, function, args);
                if (fs.left.fetchSub(1, .acq_rel) == 1) fs.done.set(io, 1, null);
            }
        };
    }

    /// Past `cap`, inline, as `Group.async` may run it.
    fn async(fs: *Futures, io: Io, function: anytype, args: std.meta.ArgsTuple(@TypeOf(function))) void {
        if (cap == 0 or fs.len == cap) return @call(.auto, function, args);
        _ = fs.left.fetchAdd(1, .monotonic);
        fs.list[fs.len] = io.async(Task(function).run, .{ fs, io, args });
        fs.len += 1;
    }

    fn concurrent(fs: *Futures, io: Io, function: anytype, args: std.meta.ArgsTuple(@TypeOf(function))) Io.ConcurrentError!void {
        if (cap == 0 or fs.len == cap) return error.ConcurrencyUnavailable;
        _ = fs.left.fetchAdd(1, .monotonic);
        fs.list[fs.len] = io.concurrent(Task(function).run, .{ fs, io, args }) catch |err| {
            _ = fs.left.fetchSub(1, .monotonic);
            return err;
        };
        fs.len += 1;
    }

    fn await(fs: *Futures, io: Io) Io.Cancelable!void {
        if (fs.left.fetchSub(1, .acq_rel) != 1) _ = try fs.done.waitWhile(io, 0, null);
        for (fs.list[0..fs.len]) |*f| f.await(io);
        fs.len = 0;
    }

    fn cancel(fs: *Futures, io: Io) void {
        for (fs.list[0..fs.len]) |*f| f.cancel(io);
        fs.len = 0;
    }
};

/// `io.checkCancel()`, or nothing on an `Io` without cancellation, whose
/// `checkCancel` is unreachable (`Io.failing`, rstd's UEFI `Io`).
pub fn checkCancel(io: Io) Io.Cancelable!void {
    if (io.vtable.checkCancel == &Io.unreachableCheckCancel) return;
    return io.checkCancel();
}

// ---- partial vectors ---------------------------------------------------------------

/// Rows `lo..hi` of a vector.
pub const Rows = struct { lo: usize, hi: usize };

/// Per-chunk partial vectors, for a level-2 product whose elements of `y` sum
/// over several chunks' columns: chunk `i` adds its part into a zeroed vector
/// of its own over the rows it reaches, and `fold` adds those to `y` in chunk
/// order. Chunks follow the size alone, so the bits do too, at any worker
/// count (OpenBLAS's `symv_thread.c` keeps one vector per thread instead).
pub fn Partials(comptime T: type) type {
    return struct {
        const Self = @This();
        /// Each vector starts on a cache line.
        const line = 64 / std.math.gcd(64, @sizeOf(T));

        buf: []align(64) T,
        rows: [chunks_max]Rows,
        at: [chunks_max]usize,
        count: usize,
        /// The rows `fold` reaches: the highest `hi`.
        len: usize,

        /// Room for `p`'s chunks, chunk `c` over `rows(ctx, c)`.
        pub fn init(allocator: Allocator, p: *const Plan, ctx: anytype, comptime rows: fn (@TypeOf(ctx), Chunk) Rows) Allocator.Error!Self {
            var s: Self = .{ .buf = &.{}, .rows = undefined, .at = undefined, .count = p.chunks.count, .len = 0 };
            var total: usize = 0;
            for (0..s.count) |i| {
                const r = rows(ctx, p.chunks.get(i));
                s.rows[i] = r;
                s.at[i] = total;
                s.len = @max(s.len, r.hi);
                total += std.mem.alignForward(usize, r.hi - r.lo, line);
            }
            s.buf = try allocator.alignedAlloc(T, .@"64", total);
            return s;
        }

        pub fn deinit(s: *const Self, allocator: Allocator) void {
            allocator.free(s.buf);
        }

        /// Chunk `i`'s vector, zeroed; element `r - rows[i].lo` is row `r`.
        pub fn take(s: *const Self, i: usize) []T {
            const v = s.buf[s.at[i]..][0 .. s.rows[i].hi - s.rows[i].lo];
            @memset(v, std.mem.zeroes(T));
            return v;
        }

        /// `y[r * incy] += ` each chunk's row `r`, chunks in order, each sum
        /// rounded to its type; rows split over workers, which cannot change
        /// a bit.
        pub fn fold(s: *const Self, io: Io, y: [*]T, incy: isize) Io.Cancelable!void {
            var work: u64 = 0;
            for (s.rows[0..s.count]) |r| work += r.hi - r.lo;
            // split as `axpy` is
            const p = plan(.axpy, T, work, s.len, .even, line);
            return forEach(io, &p, Fold{ .s = s, .y = y, .incy = incy }, Fold.rows);
        }

        const Fold = struct {
            s: *const Self,
            y: [*]T,
            incy: isize,

            fn rows(f: Fold, c: Chunk) void {
                for (f.s.rows[0..f.s.count], f.s.at[0..f.s.count]) |r, at| {
                    const lo = @max(r.lo, c.start);
                    const hi = @min(r.hi, c.end);
                    if (lo >= hi) continue;
                    for (lo..hi) |row| {
                        const d = util.at(T, f.y, @as(i64, @intCast(row)) * f.incy);
                        d.* = add(d.*, f.s.buf[at + row - r.lo]);
                    }
                }
            }
        };

        fn add(a: T, b: T) T {
            if (comptime types.isComplex(T)) {
                const F = types.Float(T);
                return .init(util.rounded(F, a.re + b.re), util.rounded(F, a.im + b.im));
            }
            return util.rounded(T, a + b);
        }
    };
}

// ---- level 3: workers that wait on each other ------------------------------------

/// Runs `job` on up to `want` workers that wait on each other, the caller
/// as worker 0: `job.plan(n)` sizes the work for the `n` that started and
/// returns how many take part, `job.work(io, i, lent)` runs worker `i`.
/// Workers hold the caller's control state (`Lent`); returns when all end.
pub fn team(io: Io, job: anytype, want: u32) Io.Cancelable!void {
    if (comptime sysinfo.single_threaded) {
        _ = job.plan(1);
        return job.work(io, 0, null);
    }
    const Job = @TypeOf(job.*);
    const S = struct {
        fn worker(j: *Job, go: *Signal, i: u32, fp: FpState, wio: Io) void {
            var lent: Lent = .{ .want = fp };
            lent.take();
            defer lent.give();
            const active = go.waitWhile(wio, 0, &lent) catch return wio.recancel();
            if (i < active) j.work(wio, i, &lent) catch wio.recancel();
        }
    };
    const fp: FpState = .capture();
    // the caller's own: a fiber `Io` may resume it on another thread too
    var mine: Lent = .{ .want = fp };
    var go: Signal = .{};
    var tasks: Tasks = .init(io);
    defer {
        // first, canceling workers still waiting on `go`
        tasks.cancel(io);
        // the caller keeps its state where it resumed
        mine.take();
    }
    var started: u32 = 1;
    while (started < want) : (started += 1) tasks.concurrent(io, S.worker, .{ job, &go, started, fp, io }) catch break;
    go.set(io, job.plan(started), &mine);
    try job.work(io, 0, &mine);
    try tasks.await(io);
}

/// A worker's hold on the caller's control state (`FpState`), in force on
/// its thread only while it runs: a fiber `Io` (`Dispatch`) may resume a
/// waiter on another thread, which must not keep it, nor the first its own.
/// The caller's starts given (its thread's own is `want`).
pub const Lent = struct {
    want: FpState,
    /// The thread's own state while `want` replaces it.
    prior: ?FpState = null,

    /// `want` in force on this thread.
    pub fn take(l: *Lent) void {
        l.prior = l.want.adopt();
    }

    /// The thread's own state back, before a wait that may move the worker.
    pub fn give(l: *Lent) void {
        if (l.prior) |p| p.apply();
        l.prior = null;
    }
};

/// Spins of a wait before it sleeps; its time per spin (`pause` on Zen 3)
/// is unmeasured, so the wait's length is too. BLIS and OpenBLAS spin without
/// bound; the bound lets a fiber `Io` run the peer a waiter waits on, and a
/// canceled waiter reach its futex.
pub const spin_limit = 1 << 12;

/// A word workers wait on to change. They spin first, as BLIS's
/// `bli_thrcomm_barrier_atomic` and OpenBLAS's `level3_thread.c` do, then
/// sleep on a futex; setters wake only when someone sleeps. A worker sleeps
/// and wakes others without its `Lent` state, and takes it again after.
pub const Signal = struct {
    value: std.atomic.Value(u32) = .init(0),
    sleepers: std.atomic.Value(u32) = .init(0),

    /// Waits until the value is other than `old`; returns it.
    pub fn waitWhile(s: *Signal, io: Io, old: u32, lent: ?*Lent) Io.Cancelable!u32 {
        for (0..spin_limit) |_| {
            const x = s.value.load(.acquire);
            if (x != old) return x;
            std.atomic.spinLoopHint();
        }
        while (true) {
            // with `set`'s order, either it sees this sleeper or this sees its value
            _ = s.sleepers.fetchAdd(1, .seq_cst);
            defer _ = s.sleepers.fetchSub(1, .monotonic);
            const x = s.value.load(.seq_cst);
            if (x != old) return x;
            if (lent) |l| l.give();
            defer if (lent) |l| l.take();
            try io.futexWait(u32, &s.value.raw, old);
        }
    }

    /// Sets the value; a wake may move the setter too (`Dispatch`).
    pub fn set(s: *Signal, io: Io, v: u32, lent: ?*Lent) void {
        s.value.store(v, .seq_cst);
        if (s.sleepers.load(.seq_cst) == 0) return;
        if (lent) |l| l.give();
        defer if (lent) |l| l.take();
        // all: `Uring` passes the count on as a C `int`, where `maxInt(u32)`
        // wakes one (issue #207)
        io.futexWake(u32, &s.value.raw, std.math.maxInt(i32));
    }
};

/// A reusable sense-reversing barrier of `parties` workers
/// (`bli_thrcomm_barrier_atomic`), waiting as `Signal` does.
pub const Barrier = struct {
    parties: u32,
    arrived: std.atomic.Value(u32) align(std.atomic.cache_line) = .init(0),
    phase: Signal align(std.atomic.cache_line) = .{},

    pub fn wait(b: *Barrier, io: Io, lent: ?*Lent) Io.Cancelable!void {
        if (b.parties <= 1) return;
        const p = b.phase.value.load(.monotonic);
        if (b.arrived.fetchAdd(1, .acq_rel) + 1 == b.parties) {
            b.arrived.store(0, .monotonic);
            b.phase.set(io, p +% 1, lent);
            return;
        }
        _ = try b.phase.waitWhile(io, p, lent);
    }
};

// ---- tests -------------------------------------------------------------------------

const X86 = coretypes.X86CoreType;
const Arm64 = coretypes.Arm64CoreType;

test "chunks cover the range once, in order, on quantum bounds" {
    const lens = [_]usize{ 0, 1, 2, 3, 10, 63, 64, 65, 100, 1000, 12345, std.math.maxInt(usize) };
    const wants = [_]usize{ 0, 1, 2, 3, 7, 64, 65, 1000 };
    const quanta = [_]usize{ 0, 1, 4, 16 };
    for (lens) |len| for (wants) |want| for (quanta) |q| inline for (.{ Shape.even, Shape.rising, Shape.falling }) |shape| {
        const c: Chunks = .init(len, want, shape, q);
        try zigstd.testExpectTrue(c.count >= 1 and c.count <= @max(@min(want, chunks_max), 1));
        try zigstd.testExpectEqual(@as(usize, 0), c.bounds[0]);
        try zigstd.testExpectEqual(len, c.bounds[c.count]);
        for (0..c.count) |i| {
            const ch = c.get(i);
            try zigstd.testExpectEqual(i, ch.index);
            if (len > 0) try zigstd.testExpectTrue(ch.start < ch.end);
            if (i > 0) try zigstd.testExpectEqual(@as(usize, 0), ch.start % @max(q, 1));
        }
    };
}

test "even chunks differ by at most a quantum; rising and falling ones balance a triangle" {
    const even: Chunks = .init(1000, 7, .even, 1);
    try zigstd.testExpectEqual(@as(usize, 7), even.count);
    for (0..even.count) |i| {
        const w = even.get(i).end - even.get(i).start;
        try zigstd.testExpectTrue(w == 142 or w == 143);
    }
    const len = 10000;
    inline for (.{ Shape.rising, Shape.falling }) |shape| {
        const c: Chunks = .init(len, 8, shape, 1);
        try zigstd.testExpectEqual(@as(usize, 8), c.count);
        for (0..c.count) |i| {
            var cost: u64 = 0;
            for (c.get(i).start..c.get(i).end) |u| cost += if (shape == .rising) u + 1 else len - u;
            const share = @as(u64, len) * (len + 1) / 2 / 8;
            // a bound is off by at most a row or two of `len`
            try zigstd.testExpectTrue(cost + 3 * len >= share and cost <= share + 3 * len);
        }
    }
}

test "rules: thresholds at their edges, tiers and the chunk cap" {
    const axpy = rule(.axpy, f64, X86.haswell);
    try zigstd.testExpectEqual(@as(usize, 1), axpy.chunkCount(10000));
    try zigstd.testExpectEqual(@as(u32, 1), axpy.sizeMax(10000));
    try zigstd.testExpectEqual(@as(usize, 2), axpy.chunkCount(10001));
    try zigstd.testExpectEqual(@as(u32, std.math.maxInt(u32)), axpy.sizeMax(10001));
    try zigstd.testExpectEqual(@as(usize, chunks_max), axpy.chunkCount(1 << 40));

    const gemm = rule(.gemm, f64, X86.haswell);
    try zigstd.testExpectEqual(@as(u32, 1), gemm.sizeMax(262144));
    try zigstd.testExpectEqual(@as(u32, 5), gemm.sizeMax(262144 * 5 + 1));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.gemm, std.math.Complex(f32), X86.haswell).sizeMax(32768));
    try zigstd.testExpectEqual(@as(u32, 2), rule(.gemm, std.math.Complex(f32), X86.haswell).sizeMax(65536));

    const ger = rule(.ger, f32, X86.haswell);
    try zigstd.testExpectEqual(@as(u32, 1), ger.sizeMax(8192));
    try zigstd.testExpectEqual(@as(u32, 3), ger.sizeMax(8193));
    try zigstd.testExpectEqual(@as(u32, 3), ger.sizeMax(81920));
    try zigstd.testExpectEqual(@as(u32, std.math.maxInt(u32)), ger.sizeMax(81921));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.ger, std.math.Complex(f64), X86.haswell).sizeMax(9216));
    try zigstd.testExpectEqual(@as(u32, std.math.maxInt(u32)), rule(.ger, std.math.Complex(f64), X86.haswell).sizeMax(9217));

    // `n < 8 MiB / sizeof` serial; complex counts its elements
    try zigstd.testExpectEqual(@as(u32, 1), rule(.swap, f64, X86.haswell).sizeMax(1048575));
    try zigstd.testExpectEqual(@as(usize, 2), rule(.swap, f64, X86.haswell).chunkCount(1048576));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.swap, std.math.Complex(f64), X86.haswell).sizeMax(524287));
    try zigstd.testExpectEqual(@as(usize, 2), rule(.swap, std.math.Complex(f64), X86.haswell).chunkCount(524288));

    const ztrmv = rule(.trmv, std.math.Complex(f64), X86.haswell);
    try zigstd.testExpectEqual(@as(u32, 1), ztrmv.sizeMax(96));
    try zigstd.testExpectEqual(@as(u32, 2), ztrmv.sizeMax(97));
    try zigstd.testExpectEqual(@as(u32, 2), ztrmv.sizeMax(127));
    try zigstd.testExpectEqual(@as(u32, std.math.maxInt(u32)), ztrmv.sizeMax(128));
    try zigstd.testExpectEqual(@as(u32, 2), rule(.trmv, f64, X86.haswell).sizeMax(499));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.symv, f64, X86.haswell).sizeMax(199));
    try zigstd.testExpectEqual(@as(usize, 2), rule(.symv, f64, X86.haswell).chunkCount(200));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.hemv, std.math.Complex(f64), X86.haswell).sizeMax(361));
    try zigstd.testExpectEqual(@as(usize, 2), rule(.hemv, std.math.Complex(f64), X86.haswell).chunkCount(362));
    try zigstd.testExpectEqual(@as(u32, 2), rule(.tpmv, std.math.Complex(f32), X86.haswell).sizeMax(499));

    try zigstd.testExpectEqual(@as(u32, 1), rule(.asum, f32, X86.haswell).sizeMax(100000));
    try zigstd.testExpectEqual(@as(u32, 3), rule(.asum, f32, X86.haswell).sizeMax(300000));
    try zigstd.testExpectEqual(@as(u32, 3), rule(.asum, std.math.Complex(f32), X86.haswell).sizeMax(30000));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.trsm, f64, X86.haswell).sizeMax(1023));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.syrk, f64, X86.haswell).sizeMax(product(&.{ 101, 100, 43 })));
    try zigstd.testExpectEqual(@as(u32, std.math.maxInt(u32)), rule(.syrk, f64, X86.haswell).sizeMax(product(&.{ 101, 100, 44 })));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.syr2k, f64, X86.haswell).sizeMax(999));
    try zigstd.testExpectEqual(@as(u32, std.math.maxInt(u32)), rule(.syr2k, std.math.Complex(f32), X86.haswell).sizeMax(1000));
}

test "rules: OpenBLAS's per-core overrides apply to their core and precision only" {
    try zigstd.testExpectEqual(@as(u32, 1), rule(.dot, f64, Arm64.neoversev1).sizeMax(64500));
    try zigstd.testExpectEqual(@as(u32, 2), rule(.dot, f64, Arm64.neoversev1).sizeMax(64501));
    try zigstd.testExpectEqual(@as(u32, 56), rule(.dot, f64, Arm64.neoversev1).sizeMax(1010000));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.dot, f32, Arm64.neoversev1).sizeMax(110000));
    try zigstd.testExpectEqual(@as(u32, 8), rule(.dot, f32, Arm64.a64fx).sizeMax(1020000));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.dot, std.math.Complex(f64), Arm64.neoversev1).sizeMax(10000));
    try zigstd.testExpectEqual(@as(u32, std.math.maxInt(u32)), rule(.dot, std.math.Complex(f64), Arm64.neoversev1).sizeMax(10001));

    try zigstd.testExpectEqual(@as(u32, 1), rule(.gemm, f32, Arm64.neoversev1).sizeMax(262143));
    try zigstd.testExpectEqual(@as(u32, 6), rule(.gemm, f32, Arm64.neoversev1).sizeMax(262144));
    try zigstd.testExpectEqual(@as(u32, 64), rule(.gemm, f32, Arm64.neoversev2).sizeMax(1073741823));
    // double gemm keeps the generic rule
    try zigstd.testExpectEqual(@as(u32, 1), rule(.gemm, f64, Arm64.neoversev1).sizeMax(262144));

    try zigstd.testExpectEqual(@as(u32, 1), rule(.gemv, f64, Arm64.neoversev1).sizeMax(8099));
    try zigstd.testExpectEqual(@as(u32, 2), rule(.gemv, f64, Arm64.neoversev1).sizeMax(8100));
    try zigstd.testExpectEqual(@as(u32, 4), rule(.gemv, f32, Arm64.neoversev1).sizeMax(25600));
    try zigstd.testExpectEqual(@as(u32, 8), rule(.gemv, f64, Arm64.neoversev2).sizeMax(24964));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.gemv, f64, X86.haswell).sizeMax(460799));
}

test "serial is exact: its floor is one chunk and one worker on every core" {
    inline for (comptime std.meta.tags(Family)) |f| {
        inline for (.{ f32, f64, std.math.Complex(f32), std.math.Complex(f64) }) |T| {
            const floor = comptime serialFloor(f, T);
            try zigstd.testExpectTrue(serial(f, T, floor));
            try zigstd.testExpectTrue(!serial(f, T, floor + 1));
            inline for (coretypes.all) |C| {
                for (std.meta.tags(C)) |c| {
                    const r = rule(f, T, c);
                    try zigstd.testExpectEqual(@as(usize, 1), r.chunkCount(floor));
                    try zigstd.testExpectEqual(@as(u32, 1), r.sizeMax(floor));
                }
            }
        }
    }
    // the per-core overrides lower the floor below a default
    try zigstd.testExpectEqual(@as(u64, 8099), comptime serialFloor(.gemv, f64));
}

test "the worker count is the least of the CPUs, the ceiling and the size" {
    var cfg: Config = .{ .cpu_count = 8 };
    try zigstd.testExpectTrue(cfg.threadCeiling() == null);
    const expect = struct {
        fn f(c: *const Config, size_max: u32, want: u32) !void {
            try zigstd.testExpectEqual(if (sysinfo.single_threaded) 1 else want, workers(c, size_max));
        }
    }.f;
    try expect(&cfg, std.math.maxInt(u32), 8);
    try expect(&cfg, 5, 5);
    try expect(&cfg, 0, 1);
    cfg.goto_num_threads = 6;
    try expect(&cfg, 7, 6);
    cfg.num_threads = 3;
    try expect(&cfg, 7, 3);
    try expect(&cfg, 2, 2);
    cfg.cpu_count = 1;
    try expect(&cfg, 7, 1);
    cfg = .{ .cpu_count = 2, .num_threads = 16 };
    try expect(&cfg, std.math.maxInt(u32), 2);
    try zigstd.testExpectTrue(cpuCount() >= 1);
}

test "chunks follow the size and the core, never the CPU count or the ceiling" {
    const cfgs = [_]Config{ .{}, .{ .cpu_count = 4 }, .{ .cpu_count = 64 }, .{ .cpu_count = 64, .num_threads = 2 }, .{ .cpu_count = 3, .goto_num_threads = 5 } };
    const r = rule(.dot, f64, X86.haswell);
    for ([_]usize{ 1, 10000, 10001, 123457, 1 << 30 }) |n| {
        const base = planWith(&cfgs[0], r, n, n, .even, 8);
        for (cfgs[1..]) |*cfg| {
            const p = planWith(cfg, r, n, n, .even, 8);
            try zigstd.testExpectEqual(base.chunks.count, p.chunks.count);
            try zigstd.testExpectEqualSlices(usize, base.chunks.bounds[0 .. base.chunks.count + 1], p.chunks.bounds[0 .. p.chunks.count + 1]);
            try zigstd.testExpectTrue(p.workers >= 1 and p.workers <= p.chunks.count);
        }
    }
}

test "single-threaded builds resolve to one worker over the same chunks" {
    if (!sysinfo.single_threaded) return error.SkipZigTest;
    _ = blas.vtable();
    try zigstd.testExpectEqual(@as(u32, 1), blas.config.cpu_count);
    const p = plan(.dot, f64, 1 << 20, 1 << 20, .even, 1);
    try zigstd.testExpectEqual(@as(u32, 1), p.workers);
    try zigstd.testExpectEqual(@as(usize, chunks_max), p.chunks.count);
}

/// Mixed magnitudes, so a sum's bits follow its order.
fn testData(x: []f64) void {
    var prng: std.Random.DefaultPrng = .init(0x7e57_da7a);
    const rand = prng.random();
    for (x) |*v| v.* = std.math.ldexp(rand.float(f64) * 2 - 1, rand.intRangeAtMost(i32, -30, 30));
}

/// Each sum rounded to `f64` (`util.rounded`), so on the x87 the bits do not
/// hang on what stays in a register.
const Sum = struct {
    x: []const f64,

    fn partial(s: Sum, c: Chunk) f64 {
        var acc: f64 = 0;
        for (s.x[c.start..c.end]) |v| acc = add(acc, v);
        return acc;
    }

    fn add(a: f64, b: f64) f64 {
        return util.rounded(f64, a + b);
    }
};

/// Every output written once: `marks[i]` counts the chunks that covered `i`.
const Marks = struct {
    marks: []u8,

    fn mark(m: Marks, c: Chunk) void {
        for (m.marks[c.start..c.end]) |*v| v.* += 1;
    }
};

/// The ordered fold, without the machinery; `reverse`d, last chunk first.
fn folded(x: []const f64, chunks: *const Chunks, reverse: bool) f64 {
    var acc: f64 = 0;
    for (0..chunks.count) |i| {
        const part = Sum.partial(.{ .x = x }, chunks.get(if (reverse) chunks.count - 1 - i else i));
        acc = if (i == 0) part else Sum.add(acc, part);
    }
    return acc;
}

/// The reduction and the element-wise loop on `io`, at every worker count and
/// in several claim orders, against `expected`.
fn checkRegions(io: Io, x: []const f64, marks: []u8, chunks: *const Chunks, expected: f64) !void {
    var orders: [5][chunks_max]u8 = undefined;
    for (&orders, 0..) |*o, seed| {
        for (o, 0..) |*v, i| v.* = @intCast(i);
        if (seed == 1) std.mem.reverse(u8, o[0..chunks.count]);
        if (seed > 1) {
            var prng: std.Random.DefaultPrng = .init(seed);
            prng.random().shuffle(u8, o[0..chunks.count]);
        }
    }
    for (test_io.counts()) |w| {
        const p: Plan = .{ .chunks = chunks.*, .workers = w };
        try zigstd.testExpectEqual(@as(u64, @bitCast(expected)), @as(u64, @bitCast(try reduce(io, &p, f64, Sum{ .x = x }, Sum.partial, Sum.add))));
        for (&orders) |*o| {
            const got = try reduceIn(io, &p, f64, Sum{ .x = x }, Sum.partial, Sum.add, o[0..chunks.count]);
            try zigstd.testExpectEqual(@as(u64, @bitCast(expected)), @as(u64, @bitCast(got)));
            @memset(marks, 0);
            try forEachIn(io, &p, Marks{ .marks = marks }, Marks.mark, o[0..chunks.count]);
            for (marks) |v| try zigstd.testExpectEqual(@as(u8, 1), v);
        }
    }
}

test "a reduction gives the same bits for any worker count, Io and claim order" {
    const n = 64123;
    const x = try zigstd.testAllocator.alloc(f64, n);
    defer zigstd.testAllocator.free(x);
    const marks = try zigstd.testAllocator.alloc(u8, n);
    defer zigstd.testAllocator.free(marks);
    testData(x);
    const r: Rule = .{ .serial_max = 1000 };
    const chunks: Chunks = .init(n, r.chunkCount(n), .even, 1);
    try zigstd.testExpectEqual(@as(usize, chunks_max), chunks.count);

    const expected = folded(x, &chunks, false);
    // the data can tell the orders apart
    try zigstd.testExpectTrue(@as(u64, @bitCast(expected)) != @as(u64, @bitCast(folded(x, &chunks, true))));

    try checkRegions(Io.failing, x, marks, &chunks, expected);
    if (comptime test_io.threaded_tests) {
        for (test_io.limits) |limit| {
            var threaded: test_io.TestIo = .init(zigstd.testAllocator, .{ .async_limit = limit });
            defer threaded.deinit();
            try checkRegions(threaded.io(), x, marks, &chunks, expected);
        }
    }
    if (comptime test_io.evented_workers) {
        var ev: Evented = undefined;
        try zigstd.testExpectTrue(test_io.initEvented(&ev, 3));
        defer ev.deinit();
        try checkRegions(ev.io(), x, marks, &chunks, expected);
    }
}

test "workers adopt the caller's rounding and flush to zero and leave the caller's state as it was" {
    if (comptime !test_io.threaded_tests or !FpState.governs) return error.SkipZigTest;
    const n = 64123;
    const x = try zigstd.testAllocator.alloc(f64, n);
    defer zigstd.testAllocator.free(x);
    const marks = try zigstd.testAllocator.alloc(u8, n);
    defer zigstd.testAllocator.free(marks);
    testData(x);
    const chunks: Chunks = .init(n, chunks_max, .even, 1);
    const own: FpState = .capture();
    defer own.apply();
    const rz = own.towardZero();
    const nearest = folded(x, &chunks, false);
    rz.apply();
    const expected = folded(x, &chunks, false);
    // the data can tell the modes apart
    try zigstd.testExpectTrue(@as(u64, @bitCast(nearest)) != @as(u64, @bitCast(expected)));
    for (test_io.limits) |limit| {
        var threaded: test_io.TestIo = .init(zigstd.testAllocator, .{ .async_limit = limit });
        defer threaded.deinit();
        // workers started under the default state, which a spawned thread
        // may inherit
        own.apply();
        try checkRegions(threaded.io(), x, marks, &chunks, nearest);
        rz.apply();
        try checkRegions(threaded.io(), x, marks, &chunks, expected);
        try zigstd.testExpectTrue(FpState.capture().eql(rz));
    }
    if (comptime test_io.evented_workers) {
        var ev: Evented = undefined;
        try zigstd.testExpectTrue(test_io.initEvented(&ev, 3));
        defer ev.deinit();
        own.apply();
        try checkRegions(ev.io(), x, marks, &chunks, nearest);
        rz.apply();
        try checkRegions(ev.io(), x, marks, &chunks, expected);
        try zigstd.testExpectTrue(FpState.capture().eql(rz));
    }
}

/// Chunks that block until their task is canceled.
const Blocking = struct {
    io: Io,
    entered: std.atomic.Value(usize) = .init(0),
    left: std.atomic.Value(usize) = .init(0),

    fn chunk(b: *Blocking, _: Chunk) void {
        _ = b.entered.fetchAdd(1, .acq_rel);
        defer _ = b.left.fetchAdd(1, .acq_rel);
        while (true) {
            b.io.checkCancel() catch {
                b.io.recancel();
                return;
            };
            std.Thread.yield() catch {};
        }
    }

    fn region(io: Io, p: *const Plan, b: *Blocking) Io.Cancelable!void {
        return forEach(io, p, b, chunk);
    }
};

/// Cancels a region of blocking chunks on `io`.
fn cancelRegion(io: Io) !void {
    var b: Blocking = .{ .io = io };
    const p: Plan = .{ .chunks = .init(chunks_max, chunks_max, .even, 1), .workers = 4 };
    var future = io.concurrent(Blocking.region, .{ io, &p, &b }) catch return error.SkipZigTest;
    while (b.entered.load(.acquire) == 0) std.Thread.yield() catch {};
    try zigstd.testExpectError(error.Canceled, future.cancel(io));
    try zigstd.testExpectEqual(b.entered.load(.acquire), b.left.load(.acquire));
    try zigstd.testExpectTrue(b.left.load(.acquire) < chunks_max);
}

test "a canceled region returns error.Canceled after its workers stop" {
    if (comptime !test_io.threaded_tests) return error.SkipZigTest;
    for (test_io.limits) |limit| {
        var threaded: test_io.TestIo = .init(zigstd.testAllocator, .{ .async_limit = limit });
        defer threaded.deinit();
        try cancelRegion(threaded.io());
    }
    if (comptime test_io.evented_workers) {
        var ev: Evented = undefined;
        try zigstd.testExpectTrue(test_io.initEvented(&ev, 3));
        defer ev.deinit();
        try cancelRegion(ev.io());
    }
}

test "worker tasks run on plain futures on Io.Dispatch only" {
    try zigstd.testExpectTrue(!usesFutures(Io.failing));
    if (comptime test_io.hosted) {
        var threaded: Io.Threaded = .init(zigstd.testAllocator, .{});
        defer threaded.deinit();
        try zigstd.testExpectTrue(!usesFutures(threaded.io()));
    }
    if (comptime test_io.hosted and Evented != void) {
        var ev: Evented = undefined;
        if (!test_io.initEvented(&ev, if (Evented == Io.Dispatch) 3 else 0)) return error.SkipZigTest;
        defer ev.deinit();
        try zigstd.testExpectEqual(test_io.evented_workers, usesFutures(ev.io()));
    }
}

threadlocal var probe: u32 = 7;

/// A task's thread and its view of `probe`, which it then sets.
const Probe = struct {
    id: std.Thread.Id = undefined,
    seen: u32 = undefined,

    fn run(p: *Probe) void {
        p.id = std.Thread.getCurrentId();
        p.seen = probe;
        probe = 9;
    }
};

test "a concurrent task runs on a thread of its own with its own thread-locals" {
    if (comptime !test_io.threaded_tests) return error.SkipZigTest;
    var threaded: test_io.TestIo = .init(zigstd.testAllocator, .{});
    defer threaded.deinit();
    const io = threaded.io();
    probe = 8;
    defer probe = 7;
    var p: Probe = .{};
    var group: Io.Group = .init;
    group.concurrent(io, Probe.run, .{&p}) catch return error.SkipZigTest;
    try group.await(io);
    try zigstd.testExpectTrue(p.id != std.Thread.getCurrentId());
    // a fresh copy of the initial value, and this thread's untouched
    try zigstd.testExpectEqual(@as(u32, 7), p.seen);
    try zigstd.testExpectEqual(@as(u32, 8), probe);
}

/// Sleeps on `s` until it leaves 0, then counts itself in `woken`.
fn sleeper(s: *Signal, io: Io, woken: *std.atomic.Value(u32)) void {
    _ = s.waitWhile(io, 0, null) catch return io.recancel();
    _ = woken.fetchAdd(1, .acq_rel);
}

test "a signal wakes every sleeper on an Io.Evented with worker threads" {
    if (comptime !test_io.hosted or Evented == void) return error.SkipZigTest;
    var ev: Evented = undefined;
    if (!test_io.initEvented(&ev, 3)) return error.SkipZigTest;
    defer ev.deinit();
    const io = ev.io();
    const n = 3;
    var s: Signal = .{};
    var woken: std.atomic.Value(u32) = .init(0);
    var group: Io.Group = .init;
    var started: u32 = 0;
    while (started < n) : (started += 1) group.concurrent(io, sleeper, .{ &s, io, &woken }) catch break;
    // all past their spins and on the futex
    while (started == n and s.sleepers.load(.acquire) < n) try io.sleep(.fromMilliseconds(1), .awake);
    try io.sleep(.fromMilliseconds(50), .awake);
    s.set(io, 1, null);
    var waited: u32 = 0;
    while (woken.load(.acquire) < started and waited < 5000) : (waited += 1) try io.sleep(.fromMilliseconds(1), .awake);
    const got = woken.load(.acquire);
    // free any left asleep, so a failure does not hang
    io.futexWake(u32, &s.value.raw, std.math.maxInt(i32));
    // let the tasks end first: an await racing a task's end crashes `Group` (issue #208)
    try io.sleep(.fromMilliseconds(50), .awake);
    try group.await(io);
    if (started < n) return error.SkipZigTest;
    try zigstd.testExpectEqual(@as(u32, n), got);
}

test "product saturates" {
    try zigstd.testExpectEqual(@as(u64, 24), product(&.{ 2, 3, 4 }));
    try zigstd.testExpectEqual(@as(u64, std.math.maxInt(u64)), product(&.{ 1 << 40, 1 << 40 }));
}
