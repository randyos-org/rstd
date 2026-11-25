//! The machinery of `threading.md`'s levels 1 and 2: thresholds, the worker
//! count, fixed chunks, the claim loop and the ordered combine. Chunks follow
//! the problem size and the core only, never the worker count, and partials
//! fold in chunk order, so the bits are the same for any worker count, `Io` or
//! scheduling. Level 3 uses the thresholds, the worker count and the grain.
//!
//! Here, the worker count and the plan; the rest is in `threading/`: the
//! thresholds (`rules.zig`), `chunks.zig`, the claim loop (`region.zig`,
//! `tasks.zig`), `partials.zig`, level 3's `team.zig` and `grain.zig`.

const std = @import("std");
const sysinfo = @import("builtin");
const blas = @import("__root__.zig");
const Config = @import("Config.zig");
const coretypes = @import("cpu/coretypes.zig");
const rules = @import("threading/rules.zig");
const chunking = @import("threading/chunks.zig");
const region = @import("threading/region.zig");
const tasks = @import("threading/tasks.zig");
const partials = @import("threading/partials.zig");
const teams = @import("threading/team.zig");
const grain = @import("threading/grain.zig");

pub const chunks_max = chunking.chunks_max;
pub const Shape = chunking.Shape;
pub const Chunk = chunking.Chunk;
pub const Chunks = chunking.Chunks;
pub const Family = rules.Family;
pub const Tier = rules.Tier;
pub const Rule = rules.Rule;
pub const rule = rules.rule;
pub const product = rules.product;
pub const team_grain = grain.team_grain;
pub const crew = grain.crew;
pub const matrix_grain = grain.matrix_grain;
pub const share = grain.share;
pub const forEach = region.forEach;
pub const reduce = region.reduce;
pub const forEachPart = region.forEachPart;
pub const checkCancel = region.checkCancel;
pub const dispatch_futures = tasks.dispatch_futures;
pub const Evented = tasks.Evented;
pub const usesFutures = tasks.usesFutures;
pub const Rows = partials.Rows;
pub const Partials = partials.Partials;
pub const team = teams.team;
pub const Lent = teams.Lent;
pub const spin_limit = teams.spin_limit;
pub const Signal = teams.Signal;
pub const Barrier = teams.Barrier;

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

/// A call's chunks and the most workers that run them.
pub const Plan = struct {
    chunks: Chunks,
    workers: u32,
};

/// The least `serial_max` of `family`'s rule for `T` on any core.
pub fn serialFloor(comptime family: Family, comptime T: type) u64 {
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
