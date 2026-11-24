//! `select.zig`'s tests: the vtable's defaults, the init handshake, `init`'s
//! environment, a test's own table, and the `q`/`x` slots no core fills.

const std = @import("std");
const sysinfo = @import("builtin");
const zigstd = @import("../../../../zigstd.zig");
const blas = @import("../__root__.zig");
const select = @import("../select.zig");
const portable = @import("../portable.zig");
const registry = @import("../registry.zig");
const cpu = @import("../cpu/__root__.zig");
const Vtable = @import("../Vtable.zig");

const init = select.init;
const shutdown = select.shutdown;
const vtable = select.vtable;
const table = select.table;
const kernel = select.kernel;
const coreType = select.coreType;
const coreName = select.coreName;
const testUse = select.testUse;
const testRestore = select.testRestore;
const LogSink = select.LogSink;

// lives here, not in `Vtable.zig`, so that file stays import-free for the
// per-core kernel objects, which use its type (see `kernels/__root__.zig`)
test "VTable default-initializes with every kernel pointer null" {
    const vt: Vtable = .{};
    try zigstd.testExpectTrue(vt.sgemm_kernel == null);
    try zigstd.testExpectTrue(vt.dgemm_kernel == null);
    try zigstd.testExpectTrue(vt.cgemm_kernel_n == null);
    try zigstd.testExpectTrue(vt.zgemm_kernel_n == null);
    try zigstd.testExpectTrue(vt.init == null);
    try zigstd.testExpectEqual(@as(i32, 0), vt.dtb_entries);
}

test "a kernel call initializes by itself, and init after that does nothing" {
    shutdown();
    defer _ = vtable();
    const before = select.inits.load(.monotonic);
    try zigstd.testExpectTrue(select.phase.load(.acquire) == .uninit);
    try zigstd.testExpectTrue(kernel(.daxpy_k) == table().daxpy_k.?);
    try zigstd.testExpectTrue(select.phase.load(.acquire) == .ready);
    try zigstd.testExpectEqual(before + 1, select.inits.load(.monotonic));
    try zigstd.testExpectTrue(table().dtb_entries > 0);
    try zigstd.testExpectTrue(coreType() != null);
    var env: std.process.Environ.Map = .init(zigstd.testAllocator);
    defer env.deinit();
    try zigstd.testExpectTrue(!init(&env));
    try zigstd.testExpectEqual(before + 1, select.inits.load(.monotonic));
}

test "init keeps its own copy of the environment, shutdown goes back to the portable kernels" {
    shutdown();
    defer _ = vtable();
    try zigstd.testExpectTrue(select._vtable.daxpy_k == select.generic_vtable.daxpy_k);
    const detected = @tagName(cpu.detect());
    // the core that runs it: a detected core not linked into this build falls back
    const chosen = cpu.choose(detected, null, LogSink{});
    const runs = @tagName(chosen.served orelse chosen.core);
    {
        var env: std.process.Environ.Map = .init(zigstd.testAllocator);
        defer env.deinit();
        try env.put("OPENBLAS_CORETYPE", detected);
        try env.put("OPENBLAS_NUM_THREADS", "3");
        try zigstd.testExpectTrue(init(&env));
    }
    // the map is gone; the copy is not
    try zigstd.testExpectEqualStrings(detected, blas.config.coreType().?);
    try zigstd.testExpectEqual(@as(i32, 3), blas.config.num_threads);
    try zigstd.testExpectEqualStrings(runs, coreName());

    shutdown();
    try zigstd.testExpectEqual(@as(i32, 0), blas.config.num_threads);
    try zigstd.testExpectTrue(select._vtable.daxpy_k == select.generic_vtable.daxpy_k);
    var env: std.process.Environ.Map = .init(zigstd.testAllocator);
    defer env.deinit();
    try env.put("OPENBLAS_CORETYPE", "no-such-core");
    try zigstd.testExpectTrue(init(&env));
    try zigstd.testExpectEqualStrings(runs, coreName());
}

test "a test's table needs no core, and goes back to the unset state" {
    shutdown();
    defer _ = vtable();
    const before = select.inits.load(.monotonic);
    const vt: Vtable = select.generic_vtable;
    const saved = testUse(&vt);
    try zigstd.testExpectTrue(select.phase.load(.acquire) == .ready);
    try zigstd.testExpectTrue(coreType() == null);
    testRestore(saved);
    try zigstd.testExpectTrue(select.phase.load(.acquire) == .uninit);
    try zigstd.testExpectEqual(before, select.inits.load(.monotonic));
}

var race_seen: [4]?cpu.CoreType = @splat(null);

fn raceKernel(i: usize) void {
    race_seen[i] = if (kernel(.daxpy_k) == select._vtable.daxpy_k.?) select.core_type else null;
}

test "threads racing the first kernel call select the core once and all see it" {
    if (sysinfo.single_threaded) return error.SkipZigTest;
    shutdown();
    defer _ = vtable();
    race_seen = @splat(null);
    const before = select.inits.load(.monotonic);
    var threads: [race_seen.len]std.Thread = undefined;
    {
        // a failed spawn still joins the threads already running
        var spawned: usize = 0;
        defer for (threads[0..spawned]) |t| t.join();
        for (&threads, 0..) |*t, i| {
            t.* = try std.Thread.spawn(.{}, raceKernel, .{i});
            spawned += 1;
        }
    }
    try zigstd.testExpectEqual(before + 1, select.inits.load(.monotonic));
    try zigstd.testExpectTrue(select.phase.load(.acquire) == .ready);
    for (race_seen) |seen| try zigstd.testExpectTrue(seen != null and seen.? == select.core_type.?);
    try zigstd.testExpectTrue(select._vtable.dtb_entries > 0);
}

/// The `q`/`x` kernel slots' names.
const qx_slots = blk: {
    @setEvalBranchQuota(100_000);
    const info = @typeInfo(Vtable).@"struct";
    var names: []const []const u8 = &.{};
    for (info.field_names, info.field_types) |f, F| {
        const qx = f[0] == 'q' or f[0] == 'x' or (f[0] == 'i' and (f[1] == 'q' or f[1] == 'x'));
        if (qx and @typeInfo(F) == .optional) names = names ++ .{f};
    }
    break :blk names;
};

test "no core puts a kernel of its own in a q or x slot" {
    @setEvalBranchQuota(100_000);
    var base: Vtable = .{};
    portable.install(&base);
    inline for (@typeInfo(cpu.CoreType).@"enum".field_names) |name| {
        if (registry.isLinked(@field(cpu.CoreType, name))) try qxUntouched(&base, @field(cpu.CoreType, name));
    }
}

// the Linux rigs match this line against the core they expect; only
// there, as `zig build` reports a passing test's stderr as a failed command
// (issue #276)
test "the core in use" {
    if (sysinfo.os.tag != .linux) return error.SkipZigTest;
    if (std.testing.environ.getPosix("RSTD_TEST_PROGRESS") == null) return error.SkipZigTest;
    std.debug.print("blas core: {s}\n", .{coreName()});
}

/// One core per call: inlined into the loop, each core's table and install
/// would keep a frame slot of their own (~400 KB on a 64-bit debug build).
fn qxUntouched(base: *const Vtable, comptime core: cpu.CoreType) !void {
    var vt = base.*;
    cpu.kernel.install(&vt, core);
    inline for (qx_slots) |f| try zigstd.testExpectTrue(@field(vt, f) == @field(base, f));
}
