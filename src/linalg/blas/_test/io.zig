//! The threaded tests' `Io` setups and worker counts.

const std = @import("std");
const sysinfo = @import("builtin");
const zigstd = @import("../../../zigstd.zig");
const blas = @import("../__root__.zig");
const threading = @import("../threading.zig");
const root = @import("root");

const Evented = threading.Evented;
const Io = std.Io;

/// Whether `Io.Threaded` runs here, for the threaded tests.
pub const hosted = !sysinfo.single_threaded and switch (sysinfo.target.os.tag) {
    .freestanding, .other, .uefi => false,
    else => true,
};

/// A multicore test kernel's cores (`test/qemu/smp.zig`); null elsewhere.
const kernel_smp: ?type = if (sysinfo.is_test and @hasDecl(root, "qemutest_smp")) root.qemutest_smp else null;

/// Whether the threaded tests run here: hosted, or on a multicore test kernel.
pub const threaded_tests = hosted or kernel_smp != null;

/// The threaded tests' `Io`: `Io.Threaded`, or on a test kernel its cores,
/// set up alike (`init(allocator, .{ .async_limit, .concurrent_limit })`).
pub const TestIo = if (hosted) Io.Threaded else if (kernel_smp) |k| k.Pool else void;

/// Starts `ev` with `extra` threads besides the caller's; false where it
/// cannot be sized so (`Dispatch` sizes its own pool) or does not start.
/// Fibers take large stacks: the page allocator neither fills nor checks them.
pub fn initEvented(ev: *Evented, extra: u32) bool {
    var opts: Evented.InitOptions = .{};
    // the page allocator is thread-safe
    if (@hasField(Evented.InitOptions, "backing_allocator_needs_mutex")) opts.backing_allocator_needs_mutex = false;
    if (@hasField(Evented.InitOptions, "thread_limit")) {
        opts.thread_limit = extra;
    } else if (@hasField(Evented.InitOptions, "n_threads")) {
        opts.n_threads = extra + 1;
    } else if (extra == 0) return false;
    ev.init(std.heap.page_allocator, opts) catch return false;
    return true;
}

/// Whether the threaded tests also run on `Io.Evented` with worker threads:
/// on `Dispatch` through `Futures` only. `Group` crashes there and on `Uring`,
/// whose futures are reused early (issues.md).
pub const evented_workers = hosted and threading.dispatch_futures and Evented == Io.Dispatch;

/// Rounds of the race stress tests (`race stress`, levels 1 to 3): enough for
/// a weak-memory rig to reorder something, few enough to stay quick.
pub const stress_rounds = 1000;

/// The worker counts the tests try: 1, 2, 3, 7, this machine's, and more
/// than there are chunks.
pub fn counts() [6]u32 {
    return .{ 1, 2, 3, 7, threading.cpuCount(), threading.chunks_max + 1 };
}

/// The `async_limit`s of the threaded tests: all inline, fewer than asked,
/// oversubscribed.
pub const limits = [_]Io.Limit{ .nothing, .limited(3), .limited(8) };

/// `f(ctx, io)` on each threaded-test `Io`: `Io.failing`; `TestIo` at each of
/// `limits` (async and concurrent); `Evented` with 3 workers.
pub fn eachTestIo(ctx: anytype, comptime f: anytype) !void {
    try f(ctx, Io.failing);
    if (comptime threaded_tests) for (limits) |limit| {
        var threaded: TestIo = .init(std.testing.allocator, .{ .async_limit = limit, .concurrent_limit = limit });
        defer threaded.deinit();
        try f(ctx, threaded.io());
    };
    if (comptime evented_workers) {
        var ev: Evented = undefined;
        try zigstd.testExpectTrue(initEvented(&ev, 3));
        defer ev.deinit();
        try f(ctx, ev.io());
    }
}

/// The worker counts of the bit-identity runs (`WorkerRuns`).
pub const worker_runs = [_]u32{ 1, 3, 8 };

/// The bit-identity runs at each of `worker_runs`: one worker on `Io.failing`,
/// the rest many chunks on an `Io.Threaded`. Hosted only; `init` after the
/// vtable's, which sets `cpu_count`; `deinit` restores it.
pub const WorkerRuns = struct {
    threaded: Io.Threaded,
    cpus: u32,

    pub fn init(allocator: std.mem.Allocator) WorkerRuns {
        return .{ .threaded = .init(allocator, .{ .async_limit = .limited(8) }), .cpus = blas.config.cpu_count };
    }

    pub fn deinit(r: *WorkerRuns) void {
        r.threaded.deinit();
        threading.test_serial_max = null;
        blas.config.cpu_count = r.cpus;
    }

    /// Sets up run `c` and returns its `Io`.
    pub fn io(r: *WorkerRuns, c: usize) Io {
        const w = worker_runs[c];
        blas.config.cpu_count = w;
        threading.test_serial_max = if (w == 1) null else 4096;
        return if (w == 1) Io.failing else r.threaded.io();
    }
};

/// Each run's bytes in `outs` against the first's.
pub fn expectSameBytes(outs: anytype) !void {
    for (1..outs.len) |c| try std.testing.expectEqualSlices(u8, std.mem.sliceAsBytes(outs[0]), std.mem.sliceAsBytes(outs[c]));
}
