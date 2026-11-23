//! A canceled threaded `gemm` returns `error.Canceled` with no worker left running.

const std = @import("std");
const zigstd = @import("../../../../../../zigstd.zig");
const Vtable = @import("../../../Vtable.zig");
const driver = @import("../../../blas3/driver.zig");
const threading = @import("../../../threading.zig");
const test_io = @import("../../io.zig");
const types = @import("../../../types.zig");
const check = @import("../../check.zig");
const th = @import("../threaded.zig");

const Io = std.Io;
const Error = driver.Error;
const fill = check.fillMixed;
const shrunk = th.shrunk;

/// A kernel that holds every call until released, counting them and the
/// threads that make them.
pub const Hold = struct {
    pub const Kernel = @typeInfo(@FieldType(Vtable, "dgemm_kernel")).optional.child;
    pub const cc = types.kernelCc(f64);
    pub const BlasLong = types.BlasLong;
    pub const Atomic = std.atomic.Value(u32);

    pub var real: Kernel = undefined;
    pub var calls: Atomic = .init(0);
    pub var inside: Atomic = .init(0);
    pub var held: std.atomic.Value(bool) = .init(true);
    pub var first: std.atomic.Value(std.Thread.Id) = .init(0);
    pub var others: Atomic = .init(0);
    pub var settled: std.atomic.Value(bool) = .init(false);

    pub fn kernel(m: BlasLong, n: BlasLong, k: BlasLong, alpha: f64, a: [*]f64, b: [*]f64, c: [*]f64, ldc: BlasLong) callconv(cc) i32 {
        _ = calls.fetchAdd(1, .acq_rel);
        _ = inside.fetchAdd(1, .acq_rel);
        const id = std.Thread.getCurrentId();
        if (first.cmpxchgStrong(0, id, .acq_rel, .acquire)) |f| {
            if (f != id) _ = others.fetchAdd(1, .acq_rel);
        }
        defer _ = inside.fetchSub(1, .acq_rel);
        while (held.load(.acquire)) std.Thread.yield() catch {};
        return real(m, n, k, alpha, a, b, c, ldc);
    }

    /// Settles once a second thread is in the kernel (or after ten seconds
    /// without one), then releases the calls 50 ms later.
    pub fn release(io: Io) void {
        while (calls.load(.acquire) == 0) std.Thread.yield() catch {};
        // workers start well after the caller on a slow runner (issue #300)
        for (0..10_000) |_| {
            if (others.load(.acquire) > 0) break;
            io.sleep(.fromMilliseconds(1), .awake) catch {};
        }
        settled.store(true, .release);
        io.sleep(.fromMilliseconds(50), .awake) catch {};
        held.store(false, .release);
    }
};

pub fn canceledGemm(io: Io, vt: *const Vtable, a: []f64, b: []f64, c: []f64, n: usize, k: usize) Error!void {
    return driver.gemm(std.testing.allocator, io, f64, vt, false, false, .n, n, n, k, 1, a.ptr, n, b.ptr, k, 1, c.ptr, n, 4);
}

test "a canceled threaded gemm returns error.Canceled with no worker left running" {
    if (comptime !test_io.hosted) return error.SkipZigTest;
    var vt = shrunk(f64);
    Hold.real = vt.dgemm_kernel.?;
    vt.dgemm_kernel = &Hold.kernel;
    const n = 64;
    const k = 2000;
    const al = std.testing.allocator;
    const a = try al.alloc(f64, n * k);
    defer al.free(a);
    const b = try al.alloc(f64, k * n);
    defer al.free(b);
    const c = try al.alloc(f64, n * n);
    defer al.free(c);
    fill(f64, a, 1);
    fill(f64, b, 2);
    for ([_]Io.Limit{ .limited(3), .limited(8) }) |limit| {
        var threaded: Io.Threaded = .init(al, .{ .async_limit = limit, .concurrent_limit = limit });
        defer threaded.deinit();
        try cancelGemm(threaded.io(), &vt, a, b, c, n, k);
    }
    if (comptime test_io.evented_workers) {
        var ev: threading.Evented = undefined;
        try zigstd.testExpectTrue(test_io.initEvented(&ev, 3));
        defer ev.deinit();
        try cancelGemm(ev.io(), &vt, a, b, c, n, k);
    }
}

/// Cancels a held `gemm` on `io` once workers are in the kernel.
pub fn cancelGemm(io: Io, vt: *const Vtable, a: []f64, b: []f64, c: []f64, n: usize, k: usize) !void {
    Hold.calls.store(0, .release);
    Hold.held.store(true, .release);
    Hold.first.store(0, .release);
    Hold.others.store(0, .release);
    Hold.settled.store(false, .release);
    // a thread of its own sleeps on a blocking `Io`; a fiber `Io` needs a fiber
    const releaser = try std.Thread.spawn(.{}, Hold.release, .{Io.Threaded.global_single_threaded.io()});
    var future = io.concurrent(canceledGemm, .{ io, vt, a, b, c, n, k }) catch return error.SkipZigTest;
    while (!Hold.settled.load(.acquire)) std.Thread.yield() catch {};
    const result = future.cancel(io);
    releaser.join();
    try zigstd.testExpectError(error.Canceled, result);
    // workers did run beside the caller
    try zigstd.testExpectTrue(Hold.others.load(.acquire) > 0);
    // nothing calls the kernel once the call has returned
    const after = Hold.calls.load(.acquire);
    try zigstd.testExpectEqual(@as(u32, 0), Hold.inside.load(.acquire));
    io.sleep(.fromMilliseconds(20), .awake) catch {};
    try zigstd.testExpectEqual(after, Hold.calls.load(.acquire));
}
