//! `threading.zig`'s tests of the regions, their tasks and level 3's signal.

const std = @import("std");
const FpState = @import("../../FpState.zig");
const util = @import("../../util.zig");
const test_io = @import("../../_test/io.zig");
const threading = @import("../../threading.zig");
const region = @import("../region.zig");
const zigstd = @import("../../../../../zigstd.zig");

const Io = std.Io;
const chunks_max = threading.chunks_max;
const Rule = threading.Rule;
const Chunk = threading.Chunk;
const Chunks = threading.Chunks;
const Plan = threading.Plan;
const forEach = threading.forEach;
const reduce = threading.reduce;
const forEachIn = region.forEachIn;
const reduceIn = region.reduceIn;
const Evented = threading.Evented;
const usesFutures = threading.usesFutures;
const Signal = threading.Signal;

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
