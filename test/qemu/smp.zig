//! The multicore rigs' cores: the boot core starts the rest (PSCI `CPU_ON`, INIT-SIPI-SIPI), and `Pool` is an `Io` handing tasks to idle
//! cores. No queue: with none idle `async` runs inline and `concurrent` fails, as `Io.Threaded` past its limits. Waits poll (`wfe`, `mwait`). Each architecture's bring-up is in `smp/`.

const std = @import("std");
const sysinfo = @import("builtin");
const root_mod = @import("root");

const Io = std.Io;
const Alignment = std.mem.Alignment;
const Allocator = std.mem.Allocator;
const Atomic = std.atomic.Value;

/// The most cores used, the boot core included.
pub const max_cores = 8;
/// Each started core's stack.
pub const stack_size = 256 << 10;
/// Room for a core's TLS block and its thread control block.
const tls_room = 4096;
/// The largest task context and future result a core takes, and their
/// alignment; larger ones run inline or fail as with no idle core.
const context_max = 256;
const result_max = 256;
const slot_align = 64;

const arch = switch (sysinfo.target.cpu.arch) {
    .x86_64 => @import("smp/x86_64.zig"),
    .aarch64 => @import("smp/aarch64.zig"),
    else => @compileError("qemutest: no SMP boot for " ++ @tagName(sysinfo.target.cpu.arch)),
};

/// The started cores' stacks: index `i`'s top is `i * stack_size` past the
/// start (index 0 is the boot core, on the linker script's stack).
export var qemutest_ap_stacks: [max_cores - 1][stack_size]u8 align(16) = @splat(@splat(0));

/// Cores up, the boot core included.
var online: Atomic(usize) = .init(1);

threadlocal var core_index: usize = 0;
/// The running task's cancellation; the boot core's never comes.
threadlocal var current: *Ctx = &root;
var root: Ctx = .{ .core = 0 };

/// This core: 0 for the boot core.
pub fn coreIndex() usize {
    return core_index;
}

/// Cores up, the boot core included.
pub fn count() usize {
    return online.load(.acquire);
}

/// Sets up the boot core's TLS, then starts the other cores and waits for
/// them. First thing at boot: nothing before may touch a thread-local.
pub fn start() void {
    setupTls(0);
    arch.initCore();
    for (&cores, 0..) |*c, i| c.slot.ctx.core = @intCast(i);
    const asked = arch.startCores(max_cores - 1);
    // a core that never comes stays offline; none has been seen not to
    var spins: usize = 0;
    while (online.load(.acquire) < 1 + asked and spins < 1 << 30) : (spins += 1) std.atomic.spinLoopHint();
    var buf: [64]u8 = undefined;
    root_mod.testWrite(std.mem.print(&buf, "{d} of {d} cores up\n", .{ count(), 1 + asked }) catch unreachable);
}

/// Where a started core lands, on its stack: TLS, then tasks forever.
export fn qemutest_ap_main(index: usize) callconv(.c) noreturn {
    setupTls(index);
    arch.initCore();
    core_index = index;
    const s = &cores[index].slot;
    current = &s.ctx;
    s.state.store(.idle, .release);
    _ = online.fetchAdd(1, .acq_rel);
    while (true) {
        wait(@intFromPtr(&s.state), s, Slot.isReady);
        s.run();
    }
}

// ---- thread-local storage ----------------------------------------------------------

fn linkerValue(comptime name: []const u8) usize {
    return @intFromPtr(@extern(*allowzero const u8, .{ .name = name }));
}

/// Core `i`'s TLS block from the image the linker script bounds (`__tls_*`),
/// laid out as lld assigns offsets: below the thread pointer on x86_64 (a
/// self pointer at it), past a 16-byte control block on aarch64. The script
/// aligns the segment to 64, so its offsets need no misalignment term.
fn setupTls(i: usize) void {
    const image = linkerValue("__tls_image");
    const size = linkerValue("__tls_size");
    const al = @max(linkerValue("__tls_align"), 1);
    const base = @intFromPtr(&cores[i].tls);
    if (al > 64) fail("TLS aligned past 64");
    const tp, const block = switch (sysinfo.target.cpu.arch) {
        .x86_64 => blk: {
            const tp = std.mem.alignForward(usize, base + size + al, al);
            break :blk .{ tp, tp - std.mem.alignForward(usize, size, al) };
        },
        .aarch64 => blk: {
            const tp = std.mem.alignForward(usize, base, al);
            break :blk .{ tp, tp + std.mem.alignForward(usize, 16, al) };
        },
        else => unreachable,
    };
    if (@max(tp + 16, block + size) > base + tls_room) fail("TLS block past its room");
    const dst: [*]u8 = @ptrFromInt(block);
    const src: [*]const u8 = @ptrFromInt(linkerValue("__tls_start"));
    @memcpy(dst[0..image], src[0..image]);
    @memset(dst[image..size], 0);
    arch.setThreadPointer(tp);
}

pub fn fail(msg: []const u8) noreturn {
    root_mod.testWrite("smp: ");
    @panic(msg);
}

// ---- cores, waiting and waking -----------------------------------------------------

const Core = struct {
    /// Bumped to wake the core (`ring`).
    doorbell: Atomic(u32) align(std.atomic.cache_line) = .init(0),
    /// The address the core waits on (`wait`), 0 when it does not.
    waiting_on: Atomic(usize) = .init(0),
    /// Its task (started cores).
    slot: Slot = .{},
    tls: [tls_room]u8 align(64) = @splat(0),
};

var cores: [max_cores]Core = @splat(.{});

/// Waits until `ready(arg)`, sleeping until a `wake` at `addr` or a ring of
/// this core. `ready` must read with `.seq_cst`: with `wake`'s order either
/// it sees the change or `wake` sees this waiter.
fn wait(addr: usize, arg: anytype, comptime ready: fn (@TypeOf(arg)) bool) void {
    const c = &cores[core_index];
    defer c.waiting_on.store(0, .release);
    while (true) {
        _ = c.waiting_on.swap(addr, .seq_cst);
        const seen = c.doorbell.load(.acquire);
        arch.arm(&c.doorbell);
        if (ready(arg)) return;
        arch.sleep(&c.doorbell, seen);
    }
}

/// Wakes every core waiting at `addr`; the caller changed what it waits for.
fn wake(addr: usize) void {
    // a full barrier between that change and reading the waiters
    _ = fence.fetchAdd(0, .seq_cst);
    for (&cores, 0..) |*c, i| if (c.waiting_on.load(.seq_cst) == addr) ring(i);
}

var fence: Atomic(u32) = .init(0);

fn ring(i: usize) void {
    _ = cores[i].doorbell.fetchAdd(1, .seq_cst);
    arch.signal();
}

// ---- tasks -------------------------------------------------------------------------

const Cancel = enum(u8) { none, requested, acknowledged };

/// A task's cancellation, shared by the tasks it runs inline.
const Ctx = struct {
    cancel: Atomic(Cancel) = .init(.none),
    protection: Io.CancelProtection = .unblocked,
    /// The core that runs it, which a cancel rings.
    core: u8,

    fn pending(ctx: *const Ctx) bool {
        return ctx.protection == .unblocked and ctx.cancel.load(.seq_cst) == .requested;
    }

    fn check(ctx: *Ctx) Io.Cancelable!void {
        if (!ctx.pending()) return;
        if (ctx.cancel.cmpxchgStrong(.requested, .acknowledged, .acq_rel, .monotonic) == null) return error.Canceled;
    }

    /// From `none` only: an acknowledged cancel stays so until `recancel`.
    fn request(ctx: *Ctx) void {
        if (ctx.cancel.cmpxchgStrong(.none, .requested, .acq_rel, .monotonic) == null) ring(ctx.core);
    }
};

const Kind = enum { async, concurrent };

/// A started core's one task: idle, claimed by a dispatcher, ready, running,
/// then idle again (group tasks) or done until awaited (futures).
const Slot = struct {
    state: Atomic(State) = .init(.offline),
    /// Guards `group` and a cancel's choice of tasks against the task's end.
    lock: std.atomic.Mutex = .unlocked,
    ctx: Ctx = .{ .core = 0 },
    kind: Kind = .async,
    pool: *Pool = undefined,
    group: ?*Io.Group = null,
    start: Start = undefined,
    context: [context_max]u8 align(slot_align) = undefined,
    result: [result_max]u8 align(slot_align) = undefined,

    const State = enum(u8) { offline, idle, claimed, ready, running, done };
    const Start = union(enum) {
        group: *const fn (*const anyopaque) void,
        future: *const fn (*const anyopaque, *anyopaque) void,
    };

    fn isReady(s: *Slot) bool {
        return s.state.load(.seq_cst) == .ready;
    }

    fn isDone(s: *Slot) bool {
        return s.state.load(.seq_cst) == .done;
    }

    fn acquire(s: *Slot) void {
        while (!s.lock.tryLock()) std.atomic.spinLoopHint();
    }

    fn run(s: *Slot) void {
        s.state.store(.running, .monotonic);
        const p = s.pool;
        switch (s.start) {
            .group => |f| {
                f(&s.context);
                s.acquire();
                const g = s.group.?;
                s.group = null;
                s.lock.unlock();
                // the pool and the group may go once their counts say so
                p.release(s.kind);
                s.state.store(.idle, .release);
                if (@atomicRmw(usize, &g.state, .Sub, 1, .acq_rel) == 1) wake(@intFromPtr(&g.state));
            },
            .future => |f| {
                f(&s.context, &s.result);
                p.release(s.kind);
                s.state.store(.done, .release);
                wake(@intFromPtr(&s.state));
            },
        }
    }
};

/// An `Io` over the started cores, set up as `Io.Threaded` is for the tests:
/// `async_limit` tasks at most run on other cores at once, the rest inline,
/// and `concurrent_limit` concurrent ones, past which they fail.
pub const Pool = struct {
    limits: std.EnumArray(Kind, Io.Limit),
    busy: std.EnumArray(Kind, Atomic(usize)) = .initFill(.init(0)),

    pub const Options = struct {
        async_limit: Io.Limit = .unlimited,
        concurrent_limit: Io.Limit = .unlimited,
    };

    pub fn init(_: Allocator, o: Options) Pool {
        return .{ .limits = .init(.{ .async = o.async_limit, .concurrent = o.concurrent_limit }) };
    }

    pub fn deinit(p: *Pool) void {
        for (p.busy.values) |b| std.debug.assert(b.load(.acquire) == 0);
    }

    pub fn io(p: *Pool) Io {
        return .{ .userdata = p, .vtable = &vtable };
    }

    /// An idle core claimed for a task within `kind`'s limit, given the
    /// context and result fit; null otherwise.
    fn launch(p: *Pool, kind: Kind, group: ?*Io.Group, result_len: usize, result_alignment: Alignment, context: []const u8, context_alignment: Alignment, start_fn: Slot.Start) ?*Slot {
        if (context.len > context_max or result_len > result_max) return null;
        if (context_alignment.toByteUnits() > slot_align or result_alignment.toByteUnits() > slot_align) return null;
        const busy = p.busy.getPtr(kind);
        if (busy.fetchAdd(1, .acq_rel) >= @backingInt(p.limits.get(kind))) {
            _ = busy.fetchSub(1, .acq_rel);
            return null;
        }
        const s = for (cores[1..]) |*c| {
            if (c.slot.state.cmpxchgStrong(.idle, .claimed, .acquire, .monotonic) == null) break &c.slot;
        } else {
            _ = busy.fetchSub(1, .acq_rel);
            return null;
        };
        if (group) |g| {
            _ = @atomicRmw(usize, &g.state, .Add, 1, .acq_rel);
            g.token.store(g, .release);
        }
        s.acquire();
        s.group = group;
        s.ctx.cancel.store(.none, .release);
        s.lock.unlock();
        s.ctx.protection = .unblocked;
        s.kind = kind;
        s.pool = p;
        s.start = start_fn;
        @memcpy(s.context[0..context.len], context);
        s.state.store(.ready, .seq_cst);
        ring(s.ctx.core);
        return s;
    }

    fn release(p: *Pool, kind: Kind) void {
        _ = p.busy.getPtr(kind).fetchSub(1, .acq_rel);
    }

    const vtable: Io.VTable = blk: {
        var vt = Io.failing.vtable.*;
        vt.async = asyncFn;
        vt.concurrent = concurrentFn;
        vt.await = awaitFn;
        vt.cancel = cancelFn;
        vt.groupAsync = groupAsync;
        vt.groupConcurrent = groupConcurrent;
        vt.groupAwait = groupAwait;
        vt.groupCancel = groupCancel;
        vt.recancel = recancel;
        vt.swapCancelProtection = swapCancelProtection;
        vt.checkCancel = checkCancel;
        vt.futexWait = futexWait;
        vt.futexWaitUncancelable = futexWaitUncancelable;
        vt.futexWake = futexWake;
        break :blk vt;
    };
};

fn poolOf(userdata: ?*anyopaque) *Pool {
    return @ptrCast(@alignCast(userdata));
}

fn asyncFn(userdata: ?*anyopaque, result: []u8, result_alignment: Alignment, context: []const u8, context_alignment: Alignment, start_fn: *const fn (*const anyopaque, *anyopaque) void) ?*Io.AnyFuture {
    if (poolOf(userdata).launch(.async, null, result.len, result_alignment, context, context_alignment, .{ .future = start_fn })) |s| return @ptrCast(s);
    start_fn(context.ptr, result.ptr);
    return null;
}

fn concurrentFn(userdata: ?*anyopaque, result_len: usize, result_alignment: Alignment, context: []const u8, context_alignment: Alignment, start_fn: *const fn (*const anyopaque, *anyopaque) void) Io.ConcurrentError!*Io.AnyFuture {
    const s = poolOf(userdata).launch(.concurrent, null, result_len, result_alignment, context, context_alignment, .{ .future = start_fn }) orelse return error.ConcurrencyUnavailable;
    return @ptrCast(s);
}

fn awaitFn(_: ?*anyopaque, any_future: *Io.AnyFuture, result: []u8, _: Alignment) void {
    const s: *Slot = @ptrCast(@alignCast(any_future));
    wait(@intFromPtr(&s.state), s, Slot.isDone);
    @memcpy(result, s.result[0..result.len]);
    s.state.store(.idle, .release);
}

fn cancelFn(userdata: ?*anyopaque, any_future: *Io.AnyFuture, result: []u8, result_alignment: Alignment) void {
    const s: *Slot = @ptrCast(@alignCast(any_future));
    s.ctx.request();
    awaitFn(userdata, any_future, result, result_alignment);
}

fn groupAsync(userdata: ?*anyopaque, g: *Io.Group, context: []const u8, context_alignment: Alignment, start_fn: *const fn (*const anyopaque) void) void {
    if (poolOf(userdata).launch(.async, g, 0, .@"1", context, context_alignment, .{ .group = start_fn }) == null) start_fn(context.ptr);
}

fn groupConcurrent(userdata: ?*anyopaque, g: *Io.Group, context: []const u8, context_alignment: Alignment, start_fn: *const fn (*const anyopaque) void) Io.ConcurrentError!void {
    _ = poolOf(userdata).launch(.concurrent, g, 0, .@"1", context, context_alignment, .{ .group = start_fn }) orelse return error.ConcurrencyUnavailable;
}

/// Requests a cancel of every task of `g` still running or about to.
fn cancelGroup(g: *Io.Group) void {
    for (cores[1..]) |*c| {
        const s = &c.slot;
        s.acquire();
        defer s.lock.unlock();
        if (s.group == g) s.ctx.request();
    }
}

const GroupWait = struct {
    g: *Io.Group,
    /// The awaiting task's, while a cancel of it would pass on to `g`.
    ctx: ?*Ctx,

    fn ready(w: GroupWait) bool {
        if (@atomicLoad(usize, &w.g.state, .seq_cst) == 0) return true;
        return if (w.ctx) |ctx| ctx.pending() else false;
    }
};

fn groupAwait(_: ?*anyopaque, g: *Io.Group, _: *anyopaque) Io.Cancelable!void {
    const ctx = current;
    var canceled = false;
    while (true) {
        wait(@intFromPtr(&g.state), GroupWait{ .g = g, .ctx = if (canceled) null else ctx }, GroupWait.ready);
        if (@atomicLoad(usize, &g.state, .seq_cst) == 0) break;
        if (ctx.check()) |_| {} else |_| {
            canceled = true;
            cancelGroup(g);
        }
    }
    g.token.store(null, .release);
    if (canceled) return error.Canceled;
}

fn groupCancel(_: ?*anyopaque, g: *Io.Group, _: *anyopaque) void {
    cancelGroup(g);
    wait(@intFromPtr(&g.state), GroupWait{ .g = g, .ctx = null }, GroupWait.ready);
    g.token.store(null, .release);
}

fn recancel(_: ?*anyopaque) void {
    current.cancel.store(.requested, .release);
}

fn swapCancelProtection(_: ?*anyopaque, new: Io.CancelProtection) Io.CancelProtection {
    const ctx = current;
    defer ctx.protection = new;
    return ctx.protection;
}

fn checkCancel(_: ?*anyopaque) Io.Cancelable!void {
    return current.check();
}

const FutexWait = struct {
    ptr: *const u32,
    expected: u32,
    ctx: ?*Ctx,

    fn ready(w: FutexWait) bool {
        if (@atomicLoad(u32, w.ptr, .seq_cst) != w.expected) return true;
        return if (w.ctx) |ctx| ctx.pending() else false;
    }
};

/// With a timeout, one look: there is no clock, and futex waits may end early.
fn futexWait(_: ?*anyopaque, ptr: *const u32, expected: u32, timeout: Io.Timeout) Io.Cancelable!void {
    const ctx = current;
    try ctx.check();
    if (timeout != .none) return std.atomic.spinLoopHint();
    wait(@intFromPtr(ptr), FutexWait{ .ptr = ptr, .expected = expected, .ctx = ctx }, FutexWait.ready);
    try ctx.check();
}

fn futexWaitUncancelable(_: ?*anyopaque, ptr: *const u32, expected: u32) void {
    wait(@intFromPtr(ptr), FutexWait{ .ptr = ptr, .expected = expected, .ctx = null }, FutexWait.ready);
}

fn futexWake(_: ?*anyopaque, ptr: *const u32, _: u32) void {
    wake(@intFromPtr(ptr));
}
