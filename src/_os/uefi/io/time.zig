//! "seconds since boot" clock for the bootloader.
//!
//! RuntimeServices.getTime() only tracks whole seconds (e.g. OVMF always
//! reports nanosecond=0) -- no sub-second data there regardless of print
//! precision.
//!
//! A cycle counter (rdtsc) would give real precision but is arch-specific;
//! UEFI isn't, so don't bake that assumption in even though we're
//! x86_64-only today.
//!
//! Instead: periodic Boot Services timer event (100us), plain API on every
//! UEFI arch.

const std = @import("std");
const uefi = std.os.uefi;
const Io = std.Io;
const log = std.log.scoped(.boottime);
// const BootServices = uefi.tables.BootServices;

const state = @import("state.zig");

/// ns elapsed since init(), updated by tick on each timer fire
var elapsed_ns: u64 = std.math.maxInt(u64); // "-1", as a sentinel for "never started"

/// granularity of getElapsedNanoseconds(); pub for io.zig's clock resolution
pub const tick_period_ns: u64 = 100_000; // 100us
/// `setTimer`'s `trigger_time` is in units of 100ns.
const tick_period_100ns: u64 = tick_period_ns / 100;

fn tick(event: uefi.Event, context: ?*anyopaque) callconv(uefi.cc) void {
    _ = event;
    _ = context;
    _ = @atomicRmw(u64, &elapsed_ns, .Add, tick_period_ns, .seq_cst);
}

/// Start periodic timer event backing getTime(). Call once, early in init.
///
/// Never needs explicit stop: `tick`'s callback only runs via Boot
/// Services' event loop (firmware code) -- once ExitBootServices succeeds
/// that loop never runs again. Event memory is Boot Services-owned,
/// reclaimable at that point (see main.zig), nothing to leak. Tearing down
/// the real HW timer is firmware's job during ExitBootServices anyway --
/// moot since we can't call setTimer/closeEvent after that point either.
pub fn init() !void {
    const bs = state.bootServices() orelse return error.BootServicesUnavailable;
    const event = try bs.createEvent(
        .{ .timer = true, .signal = true },
        .{ .tpl = .notify, .function = &tick },
    );
    try bs.setTimer(event, .periodic, tick_period_100ns);
    elapsed_ns = 0;
}

/// ns elapsed since init(), in std.Io clock units. Returns maxInt(u64)
/// sentinel (not a crash) if init() never called/failed.
pub fn getElapsedNanoseconds() u64 {
    return @atomicLoad(u64, &elapsed_ns, .seq_cst);
}

/// seconds since init() as f64, shape rstd.logging.get_time wants; -1 if
/// clock never started
pub fn getTimeSeconds() f64 {
    const ns = getElapsedNanoseconds();
    if (ns == std.math.maxInt(u64)) return -1;
    return @as(f64, @floatFromInt(ns)) / std.time.ns_per_s;
}

fn monotonicNow() Io.Timestamp {
    const ns = getElapsedNanoseconds();
    // map "never started" sentinel to zero, better than starting at end of time
    return .{ .nanoseconds = if (ns == std.math.maxInt(u64)) 0 else ns };
}

pub fn now(userdata: ?*anyopaque, clock: Io.Clock) Io.Timestamp {
    _ = userdata;
    switch (clock) {
        .real => {
            const result = uefi.system_table.runtime_services.getTime() catch return monotonicNow();
            return .{ .nanoseconds = result[0].toEpoch() };
        },
        // no other process/thread to distinguish; every other clock = tick counter
        .awake, .boot, .cpu_process, .cpu_thread => return monotonicNow(),
    }
}

pub fn clockResolution(userdata: ?*anyopaque, clock: Io.Clock) Io.Clock.ResolutionError!Io.Duration {
    _ = userdata;
    return switch (clock) {
        // RTC only reliably tracks whole seconds (see file doc comment)
        .real => .fromSeconds(1),
        .awake, .boot, .cpu_process, .cpu_thread => .fromNanoseconds(tick_period_ns),
    };
}

pub fn sleep(userdata: ?*anyopaque, timeout: Io.Timeout) Io.Cancelable!void {
    _ = userdata;
    const ns: i96 = switch (timeout) {
        // "no timeout" = block forever; would hang single-threaded world, so return
        .none => return,
        .duration => |d| d.raw.nanoseconds,
        .deadline => |ts| ts.raw.nanoseconds - now(null, ts.clock).nanoseconds,
    };
    if (ns <= 0) return;
    const bs = state.bootServices() orelse return;
    bs.stall(@intCast(@divTrunc(ns, std.time.ns_per_us))) catch {};
}
