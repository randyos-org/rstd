const std = @import("std");
const uefi = std.os.uefi;
const Io = std.Io;
const Alignment = std.mem.Alignment;

const state = @import("state.zig");

pub fn crashHandler(_: ?*anyopaque) void {
    // marks crashing task canceled so others stop; one task, nothing to notify
}

pub fn async(
    _: ?*anyopaque,
    result: []u8,
    _: Alignment,
    context: []const u8,
    _: Alignment,
    start: *const fn (context: *const anyopaque, result: *anyopaque) void,
) ?*Io.AnyFuture {
    // run eagerly to completion; null tells caller result's ready, await is no-op
    start(context.ptr, result.ptr);
    return null;
}

pub const concurrent = Io.failingConcurrent;
pub const await = Io.unreachableAwait;
pub const cancel = Io.unreachableCancel;

pub fn groupAsync(
    _: ?*anyopaque,
    _: *Io.Group,
    context: []const u8,
    _: Alignment,
    start: *const fn (context: *const anyopaque) void,
) void {
    start(context.ptr);
}

pub const groupConcurrent = Io.failingGroupConcurrent;
pub const groupAwait = Io.unreachableGroupAwait;
pub const groupCancel = Io.unreachableGroupCancel;

pub fn swapCancelProtection(userdata: ?*anyopaque, new: Io.CancelProtection) Io.CancelProtection {
    _ = userdata;
    const prev = state.cancel_protection;
    state.cancel_protection = new;
    return prev;
}

pub const recancel = Io.unreachableRecancel;
pub const checkCancel = Io.unreachableCheckCancel;

pub const futexWait = Io.noFutexWait;

pub fn futexWaitUncancelable(userdata: ?*anyopaque, _: *const u32, _: u32) void {
    _ = userdata;
    // std.debug parks a panicking task in `while (true)`; stall so it's not a busy spin
    if (state.bootServices()) |bs| bs.stall(1000) catch {};
}

pub const futexWake = Io.noFutexWake;

pub const batchAwaitAsync = Io.unreachableBatchAwaitAsync;
pub const batchAwaitConcurrent = Io.unreachableBatchAwaitConcurrent;
pub const batchCancel = Io.unreachableBatchCancel;
