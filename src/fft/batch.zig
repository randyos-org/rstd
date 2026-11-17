//! A batch's transforms through `io`: as many at once as the scratch holds,
//! each run on a scratch slot of its own, the slots claimed by the BLAS's
//! `threading.forEach` (its workers, `Io.Dispatch` workaround and
//! floating-point state). Transforms are independent, so the bits depend on
//! neither the worker count nor the order.

const std = @import("std");
const threading = @import("../linalg/blas/threading.zig");

const Io = std.Io;

/// `one(ctx, k, slot)` for every transform `k < count`, at most `slots` at
/// once; `slot` (below `slots`) is `k`'s scratch.
pub fn run(io: Io, count: usize, slots: usize, ctx: anytype, comptime one: fn (@TypeOf(ctx), usize, usize) void) Io.Cancelable!void {
    if (count == 0) return;
    const want = @min(count, slots);
    if (want <= 1) {
        for (0..count) |k| {
            try threading.checkCancel(io);
            one(ctx, k, 0);
        }
        return;
    }
    // a chunk per slot: chunk `index` owns slot `index`
    const p: threading.Plan = .{ .chunks = .init(count, want, .even, 1), .workers = @min(@as(u32, @intCast(@min(want, threading.chunks_max))), threading.cpuCount()) };
    const W = struct {
        c: @TypeOf(ctx),
        fn each(w: @This(), ch: threading.Chunk) void {
            for (ch.start..ch.end) |k| one(w.c, k, ch.index);
        }
    };
    try threading.forEach(io, &p, W{ .c = ctx }, W.each);
}
