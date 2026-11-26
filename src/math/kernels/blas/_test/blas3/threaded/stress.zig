//! Race stress: threaded level 3 keeps one worker's bits and its guard bands, round after round.

const std = @import("std");
const zigstd = @import("../../../../../../zigstd.zig");
const driver = @import("../../../blas3/driver.zig");
const threading = @import("../../../threading.zig");
const test_io = @import("../../io.zig");
const types = @import("../../../types.zig");
const check = @import("../../check.zig");
const bitwise = @import("../../bitwise.zig");
const guard = @import("../../guard.zig");
const trt = @import("routines.zig");

const Io = std.Io;
const Blocking = driver.Blocking;
const fill = check.fillMixed;
const scalar = check.val;
const sameBits = bitwise.sameBits;
const Routine = trt.Routine;
const routine = trt.routine;
const shrunkCore = trt.shrunkCore;
const tame = trt.tame;

/// `r` (`gemm` if null) on `io` at 2 to 8 workers, round after round,
/// against one worker's bits on `Io.failing`. `C` (`B`) is `m x n` at
/// `ld`; every element past it is a guard band that must keep its value.
/// Scratch comes NaN-poisoned, so a read before its packing shows.
pub fn stress(comptime T: type, io: Io, r: ?Routine, opt: [3]bool, m: usize, n: usize, k: usize, rounds: usize) !void {
    const al = std.testing.allocator;
    const vt = shrunkCore(T);
    const ld = @max(m, n, k) + 3;
    const a = try al.alloc(T, ld * ld);
    defer al.free(a);
    const b = try al.alloc(T, ld * ld);
    defer al.free(b);
    const c0 = try al.alloc(T, ld * ld);
    defer al.free(c0);
    const want = try al.alloc(T, c0.len);
    defer al.free(want);
    const got = try al.alloc(T, c0.len);
    defer al.free(got);
    fill(T, a, 41);
    fill(T, b, 42);
    fill(T, c0, 43);
    if (r == .trsm) tame(T, a, ld);
    var poison: guard.Poison = .{ .child = al };
    const pa = poison.allocator();
    const alpha = scalar(T, -0.75, 0.25);
    const beta = scalar(T, 0.5, 0);
    for (0..rounds + 1) |round| {
        const w: u32 = if (round == 0) 1 else @intCast(2 + round % 7);
        const out = if (round == 0) want else got;
        @memcpy(out, c0);
        const rio = if (round == 0) Io.failing else io;
        if (r) |rr| {
            try routine(T, &vt, pa, rio, rr, opt, m, n, k, alpha, beta, a.ptr, b.ptr, ld, out.ptr, w);
        } else {
            try driver.gemm(pa, rio, T, &vt, opt[0], opt[1], .n, m, n, k, alpha, a.ptr, ld, b.ptr, ld, beta, out.ptr, ld, w);
        }
        if (round > 0 and !sameBits(T, want, got)) {
            std.debug.print("{s} {?t} round {d} w={d}\n", .{ @typeName(T), r, round, w });
            return error.TestUnexpectedResult;
        }
    }
}

test "race stress: threaded level 3 keeps one worker's bits and its guard bands, round after round" {
    if (comptime !test_io.threaded_tests) return error.SkipZigTest;
    var threaded: test_io.TestIo = .init(std.testing.allocator, .{});
    defer threaded.deinit();
    try stresses(threaded.io());
    if (comptime test_io.evented_workers) {
        var ev: threading.Evented = undefined;
        try zigstd.testExpectTrue(test_io.initEvented(&ev, 3));
        defer ev.deinit();
        try stresses(ev.io());
    }
}

/// The race stress test's routines on `io`.
pub fn stresses(io: Io) !void {
    const rounds = test_io.stress_rounds;
    // shared packed `B` behind its barriers: `ic` and `jc` past 1, several slices
    {
        const blk = Blocking.of(f64, &shrunkCore(f64));
        try stress(f64, io, null, .{ false, true, false }, 9 * blk.um + 5, 2 * blk.r + 3, 2 * blk.q + 3, rounds);
    }
    // a triangle's balanced groups, both passes
    {
        const blk = Blocking.of(types.Complex32, &shrunkCore(types.Complex32));
        try stress(types.Complex32, io, .syr2k, .{ false, true, false }, 2 * blk.r + 3, 2 * blk.r + 3, 2 * blk.q + 3, rounds / 2);
    }
    // independent panels through `forEach`
    {
        const blk = Blocking.of(f32, &shrunkCore(f32));
        try stress(f32, io, .trsm, .{ true, false, false }, 2 * blk.r + 3, 2 * blk.r + 3, 1, rounds / 2);
    }
}
