//! Threaded level 1 (`blas1/interface.zig`) against one worker, bit for bit,
//! at any worker count, `Io` and core.

const std = @import("std");
const blas = @import("../../__root__.zig");
const threading = @import("../../threading.zig");
const test_io = @import("../io.zig");
const types = @import("../../types.zig");
const util = @import("../../util.zig");
const portable = @import("../../portable.zig");
const Vtable = @import("../../Vtable.zig");
const check = @import("../check.zig");
const bitwise = @import("../bitwise.zig");
const interface = @import("../../blas1/interface.zig");

const Io = std.Io;
const Complex32 = types.Complex32;
const Float = types.Float;
const isComplex = types.isComplex;
const raw = util.raw;
const all = check.all_types;
const quantum = interface.quantum;
const Extreme = interface.Extreme;
const asum = interface.asum;
const axpy = interface.axpy;
const axpyc = interface.axpyc;
const copy = interface.copy;
const dot = interface.dot;
const dotc = interface.dotc;
const dsdot = interface.dsdot;
const extreme = interface.extreme;
const index = interface.index;
const rot = interface.rot;
const scal = interface.scal;
const sum = interface.sum;
const swap = interface.swap;

const scalarOf = check.val;

/// A threaded test's length: `test_serial_max` makes it a dozen chunks.
const threaded_n = 17 * quantum + 517;

/// Every threaded routine's results and outputs for `x0` and `y0`, as bytes.
fn level1Bits(comptime T: type, al: std.mem.Allocator, io: Io, inc: [2]isize, x0: []const T, y0: []const T, out: *std.ArrayList(u8)) !void {
    const n: isize = threaded_n;
    const x = try al.dupe(T, x0);
    defer al.free(x);
    const y = try al.dupe(T, y0);
    defer al.free(y);
    const a = scalarOf(T, 1.25, -0.75);
    const pos: isize = @intCast(@abs(inc[0]));
    try copy(io, T, n, x.ptr, inc[0], y.ptr, inc[1]);
    try bitwise.appendBits(al, out, y);
    try swap(io, T, n, x.ptr, inc[0], y.ptr, inc[1]);
    try bitwise.appendBits(al, out, x);
    try scal(io, T, n, a, x.ptr, pos);
    try bitwise.appendBits(al, out, x);
    try rot(io, T, n, x.ptr, inc[0], y.ptr, inc[1], 0.6, -0.8);
    try bitwise.appendBits(al, out, x);
    try bitwise.appendBits(al, out, y);
    try axpy(io, T, n, a, x.ptr, inc[0], y.ptr, inc[1]);
    try bitwise.appendBits(al, out, y);
    try bitwise.appendBits(al, out, try dot(io, T, n, x.ptr, inc[0], y.ptr, inc[1]));
    if (comptime isComplex(T)) {
        try axpyc(io, T, n, a, x.ptr, inc[0], y.ptr, inc[1]);
        try bitwise.appendBits(al, out, y);
        try bitwise.appendBits(al, out, try dotc(io, T, n, x.ptr, inc[0], y.ptr, inc[1]));
    }
    if (comptime T == f32) try bitwise.appendBits(al, out, try dsdot(io, n, x.ptr, inc[0], y.ptr, inc[1]));
    try bitwise.appendBits(al, out, try asum(io, T, n, x.ptr, pos));
    try bitwise.appendBits(al, out, try sum(io, T, n, x.ptr, pos));
    try extremeBits(T, io, x.ptr, pos, out, al);
}

fn extremeBits(comptime T: type, io: Io, x: [*]T, inc: isize, out: *std.ArrayList(u8), al: std.mem.Allocator) !void {
    const whiches = if (comptime isComplex(T)) .{ Extreme.amax, Extreme.amin } else .{ Extreme.amax, Extreme.amin, Extreme.max, Extreme.min };
    inline for (whiches) |w| {
        try bitwise.appendBits(al, out, try extreme(io, T, w, threaded_n, x, inc));
        try bitwise.appendBits(al, out, try index(io, T, w, threaded_n, x, inc));
    }
}

/// The extremes over chunks against one kernel call: the fold is exact.
fn expectExtremesExact(comptime T: type, io: Io, x: [*]T, inc: isize) !void {
    const whiches = if (comptime isComplex(T)) .{ Extreme.amax, Extreme.amin } else .{ Extreme.amax, Extreme.amin, Extreme.max, Extreme.min };
    inline for (whiches) |w| {
        const v = portable.fetch(T, comptime w.op(), blas.vtable)(threaded_n, raw(T, x), inc);
        const i = portable.fetch(T, comptime w.indexOp(), blas.vtable)(threaded_n, raw(T, x), inc);
        try std.testing.expect(bitwise.same(@TypeOf(v), v, try extreme(io, T, w, threaded_n, x, inc), true));
        try std.testing.expectEqual(@as(isize, @intCast(i)), try index(io, T, w, threaded_n, x, inc));
    }
}

/// `x` with NaNs, infinities and, for complex, finite elements whose keys
/// overflow, spread over several chunks (`variant` 1 and 2).
fn specials(comptime T: type, x: []T, inc: usize, variant: usize) void {
    const F = Float(T);
    const at = [_]usize{ 2 * quantum + 3, 7 * quantum + 11, 11 * quantum + 1, 16 * quantum + 9 };
    for (at, 0..) |i, k| {
        const v: F = switch (variant) {
            1 => if (k % 2 == 0) std.math.nan(F) else std.math.inf(F),
            else => std.math.floatMax(F) * (0.5 + @as(F, @floatFromInt(k % 3)) / 8),
        };
        x[i * inc] = if (comptime isComplex(T)) .init(v, if (k == 3) -v else v) else if (k == 1) -v else v;
    }
}

fn threadedType(comptime T: type, al: std.mem.Allocator, tio: Io) !void {
    const cases = [_][2]isize{ .{ 1, 1 }, .{ 2, -1 }, .{ -3, 2 } };
    for (cases) |inc| for (0..3) |variant| {
        if (variant == 2 and !isComplex(T)) continue;
        const x0 = try al.alloc(T, threaded_n * @abs(inc[0]));
        defer al.free(x0);
        const y0 = try al.alloc(T, threaded_n * @abs(inc[1]));
        defer al.free(y0);
        check.fill(T, x0, 11 + @as(u64, @abs(inc[0])));
        check.fill(T, y0, 12 + @as(u64, @abs(inc[1])));
        if (variant > 0) specials(T, x0, @abs(inc[0]), variant);
        var want: std.ArrayList(u8) = .empty;
        defer want.deinit(al);
        blas.config.cpu_count = 1;
        try level1Bits(T, al, Io.failing, inc, x0, y0, &want);
        for ([_]u32{ 3, 8 }) |cpus| for ([_]Io{ tio, Io.failing }) |io| {
            blas.config.cpu_count = cpus;
            var got: std.ArrayList(u8) = .empty;
            defer got.deinit(al);
            try level1Bits(T, al, io, inc, x0, y0, &got);
            try std.testing.expectEqualSlices(u8, want.items, got.items);
            try expectExtremesExact(T, io, x0.ptr, @intCast(@abs(inc[0])));
        };
    };
}

fn threadedCore(vt: *const Vtable) !void {
    const saved = blas.testUse(vt);
    defer blas.testRestore(saved);
    const al = std.testing.allocator;
    var threaded: test_io.TestIo = .init(al, .{});
    defer threaded.deinit();
    inline for (all) |T| threadedType(T, al, threaded.io()) catch |e| {
        check.say("{s}: ", .{@typeName(T)});
        return e;
    };
}

test "interface: threaded level 1 gives one worker's bits at any worker count, Io and core" {
    if (comptime !test_io.threaded_tests) return error.SkipZigTest;
    _ = blas.coreName();
    const cpus = blas.config.cpu_count;
    defer blas.config.cpu_count = cpus;
    // a dozen chunks of `threaded_n`
    threading.test_serial_max = 3 * quantum / 2;
    defer threading.test_serial_max = null;
    try check.eachCore(threadedCore);
}

/// The race stress test's length: eight chunks under `test_serial_max`.
const stress_n = 8 * quantum + 517;

/// One round of the race stress test: an element-wise routine, a reduction
/// and an exact fold, strided, so the gaps between elements are guard bands.
fn stressBits(comptime T: type, al: std.mem.Allocator, io: Io, x: []T, y0: []const T, out: *std.ArrayList(u8)) !void {
    const y = try al.dupe(T, y0);
    defer al.free(y);
    try axpy(io, T, stress_n, scalarOf(T, 1.25, -0.75), x.ptr, 2, y.ptr, -3);
    try bitwise.appendBits(al, out, y);
    try bitwise.appendBits(al, out, try dot(io, T, stress_n, x.ptr, 2, y.ptr, -3));
    try bitwise.appendBits(al, out, try index(io, T, .amax, stress_n, x.ptr, 2));
}

test "interface: race stress: threaded level 1 keeps one worker's bits and its guard bands, round after round" {
    if (comptime !test_io.threaded_tests) return error.SkipZigTest;
    _ = blas.coreName();
    const cpus = blas.config.cpu_count;
    defer blas.config.cpu_count = cpus;
    threading.test_serial_max = quantum;
    defer threading.test_serial_max = null;
    const al = std.testing.allocator;
    var threaded: test_io.TestIo = .init(al, .{ .async_limit = .limited(8) });
    defer threaded.deinit();
    inline for (.{ f64, Complex32 }) |T| {
        const x = try al.alloc(T, stress_n * 2);
        defer al.free(x);
        const y0 = try al.alloc(T, stress_n * 3);
        defer al.free(y0);
        check.fill(T, x, 21);
        check.fill(T, y0, 22);
        var want: std.ArrayList(u8) = .empty;
        defer want.deinit(al);
        blas.config.cpu_count = 1;
        try stressBits(T, al, Io.failing, x, y0, &want);
        var got: std.ArrayList(u8) = .empty;
        defer got.deinit(al);
        for (0..test_io.stress_rounds) |round| {
            blas.config.cpu_count = @intCast(2 + round % 7);
            got.clearRetainingCapacity();
            try stressBits(T, al, threaded.io(), x, y0, &got);
            if (!bitwise.sameBytes(want.items, got.items)) std.testing.expectEqualSlices(u8, want.items, got.items) catch |e| {
                check.say("{s} round {d}: ", .{ @typeName(T), round });
                return e;
            };
        }
    }
}
