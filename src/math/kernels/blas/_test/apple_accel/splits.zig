//! rstd's workers around Accelerate: size-determined splits, and what does not fit `c_int`.

const std = @import("std");
const types = @import("../../types.zig");
const blas = @import("../../__root__.zig");
const l3 = @import("../../blas3/core.zig");
const l2 = @import("../../blas2/core.zig");
const l1 = @import("../../blas1/interface.zig");
const threading = @import("../../threading.zig");
const test_io = @import("../io.zig");
const accel = @import("../../apple_accel.zig");
const zigstd = @import("../../../../../zigstd.zig");
const check = @import("../check.zig");
const rt = @import("routing.zig");

const Io = std.Io;
const Complex32 = types.Complex32;
const Complex64 = types.Complex64;
const isComplex = types.isComplex;
const Float = types.Float;
const call = accel.call;
const calls = &accel.calls;
const load = accel.load;
const matrix_quantum = accel.matrix_quantum;
const routing = accel.routing;
const span = accel.span;
const unpinned = &accel.unpinned;
const vector_quantum = accel.vector_quantum;
const expectCalled = rt.expectCalled;
const val = rt.val;

// ---- rstd's workers around Accelerate --------------------------------------------

/// A split test's sizes: under a `test_serial_max` of `split_unit`, several
/// parts each, every one with a tail.
pub const split_unit = vector_quantum;
pub const split_n = 5 * vector_quantum + 37;
pub const split_m = 3 * matrix_quantum + 11;
pub const split_k = 2 * matrix_quantum + 5;

/// A split test's outputs: their bits, and their values as `T`.
pub fn Outputs(comptime T: type) type {
    return struct {
        const Self = @This();
        bits: std.ArrayList(u8) = .empty,
        vals: std.ArrayList(T) = .empty,

        fn deinit(o: *Self, al: std.mem.Allocator) void {
            o.bits.deinit(al);
            o.vals.deinit(al);
        }

        fn add(o: *Self, al: std.mem.Allocator, v: []const T) !void {
            try o.bits.appendSlice(al, std.mem.sliceAsBytes(v));
            try o.vals.appendSlice(al, v);
        }

        /// A scalar result: `T`, or a real one of any float type.
        fn one(o: *Self, al: std.mem.Allocator, v: anytype) !void {
            try o.bits.appendSlice(al, std.mem.asBytes(&v));
            const V = @TypeOf(v);
            try o.vals.append(al, if (V == T) v else val(T, @floatCast(v), 0));
        }
    };
}

/// Every routine that splits, routed, over `T`, on `io`.
pub fn splitRun(comptime T: type, al: std.mem.Allocator, io: Io, o: *Outputs(T)) !void {
    const cx = comptime isComplex(T);
    const al1 = val(T, 1.25, -0.5);
    const be1 = val(T, 0.5, 0.25);
    const nl: isize = split_n;
    const x = try al.alloc(T, 2 * split_n);
    defer al.free(x);
    const y = try al.alloc(T, 3 * split_n);
    defer al.free(y);
    for ([_][2]isize{ .{ 1, 1 }, .{ -2, 3 } }) |inc| {
        check.fill(T, x, 1);
        check.fill(T, y, 2);
        const px: isize = @intCast(@abs(inc[0]));
        try l1.copy(io, T, nl, x.ptr, inc[0], y.ptr, inc[1]);
        try o.add(al, y);
        try l1.swap(io, T, nl, x.ptr, inc[0], y.ptr, inc[1]);
        try o.add(al, x);
        try l1.scal(io, T, nl, al1, x.ptr, px);
        try o.add(al, x);
        try l1.rot(io, T, nl, x.ptr, inc[0], y.ptr, inc[1], 0.6, -0.8);
        try o.add(al, x);
        try o.add(al, y);
        try l1.axpy(io, T, nl, al1, x.ptr, inc[0], y.ptr, inc[1]);
        try o.add(al, y);
        try o.one(al, try l1.dot(io, T, nl, x.ptr, inc[0], y.ptr, inc[1]));
        if (cx) try o.one(al, try l1.dotc(io, T, nl, x.ptr, inc[0], y.ptr, inc[1]));
        if (T == f32) try o.one(al, try l1.dsdot(io, nl, x.ptr, inc[0], y.ptr, inc[1]));
        try o.one(al, try l1.asum(io, T, nl, x.ptr, px));
    }

    const ld = split_m + 1;
    const len = ld * split_m;
    const a = try al.alloc(T, len);
    defer al.free(a);
    const b0 = try al.alloc(T, len);
    defer al.free(b0);
    const c = try al.alloc(T, len);
    defer al.free(c);
    check.fill(T, a, 3);
    check.fill(T, b0, 4);
    // a dominant diagonal, for the solves
    for (0..split_m) |i| a[i * (ld + 1)] = val(T, 2 * split_m, 0);
    const l: isize = ld;
    const big: isize = split_m;
    const small: isize = split_k;
    // the long side first, then the other way about: each side splits
    for ([_][2]isize{ .{ big, small }, .{ small, big } }) |mn| {
        const m = mn[0];
        const n = mn[1];
        for ([_]types.Trans{ .n, if (cx) .c else .t }) |t| {
            check.fill(T, y, 5);
            try l2.gemv(al, io, T, t, m, n, al1, a.ptr, l, x.ptr, -2, be1, y.ptr, 3);
            try o.add(al, y);
        }
        @memcpy(c, b0);
        try l2.ger(al, io, T, m, n, al1, x.ptr, 1, y.ptr, -3, c.ptr, l);
        try o.add(al, c);
        for ([_][2]types.Trans{ .{ .n, .n }, .{ .t, if (cx) .c else .t } }) |tt| {
            @memcpy(c, b0);
            try l3.gemm(al, io, T, tt[0], tt[1], m, n, small, al1, a.ptr, l, a.ptr + 1, l, be1, c.ptr, l);
            try o.add(al, c);
        }
        for ([_]types.Side{ .l, .r }) |sd| {
            @memcpy(c, b0);
            try l3.symm(al, io, T, sd, .u, m, n, al1, a.ptr, l, a.ptr + 1, l, be1, c.ptr, l);
            try o.add(al, c);
            if (cx) {
                @memcpy(c, b0);
                try l3.hemm(al, io, T, sd, .l, m, n, al1, a.ptr, l, a.ptr + 1, l, be1, c.ptr, l);
                try o.add(al, c);
            }
            @memcpy(c, b0);
            try l3.trmm(al, io, T, sd, .u, .t, .n, m, n, al1, a.ptr, l, c.ptr, l);
            try o.add(al, c);
            @memcpy(c, b0);
            try l3.trsm(al, io, T, sd, .l, .n, .n, m, n, al1, a.ptr, l, c.ptr, l);
            try o.add(al, c);
        }
    }
}

/// `splitRun`'s outputs at one worker against 3 and 8, on a threaded `Io`
/// and an inline one; then against unsplit calls, within rounding.
pub fn splitType(comptime T: type, al: std.mem.Allocator, tio: Io) !void {
    var want: Outputs(T) = .{};
    defer want.deinit(al);
    blas.config.cpu_count = 1;
    threading.test_serial_max = split_unit;
    const before = calls.load(.monotonic);
    try splitRun(T, al, Io.failing, &want);
    const split_calls = calls.load(.monotonic) - before;
    for ([_]u32{ 3, 8 }) |cpus| for ([_]Io{ tio, Io.failing }) |io| {
        blas.config.cpu_count = cpus;
        var got: Outputs(T) = .{};
        defer got.deinit(al);
        try splitRun(T, al, io, &got);
        try zigstd.testExpectEqualSlices(u8, want.bits.items, got.bits.items);
    };
    var whole_: Outputs(T) = .{};
    defer whole_.deinit(al);
    threading.test_serial_max = std.math.maxInt(u64);
    const before_whole = calls.load(.monotonic);
    try splitRun(T, al, Io.failing, &whole_);
    // the parts were Accelerate calls of their own (on arm64 level 3, and
    // level 2 to 8 MB, make one: issues #729, #838)
    try zigstd.testExpectTrue(split_calls > calls.load(.monotonic) - before_whole);
    const tol = check.tolerance(Float(T), split_k);
    for (whole_.vals.items, want.vals.items, 0..) |w, g, i| {
        const scale = @max(1, @abs(if (comptime isComplex(T)) w.re else w) + @abs(if (comptime isComplex(T)) w.im else 0));
        check.expectClose(T, w, g, tol * scale) catch |e| {
            check.say("{s} output {d}: ", .{ @typeName(T), i });
            return e;
        };
    }
}

test "routed calls split by size on rstd's workers: one worker's bits at any worker count and Io" {
    if (comptime !test_io.hosted) return error.SkipZigTest;
    _ = blas.vtable();
    if (!routing()) return error.SkipZigTest;
    const cpus = blas.config.cpu_count;
    defer blas.config.cpu_count = cpus;
    defer threading.test_serial_max = null;
    const al = zigstd.testAllocator;
    var threaded: Io.Threaded = .init(al, .{});
    defer threaded.deinit();
    const before = unpinned.load(.monotonic);
    inline for (.{ f32, f64, Complex32, Complex64 }) |T| splitType(T, al, threaded.io()) catch |e| {
        check.say("{s}: ", .{@typeName(T)});
        return e;
    };
    // every worker ran its parts single-threaded (macOS 15+)
    try zigstd.testExpectEqual(before, unpinned.load(.monotonic));
}

test "on the 32-bit table, what does not fit c_int runs rstd's kernels" {
    // sizes past `c_int` need a 64-bit `isize`, as on every macOS
    if (comptime @bitSizeOf(isize) <= 32) return error.SkipZigTest;
    _ = blas.vtable();
    if (!routing()) return error.SkipZigTest;
    const lp = load(.{ .lp64 = true }) orelse return error.SkipZigTest;
    try zigstd.testExpectTrue(lp.blas == .lp64);
    const saved = accel.loaded;
    accel.loaded = lp;
    defer accel.loaded = saved;
    const io = check.io;
    const al = zigstd.testAllocator;
    const big: isize = std.math.maxInt(c_int) + 1;
    var x = [_]f64{ 1, -2, 3 };
    var y = [_]f64{ 0.5, 4, -1 };

    var before = calls.load(.monotonic);
    try zigstd.testExpectEqual(@as(f64, -10.5), try l1.dot(io, f64, 3, &x, 1, &y, 1));
    try expectCalled(true, before, "dot");
    try zigstd.testExpectEqual(@as(isize, 3), try l1.index(io, f64, .largest_magnitude, 3, &x, 1));
    // one element at a stride past `c_int`
    before = calls.load(.monotonic);
    try zigstd.testExpectEqual(@as(f64, 0.5), try l1.dot(io, f64, 1, &x, big, &y, 1));
    try zigstd.testExpectEqual(@as(f64, 1), try l1.asum(io, f64, 1, &x, big));
    try zigstd.testExpectEqual(@as(isize, 1), try l1.index(io, f64, .largest_magnitude, 1, &x, big));
    // one column at a leading dimension past it
    var c = [_]f64{ 1, 1, 1 };
    try l2.gemv(al, io, f64, .n, 3, 1, 2, &x, big, &y, 1, 1, &c, 1);
    try zigstd.testExpectEqualSlices(f64, &.{ 2, -1, 4 }, &c);
    try l3.gemm(al, io, f64, .n, .n, 3, 1, 1, 1, &x, 3, &y, big, 0, &c, big);
    try zigstd.testExpectEqualSlices(f64, &.{ 0.5, -1, 1.5 }, &c);
    try expectCalled(false, before, "past c_int");
    // an operand whose span is past `c_int`, though each argument fits
    const half: isize = std.math.maxInt(c_int) / 2 + 1;
    try zigstd.testExpectTrue(call(.cblas_ddot, &accel.loaded.?.blas, &.{ span(3, half), 3 }, .{ @as(isize, 3), &x, half, &y, @as(isize, 1) }) == null);
    try zigstd.testExpectEqual(@as(f64, -10.5), call(.cblas_ddot, &accel.loaded.?.blas, &.{ 3, 3 }, .{ @as(isize, 3), &x, @as(isize, 1), &y, @as(isize, 1) }).?);
}
