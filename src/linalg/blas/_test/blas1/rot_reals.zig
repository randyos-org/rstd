//! `rot`'s route to the real kernel (`blas1/rot_reals.zig`): where it
//! applies, what it gives, and that no worker count changes its bits.

const std = @import("std");
const sysinfo = @import("builtin");
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
const rot_reals = @import("../../blas1/rot_reals.zig");

const Io = std.Io;
const Float = types.Float;
const raw = util.raw;
const quantum = interface.quantum;

const x86_64 = sysinfo.cpu.arch == .x86_64;
const complex_types = .{ types.Complex32, types.Complex64 };

/// Lengths around the threshold, in whole chunks and not.
const lengths = [_]usize{ rot_reals.n_min - 1, rot_reals.n_min, rot_reals.n_min + 1, rot_reals.n_min + quantum + 5, 3 * rot_reals.n_min + 17 };
const strides = [_][2]isize{ .{ 1, 1 }, .{ 1, 2 }, .{ 2, 1 }, .{ -1, -1 }, .{ -1, 1 }, .{ 1, -1 } };

const c = 0.6;
const s = -0.8;

/// `x` and `y` rotated by one call of the kernel the reference's rule
/// names, spelled out here: the real one on `2n` reals from 4096 at unit
/// strides on x86_64, else the complex one.
fn kernelRot(comptime T: type, vt: *const Vtable, n: usize, inc: [2]isize, x: []T, y: []T) bool {
    const F = Float(T);
    const ni: isize = @intCast(n);
    const real = x86_64 and n >= 4096 and inc[0] == 1 and inc[1] == 1;
    if (real) {
        _ = portable.kernel(F, vt, .rot_k)(2 * ni, raw(T, x.ptr), 1, raw(T, y.ptr), 1, c, s);
    } else {
        _ = portable.kernel(T, vt, .rot_k)(ni, raw(T, util.first(T, x.ptr, ni, inc[0])), inc[0], raw(T, util.first(T, y.ptr, ni, inc[1])), inc[1], c, s);
    }
    return real;
}

fn routeCase(comptime T: type, al: std.mem.Allocator, vt: *const Vtable, tio: ?Io, n: usize, inc: [2]isize) !void {
    const x0 = try al.alloc(T, n * @abs(inc[0]));
    defer al.free(x0);
    const y0 = try al.alloc(T, n * @abs(inc[1]));
    defer al.free(y0);
    check.fill(T, x0, 31 + n);
    check.fill(T, y0, 32 + n);
    const bufs = try al.alloc(T, 2 * (x0.len + y0.len));
    defer al.free(bufs);
    const want_x = bufs[0..x0.len];
    const want_y = bufs[x0.len..][0..y0.len];
    const x = bufs[x0.len + y0.len ..][0..x0.len];
    const y = bufs[2 * x0.len + y0.len ..];
    @memcpy(want_x, x0);
    @memcpy(want_y, y0);
    const real = kernelRot(T, vt, n, inc, want_x, want_y);
    try std.testing.expectEqual(real, rot_reals.routed(T, @intCast(n), inc[0], inc[1]));
    if (real) {
        // the complex kernel's result, to rounding
        @memcpy(x, x0);
        @memcpy(y, y0);
        _ = portable.kernel(T, vt, .rot_k)(@intCast(n), raw(T, x.ptr), 1, raw(T, y.ptr), 1, c, s);
        const tol = check.tolerance(Float(T), 2);
        try check.expectCloseSlices(T, x, want_x, tol);
        try check.expectCloseSlices(T, y, want_y, tol);
    }
    // one worker in one call, then chunks shorter than the threshold
    for (0..10) |workers| {
        if (workers > 1 and tio == null) break;
        blas.config.cpu_count = @intCast(@max(workers, 1));
        threading.test_serial_max = if (workers == 0) null else quantum;
        @memcpy(x, x0);
        @memcpy(y, y0);
        try interface.rot(if (workers > 1) tio.? else Io.failing, T, @intCast(n), x.ptr, inc[0], y.ptr, inc[1], c, s);
        if (!bitwise.sameBits(T, want_x, x) or !bitwise.sameBits(T, want_y, y)) {
            check.say("{d} workers: ", .{workers});
            return error.TestExpectedEqual;
        }
    }
}

fn routeCore(vt: *const Vtable) !void {
    const saved = blas.testUse(vt);
    defer blas.testRestore(saved);
    const al = std.testing.allocator;
    var threaded: if (test_io.threaded_tests) test_io.TestIo else void = if (comptime test_io.threaded_tests) .init(al, .{}) else {};
    defer if (comptime test_io.threaded_tests) threaded.deinit();
    const tio: ?Io = if (comptime test_io.threaded_tests) threaded.io() else null;
    inline for (complex_types) |T| for (lengths) |n| for (strides) |inc| routeCase(T, al, vt, tio, n, inc) catch |e| {
        check.say("{s} n={d} inc={d},{d}: ", .{ @typeName(T), n, inc[0], inc[1] });
        return e;
    };
}

test "interface: complex rot takes the real kernel from n = 4096 at unit strides on x86_64 only, the same bits at 1 to 9 workers, on every core" {
    _ = blas.coreName();
    const cpus = blas.config.cpu_count;
    defer blas.config.cpu_count = cpus;
    defer threading.test_serial_max = null;
    try check.eachCore(routeCore);
}

test "rot_reals.routed: the reference's bounds, unit strides, single and double complex on x86_64" {
    inline for (complex_types) |T| {
        try std.testing.expectEqual(x86_64, rot_reals.serves(T));
        try std.testing.expect(!rot_reals.routed(T, 4095, 1, 1));
        try std.testing.expectEqual(x86_64, rot_reals.routed(T, 4096, 1, 1));
        try std.testing.expectEqual(x86_64, rot_reals.routed(T, 1 << 29, 1, 1));
        if (@bitSizeOf(isize) > 32) try std.testing.expect(!rot_reals.routed(T, (1 << 29) + 1, 1, 1));
        for ([_][2]isize{ .{ 1, 2 }, .{ 2, 1 }, .{ -1, -1 }, .{ 1, -1 }, .{ 0, 1 }, .{ 1, 0 } }) |inc| try std.testing.expect(!rot_reals.routed(T, 1 << 20, inc[0], inc[1]));
    }
    // real rot is its own kernel; the reference leaves `XDOUBLE` alone, where it is a type of its own
    inline for (.{ f32, f64 } ++ check.x_real ++ check.x_complex) |T| {
        try std.testing.expect(!rot_reals.serves(T));
        try std.testing.expect(!rot_reals.routed(T, 1 << 20, 1, 1));
    }
}
