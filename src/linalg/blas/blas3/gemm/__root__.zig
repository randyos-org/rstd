//! The gemm blocking parameters and Zig kernels (`params.zig`, `kernel.zig`),
//! and their tests.

const types = @import("../../types.zig");

const BlasLong = types.BlasLong;

pub const params = @import("params.zig");
pub const kernel = @import("kernel.zig");

/// Points the level-3 slots at the portable kernels (`kernel.zig`).
pub const installGeneric = kernel.installGeneric;

const check = @import("../../_test/blas3/check.zig");
const Vtable = @import("../../Vtable.zig");

fn cases(vt: *const Vtable) !void {
    try check.routineAll(vt, .gemm);
}

const std = @import("std");
const sysinfo = @import("builtin");
const driver = @import("../driver.zig");
const portable = @import("../../portable.zig");
const blas1check = @import("../../_test/blas1/check.zig");

test "portable gemm kernels match the wide reference through the driver" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try cases(&vt);
}

test "every linked core's gemm kernels match the wide reference, where this machine runs them" {
    try blas1check.eachCore(cases);
}

test "every routine over tiles other targets choose: NEON's, 8 x 3 complex" {
    // natively only: the rigs run their own tiles, and the freestanding ones
    // have no room for four more sets of kernels
    if (comptime sysinfo.target.os.tag == .freestanding) return error.SkipZigTest;
    if (blas1check.rig()) return error.SkipZigTest;
    var vt: Vtable = .{};
    portable.install(&vt);
    kernel.installTile(&vt, "s", f32, .{ .vl = 4, .mv = 4, .nr = 6 });
    kernel.installTile(&vt, "d", f64, .{ .vl = 2, .mv = 4, .nr = 6 });
    kernel.installTile(&vt, "c", std.math.Complex(f32), .{ .vl = 2, .mv = 4, .nr = 3 });
    kernel.installTile(&vt, "z", std.math.Complex(f64), .{ .vl = 1, .mv = 4, .nr = 3 });
    try check.coreCases(&vt);
}

test "the unrolls are the installed kernels' own after a core's blocking" {
    const K = kernel.Kernels(f32, kernel.geometry(f32));
    var vt: Vtable = .{};
    portable.install(&vt);
    // a core whose object is not linked keeps the portable kernels
    params.installCore(&vt, .x86_64, .skylakex, .{});
    portable.finish(&vt);
    try std.testing.expectEqual(K.unroll_m, vt.sgemm_unroll_m);
    try std.testing.expectEqual(K.unroll_n, vt.sgemm_unroll_n);
    try std.testing.expectEqual(@max(K.unroll_m, K.unroll_n), vt.sgemm_unroll_mn);
}

test "negative or zero blocking falls back to the generic sizes" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    vt.dgemm_p = -7;
    vt.dgemm_q = 0;
    vt.dgemm_r = -100;
    const blk = driver.Blocking.of(f64, &vt);
    try std.testing.expect(blk.p > 0 and blk.q == 128 and blk.r == 4096);
    try check.routineCases(f64, &vt, .gemm);
}

test "a packing buffer on the stack is poisoned in tests: a copy that writes nothing shows" {
    const Nop = struct {
        fn copy(_: BlasLong, _: BlasLong, _: [*]f64, _: BlasLong, _: [*]f64) callconv(.c) i32 {
            return 0;
        }
    };
    var vt: Vtable = .{};
    installGeneric(&vt);
    vt.dgemm_itcopy = &Nop.copy;
    var a: [4]f64 = @splat(1);
    var c: [4]f64 = @splat(0);
    try driver.gemm(std.testing.failing_allocator, check.io, f64, &vt, false, false, .n, 2, 2, 2, 1, &a, 2, &a, 2, 0, &c, 2, 1);
    for (c) |v| try std.testing.expect(std.math.isNan(v));
}

test {
    _ = kernel;
    _ = params;
    _ = @import("small.zig");
}
