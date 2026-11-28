//! The public level-3 routines thread above their thresholds and keep the serial bits.

const std = @import("std");
const zigstd = @import("../../../../../../zigstd.zig");
const blas = @import("../../../__root__.zig");
const driver = @import("../../../blas3/driver.zig");
const test_io = @import("../../io.zig");
const types = @import("../../../types.zig");
const apple_accel = @import("../../../apple_accel.zig");
const check = @import("../../check.zig");
const bitwise = @import("../../bitwise.zig");
const interface = @import("../../../blas3/interface.zig");
const l3 = @import("../../../blas3/threaded.zig");
const th = @import("../threaded.zig");

const Io = std.Io;
const fill = check.fillMixed;
const scalar = check.val;
const sameBits = bitwise.sameBits;
const precisions = th.precisions;

/// Takes the public routines off Accelerate (macOS), so they run rstd's
/// threading; returns whether it was routing, for `apple_accel.route`.
pub fn rstdOnly() bool {
    // after init, which turns routing on
    _ = blas.vtable();
    const was = apple_accel.routing();
    _ = apple_accel.route(false);
    return was;
}

test "the public gemm threads above the threshold and keeps the serial bits" {
    if (comptime !test_io.threaded_tests) return error.SkipZigTest;
    _ = blas.coreName();
    const saved = blas.config.cpu_count;
    defer blas.config.cpu_count = saved;
    const accel = rstdOnly();
    defer _ = apple_accel.route(accel);
    var threaded: test_io.TestIo = .init(std.testing.allocator, .{});
    defer threaded.deinit();
    const al = std.testing.allocator;
    inline for (precisions) |T| {
        // past `m * n * k` of 8 serial units for every type's rule, and of 8
        // `team_grain`s on any system
        const m = 331;
        const n = 317;
        const k = 307;
        const a = try al.alloc(T, m * k);
        defer al.free(a);
        const b = try al.alloc(T, k * n);
        defer al.free(b);
        const want = try al.alloc(T, m * n);
        defer al.free(want);
        const got = try al.alloc(T, m * n);
        defer al.free(got);
        fill(T, a, 5);
        fill(T, b, 6);
        fill(T, want, 7);
        @memcpy(got, want);
        const alpha = scalar(T, 0.5, 2);
        const beta = scalar(T, -1.5, 0.25);
        try driver.gemm(al, Io.failing, T, blas.vtable(), true, false, .n, m, n, k, alpha, a.ptr, k, b.ptr, k, beta, want.ptr, m, 1);
        for ([_]u32{ 1, 3, 8 }) |cpus| {
            blas.config.cpu_count = cpus;
            try zigstd.testExpectEqual(cpus, l3.gemmWorkers(T, m, n, k));
            const c = try al.dupe(T, got);
            defer al.free(c);
            try interface.gemm(al, threaded.io(), T, .t, .n, m, n, k, alpha, a.ptr, k, b.ptr, k, beta, c.ptr, m);
            try zigstd.testExpectTrue(sameBits(T, want, c));
        }
    }
}

test "the public symm, herk, syr2k, gemmt and trmm thread above their thresholds and keep the serial bits" {
    if (comptime !test_io.threaded_tests) return error.SkipZigTest;
    // selecting the core fills `blas.config` once; override after
    _ = blas.coreName();
    const saved = blas.config.cpu_count;
    defer blas.config.cpu_count = saved;
    blas.config.cpu_count = 5;
    const accel = rstdOnly();
    defer _ = apple_accel.route(accel);
    var threaded: test_io.TestIo = .init(std.testing.allocator, .{});
    defer threaded.deinit();
    const io = threaded.io();
    const al = std.testing.allocator;
    const T = types.Complex64;
    // five `team_grain`s or more for each call on any system
    const n = 241;
    const k = 229;
    const a = try al.alloc(T, n * n);
    defer al.free(a);
    const b = try al.alloc(T, n * n);
    defer al.free(b);
    const want = try al.alloc(T, n * n);
    defer al.free(want);
    const got = try al.alloc(T, n * n);
    defer al.free(got);
    fill(T, a, 21);
    fill(T, b, 22);
    fill(T, want, 23);
    @memcpy(got, want);
    const vt = blas.vtable();
    const alpha = scalar(T, 0.5, -1);
    try zigstd.testExpectEqual(@as(u32, 5), l3.symmWorkers(T, n, n, n));
    try driver.symm(al, Io.failing, T, vt, false, false, true, n, n, alpha, a.ptr, n, b.ptr, n, alpha, want.ptr, n, 1);
    try interface.symm(al, io, T, .r, .u, n, n, alpha, a.ptr, n, b.ptr, n, alpha, got.ptr, n);
    try zigstd.testExpectTrue(sameBits(T, want, got));
    try zigstd.testExpectEqual(@as(u32, 5), l3.syrkWorkers(T, false, n, k));
    try driver.syrk(al, Io.failing, T, vt, true, false, true, n, k, scalar(T, 0.25, 0), a.ptr, k, null, 0, scalar(T, 2, 0), want.ptr, n, 1);
    try interface.herk(al, io, T, .l, .c, n, k, 0.25, a.ptr, k, 2, got.ptr, n);
    try zigstd.testExpectTrue(sameBits(T, want, got));
    try zigstd.testExpectEqual(@as(u32, 5), l3.syrkWorkers(T, true, n, k));
    try driver.syrk(al, Io.failing, T, vt, false, true, false, n, k, alpha, a.ptr, n, b.ptr, n, alpha, want.ptr, n, 1);
    try interface.syr2k(al, io, T, .u, .n, n, k, alpha, a.ptr, n, b.ptr, n, alpha, got.ptr, n);
    try zigstd.testExpectTrue(sameBits(T, want, got));
    try zigstd.testExpectEqual(@as(u32, 5), l3.gemmtWorkers(T, n, k));
    try driver.gemmt(al, Io.failing, T, vt, false, true, false, .of(true, false), n, k, alpha, a.ptr, k, b.ptr, k, alpha, want.ptr, n, 1);
    try interface.gemmt(al, io, T, .l, .c, .n, n, k, alpha, a.ptr, k, b.ptr, k, alpha, got.ptr, n);
    try zigstd.testExpectTrue(sameBits(T, want, got));
    try zigstd.testExpectEqual(@as(u32, 5), l3.triWorkers(T, k, n, n));
    try driver.trmm(al, Io.failing, T, vt, false, false, true, true, true, k, n, alpha, a.ptr, n, want.ptr, n, 1);
    try interface.trmm(al, io, T, .r, .l, .c, .u, k, n, alpha, a.ptr, n, got.ptr, n);
    try zigstd.testExpectTrue(sameBits(T, want, got));
}
