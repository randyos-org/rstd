//! Routing: which public calls reach Accelerate, the off switch, and one Accelerate thread per caller.

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("../../types.zig");
const blas = @import("../../__root__.zig");
const l3 = @import("../../blas3/core.zig");
const l2 = @import("../../blas2/core.zig");
const l1 = @import("../../blas1/interface.zig");
const accel = @import("../../apple_accel.zig");
const zigstd = @import("../../../../../zigstd.zig");
const check = @import("../check.zig");

const Complex32 = types.Complex32;
const Complex64 = types.Complex64;
const isComplex = types.isComplex;
const Float = types.Float;
const Threading = accel.Threading;
const available = accel.available;
const calls = &accel.calls;
const route = accel.route;
const routing = accel.routing;
const x86_64 = accel.x86_64;
const sent = accel.sent;
const routes = accel.routes;

// ---- routing tests -------------------------------------------------------------

/// Fails unless Accelerate was called since `before` exactly when `want`.
pub fn expectCalled(want: bool, before: usize, what: []const u8) !void {
    if ((calls.load(.monotonic) > before) == want) return;
    check.say("{s}: Accelerate {s}called: ", .{ what, if (want) "not " else "" });
    return error.TestUnexpectedResult;
}

/// Counts calls into Accelerate routine by routine.
pub const Probe = struct {
    on: bool,
    last: usize,

    pub fn init(on_: bool) Probe {
        return .{ .on = on_, .last = calls.load(.monotonic) };
    }

    /// The routine just run called Accelerate if `routed` and `on`.
    pub fn saw(p: *Probe, routed: bool, what: []const u8) !void {
        try expectCalled(p.on and routed, p.last, what);
        p.last = calls.load(.monotonic);
    }
};

pub fn val(comptime T: type, re: Float(T), im: Float(T)) T {
    return if (comptime isComplex(T)) .init(re, im) else re;
}

/// Each public routine over `T` once, on small finite data: each one
/// Accelerate has calls it when `on`, and nothing else does.
pub fn everyRoutine(comptime T: type, on_: bool) !void {
    const al = zigstd.testAllocator;
    const io = check.io;
    const cx = comptime isComplex(T);
    const real = !cx;
    const al1 = val(T, 1.5, -0.5);
    const be1 = val(T, 0.5, 0.25);
    var a: [36]T = undefined;
    var bm: [36]T = undefined;
    var c: [36]T = undefined;
    var x: [6]T = undefined;
    var y: [6]T = undefined;
    check.fill(T, &a, 1);
    check.fill(T, &bm, 2);
    check.fill(T, &c, 3);
    check.fill(T, &x, 4);
    check.fill(T, &y, 5);
    // dominant diagonals, for the solves
    for (0..6) |i| a[i * 7] = val(T, 4, 0);
    for (0..6) |i| a[i] = val(T, 4, 0);
    const n: isize = 6;
    var p: Probe = .init(on_);

    try l1.copy(io, T, n, &x, 1, &y, 1);
    try p.saw(sent(.copy, T, true), "copy");
    try l1.swap(io, T, n, &x, 1, &y, 1);
    try p.saw(sent(.swap, T, true), "swap");
    // strided calls go by `sends`; see the test of it
    try l1.swap(io, T, n, &x, 1, &y, -1);
    try p.saw(x86_64 and sent(.swap, T, false), "swap, a negative stride");
    try l1.scal(io, T, n, al1, &x, 1);
    try p.saw(sent(.scal, T, true), "scal");
    _ = try l1.dot(io, T, n, &x, 1, &y, 1);
    // `.plain`'s complex `f32` sums in rstd's `f64` slot (issue #584)
    try p.saw(T != Complex32 and (T != f32 or (x86_64 and !accel.rosetta())) and sent(.dot, T, true), "dot");
    if (cx) {
        _ = try l1.dotc(io, T, n, &x, 1, &y, 1);
        try p.saw(T != Complex32 and sent(.dotc, T, true), "dotc");
        try l1.axpyc(io, T, n, al1, &x, 1, &y, 1);
        try p.saw(false, "axpyc");
        // `alpha` has a zero imaginary part, a quirk on arm64 only
        try l1.scalReal(io, T, n, 0.5, &x, 1);
        try p.saw(x86_64 and sent(.scal, T, true), "scalReal");
    } else {
        var param = [5]T{ -1, 0.5, 0.25, -0.25, 0.75 };
        try l1.rotm(io, T, n, &x, 1, &y, 1, &param);
        try p.saw(sent(.rotm, T, true), "rotm");
    }
    if (T == f32) {
        _ = try l1.dsdot(io, n, &x, 1, &y, 1);
        try p.saw(!accel.rosetta() and sent(.dsdot, T, true), "dsdot");
        _ = try l1.sdsdot(io, n, 0.5, &x, 1, &y, 1);
        try p.saw(!accel.rosetta() and sent(.dsdot, T, true), "sdsdot");
    }
    _ = try l1.nrm2(io, T, n, &x, 1);
    try p.saw(x86_64 and (Float(T) == f64 or !accel.rosetta()) and sent(.nrm2, T, true), "nrm2");
    _ = try l1.nrm2(io, T, n, &x, -1);
    try p.saw(false, "nrm2, a negative stride");
    _ = try l1.asum(io, T, n, &x, 1);
    try p.saw(Float(T) != f32 and sent(.asum, T, true), "asum");
    _ = try l1.sum(io, T, n, &x, 1);
    try p.saw(false, "sum");
    _ = try l1.index(io, T, .largest_magnitude, n, &x, 1);
    try p.saw(false, "index amax");
    _ = try l1.index(io, T, .smallest_magnitude, n, &x, 1);
    try p.saw(false, "index amin");
    _ = try l1.extreme(io, T, .largest_magnitude, n, &x, 1);
    try p.saw(false, "extreme amax");
    try l1.rot(io, T, n, &x, 1, &y, 1, 0.6, 0.8);
    try p.saw(sent(.rot, T, true), "rot");
    try l1.axpy(io, T, n, al1, &x, 1, &y, 1);
    try p.saw(sent(.axpy, T, true), "axpy");
    try l1.axpby(io, T, n, al1, &x, 1, be1, &y, 1);
    try p.saw(sent(.axpby, T, true), "axpby");

    try l2.gemv(al, io, T, .n, n, n, al1, &a, n, &x, 1, be1, &y, 1);
    try p.saw(routes(.gemv, T), "gemv");
    try l2.gemv(al, io, T, if (cx) .c else .t, n, n, al1, &a, n, &x, 1, be1, &y, 1);
    try p.saw(routes(.gemv, T), "gemv transposed");
    try l2.gemv(al, io, T, .r, n, n, al1, &a, n, &x, 1, be1, &y, 1);
    try p.saw(real and routes(.gemv, T), "gemv r");
    try l2.gbmv(al, io, T, .t, n, n, 1, 1, al1, &a, 3, &x, 1, be1, &y, 1);
    try p.saw(routes(.gbmv, T), "gbmv");
    try l2.ger(al, io, T, n, n, al1, &x, 1, &y, 1, &c, n);
    try p.saw(routes(.ger, T), "ger");
    if (cx) {
        try l2.gerc(al, io, T, n, n, al1, &x, 1, &y, 1, &c, n);
        try p.saw(routes(.ger, T), "gerc");
        try l2.hemv(al, io, T, .l, n, al1, &a, n, &x, 1, be1, &y, 1);
        try p.saw(true, "hemv");
        try l2.hemv(al, io, T, .v, n, al1, &a, n, &x, 1, be1, &y, 1);
        try p.saw(false, "hemv v");
        try l2.hbmv(al, io, T, .u, n, 1, al1, &a, 2, &x, 1, be1, &y, 1);
        try p.saw(true, "hbmv");
        try l2.hpmv(al, io, T, .l, n, al1, &c, &x, 1, be1, &y, 1);
        try p.saw(true, "hpmv");
        try l2.her(al, io, T, .u, n, 0.5, &x, 1, &c, n);
        try p.saw(true, "her");
        try l2.hpr(al, io, T, .l, n, 0.5, &x, 1, &c);
        try p.saw(true, "hpr");
        try l2.her2(al, io, T, .u, n, al1, &x, 1, &y, 1, &c, n);
        try p.saw(true, "her2");
        try l2.hpr2(al, io, T, .l, n, al1, &x, 1, &y, 1, &c);
        try p.saw(true, "hpr2");
    }
    try l2.symv(al, io, T, .u, n, al1, &a, n, &x, 1, be1, &y, 1);
    try p.saw(real and routes(.symv, T), "symv");
    try l2.sbmv(al, io, T, .l, n, 1, al1, &a, 2, &x, 1, be1, &y, 1);
    try p.saw(real, "sbmv");
    try l2.spmv(al, io, T, .u, n, al1, &c, &x, 1, be1, &y, 1);
    try p.saw(real, "spmv");
    try l2.syr(al, io, T, .u, n, al1, &x, 1, &c, n);
    try p.saw(real, "syr");
    try l2.spr(al, io, T, .l, n, al1, &x, 1, &c);
    try p.saw(real, "spr");
    try l2.syr2(al, io, T, .u, n, al1, &x, 1, &y, 1, &c, n);
    try p.saw(real, "syr2");
    try l2.spr2(al, io, T, .l, n, al1, &x, 1, &y, 1, &c);
    try p.saw(real, "spr2");
    try l2.trmv(al, io, T, .u, .t, .n, n, &a, n, &x, 1);
    try p.saw(true, "trmv");
    try l2.trmv(al, io, T, .u, .r, .n, n, &a, n, &x, 1);
    try p.saw(real, "trmv r");
    try l2.trsv(al, io, T, .l, .n, .u, n, &a, n, &x, -1);
    try p.saw(routes(.trsv, T), "trsv");
    try l2.tbmv(al, io, T, .u, .n, .n, n, 1, &a, 2, &x, 1);
    try p.saw(true, "tbmv");
    try l2.tbsv(al, io, T, .l, .c, .n, n, 1, &a, 2, &x, 1);
    try p.saw(true, "tbsv");
    try l2.tpmv(al, io, T, .u, .n, .u, n, &a, &x, 1);
    try p.saw(true, "tpmv");
    try l2.tpsv(al, io, T, .l, .n, .u, n, &a, &x, 1);
    try p.saw(true, "tpsv");

    try l3.gemm(al, io, T, .n, if (cx) .c else .t, n, n, n, al1, &a, n, &bm, n, be1, &c, n);
    try p.saw(routes(.gemm, T), "gemm");
    try l3.gemm(al, io, T, .r, .n, n, n, n, al1, &a, n, &bm, n, be1, &c, n);
    try p.saw(real and routes(.gemm, T), "gemm r");
    try l3.symm(al, io, T, .l, .u, n, n, al1, &a, n, &bm, n, be1, &c, n);
    try p.saw(routes(.symm, T), "symm");
    try l3.syrk(al, io, T, .u, .n, n, n, al1, &a, n, be1, &c, n);
    try p.saw(routes(.syrk, T), "syrk");
    try l3.syr2k(al, io, T, .l, .t, n, n, al1, &a, n, &bm, n, be1, &c, n);
    try p.saw(true, "syr2k");
    if (cx) {
        try l3.hemm(al, io, T, .r, .l, n, n, al1, &a, n, &bm, n, be1, &c, n);
        try p.saw(true, "hemm");
        try l3.herk(al, io, T, .l, .c, n, n, 1.5, &a, n, 0.5, &c, n);
        try p.saw(true, "herk");
        try l3.her2k(al, io, T, .u, .n, n, n, al1, &a, n, &bm, n, 0.5, &c, n);
        try p.saw(true, "her2k");
    }
    try l3.trmm(al, io, T, .r, .u, .t, .n, n, n, al1, &a, n, &c, n);
    try p.saw(routes(.trmm, T), "trmm");
    try l3.trsm(al, io, T, .l, .l, .n, .u, n, n, al1, &a, n, &c, n);
    try p.saw(routes(.trsm, T), "trsm");
}

test "the public routines call Accelerate on macOS where it has them, for s, d, c and z only" {
    _ = blas.vtable();
    const was = routing();
    defer _ = route(was);
    if (!route(true)) return zigstd.testExpectTrue(!available or accel.loaded == null);
    inline for (.{ f32, f64, Complex32, Complex64 }) |T| try everyRoutine(T, true);
    // types it has nothing for
    var h: [4]f16 = @splat(1);
    const before = calls.load(.monotonic);
    _ = try l1.dot(check.io, f16, 4, &h, 1, &h, 1);
    inline for (check.x_real) |X| {
        var q: [4]X = @splat(1);
        _ = try l1.dot(check.io, X, 4, &q, 1, &q, 1);
    }
    try expectCalled(false, before, "f16 and XDouble dot");
}

test "the off switch runs rstd's kernels, and RSTD_BLAS_ACCELERATE=0 sets it at init" {
    _ = blas.vtable();
    const was = routing();
    defer _ = route(was);
    try zigstd.testExpectTrue(!blas.useAccelerate(false));
    inline for (.{ f32, f64, Complex32, Complex64 }) |T| try everyRoutine(T, false);
    if (!blas.useAccelerate(true)) return error.SkipZigTest;
    try everyRoutine(f64, true);

    blas.shutdown();
    defer {
        blas.shutdown();
        _ = blas.vtable();
    }
    var env: std.process.Environ.Map = .init(zigstd.testAllocator);
    defer env.deinit();
    try env.put("RSTD_BLAS_ACCELERATE", "0");
    try zigstd.testExpectTrue(blas.init(&env));
    try zigstd.testExpectTrue(!routing());
    try everyRoutine(Complex32, false);
    // loaded before, and never unloaded
    try zigstd.testExpectTrue(blas.useAccelerate(true));
    try everyRoutine(Complex32, true);
}

pub fn dotOnThread(out: *?Threading, before: usize, called: *bool) void {
    var x = [_]f64{ 1, 2, 3 };
    _ = l1.dot(check.io, f64, 3, &x, 1, &x, 1) catch {};
    called.* = calls.load(.monotonic) > before;
    out.* = accel.loaded.?.threading();
}

test "each thread that calls in runs Accelerate single-threaded" {
    if (sysinfo.single_threaded) return error.SkipZigTest;
    _ = blas.vtable();
    if (!routing()) return error.SkipZigTest;
    var model: ?Threading = .multi_threaded;
    var called = false;
    const th = try std.Thread.spawn(.{}, dotOnThread, .{ &model, calls.load(.monotonic), &called });
    th.join();
    try zigstd.testExpectTrue(called);
    // null before macOS 15, where only `VECLIB_MAXIMUM_THREADS` does it
    if (model) |m| try zigstd.testExpectEqual(Threading.single_threaded, m);
}
