//! Tests of `../../blas3/trsm/kernel.zig`'s vector solves (#720).

const std = @import("std");
const gk = @import("../../blas3/gemm/kernel.zig");
const tk = @import("../../blas3/trsm/kernel.zig");
const fill = @import("trmm_diagonal.zig").fill;

/// Equal bits, a NaN for a NaN: its sign and payload are LLVM's to choose.
fn expectBits(comptime F: type, want: []const F, got: []const F) !void {
    const U = if (F == f32) u32 else u64;
    for (want, got) |x, y| {
        if (std.math.isNan(x) and std.math.isNan(y)) continue;
        try std.testing.expectEqual(@as(U, @bitCast(x)), @as(U, @bitCast(y)));
    }
}

/// The vector solves against the scalar ones, bit for bit, on whole tiles
/// of every width up to `nr`, both directions, unit or not.
fn checkGeometry(comptime F: type, comptime geo: gk.Geometry) !void {
    const K = tk.Kernels(F, geo);
    if (comptime !K.vector_solve) return error.SkipZigTest;
    const mr: usize = K.tile_rows;
    const nr: usize = geo.nr;
    const ld = mr + 3;
    var prng: std.Random.DefaultPrng = .init(0x720);
    const rand = prng.random();
    var tri: [mr * mr]F = undefined;
    var c0: [ld * nr]F = undefined;
    var c1: [ld * nr]F = undefined;
    var p0: [mr * nr]F = undefined;
    var p1: [mr * nr]F = undefined;
    var right: [nr * nr]F = undefined;
    for (0..20) |_| for (1..nr + 1) |n| for ([_]bool{ true, false }) |forward| for ([_]bool{ true, false }) |unit| {
        fill(F, rand, &tri);
        fill(F, rand, &right);
        fill(F, rand, &c0);
        c1 = c0;
        p0 = @splat(0);
        p1 = p0;
        K.solveLeft(forward, false, unit, mr, n, &tri, &p0, &c0, ld);
        K.solveLeftScalar(forward, false, unit, mr, n, &tri, &p1, &c1, ld);
        try expectBits(F, &c1, &c0);
        try expectBits(F, p1[0 .. mr * n], p0[0 .. mr * n]);

        fill(F, rand, &c0);
        c1 = c0;
        K.solveRight(forward, false, unit, mr, n, &p0, &right, &c0, ld);
        K.solveRightScalar(forward, false, unit, mr, n, &p1, &right, &c1, ld);
        try expectBits(F, &c1, &c0);
        try expectBits(F, p1[0 .. mr * n], p0[0 .. mr * n]);
    };
}

test "trsm's vector solves have the scalar solves' bits" {
    inline for (.{ f32, f64 }) |F| {
        checkGeometry(F, gk.geometry(F)) catch |e| if (e != error.SkipZigTest) return e;
        try checkGeometry(F, .{ .vl = 4, .mv = 2, .nr = 3 });
    }
}
