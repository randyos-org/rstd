//! Tests of `../../blas3/trmm/kernel.zig`'s vector diagonal block (#719).

const std = @import("std");
const gk = @import("../../blas3/gemm/kernel.zig");
const tk = @import("../../blas3/trmm/kernel.zig");

/// Mixed magnitudes with infinities, NaNs and signed zeros sprinkled in.
fn fill(comptime F: type, rand: std.Random, x: []F) void {
    for (x) |*v| v.* = switch (rand.uintLessThan(u8, 16)) {
        0 => std.math.inf(F),
        1 => -std.math.inf(F),
        2 => std.math.nan(F),
        3 => -0.0,
        4 => 0.0,
        else => std.math.ldexp(rand.float(F) * 2 - 1, rand.intRangeAtMost(i32, -20, 20)),
    };
}

/// The vector diagonal block against the scalar one, bit for bit, on whole
/// tiles for every side and direction and the offsets that cut the block.
fn checkGeometry(comptime F: type, comptime geo: gk.Geometry) !void {
    const K = tk.Kernels(F, geo);
    if (comptime !K.vector_diag) return error.SkipZigTest;
    const mr: usize = K.tile_rows;
    const wn: usize = geo.nr;
    const k: usize = mr + wn + 3;
    const ld = mr + 3;
    var prng: std.Random.DefaultPrng = .init(0x719);
    const rand = prng.random();
    var pa: [k * mr]F = undefined;
    var pb: [k * wn]F = undefined;
    var c0: [ld * wn]F = undefined;
    var c1: [ld * wn]F = undefined;
    for (0..20) |_| {
        fill(F, rand, &pa);
        fill(F, rand, &pb);
        fill(F, rand, &c0);
        c1 = c0;
        const alpha: F = if (rand.boolean()) 1 else std.math.ldexp(rand.float(F) + 0.5, rand.intRangeAtMost(i32, -3, 3));
        for ([_]bool{ true, false }) |left| for ([_]bool{ true, false }) |from| {
            var o: isize = -@as(isize, @intCast(mr));
            while (o <= @as(isize, @intCast(k))) : (o += 1) {
                const width: isize = @intCast(if (left) mr else wn);
                const d0: usize = @intCast(@max(0, @min(o, @as(isize, @intCast(k)))));
                const d1: usize = @intCast(@max(0, @min(o + width, @as(isize, @intCast(k)))));
                K.diagonal(left, from, wn, d0, d1, o, &pa, &pb, &c0, ld, alpha);
                K.diagonalScalar(left, from, false, mr, wn, d0, d1, o, &pa, &pb, &c1, ld, &.{alpha});
                const U = if (F == f32) u32 else u64;
                for (c0, c1) |x, y| try std.testing.expectEqual(@as(U, @bitCast(y)), @as(U, @bitCast(x)));
            }
        };
    }
}

test "trmm's vector diagonal block has the scalar loop's bits, infinities and NaNs too" {
    inline for (.{ f32, f64 }) |F| {
        checkGeometry(F, gk.geometry(F)) catch |e| if (e != error.SkipZigTest) return e;
        try checkGeometry(F, .{ .vl = 4, .mv = 2, .nr = 3 });
    }
}
