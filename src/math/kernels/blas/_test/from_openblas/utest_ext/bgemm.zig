//! `test_bgemm.c`: `bgemm` against `sbgemm`.

const std = @import("std");
const ob = @import("../__root__.zig");
const blas = @import("../../../__root__.zig");
const l3 = @import("../../../blas3/core.zig");
const util = @import("util.zig");

const testing = std.testing;
const transOf = util.transOf;
const transposedLetter = util.transposedLetter;

// ---- bgemm ---------------------------------------------------------------------

/// One `CTEST` of `test_bgemm.c`'s `check_bgemm`.
const BgemmCase = struct { name: []const u8, ta: u8, tb: u8, m: isize, n: isize, k: isize, alpha: f32, beta: f32, lda: isize, ldb: isize, ldc: isize };

/// `test_bgemm.c`'s `float_to_bf16`: the high half of the bits, truncated.
fn truncBf16(v: f32) blas.types.BFloat16 {
    return @fromBackingInt(@as(u16, @intCast(@as(u32, @bitCast(v)) >> 16)));
}

/// `test_bgemm.c`'s `bf16_to_float`: the bits widened, nothing flushed.
fn bitsF32(v: blas.types.BFloat16) f32 {
    return @bitCast(@as(u32, @backingInt(v)) << 16);
}

/// `check_bgemm`: `bgemm` against `sbgemm` on the same bfloat16 operands and
/// scalars, the largest `|bgemm - sbgemm| / (|sbgemm| + 1)` within `BGEMM_TOL`.
fn bgemmCase(c: BgemmCase, r: std.Random) !void {
    const al = testing.allocator;
    const BFloat16 = blas.types.BFloat16;
    const ta = transposedLetter(c.ta);
    const a_cols: usize = @intCast(if (ta) c.m else c.k);
    const b_cols: usize = @intCast(if (transposedLetter(c.tb)) c.k else c.n);
    const n: usize = @intCast(c.n);
    const m: usize = @intCast(c.m);
    const ldc: usize = @intCast(c.ldc);
    const a = try al.alloc(BFloat16, a_cols * @as(usize, @intCast(c.lda)));
    defer al.free(a);
    const b = try al.alloc(BFloat16, b_cols * @as(usize, @intCast(c.ldb)));
    defer al.free(b);
    const cb = try al.alloc(BFloat16, n * ldc);
    defer al.free(cb);
    const cs = try al.alloc(f32, n * ldc);
    defer al.free(cs);
    for (a) |*v| v.* = truncBf16(0.5 + r.float(f32));
    for (b) |*v| v.* = truncBf16(0.5 + r.float(f32));
    for (cb, cs) |*v, *w| {
        v.* = truncBf16(0.5 + r.float(f32));
        w.* = bitsF32(v.*);
    }
    const alpha = truncBf16(c.alpha);
    const beta = truncBf16(c.beta);
    try l3.gemmHalf(al, ob.io, .bf16_to_f32, transOf(c.ta), transOf(c.tb), c.m, c.n, c.k, bitsF32(alpha), a.ptr, c.lda, b.ptr, c.ldb, bitsF32(beta), cs.ptr, c.ldc);
    try l3.gemmHalf(al, ob.io, .bf16, transOf(c.ta), transOf(c.tb), c.m, c.n, c.k, alpha, a.ptr, c.lda, b.ptr, c.ldb, beta, cb.ptr, c.ldc);
    var worst: f32 = 0;
    for (0..n) |j| for (0..m) |i| {
        const want = cs[j * ldc + i];
        worst = @max(worst, @abs(bitsF32(cb[j * ldc + i]) - want) / (@abs(want) + 1));
    };
    // `BGEMM_TOL`
    if (worst <= 0.01) return;
    ob.say("bgemm {s}: relative difference {d}\n", .{ c.name, worst });
    return error.TestUnexpectedResult;
}

pub fn bgemm() !void {
    var prng = std.Random.DefaultPrng.init(0xB6);
    for (bgemm_cases) |c| try bgemmCase(c, prng.random());
}

test "utest_ext: bgemm against sbgemm on the same bfloat16 operands" {
    try bgemm();
}

/// `test_bgemm.c`.
const bgemm_cases = [_]BgemmCase{
    .{ .name = "nn_M50_N50_K50", .ta = 'N', .tb = 'N', .m = 50, .n = 50, .k = 50, .alpha = 1.5, .beta = 0, .lda = 50, .ldb = 50, .ldc = 50 },
    .{ .name = "nn_M100_N100_K100", .ta = 'N', .tb = 'N', .m = 100, .n = 100, .k = 100, .alpha = 1, .beta = 0, .lda = 100, .ldb = 100, .ldc = 100 },
    .{ .name = "nn_M16_N8_K16", .ta = 'N', .tb = 'N', .m = 16, .n = 8, .k = 16, .alpha = 1, .beta = 0, .lda = 16, .ldb = 16, .ldc = 16 },
    .{ .name = "nn_remainder_M17_N9_K7", .ta = 'N', .tb = 'N', .m = 17, .n = 9, .k = 7, .alpha = 1, .beta = 0, .lda = 17, .ldb = 7, .ldc = 17 },
    .{ .name = "nn_M1_N1_K8", .ta = 'N', .tb = 'N', .m = 1, .n = 1, .k = 8, .alpha = 1, .beta = 0, .lda = 1, .ldb = 8, .ldc = 1 },
    .{ .name = "tn_M50_N50_K50", .ta = 'T', .tb = 'N', .m = 50, .n = 50, .k = 50, .alpha = 1.5, .beta = 0, .lda = 50, .ldb = 50, .ldc = 50 },
    .{ .name = "tn_M100_N50_K50", .ta = 'T', .tb = 'N', .m = 100, .n = 50, .k = 50, .alpha = 1, .beta = 0, .lda = 50, .ldb = 50, .ldc = 100 },
    .{ .name = "nt_M50_N50_K100", .ta = 'N', .tb = 'T', .m = 50, .n = 50, .k = 100, .alpha = 1, .beta = 0, .lda = 50, .ldb = 50, .ldc = 50 },
    .{ .name = "tt_M50_N50_K50", .ta = 'T', .tb = 'T', .m = 50, .n = 50, .k = 50, .alpha = 1.5, .beta = 2, .lda = 50, .ldb = 50, .ldc = 50 },
    .{ .name = "alpha_zero_nn", .ta = 'N', .tb = 'N', .m = 50, .n = 50, .k = 50, .alpha = 0, .beta = 2, .lda = 50, .ldb = 50, .ldc = 50 },
    .{ .name = "beta_zero_nn", .ta = 'N', .tb = 'N', .m = 50, .n = 50, .k = 50, .alpha = 1, .beta = 0, .lda = 50, .ldb = 50, .ldc = 50 },
    .{ .name = "beta_one_nn", .ta = 'N', .tb = 'N', .m = 50, .n = 50, .k = 50, .alpha = 2, .beta = 1, .lda = 50, .ldb = 50, .ldc = 50 },
    .{ .name = "nn_odd_k_M32_N8_K7", .ta = 'N', .tb = 'N', .m = 32, .n = 8, .k = 7, .alpha = 1, .beta = 0, .lda = 32, .ldb = 7, .ldc = 32 },
    .{ .name = "nn_odd_k_M16_N8_K1", .ta = 'N', .tb = 'N', .m = 16, .n = 8, .k = 1, .alpha = 1, .beta = 0, .lda = 16, .ldb = 1, .ldc = 16 },
    .{ .name = "nn_n4_path_M32_N4_K16", .ta = 'N', .tb = 'N', .m = 32, .n = 4, .k = 16, .alpha = 1, .beta = 0, .lda = 32, .ldb = 16, .ldc = 32 },
    .{ .name = "nn_n2_path_M32_N2_K16", .ta = 'N', .tb = 'N', .m = 32, .n = 2, .k = 16, .alpha = 1, .beta = 0, .lda = 32, .ldb = 16, .ldc = 32 },
    .{ .name = "nn_n1_path_M32_N1_K16", .ta = 'N', .tb = 'N', .m = 32, .n = 1, .k = 16, .alpha = 1, .beta = 0, .lda = 32, .ldb = 16, .ldc = 32 },
    .{ .name = "nn_m8_path_M8_N8_K16", .ta = 'N', .tb = 'N', .m = 8, .n = 8, .k = 16, .alpha = 1, .beta = 0, .lda = 8, .ldb = 16, .ldc = 8 },
    .{ .name = "nn_m4_path_M4_N8_K16", .ta = 'N', .tb = 'N', .m = 4, .n = 8, .k = 16, .alpha = 1, .beta = 0, .lda = 4, .ldb = 16, .ldc = 4 },
    .{ .name = "nn_m2_path_M2_N8_K16", .ta = 'N', .tb = 'N', .m = 2, .n = 8, .k = 16, .alpha = 1, .beta = 0, .lda = 2, .ldb = 16, .ldc = 2 },
    .{ .name = "nn_m1_path_M1_N8_K16", .ta = 'N', .tb = 'N', .m = 1, .n = 8, .k = 16, .alpha = 1, .beta = 0, .lda = 1, .ldb = 16, .ldc = 1 },
    .{ .name = "nn_m32_N8_K16", .ta = 'N', .tb = 'N', .m = 32, .n = 8, .k = 16, .alpha = 1, .beta = 0, .lda = 32, .ldb = 16, .ldc = 32 },
    .{ .name = "nn_alpha1p5_beta2_M50_N50_K50", .ta = 'N', .tb = 'N', .m = 50, .n = 50, .k = 50, .alpha = 1.5, .beta = 2, .lda = 50, .ldb = 50, .ldc = 50 },
};
