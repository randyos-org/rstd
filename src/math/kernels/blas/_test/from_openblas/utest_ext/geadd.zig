//! `test_?geadd.c`: `geadd` against the plain sum.

const std = @import("std");
const ob = @import("../__root__.zig");
const ext = @import("../../../ext/interface.zig");
const util = @import("util.zig");

const testing = std.testing;
const isComplex = util.isComplex;
const C32 = util.C32;
const C64 = util.C64;
const single_eps = util.single_eps;
const double_eps = util.double_eps;
const tol = util.tol;
const Case = util.Case;
const val = util.val;
const fill = util.fill;
const expectNorm = util.expectNorm;
const transOf = util.transOf;
const transposed = util.transposed;
const matDiff = util.matDiff;
const transposedLetter = util.transposedLetter;

/// `geadd` against `geaddRef` for one case of `test_?geadd.c`, or its bad
/// argument.
fn geaddCase(comptime T: type, c: Case, r: std.Random) !void {
    const al = testing.allocator;
    const ta = transposedLetter(c.transa);
    const tc = transposedLetter(c.transb);
    if (c.kind == .bad) {
        var a = [_]T{val(T, .{ 1, 1 })};
        var x = [_]T{val(T, .{ 1, 1 })};
        const want: anyerror = switch (c.info) {
            1 => error.InvalidM,
            2 => error.InvalidN,
            5 => error.InvalidLda,
            else => error.InvalidLdc,
        };
        return ob.expectError(want, ext.geadd(ob.io, T, c.m, c.n, val(T, .{ 1, 1 }), &a, c.lda, val(T, .{ 1, 1 }), &x, c.ldc, transOf(c.transa), transOf(c.transb)), c.name);
    }
    const m: usize = @intCast(c.m);
    const n: usize = @intCast(c.n);
    const lda: usize = @intCast(c.lda);
    const ldc: usize = @intCast(c.ldc);
    // `op(A)` and `op(C)` are `m x n`; stored transposed, `n x m`
    const a = try al.alloc(T, lda * @max(if (ta) m else n, 1));
    defer al.free(a);
    const cm = try al.alloc(T, ldc * @max(if (tc) m else n, 1));
    defer al.free(cm);
    fill(T, a, r);
    fill(T, cm, r);
    const want = try al.dupe(T, cm);
    defer al.free(want);
    const alpha = val(T, c.alpha);
    const beta = val(T, c.beta);
    for (0..n) |j| for (0..m) |i| {
        const k = if (tc) j + i * ldc else i + j * ldc;
        const x = if (ta) a[j + i * lda] else a[i + j * lda];
        want[k] = ob.add(ob.mul(alpha, x), ob.mul(beta, want[k]));
    };
    try ext.geadd(ob.io, T, c.m, c.n, alpha, a.ptr, c.lda, beta, cm.ptr, c.ldc, transOf(c.transa), transOf(c.transb));
    try expectNorm(c.name, T, matDiff(T, if (tc) n else m, if (tc) m else n, cm, want, ldc), tol(T, single_eps, double_eps));
}

/// The two `c_api_matrix_2x2_*` cases of `test_?geadd.c`: `A^T` into `C`,
/// and into `C^T`, which is `A` again.
fn geaddFixed(comptime T: type) !void {
    const complex = comptime isComplex(T);
    var a: [4]T = undefined;
    for (&a, 0..) |*v, i| {
        const f: f64 = @floatFromInt(i);
        v.* = val(T, if (complex) .{ 2 * f + 1, 2 * f + 2 } else .{ f + 1, 0 });
    }
    inline for (.{ .{ .n, [_]usize{ 0, 2, 1, 3 } }, .{ .t, [_]usize{ 0, 1, 2, 3 } } }) |p| {
        var c: [4]T = @splat(ob.zero(T));
        try ext.geadd(ob.io, T, 2, 2, val(T, .{ 1, 0 }), &a, 2, ob.zero(T), &c, 2, .t, p[0]);
        for (c, p[1]) |got, from| try ob.expect(ob.eql(got, a[from]), "c_api_matrix_2x2_transA");
    }
}

pub fn geadd() !void {
    var prng = std.Random.DefaultPrng.init(0xC0);
    inline for (.{ .{ f32, &geadd_s }, .{ f64, &geadd_d }, .{ C32, &geadd_c }, .{ C64, &geadd_z } }) |p| {
        for (p[1]) |c| try geaddCase(p[0], c, prng.random());
        try geaddFixed(p[0]);
    }
}

test "utest_ext: geadd against the plain sum, and its bad arguments" {
    try geadd();
}

/// `test_sgeadd.c`, column-major.
const geadd_s = [_]Case{
    .{ .name = "matrix_n_100_m_100", .kind = .run, .transa = 'N', .transb = 'N', .m = 100, .n = 100, .lda = 100, .ldc = 100, .alpha = .{ 3.0, 0 }, .beta = .{ 3.0, 0 } },
    .{ .name = "matrix_n_100_m_100_alpha_zero", .kind = .run, .transa = 'N', .transb = 'N', .m = 100, .n = 100, .lda = 100, .ldc = 100, .alpha = .{ 0.0, 0 }, .beta = .{ 2.5, 0 } },
    .{ .name = "matrix_n_100_m_100_beta_zero", .kind = .run, .transa = 'N', .transb = 'N', .m = 100, .n = 100, .lda = 100, .ldc = 100, .alpha = .{ 3.0, 0 }, .beta = .{ 0.0, 0 } },
    .{ .name = "matrix_n_100_m_100_alpha_beta_zero", .kind = .run, .transa = 'N', .transb = 'N', .m = 100, .n = 100, .lda = 100, .ldc = 100, .alpha = .{ 0.0, 0 }, .beta = .{ 0.0, 0 } },
    .{ .name = "matrix_n_100_m_50", .kind = .run, .transa = 'N', .transb = 'N', .m = 50, .n = 100, .lda = 50, .ldc = 50, .alpha = .{ 1.0, 0 }, .beta = .{ 1.0, 0 } },
    .{ .name = "xerbla_n_invalid", .kind = .bad, .transa = 'N', .transb = 'N', .m = 1, .n = -1, .lda = 1, .ldc = 1, .info = 2 },
    .{ .name = "xerbla_m_invalid", .kind = .bad, .transa = 'N', .transb = 'N', .m = -1, .n = 1, .lda = 1, .ldc = 1, .info = 1 },
    .{ .name = "xerbla_lda_invalid", .kind = .bad, .transa = 'N', .transb = 'N', .m = 1, .n = 1, .lda = -1, .ldc = 1, .info = 5 },
    .{ .name = "xerbla_ldc_invalid", .kind = .bad, .transa = 'N', .transb = 'N', .m = 1, .n = 1, .lda = 1, .ldc = -1, .info = 8 },
    .{ .name = "n_zero", .kind = .run, .transa = 'N', .transb = 'N', .m = 1, .n = 0, .lda = 1, .ldc = 1, .alpha = .{ 1.0, 0 }, .beta = .{ 1.0, 0 } },
    .{ .name = "m_zero", .kind = .run, .transa = 'N', .transb = 'N', .m = 0, .n = 1, .lda = 1, .ldc = 1, .alpha = .{ 1.0, 0 }, .beta = .{ 1.0, 0 } },
    .{ .name = "c_api_matrix_n_100_m_100", .kind = .run, .transa = 'N', .transb = 'N', .m = 100, .n = 100, .lda = 100, .ldc = 100, .alpha = .{ 2.0, 0 }, .beta = .{ 3.0, 0 } },
    .{ .name = "c_api_matrix_n_100_m_100_transA", .kind = .run, .transa = 'T', .transb = 'N', .m = 100, .n = 100, .lda = 100, .ldc = 100, .alpha = .{ 2.0, 0 }, .beta = .{ 3.0, 0 } },
    .{ .name = "c_api_matrix_n_100_m_100_alpha_zero", .kind = .run, .transa = 'N', .transb = 'N', .m = 100, .n = 100, .lda = 100, .ldc = 100, .alpha = .{ 0.0, 0 }, .beta = .{ 1.0, 0 } },
    .{ .name = "c_api_matrix_n_100_m_50", .kind = .run, .transa = 'N', .transb = 'N', .m = 50, .n = 100, .lda = 50, .ldc = 50, .alpha = .{ 3.0, 0 }, .beta = .{ 4.0, 0 } },
};

/// `test_dgeadd.c`, column-major.
const geadd_d = [_]Case{
    .{ .name = "matrix_n_100_m_100", .kind = .run, .transa = 'N', .transb = 'N', .m = 100, .n = 100, .lda = 100, .ldc = 100, .alpha = .{ 3.0, 0 }, .beta = .{ 3.0, 0 } },
    .{ .name = "matrix_n_100_m_100_alpha_zero", .kind = .run, .transa = 'N', .transb = 'N', .m = 100, .n = 100, .lda = 100, .ldc = 100, .alpha = .{ 0.0, 0 }, .beta = .{ 2.5, 0 } },
    .{ .name = "matrix_n_100_m_100_beta_zero", .kind = .run, .transa = 'N', .transb = 'N', .m = 100, .n = 100, .lda = 100, .ldc = 100, .alpha = .{ 3.0, 0 }, .beta = .{ 0.0, 0 } },
    .{ .name = "matrix_n_100_m_100_alpha_beta_zero", .kind = .run, .transa = 'N', .transb = 'N', .m = 100, .n = 100, .lda = 100, .ldc = 100, .alpha = .{ 0.0, 0 }, .beta = .{ 0.0, 0 } },
    .{ .name = "matrix_n_100_m_50", .kind = .run, .transa = 'N', .transb = 'N', .m = 50, .n = 100, .lda = 50, .ldc = 50, .alpha = .{ 1.0, 0 }, .beta = .{ 1.0, 0 } },
    .{ .name = "xerbla_n_invalid", .kind = .bad, .transa = 'N', .transb = 'N', .m = 1, .n = -1, .lda = 1, .ldc = 1, .info = 2 },
    .{ .name = "xerbla_m_invalid", .kind = .bad, .transa = 'N', .transb = 'N', .m = -1, .n = 1, .lda = 1, .ldc = 1, .info = 1 },
    .{ .name = "xerbla_lda_invalid", .kind = .bad, .transa = 'N', .transb = 'N', .m = 1, .n = 1, .lda = -1, .ldc = 1, .info = 5 },
    .{ .name = "xerbla_ldc_invalid", .kind = .bad, .transa = 'N', .transb = 'N', .m = 1, .n = 1, .lda = 1, .ldc = -1, .info = 8 },
    .{ .name = "n_zero", .kind = .run, .transa = 'N', .transb = 'N', .m = 1, .n = 0, .lda = 1, .ldc = 1, .alpha = .{ 1.0, 0 }, .beta = .{ 1.0, 0 } },
    .{ .name = "m_zero", .kind = .run, .transa = 'N', .transb = 'N', .m = 0, .n = 1, .lda = 1, .ldc = 1, .alpha = .{ 1.0, 0 }, .beta = .{ 1.0, 0 } },
    .{ .name = "c_api_matrix_n_100_m_100", .kind = .run, .transa = 'N', .transb = 'N', .m = 100, .n = 100, .lda = 100, .ldc = 100, .alpha = .{ 2.0, 0 }, .beta = .{ 3.0, 0 } },
    .{ .name = "c_api_matrix_n_100_m_100_alpha_zero", .kind = .run, .transa = 'N', .transb = 'N', .m = 100, .n = 100, .lda = 100, .ldc = 100, .alpha = .{ 0.0, 0 }, .beta = .{ 1.0, 0 } },
    .{ .name = "c_api_matrix_n_100_m_50", .kind = .run, .transa = 'N', .transb = 'N', .m = 50, .n = 100, .lda = 50, .ldc = 50, .alpha = .{ 3.0, 0 }, .beta = .{ 4.0, 0 } },
    .{ .name = "c_api_matrix_n_100_m_100_transA", .kind = .run, .transa = 'T', .transb = 'N', .m = 100, .n = 100, .lda = 100, .ldc = 100, .alpha = .{ 2.0, 0 }, .beta = .{ 3.0, 0 } },
};

/// `test_cgeadd.c`, column-major.
const geadd_c = [_]Case{
    .{ .name = "matrix_n_100_m_100", .kind = .run, .transa = 'N', .transb = 'N', .m = 100, .n = 100, .lda = 100, .ldc = 100, .alpha = .{ 3.0, 2.0 }, .beta = .{ 1.0, 3.0 } },
    .{ .name = "matrix_n_100_m_100_alpha_zero", .kind = .run, .transa = 'N', .transb = 'N', .m = 100, .n = 100, .lda = 100, .ldc = 100, .alpha = .{ 0.0, 0.0 }, .beta = .{ 2.5, 1.0 } },
    .{ .name = "matrix_n_100_m_100_beta_zero", .kind = .run, .transa = 'N', .transb = 'N', .m = 100, .n = 100, .lda = 100, .ldc = 100, .alpha = .{ 3.0, 1.5 }, .beta = .{ 0.0, 0.0 } },
    .{ .name = "matrix_n_100_m_100_alpha_beta_zero", .kind = .run, .transa = 'N', .transb = 'N', .m = 100, .n = 100, .lda = 100, .ldc = 100, .alpha = .{ 0.0, 0.0 }, .beta = .{ 0.0, 0.0 } },
    .{ .name = "matrix_n_100_m_50", .kind = .run, .transa = 'N', .transb = 'N', .m = 50, .n = 100, .lda = 50, .ldc = 50, .alpha = .{ 1.0, 1.0 }, .beta = .{ 1.0, 1.0 } },
    .{ .name = "xerbla_n_invalid", .kind = .bad, .transa = 'N', .transb = 'N', .m = 1, .n = -1, .lda = 1, .ldc = 1, .info = 2 },
    .{ .name = "xerbla_m_invalid", .kind = .bad, .transa = 'N', .transb = 'N', .m = -1, .n = 1, .lda = 1, .ldc = 1, .info = 1 },
    .{ .name = "xerbla_lda_invalid", .kind = .bad, .transa = 'N', .transb = 'N', .m = 1, .n = 1, .lda = -1, .ldc = 1, .info = 5 },
    .{ .name = "xerbla_ldc_invalid", .kind = .bad, .transa = 'N', .transb = 'N', .m = 1, .n = 1, .lda = 1, .ldc = -1, .info = 8 },
    .{ .name = "n_zero", .kind = .run, .transa = 'N', .transb = 'N', .m = 1, .n = 0, .lda = 1, .ldc = 1, .alpha = .{ 1.0, 1.0 }, .beta = .{ 1.0, 1.0 } },
    .{ .name = "m_zero", .kind = .run, .transa = 'N', .transb = 'N', .m = 0, .n = 1, .lda = 1, .ldc = 1, .alpha = .{ 1.0, 1.0 }, .beta = .{ 1.0, 1.0 } },
    .{ .name = "c_api_matrix_n_100_m_100", .kind = .run, .transa = 'N', .transb = 'N', .m = 100, .n = 100, .lda = 100, .ldc = 100, .alpha = .{ 2.0, 1.0 }, .beta = .{ 1.0, 3.0 } },
    .{ .name = "c_api_matrix_n_100_m_100_alpha_zero", .kind = .run, .transa = 'N', .transb = 'N', .m = 100, .n = 100, .lda = 100, .ldc = 100, .alpha = .{ 0.0, 0.0 }, .beta = .{ 1.0, 1.0 } },
    .{ .name = "c_api_matrix_n_100_m_50", .kind = .run, .transa = 'N', .transb = 'N', .m = 50, .n = 100, .lda = 50, .ldc = 50, .alpha = .{ 2.0, 3.0 }, .beta = .{ 2.0, 4.0 } },
    .{ .name = "c_api_matrix_n_100_m_100_transA", .kind = .run, .transa = 'T', .transb = 'N', .m = 100, .n = 100, .lda = 100, .ldc = 100, .alpha = .{ 2.0, -1.0 }, .beta = .{ 1.5, 3.0 } },
};

/// `test_zgeadd.c`, column-major.
const geadd_z = [_]Case{
    .{ .name = "matrix_n_100_m_100", .kind = .run, .transa = 'N', .transb = 'N', .m = 100, .n = 100, .lda = 100, .ldc = 100, .alpha = .{ 3.0, 2.0 }, .beta = .{ 1.0, 3.0 } },
    .{ .name = "matrix_n_100_m_100_alpha_zero", .kind = .run, .transa = 'N', .transb = 'N', .m = 100, .n = 100, .lda = 100, .ldc = 100, .alpha = .{ 0.0, 0.0 }, .beta = .{ 1.0, 1.0 } },
    .{ .name = "matrix_n_100_m_100_beta_zero", .kind = .run, .transa = 'N', .transb = 'N', .m = 100, .n = 100, .lda = 100, .ldc = 100, .alpha = .{ 3.0, 1.5 }, .beta = .{ 0.0, 0.0 } },
    .{ .name = "matrix_n_100_m_100_alpha_beta_zero", .kind = .run, .transa = 'N', .transb = 'N', .m = 100, .n = 100, .lda = 100, .ldc = 100, .alpha = .{ 0.0, 0.0 }, .beta = .{ 0.0, 0.0 } },
    .{ .name = "matrix_n_100_m_50", .kind = .run, .transa = 'N', .transb = 'N', .m = 50, .n = 100, .lda = 50, .ldc = 50, .alpha = .{ 1.0, 1.0 }, .beta = .{ 1.0, 1.0 } },
    .{ .name = "xerbla_n_invalid", .kind = .bad, .transa = 'N', .transb = 'N', .m = 1, .n = -1, .lda = 1, .ldc = 1, .info = 2 },
    .{ .name = "xerbla_m_invalid", .kind = .bad, .transa = 'N', .transb = 'N', .m = -1, .n = 1, .lda = 1, .ldc = 1, .info = 1 },
    .{ .name = "xerbla_lda_invalid", .kind = .bad, .transa = 'N', .transb = 'N', .m = 1, .n = 1, .lda = -1, .ldc = 1, .info = 5 },
    .{ .name = "xerbla_ldc_invalid", .kind = .bad, .transa = 'N', .transb = 'N', .m = 1, .n = 1, .lda = 1, .ldc = -1, .info = 8 },
    .{ .name = "n_zero", .kind = .run, .transa = 'N', .transb = 'N', .m = 1, .n = 0, .lda = 1, .ldc = 1, .alpha = .{ 1.0, 1.0 }, .beta = .{ 1.0, 1.0 } },
    .{ .name = "m_zero", .kind = .run, .transa = 'N', .transb = 'N', .m = 0, .n = 1, .lda = 1, .ldc = 1, .alpha = .{ 1.0, 1.0 }, .beta = .{ 1.0, 1.0 } },
    .{ .name = "c_api_matrix_n_100_m_100", .kind = .run, .transa = 'N', .transb = 'N', .m = 100, .n = 100, .lda = 100, .ldc = 100, .alpha = .{ 2.0, 1.0 }, .beta = .{ 1.0, 3.0 } },
    .{ .name = "c_api_matrix_n_100_m_50", .kind = .run, .transa = 'N', .transb = 'N', .m = 50, .n = 100, .lda = 50, .ldc = 50, .alpha = .{ 2.0, 3.0 }, .beta = .{ 2.0, 4.0 } },
    .{ .name = "c_api_matrix_n_100_m_100_transA", .kind = .run, .transa = 'T', .transb = 'N', .m = 100, .n = 100, .lda = 100, .ldc = 100, .alpha = .{ 2.0, -1.0 }, .beta = .{ 1.5, 3.0 } },
};
