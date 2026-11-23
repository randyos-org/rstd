//! `test_[cz]gemm.c` and `test_?gemmt.c`: complex `gemm` with `R`, and `gemmt`.

const std = @import("std");
const ob = @import("../__root__.zig");
const blas = @import("../../../__root__.zig");
const ext = @import("../../../ext/interface.zig");
const l3 = @import("../../../blas3/core.zig");
const check = @import("../../check.zig");
const apple_accel = @import("../../../apple_accel.zig");
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
const diffNorm = util.diffNorm;
const expectNorm = util.expectNorm;
const transOf = util.transOf;
const matDiff = util.matDiff;
const transposedLetter = util.transposedLetter;

/// `test_[cz]gemm.c`: `gemm` with a conjugated operand (`R`, `C`) against
/// `gemm` on a copy that `imatcopy` conjugated (`R` becoming `N`). Every
/// case is 100 square; on a rig, 23.
pub fn complexGemm() !void {
    // both calls on rstd: Accelerate has no `R`, and two backends' f32 sums
    // differ past the source's 1e-4; after init, which routes
    _ = blas.vtable();
    const was = apple_accel.routing();
    _ = apple_accel.route(false);
    defer _ = apple_accel.route(was);
    var prng = std.Random.DefaultPrng.init(0xC2);
    const al = testing.allocator;
    inline for (.{ .{ C32, &gemm_c }, .{ C64, &gemm_z } }) |p| for (p[1]) |case| {
        const T = p[0];
        var c = case;
        if (check.rig()) {
            inline for (.{ "m", "n", "k", "lda", "ldb", "ldc" }) |f| @field(c, f) = 23;
        }
        const r = prng.random();
        const ta = transposedLetter(c.transa);
        const tb = transposedLetter(c.transb);
        const rows_a: isize = if (ta) c.k else c.m;
        const cols_a: isize = if (ta) c.m else c.k;
        const rows_b: isize = if (tb) c.n else c.k;
        const cols_b: isize = if (tb) c.k else c.n;
        const a = try al.alloc(T, @intCast(c.lda * cols_a));
        defer al.free(a);
        const b = try al.alloc(T, @intCast(c.ldb * cols_b));
        defer al.free(b);
        const cm = try al.alloc(T, @intCast(c.ldc * c.n));
        defer al.free(cm);
        fill(T, a, r);
        fill(T, b, r);
        fill(T, cm, r);
        const av = try al.dupe(T, a);
        defer al.free(av);
        const bv = try al.dupe(T, b);
        defer al.free(bv);
        const cv = try al.dupe(T, cm);
        defer al.free(cv);
        const one = val(T, .{ 1, 0 });
        var tav = transOf(c.transa);
        var tbv = transOf(c.transb);
        if (c.transa == 'R') {
            try ext.imatcopy(al, ob.io, T, .r, rows_a, cols_a, one, av.ptr, c.lda, c.lda);
            tav = .n;
        }
        if (c.transb == 'R') {
            try ext.imatcopy(al, ob.io, T, .r, rows_b, cols_b, one, bv.ptr, c.ldb, c.ldb);
            tbv = .n;
        }
        const alpha = val(T, c.alpha);
        const beta = val(T, c.beta);
        try l3.gemm(al, ob.io, T, tav, tbv, c.m, c.n, c.k, alpha, av.ptr, c.lda, bv.ptr, c.ldb, beta, cv.ptr, c.ldc);
        try l3.gemm(al, ob.io, T, transOf(c.transa), transOf(c.transb), c.m, c.n, c.k, alpha, a.ptr, c.lda, b.ptr, c.ldb, beta, cm.ptr, c.ldc);
        try expectNorm(c.name, T, matDiff(T, @intCast(c.m), @intCast(c.n), cm, cv, @intCast(c.ldc)), tol(T, single_eps, double_eps));
    };
}

test "utest_ext: complex gemm with R against gemm on the conjugated copy" {
    try complexGemm();
}

/// `gemmt` against `gemm`'s triangle for one case of `test_?gemmt.c`, or its
/// bad argument; the norm of differences over the whole of `C`, by its
/// floats, as the sources take it.
fn gemmtCase(comptime T: type, c: Case, r: std.Random) !void {
    const al = testing.allocator;
    const uplo = ob.letter(blas.Uplo, c.uplo);
    if (c.kind == .bad) {
        // `O` is no `Uplo`
        if (c.info == 1) return ob.expect(uplo == null, c.name);
        var x: [1]T = @splat(ob.zero(T));
        const want: anyerror = switch (c.info) {
            2, 3 => error.InvalidTrans,
            4 => error.InvalidM,
            5 => error.InvalidK,
            8 => error.InvalidLda,
            10 => error.InvalidLdb,
            else => error.InvalidLdc,
        };
        const one = val(T, .{ 1, 1 });
        return ob.expectError(want, l3.gemmt(al, ob.io, T, uplo.?, transOf(c.transa), transOf(c.transb), c.m, c.k, one, &x, c.lda, &x, c.ldb, one, &x, c.ldc), c.name);
    }
    const m: usize = @intCast(c.m);
    const ldc: usize = @intCast(c.ldc);
    const cols_a: isize = if (transposedLetter(c.transa)) c.m else c.k;
    const cols_b: isize = if (transposedLetter(c.transb)) c.k else c.m;
    const a = try al.alloc(T, @intCast(c.lda * cols_a));
    defer al.free(a);
    const b = try al.alloc(T, @intCast(c.ldb * cols_b));
    defer al.free(b);
    const cm = try al.alloc(T, m * ldc);
    defer al.free(cm);
    fill(T, a, r);
    fill(T, b, r);
    fill(T, cm, r);
    const full = try al.dupe(T, cm);
    defer al.free(full);
    const want = try al.dupe(T, cm);
    defer al.free(want);
    const alpha = val(T, c.alpha);
    const beta = val(T, c.beta);
    try l3.gemm(al, ob.io, T, transOf(c.transa), transOf(c.transb), c.m, c.m, c.k, alpha, a.ptr, c.lda, b.ptr, c.ldb, beta, full.ptr, c.ldc);
    for (0..m) |j| for (0..m) |i| {
        if (if (c.uplo == 'U') i <= j else i >= j) want[i + j * ldc] = full[i + j * ldc];
    };
    try l3.gemmt(al, ob.io, T, uplo.?, transOf(c.transa), transOf(c.transb), c.m, c.k, alpha, a.ptr, c.lda, b.ptr, c.ldb, beta, cm.ptr, c.ldc);
    const floats: f64 = @floatFromInt(m * ldc * @as(usize, if (comptime isComplex(T)) 2 else 1));
    try expectNorm(c.name, T, diffNorm(T, @intCast(m * ldc), cm, want, 1) / floats, tol(T, single_eps, double_eps));
}

pub fn gemmt() !void {
    var prng = std.Random.DefaultPrng.init(0xC3);
    inline for (.{ .{ f32, &gemmt_s }, .{ f64, &gemmt_d }, .{ C32, &gemmt_c }, .{ C64, &gemmt_z } }) |p| {
        for (p[1]) |c| try gemmtCase(p[0], c, prng.random());
    }
}

test "utest_ext: gemmt against gemm's triangle, and its bad arguments" {
    try gemmt();
}

/// `test_cgemm.c`, column-major.
const gemm_c = [_]Case{
    .{ .name = "conjtransa_conjnotransb", .kind = .run, .transa = 'C', .transb = 'R', .m = 100, .n = 100, .k = 100, .lda = 100, .ldb = 100, .ldc = 100, .alpha = .{ -2.0, 1.0 }, .beta = .{ 1.0, -1.0 } },
    .{ .name = "notransa_conjnotransb", .kind = .run, .transa = 'N', .transb = 'R', .m = 100, .n = 100, .k = 100, .lda = 100, .ldb = 100, .ldc = 100, .alpha = .{ -2.0, 1.0 }, .beta = .{ 1.0, -1.0 } },
    .{ .name = "conjnotransa_conjtransb", .kind = .run, .transa = 'R', .transb = 'C', .m = 100, .n = 100, .k = 100, .lda = 100, .ldb = 100, .ldc = 100, .alpha = .{ -2.0, 1.0 }, .beta = .{ 1.0, -1.0 } },
    .{ .name = "conjnotransa_notransb", .kind = .run, .transa = 'R', .transb = 'N', .m = 100, .n = 100, .k = 100, .lda = 100, .ldb = 100, .ldc = 100, .alpha = .{ -2.0, 1.0 }, .beta = .{ 1.0, -1.0 } },
    .{ .name = "conjnotransa_conjnotransb", .kind = .run, .transa = 'R', .transb = 'R', .m = 100, .n = 100, .k = 100, .lda = 100, .ldb = 100, .ldc = 100, .alpha = .{ -2.0, 1.0 }, .beta = .{ 1.0, -1.0 } },
    .{ .name = "conjnotransa_transb", .kind = .run, .transa = 'R', .transb = 'T', .m = 100, .n = 100, .k = 100, .lda = 100, .ldb = 100, .ldc = 100, .alpha = .{ -2.0, 1.0 }, .beta = .{ 1.0, -1.0 } },
    .{ .name = "transa_conjnotransb", .kind = .run, .transa = 'T', .transb = 'R', .m = 100, .n = 100, .k = 100, .lda = 100, .ldb = 100, .ldc = 100, .alpha = .{ -2.0, 1.0 }, .beta = .{ 1.0, -1.0 } },
};

/// `test_zgemm.c`, column-major.
const gemm_z = [_]Case{
    .{ .name = "conjtransa_conjnotransb", .kind = .run, .transa = 'C', .transb = 'R', .m = 100, .n = 100, .k = 100, .lda = 100, .ldb = 100, .ldc = 100, .alpha = .{ -2.0, 1.0 }, .beta = .{ 1.0, -1.0 } },
    .{ .name = "notransa_conjnotransb", .kind = .run, .transa = 'N', .transb = 'R', .m = 100, .n = 100, .k = 100, .lda = 100, .ldb = 100, .ldc = 100, .alpha = .{ -2.0, 1.0 }, .beta = .{ 1.0, -1.0 } },
    .{ .name = "conjnotransa_conjtransb", .kind = .run, .transa = 'R', .transb = 'C', .m = 100, .n = 100, .k = 100, .lda = 100, .ldb = 100, .ldc = 100, .alpha = .{ -2.0, 1.0 }, .beta = .{ 1.0, -1.0 } },
    .{ .name = "conjnotransa_notransb", .kind = .run, .transa = 'R', .transb = 'N', .m = 100, .n = 100, .k = 100, .lda = 100, .ldb = 100, .ldc = 100, .alpha = .{ -2.0, 1.0 }, .beta = .{ 1.0, -1.0 } },
    .{ .name = "conjnotransa_conjnotransb", .kind = .run, .transa = 'R', .transb = 'R', .m = 100, .n = 100, .k = 100, .lda = 100, .ldb = 100, .ldc = 100, .alpha = .{ -2.0, 1.0 }, .beta = .{ 1.0, -1.0 } },
    .{ .name = "conjnotransa_transb", .kind = .run, .transa = 'R', .transb = 'T', .m = 100, .n = 100, .k = 100, .lda = 100, .ldb = 100, .ldc = 100, .alpha = .{ -2.0, 1.0 }, .beta = .{ 1.0, -1.0 } },
    .{ .name = "transa_conjnotransb", .kind = .run, .transa = 'T', .transb = 'R', .m = 100, .n = 100, .k = 100, .lda = 100, .ldb = 100, .ldc = 100, .alpha = .{ -2.0, 1.0 }, .beta = .{ 1.0, -1.0 } },
};

/// `test_sgemmt.c`, column-major.
const gemmt_s = [_]Case{
    .{ .name = "upper_M_50_K_50_a_notrans_b_notrans", .kind = .run, .uplo = 'U', .transa = 'N', .transb = 'N', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 1.5, 0 }, .beta = .{ 2.0, 0 } },
    .{ .name = "upper_M_100_K_50_a_trans_b_notrans", .kind = .run, .uplo = 'U', .transa = 'T', .transb = 'N', .m = 100, .k = 50, .lda = 50, .ldb = 50, .ldc = 100, .alpha = .{ 1.0, 0 }, .beta = .{ 2.0, 0 } },
    .{ .name = "upper_M_50_K_100_a_notrans_b_trans", .kind = .run, .uplo = 'U', .transa = 'N', .transb = 'T', .m = 50, .k = 100, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 1.0, 0 }, .beta = .{ 0.0, 0 } },
    .{ .name = "upper_M_50_K_50_a_trans_b_trans", .kind = .run, .uplo = 'U', .transa = 'T', .transb = 'T', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 1.5, 0 }, .beta = .{ 2.0, 0 } },
    .{ .name = "upper_alpha_zero", .kind = .run, .uplo = 'U', .transa = 'N', .transb = 'N', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 0.0, 0 }, .beta = .{ 2.0, 0 } },
    .{ .name = "upper_beta_one", .kind = .run, .uplo = 'U', .transa = 'N', .transb = 'N', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 2.0, 0 }, .beta = .{ 1.0, 0 } },
    .{ .name = "lower_M_50_K_50_a_notrans_b_notrans", .kind = .run, .uplo = 'L', .transa = 'N', .transb = 'N', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 1.5, 0 }, .beta = .{ 2.0, 0 } },
    .{ .name = "lower_M_100_K_50_a_trans_b_notrans", .kind = .run, .uplo = 'L', .transa = 'T', .transb = 'N', .m = 100, .k = 50, .lda = 50, .ldb = 50, .ldc = 100, .alpha = .{ 1.0, 0 }, .beta = .{ 2.0, 0 } },
    .{ .name = "lower_M_50_K_100_a_notrans_b_trans", .kind = .run, .uplo = 'L', .transa = 'N', .transb = 'T', .m = 50, .k = 100, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 1.0, 0 }, .beta = .{ 0.0, 0 } },
    .{ .name = "lower_M_50_K_50_a_trans_b_trans", .kind = .run, .uplo = 'L', .transa = 'T', .transb = 'T', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 1.5, 0 }, .beta = .{ 2.0, 0 } },
    .{ .name = "lower_alpha_zero", .kind = .run, .uplo = 'L', .transa = 'N', .transb = 'N', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 0.0, 0 }, .beta = .{ 2.0, 0 } },
    .{ .name = "lower_beta_one", .kind = .run, .uplo = 'L', .transa = 'N', .transb = 'N', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 2.0, 0 }, .beta = .{ 1.0, 0 } },
    .{ .name = "xerbla_uplo_invalid", .kind = .bad, .uplo = 'O', .transa = 'N', .transb = 'N', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .info = 1 },
    .{ .name = "xerbla_transa_invalid", .kind = .bad, .uplo = 'U', .transa = 'O', .transb = 'N', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .info = 2 },
    .{ .name = "xerbla_transb_invalid", .kind = .bad, .uplo = 'U', .transa = 'N', .transb = 'O', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .info = 3 },
    .{ .name = "xerbla_m_invalid", .kind = .bad, .uplo = 'U', .transa = 'N', .transb = 'N', .m = -1, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .info = 4 },
    .{ .name = "xerbla_k_invalid", .kind = .bad, .uplo = 'U', .transa = 'N', .transb = 'N', .m = 50, .k = -1, .lda = 50, .ldb = 50, .ldc = 50, .info = 5 },
    .{ .name = "xerbla_lda_invalid", .kind = .bad, .uplo = 'U', .transa = 'T', .transb = 'N', .m = 50, .k = 100, .lda = 50, .ldb = 100, .ldc = 50, .info = 8 },
    .{ .name = "xerbla_ldb_invalid", .kind = .bad, .uplo = 'U', .transa = 'N', .transb = 'N', .m = 50, .k = 100, .lda = 50, .ldb = 50, .ldc = 50, .info = 10 },
    .{ .name = "xerbla_ldc_invalid", .kind = .bad, .uplo = 'U', .transa = 'T', .transb = 'N', .m = 100, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .info = 13 },
};

/// `test_dgemmt.c`, column-major.
const gemmt_d = [_]Case{
    .{ .name = "upper_M_50_K_50_a_notrans_b_notrans", .kind = .run, .uplo = 'U', .transa = 'N', .transb = 'N', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 1.5, 0 }, .beta = .{ 2.0, 0 } },
    .{ .name = "upper_M_100_K_50_a_trans_b_notrans", .kind = .run, .uplo = 'U', .transa = 'T', .transb = 'N', .m = 100, .k = 50, .lda = 50, .ldb = 50, .ldc = 100, .alpha = .{ 1.0, 0 }, .beta = .{ 2.0, 0 } },
    .{ .name = "upper_M_50_K_100_a_notrans_b_trans", .kind = .run, .uplo = 'U', .transa = 'N', .transb = 'T', .m = 50, .k = 100, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 1.0, 0 }, .beta = .{ 0.0, 0 } },
    .{ .name = "upper_M_50_K_50_a_trans_b_trans", .kind = .run, .uplo = 'U', .transa = 'T', .transb = 'T', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 1.5, 0 }, .beta = .{ 2.0, 0 } },
    .{ .name = "upper_alpha_zero", .kind = .run, .uplo = 'U', .transa = 'N', .transb = 'N', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 0.0, 0 }, .beta = .{ 2.0, 0 } },
    .{ .name = "upper_beta_one", .kind = .run, .uplo = 'U', .transa = 'N', .transb = 'N', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 2.0, 0 }, .beta = .{ 1.0, 0 } },
    .{ .name = "lower_M_50_K_50_a_notrans_b_notrans", .kind = .run, .uplo = 'L', .transa = 'N', .transb = 'N', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 1.5, 0 }, .beta = .{ 2.0, 0 } },
    .{ .name = "lower_M_100_K_50_a_trans_b_notrans", .kind = .run, .uplo = 'L', .transa = 'T', .transb = 'N', .m = 100, .k = 50, .lda = 50, .ldb = 50, .ldc = 100, .alpha = .{ 1.0, 0 }, .beta = .{ 2.0, 0 } },
    .{ .name = "lower_M_50_K_100_a_notrans_b_trans", .kind = .run, .uplo = 'L', .transa = 'N', .transb = 'T', .m = 50, .k = 100, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 1.0, 0 }, .beta = .{ 0.0, 0 } },
    .{ .name = "lower_M_50_K_50_a_trans_b_trans", .kind = .run, .uplo = 'L', .transa = 'T', .transb = 'T', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 1.5, 0 }, .beta = .{ 2.0, 0 } },
    .{ .name = "lower_alpha_zero", .kind = .run, .uplo = 'L', .transa = 'N', .transb = 'N', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 0.0, 0 }, .beta = .{ 2.0, 0 } },
    .{ .name = "lower_beta_one", .kind = .run, .uplo = 'L', .transa = 'N', .transb = 'N', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 2.0, 0 }, .beta = .{ 1.0, 0 } },
    .{ .name = "xerbla_uplo_invalid", .kind = .bad, .uplo = 'O', .transa = 'N', .transb = 'N', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .info = 1 },
    .{ .name = "xerbla_transa_invalid", .kind = .bad, .uplo = 'U', .transa = 'O', .transb = 'N', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .info = 2 },
    .{ .name = "xerbla_transb_invalid", .kind = .bad, .uplo = 'U', .transa = 'N', .transb = 'O', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .info = 3 },
    .{ .name = "xerbla_m_invalid", .kind = .bad, .uplo = 'U', .transa = 'N', .transb = 'N', .m = -1, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .info = 4 },
    .{ .name = "xerbla_k_invalid", .kind = .bad, .uplo = 'U', .transa = 'N', .transb = 'N', .m = 50, .k = -1, .lda = 50, .ldb = 50, .ldc = 50, .info = 5 },
    .{ .name = "xerbla_lda_invalid", .kind = .bad, .uplo = 'U', .transa = 'T', .transb = 'N', .m = 50, .k = 100, .lda = 50, .ldb = 100, .ldc = 50, .info = 8 },
    .{ .name = "xerbla_ldb_invalid", .kind = .bad, .uplo = 'U', .transa = 'N', .transb = 'N', .m = 50, .k = 100, .lda = 50, .ldb = 50, .ldc = 50, .info = 10 },
    .{ .name = "xerbla_ldc_invalid", .kind = .bad, .uplo = 'U', .transa = 'T', .transb = 'N', .m = 100, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .info = 13 },
};

/// `test_cgemmt.c`, column-major.
const gemmt_c = [_]Case{
    .{ .name = "upper_M_50_K_50_a_notrans_b_notrans", .kind = .run, .uplo = 'U', .transa = 'N', .transb = 'N', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ -1.0, 1.0 }, .beta = .{ 2.0, 1.0 } },
    .{ .name = "upper_M_50_K_25_a_trans_b_notrans", .kind = .run, .uplo = 'U', .transa = 'T', .transb = 'N', .m = 50, .k = 25, .lda = 25, .ldb = 25, .ldc = 50, .alpha = .{ 1.0, 1.0 }, .beta = .{ 2.0, 1.0 } },
    .{ .name = "upper_M_25_K_50_a_notrans_b_trans", .kind = .run, .uplo = 'U', .transa = 'N', .transb = 'T', .m = 25, .k = 50, .lda = 25, .ldb = 25, .ldc = 25, .alpha = .{ 1.0, 1.0 }, .beta = .{ 0.0, 0.0 } },
    .{ .name = "upper_M_50_K_50_a_trans_b_trans", .kind = .run, .uplo = 'U', .transa = 'T', .transb = 'T', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ -1.0, 1.0 }, .beta = .{ 2.0, 1.0 } },
    .{ .name = "upper_M_25_K_50_a_conjnotrans_b_conjnotrans", .kind = .run, .uplo = 'U', .transa = 'R', .transb = 'R', .m = 25, .k = 50, .lda = 25, .ldb = 50, .ldc = 25, .alpha = .{ 2.0, 1.0 }, .beta = .{ -1.0, 2.0 } },
    .{ .name = "upper_M_50_K_50_a_conjtrans_b_conjnotrans", .kind = .run, .uplo = 'U', .transa = 'C', .transb = 'R', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 1.0, 1.0 }, .beta = .{ 1.0, 1.0 } },
    .{ .name = "upper_M_50_K_50_a_conjnotrans_b_conjtrans", .kind = .run, .uplo = 'U', .transa = 'R', .transb = 'C', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 2.0, 1.0 }, .beta = .{ -1.0, 2.0 } },
    .{ .name = "upper_M_50_K_25_a_conjtrans_b_conjtrans", .kind = .run, .uplo = 'U', .transa = 'C', .transb = 'C', .m = 50, .k = 25, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 2.0, 1.0 }, .beta = .{ -1.0, 2.0 } },
    .{ .name = "upper_alpha_zero", .kind = .run, .uplo = 'U', .transa = 'N', .transb = 'N', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 0.0, 0.0 }, .beta = .{ 2.0, 1.0 } },
    .{ .name = "upper_beta_one", .kind = .run, .uplo = 'U', .transa = 'N', .transb = 'N', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 1.0, 1.0 }, .beta = .{ 1.0, 0.0 } },
    .{ .name = "lower_M_50_K_50_a_notrans_b_notrans", .kind = .run, .uplo = 'L', .transa = 'N', .transb = 'N', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ -1.0, 1.0 }, .beta = .{ 2.0, 1.0 } },
    .{ .name = "lower_M_50_K_25_a_trans_b_notrans", .kind = .run, .uplo = 'L', .transa = 'T', .transb = 'N', .m = 50, .k = 25, .lda = 25, .ldb = 25, .ldc = 50, .alpha = .{ 1.0, 1.0 }, .beta = .{ 2.0, 1.0 } },
    .{ .name = "lower_M_25_K_50_a_notrans_b_trans", .kind = .run, .uplo = 'L', .transa = 'N', .transb = 'T', .m = 25, .k = 50, .lda = 25, .ldb = 25, .ldc = 25, .alpha = .{ 1.0, 1.0 }, .beta = .{ 0.0, 0.0 } },
    .{ .name = "lower_M_50_K_50_a_trans_b_trans", .kind = .run, .uplo = 'L', .transa = 'T', .transb = 'T', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ -1.0, 1.0 }, .beta = .{ 2.0, 1.0 } },
    .{ .name = "lower_M_25_K_50_a_conjnotrans_b_conjnotrans", .kind = .run, .uplo = 'L', .transa = 'R', .transb = 'R', .m = 25, .k = 50, .lda = 25, .ldb = 50, .ldc = 25, .alpha = .{ 2.0, 1.0 }, .beta = .{ -1.0, 2.0 } },
    .{ .name = "lower_M_50_K_50_a_conjtrans_b_conjnotrans", .kind = .run, .uplo = 'L', .transa = 'C', .transb = 'R', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 1.0, 1.0 }, .beta = .{ 1.0, 1.0 } },
    .{ .name = "lower_M_50_K_50_a_conjnotrans_b_conjtrans", .kind = .run, .uplo = 'L', .transa = 'R', .transb = 'C', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 2.0, 1.0 }, .beta = .{ -1.0, 2.0 } },
    .{ .name = "lower_M_50_K_25_a_conjtrans_b_conjtrans", .kind = .run, .uplo = 'L', .transa = 'C', .transb = 'C', .m = 50, .k = 25, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 2.0, 1.0 }, .beta = .{ -1.0, 2.0 } },
    .{ .name = "lower_alpha_zero", .kind = .run, .uplo = 'L', .transa = 'N', .transb = 'N', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 0.0, 0.0 }, .beta = .{ 2.0, 1.0 } },
    .{ .name = "lower_beta_one", .kind = .run, .uplo = 'L', .transa = 'N', .transb = 'N', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 1.0, 1.0 }, .beta = .{ 1.0, 0.0 } },
    .{ .name = "c_api_colmajor_lower_M_50_K_50_a_conjtrans_b_conjnotrans", .kind = .run, .uplo = 'L', .transa = 'C', .transb = 'R', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 2.0, 1.0 }, .beta = .{ -1.0, 2.0 } },
    .{ .name = "xerbla_uplo_invalid", .kind = .bad, .uplo = 'O', .transa = 'N', .transb = 'N', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .info = 1 },
    .{ .name = "xerbla_transa_invalid", .kind = .bad, .uplo = 'U', .transa = 'O', .transb = 'N', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .info = 2 },
    .{ .name = "xerbla_transb_invalid", .kind = .bad, .uplo = 'U', .transa = 'N', .transb = 'O', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .info = 3 },
    .{ .name = "xerbla_m_invalid", .kind = .bad, .uplo = 'U', .transa = 'N', .transb = 'N', .m = -1, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .info = 4 },
    .{ .name = "xerbla_k_invalid", .kind = .bad, .uplo = 'U', .transa = 'N', .transb = 'N', .m = 50, .k = -1, .lda = 50, .ldb = 50, .ldc = 50, .info = 5 },
    .{ .name = "xerbla_lda_invalid", .kind = .bad, .uplo = 'U', .transa = 'T', .transb = 'N', .m = 50, .k = 100, .lda = 50, .ldb = 100, .ldc = 50, .info = 8 },
    .{ .name = "xerbla_ldb_invalid", .kind = .bad, .uplo = 'U', .transa = 'N', .transb = 'N', .m = 50, .k = 100, .lda = 50, .ldb = 50, .ldc = 50, .info = 10 },
    .{ .name = "xerbla_ldc_invalid", .kind = .bad, .uplo = 'U', .transa = 'T', .transb = 'N', .m = 100, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .info = 13 },
};

/// `test_zgemmt.c`, column-major.
const gemmt_z = [_]Case{
    .{ .name = "upper_M_50_K_50_a_notrans_b_notrans", .kind = .run, .uplo = 'U', .transa = 'N', .transb = 'N', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 1.5, 0.5 }, .beta = .{ 2.0, 1.0 } },
    .{ .name = "upper_M_50_K_25_a_trans_b_notrans", .kind = .run, .uplo = 'U', .transa = 'T', .transb = 'N', .m = 50, .k = 25, .lda = 25, .ldb = 25, .ldc = 50, .alpha = .{ 1.0, 1.0 }, .beta = .{ 2.0, 1.0 } },
    .{ .name = "upper_M_25_K_50_a_notrans_b_trans", .kind = .run, .uplo = 'U', .transa = 'N', .transb = 'T', .m = 25, .k = 50, .lda = 25, .ldb = 25, .ldc = 25, .alpha = .{ 1.0, 1.0 }, .beta = .{ 0.0, 0.0 } },
    .{ .name = "upper_M_50_K_50_a_trans_b_trans", .kind = .run, .uplo = 'U', .transa = 'T', .transb = 'T', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 1.5, 0.5 }, .beta = .{ 2.0, 1.0 } },
    .{ .name = "upper_M_25_K_50_a_conjnotrans_b_conjnotrans", .kind = .run, .uplo = 'U', .transa = 'R', .transb = 'R', .m = 25, .k = 50, .lda = 25, .ldb = 50, .ldc = 25, .alpha = .{ 2.0, 1.0 }, .beta = .{ 1.5, 2.0 } },
    .{ .name = "upper_M_50_K_50_a_conjtrans_b_conjnotrans", .kind = .run, .uplo = 'U', .transa = 'C', .transb = 'R', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 1.0, 1.0 }, .beta = .{ 1.0, 1.0 } },
    .{ .name = "upper_M_50_K_50_a_conjnotrans_b_conjtrans", .kind = .run, .uplo = 'U', .transa = 'R', .transb = 'C', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 2.0, 1.0 }, .beta = .{ 1.5, 2.0 } },
    .{ .name = "upper_M_50_K_25_a_conjtrans_b_conjtrans", .kind = .run, .uplo = 'U', .transa = 'C', .transb = 'C', .m = 50, .k = 25, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 2.0, 1.0 }, .beta = .{ 1.5, 2.0 } },
    .{ .name = "upper_alpha_zero", .kind = .run, .uplo = 'U', .transa = 'N', .transb = 'N', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 0.0, 0.0 }, .beta = .{ 2.0, 1.0 } },
    .{ .name = "upper_beta_one", .kind = .run, .uplo = 'U', .transa = 'N', .transb = 'N', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 1.0, 1.0 }, .beta = .{ 1.0, 0.0 } },
    .{ .name = "lower_M_50_K_50_a_notrans_b_notrans", .kind = .run, .uplo = 'L', .transa = 'N', .transb = 'N', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 1.5, 0.5 }, .beta = .{ 2.0, 1.0 } },
    .{ .name = "lower_M_50_K_25_a_trans_b_notrans", .kind = .run, .uplo = 'L', .transa = 'T', .transb = 'N', .m = 50, .k = 25, .lda = 25, .ldb = 25, .ldc = 50, .alpha = .{ 1.0, 1.0 }, .beta = .{ 2.0, 1.0 } },
    .{ .name = "lower_M_25_K_50_a_notrans_b_trans", .kind = .run, .uplo = 'L', .transa = 'N', .transb = 'T', .m = 25, .k = 50, .lda = 25, .ldb = 25, .ldc = 25, .alpha = .{ 1.0, 1.0 }, .beta = .{ 0.0, 0.0 } },
    .{ .name = "lower_M_50_K_50_a_trans_b_trans", .kind = .run, .uplo = 'L', .transa = 'T', .transb = 'T', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 1.5, 0.5 }, .beta = .{ 2.0, 1.0 } },
    .{ .name = "lower_M_25_K_50_a_conjnotrans_b_conjnotrans", .kind = .run, .uplo = 'L', .transa = 'R', .transb = 'R', .m = 25, .k = 50, .lda = 25, .ldb = 50, .ldc = 25, .alpha = .{ 2.0, 1.0 }, .beta = .{ 1.5, 2.0 } },
    .{ .name = "lower_M_50_K_50_a_conjtrans_b_conjnotrans", .kind = .run, .uplo = 'L', .transa = 'C', .transb = 'R', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 1.0, 1.0 }, .beta = .{ 1.0, 1.0 } },
    .{ .name = "lower_M_50_K_50_a_conjnotrans_b_conjtrans", .kind = .run, .uplo = 'L', .transa = 'R', .transb = 'C', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 2.0, 1.0 }, .beta = .{ 1.5, 2.0 } },
    .{ .name = "lower_M_50_K_25_a_conjtrans_b_conjtrans", .kind = .run, .uplo = 'L', .transa = 'C', .transb = 'C', .m = 50, .k = 25, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 2.0, 1.0 }, .beta = .{ 1.5, 2.0 } },
    .{ .name = "lower_alpha_zero", .kind = .run, .uplo = 'L', .transa = 'N', .transb = 'N', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 0.0, 0.0 }, .beta = .{ 2.0, 1.0 } },
    .{ .name = "lower_beta_one", .kind = .run, .uplo = 'L', .transa = 'N', .transb = 'N', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 1.0, 1.0 }, .beta = .{ 1.0, 0.0 } },
    .{ .name = "c_api_colmajor_lower_M_50_K_50_a_conjtrans_b_conjnotrans", .kind = .run, .uplo = 'L', .transa = 'C', .transb = 'R', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .alpha = .{ 2.0, 1.0 }, .beta = .{ 1.5, 2.0 } },
    .{ .name = "xerbla_uplo_invalid", .kind = .bad, .uplo = 'O', .transa = 'N', .transb = 'N', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .info = 1 },
    .{ .name = "xerbla_transa_invalid", .kind = .bad, .uplo = 'U', .transa = 'O', .transb = 'N', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .info = 2 },
    .{ .name = "xerbla_transb_invalid", .kind = .bad, .uplo = 'U', .transa = 'N', .transb = 'O', .m = 50, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .info = 3 },
    .{ .name = "xerbla_m_invalid", .kind = .bad, .uplo = 'U', .transa = 'N', .transb = 'N', .m = -1, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .info = 4 },
    .{ .name = "xerbla_k_invalid", .kind = .bad, .uplo = 'U', .transa = 'N', .transb = 'N', .m = 50, .k = -1, .lda = 50, .ldb = 50, .ldc = 50, .info = 5 },
    .{ .name = "xerbla_lda_invalid", .kind = .bad, .uplo = 'U', .transa = 'T', .transb = 'N', .m = 50, .k = 100, .lda = 50, .ldb = 100, .ldc = 50, .info = 8 },
    .{ .name = "xerbla_ldb_invalid", .kind = .bad, .uplo = 'U', .transa = 'N', .transb = 'N', .m = 50, .k = 100, .lda = 50, .ldb = 50, .ldc = 50, .info = 10 },
    .{ .name = "xerbla_ldc_invalid", .kind = .bad, .uplo = 'U', .transa = 'T', .transb = 'N', .m = 100, .k = 50, .lda = 50, .ldb = 50, .ldc = 50, .info = 13 },
};
