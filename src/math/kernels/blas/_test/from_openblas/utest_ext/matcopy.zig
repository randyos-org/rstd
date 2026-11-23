//! `test_?omatcopy.c` and `test_?imatcopy.c`: the copies against the definition.

const std = @import("std");
const ob = @import("../__root__.zig");
const ext = @import("../../../ext/interface.zig");
const util = @import("util.zig");

const testing = std.testing;
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
const matDiff = util.matDiff;
const transposedLetter = util.transposedLetter;

/// `omatcopy` (or, `in_place`, `imatcopy`) against the definition for one
/// case of `test_?omatcopy.c` or `test_?imatcopy.c`, or its bad argument.
/// The sources' `m` and `n` are the rows and columns of `A`.
fn matcopyCase(comptime T: type, comptime in_place: bool, c: Case, r: std.Random) !void {
    const al = testing.allocator;
    const tr = transOf(c.trans);
    if (c.kind == .bad) {
        var a: [1]T = @splat(ob.zero(T));
        var b: [1]T = @splat(ob.zero(T));
        const want: anyerror = switch (c.info) {
            2 => error.InvalidTrans,
            3 => error.InvalidRows,
            4 => error.InvalidCols,
            7 => error.InvalidLda,
            else => error.InvalidLdb,
        };
        const one = val(T, .{ 1, 1 });
        const res = if (in_place) ext.imatcopy(al, ob.io, T, tr, c.m, c.n, one, &a, c.lda, c.ldb) else ext.omatcopy(ob.io, T, tr, c.m, c.n, one, &a, c.lda, &b, c.ldb);
        return ob.expectError(want, res, c.name);
    }
    const rows: usize = @intCast(c.m);
    const cols: usize = @intCast(c.n);
    const lda: usize = @intCast(c.lda);
    const ldb: usize = @intCast(c.ldb);
    const t = transposedLetter(c.trans);
    const out_rows = if (t) cols else rows;
    const out_cols = if (t) rows else cols;
    const a = try al.alloc(T, @max(lda * cols, ldb * out_cols));
    defer al.free(a);
    const b = try al.alloc(T, ldb * out_cols);
    defer al.free(b);
    fill(T, a, r);
    @memset(b, ob.zero(T));
    const want = try al.dupe(T, b);
    defer al.free(want);
    const alpha = val(T, c.alpha);
    for (0..out_cols) |j| for (0..out_rows) |i| {
        const x = if (t) a[j + i * lda] else a[i + j * lda];
        want[i + j * ldb] = ob.mul(alpha, if (c.trans == 'R' or c.trans == 'C') ob.conj(x) else x);
    };
    if (in_place) {
        try ext.imatcopy(al, ob.io, T, tr, c.m, c.n, alpha, a.ptr, c.lda, c.ldb);
        @memcpy(b, a[0..b.len]);
    } else try ext.omatcopy(ob.io, T, tr, c.m, c.n, alpha, a.ptr, c.lda, b.ptr, c.ldb);
    try expectNorm(c.name, T, matDiff(T, out_rows, out_cols, b, want, ldb), tol(T, single_eps, double_eps));
}

pub fn matcopy() !void {
    var prng = std.Random.DefaultPrng.init(0xC1);
    inline for (.{ .{ f32, &omatcopy_s, &imatcopy_s }, .{ f64, &omatcopy_d, &imatcopy_d }, .{ C32, &omatcopy_c, &imatcopy_c }, .{ C64, &omatcopy_z, &imatcopy_z } }) |p| {
        for (p[1]) |c| try matcopyCase(p[0], false, c, prng.random());
        for (p[2]) |c| try matcopyCase(p[0], true, c, prng.random());
    }
}

test "utest_ext: omatcopy and imatcopy against the definition, and their bad arguments" {
    try matcopy();
}

/// `test_somatcopy.c`, column-major.
const omatcopy_s = [_]Case{
    .{ .name = "colmajor_trans_col_100_row_100", .kind = .run, .trans = 'T', .m = 100, .n = 100, .lda = 100, .ldb = 100, .alpha = .{ 1.0, 0 } },
    .{ .name = "colmajor_notrans_col_100_row_100", .kind = .run, .trans = 'N', .m = 100, .n = 100, .lda = 100, .ldb = 100, .alpha = .{ 1.0, 0 } },
    .{ .name = "colmajor_trans_col_50_row_100", .kind = .run, .trans = 'T', .m = 100, .n = 50, .lda = 100, .ldb = 50, .alpha = .{ 2.0, 0 } },
    .{ .name = "colmajor_notrans_col_50_row_100", .kind = .run, .trans = 'N', .m = 100, .n = 50, .lda = 100, .ldb = 100, .alpha = .{ 2.0, 0 } },
    .{ .name = "colmajor_trans_col_100_row_50", .kind = .run, .trans = 'T', .m = 50, .n = 100, .lda = 50, .ldb = 100, .alpha = .{ 0.0, 0 } },
    .{ .name = "colmajor_notrans_col_100_row_50", .kind = .run, .trans = 'N', .m = 50, .n = 100, .lda = 50, .ldb = 50, .alpha = .{ 0.0, 0 } },
    .{ .name = "xerbla_invalid_trans", .kind = .bad, .trans = 'O', .m = 100, .n = 100, .lda = 100, .ldb = 100, .info = 2 },
    .{ .name = "xerbla_colmajor_invalid_lda", .kind = .bad, .trans = 'T', .m = 100, .n = 50, .lda = 50, .ldb = 100, .info = 7 },
    .{ .name = "xerbla_colmajor_notrans_invalid_ldb", .kind = .bad, .trans = 'N', .m = 100, .n = 50, .lda = 100, .ldb = 50, .info = 9 },
    .{ .name = "xerbla_colmajor_trans_invalid_ldb", .kind = .bad, .trans = 'T', .m = 50, .n = 100, .lda = 100, .ldb = 50, .info = 9 },
};

/// `test_domatcopy.c`, column-major.
const omatcopy_d = [_]Case{
    .{ .name = "colmajor_trans_col_100_row_100", .kind = .run, .trans = 'T', .m = 100, .n = 100, .lda = 100, .ldb = 100, .alpha = .{ 1.0, 0 } },
    .{ .name = "colmajor_notrans_col_100_row_100", .kind = .run, .trans = 'N', .m = 100, .n = 100, .lda = 100, .ldb = 100, .alpha = .{ 1.0, 0 } },
    .{ .name = "colmajor_trans_col_50_row_100", .kind = .run, .trans = 'T', .m = 100, .n = 50, .lda = 100, .ldb = 50, .alpha = .{ 2.0, 0 } },
    .{ .name = "colmajor_notrans_col_50_row_100", .kind = .run, .trans = 'N', .m = 100, .n = 50, .lda = 100, .ldb = 100, .alpha = .{ 2.0, 0 } },
    .{ .name = "colmajor_trans_col_100_row_50", .kind = .run, .trans = 'T', .m = 50, .n = 100, .lda = 50, .ldb = 100, .alpha = .{ 0.0, 0 } },
    .{ .name = "colmajor_notrans_col_100_row_50", .kind = .run, .trans = 'N', .m = 50, .n = 100, .lda = 50, .ldb = 50, .alpha = .{ 0.0, 0 } },
    // large enough to thread (OpenBLAS f2ccc5045)
    .{ .name = "colmajor_trans_col_600_row_500", .kind = .run, .trans = 'T', .m = 500, .n = 600, .lda = 500, .ldb = 600, .alpha = .{ 2.0, 0 } },
    .{ .name = "colmajor_notrans_col_600_row_500", .kind = .run, .trans = 'N', .m = 500, .n = 600, .lda = 500, .ldb = 500, .alpha = .{ 2.0, 0 } },
    .{ .name = "xerbla_invalid_trans", .kind = .bad, .trans = 'O', .m = 100, .n = 100, .lda = 100, .ldb = 100, .info = 2 },
    .{ .name = "xerbla_colmajor_invalid_lda", .kind = .bad, .trans = 'T', .m = 100, .n = 50, .lda = 50, .ldb = 100, .info = 7 },
    .{ .name = "xerbla_colmajor_notrans_invalid_ldb", .kind = .bad, .trans = 'N', .m = 100, .n = 50, .lda = 100, .ldb = 50, .info = 9 },
    .{ .name = "xerbla_colmajor_trans_invalid_ldb", .kind = .bad, .trans = 'T', .m = 50, .n = 100, .lda = 100, .ldb = 50, .info = 9 },
};

/// `test_comatcopy.c`, column-major.
const omatcopy_c = [_]Case{
    .{ .name = "colmajor_notrans_col_50_row_100", .kind = .run, .trans = 'N', .m = 100, .n = 50, .lda = 100, .ldb = 100, .alpha = .{ 1.0, 2.0 } },
    .{ .name = "colmajor_trans_col_50_row_100", .kind = .run, .trans = 'T', .m = 100, .n = 50, .lda = 100, .ldb = 50, .alpha = .{ -1.0, 2.0 } },
    .{ .name = "colmajor_conj_col_100_row_100", .kind = .run, .trans = 'R', .m = 100, .n = 100, .lda = 100, .ldb = 100, .alpha = .{ 1.0, 2.0 } },
    .{ .name = "colmajor_conjtrnas_col_100_row_100", .kind = .run, .trans = 'C', .m = 100, .n = 100, .lda = 100, .ldb = 100, .alpha = .{ 2.0, 1.0 } },
    .{ .name = "xerbla_invalid_trans", .kind = .bad, .trans = 'O', .m = 100, .n = 100, .lda = 100, .ldb = 100, .info = 2 },
    .{ .name = "xerbla_colmajor_invalid_lda", .kind = .bad, .trans = 'T', .m = 100, .n = 50, .lda = 50, .ldb = 100, .info = 7 },
    .{ .name = "xerbla_colmajor_notrans_invalid_ldb", .kind = .bad, .trans = 'N', .m = 100, .n = 50, .lda = 100, .ldb = 50, .info = 9 },
    .{ .name = "xerbla_colmajor_trans_invalid_ldb", .kind = .bad, .trans = 'T', .m = 50, .n = 100, .lda = 100, .ldb = 50, .info = 9 },
    .{ .name = "xerbla_colmajor_conj_invalid_ldb", .kind = .bad, .trans = 'R', .m = 100, .n = 50, .lda = 100, .ldb = 50, .info = 9 },
    .{ .name = "xerbla_colmajor_transconj_invalid_ldb", .kind = .bad, .trans = 'C', .m = 50, .n = 100, .lda = 100, .ldb = 50, .info = 9 },
};

/// `test_zomatcopy.c`, column-major.
const omatcopy_z = [_]Case{
    .{ .name = "colmajor_notrans_col_50_row_100", .kind = .run, .trans = 'N', .m = 100, .n = 50, .lda = 100, .ldb = 100, .alpha = .{ 1.0, 2.0 } },
    .{ .name = "colmajor_trans_col_50_row_100", .kind = .run, .trans = 'T', .m = 100, .n = 50, .lda = 100, .ldb = 50, .alpha = .{ -1.0, 2.0 } },
    .{ .name = "colmajor_conj_col_100_row_100", .kind = .run, .trans = 'R', .m = 100, .n = 100, .lda = 100, .ldb = 100, .alpha = .{ 1.0, 2.0 } },
    .{ .name = "colmajor_conjtrnas_col_100_row_100", .kind = .run, .trans = 'C', .m = 100, .n = 100, .lda = 100, .ldb = 100, .alpha = .{ 2.0, 1.0 } },
    .{ .name = "xerbla_invalid_trans", .kind = .bad, .trans = 'O', .m = 100, .n = 100, .lda = 100, .ldb = 100, .info = 2 },
    .{ .name = "xerbla_colmajor_invalid_lda", .kind = .bad, .trans = 'T', .m = 100, .n = 50, .lda = 50, .ldb = 100, .info = 7 },
    .{ .name = "xerbla_colmajor_notrans_invalid_ldb", .kind = .bad, .trans = 'N', .m = 100, .n = 50, .lda = 100, .ldb = 50, .info = 9 },
    .{ .name = "xerbla_colmajor_trans_invalid_ldb", .kind = .bad, .trans = 'T', .m = 50, .n = 100, .lda = 100, .ldb = 50, .info = 9 },
    .{ .name = "xerbla_colmajor_conj_invalid_ldb", .kind = .bad, .trans = 'R', .m = 100, .n = 50, .lda = 100, .ldb = 50, .info = 9 },
    .{ .name = "xerbla_colmajor_transconj_invalid_ldb", .kind = .bad, .trans = 'C', .m = 50, .n = 100, .lda = 100, .ldb = 50, .info = 9 },
};

/// `test_simatcopy.c`, column-major.
const imatcopy_s = [_]Case{
    .{ .name = "colmajor_trans_col_100_row_100_alpha_one", .kind = .run, .trans = 'T', .m = 100, .n = 100, .lda = 100, .ldb = 100, .alpha = .{ 1.0, 0 } },
    .{ .name = "colmajor_notrans_col_100_row_100_alpha_one", .kind = .run, .trans = 'N', .m = 100, .n = 100, .lda = 100, .ldb = 100, .alpha = .{ 1.0, 0 } },
    .{ .name = "colmajor_trans_col_100_row_100_alpha_zero", .kind = .run, .trans = 'T', .m = 100, .n = 100, .lda = 100, .ldb = 100, .alpha = .{ 0.0, 0 } },
    .{ .name = "colmajor_notrans_col_100_row_100_alpha_zero", .kind = .run, .trans = 'N', .m = 100, .n = 100, .lda = 100, .ldb = 100, .alpha = .{ 0.0, 0 } },
    .{ .name = "colmajor_trans_col_100_row_100", .kind = .run, .trans = 'T', .m = 100, .n = 100, .lda = 100, .ldb = 100, .alpha = .{ 2.0, 0 } },
    .{ .name = "colmajor_notrans_col_100_row_100", .kind = .run, .trans = 'N', .m = 100, .n = 100, .lda = 100, .ldb = 100, .alpha = .{ 2.0, 0 } },
    .{ .name = "colmajor_trans_col_50_row_100_alpha_one", .kind = .run, .trans = 'T', .m = 100, .n = 50, .lda = 100, .ldb = 50, .alpha = .{ 1.0, 0 } },
    .{ .name = "colmajor_notrans_col_50_row_100_alpha_one", .kind = .run, .trans = 'N', .m = 100, .n = 50, .lda = 100, .ldb = 100, .alpha = .{ 1.0, 0 } },
    .{ .name = "colmajor_trans_col_50_row_100_alpha_zero", .kind = .run, .trans = 'T', .m = 100, .n = 50, .lda = 100, .ldb = 50, .alpha = .{ 0.0, 0 } },
    .{ .name = "colmajor_notrans_col_50_row_100_alpha_zero", .kind = .run, .trans = 'N', .m = 100, .n = 50, .lda = 100, .ldb = 100, .alpha = .{ 0.0, 0 } },
    .{ .name = "colmajor_trans_col_50_row_100", .kind = .run, .trans = 'T', .m = 100, .n = 50, .lda = 100, .ldb = 50, .alpha = .{ 2.0, 0 } },
    .{ .name = "colmajor_notrans_col_50_row_100", .kind = .run, .trans = 'N', .m = 100, .n = 50, .lda = 100, .ldb = 100, .alpha = .{ 2.0, 0 } },
    .{ .name = "xerbla_invalid_trans", .kind = .bad, .trans = 'O', .m = 100, .n = 100, .lda = 100, .ldb = 100, .info = 2 },
    .{ .name = "xerbla_colmajor_invalid_lda", .kind = .bad, .trans = 'T', .m = 100, .n = 50, .lda = 50, .ldb = 100, .info = 7 },
    .{ .name = "xerbla_colmajor_notrans_invalid_ldb", .kind = .bad, .trans = 'N', .m = 100, .n = 50, .lda = 100, .ldb = 50, .info = 8 },
    .{ .name = "xerbla_colmajor_trans_invalid_ldb", .kind = .bad, .trans = 'T', .m = 50, .n = 100, .lda = 100, .ldb = 50, .info = 8 },
};

/// `test_dimatcopy.c`, column-major.
const imatcopy_d = [_]Case{
    .{ .name = "colmajor_trans_col_100_row_100_alpha_one", .kind = .run, .trans = 'T', .m = 100, .n = 100, .lda = 100, .ldb = 100, .alpha = .{ 1.0, 0 } },
    .{ .name = "colmajor_notrans_col_100_row_100_alpha_one", .kind = .run, .trans = 'N', .m = 100, .n = 100, .lda = 100, .ldb = 100, .alpha = .{ 1.0, 0 } },
    .{ .name = "colmajor_trans_col_100_row_100_alpha_zero", .kind = .run, .trans = 'T', .m = 100, .n = 100, .lda = 100, .ldb = 100, .alpha = .{ 0.0, 0 } },
    .{ .name = "colmajor_notrans_col_100_row_100_alpha_zero", .kind = .run, .trans = 'N', .m = 100, .n = 100, .lda = 100, .ldb = 100, .alpha = .{ 0.0, 0 } },
    .{ .name = "colmajor_trans_col_100_row_100", .kind = .run, .trans = 'T', .m = 100, .n = 100, .lda = 100, .ldb = 100, .alpha = .{ 2.0, 0 } },
    .{ .name = "colmajor_notrans_col_100_row_100", .kind = .run, .trans = 'N', .m = 100, .n = 100, .lda = 100, .ldb = 100, .alpha = .{ 2.0, 0 } },
    .{ .name = "colmajor_trans_col_50_row_100_alpha_one", .kind = .run, .trans = 'T', .m = 100, .n = 50, .lda = 100, .ldb = 50, .alpha = .{ 1.0, 0 } },
    .{ .name = "colmajor_notrans_col_50_row_100_alpha_one", .kind = .run, .trans = 'N', .m = 100, .n = 50, .lda = 100, .ldb = 100, .alpha = .{ 1.0, 0 } },
    .{ .name = "colmajor_trans_col_50_row_100_alpha_zero", .kind = .run, .trans = 'T', .m = 100, .n = 50, .lda = 100, .ldb = 50, .alpha = .{ 0.0, 0 } },
    .{ .name = "colmajor_notrans_col_50_row_100_alpha_zero", .kind = .run, .trans = 'N', .m = 100, .n = 50, .lda = 100, .ldb = 100, .alpha = .{ 0.0, 0 } },
    .{ .name = "colmajor_trans_col_50_row_100", .kind = .run, .trans = 'T', .m = 100, .n = 50, .lda = 100, .ldb = 50, .alpha = .{ 2.0, 0 } },
    .{ .name = "colmajor_notrans_col_50_row_100", .kind = .run, .trans = 'N', .m = 100, .n = 50, .lda = 100, .ldb = 100, .alpha = .{ 2.0, 0 } },
    .{ .name = "xerbla_invalid_trans", .kind = .bad, .trans = 'O', .m = 100, .n = 100, .lda = 100, .ldb = 100, .info = 2 },
    .{ .name = "xerbla_colmajor_invalid_lda", .kind = .bad, .trans = 'T', .m = 100, .n = 50, .lda = 50, .ldb = 100, .info = 7 },
    .{ .name = "xerbla_colmajor_notrans_invalid_ldb", .kind = .bad, .trans = 'N', .m = 100, .n = 50, .lda = 100, .ldb = 50, .info = 8 },
    .{ .name = "xerbla_colmajor_trans_invalid_ldb", .kind = .bad, .trans = 'T', .m = 50, .n = 100, .lda = 100, .ldb = 50, .info = 8 },
};

/// `test_cimatcopy.c`, column-major.
const imatcopy_c = [_]Case{
    .{ .name = "colmajor_trans_col_100_row_100", .kind = .run, .trans = 'T', .m = 100, .n = 100, .lda = 100, .ldb = 100, .alpha = .{ 1.0, 2.0 } },
    .{ .name = "colmajor_notrans_col_100_row_100", .kind = .run, .trans = 'N', .m = 100, .n = 100, .lda = 100, .ldb = 100, .alpha = .{ -3.0, 1.0 } },
    .{ .name = "colmajor_conj_col_100_row_100", .kind = .run, .trans = 'R', .m = 100, .n = 100, .lda = 100, .ldb = 100, .alpha = .{ 1.0, 2.0 } },
    .{ .name = "colmajor_conjtrans_col_100_row_100", .kind = .run, .trans = 'C', .m = 100, .n = 100, .lda = 100, .ldb = 100, .alpha = .{ 2.0, 1.0 } },
    .{ .name = "colmajor_trans_col_50_row_100", .kind = .run, .trans = 'T', .m = 100, .n = 50, .lda = 100, .ldb = 50, .alpha = .{ 1.0, 2.0 } },
    .{ .name = "colmajor_notrans_col_100_row_50", .kind = .run, .trans = 'N', .m = 50, .n = 100, .lda = 50, .ldb = 50, .alpha = .{ 1.0, 2.0 } },
    .{ .name = "colmajor_conjtrans_col_50_row_100", .kind = .run, .trans = 'C', .m = 100, .n = 50, .lda = 100, .ldb = 50, .alpha = .{ 1.0, 1.0 } },
    .{ .name = "colmajor_conj_col_100_row_50", .kind = .run, .trans = 'R', .m = 50, .n = 100, .lda = 50, .ldb = 50, .alpha = .{ 1.0, 2.0 } },
    .{ .name = "c_api_colmajor_trans_col_100_row_100", .kind = .run, .trans = 'T', .m = 100, .n = 100, .lda = 100, .ldb = 100, .alpha = .{ 3.0, 2.0 } },
    .{ .name = "c_api_colmajor_notrans_col_100_row_100", .kind = .run, .trans = 'N', .m = 100, .n = 100, .lda = 100, .ldb = 100, .alpha = .{ 3.0, 1.5 } },
    .{ .name = "xerbla_invalid_trans", .kind = .bad, .trans = 'O', .m = 100, .n = 100, .lda = 100, .ldb = 100, .info = 2 },
    .{ .name = "xerbla_colmajor_invalid_lda", .kind = .bad, .trans = 'T', .m = 100, .n = 50, .lda = 50, .ldb = 100, .info = 7 },
    .{ .name = "xerbla_colmajor_notrans_invalid_ldb", .kind = .bad, .trans = 'N', .m = 100, .n = 50, .lda = 100, .ldb = 50, .info = 9 },
    .{ .name = "xerbla_colmajor_trans_invalid_ldb", .kind = .bad, .trans = 'T', .m = 50, .n = 100, .lda = 100, .ldb = 50, .info = 9 },
};

/// `test_zimatcopy.c`, column-major.
const imatcopy_z = [_]Case{
    .{ .name = "colmajor_trans_col_100_row_100", .kind = .run, .trans = 'T', .m = 100, .n = 100, .lda = 100, .ldb = 100, .alpha = .{ 1.0, 2.0 } },
    .{ .name = "colmajor_notrans_col_100_row_100", .kind = .run, .trans = 'N', .m = 100, .n = 100, .lda = 100, .ldb = 100, .alpha = .{ -3.0, 1.0 } },
    .{ .name = "colmajor_conj_col_100_row_100", .kind = .run, .trans = 'R', .m = 100, .n = 100, .lda = 100, .ldb = 100, .alpha = .{ 1.0, 2.0 } },
    .{ .name = "colmajor_conjtrans_col_100_row_100", .kind = .run, .trans = 'C', .m = 100, .n = 100, .lda = 100, .ldb = 100, .alpha = .{ 2.0, 1.0 } },
    .{ .name = "colmajor_trans_col_50_row_100", .kind = .run, .trans = 'T', .m = 100, .n = 50, .lda = 100, .ldb = 50, .alpha = .{ 1.0, 2.0 } },
    .{ .name = "colmajor_notrans_col_100_row_50", .kind = .run, .trans = 'N', .m = 50, .n = 100, .lda = 50, .ldb = 50, .alpha = .{ 1.0, 2.0 } },
    .{ .name = "colmajor_conjtrans_col_50_row_100", .kind = .run, .trans = 'C', .m = 100, .n = 50, .lda = 100, .ldb = 50, .alpha = .{ 1.0, 1.0 } },
    .{ .name = "colmajor_conj_col_100_row_50", .kind = .run, .trans = 'R', .m = 50, .n = 100, .lda = 50, .ldb = 50, .alpha = .{ 1.0, 2.0 } },
    .{ .name = "c_api_colmajor_trans_col_100_row_100", .kind = .run, .trans = 'T', .m = 100, .n = 100, .lda = 100, .ldb = 100, .alpha = .{ 3.0, 2.0 } },
    .{ .name = "c_api_colmajor_notrans_col_100_row_100", .kind = .run, .trans = 'N', .m = 100, .n = 100, .lda = 100, .ldb = 100, .alpha = .{ 3.0, 1.5 } },
    .{ .name = "xerbla_invalid_trans", .kind = .bad, .trans = 'O', .m = 100, .n = 100, .lda = 100, .ldb = 100, .info = 2 },
    .{ .name = "xerbla_colmajor_invalid_lda", .kind = .bad, .trans = 'T', .m = 100, .n = 50, .lda = 50, .ldb = 100, .info = 7 },
    .{ .name = "xerbla_colmajor_notrans_invalid_ldb", .kind = .bad, .trans = 'N', .m = 100, .n = 50, .lda = 100, .ldb = 50, .info = 9 },
    .{ .name = "xerbla_colmajor_trans_invalid_ldb", .kind = .bad, .trans = 'T', .m = 50, .n = 100, .lda = 100, .ldb = 50, .info = 9 },
};
