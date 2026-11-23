//! The level-2 cases: `gbmv`, `gemv`'s extra letters, complex `sbmv` and `spmv`, `trmv` and `trsv` with `R`.

const std = @import("std");
const ob = @import("../__root__.zig");
const blas = @import("../../../__root__.zig");
const l2 = @import("../../../blas2/core.zig");
const util = @import("util.zig");

const testing = std.testing;
const C32 = util.C32;
const C64 = util.C64;
const single_eps = util.single_eps;
const double_eps = util.double_eps;
const single_tol = util.single_tol;
const double_tol = util.double_tol;
const single_eps_zgemv = util.single_eps_zgemv;
const double_eps_zgemv = util.double_eps_zgemv;
const tol = util.tol;
const Case = util.Case;
const val = util.val;
const fill = util.fill;
const span = util.span;
const diffNorm = util.diffNorm;
const norm2 = util.norm2;
const expectNorm = util.expectNorm;
const conjugate = util.conjugate;
const transOf = util.transOf;
const transposed = util.transposed;

/// `y := alpha*op(A)*x + beta*y` by the definition, for `trans` any of the
/// eight letters: the sources' `?gemv_n_trusted` and its kin.
fn gemvRef(comptime T: type, trans: u8, m: usize, n: usize, alpha: T, a: []const T, lda: usize, x: []const T, incx: isize, beta: T, y: []T, incy: isize) void {
    const tr = transposed(trans);
    const conj_a = trans == 'R' or trans == 'C' or trans == 'S' or trans == 'D';
    const conj_x = trans == 'O' or trans == 'U' or trans == 'S' or trans == 'D';
    const leny = if (tr) n else m;
    const lenx = if (tr) m else n;
    const ix: usize = @abs(incx);
    const iy: usize = @abs(incy);
    for (0..leny) |r| {
        var acc = ob.zero(T);
        for (0..lenx) |q| {
            var e = if (tr) a[q + r * lda] else a[r + q * lda];
            var xv = x[(if (incx < 0) lenx - 1 - q else q) * ix];
            if (conj_a) e = ob.conj(e);
            if (conj_x) xv = ob.conj(xv);
            acc = ob.add(acc, ob.mul(e, xv));
        }
        const k = (if (incy < 0) leny - 1 - r else r) * iy;
        y[k] = ob.add(ob.mul(alpha, acc), ob.mul(beta, y[k]));
    }
}

pub fn gbmv() !void {
    var prng = std.Random.DefaultPrng.init(0xB0);
    const al = testing.allocator;
    inline for (.{ C32, C64 }) |T| for (gbmv_cases) |c| {
        const m: usize = @intCast(c.m);
        const n: usize = @intCast(c.n);
        const lda: usize = @intCast(c.lda);
        const kl: usize = @intCast(c.kl);
        const ku: usize = @intCast(c.ku);
        const tr = transposed(c.trans);
        const lenb: isize = if (tr) c.m else c.n;
        const lenc: isize = if (tr) c.n else c.m;
        // the band of a random m by n matrix, full and in band storage
        const full = try al.alloc(T, m * n);
        defer al.free(full);
        const band = try al.alloc(T, lda * n);
        defer al.free(band);
        @memset(band, ob.zero(T));
        fill(T, full, prng.random());
        for (0..n) |j| for (0..m) |i| {
            if (i > kl + j or j > ku + i) full[i + j * m] = ob.zero(T) else band[ku + i - j + j * lda] = full[i + j * m];
        };
        const b = try al.alloc(T, span(lenb, c.inc_b));
        defer al.free(b);
        const y = try al.alloc(T, span(lenc, c.inc_c));
        defer al.free(y);
        fill(T, b, prng.random());
        fill(T, y, prng.random());
        const yv = try al.dupe(T, y);
        defer al.free(yv);
        try l2.gemv(al, ob.io, T, transOf(c.trans), c.m, c.n, val(T, c.alpha), full.ptr, c.m, b.ptr, c.inc_b, val(T, c.beta), yv.ptr, c.inc_c);
        try l2.gbmv(al, ob.io, T, transOf(c.trans), c.m, c.n, c.kl, c.ku, val(T, c.alpha), band.ptr, c.lda, b.ptr, c.inc_b, val(T, c.beta), y.ptr, c.inc_c);
        try expectNorm(c.name, T, diffNorm(T, lenc, y, yv, c.inc_c), tol(T, single_tol, double_tol));
    };
}

test "utest_ext: gbmv against gemv on the band" {
    try gbmv();
}

/// `gemv` against `gemvRef` for one case of `test_[cz]gemv_[nt].c`.
fn gemvCase(comptime T: type, c: Case, r: std.Random, t: f64) !void {
    const al = testing.allocator;
    const m: usize = @intCast(c.m);
    const n: usize = @intCast(c.n);
    const lda: usize = @intCast(c.lda);
    const tr = transposed(c.trans);
    const ib = if (c.inc_b != 1 or c.inc_x == 1) c.inc_b else c.inc_x;
    const ic = if (c.inc_c != 1 or c.inc_y == 1) c.inc_c else c.inc_y;
    const lenx: isize = if (tr) c.m else c.n;
    const leny: isize = if (tr) c.n else c.m;
    const a = try al.alloc(T, lda * n);
    defer al.free(a);
    const x = try al.alloc(T, span(lenx, ib));
    defer al.free(x);
    const y = try al.alloc(T, span(leny, ic));
    defer al.free(y);
    fill(T, a, r);
    fill(T, x, r);
    fill(T, y, r);
    const yv = try al.dupe(T, y);
    defer al.free(yv);
    gemvRef(T, c.trans, m, n, val(T, c.alpha), a, lda, x, ib, val(T, c.beta), yv, ic);
    try l2.gemv(al, ob.io, T, transOf(c.trans), c.m, c.n, val(T, c.alpha), a.ptr, c.lda, x.ptr, ib, val(T, c.beta), y.ptr, ic);
    try expectNorm(c.name, T, diffNorm(T, leny, y, yv, ic), t);
}

pub fn gemvOS() !void {
    var prng = std.Random.DefaultPrng.init(0xB1);
    inline for (.{ C32, C64 }) |T| for (gemv_n_cases) |c| try gemvCase(T, c, prng.random(), tol(T, single_tol, double_tol));
}

test "utest_ext: gemv with O and S against the plain product" {
    try gemvOS();
}

pub fn gemvTCUD() !void {
    var prng = std.Random.DefaultPrng.init(0xB2);
    const al = testing.allocator;
    inline for (.{ C32, C64 }) |T| for (gemv_t_cases) |c| {
        if (c.kind == .run) {
            try gemvCase(T, c, prng.random(), tol(T, single_eps_zgemv, double_eps_zgemv));
            continue;
        }
        // `Z` is no `Trans`
        if (c.info == 1) continue;
        var a = [_]T{val(T, .{ 1, 1 })};
        var x = [_]T{val(T, .{ 1, 1 })};
        var y = [_]T{val(T, .{ 1, 1 })};
        const want: anyerror = switch (c.info) {
            2 => error.InvalidM,
            3 => error.InvalidN,
            6 => error.InvalidLda,
            8 => error.InvalidIncx,
            else => error.InvalidIncy,
        };
        try ob.expectError(want, l2.gemv(al, ob.io, T, transOf(c.trans), c.m, c.n, val(T, .{ 1, 1 }), &a, c.lda, &x, c.inc_x, val(T, .{ 1, 1 }), &y, c.inc_y), c.name);
    };
}

test "utest_ext: gemv with T, C, U, D against the plain product, and its bad arguments" {
    try gemvTCUD();
}

/// The bad-argument cases of `test_[cz]sbmv.c` and `test_[cz]spmv.c`.
fn symBad(comptime T: type, comptime pack: bool, c: Case) !void {
    const al = testing.allocator;
    if (c.info == 1) return ob.expect(ob.letter(blas.Uplo, c.uplo) == null, c.name);
    var a: [2]T = @splat(ob.zero(T));
    var x: [2]T = @splat(ob.zero(T));
    var y: [2]T = @splat(ob.zero(T));
    const n = if (pack) c.N else c.n;
    const want: anyerror = switch (c.info) {
        2 => error.InvalidN,
        3 => error.InvalidK,
        6 => if (pack) error.InvalidIncx else error.InvalidLda,
        8 => error.InvalidIncx,
        else => error.InvalidIncy,
    };
    const one = val(T, .{ 1, 1 });
    const res = if (pack)
        l2.spmv(al, ob.io, T, .u, n, one, &a, &x, c.inc_b, ob.zero(T), &y, c.inc_c)
    else
        l2.sbmv(al, ob.io, T, .u, n, c.k, one, &a, c.lda, &x, c.inc_b, ob.zero(T), &y, c.inc_c);
    try ob.expectError(want, res, c.name);
}

pub fn sbmvSpmv() !void {
    var prng = std.Random.DefaultPrng.init(0xB3);
    const al = testing.allocator;
    inline for (.{ C32, C64 }) |T| inline for (.{ false, true }) |pack| {
        const cases: []const Case = if (pack) &spmv_cases else &sbmv_cases;
        for (cases) |c| {
            if (c.kind == .bad) {
                try symBad(T, pack, c);
                continue;
            }
            const n: usize = @intCast(if (pack) c.N else c.n);
            const k: usize = if (pack) n else @intCast(c.k);
            const up = c.uplo == 'U';
            // a random symmetric matrix of k bands (all of them, packed)
            const full = try al.alloc(T, @max(n * n, 1));
            defer al.free(full);
            @memset(full, ob.zero(T));
            for (0..n) |j| for (0..n) |i| if (i <= j and j - i <= k) {
                var v: [1]T = undefined;
                fill(T, &v, prng.random());
                full[i + j * n] = v[0];
                full[j + i * n] = v[0];
            };
            const lda: usize = if (pack) 1 else @intCast(c.lda);
            const store = try al.alloc(T, @max(if (pack) n * (n + 1) / 2 else lda * n, 1));
            defer al.free(store);
            var p: usize = 0;
            for (0..n) |j| for (0..n) |i| {
                const in = if (up) i <= j and j - i <= k else i >= j and i - j <= k;
                if (!in) continue;
                if (pack) {
                    store[p] = full[i + j * n];
                    p += 1;
                } else store[(if (up) k + i - j else i - j) + j * lda] = full[i + j * n];
            };
            const b = try al.alloc(T, span(@intCast(n), c.inc_b));
            defer al.free(b);
            const y = try al.alloc(T, span(@intCast(n), c.inc_c));
            defer al.free(y);
            fill(T, b, prng.random());
            fill(T, y, prng.random());
            const yv = try al.dupe(T, y);
            defer al.free(yv);
            const in: isize = @intCast(n);
            try l2.gemv(al, ob.io, T, .n, in, in, val(T, c.alpha), full.ptr, @max(in, 1), b.ptr, c.inc_b, val(T, c.beta), yv.ptr, c.inc_c);
            if (pack) {
                try l2.spmv(al, ob.io, T, if (up) .u else .l, in, val(T, c.alpha), store.ptr, b.ptr, c.inc_b, val(T, c.beta), y.ptr, c.inc_c);
            } else {
                try l2.sbmv(al, ob.io, T, if (up) .u else .l, in, c.k, val(T, c.alpha), store.ptr, c.lda, b.ptr, c.inc_b, val(T, c.beta), y.ptr, c.inc_c);
            }
            const norm = diffNorm(T, in, y, yv, c.inc_c);
            if (c.kind == .nan) {
                try ob.expect(norm == norm, c.name);
            } else try expectNorm(c.name, T, norm, if (pack) tol(T, single_tol, double_tol) else tol(T, single_eps_zgemv, double_eps_zgemv));
        }
    };
}

test "utest_ext: complex symmetric sbmv and spmv against gemv" {
    try sbmvSpmv();
}

pub fn trmvTrsv() !void {
    var prng = std.Random.DefaultPrng.init(0xB4);
    const al = testing.allocator;
    inline for (.{ C32, C64 }) |T| inline for (.{ false, true }) |solve| {
        const cases: []const Case = if (solve) &trsv_cases else &trmv_cases;
        for (cases) |c| {
            const n: usize = @intCast(c.n);
            const lda: usize = @intCast(c.lda);
            const a = try al.alloc(T, lda * n);
            defer al.free(a);
            fill(T, a, prng.random());
            const av = try al.dupe(T, a);
            defer al.free(av);
            conjugate(T, av);
            const x = try al.alloc(T, span(c.n, c.incx));
            defer al.free(x);
            fill(T, x, prng.random());
            const xv = try al.dupe(T, x);
            defer al.free(xv);
            const u = ob.letter(blas.Uplo, c.uplo).?;
            const d = ob.letter(blas.Diag, c.diag).?;
            const f = if (solve) l2.trsv else l2.trmv;
            try f(al, ob.io, T, u, .n, d, c.n, av.ptr, c.lda, xv.ptr, c.incx);
            try f(al, ob.io, T, u, ob.letter(blas.Trans, c.trans).?, d, c.n, a.ptr, c.lda, x.ptr, c.incx);
            // relative to the result, which a solve grows to 1e32 in single
            // precision; there `test_ctrsv.c`'s DOUBLE_EPS is single's, as the
            // R and N paths round differently on the older cores' kernels
            const t = tol(T, single_eps, double_eps);
            const scale = @max(1, norm2(T, c.n, xv, c.incx));
            try expectNorm(c.name, T, diffNorm(T, c.n, x, xv, c.incx) / scale, t);
        }
    };
}

test "utest_ext: trmv and trsv with R against N on the conjugate" {
    try trmvTrsv();
}

/// `test_zgbmv.c (test_cgbmv.c)`.
const gbmv_cases = [_]Case{
    .{ .name = "trans_D", .kind = .run, .m = 50, .n = 25, .inc_b = 1, .inc_c = 1, .kl = 20, .ku = 11, .lda = 50, .trans = 'D', .alpha = .{ 7.0, 1.0 }, .beta = .{ 1.5, -1.5 } },
    .{ .name = "trans_O", .kind = .run, .m = 50, .n = 25, .inc_b = 1, .inc_c = 1, .kl = 20, .ku = 10, .lda = 50, .trans = 'O', .alpha = .{ 7.0, 1.0 }, .beta = .{ 1.5, -1.5 } },
    .{ .name = "trans_S", .kind = .run, .m = 50, .n = 25, .inc_b = 1, .inc_c = 1, .kl = 6, .ku = 9, .lda = 50, .trans = 'S', .alpha = .{ 7.0, 1.0 }, .beta = .{ 1.5, -1.5 } },
    .{ .name = "trans_U", .kind = .run, .m = 25, .n = 50, .inc_b = 1, .inc_c = 1, .kl = 7, .ku = 11, .lda = 19, .trans = 'U', .alpha = .{ 7.0, 1.0 }, .beta = .{ 1.5, -1.5 } },
    .{ .name = "trans_C", .kind = .run, .m = 50, .n = 25, .inc_b = 1, .inc_c = 1, .kl = 20, .ku = 11, .lda = 50, .trans = 'C', .alpha = .{ 7.0, 1.0 }, .beta = .{ 1.5, -1.5 } },
    .{ .name = "trans_R", .kind = .run, .m = 50, .n = 100, .inc_b = 1, .inc_c = 1, .kl = 20, .ku = 11, .lda = 50, .trans = 'R', .alpha = .{ 7.0, 1.0 }, .beta = .{ 1.5, -1.5 } },
};

/// `test_zgemv_n.c (test_cgemv_n.c)`.
const gemv_n_cases = [_]Case{
    .{ .name = "trans_o_square_matrix", .kind = .run, .n = 100, .m = 100, .lda = 100, .inc_b = 1, .inc_c = 1, .trans = 'O', .alpha = .{ 2.0, -1.0 }, .beta = .{ 1.4, 5.0 } },
    .{ .name = "trans_o_rectangular_matrix_rows_less_then_cols", .kind = .run, .n = 100, .m = 50, .lda = 50, .inc_b = 1, .inc_c = 1, .trans = 'O', .alpha = .{ 2.0, -1.0 }, .beta = .{ 1.4, 5.0 } },
    .{ .name = "trans_o_rectangular_matrix_cols_less_then_rows", .kind = .run, .n = 50, .m = 100, .lda = 100, .inc_b = 1, .inc_c = 1, .trans = 'O', .alpha = .{ 2.0, -1.0 }, .beta = .{ 1.4, 5.0 } },
    .{ .name = "trans_o_double_strides", .kind = .run, .n = 100, .m = 100, .lda = 100, .inc_b = 2, .inc_c = 2, .trans = 'O', .alpha = .{ 2.0, -1.0 }, .beta = .{ 1.4, 5.0 } },
    .{ .name = "trans_s_square_matrix", .kind = .run, .n = 100, .m = 100, .lda = 100, .inc_b = 1, .inc_c = 1, .trans = 'S', .alpha = .{ 1.0, 1.0 }, .beta = .{ 1.4, 5.0 } },
    .{ .name = "trans_s_rectangular_matrix_rows_less_then_cols", .kind = .run, .n = 100, .m = 50, .lda = 50, .inc_b = 1, .inc_c = 1, .trans = 'S', .alpha = .{ 2.0, -1.0 }, .beta = .{ 1.4, 5.0 } },
    .{ .name = "trans_s_rectangular_matrix_cols_less_then_rows", .kind = .run, .n = 50, .m = 100, .lda = 100, .inc_b = 1, .inc_c = 1, .trans = 'S', .alpha = .{ 2.0, -1.0 }, .beta = .{ 1.4, 0.0 } },
    .{ .name = "trans_s_double_strides", .kind = .run, .n = 100, .m = 100, .lda = 100, .inc_b = 2, .inc_c = 2, .trans = 'S', .alpha = .{ 2.0, -1.0 }, .beta = .{ 1.0, 5.0 } },
};

/// `test_zgemv_t.c (test_cgemv_t.c)`.
const gemv_t_cases = [_]Case{
    .{ .name = "colmajor_trans_col_100_row_100_inc_x_1_y_1", .kind = .run, .m = 100, .n = 100, .lda = 100, .inc_x = 1, .inc_y = 1, .order = 'C', .trans = 'T', .alpha = .{ 2.0, 1.0 }, .beta = .{ 1.0, 2.0 } },
    .{ .name = "colmajor_trans_col_100_row_100_inc_x_2_y_1", .kind = .run, .m = 100, .n = 100, .lda = 100, .inc_x = 2, .inc_y = 1, .order = 'C', .trans = 'T', .alpha = .{ 1.0, 1.0 }, .beta = .{ 2.0, 2.0 } },
    .{ .name = "colmajor_conjtrans_col_100_row_100_inc_x_1_y_1", .kind = .run, .m = 100, .n = 100, .lda = 100, .inc_x = 1, .inc_y = 1, .order = 'C', .trans = 'C', .alpha = .{ 2.0, 1.0 }, .beta = .{ 2.0, 1.0 } },
    .{ .name = "colmajor_conjtrans_col_100_row_100_inc_x_1_y_2", .kind = .run, .m = 100, .n = 100, .lda = 100, .inc_x = 1, .inc_y = 2, .order = 'C', .trans = 'C', .alpha = .{ 1.0, 1.0 }, .beta = .{ 1.0, 1.0 } },
    .{ .name = "colmajor_trans_x_conj_col_100_row_100_inc_x_1_y_1", .kind = .run, .m = 100, .n = 100, .lda = 100, .inc_x = 1, .inc_y = 1, .order = 'C', .trans = 'U', .alpha = .{ 1.0, 1.0 }, .beta = .{ 1.0, 1.0 } },
    .{ .name = "colmajor_trans_x_conj_col_100_row_100_inc_x_2_y_2", .kind = .run, .m = 100, .n = 100, .lda = 100, .inc_x = 2, .inc_y = 2, .order = 'C', .trans = 'U', .alpha = .{ 1.0, 2.0 }, .beta = .{ 1.0, 1.0 } },
    .{ .name = "colmajor_conjtrans_x_conj_col_100_row_100_inc_x_1_y_2", .kind = .run, .m = 100, .n = 100, .lda = 100, .inc_x = 1, .inc_y = 2, .order = 'C', .trans = 'D', .alpha = .{ 2.0, 1.0 }, .beta = .{ 1.0, 2.0 } },
    .{ .name = "xerbla_invalid_inc_y", .kind = .bad, .m = 1, .n = 1, .lda = 1, .inc_x = 1, .inc_y = 0, .order = 'C', .trans = 'T', .info = 11 },
    .{ .name = "xerbla_invalid_inc_x", .kind = .bad, .m = 1, .n = 1, .lda = 1, .inc_x = 0, .inc_y = 1, .order = 'C', .trans = 'T', .info = 8 },
    .{ .name = "xerbla_invalid_n", .kind = .bad, .m = 1, .n = -1, .lda = 1, .inc_x = 1, .inc_y = 1, .order = 'C', .trans = 'T', .info = 3 },
    .{ .name = "xerbla_invalid_m", .kind = .bad, .m = -1, .n = 1, .lda = 1, .inc_x = 1, .inc_y = 1, .order = 'C', .trans = 'T', .info = 2 },
    .{ .name = "xerbla_invalid_lda", .kind = .bad, .m = 1, .n = 1, .lda = -1, .inc_x = 1, .inc_y = 1, .order = 'C', .trans = 'T', .info = 6 },
    .{ .name = "xerbla_invalid_trans", .kind = .bad, .m = 1, .n = 1, .lda = 1, .inc_x = 1, .inc_y = 1, .order = 'C', .trans = 'Z', .info = 1 },
};

/// `test_zsbmv.c (test_csbmv.c)`.
const sbmv_cases = [_]Case{
    .{ .name = "upper_k_0_inc_b_1_inc_c_1_n_100", .kind = .run, .n = 100, .inc_b = 1, .inc_c = 1, .k = 0, .lda = 1, .ldm = 100, .uplo = 'U', .alpha = .{ 1.0, 1.0 }, .beta = .{ 1.0, 1.0 } },
    .{ .name = "upper_k_1_inc_b_1_inc_c_1_n_100", .kind = .run, .n = 100, .inc_b = 1, .inc_c = 1, .k = 1, .lda = 2, .ldm = 100, .uplo = 'U', .alpha = .{ 1.0, 1.0 }, .beta = .{ 1.0, 1.0 } },
    .{ .name = "upper_k_2_inc_b_1_inc_c_1_n_100", .kind = .run, .n = 100, .inc_b = 1, .inc_c = 1, .k = 2, .lda = 3, .ldm = 100, .uplo = 'U', .alpha = .{ 1.0, 1.0 }, .beta = .{ 1.0, 1.0 } },
    .{ .name = "upper_k_2_inc_b_2_inc_c_1_n_100", .kind = .run, .n = 100, .inc_b = 2, .inc_c = 1, .k = 2, .lda = 3, .ldm = 100, .uplo = 'U', .alpha = .{ 2.0, 1.0 }, .beta = .{ 2.0, 1.0 } },
    .{ .name = "upper_k_2_inc_b_2_inc_c_2_n_100", .kind = .run, .n = 100, .inc_b = 2, .inc_c = 2, .k = 2, .lda = 3, .ldm = 100, .uplo = 'U', .alpha = .{ 2.0, 1.0 }, .beta = .{ 2.0, 1.0 } },
    .{ .name = "lower_k_0_inc_b_1_inc_c_1_n_100", .kind = .run, .n = 100, .inc_b = 1, .inc_c = 1, .k = 0, .lda = 1, .ldm = 100, .uplo = 'L', .alpha = .{ 1.0, 1.0 }, .beta = .{ 1.0, 1.0 } },
    .{ .name = "lower_k_1_inc_b_1_inc_c_1_n_100", .kind = .run, .n = 100, .inc_b = 1, .inc_c = 1, .k = 1, .lda = 2, .ldm = 100, .uplo = 'L', .alpha = .{ 1.0, 1.0 }, .beta = .{ 1.0, 1.0 } },
    .{ .name = "lower_k_2_inc_b_1_inc_c_1_n_100", .kind = .run, .n = 100, .inc_b = 1, .inc_c = 1, .k = 2, .lda = 3, .ldm = 100, .uplo = 'L', .alpha = .{ 1.0, 1.0 }, .beta = .{ 1.0, 1.0 } },
    .{ .name = "lower_k_2_inc_b_2_inc_c_1_n_100", .kind = .run, .n = 100, .inc_b = 2, .inc_c = 1, .k = 2, .lda = 3, .ldm = 100, .uplo = 'L', .alpha = .{ 2.0, 1.0 }, .beta = .{ 2.0, 1.0 } },
    .{ .name = "lower_k_2_inc_b_2_inc_c_2_n_100", .kind = .run, .n = 100, .inc_b = 2, .inc_c = 2, .k = 2, .lda = 3, .ldm = 100, .uplo = 'L', .alpha = .{ 2.0, 1.0 }, .beta = .{ 2.0, 1.0 } },
    .{ .name = "check_for_NaN", .kind = .nan, .n = 100, .inc_b = 1, .inc_c = 1, .k = 0, .lda = 1, .ldm = 100, .uplo = 'U', .alpha = .{ 1.0, 1.0 }, .beta = .{ 1.0, 1.0 } },
    .{ .name = "xerbla_uplo_invalid", .kind = .bad, .n = 1, .inc_b = 1, .inc_c = 1, .k = 0, .lda = 1, .uplo = 'O', .info = 1 },
    .{ .name = "xerbla_n_invalid", .kind = .bad, .n = -1, .inc_b = 1, .inc_c = 1, .k = 0, .lda = 1, .uplo = 'U', .info = 2 },
    .{ .name = "check_n_zero", .kind = .run, .n = 0, .inc_b = 1, .inc_c = 1, .k = 0, .lda = 1, .ldm = 1, .uplo = 'U', .alpha = .{ 1.0, 1.0 }, .beta = .{ 0.0, 0.0 } },
    .{ .name = "xerbla_inc_b_zero", .kind = .bad, .n = 1, .inc_b = 0, .inc_c = 1, .k = 0, .lda = 1, .uplo = 'U', .info = 8 },
    .{ .name = "xerbla_inc_c_zero", .kind = .bad, .n = 1, .inc_b = 1, .inc_c = 0, .k = 0, .lda = 1, .uplo = 'U', .info = 11 },
    .{ .name = "xerbla_k_invalid", .kind = .bad, .n = 1, .inc_b = 1, .inc_c = 1, .k = -1, .lda = 1, .uplo = 'U', .info = 3 },
    .{ .name = "xerbla_lda_invalid", .kind = .bad, .n = 1, .inc_b = 1, .inc_c = 1, .k = 0, .lda = -1, .uplo = 'U', .info = 6 },
};

/// `test_zspmv.c (test_cspmv.c)`.
const spmv_cases = [_]Case{
    .{ .name = "upper_inc_b_1_inc_c_1_N_100", .kind = .run, .N = 100, .inc_b = 1, .inc_c = 1, .uplo = 'U', .alpha = .{ 1.0, 1.0 }, .beta = .{ 0.0, 0.0 } },
    .{ .name = "upper_inc_b_1_inc_c_2_N_100", .kind = .run, .N = 100, .inc_b = 1, .inc_c = 2, .uplo = 'U', .alpha = .{ 1.0, 1.0 }, .beta = .{ 0.0, 0.0 } },
    .{ .name = "upper_inc_b_2_inc_c_1_N_100", .kind = .run, .N = 100, .inc_b = 2, .inc_c = 1, .uplo = 'U', .alpha = .{ 1.0, 0.0 }, .beta = .{ 1.0, 0.0 } },
    .{ .name = "upper_inc_b_2_inc_c_2_N_100", .kind = .run, .N = 100, .inc_b = 2, .inc_c = 2, .uplo = 'U', .alpha = .{ 2.5, -2.1 }, .beta = .{ 0.0, 1.0 } },
    .{ .name = "lower_inc_b_1_inc_c_1_N_100", .kind = .run, .N = 100, .inc_b = 1, .inc_c = 1, .uplo = 'L', .alpha = .{ 1.0, 1.0 }, .beta = .{ 0.0, 0.0 } },
    .{ .name = "lower_inc_b_1_inc_c_2_N_100", .kind = .run, .N = 100, .inc_b = 1, .inc_c = 2, .uplo = 'L', .alpha = .{ 1.0, 1.0 }, .beta = .{ 0.0, 0.0 } },
    .{ .name = "lower_inc_b_2_inc_c_1_N_100", .kind = .run, .N = 100, .inc_b = 2, .inc_c = 1, .uplo = 'L', .alpha = .{ 1.0, 0.0 }, .beta = .{ 1.0, 0.0 } },
    .{ .name = "lower_inc_b_2_inc_c_2_N_100", .kind = .run, .N = 100, .inc_b = 2, .inc_c = 2, .uplo = 'L', .alpha = .{ 2.5, -2.1 }, .beta = .{ 0.0, 1.0 } },
    .{ .name = "check_for_NaN", .kind = .nan, .N = 100, .inc_b = 1, .inc_c = 1, .uplo = 'U', .alpha = .{ 1.0, 1.0 }, .beta = .{ 0.0, 0.0 } },
    .{ .name = "xerbla_uplo_invalid", .kind = .bad, .N = 100, .inc_b = 1, .inc_c = 1, .uplo = 'O', .info = 1 },
    .{ .name = "xerbla_N_invalid", .kind = .bad, .N = -1, .inc_b = 1, .inc_c = 1, .uplo = 'U', .info = 2 },
    .{ .name = "xerbla_inc_b_zero", .kind = .bad, .N = 100, .inc_b = 0, .inc_c = 1, .uplo = 'U', .info = 6 },
    .{ .name = "xerbla_inc_c_zero", .kind = .bad, .N = 100, .inc_b = 1, .inc_c = 0, .uplo = 'U', .info = 9 },
};

/// `test_ztrmv.c (test_ctrmv.c)`.
const trmv_cases = [_]Case{
    .{ .name = "conj_notrans_upper_not_unit_triangular", .kind = .run, .n = 300, .incx = 1, .lda = 300, .uplo = 'U', .diag = 'N', .trans = 'R' },
    .{ .name = "conj_notrans_upper_unit_triangular", .kind = .run, .n = 300, .incx = 1, .lda = 300, .uplo = 'U', .diag = 'U', .trans = 'R' },
    .{ .name = "conj_notrans_lower_not_triangular", .kind = .run, .n = 300, .incx = 1, .lda = 300, .uplo = 'L', .diag = 'N', .trans = 'R' },
    .{ .name = "conj_notrans_lower_unit_triangular", .kind = .run, .n = 300, .incx = 1, .lda = 300, .uplo = 'L', .diag = 'U', .trans = 'R' },
    .{ .name = "conj_notrans_upper_not_unit_triangular_incx_2", .kind = .run, .n = 300, .incx = 2, .lda = 300, .uplo = 'U', .diag = 'N', .trans = 'R' },
    .{ .name = "conj_notrans_upper_unit_triangular_incx_2", .kind = .run, .n = 300, .incx = 2, .lda = 300, .uplo = 'U', .diag = 'U', .trans = 'R' },
    .{ .name = "conj_notrans_lower_not_triangular_incx_2", .kind = .run, .n = 300, .incx = 2, .lda = 300, .uplo = 'L', .diag = 'N', .trans = 'R' },
    .{ .name = "conj_notrans_lower_unit_triangular_incx_2", .kind = .run, .n = 300, .incx = 2, .lda = 300, .uplo = 'L', .diag = 'U', .trans = 'R' },
};

/// `test_ztrsv.c (test_ctrsv.c)`.
const trsv_cases = [_]Case{
    .{ .name = "conj_notrans_upper_not_unit_triangular", .kind = .run, .n = 300, .incx = 1, .lda = 300, .uplo = 'U', .diag = 'N', .trans = 'R' },
    .{ .name = "conj_notrans_upper_unit_triangular", .kind = .run, .n = 300, .incx = 1, .lda = 300, .uplo = 'U', .diag = 'U', .trans = 'R' },
    .{ .name = "conj_notrans_lower_not_triangular", .kind = .run, .n = 300, .incx = 1, .lda = 300, .uplo = 'L', .diag = 'N', .trans = 'R' },
    .{ .name = "conj_notrans_lower_unit_triangular", .kind = .run, .n = 300, .incx = 1, .lda = 300, .uplo = 'L', .diag = 'U', .trans = 'R' },
    .{ .name = "conj_notrans_upper_not_unit_triangular_incx_2", .kind = .run, .n = 300, .incx = 2, .lda = 300, .uplo = 'U', .diag = 'N', .trans = 'R' },
    .{ .name = "conj_notrans_upper_unit_triangular_incx_2", .kind = .run, .n = 300, .incx = 2, .lda = 300, .uplo = 'U', .diag = 'U', .trans = 'R' },
    .{ .name = "conj_notrans_lower_not_triangular_incx_2", .kind = .run, .n = 300, .incx = 2, .lda = 300, .uplo = 'L', .diag = 'N', .trans = 'R' },
    .{ .name = "conj_notrans_lower_unit_triangular_incx_2", .kind = .run, .n = 300, .incx = 2, .lda = 300, .uplo = 'L', .diag = 'U', .trans = 'R' },
};
