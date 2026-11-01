//! OpenBLAS's extension unit tests, `utest/test_extensions/test_*.c`,
//! OpenBLAS v0.3.34-317: the `?amin`, `?amax`, `?sum` and `i?amin` vectors,
//! `axpby`, `axpyc`, the complex `rot`, `rotg` and `scal`, `gbmv` and
//! `gemv` with every `trans` letter, the complex symmetric `sbmv` and
//! `spmv`, `trmv` and `trsv` with `R`, `rotmg`, the column-major `geadd`,
//! `omatcopy`, `imatcopy` and `gemmt`, and complex `gemm` with `R`. The
//! tables at the end are the sources' cases (the `c` and `s` files' are the
//! `z` and `d` ones', but for the matrix extensions, whose files differ, so
//! each has its own); the `c_api_` cases are the same calls through CBLAS,
//! so not repeated, but for the column-major ones whose arguments differ.
//! Random operands come from a seeded generator, not `rand()`.
//!
//! `test_bgemm.c` is here too. Not translatable: `gemv`'s bad
//! letter `Z`, `gemmt`'s bad `UPLO` `O` (only checked to be no `Uplo`), and
//! the matrix extensions' bad `ORDER`.
//!
//! A norm of differences counts elements that agree exactly (infinity with
//! infinity, NaN with NaN) as equal, so overflow in both is not NaN; the C
//! passes any NaN norm. `trmv` and `trsv` compare relative to the result,
//! and single `trsv` with single precision's tolerance (see `trmvTrsv`).

const std = @import("std");
const ob = @import("__root__.zig");
const blas = @import("../../__root__.zig");
const apple_accel = @import("../../apple_accel.zig");

const Float = ob.Float;
const C32 = std.math.Complex(f32);
const C64 = std.math.Complex(f64);

const single_eps = 1e-4;
const double_eps = 1e-13;
const single_tol = 1e-2;
const double_tol = 1e-10;
const single_eps_zgemv = 5e-3;
const double_eps_zgemv = 1e-11;

fn tol(comptime T: type, single: f64, double: f64) f64 {
    return if (Float(T) == f32) single else double;
}

// ---- the vector cases ----------------------------------------------------------

/// `x` as raw floats (two per complex element), `n` elements at stride
/// `inc`, and the expected result.
const VecCase = struct { name: []const u8, n: isize, inc: isize, want: f64, x: []const f64 };

const Reduce = enum { amin, amax, sum, iamin };

fn vecCases(comptime T: type, comptime what: Reduce, cases: []const VecCase) !void {
    const F = Float(T);
    var buf: [512]F align(@alignOf(T)) = undefined;
    for (cases) |c| {
        for (buf[0..c.x.len], c.x) |*b, v| b.* = @floatCast(v);
        const x: [*]T = @ptrCast(&buf);
        const got: f64 = switch (what) {
            .amin => try blas.extreme(ob.io, T, .amin, c.n, x, c.inc),
            .amax => try blas.extreme(ob.io, T, .amax, c.n, x, c.inc),
            .sum => try blas.sum(ob.io, T, c.n, x, c.inc),
            .iamin => @floatFromInt(try blas.index(ob.io, T, .amin, c.n, x, c.inc)),
        };
        const t: f64 = if (what == .iamin) 0 else tol(T, single_eps, double_eps);
        if (@abs(got - c.want) <= t) continue;
        ob.say("{s} {s}: got {d}, want {d}\n", .{ @typeName(T), c.name, got, c.want });
        return error.TestUnexpectedResult;
    }
}

/// The vector cases, for `check.eachCore`; with `level2All` and `level3All`,
/// every case here. Apart, so each fits a rig image.
pub fn all() !void {
    try reductions();
    try axpby();
    try axpyc();
    try complexRot();
    try complexScal();
    try complexRotg();
    try rotmg();
}

/// The level-2 cases, for `check.eachCore`.
pub fn level2All() !void {
    try gbmv();
    try gemvOS();
    try gemvTCUD();
    try sbmvSpmv();
    try trmvTrsv();
}

/// The matrix extensions and `gemm`, for `check.eachCore`; they run with the
/// level-3 suite, whose drivers they reach.
pub fn level3All() !void {
    try geadd();
    try matcopy();
    try complexGemm();
    try gemmt();
    try bgemm();
}

fn reductions() !void {
    inline for (.{ f32, f64 }) |F| {
        try vecCases(F, .amin, &amin_real);
        try vecCases(F, .sum, &sum_real);
        try vecCases(F, .iamin, &iamin_real);
    }
    inline for (.{ C32, C64 }) |C| {
        try vecCases(C, .amax, &amax_complex);
        try vecCases(C, .amin, &amin_complex);
        try vecCases(C, .sum, &sum_complex);
        try vecCases(C, .iamin, &iamin_complex);
    }
}

const isComplex = ob.isComplex;
const testing = std.testing;

test "utest_ext: amin, amax, sum, iamin" {
    try reductions();
}

// ---- the parameter cases ---------------------------------------------------------

/// One case: the source's variables by their own names.
const Case = struct {
    name: []const u8,
    kind: enum { run, nan, bad } = .run,
    n: isize = 0,
    N: isize = 0,
    m: isize = 0,
    k: isize = 0,
    kl: isize = 0,
    ku: isize = 0,
    lda: isize = 0,
    ldb: isize = 0,
    ldc: isize = 0,
    ldm: isize = 0,
    incx: isize = 1,
    incy: isize = 1,
    inc: isize = 1,
    inc_x: isize = 1,
    inc_y: isize = 1,
    inc_b: isize = 1,
    inc_c: isize = 1,
    trans: u8 = 'N',
    transa: u8 = 'N',
    transb: u8 = 'N',
    uplo: u8 = 'U',
    diag: u8 = 'N',
    order: u8 = 'C',
    alpha: [2]f64 = .{ 0, 0 },
    beta: [2]f64 = .{ 0, 0 },
    c: [2]f64 = .{ 0, 0 },
    s: [2]f64 = .{ 0, 0 },
    info: u8 = 0,
};

fn val(comptime T: type, v: [2]f64) T {
    return ob.val(T, v[0], v[1]);
}

/// `drand_generate`: values in [0, 1).
fn fill(comptime T: type, p: []T, r: std.Random) void {
    for (p) |*v| v.* = if (comptime isComplex(T)) .init(@floatCast(r.float(f64)), @floatCast(r.float(f64))) else @floatCast(r.float(f64));
}

fn span(n: isize, inc: isize) usize {
    return @intCast(@max(n, 1) * @max(@as(isize, @intCast(@abs(inc))), 1));
}

/// The norm of `got - want` over `n` elements at stride `|inc|`; elements
/// that agree exactly count as zero.
fn diffNorm(comptime T: type, n: isize, got: []const T, want: []const T, inc: isize) f64 {
    if (n <= 0) return 0;
    var s: f64 = 0;
    const st: usize = @abs(inc);
    for (0..@intCast(n)) |i| {
        const g = got[i * st];
        const w = want[i * st];
        if (ob.eql(g, w)) continue;
        const d = ob.sub(g, w);
        s += if (comptime isComplex(T)) @as(f64, d.re) * d.re + @as(f64, d.im) * d.im else @as(f64, d) * d;
    }
    return @sqrt(s);
}

/// The norm of `v`, `n` elements at stride `|inc|`.
fn norm2(comptime T: type, n: isize, v: []const T, inc: isize) f64 {
    var s: f64 = 0;
    for (0..@intCast(@max(n, 0))) |i| s += @as(f64, ob.absv(v[i * @abs(inc)])) * ob.absv(v[i * @abs(inc)]);
    return @sqrt(s);
}

fn expectNorm(name: []const u8, comptime T: type, norm: f64, t: f64) !void {
    if (norm <= t) return;
    ob.say("{s} {s}: norm of differences {e} over {e}\n", .{ @typeName(T), name, norm, t });
    return error.TestUnexpectedResult;
}

/// `T`'s elements of `p`, conjugated: the source's `?conjugate_vector`.
fn conjugate(comptime T: type, p: []T) void {
    if (comptime isComplex(T)) for (p) |*v| {
        v.* = v.conjugate();
    };
}

fn axpby() !void {
    var prng = std.Random.DefaultPrng.init(0xA8);
    const al = testing.allocator;
    inline for (.{ f32, f64, C32, C64 }) |T| {
        const cases: []const Case = if (comptime isComplex(T)) &axpby_complex else &axpby_real;
        for (cases) |c| {
            const x = try al.alloc(T, span(c.n, c.incx));
            defer al.free(x);
            const y = try al.alloc(T, span(c.n, c.incy));
            defer al.free(y);
            fill(T, x, prng.random());
            fill(T, y, prng.random());
            const yv = try al.dupe(T, y);
            defer al.free(yv);
            // scal takes only positive strides
            try blas.scal(ob.io, T, c.n, val(T, c.beta), yv.ptr, @intCast(@abs(c.incy)));
            try blas.axpy(ob.io, T, c.n, val(T, c.alpha), x.ptr, c.incx, yv.ptr, c.incy);
            try blas.axpby(ob.io, T, c.n, val(T, c.alpha), x.ptr, c.incx, val(T, c.beta), y.ptr, c.incy);
            try expectNorm(c.name, T, diffNorm(T, c.n, y, yv, c.incy), tol(T, single_eps, double_eps));
        }
    }
}

test "utest_ext: axpby against scal and axpy" {
    try axpby();
}

fn axpyc() !void {
    var prng = std.Random.DefaultPrng.init(0xA9);
    const al = testing.allocator;
    inline for (.{ C32, C64 }) |T| for (axpyc_cases) |c| {
        const x = try al.alloc(T, span(c.n, c.incx));
        defer al.free(x);
        const y = try al.alloc(T, span(c.n, c.incy));
        defer al.free(y);
        fill(T, x, prng.random());
        fill(T, y, prng.random());
        const xv = try al.dupe(T, x);
        defer al.free(xv);
        const yv = try al.dupe(T, y);
        defer al.free(yv);
        conjugate(T, xv);
        try blas.axpy(ob.io, T, c.n, val(T, c.alpha), xv.ptr, c.incx, yv.ptr, c.incy);
        try blas.axpyc(ob.io, T, c.n, val(T, c.alpha), x.ptr, c.incx, y.ptr, c.incy);
        try expectNorm(c.name, T, diffNorm(T, c.n, y, yv, c.incy), tol(T, single_eps, double_eps));
    };
}

test "utest_ext: axpyc against conjugation and axpy" {
    try axpyc();
}

fn complexRot() !void {
    var prng = std.Random.DefaultPrng.init(0xAA);
    const al = testing.allocator;
    inline for (.{ C32, C64 }) |T| for (rot_cases) |c| {
        const x = try al.alloc(T, span(c.n, c.inc_x));
        defer al.free(x);
        const y = try al.alloc(T, span(c.n, c.inc_y));
        defer al.free(y);
        fill(T, x, prng.random());
        fill(T, y, prng.random());
        const xv = try al.dupe(T, x);
        defer al.free(xv);
        const yv = try al.dupe(T, y);
        defer al.free(yv);
        // c*x + s*y, then -conj(s)*x + c*y
        try blas.axpby(ob.io, T, c.n, val(T, c.s), y.ptr, c.inc_y, val(T, c.c), xv.ptr, c.inc_x);
        try blas.axpby(ob.io, T, c.n, val(T, .{ -c.s[0], c.s[1] }), x.ptr, c.inc_x, val(T, c.c), yv.ptr, c.inc_y);
        try blas.rot(ob.io, T, c.n, x.ptr, c.inc_x, y.ptr, c.inc_y, @floatCast(c.c[0]), @floatCast(c.s[0]));
        const norm = (diffNorm(T, c.n, x, xv, c.inc_x) + diffNorm(T, c.n, y, yv, c.inc_y)) / 2;
        try expectNorm(c.name, T, norm, tol(T, single_eps, double_eps));
    };
}

test "utest_ext: complex rot against axpby" {
    try complexRot();
}

fn complexScal() !void {
    var prng = std.Random.DefaultPrng.init(0xAB);
    const al = testing.allocator;
    inline for (.{ C32, C64 }) |T| for (scal_cases) |c| {
        const x = try al.alloc(T, span(c.N, c.inc));
        defer al.free(x);
        fill(T, x, prng.random());
        const xv = try al.dupe(T, x);
        defer al.free(xv);
        const a = val(T, c.alpha);
        for (0..@intCast(c.N)) |i| xv[i * @as(usize, @intCast(c.inc))] = a.mul(xv[i * @as(usize, @intCast(c.inc))]);
        try blas.scal(ob.io, T, c.N, a, x.ptr, c.inc);
        try expectNorm(c.name, T, diffNorm(T, c.N, x, xv, c.inc), tol(T, single_eps, double_eps));
    };
}

test "utest_ext: complex scal against the plain product" {
    try complexScal();
}

fn complexRotg() !void {
    const Rotg = struct { name: []const u8, a: [2]f64, b: [2]f64, c: f64, s: [2]f64, r: [2]f64 };
    // the double cases' expected values; the single ones' are the same to 5 places
    const cases = [_]Rotg{
        .{ .name = "zero_a", .a = .{ 0, 0 }, .b = .{ 1, 1 }, .c = 0.0, .s = .{ 0.70710678118655, -0.70710678118655 }, .r = .{ 1.41421356237310, 0.0 } },
        .{ .name = "zero_b", .a = .{ 1, 1 }, .b = .{ 0, 0 }, .c = 1.0, .s = .{ 0.0, 0.0 }, .r = .{ 1.0, 1.0 } },
        .{ .name = "zero_real", .a = .{ 0, 1 }, .b = .{ 0, 1 }, .c = 0.70710678118654, .s = .{ 0.70710678118654, 0.0 }, .r = .{ 0.0, 1.41421356237309 } },
        .{ .name = "positive_real_positive_img", .a = .{ 3, 4 }, .b = .{ 4, 6 }, .c = 0.56980288229818, .s = .{ 0.82051615050939, -0.04558423058385 }, .r = .{ 5.26497863243527, 7.01997150991369 } },
        .{ .name = "negative_real_positive_img", .a = .{ -3, 4 }, .b = .{ -4, 6 }, .c = 0.56980288229818, .s = .{ 0.82051615050939, 0.04558423058385 }, .r = .{ -5.26497863243527, 7.01997150991369 } },
        .{ .name = "positive_real_negative_img", .a = .{ 3, -4 }, .b = .{ 4, -6 }, .c = 0.56980288229818, .s = .{ 0.82051615050939, 0.04558423058385 }, .r = .{ 5.26497863243527, -7.01997150991369 } },
        .{ .name = "negative_real_negative_img", .a = .{ -3, -4 }, .b = .{ -4, -6 }, .c = 0.56980288229818, .s = .{ 0.82051615050939, -0.04558423058385 }, .r = .{ -5.26497863243527, -7.01997150991369 } },
    };
    inline for (.{ C32, C64 }) |T| for (cases) |c| {
        var a = val(T, c.a);
        var cc: Float(T) = undefined;
        var s: T = undefined;
        blas.zrotg(T, &a, val(T, c.b), &cc, &s);
        const t = tol(T, single_eps, double_eps);
        const got = [_]f64{ cc, s.re, s.im, a.re, a.im };
        const want = [_]f64{ c.c, c.s[0], c.s[1], c.r[0], c.r[1] };
        for (got, want) |g, w| if (!(@abs(g - w) <= t)) {
            ob.say("{s} {s}: got {d}, want {d}\n", .{ @typeName(T), c.name, g, w });
            return error.TestUnexpectedResult;
        };
    };
}

test "utest_ext: complex rotg" {
    try complexRotg();
}

fn rotmg() !void {
    const Rotmg = struct { name: []const u8, in: [4]f64, out: [4]f64, param: [5]f64 };
    const cases = [_]Rotmg{
        .{ .name = "y1_zero", .in = .{ 2, 2, 8, 0 }, .out = .{ 2, 2, 8, 0 }, .param = .{ -2, 0, 0, 0, 0 } },
        .{ .name = "d1_negative", .in = .{ -1, 2, 8, 8 }, .out = .{ 0, 0, 0, 8 }, .param = .{ -1, 0, 0, 0, 0 } },
        .{ .name = "d1_positive_d2_positive_x1_zero", .in = .{ 2, 2, 0, 8 }, .out = .{ 2, 2, 8, 8 }, .param = .{ 1, 0, 0, 0, 0 } },
        .{ .name = "scaled_y_greater_than_scaled_x", .in = .{ 1, -2, 8, 8 }, .out = .{ 0, 0, 0, 8 }, .param = .{ -1, 0, 0, 0, 0 } },
    };
    inline for (.{ f32, f64 }) |F| for (cases) |c| {
        var d1: F = @floatCast(c.in[0]);
        var d2: F = @floatCast(c.in[1]);
        var x1: F = @floatCast(c.in[2]);
        const y1: F = @floatCast(c.in[3]);
        var param: [5]F = @splat(0);
        blas.rotmg(F, &d1, &d2, &x1, y1, &param);
        const t = tol(F, single_eps, double_eps);
        const got = [_]f64{ d1, d2, x1, y1, param[0], param[1], param[2], param[3], param[4] };
        const want = c.out ++ c.param;
        for (got, want) |g, w| if (!(@abs(g - w) <= t)) {
            ob.say("{s} {s}: got {d}, want {d}\n", .{ @typeName(F), c.name, g, w });
            return error.TestUnexpectedResult;
        };
    };
}

test "utest_ext: rotmg" {
    try rotmg();
}

// ---- level 2 -------------------------------------------------------------------

fn transOf(c: u8) blas.Trans {
    return ob.letter(blas.Trans, c) orelse .d;
}

fn transposed(c: u8) bool {
    return c == 'T' or c == 'C' or c == 'U' or c == 'D';
}

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

fn gbmv() !void {
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
        try blas.gemv(al, ob.io, T, transOf(c.trans), c.m, c.n, val(T, c.alpha), full.ptr, c.m, b.ptr, c.inc_b, val(T, c.beta), yv.ptr, c.inc_c);
        try blas.gbmv(al, ob.io, T, transOf(c.trans), c.m, c.n, c.kl, c.ku, val(T, c.alpha), band.ptr, c.lda, b.ptr, c.inc_b, val(T, c.beta), y.ptr, c.inc_c);
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
    try blas.gemv(al, ob.io, T, transOf(c.trans), c.m, c.n, val(T, c.alpha), a.ptr, c.lda, x.ptr, ib, val(T, c.beta), y.ptr, ic);
    try expectNorm(c.name, T, diffNorm(T, leny, y, yv, ic), t);
}

fn gemvOS() !void {
    var prng = std.Random.DefaultPrng.init(0xB1);
    inline for (.{ C32, C64 }) |T| for (gemv_n_cases) |c| try gemvCase(T, c, prng.random(), tol(T, single_tol, double_tol));
}

test "utest_ext: gemv with O and S against the plain product" {
    try gemvOS();
}

fn gemvTCUD() !void {
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
        try ob.expectError(want, blas.gemv(al, ob.io, T, transOf(c.trans), c.m, c.n, val(T, .{ 1, 1 }), &a, c.lda, &x, c.inc_x, val(T, .{ 1, 1 }), &y, c.inc_y), c.name);
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
        blas.spmv(al, ob.io, T, .u, n, one, &a, &x, c.inc_b, ob.zero(T), &y, c.inc_c)
    else
        blas.sbmv(al, ob.io, T, .u, n, c.k, one, &a, c.lda, &x, c.inc_b, ob.zero(T), &y, c.inc_c);
    try ob.expectError(want, res, c.name);
}

fn sbmvSpmv() !void {
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
            try blas.gemv(al, ob.io, T, .n, in, in, val(T, c.alpha), full.ptr, @max(in, 1), b.ptr, c.inc_b, val(T, c.beta), yv.ptr, c.inc_c);
            if (pack) {
                try blas.spmv(al, ob.io, T, if (up) .u else .l, in, val(T, c.alpha), store.ptr, b.ptr, c.inc_b, val(T, c.beta), y.ptr, c.inc_c);
            } else {
                try blas.sbmv(al, ob.io, T, if (up) .u else .l, in, c.k, val(T, c.alpha), store.ptr, c.lda, b.ptr, c.inc_b, val(T, c.beta), y.ptr, c.inc_c);
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

fn trmvTrsv() !void {
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
            const f = if (solve) blas.trsv else blas.trmv;
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

// ---- the matrix extensions -----------------------------------------------------

/// The sources' `?matrix_difference`: the mean over `cols` columns of `ld`
/// of the norm of `got - want` down `rows` of each; elements that agree
/// exactly count as zero.
fn matDiff(comptime T: type, rows: usize, cols: usize, got: []const T, want: []const T, ld: usize) f64 {
    if (rows == 0 or cols == 0) return 0;
    var sum: f64 = 0;
    for (0..cols) |j| sum += diffNorm(T, @intCast(rows), got[j * ld ..], want[j * ld ..], 1);
    return sum / @as(f64, @floatFromInt(cols));
}

fn transposedLetter(c: u8) bool {
    return c == 'T' or c == 'C';
}

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
        return ob.expectError(want, blas.geadd(ob.io, T, c.m, c.n, val(T, .{ 1, 1 }), &a, c.lda, val(T, .{ 1, 1 }), &x, c.ldc, transOf(c.transa), transOf(c.transb)), c.name);
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
    try blas.geadd(ob.io, T, c.m, c.n, alpha, a.ptr, c.lda, beta, cm.ptr, c.ldc, transOf(c.transa), transOf(c.transb));
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
        try blas.geadd(ob.io, T, 2, 2, val(T, .{ 1, 0 }), &a, 2, ob.zero(T), &c, 2, .t, p[0]);
        for (c, p[1]) |got, from| try ob.expect(ob.eql(got, a[from]), "c_api_matrix_2x2_transA");
    }
}

fn geadd() !void {
    var prng = std.Random.DefaultPrng.init(0xC0);
    inline for (.{ .{ f32, &geadd_s }, .{ f64, &geadd_d }, .{ C32, &geadd_c }, .{ C64, &geadd_z } }) |p| {
        for (p[1]) |c| try geaddCase(p[0], c, prng.random());
        try geaddFixed(p[0]);
    }
}

test "utest_ext: geadd against the plain sum, and its bad arguments" {
    try geadd();
}

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
        const res = if (in_place) blas.imatcopy(al, ob.io, T, tr, c.m, c.n, one, &a, c.lda, c.ldb) else blas.omatcopy(ob.io, T, tr, c.m, c.n, one, &a, c.lda, &b, c.ldb);
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
        try blas.imatcopy(al, ob.io, T, tr, c.m, c.n, alpha, a.ptr, c.lda, c.ldb);
        @memcpy(b, a[0..b.len]);
    } else try blas.omatcopy(ob.io, T, tr, c.m, c.n, alpha, a.ptr, c.lda, b.ptr, c.ldb);
    try expectNorm(c.name, T, matDiff(T, out_rows, out_cols, b, want, ldb), tol(T, single_eps, double_eps));
}

fn matcopy() !void {
    var prng = std.Random.DefaultPrng.init(0xC1);
    inline for (.{ .{ f32, &omatcopy_s, &imatcopy_s }, .{ f64, &omatcopy_d, &imatcopy_d }, .{ C32, &omatcopy_c, &imatcopy_c }, .{ C64, &omatcopy_z, &imatcopy_z } }) |p| {
        for (p[1]) |c| try matcopyCase(p[0], false, c, prng.random());
        for (p[2]) |c| try matcopyCase(p[0], true, c, prng.random());
    }
}

test "utest_ext: omatcopy and imatcopy against the definition, and their bad arguments" {
    try matcopy();
}

/// `test_[cz]gemm.c`: `gemm` with a conjugated operand (`R`, `C`) against
/// `gemm` on a copy that `imatcopy` conjugated (`R` becoming `N`). Every
/// case is 100 square; on a rig, 23.
fn complexGemm() !void {
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
        if (ob.reduced()) {
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
            try blas.imatcopy(al, ob.io, T, .r, rows_a, cols_a, one, av.ptr, c.lda, c.lda);
            tav = .n;
        }
        if (c.transb == 'R') {
            try blas.imatcopy(al, ob.io, T, .r, rows_b, cols_b, one, bv.ptr, c.ldb, c.ldb);
            tbv = .n;
        }
        const alpha = val(T, c.alpha);
        const beta = val(T, c.beta);
        try blas.gemm(al, ob.io, T, tav, tbv, c.m, c.n, c.k, alpha, av.ptr, c.lda, bv.ptr, c.ldb, beta, cv.ptr, c.ldc);
        try blas.gemm(al, ob.io, T, transOf(c.transa), transOf(c.transb), c.m, c.n, c.k, alpha, a.ptr, c.lda, b.ptr, c.ldb, beta, cm.ptr, c.ldc);
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
        return ob.expectError(want, blas.gemmt(al, ob.io, T, uplo.?, transOf(c.transa), transOf(c.transb), c.m, c.k, one, &x, c.lda, &x, c.ldb, one, &x, c.ldc), c.name);
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
    try blas.gemm(al, ob.io, T, transOf(c.transa), transOf(c.transb), c.m, c.m, c.k, alpha, a.ptr, c.lda, b.ptr, c.ldb, beta, full.ptr, c.ldc);
    for (0..m) |j| for (0..m) |i| {
        if (if (c.uplo == 'U') i <= j else i >= j) want[i + j * ldc] = full[i + j * ldc];
    };
    try blas.gemmt(al, ob.io, T, uplo.?, transOf(c.transa), transOf(c.transb), c.m, c.k, alpha, a.ptr, c.lda, b.ptr, c.ldb, beta, cm.ptr, c.ldc);
    const floats: f64 = @floatFromInt(m * ldc * @as(usize, if (comptime isComplex(T)) 2 else 1));
    try expectNorm(c.name, T, diffNorm(T, @intCast(m * ldc), cm, want, 1) / floats, tol(T, single_eps, double_eps));
}

fn gemmt() !void {
    var prng = std.Random.DefaultPrng.init(0xC3);
    inline for (.{ .{ f32, &gemmt_s }, .{ f64, &gemmt_d }, .{ C32, &gemmt_c }, .{ C64, &gemmt_z } }) |p| {
        for (p[1]) |c| try gemmtCase(p[0], c, prng.random());
    }
}

test "utest_ext: gemmt against gemm's triangle, and its bad arguments" {
    try gemmt();
}

// ---- bgemm ---------------------------------------------------------------------

/// One `CTEST` of `test_bgemm.c`'s `check_bgemm`.
const BgemmCase = struct { name: []const u8, ta: u8, tb: u8, m: isize, n: isize, k: isize, alpha: f32, beta: f32, lda: isize, ldb: isize, ldc: isize };

/// `test_bgemm.c`'s `float_to_bf16`: the high half of the bits, truncated.
fn truncBf16(v: f32) blas.types.BF16 {
    return @fromBackingInt(@as(u16, @intCast(@as(u32, @bitCast(v)) >> 16)));
}

/// `test_bgemm.c`'s `bf16_to_float`: the bits widened, nothing flushed.
fn bitsF32(v: blas.types.BF16) f32 {
    return @bitCast(@as(u32, @backingInt(v)) << 16);
}

/// `check_bgemm`: `bgemm` against `sbgemm` on the same bfloat16 operands and
/// scalars, the largest `|bgemm - sbgemm| / (|sbgemm| + 1)` within `BGEMM_TOL`.
fn bgemmCase(c: BgemmCase, r: std.Random) !void {
    const al = testing.allocator;
    const BF16 = blas.types.BF16;
    const ta = transposedLetter(c.ta);
    const a_cols: usize = @intCast(if (ta) c.m else c.k);
    const b_cols: usize = @intCast(if (transposedLetter(c.tb)) c.k else c.n);
    const n: usize = @intCast(c.n);
    const m: usize = @intCast(c.m);
    const ldc: usize = @intCast(c.ldc);
    const a = try al.alloc(BF16, a_cols * @as(usize, @intCast(c.lda)));
    defer al.free(a);
    const b = try al.alloc(BF16, b_cols * @as(usize, @intCast(c.ldb)));
    defer al.free(b);
    const cb = try al.alloc(BF16, n * ldc);
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
    try blas.gemmHalf(al, ob.io, .sbgemm, transOf(c.ta), transOf(c.tb), c.m, c.n, c.k, bitsF32(alpha), a.ptr, c.lda, b.ptr, c.ldb, bitsF32(beta), cs.ptr, c.ldc);
    try blas.gemmHalf(al, ob.io, .bgemm, transOf(c.ta), transOf(c.tb), c.m, c.n, c.k, alpha, a.ptr, c.lda, b.ptr, c.ldb, beta, cb.ptr, c.ldc);
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

fn bgemm() !void {
    var prng = std.Random.DefaultPrng.init(0xB6);
    for (bgemm_cases) |c| try bgemmCase(c, prng.random());
}

test "utest_ext: bgemm against sbgemm on the same bfloat16 operands" {
    try bgemm();
}

// ---- data ----------------------------------------------------------------------

/// `test_damin.c (test_samin.c is the same)`.
const amin_real = [_]VecCase{
    .{ .name = "bad_args_N_0", .n = 0, .inc = 1, .want = 0, .x = &.{ 1000, 999, 998, 997, 996, 995, 994, 993, 992, 991, 990, 989, 988, 987, 986, 985, 984, 983, 982, 981, 980, 979, 978, 977, 976, 975, 974, 973, 972, 971, 970, 969, 968, 967, 966, 965, 964, 963, 962, 961, 960, 959, 958, 957, 956, 955, 954, 953, 952, 951, 950, 949, 948, 947, 946, 945, 944, 943, 942, 941, 940, 939, 938, 937, 936, 935, 934, 933, 932, 931 } },
    .{ .name = "step_zero", .n = 70, .inc = 0, .want = 1000, .x = &.{ -1000, -999, -998, -997, -996, -995, -994, -993, 0, -991, -990, -989, -988, -987, -986, -985, -984, -983, -982, -981, -980, -979, -978, -977, -976, -975, -974, -973, -972, -971, -970, -969, -968, -967, -966, -965, -964, -963, -962, -961, -960, -959, -958, -957, -956, -955, -954, -953, -952, -951, -950, -949, -948, -947, -946, -945, -944, -943, -942, -941, -940, -939, -938, -937, -936, -935, -934, -933, -932, -931 } },
    .{ .name = "positive_step_1_N_1", .n = 1, .inc = 1, .want = 1.1, .x = &.{1.1} },
    .{ .name = "negative_step_1_N_1", .n = 1, .inc = 1, .want = 1.1, .x = &.{-1.1} },
    .{ .name = "positive_step_2_N_1", .n = 1, .inc = 2, .want = 1.1, .x = &.{ 1.1, 0 } },
    .{ .name = "negative_step_2_N_1", .n = 1, .inc = 2, .want = 1.1, .x = &.{ -1.1, 0 } },
    .{ .name = "positive_step_1_N_2", .n = 2, .inc = 1, .want = 1, .x = &.{ 1.1, 1 } },
    .{ .name = "negative_step_1_N_2", .n = 2, .inc = 1, .want = 1, .x = &.{ -1.1, 1 } },
    .{ .name = "positive_step_2_N_2", .n = 2, .inc = 2, .want = 1, .x = &.{ 1.1, 0, 1, 0 } },
    .{ .name = "negative_step_2_N_2", .n = 2, .inc = 2, .want = 1, .x = &.{ -1.1, 0, 1, 0 } },
    .{ .name = "positive_step_1_N_3", .n = 3, .inc = 1, .want = 1, .x = &.{ 1.1, 1, 2.2 } },
    .{ .name = "negative_step_1_N_3", .n = 3, .inc = 1, .want = 1, .x = &.{ -1.1, 1, -2.2 } },
    .{ .name = "positive_step_2_N_3", .n = 3, .inc = 2, .want = 1, .x = &.{ 1.1, 0, 1, 0, 2.2, 0 } },
    .{ .name = "negative_step_2_N_3", .n = 3, .inc = 2, .want = 1, .x = &.{ -1.1, 0, 1, 0, -2.2, 0 } },
    .{ .name = "positive_step_1_N_4", .n = 4, .inc = 1, .want = 1, .x = &.{ 1.1, 1, 2.2, 3.3 } },
    .{ .name = "negative_step_1_N_4", .n = 4, .inc = 1, .want = 1, .x = &.{ -1.1, 1, -2.2, -3.3 } },
    .{ .name = "positive_step_2_N_4", .n = 4, .inc = 2, .want = 1, .x = &.{ 1.1, 0, 1, 0, 2.2, 0, 3.3, 0 } },
    .{ .name = "negative_step_2_N_4", .n = 4, .inc = 2, .want = 1, .x = &.{ -1.1, 0, 1, 0, -2.2, 0, -3.3, 0 } },
    .{ .name = "positive_step_1_N_5", .n = 5, .inc = 1, .want = 0, .x = &.{ 1.1, 1, 2.2, 3.3, 0 } },
    .{ .name = "negative_step_1_N_5", .n = 5, .inc = 1, .want = 0, .x = &.{ -1.1, 1, -2.2, -3.3, 0 } },
    .{ .name = "positive_step_2_N_5", .n = 5, .inc = 2, .want = 0, .x = &.{ 1.1, 0, 1, 0, 2.2, 0, 3.3, 0, 0, 0 } },
    .{ .name = "negative_step_2_N_5", .n = 5, .inc = 2, .want = 0, .x = &.{ -1.1, 0, 1, 0, -2.2, 0, -3.3, 0, 0, 0 } },
    .{ .name = "positive_step_1_N_70", .n = 70, .inc = 1, .want = 0, .x = &.{ 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 0, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1025, 1026, 1027, 1028, 1029, 1030, 1031, 1032, 1033, 1034, 1035, 1036, 1037, 1038, 1039, 1040, 1041, 1042, 1043, 1044, 1045, 1046, 1047, 1048, 1049, 1050, 1051, 1052, 1053, 1054, 1055, 1056, 1057, 1058, 1059, 1060, 1061, 1062, 1063, 1064, 1065, 1066, 1067, 1068, 1069 } },
    .{ .name = "negative_step_1_N_70", .n = 70, .inc = 1, .want = 1, .x = &.{ -1000, -1001, -1002, -1003, -1004, -1005, -1006, -1007, -1, -1009, -1010, -1011, -1012, -1013, -1014, -1015, -1016, -1017, -1018, -1019, -1020, -1021, -1022, -1023, -1024, -1025, -1026, -1027, -1028, -1029, -1030, -1031, -1032, -1033, -1034, -1035, -1036, -1037, -1038, -1039, -1040, -1041, -1042, -1043, -1044, -1045, -1046, -1047, -1048, -1049, -1050, -1051, -1052, -1053, -1054, -1055, -1056, -1057, -1058, -1059, -1060, -1061, -1062, -1063, -1064, -1065, -1066, -1067, -1068, -1069 } },
    .{ .name = "positive_step_2_N_70", .n = 70, .inc = 2, .want = 1, .x = &.{ 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1025, 1026, 1027, 1028, 1029, 1030, 1031, 1032, 1033, 1034, 1035, 1036, 1037, 1038, 1039, 1040, 1041, 1042, 1043, 1044, 1045, 1046, 1047, 1048, 1049, 1050, 1051, 1052, 1053, 1054, 1055, 1056, 1057, 1058, 1059, 1060, 1061, 1062, 1063, 1064, 1065, 1066, 1067, 1068, 1069, 1070, 1071, 1072, 1073, 1074, 1075, 1076, 1077, 1078, 1079, 1080, 1081, 1082, 1083, 1084, 1085, 1086, 1087, 1088, 1089, 1090, 1091, 1092, 1093, 1094, 1095, 1096, 1097, 1098, 1099, 1100, 1101, 1102, 1103, 1104, 1105, 1106, 1107, 1108, 1109, 1110, 1111, 1112, 1113, 1114, 1115, 1116, 1117, 1118, 1119, 1120, 1121, 1122, 1123, 1124, 1125, 1126, 1127, 1128, 1129, 1130, 1131, 1132, 1133, 1134, 1135, 1136, 1137, 1138, 1139 } },
    .{ .name = "negative_step_2_N_70", .n = 70, .inc = 2, .want = 1, .x = &.{ -1000, -1001, -1002, -1003, -1004, -1005, -1006, -1007, -1008, -1009, -1010, -1011, -1012, -1013, -1014, -1015, -1, -1017, -1018, -1019, -1020, -1021, -1022, -1023, -1024, -1025, -1026, -1027, -1028, -1029, -1030, -1031, -1032, -1033, -1034, -1035, -1036, -1037, -1038, -1039, -1040, -1041, -1042, -1043, -1044, -1045, -1046, -1047, -1048, -1049, -1050, -1051, -1052, -1053, -1054, -1055, -1056, -1057, -1058, -1059, -1060, -1061, -1062, -1063, -1064, -1065, -1066, -1067, -1068, -1069, -1070, -1071, -1072, -1073, -1074, -1075, -1076, -1077, -1078, -1079, -1080, -1081, -1082, -1083, -1084, -1085, -1086, -1087, -1088, -1089, -1090, -1091, -1092, -1093, -1094, -1095, -1096, -1097, -1098, -1099, -1100, -1101, -1102, -1103, -1104, -1105, -1106, -1107, -1108, -1109, -1110, -1111, -1112, -1113, -1114, -1115, -1116, -1117, -1118, -1119, -1120, -1121, -1122, -1123, -1124, -1125, -1126, -1127, -1128, -1129, -1130, -1131, -1132, -1133, -1134, -1135, -1136, -1137, -1138, -1139 } },
};

/// `test_dzamax.c (test_scamax.c)`.
const amax_complex = [_]VecCase{
    .{ .name = "bad_args_N_0", .n = 0, .inc = 1, .want = 0, .x = &.{ 1000, 999, 998, 997, 996, 995, 994, 993, 992, 991, 990, 989, 988, 987, 986, 985, 984, 983, 982, 981, 980, 979, 978, 977, 976, 975, 974, 973, 972, 971, 970, 969, 968, 967, 966, 965, 964, 963, 962, 961, 960, 959, 958, 957, 956, 955, 954, 953, 952, 951, 950, 949, 948, 947, 946, 945, 944, 943, 942, 941, 940, 939, 938, 937, 936, 935, 934, 933, 932, 931, 930, 929, 928, 927, 926, 925, 924, 923, 922, 921, 920, 919, 918, 917, 916, 915, 914, 913, 912, 911, 910, 909, 908, 907, 906, 905, 904, 903, 902, 901, 900, 899, 898, 897, 896, 895, 894, 893, 892, 891, 890, 889, 888, 887, 886, 885, 884, 883, 882, 881, 880, 879, 878, 877, 876, 875, 874, 873, 872, 871, 870, 869, 868, 867, 866, 865, 864, 863, 862, 861 } },
    .{ .name = "step_zero", .n = 140, .inc = 0, .want = 1999, .x = &.{ -1000, -999, -998, -997, -996, -995, -994, -993, -992, -991, -990, -989, -988, -987, -986, -985, -984, -983, -982, -981, -980, -979, -978, -977, -976, -975, -974, -973, -972, -971, -970, -969, -968, -967, -966, -965, -964, -963, -962, -961, -960, -959, -958, -957, -956, -955, -954, -953, -952, -951, -950, -949, -948, -947, -946, -945, -944, -943, -942, -941, -940, -939, -938, -937, -936, -935, -934, -933, -932, -931, -930, -929, -928, -927, -926, -925, -924, -923, -922, -921, -920, -919, -918, -917, -916, -915, -914, -913, -912, -911, -910, -909, -908, -907, -906, -905, -904, -903, -902, -901, -900, -899, -898, -897, -896, -895, -894, -893, -892, -891, -890, -889, -888, -887, -886, -885, -884, -883, -882, -881, -880, -879, -878, -877, -876, -875, -874, -873, -872, -871, -870, -869, -868, -867, -866, -865, -864, -863, -862, -861 } },
    .{ .name = "positive_step_1_N_1", .n = 1, .inc = 1, .want = 3, .x = &.{ 1, 2 } },
    .{ .name = "negative_step_1_N_1", .n = 1, .inc = 1, .want = 3, .x = &.{ -1, -2 } },
    .{ .name = "positive_step_2_N_1", .n = 1, .inc = 2, .want = 3, .x = &.{ 1, 2, 0, 0 } },
    .{ .name = "negative_step_2_N_1", .n = 1, .inc = 2, .want = 3, .x = &.{ -1, -2, 0, 0 } },
    .{ .name = "positive_step_1_N_2", .n = 2, .inc = 1, .want = 3, .x = &.{ 1, 2, 0, 0 } },
    .{ .name = "negative_step_1_N_2", .n = 2, .inc = 1, .want = 3, .x = &.{ -1, -2, 0, 0 } },
    .{ .name = "positive_step_2_N_2", .n = 2, .inc = 2, .want = 3, .x = &.{ 1, 2, 0, 0, 1, 1, 0, 0 } },
    .{ .name = "negative_step_2_N_2", .n = 2, .inc = 2, .want = 3, .x = &.{ -1, -2, 0, 0, -1, -1, 0, 0 } },
    .{ .name = "positive_step_1_N_3", .n = 3, .inc = 1, .want = 3, .x = &.{ 1, 2, 0, 0, 2, 1 } },
    .{ .name = "negative_step_1_N_3", .n = 3, .inc = 1, .want = 4, .x = &.{ -1, -2, 0, 0, -3, -1 } },
    .{ .name = "positive_step_2_N_3", .n = 3, .inc = 2, .want = 4, .x = &.{ 1, 2, 0, 0, 1, 1, 0, 0, 3, 1, 0, 0 } },
    .{ .name = "negative_step_2_N_3", .n = 3, .inc = 2, .want = 4, .x = &.{ -1, -2, 0, 0, -1, -1, 0, 0, -3, -1, 0, 0 } },
    .{ .name = "positive_step_1_N_4", .n = 4, .inc = 1, .want = 4, .x = &.{ 1, 2, 0, 0, 2, 1, -2, -2 } },
    .{ .name = "negative_step_1_N_4", .n = 4, .inc = 1, .want = 4, .x = &.{ -1, -2, 0, 0, -2, -1, -2, -2 } },
    .{ .name = "positive_step_2_N_4", .n = 4, .inc = 2, .want = 4, .x = &.{ 1, 2, 0, 0, 1, 1, 0, 0, 2, 1, 0, 0, -2, -2, 0, 0 } },
    .{ .name = "negative_step_2_N_4", .n = 4, .inc = 2, .want = 4, .x = &.{ -1, -2, 0, 0, -1, -1, 0, 0, -2, -1, 0, 0, -2, -2, 0, 0 } },
    .{ .name = "positive_step_1_N_70", .n = 70, .inc = 1, .want = 2000, .x = &.{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 1000, 1000, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139 } },
    .{ .name = "negative_step_1_N_70", .n = 70, .inc = 1, .want = 2000, .x = &.{ 0, -1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13, 1000, 1000, -16, -17, -18, -19, -20, -21, -22, -23, -24, -25, -26, -27, -28, -29, -30, -31, -32, -33, -34, -35, -36, -37, -38, -39, -40, -41, -42, -43, -44, -45, -46, -47, -48, -49, -50, -51, -52, -53, -54, -55, -56, -57, -58, -59, -60, -61, -62, -63, -64, -65, -66, -67, -68, -69, -70, -71, -72, -73, -74, -75, -76, -77, -78, -79, -80, -81, -82, -83, -84, -85, -86, -87, -88, -89, -90, -91, -92, -93, -94, -95, -96, -97, -98, -99, -100, -101, -102, -103, -104, -105, -106, -107, -108, -109, -110, -111, -112, -113, -114, -115, -116, -117, -118, -119, -120, -121, -122, -123, -124, -125, -126, -127, -128, -129, -130, -131, -132, -133, -134, -135, -136, -137, -138, -139 } },
    .{ .name = "positive_step_2_N_70", .n = 70, .inc = 2, .want = 2000, .x = &.{ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 1000, 1000, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279 } },
    .{ .name = "negative_step_2_N_70", .n = 70, .inc = 2, .want = 2000, .x = &.{ 0, -1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13, -14, -15, -16, -17, -18, -19, -20, -21, -22, -23, -24, -25, -26, -27, 1000, 1000, -30, -31, -32, -33, -34, -35, -36, -37, -38, -39, -40, -41, -42, -43, -44, -45, -46, -47, -48, -49, -50, -51, -52, -53, -54, -55, -56, -57, -58, -59, -60, -61, -62, -63, -64, -65, -66, -67, -68, -69, -70, -71, -72, -73, -74, -75, -76, -77, -78, -79, -80, -81, -82, -83, -84, -85, -86, -87, -88, -89, -90, -91, -92, -93, -94, -95, -96, -97, -98, -99, -100, -101, -102, -103, -104, -105, -106, -107, -108, -109, -110, -111, -112, -113, -114, -115, -116, -117, -118, -119, -120, -121, -122, -123, -124, -125, -126, -127, -128, -129, -130, -131, -132, -133, -134, -135, -136, -137, -138, -139, -140, -141, -142, -143, -144, -145, -146, -147, -148, -149, -150, -151, -152, -153, -154, -155, -156, -157, -158, -159, -160, -161, -162, -163, -164, -165, -166, -167, -168, -169, -170, -171, -172, -173, -174, -175, -176, -177, -178, -179, -180, -181, -182, -183, -184, -185, -186, -187, -188, -189, -190, -191, -192, -193, -194, -195, -196, -197, -198, -199, -200, -201, -202, -203, -204, -205, -206, -207, -208, -209, -210, -211, -212, -213, -214, -215, -216, -217, -218, -219, -220, -221, -222, -223, -224, -225, -226, -227, -228, -229, -230, -231, -232, -233, -234, -235, -236, -237, -238, -239, -240, -241, -242, -243, -244, -245, -246, -247, -248, -249, -250, -251, -252, -253, -254, -255, -256, -257, -258, -259, -260, -261, -262, -263, -264, -265, -266, -267, -268, -269, -270, -271, -272, -273, -274, -275, -276, -277, -278, -279 } },
};

/// `test_dzamin.c (test_scamin.c)`.
const amin_complex = [_]VecCase{
    .{ .name = "bad_args_N_0", .n = 0, .inc = 1, .want = 0, .x = &.{ 1000, 999, 998, 997, 996, 995, 994, 993, 992, 991, 990, 989, 988, 987, 986, 985, 984, 983, 982, 981, 980, 979, 978, 977, 976, 975, 974, 973, 972, 971, 970, 969, 968, 967, 966, 965, 964, 963, 962, 961, 960, 959, 958, 957, 956, 955, 954, 953, 952, 951, 950, 949, 948, 947, 946, 945, 944, 943, 942, 941, 940, 939, 938, 937, 936, 935, 934, 933, 932, 931, 930, 929, 928, 927, 926, 925, 924, 923, 922, 921, 920, 919, 918, 917, 916, 915, 914, 913, 912, 911, 910, 909, 908, 907, 906, 905, 904, 903, 902, 901, 900, 899, 898, 897, 896, 895, 894, 893, 892, 891, 890, 889, 888, 887, 886, 885, 884, 883, 882, 881, 880, 879, 878, 877, 876, 875, 874, 873, 872, 871, 870, 869, 868, 867, 866, 865, 864, 863, 862, 861 } },
    .{ .name = "step_zero", .n = 140, .inc = 0, .want = 1999, .x = &.{ -1000, -999, -998, -997, -996, -995, -994, -993, 0, -991, -990, -989, -988, -987, -986, -985, -984, -983, -982, -981, -980, -979, -978, -977, -976, -975, -974, -973, -972, -971, -970, -969, -968, -967, -966, -965, -964, -963, -962, -961, -960, -959, -958, -957, -956, -955, -954, -953, -952, -951, -950, -949, -948, -947, -946, -945, -944, -943, -942, -941, -940, -939, -938, -937, -936, -935, -934, -933, -932, -931, -930, -929, -928, -927, -926, -925, -924, -923, -922, -921, -920, -919, -918, -917, -916, -915, -914, -913, -912, -911, -910, -909, -908, -907, -906, -905, -904, -903, -902, -901, -900, -899, -898, -897, -896, -895, -894, -893, -892, -891, -890, -889, -888, -887, -886, -885, -884, -883, -882, -881, -880, -879, -878, -877, -876, -875, -874, -873, -872, -871, -870, -869, -868, -867, -866, -865, -864, -863, -862, -861 } },
    .{ .name = "positive_step_1_N_1", .n = 1, .inc = 1, .want = 3, .x = &.{ 1, 2 } },
    .{ .name = "negative_step_1_N_1", .n = 1, .inc = 1, .want = 3, .x = &.{ -1, -2 } },
    .{ .name = "positive_step_2_N_1", .n = 1, .inc = 2, .want = 3, .x = &.{ 1, 2, 0, 0 } },
    .{ .name = "negative_step_2_N_1", .n = 1, .inc = 2, .want = 3, .x = &.{ -1, -2, 0, 0 } },
    .{ .name = "positive_step_1_N_2", .n = 2, .inc = 1, .want = 0, .x = &.{ 1, 2, 0, 0 } },
    .{ .name = "negative_step_1_N_2", .n = 2, .inc = 1, .want = 0, .x = &.{ -1, -2, 0, 0 } },
    .{ .name = "positive_step_2_N_2", .n = 2, .inc = 2, .want = 2, .x = &.{ 1, 2, 0, 0, 1, 1, 0, 0 } },
    .{ .name = "negative_step_2_N_2", .n = 2, .inc = 2, .want = 2, .x = &.{ -1, -2, 0, 0, -1, -1, 0, 0 } },
    .{ .name = "positive_step_1_N_3", .n = 3, .inc = 1, .want = 0, .x = &.{ 1, 2, 0, 0, 2, 1 } },
    .{ .name = "negative_step_1_N_3", .n = 3, .inc = 1, .want = 0, .x = &.{ -1, -2, 0, 0, -2, -1 } },
    .{ .name = "positive_step_2_N_3", .n = 3, .inc = 2, .want = 2, .x = &.{ 1, 2, 0, 0, 1, 1, 0, 0, 2, 1, 0, 0 } },
    .{ .name = "negative_step_2_N_3", .n = 3, .inc = 2, .want = 2, .x = &.{ -1, -2, 0, 0, -1, -1, 0, 0, -2, -1, 0, 0 } },
    .{ .name = "positive_step_1_N_4", .n = 4, .inc = 1, .want = 0, .x = &.{ 1, 2, 0, 0, 2, 1, -2, -2 } },
    .{ .name = "negative_step_1_N_4", .n = 4, .inc = 1, .want = 0, .x = &.{ -1, -2, 0, 0, -2, -1, -2, -2 } },
    .{ .name = "positive_step_2_N_4", .n = 4, .inc = 2, .want = 2, .x = &.{ 1, 2, 0, 0, 1, 1, 0, 0, 2, 1, 0, 0, -2, -2, 0, 0 } },
    .{ .name = "negative_step_2_N_4", .n = 4, .inc = 2, .want = 2, .x = &.{ -1, -2, 0, 0, -1, -1, 0, 0, -2, -1, 0, 0, -2, -2, 0, 0 } },
    .{ .name = "positive_step_1_N_70", .n = 70, .inc = 1, .want = 0, .x = &.{ 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 0, 0, 1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1025, 1026, 1027, 1028, 1029, 1030, 1031, 1032, 1033, 1034, 1035, 1036, 1037, 1038, 1039, 1040, 1041, 1042, 1043, 1044, 1045, 1046, 1047, 1048, 1049, 1050, 1051, 1052, 1053, 1054, 1055, 1056, 1057, 1058, 1059, 1060, 1061, 1062, 1063, 1064, 1065, 1066, 1067, 1068, 1069, 1070, 1071, 1072, 1073, 1074, 1075, 1076, 1077, 1078, 1079, 1080, 1081, 1082, 1083, 1084, 1085, 1086, 1087, 1088, 1089, 1090, 1091, 1092, 1093, 1094, 1095, 1096, 1097, 1098, 1099, 1100, 1101, 1102, 1103, 1104, 1105, 1106, 1107, 1108, 1109, 1110, 1111, 1112, 1113, 1114, 1115, 1116, 1117, 1118, 1119, 1120, 1121, 1122, 1123, 1124, 1125, 1126, 1127, 1128, 1129, 1130, 1131, 1132, 1133, 1134, 1135, 1136, 1137, 1138, 1139 } },
    .{ .name = "negative_step_1_N_70", .n = 70, .inc = 1, .want = 0, .x = &.{ -1000, -1001, -1002, -1003, -1004, -1005, -1006, -1007, -1008, -1009, -1010, -1011, -1012, -1013, 0, 0, -1016, -1017, -1018, -1019, -1020, -1021, -1022, -1023, -1024, -1025, -1026, -1027, -1028, -1029, -1030, -1031, -1032, -1033, -1034, -1035, -1036, -1037, -1038, -1039, -1040, -1041, -1042, -1043, -1044, -1045, -1046, -1047, -1048, -1049, -1050, -1051, -1052, -1053, -1054, -1055, -1056, -1057, -1058, -1059, -1060, -1061, -1062, -1063, -1064, -1065, -1066, -1067, -1068, -1069, -1070, -1071, -1072, -1073, -1074, -1075, -1076, -1077, -1078, -1079, -1080, -1081, -1082, -1083, -1084, -1085, -1086, -1087, -1088, -1089, -1090, -1091, -1092, -1093, -1094, -1095, -1096, -1097, -1098, -1099, -1100, -1101, -1102, -1103, -1104, -1105, -1106, -1107, -1108, -1109, -1110, -1111, -1112, -1113, -1114, -1115, -1116, -1117, -1118, -1119, -1120, -1121, -1122, -1123, -1124, -1125, -1126, -1127, -1128, -1129, -1130, -1131, -1132, -1133, -1134, -1135, -1136, -1137, -1138, -1139 } },
    .{ .name = "positive_step_2_N_70", .n = 70, .inc = 2, .want = 0, .x = &.{ 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1025, 1026, 1027, 0, 0, 1030, 1031, 1032, 1033, 1034, 1035, 1036, 1037, 1038, 1039, 1040, 1041, 1042, 1043, 1044, 1045, 1046, 1047, 1048, 1049, 1050, 1051, 1052, 1053, 1054, 1055, 1056, 1057, 1058, 1059, 1060, 1061, 1062, 1063, 1064, 1065, 1066, 1067, 1068, 1069, 1070, 1071, 1072, 1073, 1074, 1075, 1076, 1077, 1078, 1079, 1080, 1081, 1082, 1083, 1084, 1085, 1086, 1087, 1088, 1089, 1090, 1091, 1092, 1093, 1094, 1095, 1096, 1097, 1098, 1099, 1100, 1101, 1102, 1103, 1104, 1105, 1106, 1107, 1108, 1109, 1110, 1111, 1112, 1113, 1114, 1115, 1116, 1117, 1118, 1119, 1120, 1121, 1122, 1123, 1124, 1125, 1126, 1127, 1128, 1129, 1130, 1131, 1132, 1133, 1134, 1135, 1136, 1137, 1138, 1139, 1140, 1141, 1142, 1143, 1144, 1145, 1146, 1147, 1148, 1149, 1150, 1151, 1152, 1153, 1154, 1155, 1156, 1157, 1158, 1159, 1160, 1161, 1162, 1163, 1164, 1165, 1166, 1167, 1168, 1169, 1170, 1171, 1172, 1173, 1174, 1175, 1176, 1177, 1178, 1179, 1180, 1181, 1182, 1183, 1184, 1185, 1186, 1187, 1188, 1189, 1190, 1191, 1192, 1193, 1194, 1195, 1196, 1197, 1198, 1199, 1200, 1201, 1202, 1203, 1204, 1205, 1206, 1207, 1208, 1209, 1210, 1211, 1212, 1213, 1214, 1215, 1216, 1217, 1218, 1219, 1220, 1221, 1222, 1223, 1224, 1225, 1226, 1227, 1228, 1229, 1230, 1231, 1232, 1233, 1234, 1235, 1236, 1237, 1238, 1239, 1240, 1241, 1242, 1243, 1244, 1245, 1246, 1247, 1248, 1249, 1250, 1251, 1252, 1253, 1254, 1255, 1256, 1257, 1258, 1259, 1260, 1261, 1262, 1263, 1264, 1265, 1266, 1267, 1268, 1269, 1270, 1271, 1272, 1273, 1274, 1275, 1276, 1277, 1278, 1279 } },
    .{ .name = "negative_step_2_N_70", .n = 70, .inc = 2, .want = 0, .x = &.{ -1000, -1001, -1002, -1003, -1004, -1005, -1006, -1007, -1008, -1009, -1010, -1011, -1012, -1013, -1014, -1015, -1016, -1017, -1018, -1019, -1020, -1021, -1022, -1023, -1024, -1025, -1026, -1027, 0, 0, -1030, -1031, -1032, -1033, -1034, -1035, -1036, -1037, -1038, -1039, -1040, -1041, -1042, -1043, -1044, -1045, -1046, -1047, -1048, -1049, -1050, -1051, -1052, -1053, -1054, -1055, -1056, -1057, -1058, -1059, -1060, -1061, -1062, -1063, -1064, -1065, -1066, -1067, -1068, -1069, -1070, -1071, -1072, -1073, -1074, -1075, -1076, -1077, -1078, -1079, -1080, -1081, -1082, -1083, -1084, -1085, -1086, -1087, -1088, -1089, -1090, -1091, -1092, -1093, -1094, -1095, -1096, -1097, -1098, -1099, -1100, -1101, -1102, -1103, -1104, -1105, -1106, -1107, -1108, -1109, -1110, -1111, -1112, -1113, -1114, -1115, -1116, -1117, -1118, -1119, -1120, -1121, -1122, -1123, -1124, -1125, -1126, -1127, -1128, -1129, -1130, -1131, -1132, -1133, -1134, -1135, -1136, -1137, -1138, -1139, -1140, -1141, -1142, -1143, -1144, -1145, -1146, -1147, -1148, -1149, -1150, -1151, -1152, -1153, -1154, -1155, -1156, -1157, -1158, -1159, -1160, -1161, -1162, -1163, -1164, -1165, -1166, -1167, -1168, -1169, -1170, -1171, -1172, -1173, -1174, -1175, -1176, -1177, -1178, -1179, -1180, -1181, -1182, -1183, -1184, -1185, -1186, -1187, -1188, -1189, -1190, -1191, -1192, -1193, -1194, -1195, -1196, -1197, -1198, -1199, -1200, -1201, -1202, -1203, -1204, -1205, -1206, -1207, -1208, -1209, -1210, -1211, -1212, -1213, -1214, -1215, -1216, -1217, -1218, -1219, -1220, -1221, -1222, -1223, -1224, -1225, -1226, -1227, -1228, -1229, -1230, -1231, -1232, -1233, -1234, -1235, -1236, -1237, -1238, -1239, -1240, -1241, -1242, -1243, -1244, -1245, -1246, -1247, -1248, -1249, -1250, -1251, -1252, -1253, -1254, -1255, -1256, -1257, -1258, -1259, -1260, -1261, -1262, -1263, -1264, -1265, -1266, -1267, -1268, -1269, -1270, -1271, -1272, -1273, -1274, -1275, -1276, -1277, -1278, -1279 } },
};

/// `test_dsum.c (test_ssum.c)`.
const sum_real = [_]VecCase{
    .{ .name = "bad_args_N_0", .n = 0, .inc = 1, .want = 0, .x = &.{ 1000, 999, 998, 997, 996, 995, 994, 993, 992, 991, 990, 989, 988, 987, 986, 985, 984, 983, 982, 981, 980, 979, 978, 977, 976, 975, 974, 973, 972, 971, 970, 969, 968, 967, 966, 965, 964, 963, 962, 961, 960, 959, 958, 957, 956, 955, 954, 953, 952, 951 } },
    .{ .name = "step_zero", .n = 50, .inc = 0, .want = 0, .x = &.{ 0, 777, 777, 777, 777, 777, 777, 777, 0, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777 } },
    .{ .name = "step_1_N_1", .n = 1, .inc = 1, .want = 1.1, .x = &.{1.1} },
    .{ .name = "step_2_N_1", .n = 1, .inc = 2, .want = 1.1, .x = &.{ 1.1, 0 } },
    .{ .name = "step_1_N_2", .n = 2, .inc = 1, .want = 0.1, .x = &.{ 1.1, -1 } },
    .{ .name = "step_2_N_2", .n = 2, .inc = 2, .want = 2.1, .x = &.{ 1.1, -1.5, 1, 1 } },
    .{ .name = "step_1_N_3", .n = 3, .inc = 1, .want = 4.3, .x = &.{ 1.1, 1, 2.2 } },
    .{ .name = "step_2_N_3", .n = 3, .inc = 2, .want = 2.3, .x = &.{ 1.1, 0, -1, -3, 2.2, 3 } },
    .{ .name = "step_1_N_4", .n = 4, .inc = 1, .want = 3.2, .x = &.{ 1.1, 1, -2.2, 3.3 } },
    .{ .name = "step_2_N_4", .n = 4, .inc = 2, .want = 1, .x = &.{ 1.1, 0, 1, 2, 2.2, 2.7, -3.3, -5.9 } },
    .{ .name = "step_1_N_5", .n = 5, .inc = 1, .want = 6.5, .x = &.{ 0, 1, 2.2, 3.3, 0 } },
    .{ .name = "step_2_N_5", .n = 5, .inc = 2, .want = 6.5, .x = &.{ 0, 3, 1, -2.2, 2.2, -1.7, 3.3, 14.5, 0, -9 } },
    .{ .name = "step_1_N_50", .n = 50, .inc = 1, .want = 0, .x = &.{ 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1 } },
    .{ .name = "step_2_N_50", .n = 50, .inc = 2, .want = 50, .x = &.{ 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1 } },
};

/// `test_dzsum.c (test_scsum.c)`.
const sum_complex = [_]VecCase{
    .{ .name = "bad_args_N_0", .n = 0, .inc = 1, .want = 0, .x = &.{ 1000, 999, 998, 997, 996, 995, 994, 993, 992, 991, 990, 989, 988, 987, 986, 985, 984, 983, 982, 981, 980, 979, 978, 977, 976, 975, 974, 973, 972, 971, 970, 969, 968, 967, 966, 965, 964, 963, 962, 961, 960, 959, 958, 957, 956, 955, 954, 953, 952, 951, 950, 949, 948, 947, 946, 945, 944, 943, 942, 941, 940, 939, 938, 937, 936, 935, 934, 933, 932, 931, 930, 929, 928, 927, 926, 925, 924, 923, 922, 921, 920, 919, 918, 917, 916, 915, 914, 913, 912, 911, 910, 909, 908, 907, 906, 905, 904, 903, 902, 901 } },
    .{ .name = "step_zero", .n = 50, .inc = 0, .want = 0, .x = &.{ 0, 777, 777, 777, 777, 777, 777, 777, 0, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777, 777 } },
    .{ .name = "step_1_N_1", .n = 1, .inc = 1, .want = 0.1, .x = &.{ 1.1, -1 } },
    .{ .name = "step_2_N_1", .n = 1, .inc = 2, .want = 1.1, .x = &.{ 1.1, 0, 2.3, -1 } },
    .{ .name = "step_1_N_2", .n = 2, .inc = 1, .want = 1.4, .x = &.{ 1.1, -1, 2.3, -1 } },
    .{ .name = "step_2_N_2", .n = 2, .inc = 2, .want = 1.6, .x = &.{ 1.1, -1.5, 1.1, -1, 1, 1, 1.1, -1 } },
    .{ .name = "step_1_N_3", .n = 3, .inc = 1, .want = 4.4, .x = &.{ 1.1, 1, 2.2, 1.1, -1, 0 } },
    .{ .name = "step_2_N_3", .n = 3, .inc = 2, .want = 2.3, .x = &.{ 1.1, 0, -1, 0, -1, -3, -1, 0, 2.2, 3, -1, 0 } },
    .{ .name = "step_1_N_4", .n = 4, .inc = 1, .want = 6.4, .x = &.{ 1.1, 1, -2.2, 3.3, 1.1, 1, -2.2, 3.3 } },
    .{ .name = "step_2_N_4", .n = 4, .inc = 2, .want = -0.2, .x = &.{ 1.1, 0, 1.1, 1, 1, 2, 1.1, 1, 2.2, 2.7, 1.1, 1, -3.3, -5.9 } },
    .{ .name = "step_1_N_5", .n = 5, .inc = 1, .want = 13, .x = &.{ 0, 1, 2.2, 3.3, 0, 0, 1, 2.2, 3.3, 0 } },
    .{ .name = "step_2_N_5", .n = 5, .inc = 2, .want = 11.1, .x = &.{ 0, 3, 1, 2.2, 1, -2.2, 1, 2.2, 2.2, -1.7, 1, 2.2, 3.3, 14.5, 1, 2.2, 0, -9, 1, 2.2 } },
    .{ .name = "step_1_N_50", .n = 50, .inc = 1, .want = 0, .x = &.{ 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1 } },
    .{ .name = "step_2_N_50", .n = 50, .inc = 2, .want = 0, .x = &.{ 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1 } },
};

/// `test_idamin.c (test_isamin.c)`.
const iamin_real = [_]VecCase{
    .{ .name = "bad_args_N_0", .n = 0, .inc = 1, .want = 0, .x = &.{ 1000, 999, 998, 997, 996, 995, 994, 993, 992, 991, 990, 989, 988, 987, 986, 985, 984, 983, 982, 981, 980, 979, 978, 977, 976, 975, 974, 973, 972, 971, 970, 969, 968, 967, 966, 965, 964, 963, 962, 961, 960, 959, 958, 957, 956, 955, 954, 953, 952, 951 } },
    .{ .name = "step_zero", .n = 50, .inc = 0, .want = 0, .x = &.{ 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 0, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1025, 1026, 1027, 1028, 1029, 1030, 1031, 1032, 1033, 1034, 1035, 1036, 1037, 1038, 1039, 1040, 1041, 1042, 1043, 1044, 1045, 1046, 1047, 1048, 1049 } },
    .{ .name = "positive_step_1_N_1", .n = 1, .inc = 1, .want = 1, .x = &.{1.1} },
    .{ .name = "negative_step_1_N_1", .n = 1, .inc = 1, .want = 1, .x = &.{-1.1} },
    .{ .name = "positive_step_2_N_1", .n = 1, .inc = 2, .want = 1, .x = &.{ 1.1, 0 } },
    .{ .name = "negative_step_2_N_1", .n = 1, .inc = 2, .want = 1, .x = &.{ -1.1, 0 } },
    .{ .name = "positive_step_1_N_2", .n = 2, .inc = 1, .want = 2, .x = &.{ 1.1, 1 } },
    .{ .name = "negative_step_1_N_2", .n = 2, .inc = 1, .want = 2, .x = &.{ -1.1, 1 } },
    .{ .name = "positive_step_2_N_2", .n = 2, .inc = 2, .want = 2, .x = &.{ 1.1, 0, 1, 0 } },
    .{ .name = "negative_step_2_N_2", .n = 2, .inc = 2, .want = 2, .x = &.{ -1.1, 0, 1, 0 } },
    .{ .name = "positive_step_1_N_3", .n = 3, .inc = 1, .want = 2, .x = &.{ 1.1, 1, 2.2 } },
    .{ .name = "negative_step_1_N_3", .n = 3, .inc = 1, .want = 2, .x = &.{ -1.1, 1, -2.2 } },
    .{ .name = "positive_step_2_N_3", .n = 3, .inc = 2, .want = 2, .x = &.{ 1.1, 0, 1, 0, 2.2, 0 } },
    .{ .name = "negative_step_2_N_3", .n = 3, .inc = 2, .want = 2, .x = &.{ -1.1, 0, 1, 0, -2.2, 0 } },
    .{ .name = "positive_step_1_N_4", .n = 4, .inc = 1, .want = 2, .x = &.{ 1.1, 1, 2.2, 3.3 } },
    .{ .name = "negative_step_1_N_4", .n = 4, .inc = 1, .want = 2, .x = &.{ -1.1, 1, -2.2, -3.3 } },
    .{ .name = "positive_step_2_N_4", .n = 4, .inc = 2, .want = 2, .x = &.{ 1.1, 0, 1, 0, 2.2, 0, 3.3, 0 } },
    .{ .name = "negative_step_2_N_4", .n = 4, .inc = 2, .want = 2, .x = &.{ -1.1, 0, 1, 0, -2.2, 0, -3.3, 0 } },
    .{ .name = "positive_step_1_N_5", .n = 5, .inc = 1, .want = 5, .x = &.{ 1.1, 1, 2.2, 3.3, 0 } },
    .{ .name = "negative_step_1_N_5", .n = 5, .inc = 1, .want = 5, .x = &.{ -1.1, 1, -2.2, -3.3, 0 } },
    .{ .name = "positive_step_2_N_5", .n = 5, .inc = 2, .want = 5, .x = &.{ 1.1, 0, 1, 0, 2.2, 0, 3.3, 0, 0, 0 } },
    .{ .name = "negative_step_2_N_5", .n = 5, .inc = 2, .want = 5, .x = &.{ -1.1, 0, 1, 0, -2.2, 0, -3.3, 0, 0, 0 } },
    .{ .name = "positive_step_1_N_50", .n = 50, .inc = 1, .want = 9, .x = &.{ 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 0, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1025, 1026, 1027, 1028, 1029, 1030, 1031, 1032, 1033, 1034, 1035, 1036, 1037, 1038, 1039, 1040, 1041, 1042, 1043, 1044, 1045, 1046, 1047, 1048, 1049 } },
    .{ .name = "negative_step_1_N_50", .n = 50, .inc = 1, .want = 9, .x = &.{ -1000, -1001, -1002, -1003, -1004, -1005, -1006, -1007, -1, -1009, -1010, -1011, -1012, -1013, -1014, -1015, -1016, -1017, -1018, -1019, -1020, -1021, -1022, -1023, -1024, -1025, -1026, -1027, -1028, -1029, -1030, -1031, -1032, -1033, -1034, -1035, -1036, -1037, -1038, -1039, -1040, -1041, -1042, -1043, -1044, -1045, -1046, -1047, -1048, -1049 } },
    .{ .name = "positive_step_2_N_50", .n = 50, .inc = 2, .want = 9, .x = &.{ 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 0, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1025, 1026, 1027, 1028, 1029, 1030, 1031, 1032, 1033, 1034, 1035, 1036, 1037, 1038, 1039, 1040, 1041, 1042, 1043, 1044, 1045, 1046, 1047, 1048, 1049, 1050, 1051, 1052, 1053, 1054, 1055, 1056, 1057, 1058, 1059, 1060, 1061, 1062, 1063, 1064, 1065, 1066, 1067, 1068, 1069, 1070, 1071, 1072, 1073, 1074, 1075, 1076, 1077, 1078, 1079, 1080, 1081, 1082, 1083, 1084, 1085, 1086, 1087, 1088, 1089, 1090, 1091, 1092, 1093, 1094, 1095, 1096, 1097, 1098, 1099 } },
    .{ .name = "negative_step_2_N_50", .n = 50, .inc = 2, .want = 9, .x = &.{ -1000, -1001, -1002, -1003, -1004, -1005, -1006, -1007, -1008, -1009, -1010, -1011, -1012, -1013, -1014, -1015, -1, -1017, -1018, -1019, -1020, -1021, -1022, -1023, -1024, -1025, -1026, -1027, -1028, -1029, -1030, -1031, -1032, -1033, -1034, -1035, -1036, -1037, -1038, -1039, -1040, -1041, -1042, -1043, -1044, -1045, -1046, -1047, -1048, -1049, -1050, -1051, -1052, -1053, -1054, -1055, -1056, -1057, -1058, -1059, -1060, -1061, -1062, -1063, -1064, -1065, -1066, -1067, -1068, -1069, -1070, -1071, -1072, -1073, -1074, -1075, -1076, -1077, -1078, -1079, -1080, -1081, -1082, -1083, -1084, -1085, -1086, -1087, -1088, -1089, -1090, -1091, -1092, -1093, -1094, -1095, -1096, -1097, -1098, -1099 } },
    .{ .name = "min_idx_in_vec_tail", .n = 50, .inc = 2, .want = 50, .x = &.{ 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1025, 1026, 1027, 1028, 1029, 1030, 1031, 1032, 1033, 1034, 1035, 1036, 1037, 1038, 1039, 1040, 1041, 1042, 1043, 1044, 1045, 1046, 1047, 1048, 1049, 1050, 1051, 1052, 1053, 1054, 1055, 1056, 1057, 1058, 1059, 1060, 1061, 1062, 1063, 1064, 1065, 1066, 1067, 1068, 1069, 1070, 1071, 1072, 1073, 1074, 1075, 1076, 1077, 1078, 1079, 1080, 1081, 1082, 1083, 1084, 1085, 1086, 1087, 1088, 1089, 1090, 1091, 1092, 1093, 1094, 1095, 1096, 1097, 0, 1099 } },
    .{ .name = "min_idx_in_vec_tail_inc_1", .n = 50, .inc = 1, .want = 50, .x = &.{ 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1025, 1026, 1027, 1028, 1029, 1030, 1031, 1032, 1033, 1034, 1035, 1036, 1037, 1038, 1039, 1040, 1041, 1042, 1043, 1044, 1045, 1046, 1047, 1048, 0 } },
};

/// `test_izamin.c (test_icamin.c)`.
const iamin_complex = [_]VecCase{
    .{ .name = "bad_args_N_0", .n = 0, .inc = 1, .want = 0, .x = &.{ 1000, 999, 998, 997, 996, 995, 994, 993, 992, 991, 990, 989, 988, 987, 986, 985, 984, 983, 982, 981, 980, 979, 978, 977, 976, 975, 974, 973, 972, 971, 970, 969, 968, 967, 966, 965, 964, 963, 962, 961, 960, 959, 958, 957, 956, 955, 954, 953, 952, 951, 950, 949, 948, 947, 946, 945, 944, 943, 942, 941, 940, 939, 938, 937, 936, 935, 934, 933, 932, 931, 930, 929, 928, 927, 926, 925, 924, 923, 922, 921, 920, 919, 918, 917, 916, 915, 914, 913, 912, 911, 910, 909, 908, 907, 906, 905, 904, 903, 902, 901 } },
    .{ .name = "step_zero", .n = 50, .inc = 0, .want = 0, .x = &.{ -1000, -999, -998, -997, -996, -995, -994, -993, -992, -991, -990, -989, -988, -987, -986, -985, -984, -983, -982, -981, -980, -979, -978, -977, -976, -975, -974, -973, -972, -971, -970, -969, -968, -967, -966, -965, -964, -963, -962, -961, -960, -959, -958, -957, -956, -955, -954, -953, -952, -951, -950, -949, -948, -947, -946, -945, -944, -943, -942, -941, -940, -939, -938, -937, -936, -935, -934, -933, -932, -931, -930, -929, -928, -927, -926, -925, -924, -923, -922, -921, -920, -919, -918, -917, -916, -915, -914, -913, -912, -911, -910, -909, -908, -907, -906, -905, -904, -903, -902, -901 } },
    .{ .name = "positive_step_1_N_1", .n = 1, .inc = 1, .want = 1, .x = &.{ 1, 2 } },
    .{ .name = "negative_step_1_N_1", .n = 1, .inc = 1, .want = 1, .x = &.{ -1, -2 } },
    .{ .name = "positive_step_2_N_1", .n = 1, .inc = 2, .want = 1, .x = &.{ 1, 2, 0, 0 } },
    .{ .name = "negative_step_2_N_1", .n = 1, .inc = 2, .want = 1, .x = &.{ -1, -2, 0, 0 } },
    .{ .name = "positive_step_1_N_2", .n = 2, .inc = 1, .want = 2, .x = &.{ 1, 2, 0, 0 } },
    .{ .name = "negative_step_1_N_2", .n = 2, .inc = 1, .want = 2, .x = &.{ -1, -2, 0, 0 } },
    .{ .name = "positive_step_2_N_2", .n = 2, .inc = 2, .want = 2, .x = &.{ 1, 2, 0, 0, 1, 1, 0, 0 } },
    .{ .name = "negative_step_2_N_2", .n = 2, .inc = 2, .want = 2, .x = &.{ -1, -2, 0, 0, -1, -1, 0, 0 } },
    .{ .name = "positive_step_1_N_3", .n = 3, .inc = 1, .want = 2, .x = &.{ 1, 2, 0, 0, 2, 1 } },
    .{ .name = "negative_step_1_N_3", .n = 3, .inc = 1, .want = 2, .x = &.{ -1, -2, 0, 0, -2, -1 } },
    .{ .name = "positive_step_2_N_3", .n = 3, .inc = 2, .want = 2, .x = &.{ 1, 2, 0, 0, 1, 1, 0, 0, 2, 1, 0, 0 } },
    .{ .name = "negative_step_2_N_3", .n = 3, .inc = 2, .want = 2, .x = &.{ -1, -2, 0, 0, -1, -1, 0, 0, -2, -1, 0, 0 } },
    .{ .name = "positive_step_1_N_4", .n = 4, .inc = 1, .want = 2, .x = &.{ 1, 2, 0, 0, 2, 1, -2, -2 } },
    .{ .name = "negative_step_1_N_4", .n = 4, .inc = 1, .want = 2, .x = &.{ -1, -2, 0, 0, -2, -1, -2, -2 } },
    .{ .name = "positive_step_2_N_4", .n = 4, .inc = 2, .want = 2, .x = &.{ 1, 2, 0, 0, 1, 1, 0, 0, 2, 1, 0, 0, -2, -2, 0, 0 } },
    .{ .name = "negative_step_2_N_4", .n = 4, .inc = 2, .want = 2, .x = &.{ -1, -2, 0, 0, -1, -1, 0, 0, -2, -1, 0, 0, -2, -2, 0, 0 } },
    .{ .name = "positive_step_1_N_50", .n = 50, .inc = 1, .want = 8, .x = &.{ 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 0, 0, 1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1025, 1026, 1027, 1028, 1029, 1030, 1031, 1032, 1033, 1034, 1035, 1036, 1037, 1038, 1039, 1040, 1041, 1042, 1043, 1044, 1045, 1046, 1047, 1048, 1049, 1050, 1051, 1052, 1053, 1054, 1055, 1056, 1057, 1058, 1059, 1060, 1061, 1062, 1063, 1064, 1065, 1066, 1067, 1068, 1069, 1070, 1071, 1072, 1073, 1074, 1075, 1076, 1077, 1078, 1079, 1080, 1081, 1082, 1083, 1084, 1085, 1086, 1087, 1088, 1089, 1090, 1091, 1092, 1093, 1094, 1095, 1096, 1097, 1098, 1099 } },
    .{ .name = "negative_step_1_N_50", .n = 50, .inc = 1, .want = 8, .x = &.{ -1000, -1001, -1002, -1003, -1004, -1005, -1006, -1007, -1008, -1009, -1010, -1011, -1012, -1013, 0, 0, -1016, -1017, -1018, -1019, -1020, -1021, -1022, -1023, -1024, -1025, -1026, -1027, -1028, -1029, -1030, -1031, -1032, -1033, -1034, -1035, -1036, -1037, -1038, -1039, -1040, -1041, -1042, -1043, -1044, -1045, -1046, -1047, -1048, -1049, -1050, -1051, -1052, -1053, -1054, -1055, -1056, -1057, -1058, -1059, -1060, -1061, -1062, -1063, -1064, -1065, -1066, -1067, -1068, -1069, -1070, -1071, -1072, -1073, -1074, -1075, -1076, -1077, -1078, -1079, -1080, -1081, -1082, -1083, -1084, -1085, -1086, -1087, -1088, -1089, -1090, -1091, -1092, -1093, -1094, -1095, -1096, -1097, -1098, -1099 } },
    .{ .name = "positive_step_2_N_50", .n = 50, .inc = 2, .want = 8, .x = &.{ 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1025, 1026, 1027, 0, 0, 1030, 1031, 1032, 1033, 1034, 1035, 1036, 1037, 1038, 1039, 1040, 1041, 1042, 1043, 1044, 1045, 1046, 1047, 1048, 1049, 1050, 1051, 1052, 1053, 1054, 1055, 1056, 1057, 1058, 1059, 1060, 1061, 1062, 1063, 1064, 1065, 1066, 1067, 1068, 1069, 1070, 1071, 1072, 1073, 1074, 1075, 1076, 1077, 1078, 1079, 1080, 1081, 1082, 1083, 1084, 1085, 1086, 1087, 1088, 1089, 1090, 1091, 1092, 1093, 1094, 1095, 1096, 1097, 1098, 1099, 1100, 1101, 1102, 1103, 1104, 1105, 1106, 1107, 1108, 1109, 1110, 1111, 1112, 1113, 1114, 1115, 1116, 1117, 1118, 1119, 1120, 1121, 1122, 1123, 1124, 1125, 1126, 1127, 1128, 1129, 1130, 1131, 1132, 1133, 1134, 1135, 1136, 1137, 1138, 1139, 1140, 1141, 1142, 1143, 1144, 1145, 1146, 1147, 1148, 1149, 1150, 1151, 1152, 1153, 1154, 1155, 1156, 1157, 1158, 1159, 1160, 1161, 1162, 1163, 1164, 1165, 1166, 1167, 1168, 1169, 1170, 1171, 1172, 1173, 1174, 1175, 1176, 1177, 1178, 1179, 1180, 1181, 1182, 1183, 1184, 1185, 1186, 1187, 1188, 1189, 1190, 1191, 1192, 1193, 1194, 1195, 1196, 1197, 1198, 1199 } },
    .{ .name = "negative_step_2_N_50", .n = 50, .inc = 2, .want = 8, .x = &.{ -1000, -1001, -1002, -1003, -1004, -1005, -1006, -1007, -1008, -1009, -1010, -1011, -1012, -1013, -1014, -1015, -1016, -1017, -1018, -1019, -1020, -1021, -1022, -1023, -1024, -1025, -1026, -1027, 0, 0, -1030, -1031, -1032, -1033, -1034, -1035, -1036, -1037, -1038, -1039, -1040, -1041, -1042, -1043, -1044, -1045, -1046, -1047, -1048, -1049, -1050, -1051, -1052, -1053, -1054, -1055, -1056, -1057, -1058, -1059, -1060, -1061, -1062, -1063, -1064, -1065, -1066, -1067, -1068, -1069, -1070, -1071, -1072, -1073, -1074, -1075, -1076, -1077, -1078, -1079, -1080, -1081, -1082, -1083, -1084, -1085, -1086, -1087, -1088, -1089, -1090, -1091, -1092, -1093, -1094, -1095, -1096, -1097, -1098, -1099, -1100, -1101, -1102, -1103, -1104, -1105, -1106, -1107, -1108, -1109, -1110, -1111, -1112, -1113, -1114, -1115, -1116, -1117, -1118, -1119, -1120, -1121, -1122, -1123, -1124, -1125, -1126, -1127, -1128, -1129, -1130, -1131, -1132, -1133, -1134, -1135, -1136, -1137, -1138, -1139, -1140, -1141, -1142, -1143, -1144, -1145, -1146, -1147, -1148, -1149, -1150, -1151, -1152, -1153, -1154, -1155, -1156, -1157, -1158, -1159, -1160, -1161, -1162, -1163, -1164, -1165, -1166, -1167, -1168, -1169, -1170, -1171, -1172, -1173, -1174, -1175, -1176, -1177, -1178, -1179, -1180, -1181, -1182, -1183, -1184, -1185, -1186, -1187, -1188, -1189, -1190, -1191, -1192, -1193, -1194, -1195, -1196, -1197, -1198, -1199 } },
    .{ .name = "min_idx_in_vec_tail", .n = 50, .inc = 2, .want = 50, .x = &.{ 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1025, 1026, 1027, 1028, 1029, 1030, 1031, 1032, 1033, 1034, 1035, 1036, 1037, 1038, 1039, 1040, 1041, 1042, 1043, 1044, 1045, 1046, 1047, 1048, 1049, 1050, 1051, 1052, 1053, 1054, 1055, 1056, 1057, 1058, 1059, 1060, 1061, 1062, 1063, 1064, 1065, 1066, 1067, 1068, 1069, 1070, 1071, 1072, 1073, 1074, 1075, 1076, 1077, 1078, 1079, 1080, 1081, 1082, 1083, 1084, 1085, 1086, 1087, 1088, 1089, 1090, 1091, 1092, 1093, 1094, 1095, 1096, 1097, 1098, 1099, 1100, 1101, 1102, 1103, 1104, 1105, 1106, 1107, 1108, 1109, 1110, 1111, 1112, 1113, 1114, 1115, 1116, 1117, 1118, 1119, 1120, 1121, 1122, 1123, 1124, 1125, 1126, 1127, 1128, 1129, 1130, 1131, 1132, 1133, 1134, 1135, 1136, 1137, 1138, 1139, 1140, 1141, 1142, 1143, 1144, 1145, 1146, 1147, 1148, 1149, 1150, 1151, 1152, 1153, 1154, 1155, 1156, 1157, 1158, 1159, 1160, 1161, 1162, 1163, 1164, 1165, 1166, 1167, 1168, 1169, 1170, 1171, 1172, 1173, 1174, 1175, 1176, 1177, 1178, 1179, 1180, 1181, 1182, 1183, 1184, 1185, 1186, 1187, 1188, 1189, 1190, 1191, 1192, 1193, 1194, 1195, 0, 0, 1198, 1199 } },
};

/// `test_daxpby.c (test_saxpby.c)`.
const axpby_real = [_]Case{
    .{ .name = "inc_x_1_inc_y_1_N_100", .kind = .run, .n = 100, .incx = 1, .incy = 1, .alpha = .{ 1.0, 0 }, .beta = .{ 1.0, 0 } },
    .{ .name = "inc_x_2_inc_y_1_N_100", .kind = .run, .n = 100, .incx = 2, .incy = 1, .alpha = .{ 2.0, 0 }, .beta = .{ 1.0, 0 } },
    .{ .name = "inc_x_1_inc_y_2_N_100", .kind = .run, .n = 100, .incx = 1, .incy = 2, .alpha = .{ 1.0, 0 }, .beta = .{ 2.0, 0 } },
    .{ .name = "inc_x_2_inc_y_2_N_100", .kind = .run, .n = 100, .incx = 2, .incy = 2, .alpha = .{ 3.0, 0 }, .beta = .{ 4.0, 0 } },
    .{ .name = "inc_x_neg_1_inc_y_2_N_100", .kind = .run, .n = 100, .incx = -1, .incy = 2, .alpha = .{ 5.0, 0 }, .beta = .{ 4.0, 0 } },
    .{ .name = "inc_x_2_inc_y_neg_1_N_100", .kind = .run, .n = 100, .incx = 2, .incy = -1, .alpha = .{ 1.0, 0 }, .beta = .{ 6.0, 0 } },
    .{ .name = "inc_x_neg_2_inc_y_neg_1_N_100", .kind = .run, .n = 100, .incx = -2, .incy = -1, .alpha = .{ 7.0, 0 }, .beta = .{ 3.5, 0 } },
    .{ .name = "inc_x_1_inc_y_1_N_100_alpha_zero", .kind = .run, .n = 100, .incx = 1, .incy = 1, .alpha = .{ 0.0, 0 }, .beta = .{ 1.0, 0 } },
    .{ .name = "inc_x_1_inc_y_2_N_100_alpha_zero", .kind = .run, .n = 100, .incx = 1, .incy = 2, .alpha = .{ 0.0, 0 }, .beta = .{ 1.0, 0 } },
    .{ .name = "inc_x_1_inc_y_1_N_100_beta_zero", .kind = .run, .n = 100, .incx = 1, .incy = 1, .alpha = .{ 1.0, 0 }, .beta = .{ 0.0, 0 } },
    .{ .name = "inc_x_2_inc_y_1_N_100_beta_zero", .kind = .run, .n = 100, .incx = 2, .incy = 1, .alpha = .{ 1.0, 0 }, .beta = .{ 0.0, 0 } },
    .{ .name = "inc_x_1_inc_y_2_N_100_beta_zero", .kind = .run, .n = 100, .incx = 1, .incy = 2, .alpha = .{ 1.0, 0 }, .beta = .{ 0.0, 0 } },
    .{ .name = "inc_x_2_inc_y_2_N_100_beta_zero", .kind = .run, .n = 100, .incx = 2, .incy = 2, .alpha = .{ 1.0, 0 }, .beta = .{ 0.0, 0 } },
    .{ .name = "inc_x_1_inc_y_1_N_100_alpha_beta_zero", .kind = .run, .n = 100, .incx = 1, .incy = 1, .alpha = .{ 0.0, 0 }, .beta = .{ 0.0, 0 } },
    .{ .name = "inc_x_1_inc_y_2_N_100_alpha_beta_zero", .kind = .run, .n = 100, .incx = 1, .incy = 2, .alpha = .{ 0.0, 0 }, .beta = .{ 0.0, 0 } },
    .{ .name = "check_n_zero", .kind = .run, .n = 0, .incx = 1, .incy = 1, .alpha = .{ 1.0, 0 }, .beta = .{ 1.0, 0 } },
};

/// `test_zaxpby.c (test_caxpby.c)`.
const axpby_complex = [_]Case{
    .{ .name = "inc_x_1_inc_y_1_N_100", .kind = .run, .n = 100, .incx = 1, .incy = 1, .alpha = .{ 1.0, 1.0 }, .beta = .{ 1.0, 1.0 } },
    .{ .name = "inc_x_2_inc_y_1_N_100", .kind = .run, .n = 100, .incx = 2, .incy = 1, .alpha = .{ 2.0, 1.0 }, .beta = .{ 1.0, 1.0 } },
    .{ .name = "inc_x_1_inc_y_2_N_100", .kind = .run, .n = 100, .incx = 1, .incy = 2, .alpha = .{ 1.0, 1.0 }, .beta = .{ 2.0, 1.0 } },
    .{ .name = "inc_x_2_inc_y_2_N_100", .kind = .run, .n = 100, .incx = 2, .incy = 2, .alpha = .{ 3.0, 1.0 }, .beta = .{ 4.0, 3.0 } },
    .{ .name = "inc_x_neg_1_inc_y_2_N_100", .kind = .run, .n = 100, .incx = -1, .incy = 2, .alpha = .{ 5.0, 2.2 }, .beta = .{ 4.0, 5.0 } },
    .{ .name = "inc_x_2_inc_y_neg_1_N_100", .kind = .run, .n = 100, .incx = 2, .incy = -1, .alpha = .{ 1.0, 1.0 }, .beta = .{ 6.0, 3.0 } },
    .{ .name = "inc_x_neg_2_inc_y_neg_1_N_100", .kind = .run, .n = 100, .incx = -2, .incy = -1, .alpha = .{ 7.0, 2.0 }, .beta = .{ 3.5, 1.3 } },
    .{ .name = "inc_x_1_inc_y_1_N_100_alpha_zero", .kind = .run, .n = 100, .incx = 1, .incy = 1, .alpha = .{ 0.0, 0.0 }, .beta = .{ 1.0, 1.0 } },
    .{ .name = "inc_x_1_inc_y_1_N_100_beta_zero", .kind = .run, .n = 100, .incx = 1, .incy = 1, .alpha = .{ 1.0, 1.0 }, .beta = .{ 0.0, 0.0 } },
    .{ .name = "inc_x_1_inc_y_1_N_100_alpha_beta_zero", .kind = .run, .n = 100, .incx = 1, .incy = 1, .alpha = .{ 0.0, 0.0 }, .beta = .{ 0.0, 0.0 } },
    .{ .name = "inc_x_1_inc_y_2_N_100_alpha_beta_zero", .kind = .run, .n = 100, .incx = 1, .incy = 2, .alpha = .{ 0.0, 0.0 }, .beta = .{ 0.0, 0.0 } },
    .{ .name = "check_n_zero", .kind = .run, .n = 0, .incx = 1, .incy = 1, .alpha = .{ 1.0, 1.0 }, .beta = .{ 1.0, 1.0 } },
};

/// `test_zaxpyc.c (test_caxpyc.c)`.
const axpyc_cases = [_]Case{
    .{ .name = "conj_strides_one", .kind = .run, .n = 100, .incx = 1, .incy = 1, .alpha = .{ 5.0, 2.2 } },
    .{ .name = "conj_incx_one", .kind = .run, .n = 100, .incx = 1, .incy = 2, .alpha = .{ 5.0, 2.2 } },
    .{ .name = "conj_incy_one", .kind = .run, .n = 100, .incx = 2, .incy = 1, .alpha = .{ 5.0, 2.2 } },
    .{ .name = "conj_strides_two", .kind = .run, .n = 100, .incx = 2, .incy = 2, .alpha = .{ 5.0, 2.2 } },
};

/// `test_zrot.c (test_crot.c); `inc_x_0_inc_y_0` is disabled there`.
const rot_cases = [_]Case{
    .{ .name = "inc_x_1_inc_y_1", .kind = .run, .n = 100, .inc_x = 1, .inc_y = 1, .c = .{ 1.0, 0.0 }, .s = .{ 1.0, 0.0 } },
    .{ .name = "inc_x_neg_1_inc_y_neg_1", .kind = .run, .n = 100, .inc_x = -1, .inc_y = -1, .c = .{ 1.0, 0.0 }, .s = .{ 1.0, 0.0 } },
    .{ .name = "inc_x_2_inc_y_1", .kind = .run, .n = 100, .inc_x = 2, .inc_y = 1, .c = .{ 3.0, 0.0 }, .s = .{ 2.0, 0.0 } },
    .{ .name = "inc_x_neg_2_inc_y_1", .kind = .run, .n = 100, .inc_x = -2, .inc_y = 1, .c = .{ 1.0, 0.0 }, .s = .{ 1.0, 0.0 } },
    .{ .name = "inc_x_1_inc_y_2", .kind = .run, .n = 100, .inc_x = 1, .inc_y = 2, .c = .{ 1.0, 0.0 }, .s = .{ 1.0, 0.0 } },
    .{ .name = "inc_x_1_inc_y_neg_2", .kind = .run, .n = 100, .inc_x = 1, .inc_y = -2, .c = .{ 2.0, 0.0 }, .s = .{ 1.0, 0.0 } },
    .{ .name = "inc_x_2_inc_y_2", .kind = .run, .n = 100, .inc_x = 2, .inc_y = 2, .c = .{ 1.0, 0.0 }, .s = .{ 2.0, 0.0 } },
    .{ .name = "inc_x_neg_2_inc_y_neg_2", .kind = .run, .n = 100, .inc_x = -2, .inc_y = -2, .c = .{ 1.0, 0.0 }, .s = .{ 1.0, 0.0 } },
    .{ .name = "inc_x_2_inc_y_2_c_zero", .kind = .run, .n = 100, .inc_x = 2, .inc_y = 2, .c = .{ 0.0, 0.0 }, .s = .{ 1.0, 0.0 } },
    .{ .name = "inc_x_2_inc_y_2_s_zero", .kind = .run, .n = 100, .inc_x = 2, .inc_y = 2, .c = .{ 1.0, 0.0 }, .s = .{ 0.0, 0.0 } },
    .{ .name = "check_n_zero", .kind = .run, .n = 0, .inc_x = 1, .inc_y = 1, .c = .{ 1.0, 0.0 }, .s = .{ 1.0, 0.0 } },
};

/// `test_zscal.c (test_cscal.c)`.
const scal_cases = [_]Case{
    .{ .name = "alpha_r_zero_alpha_i_not_zero", .kind = .run, .N = 100, .inc = 1, .alpha = .{ 0.0, 1.0 } },
    .{ .name = "alpha_r_zero_alpha_i_zero_inc_2", .kind = .run, .N = 100, .inc = 2, .alpha = .{ 0.0, 0.0 } },
};

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
