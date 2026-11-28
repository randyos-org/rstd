//! Tests: `gesvx` by Reference-LAPACK 3.12.0's `TESTING/LIN/ddrvge.f`: every
//! `fact`, `trans` and given equilibration, with `dget01.f`'s and `dget02.f`'s
//! residuals, `gerfs`'s `dget04.f` and `dget07.f` ratios, `dget06.f`'s on
//! `rcond`, the pivot growth, and the zero columns of its types 5 to 7.

const std = @import("std");
const types = @import("../../blas/types.zig");
const check = @import("../../blas/_test/check.zig");
const lcheck = @import("check.zig");
const refine_tests = @import("gerfs.zig");
const ops = @import("../ops.zig");
const interface = @import("../interface.zig");
const norm = @import("../norm.zig");
const lu = @import("../lu.zig");
const equilibrate = @import("../equilibrate.zig");
const expert = @import("../gesvx.zig");

const testing = std.testing;
const Allocator = std.mem.Allocator;
const Float = types.Float;
const Trans = types.Trans;
const Equilibration = interface.Equilibration;
const FactorizationInput = interface.FactorizationInput;
const abs1 = ops.abs1;
const opAt = refine_tests.opAt;
const thresh = lcheck.thresh;

/// The matrices: `dlatb4.f`'s well-conditioned type 1, rows and columns
/// scaled by powers of 2 for `geequ` to undo, two nearly equal rows (`rcond`
/// below `eps`), and types 5 to 7's zero first, last and trailing columns.
const Kind = enum { plain, scaled, near_singular, zero_first, zero_last, zero_trailing };

/// What one `gesvx` call tries: `fact`, and for `.factored` the scaling the
/// given factor was made with, as `ddrvge.f`'s `EQUEDS`.
const Variant = struct { fact: FactorizationInput, equed: Equilibration = .none };

const variants = [_]Variant{
    .{ .fact = .factor },
    .{ .fact = .equilibrate_and_factor },
    .{ .fact = .factored },
    .{ .fact = .factored, .equed = .rows },
    .{ .fact = .factored, .equed = .columns },
    .{ .fact = .factored, .equed = .both },
};

fn matrix(comptime T: type, kind: Kind, n: usize, a: []T) ?usize {
    check.fill(T, a, n * 31 + 5);
    for (0..n) |i| a[i + i * n] = ops.add(T, a[i + i * n], ops.scalar(T, 2));
    switch (kind) {
        .plain => {},
        // rows by 2^-3 to 2^3, columns by 2^-3, 1 and 2^3
        .scaled => for (0..n) |j| for (0..n) |i| {
            const e: i32 = @as(i32, @intCast(i % 4)) * 2 - 3 + @as(i32, @intCast(j % 3)) * 3 - 3;
            a[i + j * n] = ops.mulReal(T, a[i + j * n], std.math.ldexp(@as(Float(T), 1), e));
        },
        // row 1 is row 0 but for `2^-18 eps` in the last column, whose largest
        // entry is row 2's, so scaling cannot undo it; column 0 is zero below
        // them, so row 0 is the first pivot and `U`'s last entry that difference
        .near_singular => if (n > 2) {
            const eps = std.math.floatEps(Float(T));
            const s = std.math.ldexp(@as(Float(T), 1), -20);
            const l = (n - 1) * n;
            @memset(a[2..n], ops.zero(T));
            a[0] = ops.scalar(T, 8);
            a[2 + l] = ops.scalar(T, 2);
            for (0..n) |j| a[1 + j * n] = a[j * n];
            a[l] = ops.scalar(T, s);
            a[1 + l] = ops.scalar(T, s * (1 + 4 * eps));
        },
        .zero_first, .zero_last, .zero_trailing => {
            const z = switch (kind) {
                .zero_first => 0,
                .zero_last => n - 1,
                else => n / 2,
            };
            const last = if (kind == .zero_trailing) n else z + 1;
            @memset(a[z * n .. last * n], ops.zero(T));
            return z;
        },
    }
    return null;
}

/// `dget01.f`: `‖P L U - A‖₁ / (n ‖A‖₁ eps)`, `eps` the unit roundoff
/// (`DLAMCH('Epsilon')`, half of `floatEps`), as in the ratios below.
fn factorRatio(comptime T: type, al: Allocator, n: usize, a: []const T, af: []const T, ipiv: []const usize) !Float(T) {
    const F = Float(T);
    const lu_ = try al.alloc(T, n * n);
    defer al.free(lu_);
    for (0..n) |j| for (0..n) |i| {
        var s = ops.zero(T);
        for (0..@min(i, j) + 1) |k| {
            const l = if (k == i) ops.one(T) else af[i + k * n];
            s = ops.add(T, s, ops.mul(T, l, af[k + j * n]));
        }
        lu_[i + j * n] = s;
    };
    var k = n;
    while (k > 0) {
        k -= 1;
        if (ipiv[k] != k) for (0..n) |j| std.mem.swap(T, &lu_[k + j * n], &lu_[ipiv[k] + j * n]);
    }
    for (lu_, a) |*e, v| e.* = ops.sub(T, e.*, v);
    const anorm = try norm.lange(al, check.io, T, .one, n, n, a, n);
    const resid = try norm.lange(al, check.io, T, .one, n, n, lu_, n);
    const eps = std.math.floatEps(F) / 2;
    if (anorm == 0) return if (resid == 0) 0 else 1 / eps;
    return resid / @as(F, @floatFromInt(n)) / anorm / eps;
}

/// `dget02.f`: the largest `‖B - op(A) X‖₁ / (‖op(A)‖₁ ‖X‖₁ eps)` over the
/// columns.
fn residualRatio(comptime T: type, al: Allocator, trans: Trans, n: usize, nrhs: usize, a: []const T, b: []const T, x: []const T) !Float(T) {
    const F = Float(T);
    const eps = std.math.floatEps(F) / 2;
    const anorm = try norm.lange(al, check.io, T, if (trans == .n) .one else .infinity, n, n, a, n);
    var out: F = 0;
    for (0..nrhs) |j| {
        var rnorm: F = 0;
        var xnorm: F = 0;
        for (0..n) |i| {
            var t = b[i + j * n];
            for (0..n) |k| t = ops.sub(T, t, ops.mul(T, opAt(T, trans, a, n, i, k), x[k + j * n]));
            rnorm += abs1(T, t);
            xnorm += abs1(T, x[i + j * n]);
        }
        out = @max(out, if (xnorm == 0) 1 / eps else rnorm / anorm / xnorm / eps);
    }
    return out;
}

/// `1 / (‖A‖ ‖A⁻¹‖)` in `kind`, from the inverse, as `ddrvge.f`'s `RCONDC`.
fn trueRcond(comptime T: type, al: Allocator, kind: norm.NormType, n: usize, a: []const T) !Float(T) {
    const inv = try al.dupe(T, a);
    defer al.free(inv);
    const ipiv = try al.alloc(usize, n);
    defer al.free(ipiv);
    _ = try lu.getrf(al, check.io, T, n, n, inv, n, ipiv);
    try lu.getri(al, check.io, T, n, inv, n, ipiv);
    return 1 / (try norm.lange(al, check.io, T, kind, n, n, a, n) * try norm.lange(al, check.io, T, kind, n, n, inv, n));
}

fn fail(comptime T: type, v: Variant, trans: Trans, kind: Kind, n: usize, what: []const u8, value: anytype) error{TestUnexpectedResult} {
    check.say("gesvx {s} {t} {t} trans {t} {t} n {d}: {s} {any}\n", .{ @typeName(T), v.fact, v.equed, trans, kind, n, what, value });
    return error.TestUnexpectedResult;
}

fn case(comptime T: type, v: Variant, trans: Trans, kind: Kind, n: usize) !void {
    const F = Float(T);
    const al = testing.allocator;
    const nrhs = 2;
    const asav = try al.alloc(T, n * n);
    defer al.free(asav);
    const zero_at = matrix(T, kind, n, asav);
    if (zero_at != null and v.fact == .factored) return;
    const xact = try al.alloc(T, n * nrhs);
    defer al.free(xact);
    check.fill(T, xact, n * 7 + 3);
    const bsav = try al.alloc(T, n * nrhs);
    defer al.free(bsav);
    for (0..nrhs) |j| for (0..n) |i| {
        var t = ops.zero(T);
        for (0..n) |k| t = ops.add(T, t, ops.mul(T, opAt(T, trans, asav, n, i, k), xact[k + j * n]));
        bsav[i + j * n] = t;
    };
    const a = try al.dupe(T, asav);
    defer al.free(a);
    const af = try al.alloc(T, n * n);
    defer al.free(af);
    @memset(af, ops.zero(T));
    const ipiv = try al.alloc(usize, n);
    defer al.free(ipiv);
    const s = try al.alloc(F, 2 * n);
    defer al.free(s);
    @memset(s, 1);
    var equed: Equilibration = .none;
    if (v.fact == .factored) {
        // `ddrvge.f`: `laqge` forced to the scaling under test, then `getrf`
        if (v.equed != .none) {
            var sc = try equilibrate.geequ(T, n, n, a, n, s[0..n], s[n..]);
            sc.rowcnd = if (v.equed == .columns) 1 else 0;
            sc.colcnd = if (v.equed == .rows) 1 else 0;
            equed = try equilibrate.laqge(T, n, n, a, n, s[0..n], s[n..], sc);
            if (equed != v.equed) return fail(T, v, trans, kind, n, "laqge scaled", equed);
        }
        @memcpy(af, a);
        try testing.expectEqual(null, try lu.getrf(al, check.io, T, n, n, af, n, ipiv));
    }
    const b = try al.dupe(T, bsav);
    defer al.free(b);
    const x = try al.alloc(T, n * nrhs);
    defer al.free(x);
    @memset(x, ops.zero(T));
    var ferr: [nrhs]F = undefined;
    var berr: [nrhs]F = undefined;
    const got = try expert.gesvx(al, check.io, T, v.fact, trans, n, nrhs, a, n, af, n, ipiv, &equed, s[0..n], s[n..], b, n, x, n, &ferr, &berr);
    if (!std.meta.eql(got.singular, zero_at)) return fail(T, v, trans, kind, n, "singular at", got.singular);

    // the pivot growth, recomputed (`ddrvge.f`'s test 7)
    const cols = if (zero_at) |z| z + 1 else n;
    var amax: F = 0;
    var umax: F = 0;
    for (0..cols) |j| for (0..n) |i| {
        amax = @max(amax, ops.abs(T, a[i + j * n]));
        if (i <= j and i < cols) umax = @max(umax, ops.abs(T, af[i + j * n]));
    };
    const growth = if (umax == 0) 1 else amax / umax;
    if (!(@abs(growth - got.pivot_growth) / @max(growth, got.pivot_growth) / (std.math.floatEps(F) / 2) < thresh)) return fail(T, v, trans, kind, n, "pivot growth", got.pivot_growth);
    if (v.fact != .factored) {
        const r1 = try factorRatio(T, al, n, a, af, ipiv);
        if (!(r1 < thresh)) return fail(T, v, trans, kind, n, "dget01", r1);
    }
    if (zero_at != null) {
        if (got.rcond != 0) return fail(T, v, trans, kind, n, "rcond", got.rcond);
        return;
    }
    if (kind == .near_singular and n > 2) {
        if (!got.ill_conditioned) return fail(T, v, trans, kind, n, "ill_conditioned", got.rcond);
        return;
    }
    if (got.ill_conditioned) return fail(T, v, trans, kind, n, "ill_conditioned", got.rcond);
    const r2 = try residualRatio(T, al, trans, n, nrhs, asav, bsav, x);
    if (!(r2 < thresh)) return fail(T, v, trans, kind, n, "dget02", r2);
    const r = try refine_tests.ratios(T, al, trans, n, nrhs, asav, bsav, x, xact, &ferr, &berr);
    for (r) |e| if (!(e < thresh)) return fail(T, v, trans, kind, n, "dget07 and dget04", r);
    // `dget06` against `A` as factored, equilibrated
    const rc = try trueRcond(T, al, if (trans == .n) .one else .infinity, n, a);
    const r6 = lcheck.ratio(F, got.rcond, rc);
    if (!(r6 < thresh)) return fail(T, v, trans, kind, n, "dget06", r6);
    // `B` comes back scaled as `A` was, `diag(r)` or `diag(c)` by `trans`
    const sb = if (trans == .n) equed == .rows or equed == .both else equed == .columns or equed == .both;
    for (0..nrhs) |j| for (0..n) |i| {
        const f = if (!sb) 1 else if (trans == .n) s[i] else s[n + i];
        if (!lcheck.sameBits(T, b[i + j * n], ops.mulReal(T, bsav[i + j * n], f))) return fail(T, v, trans, kind, n, "B scaled at", .{ i, j });
    };
}

test "gesvx: ddrvge's ratios, every fact, trans and given scaling, every type" {
    inline for (lcheck.test_types) |T| for (variants) |v| for ([_]Trans{ .n, .t, .c }) |trans| for (std.meta.tags(Kind)) |kind| for ([_]usize{ 1, 2, 5, 23 }) |n| {
        try case(T, v, trans, kind, n);
    };
}

test "gesvx: equilibrate_and_factor scales a badly scaled A, factor leaves it" {
    const al = testing.allocator;
    const n = 5;
    var a: [n * n]f64 = undefined;
    _ = matrix(f64, .scaled, n, &a);
    var af: [n * n]f64 = undefined;
    var ipiv: [n]usize = undefined;
    var s: [2 * n]f64 = undefined;
    var b: [n]f64 = @splat(1);
    var x: [n]f64 = undefined;
    var ferr: [1]f64 = undefined;
    var berr: [1]f64 = undefined;
    var equed: Equilibration = .both;
    const a0 = a;
    _ = try expert.gesvx(al, check.io, f64, .factor, .n, n, 1, &a, n, &af, n, &ipiv, &equed, s[0..n], s[n..], &b, n, &x, n, &ferr, &berr);
    try testing.expectEqual(Equilibration.none, equed);
    try testing.expectEqualSlices(f64, &a0, &a);
    _ = try expert.gesvx(al, check.io, f64, .equilibrate_and_factor, .n, n, 1, &a, n, &af, n, &ipiv, &equed, s[0..n], s[n..], &b, n, &x, n, &ferr, &berr);
    try testing.expectEqual(Equilibration.both, equed);
}

test "gesvx: empty sizes, bad arguments" {
    const al = testing.allocator;
    var a = [_]f64{ 2, 1, 1, 3 };
    var af = a;
    var ipiv = [_]usize{ 0, 1 };
    var s = [_]f64{ 1, 1, 1, 1 };
    var b = [_]f64{ 1, 2 };
    var x = [_]f64{ 1, 2 };
    var ferr = [_]f64{ 7, 7 };
    var berr = [_]f64{ 7, 7 };
    var equed: Equilibration = .none;
    const e0 = try expert.gesvx(al, check.io, f64, .factor, .n, 0, 2, a[0..0], 1, af[0..0], 1, ipiv[0..0], &equed, s[0..0], s[0..0], b[0..0], 1, x[0..0], 1, &ferr, &berr);
    try testing.expectEqual(@as(f64, 1), e0.rcond);
    try testing.expectEqual(@as(f64, 1), e0.pivot_growth);
    try testing.expectEqualSlices(f64, &.{ 0, 0 }, &ferr);
    const e = testing.expectError;
    try e(error.UnsupportedTranspose, expert.gesvx(al, check.io, f64, .factor, .r, 2, 1, &a, 2, &af, 2, &ipiv, &equed, s[0..2], s[2..], &b, 2, &x, 2, &ferr, &berr));
    try e(error.LeadingDimensionTooSmall, expert.gesvx(al, check.io, f64, .factor, .n, 2, 1, &a, 2, &af, 1, &ipiv, &equed, s[0..2], s[2..], &b, 2, &x, 2, &ferr, &berr));
    try e(error.SliceTooShort, expert.gesvx(al, check.io, f64, .factor, .n, 2, 1, &a, 2, &af, 2, &ipiv, &equed, s[0..1], s[2..], &b, 2, &x, 2, &ferr, &berr));
    var bad = [_]usize{ 0, 2 };
    try e(error.PivotOutOfRange, expert.gesvx(al, check.io, f64, .factored, .n, 2, 1, &a, 2, &af, 2, &bad, &equed, s[0..2], s[2..], &b, 2, &x, 2, &ferr, &berr));
    equed = .rows;
    s[1] = 0;
    try e(error.ScaleFactorNotPositive, expert.gesvx(al, check.io, f64, .factored, .n, 2, 1, &a, 2, &af, 2, &ipiv, &equed, s[0..2], s[2..], &b, 2, &x, 2, &ferr, &berr));
    equed = .columns;
    s[1] = 1;
    s[3] = -1;
    try e(error.ScaleFactorNotPositive, expert.gesvx(al, check.io, f64, .factored, .n, 2, 1, &a, 2, &af, 2, &ipiv, &equed, s[0..2], s[2..], &b, 2, &x, 2, &ferr, &berr));
}
