//! Tests: `ormtr` by Reference-LAPACK 3.12.0's `TESTING/EIG/dsyt21.f`
//! ratios on `sytrd`'s reduction, `‖Qᴴ A Q - T‖₁ / (n ‖A‖₁ ulp)` and
//! `‖I - Qᴴ Q‖₁ / (n ulp)` under `sep.in`'s threshold (60), with `Q` formed
//! by `ormtr` from the identity, and every side and `trans` against it.

const std = @import("std");
const types = @import("../../blas/types.zig");
const check = @import("../../blas/_test/check.zig");
const lcheck = @import("check.zig");
const ops = @import("../ops.zig");
const norm = @import("../norm.zig");
const reduce = @import("../sytrd.zig");
const apply = @import("../ormtr.zig");

const testing = std.testing;
const Float = types.Float;
const Mat = ops.Mat;
const Side = types.Side;
const Trans = types.Trans;
const Uplo = types.Uplo;

/// `sep.in`'s `THRESH`.
const thresh = 60;

fn fail(comptime T: type, uplo: Uplo, n: usize, what: []const u8, value: anytype) error{TestUnexpectedResult} {
    check.say("ormtr {s} uplo {t} n {d}: {s} {any}\n", .{ @typeName(T), uplo, n, what, value });
    return error.TestUnexpectedResult;
}

fn oneNorm(comptime T: type, m: Mat(T)) !Float(T) {
    return norm.lange(testing.allocator, check.io, T, .one, m.rows, m.cols, m.slice(), m.ld);
}

fn case(comptime T: type, uplo: Uplo, n: usize) !void {
    const F = Float(T);
    const al = testing.allocator;
    const a = try ops.alloc(T, al, n, n);
    defer ops.free(T, al, a);
    lcheck.hermitian(T, a, n * 29 + 7);
    const af = try ops.alloc(T, al, n, n);
    defer ops.free(T, al, af);
    af.copyFrom(a);
    const d = try al.alloc(F, n);
    defer al.free(d);
    const e = try al.alloc(F, n -| 1);
    defer al.free(e);
    const tau = try al.alloc(T, n -| 1);
    defer al.free(tau);
    try reduce.sytrd(al, check.io, T, uplo, n, af.slice(), af.ld, d, e, tau);

    const q = try ops.alloc(T, al, n, n);
    defer ops.free(T, al, q);
    q.identity();
    try apply.ormtr(al, check.io, T, .l, uplo, .n, n, n, af.slice(), af.ld, tau, q.slice(), q.ld);
    const ulp = std.math.floatEps(F);
    const nf: F = @floatFromInt(n);
    // `Qᴴ A Q - T`
    const aq = try ops.alloc(T, al, n, n);
    defer ops.free(T, al, aq);
    lcheck.mul(T, false, a, false, q, aq);
    const r = try ops.alloc(T, al, n, n);
    defer ops.free(T, al, r);
    lcheck.mul(T, true, q, false, aq, r);
    for (0..n) |i| {
        r.at(i, i).* = ops.sub(T, r.get(i, i), ops.scalar(T, d[i]));
        if (i + 1 < n) {
            r.at(i + 1, i).* = ops.sub(T, r.get(i + 1, i), ops.scalar(T, e[i]));
            r.at(i, i + 1).* = ops.sub(T, r.get(i, i + 1), ops.scalar(T, e[i]));
        }
    }
    const anorm = @max(try oneNorm(T, a), std.math.floatMin(F));
    const r1 = try oneNorm(T, r) / (nf * anorm * ulp);
    if (!(r1 < thresh)) return fail(T, uplo, n, "Qᴴ A Q - T", r1);
    // `I - Qᴴ Q`
    lcheck.mul(T, true, q, false, q, r);
    for (0..n) |i| r.at(i, i).* = ops.sub(T, r.get(i, i), ops.one(T));
    const r2 = try oneNorm(T, r) / (nf * ulp);
    if (!(r2 < thresh)) return fail(T, uplo, n, "I - Qᴴ Q", r2);

    // each side and `trans` on a random `C`, against `q`
    const p = 3;
    const adj: Trans = if (comptime types.isComplex(T)) .c else .t;
    inline for (.{ Side.l, Side.r }) |side| for ([_]Trans{ .n, adj }) |trans| {
        const rows = if (side == .l) n else p;
        const cols = if (side == .l) p else n;
        const c = try ops.alloc(T, al, rows, cols);
        defer ops.free(T, al, c);
        check.fill(T, c.slice(), n * 3 + rows);
        const want = try ops.alloc(T, al, rows, cols);
        defer ops.free(T, al, want);
        if (side == .l) lcheck.mul(T, trans != .n, q, false, c, want) else lcheck.mul(T, false, c, trans != .n, q, want);
        try apply.ormtr(al, check.io, T, side, uplo, trans, rows, cols, af.slice(), af.ld, tau, c.slice(), c.ld);
        try lcheck.expectNear(T, "ormtr against Q", c, want, thresh * nf * ulp * 4);
    };
}

test "ormtr: dsyt21's ratios on sytrd's Q, every side, trans and uplo, every type" {
    inline for (lcheck.test_types) |T| for ([_]Uplo{ .u, .l }) |uplo| for ([_]usize{ 1, 2, 5, 23, 40 }) |n| {
        try case(T, uplo, n);
    };
}

test "ormtr: empty sizes, bad arguments" {
    const al = testing.allocator;
    var a = [_]f64{ 2, 1, 1, 3 };
    const tau = [_]f64{0};
    var c = [_]f64{ 1, 2, 3, 4 };
    try apply.ormtr(al, check.io, f64, .l, .u, .n, 0, 2, a[0..0], 1, tau[0..0], c[0..0], 1);
    const ee = testing.expectError;
    try ee(error.UnsupportedUplo, apply.ormtr(al, check.io, f64, .l, .v, .n, 2, 2, &a, 2, &tau, &c, 2));
    try ee(error.UnsupportedTranspose, apply.ormtr(al, check.io, f64, .l, .u, .r, 2, 2, &a, 2, &tau, &c, 2));
    var z = [_]types.Complex64{ .init(1, 0), .init(0, 0), .init(0, 0), .init(1, 0) };
    try ee(error.UnsupportedTranspose, apply.ormtr(al, check.io, types.Complex64, .l, .l, .t, 2, 2, &z, 2, z[0..1], &z, 2));
    try ee(error.LeadingDimensionTooSmall, apply.ormtr(al, check.io, f64, .l, .l, .n, 2, 2, &a, 1, &tau, &c, 2));
    try ee(error.SliceTooShort, apply.ormtr(al, check.io, f64, .r, .l, .n, 2, 2, &a, 2, tau[0..0], &c, 2));
}
