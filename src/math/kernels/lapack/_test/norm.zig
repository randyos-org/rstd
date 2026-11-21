//! Tests of `../norm.zig`, `../lacn2.zig`, `../latrs.zig` and `../gecon.zig`.

const std = @import("std");
const types = @import("../../blas/types.zig");
const ops = @import("../ops.zig");
const auxiliary = @import("../auxiliary.zig");
const norm = @import("../norm.zig");
const gecon_estimate = @import("../gecon.zig");
const lacn2_estimate = @import("../lacn2.zig");
const latrs_solve = @import("../latrs.zig");

const Allocator = std.mem.Allocator;
const Float = types.Float;
const isComplex = types.isComplex;
const Trans = types.Trans;
const Norm = norm.Norm;
const gecon = gecon_estimate.gecon;
const lacn2 = lacn2_estimate.lacn2;
const lange = norm.lange;
const lassq = norm.lassq;
const latrs = latrs_solve.latrs;
const rscl = latrs_solve.rscl;

const testing = std.testing;
const check = @import("../../blas/_test/check.zig");
const lcheck = @import("check.zig");
const lu = @import("../lu.zig");

const test_types = lcheck.test_types;

test "lapack norm: lange's four norms on a known 2 x 3, NaN, empty, bad arguments" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        // lda 3, columns (1, -2), (3 + 4i, 0), (-6, 0); complex (3, 4) has modulus 5
        const a = [_]T{ lcheck.cplx(T, 1, 0), lcheck.cplx(T, -2, 0), lcheck.cplx(T, 9, 9), lcheck.cplx(T, 3, if (isComplex(T)) 4 else 0), lcheck.cplx(T, 0, 0), lcheck.cplx(T, 9, 9), lcheck.cplx(T, -6, 0), lcheck.cplx(T, 0, 0) };
        const c2: F = if (comptime isComplex(T)) 5 else 3;
        try testing.expectEqual(@as(F, 6), try lange(al, check.io, T, .m, 2, 3, &a, 3));
        try testing.expectEqual(@as(F, 6), try lange(al, check.io, T, .o, 2, 3, &a, 3));
        try testing.expectEqual(1 + c2 + 6, try lange(al, check.io, T, .i, 2, 3, &a, 3));
        const fro = try lange(al, check.io, T, .f, 2, 3, &a, 3);
        try testing.expectApproxEqRel(@sqrt(41 + c2 * c2), fro, 4 * std.math.floatEps(F));
        var b = a;
        b[1] = lcheck.poisoned(T);
        for ([_]Norm{ .m, .o, .i, .f }) |nm| try testing.expect(std.math.isNan(try lange(al, check.io, T, nm, 2, 3, &b, 3)));
        try testing.expectEqual(@as(F, 0), try lange(testing.failing_allocator, check.io, T, .i, 0, 3, &b, 1));
        try testing.expectEqual(@as(F, 0), try lange(testing.failing_allocator, check.io, T, .f, 2, 0, &b, 2));
    }
    var x = [_]f64{ 0, 0, 0 };
    try testing.expectError(error.LeadingDimensionTooSmall, lange(al, check.io, f64, .m, 2, 1, &x, 1));
    try testing.expectError(error.SliceTooShort, lange(al, check.io, f64, .m, 2, 2, &x, 2));
}

test "lapack norm: lange's Frobenius norm neither overflows nor underflows (Blue's lassq)" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        // sqrt(20) big stays finite
        const big = std.math.floatMax(F) / 16;
        const tiny = std.math.floatTrueMin(F) * 16;
        for ([_]F{ big, tiny, 1 }) |s| {
            var a: [8]T = undefined;
            for (&a, 0..) |*v, i| v.* = ops.scalar(T, s * @as(F, @floatFromInt(i % 2 + 1)));
            // 4 (1 + 4) s squared
            const want = @sqrt(@as(F, 20)) * s;
            const got = try lange(al, check.io, T, .f, 4, 2, &a, 4);
            try testing.expect(std.math.isFinite(got) and got > 0);
            // denormal entries carry few bits
            const tol: F = if (s == tiny) 0.15 else 8 * std.math.floatEps(F);
            try testing.expectApproxEqRel(want, got, tol);
        }
        // a big and a small entry together
        var c = [_]T{ ops.scalar(T, big), ops.scalar(T, tiny), ops.scalar(T, 1) };
        try testing.expectEqual(big, try lange(al, check.io, T, .f, 3, 1, &c, 3));
    }
}

/// `gecon`'s estimate against `1 / (‖A‖ ‖A⁻¹‖)` from `getrs`'s inverse,
/// as the reference's `dchkge.f` checks it, in both norms.
fn geconCase(comptime T: type, al: Allocator, n: usize, grade: Float(T), seed: u64) !void {
    const F = Float(T);
    const a = try al.alloc(T, n * n);
    defer al.free(a);
    check.fill(T, a, seed);
    // row i scaled by grade^i
    for (0..n) |i| {
        const s: F = @floatCast(std.math.pow(f64, @floatCast(grade), @floatFromInt(i)));
        for (0..n) |j| a[i + j * n] = ops.mulReal(T, a[i + j * n], s);
    }
    const ipiv = try al.alloc(usize, n);
    defer al.free(ipiv);
    const inv = try al.alloc(T, n * n);
    defer al.free(inv);
    const f = try al.alloc(T, n * n);
    defer al.free(f);
    @memcpy(f, a);
    _ = try lu.getrf(al, check.io, T, n, n, f, n, ipiv);
    @memset(inv, ops.zero(T));
    for (0..n) |i| inv[i * (n + 1)] = ops.one(T);
    try lu.getrs(al, check.io, T, .n, n, n, f, n, ipiv, inv, n);
    for ([_]Norm{ .o, .i }) |nm| {
        const anorm = try lange(al, check.io, T, nm, n, n, a, n);
        const ainvnm = try lange(al, check.io, T, nm, n, n, inv, n);
        const rcondc = (1 / anorm) / ainvnm;
        const rcond = try gecon(al, check.io, T, nm, n, f, n, anorm);
        const r = lcheck.ratio(F, rcond, rcondc);
        if (!(r < lcheck.thresh)) {
            check.say("gecon {s} n {d} grade {d} norm {s}: rcond {e} against {e}, ratio {d}\n", .{ @typeName(T), n, @as(f64, @floatCast(grade)), @tagName(nm), @as(f64, @floatCast(rcond)), @as(f64, @floatCast(rcondc)), @as(f64, @floatCast(r)) });
            return error.TestUnexpectedResult;
        }
    }
}

test "lapack gecon: estimates within the reference's threshold of 1 / (|A| |inv(A)|), both norms" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        for ([_]usize{ 1, 2, 3, 5, 10, 50 }) |n| {
            // the last grade reaches a condition near 1 / sqrt(eps)
            const steep: F = @floatCast(std.math.pow(f64, @floatCast(@sqrt(std.math.floatEps(F))), 1 / @as(f64, @floatFromInt(@max(n, 2) - 1))));
            for ([_]F{ 1, 0.5, steep }, 0..) |g, k| try geconCase(T, al, n, g, n * 7 + k);
        }
    }
}

test "lapack gecon: n = 0, anorm = 0, 1 x 1, a zero pivot, NaN and infinite anorm, an overflowing estimate, bad arguments" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        const a = [_]T{ lcheck.cplx(T, 2, 0), lcheck.cplx(T, 1, 0), lcheck.cplx(T, 0, 0), lcheck.cplx(T, 3, 0) };
        try testing.expectEqual(@as(F, 1), try gecon(al, check.io, T, .o, 0, &a, 1, 1));
        try testing.expectEqual(@as(F, 0), try gecon(al, check.io, T, .o, 2, &a, 2, 0));
        // 1 x 1: exact
        try testing.expectEqual(@as(F, 1), try gecon(al, check.io, T, .i, 1, &a, 1, 2));
        // U singular: the overflowing solve stops at 0
        const s = [_]T{ lcheck.cplx(T, 2, 0), lcheck.cplx(T, 0.5, 0), lcheck.cplx(T, 1, 0), lcheck.cplx(T, 0, 0) };
        for ([_]Norm{ .o, .i }) |nm| try testing.expectEqual(@as(F, 0), try gecon(al, check.io, T, nm, 2, &s, 2, 3));
        try testing.expectError(error.NormOutOfRange, gecon(al, check.io, T, .o, 2, &a, 2, std.math.nan(F)));
        try testing.expectError(error.NormOutOfRange, gecon(al, check.io, T, .o, 2, &a, 2, std.math.inf(F)));
        // an anorm so small the estimate overflows: the reference's INFO = 1
        try testing.expectError(error.EstimateFailed, gecon(al, check.io, T, .o, 1, &a, 1, std.math.floatTrueMin(F)));
    }
    const x = [_]f64{ 1, 0, 0 };
    try testing.expectError(error.UnsupportedNorm, gecon(al, check.io, f64, .f, 1, &x, 1, 1));
    try testing.expectError(error.LeadingDimensionTooSmall, gecon(al, check.io, f64, .o, 2, &x, 1, 1));
    try testing.expectError(error.SliceTooShort, gecon(al, check.io, f64, .o, 2, &x, 2, 1));
    try testing.expectError(error.NormOutOfRange, gecon(al, check.io, f64, .o, 1, &x, 1, -1));
}

/// `latrs` on an `n x n` triangle with diagonal `d`: `scale` in [0, 1], `x`
/// finite, and `‖op(A) x - scale b‖ / (‖A‖ ‖x‖ n eps)` below the threshold.
fn latrsCase(comptime T: type, al: Allocator, upper: bool, op: Trans, unit: bool, n: usize, d: Float(T), seed: u64) !void {
    const F = Float(T);
    const a = try al.alloc(T, n * n);
    defer al.free(a);
    check.fill(T, a, seed);
    for (0..n) |j| for (0..n) |i| {
        if (if (upper) i > j else i < j) a[i + j * n] = lcheck.poisoned(T);
    };
    for (0..n) |i| a[i * (n + 1)] = if (unit) lcheck.poisoned(T) else ops.scalar(T, d);
    const b = try al.alloc(T, n);
    defer al.free(b);
    check.fill(T, b, seed + 1);
    const x = try al.alloc(T, n);
    defer al.free(x);
    @memcpy(x, b);
    const cnorm = try al.alloc(F, n);
    defer al.free(cnorm);
    const scale = try latrs(al, check.io, T, upper, op, unit, false, n, a, n, x, cnorm);
    try testing.expect(scale >= 0 and scale <= 1);
    var xmax: F = 0;
    for (x) |v| xmax = @max(xmax, ops.abs1(T, v));
    try testing.expect(std.math.isFinite(xmax));
    var amax: F = 0;
    var rmax: F = 0;
    for (0..n) |i| {
        var s = ops.mulReal(T, b[i], -scale);
        for (0..n) |j| {
            const r, const c = if (op == .n) .{ i, j } else .{ j, i };
            if (if (upper) r > c else r < c) continue;
            var e = if (r == c) (if (unit) ops.one(T) else ops.scalar(T, d)) else a[r + c * n];
            if (op == .c) e = ops.conj(T, e);
            amax = @max(amax, ops.abs1(T, e));
            s = ops.add(T, s, ops.mul(T, e, x[j]));
        }
        rmax = @max(rmax, ops.abs1(T, s));
    }
    const eps = std.math.floatEps(F);
    const nf: F = @floatFromInt(n);
    const denom = amax * xmax * nf * eps;
    if (!(rmax <= lcheck.thresh * denom)) {
        check.say("latrs {s} n {d} upper {} op {s} unit {} d {e}: residual {e}, bound {e}\n", .{ @typeName(T), n, upper, @tagName(op), unit, @as(f64, @floatCast(d)), @as(f64, @floatCast(rmax)), @as(f64, @floatCast(denom)) });
        return error.TestUnexpectedResult;
    }
}

test "lapack norm: latrs solves with scaling, ordinary and near-overflow triangles, every form" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        const op_list: []const Trans = if (comptime isComplex(T)) &.{ .n, .t, .c } else &.{ .n, .t };
        for ([_]bool{ true, false }) |upper| for (op_list) |op| for ([_]bool{ false, true }) |unit| {
            for ([_]usize{ 1, 2, 7, 40 }, 0..) |n, k| {
                // a diagonal of 3 (trsv), and one so small the solution would overflow
                for ([_]F{ 3, std.math.floatMin(F) * 64 }) |d| try latrsCase(T, al, upper, op, unit, n, d, n * 5 + k);
            }
        };
    }
}

test "lapack norm: laset fills the diagonal and the chosen strict triangle, wide and tall" {
    inline for (.{ f64, types.Complex64 }) |T| {
        const x = ops.scalar(T, -1);
        const al_ = ops.scalar(T, 2);
        const be = ops.scalar(T, 5);
        for ([_][2]usize{ .{ 3, 3 }, .{ 2, 4 }, .{ 4, 2 } }) |mn| for ([_]?types.Uplo{ null, .u, .l }) |uplo| {
            const m = mn[0];
            const n = mn[1];
            var a: [16]T = @splat(x);
            try auxiliary.laset(T, uplo, m, n, al_, be, &a, m);
            for (0..n) |j| for (0..m) |i| {
                const inside = uplo == null or (uplo == .u and i <= j) or (uplo == .l and i >= j);
                const want = if (!inside) x else if (i == j) be else al_;
                try testing.expectEqual(want, a[i + j * m]);
            };
        };
    }
    var a: [4]f64 = @splat(0);
    try testing.expectError(error.UnsupportedUplo, auxiliary.laset(f64, .v, 2, 2, 0, 0, &a, 2));
    try testing.expectError(error.LeadingDimensionTooSmall, auxiliary.laset(f64, null, 2, 2, 0, 0, &a, 1));
    try testing.expectError(error.SliceTooShort, auxiliary.laset(f64, null, 2, 3, 0, 0, &a, 2));
}

test "lapack norm: lassq, rscl and larfg on strided views; latrs and lacn2 checks" {
    // lassq: |(3, 4)| = 5 at stride 2, either sign
    const v = [_]f64{ 3, 99, 4 };
    for ([_]isize{ 2, -2 }) |st| {
        var scale: f64 = 0;
        var sumsq: f64 = 1;
        lassq(f64, try .init(&v, 2, st), &scale, &sumsq);
        try testing.expectApproxEqRel(@as(f64, 5), scale * @sqrt(sumsq), 4 * std.math.floatEps(f64));
    }
    // rscl: x / 4 at stride 2, the gap untouched
    var w = [_]f64{ 8, 7, 4 };
    try rscl(check.io, f64, 4, try .init(&w, 2, 2));
    try testing.expectEqualSlices(f64, &.{ 2, 7, 1 }, &w);
    // larfg: H [3; 4] = [-5; 0], v = [1; 0.5], tau = 1.6
    var alpha: f64 = 3;
    var x = [_]f64{4};
    const tau = try auxiliary.larfg(check.io, f64, &alpha, .of(&x));
    try testing.expectApproxEqAbs(@as(f64, -5), alpha, 1e-15);
    try testing.expectApproxEqAbs(@as(f64, 1.6), tau, 1e-15);
    try testing.expectApproxEqAbs(@as(f64, 0.5), x[0], 1e-15);
    const al = testing.allocator;
    var a = [_]f64{ 2, 0, 1, 2 };
    var b = [_]f64{ 1, 1 };
    var cn: [2]f64 = undefined;
    try testing.expectError(error.UnsupportedTranspose, latrs(al, check.io, f64, true, .o, false, false, 2, &a, 2, &b, &cn));
    try testing.expectError(error.LeadingDimensionTooSmall, latrs(al, check.io, f64, true, .n, false, false, 2, &a, 1, &b, &cn));
    try testing.expectError(error.SliceTooShort, latrs(al, check.io, f64, true, .n, false, false, 2, &a, 2, b[0..1], &cn));
    const Id = struct {
        fn apply(_: void, _: bool, _: [*]f64) error{}!void {}
    };
    try testing.expectError(error.SliceTooShort, lacn2(al, f64, error{}, 2, b[0..1], &cn, {}, Id.apply));
    // the identity's one norm
    try testing.expectEqual(@as(f64, 1), try lacn2(al, f64, error{}, 2, &b, &cn, {}, Id.apply));
}
