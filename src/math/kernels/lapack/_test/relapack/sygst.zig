//! Tests of ReLAPACK's `sygst` and `hegst`.

const std = @import("std");
const types = @import("../../../blas/types.zig");
const check = @import("../../../blas/_test/check.zig");
const reference = @import("../../../blas/_test/reference.zig");
const lcheck = @import("../check.zig");
const fixtures = @import("../fixtures.zig");
const symmetric = @import("../symmetric.zig");
const relapack = @import("../../relapack.zig");

const testing = std.testing;
const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const isComplex = types.isComplex;
const test_types = lcheck.test_types;
const isNanT = lcheck.isNanT;
const W = lcheck.W;

// ---- sygst tests ------------------------------------------------------------------

pub fn gstCase(comptime T: type, al: Allocator, problem: relapack.GeneralizedProblem, upper: bool, n: usize, seed: u64) !void {
    // `B` from a known Cholesky factor, `A` Hermitian, both in `uplo`'s triangle
    var p = try fixtures.CholProblem(T).init(al, n, upper, seed);
    defer p.deinit(al);
    var h = try symmetric.Indefinite(T, isComplex(T)).init(al, n, upper, false, seed + 1);
    defer h.deinit(al);
    const ldb = n + 1;
    const b = try al.alloc(T, @max(1, ldb * n));
    defer al.free(b);
    @memset(b, lcheck.poisoned(T));
    // `L`, or `U = Lᴴ`, whole in `W(T)`
    const f = try al.alloc(W(T), n * n);
    defer al.free(f);
    for (0..n) |j| for (0..n) |i| {
        const l = reference.wide(T, p.l[i + j * n]);
        f[i + j * n] = if (upper) reference.conjW(reference.wide(T, p.l[j + i * n])) else l;
        if (if (upper) i <= j else i >= j) b[i + j * ldb] = lcheck.toT(T, f[i + j * n]);
    };
    const routine = if (comptime isComplex(T)) relapack.hegst else relapack.sygst;
    try routine(al, check.io, T, problem, lcheck.uploOf(upper), n, h.a, h.lda, b, ldb);
    // the result, whole
    const c = try al.alloc(W(T), n * n);
    defer al.free(c);
    for (0..n) |j| for (0..n) |i| {
        const s = if (upper) i <= j else i >= j;
        const v = reference.wide(T, if (s) h.a[i + j * h.lda] else h.a[j + i * h.lda]);
        c[i + j * n] = if (s) v else reference.conjW(v);
    };
    // `M X Mᴴ` against the other side: `L C Lᴴ` or `Uᴴ C U` against `A`
    // (`ax_lbx`); `C` against `Lᴴ A L` or `U A Uᴴ`
    const m = try al.alloc(W(T), n * n);
    defer al.free(m);
    for (0..n) |j| for (0..n) |i| {
        const flip = (problem == .ax_equals_lambda_bx) == upper;
        m[i + j * n] = if (flip) reference.conjW(f[j + i * n]) else f[i + j * n];
    };
    const r = try al.alloc(W(T), n * n);
    defer al.free(r);
    const mag = try al.alloc(f128, n * n);
    defer al.free(mag);
    try symmetric.sandwich(T, al, n, m, if (problem == .ax_equals_lambda_bx) c else h.w, r, mag);
    const other = if (problem == .ax_equals_lambda_bx) h.w else c;
    const eps: f128 = std.math.floatEps(Float(T));
    const ulps: f128 = @floatFromInt(32 * (n + 2));
    for (0..n) |j| for (0..n) |i| {
        const d = reference.abs1(reference.subW(r[i + j * n], other[i + j * n]));
        const allowed = ulps * eps * (mag[i + j * n] + reference.abs1(other[i + j * n])) + std.math.floatMin(Float(T));
        if (!(d <= allowed)) {
            check.say("sygst {s} {t} {s} n={d}: ({d},{d}) off by {e}, allowed {e}\n", .{ @typeName(T), problem, if (upper) "u" else "l", n, i, j, @as(f64, @floatCast(d)), @as(f64, @floatCast(allowed)) });
            return error.TestUnexpectedResult;
        }
    };
    for (0..n) |j| for (0..h.lda) |i| {
        const x = h.a[i + j * h.lda];
        const inside = i < n and (if (upper) i <= j else i >= j);
        if (!inside and !isNanT(T, x)) {
            check.say("sygst {s} n={d}: wrote ({d},{d}), outside its triangle\n", .{ @typeName(T), n, i, j });
            return error.TestUnexpectedResult;
        }
        if (comptime isComplex(T)) if (i == j and i < n and x.im != 0) {
            check.say("hegst {s} n={d}: ({d},{d}) has imaginary part {e}\n", .{ @typeName(T), n, i, i, x.im });
            return error.TestUnexpectedResult;
        };
    };
}

test "relapack sygst, hegst: every problem, both triangles, against the products, every type" {
    const al = testing.allocator;
    inline for (test_types) |T| for (std.enums.values(relapack.GeneralizedProblem), 1..) |problem, pi| for ([_]bool{ true, false }) |upper| {
        for ([_]usize{ 1, 2, 24, 25, 40, 57 }) |n| {
            if (check.rig() and n > 40) continue;
            gstCase(T, al, problem, upper, n, n * 5 + pi) catch |e| {
                check.say("sygst {s} {t} n={d}\n", .{ @typeName(T), problem, n });
                return e;
            };
        }
    };
}

test "relapack sygst, hegst: past the crossover twice over, f64 and c64" {
    if (check.rig()) return error.SkipZigTest;
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64 }) |T| for ([_]relapack.GeneralizedProblem{ .ax_equals_lambda_bx, .abx_equals_lambda_x }) |problem| for ([_]bool{ true, false }) |upper| {
        try gstCase(T, al, problem, upper, 100, 7);
    };
}

test "relapack sygst, hegst: a reversed uplo, short slices, small leading dimensions; nothing read for n = 0" {
    const al = testing.allocator;
    var a: [4]f64 = @splat(lcheck.poisoned(f64));
    const b: [4]f64 = @splat(lcheck.poisoned(f64));
    try testing.expectError(error.UnsupportedUplo, relapack.sygst(al, check.io, f64, .ax_equals_lambda_bx, .m, 2, &a, 2, &b, 2));
    try testing.expectError(error.LeadingDimensionTooSmall, relapack.sygst(al, check.io, f64, .abx_equals_lambda_x, .u, 2, &a, 1, &b, 2));
    try testing.expectError(error.LeadingDimensionTooSmall, relapack.sygst(al, check.io, f64, .bax_equals_lambda_x, .l, 2, &a, 2, &b, 1));
    try testing.expectError(error.SliceTooShort, relapack.sygst(al, check.io, f64, .ax_equals_lambda_bx, .l, 3, &a, 3, &b, 3));
    try relapack.sygst(al, check.io, f64, .ax_equals_lambda_bx, .l, 0, &a, 1, &b, 1);
    for (a) |x| try testing.expect(std.math.isNan(x));
    var z: [1]types.Complex32 = @splat(lcheck.poisoned(types.Complex32));
    try testing.expectError(error.LeadingDimensionTooSmall, relapack.hegst(al, check.io, types.Complex32, .ax_equals_lambda_bx, .u, 1, &z, 1, &z, 0));
}

pub const BigGst = struct {
    a: []const f64,
    b: []const f64,

    pub fn run(self: BigGst, io: Io, out: []f64) anyerror!void {
        @memcpy(out, self.a);
        try relapack.sygst(testing.allocator, io, f64, .ax_equals_lambda_bx, .l, lcheck.big_n, out, lcheck.big_n, self.b, lcheck.big_n);
    }
};

test "relapack sygst: the same bits inline and threaded" {
    const al = testing.allocator;
    const a = try fixtures.dominant(al, lcheck.big_n, 5);
    defer al.free(a);
    const b = try fixtures.dominant(al, lcheck.big_n, 6);
    defer al.free(b);
    try relapack.potrf(al, check.io, f64, .l, lcheck.big_n, b, lcheck.big_n);
    try lcheck.expectIoBits(f64, BigGst{ .a = a, .b = b }, BigGst.run, lcheck.big_n * lcheck.big_n, "sygst");
}
