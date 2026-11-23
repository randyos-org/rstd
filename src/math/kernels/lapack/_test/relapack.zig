//! Tests of `relapack.zig`'s routines: `potrf`, `getrf`, `trtri` and
//! `lauum` (`-Drelapack`), `sygst`, `hegst` and `pbtrf`.

const std = @import("std");
const types = @import("../../blas/types.zig");
const check = @import("../../blas/_test/check.zig");
const guard = @import("../../blas/_test/guard.zig");
const reference = @import("../../blas/_test/reference.zig");
const Vtable = @import("../../blas/Vtable.zig");
const options = @import("lapack_options");
const lcheck = @import("check.zig");
const fixtures = @import("fixtures.zig");
const ops = @import("../ops.zig");
const chol = @import("../cholesky.zig");
const lu = @import("../lu.zig");
const tri = @import("../triangular.zig");
const relapack = @import("../relapack.zig");

const testing = std.testing;
const Io = std.Io;
const scalar = ops.scalar;
const test_types = lcheck.test_types;

/// Around the crossover and `split`'s rounding for every element size.
pub const test_sizes = [_]usize{ 1, 2, 24, 25, 33, 49, 70, 100 };

pub fn potrfCases(vt: *const Vtable) anyerror!void {
    const al = testing.allocator;
    guard.use(vt);
    defer guard.restore();
    inline for (test_types) |T| for ([_]bool{ true, false }) |upper| for (test_sizes) |n| {
        if (check.rig() and n > 70) continue;
        var p = try fixtures.CholProblem(T).init(al, n, upper, n * 37 + @intFromBool(upper));
        defer p.deinit(al);
        relapack.potrf(al, check.io, T, lcheck.uploOf(upper), n, p.a, p.lda) catch |e| {
            check.say("relapack potrf {s} n={d} upper={}\n", .{ @typeName(T), n, upper });
            return e;
        };
        p.expectFactor() catch |e| {
            check.say("relapack potrf {s} n={d}\n", .{ @typeName(T), n });
            return e;
        };
    };
}

test "relapack potrf: A = L Lᴴ from a known L, both triangles, on every linked core" {
    try check.eachCore(potrfCases);
}

test "relapack potrf: the first pivot that is not positive, either half" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex32 }) |T| for ([_]bool{ true, false }) |upper| {
        const n = 49;
        var p = try fixtures.CholProblem(T).init(al, n, upper, 5);
        defer p.deinit(al);
        for ([_]usize{ 0, 23, 24, 25, 40, 48 }) |k| {
            p.store();
            p.a[k * (p.lda + 1)] = scalar(T, -100);
            try testing.expectError(error.NotPositiveDefinite, relapack.potrf(al, check.io, T, lcheck.uploOf(upper), n, p.a, p.lda));
            p.store();
            p.a[k * (p.lda + 1)] = scalar(T, -100);
            try testing.expectEqual(@as(?usize, k), try relapack.potrfRec(al, check.io, T, lcheck.uploOf(upper), n, p.a.ptr, p.lda));
        }
    };
}

pub const getrf_shapes = [_][2]usize{ .{ 1, 1 }, .{ 24, 24 }, .{ 25, 25 }, .{ 50, 50 }, .{ 70, 33 }, .{ 33, 70 }, .{ 100, 60 }, .{ 26, 100 } };

pub fn getrfCases(vt: *const Vtable) anyerror!void {
    const al = testing.allocator;
    guard.use(vt);
    defer guard.restore();
    inline for (test_types) |T| for (getrf_shapes) |s| {
        if (check.rig() and s[0] * s[1] > 3000) continue;
        var p = try fixtures.LuProblem(T).init(al, s[0], s[1], s[0] * 101 + s[1]);
        defer p.deinit(al);
        try testing.expectEqual(@as(?usize, null), try relapack.getrf(al, check.io, T, s[0], s[1], p.a, p.lda, p.ipiv));
        p.expectFactor(al) catch |e| {
            check.say("relapack getrf {s}\n", .{@typeName(T)});
            return e;
        };
    };
}

test "relapack getrf: P L U = A, square, tall and wide, on every linked core" {
    try check.eachCore(getrfCases);
}

test "relapack getrf: a zero column gives that pivot, the factorization completed, either half" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64 }) |T| for ([_]usize{ 0, 10, 30, 45 }) |k| {
        var p = try fixtures.LuProblem(T).init(al, 50, 50, 9);
        defer p.deinit(al);
        for (0..50) |i| p.aw[i + k * 50] = reference.wide(T, scalar(T, 0));
        p.store();
        try testing.expectEqual(@as(?usize, k), try relapack.getrf(al, check.io, T, 50, 50, p.a, p.lda, p.ipiv));
        try p.expectFactor(al);
    };
}

pub fn triCases(vt: *const Vtable) anyerror!void {
    const al = testing.allocator;
    guard.use(vt);
    defer guard.restore();
    inline for (test_types) |T| for ([_]bool{ true, false }) |upper| for ([_]bool{ false, true }) |unit| for (test_sizes) |n| {
        if (check.rig() and n > 70) continue;
        var t = try fixtures.Tri(T).init(al, n, upper, unit, n * 13 + @intFromBool(upper));
        defer t.deinit(al);
        try relapack.trtri(al, check.io, T, lcheck.uploOf(upper), if (unit) .u else .n, n, t.a, t.lda);
        try t.expectInverse();
        if (unit) continue;
        t.store();
        try relapack.lauum(al, check.io, T, lcheck.uploOf(upper), n, t.a, t.lda);
        try t.expectProduct();
    };
}

test "relapack trtri and lauum: A X = I and U Uᴴ, Lᴴ L, both triangles, on every linked core" {
    try check.eachCore(triCases);
}

test "relapack trtri: a zero first diagonal element is singular, nothing written" {
    const al = testing.allocator;
    var t = try fixtures.Tri(f64).init(al, 30, true, false, 3);
    defer t.deinit(al);
    t.a[0] = 0;
    const before = try al.dupe(f64, t.a);
    defer al.free(before);
    try testing.expectError(error.Singular, relapack.trtri(al, check.io, f64, .u, .n, 30, t.a, t.lda));
    for (before, t.a) |x, y| try testing.expect(lcheck.sameBits(f64, x, y));
}

test "relapack potrf, getrf, trtri, lauum: the first bad argument, nothing read when empty" {
    const al = testing.allocator;
    var a: [4]f64 = @splat(lcheck.poisoned(f64));
    var ipiv: [2]usize = @splat(std.math.maxInt(usize));
    try testing.expectError(error.UnsupportedUplo, relapack.potrf(al, check.io, f64, .v, 2, &a, 2));
    try testing.expectError(error.SliceTooShort, relapack.potrf(al, check.io, f64, .u, 3, &a, 3));
    try testing.expectError(error.LeadingDimensionTooSmall, relapack.potrf(al, check.io, f64, .u, 2, &a, 1));
    try testing.expectError(error.LeadingDimensionTooSmall, relapack.getrf(al, check.io, f64, 2, 2, &a, 1, &ipiv));
    try testing.expectError(error.SliceTooShort, relapack.getrf(al, check.io, f64, 2, 3, &a, 2, &ipiv));
    try testing.expectError(error.SliceTooShort, relapack.getrf(al, check.io, f64, 2, 2, &a, 2, ipiv[0..1]));
    try testing.expectError(error.UnsupportedUplo, relapack.trtri(al, check.io, f64, .m, .n, 2, &a, 2));
    try testing.expectError(error.SliceTooShort, relapack.trtri(al, check.io, f64, .u, .n, 3, &a, 3));
    try testing.expectError(error.LeadingDimensionTooSmall, relapack.lauum(al, check.io, f64, .l, 2, &a, 1));
    try relapack.potrf(al, check.io, f64, .u, 0, &a, 1);
    try testing.expectEqual(@as(?usize, null), try relapack.getrf(al, check.io, f64, 0, 2, &a, 1, &ipiv));
    try relapack.trtri(al, check.io, f64, .l, .n, 0, &a, 1);
    try relapack.lauum(al, check.io, f64, .l, 0, &a, 1);
    for (a) |x| try testing.expect(std.math.isNan(x));
    try testing.expectEqual(std.math.maxInt(usize), ipiv[0]);
}

test "with -Drelapack, ReLAPACK's potrf, getrf, trtri and lauum agree with OpenBLAS's" {
    if (!options.relapack) return error.SkipZigTest;
    const root = @import("../__root__.zig");
    const al = testing.allocator;
    const io = check.io;
    inline for (.{ f32, f64, types.Complex32, types.Complex64 }) |T| for ([_]bool{ true, false }) |upper| {
        const n = 100;
        const uplo = lcheck.uploOf(upper);
        var p = try fixtures.CholProblem(T).init(al, n, upper, 77);
        defer p.deinit(al);
        const x = try al.dupe(T, p.a);
        defer al.free(x);
        try root.potrf(al, io, T, uplo, n, x, p.lda);
        try chol.potrf(al, io, T, uplo, n, p.a, p.lda);
        try fixtures.expectAgree(T, "potrf", n, x, p.a);

        var t = try fixtures.Tri(T).init(al, n, upper, false, 78);
        defer t.deinit(al);
        const z = try al.dupe(T, t.a);
        defer al.free(z);
        try root.trtri(al, io, T, uplo, .n, n, z, t.lda);
        try tri.trtri(al, io, T, uplo, .n, n, t.a, t.lda);
        try fixtures.expectAgree(T, "trtri", n, z, t.a);
        t.store();
        @memcpy(z, t.a);
        try root.lauum(al, io, T, uplo, n, z, t.lda);
        try tri.lauum(al, io, T, uplo, n, t.a, t.lda);
        try fixtures.expectAgree(T, "lauum", n, z, t.a);

        var g = try fixtures.LuProblem(T).init(al, n, n, 79 + @as(u64, @intFromBool(upper)));
        defer g.deinit(al);
        const y = try al.dupe(T, g.a);
        defer al.free(y);
        const ip = try al.alloc(usize, n);
        defer al.free(ip);
        _ = try root.getrf(al, io, T, n, n, y, g.lda, ip);
        _ = try lu.getrf(al, io, T, n, n, g.a, g.lda, g.ipiv);
        try testing.expectEqualSlices(usize, g.ipiv, ip);
        try fixtures.expectAgree(T, "getrf", n, y, g.a);
    };
}

pub const Big = struct {
    a: []const f64,

    pub fn potrfRun(self: Big, io: Io, out: []f64) anyerror!void {
        @memcpy(out, self.a);
        try relapack.potrf(testing.allocator, io, f64, .l, lcheck.big_n, out, lcheck.big_n);
    }

    pub fn getrfRun(self: Big, io: Io, out: []f64) anyerror!void {
        @memcpy(out, self.a);
        var ipiv: [lcheck.big_n]usize = undefined;
        _ = try relapack.getrf(testing.allocator, io, f64, lcheck.big_n, lcheck.big_n, out, lcheck.big_n, &ipiv);
    }
};

test "relapack potrf and getrf: the same bits inline and threaded" {
    const al = testing.allocator;
    const a = try fixtures.dominant(al, lcheck.big_n, 1);
    defer al.free(a);
    try lcheck.expectIoBits(f64, Big{ .a = a }, Big.potrfRun, lcheck.big_n * lcheck.big_n, "potrf");
    try lcheck.expectIoBits(f64, Big{ .a = a }, Big.getrfRun, lcheck.big_n * lcheck.big_n, "getrf");
}
