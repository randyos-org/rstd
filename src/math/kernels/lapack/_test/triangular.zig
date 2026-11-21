//! Tests of `../triangular.zig`.

const std = @import("std");
const types = @import("../../blas/types.zig");
const Vtable = @import("../../blas/Vtable.zig");
const ops = @import("../ops.zig");
const lapack_tri = @import("../triangular.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Uplo = types.Uplo;
const Diag = types.Diag;
const Trans = types.Trans;
const Float = types.Float;
const isComplex = types.isComplex;
const scalar = ops.scalar;
const firstZero = lapack_tri.firstZero;
const inverse = lapack_tri.inverse;
const lauu2 = lapack_tri.lauu2;
const lauum = lapack_tri.lauum;
const product = lapack_tri.product;
const trti2 = lapack_tri.trti2;
const trtri = lapack_tri.trtri;
const trtrs = lapack_tri.trtrs;

const testing = std.testing;
const blas = @import("../../blas/__root__.zig");
const check = @import("../../blas/_test/check.zig");
const guard = @import("../../blas/_test/guard.zig");
const reference = @import("../../blas/_test/reference.zig");
const lcheck = @import("check.zig");
const fixtures = @import("fixtures.zig");

const W = lcheck.W;
const test_types = lcheck.test_types;
const poisoned = lcheck.poisoned;
const isNanT = lcheck.isNanT;

const nan = std.math.nan(f64);

/// Sizes around every block edge: on `lcheck.shrink`'s blocking (`Q` 3,
/// `dtb_entries` 4) and on a core's own.
const shrunk_sizes = [_]usize{ 1, 2, 3, 4, 5, 8, 9, 12, 13, 14, 23 };
const own_sizes = [_]usize{ 1, 33, 70, 130 };

/// `f` over every test type, triangle, blocking and size, with the single
/// algorithm (`want` null) and the parallel one on one worker.
fn eachTri(vt: *const Vtable, comptime unit_too: bool, comptime f: anytype) anyerror!void {
    const al = testing.allocator;
    inline for (test_types) |T| for ([_]*const Vtable{ lcheck.shrink(vt), vt }, 0..) |v, vi| {
        guard.use(v);
        defer guard.restore();
        const sizes: []const usize = if (vi == 0) &shrunk_sizes else &own_sizes;
        for (sizes) |n| for ([_]bool{ true, false }) |upper| for ([_]bool{ false, true }) |unit| {
            if (unit and !unit_too) continue;
            if (check.rig() and n > 70) continue;
            for ([_]?u32{ null, 1 }) |want| {
                var t = try fixtures.Tri(T).init(al, n, upper, unit, n * 13 + @as(u64, @intFromBool(upper)) * 2 + @intFromBool(unit));
                defer t.deinit(al);
                f(T, &t, want) catch |e| {
                    check.say("{s}, blocking {s}, {s} algorithm\n", .{ @typeName(T), if (vi == 0) "shrunk" else "own", if (want == null) "single" else "parallel" });
                    return e;
                };
            }
        };
    };
}

fn inverseCase(comptime T: type, t: *fixtures.Tri(T), want: ?u32) !void {
    var poison: guard.Poison = .{ .child = testing.allocator };
    try inverse(poison.allocator(), check.io, T, t.upper, t.unit, t.n, t.a.ptr, t.lda, want);
    try t.expectInverse();
}

fn inverseCases(vt: *const Vtable) anyerror!void {
    return eachTri(vt, true, inverseCase);
}

test "trtri: A X = I, both triangles and diagonals, both algorithms, on every linked core" {
    try check.eachCore(inverseCases);
}

fn productCase(comptime T: type, t: *fixtures.Tri(T), want: ?u32) !void {
    var poison: guard.Poison = .{ .child = testing.allocator };
    try product(poison.allocator(), check.io, T, t.upper, t.n, t.a.ptr, t.lda, want);
    try t.expectProduct();
}

fn productCases(vt: *const Vtable) anyerror!void {
    return eachTri(vt, false, productCase);
}

test "lauum: U Uᴴ and Lᴴ L, both algorithms, on every linked core" {
    try check.eachCore(productCases);
}

test "trti2 and lauu2 are the unblocked forms, and every size picks an algorithm" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex32 }) |T| for ([_]bool{ true, false }) |upper| {
        var t = try fixtures.Tri(T).init(al, 9, upper, false, 3);
        defer t.deinit(al);
        try trti2(T, if (upper) .u else .l, .n, 9, t.a, t.lda);
        try t.expectInverse();
        t.store();
        try lauu2(T, if (upper) .u else .l, 9, t.a, t.lda);
        try t.expectProduct();
        // the public routines, past upstream's thresholds too
        for ([_]usize{ 40, 200 }) |n| {
            if (check.rig() and n > 70) continue;
            var p = try fixtures.Tri(T).init(al, n, upper, false, n);
            defer p.deinit(al);
            try trtri(al, check.io, T, if (upper) .u else .l, .n, n, p.a, p.lda);
            try p.expectInverse();
            p.store();
            try lauum(al, check.io, T, if (upper) .u else .l, n, p.a, p.lda);
            try p.expectProduct();
        }
    };
}

test "trtri: a zero on the diagonal is Singular, nothing written; a unit diagonal is not read" {
    const al = testing.failing_allocator;
    inline for (.{ f64, types.Complex64 }) |T| for ([_]Uplo{ .u, .l }) |uplo| {
        var a: [16]T = @splat(scalar(T, 0));
        for (0..4) |i| a[i * 5] = scalar(T, 2);
        a[5] = scalar(T, 0);
        a[15] = scalar(T, 0);
        // a part only zero is not zero
        if (comptime isComplex(T)) a[0] = .init(0, 1);
        const before = a;
        try testing.expectEqual(@as(?usize, 1), firstZero(T, 4, &a, 4));
        try testing.expectError(error.Singular, trtri(al, check.io, T, uplo, .n, 4, &a, 4));
        try testing.expectEqualSlices(T, &before, &a);
        // `u`: the zeros are ones
        try trtri(testing.allocator, check.io, T, uplo, .u, 4, &a, 4);
        for (0..4) |i| try testing.expect(lcheck.sameBits(T, before[i * 5], a[i * 5]));
    };
}

test "trtri, trti2, lauum, lauu2: a reversed uplo, a short slice, a small leading dimension; nothing read for n = 0" {
    const fa = testing.failing_allocator;
    var a: [4]f64 = @splat(nan);
    try testing.expectError(error.UnsupportedUplo, trtri(fa, check.io, f64, .v, .n, 2, &a, 2));
    try testing.expectError(error.LeadingDimensionTooSmall, trtri(fa, check.io, f64, .l, .u, 2, &a, 1));
    try testing.expectError(error.SliceTooShort, trtri(fa, check.io, f64, .l, .u, 2, &a, std.math.maxInt(usize) / 2));
    try testing.expectError(error.SliceTooShort, trtri(fa, check.io, f64, .l, .u, 3, &a, 3));
    try testing.expectError(error.UnsupportedUplo, trti2(f64, .m, .n, 2, &a, 2));
    try testing.expectError(error.LeadingDimensionTooSmall, trti2(f64, .u, .n, 2, &a, 1));
    try testing.expectError(error.UnsupportedUplo, lauum(fa, check.io, f64, .v, 2, &a, 2));
    try testing.expectError(error.LeadingDimensionTooSmall, lauum(fa, check.io, f64, .l, 2, &a, 1));
    try testing.expectError(error.UnsupportedUplo, lauu2(f64, .m, 2, &a, 2));
    try testing.expectError(error.LeadingDimensionTooSmall, lauu2(f64, .u, 0, &a, 0));
    try trtri(fa, check.io, f64, .u, .n, 0, &a, 1);
    try trti2(f64, .l, .n, 0, &a, 1);
    try lauum(fa, check.io, f64, .u, 0, &a, 1);
    try lauu2(f64, .l, 0, &a, 1);
    for (a) |v| try testing.expect(std.math.isNan(v));
}

/// `trtrs` on `t` for `nrhs` right-hand sides at `ldb = n + 1` (the padding
/// row NaN): `op(A) X - B` within `8 (n + 2)` ulps of `|op(A)| |X|`, the
/// padding and `A` untouched.
fn solveCase(comptime T: type, al: Allocator, t: *const fixtures.Tri(T), trans: Trans, nrhs: usize, seed: u64) !void {
    const n = t.n;
    const ldb = n + 1;
    const b = try al.alloc(T, ldb * nrhs);
    defer al.free(b);
    const b0 = try al.alloc(T, n * nrhs);
    defer al.free(b0);
    check.fill(T, b0, seed);
    for (0..nrhs) |j| {
        @memcpy(b[j * ldb ..][0..n], b0[j * n ..][0..n]);
        b[j * ldb + n] = poisoned(T);
    }
    var poison: guard.Poison = .{ .child = al };
    const uplo: Uplo = if (t.upper) .u else .l;
    const diag: Diag = if (t.unit) .u else .n;
    try trtrs(poison.allocator(), check.io, T, uplo, trans, diag, n, nrhs, t.a, t.lda, b, ldb);
    const tr = trans.transposed();
    const cj = isComplex(T) and trans.conjA();
    const eps: f128 = std.math.floatEps(Float(T));
    const ulps: f128 = @floatFromInt(8 * (n + 2));
    for (0..nrhs) |j| {
        try testing.expect(isNanT(T, b[j * ldb + n]));
        for (0..n) |i| {
            var s: W(T) = reference.wide(T, b0[i + j * n]);
            var mag: f128 = 0;
            for (0..n) |k| {
                const v = if (tr) t.w[k + i * n] else t.w[i + k * n];
                const x = reference.wide(T, b[k + j * ldb]);
                s = reference.subW(s, reference.mulW(if (cj) reference.conjW(v) else v, x));
                mag += reference.abs1(v) * reference.abs1(x);
            }
            const d = reference.abs1(s);
            if (!(d <= ulps * eps * mag)) {
                check.say("trtrs {s} {s}{s} {s} n={d} nrhs={d}: residual ({d},{d}) {e}, allowed {e}\n", .{ @typeName(T), @tagName(uplo), @tagName(diag), @tagName(trans), n, nrhs, i, j, @as(f64, @floatCast(d)), @as(f64, @floatCast(ulps * eps * mag)) });
                return error.TestUnexpectedResult;
            }
        }
    }
    try t.expectUntouched("trtrs");
}

fn solveCases(vt: *const Vtable) anyerror!void {
    const al = testing.allocator;
    inline for (test_types) |T| {
        guard.use(vt);
        defer guard.restore();
        for ([_]usize{ 1, 2, 7, 33, 70, 130 }) |n| for ([_]bool{ true, false }) |upper| for ([_]bool{ false, true }) |unit| {
            if (check.rig() and n > 70) continue;
            var t = try fixtures.Tri(T).init(al, n, upper, unit, n * 7 + @intFromBool(unit));
            defer t.deinit(al);
            for ([_]Trans{ .n, .t, .r, .c }) |trans| for ([_]usize{ 1, 3 }) |nrhs| {
                try solveCase(T, al, &t, trans, nrhs, n + nrhs);
            };
        };
    }
}

test "trtrs: op(A) X = B, every op, both triangles and diagonals, one and several right-hand sides, on every linked core" {
    try check.eachCore(solveCases);
}

test "trtrs: a zero on the diagonal before B is touched, even with no right-hand side; bad arguments" {
    const fa = testing.failing_allocator;
    var a: [4]f64 = .{ 1, 2, 3, 0 };
    var b: [4]f64 = @splat(nan);
    try testing.expectError(error.Singular, trtrs(fa, check.io, f64, .u, .n, .n, 2, 2, &a, 2, &b, 2));
    try testing.expectError(error.Singular, trtrs(fa, check.io, f64, .l, .t, .n, 2, 0, &a, 2, &b, 2));
    // a unit diagonal is not read
    try trtrs(fa, check.io, f64, .u, .n, .u, 2, 0, &a, 2, &b, 2);
    for (b) |v| try testing.expect(std.math.isNan(v));
    try testing.expectError(error.UnsupportedUplo, trtrs(fa, check.io, f64, .v, .n, .n, 2, 2, &a, 2, &b, 2));
    try testing.expectError(error.UnsupportedTranspose, trtrs(fa, check.io, f64, .u, .o, .n, 2, 2, &a, 2, &b, 2));
    try testing.expectError(error.LeadingDimensionTooSmall, trtrs(fa, check.io, f64, .u, .n, .n, 2, 1, &a, 1, &b, 2));
    try testing.expectError(error.LeadingDimensionTooSmall, trtrs(fa, check.io, f64, .u, .n, .n, 2, 1, &a, 2, &b, 1));
    try testing.expectError(error.SliceTooShort, trtrs(fa, check.io, f64, .u, .n, .n, 2, 2, &a, 2, &b, std.math.maxInt(usize) / 2));
    try testing.expectError(error.SliceTooShort, trtrs(fa, check.io, f64, .u, .n, .n, 2, 3, &a, 2, &b, 2));
    try trtrs(fa, check.io, f64, .u, .n, .n, 0, 2, &a, 1, &b, 1);
}

/// Inputs of a threaded-bits case: a triangle, copied into each run's `out`.
fn Threaded(comptime T: type) type {
    return struct {
        t: *const fixtures.Tri(T),

        pub fn inv(c: @This(), io: Io, want: u32, out: []T) anyerror!void {
            @memcpy(out, c.t.a);
            try inverse(testing.allocator, io, T, c.t.upper, c.t.unit, c.t.n, out.ptr, c.t.lda, want);
        }

        pub fn prod(c: @This(), io: Io, want: u32, out: []T) anyerror!void {
            @memcpy(out, c.t.a);
            try product(testing.allocator, io, T, c.t.upper, c.t.n, out.ptr, c.t.lda, want);
        }

        pub fn solve(c: @This(), io: Io, want: u32, out: []T) anyerror!void {
            // `trtrs` threads by `blas.config`'s workers
            const saved = blas.config.cpu_count;
            defer blas.config.cpu_count = saved;
            blas.config.cpu_count = want;
            const n = c.t.n;
            check.fill(T, out, 9);
            try trtrs(testing.allocator, io, T, if (c.t.upper) .u else .l, .c, .n, n, out.len / n, c.t.a, c.t.lda, out, n);
        }
    };
}

test "trtri's and lauum's parallel algorithms and trtrs give one worker's bits at any worker count" {
    _ = blas.coreName();
    const al = testing.allocator;
    inline for (.{ f32, f64, types.Complex32, types.Complex64 }) |T| for ([_]bool{ true, false }) |upper| {
        const n: usize = if (check.rig()) 70 else 130;
        var t = try fixtures.Tri(T).init(al, n, upper, false, 21);
        defer t.deinit(al);
        const c: Threaded(T) = .{ .t = &t };
        try lcheck.expectSameBits(T, c, Threaded(T).inv, t.a.len, "trtri");
        try lcheck.expectSameBits(T, c, Threaded(T).prod, t.a.len, "lauum");
        try lcheck.expectSameBits(T, c, Threaded(T).solve, n * 37, "trtrs");
    };
}
