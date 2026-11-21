//! Tests of `../lu.zig`'s `getrs`, `gesv` and `getri`, and of `getrf`'s lookahead.

const std = @import("std");
const mem = @import("../../../../mem.zig");
const types = @import("../../blas/types.zig");
const driver = @import("../../blas/blas3/driver.zig");
const Vtable = @import("../../blas/Vtable.zig");
const ops = @import("../ops.zig");
const interface = @import("../interface.zig");
const blas = @import("../../blas/__root__.zig");
const check = @import("../../blas/_test/check.zig");
const guard = @import("../../blas/_test/guard.zig");
const reference = @import("../../blas/_test/reference.zig");
const lcheck = @import("check.zig");
const fixtures = @import("fixtures.zig");
const lu = @import("../lu.zig");
const lu_tests = @import("lu.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Trans = types.Trans;
const Float = types.Float;
const isComplex = types.isComplex;
const scalar = ops.scalar;
const add = ops.add;
const mul = ops.mul;
const testing = std.testing;
const W = lcheck.W;
const test_types = lcheck.test_types;
const poisoned = lcheck.poisoned;
const isNanT = lcheck.isNanT;
const fromF64 = lcheck.fromF64;
const factor = lu.factor;
const factorAhead = lu.factorAhead;
const gesv = lu.gesv;
const getrf = lu.getrf;
const getri = lu.getri;
const getrs = lu.getrs;
const nan = lu_tests.nan;
const shrink = lu_tests.shrink;
const solve = lu.solve;
const unset = lu_tests.unset;

/// `getrs` on `p`'s factor, `nrhs` right-hand sides with `ldb = n + 1` (the
/// padding row NaN), then the residual in `W(T)` within `8 (n + 2)` ulps of
/// `L` and `U` times `X` in magnitude: `P (A X - B)` against
/// `|L| |U| |X|`, or for `op` transposed `op(A) X - B` against
/// `(|L| |U|)ᵀ P |X|`. With `nan_at`, `B(i, j)` is NaN: column `j` must come
/// out all NaN, the others as without it, bit for bit.
fn solveCase(comptime T: type, al: Allocator, p: *const fixtures.LuProblem(T), trans: Trans, nrhs: usize, seed: u64, nan_at: ?[2]usize) !void {
    const n = p.n;
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
    try getrs(poison.allocator(), check.io, T, trans, n, nrhs, p.a, p.lda, p.ipiv, b, ldb);
    // `|L| |U|`
    const mag = try al.alloc(f128, n * n);
    defer al.free(mag);
    for (0..n) |j| for (0..n) |i| {
        var s: f128 = 0;
        for (0..@min(i, j) + 1) |k| s += reference.abs1(p.l(i, k)) * reference.abs1(p.u(k, j));
        mag[i + j * n] = s;
    };
    const t = trans.transposed();
    const cj = isComplex(T) and trans.conjA();
    const F = Float(T);
    const eps: f128 = std.math.floatEps(F);
    const ulps: f128 = @floatFromInt(8 * (n + 2));
    const res = try al.alloc(W(T), n);
    defer al.free(res);
    const xm = try al.alloc(f128, n);
    defer al.free(xm);
    for (0..nrhs) |j| {
        const x = b[j * ldb ..][0..n];
        try testing.expect(isNanT(T, b[j * ldb + n]));
        for (0..n) |i| {
            var s: W(T) = reference.wide(T, b0[i + j * n]);
            for (0..n) |k| {
                const v = if (t) p.aw[k + i * n] else p.aw[i + k * n];
                s = reference.subW(s, reference.mulW(if (cj) reference.conjW(v) else v, reference.wide(T, x[k])));
            }
            res[i] = s;
            xm[i] = reference.abs1(reference.wide(T, x[i]));
        }
        // `P` on the residual, or on `|X|`
        for (0..n) |k| {
            const q = p.ipiv[k];
            if (t) mem.swap(f128, &xm[k], &xm[q]) else mem.swap(W(T), &res[k], &res[q]);
        }
        for (0..n) |i| {
            var bound: f128 = 0;
            for (0..n) |k| bound += (if (t) mag[k + i * n] else mag[i + k * n]) * xm[k];
            const d = reference.abs1(res[i]);
            if (!(d <= ulps * eps * bound + std.math.floatMin(F))) {
                check.say("getrs {s} {s} n={d} nrhs={d}: residual ({d},{d}) {e}, allowed {e}\n", .{ @typeName(T), @tagName(trans), n, nrhs, i, j, @as(f64, @floatCast(d)), @as(f64, @floatCast(ulps * eps * bound)) });
                return error.TestUnexpectedResult;
            }
        }
    }
    const at = nan_at orelse return;
    const bn = try al.alloc(T, ldb * nrhs);
    defer al.free(bn);
    for (0..nrhs) |j| @memcpy(bn[j * ldb ..][0..n], b0[j * n ..][0..n]);
    bn[at[0] + at[1] * ldb] = poisoned(T);
    try getrs(poison.allocator(), check.io, T, trans, n, nrhs, p.a, p.lda, p.ipiv, bn, ldb);
    try lcheck.expectNanColumn(T, "getrs", @tagName(trans), n, nrhs, ldb, at, bn, b);
}

const solve_sizes = [_]usize{ 1, 2, 3, 7, 17, 33, 50, 70, 130 };

/// The real or the complex test types: a test each, so sparc64-debug's images
/// fit the firmware's window.
fn solveCases(comptime complex: bool, vt: *const Vtable) anyerror!void {
    return lcheck.eachBlocking(vt, shrink, test_types, SolveOn(complex).f);
}

fn SolveOn(comptime complex: bool) type {
    return struct {
        fn f(comptime T: type, v: *const Vtable, vi: usize) anyerror!void {
            if (comptime isComplex(T) != complex) return;
            _ = v;
            const al = testing.allocator;
            for (solve_sizes) |n| {
                if (check.rig() and n > 70) continue;
                var p = try fixtures.LuProblem(T).init(al, n, n, n * 17 + vi);
                defer p.deinit(al);
                _ = try p.factorIt(al);
                for ([_]Trans{ .n, .t, .r, .c }) |trans| for ([_]usize{ 1, 3 }) |nrhs| {
                    const nan_at: [2]usize = .{ (n * 7 + nrhs) % n, nrhs / 2 };
                    solveCase(T, al, &p, trans, nrhs, n + nrhs, nan_at) catch |e| {
                        check.say("{s}, blocking {s}\n", .{ @typeName(T), if (vi == 0) "shrunk" else "own" });
                        return e;
                    };
                };
            }
        }
    };
}

fn solveReal(vt: *const Vtable) anyerror!void {
    return solveCases(false, vt);
}

fn solveComplex(vt: *const Vtable) anyerror!void {
    return solveCases(true, vt);
}

test "getrs: op(A) X = B on getrf's factor, every op, NaN in one column of B, real, on every linked core" {
    try check.eachCore(solveReal);
}

test "getrs: op(A) X = B on getrf's factor, every op, NaN in one column of B, complex, on every linked core" {
    try check.eachCore(solveComplex);
}

test "getrs: a known 2 x 2 solve, A and its transpose" {
    inline for (.{ f64, types.Complex64 }) |T| {
        // [[1, 2], [3, 4]] x = [5, 6]: x = [-4, 9/2]; transposed, [-1, 2]
        for ([_]struct { Trans, [2]f64 }{ .{ .n, .{ -4, 4.5 } }, .{ .t, .{ -1, 2 } } }) |c| {
            var a: [4]T = .{ scalar(T, 1), scalar(T, 3), scalar(T, 2), scalar(T, 4) };
            var ipiv: [2]usize = undefined;
            _ = try getrf(testing.failing_allocator, check.io, T, 2, 2, &a, 2, &ipiv);
            var b: [2]T = .{ scalar(T, 5), scalar(T, 6) };
            try getrs(testing.allocator, check.io, T, c[0], 2, 1, &a, 2, &ipiv, &b, 2);
            for (b, c[1]) |x, w| try check.expectClose(T, scalar(T, @floatCast(w)), x, 1e-14);
        }
    }
}

test "getrs: a trans it does not take, short slices, small leading dimensions, a pivot outside A; nothing touched when empty" {
    const fa = testing.failing_allocator;
    var a: [4]f64 = @splat(nan);
    var b: [4]f64 = @splat(nan);
    const ipiv: [2]usize = .{ 0, 1 };
    try testing.expectError(error.UnsupportedTranspose, getrs(fa, check.io, f64, .o, 2, 2, &a, 2, &ipiv, &b, 2));
    try testing.expectError(error.LeadingDimensionTooSmall, getrs(fa, check.io, f64, .n, 2, 1, &a, 1, &ipiv, &b, 2));
    try testing.expectError(error.LeadingDimensionTooSmall, getrs(fa, check.io, f64, .n, 2, 1, &a, 2, &ipiv, &b, 1));
    try testing.expectError(error.LeadingDimensionTooSmall, getrs(fa, check.io, f64, .n, 0, 1, &a, 1, &ipiv, &b, 0));
    try testing.expectError(error.SliceTooShort, getrs(fa, check.io, f64, .n, 2, 2, &a, std.math.maxInt(usize) / 2, &ipiv, &b, 2));
    try testing.expectError(error.SliceTooShort, getrs(fa, check.io, f64, .t, 2, 3, &a, 2, &ipiv, &b, 2));
    try testing.expectError(error.SliceTooShort, getrs(fa, check.io, f64, .t, 2, 2, &a, 2, ipiv[0..1], &b, 2));
    for ([_][2]usize{ .{ 2, 1 }, .{ 0, 2 }, .{ unset, 0 } }) |bad| {
        try testing.expectError(error.PivotOutOfRange, getrs(fa, check.io, f64, .t, 2, 2, &a, 2, &bad, &b, 2));
    }
    try getrs(fa, check.io, f64, .n, 0, 2, &a, 1, &ipiv, &b, 1);
    try getrs(fa, check.io, f64, .n, 2, 0, &a, 2, &ipiv, &b, 2);
    for (b) |x| try testing.expect(std.math.isNan(x));
}

/// Inputs of a threaded-bits case: `A`, copied into each run's `out`, then
/// `ipiv` and the first zero pivot after it.
fn Threaded(comptime T: type) type {
    return struct {
        p: *const fixtures.LuProblem(T),

        pub fn len(c: @This()) usize {
            return c.p.a.len + @min(c.p.m, c.p.n) + 1;
        }

        /// `factor` at one worker, else `factorAhead`.
        pub fn run(c: @This(), io: Io, want: u32, out: []T) anyerror!void {
            const p = c.p;
            const al = testing.allocator;
            const ipiv = try al.alloc(usize, @min(p.m, p.n));
            defer al.free(ipiv);
            @memcpy(out[0..p.a.len], p.a);
            const vt = interface.table(T);
            const blk = driver.Blocking.of(T, vt);
            const k = if (want == 1)
                try factor(al, io, T, vt, blk.un, blk.q, p.m, p.n, out.ptr, p.lda, ipiv.ptr)
            else
                try factorAhead(al, io, T, vt, blk, p.m, p.n, out.ptr, p.lda, ipiv.ptr, want);
            for (ipiv, out[p.a.len..][0..ipiv.len]) |v, *o| o.* = fromF64(T, @floatFromInt(v));
            out[out.len - 1] = fromF64(T, if (k) |z| @floatFromInt(z) else -1);
        }

        /// The public `getrs` and `gesv` (their solutions after their
        /// factors), threading by `blas.config`'s workers.
        pub fn solve(c: @This(), io: Io, want: u32, out: []T) anyerror!void {
            const saved = blas.config.cpu_count;
            defer blas.config.cpu_count = saved;
            blas.config.cpu_count = want;
            const p = c.p;
            const n = p.n;
            const al = testing.allocator;
            const ipiv = try al.alloc(usize, n);
            defer al.free(ipiv);
            const b = out[p.a.len..];
            const nrhs = b.len / (2 * n);
            check.fill(T, b, 4);
            @memcpy(out[0..p.a.len], p.a);
            try gesv(al, io, T, n, nrhs, out[0..p.a.len], p.lda, ipiv, b[0 .. n * nrhs], n);
            try getrs(al, io, T, .c, n, nrhs, out[0..p.a.len], p.lda, ipiv, b[n * nrhs ..], n);
        }
    };
}

fn aheadCases(vt: *const Vtable) anyerror!void {
    _ = blas.coreName();
    return lcheck.eachBlocking(vt, shrink, .{ f32, f64, types.Complex32, types.Complex64 }, aheadOn);
}

fn aheadOn(comptime T: type, v: *const Vtable, vi: usize) anyerror!void {
    _ = v;
    _ = vi;
    const al = testing.allocator;
    for ([_][2]usize{ .{ 130, 130 }, .{ 150, 90 }, .{ 90, 150 }, .{ 41, 41 } }) |s| {
        if (check.rig() and s[0] * s[1] > 70 * 70) continue;
        var p = try fixtures.LuProblem(T).init(al, s[0], s[1], s[0] + 7 * s[1]);
        defer p.deinit(al);
        try lcheck.expectSameBits(T, Threaded(T){ .p = &p }, Threaded(T).run, (Threaded(T){ .p = &p }).len(), "getrf");
        // two zero columns: the first one's pivot, the factorization completed
        for (0..s[0]) |i| for ([_]usize{ 17, 30 }) |j| {
            p.a[i + j * p.lda] = scalar(T, 0);
        };
        try lcheck.expectSameBits(T, Threaded(T){ .p = &p }, Threaded(T).run, (Threaded(T){ .p = &p }).len(), "getrf, singular");
    }
}

test "getrf's lookahead gives the serial bits at any worker count, singular too, on every linked core" {
    try check.eachCore(aheadCases);
}

test "gesv and getrs give one worker's bits at any worker count" {
    _ = blas.coreName();
    const al = testing.allocator;
    inline for (.{ f32, f64, types.Complex32, types.Complex64 }) |T| {
        // past every type's `gesv` threshold
        const n: usize = if (check.rig()) 70 else 210;
        var p = try fixtures.LuProblem(T).init(al, n, n, 77);
        defer p.deinit(al);
        try lcheck.expectSameBits(T, Threaded(T){ .p = &p }, Threaded(T).solve, p.a.len + 2 * n * 5, "gesv, getrs");
    }
}

test "gesv: getrf then getrs, bit for bit; a singular U before B is touched; the first bad argument" {
    const al = testing.allocator;
    inline for (test_types) |T| for ([_]usize{ 1, 7, 33, 70, 130 }) |n| {
        if (check.rig() and n > 70) continue;
        var p = try fixtures.LuProblem(T).init(al, n, n, n + 3);
        defer p.deinit(al);
        const nrhs = 3;
        const b = try al.alloc(T, n * nrhs);
        defer al.free(b);
        check.fill(T, b, n);
        const b2 = try al.dupe(T, b);
        defer al.free(b2);
        const a2 = try al.dupe(T, p.a);
        defer al.free(a2);
        const ipiv2 = try al.alloc(usize, n);
        defer al.free(ipiv2);
        _ = try p.factorIt(al);
        try getrs(al, check.io, T, .n, n, nrhs, p.a, p.lda, p.ipiv, b, n);
        var poison: guard.Poison = .{ .child = al };
        try gesv(poison.allocator(), check.io, T, n, nrhs, a2, p.lda, ipiv2, b2, n);
        try testing.expectEqualSlices(usize, p.ipiv, ipiv2);
        for (p.a, a2) |x, y| try testing.expect(lcheck.sameBits(T, x, y));
        for (b, b2) |x, y| try testing.expect(lcheck.sameBits(T, x, y));
    };
    inline for (.{ f64, types.Complex64 }) |T| {
        // [[1, 2], [2, 4]]: rows swapped, then U(1, 1) = 0
        var a: [4]T = .{ scalar(T, 1), scalar(T, 2), scalar(T, 2), scalar(T, 4) };
        var b: [2]T = .{ scalar(T, 5), scalar(T, 6) };
        var ipiv: [2]usize = @splat(unset);
        try testing.expectError(error.Singular, gesv(al, check.io, T, 2, 1, &a, 2, &ipiv, &b, 2));
        try testing.expectEqual([4]T{ scalar(T, 2), scalar(T, 0.5), scalar(T, 4), scalar(T, 0) }, a);
        try testing.expectEqual([2]T{ scalar(T, 5), scalar(T, 6) }, b);
        // no right-hand side: factored all the same
        a = .{ scalar(T, 1), scalar(T, 3), scalar(T, 2), scalar(T, 4) };
        try gesv(al, check.io, T, 2, 0, &a, 2, &ipiv, &b, 2);
        try testing.expectEqual([2]usize{ 1, 1 }, ipiv);
    }
    const fa = testing.failing_allocator;
    var a: [4]f64 = @splat(nan);
    var b: [4]f64 = @splat(nan);
    var ipiv: [2]usize = @splat(unset);
    try testing.expectError(error.LeadingDimensionTooSmall, gesv(fa, check.io, f64, 2, 1, &a, 1, &ipiv, &b, 2));
    try testing.expectError(error.LeadingDimensionTooSmall, gesv(fa, check.io, f64, 2, 1, &a, 2, &ipiv, &b, 1));
    try testing.expectError(error.SliceTooShort, gesv(fa, check.io, f64, 2, 2, &a, 2, &ipiv, &b, std.math.maxInt(usize) / 2));
    try testing.expectError(error.SliceTooShort, gesv(fa, check.io, f64, 2, 1, &a, 2, ipiv[0..1], &b, 2));
    try gesv(fa, check.io, f64, 0, 2, &a, 1, &ipiv, &b, 1);
    for (a) |v| try testing.expect(std.math.isNan(v));
    try testing.expectEqual([2]usize{ unset, unset }, ipiv);
}

test "getri: A A⁻¹ = I on getrf's factor, every type, sizes past the trmm blocks" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        for ([_]usize{ 1, 2, 3, 5, 16, 33 }) |n| {
            // diagonally dominant, so the inverse's error is a few eps
            const a0 = try al.alloc(T, n * n);
            defer al.free(a0);
            check.fill(T, a0, 7 * n + 1);
            for (0..n) |i| a0[i + i * n] = add(T, a0[i + i * n], scalar(T, @floatFromInt(2 * n)));
            const a = try al.dupe(T, a0);
            defer al.free(a);
            const ipiv = try al.alloc(usize, n);
            defer al.free(ipiv);
            _ = try getrf(al, check.io, T, n, n, a, n, ipiv);
            try getri(al, check.io, T, n, a, n, ipiv);
            const tol = 64 * @as(F, @floatFromInt(n)) * std.math.floatEps(F);
            for (0..n) |j| for (0..n) |i| {
                var s = scalar(T, 0);
                for (0..n) |k| s = add(T, s, mul(T, a0[i + k * n], a[k + j * n]));
                try check.expectClose(T, scalar(T, if (i == j) 1 else 0), s, tol);
            };
        }
    }
}

test "getri: a zero pivot is Singular with A as getrf left it; bad arguments; empty" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64 }) |T| {
        // [[1, 2], [2, 4]]: U(1, 1) = 0
        var a: [4]T = .{ scalar(T, 1), scalar(T, 2), scalar(T, 2), scalar(T, 4) };
        var ipiv: [2]usize = undefined;
        try testing.expectEqual(@as(?usize, 1), try getrf(al, check.io, T, 2, 2, &a, 2, &ipiv));
        const f = a;
        try testing.expectError(error.Singular, getri(al, check.io, T, 2, &a, 2, &ipiv));
        try testing.expectEqual(f, a);
    }
    const fa = testing.failing_allocator;
    var a: [4]f64 = @splat(nan);
    const ipiv: [2]usize = .{ 0, 2 };
    try testing.expectError(error.LeadingDimensionTooSmall, getri(fa, check.io, f64, 2, &a, 1, &ipiv));
    try testing.expectError(error.SliceTooShort, getri(fa, check.io, f64, 3, &a, 3, &ipiv));
    try testing.expectError(error.PivotOutOfRange, getri(fa, check.io, f64, 2, &a, 2, &ipiv));
    try getri(fa, check.io, f64, 0, &a, 1, &ipiv);
    for (a) |v| try testing.expect(std.math.isNan(v));
}
