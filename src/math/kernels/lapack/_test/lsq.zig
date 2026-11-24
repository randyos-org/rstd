//! Tests of `../lsq.zig`.

const std = @import("std");
const blas = @import("../../blas/__root__.zig");
const types = @import("../../blas/types.zig");
const ops = @import("../ops.zig");
const lroot = @import("../lsq.zig");

const Allocator = std.mem.Allocator;
const Mat = ops.Mat;
const Float = types.Float;
const gelsy = lroot.gelsy;

// ---- tests -----------------------------------------------------------------------

const testing = std.testing;
const check = @import("../../blas/_test/check.zig");
const lcheck = @import("check.zig");
const test_io = @import("../../blas/_test/io.zig");
const svd = @import("../svd.zig");

const test_types = lcheck.test_types;

const thresh = lcheck.thresh;

/// A random `m x n` of rank `r` scaled by `s`, and `B` (`max(m, n) x nrhs`)
/// scaled by `s` too, so `X` does not move.
fn problem(comptime T: type, al: Allocator, m: usize, n: usize, r: usize, nrhs: usize, s: Float(T), seed: u64) !struct { a: Mat(T), b: Mat(T) } {
    const a = try ops.alloc(T, al, m, n);
    errdefer ops.free(T, al, a);
    try lcheck.rankProduct(T, al, a, r, seed);
    for (0..n) |j| for (0..m) |i| {
        a.at(i, j).* = ops.mulReal(T, a.get(i, j), s);
    };
    const b = try ops.alloc(T, al, @max(m, n), nrhs);
    check.fill(T, b.ptr[0 .. b.ld * nrhs], seed + 2);
    for (0..nrhs) |j| for (b.col(j)[0..b.rows]) |*v| {
        v.* = ops.mulReal(T, v.*, s);
    };
    return .{ .a = a, .b = b };
}

/// `gelsy` against the SVD's minimum-norm solution (`svd.zig`) cut at
/// rank `r`: the rank exactly, `X` within `thresh n eps cond(A) ‖X‖`.
fn gelsyCase(comptime T: type, al: Allocator, m: usize, n: usize, r: usize, nrhs: usize, s: Float(T), seed: u64) !void {
    const F = Float(T);
    const eps = std.math.floatEps(F);
    const p = try problem(T, al, m, n, r, nrhs, s, seed);
    defer ops.free(T, al, p.a);
    defer ops.free(T, al, p.b);
    // the reference solution from the unscaled SVD
    const k = @min(m, n);
    const sa = try ops.alloc(T, al, m, n);
    defer ops.free(T, al, sa);
    for (0..n) |j| for (0..m) |i| {
        sa.at(i, j).* = ops.mulReal(T, p.a.get(i, j), 1 / s);
    };
    const sv = try al.alloc(F, k);
    defer al.free(sv);
    const u = try ops.alloc(T, al, m, m);
    defer ops.free(T, al, u);
    const vt = try ops.alloc(T, al, n, n);
    defer ops.free(T, al, vt);
    try svd.gesvd(al, check.io, T, .all, .all, m, n, sa.slice(), m, sv, u.slice(), m, vt.slice(), n, null);
    const want = try ops.alloc(T, al, n, nrhs);
    defer ops.free(T, al, want);
    want.fill(ops.zero(T));
    var wmax: F = 0;
    for (0..nrhs) |c| {
        for (0..r) |q| {
            // (u_qᴴ b / s_q) v_q, b unscaled
            var d = ops.zero(T);
            for (0..m) |i| d = ops.add(T, d, ops.mul(T, ops.conj(T, u.get(i, q)), ops.mulReal(T, p.b.get(i, c), 1 / s)));
            d = ops.divReal(T, d, sv[q]);
            for (0..n) |j| want.at(j, c).* = ops.add(T, want.get(j, c), ops.mul(T, d, ops.conj(T, vt.get(q, j))));
        }
        for (0..n) |j| wmax = @max(wmax, ops.abs(T, want.get(j, c)));
    }
    const jpvt = try al.alloc(usize, n);
    defer al.free(jpvt);
    @memset(jpvt, 0);
    const rcond = 1000 * eps;
    const rank = try gelsy(al, check.io, T, m, n, nrhs, p.a.slice(), p.a.ld, p.b.slice(), p.b.ld, jpvt, rcond);
    try testing.expectEqual(r, rank);
    const cond = if (r > 0) sv[0] / sv[r - 1] else 1;
    const tol = thresh * @as(F, @floatFromInt(@max(m, n))) * eps * cond * @max(wmax, 1);
    for (0..nrhs) |c| for (0..n) |j| {
        const d = ops.abs(T, ops.sub(T, p.b.get(j, c), want.get(j, c)));
        if (!(d <= tol)) {
            check.say("gelsy {s} {d} x {d} rank {d} scale {e}: X({d}, {d}) off by {e}, tolerance {e}\n", .{ @typeName(T), m, n, r, @as(f64, @floatCast(s)), j, c, @as(f64, @floatCast(d)), @as(f64, @floatCast(tol)) });
            return error.TestUnexpectedResult;
        }
    };
}

test "lapack gelsy: the minimum-norm solution and rank, full and deficient rank, tall and wide, scaled near the range's ends" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        const cases = [_][4]usize{ .{ 1, 1, 1, 1 }, .{ 6, 4, 4, 2 }, .{ 4, 6, 4, 1 }, .{ 8, 6, 3, 2 }, .{ 6, 8, 3, 3 }, .{ 20, 15, 10, 4 }, .{ 40, 40, 25, 2 } };
        for (cases, 0..) |c, i| {
            try gelsyCase(T, al, c[0], c[1], c[2], c[3], 1, i * 11 + 1);
        }
        // past smlnum and bignum, so A and B are scaled and back; normal
        // entries, so the rank survives. Not at an extended type, whose
        // `nrm2` is `nrm2.S`'s unscaled sum of squares and underflows there.
        if (comptime !types.isExtended(F)) {
            try gelsyCase(T, al, 8, 6, 3, 2, std.math.floatMin(F) / std.math.floatEps(F) / 64, 3);
            try gelsyCase(T, al, 6, 8, 4, 2, std.math.floatMax(F) / 4096, 5);
        }
    }
}

test "lapack gelsy: A and B both past BIGNUM, a small column of X is not flushed" {
    const al = testing.allocator;
    // `nrm2` at an extended type is unscaled, as above
    inline for (test_types) |T| if (comptime !types.isExtended(Float(T))) {
        const F = Float(T);
        const tm = std.math.floatTrueMin(F);
        const bignum = std.math.floatEps(F) / std.math.floatMin(F);
        // X = B / a: 2^12, and 2^16 tm where the scaled solve has 16 tm,
        // which A's step back alone (2^-8) flushes
        const av = bignum * 256;
        const bv = av * 4096;
        const tv = 16 * tm * bv;
        var a = [_]T{ ops.scalar(T, av), ops.zero(T), ops.zero(T), ops.zero(T), ops.scalar(T, av), ops.zero(T) };
        var b = [_]T{ ops.scalar(T, bv), ops.scalar(T, bv), ops.zero(T), ops.scalar(T, tv), ops.scalar(T, tv), ops.zero(T) };
        var jpvt = [_]usize{ 0, 0 };
        try testing.expectEqual(2, try gelsy(al, check.io, T, 3, 2, 2, &a, 3, &b, 3, &jpvt, 0.01));
        for ([_]usize{ 0, 1 }) |i| {
            try testing.expect(ops.abs(T, ops.sub(T, b[i], ops.scalar(T, 4096))) <= 4096 * thresh * std.math.floatEps(F));
            try testing.expect(ops.abs(T, ops.sub(T, b[3 + i], ops.scalar(T, 65536 * tm))) <= 256 * tm);
            // R, scaled back
            try testing.expect(@abs(ops.abs(T, a[i * 4]) - av) <= av * thresh * std.math.floatEps(F));
        }
    };
}

test "lapack gelsy: a zero A, a fixed column, empty, bad arguments" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        var a: [6]T = @splat(ops.zero(T));
        var b: [3]T = @splat(ops.one(T));
        var jpvt = [_]usize{ 0, 0 };
        try testing.expectEqual(0, try gelsy(al, check.io, T, 3, 2, 1, &a, 3, &b, 3, &jpvt, 0.01));
        for (b) |v| try testing.expect(ops.isZero(T, v));
        // column 1 fixed first; the solution of a square, well-conditioned A
        var c = [_]T{ ops.scalar(T, 4), ops.scalar(T, 1), ops.scalar(T, 1), ops.scalar(T, 3) };
        var d = [_]T{ ops.scalar(T, 5), ops.scalar(T, 4) };
        jpvt = .{ 0, 1 };
        try testing.expectEqual(2, try gelsy(al, check.io, T, 2, 2, 1, &c, 2, &d, 2, &jpvt, 0.01));
        try testing.expectEqualSlices(usize, &.{ 1, 0 }, &jpvt);
        const tol = 16 * std.math.floatEps(Float(T));
        try testing.expect(ops.abs(T, ops.sub(T, d[0], ops.one(T))) <= tol);
        try testing.expect(ops.abs(T, ops.sub(T, d[1], ops.one(T))) <= tol);
        try testing.expectEqual(0, try gelsy(al, check.io, T, 0, 2, 1, &a, 1, &b, 2, &jpvt, 0.01));
    }
    var a: [4]f64 = undefined;
    var jp: [2]usize = undefined;
    try testing.expectError(error.LeadingDimensionTooSmall, gelsy(al, check.io, f64, 2, 1, 1, &a, 1, &a, 2, &jp, 0));
    try testing.expectError(error.LeadingDimensionTooSmall, gelsy(al, check.io, f64, 1, 2, 1, &a, 1, &a, 1, &jp, 0));
    try testing.expectError(error.SliceTooShort, gelsy(al, check.io, f64, 2, 3, 1, &a, 2, &a, 3, &jp, 0));
    try testing.expectError(error.SliceTooShort, gelsy(al, check.io, f64, 1, 3, 1, &a, 1, &a, 3, &jp, 0));
}

test "lapack gelsy: one worker's bits at any worker count" {
    if (comptime !test_io.hosted) return error.SkipZigTest;
    _ = blas.coreName();
    const al = testing.allocator;
    var runs: test_io.WorkerRuns = .init(al);
    defer runs.deinit();
    inline for (.{ f64, types.Complex32 }) |T| {
        const m = 300;
        const n = 90;
        const nrhs = 70;
        const p = try problem(T, al, m, n, 60, nrhs, 1, 21);
        defer ops.free(T, al, p.a);
        defer ops.free(T, al, p.b);
        var outs: [3][]T = undefined;
        for (&outs) |*o| o.* = try al.alloc(T, m * n + m * nrhs);
        defer for (outs) |o| al.free(o);
        for (0..test_io.worker_runs.len) |c| {
            const io = runs.io(c);
            @memcpy(outs[c][0 .. m * n], p.a.ptr[0 .. m * n]);
            @memcpy(outs[c][m * n ..], p.b.ptr[0 .. m * nrhs]);
            var jpvt: [n]usize = @splat(0);
            _ = try gelsy(al, io, T, m, n, nrhs, outs[c][0 .. m * n], m, outs[c][m * n ..], m, &jpvt, 1e-5);
        }
        try test_io.expectSameBytes(outs);
    }
}
