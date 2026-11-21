//! `mulVec` (CSC and CSR) and `mulDense` against dense loops, and their bits
//! at any worker count.

const std = @import("std");
const matrix = @import("../../../linalg/matrix.zig");
const types = @import("../../blas/types.zig");
const ops = @import("../../lapack/ops.zig");
const blas = @import("../../blas/__root__.zig");
const threading = @import("../../blas/threading.zig");
const check = @import("../../blas/_test/check.zig");
const test_io = @import("../../blas/_test/io.zig");
const sparse = @import("../__root__.zig");
const util = @import("util.zig");

const testing = std.testing;
const Io = std.Io;
const Trans = types.Trans;

const scalars = [_][2][2]f64{
    .{ .{ 1, 0 }, .{ 0, 0 } },
    .{ .{ -0.5, 1.5 }, .{ 2, -1 } },
    .{ .{ 0, 0 }, .{ 0.25, 0.5 } },
    .{ .{ 0, 0 }, .{ 0, 0 } },
};

/// `mulVec` on CSC and CSR for every `trans` and each of `scalars`, on `io`.
fn checkMulVec(comptime T: type, io: Io, p: util.Pair(T), seed: u64) !void {
    const al = testing.allocator;
    const m = p.d.rows;
    const n = p.d.cols;
    const x = try util.vec(T, al, @max(m, n), seed);
    defer al.free(x);
    const y0 = try util.vec(T, al, @max(m, n) + 1, seed + 1);
    defer al.free(y0);
    // NaN where beta zero must not read
    const nan = check.val(T, std.math.nan(f64), 0);
    var r = try p.s.toCompressedRows(al);
    defer r.deinit(al);
    const want = try al.alloc(T, y0.len);
    defer al.free(want);
    const got = try al.alloc(T, y0.len);
    defer al.free(got);
    const got_r = try al.alloc(T, y0.len);
    defer al.free(got_r);
    for (std.meta.tags(Trans)) |t| for (scalars) |s| {
        const alpha = check.val(T, s[0][0], s[0][1]);
        const beta = check.val(T, s[1][0], s[1][1]);
        @memcpy(want, y0);
        if (ops.isZero(T, beta)) @memset(want[0 .. want.len - 1], nan);
        @memcpy(got, want);
        @memcpy(got_r, want);
        util.mulVec(T, p.d, t, alpha, x, beta, want);
        try p.s.mulVec(al, io, t, alpha, x, beta, got);
        try r.mulVec(al, io, t, alpha, x, beta, got_r);
        const tol = util.tol(T, if (t.transposed()) m else n, 2);
        check.expectCloseSlices(T, want, got, tol) catch |e| {
            check.say("csc {s} {d}x{d} trans {s}\n", .{ @typeName(T), m, n, @tagName(t) });
            return e;
        };
        check.expectCloseSlices(T, want, got_r, tol) catch |e| {
            check.say("csr {s} {d}x{d} trans {s}\n", .{ @typeName(T), m, n, @tagName(t) });
            return e;
        };
    };
}

test "mulVec: every trans, alpha and beta against dense loops" {
    const al = testing.allocator;
    inline for (util.all_types) |T| {
        for ([_][2]usize{ .{ 1, 1 }, .{ 17, 5 }, .{ 6, 23 }, .{ 40, 40 } }, 0..) |s, k| {
            var p = try util.random(T, al, s[0], s[1], 350, k + 20);
            defer p.deinit(al);
            try checkMulVec(T, check.io, p, k);
        }
    }
}

test "mulVec: threaded, chunked, on every test Io" {
    const al = testing.allocator;
    _ = blas.coreName();
    const saved = threading.test_serial_max;
    defer threading.test_serial_max = saved;
    threading.test_serial_max = 64;
    var p = try util.random(f64, al, 60, 300, 200, 4);
    defer p.deinit(al);
    var cp = try util.random(types.Complex32, al, 50, 200, 250, 5);
    defer cp.deinit(al);
    const S = struct {
        p: *const util.Pair(f64),
        cp: *const util.Pair(types.Complex32),

        fn run(s: @This(), io: Io) !void {
            try checkMulVec(f64, io, s.p.*, 1);
            try checkMulVec(types.Complex32, io, s.cp.*, 2);
        }
    };
    try test_io.eachTestIo(S{ .p = &p, .cp = &cp }, S.run);
}

test "mulVec: dimension checks and empty shapes" {
    const al = testing.allocator;
    var p = try util.random(f64, al, 4, 3, 500, 1);
    defer p.deinit(al);
    var x: [4]f64 = @splat(1);
    var y: [4]f64 = @splat(1);
    try testing.expectError(error.DimensionMismatch, p.s.mulVec(al, check.io, .n, 1, x[0..2], 0, &y));
    try testing.expectError(error.DimensionMismatch, p.s.mulVec(al, check.io, .n, 1, &x, 0, y[0..3]));
    try testing.expectError(error.DimensionMismatch, p.s.mulVec(al, check.io, .t, 1, x[0..3], 0, &y));
    try testing.expectError(error.DimensionMismatch, p.s.mulVec(al, check.io, .c, 1, &x, 0, y[0..2]));
    // no columns: y := beta y; no rows: nothing
    var e = try sparse.Triplets(f64).init(4, 0);
    defer e.deinit(al);
    var a = try e.toCompressedColumns(al);
    defer a.deinit(al);
    try a.mulVec(al, check.io, .n, 1, x[0..0], 2, &y);
    try testing.expectEqualSlices(f64, &.{ 2, 2, 2, 2 }, &y);
    try a.mulVec(al, check.io, .t, 1, &x, 2, y[0..0]);
    var at = try a.transpose(al);
    defer at.deinit(al);
    try at.mulVec(al, check.io, .n, 1, &x, 0, y[0..0]);
    try at.mulVec(al, check.io, .t, 3, x[0..0], 0, &y);
    try testing.expectEqualSlices(f64, &.{ 0, 0, 0, 0 }, &y);
}

test "mulDense against dense loops" {
    const al = testing.allocator;
    inline for (util.all_types) |T| {
        for ([_][3]usize{ .{ 1, 1, 1 }, .{ 12, 7, 5 }, .{ 5, 9, 0 }, .{ 0, 4, 3 }, .{ 30, 25, 4 } }, 0..) |s, k| {
            var p = try util.random(T, al, s[0], s[1], 300, k + 40);
            defer p.deinit(al);
            const b: matrix.Mat(T) = try .init(al, s[1], s[2]);
            defer b.deinit(al);
            check.fill(T, b.data, k);
            const want = try util.mulDense(T, al, p.d, b);
            defer want.deinit(al);
            // C at a leading dimension one past its rows, stale values in it
            const ldc = s[0] + 1;
            const got = try util.vec(T, al, ldc * s[2], k + 90);
            defer al.free(got);
            const pad = try al.dupe(T, got);
            defer al.free(pad);
            try p.s.mulDense(al, check.io, b.data, @max(s[1], 1), s[2], got, ldc);
            for (0..s[2]) |j| {
                try check.expectCloseSlices(T, want.data[j * s[0] ..][0..s[0]], got[j * ldc ..][0..s[0]], util.tol(T, s[1], 1));
                try testing.expectEqual(pad[j * ldc + s[0]], got[j * ldc + s[0]]);
            }
        }
        var p = try util.random(T, al, 3, 4, 500, 9);
        defer p.deinit(al);
        var b: [8]T = undefined;
        var c: [6]T = undefined;
        try testing.expectError(error.DimensionMismatch, p.s.mulDense(al, check.io, &b, 3, 2, &c, 3));
        try testing.expectError(error.DimensionMismatch, p.s.mulDense(al, check.io, &b, 4, 2, c[0..5], 3));
        try testing.expectError(error.DimensionMismatch, p.s.mulDense(al, check.io, &b, 4, 2, &c, 2));
    }
}

test "mulVec and mulDense: one worker's bits at any worker count" {
    if (comptime !test_io.hosted) return error.SkipZigTest;
    _ = blas.coreName();
    const al = testing.allocator;
    var runs: test_io.WorkerRuns = .init(al);
    defer runs.deinit();
    inline for (.{ f64, types.Complex32 }) |T| {
        const m = 70;
        const n = 900;
        var p = try util.random(T, al, m, n, 150, 77);
        defer p.deinit(al);
        const x = try util.vec(T, al, n, 6);
        defer al.free(x);
        const y0 = try util.vec(T, al, n, 7);
        defer al.free(y0);
        const b: matrix.Mat(T) = try .init(al, n, 3);
        defer b.deinit(al);
        check.fill(T, b.data, 8);
        const alpha = check.val(T, 0.75, 0.5);
        const beta = check.val(T, -1.5, 0.25);
        var outs: [test_io.worker_runs.len][]T = undefined;
        for (&outs) |*o| o.* = try al.alloc(T, m + n + m * 3);
        defer for (outs) |o| al.free(o);
        for (0..test_io.worker_runs.len) |c| {
            const io = runs.io(c);
            // the same chunks every run, many of them
            threading.test_serial_max = 64;
            const o = outs[c];
            @memcpy(o[0..m], y0[0..m]);
            try p.s.mulVec(al, io, .n, alpha, x, beta, o[0..m]);
            @memcpy(o[m..][0..n], y0);
            try p.s.mulVec(al, io, .c, alpha, x, beta, o[m..][0..n]);
            try p.s.mulDense(al, io, b.data, n, 3, o[m + n ..], m);
        }
        try test_io.expectSameBytes(outs);
        // and right
        const want = try al.dupe(T, y0[0..m]);
        defer al.free(want);
        util.mulVec(T, p.d, .n, alpha, x, beta, want);
        try check.expectCloseSlices(T, want, outs[0][0..m], util.tol(T, n, 2));
    }
}
