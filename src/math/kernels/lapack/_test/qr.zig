//! Tests of `../qr.zig`.

const std = @import("std");
const ops = @import("../ops.zig");
const blas_calls = @import("../blas_calls.zig");
const auxiliary = @import("../auxiliary.zig");
const hh = @import("../householder.zig");
const types = @import("../../blas/types.zig");
const qroot = @import("../qr.zig");

const Allocator = std.mem.Allocator;
const Mat = ops.Mat;
const Float = ops.Float;
const Trans = types.Trans;
const Side = types.Side;
const blockSize = qroot.blockSize;
const factor = qroot.factor;
const geqrf = qroot.geqrf;
const orgqr = qroot.orgqr;
const ormqr = qroot.ormqr;

// ---- tests -----------------------------------------------------------------------

const testing = std.testing;
const check = @import("../../blas/_test/check.zig");
const lcheck = @import("check.zig");
const blas = @import("../../blas/__root__.zig");
const test_io = @import("../../blas/_test/io.zig");

/// `Q R` against `A`, within `n eps ‖A‖` per element (faer's `test_qr` takes
/// `1e-10` at `c64`).
fn expectQr(comptime T: type, al: Allocator, a0: Mat(T), qr: Mat(T), h: Mat(T)) !void {
    const F = Float(T);
    const m = a0.rows;
    const n = a0.cols;
    const q = try ops.alloc(T, al, m, m);
    defer ops.free(T, al, q);
    q.identity();
    try hh.applySequence(al, check.io, T, qr, h, false, q);
    var amax: F = 0;
    for (0..n) |j| for (0..m) |i| {
        amax = @max(amax, ops.abs(T, a0.get(i, j)));
    };
    const tol = @as(F, @floatFromInt(4 * (m + n))) * std.math.floatEps(F) * @max(amax, 1);
    for (0..n) |j| for (0..m) |i| {
        var s = ops.zero(T);
        for (0..@min(j + 1, m)) |k| s = ops.add(T, s, ops.mul(T, q.get(i, k), qr.get(k, j)));
        const err = ops.abs(T, ops.sub(T, s, a0.get(i, j)));
        if (!(err <= tol)) {
            check.say("qr {s} {d}x{d}: ({d},{d}) off by {e}, tol {e}\n", .{ @typeName(T), m, n, i, j, err, tol });
            return error.TestUnexpectedResult;
        }
    };
}

test "faer qr: Q R = A, square, tall, wide, blocked and rank deficient" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64, f32, types.Complex32 }) |T| {
        for ([_][3]usize{ .{ 1, 1, 0 }, .{ 2, 2, 0 }, .{ 8, 8, 0 }, .{ 40, 17, 0 }, .{ 17, 40, 0 }, .{ 70, 70, 0 }, .{ 130, 60, 0 }, .{ 64, 64, 3 }, .{ 80, 50, 1 } }) |c| {
            const m = c[0];
            const n = c[1];
            const a0 = try ops.alloc(T, al, m, n);
            defer ops.free(T, al, a0);
            if (c[2] == 0) {
                check.fill(T, a0.ptr[0 .. m * n], m * 31 + n);
            } else try lcheck.rankProduct(T, al, a0, c[2], 3);
            const qr = try ops.alloc(T, al, m, n);
            defer ops.free(T, al, qr);
            qr.copyFrom(a0);
            const h = try ops.alloc(T, al, blockSize(m, n), @min(m, n));
            defer ops.free(T, al, h);
            const rank = try factor(al, check.io, T, qr, h);
            try testing.expect(rank >= @min(c[2], @min(m, n)));
            try expectQr(T, al, a0, qr, h);
        }
    }
}

const test_types = .{ f32, f64, types.Complex32, types.Complex64 };

/// `expectQr`'s tolerance, `4 (m + n) eps max(1, ‖A‖max)` per element.
fn qrTol(comptime F: type, m: usize, n: usize, amax: F) F {
    return @as(F, @floatFromInt(4 * (m + n))) * std.math.floatEps(F) * @max(amax, 1);
}

test "faer qr: geqrf and orgqr, Q R = A and Qᴴ Q = I in LAPACK's layout" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        // m, n, rank (0 full, else a product; n + 1 a zero column)
        for ([_][3]usize{ .{ 1, 1, 0 }, .{ 5, 3, 0 }, .{ 3, 5, 0 }, .{ 40, 17, 0 }, .{ 17, 40, 0 }, .{ 70, 70, 0 }, .{ 130, 60, 0 }, .{ 64, 64, 3 }, .{ 30, 20, 21 } }) |c| {
            const m = c[0];
            const n = c[1];
            const k = @min(m, n);
            const a0 = try ops.alloc(T, al, m, n);
            defer ops.free(T, al, a0);
            if (c[2] == n + 1) {
                check.fill(T, a0.ptr[0 .. m * n], 9);
                @memset(a0.col(4)[0..m], ops.zero(T));
            } else if (c[2] == 0) {
                check.fill(T, a0.ptr[0 .. m * n], m + 7 * n);
            } else try lcheck.rankProduct(T, al, a0, c[2], m + 7 * n);
            const lda = m + 1;
            const a = try ops.alloc(T, al, lda, n);
            defer ops.free(T, al, a);
            a.fill(lcheck.poisoned(T));
            a.sub(0, 0, m, n).copyFrom(a0);
            const tau = try al.alloc(T, k);
            defer al.free(tau);
            try geqrf(al, check.io, T, m, n, a.slice(), lda, tau);
            for (0..n) |j| try testing.expect(lcheck.isNanT(T, a.get(m, j)));
            // R, then Q's first k columns
            const r = try ops.alloc(T, al, k, n);
            defer ops.free(T, al, r);
            for (0..n) |j| for (0..k) |i| {
                r.at(i, j).* = if (i <= j) a.get(i, j) else ops.zero(T);
            };
            try orgqr(al, check.io, T, m, k, k, a.slice(), lda, tau);
            const q = a.sub(0, 0, m, k);
            const qr_ = try ops.alloc(T, al, m, n);
            defer ops.free(T, al, qr_);
            try blas_calls.gemm(al, check.io, T, .n, .n, ops.one(T), q, r, ops.zero(T), qr_);
            try lcheck.expectNear(T, "geqrf Q R", qr_, a0, qrTol(F, m, n, auxiliary.normMax(T, a0)));
            const g = try ops.alloc(T, al, k, k);
            defer ops.free(T, al, g);
            try blas_calls.gemm(al, check.io, T, .c, .n, ops.one(T), q, q, ops.zero(T), g);
            const id = try ops.alloc(T, al, k, k);
            defer ops.free(T, al, id);
            id.identity();
            try lcheck.expectNear(T, "orgqr Qh Q", g, id, qrTol(F, m, n, 1));
            if (c[2] == n + 1) try testing.expect(ops.isZero(T, r.get(4, 4)));
        }
    }
}

/// `C := op(H(0) ... H(k-1)) C` or `C op(...)`, `H = I - tau v vᴴ`, one
/// reflection at a time, `v` unit lower in `v`'s columns.
fn refApply(comptime T: type, side: Side, adjoint: bool, v: Mat(T), tau: []const T, c: Mat(T)) void {
    const k = tau.len;
    const left = side == .l;
    // left: Q C applies the last first; Qᴴ C the first first; right mirrors
    for (0..k) |s| {
        const j = if (left != adjoint) k - 1 - s else s;
        const t = if (adjoint) ops.conj(T, tau[j]) else tau[j];
        const nq = v.rows;
        const others = if (left) c.cols else c.rows;
        for (0..others) |o| {
            // d := vᴴ x (left, x a column) or x v (right, x a row)
            var d = ops.zero(T);
            for (j..nq) |i| {
                const vi = if (i == j) ops.one(T) else v.get(i, j);
                d = if (left) ops.add(T, d, ops.mul(T, ops.conj(T, vi), c.get(i, o))) else ops.add(T, d, ops.mul(T, c.get(o, i), vi));
            }
            d = ops.mul(T, t, d);
            for (j..nq) |i| {
                const vi = if (i == j) ops.one(T) else v.get(i, j);
                if (left) {
                    c.at(i, o).* = ops.sub(T, c.get(i, o), ops.mul(T, vi, d));
                } else {
                    c.at(o, i).* = ops.sub(T, c.get(o, i), ops.mul(T, d, ops.conj(T, vi)));
                }
            }
        }
    }
}

test "faer qr: ormqr and orgqr against one reflection at a time, faer's tau and the reference's" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        const nq = 37;
        const k = 20;
        const other = 9;
        for (0..2) |kind| {
            const v = try ops.alloc(T, al, nq, nq);
            defer ops.free(T, al, v);
            const tau = try al.alloc(T, k);
            defer al.free(tau);
            if (kind == 0) {
                // geqrf's
                check.fill(T, v.ptr[0 .. nq * nq], 21);
                try geqrf(al, check.io, T, nq, k, v.slice(), nq, tau);
            } else {
                // the reference's form: `tau = (1 - e^(i phi)) / ‖v‖²`
                // (`xLARFG`'s `(beta - alpha) / beta`), real 2 / ‖v‖²; one 0
                check.fill(T, v.ptr[0 .. nq * nq], 22);
                var prng = std.Random.DefaultPrng.init(23);
                for (tau, 0..) |*t, j| {
                    var nv: F = 1;
                    for (j + 1..nq) |i| nv += ops.abs2(T, v.get(i, j));
                    const phi = prng.random().float(F) * 2 * std.math.pi;
                    t.* = if (comptime ops.isComplex(T)) .init((1 - @cos(phi)) / nv, -@sin(phi) / nv) else 2 / nv;
                }
                tau[5] = ops.zero(T);
            }
            for ([_]Side{ .l, .r }) |side| for ([_]bool{ false, true }) |adj| {
                const rows: usize = if (side == .l) nq else other;
                const cols: usize = if (side == .l) other else nq;
                const c = try ops.alloc(T, al, rows, cols);
                defer ops.free(T, al, c);
                check.fill(T, c.ptr[0 .. rows * cols], 24);
                const want = try ops.alloc(T, al, rows, cols);
                defer ops.free(T, al, want);
                want.copyFrom(c);
                refApply(T, side, adj, v.sub(0, 0, nq, k), tau, want);
                const tr: Trans = if (adj) .c else .n;
                try ormqr(al, check.io, T, side, tr, rows, cols, k, v.slice(), nq, tau, c.slice(), rows);
                try lcheck.expectNear(T, "ormqr", c, want, qrTol(F, nq, nq, auxiliary.normMax(T, want)));
            };
            // orgqr: Q's first n columns
            const n = 30;
            const q = try ops.alloc(T, al, nq, nq);
            defer ops.free(T, al, q);
            q.copyFrom(v);
            try orgqr(al, check.io, T, nq, n, k, q.slice(), nq, tau);
            const want = try ops.alloc(T, al, nq, n);
            defer ops.free(T, al, want);
            want.identity();
            refApply(T, .l, false, v.sub(0, 0, nq, k), tau, want);
            try lcheck.expectNear(T, "orgqr", q.sub(0, 0, nq, n), want, qrTol(F, nq, nq, 1));
        }
    }
}

test "faer qr: geqrf, orgqr and ormqr report the first bad argument and read nothing when empty" {
    const fa = testing.failing_allocator;
    const nan = std.math.nan(f64);
    var a: [4]f64 = @splat(nan);
    var tau: [2]f64 = @splat(nan);
    var c: [4]f64 = @splat(nan);
    try testing.expectError(error.LeadingDimensionTooSmall, geqrf(fa, check.io, f64, 2, 2, &a, 1, &tau));
    try testing.expectError(error.SliceTooShort, geqrf(fa, check.io, f64, 3, 2, &a, 3, &tau));
    try testing.expectError(error.SliceTooShort, geqrf(fa, check.io, f64, 2, 2, &a, 2, tau[0..1]));
    try testing.expectError(error.TooManyColumns, orgqr(fa, check.io, f64, 1, 2, 0, &a, 1, &tau));
    try testing.expectError(error.TooManyReflections, orgqr(fa, check.io, f64, 2, 1, 2, &a, 2, &tau));
    try testing.expectError(error.LeadingDimensionTooSmall, orgqr(fa, check.io, f64, 2, 2, 2, &a, 1, &tau));
    try testing.expectError(error.SliceTooShort, orgqr(fa, check.io, f64, 2, 2, 2, &a, 2, tau[0..1]));
    try testing.expectError(error.UnsupportedTranspose, ormqr(fa, check.io, f64, .l, .r, 2, 2, 1, &a, 2, &tau, &c, 2));
    try testing.expectError(error.UnsupportedTranspose, ormqr(fa, check.io, types.Complex64, .l, .t, 1, 1, 1, @ptrCast(&a), 1, @ptrCast(&tau), @ptrCast(&c), 1));
    try testing.expectError(error.TooManyReflections, ormqr(fa, check.io, f64, .r, .n, 2, 1, 2, &a, 1, &tau, &c, 2));
    try testing.expectError(error.LeadingDimensionTooSmall, ormqr(fa, check.io, f64, .l, .n, 2, 1, 1, &a, 1, &tau, &c, 2));
    try testing.expectError(error.LeadingDimensionTooSmall, ormqr(fa, check.io, f64, .r, .n, 2, 1, 1, &a, 1, &tau, &c, 1));
    try testing.expectError(error.SliceTooShort, ormqr(fa, check.io, f64, .l, .n, 2, 3, 1, &a, 2, &tau, &c, 2));
    // empty: nothing read or written
    try geqrf(fa, check.io, f64, 0, 2, &a, 1, &tau);
    try orgqr(fa, check.io, f64, 2, 0, 0, &a, 2, &tau);
    try ormqr(fa, check.io, f64, .l, .n, 2, 2, 0, &a, 2, &tau, &c, 2);
    for ([_][]const f64{ &a, &tau, &c }) |x| for (x) |e| try testing.expect(std.math.isNan(e));
}

test "faer qr: geqrf and ormqr give one worker's bits at any worker count" {
    if (comptime !test_io.hosted) return error.SkipZigTest;
    _ = blas.coreName();
    const al = testing.allocator;
    var runs: test_io.WorkerRuns = .init(al);
    defer runs.deinit();
    inline for (.{ f64, types.Complex32 }) |T| {
        const m = 300;
        const n = 140;
        const a0 = try ops.alloc(T, al, m, n);
        defer ops.free(T, al, a0);
        check.fill(T, a0.ptr[0 .. m * n], 77);
        var outs: [3][]T = undefined;
        for (&outs) |*o| o.* = try al.alloc(T, m * n + n + m * 40);
        defer for (outs) |o| al.free(o);
        for (0..test_io.worker_runs.len) |ci| {
            const io = runs.io(ci);
            const o = outs[ci];
            const a: Mat(T) = .{ .ptr = o.ptr, .rows = m, .cols = n, .ld = m };
            a.copyFrom(a0);
            try geqrf(al, io, T, m, n, a.slice(), m, o[m * n ..][0..n]);
            const c: Mat(T) = .{ .ptr = o.ptr + m * n + n, .rows = m, .cols = 40, .ld = m };
            check.fill(T, c.ptr[0 .. m * 40], 78);
            try ormqr(al, io, T, .l, .c, m, 40, n, a.slice(), m, o[m * n ..][0..n], c.slice(), m);
        }
        try test_io.expectSameBytes(outs);
    }
}
