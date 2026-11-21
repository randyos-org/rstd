//! Tests of `evd.zig`: `syevd`, `syev`, `geev`, `hseqr` and `gees`.

const std = @import("std");
const types = @import("../../blas/types.zig");
const blas = @import("../../blas/__root__.zig");
const test_io = @import("../../blas/_test/io.zig");
const check = @import("../../blas/_test/check.zig");
const lcheck = @import("check.zig");
const ops = @import("../ops.zig");
const blas_calls = @import("../blas_calls.zig");
const auxiliary = @import("../auxiliary.zig");
const evd = @import("../evd.zig");

const testing = std.testing;
const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = ops.Mat;
const Float = ops.Float;
const isComplex = ops.isComplex;
const Uplo = types.Uplo;
const VectorJob = evd.VectorJob;

pub const test_types = .{ f32, f64, types.Complex32, types.Complex64 };

/// Faer's `test_self_adjoint_evd` tolerance: `ApproxEq::eps()` (`128 eps`)
/// times `sqrt(8 n)`, against `max(1, ‖A‖max)` per element.
pub fn tolerance(comptime F: type, n: usize, amax: F) F {
    return 128 * std.math.floatEps(F) * @sqrt(8 * @as(F, @floatFromInt(n))) * @max(amax, 1);
}

/// `‖A V - V Λ‖`, `‖Vᴴ V - I‖` per element within `tolerance`, `w`
/// nondecreasing.
fn expectSym(comptime T: type, al: Allocator, a0: Mat(T), w: []const Float(T), v: Mat(T), what: []const u8) !void {
    const F = Float(T);
    const n = a0.rows;
    for (1..n) |i| if (!(w[i - 1] <= w[i])) {
        check.say("{s} {s} n={d}: w out of order at {d}\n", .{ what, @typeName(T), n, i });
        return error.TestUnexpectedResult;
    };
    const r = try ops.alloc(T, al, n, n);
    defer ops.free(T, al, r);
    for (0..n) |j| for (0..n) |i| {
        r.at(i, j).* = ops.mulReal(T, v.get(i, j), w[j]);
    };
    try blas_calls.gemm(al, check.io, T, .n, .n, ops.neg(T, ops.one(T)), a0, v, ops.one(T), r);
    const tol = tolerance(F, n, auxiliary.normMax(T, a0));
    if (!(auxiliary.normMax(T, r) <= tol)) {
        check.say("{s} {s} n={d}: A V - V L = {e}, tol {e}\n", .{ what, @typeName(T), n, auxiliary.normMax(T, r), tol });
        return error.TestUnexpectedResult;
    }
    try blas_calls.gemm(al, check.io, T, .c, .n, ops.one(T), v, v, ops.zero(T), r);
    for (0..n) |i| r.at(i, i).* = ops.sub(T, r.get(i, i), ops.one(T));
    if (!(auxiliary.normMax(T, r) <= tolerance(F, n, 1))) {
        check.say("{s} {s} n={d}: Vh V - I = {e}\n", .{ what, @typeName(T), n, auxiliary.normMax(T, r) });
        return error.TestUnexpectedResult;
    }
}

const Routine = enum { syevd, syev };

/// The routine on a copy of `a0` (padded `lda`), both triangles hidden
/// behind NaN but the one `uplo` names; the vectors into `v`.
fn decompose(comptime T: type, al: Allocator, io: Io, r: Routine, jobz: VectorJob, uplo: Uplo, a0: Mat(T), w: []Float(T), v: ?Mat(T)) !void {
    const n = a0.rows;
    const lda = n + 1;
    const a = try ops.alloc(T, al, lda, n);
    defer ops.free(T, al, a);
    a.fill(lcheck.poisoned(T));
    for (0..n) |j| for (0..n) |i| {
        if (if (uplo == .l) i >= j else i <= j) a.at(i, j).* = a0.get(i, j);
    };
    switch (r) {
        .syevd => try evd.syevd(al, io, T, jobz, uplo, n, a.slice(), lda, w),
        .syev => try evd.syev(al, io, T, jobz, uplo, n, a.slice(), lda, w),
    }
    for (0..n) |j| try testing.expect(lcheck.isNanT(T, a.get(n, j)));
    if (v) |x| x.copyFrom(a.sub(0, 0, n, n));
}

test "faer evd: syevd and syev, A V = V L with V orthonormal, both triangles, faer's sizes and past 128" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        for ([_]usize{ 1, 2, 4, 10, 15, 20, 50, 130, 150 }) |n| {
            const a0 = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, a0);
            lcheck.hermitian(T, a0, 7 * n + 1);
            const w = try al.alloc(F, n);
            defer al.free(w);
            const w2 = try al.alloc(F, n);
            defer al.free(w2);
            const v = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, v);
            for ([_]Routine{ .syevd, .syev }) |r| for ([_]Uplo{ .l, .u }) |uplo| {
                try decompose(T, al, check.io, r, .compute, uplo, a0, w, v);
                try expectSym(T, al, a0, w, v, @tagName(r));
                try decompose(T, al, check.io, r, .none, uplo, a0, w2, null);
                for (w, w2) |x, y| try testing.expect(@abs(x - y) <= tolerance(F, n, auxiliary.normMax(T, a0)));
            };
        }
    }
}

test "faer evd: special matrices, faer's test_special: zero, ones, identity; rank one, clustered, graded" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        for ([_]usize{ 1, 2, 4, 15, 50, 150 }) |n| for (0..6) |kind| {
            const a0 = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, a0);
            switch (kind) {
                0 => a0.fill(ops.zero(T)),
                1 => a0.fill(ops.one(T)),
                2 => a0.identity(),
                // clustered: I + eps-sized noise
                3 => {
                    lcheck.hermitian(T, a0, n);
                    for (0..n) |j| for (a0.col(j)[0..n]) |*x| {
                        x.* = ops.mulReal(T, x.*, std.math.floatEps(F));
                    };
                    for (0..n) |i| a0.at(i, i).* = ops.add(T, a0.get(i, i), ops.one(T));
                },
                // graded: D H D, D from 1 down by 10 per two rows
                4 => {
                    lcheck.hermitian(T, a0, n + 1);
                    for (0..n) |j| for (0..n) |i| {
                        const s = std.math.pow(F, 10, -@as(F, @floatFromInt(((i + j) / 4) % 8)));
                        a0.at(i, j).* = ops.mulReal(T, a0.get(i, j), s);
                    };
                },
                // tridiagonal Wilkinson W+ (a close pair at the top)
                else => {
                    a0.fill(ops.zero(T));
                    for (0..n) |i| {
                        const c = @as(F, @floatFromInt(n / 2)) - @as(F, @floatFromInt(i));
                        a0.at(i, i).* = ops.scalar(T, @abs(c));
                        if (i + 1 < n) {
                            a0.at(i + 1, i).* = ops.one(T);
                            a0.at(i, i + 1).* = ops.one(T);
                        }
                    }
                },
            }
            const w = try al.alloc(F, n);
            defer al.free(w);
            const v = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, v);
            for ([_]Routine{ .syevd, .syev }) |r| {
                decompose(T, al, check.io, r, .compute, .l, a0, w, v) catch |e| {
                    check.say("{s} {s} n={d} kind {d}: {t}\n", .{ @tagName(r), @typeName(T), n, kind, e });
                    return e;
                };
                try expectSym(T, al, a0, w, v, @tagName(r));
            }
        };
    }
}

test "faer evd: OpenBLAS's test_laed3_nan.c: dsyevd on I + 10 v vᵀ, nearly degenerate, no NaN" {
    // OpenBLAS `test/test_laed3_nan.c` (71f22601): n = 8, v = (1, 1 + 1e-8,
    // 1, ...), jobz 'V', uplo 'U'; here through `syevd` (n below 128: QR
    // iteration), through divide and conquer at faer's test threshold 8 and
    // at 4, and at n = 200, where `syevd` divides
    const al = testing.allocator;
    for ([_]usize{ 8, 200 }) |n| {
        const a0 = try ops.alloc(f64, al, n, n);
        defer ops.free(f64, al, a0);
        const v = try al.alloc(f64, n);
        defer al.free(v);
        @memset(v, 1);
        v[1] = 1 + 1e-8;
        for (0..n) |j| for (0..n) |i| {
            a0.at(i, j).* = (if (i == j) @as(f64, 1) else 0) + 10 * v[i] * v[j];
        };
        const w = try al.alloc(f64, n);
        defer al.free(w);
        const a = try ops.alloc(f64, al, n, n);
        defer ops.free(f64, al, a);
        a.copyFrom(a0);
        try evd.syevd(al, check.io, f64, .compute, .u, n, a.slice(), n, w);
        for (w) |x| try testing.expect(!std.math.isNan(x));
        for (a.ptr[0 .. n * n]) |x| try testing.expect(!std.math.isNan(x));
        try expectSym(f64, al, a0, w, a, "laed3");
        for ([_]usize{ 4, 8 }) |thr| {
            a.copyFrom(a0);
            try evd.symmetric(al, check.io, f64, true, false, a, w, thr);
            for (a.ptr[0 .. n * n]) |x| try testing.expect(!std.math.isNan(x));
            try expectSym(f64, al, a0, w, a, "laed3 dc");
        }
    }
}

test "faer evd: NaN and infinity are NoConvergence; scaled past the bounds; 1 x 1 and empty" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex32 }) |T| {
        const F = Float(T);
        for ([_]F{ std.math.nan(F), std.math.inf(F) }) |bad| {
            var a: [9]T = @splat(ops.one(T));
            a[1] = ops.scalar(T, bad);
            var w: [3]F = undefined;
            try testing.expectError(error.NoConvergence, evd.syevd(al, check.io, T, .compute, .l, 3, &a, 3, &w));
            // the other triangle is not read
            a[1] = ops.one(T);
            a[3] = ops.scalar(T, bad);
            try evd.syevd(al, check.io, T, .none, .l, 3, &a, 3, &w);
        }
        const n = 20;
        for ([_]F{ std.math.floatMin(F) * 4, std.math.floatMax(F) / 4 }) |sc| {
            const a0 = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, a0);
            lcheck.hermitian(T, a0, 3);
            var w0: [n]F = undefined;
            var w1: [n]F = undefined;
            const a = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, a);
            a.copyFrom(a0);
            try evd.syevd(al, check.io, T, .none, .l, n, a.slice(), n, &w0);
            for (0..n) |j| for (a0.col(j)[0..n]) |*x| {
                x.* = ops.mulReal(T, x.*, sc / 8);
            };
            try evd.syevd(al, check.io, T, .none, .l, n, a0.slice(), n, &w1);
            const tol = tolerance(F, n, 1) * (sc / 8) * 8;
            for (w0, w1) |x, y| try testing.expect(@abs(x * (sc / 8) - y) <= tol);
        }
        var one = [_]T{ops.scalar(T, -3)};
        var w: [1]F = undefined;
        try evd.syevd(al, check.io, T, .compute, .u, 1, &one, 1, &w);
        try testing.expectEqual(@as(F, -3), w[0]);
        try testing.expectEqual(@as(F, 1), ops.re(T, one[0]));
    }
    const fa = testing.failing_allocator;
    var a: [4]f64 = @splat(std.math.nan(f64));
    var w: [2]f64 = @splat(std.math.nan(f64));
    try testing.expectError(error.UnsupportedUplo, evd.syevd(fa, check.io, f64, .compute, .v, 2, &a, 2, &w));
    try testing.expectError(error.LeadingDimensionTooSmall, evd.syev(fa, check.io, f64, .compute, .u, 2, &a, 1, &w));
    try testing.expectError(error.SliceTooShort, evd.syevd(fa, check.io, f64, .compute, .l, 2, a[0..3], 2, &w));
    try testing.expectError(error.SliceTooShort, evd.syev(fa, check.io, f64, .none, .l, 2, &a, 2, w[0..1]));
    try evd.syevd(fa, check.io, f64, .compute, .l, 0, &a, 1, &w);
    for (a ++ w) |x| try testing.expect(std.math.isNan(x));
}

test "faer evd: syevd gives one worker's bits at any worker count" {
    if (comptime !test_io.hosted) return error.SkipZigTest;
    _ = blas.coreName();
    const al = testing.allocator;
    var runs: test_io.WorkerRuns = .init(al);
    defer runs.deinit();
    inline for (.{ f64, types.Complex32 }) |T| {
        const n = 260;
        const a0 = try ops.alloc(T, al, n, n);
        defer ops.free(T, al, a0);
        lcheck.hermitian(T, a0, 5);
        var outs: [3][]T = undefined;
        var ws: [3][]Float(T) = undefined;
        for (0..3) |c| {
            outs[c] = try al.alloc(T, n * n);
            ws[c] = try al.alloc(Float(T), n);
        }
        defer for (0..3) |c| {
            al.free(outs[c]);
            al.free(ws[c]);
        };
        for (0..test_io.worker_runs.len) |c| {
            const io = runs.io(c);
            @memcpy(outs[c], a0.ptr[0 .. n * n]);
            try evd.syevd(al, io, T, .compute, .l, n, outs[c], n, ws[c]);
        }
        try test_io.expectSameBytes(outs);
        try test_io.expectSameBytes(ws);
    }
}

test "syevd, syev: graded large at the bottom right given as .u, and at the top left as .l, keep the small eigenvalues (T125)" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64 }) |T| {
        const F = Float(T);
        const W = if (comptime isComplex(T)) std.math.Complex(f128) else f128;
        for ([_]usize{ 4, 6, 10 }) |n| for ([_]bool{ true, false }) |increasing| for ([_]bool{ false, true }) |use_qr| {
            // M1's E2: D B D, B Hermitian with a dominant diagonal, D geometric
            const ratio: f64 = if (n == 4) 1e-4 else 1e-2;
            var prng = std.Random.DefaultPrng.init(n * 7 + @intFromBool(increasing));
            const r = prng.random();
            const a = try al.alloc(T, n * n);
            defer al.free(a);
            const a128 = try al.alloc(W, n * n);
            defer al.free(a128);
            for (0..n) |j| for (0..j + 1) |i| {
                const di = std.math.pow(f64, ratio, @floatFromInt(if (increasing) n - 1 - i else i));
                const dj = std.math.pow(f64, ratio, @floatFromInt(if (increasing) n - 1 - j else j));
                const br = if (i == j) 2 * @as(f64, @floatFromInt(n)) + r.float(f64) else r.float(f64) - 0.5;
                const bi = if (i == j or !isComplex(T)) 0 else r.float(f64) - 0.5;
                const re = di * br * dj;
                const im = di * bi * dj;
                if (comptime isComplex(T)) {
                    a[i + j * n] = .init(re, im);
                    a[j + i * n] = .init(re, -im);
                    a128[i + j * n] = .init(re, im);
                    a128[j + i * n] = .init(re, -im);
                } else {
                    a[i + j * n] = re;
                    a[j + i * n] = re;
                    a128[i + j * n] = re;
                    a128[j + i * n] = re;
                }
            };
            const uplo: Uplo = if (increasing) .u else .l;
            const w = try al.alloc(F, n);
            defer al.free(w);
            const w128 = try al.alloc(f128, n);
            defer al.free(w128);
            if (use_qr) try evd.syev(al, check.io, T, .none, uplo, n, a, n, w) else try evd.syevd(al, check.io, T, .none, uplo, n, a, n, w);
            try evd.syevd(al, check.io, W, .none, uplo, n, a128, n, w128);
            // graded small eigenvalues hold to a few hundred eps (matlab88's
            // bottom-up reduction reaches 1.7e-13 here, M1's E2); reduced
            // from the wrong end they were off by 0.3 to 3e8
            for (w, w128, 0..) |x, ref, i| {
                const err: f128 = @abs(@as(f128, x) - ref) / @abs(ref);
                if (!(err <= 1024 * @as(f128, @floatFromInt(n)) * std.math.floatEps(F))) {
                    check.say("{s} {s} n={d} increasing={}: eigenvalue {d} relative error {e}\n", .{ if (use_qr) "syev" else "syevd", @typeName(T), n, increasing, i, @as(f64, @floatCast(err)) });
                    return error.TestUnexpectedResult;
                }
            }
        };
    }
}
