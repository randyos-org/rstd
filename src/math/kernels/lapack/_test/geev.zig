//! Tests of `geev` (`../evd.zig`).

const std = @import("std");
const types = @import("../../blas/types.zig");
const blas = @import("../../blas/__root__.zig");
const apple_accel = @import("../../blas/apple_accel.zig");
const test_io = @import("../../blas/_test/io.zig");
const check = @import("../../blas/_test/check.zig");
const lcheck = @import("check.zig");
const ops = @import("../ops.zig");
const auxiliary = @import("../auxiliary.zig");
const evd = @import("../evd.zig");
const schur = @import("../schur.zig");
const evd_tests = @import("evd.zig");

const testing = std.testing;
const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = ops.Mat;
const Float = ops.Float;
const isComplex = ops.isComplex;
const Jobz = evd.Jobz;
const Eigenvalues = schur.Eigenvalues;
const test_types = evd_tests.test_types;
const tolerance = evd_tests.tolerance;

/// Column `j` of `geev`'s vectors as complex, a real pair unpacked.
fn vectorOf(comptime T: type, comptime C: type, wi: ?[]const Float(T), v: Mat(T), j: usize, out: []C) void {
    for (out, 0..) |*o, i| {
        if (comptime isComplex(T)) {
            o.* = v.get(i, j);
        } else if (wi.?[j] == 0) {
            o.* = .init(v.get(i, j), 0);
        } else if (wi.?[j] > 0) {
            o.* = .init(v.get(i, j), v.get(i, j + 1));
        } else {
            o.* = .init(v.get(i, j - 1), -v.get(i, j));
        }
    }
}

/// `‖A v - λ v‖` and `‖uᴴ A - λ uᴴ‖` per vector within `tolerance`, each
/// of unit norm with its largest entry real.
fn expectGeev(comptime T: type, al: Allocator, a0: Mat(T), wr: []const Float(T), wi: ?[]const Float(T), w: ?[]const T, vl: Mat(T), vr: Mat(T)) !void {
    const F = Float(T);
    const C = std.math.Complex(F);
    const n = a0.rows;
    const tol = tolerance(F, n, auxiliary.normMax(T, a0));
    const x = try al.alloc(C, n);
    defer al.free(x);
    for (0..n) |j| {
        const lam: C = if (comptime isComplex(T)) w.?[j] else .init(wr[j], wi.?[j]);
        for ([_]Mat(T){ vr, vl }, 0..) |v, side| {
            vectorOf(T, C, wi, v, j, x);
            // unit norm, largest entry real
            var nrm: F = 0;
            var big: F = 0;
            for (x) |e| {
                const e2 = e.re * e.re + e.im * e.im;
                nrm += e2;
                big = @max(big, e2);
            }
            try testing.expect(@abs(@sqrt(nrm) - 1) <= tolerance(F, n, 1));
            // entries tied for largest (a rotation's) may round either way first
            var real_big = false;
            for (x) |e| {
                const e2 = e.re * e.re + e.im * e.im;
                if (e2 >= big * (1 - tolerance(F, n, 1)) and @abs(e.im) <= tolerance(F, n, 1)) real_big = true;
            }
            try testing.expect(real_big);
            // right: (A x)(i) - λ x(i); left: (xᴴ A)(i) - λ conj(x(i))
            var worst: F = 0;
            for (0..n) |i| {
                var s: C = .init(0, 0);
                for (0..n) |k| {
                    const e = if (side == 0) a0.get(i, k) else a0.get(k, i);
                    const aik: C = if (comptime isComplex(T)) e else .init(e, 0);
                    const xk = if (side == 0) x[k] else x[k].conjugate();
                    s = s.add(aik.mul(xk));
                }
                const xi = if (side == 0) x[i] else x[i].conjugate();
                const r = s.sub(lam.mul(xi));
                worst = @max(worst, @sqrt(r.re * r.re + r.im * r.im));
            }
            if (!(worst <= tol)) {
                check.say("geev {s} n={d}: {s} vector {d} residual {e}, tol {e}\n", .{ @typeName(T), n, if (side == 0) "right" else "left", j, worst, tol });
                return error.TestUnexpectedResult;
            }
        }
    }
}

/// `geev` on a copy of `a0` with both vector sets, then values alone, which
/// must agree.
fn geevCase(comptime T: type, al: Allocator, io: Io, a0: Mat(T)) !void {
    const F = Float(T);
    const n = a0.rows;
    const a = try ops.alloc(T, al, n + 1, n);
    defer ops.free(T, al, a);
    const vl = try ops.alloc(T, al, n, n);
    defer ops.free(T, al, vl);
    const vr = try ops.alloc(T, al, n, n);
    defer ops.free(T, al, vr);
    const ev = try al.alloc(F, 4 * n);
    defer al.free(ev);
    const cw = try al.alloc(T, 2 * n);
    defer al.free(cw);
    for (0..2) |pass| {
        a.fill(lcheck.poisoned(T));
        a.sub(0, 0, n, n).copyFrom(a0);
        const jv: Jobz = if (pass == 0) .v else .n;
        const off = pass * 2 * n;
        const w: Eigenvalues(T) = if (comptime isComplex(T)) cw[pass * n ..][0..n] else .{ .re = ev[off..][0..n], .im = ev[off + n ..][0..n] };
        try evd.geev(al, io, T, jv, jv, n, a.slice(), n + 1, w, vl.slice(), n, vr.slice(), n, null);
        for (0..n) |j| try testing.expect(lcheck.isNanT(T, a.get(n, j)));
    }
    const etol = tolerance(F, n, auxiliary.normMax(T, a0));
    if (comptime isComplex(T)) {
        try expectGeev(T, al, a0, &.{}, null, cw[0..n], vl, vr);
        for (cw[0..n], cw[n..]) |x, y| try testing.expect(ops.abs(T, ops.sub(T, x, y)) <= etol);
    } else {
        try expectGeev(T, al, a0, ev[0..n], ev[n .. 2 * n], null, vl, vr);
        for (ev[0 .. 2 * n], ev[2 * n ..]) |x, y| try testing.expect(@abs(x - y) <= etol);
        // pairs: positive imaginary part first
        var j: usize = 0;
        while (j < n) : (j += 1) if (ev[n + j] != 0) {
            try testing.expect(ev[n + j] > 0 and ev[n + j + 1] == -ev[n + j] and ev[j] == ev[j + 1]);
            j += 1;
        };
    }
}

test "faer evd: geev, A v = λ v and uᴴ A = λ uᴴ, unit vectors, faer's sizes" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        for ([_]usize{ 1, 2, 3, 4, 10, 15, 20, 50, 100 }) |n| {
            const a0 = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, a0);
            check.fill(T, a0.ptr[0 .. n * n], 3 * n + 2);
            geevCase(T, al, check.io, a0) catch |e| {
                check.say("geev {s} n={d}\n", .{ @typeName(T), n });
                return e;
            };
        }
    }
}

test "faer evd: geev in the extended types, real and complex (T124)" {
    const al = testing.allocator;
    inline for (check.x_real ++ check.x_complex ++ check.f80_unslotted ++ .{ f128, std.math.Complex(f128) }) |T| {
        for ([_]usize{ 1, 2, 3, 10 }) |n| {
            const a0 = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, a0);
            check.fill(T, a0.ptr[0 .. n * n], 3 * n + 2);
            geevCase(T, al, check.io, a0) catch |e| {
                check.say("geev {s} n={d}\n", .{ @typeName(T), n });
                return e;
            };
        }
    }
}

test "faer evd: geev on special matrices: zero, identity, Jordan, rotation pairs, companion, triangular" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        for ([_]usize{ 1, 2, 5, 12 }) |n| for (0..6) |kind| {
            // a companion's roots are ill-conditioned past n = 5 or so; a
            // Jordan block's vectors are only as good as eps^(1/n)
            if ((kind == 4 and n > 5) or (kind == 2 and n > 2)) continue;
            const a0 = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, a0);
            a0.fill(ops.zero(T));
            switch (kind) {
                0 => {},
                1 => a0.identity(),
                2 => for (0..n) |i| {
                    a0.at(i, i).* = ops.scalar(T, 2);
                    if (i + 1 < n) a0.at(i, i + 1).* = ops.one(T);
                },
                // rotations: eigenvalues ±i, and 1 for odd n
                3 => {
                    var i: usize = 0;
                    while (i + 1 < n) : (i += 2) {
                        a0.at(i + 1, i).* = ops.one(T);
                        a0.at(i, i + 1).* = ops.scalar(T, -1);
                    }
                    if (n % 2 == 1) a0.at(n - 1, n - 1).* = ops.one(T);
                },
                // the companion matrix of (x - 1)(x - 2)...(x - n)
                4 => {
                    var c: [13]F = @splat(0);
                    c[0] = 1;
                    for (1..n + 1) |r| {
                        var k = r;
                        while (k > 0) : (k -= 1) c[k] = c[k] - @as(F, @floatFromInt(r)) * c[k - 1];
                    }
                    for (0..n) |j| a0.at(0, j).* = ops.scalar(T, -c[j + 1]);
                    for (1..n) |i| a0.at(i, i - 1).* = ops.one(T);
                },
                else => {
                    check.fill(T, a0.ptr[0 .. n * n], n);
                    for (0..n) |j| for (j + 1..n) |i| {
                        a0.at(i, j).* = ops.zero(T);
                    };
                },
            }
            geevCase(T, al, check.io, a0) catch |e| {
                check.say("geev {s} n={d} kind {d}\n", .{ @typeName(T), n, kind });
                return e;
            };
        };
    }
}

test "faer evd: geev's eigenvectors that overflow unscaled (xTREVC3's scaling)" {
    // upper triangular, T(i, i) = i h, ones above: back substitution for
    // the last (first, for left vectors) eigenvalue grows like
    // 1 / (k! h^k), past the largest float by about 2^500 for f64 (2^60
    // for f32), so only a scaled solve stays finite
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        const n = 20;
        const e = (3 * std.math.floatExponentMax(F)) / (2 * (n - 1)) + 1;
        const h = std.math.ldexp(@as(F, 1), -e);
        const a0 = try ops.alloc(T, al, n, n);
        defer ops.free(T, al, a0);
        a0.fill(ops.zero(T));
        for (0..n) |i| {
            a0.at(i, i).* = ops.scalar(T, @as(F, @floatFromInt(i)) * h);
            if (i + 1 < n) a0.at(i, i + 1).* = ops.one(T);
        }
        geevCase(T, al, check.io, a0) catch |err| {
            check.say("geev overflow {s}\n", .{@typeName(T)});
            return err;
        };
    }
}

test "faer evd: geev's NaN is NanInA, infinity NoConvergence; scaled past the bounds; arguments" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex32 }) |T| {
        const F = Float(T);
        var a: [9]T = undefined;
        check.fill(T, &a, 1);
        var e: [6]F = undefined;
        var ce: [3]T = undefined;
        const w: Eigenvalues(T) = if (comptime isComplex(T)) &ce else .{ .re = e[0..3], .im = e[3..] };
        var v: [9]T = undefined;
        a[4] = ops.scalar(T, std.math.nan(F));
        try testing.expectError(error.NanInA, evd.geev(al, check.io, T, .v, .n, 3, &a, 3, w, &v, 3, v[0..0], 1, null));
        a[4] = ops.scalar(T, std.math.inf(F));
        try testing.expectError(error.NoConvergence, evd.geev(al, check.io, T, .n, .n, 3, &a, 3, w, v[0..0], 1, v[0..0], 1, null));
        // scaled past SMLNUM and BIGNUM: the unscaled matrix's vectors
        // (residuals on it, as they would overflow on the scaled one) and
        // its values times the scale
        // powers of 2, so the scaled matrix is exact: normal, past SMLNUM;
        // and below the overflow threshold, past BIGNUM
        const exps = .{ std.math.floatExponentMin(F) + 64, std.math.floatExponentMax(F) - 8 };
        inline for (exps) |ex| {
            const sc = 8 * std.math.ldexp(@as(F, 1), ex);
            const n = 12;
            const a0 = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, a0);
            check.fill(T, a0.ptr[0 .. n * n], 9);
            const as = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, as);
            for (0..n) |j| for (0..n) |i| {
                as.at(i, j).* = ops.mulReal(T, a0.get(i, j), sc / 8);
            };
            const vl = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, vl);
            const vr = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, vr);
            var ev: [2 * n]F = undefined;
            var cw: [n]T = undefined;
            const ws: Eigenvalues(T) = if (comptime isComplex(T)) &cw else .{ .re = ev[0..n], .im = ev[n..] };
            try evd.geev(al, check.io, T, .v, .v, n, as.slice(), n, ws, vl.slice(), n, vr.slice(), n, null);
            for (&ev) |*x| x.* /= sc / 8;
            for (&cw) |*x| x.* = ops.mulReal(T, x.*, 1 / (sc / 8));
            if (comptime isComplex(T)) {
                try expectGeev(T, al, a0, &.{}, null, &cw, vl, vr);
            } else {
                try expectGeev(T, al, a0, ev[0..n], ev[n..], null, vl, vr);
            }
        }
    }
    const fa = testing.failing_allocator;
    var a: [4]f64 = @splat(std.math.nan(f64));
    var wr: [2]f64 = @splat(std.math.nan(f64));
    var wi: [2]f64 = @splat(std.math.nan(f64));
    const w: Eigenvalues(f64) = .{ .re = &wr, .im = &wi };
    try testing.expectError(error.LeadingDimensionTooSmall, evd.geev(fa, check.io, f64, .v, .v, 2, &a, 1, w, &a, 2, &a, 2, null));
    try testing.expectError(error.LeadingDimensionTooSmall, evd.geev(fa, check.io, f64, .v, .v, 2, &a, 2, w, &a, 1, &a, 2, null));
    try testing.expectError(error.LeadingDimensionTooSmall, evd.geev(fa, check.io, f64, .n, .v, 2, &a, 2, w, a[0..0], 0, &a, 1, null));
    try testing.expectError(error.SliceTooShort, evd.geev(fa, check.io, f64, .n, .n, 2, a[0..3], 2, w, a[0..0], 0, a[0..0], 0, null));
    try testing.expectError(error.SliceTooShort, evd.geev(fa, check.io, f64, .n, .v, 2, &a, 2, w, a[0..0], 0, a[0..3], 2, null));
    try testing.expectError(error.SliceTooShort, evd.geev(fa, check.io, f64, .n, .n, 2, &a, 2, .{ .re = &wr, .im = wi[0..1] }, a[0..0], 0, a[0..0], 0, null));
    try evd.geev(fa, check.io, f64, .v, .v, 0, &a, 1, w, &a, 1, &a, 1, null);
    for (a ++ wr ++ wi) |x| try testing.expect(std.math.isNan(x));
}

test "faer evd: geev's INFO > 0: w[i..n] converged and scaled back, the count in unconverged" {
    const al = testing.allocator;
    defer ops.test_iter_limit = null;
    inline for (test_types) |T| {
        const F = Float(T);
        const n = 12;
        var a0: [n * n]T = undefined;
        check.fill(T, &a0, 31);
        var v: [1]T = undefined;
        // the true eigenvalues, then capped runs: unscaled and past SMLNUM
        var ev: [2 * n]F = undefined;
        var cw: [n]T = undefined;
        var a = a0;
        const w0: Eigenvalues(T) = if (comptime isComplex(T)) &cw else .{ .re = ev[0..n], .im = ev[n..] };
        try evd.geev(al, check.io, T, .n, .n, n, &a, n, w0, &v, 1, &v, 1, null);
        for ([_]F{ 1, std.math.ldexp(@as(F, 1), std.math.floatExponentMin(F) + 8) }) |sc| {
            ops.test_iter_limit = 10;
            defer ops.test_iter_limit = null;
            for (&a, a0) |*x, y| x.* = ops.mulReal(T, y, sc);
            var ev1: [2 * n]F = undefined;
            var cw1: [n]T = undefined;
            const w1: Eigenvalues(T) = if (comptime isComplex(T)) &cw1 else .{ .re = ev1[0..n], .im = ev1[n..] };
            var info: usize = 0;
            var vr: [n * n]T = undefined;
            try testing.expectError(error.NoConvergence, evd.geev(al, check.io, T, .n, .v, n, &a, n, w1, &v, 1, &vr, n, &info));
            try testing.expect(info > 0 and info < n);
            const tol = tolerance(F, n, auxiliary.normMax(T, .{ .ptr = &a0, .rows = n, .cols = n, .ld = n }));
            for (info..n) |i| {
                var best = std.math.inf(F);
                for (0..n) |j| {
                    const d = if (comptime isComplex(T)) ops.abs(T, ops.sub(T, ops.mulReal(T, cw1[i], 1 / sc), cw[j])) else std.math.hypot(ev1[i] / sc - ev[j], ev1[n + i] / sc - ev[n + j]);
                    best = @min(best, d);
                }
                try testing.expect(best <= tol);
            }
        }
        // an infinity: none converged
        a = a0;
        a[5] = ops.scalar(T, std.math.inf(F));
        var info: usize = 0;
        try testing.expectError(error.NoConvergence, evd.geev(al, check.io, T, .n, .n, n, &a, n, w0, &v, 1, &v, 1, &info));
        try testing.expectEqual(n, info);
    }
}

test "faer evd: geev gives one worker's bits at any worker count" {
    if (comptime !test_io.hosted) return error.SkipZigTest;
    _ = blas.coreName();
    // rstd's kernels: one worker's run is one chunk and the others many, and
    // Accelerate's bits follow a part's shape (its own split keeps the chunks)
    const was = apple_accel.routing();
    _ = apple_accel.route(false);
    defer _ = apple_accel.route(was);
    const al = testing.allocator;
    var runs: test_io.WorkerRuns = .init(al);
    defer runs.deinit();
    inline for (.{ f64, types.Complex32 }) |T| {
        const n = 150;
        const a0 = try ops.alloc(T, al, n, n);
        defer ops.free(T, al, a0);
        check.fill(T, a0.ptr[0 .. n * n], 11);
        var outs: [3][]T = undefined;
        var evs: [3][]T = undefined;
        for (0..3) |c| {
            outs[c] = try al.alloc(T, 3 * n * n);
            evs[c] = try al.alloc(T, 2 * n);
        }
        defer for (0..3) |c| {
            al.free(outs[c]);
            al.free(evs[c]);
        };
        for (0..test_io.worker_runs.len) |c| {
            const io = runs.io(c);
            const o = outs[c];
            @memcpy(o[0 .. n * n], a0.ptr[0 .. n * n]);
            const w: Eigenvalues(T) = if (comptime isComplex(T)) evs[c][0..n] else .{ .re = evs[c][0..n], .im = evs[c][n..] };
            try evd.geev(al, io, T, .v, .v, n, o[0 .. n * n], n, w, o[n * n ..][0 .. n * n], n, o[2 * n * n ..], n, null);
        }
        try test_io.expectSameBytes(outs);
        try test_io.expectSameBytes(evs);
    }
}

test "faer evd: geev balances, so D M D⁻¹ keeps M's eigenvalues to M's scale (dgebal.f)" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        const n = 16;
        const m = try ops.alloc(T, al, n, n);
        defer ops.free(T, al, m);
        check.fill(T, m.ptr[0 .. n * n], 21);
        const a = try ops.alloc(T, al, n, n);
        defer ops.free(T, al, a);
        // D = diag(2^(3i)): entries from 2^-45 to 2^45 times M's
        for (0..n) |j| for (0..n) |i| {
            const d = std.math.pow(F, 2, 3 * (@as(F, @floatFromInt(i)) - @as(F, @floatFromInt(j))));
            a.at(i, j).* = ops.mulReal(T, m.get(i, j), d);
        };
        var ev: [4 * n]F = undefined;
        var cw: [2 * n]T = undefined;
        var v: [1]T = undefined;
        for ([_]Mat(T){ m, a }, 0..) |x, p| {
            const w: Eigenvalues(T) = if (comptime isComplex(T)) cw[p * n ..][0..n] else .{ .re = ev[2 * p * n ..][0..n], .im = ev[2 * p * n + n ..][0..n] };
            try evd.geev(al, check.io, T, .n, .n, n, x.slice(), n, w, &v, 1, &v, 1, null);
        }
        // each of A's within tolerance of one of M's
        const tol = tolerance(F, n, 4);
        for (0..n) |i| {
            var best = std.math.inf(F);
            for (0..n) |j| {
                const d = if (comptime isComplex(T)) ops.abs(T, ops.sub(T, cw[n + i], cw[j])) else std.math.hypot(ev[2 * n + i] - ev[j], ev[3 * n + i] - ev[n + j]);
                best = @min(best, d);
            }
            try testing.expect(best <= tol);
        }
    }
}
