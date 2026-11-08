//! Tests of `evd.zig`: `syevd`, `syev`, `geev`, `hseqr` and `gees`.

const std = @import("std");
const types = @import("../../blas/types.zig");
const blas = @import("../../blas/__root__.zig");
const apple_accel = @import("../../blas/apple_accel.zig");
const test_io = @import("../../blas/_test/io.zig");
const check = @import("../../blas/_test/check.zig");
const lcheck = @import("check.zig");
const ops = @import("../ops.zig");
const blas_calls = @import("../blas_calls.zig");
const auxiliary = @import("../auxiliary.zig");
const hessenberg = @import("../hessenberg.zig");
const evd = @import("../evd.zig");

const testing = std.testing;
const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = ops.Mat;
const Float = ops.Float;
const isComplex = ops.isComplex;
const Uplo = types.Uplo;
const Jobz = evd.Jobz;
const Eigenvalues = evd.Eigenvalues;
const GeesSelect = evd.GeesSelect;

const test_types = .{ f32, f64, types.Complex32, types.Complex64 };

/// Faer's `test_self_adjoint_evd` tolerance: `ApproxEq::eps()` (`128 eps`)
/// times `sqrt(8 n)`, against `max(1, ‖A‖max)` per element.
fn tolerance(comptime F: type, n: usize, amax: F) F {
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
fn decompose(comptime T: type, al: Allocator, io: Io, r: Routine, jobz: Jobz, uplo: Uplo, a0: Mat(T), w: []Float(T), v: ?Mat(T)) !void {
    const n = a0.rows;
    const lda = n + 1;
    const a = try ops.alloc(T, al, lda, n);
    defer ops.free(T, al, a);
    a.fill(lcheck.poisoned(T));
    for (0..n) |j| for (0..n) |i| {
        if (if (uplo == .l) i >= j else i <= j) a.at(i, j).* = a0.get(i, j);
    };
    switch (r) {
        .syevd => try evd.syevd(al, io, T, jobz, uplo, @intCast(n), a.ptr, @intCast(lda), w.ptr),
        .syev => try evd.syev(al, io, T, jobz, uplo, @intCast(n), a.ptr, @intCast(lda), w.ptr),
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
                try decompose(T, al, check.io, r, .v, uplo, a0, w, v);
                try expectSym(T, al, a0, w, v, @tagName(r));
                try decompose(T, al, check.io, r, .n, uplo, a0, w2, null);
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
                decompose(T, al, check.io, r, .v, .l, a0, w, v) catch |e| {
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
        try evd.syevd(al, check.io, f64, .v, .u, @intCast(n), a.ptr, @intCast(n), w.ptr);
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
            try testing.expectError(error.NoConvergence, evd.syevd(al, check.io, T, .v, .l, 3, &a, 3, &w));
            // the other triangle is not read
            a[1] = ops.one(T);
            a[3] = ops.scalar(T, bad);
            try evd.syevd(al, check.io, T, .n, .l, 3, &a, 3, &w);
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
            try evd.syevd(al, check.io, T, .n, .l, n, a.ptr, n, &w0);
            for (0..n) |j| for (a0.col(j)[0..n]) |*x| {
                x.* = ops.mulReal(T, x.*, sc / 8);
            };
            try evd.syevd(al, check.io, T, .n, .l, n, a0.ptr, n, &w1);
            const tol = tolerance(F, n, 1) * (sc / 8) * 8;
            for (w0, w1) |x, y| try testing.expect(@abs(x * (sc / 8) - y) <= tol);
        }
        var one = [_]T{ops.scalar(T, -3)};
        var w: [1]F = undefined;
        try evd.syevd(al, check.io, T, .v, .u, 1, &one, 1, &w);
        try testing.expectEqual(@as(F, -3), w[0]);
        try testing.expectEqual(@as(F, 1), ops.re(T, one[0]));
    }
    const fa = testing.failing_allocator;
    var a: [4]f64 = @splat(std.math.nan(f64));
    var w: [2]f64 = @splat(std.math.nan(f64));
    try testing.expectError(error.InvalidUplo, evd.syevd(fa, check.io, f64, .v, .v, -1, &a, 0, &w));
    try testing.expectError(error.InvalidN, evd.syevd(fa, check.io, f64, .v, .l, -1, &a, 0, &w));
    try testing.expectError(error.InvalidLda, evd.syev(fa, check.io, f64, .v, .u, 2, &a, 1, &w));
    try evd.syevd(fa, check.io, f64, .v, .l, 0, &a, 1, &w);
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
            try evd.syevd(al, io, T, .v, .l, n, outs[c].ptr, n, ws[c].ptr);
        }
        try test_io.expectSameBytes(outs);
        try test_io.expectSameBytes(ws);
    }
}

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
            var big: F = -1;
            var kbig: usize = 0;
            for (x, 0..) |e, k| {
                const e2 = e.re * e.re + e.im * e.im;
                nrm += e2;
                if (e2 > big) {
                    big = e2;
                    kbig = k;
                }
            }
            try testing.expect(@abs(@sqrt(nrm) - 1) <= tolerance(F, n, 1));
            try testing.expect(@abs(x[kbig].im) <= tolerance(F, n, 1));
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
        const w: Eigenvalues(T) = if (comptime isComplex(T)) cw.ptr + pass * n else .{ .re = ev.ptr + off, .im = ev.ptr + off + n };
        try evd.geev(al, io, T, jv, jv, @intCast(n), a.ptr, @intCast(n + 1), w, vl.ptr, @intCast(n), vr.ptr, @intCast(n), null);
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

test "faer evd: geev's NaN is InvalidA, infinity NoConvergence; scaled past the bounds; arguments" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex32 }) |T| {
        const F = Float(T);
        var a: [9]T = undefined;
        check.fill(T, &a, 1);
        var e: [6]F = undefined;
        var ce: [3]T = undefined;
        const w: Eigenvalues(T) = if (comptime isComplex(T)) &ce else .{ .re = &e, .im = e[3..].ptr };
        var v: [9]T = undefined;
        a[4] = ops.scalar(T, std.math.nan(F));
        try testing.expectError(error.InvalidA, evd.geev(al, check.io, T, .v, .n, 3, &a, 3, w, &v, 3, &v, 1, null));
        a[4] = ops.scalar(T, std.math.inf(F));
        try testing.expectError(error.NoConvergence, evd.geev(al, check.io, T, .n, .n, 3, &a, 3, w, &v, 1, &v, 1, null));
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
            const ws: Eigenvalues(T) = if (comptime isComplex(T)) &cw else .{ .re = &ev, .im = ev[n..].ptr };
            try evd.geev(al, check.io, T, .v, .v, n, as.ptr, n, ws, vl.ptr, n, vr.ptr, n, null);
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
    try testing.expectError(error.InvalidN, evd.geev(fa, check.io, f64, .v, .v, -1, &a, 0, w, &a, 0, &a, 0, null));
    try testing.expectError(error.InvalidLda, evd.geev(fa, check.io, f64, .v, .v, 2, &a, 1, w, &a, 0, &a, 0, null));
    try testing.expectError(error.InvalidLdvl, evd.geev(fa, check.io, f64, .v, .v, 2, &a, 2, w, &a, 1, &a, 0, null));
    try testing.expectError(error.InvalidLdvl, evd.geev(fa, check.io, f64, .n, .v, 2, &a, 2, w, &a, 0, &a, 2, null));
    try testing.expectError(error.InvalidLdvr, evd.geev(fa, check.io, f64, .n, .v, 2, &a, 2, w, &a, 1, &a, 1, null));
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
        const w0: Eigenvalues(T) = if (comptime isComplex(T)) &cw else .{ .re = &ev, .im = ev[n..].ptr };
        try evd.geev(al, check.io, T, .n, .n, n, &a, n, w0, &v, 1, &v, 1, null);
        for ([_]F{ 1, std.math.ldexp(@as(F, 1), std.math.floatExponentMin(F) + 8) }) |sc| {
            ops.test_iter_limit = 10;
            defer ops.test_iter_limit = null;
            for (&a, a0) |*x, y| x.* = ops.mulReal(T, y, sc);
            var ev1: [2 * n]F = undefined;
            var cw1: [n]T = undefined;
            const w1: Eigenvalues(T) = if (comptime isComplex(T)) &cw1 else .{ .re = &ev1, .im = ev1[n..].ptr };
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
            const w: Eigenvalues(T) = if (comptime isComplex(T)) evs[c].ptr else .{ .re = evs[c].ptr, .im = evs[c].ptr + n };
            try evd.geev(al, io, T, .v, .v, n, o.ptr, n, w, o.ptr + n * n, n, o.ptr + 2 * n * n, n, null);
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
            const w: Eigenvalues(T) = if (comptime isComplex(T)) cw[p * n ..].ptr else .{ .re = ev[2 * p * n ..].ptr, .im = ev[2 * p * n + n ..].ptr };
            try evd.geev(al, check.io, T, .n, .n, n, x.ptr, n, w, &v, 1, &v, 1, null);
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

test "hseqr: Zᴴ A Z = T after gehrd and orghr, the ilo..ihi block, both jobs, real and complex" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64, f32, types.Complex32 }) |T| {
        const F = Float(T);
        for ([_][3]usize{ .{ 1, 1, 1 }, .{ 2, 1, 2 }, .{ 6, 1, 6 }, .{ 9, 3, 7 }, .{ 30, 4, 27 } }) |c| {
            const n = c[0];
            const lo = c[1] - 1;
            const hi = c[2];
            const a0 = try al.alloc(T, n * n);
            defer al.free(a0);
            check.fill(T, a0, 11 * n + lo);
            for (0..n) |j| for (j + 1..n) |i| {
                if (j < lo or i >= hi) a0[i + j * n] = ops.zero(T);
            };
            const t = try al.dupe(T, a0);
            defer al.free(t);
            const tau = try al.alloc(T, @max(n, 2) - 1);
            defer al.free(tau);
            const ni: isize = @intCast(n);
            try hessenberg.gehrd(al, check.io, T, ni, @intCast(lo + 1), @intCast(hi), t.ptr, ni, tau.ptr);
            const z = try al.dupe(T, t);
            defer al.free(z);
            try hessenberg.orghr(al, check.io, T, ni, @intCast(lo + 1), @intCast(hi), z.ptr, ni, tau.ptr);
            const h2 = try al.dupe(T, t);
            defer al.free(h2);
            const ev = try al.alloc(F, 4 * n);
            defer al.free(ev);
            const cw = try al.alloc(T, 2 * n);
            defer al.free(cw);
            const w: Eigenvalues(T) = if (comptime isComplex(T)) cw.ptr else .{ .re = ev.ptr, .im = ev.ptr + n };
            const w2: Eigenvalues(T) = if (comptime isComplex(T)) cw.ptr + n else .{ .re = ev.ptr + 2 * n, .im = ev.ptr + 3 * n };
            try evd.hseqr(check.io, T, .s, .v, ni, @intCast(lo + 1), @intCast(hi), t.ptr, ni, w, z.ptr, ni, null);
            var dummy: [1]T = undefined;
            try evd.hseqr(check.io, T, .e, .n, ni, @intCast(lo + 1), @intCast(hi), h2.ptr, ni, w2, &dummy, 1, null);
            const nrm = @max(auxiliary.normMax(T, .{ .ptr = a0.ptr, .rows = n, .cols = n, .ld = n }), 1);
            const tol = 128 * std.math.floatEps(F) * @as(F, @floatFromInt(n + 1)) * nrm;
            for (0..n) |j| for (0..n) |i| {
                var s = ops.zero(T);
                for (0..n) |k| for (0..n) |l| {
                    s = ops.add(T, s, ops.mul(T, ops.mul(T, ops.conj(T, z[k + i * n]), a0[k + l * n]), z[l + j * n]));
                };
                // T upper triangular; real: a subdiagonal entry only inside a complex pair's block
                const pair = !isComplex(T) and i == j + 1 and (if (comptime isComplex(T)) false else ev[n + j] > 0);
                const want = if (i <= j or pair) t[i + j * n] else ops.zero(T);
                if (i > j and !pair and !ops.isZero(T, t[i + j * n])) return error.TestUnexpectedResult;
                var g = ops.zero(T);
                for (0..n) |k| g = ops.add(T, g, ops.mul(T, ops.conj(T, z[k + i * n]), z[k + j * n]));
                const err = ops.abs(T, ops.sub(T, s, want));
                const gerr = ops.abs(T, ops.sub(T, g, if (i == j) ops.one(T) else ops.zero(T)));
                if (!(err <= tol) or !(gerr <= tol / nrm)) {
                    check.say("hseqr {s} n={d} ilo={d} ihi={d}: ({d},{d}) T off by {e}, ZᴴZ by {e}\n", .{ @typeName(T), n, lo + 1, hi, i, j, err, gerr });
                    return error.TestUnexpectedResult;
                }
            };
            // the eigenvalues alone agree with the Schur form's
            if (comptime isComplex(T)) {
                for (cw[0..n], cw[n..]) |x, y| try testing.expect(ops.abs(T, ops.sub(T, x, y)) <= tol);
            } else for (0..n) |i| {
                try testing.expect(@abs(ev[i] - ev[2 * n + i]) <= tol and @abs(ev[n + i] - ev[3 * n + i]) <= tol);
            }
        }
    }
    var a: [4]f64 = @splat(0);
    var wr: [2]f64 = undefined;
    var wi: [2]f64 = undefined;
    const w: Eigenvalues(f64) = .{ .re = &wr, .im = &wi };
    try testing.expectError(error.InvalidN, evd.hseqr(check.io, f64, .e, .n, -1, 1, 0, &a, 1, w, &a, 1, null));
    try testing.expectError(error.InvalidIlo, evd.hseqr(check.io, f64, .e, .n, 2, 3, 2, &a, 2, w, &a, 1, null));
    try testing.expectError(error.InvalidIhi, evd.hseqr(check.io, f64, .e, .n, 2, 2, 1, &a, 2, w, &a, 1, null));
    try testing.expectError(error.InvalidLdh, evd.hseqr(check.io, f64, .e, .n, 2, 1, 2, &a, 1, w, &a, 1, null));
    try testing.expectError(error.InvalidLdz, evd.hseqr(check.io, f64, .s, .i, 2, 1, 2, &a, 2, w, &a, 1, null));
    try evd.hseqr(check.io, f64, .e, .n, 0, 1, 0, &a, 1, w, &a, 1, null);
}

/// `gees`'s test sort: the left half plane first.
fn LeftHalf(comptime T: type) type {
    return struct {
        fn real(re: Float(T), im: Float(T)) bool {
            _ = im;
            return re < 0;
        }
        fn cplx(x: T) bool {
            return x.re < 0;
        }
    };
}

test "gees: A = Z T Zᴴ, unsorted and the left half plane first, tiny scale, NaN; real and complex" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64, f32, types.Complex32 }) |T| {
        const F = Float(T);
        const pick: GeesSelect(T) = if (comptime isComplex(T)) &LeftHalf(T).cplx else &LeftHalf(T).real;
        for ([_]usize{ 1, 2, 6, 17 }) |n| for ([_]bool{ false, true }) |sort| for ([_]F{ 1, 1e-30 }) |sc| {
            const ni: isize = @intCast(n);
            const a0 = try al.alloc(T, n * n);
            defer al.free(a0);
            check.fill(T, a0, 17 * n + @intFromBool(sort));
            for (a0) |*e| e.* = ops.mulReal(T, e.*, sc);
            const t = try al.dupe(T, a0);
            defer al.free(t);
            const z = try al.alloc(T, n * n);
            defer al.free(z);
            const ev = try al.alloc(F, 2 * n);
            defer al.free(ev);
            const cw = try al.alloc(T, n);
            defer al.free(cw);
            const w: Eigenvalues(T) = if (comptime isComplex(T)) cw.ptr else .{ .re = ev.ptr, .im = ev.ptr + n };
            var sdim: isize = -1;
            try evd.gees(al, check.io, T, .v, if (sort) pick else null, ni, t.ptr, ni, &sdim, w, z.ptr, ni, null);
            const nrm = @max(auxiliary.normMax(T, .{ .ptr = a0.ptr, .rows = n, .cols = n, .ld = n }), std.math.floatMin(F));
            const tol = 256 * std.math.floatEps(F) * @as(F, @floatFromInt(n + 1)) * nrm;
            for (0..n) |j| for (0..n) |i| {
                var s = ops.zero(T);
                for (0..n) |k| for (0..n) |l| {
                    s = ops.add(T, s, ops.mul(T, ops.mul(T, z[i + k * n], t[k + l * n]), ops.conj(T, z[j + l * n])));
                };
                const pair = !isComplex(T) and i == j + 1 and (if (comptime isComplex(T)) false else ev[n + j] > 0);
                if (i > j and !pair and !ops.isZero(T, t[i + j * n])) return error.TestUnexpectedResult;
                const err = ops.abs(T, ops.sub(T, s, a0[i + j * n]));
                if (!(err <= tol)) {
                    check.say("gees {s} n={d} sort={} scale={e}: ({d},{d}) Z T Zᴴ off by {e}\n", .{ @typeName(T), n, sort, sc, i, j, err });
                    return error.TestUnexpectedResult;
                }
            };
            if (!sort) {
                try testing.expectEqual(@as(isize, 0), sdim);
                continue;
            }
            var neg: isize = 0;
            for (0..n) |i| {
                const re = if (comptime isComplex(T)) cw[i].re else ev[i];
                if (re < 0) neg += 1;
                try testing.expectEqual(re < 0, @as(isize, @intCast(i)) < sdim);
            }
            try testing.expectEqual(neg, sdim);
        };
    }
    var a: [4]f64 = .{ 1, std.math.nan(f64), 0, 1 };
    var z: [4]f64 = undefined;
    var wr: [2]f64 = undefined;
    var wi: [2]f64 = undefined;
    var sdim: isize = 0;
    const w: Eigenvalues(f64) = .{ .re = &wr, .im = &wi };
    try testing.expectError(error.InvalidA, evd.gees(al, check.io, f64, .v, null, 2, &a, 2, &sdim, w, &z, 2, null));
    try testing.expectError(error.InvalidN, evd.gees(al, check.io, f64, .v, null, -1, &a, 1, &sdim, w, &z, 1, null));
    try testing.expectError(error.InvalidLda, evd.gees(al, check.io, f64, .n, null, 2, &a, 1, &sdim, w, &z, 1, null));
    try testing.expectError(error.InvalidLdvs, evd.gees(al, check.io, f64, .v, null, 2, &a, 2, &sdim, w, &z, 1, null));
    try evd.gees(al, check.io, f64, .n, null, 0, &a, 1, &sdim, w, &z, 1, null);
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
            const ni: isize = @intCast(n);
            if (use_qr) try evd.syev(al, check.io, T, .n, uplo, ni, a.ptr, ni, w.ptr) else try evd.syevd(al, check.io, T, .n, uplo, ni, a.ptr, ni, w.ptr);
            try evd.syevd(al, check.io, W, .n, uplo, ni, a128.ptr, ni, w128.ptr);
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
