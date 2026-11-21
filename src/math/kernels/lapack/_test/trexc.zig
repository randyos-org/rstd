//! Tests of `../trexc.zig` and `../trsen.zig`.

const std = @import("std");
const types = @import("../../blas/types.zig");
const ops = @import("../ops.zig");
const schur = @import("../schur.zig");
const trexc_move = @import("../trexc.zig");
const trsen_reorder = @import("../trsen.zig");

const Allocator = std.mem.Allocator;
const Float = types.Float;
const isComplex = types.isComplex;
const Eigenvalues = schur.Eigenvalues;
const trexc = trexc_move.trexc;
const trsen = trsen_reorder.trsen;

const testing = std.testing;
const check = @import("../../blas/_test/check.zig");
const lcheck = @import("check.zig");
const lu = @import("../lu.zig");
const svd = @import("../svd.zig");

const test_types = lcheck.test_types;

/// `THRESH` of the reference's eigenvalue-condition tests
/// (`TESTING/dec.in`, `dchkec.f`'s `dget36` and `dget38`).
const thresh = 20;

/// An `n x n` Schur form: upper triangular with diagonal `1, 2, ...`
/// (complex: plus `i` times as much), real ones with a standardized 2 x 2
/// block `[k b; c k]`, `b c < 0`, wherever `pairs` has a bit set and room
/// remains; random above.
fn schurForm(comptime T: type, al: Allocator, n: usize, pairs: u32, seed: u64) ![]T {
    const F = Float(T);
    const t = try al.alloc(T, n * n);
    check.fill(T, t, seed);
    for (0..n) |j| for (j..n) |i| {
        t[i + j * n] = ops.zero(T);
    };
    var k: usize = 0;
    while (k < n) {
        const v: F = @floatFromInt(k + 1);
        if (!isComplex(T) and k + 1 < n and (pairs >> @intCast(k % 32)) & 1 != 0) {
            t[k + k * n] = ops.scalar(T, v);
            t[k + 1 + (k + 1) * n] = ops.scalar(T, v);
            t[k + (k + 1) * n] = ops.scalar(T, 1.5);
            t[k + 1 + k * n] = ops.scalar(T, -0.75);
            k += 2;
        } else {
            t[k + k * n] = if (comptime isComplex(T)) .init(v, v / 2) else v;
            k += 1;
        }
    }
    return t;
}

/// `Q`'s orthonormality and `Q T Qᴴ = T0` within the threshold of `n eps`
/// (`dget36.f`'s two residuals), and `T` still in Schur form.
fn expectSimilar(comptime T: type, n: usize, t0: []const T, t: []const T, q: []const T) !void {
    const F = Float(T);
    const eps = std.math.floatEps(F);
    const nf: F = @floatFromInt(@max(n, 1));
    var tnorm: F = 0;
    for (t0) |v| tnorm = @max(tnorm, ops.abs(T, v));
    for (0..n) |i| for (0..n) |j| {
        var g = ops.zero(T);
        var r = ops.zero(T);
        for (0..n) |k| {
            g = ops.add(T, g, ops.mul(T, ops.conj(T, q[k + i * n]), q[k + j * n]));
            var s = ops.zero(T);
            for (0..n) |l| s = ops.add(T, s, ops.mul(T, t[k + l * n], ops.conj(T, q[j + l * n])));
            r = ops.add(T, r, ops.mul(T, q[i + k * n], s));
        }
        if (i == j) g = ops.sub(T, g, ops.one(T));
        try testing.expect(ops.abs(T, g) <= thresh * nf * eps);
        try testing.expect(ops.abs(T, ops.sub(T, r, t0[i + j * n])) <= thresh * nf * eps * tnorm);
    };
    for (0..n) |j| for (j + 1..n) |i| {
        const v = t[i + j * n];
        if (comptime isComplex(T)) {
            try testing.expect(ops.isZero(T, v));
        } else if (i > j + 1) {
            try testing.expect(ops.isZero(T, v));
        } else if (!ops.isZero(T, v)) {
            // a standardized block, alone on the subdiagonal
            try testing.expect(t[j + j * n] == t[i + i * n]);
            try testing.expect(std.math.signbit(ops.re(T, v)) != std.math.signbit(ops.re(T, t[j + i * n])));
            if (j > 0) try testing.expect(ops.isZero(T, t[j + (j - 1) * n]));
        }
    };
}

/// The eigenvalue (real: of the block) starting at 0-based `k`.
fn eigAt(comptime T: type, n: usize, t: []const T, k: usize) [2]f64 {
    const v = t[k + k * n];
    if (comptime isComplex(T)) return .{ @floatCast(v.re), @floatCast(v.im) };
    if (k + 1 < n and t[k + 1 + k * n] != 0) return .{ @floatCast(v), @floatCast(@sqrt(@abs(t[k + (k + 1) * n])) * @sqrt(@abs(t[k + 1 + k * n]))) };
    return .{ @floatCast(v), 0 };
}

test "lapack trexc: moves every block to every place, Q T Q' = T0, Schur form kept" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        for ([_]usize{ 1, 2, 3, 5, 8 }) |n| for ([_]u32{ 0, 0b1, 0b10, 0b1001, 0b100101 }) |pairs| {
            const t0 = try schurForm(T, al, n, pairs, n * 3 + pairs);
            defer al.free(t0);
            const t = try al.alloc(T, n * n);
            defer al.free(t);
            const q = try al.alloc(T, n * n);
            defer al.free(q);
            for (0..n) |f| for (0..n) |l| {
                @memcpy(t, t0);
                @memset(q, ops.zero(T));
                for (0..n) |i| q[i * (n + 1)] = ops.one(T);
                var ifst = f;
                var ilst = l;
                try trexc(check.io, T, .v, n, t, n, q, n, &ifst, &ilst);
                try expectSimilar(T, n, t0, t, q);
                // the moved eigenvalue sits at ilst
                const want = eigAt(T, n, t0, ifst);
                const got = eigAt(T, n, t, ilst);
                const tol = 100 * @as(f64, @floatCast(std.math.floatEps(F))) * @as(f64, @floatFromInt(n));
                try testing.expectApproxEqAbs(want[0], got[0], tol);
                try testing.expectApproxEqAbs(@abs(want[1]), @abs(got[1]), tol);
            };
        };
    }
}

test "lapack trexc: n = 0 and 1, compq n, bad arguments" {
    inline for (test_types) |T| {
        var t = [_]T{ops.scalar(T, 3)};
        var q = [_]T{ops.scalar(T, 7)};
        var f: usize = 0;
        var l: usize = 0;
        try trexc(check.io, T, .v, 1, &t, 1, &q, 1, &f, &l);
        try trexc(check.io, T, .n, 0, &t, 1, &q, 1, &f, &l);
        try testing.expect(ops.re(T, t[0]) == 3 and ops.re(T, q[0]) == 7);
    }
    var t: [4]f64 = .{ 1, 0, 2, 3 };
    var q: [4]f64 = undefined;
    var f: usize = 0;
    var l: usize = 1;
    try trexc(check.io, f64, .n, 2, &t, 2, &q, 1, &f, &l);
    try testing.expect(t[0] == 3 and t[3] == 1 and @abs(t[2]) == 2 and t[1] == 0);
    try testing.expectError(error.LeadingDimensionTooSmall, trexc(check.io, f64, .n, 2, &t, 1, &q, 1, &f, &l));
    try testing.expectError(error.SliceTooShort, trexc(check.io, f64, .n, 2, t[0..3], 2, &q, 1, &f, &l));
    try testing.expectError(error.LeadingDimensionTooSmall, trexc(check.io, f64, .v, 2, &t, 2, &q, 1, &f, &l));
    f = 2;
    try testing.expectError(error.RowOutOfRange, trexc(check.io, f64, .n, 2, &t, 2, &q, 1, &f, &l));
    f = 0;
    l = 2;
    try testing.expectError(error.RowOutOfRange, trexc(check.io, f64, .n, 2, &t, 2, &q, 1, &f, &l));
}

/// `I ⊗ T11 - T22ᵀ ⊗ I`, the Sylvester operator `trsen` measures, on
/// column-major `vec(X)`.
fn kron(comptime T: type, al: Allocator, n: usize, t: []const T, n1: usize) ![]T {
    const n2 = n - n1;
    const nn = n1 * n2;
    const k = try al.alloc(T, nn * nn);
    @memset(k, ops.zero(T));
    for (0..n2) |j| for (0..n1) |i| {
        const row = i + j * n1;
        for (0..n1) |p| k[row + (p + j * n1) * nn] = ops.add(T, k[row + (p + j * n1) * nn], t[i + p * n]);
        for (0..n2) |r| k[row + (i + r * n1) * nn] = ops.sub(T, k[row + (i + r * n1) * nn], t[n1 + r + (n1 + j) * n]);
    };
    return k;
}

test "lapack trsen: the selected cluster leads, s against the Kronecker solve, sep against its smallest singular value" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        const eps = std.math.floatEps(F);
        for ([_]usize{ 2, 4, 7 }) |n| for ([_]u32{ 0, 0b101 }) |pairs| for ([_]u32{ 0b1, 0b10, 0b110, 0b1010010 }) |sel| {
            const t0 = try schurForm(T, al, n, pairs, n * 5 + sel);
            defer al.free(t0);
            const t = try al.alloc(T, n * n);
            defer al.free(t);
            @memcpy(t, t0);
            const q = try al.alloc(T, n * n);
            defer al.free(q);
            @memset(q, ops.zero(T));
            for (0..n) |i| q[i * (n + 1)] = ops.one(T);
            var select: [8]bool = undefined;
            for (0..n) |i| select[i] = (sel >> @intCast(i)) & 1 != 0;
            var wbuf: [16]T = undefined;
            const w: Eigenvalues(T) = if (comptime isComplex(T)) &wbuf else .{ .re = wbuf[0..8], .im = wbuf[8..] };
            const r = try trsen(al, check.io, T, .b, .v, select[0..n], n, t, n, q, n, w);
            const s = r.s.?;
            const sep = r.sep.?;
            try expectSimilar(T, n, t0, t, q);
            const n1 = r.m;
            // the selected eigenvalues, in their old order, lead
            var k: usize = 0;
            var got: usize = 0;
            while (k < n) {
                const e = eigAt(T, n, t0, k);
                const size: usize = if (!isComplex(T) and e[1] != 0) 2 else 1;
                if (select[k] or (size == 2 and select[k + 1])) {
                    const g = eigAt(T, n, t, got);
                    try testing.expectApproxEqAbs(e[0], g[0], 1e3 * @as(f64, @floatCast(eps)));
                    got += size;
                }
                k += size;
            }
            try testing.expectEqual(got, n1);
            if (n1 == 0 or n1 == n) {
                try testing.expect(s == 1);
                continue;
            }
            const n2 = n - n1;
            const nn = n1 * n2;
            const kk = try kron(T, al, n, t, n1);
            defer al.free(kk);
            // s = 1 / sqrt(1 + |X|_F^2) for K vec(X) = vec(T12)
            const x = try al.alloc(T, nn);
            defer al.free(x);
            for (0..n2) |j| for (0..n1) |i| {
                x[i + j * n1] = t[i + (n1 + j) * n];
            };
            const kf = try al.alloc(T, nn * nn);
            defer al.free(kf);
            @memcpy(kf, kk);
            const ipiv = try al.alloc(usize, nn);
            defer al.free(ipiv);
            try lu.gesv(al, check.io, T, nn, 1, kf, nn, ipiv, x, nn);
            var xf: F = 0;
            for (x) |v| xf += ops.abs2(T, v);
            const s_true = 1 / @sqrt(1 + xf);
            try testing.expectApproxEqRel(s_true, s, 1e4 * eps);
            // sep within sqrt(nn) and the threshold of sigma_min(K)
            @memcpy(kf, kk);
            const sv = try al.alloc(F, nn);
            defer al.free(sv);
            try svd.gesvd(al, check.io, T, .n, .n, nn, nn, kf, nn, sv, kf[0..0], 1, kf[0..0], 1, null);
            const smin = sv[nn - 1];
            const root = @sqrt(@as(F, @floatFromInt(nn)));
            if (!(sep >= smin / root * (1 - 1e3 * eps) and sep <= thresh * root * smin)) {
                check.say("trsen {s} n {d} select {b}: sep {e}, sigma_min {e}\n", .{ @typeName(T), n, sel, @as(f64, @floatCast(sep)), @as(f64, @floatCast(smin)) });
                return error.TestUnexpectedResult;
            }
        };
    }
}

test "lapack trsen: none or all selected, a block-diagonal T, n = 0, bad arguments" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        // diag(1, 4): s = 1 (T12 = 0), sep = |1 - 4| exactly
        var t = [_]T{ ops.scalar(T, 1), ops.zero(T), ops.zero(T), ops.scalar(T, 4) };
        var q: [4]T = undefined;
        var wbuf: [4]T = undefined;
        const w: Eigenvalues(T) = if (comptime isComplex(T)) &wbuf else .{ .re = wbuf[0..2], .im = wbuf[2..] };
        var r = try trsen(al, check.io, T, .b, .n, &.{ false, true }, 2, &t, 2, &q, 1, w);
        try testing.expect(r.m == 1 and r.s.? == 1 and r.sep.? == 3);
        try testing.expect(ops.re(T, t[0]) == 4 and ops.re(T, t[3]) == 1);
        r = try trsen(al, check.io, T, .b, .n, &.{ false, false }, 2, &t, 2, &q, 1, w);
        try testing.expect(r.m == 0 and r.s.? == 1 and r.sep.? == 4);
        r = try trsen(al, check.io, T, .v, .n, &.{ true, true }, 2, &t, 2, &q, 1, w);
        try testing.expect(r.m == 2 and r.s == null and r.sep.? == 4);
        r = try trsen(al, check.io, T, .b, .n, &.{}, 0, &t, 1, &q, 1, w);
        try testing.expect(r.m == 0 and r.s.? == 1);
        r = try trsen(al, check.io, T, .n, .n, &.{ true, false }, 2, &t, 2, &q, 1, w);
        try testing.expect(r.m == 1 and r.s == null and r.sep == null);
    }
    var t: [4]f64 = undefined;
    var q: [4]f64 = undefined;
    var wr: [2]f64 = undefined;
    var wi: [2]f64 = undefined;
    const sel = [_]bool{ true, false };
    const w: Eigenvalues(f64) = .{ .re = &wr, .im = &wi };
    try testing.expectError(error.LeadingDimensionTooSmall, trsen(al, check.io, f64, .n, .n, &sel, 2, &t, 1, &q, 1, w));
    try testing.expectError(error.SliceTooShort, trsen(al, check.io, f64, .n, .n, &sel, 2, t[0..3], 2, &q, 1, w));
    try testing.expectError(error.LeadingDimensionTooSmall, trsen(al, check.io, f64, .n, .v, &sel, 2, &t, 2, &q, 1, w));
    try testing.expectError(error.SliceTooShort, trsen(al, check.io, f64, .n, .n, sel[0..1], 2, &t, 2, &q, 1, w));
    try testing.expectError(error.SliceTooShort, trsen(al, check.io, f64, .n, .n, &sel, 2, &t, 2, &q, 1, .{ .re = &wr, .im = wi[0..1] }));
}
