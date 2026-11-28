//! Tests: `gerfs` by Reference-LAPACK 3.12.0's ratios (`TESTING/LIN/dget07.f`,
//! the error bound and the backward error; `dget04.f`, the forward error), against
//! an exact solution, `dchkge.f`'s threshold.

const std = @import("std");
const types = @import("../../blas/types.zig");
const check = @import("../../blas/_test/check.zig");
const lcheck = @import("check.zig");
const ops = @import("../ops.zig");
const norm = @import("../norm.zig");
const lu = @import("../lu.zig");
const refine = @import("../gerfs.zig");

const testing = std.testing;
const Float = types.Float;
const Trans = types.Trans;
const abs1 = ops.abs1;

/// `dchkge.f`'s `THRESH` in its shipped input (`dtest.in`).
const thresh = 30;

/// `op(A)(i, k)` of the `n x n` `a`.
pub fn opAt(comptime T: type, trans: Trans, a: []const T, n: usize, i: usize, k: usize) T {
    return switch (trans) {
        .n => a[i + k * n],
        .t => a[k + i * n],
        else => ops.conj(T, a[k + i * n]),
    };
}

/// The ratios of one solve: `dget07`'s two and `dget04`'s, each to stay
/// under `thresh`.
pub fn ratios(comptime T: type, al: std.mem.Allocator, trans: Trans, n: usize, nrhs: usize, a: []const T, b: []const T, x: []const T, xact: []const T, ferr: []const Float(T), berr: []const Float(T)) ![3]Float(T) {
    const F = Float(T);
    const eps = std.math.floatEps(F);
    const unfl = std.math.floatMin(F);
    // `dget04`'s `rcond`: from the inverse, as `dchkge.f`'s `RCONDC`
    const inv = try al.dupe(T, a);
    defer al.free(inv);
    const ipiv = try al.alloc(usize, n);
    defer al.free(ipiv);
    _ = try lu.getrf(al, check.io, T, n, n, inv, n, ipiv);
    try lu.getri(al, check.io, T, n, inv, n, ipiv);
    const kind: norm.NormType = if (trans == .n) .one else .infinity;
    const rcond = 1 / (try norm.lange(al, check.io, T, kind, n, n, a, n) * try norm.lange(al, check.io, T, kind, n, n, inv, n));
    var out: [3]F = @splat(0);
    for (0..nrhs) |j| {
        const xj = x[j * n ..][0..n];
        const ej = xact[j * n ..][0..n];
        var xnorm: F = unfl;
        var enorm: F = 0;
        var diff: F = 0;
        for (xj, ej) |u, v| {
            xnorm = @max(xnorm, abs1(T, u));
            enorm = @max(enorm, abs1(T, v));
            diff = @max(diff, abs1(T, ops.sub(T, u, v)));
        }
        out[0] = @max(out[0], if (diff / xnorm <= ferr[j]) (diff / xnorm) / ferr[j] else 1 / eps);
        var axbi: F = std.math.inf(F);
        for (0..n) |i| {
            var t = abs1(T, b[i + j * n]);
            for (0..n) |k| t += abs1(T, opAt(T, trans, a, n, i, k)) * abs1(T, xj[k]);
            axbi = @min(axbi, t);
        }
        const nz: F = @floatFromInt(n + 1);
        out[1] = @max(out[1], berr[j] / (nz * eps + nz * unfl / @max(axbi, nz * unfl)));
        out[2] = @max(out[2], diff / enorm * rcond / eps);
    }
    return out;
}

/// A diagonally dominant `A`, an exact `X`, `B = op(A) X`, and `getrs`'s
/// solution, perturbed by `perturb` relative to each entry when nonzero,
/// refined by `gerfs`.
fn case(comptime T: type, trans: Trans, n: usize, nrhs: usize, perturb: Float(T)) !void {
    const al = testing.allocator;
    const a = try al.alloc(T, n * n);
    defer al.free(a);
    check.fill(T, a, n * 31 + nrhs);
    for (0..n) |i| a[i + i * n] = ops.add(T, a[i + i * n], ops.scalar(T, 2));
    const xact = try al.alloc(T, n * nrhs);
    defer al.free(xact);
    check.fill(T, xact, n * 7 + 3);
    const b = try al.alloc(T, n * nrhs);
    defer al.free(b);
    for (0..nrhs) |j| for (0..n) |i| {
        var t = ops.zero(T);
        for (0..n) |k| t = ops.add(T, t, ops.mul(T, opAt(T, trans, a, n, i, k), xact[k + j * n]));
        b[i + j * n] = t;
    };
    const af = try al.dupe(T, a);
    defer al.free(af);
    const ipiv = try al.alloc(usize, n);
    defer al.free(ipiv);
    try testing.expectEqual(null, try lu.getrf(al, check.io, T, n, n, af, n, ipiv));
    const x = try al.dupe(T, b);
    defer al.free(x);
    try lu.getrs(al, check.io, T, trans, n, nrhs, af, n, ipiv, x, n);
    for (x) |*e| e.* = ops.mulReal(T, e.*, 1 + perturb);
    const ferr = try al.alloc(Float(T), nrhs);
    defer al.free(ferr);
    const berr = try al.alloc(Float(T), nrhs);
    defer al.free(berr);
    try refine.gerfs(al, check.io, T, trans, n, nrhs, a, n, af, n, ipiv, b, n, x, n, ferr, berr);
    const r = try ratios(T, al, trans, n, nrhs, a, b, x, xact, ferr, berr);
    for (r, 0..) |v, k| if (!(v < thresh)) {
        check.say("gerfs {s} trans {t} n {d}: ratio {d} is {e}\n", .{ @typeName(T), trans, n, k, v });
        return error.TestUnexpectedResult;
    };
}

test "gerfs: dget07's and dget04's ratios, every type and trans, refined from getrs and from a perturbed solution" {
    inline for (lcheck.test_types) |T| for ([_]Trans{ .n, .t, .c }) |trans| for ([_]usize{ 1, 2, 5, 23 }) |n| {
        try case(T, trans, n, 3, 0);
        try case(T, trans, n, 2, 1e-3);
    };
}

test "gerfs: empty sizes, bad arguments" {
    const al = testing.allocator;
    var a = [_]f64{ 2, 1, 1, 3 };
    const ipiv = [_]usize{ 0, 1 };
    var b = [_]f64{ 1, 2 };
    var x = [_]f64{ 1, 2 };
    var ferr = [_]f64{ 7, 7 };
    var berr = [_]f64{ 7, 7 };
    try refine.gerfs(al, check.io, f64, .n, 0, 2, a[0..0], 1, a[0..0], 1, ipiv[0..0], b[0..0], 1, x[0..0], 1, &ferr, &berr);
    try testing.expectEqualSlices(f64, &.{ 0, 0 }, &ferr);
    try testing.expectEqualSlices(f64, &.{ 0, 0 }, &berr);
    try testing.expectError(error.UnsupportedTranspose, refine.gerfs(al, check.io, f64, .r, 2, 1, &a, 2, &a, 2, &ipiv, &b, 2, &x, 2, &ferr, &berr));
    try testing.expectError(error.PivotOutOfRange, refine.gerfs(al, check.io, f64, .n, 2, 1, &a, 2, &a, 2, &.{ 0, 2 }, &b, 2, &x, 2, &ferr, &berr));
    try testing.expectError(error.LeadingDimensionTooSmall, refine.gerfs(al, check.io, f64, .n, 2, 1, &a, 2, &a, 1, &ipiv, &b, 2, &x, 2, &ferr, &berr));
    try testing.expectError(error.SliceTooShort, refine.gerfs(al, check.io, f64, .n, 2, 1, &a, 2, &a, 2, &ipiv, &b, 2, &x, 2, ferr[0..0], &berr));
}
