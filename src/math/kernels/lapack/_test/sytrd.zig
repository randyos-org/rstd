//! Tests of `sytrd.zig`: `Qᴴ A Q = T` rebuilt from the reflectors, and the
//! eigenvalues against faer's reduction, around the block and crossover.

const std = @import("std");
const check = @import("../../blas/_test/check.zig");
const types = @import("../../blas/types.zig");
const lcheck = @import("check.zig");
const ops = @import("../ops.zig");
const sytrd = @import("../sytrd.zig");
const evd = @import("../evd.zig");

const testing = std.testing;
const Mat = ops.Mat;
const Float = ops.Float;

const sizes = [_]usize{ 1, 2, 3, 5, 31, 32, 33, 34, 63, 64, 65, 97, 130 };

/// `H(i) = I - tau(i) v vᴴ` applied from both sides to `x` (`n x n`) in
/// order `i = 0, 1, ...`: `x := H(i)ᴴ x H(i)`, so `x` goes from `A` to `T`.
fn applyBoth(comptime T: type, x: Mat(T), a: Mat(T), tau: []const T) void {
    const n = x.rows;
    const v = std.heap.page_allocator.alloc(T, n) catch unreachable;
    defer std.heap.page_allocator.free(v);
    for (0..n -| 1) |i| {
        @memset(v, ops.zero(T));
        v[i + 1] = ops.one(T);
        for (i + 2..n) |r| v[r] = a.get(r, i);
        const t = tau[i];
        // x := (I - conj(t) v vᴴ) x
        for (0..n) |j| {
            var s = ops.zero(T);
            for (0..n) |r| s = ops.add(T, s, ops.mul(T, ops.conj(T, v[r]), x.get(r, j)));
            const f = ops.mul(T, ops.conj(T, t), s);
            for (0..n) |r| x.at(r, j).* = ops.sub(T, x.get(r, j), ops.mul(T, v[r], f));
        }
        // x := x (I - t v vᴴ)
        for (0..n) |r| {
            var s = ops.zero(T);
            for (0..n) |j| s = ops.add(T, s, ops.mul(T, x.get(r, j), v[j]));
            const f = ops.mul(T, s, t);
            for (0..n) |j| x.at(r, j).* = ops.sub(T, x.get(r, j), ops.mul(T, f, ops.conj(T, v[j])));
        }
    }
}

test "sytrd: Qᴴ A Q is the tridiagonal of d and e, every type, around the block" {
    const al = testing.allocator;
    inline for (.{ f64, f32, types.Complex64, types.Complex32 }) |T| {
        const F = Float(T);
        for (sizes) |n| {
            const a0 = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, a0);
            lcheck.hermitian(T, a0, 7 * n + 1);
            const a = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, a);
            a.copyFrom(a0);
            const d = try al.alloc(F, n);
            defer al.free(d);
            const e = try al.alloc(F, @max(n, 1) - 1);
            defer al.free(e);
            const tau = try al.alloc(T, @max(n, 1) - 1);
            defer al.free(tau);
            try sytrd.sytrd(al, check.io, T, a, d, e, tau);
            const x = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, x);
            x.copyFrom(a0);
            // the full Hermitian `A`, from its lower triangle
            for (0..n) |j| for (0..j) |i| {
                x.at(i, j).* = ops.conj(T, x.get(j, i));
            };
            applyBoth(T, x, a, tau);
            const tol = 64 * std.math.floatEps(F) * @as(F, @floatFromInt(n + 1)) * @max(normLower(T, a0), 1);
            for (0..n) |j| for (0..n) |i| {
                const want: T = if (i == j) ops.scalar(T, d[i]) else if (i == j + 1) ops.scalar(T, e[j]) else if (j == i + 1) ops.scalar(T, e[i]) else ops.zero(T);
                const err = ops.abs(T, ops.sub(T, x.get(i, j), want));
                if (!(err <= tol)) {
                    check.say("sytrd {s} n={d}: (Qᴴ A Q)({d},{d}) off by {e}, allowed {e}\n", .{ @typeName(T), n, i, j, err, tol });
                    return error.TestUnexpectedResult;
                }
            };
        }
    }
}

fn normLower(comptime T: type, a: Mat(T)) Float(T) {
    var m: Float(T) = 0;
    for (0..a.rows) |j| for (j..a.rows) |i| {
        m = @max(m, ops.abs(T, a.get(i, j)));
    };
    return m;
}

test "sytrd: syevd's values alone agree with its values beside the vectors" {
    const al = testing.allocator;
    inline for (.{ f64, f32, types.Complex64, types.Complex32 }) |T| {
        const F = Float(T);
        for ([_]usize{ 3, 33, 64, 130, 260 }) |n| {
            const a0 = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, a0);
            lcheck.hermitian(T, a0, 11 * n + 3);
            const a = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, a);
            const w_vals = try al.alloc(F, n);
            defer al.free(w_vals);
            const w_vecs = try al.alloc(F, n);
            defer al.free(w_vecs);
            a.copyFrom(a0);
            try evd.symmetric(al, check.io, T, false, true, a, w_vals, 25);
            a.copyFrom(a0);
            try evd.symmetric(al, check.io, T, true, true, a, w_vecs, 25);
            const tol = 64 * std.math.floatEps(F) * @as(F, @floatFromInt(n)) * @max(normLower(T, a0), 1);
            for (w_vals, w_vecs, 0..) |x, y, i| {
                if (!(@abs(x - y) <= tol)) {
                    check.say("syevd {s} n={d}: value {d} {e} alone, {e} with vectors\n", .{ @typeName(T), n, i, x, y });
                    return error.TestUnexpectedResult;
                }
            }
        }
    }
}
