//! Half-precision `GEMV` kernels: `y = alpha*op(A)*x + beta*y`, `A` column-major
//! `m` by `n`. `shgemv` is `f16` in and `f32` out, `sbgemv` is `bf16` in and `f32`
//! out, `bgemv` is `bf16` throughout.
//!
//! The kernels are `ext/bf16/kernel.zig`'s `Gemv`, the same on every core;
//! this file holds the portable install and the tests. See the readme.

const types = @import("../../types.zig");
const kernel = @import("../../ext/bf16/kernel.zig");

const BlasLong = types.BlasLong;
const HF16 = types.HF16;

pub const installGeneric = kernel.installGenericGemv;

// ---- tests -----------------------------------------------------------------

const std = @import("std");
const Vtable = @import("../../Vtable.zig");
const check = @import("../../_test/blas1/check.zig");
const fma = @import("../../fma.zig");
const BF16 = types.BF16;
const families = kernel.families;
const wide = kernel.wide;
const narrow = kernel.narrow;

const sizes = [_][2]usize{ .{ 1, 1 }, .{ 2, 3 }, .{ 3, 2 }, .{ 4, 4 }, .{ 7, 5 }, .{ 8, 8 }, .{ 9, 16 }, .{ 16, 9 }, .{ 17, 17 }, .{ 33, 4 }, .{ 4, 33 } };
const strides = [_][2]usize{ .{ 1, 1 }, .{ 2, 3 }, .{ 3, 2 } };
const alphas = [_]f32{ 1.25, -0.5 };
const betas = [_]f32{ 0, 1, 0.5 };
const cap = 33 * 3 + 2; // a vector's span at stride 3

fn fillAs(comptime T: type, buf: []T, seed: u64) void {
    var f: [400]f32 = undefined;
    check.fill(f32, f[0..buf.len], seed);
    for (buf, f[0..buf.len]) |*d, s| d.* = narrow(T, s);
}

fn checkHalf(comptime E: type, comptime Y: type, comptime trans: bool, kern: anytype) !void {
    for (sizes) |mn| for (strides) |inc| for (alphas) |alpha| for (betas) |beta| {
        const m = mn[0];
        const n = mn[1];
        const lda = m + 1;
        var a: [400]E = undefined; // lda * n at most 18 * 17
        var x: [cap]E = undefined;
        var y: [cap]Y = undefined;
        fillAs(E, &a, m * 7 + n);
        fillAs(E, &x, m + 3);
        fillAs(Y, &y, n + 5);
        const rows = if (trans) n else m;
        const inner = if (trans) m else n;
        // With beta 0, y must never be read: poison it.
        if (beta == 0) for (&y) |*v| {
            v.* = narrow(Y, std.math.nan(f32));
        };
        const orig = y;
        _ = kern(@intCast(m), @intCast(n), narrow(Y, alpha), &a, @intCast(lda), &x, @intCast(inc[0]), narrow(Y, beta), &y, @intCast(inc[1]));
        var touched: [cap]bool = @splat(false);
        for (0..rows) |r| {
            var acc: f64 = 0;
            for (0..inner) |q| {
                const e = if (trans) a[q + r * lda] else a[r + q * lda];
                acc += @as(f64, wide(E, e)) * wide(E, x[q * inc[0]]);
            }
            const yi = r * inc[1];
            touched[yi] = true;
            const prior: f64 = if (beta == 0) 0 else wide(Y, orig[yi]);
            const want = @as(f64, alpha) * acc + @as(f64, beta) * prior;
            const got: f64 = wide(Y, y[yi]);
            const tol: f64 = if (Y == BF16)
                0.005 * @abs(want) + 1e-4
            else
                check.tolerance(f32, 4 * inner + 4);
            if (!(@abs(want - got) <= tol)) {
                check.say("got {e}, want {e} within {e}: m={d} n={d} trans={} alpha={d} beta={d} row={d}: ", .{ got, want, tol, m, n, trans, alpha, beta, r });
                return check.case(error.TestExpectedApproxEqAbs, "half gemv", E, rows, inc);
            }
        }
        for (0..cap) |i| {
            if (touched[i]) continue;
            if (!std.mem.eql(u8, std.mem.asBytes(&y[i]), std.mem.asBytes(&orig[i]))) {
                check.say("untouched y[{d}] changed: m={d} n={d} trans={}: ", .{ i, m, n, trans });
                return check.case(error.TestUnexpectedResult, "half gemv", E, rows, inc);
            }
        }
    };
}

fn checkAll(vt: *const Vtable) !void {
    inline for (families) |f| {
        try checkHalf(f[1], f[2], false, @field(vt, f[0] ++ "_n").?);
        try checkHalf(f[1], f[2], true, @field(vt, f[0] ++ "_t").?);
    }
}

test "portable half-precision gemv kernels match plain arithmetic" {
    var vt: Vtable = .{};
    installGeneric(&vt);
    try checkAll(&vt);
}

test "every linked core's half-precision gemv kernels match plain arithmetic, where this machine runs them" {
    try check.eachCore(checkAll);
}

/// `y` after the kernel, by a loop of a different shape (columns outermost for
/// `n`, every row's accumulator kept in an array), bit for bit.
fn model(comptime E: type, comptime Y: type, comptime trans: bool, comptime fuse: bool, m: usize, n: usize, alpha: Y, a: []const E, lda: usize, x: []const E, incx: usize, beta: Y, y: []Y, incy: usize) void {
    const rows = if (trans) n else m;
    const inner = if (trans) m else n;
    var acc: [40]f32 = @splat(0);
    for (0..inner) |q| for (0..rows) |r| {
        const av = wide(E, if (trans) a[q + r * lda] else a[r + q * lda]);
        const xv = wide(E, x[q * incx]);
        acc[r] = if (fuse) @mulAdd(f32, av, xv, acc[r]) else av * xv + acc[r];
    };
    for (0..rows) |r| {
        const af = wide(Y, alpha);
        const bf = wide(Y, beta);
        const yi = r * incy;
        y[yi] = narrow(Y, if (bf == 0) af * acc[r] else if (fuse) @mulAdd(f32, af, acc[r], bf * wide(Y, y[yi])) else af * acc[r] + bf * wide(Y, y[yi]));
    }
}

fn checkBits(comptime E: type, comptime Y: type, comptime trans: bool, comptime fuse: bool) !void {
    const kern = kernel.Gemv(E, Y, trans, fuse).k;
    for (sizes) |mn| for (strides) |inc| for ([_]f32{ 1.25, -0.5, 0 }) |alpha| for (betas) |beta| {
        const m = mn[0];
        const n = mn[1];
        const lda = m + 2;
        var a: [400]E = undefined;
        var x: [cap]E = undefined;
        var y: [cap]Y = undefined;
        fillAs(E, &a, m * 5 + n);
        fillAs(E, &x, m + 9);
        fillAs(Y, &y, n + 1);
        var want = y;
        model(E, Y, trans, fuse, m, n, narrow(Y, alpha), &a, lda, &x, inc[0], narrow(Y, beta), &want, inc[1]);
        _ = kern(@intCast(m), @intCast(n), narrow(Y, alpha), &a, @intCast(lda), &x, @intCast(inc[0]), narrow(Y, beta), &y, @intCast(inc[1]));
        if (!std.mem.eql(u8, std.mem.asBytes(&want), std.mem.asBytes(&y))) {
            check.say("m={d} n={d} trans={} fuse={} alpha={d} beta={d}: ", .{ m, n, trans, fuse, alpha, beta });
            return check.case(error.TestUnexpectedResult, "half gemv bits", E, m, inc);
        }
    };
}

test "the half-precision gemv kernels are the model, bit for bit, fused or not" {
    inline for (families) |f| inline for (.{ false, true }) |trans| {
        try checkBits(f[1], f[2], trans, false);
        if (comptime fma.hardware) try checkBits(f[1], f[2], trans, true);
    };
}

test "an empty operand leaves beta*y" {
    inline for (.{ false, true }) |trans| {
        const gen = kernel.Gemv(BF16, f32, trans, false).k;
        var a: [4]BF16 = @splat(narrow(BF16, 1));
        var x: [4]BF16 = @splat(narrow(BF16, 1));
        var y0: [4]f32 = .{ 2, 4, 6, 8 };
        _ = gen(if (trans) 0 else 4, if (trans) 4 else 0, 3, &a, 4, &x, 1, 0.5, &y0, 1);
        try std.testing.expectEqualSlices(f32, &.{ 1, 2, 3, 4 }, &y0);
    }
}

test "the row blocks of the non-transposed kernel sum each row as a row-by-row loop does" {
    const al = std.testing.allocator;
    inline for (families) |f| {
        const E = f[1];
        const Y = f[2];
        for ([_][2]usize{ .{ 255, 3 }, .{ 256, 2 }, .{ 257, 5 }, .{ 600, 7 } }) |mn| for ([_][2]usize{ .{ 1, 1 }, .{ 3, 2 } }) |inc| {
            const m = mn[0];
            const n = mn[1];
            const lda = m + 1;
            const a = try al.alloc(E, lda * n);
            defer al.free(a);
            const x = try al.alloc(E, n * inc[0]);
            defer al.free(x);
            const y = try al.alloc(Y, m * inc[1]);
            defer al.free(y);
            const want = try al.alloc(Y, m * inc[1]);
            defer al.free(want);
            var prng = std.Random.DefaultPrng.init(m + n);
            const r = prng.random();
            for (a) |*v| v.* = narrow(E, r.float(f32) * 4 - 2);
            for (x) |*v| v.* = narrow(E, r.float(f32) * 4 - 2);
            for (y, want) |*v, *w| {
                v.* = narrow(Y, r.float(f32) * 4 - 2);
                w.* = v.*;
            }
            const alpha: f32 = 1.25;
            const beta: f32 = 0.5;
            for (0..m) |i| {
                var acc: f32 = 0;
                for (0..n) |j| acc = wide(E, a[i + j * lda]) * wide(E, x[j * inc[0]]) + acc;
                const w = &want[i * inc[1]];
                w.* = narrow(Y, wide(Y, narrow(Y, alpha)) * acc + wide(Y, narrow(Y, beta)) * wide(Y, w.*));
            }
            _ = kernel.Gemv(E, Y, false, false).k(@intCast(m), @intCast(n), narrow(Y, alpha), a.ptr, @intCast(lda), x.ptr, @intCast(inc[0]), narrow(Y, beta), y.ptr, @intCast(inc[1]));
            try std.testing.expect(std.mem.eql(u8, std.mem.sliceAsBytes(want), std.mem.sliceAsBytes(y)));
        };
    }
}
