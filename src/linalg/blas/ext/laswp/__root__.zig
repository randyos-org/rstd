//! Tests of `kernel.zig`'s row interchanges and `laswp_ncopy`, which the gemm
//! kernels' sets install (`gemm/kernel.zig`).

const std = @import("std");
const types = @import("../../types.zig");
const Vtable = @import("../../Vtable.zig");
const portable = @import("../../portable.zig");
const check = @import("../../_test/check.zig");
const bitwise = @import("../../_test/bitwise.zig");
const pattern = @import("../../_test/pattern.zig");
const kernel = @import("kernel.zig");

const XDouble = types.XDouble;
const Complex32 = types.Complex32;
const Complex64 = types.Complex64;
const ComplexXD = types.ComplexXD;

// ---- tests -----------------------------------------------------------------

/// Reference-LAPACK's `dlaswp.f` as written: `IX0`, `I1`, `I2` and `INC`.
pub fn laswpRef(comptime T: type, n: usize, a: []T, lda: usize, k1: isize, k2: isize, ipiv: []const isize, incx: isize) void {
    if (incx == 0) return;
    const ix0 = if (incx > 0) k1 else k1 + (k1 - k2) * incx;
    const inc: isize = if (incx > 0) 1 else -1;
    for (0..n) |j| {
        var ix = ix0;
        var i = if (incx > 0) k1 else k2;
        while (if (incx > 0) i <= k2 else i >= k1) : (i += inc) {
            const ip = ipiv[@intCast(ix - 1)];
            if (ip != i) std.mem.swap(T, &a[@as(usize, @intCast(i - 1)) + j * lda], &a[@as(usize, @intCast(ip - 1)) + j * lda]);
            ix += incx;
        }
    }
}

/// Pivots as `getrf` makes them (`ipiv(k) >= k`), with repeats and rows
/// swapped with themselves, for rows `1..m`.
pub fn pivots(ipiv: []isize, m: usize, seed: u64) void {
    var prng = std.Random.DefaultPrng.init(seed);
    const r = prng.random();
    for (ipiv, 0..) |*p, i| {
        const k = i % m;
        p.* = @intCast(k + 1 + r.uintLessThan(usize, m - k));
    }
}

/// The precisions with a `laswp_ncopy` slot; `q` and `x` where they are
/// types of their own.
const precisions = .{ .{ "s", f32 }, .{ "d", f64 }, .{ "c", Complex32 }, .{ "z", Complex64 } } ++
    (if (XDouble == f64) .{} else .{ .{ "q", XDouble }, .{ "x", ComplexXD } });

/// `vt`'s `laswp_ncopy` against the interchanges by `laswpRef` and then its
/// `gemm_oncopy` of rows `k1` to `k2`, bit for bit, and `A` swapped as the
/// reference swaps it.
fn checkAll(vt: *const Vtable) !void {
    const al = std.testing.allocator;
    inline for (precisions) |p| {
        const T = p[1];
        const F = types.Float(T);
        const ncopy = portable.slotPtr(vt, .from(p[0]), .laswp_ncopy).*.?;
        const oncopy = portable.slotPtr(vt, .from(p[0]), .gemm_oncopy).*.?;
        for ([_]usize{ 1, 2, 3, 5, 13, 40 }) |n| for ([_][2]usize{ .{ 1, 1 }, .{ 1, 4 }, .{ 3, 9 }, .{ 2, 2 } }) |ks| {
            const m = 11;
            const lda = m + 2;
            const k1 = ks[0];
            const k2 = ks[1];
            const a = try al.alloc(T, lda * n);
            defer al.free(a);
            check.fill(T, a, n * 3 + k2);
            for (a, 0..) |*v, i| v.* = pattern.shape(T, v.*, i, .mixed);
            var ipiv: [m]isize = undefined;
            pivots(&ipiv, m, n + k1);
            const want = try al.dupe(T, a);
            defer al.free(want);
            laswpRef(T, n, want, lda, @intCast(k1), @intCast(k2), &ipiv, 1);
            const rows = k2 - k1 + 1;
            const want_b = try al.alloc(T, rows * n);
            defer al.free(want_b);
            const got_b = try al.alloc(T, rows * n);
            defer al.free(got_b);
            _ = oncopy(@intCast(rows), @intCast(n), @ptrCast(want.ptr + (k1 - 1)), lda, @ptrCast(want_b.ptr));
            _ = ncopy(@intCast(n), @intCast(k1), @intCast(k2), @as([*]F, @ptrCast(a.ptr)), lda, &ipiv, @ptrCast(got_b.ptr));
            for (want, a) |w, g| if (!bitwise.same(T, w, g, true)) {
                check.say("{s}laswp_ncopy n={d} k1={d} k2={d}, A: ", .{ p[0], n, k1, k2 });
                return error.TestUnexpectedResult;
            };
            for (want_b, got_b) |w, g| if (!bitwise.same(T, w, g, true)) {
                check.say("{s}laswp_ncopy n={d} k1={d} k2={d}, buffer: ", .{ p[0], n, k1, k2 });
                return error.TestUnexpectedResult;
            };
        };
    }
}

test "portable laswp_ncopy is the interchanges and gemm_oncopy" {
    var vt: Vtable = .{};
    portable.install(&vt);
    try checkAll(&vt);
}

test "every linked core's laswp_ncopy is the interchanges and its gemm_oncopy, where this machine runs it" {
    try check.eachCore(checkAll);
}

test "the interchanges follow dlaswp.f for either sign of incx, pivots repeated" {
    const al = std.testing.allocator;
    inline for (.{ f64, Complex32, f16 }) |T| {
        for ([_]isize{ 1, 2, -1, -3 }) |incx| for ([_][2]isize{ .{ 1, 7 }, .{ 3, 5 }, .{ 4, 4 } }) |ks| {
            const m = 7;
            const n = 3;
            const lda = m + 1;
            const step: usize = @abs(incx);
            var ipiv: [m * 3]isize = undefined;
            pivots(&ipiv, m, step);
            const a = try al.alloc(T, lda * n);
            defer al.free(a);
            check.fill(T, a, 9);
            const want = try al.dupe(T, a);
            defer al.free(want);
            laswpRef(T, n, want, lda, ks[0], ks[1], &ipiv, incx);
            const k1: usize = @intCast(ks[0]);
            kernel.swapRows(T, n, a.ptr, lda, k1, @intCast(ks[1]), ipiv[k1 - 1 ..].ptr, incx);
            for (want, a) |w, g| try std.testing.expect(bitwise.same(T, w, g, true));
        };
    }
}
