//! Tests of `kernel.zig`'s `neg` kernels (`neg_tcopy`), which the gemm
//! kernels' sets install (`gemm/kernel.zig`).

const std = @import("std");
const blas = @import("../../__root__.zig");
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

/// The precisions with a `neg_tcopy` slot; `q` and `x` where they are types
/// of their own.
const precisions = .{ .{ "s", f32 }, .{ "d", f64 }, .{ "c", Complex32 }, .{ "z", Complex64 } } ++
    (if (XDouble == f64) .{} else .{ .{ "q", XDouble }, .{ "x", ComplexXD } });

/// `vt`'s `neg_tcopy` against its `gemm_itcopy` negated, bit for bit, NaN and
/// infinities included: the layout is the gemm kernels', whatever the tile.
fn checkAll(vt: *const Vtable) !void {
    const al = std.testing.allocator;
    inline for (precisions) |p| {
        const T = p[1];
        const F = types.Float(T);
        const neg = portable.slotPtr(vt, .from(p[0]), .neg_tcopy).*.?;
        const tcopy = portable.slotPtr(vt, .from(p[0]), .gemm_itcopy).*.?;
        for ([_]usize{ 1, 2, 3, 7, 17 }) |m| for ([_]usize{ 1, 2, 3, 5, 13, 33, 70 }) |n| {
            const lda = n + 2;
            const a = try al.alloc(T, lda * m);
            defer al.free(a);
            const want = try al.alloc(T, m * n);
            defer al.free(want);
            const got = try al.alloc(T, m * n);
            defer al.free(got);
            check.fill(T, a, m * 101 + n);
            // signed zeros, denormals, huge values, infinities and NaN
            for (a, 0..) |*v, i| v.* = pattern.shape(T, v.*, i, .mixed);
            const pa: [*]F = @ptrCast(a.ptr);
            _ = tcopy(@intCast(m), @intCast(n), pa, @intCast(lda), @ptrCast(want.ptr));
            _ = neg(@intCast(m), @intCast(n), pa, @intCast(lda), @ptrCast(got.ptr));
            for (want, got) |w, g| {
                const nw = if (comptime types.isComplex(T)) T.init(-w.re, -w.im) else -w;
                if (!bitwise.same(T, nw, g, true)) {
                    check.say("{s}neg_tcopy m={d} n={d}: ", .{ p[0], m, n });
                    return error.TestUnexpectedResult;
                }
            }
        };
    }
}

test "portable neg_tcopy is gemm_itcopy negated" {
    var vt: Vtable = .{};
    portable.install(&vt);
    try checkAll(&vt);
}

test "every linked core's neg_tcopy is its gemm_itcopy negated, where this machine runs it" {
    try check.eachCore(checkAll);
}

test "neg_tcopy at unroll 1 is neg_tcopy_1.c" {
    inline for (.{ f64, Complex64 }) |T| {
        const K = kernel.Kernels(T, .{ .vl = 1, .mv = 1, .nr = 1, .split = types.isComplex(T) });
        const m = 3;
        const n = 4;
        const lda = 6;
        var a: [lda * m]T = undefined;
        check.fill(T, &a, 5);
        var b: [m * n]T = undefined;
        _ = K.neg_tcopy(m, n, @ptrCast(&a), lda, @ptrCast(&b));
        // b[i + j*m] = -a[j + i*lda]
        for (0..m) |i| for (0..n) |j| {
            const v = a[j + i * lda];
            const want: T = if (comptime types.isComplex(T)) T.init(-v.re, -v.im) else -v;
            try std.testing.expect(bitwise.same(T, want, b[i + j * m], true));
        };
    }
}

test "sbneg_tcopy is sneg_tcopy" {
    try std.testing.expect(blas.vtable().sbneg_tcopy.? == blas.vtable().sneg_tcopy.?);
}
