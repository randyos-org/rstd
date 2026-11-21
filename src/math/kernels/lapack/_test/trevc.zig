//! Tests of `../trevc.zig`.

const std = @import("std");
const ops = @import("../ops.zig");
const auxiliary = @import("../auxiliary.zig");
const trevc = @import("../trevc.zig");

const Float = ops.Float;
const trevc3 = trevc.trevc3;

const testing = std.testing;
const check = @import("../../blas/_test/check.zig");
const types = @import("../../blas/types.zig");
const schur = @import("../schur.zig");
const hessenberg = @import("../hessenberg.zig");

/// `‖A v - λ v‖` (right) or `‖vᴴ A - λ vᴴ‖` (left) for column `j` of `v`
/// (a real pair's two columns as one complex vector), and its largest
/// `|re| + |im|`.
fn eigResidual(comptime T: type, n: usize, a: []const T, v: []const T, j: usize, pair: bool, lam: std.math.Complex(Float(T)), right: bool) [2]Float(T) {
    const F = Float(T);
    const C = std.math.Complex(F);
    var worst: F = 0;
    var big: F = 0;
    for (0..n) |k| {
        const x: C = if (comptime types.isComplex(T)) v[k + j * n] else .init(v[k + j * n], if (pair) v[k + (j + 1) * n] else 0);
        big = @max(big, @abs(x.re) + @abs(x.im));
    }
    for (0..n) |i| {
        var s: C = .init(0, 0);
        for (0..n) |k| {
            const e = if (right) a[i + k * n] else a[k + i * n];
            const aik: C = if (comptime types.isComplex(T)) e else .init(e, 0);
            const xk_: C = if (comptime types.isComplex(T)) v[k + j * n] else .init(v[k + j * n], if (pair) v[k + (j + 1) * n] else 0);
            s = s.add(aik.mul(if (right) xk_ else xk_.conjugate()));
        }
        const xi_: C = if (comptime types.isComplex(T)) v[i + j * n] else .init(v[i + j * n], if (pair) v[i + (j + 1) * n] else 0);
        const r = s.sub(if (right) lam.mul(xi_) else lam.mul(xi_.conjugate()));
        worst = @max(worst, @sqrt(r.re * r.re + r.im * r.im));
    }
    return .{ worst, big };
}

test "trevc3: A's eigenvectors after gehrd, orghr and hseqr, back-transformed and selected, real and complex" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64, f32, types.Complex32 }) |T| {
        const F = Float(T);
        const C = std.math.Complex(F);
        for ([_]usize{ 1, 2, 5, 12, 31 }) |n| {
            const a0 = try al.alloc(T, n * n);
            defer al.free(a0);
            check.fill(T, a0, 13 * n + 1);
            const t = try al.dupe(T, a0);
            defer al.free(t);
            const tau = try al.alloc(T, @max(n, 2) - 1);
            defer al.free(tau);
            try hessenberg.gehrd(al, check.io, T, n, 0, n, t, n, tau);
            const z = try al.dupe(T, t);
            defer al.free(z);
            try hessenberg.orghr(al, check.io, T, n, 0, n, z, n, tau);
            const ev = try al.alloc(F, 2 * n);
            defer al.free(ev);
            const cw = try al.alloc(T, n);
            defer al.free(cw);
            const w: schur.Eigenvalues(T) = if (comptime types.isComplex(T)) cw else .{ .re = ev[0..n], .im = ev[n..] };
            try schur.hseqr(check.io, T, .s, .v, n, 0, n, t, n, w, z, n, null);
            const vl = try al.dupe(T, z);
            defer al.free(vl);
            const vr = try al.dupe(T, z);
            defer al.free(vr);
            var m = try trevc3(al, check.io, T, .b, .b, null, n, t, n, vl, n, vr, n, n);
            try testing.expectEqual(n, m);
            const tol = 256 * std.math.floatEps(F) * @as(F, @floatFromInt(n + 1)) * @max(auxiliary.normMax(T, .{ .ptr = a0.ptr, .rows = n, .cols = n, .ld = n }), 1);
            var j: usize = 0;
            while (j < n) : (j += 1) {
                const pair = !types.isComplex(T) and ev[n + j] > 0;
                const lam: C = if (comptime types.isComplex(T)) cw[j] else .init(ev[j], ev[n + j]);
                for ([_]bool{ true, false }) |right| {
                    const r = eigResidual(T, n, a0, if (right) vr else vl, j, pair, lam, right);
                    if (!(r[0] <= tol) or !(@abs(r[1] - 1) <= tol)) {
                        check.say("trevc3 {s} n={d} vector {d} {s}: residual {e}, largest {e}\n", .{ @typeName(T), n, j, if (right) "right" else "left", r[0], r[1] });
                        return error.TestUnexpectedResult;
                    }
                }
                if (pair) j += 1;
            }
            // selected: every other eigenvalue, a pair picked by its second entry
            const sel = try al.alloc(bool, n);
            defer al.free(sel);
            for (sel, 0..) |*s, i| s.* = i % 2 == 1;
            const all = try al.alloc(T, n * n);
            defer al.free(all);
            const some = try al.alloc(T, n * n);
            defer al.free(some);
            m = try trevc3(al, check.io, T, .r, .a, null, n, t, n, vl, n, all, n, n);
            try testing.expectEqual(n, m);
            m = try trevc3(al, check.io, T, .r, .s, sel, n, t, n, vl, n, some, n, n);
            var c: usize = 0;
            j = 0;
            while (j < n) : (j += 1) {
                const pair = !types.isComplex(T) and ev[n + j] > 0;
                if (sel[j]) {
                    for (0..if (pair) 2 else 1) |o| {
                        try testing.expectEqualSlices(T, all[(j + o) * n ..][0..n], some[(c + o) * n ..][0..n]);
                    }
                    c += if (pair) 2 else 1;
                }
                if (pair) {
                    try testing.expect(!sel[j + 1]);
                    j += 1;
                }
            }
            try testing.expectEqual(c, m);
        }
    }
    var a: [4]f64 = @splat(0);
    var sel = [_]bool{ true, false };
    try testing.expectError(error.MissingSelect, trevc3(al, check.io, f64, .r, .s, null, 2, &a, 2, &a, 1, &a, 2, 2));
    try testing.expectError(error.SliceTooShort, trevc3(al, check.io, f64, .r, .s, sel[0..1], 2, &a, 2, &a, 1, &a, 2, 2));
    try testing.expectError(error.LeadingDimensionTooSmall, trevc3(al, check.io, f64, .r, .a, null, 2, &a, 2, &a, 1, &a, 1, 2));
    try testing.expectError(error.LeadingDimensionTooSmall, trevc3(al, check.io, f64, .l, .a, null, 2, &a, 2, &a, 1, &a, 1, 2));
    try testing.expectError(error.SliceTooShort, trevc3(al, check.io, f64, .r, .a, null, 2, &a, 2, &a, 1, a[0..3], 2, 2));
    try testing.expectError(error.SliceTooShort, trevc3(al, check.io, f64, .r, .a, null, 2, a[0..3], 2, &a, 1, &a, 2, 2));
    try testing.expectError(error.TooFewColumns, trevc3(al, check.io, f64, .r, .a, null, 2, &a, 2, &a, 1, &a, 2, 1));
    try testing.expectEqual(1, try trevc3(al, check.io, f64, .r, .s, &sel, 2, &a, 2, &a, 1, &a, 2, 1));
}

test "trevc3: a zero Q back-transforms to zero columns, not NaN, real with a pair and complex" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64 }) |T| {
        const n = 3;
        // real: a standardized 2 x 2 pair, then a single; complex: upper triangular
        const t: [n * n]T = if (comptime types.isComplex(T)) .{
            .init(1, 0),   .init(0, 0),    .init(0, 0),
            .init(0.5, 1), .init(2, 1),    .init(0, 0),
            .init(1, 0),   .init(0.25, 0), .init(3, -1),
        } else .{ 1, -0.5, 0, 2, 1, 0, 0.5, 0.25, 3 };
        var tm = t;
        var vl: [n * n]T = @splat(ops.zero(T));
        var vr: [n * n]T = @splat(ops.zero(T));
        try testing.expectEqual(n, try trevc3(al, check.io, T, .b, .b, null, n, &tm, n, &vl, n, &vr, n, n));
        for (vl, vr) |l, r| try testing.expect(ops.abs1(T, l) == 0 and ops.abs1(T, r) == 0);
    }
}
