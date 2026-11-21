//! Tests of `hseqr` and `gees` (`../schur.zig`).

const std = @import("std");
const types = @import("../../blas/types.zig");
const check = @import("../../blas/_test/check.zig");
const ops = @import("../ops.zig");
const auxiliary = @import("../auxiliary.zig");
const hessenberg = @import("../hessenberg.zig");
const schur = @import("../schur.zig");

const testing = std.testing;
const Float = ops.Float;
const isComplex = ops.isComplex;
const Eigenvalues = schur.Eigenvalues;
const GeesSelect = schur.GeesSelect;

test "hseqr: Zᴴ A Z = T after gehrd and orghr, the lo..hi block, both jobs, real and complex" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64, f32, types.Complex32 }) |T| {
        const F = Float(T);
        // n, lo, hi; an empty block leaves A triangular
        for ([_][3]usize{ .{ 1, 0, 1 }, .{ 2, 0, 2 }, .{ 6, 0, 6 }, .{ 9, 2, 7 }, .{ 30, 3, 27 }, .{ 5, 2, 2 } }) |c| {
            const n = c[0];
            const lo = c[1];
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
            try hessenberg.gehrd(al, check.io, T, n, lo, hi, t, n, tau);
            const z = try al.dupe(T, t);
            defer al.free(z);
            try hessenberg.orghr(al, check.io, T, n, lo, hi, z, n, tau);
            const h2 = try al.dupe(T, t);
            defer al.free(h2);
            const ev = try al.alloc(F, 4 * n);
            defer al.free(ev);
            const cw = try al.alloc(T, 2 * n);
            defer al.free(cw);
            const w: Eigenvalues(T) = if (comptime isComplex(T)) cw[0..n] else .{ .re = ev[0..n], .im = ev[n..][0..n] };
            const w2: Eigenvalues(T) = if (comptime isComplex(T)) cw[n..] else .{ .re = ev[2 * n ..][0..n], .im = ev[3 * n ..] };
            try schur.hseqr(check.io, T, .s, .v, n, lo, hi, t, n, w, z, n, null);
            try schur.hseqr(check.io, T, .e, .n, n, lo, hi, h2, n, w2, h2[0..0], 0, null);
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
                    check.say("hseqr {s} n={d} lo={d} hi={d}: ({d},{d}) T off by {e}, ZᴴZ by {e}\n", .{ @typeName(T), n, lo, hi, i, j, err, gerr });
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
    try testing.expectError(error.BlockOutOfRange, schur.hseqr(check.io, f64, .e, .n, 2, 2, 1, &a, 2, w, &a, 1, null));
    try testing.expectError(error.BlockOutOfRange, schur.hseqr(check.io, f64, .e, .n, 2, 0, 3, &a, 2, w, &a, 1, null));
    try testing.expectError(error.LeadingDimensionTooSmall, schur.hseqr(check.io, f64, .e, .n, 2, 0, 2, &a, 1, w, &a, 1, null));
    try testing.expectError(error.LeadingDimensionTooSmall, schur.hseqr(check.io, f64, .s, .i, 2, 0, 2, &a, 2, w, &a, 1, null));
    try testing.expectError(error.SliceTooShort, schur.hseqr(check.io, f64, .s, .i, 2, 0, 2, &a, 2, w, a[0..3], 2, null));
    try testing.expectError(error.SliceTooShort, schur.hseqr(check.io, f64, .e, .n, 2, 0, 2, &a, 2, .{ .re = wr[0..1], .im = &wi }, &a, 1, null));
    try schur.hseqr(check.io, f64, .e, .n, 0, 0, 0, &a, 1, w, &a, 1, null);
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
            const w: Eigenvalues(T) = if (comptime isComplex(T)) cw else .{ .re = ev[0..n], .im = ev[n..] };
            const sdim = try schur.gees(al, check.io, T, .v, if (sort) pick else null, n, t, n, w, z, n, null);
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
                try testing.expectEqual(0, sdim);
                continue;
            }
            var neg: usize = 0;
            for (0..n) |i| {
                const re = if (comptime isComplex(T)) cw[i].re else ev[i];
                if (re < 0) neg += 1;
                try testing.expectEqual(re < 0, i < sdim);
            }
            try testing.expectEqual(neg, sdim);
        };
    }
    var a: [4]f64 = .{ 1, std.math.nan(f64), 0, 1 };
    var z: [4]f64 = undefined;
    var wr: [2]f64 = undefined;
    var wi: [2]f64 = undefined;
    const w: Eigenvalues(f64) = .{ .re = &wr, .im = &wi };
    try testing.expectError(error.NanInA, schur.gees(al, check.io, f64, .v, null, 2, &a, 2, w, &z, 2, null));
    try testing.expectError(error.LeadingDimensionTooSmall, schur.gees(al, check.io, f64, .n, null, 2, &a, 1, w, &z, 1, null));
    try testing.expectError(error.LeadingDimensionTooSmall, schur.gees(al, check.io, f64, .v, null, 2, &a, 2, w, &z, 1, null));
    try testing.expectError(error.SliceTooShort, schur.gees(al, check.io, f64, .v, null, 2, &a, 2, w, z[0..3], 2, null));
    try testing.expectError(error.SliceTooShort, schur.gees(al, check.io, f64, .n, null, 2, &a, 2, .{ .re = &wr, .im = wi[0..1] }, z[0..0], 0, null));
    try testing.expectEqual(0, try schur.gees(al, check.io, f64, .n, null, 0, &a, 1, w, z[0..0], 0, null));
}
