//! Tests of the `trsm` kernels.

const blas = @import("../../__root__.zig");
const apple_accel = @import("../../apple_accel.zig");
const types = @import("../../types.zig");

const Complex32 = types.Complex32;
const Complex64 = types.Complex64;

const std = @import("std");
const check = @import("../../_test/blas3/check.zig");
const Vtable = @import("../../Vtable.zig");
const driver = @import("../driver.zig");
const iface = @import("../interface.zig");
const triangular = @import("../../blas2/triangular.zig");
const corespec = @import("../../cpu/_corespec.zig");
const registry = @import("../../asm_registry.zig");
const util = @import("../../util.zig");

fn cases(vt: *const Vtable) !void {
    try check.routineAll(vt, .trsm);
}

const blas1check = @import("../../_test/blas1/check.zig");
const portable = @import("../../portable.zig");

test "portable trsm matches the wide reference through the driver" {
    var vt: Vtable = .{};
    portable.install(&vt);
    try cases(&vt);
}

test "every linked core's trsm matches the wide reference, where this machine runs it" {
    try blas1check.eachCore(cases);
}

const testing = std.testing;
const precisions = .{ f32, f64, Complex32, Complex64 };

/// A small integer, complex or not, as `[re, im]`.
const Int = [2]i64;

fn imul(x: Int, y: Int) Int {
    return .{ x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0] };
}

/// `t` times `v`, exact for a small `v`.
fn scaled(comptime T: type, v: Int, t: types.Float(T)) T {
    const re = @as(types.Float(T), @floatFromInt(v[0])) * t;
    const im = @as(types.Float(T), @floatFromInt(v[1])) * t;
    return if (comptime types.isComplex(T)) .init(re, im) else re;
}

/// `trsm` where `1 / A(i, i)` overflows: `A` is small integers times `t`, a
/// denormal whose reciprocal is infinite, the solution `X` small integers,
/// and `B` is `op(A) X` or `X op(A)`, so every step is exact and the wide
/// reference is `X` itself. Every side, triangle and transpose.
fn tinyDiagonal(comptime T: type, vt: *const Vtable) !void {
    const complex = comptime types.isComplex(T);
    const F = types.Float(T);
    const t = std.math.ldexp(@as(F, 1), std.math.floatExponentMin(F) - 4);
    try testing.expect(std.math.isInf(1 / t));
    const m = 9;
    const n = 7;
    const nan = std.math.nan(F);
    const ts: []const check.Trans = if (complex) &.{ .n, .t, .r, .c } else &.{ .n, .t };
    for ([_]bool{ true, false }) |left| for ([_]bool{ true, false }) |upper| for (ts) |tr| {
        const ka: usize = if (left) m else n;
        // `A` as integers, then stored scaled, NaN where it is not stored
        var ai: [m * m]Int = undefined;
        var a: [m * m]T = undefined;
        for (0..ka) |c| for (0..ka) |r| {
            const v: Int = if (r == c)
                (if (r % 2 == 0) (if (complex) .{ 1, 1 } else .{ 1, 0 }) else (if (complex) .{ 2, -1 } else .{ -2, 0 }))
            else
                .{ @as(i64, @intCast((r + 2 * c) % 3)) - 1, if (complex) @intCast((r + c) % 2) else 0 };
            const kept = if (upper) r <= c else r >= c;
            ai[r + c * ka] = if (kept) v else .{ 0, 0 };
            a[r + c * ka] = if (kept) scaled(T, v, t) else if (complex) .init(nan, nan) else nan;
        };
        const opA = struct {
            fn f(p: []const Int, k: usize, trc: check.Trans, i: usize, l: usize) Int {
                const v = if (trc.transposed()) p[l + i * k] else p[i + l * k];
                return if (trc.conjugated()) .{ v[0], -v[1] } else v;
            }
        }.f;
        var x: [m * n]Int = undefined;
        for (0..n) |j| for (0..m) |i| {
            x[i + j * m] = .{ @as(i64, @intCast((3 * i + j) % 5)) - 2, if (complex) @as(i64, @intCast((i + 2 * j) % 3)) - 1 else 0 };
        };
        var b: [m * n]T = undefined;
        for (0..n) |j| for (0..m) |i| {
            var sum: Int = .{ 0, 0 };
            for (0..ka) |l| {
                const p = if (left) imul(opA(&ai, ka, tr, i, l), x[l + j * m]) else imul(x[i + l * m], opA(&ai, ka, tr, l, j));
                sum = .{ sum[0] + p[0], sum[1] + p[1] };
            }
            b[i + j * m] = scaled(T, sum, t);
        };
        const one: T = if (complex) .init(1, 0) else 1;
        try driver.trsm(testing.allocator, check.io, T, vt, left, upper, tr.transposed(), tr.conjugated(), false, m, n, one, &a, ka, &b, m, 1);
        for (0..n) |j| for (0..m) |i| {
            const want = scaled(T, x[i + j * m], 1);
            const got = b[i + j * m];
            const ok = if (complex) got.re == want.re and got.im == want.im else got == want;
            if (!ok) {
                std.debug.print("{s} left={} upper={} trans={t} at ({d},{d}): got {any}, want {any}\n", .{ @typeName(T), left, upper, tr, i, j, got, want });
                return error.TestUnexpectedResult;
            }
        };
    };
}

fn tinyCases(vt: *const Vtable) !void {
    inline for (precisions) |T| try tinyDiagonal(T, vt);
}

test "trsm divides by a diagonal whose reciprocal overflows, and stays finite and exact" {
    var vt: Vtable = .{};
    portable.install(&vt);
    try tinyCases(&vt);
}

test "every linked core's trsm divides by a tiny diagonal, where this machine runs it" {
    try blas1check.eachCore(tinyCases);
}

/// `trsv` and a one-column `trsm` on the same system, every option, through
/// the interfaces: the same bits. `A` is diagonal, its other stored elements
/// zero, so only the divisions round (the updates' orders differ); every
/// third diagonal element is so small its reciprocal overflows.
fn sameAsTrsv(comptime T: type) !void {
    const complex = comptime types.isComplex(T);
    const F = types.Float(T);
    const n = 37;
    const t = std.math.ldexp(@as(F, 1), std.math.floatExponentMin(F) - 4);
    const zero: T = if (complex) .init(0, 0) else 0;
    const one: T = if (complex) .init(1, 0) else 1;
    const nan: T = if (complex) .init(std.math.nan(F), std.math.nan(F)) else std.math.nan(F);
    var d: [n]T = undefined;
    var rhs: [n]T = undefined;
    blas1check.fill(T, &d, 5);
    blas1check.fill(T, &rhs, 6);
    for (0..n) |i| {
        // a third tiny, all at least 1 in magnitude before scaling
        const s: F = if (i % 3 == 0) t else 1;
        if (complex) {
            d[i] = .init((d[i].re + std.math.copysign(@as(F, 1), d[i].re)) * s, d[i].im * s);
            rhs[i] = .init(rhs[i].re * s, rhs[i].im * s);
        } else {
            d[i] = (d[i] + std.math.copysign(@as(F, 1), d[i])) * s;
            rhs[i] *= s;
        }
    }
    const ts: []const types.Trans = if (complex) &.{ .n, .t, .r, .c } else &.{ .n, .t };
    var a: [n * n]T = undefined;
    for ([_]bool{ true, false }) |upper| for (ts) |tr| for ([_]bool{ false, true }) |unit| {
        for (0..n) |c| for (0..n) |r| {
            a[r + c * n] = if (r == c) d[r] else if (if (upper) r < c else r > c) zero else nan;
        };
        var x = rhs;
        var y = rhs;
        const uplo: types.Uplo = if (upper) .u else .l;
        const diag: types.Diag = if (unit) .u else .n;
        try triangular.trsv(testing.allocator, check.io, T, uplo, tr, diag, n, &a, n, &x, 1);
        try iface.trsm(testing.allocator, check.io, T, .l, uplo, tr, diag, n, 1, one, &a, n, &y, n);
        for (x) |u| {
            const fin = if (complex) std.math.isFinite(u.re) and std.math.isFinite(u.im) else std.math.isFinite(u);
            try testing.expect(fin);
        }
        try testing.expectEqualSlices(u8, std.mem.asBytes(&x), std.mem.asBytes(&y));
    };
}

/// Whether the selected core's `trsm` computes `F` where `trsv` does. An
/// x87 build's `trsv` divides on the x87, and the objects of its cores with
/// SSE2 in SSE2 (`corespec.objectQuery`); every feature an x86 spec lists
/// implies SSE2. An unlinked core runs the target's own object. Where both
/// are on the x87, `util.cdiv` rounds each step, so the bits must match.
fn sameUnit(comptime F: type) bool {
    if (comptime !util.onX87(F)) return true;
    const core = blas.coreType() orelse return true;
    if (!registry.isLinked(core)) return true;
    const spec = corespec.of(core);
    return spec.model == null and spec.features.len == 0;
}

test "trsv and trsm give the same bits on the same one-column system" {
    // rstd's own: Accelerate takes `trsv` here but not `trsm` (its tiny
    // diagonal), and the two round otherwise; after init, which routes
    _ = blas.vtable();
    const was = apple_accel.routing();
    _ = apple_accel.route(false);
    defer _ = apple_accel.route(was);
    var skipped = false;
    inline for (precisions) |T| {
        if (sameUnit(types.Float(T))) try sameAsTrsv(T) else {
            std.debug.print("{s} skipped: the core's trsm divides in SSE2, trsv on the x87\n", .{@typeName(T)});
            skipped = true;
        }
    }
    if (skipped) return error.SkipZigTest;
}

/// `trsm` with a unit diagonal, which skips the division, against the same
/// solve dividing by a stored diagonal of ones: the same bits, every side,
/// triangle and transpose. The unit run's diagonal is NaN, so never read.
fn unitCases(vt: *const Vtable) !void {
    inline for (precisions) |T| {
        const complex = comptime types.isComplex(T);
        const F = types.Float(T);
        const m = 23;
        const n = 17;
        const nan: T = if (complex) .init(std.math.nan(F), std.math.nan(F)) else std.math.nan(F);
        const one: T = if (complex) .init(1, 0) else 1;
        var a: [m * m]T = undefined;
        var rhs: [m * n]T = undefined;
        blas1check.fill(T, &a, 7);
        blas1check.fill(T, &rhs, 8);
        const ts: []const check.Trans = if (complex) &.{ .n, .t, .r, .c } else &.{ .n, .t };
        for ([_]bool{ true, false }) |left| for ([_]bool{ true, false }) |upper| for (ts) |tr| {
            const ka: usize = if (left) m else n;
            var au = a;
            var an = a;
            for (0..ka) |i| {
                au[i + i * ka] = nan;
                an[i + i * ka] = one;
            }
            var bu = rhs;
            var bn = rhs;
            try driver.trsm(testing.allocator, check.io, T, vt, left, upper, tr.transposed(), tr.conjugated(), true, m, n, one, &au, ka, &bu, m, 1);
            try driver.trsm(testing.allocator, check.io, T, vt, left, upper, tr.transposed(), tr.conjugated(), false, m, n, one, &an, ka, &bn, m, 1);
            testing.expectEqualSlices(u8, std.mem.asBytes(&bn), std.mem.asBytes(&bu)) catch |e| {
                std.debug.print("{s} left={} upper={} trans={t}\n", .{ @typeName(T), left, upper, tr });
                return e;
            };
        };
    }
}

test "trsm on a unit diagonal gives the bits of dividing by ones" {
    var vt: Vtable = .{};
    portable.install(&vt);
    try unitCases(&vt);
}

test "every linked core's trsm on a unit diagonal gives the bits of dividing by ones" {
    try blas1check.eachCore(unitCases);
}
