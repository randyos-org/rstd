//! `GER` kernels: the rank-1 update `A += alpha*x*y^T` (conjugated for the
//! complex `c` form). The portable kernels and tests of every core's kernels
//! against plain arithmetic.

const blas = @import("../../__root__.zig");
const types = @import("../../types.zig");

const BlasLong = types.BlasLong;
const XDouble = types.XDouble;

// ---- portable kernels ------------------------------------------------------

const std = @import("std");
const util = @import("../../util.zig");
const Vtable = @import("../../Vtable.zig");
const portable = @import("../../portable.zig");

/// `kernel/generic/ger.c` and `zger.c`: the rank-1 update `A += alpha*x*y^T`,
/// `A` column-major `m` by `n`, as a copy of `x` where its stride is not 1
/// (into `buffer`, `m` elements) and one axpy of the core per column. `y` is
/// conjugated for `conj` (`gerc`). `T` is the element type, whose kernels it
/// calls. Upstream's `XCONJ` build (`gerv`, `gerd`), for row-major CBLAS
/// `zgerc` only, also conjugates `x` by taking `axpyc`.
fn Ger(comptime T: type, comptime conj: bool) type {
    const F = types.Float(T);
    const complex = comptime types.isComplex(T);
    const cs: BlasLong = if (complex) 2 else 1;
    return struct {
        fn run(m: BlasLong, n: BlasLong, alpha_r: F, alpha_i: F, x: [*]F, inc_x: BlasLong, y0: [*]F, inc_y: BlasLong, a0: [*]F, lda: BlasLong, buffer: [*]F) void {
            const axpy = portable.fetch(T, .axpy_k, blas.vtable);
            var xs = x;
            if (inc_x != 1) {
                xs = buffer;
                _ = portable.fetch(T, .copy_k, blas.vtable)(m, x, inc_x, xs, 1);
            }
            var a = a0;
            var y = y0;
            var j = n;
            while (j > 0) : (j -= 1) {
                // the kernels never read the dummies
                if (comptime complex) {
                    const br = y[0];
                    const bi = y[1];
                    const tr, const ti = if (conj)
                        .{ alpha_r * br + alpha_i * bi, -alpha_r * bi + alpha_i * br }
                    else
                        .{ alpha_r * br - alpha_i * bi, alpha_r * bi + alpha_i * br };
                    _ = axpy(m, 0, 0, tr, ti, xs, 1, a, 1, a, 0);
                } else {
                    _ = axpy(m, 0, 0, alpha_r * y[0], xs, 1, a, 1, a, 0);
                }
                a += @as(usize, @intCast(cs * lda));
                y = @ptrCast(util.at(F, y, cs * inc_y));
            }
        }
        fn real(m: BlasLong, n: BlasLong, _: BlasLong, alpha: F, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, a: [*]F, lda: BlasLong, buffer: [*]F) callconv(types.kernelCc(F)) i32 {
            run(m, n, alpha, 0, x, inc_x, y, inc_y, a, lda, buffer);
            return 0;
        }
        fn cplx(m: BlasLong, n: BlasLong, _: BlasLong, alpha_r: F, alpha_i: F, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, a: [*]F, lda: BlasLong, buffer: [*]F) callconv(types.kernelCc(F)) i32 {
            run(m, n, alpha_r, alpha_i, x, inc_x, y, inc_y, a, lda, buffer);
            return 0;
        }
        const k = if (complex) cplx else real;
    };
}

/// The portable kernels over `T`, any float or its `Complex`, by the slot
/// they would fill (`portable.kernel`).
pub fn Portable(comptime T: type) type {
    return struct {
        pub const ger_k = &Ger(T, false).k;
        pub const geru_k = &Ger(T, false).k;
        pub const gerc_k = &Ger(T, true).k;
    };
}

/// The complex variants by suffix, and whether `y` is conjugated.
const variants = .{ .{ "u", false }, .{ "c", true } };

/// Points every ger slot at the portable kernels.
pub fn installGeneric(vt: *Vtable) void {
    inline for (.{ .{ "s", f32 }, .{ "d", f64 }, .{ "q", XDouble } }) |p| {
        @field(vt, p[0] ++ "ger_k") = &Ger(p[1], false).k;
    }
    inline for (.{ .{ "c", f32 }, .{ "z", f64 }, .{ "x", XDouble } }) |p| {
        inline for (variants) |v| @field(vt, p[0] ++ "ger" ++ v[0] ++ "_k") = &Ger(std.math.Complex(p[1]), v[1]).k;
    }
}

// ---- tests -----------------------------------------------------------------

const check = @import("../../_test/blas1/check.zig");
const guard = @import("../../_test/blas2/guard.zig");

const C = std.math.Complex;
/// Shapes across the kernels' blocks, then bigger ones (static storage:
/// `guard.operand`).
const sizes = [_][2]usize{
    .{ 1, 1 },   .{ 2, 3 },    .{ 3, 2 },    .{ 4, 4 },   .{ 7, 5 },   .{ 8, 8 }, .{ 9, 16 }, .{ 16, 9 }, .{ 17, 17 }, .{ 33, 4 }, .{ 4, 33 },
    .{ 65, 65 }, .{ 100, 37 }, .{ 37, 100 }, .{ 129, 3 }, .{ 3, 129 },
};
const strides = [_][2]usize{ .{ 1, 1 }, .{ 2, 3 }, .{ 1, 2 }, .{ 3, 1 } };

/// `kern` against plain arithmetic; `T` is the element type and the flags
/// are the variant's. The buffer is the contract's: `m` elements for a
/// strided `x`, nothing for a unit one, 32-byte aligned as a stack one is,
/// and nothing written outside it. `x` and `y` are left alone.
fn checkGer(comptime T: type, comptime conj: bool, kern: anytype) !void {
    const F = check.FloatOf(T);
    const complex = comptime check.isComplex(T);
    const cs: usize = if (complex) 2 else 1;
    const alpha: T = if (complex) .init(1.25, -0.5) else 1.25;
    for (sizes) |mn| for (strides) |inc| {
        const m = mn[0];
        const n = mn[1];
        const lda = m + 1;
        const a = guard.operand(T, .a, lda * n);
        const x = guard.operand(T, .x, guard.span(m, inc[0]));
        const y = guard.operand(T, .y, guard.span(n, inc[1]));
        const want = guard.operand(T, .want, a.len);
        check.fill(T, a, m * 7 + n);
        check.fill(T, x, m + 3);
        check.fill(T, y, n + 5);
        @memcpy(want, a);
        // `x` and `y` as they were, for the comparison afterwards
        const x0 = guard.operand(T, .x0, x.len);
        const y0 = guard.operand(T, .y0, y.len);
        @memcpy(x0, x);
        @memcpy(y0, y);
        for (0..n) |j| for (0..m) |i| {
            const xv = x[i * inc[0]];
            const yv = y[j * inc[1]];
            if (complex) {
                const ye = if (conj) yv.conjugate() else yv;
                want[i + j * lda] = want[i + j * lda].add(alpha.mul(xv).mul(ye));
            } else want[i + j * lda] += alpha * xv * yv;
        };
        const bytes = if (inc[0] == 1) 0 else cs * m * @sizeOf(F);
        const buf = guard.open(F, bytes, 32, .nan);
        if (complex) {
            _ = kern(@intCast(m), @intCast(n), 0, alpha.re, alpha.im, @ptrCast(x.ptr), @intCast(inc[0]), @ptrCast(y.ptr), @intCast(inc[1]), @ptrCast(a.ptr), @intCast(lda), buf);
        } else {
            _ = kern(@intCast(m), @intCast(n), 0, alpha, x.ptr, @intCast(inc[0]), y.ptr, @intCast(inc[1]), a.ptr, @intCast(lda), buf);
        }
        check.expectCloseSlices(T, want, a, check.tolerance(F, 4)) catch |e| {
            check.say("m={d} n={d}: ", .{ m, n });
            return check.case(e, "ger", T, m, inc);
        };
        guard.intact(bytes, 32) catch |e| {
            check.say("m={d} n={d}: ", .{ m, n });
            return check.case(e, "ger buffer", T, m, inc);
        };
        if (!std.mem.eql(u8, std.mem.sliceAsBytes(x0), std.mem.sliceAsBytes(x)) or !std.mem.eql(u8, std.mem.sliceAsBytes(y0), std.mem.sliceAsBytes(y))) {
            check.say("m={d} n={d}: ", .{ m, n });
            return check.case(error.TestUnexpectedResult, "ger changed x or y", T, m, inc);
        }
    };
}

fn checkAll(vt: *const Vtable) !void {
    guard.use(vt);
    defer guard.restore();
    try checkSlots(vt);
}

fn checkSlots(vt: *const Vtable) !void {
    try checkGer(f32, false, vt.sger_k.?);
    try checkGer(f64, false, vt.dger_k.?);
    inline for (.{ .{ "c", f32 }, .{ "z", f64 } }) |p| {
        inline for (variants) |v| try checkGer(C(p[1]), v[1], @field(vt, p[0] ++ "ger" ++ v[0] ++ "_k").?);
    }
}

test "portable ger kernels match plain arithmetic" {
    const vt = guard.portableVtable();
    guard.use(vt);
    defer guard.restore();
    try checkSlots(vt);
    try checkGer(XDouble, false, vt.qger_k.?);
    inline for (variants) |v| try checkGer(C(XDouble), v[1], @field(vt, "xger" ++ v[0] ++ "_k").?);
}

test "every linked core's ger kernels match plain arithmetic, where this machine runs them" {
    try check.eachCore(checkAll);
}
