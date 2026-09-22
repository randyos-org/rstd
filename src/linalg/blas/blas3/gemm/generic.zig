//! Portable `gemm` building blocks: the beta scaling, the packing copies and the
//! packed-panel kernel of OpenBLAS `kernel/generic/`
//! (`gemm_beta.c`, `gemm_{n,t}copy_<U>.c`, `gemmkernel_{2x2,4x4}.c`), for real
//! elements.
//!
//! Upstream has one file per unroll size; the layout they write is one rule.
//! The packed operand is a run of panels: `n / U` of width `U`, then one of each
//! power of two below `U` that is set in the remainder. Panel `p`, starting at
//! column `c0` with width `w`, is stored at `m * c0` as `m` groups of `w` values
//! (row-major within the panel for the `t` copies, taken down each column for the
//! `n` copies). `gemm_{n,t}copy_6.c` packs as unroll 4 (checked against the C for
//! every `m`, `n` up to 39). The kernel reads the same layout, so a kernel of unroll `UM x UN`
//! needs its operands packed with `UM` and `UN`.
//!
//! Every function here is called directly (`generic.ncopy(f64, 4, m, n, a, lda,
//! b)`), not returned as a value: nothing currently needs `ncopy(f64, 4)` on its
//! own as a first-class `fn`. A helper that a function passes to `panels` as a
//! value (`one`, `tile`, `column`, `fma_`) still goes through
//! `struct { fn call(...) {...} }.call`, because Zig has no nested named `fn`
//! at block scope; that indirection is not optional there.

const std = @import("std");
const types = @import("../../types.zig");
const fma = @import("../../fma.zig");

const BlasLong = types.BlasLong;

/// Calls `visit(ctx, start, width)` for each panel that `n` columns pack into
/// for unroll `u`, in storage order.
fn panels(comptime u: comptime_int, n: usize, comptime visit: anytype, ctx: anytype) void {
    var c0: usize = 0;
    while (c0 + u <= n) : (c0 += u) visit(ctx, c0, u);
    comptime var w = u / 2;
    inline while (w >= 1) : (w /= 2) {
        if ((n % u) & w != 0) {
            visit(ctx, c0, w);
            c0 += w;
        }
    }
}

/// `gemm_beta.c`: `C := beta * C`, or zero where `beta` is zero (`C` is not
/// read, so it may hold anything).
pub fn beta(comptime T: type, m: BlasLong, n: BlasLong, b: T, c: [*]T, ldc: BlasLong) void {
    const rows: usize = @intCast(m);
    const cols: usize = @intCast(n);
    const ld: usize = @intCast(ldc);
    for (0..cols) |j| {
        const col = c + j * ld;
        if (b == 0) {
            @memset(col[0..rows], 0);
        } else {
            for (col[0..rows]) |*v| v.* *= b;
        }
    }
}

/// `gemm_ncopy_<U>.c`: packs `n` columns of the column-major `m x n` matrix
/// `a` in panels of `u`, taking each panel down its columns.
pub fn ncopy(comptime T: type, comptime u: comptime_int, m: BlasLong, n: BlasLong, a: [*]const T, lda: BlasLong, b: [*]T) void {
    const Ctx = struct { m: usize, a: [*]const T, lda: usize, b: [*]T };
    const one = struct {
        fn call(ctx: Ctx, c0: usize, comptime w: comptime_int) void {
            const dst = ctx.b + ctx.m * c0;
            for (0..ctx.m) |i| {
                inline for (0..w) |c| dst[i * w + c] = ctx.a[(c0 + c) * ctx.lda + i];
            }
        }
    }.call;
    const ctx: Ctx = .{ .m = @intCast(m), .a = a, .lda = @intCast(lda), .b = b };
    panels(u, @intCast(n), one, ctx);
}

/// `gemm_tcopy_<U>.c`: packs the columns of the `m x n` matrix `a` whose rows
/// are `lda` apart, in panels of `u`, taking each panel along its rows.
pub fn tcopy(comptime T: type, comptime u: comptime_int, m: BlasLong, n: BlasLong, a: [*]const T, lda: BlasLong, b: [*]T) void {
    const Ctx = struct { m: usize, a: [*]const T, lda: usize, b: [*]T };
    const one = struct {
        fn call(ctx: Ctx, c0: usize, comptime w: comptime_int) void {
            const dst = ctx.b + ctx.m * c0;
            for (0..ctx.m) |i| {
                inline for (0..w) |c| dst[i * w + c] = ctx.a[i * ctx.lda + c0 + c];
            }
        }
    }.call;
    const ctx: Ctx = .{ .m = @intCast(m), .a = a, .lda = @intCast(lda), .b = b };
    panels(u, @intCast(n), one, ctx);
}

/// `gemmkernel_<UM>x<UN>.c`: `C += alpha * A * B` for the packed panels `ba`
/// (`bm x bk`, packed with `um`) and `bb` (`bk x bn`, packed with `un`).
pub fn kernel(comptime T: type, comptime um: comptime_int, comptime un: comptime_int, bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha: T, ba: [*]const T, bb: [*]const T, c: [*]T, ldc: BlasLong) void {
    const Ctx = struct { bm: usize, bk: usize, alpha: T, ba: [*]const T, bb: [*]const T, c: [*]T, ldc: usize };
    const JCtx = struct { ctx: Ctx, q0: usize };

    const fma_ = struct {
        inline fn call(acc: T, a: T, b: T) T {
            return if (fma.fuse) @mulAdd(T, a, b, acc) else acc + a * b;
        }
    }.call;

    const tile = struct {
        fn call(jc: JCtx, r0: usize, comptime wm: comptime_int, comptime wn: comptime_int) void {
            const ctx = jc.ctx;
            const pa = ctx.ba + r0 * ctx.bk;
            const pb = ctx.bb + jc.q0 * ctx.bk;
            var res: [wm][wn]T = undefined;
            inline for (0..wm) |r| inline for (0..wn) |q| {
                res[r][q] = 0;
            };
            for (0..ctx.bk) |k| {
                inline for (0..wn) |q| inline for (0..wm) |r| {
                    res[r][q] = fma_(res[r][q], pa[k * wm + r], pb[k * wn + q]);
                };
            }
            inline for (0..wn) |q| inline for (0..wm) |r| {
                const p = &ctx.c[(r0 + r) + (jc.q0 + q) * ctx.ldc];
                p.* = if (fma.fuse) @mulAdd(T, res[r][q], ctx.alpha, p.*) else p.* + res[r][q] * ctx.alpha;
            };
        }
    }.call;

    const column = struct {
        fn call(ctx: Ctx, q0: usize, comptime wn: comptime_int) void {
            const jc: JCtx = .{ .ctx = ctx, .q0 = q0 };
            var r0: usize = 0;
            while (r0 + um <= ctx.bm) : (r0 += um) tile(jc, r0, um, wn);
            comptime var w = um / 2;
            inline while (w >= 1) : (w /= 2) {
                if ((ctx.bm % um) & w != 0) {
                    tile(jc, r0, w, wn);
                    r0 += w;
                }
            }
        }
    }.call;

    const ctx: Ctx = .{ .bm = @intCast(bm), .bk = @intCast(bk), .alpha = alpha, .ba = ba, .bb = bb, .c = c, .ldc = @intCast(ldc) };
    const n: usize = @intCast(bn);
    var q0: usize = 0;
    while (q0 + un <= n) : (q0 += un) column(ctx, q0, un);
    comptime var w = un / 2;
    inline while (w >= 1) : (w /= 2) {
        if ((n % un) & w != 0) {
            column(ctx, q0, w);
            q0 += w;
        }
    }
}

// ---- tests -----------------------------------------------------------------

fn fillSeq(comptime T: type, buf: []T, seed: usize) void {
    for (buf, 0..) |*v, i| v.* = @floatFromInt(@as(i64, @intCast((i * 7 + seed * 13) % 23)) - 11);
}

/// `C += alpha * A * B` on column-major matrices, the plain way.
fn refGemm(comptime T: type, m: usize, n: usize, k: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, c: []T, ldc: usize) void {
    for (0..n) |j| for (0..m) |i| {
        var s: T = 0;
        for (0..k) |l| s += a[i + l * lda] * b[l + j * ldb];
        c[i + j * ldc] += alpha * s;
    };
}

fn checkKernel(comptime T: type, comptime um: comptime_int, comptime un: comptime_int) !void {
    const sizes = [_]usize{ 1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17 };
    for (sizes) |m| for (sizes) |n| for ([_]usize{ 1, 3, 8, 13 }) |k| {
        var a: [20 * 16]T = undefined;
        var b: [20 * 20]T = undefined;
        var c: [20 * 20]T = undefined;
        var want: [20 * 20]T = undefined;
        fillSeq(T, &a, 1);
        fillSeq(T, &b, 2);
        fillSeq(T, &c, 3);
        want = c;
        refGemm(T, m, n, k, 0.5, &a, 20, &b, 20, &want, 20);

        var pa: [20 * 16]T = undefined;
        var pb: [20 * 20]T = undefined;
        // A (m x k, column-major) goes in panels of its rows, which is `tcopy` of
        // the k x m view with the same stride; B (k x n) in panels of its columns.
        tcopy(T, um, @intCast(k), @intCast(m), &a, 20, &pa);
        ncopy(T, un, @intCast(k), @intCast(n), &b, 20, &pb);
        kernel(T, um, un, @intCast(m), @intCast(n), @intCast(k), 0.5, &pa, &pb, &c, 20);
        for (0..n) |j| for (0..m) |i| try std.testing.expectApproxEqAbs(want[i + j * 20], c[i + j * 20], 1e-9);
    };
}

test "packed kernels match a plain gemm for every panel remainder" {
    try checkKernel(f64, 2, 2);
    try checkKernel(f64, 4, 4);
    try checkKernel(f64, 8, 4);
    try checkKernel(f32, 4, 2);
}

test "beta scales the m x n block only, and zero does not read c" {
    var c = [_]f64{ std.math.nan(f64), 2, 3, 4, 5, 6, 7, 8 };
    beta(f64, 2, 2, 0, &c, 4);
    try std.testing.expectEqual(@as(f64, 0), c[0]);
    try std.testing.expectEqual(@as(f64, 0), c[1]);
    try std.testing.expectEqual(@as(f64, 3), c[2]);
    try std.testing.expectEqual(@as(f64, 0), c[4]);
    c[4] = 5;
    beta(f64, 2, 2, 2, &c, 4);
    try std.testing.expectEqual(@as(f64, 10), c[4]);
    try std.testing.expectEqual(@as(f64, 3), c[2]);
    try std.testing.expectEqual(@as(f64, 14), c[6] * 2);
}
