//! The kernels of `trsyl` (`trsyl.zig`) and `tgsyl` (`tgsyl.zig`):
//! Reference-LAPACK's `xTRSYL` (`trsy2`, here), `xLALN2` (`laln2.zig`),
//! `xLASY2` (`lasy2.zig`), and `xTGSY2` with `xGETC2`, `xGESC2` and `xLATDF`
//! (`tgsy2.zig`; `lapack-netlib/SRC/dtrsyl.f`, `dlaln2.f`, `dlasy2.f`,
//! `dtgsy2.f` and their forms, OpenBLAS v0.3.34-317). Complex division is
//! Smith's (`util.cdiv`) where the reference takes `ZLADIV` or `DLADIV`.

const std = @import("std");
const types = @import("../blas/types.zig");
const ops = @import("ops.zig");
const driver = @import("../blas/blas3/driver.zig");
const laln2_solve = @import("laln2.zig");
const lasy2_solve = @import("lasy2.zig");
const tgsy2_kernel = @import("tgsy2.zig");

const Float = types.Float;
const isComplex = types.isComplex;
const Trans = types.Trans;
const scalar = ops.scalar;
const add = ops.add;
const sub = ops.sub;
const mul = ops.mul;

const Error = driver.Error;

pub const Small = laln2_solve.Small;
const laln2Small = laln2_solve.laln2Small;
pub const Laln2 = laln2_solve.Laln2;
pub const laln2 = laln2_solve.laln2;
pub const lasy2 = lasy2_solve.lasy2;
const lasy2Square = lasy2_solve.lasy2Square;
pub const Dif = tgsy2_kernel.Dif;
pub const tgsy2 = tgsy2_kernel.tgsy2;

/// What `trsy2` returns: `scale` and whether a near-singular system was
/// perturbed (the reference's `INFO = 1`).
pub fn Solved(comptime F: type) type {
    return struct { scale: F, perturbed: bool };
}

/// `xLANGE('M')`: the largest `|a(i, j)|` of the `n x n` `a`, the modulus
/// for complex.
fn maxAbs(comptime T: type, n: usize, a: [*]const T, lda: usize) Float(T) {
    var m: Float(T) = 0;
    for (0..n) |j| for (0..n) |i| {
        const v = a[i + j * lda];
        const x = if (comptime isComplex(T)) std.math.hypot(v.re, v.im) else @abs(v);
        if (x > m or std.math.isNan(x)) m = x;
    };
    return m;
}

/// The diagonal blocks of the `n x n` upper quasi-triangular `a`: their
/// first indices, 1 x 1 unless a subdiagonal is nonzero (real only).
pub fn blockOf(comptime T: type, n: usize, a: [*]const T, lda: usize, i: usize) usize {
    if (comptime isComplex(T)) return 1;
    return if (i + 1 < n and a[i + 1 + i * lda] != 0) 2 else 1;
}

/// `xTRSYL`: `ta`, `tb` transpose (complex: conjugate-transpose) `A`, `B`.
pub fn trsy2(comptime T: type, ta: bool, tb: bool, isgn: isize, m: usize, n: usize, a: [*]const T, lda: usize, b: [*]const T, ldb: usize, c: [*]T, ldc: usize) Solved(Float(T)) {
    const F = Float(T);
    var out: Solved(F) = .{ .scale = 1, .perturbed = false };
    if (m == 0 or n == 0) return out;
    const eps = std.math.floatEps(F);
    const smlnum = std.math.floatMin(F) * @as(F, @floatFromInt(m * n)) / eps;
    const bignum = 1 / smlnum;
    const smin = @max(smlnum, eps * maxAbs(T, m, a, lda), eps * maxAbs(T, n, b, ldb));
    const sgn: F = @floatFromInt(isgn);
    // the blocks of `B` (columns of `C`) and `A` (rows), in order
    var l: usize = if (tb) n else 0;
    while (if (tb) l > 0 else l < n) {
        // backward from the last block for `tb`, forward otherwise
        const l1, const l2 = blk: {
            if (!tb) {
                const s = blockOf(T, n, b, ldb, l);
                break :blk .{ l, l + s - 1 };
            }
            const s: usize = if (comptime !isComplex(T)) (if (l >= 2 and b[l - 1 + (l - 2) * ldb] != 0) 2 else 1) else 1;
            break :blk .{ l - s, l - 1 };
        };
        var k: usize = if (ta) 0 else m;
        while (if (ta) k < m else k > 0) {
            const k1, const k2 = blk: {
                if (ta) {
                    const s = blockOf(T, m, a, lda, k);
                    break :blk .{ k, k + s - 1 };
                }
                const s: usize = if (comptime !isComplex(T)) (if (k >= 2 and a[k - 1 + (k - 2) * lda] != 0) 2 else 1) else 1;
                break :blk .{ k - s, k - 1 };
            };
            const ctx: Ctx(T) = .{ .ta = ta, .tb = tb, .m = m, .n = n, .a = a, .lda = lda, .b = b, .ldb = ldb, .c = c, .ldc = ldc, .k1 = k1, .k2 = k2, .l1 = l1, .l2 = l2, .sgn = sgn };
            var x: [4]T = undefined;
            var scaloc: F = 1;
            if (l1 == l2 and k1 == k2) {
                const vec = ctx.rhs(k1, l1);
                var a11 = ctx.diag(k1, l1);
                var da11 = ops.abs1(T, a11);
                if (da11 <= smin) {
                    a11 = scalar(T, smin);
                    da11 = smin;
                    out.perturbed = true;
                }
                const db = ops.abs1(T, vec);
                if (da11 < 1 and db > 1 and db > bignum * da11) scaloc = 1 / db;
                x[0] = ops.div(T, ops.mulReal(T, vec, scaloc), a11);
            } else if (comptime !isComplex(T)) {
                const r = if (l1 == l2) blk: {
                    // a 2 x 1 block of `A`
                    const v = [2]F{ ctx.rhs(k1, l1), ctx.rhs(k2, l1) };
                    break :blk laln2Small(F, ta, smin, a + k1 + k1 * lda, lda, v, -sgn * b[l1 + l1 * ldb]);
                } else if (k1 == k2) blk: {
                    // a 1 x 2 block of `B`, from the right
                    const v = [2]F{ sgn * ctx.rhs(k1, l1), sgn * ctx.rhs(k1, l2) };
                    break :blk laln2Small(F, !tb, smin, b + l1 + l1 * ldb, ldb, v, -sgn * a[k1 + k1 * lda]);
                } else blk: {
                    const v = [4]F{ ctx.rhs(k1, l1), ctx.rhs(k2, l1), ctx.rhs(k1, l2), ctx.rhs(k2, l2) };
                    break :blk lasy2Square(F, ta, tb, sgn, a + k1 + k1 * lda, lda, b + l1 + l1 * ldb, ldb, v);
                };
                if (r.perturbed) out.perturbed = true;
                scaloc = r.scale;
                x = r.x;
            }
            if (scaloc != 1) {
                for (0..n) |j| for (0..m) |i| {
                    c[i + j * ldc] = ops.mulReal(T, c[i + j * ldc], scaloc);
                };
                out.scale *= scaloc;
            }
            if (l1 == l2 and k1 == k2) {
                c[k1 + l1 * ldc] = x[0];
            } else if (l1 == l2) {
                c[k1 + l1 * ldc] = x[0];
                c[k2 + l1 * ldc] = x[1];
            } else if (k1 == k2) {
                c[k1 + l1 * ldc] = x[0];
                c[k1 + l2 * ldc] = x[1];
            } else {
                c[k1 + l1 * ldc] = x[0];
                c[k2 + l1 * ldc] = x[1];
                c[k1 + l2 * ldc] = x[2];
                c[k2 + l2 * ldc] = x[3];
            }
            k = if (ta) k2 + 1 else k1;
        }
        l = if (tb) l1 else l2 + 1;
    }
    return out;
}

/// One block's view of the problem.
fn Ctx(comptime T: type) type {
    return struct {
        ta: bool,
        tb: bool,
        m: usize,
        n: usize,
        a: [*]const T,
        lda: usize,
        b: [*]const T,
        ldb: usize,
        c: [*]T,
        ldc: usize,
        k1: usize,
        k2: usize,
        l1: usize,
        l2: usize,
        sgn: Float(T),

        const Self = @This();

        fn cj(v: T) T {
            return if (comptime isComplex(T)) v.conjugate() else v;
        }

        /// `C(i, j) - (suml + sgn sumr)`, the solved rows (`op(A)`) and
        /// columns (`op(B)`) summed as the reference's `?DOT`s.
        fn rhs(s: Self, i: usize, j: usize) T {
            var suml = scalar(T, 0);
            if (s.ta) {
                // `A(p, i)` conjugated, rows before the block
                for (0..s.k1) |p| suml = add(T, suml, mul(T, cj(s.a[p + i * s.lda]), s.c[p + j * s.ldc]));
            } else {
                for (s.k2 + 1..s.m) |p| suml = add(T, suml, mul(T, s.a[i + p * s.lda], s.c[p + j * s.ldc]));
            }
            var sumr = scalar(T, 0);
            if (s.tb) {
                // `ZDOTC(C(i, l2+1:), B(j, l2+1:))` conjugated back
                for (s.l2 + 1..s.n) |q| sumr = add(T, sumr, mul(T, cj(s.c[i + q * s.ldc]), s.b[j + q * s.ldb]));
                sumr = cj(sumr);
            } else {
                for (0..s.l1) |q| sumr = add(T, sumr, mul(T, s.c[i + q * s.ldc], s.b[q + j * s.ldb]));
            }
            return sub(T, s.c[i + j * s.ldc], add(T, suml, ops.mulReal(T, sumr, s.sgn)));
        }

        /// `op(A)(i, i) + sgn op(B)(j, j)`, as each case of `ZTRSYL` forms it.
        fn diag(s: Self, i: usize, j: usize) T {
            const aii = s.a[i + i * s.lda];
            const bjj = s.b[j + j * s.ldb];
            if (s.ta and s.tb) return cj(add(T, aii, ops.mulReal(T, bjj, s.sgn)));
            return add(T, if (s.ta) cj(aii) else aii, ops.mulReal(T, if (s.tb) cj(bjj) else bjj, s.sgn));
        }
    };
}

/// `?LASCL('G', 1, s)`: `C *= s`, part by part; `trsyl`'s and `tgsyl`'s.
pub fn scaleBlock(comptime T: type, m: usize, n: usize, s: Float(T), c: [*]T, ldc: usize) void {
    for (0..n) |j| for (0..m) |i| {
        c[i + j * ldc] = ops.mulReal(T, c[i + j * ldc], s);
    };
}
