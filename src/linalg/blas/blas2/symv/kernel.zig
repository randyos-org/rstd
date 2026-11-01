//! The drivers of OpenBLAS's per-core real `symv` kernels, around their 4x4
//! blocks written on `@Vector`s: x86_64's `ssymv_L.c`, `ssymv_U.c`,
//! `dsymv_L.c` and `dsymv_U.c` with their `*_microk_*-2.c` blocks, and
//! arm64's `symv_L_asimd_4x4.c` and `symv_U_asimd_4x4.c` (`symv_microk_asimd_4x4.c`).
//! A block's `Shape` says how its microkernel computes; the drivers' C
//! contracts as the core's flags say (`fma.fuse`). The arches' `kernel.zig`
//! name the shapes per core. Nothing here depends on the BLAS root, so the
//! per-core objects (`_kernel.zig`) can use it. Tests live in `__root__.zig`.

const std = @import("std");
const types = @import("../../types.zig");
const fma = @import("../../fma.zig");
const util = @import("../../util.zig");

const BlasLong = types.BlasLong;
const Slot = util.Slot;

/// How a microkernel computes its four columns.
pub const Shape = struct {
    /// Elements a vector.
    lanes: comptime_int,
    /// Whether its multiply-adds fuse.
    fused: bool,
    /// `temp2 = sum` rather than `temp2 += sum`: the U kernels that store
    /// their sums (`movsd`), and `ssymv_U`'s Sandy Bridge and Haswell ones,
    /// whose `vaddsd` of two floats as a double leaves the low one, the sum,
    /// as it is (their `temp2` is 0 on entry).
    store: bool = false,
    /// `y += t0*a0 + t1*a1 + t2*a2 + t3*a3` as one contracted expression
    /// (intrinsics), not four multiply-adds into `y`.
    sum_y: bool = false,
    /// A pass at twice the lanes first, its sums folded in halves into the
    /// narrow ones (`dsymv_L_microk_skylakex-2.c`'s 512 bits).
    wide: bool = false,
    /// Half the lanes for the whole range when its length has bit
    /// `lanes / 2` set (`ssymv_L_microk_sandy-2.c`).
    narrow_odd: bool = false,
};

/// `symv_kernel_4x4(from, to, ...)` of shape `s`: rows `from..to` of four
/// columns `a` into `y` (times `t`) and their dots with `x` into `t2`.
pub fn Block(comptime F: type, comptime s: Shape) type {
    return struct {
        fn pass(comptime L: comptime_int, i: usize, a: *const [4][*]F, x: [*]F, y: [*]F, t: [4]F, acc: *[4]@Vector(L, F)) void {
            const V = @Vector(L, F);
            const xv = util.load(V, x + i);
            var yv = util.load(V, y + i);
            var av: [4]V = undefined;
            inline for (0..4) |c| av[c] = util.load(V, a[c] + i);
            if (s.sum_y) {
                var sum = fma.madd(s.fused, V, @splat(t[0]), av[0], @as(V, @splat(t[1])) * av[1]);
                inline for (2..4) |c| sum = fma.madd(s.fused, V, @splat(t[c]), av[c], sum);
                yv += sum;
            } else inline for (0..4) |c| {
                yv = fma.madd(s.fused, V, @splat(t[c]), av[c], yv);
            }
            inline for (0..4) |c| acc[c] = fma.madd(s.fused, V, xv, av[c], acc[c]);
            util.store(V, y + i, yv);
        }

        fn sweep(comptime L: comptime_int, from: usize, to: usize, a: *const [4][*]F, x: [*]F, y: [*]F, t: [4]F, acc: *[4]@Vector(L, F)) void {
            var i = from;
            while (i + L <= to) : (i += L) pass(L, i, a, x, y, t, acc);
        }

        fn finish(comptime L: comptime_int, acc: [4]@Vector(L, F), t2: *[4]F) void {
            inline for (0..4) |c| {
                const r = util.hsum(acc[c]);
                t2[c] = if (s.store) r else r + t2[c];
            }
        }

        pub fn run(from: usize, to: usize, a: *const [4][*]F, x: [*]F, y: [*]F, t: [4]F, t2: *[4]F) void {
            const L = s.lanes;
            if (s.narrow_odd and (to - from) & (L / 2) != 0) {
                var acc: [4]@Vector(L / 2, F) = @splat(@splat(0));
                sweep(L / 2, from, to, a, x, y, t, &acc);
                return finish(L / 2, acc, t2);
            }
            var acc: [4]@Vector(L, F) = @splat(@splat(0));
            var i = from;
            if (s.wide) {
                const W = 2 * L;
                var wacc: [4]@Vector(W, F) = @splat(@splat(0));
                const end = from + ((to - from) & ~@as(usize, W - 1));
                sweep(W, i, end, a, x, y, t, &wacc);
                i = end;
                const lo = std.simd.iota(i32, L);
                const hi = lo + @as(@Vector(L, i32), @splat(L));
                inline for (0..4) |c| acc[c] = @shuffle(F, wacc[c], undefined, lo) + @shuffle(F, wacc[c], undefined, hi);
            }
            sweep(L, i, to, a, x, y, t, &acc);
            finish(L, acc, t2);
        }
    };
}

/// `a*b + c` as the drivers' C contracts it.
inline fn md(comptime F: type, a: F, b: F, c: F) F {
    return fma.madd(fma.fuse, F, a, b, c);
}

/// `a[0]*b[0] + a[1]*b[1] + ...` contracted as the C's chain: the first
/// product fused with the second, each later one into the sum.
inline fn chain(comptime F: type, a: [4]F, b: [4]F) F {
    var r = md(F, a[0], b[0], a[1] * b[1]);
    inline for (2..4) |k| r = md(F, a[k], b[k], r);
    return r;
}

fn Kernel(comptime F: type) type {
    return fn (BlasLong, BlasLong, F, [*]F, BlasLong, [*]F, BlasLong, [*]F, BlasLong, [*]F) callconv(.c) i32;
}

/// Column `j` of `a`.
inline fn col(comptime F: type, a: [*]F, lda: usize, j: usize) [*]F {
    return a + j * lda;
}

/// x86_64 `{s,d}symv_L.c` around block `B`.
pub fn X86L(comptime F: type, comptime B: type) type {
    return struct {
        pub fn k(m: BlasLong, offset: BlasLong, alpha: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, _: [*]F) callconv(.c) i32 {
            if (m <= 0 or offset <= 0) return 0;
            const ld: usize = @intCast(lda);
            if (inc_x != 1 or inc_y != 1) {
                var jx: BlasLong = 0;
                var jy: BlasLong = 0;
                var j: BlasLong = 0;
                while (j < offset) : (j += 1) {
                    const c = col(F, a, ld, @intCast(j));
                    const temp1 = alpha * util.at(F, x, jx).*;
                    var temp2: F = 0;
                    const yj = util.at(F, y, jy);
                    yj.* = md(F, temp1, c[@intCast(j)], yj.*);
                    var ix = jx;
                    var iy = jy;
                    var i = j + 1;
                    while (i < m) : (i += 1) {
                        ix += inc_x;
                        iy += inc_y;
                        const e = c[@intCast(i)];
                        const yi = util.at(F, y, iy);
                        yi.* = md(F, temp1, e, yi.*);
                        temp2 = md(F, e, util.at(F, x, ix).*, temp2);
                    }
                    yj.* = md(F, alpha, temp2, yj.*);
                    jx += inc_x;
                    jy += inc_y;
                }
                return 0;
            }
            const mu: usize = @intCast(m);
            const off: usize = @intCast(offset);
            const offset1 = off / 4 * 4;
            var j: usize = 0;
            while (j < offset1) : (j += 4) {
                var tmp1: [4]F = undefined;
                var tmp2: [4]F = @splat(0);
                var ap: [4][*]F = undefined;
                inline for (0..4) |c| {
                    tmp1[c] = alpha * x[j + c];
                    ap[c] = col(F, a, ld, j + c);
                }
                inline for (0..4) |c| y[j + c] = md(F, tmp1[c], ap[c][j + c], y[j + c]);
                // the triangle below the diagonal block's diagonal
                inline for (0..3) |c| {
                    for (j + c + 1..j + 4) |i| {
                        y[i] = md(F, tmp1[c], ap[c][i], y[i]);
                        tmp2[c] = md(F, ap[c][i], x[i], tmp2[c]);
                    }
                }
                var i = j + 4;
                if (mu - (j + 1) >= 12) {
                    const m2 = mu / 4 * 4;
                    if (m2 > j + 4) B.run(j + 4, m2, &ap, x, y, tmp1, &tmp2);
                    i = m2;
                }
                while (i < mu) : (i += 1) {
                    inline for (0..4) |c| {
                        y[i] = md(F, tmp1[c], ap[c][i], y[i]);
                        tmp2[c] = md(F, ap[c][i], x[i], tmp2[c]);
                    }
                }
                inline for (0..4) |c| y[j + c] = md(F, alpha, tmp2[c], y[j + c]);
            }
            // the C's three ranges have one body
            while (j < off) : (j += 1) {
                const c = col(F, a, ld, j);
                const temp1 = alpha * x[j];
                var temp2: F = 0;
                y[j] = md(F, temp1, c[j], y[j]);
                for (j + 1..mu) |i| {
                    y[i] = md(F, temp1, c[i], y[i]);
                    temp2 = md(F, c[i], x[i], temp2);
                }
                y[j] = md(F, alpha, temp2, y[j]);
            }
            return 0;
        }
    };
}

/// x86_64 `{s,d}symv_U.c` around block `B`, with its C `kernel_1x4` and
/// `kernel_8x1`.
pub fn X86U(comptime F: type, comptime B: type) type {
    return struct {
        fn k1x4(from: usize, to: usize, ap: *const [4][*]F, x: [*]F, y: [*]F, tp: [4]F, temp2: *[4]F) void {
            var tmp2: [4]F = @splat(0);
            for (from..to) |i| {
                const at = [4]F{ ap[0][i], ap[1][i], ap[2][i], ap[3][i] };
                y[i] += chain(F, tp, at);
                inline for (0..4) |c| tmp2[c] = md(F, at[c], x[i], tmp2[c]);
            }
            inline for (0..4) |c| temp2[c] += tmp2[c];
        }

        fn k8x1(n: usize, a0: [*]F, x: [*]F, y: [*]F, t1: F, temp2: *F) void {
            var temp: F = 0;
            for (0..n / 4 * 4) |i| {
                y[i] = md(F, t1, a0[i], y[i]);
                temp = md(F, a0[i], x[i], temp);
            }
            temp2.* = temp;
        }

        pub fn k(m: BlasLong, offset: BlasLong, alpha: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, _: [*]F) callconv(.c) i32 {
            if (m <= 0 or offset <= 0) return 0;
            const ld: usize = @intCast(lda);
            const m1 = m - offset;
            const mrange = m - m1;
            if (inc_x != 1 or inc_y != 1 or mrange < 16) {
                var jx = m1 * inc_x;
                var jy = m1 * inc_y;
                var j = m1;
                while (j < m) : (j += 1) {
                    const c = col(F, a, ld, @intCast(j));
                    const temp1 = alpha * util.at(F, x, jx).*;
                    var temp2: F = 0;
                    var ix: BlasLong = 0;
                    var iy: BlasLong = 0;
                    var i: BlasLong = 0;
                    while (i < j) : (i += 1) {
                        const e = c[@intCast(i)];
                        const yi = util.at(F, y, iy);
                        yi.* = md(F, temp1, e, yi.*);
                        temp2 = md(F, e, util.at(F, x, ix).*, temp2);
                        ix += inc_x;
                        iy += inc_y;
                    }
                    const yj = util.at(F, y, jy);
                    yj.* += md(F, temp1, c[@intCast(j)], alpha * temp2);
                    jx += inc_x;
                    jy += inc_y;
                }
                return 0;
            }
            const mu: usize = @intCast(m);
            const m2: usize = @intCast(m - @rem(mrange, 4));
            var j: usize = @intCast(m1);
            while (j < m2) : (j += 4) {
                var tmp1: [4]F = undefined;
                var tmp2: [4]F = @splat(0);
                var ap: [4][*]F = undefined;
                inline for (0..4) |c| {
                    tmp1[c] = alpha * x[j + c];
                    ap[c] = col(F, a, ld, j + c);
                }
                const j1 = j / 8 * 8;
                if (j1 != 0) B.run(0, j1, &ap, x, y, tmp1, &tmp2);
                if (j1 < j) k1x4(j1, j, &ap, x, y, tmp1, &tmp2);
                inline for (0..4) |c| {
                    const temp1 = tmp1[c];
                    var temp2 = tmp2[c];
                    for (j..j + c) |i| {
                        y[i] = md(F, temp1, ap[c][i], y[i]);
                        temp2 = md(F, ap[c][i], x[i], temp2);
                    }
                    y[j + c] += md(F, temp1, ap[c][j + c], alpha * temp2);
                }
            }
            while (j < mu) : (j += 1) {
                const c = col(F, a, ld, j);
                const temp1 = alpha * x[j];
                var temp2: F = 0;
                const j1 = j / 8 * 8;
                if (j1 != 0) k8x1(j1, c, x, y, temp1, &temp2);
                for (j1..j) |i| {
                    y[i] = md(F, temp1, c[i], y[i]);
                    temp2 = md(F, c[i], x[i], temp2);
                }
                y[j] += md(F, temp1, c[j], alpha * temp2);
            }
            return 0;
        }
    };
}

/// `COPY_K(n, src, inc, dst, 1)`.
fn gather(comptime F: type, n: usize, src: [*]F, inc: BlasLong, dst: [*]F) void {
    for (0..n) |i| dst[i] = util.at(F, src, @as(BlasLong, @intCast(i)) * inc).*;
}

/// `COPY_K(n, src, 1, dst, inc)`.
fn scatter(comptime F: type, n: usize, src: [*]F, dst: [*]F, inc: BlasLong) void {
    for (0..n) |i| util.at(F, dst, @as(BlasLong, @intCast(i)) * inc).* = src[i];
}

/// arm64's drivers copy a strided `y`, then a strided `x`, to the buffer
/// (`2m` elements at most) and work at unit stride.
const Unit = struct {
    fn of(comptime F: type, m: usize, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, buffer: [*]F) [2][*]F {
        var yy = y;
        var xx = x;
        if (inc_y != 1) {
            yy = buffer;
            gather(F, m, y, inc_y, yy);
        }
        if (inc_x != 1) {
            xx = if (inc_y != 1) yy + m else buffer;
            gather(F, m, x, inc_x, xx);
        }
        return .{ xx, yy };
    }
};

/// arm64 `symv_L_asimd_4x4.c` around block `B`.
pub fn ArmL(comptime F: type, comptime B: type) type {
    return struct {
        pub fn k(m: BlasLong, offset: BlasLong, alpha: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, buffer: [*]F) callconv(.c) i32 {
            if (m <= 0) return 0;
            const mu: usize = @intCast(m);
            const off: usize = @intCast(@max(offset, 0));
            const ld: usize = @intCast(lda);
            const xx, const yy = Unit.of(F, mu, x, inc_x, y, inc_y, buffer);
            const offset1 = off / 4 * 4;
            var j: usize = 0;
            while (j < offset1) : (j += 4) {
                var ap: [4][*]F = undefined;
                var xs: [4]F = undefined;
                inline for (0..4) |c| {
                    ap[c] = col(F, a, ld, j + c);
                    xs[c] = xx[j + c];
                }
                // the diagonal block's rows, from the stored triangle
                var tmp2: [4]F = undefined;
                inline for (0..4) |r| {
                    var row: [4]F = undefined;
                    inline for (0..4) |c| row[c] = if (c <= r) ap[c][j + r] else ap[r][j + c];
                    tmp2[r] = chain(F, row, xs);
                }
                var tmp1: [4]F = undefined;
                inline for (0..4) |c| tmp1[c] = alpha * xs[c];
                const m2 = mu / 4 * 4;
                if (m2 > j + 4) B.run(j + 4, m2, &ap, xx, yy, tmp1, &tmp2);
                var i = @max(m2, j + 4);
                while (i < mu) : (i += 1) {
                    inline for (0..4) |c| {
                        yy[i] = md(F, tmp1[c], ap[c][i], yy[i]);
                        tmp2[c] = md(F, ap[c][i], xx[i], tmp2[c]);
                    }
                }
                inline for (0..4) |c| yy[j + c] = md(F, alpha, tmp2[c], yy[j + c]);
            }
            while (j < off) : (j += 1) {
                const c = col(F, a, ld, j);
                const temp1 = alpha * xx[j];
                var temp2: F = 0;
                yy[j] = md(F, temp1, c[j], yy[j]);
                for (j + 1..mu) |i| {
                    yy[i] = md(F, temp1, c[i], yy[i]);
                    temp2 = md(F, c[i], xx[i], temp2);
                }
                yy[j] = md(F, alpha, temp2, yy[j]);
            }
            if (inc_y != 1) scatter(F, mu, yy, y, inc_y);
            return 0;
        }
    };
}

/// arm64 `symv_U_asimd_4x4.c` around block `B`. Where the first column is
/// not a multiple of 4 it drops the rows from the last multiple of 4 up to
/// it (`j1 = (j/4)*4`); they are taken here, as `symv_U_sve_v1x4.c`'s block
/// on `0..j` takes them (`issues.md`).
pub fn ArmU(comptime F: type, comptime B: type) type {
    return struct {
        pub fn k(m: BlasLong, offset: BlasLong, alpha: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, buffer: [*]F) callconv(.c) i32 {
            if (m <= 0) return 0;
            const mu: usize = @intCast(m);
            const ld: usize = @intCast(lda);
            const m1: usize = @intCast(@max(m - offset, 0));
            const off = mu - m1;
            const xx, const yy = Unit.of(F, mu, x, inc_x, y, inc_y, buffer);
            const m2 = mu - off % 4;
            var j = m1;
            while (j < m2) : (j += 4) {
                var tmp1: [4]F = undefined;
                var tmp2: [4]F = @splat(0);
                var ap: [4][*]F = undefined;
                inline for (0..4) |c| {
                    tmp1[c] = alpha * xx[j + c];
                    ap[c] = col(F, a, ld, j + c);
                }
                const j1 = j / 4 * 4;
                if (j1 != 0) B.run(0, j1, &ap, xx, yy, tmp1, &tmp2);
                for (j1..j) |i| {
                    inline for (0..4) |c| {
                        yy[i] = md(F, tmp1[c], ap[c][i], yy[i]);
                        tmp2[c] = md(F, ap[c][i], xx[i], tmp2[c]);
                    }
                }
                inline for (0..4) |c| {
                    const temp1 = tmp1[c];
                    var temp2 = tmp2[c];
                    for (j..j + c) |i| {
                        yy[i] = md(F, temp1, ap[c][i], yy[i]);
                        temp2 = md(F, ap[c][i], xx[i], temp2);
                    }
                    yy[j + c] += md(F, temp1, ap[c][j + c], alpha * temp2);
                }
            }
            while (j < mu) : (j += 1) {
                const c = col(F, a, ld, j);
                const temp1 = alpha * xx[j];
                var temp2: F = 0;
                for (0..j) |i| {
                    yy[i] = md(F, temp1, c[i], yy[i]);
                    temp2 = md(F, c[i], xx[i], temp2);
                }
                yy[j] += md(F, temp1, c[j], alpha * temp2);
            }
            if (inc_y != 1) scatter(F, mu, yy, y, inc_y);
            return 0;
        }
    };
}

/// The element type of a real `symv` slot, and whether it is the lower one.
pub fn form(comptime slot: Slot) struct { F: type, lower: bool } {
    return switch (slot) {
        .ssymv_L => .{ .F = f32, .lower = true },
        .ssymv_U => .{ .F = f32, .lower = false },
        .dsymv_L => .{ .F = f64, .lower = true },
        .dsymv_U => .{ .F = f64, .lower = false },
        else => @compileError(@tagName(slot) ++ " is not a real symv slot"),
    };
}

/// The real `symv` slots.
pub const slots = [_]Slot{ .ssymv_U, .ssymv_L, .dsymv_U, .dsymv_L };

/// The x86_64 kernel of `slot` with block shape `s`.
pub fn x86(comptime slot: Slot, comptime s: Shape) *const Kernel(form(slot).F) {
    const f = form(slot);
    const B = Block(f.F, s);
    return if (f.lower) &X86L(f.F, B).k else &X86U(f.F, B).k;
}

/// The arm64 kernel of `slot`: the asimd driver and block, NEON-wide and fused.
pub fn arm64(comptime slot: Slot) *const Kernel(form(slot).F) {
    const f = form(slot);
    const B = Block(f.F, .{ .lanes = 16 / @sizeOf(f.F), .fused = true });
    return if (f.lower) &ArmL(f.F, B).k else &ArmU(f.F, B).k;
}
