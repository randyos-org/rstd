//! The drivers of `kernel/x86_64/{s,d,c,z}gemv_n_4.c` and `gemv_t_4.c` around
//! the blocks `uk/index_n.zig` and `uk/index_t.zig` name (microkernel, the
//! driver's assembly, else its C, as `HAVE_KERNEL_*` picks). The C contracts
//! by the core's flags (`fma.madd`); callers pass `m`, `n` of at least 1. No
//! BLAS-root imports (`_kernel.zig`).

const std = @import("std");
const sysinfo = @import("builtin");
const util = @import("../../util.zig");
const types = @import("../../types.zig");
const fma = @import("../../fma.zig");
const uk = @import("../../uk_block.zig");
const coretypes = @import("../../cpu/coretypes.zig");

const BlasLong = types.BlasLong;
const Slot = util.Slot;
const Core = coretypes.X86CoreType;

/// The kernels an architecture's `<arch>/kernel.zig` iterates over.
pub const slots = [_]Slot{
    .sgemv_n, .sgemv_t, .dgemv_n, .dgemv_t,
    .cgemv_n, .cgemv_t, .cgemv_r, .cgemv_c,
    .cgemv_o, .cgemv_u, .cgemv_s, .cgemv_d,
    .zgemv_n, .zgemv_t, .zgemv_r, .zgemv_c,
    .zgemv_o, .zgemv_u, .zgemv_s, .zgemv_d,
};

/// The element type of `slot`, one of the real `slots`.
pub fn Elem(comptime slot: Slot) type {
    return switch (slot) {
        .sgemv_n, .sgemv_t => f32,
        .dgemv_n, .dgemv_t => f64,
        else => @compileError(@tagName(slot) ++ " is not a real gemv kernel"),
    };
}

/// Vtable slot type of `slot`.
pub fn Kernel(comptime slot: Slot) type {
    return util.SlotFn(slot);
}

/// `NBMAX`, the rows a block takes: `sgemv` halves it on Steamroller and
/// Excavator.
fn nbmax(comptime slot: Slot, comptime core: Core) BlasLong {
    if (Elem(slot) == f64) return 2048;
    return switch (core) {
        .steamroller, .excavator => 2048,
        else => 4096,
    };
}

/// `p` moved by `k` elements, either way.
inline fn off(comptime F: type, p: [*]F, k: BlasLong) [*]F {
    return @ptrCast(util.at(F, p, k));
}

/// The block types, as the C declares them.
fn N4x8(comptime F: type) type {
    return *const fn (BlasLong, [*]const [*]F, [*]F, [*]F, BlasLong, [*]const F) void;
}
fn N4(comptime F: type) type {
    return *const fn (BlasLong, [*]const [*]F, [*]F, [*]F, [*]const F) void;
}
fn N1(comptime F: type) type {
    return *const fn (BlasLong, [*]F, [*]F, [*]F, [*]const F) void;
}
fn T4(comptime F: type) type {
    return *const fn (BlasLong, [*]const [*]F, [*]F, [*]F) void;
}
fn T2(comptime F: type) type {
    return *const fn (BlasLong, [*]F, [*]F, [*]F, [*]F) void;
}
fn T1(comptime F: type) type {
    return *const fn (BlasLong, [*]F, [*]F, [*]F) void;
}

/// `<p>gemv_kernel_<shape>` from the microkernel `M`, else from the driver's
/// own assembly `O`, else `fallback`, the driver's C.
fn choose(comptime Fn: type, comptime M: type, comptime O: type, comptime name: []const u8, comptime fallback: ?Fn) Fn {
    if (@hasDecl(M, name)) return &@field(M, name);
    if (O != void and @hasDecl(O, name)) return &@field(O, name);
    return fallback orelse @compileError("no block " ++ name);
}

/// A chain `p0 + p1 + ... ` of the products `a[k]*b[k]`, as clang contracts
/// it: the first product fused with the second, each later one fused into
/// the running sum.
inline fn chain(comptime F: type, comptime fuse: bool, a: anytype, b: anytype) F {
    var r = fma.madd(fuse, F, a[0], b[0], a[1] * b[1]);
    inline for (2..a.len) |k| r = fma.madd(fuse, F, a[k], b[k], r);
    return r;
}

/// The C blocks of the drivers (`#ifndef HAVE_KERNEL_4x4`), and on a core
/// without SSE3 the `t` drivers' own `4x2` and `4x1`, whose assembly needs
/// it (`haddpd`), as their SSE2 arithmetic.
pub fn Fallback(comptime F: type, comptime fuse: bool) type {
    return struct {
        /// `gemv_n_4.c`'s `kernel_4x4`: `x` times `alpha` first.
        fn n4x4(n: BlasLong, ap: [*]const [*]F, xo: [*]F, y: [*]F, alpha: [*]const F) void {
            var x: [4]F = undefined;
            for (0..4) |i| x[i] = xo[i] * alpha[0];
            for (0..@intCast(n)) |i| {
                y[i] += chain(F, fuse, [4]F{ ap[0][i], ap[1][i], ap[2][i], ap[3][i] }, x);
            }
        }

        /// `sgemv_n_4.c`'s `kernel_4x8`: the second four columns `lda4` on.
        fn n4x8(n: BlasLong, ap: [*]const [*]F, xo: [*]F, y: [*]F, lda4: BlasLong, alpha: [*]const F) void {
            var x: [8]F = undefined;
            for (0..8) |i| x[i] = xo[i] * alpha[0];
            const l: usize = @intCast(lda4);
            for (0..@intCast(n)) |i| {
                y[i] += chain(F, fuse, [4]F{ ap[0][i], ap[1][i], ap[2][i], ap[3][i] }, x[0..4].*);
                y[i] += chain(F, fuse, [4]F{ ap[0][l + i], ap[1][l + i], ap[2][l + i], ap[3][l + i] }, x[4..8].*);
            }
        }

        /// `gemv_t_4.c`'s `kernel_4x4`: four dot products, four at a time.
        fn t4x4(n: BlasLong, ap: [*]const [*]F, x: [*]F, y: [*]F) void {
            var temp: [4]F = @splat(0);
            var i: usize = 0;
            while (i < n) : (i += 4) {
                inline for (0..4) |c| temp[c] += chain(F, fuse, ap[c][i..][0..4].*, x[i..][0..4].*);
            }
            y[0..4].* = temp;
        }

        const V = @Vector(2, F);

        /// `dgemv_t_4.c`'s `kernel_4x2` in SSE2: two lanes per column, a pair
        /// at a time, then the halves added.
        pub fn t4x2(n: BlasLong, ap0: [*]F, ap1: [*]F, x: [*]F, y: [*]F) void {
            var acc: [2]V = @splat(@splat(0));
            const len: usize = @intCast(n);
            var i: usize = 0;
            if (len & 2 != 0) {
                const xv = util.load(V, x);
                acc[0] += util.load(V, ap0) * xv;
                acc[1] += util.load(V, ap1) * xv;
                i = 2;
            }
            // the loop's two pairs go to the same sums in turn
            while (i < len) : (i += 2) {
                const xv = util.load(V, x + i);
                acc[0] += util.load(V, ap0 + i) * xv;
                acc[1] += util.load(V, ap1 + i) * xv;
            }
            y[0] = acc[0][0] + acc[0][1];
            y[1] = acc[1][0] + acc[1][1];
        }

        /// `dgemv_t_4.c`'s `kernel_4x1` in SSE2: two accumulators of two lanes.
        pub fn t4x1(n: BlasLong, ap: [*]F, x: [*]F, y: [*]F) void {
            var lo: V = @splat(0);
            var hi: V = @splat(0);
            const len: usize = @intCast(n);
            var i: usize = 0;
            if (len & 2 != 0) {
                lo += util.load(V, ap) * util.load(V, x);
                i = 2;
            }
            while (i < len) : (i += 4) {
                lo += util.load(V, ap + i) * util.load(V, x + i);
                hi += util.load(V, ap + i + 2) * util.load(V, x + i + 2);
            }
            const s = lo + hi;
            y[0] = s[0] + s[1];
        }
    };
}

/// Whether the core being built has SSE3, which the `t` drivers' own
/// assembly needs.
const has_sse3 = util.x86Has(.sse3);

/// `gemv_n_4.c`'s `CNAME` of `slot` on `core`, around blocks `M` (the
/// microkernel's) and `O` (the driver's own), and for `sgemv` on Skylake-X
/// `S`'s (`sgemv_n_microk_skylakex-8.c`: `n128`, `n64`), else void.
fn GemvN(comptime slot: Slot, comptime core: Core, comptime M: type, comptime O: type, comptime fuse: bool, comptime S: type) type {
    const F = Elem(slot);
    const p = if (F == f32) "s" else "d";
    const C = Fallback(F, fuse);
    const k4x8 = if (F == f32) choose(N4x8(F), M, O, p ++ "gemv_kernel_4x8", &C.n4x8) else undefined;
    const k4x4 = choose(N4(F), M, O, p ++ "gemv_kernel_4x4", &C.n4x4);
    const k4x2 = if (S != void) undefined else choose(N4(F), M, O, p ++ "gemv_kernel_4x2", null);
    const k4x1 = choose(N1(F), M, O, p ++ "gemv_kernel_4x1", null);
    const nb_max = nbmax(slot, core);
    return struct {
        /// `y` (stride `inc`) plus `src`, as `add_y` does where `inc` is not 1.
        fn addY(n: BlasLong, src: [*]F, dest: [*]F, inc_dest: BlasLong) void {
            var d: BlasLong = 0;
            for (0..@intCast(n)) |i| {
                util.at(F, dest, d).* += src[i];
                d += inc_dest;
            }
        }

        /// Skylake-X's whole product for few columns (`HAVE_SGEMV_N_SKYLAKE_KERNEL`):
        /// `x` copied to the buffer, `y` after it, where strided.
        fn small(m: BlasLong, n: BlasLong, alpha: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, buffer: [*]F) void {
            const nu: usize = @intCast(n);
            const mu: usize = @intCast(m);
            var xb = x;
            var yb = y;
            if (inc_x != 1) {
                xb = buffer;
                for (0..nu) |i| xb[i] = util.at(F, x, @as(BlasLong, @intCast(i)) * inc_x).*;
            }
            if (inc_y != 1) {
                yb = buffer + nu;
                for (0..mu) |i| yb[i] = util.at(F, y, @as(BlasLong, @intCast(i)) * inc_y).*;
            }
            S.n128(m, n, alpha, a, lda, xb, yb);
            if (inc_y != 1) for (0..mu) |i| {
                util.at(F, y, @as(BlasLong, @intCast(i)) * inc_y).* = yb[i];
            };
        }

        pub fn call(m: BlasLong, n: BlasLong, _: BlasLong, alpha_in: F, a_in: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, buffer: [*]F) callconv(.c) i32 {
            if (m < 1 or n < 1) return 0;
            if (S != void and m <= 16384 and n <= 48 and n != 4) {
                small(m, n, alpha_in, a_in, lda, x, inc_x, y, inc_y, buffer);
                return 0;
            }
            var alpha = alpha_in;
            const pa: [*]const F = @ptrCast(&alpha);
            var a = a_in;
            var ybuffer = buffer;
            const lda4 = lda << 2;
            const lda8 = lda << 3;
            var xbuffer: [8]F = undefined;

            // `sgemv` takes eight columns a step from a contiguous `x`
            const eight = F == f32 and inc_x == 1;
            const n1 = if (eight) n >> 3 else n >> 2;
            const n2 = if (eight) n & 7 else n & 3;
            const m3 = m & 3;
            var m1 = m & -4;
            const m2 = (m & (nb_max - 1)) - m3;

            var y_ptr = y;
            var nb = nb_max;
            while (nb == nb_max) {
                m1 -= nb;
                if (m1 < 0) {
                    if (m2 == 0) break;
                    nb = m2;
                }
                var a_ptr = a;
                var x_ptr = x;
                var ap: [4][*]F = .{ a_ptr, a_ptr + @as(usize, @intCast(lda)), undefined, undefined };
                ap[2] = ap[1] + @as(usize, @intCast(lda));
                ap[3] = ap[2] + @as(usize, @intCast(lda));
                if (inc_y != 1) @memset(ybuffer[0..@intCast(nb)], 0) else ybuffer = y_ptr;

                if (inc_x == 1) {
                    const step: BlasLong = if (F == f32) lda8 else lda4;
                    var i: BlasLong = 0;
                    while (i < n1) : (i += 1) {
                        if (F == f32) k4x8(nb, &ap, x_ptr, ybuffer, lda4, pa) else k4x4(nb, &ap, x_ptr, ybuffer, pa);
                        for (&ap) |*q| q.* += @intCast(step);
                        a_ptr += @intCast(step);
                        x_ptr += if (F == f32) 8 else 4;
                    }
                    if (F == f32 and n2 & 4 != 0) {
                        k4x4(nb, &ap, x_ptr, ybuffer, pa);
                        for (&ap) |*q| q.* += @intCast(lda4);
                        a_ptr += @intCast(lda4);
                        x_ptr += 4;
                    }
                    if (n2 & 2 != 0) {
                        if (S != void) S.n64(nb, 2, alpha, a_ptr, lda, x_ptr, ybuffer) else k4x2(nb, &ap, x_ptr, ybuffer, pa);
                        a_ptr += @intCast(lda * 2);
                        x_ptr += 2;
                    }
                    if (n2 & 1 != 0) {
                        if (S != void) S.n64(nb, 1, alpha, a_ptr, lda, x_ptr, ybuffer) else k4x1(nb, a_ptr, x_ptr, ybuffer, pa);
                    }
                } else {
                    var i: BlasLong = 0;
                    while (i < n1) : (i += 1) {
                        for (0..4) |k| {
                            xbuffer[k] = x_ptr[0];
                            x_ptr = off(F, x_ptr, inc_x);
                        }
                        k4x4(nb, &ap, &xbuffer, ybuffer, pa);
                        for (&ap) |*q| q.* += @intCast(lda4);
                        a_ptr += @intCast(lda4);
                    }
                    i = 0;
                    while (i < n2) : (i += 1) {
                        xbuffer[0] = x_ptr[0];
                        x_ptr = off(F, x_ptr, inc_x);
                        k4x1(nb, a_ptr, &xbuffer, ybuffer, pa);
                        a_ptr += @intCast(lda);
                    }
                }

                a += @intCast(nb);
                if (inc_y != 1) {
                    addY(nb, ybuffer, y_ptr, inc_y);
                    y_ptr = off(F, y_ptr, nb * inc_y);
                } else y_ptr += @intCast(nb);
            }

            if (m3 == 0) return 0;
            // the last rows, each a dot product of its row with `x`
            const rows: usize = @intCast(m3);
            var temp: [3]F = @splat(0);
            var a_ptr = a;
            var x_ptr = x;
            if (lda == m3 and inc_x == 1) {
                const nn: usize = @intCast(n);
                var i: usize = 0;
                if (rows == 1) {
                    while (i < (nn & ~@as(usize, 3))) : (i += 4) temp[0] += chain(F, fuse, a_ptr[i..][0..4].*, x_ptr[i..][0..4].*);
                    while (i < nn) : (i += 1) temp[0] = fma.madd(fuse, F, a_ptr[i], x_ptr[i], temp[0]);
                } else {
                    // two columns a statement, as the C pairs them
                    while (i < (nn & ~@as(usize, 3))) : (i += 4) {
                        for (0..2) |h| {
                            for (0..rows) |r| {
                                const c = 2 * h;
                                temp[r] += fma.madd(fuse, F, a_ptr[c * rows + r], x_ptr[c], a_ptr[(c + 1) * rows + r] * x_ptr[c + 1]);
                            }
                        }
                        a_ptr += 4 * rows;
                        x_ptr += 4;
                    }
                    while (i < nn) : (i += 1) {
                        for (0..rows) |r| temp[r] = fma.madd(fuse, F, a_ptr[r], x_ptr[0], temp[r]);
                        a_ptr += rows;
                        x_ptr += 1;
                    }
                }
            } else {
                var i: BlasLong = 0;
                while (i < n) : (i += 1) {
                    for (0..rows) |r| temp[r] = fma.madd(fuse, F, a_ptr[r], x_ptr[0], temp[r]);
                    a_ptr += @intCast(lda);
                    x_ptr = off(F, x_ptr, inc_x);
                }
            }
            for (0..rows) |r| {
                y_ptr[0] = fma.madd(fuse, F, alpha, temp[r], y_ptr[0]);
                if (r + 1 < rows) y_ptr = off(F, y_ptr, inc_y);
            }
            return 0;
        }
    };
}

/// `gemv_t_4.c`'s `CNAME` of `slot` on `core`, around blocks `M` (the
/// microkernel's) and `O` (the driver's own), and for `sgemv` on Skylake-X
/// `S`'s (`sgemv_t_microk_skylakex.c`: `t`), else void.
fn GemvT(comptime slot: Slot, comptime core: Core, comptime M: type, comptime O: type, comptime fuse: bool, comptime S: type) type {
    const F = Elem(slot);
    const p = if (F == f32) "s" else "d";
    const C = Fallback(F, fuse);
    const own = if (has_sse3) O else void;
    const k4x4 = choose(T4(F), M, own, p ++ "gemv_kernel_4x4", &C.t4x4);
    const k4x2 = choose(T2(F), M, own, p ++ "gemv_kernel_4x2", if (F == f64) &C.t4x2 else null);
    const k4x1 = choose(T1(F), M, own, p ++ "gemv_kernel_4x1", if (F == f64) &C.t4x1 else null);
    const add_asm = own != void and @hasDecl(own, "add_y_asm");
    const nb_max = nbmax(slot, core);
    return struct {
        /// `y` (stride `inc_dest`) plus `da` times `src`: the C where the
        /// stride is not 1, else the assembly's (`add_y_asm`) unfused SSE
        /// products, the same arithmetic in Zig where it is left out.
        fn addY(n: BlasLong, da: F, src: [*]F, dest: [*]F, inc_dest: BlasLong) void {
            if (inc_dest != 1) {
                var d: BlasLong = 0;
                for (0..@intCast(n)) |i| {
                    const v = util.at(F, dest, d);
                    v.* = fma.madd(fuse, F, src[i], da, v.*);
                    d += inc_dest;
                }
                return;
            }
            if (add_asm) return own.add_y_asm(n, da, src, dest, inc_dest);
            for (0..@intCast(n)) |i| dest[i] += src[i] * da;
        }

        /// Skylake-X's whole product for a contiguous `A` of up to 8 rows
        /// (`HAVE_SGEMV_T_SKYLAKE_KERNEL`): `x` copied to the buffer, `y`
        /// after it, where strided.
        fn small(m: BlasLong, n: BlasLong, alpha: F, a: [*]F, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, buffer: [*]F) void {
            const mu: usize = @intCast(m);
            const nu: usize = @intCast(n);
            var xb = x;
            var yb = y;
            if (inc_x != 1) {
                xb = buffer;
                for (0..mu) |i| xb[i] = util.at(F, x, @as(BlasLong, @intCast(i)) * inc_x).*;
            }
            if (inc_y != 1) {
                yb = buffer + mu;
                for (0..nu) |i| yb[i] = util.at(F, y, @as(BlasLong, @intCast(i)) * inc_y).*;
            }
            S.t(mu, nu, alpha, a, xb, yb);
            if (inc_y != 1) for (0..nu) |i| {
                util.at(F, y, @as(BlasLong, @intCast(i)) * inc_y).* = yb[i];
            };
        }

        pub fn call(m: BlasLong, n: BlasLong, _: BlasLong, alpha: F, a_in: [*]F, lda: BlasLong, x_in: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, buffer: [*]F) callconv(.c) i32 {
            if (m < 1 or n < 1) return 0;
            if (S != void and lda == m and n <= 16384 and m <= 8) {
                small(m, n, alpha, a_in, x_in, inc_x, y, inc_y, buffer);
                return 0;
            }
            var a = a_in;
            var x = x_in;
            var xbuffer = buffer;
            const ytemp = buffer + @as(usize, @intCast(@min(m, nb_max)));
            var ybuffer: [4]F = undefined;

            const n0 = @divTrunc(n, nb_max);
            const n1 = @rem(n, nb_max) >> 2;
            const n2 = n & 3;
            const m3 = m & 3;
            var m1 = m & -4;
            const m2 = (m & (nb_max - 1)) - m3;
            const lda4 = 4 * lda;

            var nb = nb_max;
            while (nb == nb_max) {
                m1 -= nb;
                if (m1 < 0) {
                    if (m2 == 0) break;
                    nb = m2;
                }
                var y_ptr = y;
                var a_ptr = a;
                if (inc_x == 1) xbuffer = x else {
                    var s: BlasLong = 0;
                    for (0..@intCast(nb)) |i| {
                        xbuffer[i] = util.at(F, x, s).*;
                        s += inc_x;
                    }
                }
                var ap: [4][*]F = .{ a_ptr, a_ptr + @as(usize, @intCast(lda)), undefined, undefined };
                ap[2] = ap[1] + @as(usize, @intCast(lda));
                ap[3] = ap[2] + @as(usize, @intCast(lda));

                if (n0 > 0) {
                    const nb1 = @divTrunc(nb_max, 4);
                    var j: BlasLong = 0;
                    while (j < n0) : (j += 1) {
                        var yp = ytemp;
                        var i: BlasLong = 0;
                        while (i < nb1) : (i += 1) {
                            k4x4(nb, &ap, xbuffer, yp);
                            for (&ap) |*q| q.* += @intCast(lda4);
                            yp += 4;
                        }
                        addY(nb1 * 4, alpha, ytemp, y_ptr, inc_y);
                        y_ptr = off(F, y_ptr, nb1 * inc_y * 4);
                        a_ptr += @intCast(nb1 * lda4);
                    }
                }

                var yp = ytemp;
                var i: BlasLong = 0;
                while (i < n1) : (i += 1) {
                    k4x4(nb, &ap, xbuffer, yp);
                    for (&ap) |*q| q.* += @intCast(lda4);
                    yp += 4;
                }
                if (n1 > 0) {
                    addY(n1 * 4, alpha, ytemp, y_ptr, inc_y);
                    y_ptr = off(F, y_ptr, n1 * inc_y * 4);
                    a_ptr += @intCast(n1 * lda4);
                }
                if (n2 & 2 != 0) {
                    k4x2(nb, ap[0], ap[1], xbuffer, &ybuffer);
                    a_ptr += @intCast(lda * 2);
                    y_ptr[0] = fma.madd(fuse, F, ybuffer[0], alpha, y_ptr[0]);
                    y_ptr = off(F, y_ptr, inc_y);
                    y_ptr[0] = fma.madd(fuse, F, ybuffer[1], alpha, y_ptr[0]);
                    y_ptr = off(F, y_ptr, inc_y);
                }
                if (n2 & 1 != 0) {
                    k4x1(nb, a_ptr, xbuffer, &ybuffer);
                    y_ptr[0] = fma.madd(fuse, F, ybuffer[0], alpha, y_ptr[0]);
                }
                a += @intCast(nb);
                x = off(F, x, nb * inc_x);
            }

            if (m3 == 0) return 0;
            // the last rows of `A`: every column's dot product with those of `x`
            const rows: usize = @intCast(m3);
            var xtemp: [3]F = undefined;
            var x_ptr = x;
            for (0..rows) |r| {
                xtemp[r] = x_ptr[0] * alpha;
                if (r + 1 < rows) x_ptr = off(F, x_ptr, inc_x);
            }
            // every C loop here (packed, `lda` apart, strided `y`) does the
            // same per column
            var aj = a;
            var y_ptr = y;
            for (0..@intCast(n)) |_| {
                y_ptr[0] = column(rows, aj, &xtemp, y_ptr[0]);
                y_ptr = off(F, y_ptr, inc_y);
                aj += @intCast(lda);
            }
            return 0;
        }

        /// `y += aj[0]*xtemp0 + ...` over the last rows, as clang contracts
        /// it: one row fuses into `y`, more make a chain added to it.
        inline fn column(rows: usize, aj: [*]F, xtemp: *const [3]F, yv: F) F {
            return switch (rows) {
                3 => yv + chain(F, fuse, aj[0..3].*, xtemp.*),
                2 => yv + chain(F, fuse, aj[0..2].*, xtemp[0..2].*),
                else => fma.madd(fuse, F, aj[0], xtemp[0], yv),
            };
        }
    };
}

// ---- complex -----------------------------------------------------------------

/// What a complex slot is built as: `trans` is `TRANS`, `conj` `CONJ`
/// (`A` conjugated), `xconj` `XCONJ` (`x` conjugated, written as the product
/// conjugated, then `alpha` applied to its conjugate).
const Cform = struct { trans: bool, conj: bool, xconj: bool };

fn cform(comptime slot: Slot) Cform {
    const t: types.Trans = switch (slot) {
        .cgemv_n, .zgemv_n => .n,
        .cgemv_t, .zgemv_t => .t,
        .cgemv_r, .zgemv_r => .r,
        .cgemv_c, .zgemv_c => .c,
        .cgemv_o, .zgemv_o => .o,
        .cgemv_u, .zgemv_u => .u,
        .cgemv_s, .zgemv_s => .s,
        .cgemv_d, .zgemv_d => .d,
        else => @compileError(@tagName(slot) ++ " is not a complex gemv kernel"),
    };
    return .{ .trans = t.transposed(), .conj = t.conjA(), .xconj = t.conjX() };
}

/// The float type of complex `slot`.
fn Cfloat(comptime slot: Slot) type {
    return if (@tagName(slot)[0] == 'c') f32 else f64;
}

/// The suffix of the converted blocks' form (`uk.zig`'s `Form.suffix`).
fn suffix(comptime f: Cform) []const u8 {
    return if (f.conj and f.xconj) "_cx" else if (f.conj) "_c" else if (f.xconj) "_x" else "";
}

/// `a*b + c*d` and `a*b - c*d` as clang contracts them: the first product
/// fused, the second rounded.
inline fn sum2(comptime F: type, comptime fuse: bool, a: F, b: F, c: F, d: F) F {
    return fma.madd(fuse, F, a, b, c * d);
}
inline fn diff2(comptime F: type, comptime fuse: bool, a: F, b: F, c: F, d: F) F {
    return fma.madd(fuse, F, a, b, -(c * d));
}

/// One element's product terms as the C writes them: `(ar*xr - ai*xi,
/// ar*xi + ai*xr)` where `CONJ` and `XCONJ` agree, else `(ar*xr + ai*xi,
/// ar*xi - ai*xr)`.
inline fn cterm(comptime F: type, comptime fuse: bool, comptime same: bool, ar: F, ai: F, xr: F, xi: F) [2]F {
    return if (same)
        .{ diff2(F, fuse, ar, xr, ai, xi), sum2(F, fuse, ar, xi, ai, xr) }
    else
        .{ sum2(F, fuse, ar, xr, ai, xi), diff2(F, fuse, ar, xi, ai, xr) };
}

/// `y += alpha * t` as the C writes it: plain, or with `XCONJ`
/// `y_r += ar*tr + ai*ti`, `y_i -= ar*ti - ai*tr`.
inline fn cupdate(comptime F: type, comptime fuse: bool, comptime xconj: bool, y: [*]F, ar: F, ai: F, tr: F, ti: F) void {
    if (xconj) {
        y[0] += sum2(F, fuse, ar, tr, ai, ti);
        y[1] -= diff2(F, fuse, ar, ti, ai, tr);
    } else {
        y[0] += diff2(F, fuse, ar, tr, ai, ti);
        y[1] += sum2(F, fuse, ar, ti, ai, tr);
    }
}

/// The C blocks of `{c,z}gemv_n_4.c` and `{c,z}gemv_t_4.c` and `n`'s `add_y`.
fn CFallback(comptime F: type, comptime f: Cform, comptime fuse: bool) type {
    const same = f.conj == f.xconj;
    return struct {
        /// `kernel_4x<cols>` of `n`: `y += A x` over `cols` columns, column by
        /// column for each row.
        fn nBlock(comptime cols: comptime_int, rows: BlasLong, ap: [*]const [*]F, x: [*]F, y: [*]F) void {
            for (0..@intCast(rows)) |r| {
                const i = 2 * r;
                inline for (0..cols) |c| {
                    const t = cterm(F, fuse, same, ap[c][i], ap[c][i + 1], x[2 * c], x[2 * c + 1]);
                    y[i] += t[0];
                    y[i + 1] += t[1];
                }
            }
        }
        fn n4x4(rows: BlasLong, ap: [*]const [*]F, x: [*]F, y: [*]F) void {
            nBlock(4, rows, ap, x, y);
        }
        fn n4x2(rows: BlasLong, ap: [*]const [*]F, x: [*]F, y: [*]F) void {
            nBlock(2, rows, ap, x, y);
        }
        fn n4x1(rows: BlasLong, ap: [*]F, x: [*]F, y: [*]F) void {
            const one = [1][*]F{ap};
            nBlock(1, rows, &one, x, y);
        }

        /// `kernel_4x<cols>` of `t`: `cols` dot products, then `alpha` into `y`.
        fn tBlock(comptime cols: comptime_int, rows: BlasLong, ap: [*]const [*]F, x: [*]F, y: [*]F, alpha: [*]const F) void {
            var tr: [cols]F = @splat(0);
            var ti: [cols]F = @splat(0);
            for (0..@intCast(rows)) |r| {
                const i = 2 * r;
                inline for (0..cols) |c| {
                    const v = cterm(F, fuse, same, ap[c][i], ap[c][i + 1], x[i], x[i + 1]);
                    tr[c] += v[0];
                    ti[c] += v[1];
                }
            }
            inline for (0..cols) |c| cupdate(F, fuse, f.xconj, y + 2 * c, alpha[0], alpha[1], tr[c], ti[c]);
        }
        fn t4x4(rows: BlasLong, ap: [*]const [*]F, x: [*]F, y: [*]F, alpha: [*]const F) void {
            tBlock(4, rows, ap, x, y, alpha);
        }
        fn t4x2(rows: BlasLong, ap: [*]const [*]F, x: [*]F, y: [*]F, alpha: [*]const F) void {
            tBlock(2, rows, ap, x, y, alpha);
        }
        fn t4x1(rows: BlasLong, ap: [*]F, x: [*]F, y: [*]F, alpha: [*]const F) void {
            const one = [1][*]F{ap};
            tBlock(1, rows, &one, x, y, alpha);
        }

        /// `n`'s `add_y`: `dest += alpha * src`, or with `XCONJ` the product's
        /// conjugate, `(ar*sr + ai*si, -ar*si + ai*sr)`.
        fn addY(count: BlasLong, src: [*]F, dest: [*]F, inc_dest: BlasLong, ar: F, ai: F) void {
            var d: BlasLong = 0;
            for (0..@intCast(count)) |k| {
                const sr = src[2 * k];
                const si = src[2 * k + 1];
                const tr = if (f.xconj) sum2(F, fuse, ar, sr, ai, si) else diff2(F, fuse, ar, sr, ai, si);
                const ti = if (f.xconj) sum2(F, fuse, -ar, si, ai, sr) else sum2(F, fuse, ar, si, ai, sr);
                util.at(F, dest, d).* += tr;
                util.at(F, dest, d + 1).* += ti;
                d += inc_dest;
            }
        }
    };
}

fn CN4(comptime F: type) type {
    return *const fn (BlasLong, [*]const [*]F, [*]F, [*]F) void;
}
fn CN1(comptime F: type) type {
    return *const fn (BlasLong, [*]F, [*]F, [*]F) void;
}
fn CT4(comptime F: type) type {
    return *const fn (BlasLong, [*]const [*]F, [*]F, [*]F, [*]const F) void;
}
fn CT1(comptime F: type) type {
    return *const fn (BlasLong, [*]F, [*]F, [*]F, [*]const F) void;
}
fn AddY(comptime F: type) type {
    return *const fn (BlasLong, [*]F, [*]F, BlasLong, F, F) void;
}

/// `{c,z}gemv_n_4.c`'s `CNAME` of `slot`, around the microkernel `M`.
fn ZgemvN(comptime slot: Slot, comptime M: type, comptime fuse: bool) type {
    const F = Cfloat(slot);
    const f = cform(slot);
    const p = @tagName(slot)[0..1];
    const sx = suffix(f);
    const C = CFallback(F, f, fuse);
    const k4x4 = choose(CN4(F), M, void, p ++ "gemv_kernel_4x4" ++ sx, &C.n4x4);
    const k4x2 = choose(CN4(F), M, void, p ++ "gemv_kernel_4x2" ++ sx, &C.n4x2);
    const k4x1 = choose(CN1(F), M, void, p ++ "gemv_kernel_4x1" ++ sx, &C.n4x1);
    // the microkernel's `add_y` past its strided C, which is `C.addY`
    const add_asm: ?AddY(F) = if (@hasDecl(M, "add_y" ++ sx ++ "_asm")) &@field(M, "add_y" ++ sx ++ "_asm") else null;
    const nb_max: BlasLong = if (F == f32) 2048 else 1024;
    return struct {
        fn addY(count: BlasLong, src: [*]F, dest: [*]F, inc_dest: BlasLong, ar: F, ai: F) void {
            if (add_asm) |k| if (inc_dest == 2) return k(count, src, dest, inc_dest, ar, ai);
            C.addY(count, src, dest, inc_dest, ar, ai);
        }

        pub fn call(m: BlasLong, n: BlasLong, _: BlasLong, alpha_r: F, alpha_i: F, a_in: [*]F, lda_in: BlasLong, x: [*]F, inc_x_in: BlasLong, y: [*]F, inc_y_in: BlasLong, buffer: [*]F) callconv(.c) i32 {
            if (m < 1 or n < 1) return 0;
            const ybuffer = buffer;
            const inc_x = inc_x_in * 2;
            const inc_y = inc_y_in * 2;
            const lda = lda_in * 2;
            const lda4 = 4 * lda;
            var a = a_in;
            var xbuffer: [8]F = undefined;

            const n1 = @divTrunc(n, 4);
            const n2 = @rem(n, 4);
            const m3 = @rem(m, 4);
            var m1 = m - m3;
            const m2 = @rem(m, nb_max) - m3;

            var y_ptr = y;
            var nb = nb_max;
            while (nb == nb_max) {
                m1 -= nb;
                if (m1 < 0) {
                    if (m2 == 0) break;
                    nb = m2;
                }
                var a_ptr = a;
                var ap: [4][*]F = undefined;
                for (&ap, 0..) |*q, k| q.* = a_ptr + k * @as(usize, @intCast(lda));
                var x_ptr = x;
                @memset(ybuffer[0..@intCast(2 * nb)], 0);

                if (inc_x == 2) {
                    var i: BlasLong = 0;
                    while (i < n1) : (i += 1) {
                        k4x4(nb, &ap, x_ptr, ybuffer);
                        for (&ap) |*q| q.* += @intCast(lda4);
                        a_ptr += @intCast(lda4);
                        x_ptr += 8;
                    }
                    if (n2 & 2 != 0) {
                        k4x2(nb, &ap, x_ptr, ybuffer);
                        x_ptr += 4;
                        a_ptr += @intCast(2 * lda);
                    }
                    if (n2 & 1 != 0) k4x1(nb, a_ptr, x_ptr, ybuffer);
                } else {
                    var i: BlasLong = 0;
                    while (i < n1) : (i += 1) {
                        for (0..4) |k| {
                            xbuffer[2 * k] = x_ptr[0];
                            xbuffer[2 * k + 1] = x_ptr[1];
                            x_ptr = off(F, x_ptr, inc_x);
                        }
                        k4x4(nb, &ap, &xbuffer, ybuffer);
                        for (&ap) |*q| q.* += @intCast(lda4);
                        a_ptr += @intCast(lda4);
                    }
                    i = 0;
                    while (i < n2) : (i += 1) {
                        xbuffer[0] = x_ptr[0];
                        xbuffer[1] = x_ptr[1];
                        x_ptr = off(F, x_ptr, inc_x);
                        k4x1(nb, a_ptr, &xbuffer, ybuffer);
                        a_ptr += @intCast(lda);
                    }
                }
                addY(nb, ybuffer, y_ptr, inc_y, alpha_r, alpha_i);
                a += @intCast(2 * nb);
                y_ptr = off(F, y_ptr, nb * inc_y);
            }

            if (m3 == 0) return 0;
            // the last rows: each a sum over the columns from 0, in the
            // C's order whichever of its loops (packed or strided) runs
            const rows: usize = @intCast(m3);
            var tr: [3]F = @splat(0);
            var ti: [3]F = @splat(0);
            var a_ptr = a;
            var x_ptr = x;
            for (0..@intCast(n)) |_| {
                for (0..rows) |r| {
                    const v = cterm(F, fuse, f.conj == f.xconj, a_ptr[2 * r], a_ptr[2 * r + 1], x_ptr[0], x_ptr[1]);
                    tr[r] += v[0];
                    ti[r] += v[1];
                }
                a_ptr = off(F, a_ptr, lda);
                x_ptr = off(F, x_ptr, inc_x);
            }
            for (0..rows) |r| {
                cupdate(F, fuse, f.xconj, y_ptr, alpha_r, alpha_i, tr[r], ti[r]);
                if (r + 1 < rows) y_ptr = off(F, y_ptr, inc_y);
            }
            return 0;
        }
    };
}

/// `{c,z}gemv_t_4.c`'s `CNAME` of `slot`, around the microkernel `M`.
fn ZgemvT(comptime slot: Slot, comptime M: type, comptime fuse: bool) type {
    const F = Cfloat(slot);
    const f = cform(slot);
    const p = @tagName(slot)[0..1];
    const sx = suffix(f);
    const C = CFallback(F, f, fuse);
    const k4x4 = choose(CT4(F), M, void, p ++ "gemv_kernel_4x4" ++ sx, &C.t4x4);
    const k4x2 = choose(CT4(F), M, void, p ++ "gemv_kernel_4x2" ++ sx, &C.t4x2);
    const k4x1 = choose(CT1(F), M, void, p ++ "gemv_kernel_4x1" ++ sx, &C.t4x1);
    const nb_max: BlasLong = if (F == f32) 2048 else 1024;
    return struct {
        pub fn call(m: BlasLong, n: BlasLong, _: BlasLong, alpha_r: F, alpha_i: F, a_in: [*]F, lda_in: BlasLong, x_in: [*]F, inc_x_in: BlasLong, y: [*]F, inc_y_in: BlasLong, buffer: [*]F) callconv(.c) i32 {
            if (m < 1 or n < 1) return 0;
            const inc_x = inc_x_in * 2;
            const inc_y = inc_y_in * 2;
            const lda = lda_in * 2;
            const lda4 = lda * 4;
            var a = a_in;
            var x = x_in;
            var xbuffer = buffer;
            var ybuffer: [8]F = undefined;
            const alpha = [2]F{ alpha_r, alpha_i };

            const n1 = n >> 2;
            const n2 = n & 3;
            const m3 = m & 3;
            var m1 = m - m3;
            const m2 = (m & (nb_max - 1)) - m3;

            var nb = nb_max;
            while (nb == nb_max) {
                m1 -= nb;
                if (m1 < 0) {
                    if (m2 == 0) break;
                    nb = m2;
                }
                var y_ptr = y;
                var a_ptr = a;
                var ap: [4][*]F = undefined;
                for (&ap, 0..) |*q, k| q.* = a_ptr + k * @as(usize, @intCast(lda));
                if (inc_x != 2) {
                    var s: BlasLong = 0;
                    for (0..@intCast(nb)) |k| {
                        xbuffer[2 * k] = util.at(F, x, s).*;
                        xbuffer[2 * k + 1] = util.at(F, x, s + 1).*;
                        s += inc_x;
                    }
                } else xbuffer = x;

                if (inc_y == 2) {
                    var i: BlasLong = 0;
                    while (i < n1) : (i += 1) {
                        k4x4(nb, &ap, xbuffer, y_ptr, &alpha);
                        for (&ap) |*q| q.* += @intCast(lda4);
                        a_ptr += @intCast(lda4);
                        y_ptr += 8;
                    }
                    if (n2 & 2 != 0) {
                        k4x2(nb, &ap, xbuffer, y_ptr, &alpha);
                        a_ptr += @intCast(lda * 2);
                        y_ptr += 4;
                    }
                    if (n2 & 1 != 0) k4x1(nb, a_ptr, xbuffer, y_ptr, &alpha);
                } else {
                    var i: BlasLong = 0;
                    while (i < n1) : (i += 1) {
                        ybuffer = @splat(0);
                        k4x4(nb, &ap, xbuffer, &ybuffer, &alpha);
                        for (&ap) |*q| q.* += @intCast(lda4);
                        a_ptr += @intCast(lda4);
                        for (0..4) |k| {
                            y_ptr[0] += ybuffer[2 * k];
                            y_ptr[1] += ybuffer[2 * k + 1];
                            y_ptr = off(F, y_ptr, inc_y);
                        }
                    }
                    i = 0;
                    while (i < n2) : (i += 1) {
                        ybuffer = @splat(0);
                        k4x1(nb, a_ptr, xbuffer, &ybuffer, &alpha);
                        a_ptr += @intCast(lda);
                        y_ptr[0] += ybuffer[0];
                        y_ptr[1] += ybuffer[1];
                        y_ptr = off(F, y_ptr, inc_y);
                    }
                }
                a += @intCast(2 * nb);
                x = off(F, x, nb * inc_x);
            }

            if (m3 == 0) return 0;
            // the last rows: per column, the first row's term, the others added
            const rows: usize = @intCast(m3);
            var xs: [3][2]F = undefined;
            var x_ptr = x;
            for (0..rows) |r| {
                xs[r] = .{ x_ptr[0], x_ptr[1] };
                if (r + 1 < rows) x_ptr = off(F, x_ptr, inc_x);
            }
            var a_ptr = a;
            var y_ptr = y;
            for (0..@intCast(n)) |_| {
                var t = cterm(F, fuse, f.conj == f.xconj, a_ptr[0], a_ptr[1], xs[0][0], xs[0][1]);
                for (1..rows) |r| {
                    const v = cterm(F, fuse, f.conj == f.xconj, a_ptr[2 * r], a_ptr[2 * r + 1], xs[r][0], xs[r][1]);
                    t[0] += v[0];
                    t[1] += v[1];
                }
                cupdate(F, fuse, f.xconj, y_ptr, alpha_r, alpha_i, t[0], t[1]);
                a_ptr = off(F, a_ptr, lda);
                y_ptr = off(F, y_ptr, inc_y);
            }
            return 0;
        }
    };
}

/// The kernel `slot` of `core`, whose index entry is `e` (`uk_block.Kernels`).
pub fn entryOf(comptime slot: Slot, comptime core: Core, comptime e: uk.Entry, comptime Vector: fn (comptime Slot) type) Kernel(slot) {
    const M = switch (e.block) {
        .assembly => e.module,
        .vector => Vector(slot),
        .none => struct {},
    };
    return switch (slot) {
        .sgemv_n, .dgemv_n => &GemvN(slot, core, M, e.own, fma.fuse, void).call,
        .sgemv_t, .dgemv_t => &GemvT(slot, core, M, e.own, fma.fuse, void).call,
        .cgemv_n, .cgemv_r, .cgemv_o, .cgemv_s, .zgemv_n, .zgemv_r, .zgemv_o, .zgemv_s => &ZgemvN(slot, M, fma.fuse).call,
        .cgemv_t, .cgemv_c, .cgemv_u, .cgemv_d, .zgemv_t, .zgemv_c, .zgemv_u, .zgemv_d => &ZgemvT(slot, M, fma.fuse).call,
        else => @compileError(@tagName(slot) ++ " is not driven here"),
    };
}

/// Skylake-X's `sgemv_n`: `sgemv_n_4.c` around Haswell's block `e` (its
/// index entry) and `S`, the intrinsics block (`sgemv_n_microk_skylakex-8.c`).
pub fn skylakeN(comptime e: uk.Entry, comptime S: type) Kernel(.sgemv_n) {
    return &GemvN(.sgemv_n, .skylakex, e.module, e.own, fma.fuse, S).call;
}

/// Skylake-X's `sgemv_t`: `sgemv_t_4.c` around Haswell's block `e` and `S`
/// (`sgemv_t_microk_skylakex.c`).
pub fn skylakeT(comptime e: uk.Entry, comptime S: type) Kernel(.sgemv_t) {
    return &GemvT(.sgemv_t, .skylakex, e.module, e.own, fma.fuse, S).call;
}
