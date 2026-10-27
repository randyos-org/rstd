//! The drivers of `kernel/x86_64/{s,d}gemv_n_4.c` and `{s,d}gemv_t_4.c`,
//! around the blocks a core's `uk/index_n.zig` and `uk/index_t.zig` name.
//! Each block function is the microkernel's where it has one, else the
//! driver's own assembly (`uk/<p>gemv_<n|t>_4.zig`), else the driver's C, as
//! the C's `HAVE_KERNEL_*` tests pick. Nothing here depends on the BLAS root,
//! so the per-core objects (`_kernel.zig`) can use it.
//!
//! The C statements contract as the core's flags say (`fma.madd`): in
//! `t += a0*x0 + a1*x1` the sum fuses its first product and adds to `t`
//! unfused; `y += a*x` fuses. Called with `m` and `n` at least 1, as the
//! interface does; the reference's quick return stays.

const std = @import("std");
const sysinfo = @import("builtin");
const util = @import("../../util.zig");
const types = @import("../../types.zig");
const fma = @import("../../fma.zig");
const uk = @import("../../blas1/uk_block.zig");

const BlasLong = types.BlasLong;
const Slot = util.Slot;

/// The kernels an architecture's `<arch>/kernel.zig` iterates over: the real
/// ones driven here, and the complex ones only converted assembly fills.
pub const slots = [_]Slot{ .sgemv_n, .sgemv_t, .dgemv_n, .dgemv_t, .zgemv_n, .zgemv_t, .zgemv_r, .zgemv_c, .zgemv_o, .zgemv_u, .zgemv_s, .zgemv_d };

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
fn nbmax(comptime slot: Slot, comptime core: []const u8) BlasLong {
    if (Elem(slot) == f64) return 2048;
    const short = std.mem.eql(u8, core, "steamroller") or std.mem.eql(u8, core, "excavator");
    return if (short) 2048 else 4096;
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
const has_sse3 = sysinfo.target.cpu.arch.isX86() and std.Target.x86.featureSetHas(sysinfo.target.cpu.features, .sse3);

/// `gemv_n_4.c`'s `CNAME` of `slot` on `core`, around blocks `M` (the
/// microkernel's) and `O` (the driver's own).
fn GemvN(comptime slot: Slot, comptime core: []const u8, comptime M: type, comptime O: type, comptime fuse: bool) type {
    const F = Elem(slot);
    const p = if (F == f32) "s" else "d";
    const C = Fallback(F, fuse);
    const k4x8 = if (F == f32) choose(N4x8(F), M, O, p ++ "gemv_kernel_4x8", &C.n4x8) else undefined;
    const k4x4 = choose(N4(F), M, O, p ++ "gemv_kernel_4x4", &C.n4x4);
    const k4x2 = choose(N4(F), M, O, p ++ "gemv_kernel_4x2", null);
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

        pub fn call(m: BlasLong, n: BlasLong, _: BlasLong, alpha_in: F, a_in: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, buffer: [*]F) callconv(.c) i32 {
            if (m < 1 or n < 1) return 0;
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
                        k4x2(nb, &ap, x_ptr, ybuffer, pa);
                        a_ptr += @intCast(lda * 2);
                        x_ptr += 2;
                    }
                    if (n2 & 1 != 0) k4x1(nb, a_ptr, x_ptr, ybuffer, pa);
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
/// microkernel's) and `O` (the driver's own).
fn GemvT(comptime slot: Slot, comptime core: []const u8, comptime M: type, comptime O: type, comptime fuse: bool) type {
    const F = Elem(slot);
    const p = if (F == f32) "s" else "d";
    const C = Fallback(F, fuse);
    const own = if (has_sse3) O else void;
    const k4x4 = choose(T4(F), M, own, p ++ "gemv_kernel_4x4", &C.t4x4);
    const k4x2 = choose(T2(F), M, own, p ++ "gemv_kernel_4x2", if (F == f64) &C.t4x2 else null);
    const k4x1 = choose(T1(F), M, own, p ++ "gemv_kernel_4x1", if (F == f64) &C.t4x1 else null);
    const nb_max = nbmax(slot, core);
    return struct {
        /// `y` (stride `inc_dest`) plus `da` times `src`: the C where the
        /// stride is not 1, else the assembly's unfused SSE products.
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
            for (0..@intCast(n)) |i| dest[i] += src[i] * da;
        }

        pub fn call(m: BlasLong, n: BlasLong, _: BlasLong, alpha: F, a_in: [*]F, lda: BlasLong, x_in: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, buffer: [*]F) callconv(.c) i32 {
            if (m < 1 or n < 1) return 0;
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

/// The kernel `slot` of `core`, whose index entry is `e` (`uk_block.Kernels`).
pub fn entryOf(comptime slot: Slot, comptime core: []const u8, comptime e: uk.Entry, comptime Vector: fn (comptime Slot) type) Kernel(slot) {
    const M = switch (e.block) {
        .assembly => e.module,
        .vector => Vector(slot),
        .none => struct {},
    };
    return switch (slot) {
        .sgemv_n, .dgemv_n => &GemvN(slot, core, M, e.own, fma.fuse).call,
        .sgemv_t, .dgemv_t => &GemvT(slot, core, M, e.own, fma.fuse).call,
        else => @compileError(@tagName(slot) ++ " is not driven here"),
    };
}
