//! MIPS `axpy` on the SIMD architecture: OpenBLAS `kernel/mips/{s,d,c,z}axpy_msa.c`,
//! for the cores whose `KERNEL` file names them.
//!
//! The C is GCC vector-extension arithmetic (`y0 += da_vec * x0` on `v2f64` and
//! `v4f32`) with a few MSA builtins for the de-interleaving of complex data and
//! the strided gathers. Those are Zig `@Vector` operations here, `@shuffle` for
//! the pack and interleave builtins. The branches are the reference's: both
//! strides one, `y` alone unit, and everything else scalar.
//!
//! The unrolling into eight vectors and the prefetch distances are kept as
//! the reference sets them for the contiguous case.

const std = @import("std");
const types = @import("../../../types.zig");
const util = @import("../../../util.zig");
const fma = @import("../../../fma.zig");
const kernel = @import("../kernel.zig");

const BlasLong = types.BlasLong;

/// `y := alpha*x + y`, the reference's kernel ABI.
pub fn Real(comptime T: type) type {
    // a 128-bit vector
    const W = 16 / @sizeOf(T);
    const V = @Vector(W, T);
    return struct {
        inline fn load(p: [*]const T) V {
            return p[0..W].*;
        }
        inline fn store(p: [*]T, v: V) void {
            p[0..W].* = v;
        }
        /// `y + da * x`, contracted as the build's policy says.
        inline fn fmadd(x: V, da: V, y: V) V {
            return kernel.fmadd(fma.fuse, V, x, da, y);
        }

        pub fn call(
            n: BlasLong,
            _: BlasLong,
            _: BlasLong,
            da: T,
            x: [*]T,
            inc_x: BlasLong,
            y: [*]T,
            inc_y: BlasLong,
            _: [*]T,
            _: BlasLong,
        ) callconv(.c) i32 {
            if (n < 0 or da == 0.0) return 0;
            const len: usize = @intCast(n);
            const dv: V = @splat(da);

            if (inc_x == 1 and inc_y == 1) {
                var i: usize = 0;
                // eight vectors at a time, the prefetch a line ahead of the block
                const block = 8 * W;
                while (i + block <= len) : (i += block) {
                    @prefetch(x + i + 32 * 4 / @sizeOf(T), .{});
                    @prefetch(y + i + 32 * 4 / @sizeOf(T), .{ .rw = .write });
                    inline for (0..8) |k| {
                        const at = i + k * W;
                        store(y + at, fmadd(load(x + at), dv, load(y + at)));
                    }
                }
                // then four, two and one vector, as the tail bits of `n` say
                inline for ([_]usize{ 4, 2, 1 }) |vecs| {
                    if (len & (vecs * W) != 0 and len - i >= vecs * W) {
                        inline for (0..vecs) |k| {
                            const at = i + k * W;
                            store(y + at, fmadd(load(x + at), dv, load(y + at)));
                        }
                        i += vecs * W;
                    }
                }
                while (i < len) : (i += 1) y[i] = kernel.fmadd(fma.fuse, T, da, x[i], y[i]);
                return 0;
            }

            if (inc_y == 1) {
                var i: usize = 0;
                var ix: BlasLong = 0;
                // `y` is unit: gather `x` a vector at a time, four vectors
                // while they last and one after
                const group = 4 * W;
                while (i + group <= len) : (i += group) {
                    inline for (0..4) |k| {
                        var xv: V = undefined;
                        inline for (0..W) |lane| {
                            xv[lane] = util.at(T, x, ix + @as(BlasLong, @intCast(lane)) * inc_x).*;
                        }
                        ix += @as(BlasLong, W) * inc_x;
                        const at = i + k * W;
                        store(y + at, fmadd(xv, dv, load(y + at)));
                    }
                }
                while (i + W <= len) : (i += W) {
                    var xv: V = undefined;
                    inline for (0..W) |lane| {
                        xv[lane] = util.at(T, x, ix + @as(BlasLong, @intCast(lane)) * inc_x).*;
                    }
                    ix += @as(BlasLong, W) * inc_x;
                    store(y + i, fmadd(xv, dv, load(y + i)));
                }
                while (i < len) : (i += 1) {
                    y[i] = kernel.fmadd(fma.fuse, T, da, util.at(T, x, ix).*, y[i]);
                    ix += inc_x;
                }
                return 0;
            }

            // both strided: scalars, four at a time as the reference does
            var ix: BlasLong = 0;
            var iy: BlasLong = 0;
            var i: usize = 0;
            while (i + 4 <= len) : (i += 4) {
                var xs: [4]T = undefined;
                var ys: [4]T = undefined;
                inline for (0..4) |k| {
                    xs[k] = util.at(T, x, ix).*;
                    ix += inc_x;
                }
                inline for (0..4) |k| {
                    ys[k] = util.at(T, y, iy + @as(BlasLong, k) * inc_y).*;
                }
                inline for (0..4) |k| ys[k] = kernel.fmadd(fma.fuse, T, da, xs[k], ys[k]);
                inline for (0..4) |k| {
                    util.at(T, y, iy).* = ys[k];
                    iy += inc_y;
                }
            }
            while (i < len) : (i += 1) {
                util.at(T, y, iy).* = kernel.fmadd(fma.fuse, T, da, util.at(T, x, ix).*, util.at(T, y, iy).*);
                ix += inc_x;
                iy += inc_y;
            }
            return 0;
        }
    };
}

/// The complex kernel over interleaved pairs; `conj` selects the `axpyc`
/// formula, which the reference writes as the signs of three `+=`s.
pub fn Cplx(comptime F: type, comptime conj: bool) type {
    const W = 16 / @sizeOf(F);
    const V = @Vector(W, F);
    return struct {
        inline fn load(p: [*]const F) V {
            return p[0..W].*;
        }
        inline fn store(p: [*]F, v: V) void {
            p[0..W].* = v;
        }

        /// The even and the odd lanes of two vectors, the second vector's first: the
        /// packs (`pckev`/`pckod`, in their argument order) that split
        /// interleaved complex into parts.
        inline fn evens(a: V, b: V) V {
            comptime var mask: [W]i32 = undefined;
            comptime for (0..W) |i| {
                mask[i] = if (i < W / 2) @intCast(2 * i) else ~@as(i32, @intCast(2 * (i - W / 2)));
            };
            return @shuffle(F, b, a, mask);
        }
        inline fn odds(a: V, b: V) V {
            comptime var mask: [W]i32 = undefined;
            comptime for (0..W) |i| {
                mask[i] = if (i < W / 2) @intCast(2 * i + 1) else ~@as(i32, @intCast(2 * (i - W / 2) + 1));
            };
            return @shuffle(F, b, a, mask);
        }
        /// The interleave of the low and of the high halves (`ilvr`/`ilvl`) of
        /// real and imaginary parts, back into pairs.
        inline fn zipLow(re: V, im: V) V {
            comptime var mask: [W]i32 = undefined;
            comptime for (0..W) |i| {
                mask[i] = if (i % 2 == 0) @intCast(i / 2) else ~@as(i32, @intCast(i / 2));
            };
            return @shuffle(F, re, im, mask);
        }
        inline fn zipHigh(re: V, im: V) V {
            comptime var mask: [W]i32 = undefined;
            comptime for (0..W) |i| {
                mask[i] = if (i % 2 == 0) @intCast(W / 2 + i / 2) else ~@as(i32, @intCast(W / 2 + i / 2));
            };
            return @shuffle(F, re, im, mask);
        }

        /// The reference's `OP0`, `OP1`, `OP2`: the signs of the three
        /// accumulations after the first. Its conjugate form has `-=` for the
        /// third, which makes `y.im -= ai*x.re` where its own generic kernel
        /// and the definition, `alpha * conj(x)`, have `+=`; it is `+=` here.
        const op0: F = if (conj) -1 else 1;
        const op1: F = if (conj) 1 else -1;
        const op2: F = 1;

        /// One vector's worth of complex numbers, split into parts: the
        /// reference's four statements, in its order.
        inline fn madd(comptime fuse: bool, dar: V, dai: V, xr: V, xi: V, yr: *V, yi: *V) void {
            yr.* = kernel.fmadd(fuse, V, xr, dar, yr.*);
            yi.* = kernel.fmadd(fuse, V, xi, dar * @as(V, @splat(op0)), yi.*);
            yr.* = kernel.fmadd(fuse, V, xi, dai * @as(V, @splat(op1)), yr.*);
            yi.* = kernel.fmadd(fuse, V, xr, dai * @as(V, @splat(op2)), yi.*);
        }

        /// The same for one number, the scalar tail.
        inline fn maddScalar(da_r: F, da_i: F, xr: F, xi: F, yr: *F, yi: *F) void {
            yr.* = kernel.fmadd(fma.fuse, F, da_r, xr, yr.*);
            yi.* = kernel.fmadd(fma.fuse, F, da_r * op0, xi, yi.*);
            yr.* = kernel.fmadd(fma.fuse, F, da_i * op1, xi, yr.*);
            yi.* = kernel.fmadd(fma.fuse, F, da_i * op2, xr, yi.*);
        }

        pub fn call(
            n: BlasLong,
            _: BlasLong,
            _: BlasLong,
            da_r: F,
            da_i: F,
            x: [*]F,
            inc_x: BlasLong,
            y: [*]F,
            inc_y: BlasLong,
            _: [*]F,
            _: BlasLong,
        ) callconv(.c) i32 {
            if (n < 0) return 0;
            if (da_r == 0.0 and da_i == 0.0) return 0;
            const len: usize = @intCast(n);
            const dar: V = @splat(da_r);
            const dai: V = @splat(da_i);
            // complex numbers a vector holds
            const per = W / 2;
            // complex numbers per pair of vectors, which is what `madd` takes
            const chunk = 2 * per;

            if (inc_x == 1 and inc_y == 1) {
                var i: usize = 0;
                // four chunks at a time, then one, then the scalars
                while (i + 4 * chunk <= len) : (i += 4 * chunk) {
                    @prefetch(x + 2 * i + 64, .{});
                    @prefetch(y + 2 * i + 64, .{ .rw = .write });
                    inline for (0..4) |k| {
                        const at = 2 * (i + k * chunk);
                        const x0 = load(x + at);
                        const x1 = load(x + at + W);
                        var y0 = load(y + at);
                        var y1 = load(y + at + W);
                        var yr = evens(y1, y0);
                        var yi = odds(y1, y0);
                        madd(fma.fuse, dar, dai, evens(x1, x0), odds(x1, x0), &yr, &yi);
                        y0 = zipLow(yr, yi);
                        y1 = zipHigh(yr, yi);
                        store(y + at, y0);
                        store(y + at + W, y1);
                    }
                }
                while (i + chunk <= len) : (i += chunk) {
                    const at = 2 * i;
                    const x0 = load(x + at);
                    const x1 = load(x + at + W);
                    var yr = evens(load(y + at + W), load(y + at));
                    var yi = odds(load(y + at + W), load(y + at));
                    madd(fma.fuse, dar, dai, evens(x1, x0), odds(x1, x0), &yr, &yi);
                    store(y + at, zipLow(yr, yi));
                    store(y + at + W, zipHigh(yr, yi));
                }
                while (i < len) : (i += 1) {
                    var yr = y[2 * i];
                    var yi = y[2 * i + 1];
                    maddScalar(da_r, da_i, x[2 * i], x[2 * i + 1], &yr, &yi);
                    y[2 * i] = yr;
                    y[2 * i + 1] = yi;
                }
                return 0;
            }

            if (inc_y == 1) {
                // `y` is unit: gather `x` a chunk at a time
                var i: usize = 0;
                var ix: BlasLong = 0;
                while (i + chunk <= len) : (i += chunk) {
                    var xs: [2]V = undefined;
                    inline for (0..2) |v| {
                        inline for (0..per) |c| {
                            xs[v][2 * c] = util.at(F, x, 2 * ix).*;
                            xs[v][2 * c + 1] = util.at(F, x, 2 * ix + 1).*;
                            ix += inc_x;
                        }
                    }
                    const at = 2 * i;
                    var yr = evens(load(y + at + W), load(y + at));
                    var yi = odds(load(y + at + W), load(y + at));
                    madd(fma.fuse, dar, dai, evens(xs[1], xs[0]), odds(xs[1], xs[0]), &yr, &yi);
                    store(y + at, zipLow(yr, yi));
                    store(y + at + W, zipHigh(yr, yi));
                }
                while (i < len) : (i += 1) {
                    var yr = y[2 * i];
                    var yi = y[2 * i + 1];
                    maddScalar(da_r, da_i, util.at(F, x, 2 * ix).*, util.at(F, x, 2 * ix + 1).*, &yr, &yi);
                    y[2 * i] = yr;
                    y[2 * i + 1] = yi;
                    ix += inc_x;
                }
                return 0;
            }

            // everything else: scalars
            var ix: BlasLong = 0;
            var iy: BlasLong = 0;
            var i: usize = 0;
            while (i < len) : (i += 1) {
                var yr = util.at(F, y, 2 * iy).*;
                var yi = util.at(F, y, 2 * iy + 1).*;
                maddScalar(da_r, da_i, util.at(F, x, 2 * ix).*, util.at(F, x, 2 * ix + 1).*, &yr, &yi);
                util.at(F, y, 2 * iy).* = yr;
                util.at(F, y, 2 * iy + 1).* = yi;
                ix += inc_x;
                iy += inc_y;
            }
            return 0;
        }
    };
}
