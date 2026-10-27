//! `@Vector` kernels for level-1 and level-2 slots, in the kernel ABI: plain
//! Zig, `U` native vectors a step. At `U = 4` they match or beat OpenBLAS's
//! arm64 kernels on Apple's cores (`tools/blasbench/results-m6.md`). Unit
//! strides take the vectors; other strides, zero and negative ones too, a
//! scalar loop. NaN, infinities, `scal`'s flag, `i?amax`'s first index and
//! which strides a kernel refuses are the portable kernels'; `nrm2` sums
//! plain squares, in `f64` for `f32`, and for `f64` takes Blue's algorithm
//! where they could have overflowed or underflowed. `Registry` exports them
//! from the per-core objects for the cores and slots an architecture's
//! registry picks (`kernels/<arch>.zig`), among them every `routed` slot.
//! Imports stay inside `linalg/blas/`: the per-core objects compile this
//! file on their own.

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("types.zig");
const util = @import("util.zig");
const Vtable = @import("Vtable.zig");
const fma = @import("fma.zig");
const extremum = @import("blas1/extremum.zig");
const coretypes = @import("cpu/_coretypes.zig");

const BlasLong = types.BlasLong;
const Find = extremum.Find;
const Slot = util.Slot;

/// Lanes of one native vector of `F`: a NEON register on AArch64, where
/// `std.simd` counts SVE as 256 bits, which LLVM splits into NEON pairs.
pub fn lanes(comptime F: type) comptime_int {
    if (comptime sysinfo.target.cpu.arch.isAarch64()) return 16 / @sizeOf(F);
    return std.simd.suggestVectorLength(F) orelse 1;
}

inline fn at(comptime F: type, p: [*]F, i: BlasLong) *F {
    return util.at(F, p, i);
}

/// `a*b + c`, fused where the hardware has an FMA unit (and `fma.compiles`);
/// elsewhere rounded twice, as `@mulAdd` there is a software `fma` call per
/// element.
inline fn mulAdd(comptime T: type, a: T, b: T, c: T) T {
    return if (comptime fma.hardware and fma.compiles) @mulAdd(T, a, b, c) else a * b + c;
}

/// The real kernels over `F`; reductions keep two vectors of sums.
pub fn Real(comptime F: type, comptime U: comptime_int) type {
    return struct {
        pub const W = U * lanes(F);
        const V = @Vector(W, F);
        /// Four strided elements, for four sums.
        const G = @Vector(4, F);

        inline fn ld(p: [*]const F) V {
            return p[0..W].*;
        }
        inline fn st(p: [*]F, v: V) void {
            p[0..W].* = v;
        }
        inline fn gather(p: [*]F, inc: BlasLong, i: BlasLong) G {
            return .{ at(F, p, i * inc).*, at(F, p, (i + 1) * inc).*, at(F, p, (i + 2) * inc).*, at(F, p, (i + 3) * inc).* };
        }

        pub fn axpy(n: BlasLong, _: BlasLong, _: BlasLong, a: F, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong, _: [*]F, _: BlasLong) callconv(.c) i32 {
            if (n <= 0) return 0;
            if (incx == 1 and incy == 1) {
                const len: usize = @intCast(n);
                const av: V = @splat(a);
                var i: usize = 0;
                while (i + W <= len) : (i += W) st(y + i, mulAdd(V, av, ld(x + i), ld(y + i)));
                while (i < len) : (i += 1) y[i] = mulAdd(F, a, x[i], y[i]);
                return 0;
            }
            var i: BlasLong = 0;
            while (i < n) : (i += 1) {
                const yp = at(F, y, i * incy);
                yp.* = mulAdd(F, a, at(F, x, i * incx).*, yp.*);
            }
            return 0;
        }

        /// `dot`'s and `nrm2`'s sums: `f64` for `f32`, as the portable kernels
        /// (products of `f32` are exact there), else `F`.
        pub const D = if (F == f32) f64 else F;
        const DV = @Vector(W, D);
        const DG = @Vector(4, D);

        pub fn dot(n: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong) callconv(.c) F {
            if (n <= 0) return 0;
            if (incx == 1 and incy == 1) {
                const len: usize = @intCast(n);
                var acc: [2]DV = @splat(@splat(0));
                var i: usize = 0;
                while (i + 2 * W <= len) : (i += 2 * W) {
                    inline for (0..2) |q| acc[q] = mulAdd(DV, ld(x + i + q * W), ld(y + i + q * W), acc[q]);
                }
                if (i + W <= len) {
                    acc[0] = mulAdd(DV, ld(x + i), ld(y + i), acc[0]);
                    i += W;
                }
                var s = @reduce(.Add, acc[0] + acc[1]);
                while (i < len) : (i += 1) s = mulAdd(D, x[i], y[i], s);
                return @floatCast(s);
            }
            var acc: DG = @splat(0);
            var i: BlasLong = 0;
            while (i + 4 <= n) : (i += 4) acc = mulAdd(DG, gather(x, incx, i), gather(y, incy, i), acc);
            var s = @reduce(.Add, acc);
            while (i < n) : (i += 1) s = mulAdd(D, at(F, x, i * incx).*, at(F, y, i * incy).*, s);
            return @floatCast(s);
        }

        const Scale = enum { mul, zero, keep };

        /// `a*v`; for a zero `a`, +0, or NaN for NaN and infinity (`keep`).
        inline fn scaled(comptime T: type, comptime how: Scale, a: F, v: T) T {
            if (comptime @typeInfo(T) == .vector) return switch (how) {
                .mul => @as(T, @splat(a)) * v,
                .zero => @splat(0),
                .keep => @select(F, @abs(v) <= @as(T, @splat(std.math.floatMax(F))), @as(T, @splat(0)), @as(T, @splat(std.math.nan(F)))),
            };
            return switch (how) {
                .mul => a * v,
                .zero => 0,
                .keep => if (std.math.isFinite(v)) 0 else std.math.nan(F),
            };
        }

        fn scalBy(comptime how: Scale, a: F, len: usize, x: [*]F, sx: usize) void {
            var i: usize = 0;
            if (sx == 1) {
                while (i + W <= len) : (i += W) st(x + i, scaled(V, how, a, ld(x + i)));
            }
            while (i < len) : (i += 1) x[i * sx] = scaled(F, how, a, x[i * sx]);
        }

        /// A zero `alpha` stores +0, or with `flag` 1 NaN for NaN and infinity.
        pub fn scal(n: BlasLong, _: BlasLong, _: BlasLong, a: F, x: [*]F, incx: BlasLong, _: [*]F, _: BlasLong, _: [*]F, flag: BlasLong) callconv(.c) i32 {
            if (n <= 0 or incx <= 0) return 0;
            const how: Scale = if (a != 0) .mul else if (flag == 0) .zero else .keep;
            switch (how) {
                inline else => |h| scalBy(h, a, @intCast(n), x, @intCast(incx)),
            }
            return 0;
        }

        /// The plain sum of squares of `x[0..len]`, in `D`.
        pub fn squares(len: usize, x: [*]F) D {
            var acc: [2]DV = @splat(@splat(0));
            var i: usize = 0;
            while (i + 2 * W <= len) : (i += 2 * W) {
                inline for (0..2) |q| {
                    const v: DV = ld(x + i + q * W);
                    acc[q] = mulAdd(DV, v, v, acc[q]);
                }
            }
            var s = @reduce(.Add, acc[0] + acc[1]);
            while (i < len) : (i += 1) s = mulAdd(D, x[i], x[i], s);
            return s;
        }

        pub fn nrm2(n: BlasLong, x: [*]F, incx: BlasLong) callconv(.c) F {
            if (n <= 0 or incx == 0) return 0;
            if (n == 1) return @abs(x[0]);
            if (incx == 1) return root(F, 1, squares(@intCast(n), x), n, x, incx);
            var acc: DG = @splat(0);
            var i: BlasLong = 0;
            while (i + 4 <= n) : (i += 4) {
                const v: DG = gather(x, incx, i);
                acc = mulAdd(DG, v, v, acc);
            }
            var s = @reduce(.Add, acc);
            while (i < n) : (i += 1) {
                const v: D = at(F, x, i * incx).*;
                s = mulAdd(D, v, v, s);
            }
            return root(F, 1, s, n, x, incx);
        }

        /// The sum of `|x[0..len]|`.
        pub fn absSum(len: usize, x: [*]F) F {
            var acc: [2]V = @splat(@splat(0));
            var i: usize = 0;
            while (i + 2 * W <= len) : (i += 2 * W) {
                inline for (0..2) |q| acc[q] += @abs(ld(x + i + q * W));
            }
            var s = @reduce(.Add, acc[0] + acc[1]);
            while (i < len) : (i += 1) s += @abs(x[i]);
            return s;
        }

        pub fn asum(n: BlasLong, x: [*]F, incx: BlasLong) callconv(.c) F {
            if (n <= 0 or incx <= 0) return 0;
            if (incx == 1) return absSum(@intCast(n), x);
            var acc: G = @splat(0);
            var i: BlasLong = 0;
            while (i + 4 <= n) : (i += 4) acc += @abs(gather(x, incx, i));
            var s = @reduce(.Add, acc);
            while (i < n) : (i += 1) s += @abs(at(F, x, i * incx).*);
            return s;
        }

        /// The extremes by `find`, value and 1-based index, as the portable
        /// kernels (`extremum.zig`), to which they fall back for a NaN.
        fn Extreme(comptime find: Find) type {
            return struct {
                const P = extremum.Extremum(F, 1, find);
                const Key = struct {
                    inline fn key(p: [*]const F, s: usize, i: usize) F {
                        return if (find.abs) @abs(p[i * s]) else p[i * s];
                    }
                    inline fn block(p: [*]const F, s: usize, i: usize, chk: *V) F {
                        var m: [4]V = undefined;
                        inline for (0..4) |q| {
                            if (s == 1) {
                                m[q] = ld(p + i + q * W);
                            } else {
                                var r: [W]F = undefined;
                                inline for (0..W) |k| r[k] = p[(i + q * W + k) * s];
                                m[q] = r;
                            }
                            if (find.abs) m[q] = @abs(m[q]);
                            chk.* += if (find.abs) m[q] else @abs(m[q]);
                        }
                        return pick(find, m);
                    }
                };

                pub fn value(n: BlasLong, x: [*]F, incx: BlasLong) callconv(.c) F {
                    if (n <= 0 or incx <= 0) return 0;
                    const r = search(@intCast(n), x, @intCast(incx), find, Key) orelse return P.value(n, x, incx);
                    return r[1];
                }

                pub fn index(n: BlasLong, x: [*]F, incx: BlasLong) callconv(.c) BlasLong {
                    if (n <= 0 or incx <= 0) return 0;
                    const r = search(@intCast(n), x, @intCast(incx), find, Key) orelse return P.index(n, x, incx);
                    return @intCast(r[0] + 1);
                }
            };
        }

        pub const amax = Extreme(.{ .abs = true, .max = true }).value;
        pub const iamax = Extreme(.{ .abs = true, .max = true }).index;
        pub const amin = Extreme(.{ .abs = true, .max = false }).value;
        pub const iamin = Extreme(.{ .abs = true, .max = false }).index;
        pub const max = Extreme(.{ .abs = false, .max = true }).value;
        pub const imax = Extreme(.{ .abs = false, .max = true }).index;
        pub const min = Extreme(.{ .abs = false, .max = false }).value;
        pub const imin = Extreme(.{ .abs = false, .max = false }).index;

        /// The extreme of four vectors of keys.
        pub inline fn pick(comptime find: Find, m: [4]V) F {
            if (find.max) return @reduce(.Max, @max(@max(m[0], m[1]), @max(m[2], m[3])));
            return @reduce(.Min, @min(@min(m[0], m[1]), @min(m[2], m[3])));
        }

        inline fn better(comptime find: Find, k: F, best: F) bool {
            return if (find.max) k > best else k < best;
        }

        /// `find`'s first extreme of `len` keys at a stride of `s` floats,
        /// `Key.key` giving one and `Key.block` the extreme of the `4 * W`
        /// from `i`, adding their magnitudes to `chk`, which only a NaN key
        /// makes NaN. One pass keeps the first block with the extreme, no
        /// index kept per element; only that block is searched for the
        /// index. Its index from 0 and its key, or null for a NaN key.
        pub fn search(len: usize, x: [*]const F, s: usize, comptime find: Find, comptime Key: type) ?struct { usize, F } {
            const B = 4 * W;
            var best = Key.key(x, s, 0);
            if (best != best) return null;
            var pos: usize = 0;
            var i: usize = 0;
            if (len >= B) {
                var chk: V = @splat(0);
                var blk: usize = 0;
                while (i + B <= len) : (i += B) {
                    const m = Key.block(x, s, i, &chk);
                    if (better(find, m, best)) {
                        best = m;
                        blk = i;
                    }
                }
                if (std.math.isNan(@reduce(.Add, chk))) return null;
                pos = blk;
                while (Key.key(x, s, pos) != best) pos += 1;
            }
            while (i < len) : (i += 1) {
                const k = Key.key(x, s, i);
                if (k != k) return null;
                if (better(find, k, best)) {
                    best = k;
                    pos = i;
                }
            }
            return .{ pos, Key.key(x, s, pos) };
        }

        pub fn copy(n: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong) callconv(.c) i32 {
            if (n <= 0) return 0;
            if (incx == 1 and incy == 1) {
                const len: usize = @intCast(n);
                var i: usize = 0;
                while (i + W <= len) : (i += W) st(y + i, ld(x + i));
                while (i < len) : (i += 1) y[i] = x[i];
                return 0;
            }
            var i: BlasLong = 0;
            while (i < n) : (i += 1) at(F, y, i * incy).* = at(F, x, i * incx).*;
            return 0;
        }

        pub fn swap(n: BlasLong, _: BlasLong, _: BlasLong, _: F, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong, _: [*]F, _: BlasLong) callconv(.c) i32 {
            if (n <= 0) return 0;
            if (incx == 1 and incy == 1) {
                const len: usize = @intCast(n);
                var i: usize = 0;
                while (i + W <= len) : (i += W) {
                    const t = ld(x + i);
                    st(x + i, ld(y + i));
                    st(y + i, t);
                }
                while (i < len) : (i += 1) std.mem.swap(F, &x[i], &y[i]);
                return 0;
            }
            var i: BlasLong = 0;
            while (i < n) : (i += 1) std.mem.swap(F, at(F, x, i * incx), at(F, y, i * incy));
            return 0;
        }

        /// `x := c*x + s*y`, `y := c*y - s*x`.
        pub fn rot(n: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong, c: F, s: F) callconv(.c) i32 {
            if (n <= 0) return 0;
            if (incx == 1 and incy == 1) {
                const len: usize = @intCast(n);
                const cv: V = @splat(c);
                const sv: V = @splat(s);
                var i: usize = 0;
                while (i + W <= len) : (i += W) {
                    const xv = ld(x + i);
                    const yv = ld(y + i);
                    st(x + i, mulAdd(V, cv, xv, sv * yv));
                    st(y + i, mulAdd(V, cv, yv, -(sv * xv)));
                }
                while (i < len) : (i += 1) rot1(&x[i], &y[i], c, s);
                return 0;
            }
            var i: BlasLong = 0;
            while (i < n) : (i += 1) rot1(at(F, x, i * incx), at(F, y, i * incy), c, s);
            return 0;
        }

        inline fn rot1(xp: *F, yp: *F, c: F, s: F) void {
            const xv = xp.*;
            const yv = yp.*;
            xp.* = mulAdd(F, c, xv, s * yv);
            yp.* = mulAdd(F, c, yv, -(s * xv));
        }

        /// `y += alpha*A*x`, four columns a pass over a unit-stride `y`.
        pub fn gemv_n(m: BlasLong, n: BlasLong, _: BlasLong, alpha: F, a: [*]F, lda: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong, _: [*]F) callconv(.c) i32 {
            if (m <= 0 or n <= 0) return 0;
            const rows: usize = @intCast(m);
            const cols: usize = @intCast(n);
            const ld_a: usize = @intCast(lda);
            var j: usize = 0;
            if (incy == 1) {
                while (j + 4 <= cols) : (j += 4) {
                    var t: [4]F = undefined;
                    var c: [4][*]F = undefined;
                    inline for (0..4) |q| {
                        t[q] = alpha * at(F, x, @as(BlasLong, @intCast(j + q)) * incx).*;
                        c[q] = a + (j + q) * ld_a;
                    }
                    var i: usize = 0;
                    while (i + W <= rows) : (i += W) {
                        var yv = ld(y + i);
                        inline for (0..4) |q| yv = mulAdd(V, ld(c[q] + i), @splat(t[q]), yv);
                        st(y + i, yv);
                    }
                    while (i < rows) : (i += 1) {
                        inline for (0..4) |q| y[i] = mulAdd(F, c[q][i], t[q], y[i]);
                    }
                }
            }
            while (j < cols) : (j += 1) {
                const t = alpha * at(F, x, @as(BlasLong, @intCast(j)) * incx).*;
                const c = a + j * ld_a;
                var i: usize = 0;
                if (incy == 1) {
                    while (i + W <= rows) : (i += W) st(y + i, mulAdd(V, ld(c + i), @splat(t), ld(y + i)));
                }
                while (i < rows) : (i += 1) {
                    const yp = at(F, y, @as(BlasLong, @intCast(i)) * incy);
                    yp.* = mulAdd(F, c[i], t, yp.*);
                }
            }
            return 0;
        }

        /// `y += alpha*A'*x`, four columns' sums a pass over a unit-stride
        /// `x`. With no rows, `y += alpha*0`, as the portable kernel.
        pub fn gemv_t(m: BlasLong, n: BlasLong, _: BlasLong, alpha: F, a: [*]F, lda: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong, _: [*]F) callconv(.c) i32 {
            if (n <= 0) return 0;
            const rows: usize = @intCast(@max(m, 0));
            const cols: usize = @intCast(n);
            const ld_a: usize = @intCast(lda);
            var j: usize = 0;
            if (incx == 1) {
                while (j + 4 <= cols) : (j += 4) {
                    var c: [4][*]F = undefined;
                    var acc: [4]V = @splat(@splat(0));
                    inline for (0..4) |q| c[q] = a + (j + q) * ld_a;
                    var i: usize = 0;
                    while (i + W <= rows) : (i += W) {
                        const xv = ld(x + i);
                        inline for (0..4) |q| acc[q] = mulAdd(V, ld(c[q] + i), xv, acc[q]);
                    }
                    inline for (0..4) |q| {
                        var s = @reduce(.Add, acc[q]);
                        var r = i;
                        while (r < rows) : (r += 1) s = mulAdd(F, c[q][r], x[r], s);
                        const yp = at(F, y, @as(BlasLong, @intCast(j + q)) * incy);
                        yp.* = mulAdd(F, alpha, s, yp.*);
                    }
                }
            }
            while (j < cols) : (j += 1) {
                const c = a + j * ld_a;
                var s: F = 0;
                var i: usize = 0;
                if (incx == 1) {
                    var acc: V = @splat(0);
                    while (i + W <= rows) : (i += W) acc = mulAdd(V, ld(c + i), ld(x + i), acc);
                    s = @reduce(.Add, acc);
                }
                while (i < rows) : (i += 1) s = mulAdd(F, c[i], at(F, x, @as(BlasLong, @intCast(i)) * incx).*, s);
                const yp = at(F, y, @as(BlasLong, @intCast(j)) * incy);
                yp.* = mulAdd(F, alpha, s, yp.*);
            }
            return 0;
        }
    };
}

/// `nrm2` from `s`, the plain sum of squares of `n` elements of `w` floats.
/// Summed in `f64` for `f32`, it is exact enough as it is. Else `@sqrt(s)`
/// where it neither overflowed nor could have lost more than `eps^2` of
/// itself to underflow, and Blue's algorithm (`util.blue`) where it could have.
fn root(comptime F: type, comptime w: comptime_int, s: anytype, n: BlasLong, x: [*]F, incx: BlasLong) F {
    if (comptime @TypeOf(s) != F) return @floatCast(@sqrt(s));
    const lo = std.math.floatMin(F) / std.math.floatEps(F) * @as(F, @floatFromInt(w * n));
    if (s >= lo and s <= std.math.floatMax(F)) return @sqrt(s);
    return util.blue(F, w, n, x, incx);
}

/// The complex kernels over interleaved `F` pairs.
pub fn Complex(comptime F: type, comptime U: comptime_int) type {
    return struct {
        const R = Real(F, U);
        /// Floats a vector; whole elements.
        pub const W = if (R.W % 2 == 0) R.W else 2 * R.W;
        const V = @Vector(W, F);
        /// One element.
        const V2 = @Vector(2, F);
        const Result = types.CComplex(F);

        inline fn ld(p: [*]const F) V {
            return p[0..W].*;
        }
        inline fn st(p: [*]F, v: V) void {
            p[0..W].* = v;
        }
        inline fn ld2(p: [*]F, i: BlasLong) V2 {
            return @as(*[2]F, @ptrCast(at(F, p, 2 * i))).*;
        }
        inline fn st2(p: [*]F, i: BlasLong, v: V2) void {
            @as(*[2]F, @ptrCast(at(F, p, 2 * i))).* = v;
        }

        /// `(re, im, re, im, ...)` as a `T`.
        inline fn pair(comptime T: type, re: F, im: F) T {
            var r: [@typeInfo(T).vector.len]F = undefined;
            inline for (0..r.len) |k| r[k] = if (k % 2 == 0) re else im;
            return r;
        }

        /// Each element's parts exchanged.
        inline fn swapped(comptime T: type, v: T) T {
            const len = @typeInfo(T).vector.len;
            const mask = comptime blk: {
                var m: [len]i32 = undefined;
                for (0..len) |k| m[k] = @intCast(k ^ 1);
                break :blk m;
            };
            return @shuffle(F, v, undefined, mask);
        }

        /// Multiplying by `(tr, ti)`, the operand conjugated for `conj`: lane
        /// factors on it and on it with its parts exchanged.
        fn By(comptime T: type) type {
            return struct {
                r: T,
                i: T,

                inline fn init(tr: F, ti: F, comptime conj: bool) @This() {
                    return .{ .r = pair(T, tr, if (conj) -tr else tr), .i = pair(T, if (conj) ti else -ti, ti) };
                }
                /// `y + t*v`.
                inline fn madd(b: @This(), v: T, y: T) T {
                    return mulAdd(T, b.r, v, mulAdd(T, b.i, swapped(T, v), y));
                }
                /// `t*v`.
                inline fn mul(b: @This(), v: T) T {
                    return mulAdd(T, b.i, swapped(T, v), b.r * v);
                }
            };
        }

        /// The sums of `a*x` lane by lane and of `a*x` with `x`'s parts
        /// exchanged, as `Dots.finish` makes `sum(op(a)*op(x))` of them.
        fn Dots(comptime ca: bool, comptime cx: bool) type {
            return struct {
                rr: V = @splat(0),
                ri: V = @splat(0),
                rr2: V2 = @splat(0),
                ri2: V2 = @splat(0),

                inline fn add(d: *@This(), av: V, xv: V) void {
                    d.rr = mulAdd(V, av, xv, d.rr);
                    d.ri = mulAdd(V, av, swapped(V, xv), d.ri);
                }
                inline fn add2(d: *@This(), av: V2, xv: V2) void {
                    d.rr2 = mulAdd(V2, av, xv, d.rr2);
                    d.ri2 = mulAdd(V2, av, swapped(V2, xv), d.ri2);
                }
                /// `re = ar*xr - sa*sx*ai*xi`, `im = sx*ar*xi + sa*ai*xr`,
                /// `sa` and `sx` the signs conjugation gives the imaginary parts.
                fn finish(d: @This()) V2 {
                    const sa: F = if (ca) -1 else 1;
                    const sx: F = if (cx) -1 else 1;
                    const re = @reduce(.Add, d.rr * pair(V, 1, -sa * sx)) + @reduce(.Add, d.rr2 * pair(V2, 1, -sa * sx));
                    const im = @reduce(.Add, d.ri * pair(V, sx, sa)) + @reduce(.Add, d.ri2 * pair(V2, sx, sa));
                    return .{ re, im };
                }
            };
        }

        fn Axpy(comptime conj: bool) type {
            return struct {
                fn k(n: BlasLong, _: BlasLong, _: BlasLong, ar: F, ai: F, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong, _: [*]F, _: BlasLong) callconv(.c) i32 {
                    if (n <= 0) return 0;
                    var i: BlasLong = 0;
                    if (incx == 1 and incy == 1) {
                        const by = By(V).init(ar, ai, conj);
                        const len: usize = @intCast(2 * n);
                        var f: usize = 0;
                        while (f + W <= len) : (f += W) st(y + f, by.madd(ld(x + f), ld(y + f)));
                        i = @intCast(f / 2);
                    }
                    const by = By(V2).init(ar, ai, conj);
                    while (i < n) : (i += 1) st2(y, i * incy, by.madd(ld2(x, i * incx), ld2(y, i * incy)));
                    return 0;
                }
            };
        }

        /// `y += alpha*x`.
        pub const axpy = Axpy(false).k;
        /// `y += alpha*conj(x)`.
        pub const axpyc = Axpy(true).k;

        /// `sum(op(x)*y)`, `x` conjugated for `conj`.
        fn dot(comptime conj: bool, n: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong) Result {
            var d: [2]Dots(conj, false) = .{ .{}, .{} };
            var i: BlasLong = 0;
            if (n > 0 and incx == 1 and incy == 1) {
                const len: usize = @intCast(2 * n);
                var f: usize = 0;
                while (f + 2 * W <= len) : (f += 2 * W) {
                    inline for (0..2) |q| d[q].add(ld(x + f + q * W), ld(y + f + q * W));
                }
                i = @intCast(f / 2);
            }
            while (i + 2 <= n) : (i += 2) {
                inline for (0..2) |q| d[q].add2(ld2(x, (i + q) * incx), ld2(y, (i + q) * incy));
            }
            if (i < n) d[0].add2(ld2(x, i * incx), ld2(y, i * incy));
            const s = d[0].finish() + d[1].finish();
            return .{ .re = s[0], .im = s[1] };
        }

        pub fn dotu(n: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong) callconv(.c) Result {
            return dot(false, n, x, incx, y, incy);
        }
        pub fn dotc(n: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong) callconv(.c) Result {
            return dot(true, n, x, incx, y, incy);
        }

        /// A zero `alpha` with `flag` 0 stores zeros; any other multiplies.
        pub fn scal(n: BlasLong, _: BlasLong, _: BlasLong, ar: F, ai: F, x: [*]F, incx: BlasLong, _: [*]F, _: BlasLong, _: [*]F, flag: BlasLong) callconv(.c) i32 {
            if (n <= 0 or incx <= 0) return 0;
            if (flag == 0 and ar == 0 and ai == 0) {
                var i: BlasLong = 0;
                while (i < n) : (i += 1) st2(x, i * incx, @splat(0));
                return 0;
            }
            var i: BlasLong = 0;
            if (incx == 1) {
                const by = By(V).init(ar, ai, false);
                const len: usize = @intCast(2 * n);
                var f: usize = 0;
                while (f + W <= len) : (f += W) st(x + f, by.mul(ld(x + f)));
                i = @intCast(f / 2);
            }
            const by = By(V2).init(ar, ai, false);
            while (i < n) : (i += 1) st2(x, i * incx, by.mul(ld2(x, i * incx)));
            return 0;
        }

        pub fn nrm2(n: BlasLong, x: [*]F, incx: BlasLong) callconv(.c) F {
            if (n <= 0 or incx == 0) return 0;
            if (incx == 1) return root(F, 2, R.squares(@intCast(2 * n), x), n, x, incx);
            const D2 = @Vector(2, R.D);
            var acc: [2]D2 = @splat(@splat(0));
            var i: BlasLong = 0;
            while (i + 2 <= n) : (i += 2) {
                inline for (0..2) |q| {
                    const v: D2 = ld2(x, (i + q) * incx);
                    acc[q] = mulAdd(D2, v, v, acc[q]);
                }
            }
            if (i < n) {
                const v: D2 = ld2(x, i * incx);
                acc[0] = mulAdd(D2, v, v, acc[0]);
            }
            return root(F, 2, @reduce(.Add, acc[0] + acc[1]), n, x, incx);
        }

        /// The sum of `|re| + |im|`.
        pub fn asum(n: BlasLong, x: [*]F, incx: BlasLong) callconv(.c) F {
            if (n <= 0 or incx <= 0) return 0;
            if (incx == 1) return R.absSum(@intCast(2 * n), x);
            var acc: [2]V2 = @splat(@splat(0));
            var i: BlasLong = 0;
            while (i + 2 <= n) : (i += 2) {
                inline for (0..2) |q| acc[q] += @abs(ld2(x, (i + q) * incx));
            }
            if (i < n) acc[0] += @abs(ld2(x, i * incx));
            return @reduce(.Add, acc[0] + acc[1]);
        }

        /// The even (`odd` false) or odd floats of `a`, then of `b`.
        inline fn parts(a: V, b: V, comptime odd: bool) V {
            const mask = comptime blk: {
                var m: [W]i32 = undefined;
                for (0..W) |k| {
                    const s = 2 * k + @intFromBool(odd);
                    m[k] = if (s < W) @intCast(s) else ~@as(i32, @intCast(s - W));
                }
                break :blk m;
            };
            return @shuffle(F, a, b, mask);
        }

        /// The extremes of `|re| + |im|`, as `Real`'s. The index falls back
        /// to the portable kernel for an infinite key too, which may be a
        /// finite element's that overflowed.
        fn Extreme(comptime find: Find) type {
            return struct {
                const P = extremum.Extremum(F, 2, find);
                const S = Real(F, W / lanes(F));
                const Key = struct {
                    inline fn key(p: [*]const F, s: usize, i: usize) F {
                        return @abs(p[i * s]) + @abs(p[i * s + 1]);
                    }
                    inline fn block(p: [*]const F, s: usize, i: usize, chk: *V) F {
                        var m: [4]V = undefined;
                        inline for (0..4) |q| {
                            if (s == 2) {
                                const a = @abs(ld(p + 2 * (i + q * W)));
                                const b = @abs(ld(p + 2 * (i + q * W) + W));
                                m[q] = parts(a, b, false) + parts(a, b, true);
                            } else {
                                var r: [W]F = undefined;
                                inline for (0..W) |k| r[k] = key(p, s, i + q * W + k);
                                m[q] = r;
                            }
                            chk.* += m[q];
                        }
                        return S.pick(find, m);
                    }
                };

                fn search(n: BlasLong, x: [*]F, incx: BlasLong) ?struct { usize, F } {
                    return S.search(@intCast(n), x, 2 * @as(usize, @intCast(incx)), find, Key);
                }

                pub fn value(n: BlasLong, x: [*]F, incx: BlasLong) callconv(.c) F {
                    if (n <= 0 or incx <= 0) return 0;
                    const r = search(n, x, incx) orelse return P.value(n, x, incx);
                    return r[1];
                }

                pub fn index(n: BlasLong, x: [*]F, incx: BlasLong) callconv(.c) BlasLong {
                    if (n <= 0 or incx <= 0) return 0;
                    const r = search(n, x, incx) orelse return P.index(n, x, incx);
                    if (std.math.isInf(r[1])) return P.index(n, x, incx);
                    return @intCast(r[0] + 1);
                }
            };
        }

        pub const amax = Extreme(.{ .abs = true, .max = true }).value;
        pub const iamax = Extreme(.{ .abs = true, .max = true }).index;
        pub const amin = Extreme(.{ .abs = true, .max = false }).value;
        pub const iamin = Extreme(.{ .abs = true, .max = false }).index;

        pub fn copy(n: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong) callconv(.c) i32 {
            if (n <= 0) return 0;
            if (incx == 1 and incy == 1) return R.copy(2 * n, x, 1, y, 1);
            var i: BlasLong = 0;
            while (i < n) : (i += 1) st2(y, i * incy, ld2(x, i * incx));
            return 0;
        }

        pub fn swap(n: BlasLong, _: BlasLong, _: BlasLong, _: F, _: F, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong, _: [*]F, _: BlasLong) callconv(.c) i32 {
            if (n <= 0) return 0;
            if (incx == 1 and incy == 1) return R.swap(2 * n, 0, 0, 0, x, 1, y, 1, x, 0);
            var i: BlasLong = 0;
            while (i < n) : (i += 1) {
                const t = ld2(x, i * incx);
                st2(x, i * incx, ld2(y, i * incy));
                st2(y, i * incy, t);
            }
            return 0;
        }

        /// `csrot`, `zdrot`: both parts rotated by the real `c` and `s`.
        pub fn rot(n: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong, c: F, s: F) callconv(.c) i32 {
            if (n <= 0) return 0;
            if (incx == 1 and incy == 1) return R.rot(2 * n, x, 1, y, 1, c, s);
            const cv: V2 = @splat(c);
            const sv: V2 = @splat(s);
            var i: BlasLong = 0;
            while (i < n) : (i += 1) {
                const xv = ld2(x, i * incx);
                const yv = ld2(y, i * incy);
                st2(x, i * incx, mulAdd(V2, cv, xv, sv * yv));
                st2(y, i * incy, mulAdd(V2, cv, yv, -(sv * xv)));
            }
            return 0;
        }

        /// `y += alpha*op(A)*x`: transposed for `trans`, `A` conjugated for
        /// `ca`, `x` for `cx`, as `gemv/`'s portable `Zgemv`. Without a
        /// transpose, four columns a pass over a unit-stride `y`; with one,
        /// four columns' sums a pass over a unit-stride `x`.
        fn Gemv(comptime trans: bool, comptime ca: bool, comptime cx: bool) type {
            return struct {
                fn k(m: BlasLong, n: BlasLong, _: BlasLong, ar: F, ai: F, a: [*]F, lda: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong, _: [*]F) callconv(.c) i32 {
                    if (n <= 0 or (!trans and m <= 0)) return 0;
                    const rows: BlasLong = @max(m, 0);
                    const ld_a: usize = 2 * @as(usize, @intCast(lda));
                    const alpha = By(V2).init(ar, ai, false);
                    // the vectors' stride: `y`'s, or with a transpose `x`'s
                    const unit = (if (trans) incx else incy) == 1;
                    var j: BlasLong = 0;
                    while (j < n) {
                        const g: BlasLong = if (unit and j + 4 <= n) 4 else 1;
                        var c: [4][*]F = undefined;
                        for (0..@intCast(g)) |q| c[q] = a + (@as(usize, @intCast(j)) + q) * ld_a;
                        if (trans) {
                            var d: [4]Dots(ca, cx) = @splat(.{});
                            var i: BlasLong = 0;
                            if (g == 4) {
                                const len: usize = @intCast(2 * rows);
                                var f: usize = 0;
                                while (f + W <= len) : (f += W) {
                                    const xv = ld(x + f);
                                    inline for (0..4) |q| d[q].add(ld(c[q] + f), xv);
                                }
                                i = @intCast(f / 2);
                            }
                            while (i < rows) : (i += 1) {
                                const xv = ld2(x, i * incx);
                                for (0..@intCast(g)) |q| d[q].add2(ld2(c[q], i), xv);
                            }
                            for (0..@intCast(g)) |q| {
                                const at_y = (j + @as(BlasLong, @intCast(q))) * incy;
                                st2(y, at_y, alpha.madd(d[q].finish(), ld2(y, at_y)));
                            }
                        } else {
                            var t: [4]By(V) = undefined;
                            var t2: [4]By(V2) = undefined;
                            for (0..@intCast(g)) |q| {
                                const tv = By(V2).init(ar, ai, cx).mul(ld2(x, (j + @as(BlasLong, @intCast(q))) * incx));
                                t[q] = .init(tv[0], tv[1], ca);
                                t2[q] = .init(tv[0], tv[1], ca);
                            }
                            var i: BlasLong = 0;
                            if (g == 4) {
                                const len: usize = @intCast(2 * rows);
                                var f: usize = 0;
                                while (f + W <= len) : (f += W) {
                                    var yv = ld(y + f);
                                    inline for (0..4) |q| yv = t[q].madd(ld(c[q] + f), yv);
                                    st(y + f, yv);
                                }
                                i = @intCast(f / 2);
                            }
                            while (i < rows) : (i += 1) {
                                var yv = ld2(y, i * incy);
                                for (0..@intCast(g)) |q| yv = t2[q].madd(ld2(c[q], i), yv);
                                st2(y, i * incy, yv);
                            }
                        }
                        j += g;
                    }
                    return 0;
                }
            };
        }

        pub const gemv_n = Gemv(false, false, false).k;
        pub const gemv_t = Gemv(true, false, false).k;
        pub const gemv_r = Gemv(false, true, false).k;
        pub const gemv_c = Gemv(true, true, false).k;
        pub const gemv_o = Gemv(false, false, true).k;
        pub const gemv_u = Gemv(true, false, true).k;
        pub const gemv_s = Gemv(false, true, true).k;
        pub const gemv_d = Gemv(true, true, true).k;
    };
}

// ---- the slots they fill -----------------------------------------------------

/// Whether `slot`'s converted kernels give way to the portable or these on
/// every core, so its results follow one contract whatever the core
/// (`readme.md`, "One contract"): `nrm2` sums in `f64` for `s` and `c` and
/// is Blue's for `d` and `z`; in the extremes the first NaN wins; `sdot`
/// sums in `f64`; and no result depends on an operand's address (`peels`).
/// `tools/blas_asm` does not convert their assembly (`roles.zig`'s `skip`,
/// kept in step with this); an architecture's `Registry` puts these in where
/// its vector unit is IEEE.
pub fn routed(comptime slot: Slot) bool {
    return switch (slot) {
        // the extremes, `nrm2` and `sdot`
        .samax_k,
        .damax_k,
        .camax_k,
        .zamax_k,
        .samin_k,
        .damin_k,
        .camin_k,
        .zamin_k,
        .smax_k,
        .dmax_k,
        .smin_k,
        .dmin_k,
        .isamax_k,
        .idamax_k,
        .icamax_k,
        .izamax_k,
        .isamin_k,
        .idamin_k,
        .icamin_k,
        .izamin_k,
        .ismax_k,
        .idmax_k,
        .ismin_k,
        .idmin_k,
        .snrm2_k,
        .dnrm2_k,
        .cnrm2_k,
        .znrm2_k,
        .sdot_k,
        => true,
        else => peels(slot),
    };
}

/// The x86_64 reductions whose converted kernels peel elements until an
/// operand is aligned, so their bits follow its address: `asum_sse.S`,
/// `asum_sse2.S`, `zasum_sse*.S`, `zsum_sse*.S`, `dot_sse2.S`,
/// `ddot_bulldozer.S`, `zdot_sse*.S`, `dgemv_n.S` and their `atom` and
/// `bulldozer` kin. Their C kernels on other cores stay.
pub fn peels(comptime slot: Slot) bool {
    if (sysinfo.target.cpu.arch != .x86_64) return false;
    return switch (slot) {
        .sasum_k, .dasum_k, .casum_k, .zasum_k, .csum_k, .zsum_k, .ddot_k, .cdotu_k, .cdotc_k, .zdotu_k, .zdotc_k, .dgemv_n => true,
        else => false,
    };
}

/// Whether these kernels take `slot` from `core`'s converted ones, beyond
/// `routed`, so its object leaves those out: every slot of theirs but `rot`
/// on arm64's `vortex` and `vortexm4` off Apple, where they match or beat
/// the converted kernels (`tools/blasbench/results-m6.md`). On Apple's
/// systems Accelerate is to serve instead, and the converted kernels stay.
pub fn displaces(comptime core: anytype, comptime slot: Slot) bool {
    if (sysinfo.target.os.tag.isDarwin()) return false;
    return switch (@TypeOf(core)) {
        coretypes.Arm64CoreType => (core == .vortex or core == .vortexm4) and switch (slot) {
            .srot_k, .drot_k, .csrot_k, .zdrot_k => false,
            else => fills(slot),
        },
        else => false,
    };
}

/// Whether one of these kernels fills `slot` on this target.
pub fn fills(comptime slot: Slot) bool {
    inline for (kernels) |k| {
        if (k.slot == slot) return true;
    }
    return false;
}

/// The slot each function of `Real` fills, by its name there.
const RealSlots = struct {
    axpy: Slot,
    dot: Slot,
    scal: Slot,
    nrm2: Slot,
    asum: Slot,
    copy: Slot,
    swap: Slot,
    rot: Slot,
    gemv_n: Slot,
    gemv_t: Slot,
    amax: Slot,
    iamax: Slot,
    amin: Slot,
    iamin: Slot,
    max: Slot,
    imax: Slot,
    min: Slot,
    imin: Slot,
};

/// The slot each function of `Complex` fills, by its name there.
const ComplexSlots = struct {
    axpy: Slot,
    axpyc: Slot,
    dotu: Slot,
    dotc: Slot,
    scal: Slot,
    nrm2: Slot,
    asum: Slot,
    copy: Slot,
    swap: Slot,
    amax: Slot,
    iamax: Slot,
    amin: Slot,
    iamin: Slot,
    rot: Slot,
    gemv_n: Slot,
    gemv_t: Slot,
    gemv_r: Slot,
    gemv_c: Slot,
    gemv_o: Slot,
    gemv_u: Slot,
    gemv_s: Slot,
    gemv_d: Slot,
};

const real_slots = .{
    .{ f32, RealSlots{
        .axpy = .saxpy_k,
        .dot = .sdot_k,
        .scal = .sscal_k,
        .nrm2 = .snrm2_k,
        .asum = .sasum_k,
        .copy = .scopy_k,
        .swap = .sswap_k,
        .rot = .srot_k,
        .gemv_n = .sgemv_n,
        .gemv_t = .sgemv_t,
        .amax = .samax_k,
        .iamax = .isamax_k,
        .amin = .samin_k,
        .iamin = .isamin_k,
        .max = .smax_k,
        .imax = .ismax_k,
        .min = .smin_k,
        .imin = .ismin_k,
    } },
    .{ f64, RealSlots{
        .axpy = .daxpy_k,
        .dot = .ddot_k,
        .scal = .dscal_k,
        .nrm2 = .dnrm2_k,
        .asum = .dasum_k,
        .copy = .dcopy_k,
        .swap = .dswap_k,
        .rot = .drot_k,
        .gemv_n = .dgemv_n,
        .gemv_t = .dgemv_t,
        .amax = .damax_k,
        .iamax = .idamax_k,
        .amin = .damin_k,
        .iamin = .idamin_k,
        .max = .dmax_k,
        .imax = .idmax_k,
        .min = .dmin_k,
        .imin = .idmin_k,
    } },
};

const complex_slots = .{
    .{ f32, ComplexSlots{
        .axpy = .caxpy_k,
        .axpyc = .caxpyc_k,
        .dotu = .cdotu_k,
        .dotc = .cdotc_k,
        .scal = .cscal_k,
        .nrm2 = .cnrm2_k,
        .asum = .casum_k,
        .copy = .ccopy_k,
        .swap = .cswap_k,
        .amax = .camax_k,
        .iamax = .icamax_k,
        .amin = .camin_k,
        .iamin = .icamin_k,
        .rot = .csrot_k,
        .gemv_n = .cgemv_n,
        .gemv_t = .cgemv_t,
        .gemv_r = .cgemv_r,
        .gemv_c = .cgemv_c,
        .gemv_o = .cgemv_o,
        .gemv_u = .cgemv_u,
        .gemv_s = .cgemv_s,
        .gemv_d = .cgemv_d,
    } },
    .{ f64, ComplexSlots{
        .axpy = .zaxpy_k,
        .axpyc = .zaxpyc_k,
        .dotu = .zdotu_k,
        .dotc = .zdotc_k,
        .scal = .zscal_k,
        .nrm2 = .znrm2_k,
        .asum = .zasum_k,
        .copy = .zcopy_k,
        .swap = .zswap_k,
        .amax = .zamax_k,
        .iamax = .izamax_k,
        .amin = .zamin_k,
        .iamin = .izamin_k,
        .rot = .zdrot_k,
        .gemv_n = .zgemv_n,
        .gemv_t = .zgemv_t,
        .gemv_r = .zgemv_r,
        .gemv_c = .zgemv_c,
        .gemv_o = .zgemv_o,
        .gemv_u = .zgemv_u,
        .gemv_s = .zgemv_s,
        .gemv_d = .zgemv_d,
    } },
};

/// A kernel: the slot it fills, and the function of `Real` or `Complex` over
/// `F` that is it.
const Kernel = struct { slot: Slot, F: type, complex: bool, f: [:0]const u8 };

/// Every kernel here: none in a type whose lanes flush denormals
/// (`util.vectorIeee`), whose slots keep the portable kernels.
pub const kernels: []const Kernel = blk: {
    var r: []const Kernel = &.{};
    for (real_slots) |p| {
        if (!util.vectorIeee(p[0])) continue;
        for (@typeInfo(RealSlots).@"struct".field_names) |name| {
            r = r ++ &[_]Kernel{.{ .slot = @field(p[1], name), .F = p[0], .complex = false, .f = name }};
        }
    }
    for (complex_slots) |p| {
        if (!util.vectorIeee(p[0])) continue;
        for (@typeInfo(ComplexSlots).@"struct".field_names) |name| {
            r = r ++ &[_]Kernel{.{ .slot = @field(p[1], name), .F = p[0], .complex = true, .f = name }};
        }
    }
    break :blk r;
};

/// The kernel of `k`, `U` vectors a step.
fn kernel(comptime k: Kernel, comptime U: comptime_int) util.SlotFn(k.slot) {
    return &@field(if (k.complex) Complex(k.F, U) else Real(k.F, U), k.f);
}

/// The symbol of `slot`'s kernel here for `core`, apart from a converted one's.
fn symbol(comptime slot: Slot, comptime core: anytype) []const u8 {
    return util.kernelSymbol(slot, core) ++ "_vec";
}

/// Points every slot these kernels fill at this build's own, four vectors a
/// step, for the tests.
pub fn installAll(vt: *Vtable) void {
    inline for (kernels) |k| @field(vt, @tagName(k.slot)) = kernel(k, 4);
}

/// `Registry`'s `uses` for an architecture whose vector unit is IEEE on
/// every core: the `routed` slots.
pub fn usesRouted(comptime core: anytype, comptime slot: Slot) bool {
    _ = core;
    return routed(slot);
}

/// The kernels four vectors a step, in a per-core object for the cores and
/// slots `uses(core, slot)` picks. Like the gemm kernels, they never make a
/// core need an object of its own.
pub fn Registry(comptime Core: type, comptime uses: anytype) type {
    return struct {
        pub fn has(comptime core: Core) bool {
            _ = core;
            return false;
        }

        pub fn exportCore(comptime core: Core) void {
            inline for (kernels) |k| {
                if (comptime uses(core, k.slot)) @export(kernel(k, 4), .{ .name = symbol(k.slot, core), .visibility = .hidden });
            }
        }

        pub fn install(vt: *Vtable, comptime core: Core) void {
            inline for (kernels) |k| {
                if (comptime uses(core, k.slot)) {
                    @field(vt, @tagName(k.slot)) = @extern(util.SlotFn(k.slot), .{ .name = symbol(k.slot, core), .visibility = .hidden });
                }
            }
        }
    };
}
