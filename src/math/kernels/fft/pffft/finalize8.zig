//! The 8-lane steps around pffft's passes (N35), derived from pffft's 4x4
//! ones (`finalize.zig`, marton78/pffft efc1dd2d, UCAR, LICENSE) for AVX.
//! Each lane `j` transforms `x[8 m + j]`; a block of 8 bins `q` is
//! transposed, twiddled by `W_n^(j q)` and joined by a radix-8 butterfly
//! across `j`. Complex blocks hold `X[q + p M]`, `p` 0 to 7 (`M = n / 8`).
//! Real blocks hold `X[q + p M]` for `p` below 4 and, between them,
//! `conj X[q + (7 - p) M]`, which is `X[(p + 1) M - q]`; block 0's lane 0
//! holds `X[0]` and `X[n / 2]`, then `X[M / 2 + p M]` in the mirrored slots.

const std = @import("std");
const simd_file = @import("simd.zig");

const simd = simd_file.For(8);
const V = simd.V;
const Ptr = simd.Ptr;
const ConstPtr = simd.ConstPtr;
const cplxMul = simd.cplxMul;
const cplxMulConj = simd.cplxMulConj;
const s = simd.lit(simd.m_sqrt1_2);

const C = struct {
    re: V,
    im: V,

    fn add(a: C, b: C) C {
        return .{ .re = a.re + b.re, .im = a.im + b.im };
    }

    fn sub(a: C, b: C) C {
        return .{ .re = a.re - b.re, .im = a.im - b.im };
    }

    fn conj(a: C) C {
        return .{ .re = a.re, .im = -a.im };
    }

    /// `a` times `-i` (forward) or `i`.
    fn rot(a: C, comptime inverse: bool) C {
        return if (inverse) .{ .re = -a.im, .im = a.re } else .{ .re = a.im, .im = -a.re };
    }
};

/// `Y_q = sum_m W_4^(±m q) a_m`.
inline fn dft4(a: [4]C, comptime inverse: bool) [4]C {
    const s02 = a[0].add(a[2]);
    const d02 = a[0].sub(a[2]);
    const s13 = a[1].add(a[3]);
    const r13 = a[1].sub(a[3]).rot(inverse);
    return .{ s02.add(s13), d02.add(r13), s02.sub(s13), d02.sub(r13) };
}

/// `X_p = sum_j W_8^(±j p) t_j`, as two radix-4 halves joined.
inline fn dft8(t: [8]C, comptime inverse: bool) [8]C {
    const e = dft4(.{ t[0], t[2], t[4], t[6] }, inverse);
    const o = dft4(.{ t[1], t[3], t[5], t[7] }, inverse);
    const sv: V = @splat(s);
    // W_8^p o_p, p 0 to 3; conjugate twiddles inverse
    const w = [4]C{
        o[0],
        if (inverse)
            .{ .re = sv * (o[1].re - o[1].im), .im = sv * (o[1].re + o[1].im) }
        else
            .{ .re = sv * (o[1].re + o[1].im), .im = sv * (o[1].im - o[1].re) },
        o[2].rot(inverse),
        if (inverse)
            .{ .re = -sv * (o[3].re + o[3].im), .im = sv * (o[3].re - o[3].im) }
        else
            .{ .re = sv * (o[3].im - o[3].re), .im = -sv * (o[3].re + o[3].im) },
    };
    var x: [8]C = undefined;
    inline for (0..4) |p| {
        x[p] = e[p].add(w[p]);
        x[p + 4] = e[p].sub(w[p]);
    }
    return x;
}

/// Twiddles block `k`'s lanes `j` 1 to 7 by `e` (`conj` undoes it).
inline fn twiddle(t: *[8]C, e: ConstPtr, k: usize, comptime conj: bool) void {
    inline for (1..8) |j| {
        const er = e[k * 14 + 2 * (j - 1)];
        const ei = e[k * 14 + 2 * (j - 1) + 1];
        if (conj) cplxMulConj(&t[j].re, &t[j].im, er, ei) else cplxMul(&t[j].re, &t[j].im, er, ei);
    }
}

/// Bins as rows to lanes as rows, or back.
inline fn transpose(t: *[8]C) void {
    var re: [8]V = undefined;
    var im: [8]V = undefined;
    inline for (0..8) |i| {
        re[i] = t[i].re;
        im[i] = t[i].im;
    }
    simd.transpose(&re);
    simd.transpose(&im);
    inline for (0..8) |i| t[i] = .{ .re = re[i], .im = im[i] };
}

pub fn cplxFinalize(ncvec: usize, in: ConstPtr, out: Ptr, e: ConstPtr) void {
    std.debug.assert(in != out);
    for (0..ncvec / 8) |k| {
        var t: [8]C = undefined;
        inline for (0..8) |i| t[i] = .{ .re = in[16 * k + 2 * i], .im = in[16 * k + 2 * i + 1] };
        transpose(&t);
        twiddle(&t, e, k, false);
        const x = dft8(t, false);
        inline for (0..8) |p| out[16 * k + 2 * p ..][0..2].* = .{ x[p].re, x[p].im };
    }
}

pub fn cplxPreprocess(ncvec: usize, in: ConstPtr, out: Ptr, e: ConstPtr) void {
    std.debug.assert(in != out);
    for (0..ncvec / 8) |k| {
        var x: [8]C = undefined;
        inline for (0..8) |p| x[p] = .{ .re = in[16 * k + 2 * p], .im = in[16 * k + 2 * p + 1] };
        var t = dft8(x, true);
        twiddle(&t, e, k, true);
        transpose(&t);
        inline for (0..8) |i| out[16 * k + 2 * i ..][0..2].* = .{ t[i].re, t[i].im };
    }
}

/// `cos` and `sin` of `pi m / 8`, `m` 0 to 15.
const turn = blk: {
    var c: [16]f32 = undefined;
    var sn: [16]f32 = undefined;
    for (0..16) |m| {
        const a = @as(f64, std.math.pi) * @as(f64, @floatFromInt(m)) / 8;
        c[m] = @floatCast(@cos(a));
        sn[m] = @floatCast(@sin(a));
    }
    break :blk .{ c, sn };
};

/// The lanes' halfcomplex spectra (FFTPACK's order, `M` vectors) to the
/// real blocks.
pub noinline fn realFinalize(ncvec: usize, in: ConstPtr, out: Ptr, e: ConstPtr) void {
    std.debug.assert(in != out);
    const zero: V = @splat(0);
    for (0..ncvec / 8) |k| {
        var t: [8]C = undefined;
        inline for (0..8) |i| {
            // `Y[q]` at `in[2 q - 1]`, `in[2 q]`; bin 0 is block 0's special
            t[i] = if (k == 0 and i == 0) .{ .re = zero, .im = zero } else .{ .re = in[16 * k + 2 * i - 1], .im = in[16 * k + 2 * i] };
        }
        transpose(&t);
        twiddle(&t, e, k, false);
        const x = dft8(t, false);
        inline for (0..4) |p| {
            out[16 * k + 4 * p ..][0..2].* = .{ x[p].re, x[p].im };
            const m = x[7 - p].conj();
            out[16 * k + 4 * p + 2 ..][0..2].* = .{ m.re, m.im };
        }
    }
    // block 0's lane 0, from the lanes' `Y[0]` and `Y[M / 2]`, both real
    const c0: [8]f32 = in[0];
    const c1: [8]f32 = in[2 * ncvec - 1];
    const u = simd.floats(out);
    var dc: f32 = 0;
    var ny: f32 = 0;
    for (0..8) |j| {
        dc += c0[j];
        ny += if (j % 2 == 0) c0[j] else -c0[j];
    }
    u[0] = dc;
    u[8] = ny;
    for (1..4) |p| {
        // X[p M] = sum_j W_8^(j p) c0_j, in slot 2 p
        var re: f32 = 0;
        var im: f32 = 0;
        for (0..8) |j| {
            const m = (2 * j * p) % 16;
            re += turn[0][m] * c0[j];
            im -= turn[1][m] * c0[j];
        }
        u[(4 * p) * 8] = re;
        u[(4 * p + 1) * 8] = im;
    }
    for (0..4) |p| {
        // X[M / 2 + p M] = sum_j W_16^(j (2 p + 1)) c1_j, in slot 2 p + 1
        var re: f32 = 0;
        var im: f32 = 0;
        for (0..8) |j| {
            const m = (j * (2 * p + 1)) % 16;
            re += turn[0][m] * c1[j];
            im -= turn[1][m] * c1[j];
        }
        u[(4 * p + 2) * 8] = re;
        u[(4 * p + 3) * 8] = im;
    }
}

/// `realFinalize` undone, unscaled: the lanes' halfcomplex spectra times 8.
pub noinline fn realPreprocess(ncvec: usize, in: ConstPtr, out: Ptr, e: ConstPtr) void {
    std.debug.assert(in != out);
    for (0..ncvec / 8) |k| {
        var x: [8]C = undefined;
        inline for (0..4) |p| {
            x[p] = .{ .re = in[16 * k + 4 * p], .im = in[16 * k + 4 * p + 1] };
            x[7 - p] = (C{ .re = in[16 * k + 4 * p + 2], .im = in[16 * k + 4 * p + 3] }).conj();
        }
        var t = dft8(x, true);
        twiddle(&t, e, k, true);
        transpose(&t);
        inline for (0..8) |i| {
            if (k != 0 or i != 0) out[16 * k + 2 * i - 1 ..][0..2].* = .{ t[i].re, t[i].im };
        }
    }
    const f = simd.constFloats(in);
    var xr: [8]f32 = undefined;
    var xi: [8]f32 = undefined;
    for (0..8) |slot| {
        xr[slot] = f[(2 * slot) * 8];
        xi[slot] = f[(2 * slot + 1) * 8];
    }
    var c0: [8]f32 = undefined;
    var c1: [8]f32 = undefined;
    for (0..8) |j| {
        // 8 Y_j[0] = X[0] + (-1)^j X[n/2] + 2 Re sum_p W_8^(-j p) X[p M]
        var a: f32 = 0;
        for (1..4) |p| {
            const m = (2 * j * p) % 16;
            a += turn[0][m] * xr[2 * p] - turn[1][m] * xi[2 * p];
        }
        c0[j] = (xr[0] + (if (j % 2 == 0) xi[0] else -xi[0])) + 2 * a;
        // 8 Y_j[M/2] = 2 Re sum_p W_16^(-j (2 p + 1)) X[M/2 + p M]
        var b: f32 = 0;
        for (0..4) |p| {
            const m = (j * (2 * p + 1)) % 16;
            b += turn[0][m] * xr[2 * p + 1] - turn[1][m] * xi[2 * p + 1];
        }
        c1[j] = 2 * b;
    }
    out[0] = c0;
    out[2 * ncvec - 1] = c1;
}

fn reverse(v: V) V {
    return @shuffle(f32, v, undefined, [8]i32{ 7, 6, 5, 4, 3, 2, 1, 0 });
}

inline fn at(p: [*]f32, i: usize) Ptr {
    return @ptrCast(p + i);
}

inline fn atConst(p: [*]const f32, i: usize) ConstPtr {
    return @ptrCast(p + i);
}

/// pffft's order to interleaved complex values (forward), or back; real
/// spectra as pffft's ordered form, `X[0]`, `X[n/2]`, then `X[1]` on.
/// `in` and `out` must not alias.
pub fn zreorder(n: usize, ncvec: usize, real: bool, in: [*]const f32, out: [*]f32, forward: bool) void {
    const vin: ConstPtr = @ptrCast(in);
    const vout: Ptr = @ptrCast(out);
    std.debug.assert(@intFromPtr(in) != @intFromPtr(out));
    if (!real) {
        for (0..ncvec) |k| {
            const kk = (k / 8) + (k % 8) * (ncvec / 8);
            if (forward) {
                vout[kk * 2 ..][0..2].* = simd.interleave2(vin[k * 2], vin[k * 2 + 1]);
            } else {
                vout[k * 2 ..][0..2].* = simd.uninterleave2(vin[kk * 2], vin[kk * 2 + 1]);
            }
        }
        return;
    }
    const m = n / 8;
    for (0..ncvec / 8) |k| {
        for (0..4) |p| {
            // the direct slot: `X[8 k + l + p M]`, rising with `l`
            const d = 2 * (8 * k + p * m);
            // the mirrored slot: `X[(p + 1) M - 8 k - l]`, falling
            const r = 2 * ((p + 1) * m - 8 * k - 7);
            const v = 16 * k + 4 * p;
            if (forward) {
                at(out, d)[0..2].* = simd.interleave2(vin[v], vin[v + 1]);
            } else {
                vout[v..][0..2].* = simd.uninterleave2(atConst(in, d)[0], atConst(in, d)[1]);
            }
            if (k == 0) {
                // lane 0 is `X[M / 2 + p M]`; the rest one at a time
                const fv = 8 * (v + 2);
                for (0..8) |l| {
                    const pos = if (l == 0) m / 2 + p * m else (p + 1) * m - l;
                    if (forward) {
                        out[2 * pos] = in[fv + l];
                        out[2 * pos + 1] = in[fv + 8 + l];
                    } else {
                        out[fv + l] = in[2 * pos];
                        out[fv + 8 + l] = in[2 * pos + 1];
                    }
                }
            } else if (forward) {
                at(out, r)[0..2].* = simd.interleave2(reverse(vin[v + 2]), reverse(vin[v + 3]));
            } else {
                const g = simd.uninterleave2(atConst(in, r)[0], atConst(in, r)[1]);
                vout[v + 2 ..][0..2].* = .{ reverse(g[0]), reverse(g[1]) };
            }
        }
    }
}
