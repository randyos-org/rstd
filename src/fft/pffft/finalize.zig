//! pffft's 4x4 block steps around the passes, `pffft_cplx_finalize`,
//! `pffft_cplx_preprocess`, `pffft_real_finalize` and
//! `pffft_real_preprocess`, and the reordering, `pffft_zreorder`
//! (`src/pffft_priv_impl.h`, marton78/pffft efc1dd2d, UCAR, LICENSE), in
//! its operation order. The 4-lane steps; `finalize8.zig` has the 8-lane ones.

const std = @import("std");
const simd = @import("simd.zig").For(4);

const V = simd.V;
const Ptr = simd.Ptr;
const ConstPtr = simd.ConstPtr;
const lit = simd.lit;
const cplxMul = simd.cplxMul;
const cplxMulConj = simd.cplxMulConj;
const transpose4 = simd.transpose4;
const swapHL = simd.swapHL;

pub fn cplxFinalize(ncvec: usize, in: ConstPtr, out_: Ptr, e: ConstPtr) void {
    var out = out_;
    const dk = ncvec / 4; // 4x4 blocks
    std.debug.assert(in != out);
    for (0..dk) |k| {
        var r0 = in[8 * k + 0];
        var im0 = in[8 * k + 1];
        var r1 = in[8 * k + 2];
        var im1 = in[8 * k + 3];
        var r2 = in[8 * k + 4];
        var im2 = in[8 * k + 5];
        var r3 = in[8 * k + 6];
        var im3 = in[8 * k + 7];
        transpose4(&r0, &r1, &r2, &r3);
        transpose4(&im0, &im1, &im2, &im3);
        cplxMul(&r1, &im1, e[k * 6 + 0], e[k * 6 + 1]);
        cplxMul(&r2, &im2, e[k * 6 + 2], e[k * 6 + 3]);
        cplxMul(&r3, &im3, e[k * 6 + 4], e[k * 6 + 5]);

        const sr0 = r0 + r2;
        const dr0 = r0 - r2;
        const sr1 = r1 + r3;
        const dr1 = r1 - r3;
        const si0 = im0 + im2;
        const di0 = im0 - im2;
        const si1 = im1 + im3;
        const di1 = im1 - im3;

        r0 = sr0 + sr1;
        im0 = si0 + si1;
        r1 = dr0 + di1;
        im1 = di0 - dr1;
        r2 = sr0 - sr1;
        im2 = si0 - si1;
        r3 = dr0 - di1;
        im3 = di0 + dr1;

        out[0..8].* = .{ r0, im0, r1, im1, r2, im2, r3, im3 };
        out += 8;
    }
}

pub fn cplxPreprocess(ncvec: usize, in: ConstPtr, out_: Ptr, e: ConstPtr) void {
    var out = out_;
    const dk = ncvec / 4;
    std.debug.assert(in != out);
    for (0..dk) |k| {
        var r0 = in[8 * k + 0];
        var im0 = in[8 * k + 1];
        var r1 = in[8 * k + 2];
        var im1 = in[8 * k + 3];
        var r2 = in[8 * k + 4];
        var im2 = in[8 * k + 5];
        var r3 = in[8 * k + 6];
        var im3 = in[8 * k + 7];

        const sr0 = r0 + r2;
        const dr0 = r0 - r2;
        const sr1 = r1 + r3;
        const dr1 = r1 - r3;
        const si0 = im0 + im2;
        const di0 = im0 - im2;
        const si1 = im1 + im3;
        const di1 = im1 - im3;

        r0 = sr0 + sr1;
        im0 = si0 + si1;
        r1 = dr0 - di1;
        im1 = di0 + dr1;
        r2 = sr0 - sr1;
        im2 = si0 - si1;
        r3 = dr0 + di1;
        im3 = di0 - dr1;

        cplxMulConj(&r1, &im1, e[k * 6 + 0], e[k * 6 + 1]);
        cplxMulConj(&r2, &im2, e[k * 6 + 2], e[k * 6 + 3]);
        cplxMulConj(&r3, &im3, e[k * 6 + 4], e[k * 6 + 5]);

        transpose4(&r0, &r1, &r2, &r3);
        transpose4(&im0, &im1, &im2, &im3);

        out[0..8].* = .{ r0, im0, r1, im1, r2, im2, r3, im3 };
        out += 8;
    }
}

inline fn realFinalize4x4(in0: V, in1: V, in: ConstPtr, e: ConstPtr, out: Ptr) void {
    var r0 = in0;
    var im0 = in1;
    var r1 = in[0];
    var im1 = in[1];
    var r2 = in[2];
    var im2 = in[3];
    var r3 = in[4];
    var im3 = in[5];
    transpose4(&r0, &r1, &r2, &r3);
    transpose4(&im0, &im1, &im2, &im3);

    cplxMul(&r1, &im1, e[0], e[1]);
    cplxMul(&r2, &im2, e[2], e[3]);
    cplxMul(&r3, &im3, e[4], e[5]);

    const sr0 = r0 + r2;
    const dr0 = r0 - r2;
    const sr1 = r1 + r3;
    const dr1 = r3 - r1;
    const si0 = im0 + im2;
    const di0 = im0 - im2;
    const si1 = im1 + im3;
    const di1 = im3 - im1;

    r0 = sr0 + sr1;
    r3 = sr0 - sr1;
    im0 = si0 + si1;
    im3 = si1 - si0;
    r1 = dr0 + di1;
    r2 = dr0 - di1;
    im1 = dr1 - di0;
    im2 = dr1 + di0;

    out[0..8].* = .{ r0, im0, r1, im1, r2, im2, r3, im3 };
}

/// FFTPACK's order to pffft's: `f0r f1r f1i f2r f2i ... f(n-1)r f(n-1)i f(n)r`.
pub noinline fn realFinalize(ncvec: usize, in: ConstPtr, out: Ptr, e: ConstPtr) void {
    const dk = ncvec / 4;
    const uout = simd.floats(out);
    var save = in[7];
    const zero: V = @splat(0);
    const s = lit(simd.m_sqrt1_2);

    const cr = in[0];
    const ci = in[ncvec * 2 - 1];
    std.debug.assert(in != out);
    realFinalize4x4(zero, zero, in + 1, e, out);

    const xr0 = (cr[0] + cr[2]) + (cr[1] + cr[3]);
    uout[0 * 4] = xr0;
    const xi0 = (cr[0] + cr[2]) - (cr[1] + cr[3]);
    uout[1 * 4] = xi0;
    const xr2 = (cr[0] - cr[2]);
    uout[4 * 4] = xr2;
    const xi2 = (cr[3] - cr[1]);
    uout[5 * 4] = xi2;
    const xr1 = ci[0] + s * (ci[1] - ci[3]);
    uout[2 * 4] = xr1;
    const xi1 = -ci[2] - s * (ci[1] + ci[3]);
    uout[3 * 4] = xi1;
    const xr3 = ci[0] - s * (ci[1] - ci[3]);
    uout[6 * 4] = xr3;
    const xi3 = ci[2] - s * (ci[1] + ci[3]);
    uout[7 * 4] = xi3;

    for (1..dk) |k| {
        const save_next = in[8 * k + 7];
        realFinalize4x4(save, in[8 * k + 0], in + 8 * k + 1, e + k * 6, out + k * 8);
        save = save_next;
    }
}

inline fn realPreprocess4x4(in: ConstPtr, e: ConstPtr, out_: Ptr, comptime first: bool) void {
    var out = out_;
    var r0 = in[0];
    var im0 = in[1];
    var r1 = in[2];
    var im1 = in[3];
    var r2 = in[4];
    var im2 = in[5];
    var r3 = in[6];
    var im3 = in[7];

    const sr0 = r0 + r3;
    const dr0 = r0 - r3;
    const sr1 = r1 + r2;
    const dr1 = r1 - r2;
    const si0 = im0 + im3;
    const di0 = im0 - im3;
    const si1 = im1 + im2;
    const di1 = im1 - im2;

    r0 = sr0 + sr1;
    r2 = sr0 - sr1;
    r1 = dr0 - si1;
    r3 = dr0 + si1;
    im0 = di0 - di1;
    im2 = di0 + di1;
    im1 = si0 - dr1;
    im3 = si0 + dr1;

    cplxMulConj(&r1, &im1, e[0], e[1]);
    cplxMulConj(&r2, &im2, e[2], e[3]);
    cplxMulConj(&r3, &im3, e[4], e[5]);

    transpose4(&r0, &r1, &r2, &r3);
    transpose4(&im0, &im1, &im2, &im3);

    if (!first) {
        out[0..2].* = .{ r0, im0 };
        out += 2;
    }
    out[0..6].* = .{ r1, im1, r2, im2, r3, im3 };
}

/// pffft's order back to FFTPACK's, `realFinalize` undone.
pub noinline fn realPreprocess(ncvec: usize, in: ConstPtr, out: Ptr, e: ConstPtr) void {
    const dk = ncvec / 4;
    const uout = simd.floats(out);
    const fin = simd.constFloats(in);
    const s = lit(simd.m_sqrt2);
    std.debug.assert(in != out);
    var xr: [4]f32 = undefined;
    var xi: [4]f32 = undefined;
    for (0..4) |k| {
        xr[k] = fin[8 * k];
        xi[k] = fin[8 * k + 4];
    }

    realPreprocess4x4(in, e, out + 1, true); // writes 6 vectors

    for (1..dk) |k| {
        realPreprocess4x4(in + 8 * k, e + k * 6, out - 1 + k * 8, false);
    }

    const last = (2 * ncvec - 1) * 4;
    const cr0 = (xr[0] + xi[0]) + 2 * xr[2];
    uout[0] = cr0;
    const cr1 = (xr[0] - xi[0]) - 2 * xi[2];
    uout[1] = cr1;
    const cr2 = (xr[0] + xi[0]) - 2 * xr[2];
    uout[2] = cr2;
    const cr3 = (xr[0] - xi[0]) + 2 * xi[2];
    uout[3] = cr3;
    const ci0 = 2 * (xr[1] + xr[3]);
    uout[last + 0] = ci0;
    const ci1 = s * (xr[1] - xr[3]) - s * (xi[1] + xi[3]);
    uout[last + 1] = ci1;
    const ci2 = 2 * (xi[3] - xi[1]);
    uout[last + 2] = ci2;
    const ci3 = -s * (xr[1] - xr[3]) - s * (xi[1] + xi[3]);
    uout[last + 3] = ci3;
}

/// `[0 0 1 2 3 4 5 6 7 8] -> [0 8 7 6 5 4 3 2 1]`; `out` is one past the end.
fn reversedCopy(n: usize, in_: ConstPtr, in_stride: usize, out_: Ptr) void {
    var in = in_;
    var out = out_;
    var g = simd.interleave2(in[0], in[1]);
    in += in_stride;

    out -= 1;
    out[0] = swapHL(g[0], g[1]); // [g0l, g0h], [g1l g1h] -> [g1l, g0h]
    for (1..n) |_| {
        const h = simd.interleave2(in[0], in[1]);
        in += in_stride;
        out -= 1;
        out[0] = swapHL(g[1], h[0]);
        out -= 1;
        out[0] = swapHL(h[0], h[1]);
        g[1] = h[1];
    }
    out -= 1;
    out[0] = swapHL(g[1], g[0]);
}

/// `reversedCopy` undone; `out` steps back 8 vectors at a time.
fn unreversedCopy(n: usize, in_: ConstPtr, out_: Ptr) void {
    var in = in_;
    var out = out_;
    const g0 = in[0];
    var g1 = g0;
    in += 1;
    for (1..n) |_| {
        var h0 = in[0];
        const h1 = in[1];
        in += 2;
        g1 = swapHL(g1, h0);
        h0 = swapHL(h0, h1);
        const u = simd.uninterleave2(h0, g1);
        out[0] = u[0];
        out[1] = u[1];
        out -= 8;
        g1 = h1;
    }
    var h0 = in[0];
    const h1 = g0;
    g1 = swapHL(g1, h0);
    h0 = swapHL(h0, h1);
    const u = simd.uninterleave2(h0, g1);
    out[0] = u[0];
    out[1] = u[1];
}

/// `pffft_zreorder`: pffft's order to interleaved complex values
/// (forward), or back. `in` and `out` must not alias.
pub fn zreorder(n: usize, ncvec: usize, real: bool, in: [*]const f32, out: [*]f32, forward: bool) void {
    const vin: ConstPtr = @ptrCast(in);
    const vout: Ptr = @ptrCast(out);
    std.debug.assert(@intFromPtr(in) != @intFromPtr(out));
    if (real) {
        const dk = n / 32;
        if (forward) {
            for (0..dk) |k| {
                vout[2 * (0 * dk + k) ..][0..2].* = simd.interleave2(vin[k * 8 + 0], vin[k * 8 + 1]);
                vout[2 * (2 * dk + k) ..][0..2].* = simd.interleave2(vin[k * 8 + 4], vin[k * 8 + 5]);
            }
            reversedCopy(dk, vin + 2, 8, @ptrCast(out + n / 2));
            reversedCopy(dk, vin + 6, 8, @ptrCast(out + n));
        } else {
            for (0..dk) |k| {
                vout[k * 8 + 0 ..][0..2].* = simd.uninterleave2(vin[2 * (0 * dk + k) + 0], vin[2 * (0 * dk + k) + 1]);
                vout[k * 8 + 4 ..][0..2].* = simd.uninterleave2(vin[2 * (2 * dk + k) + 0], vin[2 * (2 * dk + k) + 1]);
            }
            unreversedCopy(dk, @ptrCast(in + n / 4), @ptrCast(out + n - 6 * 4));
            unreversedCopy(dk, @ptrCast(in + 3 * n / 4), @ptrCast(out + n - 2 * 4));
        }
    } else {
        if (forward) {
            for (0..ncvec) |k| {
                const kk = (k / 4) + (k % 4) * (ncvec / 4);
                vout[kk * 2 ..][0..2].* = simd.interleave2(vin[k * 2], vin[k * 2 + 1]);
            }
        } else {
            for (0..ncvec) |k| {
                const kk = (k / 4) + (k % 4) * (ncvec / 4);
                vout[k * 2 ..][0..2].* = simd.uninterleave2(vin[kk * 2], vin[kk * 2 + 1]);
            }
        }
    }
}
