//! pffft's real backward passes, `radb2_ps` to `radb5_ps` and `rfftb1_ps`
//! (`src/pffft_priv_impl.h`, marton78/pffft efc1dd2d, UCAR, LICENSE), in
//! its operation order.

const std = @import("std");
const simd = @import("simd.zig");

const V = simd.V;
const Ptr = simd.Ptr;
const ConstPtr = simd.ConstPtr;
const ld = simd.ld;
const lit = simd.lit;
const smul = simd.smul;
const smadd = simd.smadd;
const smsub = simd.smsub;
const cplxMul = simd.cplxMul;

noinline fn radb2(ido: usize, l1: usize, cc: ConstPtr, ch: Ptr, wa1: [*]const f32) void {
    const minus_two: f32 = -2;
    const l1ido = l1 * ido;
    var k: usize = 0;
    while (k < l1ido) : (k += ido) {
        const a = cc[2 * k];
        const b = cc[2 * (k + ido) - 1];
        ch[k] = a + b;
        ch[k + l1ido] = a - b;
    }
    if (ido < 2) return;
    if (ido != 2) {
        k = 0;
        while (k < l1ido) : (k += ido) {
            var i: usize = 2;
            while (i < ido) : (i += 2) {
                const a = cc[i - 1 + 2 * k];
                const b = cc[2 * (k + ido) - i - 1];
                const c = cc[i + 0 + 2 * k];
                const d = cc[2 * (k + ido) - i + 0];
                ch[i - 1 + k] = a + b;
                var tr2 = a - b;
                ch[i + 0 + k] = c - d;
                var ti2 = c + d;
                cplxMul(&tr2, &ti2, ld(wa1[i - 2]), ld(wa1[i - 1]));
                ch[i - 1 + k + l1ido] = tr2;
                ch[i + 0 + k + l1ido] = ti2;
            }
        }
        if (ido % 2 == 1) return;
    }
    k = 0;
    while (k < l1ido) : (k += ido) {
        const a = cc[2 * k + ido - 1];
        const b = cc[2 * k + ido];
        ch[k + ido - 1] = a + a;
        ch[k + ido - 1 + l1ido] = smul(minus_two, b);
    }
}

fn radb3(ido: usize, l1: usize, cc: ConstPtr, ch: Ptr, wa1: [*]const f32, wa2: [*]const f32) void {
    const taur: f32 = -0.5; // cos(2 pi / 3)
    const taui = lit(0.86602540378443864676); // sin(2 pi / 3)
    const taui_2 = lit(1.7320508075688772935); // 2 sin(2 pi / 3)
    for (0..l1) |k| {
        var tr2 = cc[ido - 1 + (3 * k + 1) * ido];
        tr2 = tr2 + tr2;
        const cr2 = smadd(taur, tr2, cc[3 * k * ido]);
        ch[k * ido] = cc[3 * k * ido] + tr2;
        const ci3 = smul(taui_2, cc[(3 * k + 2) * ido]);
        ch[(k + l1) * ido] = cr2 - ci3;
        ch[(k + 2 * l1) * ido] = cr2 + ci3;
    }
    if (ido == 1) return;
    if (ido == 3) {
        // the inner loop's one pass, i = 2 and ic = 1
        const wa1_r = ld(wa1[0]);
        const wa1_i = ld(wa1[1]);
        const wa2_r = ld(wa2[0]);
        const wa2_i = ld(wa2[1]);
        for (0..l1) |k| {
            const tr2 = cc[1 + (3 * k + 2) * ido] + cc[0 + (3 * k + 1) * ido];
            const cr2 = smadd(taur, tr2, cc[1 + 3 * k * ido]);
            ch[1 + k * ido] = cc[1 + 3 * k * ido] + tr2;
            const ti2 = cc[2 + (3 * k + 2) * ido] - cc[1 + (3 * k + 1) * ido];
            const ci2 = smadd(taur, ti2, cc[2 + 3 * k * ido]);
            ch[2 + k * ido] = cc[2 + 3 * k * ido] + ti2;
            const cr3 = smul(taui, cc[1 + (3 * k + 2) * ido] - cc[0 + (3 * k + 1) * ido]);
            const ci3 = smul(taui, cc[2 + (3 * k + 2) * ido] + cc[1 + (3 * k + 1) * ido]);
            var dr2 = cr2 - ci3;
            var dr3 = cr2 + ci3;
            var di2 = ci2 + cr3;
            var di3 = ci2 - cr3;
            cplxMul(&dr2, &di2, wa1_r, wa1_i);
            ch[1 + (k + l1) * ido] = dr2;
            ch[2 + (k + l1) * ido] = di2;
            cplxMul(&dr3, &di3, wa2_r, wa2_i);
            ch[1 + (k + 2 * l1) * ido] = dr3;
            ch[2 + (k + 2 * l1) * ido] = di3;
        }
        return;
    }
    for (0..l1) |k| {
        var i: usize = 2;
        while (i < ido) : (i += 2) {
            const ic = ido - i;
            const tr2 = cc[i - 1 + (3 * k + 2) * ido] + cc[ic - 1 + (3 * k + 1) * ido];
            const cr2 = smadd(taur, tr2, cc[i - 1 + 3 * k * ido]);
            ch[i - 1 + k * ido] = cc[i - 1 + 3 * k * ido] + tr2;
            const ti2 = cc[i + (3 * k + 2) * ido] - cc[ic + (3 * k + 1) * ido];
            const ci2 = smadd(taur, ti2, cc[i + 3 * k * ido]);
            ch[i + k * ido] = cc[i + 3 * k * ido] + ti2;
            const cr3 = smul(taui, cc[i - 1 + (3 * k + 2) * ido] - cc[ic - 1 + (3 * k + 1) * ido]);
            const ci3 = smul(taui, cc[i + (3 * k + 2) * ido] + cc[ic + (3 * k + 1) * ido]);
            var dr2 = cr2 - ci3;
            var dr3 = cr2 + ci3;
            var di2 = ci2 + cr3;
            var di3 = ci2 - cr3;
            cplxMul(&dr2, &di2, ld(wa1[i - 2]), ld(wa1[i - 1]));
            ch[i - 1 + (k + l1) * ido] = dr2;
            ch[i + (k + l1) * ido] = di2;
            cplxMul(&dr3, &di3, ld(wa2[i - 2]), ld(wa2[i - 1]));
            ch[i - 1 + (k + 2 * l1) * ido] = dr3;
            ch[i + (k + 2 * l1) * ido] = di3;
        }
    }
}

noinline fn radb4(ido: usize, l1: usize, cc: ConstPtr, ch: Ptr, wa1: [*]const f32, wa2: [*]const f32, wa3: [*]const f32) void {
    const minus_sqrt2 = lit(-simd.m_sqrt2); // -2 sin(pi / 4)
    const two: f32 = 2;
    const l1ido = l1 * ido;
    {
        var c = cc;
        var h = ch;
        const end = ch + l1ido;
        while (@intFromPtr(h) < @intFromPtr(end)) {
            const a = c[0];
            const b = c[4 * ido - 1];
            const cv = c[2 * ido];
            const d = c[2 * ido - 1];
            const tr3 = smul(two, d);
            const tr2 = a + b;
            const tr1 = a - b;
            const tr4 = smul(two, cv);
            h[0 * l1ido] = tr2 + tr3;
            h[2 * l1ido] = tr2 - tr3;
            h[1 * l1ido] = tr1 - tr4;
            h[3 * l1ido] = tr1 + tr4;
            c += 4 * ido;
            h += ido;
        }
    }
    if (ido < 2) return;
    var k: usize = 0;
    if (ido == 3) {
        // the inner loop's one pass, i = 2 and ic = 1
        const wa1_r = ld(wa1[0]);
        const wa1_i = ld(wa1[1]);
        const wa2_r = ld(wa2[0]);
        const wa2_i = ld(wa2[1]);
        const wa3_r = ld(wa3[0]);
        const wa3_i = ld(wa3[1]);
        while (k < l1ido) : (k += ido) {
            const tr1 = cc[1 + 4 * k] - cc[9 + 4 * k];
            const tr2 = cc[1 + 4 * k] + cc[9 + 4 * k];
            const ti4 = cc[7 + 4 * k] - cc[3 + 4 * k];
            const tr3 = cc[7 + 4 * k] + cc[3 + 4 * k];
            ch[k + 1] = tr2 + tr3;
            var cr3 = tr2 - tr3;

            const ti3 = cc[8 + 4 * k] - cc[4 + 4 * k];
            const tr4 = cc[8 + 4 * k] + cc[4 + 4 * k];
            var cr2 = tr1 - tr4;
            var cr4 = tr1 + tr4;

            const ti1 = cc[2 + 4 * k] + cc[10 + 4 * k];
            const ti2 = cc[2 + 4 * k] - cc[10 + 4 * k];

            ch[k + 2] = ti2 + ti3;
            var ci3 = ti2 - ti3;
            var ci2 = ti1 + ti4;
            var ci4 = ti1 - ti4;
            cplxMul(&cr2, &ci2, wa1_r, wa1_i);
            ch[k + 1 + l1ido] = cr2;
            ch[k + 2 + l1ido] = ci2;
            cplxMul(&cr3, &ci3, wa2_r, wa2_i);
            ch[k + 1 + 2 * l1ido] = cr3;
            ch[k + 2 + 2 * l1ido] = ci3;
            cplxMul(&cr4, &ci4, wa3_r, wa3_i);
            ch[k + 1 + 3 * l1ido] = cr4;
            ch[k + 2 + 3 * l1ido] = ci4;
        }
        return;
    }
    if (ido != 2) {
        while (k < l1ido) : (k += ido) {
            // pffft's `cc - 1 + 4 k`, one before `cc`'s start at k = 0
            const pc = cc + 4 * k;
            var ph = ch + k + 1;
            var i: usize = 2;
            while (i < ido) : (i += 2) {
                const tr1 = pc[i - 1] - pc[4 * ido - i - 1];
                const tr2 = pc[i - 1] + pc[4 * ido - i - 1];
                const ti4 = pc[2 * ido + i - 1] - pc[2 * ido - i - 1];
                const tr3 = pc[2 * ido + i - 1] + pc[2 * ido - i - 1];
                ph[0] = tr2 + tr3;
                var cr3 = tr2 - tr3;

                const ti3 = pc[2 * ido + i] - pc[2 * ido - i];
                const tr4 = pc[2 * ido + i] + pc[2 * ido - i];
                var cr2 = tr1 - tr4;
                var cr4 = tr1 + tr4;

                const ti1 = pc[i] + pc[4 * ido - i];
                const ti2 = pc[i] - pc[4 * ido - i];

                ph[1] = ti2 + ti3;
                ph += l1ido;
                var ci3 = ti2 - ti3;
                var ci2 = ti1 + ti4;
                var ci4 = ti1 - ti4;
                cplxMul(&cr2, &ci2, ld(wa1[i - 2]), ld(wa1[i - 1]));
                ph[0] = cr2;
                ph[1] = ci2;
                ph += l1ido;
                cplxMul(&cr3, &ci3, ld(wa2[i - 2]), ld(wa2[i - 1]));
                ph[0] = cr3;
                ph[1] = ci3;
                ph += l1ido;
                cplxMul(&cr4, &ci4, ld(wa3[i - 2]), ld(wa3[i - 1]));
                ph[0] = cr4;
                ph[1] = ci4;
                ph = ph - 3 * l1ido + 2;
            }
        }
        if (ido % 2 == 1) return;
    }
    k = 0;
    while (k < l1ido) : (k += ido) {
        const ia = 4 * k + ido;
        const c = cc[ia - 1];
        const d = cc[ia + 2 * ido - 1];
        const a = cc[ia + 0];
        const b = cc[ia + 2 * ido + 0];
        const tr1 = c - d;
        const tr2 = c + d;
        const ti1 = b + a;
        const ti2 = b - a;
        ch[ido - 1 + k + 0 * l1ido] = tr2 + tr2;
        ch[ido - 1 + k + 1 * l1ido] = smul(minus_sqrt2, ti1 - tr1);
        ch[ido - 1 + k + 2 * l1ido] = ti2 + ti2;
        ch[ido - 1 + k + 3 * l1ido] = smul(minus_sqrt2, ti1 + tr1);
    }
}

fn radb5(ido: usize, l1: usize, cc: ConstPtr, ch: Ptr, wa1: [*]const f32, wa2: [*]const f32, wa3: [*]const f32, wa4: [*]const f32) void {
    const tr11 = lit(0.30901699437494742410); // cos(2 pi / 5)
    const ti11 = lit(0.95105651629515357212); // sin(2 pi / 5)
    const tr12 = lit(-0.80901699437494742410); // cos(4 pi / 5)
    const ti12 = lit(0.58778525229247312917); // sin(4 pi / 5)

    // pffft's 1-based `cc_ref(a_1, a_2, a_3)` and `ch_ref`
    const R = struct {
        ido: usize,
        l1: usize,
        inline fn c(r: @This(), a1: usize, a2: usize, a3: usize) usize {
            return ((a3 - 1) * 5 + (a2 - 1)) * r.ido + (a1 - 1);
        }
        inline fn h(r: @This(), a1: usize, a2: usize, a3: usize) usize {
            return ((a3 - 1) * r.l1 + (a2 - 1)) * r.ido + (a1 - 1);
        }
    };
    const r: R = .{ .ido = ido, .l1 = l1 };

    for (1..l1 + 1) |k| {
        const ti5 = cc[r.c(1, 3, k)] + cc[r.c(1, 3, k)];
        const ti4 = cc[r.c(1, 5, k)] + cc[r.c(1, 5, k)];
        const tr2 = cc[r.c(ido, 2, k)] + cc[r.c(ido, 2, k)];
        const tr3 = cc[r.c(ido, 4, k)] + cc[r.c(ido, 4, k)];
        ch[r.h(1, k, 1)] = cc[r.c(1, 1, k)] + (tr2 + tr3);
        const cr2 = cc[r.c(1, 1, k)] + smadd(tr11, tr2, smul(tr12, tr3));
        const cr3 = cc[r.c(1, 1, k)] + smadd(tr12, tr2, smul(tr11, tr3));
        const ci5 = smadd(ti11, ti5, smul(ti12, ti4));
        const ci4 = smsub(ti11, ti4, smul(ti12, ti5));
        ch[r.h(1, k, 2)] = cr2 - ci5;
        ch[r.h(1, k, 3)] = cr3 - ci4;
        ch[r.h(1, k, 4)] = cr3 + ci4;
        ch[r.h(1, k, 5)] = cr2 + ci5;
    }
    if (ido == 1) return;
    const idp2 = ido + 2;
    for (1..l1 + 1) |k| {
        var i: usize = 3;
        while (i <= ido) : (i += 2) {
            const ic = idp2 - i;
            const ti5 = cc[r.c(i, 3, k)] + cc[r.c(ic, 2, k)];
            const ti2 = cc[r.c(i, 3, k)] - cc[r.c(ic, 2, k)];
            const ti4 = cc[r.c(i, 5, k)] + cc[r.c(ic, 4, k)];
            const ti3 = cc[r.c(i, 5, k)] - cc[r.c(ic, 4, k)];
            const tr5 = cc[r.c(i - 1, 3, k)] - cc[r.c(ic - 1, 2, k)];
            const tr2 = cc[r.c(i - 1, 3, k)] + cc[r.c(ic - 1, 2, k)];
            const tr4 = cc[r.c(i - 1, 5, k)] - cc[r.c(ic - 1, 4, k)];
            const tr3 = cc[r.c(i - 1, 5, k)] + cc[r.c(ic - 1, 4, k)];
            ch[r.h(i - 1, k, 1)] = cc[r.c(i - 1, 1, k)] + (tr2 + tr3);
            ch[r.h(i, k, 1)] = cc[r.c(i, 1, k)] + (ti2 + ti3);
            const cr2 = cc[r.c(i - 1, 1, k)] + smadd(tr11, tr2, smul(tr12, tr3));
            const ci2 = cc[r.c(i, 1, k)] + smadd(tr11, ti2, smul(tr12, ti3));
            const cr3 = cc[r.c(i - 1, 1, k)] + smadd(tr12, tr2, smul(tr11, tr3));
            const ci3 = cc[r.c(i, 1, k)] + smadd(tr12, ti2, smul(tr11, ti3));
            const cr5 = smadd(ti11, tr5, smul(ti12, tr4));
            const ci5 = smadd(ti11, ti5, smul(ti12, ti4));
            const cr4 = smsub(ti11, tr4, smul(ti12, tr5));
            const ci4 = smsub(ti11, ti4, smul(ti12, ti5));
            var dr3 = cr3 - ci4;
            var dr4 = cr3 + ci4;
            var di3 = ci3 + cr4;
            var di4 = ci3 - cr4;
            var dr5 = cr2 + ci5;
            var dr2 = cr2 - ci5;
            var di5 = ci2 - cr5;
            var di2 = ci2 + cr5;
            cplxMul(&dr2, &di2, ld(wa1[i - 3]), ld(wa1[i - 2]));
            cplxMul(&dr3, &di3, ld(wa2[i - 3]), ld(wa2[i - 2]));
            cplxMul(&dr4, &di4, ld(wa3[i - 3]), ld(wa3[i - 2]));
            cplxMul(&dr5, &di5, ld(wa4[i - 3]), ld(wa4[i - 2]));

            ch[r.h(i - 1, k, 2)] = dr2;
            ch[r.h(i, k, 2)] = di2;
            ch[r.h(i - 1, k, 3)] = dr3;
            ch[r.h(i, k, 3)] = di3;
            ch[r.h(i - 1, k, 4)] = dr4;
            ch[r.h(i, k, 4)] = di4;
            ch[r.h(i - 1, k, 5)] = dr5;
            ch[r.h(i, k, 5)] = di5;
        }
    }
}

/// `rfftb1_ps`: the backward passes over `n` vectors, ping-ponging between
/// `work1` and `work2`; returns the buffer holding the result.
pub noinline fn rfftb1(n: usize, input: ConstPtr, work1: Ptr, work2: Ptr, wa: [*]const f32, ifac: []const usize) ConstPtr {
    var in: ConstPtr = input;
    var out: Ptr = if (in == work2) work1 else work2;
    const nf = ifac[1];
    var l1: usize = 1;
    var iw: usize = 0;
    std.debug.assert(in != out);
    for (1..nf + 1) |k1| {
        const ip = ifac[k1 + 1];
        const l2 = ip * l1;
        const ido = n / l2;
        switch (ip) {
            5 => {
                const ix2 = iw + ido;
                const ix3 = ix2 + ido;
                const ix4 = ix3 + ido;
                radb5(ido, l1, in, out, wa + iw, wa + ix2, wa + ix3, wa + ix4);
            },
            4 => {
                const ix2 = iw + ido;
                const ix3 = ix2 + ido;
                radb4(ido, l1, in, out, wa + iw, wa + ix2, wa + ix3);
            },
            3 => {
                const ix2 = iw + ido;
                radb3(ido, l1, in, out, wa + iw, wa + ix2);
            },
            2 => radb2(ido, l1, in, out, wa + iw),
            else => unreachable,
        }
        l1 = l2;
        iw += (ip - 1) * ido;
        if (out == work2) {
            out = work1;
            in = work2;
        } else {
            out = work2;
            in = work1;
        }
    }
    return in;
}
