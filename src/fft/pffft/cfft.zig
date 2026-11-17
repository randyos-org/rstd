//! pffft's complex passes, `passf2_ps` to `passf8_ps` and `cfftf1_ps`
//! (`src/pffft_priv_impl.h`, marton78/pffft efc1dd2d, UCAR, LICENSE), in
//! its operation order. `fsign` is -1 forward, +1 backward.

const std = @import("std");
const simd = @import("simd.zig");

const V = simd.V;
const Ptr = simd.Ptr;
const ConstPtr = simd.ConstPtr;
const ld = simd.ld;
const smul = simd.smul;
const smadd = simd.smadd;
const smsub = simd.smsub;
const cplxMul = simd.cplxMul;
const lit = simd.lit;
const m_sqrt1_2 = simd.m_sqrt1_2;

noinline fn passf2(ido: usize, l1: usize, cc_: ConstPtr, ch_: Ptr, wa1: [*]const f32, fsign: f32) void {
    var cc = cc_;
    var ch = ch_;
    const l1ido = l1 * ido;
    var k: usize = 0;
    if (ido <= 2) {
        while (k < l1ido) : ({
            k += ido;
            ch += ido;
            cc += 2 * ido;
        }) {
            ch[0] = cc[0] + cc[ido + 0];
            ch[l1ido] = cc[0] - cc[ido + 0];
            ch[1] = cc[1] + cc[ido + 1];
            ch[l1ido + 1] = cc[1] - cc[ido + 1];
        }
    } else {
        while (k < l1ido) : ({
            k += ido;
            ch += ido;
            cc += 2 * ido;
        }) {
            var i: usize = 0;
            while (i < ido - 1) : (i += 2) {
                var tr2 = cc[i + 0] - cc[i + ido + 0];
                var ti2 = cc[i + 1] - cc[i + ido + 1];
                const wr = ld(wa1[i]);
                const wi = ld(fsign) * ld(wa1[i + 1]);
                ch[i] = cc[i + 0] + cc[i + ido + 0];
                ch[i + 1] = cc[i + 1] + cc[i + ido + 1];
                cplxMul(&tr2, &ti2, wr, wi);
                ch[i + l1ido] = tr2;
                ch[i + l1ido + 1] = ti2;
            }
        }
    }
}

noinline fn passf3(ido: usize, l1: usize, cc_: ConstPtr, ch_: Ptr, wa1: [*]const f32, wa2: [*]const f32, fsign: f32) void {
    var cc = cc_;
    var ch = ch_;
    const taur: f32 = -0.5; // cos(2 pi / 3)
    const taui: f32 = lit(0.86602540378443864676) * fsign; // sin(2 pi / 3)
    const l1ido = l1 * ido;
    std.debug.assert(ido > 2);
    var k: usize = 0;
    while (k < l1ido) : ({
        k += ido;
        cc += 3 * ido;
        ch += ido;
    }) {
        var i: usize = 0;
        while (i < ido - 1) : (i += 2) {
            const tr2 = cc[i + ido] + cc[i + 2 * ido];
            const cr2 = smadd(taur, tr2, cc[i]);
            ch[i] = cc[i] + tr2;
            const ti2 = cc[i + ido + 1] + cc[i + 2 * ido + 1];
            const ci2 = smadd(taur, ti2, cc[i + 1]);
            ch[i + 1] = cc[i + 1] + ti2;
            const cr3 = smul(taui, cc[i + ido] - cc[i + 2 * ido]);
            const ci3 = smul(taui, cc[i + ido + 1] - cc[i + 2 * ido + 1]);
            var dr2 = cr2 - ci3;
            var dr3 = cr2 + ci3;
            var di2 = ci2 + cr3;
            var di3 = ci2 - cr3;
            const wr1 = wa1[i];
            const wi1 = fsign * wa1[i + 1];
            const wr2 = wa2[i];
            const wi2 = fsign * wa2[i + 1];
            cplxMul(&dr2, &di2, ld(wr1), ld(wi1));
            ch[i + l1ido] = dr2;
            ch[i + l1ido + 1] = di2;
            cplxMul(&dr3, &di3, ld(wr2), ld(wi2));
            ch[i + 2 * l1ido] = dr3;
            ch[i + 2 * l1ido + 1] = di3;
        }
    }
}

noinline fn passf4(ido: usize, l1: usize, cc_: ConstPtr, ch_: Ptr, wa1: [*]const f32, wa2: [*]const f32, wa3: [*]const f32, fsign: f32) void {
    var cc = cc_;
    var ch = ch_;
    const l1ido = l1 * ido;
    var k: usize = 0;
    if (ido == 2) {
        while (k < l1ido) : ({
            k += ido;
            ch += ido;
            cc += 4 * ido;
        }) {
            const tr1 = cc[0] - cc[2 * ido + 0];
            const tr2 = cc[0] + cc[2 * ido + 0];
            const ti1 = cc[1] - cc[2 * ido + 1];
            const ti2 = cc[1] + cc[2 * ido + 1];
            const ti4 = (cc[1 * ido + 0] - cc[3 * ido + 0]) * ld(fsign);
            const tr4 = (cc[3 * ido + 1] - cc[1 * ido + 1]) * ld(fsign);
            const tr3 = cc[ido + 0] + cc[3 * ido + 0];
            const ti3 = cc[ido + 1] + cc[3 * ido + 1];

            ch[0 * l1ido + 0] = tr2 + tr3;
            ch[0 * l1ido + 1] = ti2 + ti3;
            ch[1 * l1ido + 0] = tr1 + tr4;
            ch[1 * l1ido + 1] = ti1 + ti4;
            ch[2 * l1ido + 0] = tr2 - tr3;
            ch[2 * l1ido + 1] = ti2 - ti3;
            ch[3 * l1ido + 0] = tr1 - tr4;
            ch[3 * l1ido + 1] = ti1 - ti4;
        }
    } else {
        while (k < l1ido) : ({
            k += ido;
            ch += ido;
            cc += 4 * ido;
        }) {
            var i: usize = 0;
            while (i < ido - 1) : (i += 2) {
                const tr1 = cc[i + 0] - cc[i + 2 * ido + 0];
                const tr2 = cc[i + 0] + cc[i + 2 * ido + 0];
                const ti1 = cc[i + 1] - cc[i + 2 * ido + 1];
                const ti2 = cc[i + 1] + cc[i + 2 * ido + 1];
                const tr4 = (cc[i + 3 * ido + 1] - cc[i + 1 * ido + 1]) * ld(fsign);
                const ti4 = (cc[i + 1 * ido + 0] - cc[i + 3 * ido + 0]) * ld(fsign);
                const tr3 = cc[i + ido + 0] + cc[i + 3 * ido + 0];
                const ti3 = cc[i + ido + 1] + cc[i + 3 * ido + 1];

                ch[i] = tr2 + tr3;
                var cr3 = tr2 - tr3;
                ch[i + 1] = ti2 + ti3;
                var ci3 = ti2 - ti3;

                var cr2 = tr1 + tr4;
                var cr4 = tr1 - tr4;
                var ci2 = ti1 + ti4;
                var ci4 = ti1 - ti4;
                const wr1 = wa1[i];
                const wi1 = fsign * wa1[i + 1];
                cplxMul(&cr2, &ci2, ld(wr1), ld(wi1));
                const wr2 = wa2[i];
                const wi2 = fsign * wa2[i + 1];
                ch[i + l1ido] = cr2;
                ch[i + l1ido + 1] = ci2;

                cplxMul(&cr3, &ci3, ld(wr2), ld(wi2));
                const wr3 = wa3[i];
                const wi3 = fsign * wa3[i + 1];
                ch[i + 2 * l1ido] = cr3;
                ch[i + 2 * l1ido + 1] = ci3;

                cplxMul(&cr4, &ci4, ld(wr3), ld(wi3));
                ch[i + 3 * l1ido] = cr4;
                ch[i + 3 * l1ido + 1] = ci4;
            }
        }
    }
}

noinline fn passf5(ido: usize, l1: usize, cc_: ConstPtr, ch_: Ptr, wa1: [*]const f32, wa2: [*]const f32, wa3: [*]const f32, wa4: [*]const f32, fsign: f32) void {
    var cc = cc_;
    var ch = ch_;
    const tr11 = lit(0.30901699437494742410); // cos(2 pi / 5)
    const ti11 = lit(0.95105651629515357212) * fsign; // sin(2 pi / 5)
    const tr12 = lit(-0.80901699437494742410); // cos(4 pi / 5)
    const ti12 = lit(0.58778525229247312917) * fsign; // sin(4 pi / 5)
    const l1ido = l1 * ido;

    std.debug.assert(ido > 2);
    for (0..l1) |_| {
        var i: usize = 0;
        while (i < ido - 1) : (i += 2) {
            // `cc_ref(i - 1, a)` is `cc[(a - 1) ido + i]`, `cc_ref(i, a)` the next; `ch_ref` too
            const ti5 = cc[ido + i + 1] - cc[4 * ido + i + 1];
            const ti2 = cc[ido + i + 1] + cc[4 * ido + i + 1];
            const ti4 = cc[2 * ido + i + 1] - cc[3 * ido + i + 1];
            const ti3 = cc[2 * ido + i + 1] + cc[3 * ido + i + 1];
            const tr5 = cc[ido + i] - cc[4 * ido + i];
            const tr2 = cc[ido + i] + cc[4 * ido + i];
            const tr4 = cc[2 * ido + i] - cc[3 * ido + i];
            const tr3 = cc[2 * ido + i] + cc[3 * ido + i];
            ch[i] = cc[i] + (tr2 + tr3);
            ch[i + 1] = cc[i + 1] + (ti2 + ti3);
            const cr2 = cc[i] + smadd(tr11, tr2, smul(tr12, tr3));
            const ci2 = cc[i + 1] + smadd(tr11, ti2, smul(tr12, ti3));
            const cr3 = cc[i] + smadd(tr12, tr2, smul(tr11, tr3));
            const ci3 = cc[i + 1] + smadd(tr12, ti2, smul(tr11, ti3));
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
            const wr1 = wa1[i];
            const wi1 = fsign * wa1[i + 1];
            const wr2 = wa2[i];
            const wi2 = fsign * wa2[i + 1];
            const wr3 = wa3[i];
            const wi3 = fsign * wa3[i + 1];
            const wr4 = wa4[i];
            const wi4 = fsign * wa4[i + 1];
            cplxMul(&dr2, &di2, ld(wr1), ld(wi1));
            ch[l1ido + i] = dr2;
            ch[l1ido + i + 1] = di2;
            cplxMul(&dr3, &di3, ld(wr2), ld(wi2));
            ch[2 * l1ido + i] = dr3;
            ch[2 * l1ido + i + 1] = di3;
            cplxMul(&dr4, &di4, ld(wr3), ld(wi3));
            ch[3 * l1ido + i] = dr4;
            ch[3 * l1ido + i + 1] = di4;
            cplxMul(&dr5, &di5, ld(wr4), ld(wi4));
            ch[4 * l1ido + i] = dr5;
            ch[4 * l1ido + i + 1] = di5;
        }
        cc += 5 * ido;
        ch += ido;
    }
}

/// Only `n == 8` takes radix 8 (`cffti1_ps`), so `ido` is 2 and `l1` 1.
noinline fn passf8(ido: usize, l1: usize, cc: ConstPtr, ch: Ptr, fsign: f32) void {
    const k8 = lit(m_sqrt1_2); // cos(pi / 4) = sin(pi / 4)
    const k8s = fsign * k8;
    const nk8 = -k8;
    const nfsign = -fsign;

    std.debug.assert(ido == 2 and l1 == 1);

    const x0r = cc[0 * ido + 0];
    const x0i = cc[0 * ido + 1];
    const x1r = cc[1 * ido + 0];
    const x1i = cc[1 * ido + 1];
    const x2r = cc[2 * ido + 0];
    const x2i = cc[2 * ido + 1];
    const x3r = cc[3 * ido + 0];
    const x3i = cc[3 * ido + 1];
    const x4r = cc[4 * ido + 0];
    const x4i = cc[4 * ido + 1];
    const x5r = cc[5 * ido + 0];
    const x5i = cc[5 * ido + 1];
    const x6r = cc[6 * ido + 0];
    const x6i = cc[6 * ido + 1];
    const x7r = cc[7 * ido + 0];
    const x7i = cc[7 * ido + 1];

    // sums and differences at stride 4
    const s0r = x0r + x4r;
    const s0i = x0i + x4i;
    const s1r = x1r + x5r;
    const s1i = x1i + x5i;
    const s2r = x2r + x6r;
    const s2i = x2i + x6i;
    const s3r = x3r + x7r;
    const s3i = x3i + x7i;
    const d0r = x0r - x4r;
    const d0i = x0i - x4i;
    var d1r = x1r - x5r;
    var d1i = x1i - x5i;
    var d2r = x2r - x6r;
    var d2i = x2i - x6i;
    var d3r = x3r - x7r;
    var d3i = x3i - x7i;

    // the differences by W8^k
    cplxMul(&d1r, &d1i, ld(k8), ld(k8s));
    {
        const tr = ld(nfsign) * d2i;
        const ti = ld(fsign) * d2r;
        d2r = tr;
        d2i = ti;
    }
    cplxMul(&d3r, &d3i, ld(nk8), ld(k8s));

    // a DFT-4 of the sums: the even outputs
    {
        const t1r = s0r - s2r;
        const t1i = s0i - s2i;
        const t2r = s0r + s2r;
        const t2i = s0i + s2i;
        const t4i = (s1r - s3r) * ld(fsign);
        const t4r = (s3i - s1i) * ld(fsign);
        const t3r = s1r + s3r;
        const t3i = s1i + s3i;

        ch[0 * ido + 0] = t2r + t3r;
        ch[0 * ido + 1] = t2i + t3i;
        ch[2 * ido + 0] = t1r + t4r;
        ch[2 * ido + 1] = t1i + t4i;
        ch[4 * ido + 0] = t2r - t3r;
        ch[4 * ido + 1] = t2i - t3i;
        ch[6 * ido + 0] = t1r - t4r;
        ch[6 * ido + 1] = t1i - t4i;
    }

    // a DFT-4 of the differences: the odd outputs
    {
        const t1r = d0r - d2r;
        const t1i = d0i - d2i;
        const t2r = d0r + d2r;
        const t2i = d0i + d2i;
        const t4i = (d1r - d3r) * ld(fsign);
        const t4r = (d3i - d1i) * ld(fsign);
        const t3r = d1r + d3r;
        const t3i = d1i + d3i;

        ch[1 * ido + 0] = t2r + t3r;
        ch[1 * ido + 1] = t2i + t3i;
        ch[3 * ido + 0] = t1r + t4r;
        ch[3 * ido + 1] = t1i + t4i;
        ch[5 * ido + 0] = t2r - t3r;
        ch[5 * ido + 1] = t2i - t3i;
        ch[7 * ido + 0] = t1r - t4r;
        ch[7 * ido + 1] = t1i - t4i;
    }
}

/// `cfftf1_ps`: the passes over `n` complex vectors, ping-ponging between
/// `work1` and `work2`; returns the buffer holding the result.
pub fn cfftf1(n: usize, input: ConstPtr, work1: Ptr, work2: Ptr, wa: [*]const f32, ifac: []const usize, isign: f32) ConstPtr {
    var in: ConstPtr = input;
    var out: Ptr = if (in == work2) work1 else work2;
    const nf = ifac[1];
    var l1: usize = 1;
    var iw: usize = 0;
    std.debug.assert(in != out and work1 != work2);
    for (2..nf + 2) |k1| {
        const ip = ifac[k1];
        const l2 = ip * l1;
        const ido = n / l2;
        const idot = ido + ido;
        switch (ip) {
            8 => passf8(idot, l1, in, out, isign),
            5 => {
                const ix2 = iw + idot;
                const ix3 = ix2 + idot;
                const ix4 = ix3 + idot;
                passf5(idot, l1, in, out, wa + iw, wa + ix2, wa + ix3, wa + ix4, isign);
            },
            4 => {
                const ix2 = iw + idot;
                const ix3 = ix2 + idot;
                passf4(idot, l1, in, out, wa + iw, wa + ix2, wa + ix3, isign);
            },
            2 => passf2(idot, l1, in, out, wa + iw, isign),
            3 => {
                const ix2 = iw + idot;
                passf3(idot, l1, in, out, wa + iw, wa + ix2, isign);
            },
            else => unreachable,
        }
        l1 = l2;
        iw += (ip - 1) * idot;
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
