//! pffft's real forward passes, `radf2_ps` to `radf5_ps` and `rfftf1_ps`
//! (`src/pffft_priv_impl.h`, marton78/pffft efc1dd2d, UCAR, LICENSE), in
//! its operation order.

const std = @import("std");
const simd_file = @import("simd.zig");

/// The steps over `lanes`-wide vectors.
pub fn For(comptime lanes: comptime_int) type {
    return struct {
        const simd = simd_file.For(lanes);

        const V = simd.V;
        const Ptr = simd.Ptr;
        const ConstPtr = simd.ConstPtr;
        const ld = simd.ld;
        const lit = simd.lit;
        const smul = simd.smul;
        const smadd = simd.smadd;
        const smsub = simd.smsub;
        const cplxMulConj = simd.cplxMulConj;

        noinline fn radf2(ido: usize, l1: usize, cc: ConstPtr, ch: Ptr, wa1: [*]const f32) void {
            const minus_one: f32 = -1;
            const l1ido = l1 * ido;
            var k: usize = 0;
            while (k < l1ido) : (k += ido) {
                const a = cc[k];
                const b = cc[k + l1ido];
                ch[2 * k] = a + b;
                ch[2 * (k + ido) - 1] = a - b;
            }
            if (ido < 2) return;
            if (ido != 2) {
                k = 0;
                while (k < l1ido) : (k += ido) {
                    var i: usize = 2;
                    while (i < ido) : (i += 2) {
                        var tr2 = cc[i - 1 + k + l1ido];
                        var ti2 = cc[i + k + l1ido];
                        const br = cc[i - 1 + k];
                        const bi = cc[i + k];
                        cplxMulConj(&tr2, &ti2, ld(wa1[i - 2]), ld(wa1[i - 1]));
                        ch[i + 2 * k] = bi + ti2;
                        ch[2 * (k + ido) - i] = ti2 - bi;
                        ch[i - 1 + 2 * k] = br + tr2;
                        ch[2 * (k + ido) - i - 1] = br - tr2;
                    }
                }
                if (ido % 2 == 1) return;
            }
            k = 0;
            while (k < l1ido) : (k += ido) {
                ch[2 * k + ido] = smul(minus_one, cc[ido - 1 + k + l1ido]);
                ch[2 * k + ido - 1] = cc[k + ido - 1];
            }
        }

        fn radf3(ido: usize, l1: usize, cc: ConstPtr, ch: Ptr, wa1: [*]const f32, wa2: [*]const f32) void {
            const taur: f32 = -0.5; // cos(2 pi / 3)
            const taui = lit(0.86602540378443864676); // sin(2 pi / 3)
            for (0..l1) |k| {
                const cr2 = cc[(k + l1) * ido] + cc[(k + 2 * l1) * ido];
                ch[3 * k * ido] = cc[k * ido] + cr2;
                ch[(3 * k + 2) * ido] = smul(taui, cc[(k + l1 * 2) * ido] - cc[(k + l1) * ido]);
                ch[ido - 1 + (3 * k + 1) * ido] = smadd(taur, cr2, cc[k * ido]);
            }
            if (ido == 1) return;
            if (ido == 3) {
                // the inner loop's one pass, i = 2 and ic = 1
                const wr1 = ld(wa1[0]);
                const wi1 = ld(wa1[1]);
                const wr2 = ld(wa2[0]);
                const wi2 = ld(wa2[1]);
                for (0..l1) |k| {
                    var dr2 = cc[1 + (k + l1) * ido];
                    var di2 = cc[2 + (k + l1) * ido];
                    cplxMulConj(&dr2, &di2, wr1, wi1);
                    var dr3 = cc[1 + (k + l1 * 2) * ido];
                    var di3 = cc[2 + (k + l1 * 2) * ido];
                    cplxMulConj(&dr3, &di3, wr2, wi2);
                    const cr2 = dr2 + dr3;
                    const ci2 = di2 + di3;
                    ch[1 + 3 * k * ido] = cc[1 + k * ido] + cr2;
                    ch[2 + 3 * k * ido] = cc[2 + k * ido] + ci2;
                    const tr2 = smadd(taur, cr2, cc[1 + k * ido]);
                    const ti2 = smadd(taur, ci2, cc[2 + k * ido]);
                    const tr3 = smul(taui, di2 - di3);
                    const ti3 = smul(taui, dr3 - dr2);
                    ch[1 + (3 * k + 2) * ido] = tr2 + tr3;
                    ch[0 + (3 * k + 1) * ido] = tr2 - tr3;
                    ch[2 + (3 * k + 2) * ido] = ti2 + ti3;
                    ch[1 + (3 * k + 1) * ido] = ti3 - ti2;
                }
                return;
            }
            for (0..l1) |k| {
                var i: usize = 2;
                while (i < ido) : (i += 2) {
                    const ic = ido - i;
                    const wr1 = ld(wa1[i - 2]);
                    const wi1 = ld(wa1[i - 1]);
                    var dr2 = cc[i - 1 + (k + l1) * ido];
                    var di2 = cc[i + (k + l1) * ido];
                    cplxMulConj(&dr2, &di2, wr1, wi1);

                    const wr2 = ld(wa2[i - 2]);
                    const wi2 = ld(wa2[i - 1]);
                    var dr3 = cc[i - 1 + (k + l1 * 2) * ido];
                    var di3 = cc[i + (k + l1 * 2) * ido];
                    cplxMulConj(&dr3, &di3, wr2, wi2);

                    const cr2 = dr2 + dr3;
                    const ci2 = di2 + di3;
                    ch[i - 1 + 3 * k * ido] = cc[i - 1 + k * ido] + cr2;
                    ch[i + 3 * k * ido] = cc[i + k * ido] + ci2;
                    const tr2 = smadd(taur, cr2, cc[i - 1 + k * ido]);
                    const ti2 = smadd(taur, ci2, cc[i + k * ido]);
                    const tr3 = smul(taui, di2 - di3);
                    const ti3 = smul(taui, dr3 - dr2);
                    ch[i - 1 + (3 * k + 2) * ido] = tr2 + tr3;
                    ch[ic - 1 + (3 * k + 1) * ido] = tr2 - tr3;
                    ch[i + (3 * k + 2) * ido] = ti2 + ti3;
                    ch[ic + (3 * k + 1) * ido] = ti3 - ti2;
                }
            }
        }

        noinline fn radf4(ido: usize, l1: usize, cc: ConstPtr, ch: Ptr, wa1: [*]const f32, wa2: [*]const f32, wa3: [*]const f32) void {
            const minus_hsqt2 = lit(-simd.m_sqrt1_2); // -sin(pi / 4)
            const l1ido = l1 * ido;
            {
                var c = cc;
                var h = ch;
                const end = cc + l1ido;
                while (@intFromPtr(c) < @intFromPtr(end)) {
                    // between 25% and 40% of radf4's cost (pffft)
                    const a0 = c[0];
                    const a1 = c[l1ido];
                    const a2 = c[2 * l1ido];
                    const a3 = c[3 * l1ido];
                    const tr1 = a1 + a3;
                    const tr2 = a0 + a2;
                    h[2 * ido - 1] = a0 - a2;
                    h[2 * ido] = a3 - a1;
                    h[0] = tr1 + tr2;
                    h[4 * ido - 1] = tr2 - tr1;
                    c += ido;
                    h += 4 * ido;
                }
            }
            if (ido < 2) return;
            var k: usize = 0;
            if (ido == 3) {
                // the inner loop's one pass, i = 2 and ic = 1
                while (k < l1ido) : (k += ido) {
                    var cr2 = cc[1 + k + 1 * l1ido];
                    var ci2 = cc[1 + k + 1 * l1ido + 1];
                    var wr = ld(wa1[0]);
                    var wi = ld(wa1[1]);
                    cplxMulConj(&cr2, &ci2, wr, wi);

                    var cr3 = cc[1 + k + 2 * l1ido];
                    var ci3 = cc[1 + k + 2 * l1ido + 1];
                    wr = ld(wa2[0]);
                    wi = ld(wa2[1]);
                    cplxMulConj(&cr3, &ci3, wr, wi);

                    var cr4 = cc[1 + k + 3 * l1ido];
                    var ci4 = cc[1 + k + 3 * l1ido + 1];
                    wr = ld(wa3[0]);
                    wi = ld(wa3[1]);
                    cplxMulConj(&cr4, &ci4, wr, wi);

                    const tr1 = cr2 + cr4;
                    const tr4 = cr4 - cr2;
                    const tr2 = cc[1 + k] + cr3;
                    const tr3 = cc[1 + k] - cr3;
                    ch[1 + 4 * k] = tr1 + tr2;
                    ch[0 + 4 * k + 3 * ido] = tr2 - tr1;
                    const ti1 = ci2 + ci4;
                    const ti4 = ci2 - ci4;
                    ch[1 + 4 * k + 2 * ido] = ti4 + tr3;
                    ch[0 + 4 * k + 1 * ido] = tr3 - ti4;
                    const ti2 = cc[1 + k + 1] + ci3;
                    const ti3 = cc[1 + k + 1] - ci3;
                    ch[2 + 4 * k] = ti1 + ti2;
                    ch[1 + 4 * k + 3 * ido] = ti1 - ti2;
                    ch[2 + 4 * k + 2 * ido] = tr4 + ti3;
                    ch[1 + 4 * k + 1 * ido] = tr4 - ti3;
                }
                return;
            }
            if (ido != 2) {
                while (k < l1ido) : (k += ido) {
                    var pc = cc + 1 + k;
                    var i: usize = 2;
                    while (i < ido) : ({
                        i += 2;
                        pc += 2;
                    }) {
                        const ic = ido - i;

                        var cr2 = pc[1 * l1ido + 0];
                        var ci2 = pc[1 * l1ido + 1];
                        var wr = ld(wa1[i - 2]);
                        var wi = ld(wa1[i - 1]);
                        cplxMulConj(&cr2, &ci2, wr, wi);

                        var cr3 = pc[2 * l1ido + 0];
                        var ci3 = pc[2 * l1ido + 1];
                        wr = ld(wa2[i - 2]);
                        wi = ld(wa2[i - 1]);
                        cplxMulConj(&cr3, &ci3, wr, wi);

                        var cr4 = pc[3 * l1ido];
                        var ci4 = pc[3 * l1ido + 1];
                        wr = ld(wa3[i - 2]);
                        wi = ld(wa3[i - 1]);
                        cplxMulConj(&cr4, &ci4, wr, wi);

                        const tr1 = cr2 + cr4;
                        const tr4 = cr4 - cr2;
                        const tr2 = pc[0] + cr3;
                        const tr3 = pc[0] - cr3;
                        ch[i - 1 + 4 * k] = tr1 + tr2;
                        ch[ic - 1 + 4 * k + 3 * ido] = tr2 - tr1;
                        const ti1 = ci2 + ci4;
                        const ti4 = ci2 - ci4;
                        ch[i - 1 + 4 * k + 2 * ido] = ti4 + tr3;
                        ch[ic - 1 + 4 * k + 1 * ido] = tr3 - ti4;
                        const ti2 = pc[1] + ci3;
                        const ti3 = pc[1] - ci3;
                        ch[i + 4 * k] = ti1 + ti2;
                        ch[ic + 4 * k + 3 * ido] = ti1 - ti2;
                        ch[i + 4 * k + 2 * ido] = tr4 + ti3;
                        ch[ic + 4 * k + 1 * ido] = tr4 - ti3;
                    }
                }
                if (ido % 2 == 1) return;
            }
            k = 0;
            while (k < l1ido) : (k += ido) {
                const a = cc[ido - 1 + k + l1ido];
                const b = cc[ido - 1 + k + 3 * l1ido];
                const c = cc[ido - 1 + k];
                const d = cc[ido - 1 + k + 2 * l1ido];
                const ti1 = smul(minus_hsqt2, a + b);
                const tr1 = smul(minus_hsqt2, b - a);
                ch[ido - 1 + 4 * k] = tr1 + c;
                ch[ido - 1 + 4 * k + 2 * ido] = c - tr1;
                ch[4 * k + 1 * ido] = ti1 - d;
                ch[4 * k + 3 * ido] = ti1 + d;
            }
        }

        fn radf5(ido: usize, l1: usize, cc: ConstPtr, ch: Ptr, wa1: [*]const f32, wa2: [*]const f32, wa3: [*]const f32, wa4: [*]const f32) void {
            const tr11 = lit(0.30901699437494742410); // cos(2 pi / 5)
            const ti11 = lit(0.95105651629515357212); // sin(2 pi / 5)
            const tr12 = lit(-0.80901699437494742410); // cos(4 pi / 5)
            const ti12 = lit(0.58778525229247312917); // sin(4 pi / 5)

            // pffft's 1-based `cc_ref(a_1, a_2, a_3)` and `ch_ref`
            const R = struct {
                ido: usize,
                l1: usize,
                inline fn c(r: @This(), a1: usize, a2: usize, a3: usize) usize {
                    return ((a3 - 1) * r.l1 + (a2 - 1)) * r.ido + (a1 - 1);
                }
                inline fn h(r: @This(), a1: usize, a2: usize, a3: usize) usize {
                    return ((a3 - 1) * 5 + (a2 - 1)) * r.ido + (a1 - 1);
                }
            };
            const r: R = .{ .ido = ido, .l1 = l1 };

            for (1..l1 + 1) |k| {
                const cr2 = cc[r.c(1, k, 5)] + cc[r.c(1, k, 2)];
                const ci5 = cc[r.c(1, k, 5)] - cc[r.c(1, k, 2)];
                const cr3 = cc[r.c(1, k, 4)] + cc[r.c(1, k, 3)];
                const ci4 = cc[r.c(1, k, 4)] - cc[r.c(1, k, 3)];
                ch[r.h(1, 1, k)] = cc[r.c(1, k, 1)] + (cr2 + cr3);
                ch[r.h(ido, 2, k)] = cc[r.c(1, k, 1)] + smadd(tr11, cr2, smul(tr12, cr3));
                ch[r.h(1, 3, k)] = smadd(ti11, ci5, smul(ti12, ci4));
                ch[r.h(ido, 4, k)] = cc[r.c(1, k, 1)] + smadd(tr12, cr2, smul(tr11, cr3));
                ch[r.h(1, 5, k)] = smsub(ti11, ci4, smul(ti12, ci5));
            }
            if (ido == 1) return;
            const idp2 = ido + 2;
            for (1..l1 + 1) |k| {
                var i: usize = 3;
                while (i <= ido) : (i += 2) {
                    const ic = idp2 - i;
                    var dr2 = ld(wa1[i - 3]);
                    var di2 = ld(wa1[i - 2]);
                    var dr3 = ld(wa2[i - 3]);
                    var di3 = ld(wa2[i - 2]);
                    var dr4 = ld(wa3[i - 3]);
                    var di4 = ld(wa3[i - 2]);
                    var dr5 = ld(wa4[i - 3]);
                    var di5 = ld(wa4[i - 2]);
                    cplxMulConj(&dr2, &di2, cc[r.c(i - 1, k, 2)], cc[r.c(i, k, 2)]);
                    cplxMulConj(&dr3, &di3, cc[r.c(i - 1, k, 3)], cc[r.c(i, k, 3)]);
                    cplxMulConj(&dr4, &di4, cc[r.c(i - 1, k, 4)], cc[r.c(i, k, 4)]);
                    cplxMulConj(&dr5, &di5, cc[r.c(i - 1, k, 5)], cc[r.c(i, k, 5)]);
                    const cr2 = dr2 + dr5;
                    const ci5 = dr5 - dr2;
                    const cr5 = di2 - di5;
                    const ci2 = di2 + di5;
                    const cr3 = dr3 + dr4;
                    const ci4 = dr4 - dr3;
                    const cr4 = di3 - di4;
                    const ci3 = di3 + di4;
                    ch[r.h(i - 1, 1, k)] = cc[r.c(i - 1, k, 1)] + (cr2 + cr3);
                    ch[r.h(i, 1, k)] = cc[r.c(i, k, 1)] - (ci2 + ci3);
                    const tr2 = cc[r.c(i - 1, k, 1)] + smadd(tr11, cr2, smul(tr12, cr3));
                    const ti2 = cc[r.c(i, k, 1)] - smadd(tr11, ci2, smul(tr12, ci3));
                    const tr3 = cc[r.c(i - 1, k, 1)] + smadd(tr12, cr2, smul(tr11, cr3));
                    const ti3 = cc[r.c(i, k, 1)] - smadd(tr12, ci2, smul(tr11, ci3));
                    const tr5 = smadd(ti11, cr5, smul(ti12, cr4));
                    const ti5 = smadd(ti11, ci5, smul(ti12, ci4));
                    const tr4 = smsub(ti11, cr4, smul(ti12, cr5));
                    const ti4 = smsub(ti11, ci4, smul(ti12, ci5));
                    ch[r.h(i - 1, 3, k)] = tr2 - tr5;
                    ch[r.h(ic - 1, 2, k)] = tr2 + tr5;
                    ch[r.h(i, 3, k)] = ti2 + ti5;
                    ch[r.h(ic, 2, k)] = ti5 - ti2;
                    ch[r.h(i - 1, 5, k)] = tr3 - tr4;
                    ch[r.h(ic - 1, 4, k)] = tr3 + tr4;
                    ch[r.h(i, 5, k)] = ti3 + ti4;
                    ch[r.h(ic, 4, k)] = ti4 - ti3;
                }
            }
        }

        /// `rfftf1_ps`: the forward passes over `n` vectors, ping-ponging between
        /// `work1` and `work2`; returns the buffer holding the result.
        pub noinline fn rfftf1(n: usize, input: ConstPtr, work1: Ptr, work2: Ptr, wa: [*]const f32, ifac: []const usize) ConstPtr {
            var in: ConstPtr = input;
            var out: Ptr = if (in == work2) work1 else work2;
            const nf = ifac[1];
            var l2 = n;
            var iw = n - 1;
            std.debug.assert(in != out and work1 != work2);
            for (1..nf + 1) |k1| {
                const kh = nf - k1;
                const ip = ifac[kh + 2];
                const l1 = l2 / ip;
                const ido = n / l2;
                iw -= (ip - 1) * ido;
                switch (ip) {
                    5 => {
                        const ix2 = iw + ido;
                        const ix3 = ix2 + ido;
                        const ix4 = ix3 + ido;
                        radf5(ido, l1, in, out, wa + iw, wa + ix2, wa + ix3, wa + ix4);
                    },
                    4 => {
                        const ix2 = iw + ido;
                        const ix3 = ix2 + ido;
                        radf4(ido, l1, in, out, wa + iw, wa + ix2, wa + ix3);
                    },
                    3 => {
                        const ix2 = iw + ido;
                        radf3(ido, l1, in, out, wa + iw, wa + ix2);
                    },
                    2 => radf2(ido, l1, in, out, wa + iw),
                    else => unreachable,
                }
                l2 = l1;
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
    };
}
