//! `rfftp`'s forward (real to halfcomplex) passes from PocketFFT
//! (`pocketfft_hdronly.h`, cpp branch, commit 4d865cba, BSD-3, LICENSE):
//! `radf2`, `radf3`, `radf4`, `radf5`. `cc` in, `ch` out, `wa` the twiddles
//! as interleaved real and imaginary parts.

const ix = @import("ix.zig");

pub fn Forward(comptime T: type) type {
    return struct {
        const taur: T = -0.5;
        const taui: T = 0.8660254037844386467637231707529362;
        const hsqt2: T = 0.707106781186547524400844362104849;
        const tr11: T = 0.3090169943749474241022934171828191;
        const ti11: T = 0.9510565162951535721164393333793821;
        const tr12: T = -0.8090169943749474241022934171828191;
        const ti12: T = 0.5877852522924731291687059546390728;

        pub fn radf2(ido: usize, l1: usize, cc: [*]const T, ch: [*]T, wa: [*]const T) void {
            const x: ix.Forward(T) = .{ .ido = ido, .l1 = l1, .cdim = 2, .cc = cc, .ch = ch, .wa = wa };
            for (0..l1) |k| ix.pm(T, x.CH(0, 0, k), x.CH(ido - 1, 1, k), x.CC(0, k, 0), x.CC(0, k, 1));
            if (ido & 1 == 0) for (0..l1) |k| {
                x.CH(0, 1, k).* = -x.CC(ido - 1, k, 1);
                x.CH(ido - 1, 0, k).* = x.CC(ido - 1, k, 0);
            };
            if (ido <= 2) return;
            for (0..l1) |k| {
                var i: usize = 2;
                while (i < ido) : (i += 2) {
                    const ic = ido - i;
                    const m = ix.mulpm(T, x.WA(0, i - 2), x.WA(0, i - 1), x.CC(i - 1, k, 1), x.CC(i, k, 1));
                    ix.pm(T, x.CH(i - 1, 0, k), x.CH(ic - 1, 1, k), x.CC(i - 1, k, 0), m.a);
                    ix.pm(T, x.CH(i, 0, k), x.CH(ic, 1, k), m.b, x.CC(i, k, 0));
                }
            }
        }

        pub fn radf3(ido: usize, l1: usize, cc: [*]const T, ch: [*]T, wa: [*]const T) void {
            const x: ix.Forward(T) = .{ .ido = ido, .l1 = l1, .cdim = 3, .cc = cc, .ch = ch, .wa = wa };
            for (0..l1) |k| {
                const cr2 = x.CC(0, k, 1) + x.CC(0, k, 2);
                x.CH(0, 0, k).* = x.CC(0, k, 0) + cr2;
                x.CH(0, 2, k).* = taui * (x.CC(0, k, 2) - x.CC(0, k, 1));
                x.CH(ido - 1, 1, k).* = x.CC(0, k, 0) + taur * cr2;
            }
            if (ido == 1) return;
            for (0..l1) |k| {
                var i: usize = 2;
                while (i < ido) : (i += 2) {
                    const ic = ido - i;
                    const d2 = ix.mulpm(T, x.WA(0, i - 2), x.WA(0, i - 1), x.CC(i - 1, k, 1), x.CC(i, k, 1));
                    const d3 = ix.mulpm(T, x.WA(1, i - 2), x.WA(1, i - 1), x.CC(i - 1, k, 2), x.CC(i, k, 2));
                    const r = ix.rearrange(T, d2.a, d2.b, d3.a, d3.b);
                    x.CH(i - 1, 0, k).* = x.CC(i - 1, k, 0) + r.rx;
                    x.CH(i, 0, k).* = x.CC(i, k, 0) + r.ix;
                    const tr2 = x.CC(i - 1, k, 0) + taur * r.rx;
                    const ti2 = x.CC(i, k, 0) + taur * r.ix;
                    const tr3 = taui * r.ry;
                    const ti3 = taui * r.iy;
                    ix.pm(T, x.CH(i - 1, 2, k), x.CH(ic - 1, 1, k), tr2, tr3);
                    ix.pm(T, x.CH(i, 2, k), x.CH(ic, 1, k), ti3, ti2);
                }
            }
        }

        pub fn radf4(ido: usize, l1: usize, cc: [*]const T, ch: [*]T, wa: [*]const T) void {
            const x: ix.Forward(T) = .{ .ido = ido, .l1 = l1, .cdim = 4, .cc = cc, .ch = ch, .wa = wa };
            for (0..l1) |k| {
                var tr1: T = undefined;
                var tr2: T = undefined;
                ix.pm(T, &tr1, x.CH(0, 2, k), x.CC(0, k, 3), x.CC(0, k, 1));
                ix.pm(T, &tr2, x.CH(ido - 1, 1, k), x.CC(0, k, 0), x.CC(0, k, 2));
                ix.pm(T, x.CH(0, 0, k), x.CH(ido - 1, 3, k), tr2, tr1);
            }
            if (ido & 1 == 0) for (0..l1) |k| {
                const ti1 = -hsqt2 * (x.CC(ido - 1, k, 1) + x.CC(ido - 1, k, 3));
                const tr1 = hsqt2 * (x.CC(ido - 1, k, 1) - x.CC(ido - 1, k, 3));
                ix.pm(T, x.CH(ido - 1, 0, k), x.CH(ido - 1, 2, k), x.CC(ido - 1, k, 0), tr1);
                ix.pm(T, x.CH(0, 3, k), x.CH(0, 1, k), ti1, x.CC(ido - 1, k, 2));
            };
            if (ido <= 2) return;
            for (0..l1) |k| {
                var i: usize = 2;
                while (i < ido) : (i += 2) {
                    const ic = ido - i;
                    const c2 = ix.mulpm(T, x.WA(0, i - 2), x.WA(0, i - 1), x.CC(i - 1, k, 1), x.CC(i, k, 1));
                    const c3 = ix.mulpm(T, x.WA(1, i - 2), x.WA(1, i - 1), x.CC(i - 1, k, 2), x.CC(i, k, 2));
                    const c4 = ix.mulpm(T, x.WA(2, i - 2), x.WA(2, i - 1), x.CC(i - 1, k, 3), x.CC(i, k, 3));
                    const tr1 = c4.a + c2.a;
                    const tr4 = c4.a - c2.a;
                    const ti1 = c2.b + c4.b;
                    const ti4 = c2.b - c4.b;
                    const tr2 = x.CC(i - 1, k, 0) + c3.a;
                    const tr3 = x.CC(i - 1, k, 0) - c3.a;
                    const ti2 = x.CC(i, k, 0) + c3.b;
                    const ti3 = x.CC(i, k, 0) - c3.b;
                    ix.pm(T, x.CH(i - 1, 0, k), x.CH(ic - 1, 3, k), tr2, tr1);
                    ix.pm(T, x.CH(i, 0, k), x.CH(ic, 3, k), ti1, ti2);
                    ix.pm(T, x.CH(i - 1, 2, k), x.CH(ic - 1, 1, k), tr3, ti4);
                    ix.pm(T, x.CH(i, 2, k), x.CH(ic, 1, k), tr4, ti3);
                }
            }
        }

        pub fn radf5(ido: usize, l1: usize, cc: [*]const T, ch: [*]T, wa: [*]const T) void {
            const x: ix.Forward(T) = .{ .ido = ido, .l1 = l1, .cdim = 5, .cc = cc, .ch = ch, .wa = wa };
            for (0..l1) |k| {
                const cr2 = x.CC(0, k, 4) + x.CC(0, k, 1);
                const ci5 = x.CC(0, k, 4) - x.CC(0, k, 1);
                const cr3 = x.CC(0, k, 3) + x.CC(0, k, 2);
                const ci4 = x.CC(0, k, 3) - x.CC(0, k, 2);
                x.CH(0, 0, k).* = x.CC(0, k, 0) + cr2 + cr3;
                x.CH(ido - 1, 1, k).* = x.CC(0, k, 0) + tr11 * cr2 + tr12 * cr3;
                x.CH(0, 2, k).* = ti11 * ci5 + ti12 * ci4;
                x.CH(ido - 1, 3, k).* = x.CC(0, k, 0) + tr12 * cr2 + tr11 * cr3;
                x.CH(0, 4, k).* = ti12 * ci5 - ti11 * ci4;
            }
            if (ido == 1) return;
            for (0..l1) |k| {
                var i: usize = 2;
                while (i < ido) : (i += 2) {
                    const ic = ido - i;
                    const d2 = ix.mulpm(T, x.WA(0, i - 2), x.WA(0, i - 1), x.CC(i - 1, k, 1), x.CC(i, k, 1));
                    const d3 = ix.mulpm(T, x.WA(1, i - 2), x.WA(1, i - 1), x.CC(i - 1, k, 2), x.CC(i, k, 2));
                    const d4 = ix.mulpm(T, x.WA(2, i - 2), x.WA(2, i - 1), x.CC(i - 1, k, 3), x.CC(i, k, 3));
                    const d5 = ix.mulpm(T, x.WA(3, i - 2), x.WA(3, i - 1), x.CC(i - 1, k, 4), x.CC(i, k, 4));
                    const r25 = ix.rearrange(T, d2.a, d2.b, d5.a, d5.b);
                    const r34 = ix.rearrange(T, d3.a, d3.b, d4.a, d4.b);
                    const dr2 = r25.rx;
                    const di2 = r25.ix;
                    const dr5 = r25.ry;
                    const di5 = r25.iy;
                    const dr3 = r34.rx;
                    const di3 = r34.ix;
                    const dr4 = r34.ry;
                    const di4 = r34.iy;
                    x.CH(i - 1, 0, k).* = x.CC(i - 1, k, 0) + dr2 + dr3;
                    x.CH(i, 0, k).* = x.CC(i, k, 0) + di2 + di3;
                    const tr2 = x.CC(i - 1, k, 0) + tr11 * dr2 + tr12 * dr3;
                    const ti2 = x.CC(i, k, 0) + tr11 * di2 + tr12 * di3;
                    const tr3 = x.CC(i - 1, k, 0) + tr12 * dr2 + tr11 * dr3;
                    const ti3 = x.CC(i, k, 0) + tr12 * di2 + tr11 * di3;
                    const tr5 = ti11 * dr5 + ti12 * dr4;
                    const ti5 = ti11 * di5 + ti12 * di4;
                    const tr4 = ti12 * dr5 - ti11 * dr4;
                    const ti4 = ti12 * di5 - ti11 * di4;
                    ix.pm(T, x.CH(i - 1, 2, k), x.CH(ic - 1, 1, k), tr2, tr5);
                    ix.pm(T, x.CH(i, 2, k), x.CH(ic, 1, k), ti5, ti2);
                    ix.pm(T, x.CH(i - 1, 4, k), x.CH(ic - 1, 3, k), tr3, tr4);
                    ix.pm(T, x.CH(i, 4, k), x.CH(ic, 3, k), ti4, ti3);
                }
            }
        }
    };
}
