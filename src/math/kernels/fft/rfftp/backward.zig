//! `rfftp`'s backward (halfcomplex to real) passes from PocketFFT
//! (`pocketfft_hdronly.h`, cpp branch, commit 4d865cba, BSD-3, LICENSE):
//! `radb2`, `radb3`, `radb4`, `radb5`.

const ix = @import("ix.zig");

pub fn Backward(comptime T: type) type {
    return struct {
        const taur: T = -0.5;
        const taui: T = 0.8660254037844386467637231707529362;
        const sqrt2: T = 1.414213562373095048801688724209698;
        const tr11: T = 0.3090169943749474241022934171828191;
        const ti11: T = 0.9510565162951535721164393333793821;
        const tr12: T = -0.8090169943749474241022934171828191;
        const ti12: T = 0.5877852522924731291687059546390728;

        /// `MULPM` into `CH(i, k, j)` and `CH(i - 1, k, j)`.
        inline fn twiddleOut(x: ix.Backward(T), i: usize, k: usize, j: usize, di: T, dr: T) void {
            const m = ix.mulpm(T, x.WA(j - 1, i - 2), x.WA(j - 1, i - 1), di, dr);
            x.CH(i, k, j).* = m.a;
            x.CH(i - 1, k, j).* = m.b;
        }

        pub fn radb2(ido: usize, l1: usize, cc: [*]const T, ch: [*]T, wa: [*]const T) void {
            const x: ix.Backward(T) = .{ .ido = ido, .l1 = l1, .cdim = 2, .cc = cc, .ch = ch, .wa = wa };
            for (0..l1) |k| ix.pm(T, x.CH(0, k, 0), x.CH(0, k, 1), x.CC(0, 0, k), x.CC(ido - 1, 1, k));
            if (ido & 1 == 0) for (0..l1) |k| {
                x.CH(ido - 1, k, 0).* = 2 * x.CC(ido - 1, 0, k);
                x.CH(ido - 1, k, 1).* = -2 * x.CC(0, 1, k);
            };
            if (ido <= 2) return;
            for (0..l1) |k| {
                var i: usize = 2;
                while (i < ido) : (i += 2) {
                    const ic = ido - i;
                    x.CH(i - 1, k, 0).* = x.CC(i - 1, 0, k) + x.CC(ic - 1, 1, k);
                    const tr2 = x.CC(i - 1, 0, k) - x.CC(ic - 1, 1, k);
                    const ti2 = x.CC(i, 0, k) + x.CC(ic, 1, k);
                    x.CH(i, k, 0).* = x.CC(i, 0, k) - x.CC(ic, 1, k);
                    twiddleOut(x, i, k, 1, ti2, tr2);
                }
            }
        }

        pub fn radb3(ido: usize, l1: usize, cc: [*]const T, ch: [*]T, wa: [*]const T) void {
            const x: ix.Backward(T) = .{ .ido = ido, .l1 = l1, .cdim = 3, .cc = cc, .ch = ch, .wa = wa };
            for (0..l1) |k| {
                const tr2 = 2 * x.CC(ido - 1, 1, k);
                const cr2 = x.CC(0, 0, k) + taur * tr2;
                x.CH(0, k, 0).* = x.CC(0, 0, k) + tr2;
                const ci3 = 2 * taui * x.CC(0, 2, k);
                ix.pm(T, x.CH(0, k, 2), x.CH(0, k, 1), cr2, ci3);
            }
            if (ido == 1) return;
            for (0..l1) |k| {
                var i: usize = 2;
                while (i < ido) : (i += 2) {
                    const ic = ido - i;
                    const tr2 = x.CC(i - 1, 2, k) + x.CC(ic - 1, 1, k);
                    const ti2 = x.CC(i, 2, k) - x.CC(ic, 1, k);
                    const cr2 = x.CC(i - 1, 0, k) + taur * tr2;
                    const ci2 = x.CC(i, 0, k) + taur * ti2;
                    x.CH(i - 1, k, 0).* = x.CC(i - 1, 0, k) + tr2;
                    x.CH(i, k, 0).* = x.CC(i, 0, k) + ti2;
                    const cr3 = taui * (x.CC(i - 1, 2, k) - x.CC(ic - 1, 1, k));
                    const ci3 = taui * (x.CC(i, 2, k) + x.CC(ic, 1, k));
                    const dr3 = cr2 + ci3;
                    const dr2 = cr2 - ci3;
                    const di2 = ci2 + cr3;
                    const di3 = ci2 - cr3;
                    twiddleOut(x, i, k, 1, di2, dr2);
                    twiddleOut(x, i, k, 2, di3, dr3);
                }
            }
        }

        pub fn radb4(ido: usize, l1: usize, cc: [*]const T, ch: [*]T, wa: [*]const T) void {
            const x: ix.Backward(T) = .{ .ido = ido, .l1 = l1, .cdim = 4, .cc = cc, .ch = ch, .wa = wa };
            for (0..l1) |k| {
                const tr2 = x.CC(0, 0, k) + x.CC(ido - 1, 3, k);
                const tr1 = x.CC(0, 0, k) - x.CC(ido - 1, 3, k);
                const tr3 = 2 * x.CC(ido - 1, 1, k);
                const tr4 = 2 * x.CC(0, 2, k);
                ix.pm(T, x.CH(0, k, 0), x.CH(0, k, 2), tr2, tr3);
                ix.pm(T, x.CH(0, k, 3), x.CH(0, k, 1), tr1, tr4);
            }
            if (ido & 1 == 0) for (0..l1) |k| {
                const ti1 = x.CC(0, 3, k) + x.CC(0, 1, k);
                const ti2 = x.CC(0, 3, k) - x.CC(0, 1, k);
                const tr2 = x.CC(ido - 1, 0, k) + x.CC(ido - 1, 2, k);
                const tr1 = x.CC(ido - 1, 0, k) - x.CC(ido - 1, 2, k);
                x.CH(ido - 1, k, 0).* = tr2 + tr2;
                x.CH(ido - 1, k, 1).* = sqrt2 * (tr1 - ti1);
                x.CH(ido - 1, k, 2).* = ti2 + ti2;
                x.CH(ido - 1, k, 3).* = -sqrt2 * (tr1 + ti1);
            };
            if (ido <= 2) return;
            for (0..l1) |k| {
                var i: usize = 2;
                while (i < ido) : (i += 2) {
                    const ic = ido - i;
                    const tr2 = x.CC(i - 1, 0, k) + x.CC(ic - 1, 3, k);
                    const tr1 = x.CC(i - 1, 0, k) - x.CC(ic - 1, 3, k);
                    const ti1 = x.CC(i, 0, k) + x.CC(ic, 3, k);
                    const ti2 = x.CC(i, 0, k) - x.CC(ic, 3, k);
                    const tr4 = x.CC(i, 2, k) + x.CC(ic, 1, k);
                    const ti3 = x.CC(i, 2, k) - x.CC(ic, 1, k);
                    const tr3 = x.CC(i - 1, 2, k) + x.CC(ic - 1, 1, k);
                    const ti4 = x.CC(i - 1, 2, k) - x.CC(ic - 1, 1, k);
                    x.CH(i - 1, k, 0).* = tr2 + tr3;
                    const cr3 = tr2 - tr3;
                    x.CH(i, k, 0).* = ti2 + ti3;
                    const ci3 = ti2 - ti3;
                    const cr4 = tr1 + tr4;
                    const cr2 = tr1 - tr4;
                    const ci2 = ti1 + ti4;
                    const ci4 = ti1 - ti4;
                    twiddleOut(x, i, k, 1, ci2, cr2);
                    twiddleOut(x, i, k, 2, ci3, cr3);
                    twiddleOut(x, i, k, 3, ci4, cr4);
                }
            }
        }

        pub fn radb5(ido: usize, l1: usize, cc: [*]const T, ch: [*]T, wa: [*]const T) void {
            const x: ix.Backward(T) = .{ .ido = ido, .l1 = l1, .cdim = 5, .cc = cc, .ch = ch, .wa = wa };
            for (0..l1) |k| {
                const ti5 = x.CC(0, 2, k) + x.CC(0, 2, k);
                const ti4 = x.CC(0, 4, k) + x.CC(0, 4, k);
                const tr2 = x.CC(ido - 1, 1, k) + x.CC(ido - 1, 1, k);
                const tr3 = x.CC(ido - 1, 3, k) + x.CC(ido - 1, 3, k);
                x.CH(0, k, 0).* = x.CC(0, 0, k) + tr2 + tr3;
                const cr2 = x.CC(0, 0, k) + tr11 * tr2 + tr12 * tr3;
                const cr3 = x.CC(0, 0, k) + tr12 * tr2 + tr11 * tr3;
                const c = ix.mulpm(T, ti5, ti4, ti11, ti12);
                ix.pm(T, x.CH(0, k, 4), x.CH(0, k, 1), cr2, c.a);
                ix.pm(T, x.CH(0, k, 3), x.CH(0, k, 2), cr3, c.b);
            }
            if (ido == 1) return;
            for (0..l1) |k| {
                var i: usize = 2;
                while (i < ido) : (i += 2) {
                    const ic = ido - i;
                    const tr2 = x.CC(i - 1, 2, k) + x.CC(ic - 1, 1, k);
                    const tr5 = x.CC(i - 1, 2, k) - x.CC(ic - 1, 1, k);
                    const ti5 = x.CC(i, 2, k) + x.CC(ic, 1, k);
                    const ti2 = x.CC(i, 2, k) - x.CC(ic, 1, k);
                    const tr3 = x.CC(i - 1, 4, k) + x.CC(ic - 1, 3, k);
                    const tr4 = x.CC(i - 1, 4, k) - x.CC(ic - 1, 3, k);
                    const ti4 = x.CC(i, 4, k) + x.CC(ic, 3, k);
                    const ti3 = x.CC(i, 4, k) - x.CC(ic, 3, k);
                    x.CH(i - 1, k, 0).* = x.CC(i - 1, 0, k) + tr2 + tr3;
                    x.CH(i, k, 0).* = x.CC(i, 0, k) + ti2 + ti3;
                    const cr2 = x.CC(i - 1, 0, k) + tr11 * tr2 + tr12 * tr3;
                    const ci2 = x.CC(i, 0, k) + tr11 * ti2 + tr12 * ti3;
                    const cr3 = x.CC(i - 1, 0, k) + tr12 * tr2 + tr11 * tr3;
                    const ci3 = x.CC(i, 0, k) + tr12 * ti2 + tr11 * ti3;
                    const r = ix.mulpm(T, tr5, tr4, ti11, ti12);
                    const s = ix.mulpm(T, ti5, ti4, ti11, ti12);
                    const cr5 = r.a;
                    const cr4 = r.b;
                    const ci5 = s.a;
                    const ci4 = s.b;
                    const dr4 = cr3 + ci4;
                    const dr3 = cr3 - ci4;
                    const di3 = ci3 + cr4;
                    const di4 = ci3 - cr4;
                    const dr5 = cr2 + ci5;
                    const dr2 = cr2 - ci5;
                    const di2 = ci2 + cr5;
                    const di5 = ci2 - cr5;
                    twiddleOut(x, i, k, 1, di2, dr2);
                    twiddleOut(x, i, k, 2, di3, dr3);
                    twiddleOut(x, i, k, 3, di4, dr4);
                    twiddleOut(x, i, k, 4, di5, dr5);
                }
            }
        }
    };
}
