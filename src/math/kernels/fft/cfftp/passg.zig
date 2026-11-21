//! `cfftp::passg` from PocketFFT (`pocketfft_hdronly.h`, cpp branch, commit
//! 4d865cba, BSD-3, LICENSE): one pass of a prime factor past 11. The
//! result lands in `cc`; `ch` is scratch. `wal` is read from `csarr` on the
//! fly rather than copied per call, so a transform allocates nothing.

const cmplx = @import("../cmplx.zig");

const C = cmplx.C;

pub fn passg(comptime T: type, comptime fwd: bool, ido: usize, ip: usize, l1: usize, cc: [*]C(T), ch: [*]C(T), wa: [*]const C(T), csarr: [*]const C(T)) void {
    const add = cmplx.add;
    const sub = cmplx.sub;
    const cdim = ip;
    const ipph = (ip + 1) / 2;
    const idl1 = ido * l1;
    const g: G(T, fwd) = .{ .ido = ido, .l1 = l1, .cdim = cdim, .idl1 = idl1, .cc = cc, .ch = ch, .csarr = csarr };

    for (0..l1) |k| for (0..ido) |i| {
        g.CH(i, k, 0).* = g.CC(i, 0, k);
    };
    {
        var j: usize = 1;
        var jc: usize = ip - 1;
        while (j < ipph) : ({
            j += 1;
            jc -= 1;
        }) for (0..l1) |k| for (0..ido) |i| {
            const a = g.CC(i, j, k);
            const b = g.CC(i, jc, k);
            g.CH(i, k, j).* = add(T, a, b);
            g.CH(i, k, jc).* = sub(T, a, b);
        };
    }
    for (0..l1) |k| for (0..ido) |i| {
        var tmp = g.CH(i, k, 0).*;
        for (1..ipph) |j| tmp = add(T, tmp, g.CH(i, k, j).*);
        g.CX(i, k, 0).* = tmp;
    };
    {
        var l: usize = 1;
        var lc: usize = ip - 1;
        while (l < ipph) : ({
            l += 1;
            lc -= 1;
        }) {
            // j = 0
            const w1 = g.wal(l);
            const w2 = g.wal(2 * l);
            for (0..idl1) |ik| {
                const h0 = g.CH2(ik, 0);
                const h1 = g.CH2(ik, 1);
                const h2 = g.CH2(ik, 2);
                const hp1 = g.CH2(ik, ip - 1);
                const hp2 = g.CH2(ik, ip - 2);
                g.CX2(ik, l).* = .{ .re = h0.re + w1.re * h1.re + w2.re * h2.re, .im = h0.im + w1.re * h1.im + w2.re * h2.im };
                g.CX2(ik, lc).* = .{ .re = -w1.im * hp1.im - w2.im * hp2.im, .im = w1.im * hp1.re + w2.im * hp2.re };
            }
            var iwal = 2 * l;
            var j: usize = 3;
            var jc: usize = ip - 3;
            while (j < ipph - 1) : ({
                j += 2;
                jc -= 2;
            }) {
                iwal += l;
                if (iwal > ip) iwal -= ip;
                const xw = g.wal(iwal);
                iwal += l;
                if (iwal > ip) iwal -= ip;
                const xw2 = g.wal(iwal);
                for (0..idl1) |ik| {
                    const a = g.CX2(ik, l);
                    const b = g.CX2(ik, lc);
                    const hj = g.CH2(ik, j);
                    const hj1 = g.CH2(ik, j + 1);
                    const hc = g.CH2(ik, jc);
                    const hc1 = g.CH2(ik, jc - 1);
                    a.re += hj.re * xw.re + hj1.re * xw2.re;
                    a.im += hj.im * xw.re + hj1.im * xw2.re;
                    b.re -= hc.im * xw.im + hc1.im * xw2.im;
                    b.im += hc.re * xw.im + hc1.re * xw2.im;
                }
            }
            while (j < ipph) : ({
                j += 1;
                jc -= 1;
            }) {
                iwal += l;
                if (iwal > ip) iwal -= ip;
                const xw = g.wal(iwal);
                for (0..idl1) |ik| {
                    const a = g.CX2(ik, l);
                    const b = g.CX2(ik, lc);
                    const hj = g.CH2(ik, j);
                    const hc = g.CH2(ik, jc);
                    a.re += hj.re * xw.re;
                    a.im += hj.im * xw.re;
                    b.re -= hc.im * xw.im;
                    b.im += hc.re * xw.im;
                }
            }
        }
    }

    // shuffling and twiddling
    var j: usize = 1;
    var jc: usize = ip - 1;
    while (j < ipph) : ({
        j += 1;
        jc -= 1;
    }) {
        if (ido == 1) {
            for (0..idl1) |ik| {
                const t1 = g.CX2(ik, j).*;
                const t2 = g.CX2(ik, jc).*;
                g.CX2(ik, j).* = add(T, t1, t2);
                g.CX2(ik, jc).* = sub(T, t1, t2);
            }
        } else for (0..l1) |k| {
            const t1 = g.CX(0, k, j).*;
            const t2 = g.CX(0, k, jc).*;
            g.CX(0, k, j).* = add(T, t1, t2);
            g.CX(0, k, jc).* = sub(T, t1, t2);
            for (1..ido) |i| {
                const x1 = add(T, g.CX(i, k, j).*, g.CX(i, k, jc).*);
                const x2 = sub(T, g.CX(i, k, j).*, g.CX(i, k, jc).*);
                g.CX(i, k, j).* = cmplx.twiddle(T, fwd, x1, wa[(j - 1) * (ido - 1) + i - 1]);
                g.CX(i, k, jc).* = cmplx.twiddle(T, fwd, x2, wa[(jc - 1) * (ido - 1) + i - 1]);
            }
        }
    }
}

/// `passg`'s accessors, named as PocketFFT's lambdas.
fn G(comptime T: type, comptime fwd: bool) type {
    return struct {
        ido: usize,
        l1: usize,
        cdim: usize,
        idl1: usize,
        cc: [*]C(T),
        ch: [*]C(T),
        csarr: [*]const C(T),

        const Self = @This();

        inline fn CH(s: Self, a: usize, b: usize, c: usize) *C(T) {
            return &s.ch[a + s.ido * (b + s.l1 * c)];
        }
        inline fn CC(s: Self, a: usize, b: usize, c: usize) C(T) {
            return s.cc[a + s.ido * (b + s.cdim * c)];
        }
        inline fn CX(s: Self, a: usize, b: usize, c: usize) *C(T) {
            return &s.cc[a + s.ido * (b + s.l1 * c)];
        }
        inline fn CX2(s: Self, a: usize, b: usize) *C(T) {
            return &s.cc[a + s.idl1 * b];
        }
        inline fn CH2(s: Self, a: usize, b: usize) C(T) {
            return s.ch[a + s.idl1 * b];
        }
        /// `wal[i]`: `csarr[i]`, conjugated forward; `wal[0] = 1`.
        inline fn wal(s: Self, i: usize) C(T) {
            if (i == 0) return .{ .re = 1, .im = 0 };
            const w = s.csarr[i];
            return .{ .re = w.re, .im = if (fwd) -w.im else w.im };
        }
    };
}
