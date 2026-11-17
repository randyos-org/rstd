//! `rfftp`'s passes for a prime factor past 5 from PocketFFT
//! (`pocketfft_hdronly.h`, cpp branch, commit 4d865cba, BSD-3, LICENSE):
//! `radfg` (its result in `cc`) and `radbg` (in `ch`). `csarr` holds the
//! factor's `cos`, `sin` pairs.

/// The accessors over both buffers, named as PocketFFT's lambdas.
fn G(comptime T: type) type {
    return struct {
        ido: usize,
        l1: usize,
        cdim: usize,
        idl1: usize,
        cc: [*]T,
        ch: [*]T,

        const Self = @This();

        /// `cc` by the radix.
        inline fn CC(s: Self, a: usize, b: usize, c: usize) *T {
            return &s.cc[a + s.ido * (b + s.cdim * c)];
        }
        /// `cc` by `l1`.
        inline fn C1(s: Self, a: usize, b: usize, c: usize) *T {
            return &s.cc[a + s.ido * (b + s.l1 * c)];
        }
        inline fn C2(s: Self, a: usize, b: usize) *T {
            return &s.cc[a + s.idl1 * b];
        }
        inline fn CH(s: Self, a: usize, b: usize, c: usize) *T {
            return &s.ch[a + s.ido * (b + s.l1 * c)];
        }
        inline fn CH2(s: Self, a: usize, b: usize) *T {
            return &s.ch[a + s.idl1 * b];
        }
    };
}

pub fn radfg(comptime T: type, ido: usize, ip: usize, l1: usize, cc: [*]T, ch: [*]T, wa: [*]const T, csarr: [*]const T) void {
    const ipph = (ip + 1) / 2;
    const idl1 = ido * l1;
    const g: G(T) = .{ .ido = ido, .l1 = l1, .cdim = ip, .idl1 = idl1, .cc = cc, .ch = ch };

    if (ido > 1) {
        var j: usize = 1;
        var jc: usize = ip - 1;
        while (j < ipph) : ({
            j += 1;
            jc -= 1;
        }) {
            const is = (j - 1) * (ido - 1);
            const is2 = (jc - 1) * (ido - 1);
            for (0..l1) |k| {
                var idij = is;
                var idij2 = is2;
                var i: usize = 1;
                while (i <= ido - 2) : (i += 2) {
                    const t1 = g.C1(i, k, j).*;
                    const t2 = g.C1(i + 1, k, j).*;
                    const t3 = g.C1(i, k, jc).*;
                    const t4 = g.C1(i + 1, k, jc).*;
                    const x1 = wa[idij] * t1 + wa[idij + 1] * t2;
                    const x2 = wa[idij] * t2 - wa[idij + 1] * t1;
                    const x3 = wa[idij2] * t3 + wa[idij2 + 1] * t4;
                    const x4 = wa[idij2] * t4 - wa[idij2 + 1] * t3;
                    g.C1(i, k, j).* = x3 + x1;
                    g.C1(i + 1, k, jc).* = x3 - x1;
                    g.C1(i + 1, k, j).* = x2 + x4;
                    g.C1(i, k, jc).* = x2 - x4;
                    idij += 2;
                    idij2 += 2;
                }
            }
        }
    }
    {
        var j: usize = 1;
        var jc: usize = ip - 1;
        while (j < ipph) : ({
            j += 1;
            jc -= 1;
        }) for (0..l1) |k| {
            // MPINPLACE(C1(0, k, jc), C1(0, k, j))
            const a = g.C1(0, k, jc);
            const b = g.C1(0, k, j);
            const t = a.*;
            a.* -= b.*;
            b.* = t + b.*;
        };
    }
    {
        var l: usize = 1;
        var lc: usize = ip - 1;
        while (l < ipph) : ({
            l += 1;
            lc -= 1;
        }) {
            for (0..idl1) |ik| {
                g.CH2(ik, l).* = g.C2(ik, 0).* + csarr[2 * l] * g.C2(ik, 1).* + csarr[4 * l] * g.C2(ik, 2).*;
                g.CH2(ik, lc).* = csarr[2 * l + 1] * g.C2(ik, ip - 1).* + csarr[4 * l + 1] * g.C2(ik, ip - 2).*;
            }
            var iang = 2 * l;
            var j: usize = 3;
            var jc: usize = ip - 3;
            while (j + 3 < ipph) : ({
                j += 4;
                jc -= 4;
            }) {
                var ar: [4]T = undefined;
                var ai: [4]T = undefined;
                for (0..4) |q| {
                    iang += l;
                    if (iang >= ip) iang -= ip;
                    ar[q] = csarr[2 * iang];
                    ai[q] = csarr[2 * iang + 1];
                }
                for (0..idl1) |ik| {
                    g.CH2(ik, l).* += ar[0] * g.C2(ik, j).* + ar[1] * g.C2(ik, j + 1).* + ar[2] * g.C2(ik, j + 2).* + ar[3] * g.C2(ik, j + 3).*;
                    g.CH2(ik, lc).* += ai[0] * g.C2(ik, jc).* + ai[1] * g.C2(ik, jc - 1).* + ai[2] * g.C2(ik, jc - 2).* + ai[3] * g.C2(ik, jc - 3).*;
                }
            }
            while (j + 1 < ipph) : ({
                j += 2;
                jc -= 2;
            }) {
                iang += l;
                if (iang >= ip) iang -= ip;
                const ar1 = csarr[2 * iang];
                const ai1 = csarr[2 * iang + 1];
                iang += l;
                if (iang >= ip) iang -= ip;
                const ar2 = csarr[2 * iang];
                const ai2 = csarr[2 * iang + 1];
                for (0..idl1) |ik| {
                    g.CH2(ik, l).* += ar1 * g.C2(ik, j).* + ar2 * g.C2(ik, j + 1).*;
                    g.CH2(ik, lc).* += ai1 * g.C2(ik, jc).* + ai2 * g.C2(ik, jc - 1).*;
                }
            }
            while (j < ipph) : ({
                j += 1;
                jc -= 1;
            }) {
                iang += l;
                if (iang >= ip) iang -= ip;
                const ar = csarr[2 * iang];
                const ai = csarr[2 * iang + 1];
                for (0..idl1) |ik| {
                    g.CH2(ik, l).* += ar * g.C2(ik, j).*;
                    g.CH2(ik, lc).* += ai * g.C2(ik, jc).*;
                }
            }
        }
    }
    for (0..idl1) |ik| g.CH2(ik, 0).* = g.C2(ik, 0).*;
    for (1..ipph) |j| for (0..idl1) |ik| {
        g.CH2(ik, 0).* += g.C2(ik, j).*;
    };
    for (0..l1) |k| for (0..ido) |i| {
        g.CC(i, 0, k).* = g.CH(i, k, 0).*;
    };
    {
        var j: usize = 1;
        var jc: usize = ip - 1;
        while (j < ipph) : ({
            j += 1;
            jc -= 1;
        }) {
            const j2 = 2 * j - 1;
            for (0..l1) |k| {
                g.CC(ido - 1, j2, k).* = g.CH(0, k, j).*;
                g.CC(0, j2 + 1, k).* = g.CH(0, k, jc).*;
            }
        }
    }
    if (ido == 1) return;
    var j: usize = 1;
    var jc: usize = ip - 1;
    while (j < ipph) : ({
        j += 1;
        jc -= 1;
    }) {
        const j2 = 2 * j - 1;
        for (0..l1) |k| {
            var i: usize = 1;
            var ic: usize = ido -% 3;
            while (i <= ido - 2) : ({
                i += 2;
                ic -%= 2;
            }) {
                g.CC(i, j2 + 1, k).* = g.CH(i, k, j).* + g.CH(i, k, jc).*;
                g.CC(ic, j2, k).* = g.CH(i, k, j).* - g.CH(i, k, jc).*;
                g.CC(i + 1, j2 + 1, k).* = g.CH(i + 1, k, j).* + g.CH(i + 1, k, jc).*;
                g.CC(ic + 1, j2, k).* = g.CH(i + 1, k, jc).* - g.CH(i + 1, k, j).*;
            }
        }
    }
}

pub fn radbg(comptime T: type, ido: usize, ip: usize, l1: usize, cc: [*]T, ch: [*]T, wa: [*]const T, csarr: [*]const T) void {
    const ipph = (ip + 1) / 2;
    const idl1 = ido * l1;
    const g: G(T) = .{ .ido = ido, .l1 = l1, .cdim = ip, .idl1 = idl1, .cc = cc, .ch = ch };

    for (0..l1) |k| for (0..ido) |i| {
        g.CH(i, k, 0).* = g.CC(i, 0, k).*;
    };
    {
        var j: usize = 1;
        var jc: usize = ip - 1;
        while (j < ipph) : ({
            j += 1;
            jc -= 1;
        }) {
            const j2 = 2 * j - 1;
            for (0..l1) |k| {
                g.CH(0, k, j).* = 2 * g.CC(ido - 1, j2, k).*;
                g.CH(0, k, jc).* = 2 * g.CC(0, j2 + 1, k).*;
            }
        }
    }
    if (ido != 1) {
        var j: usize = 1;
        var jc: usize = ip - 1;
        while (j < ipph) : ({
            j += 1;
            jc -= 1;
        }) {
            const j2 = 2 * j - 1;
            for (0..l1) |k| {
                var i: usize = 1;
                var ic: usize = ido -% 3;
                while (i <= ido - 2) : ({
                    i += 2;
                    ic -%= 2;
                }) {
                    g.CH(i, k, j).* = g.CC(i, j2 + 1, k).* + g.CC(ic, j2, k).*;
                    g.CH(i, k, jc).* = g.CC(i, j2 + 1, k).* - g.CC(ic, j2, k).*;
                    g.CH(i + 1, k, j).* = g.CC(i + 1, j2 + 1, k).* - g.CC(ic + 1, j2, k).*;
                    g.CH(i + 1, k, jc).* = g.CC(i + 1, j2 + 1, k).* + g.CC(ic + 1, j2, k).*;
                }
            }
        }
    }
    {
        var l: usize = 1;
        var lc: usize = ip - 1;
        while (l < ipph) : ({
            l += 1;
            lc -= 1;
        }) {
            for (0..idl1) |ik| {
                g.C2(ik, l).* = g.CH2(ik, 0).* + csarr[2 * l] * g.CH2(ik, 1).* + csarr[4 * l] * g.CH2(ik, 2).*;
                g.C2(ik, lc).* = csarr[2 * l + 1] * g.CH2(ik, ip - 1).* + csarr[4 * l + 1] * g.CH2(ik, ip - 2).*;
            }
            var iang = 2 * l;
            var j: usize = 3;
            var jc: usize = ip - 3;
            while (j + 3 < ipph) : ({
                j += 4;
                jc -= 4;
            }) {
                var ar: [4]T = undefined;
                var ai: [4]T = undefined;
                for (0..4) |q| {
                    iang += l;
                    if (iang > ip) iang -= ip;
                    ar[q] = csarr[2 * iang];
                    ai[q] = csarr[2 * iang + 1];
                }
                for (0..idl1) |ik| {
                    g.C2(ik, l).* += ar[0] * g.CH2(ik, j).* + ar[1] * g.CH2(ik, j + 1).* + ar[2] * g.CH2(ik, j + 2).* + ar[3] * g.CH2(ik, j + 3).*;
                    g.C2(ik, lc).* += ai[0] * g.CH2(ik, jc).* + ai[1] * g.CH2(ik, jc - 1).* + ai[2] * g.CH2(ik, jc - 2).* + ai[3] * g.CH2(ik, jc - 3).*;
                }
            }
            while (j + 1 < ipph) : ({
                j += 2;
                jc -= 2;
            }) {
                iang += l;
                if (iang > ip) iang -= ip;
                const ar1 = csarr[2 * iang];
                const ai1 = csarr[2 * iang + 1];
                iang += l;
                if (iang > ip) iang -= ip;
                const ar2 = csarr[2 * iang];
                const ai2 = csarr[2 * iang + 1];
                for (0..idl1) |ik| {
                    g.C2(ik, l).* += ar1 * g.CH2(ik, j).* + ar2 * g.CH2(ik, j + 1).*;
                    g.C2(ik, lc).* += ai1 * g.CH2(ik, jc).* + ai2 * g.CH2(ik, jc - 1).*;
                }
            }
            while (j < ipph) : ({
                j += 1;
                jc -= 1;
            }) {
                iang += l;
                if (iang > ip) iang -= ip;
                const war = csarr[2 * iang];
                const wai = csarr[2 * iang + 1];
                for (0..idl1) |ik| {
                    g.C2(ik, l).* += war * g.CH2(ik, j).*;
                    g.C2(ik, lc).* += wai * g.CH2(ik, jc).*;
                }
            }
        }
    }
    for (1..ipph) |j| for (0..idl1) |ik| {
        g.CH2(ik, 0).* += g.CH2(ik, j).*;
    };
    {
        var j: usize = 1;
        var jc: usize = ip - 1;
        while (j < ipph) : ({
            j += 1;
            jc -= 1;
        }) for (0..l1) |k| {
            const a = g.C1(0, k, j).*;
            const b = g.C1(0, k, jc).*;
            g.CH(0, k, jc).* = a + b;
            g.CH(0, k, j).* = a - b;
        };
    }
    if (ido == 1) return;
    {
        var j: usize = 1;
        var jc: usize = ip - 1;
        while (j < ipph) : ({
            j += 1;
            jc -= 1;
        }) for (0..l1) |k| {
            var i: usize = 1;
            while (i <= ido - 2) : (i += 2) {
                g.CH(i, k, j).* = g.C1(i, k, j).* - g.C1(i + 1, k, jc).*;
                g.CH(i, k, jc).* = g.C1(i, k, j).* + g.C1(i + 1, k, jc).*;
                g.CH(i + 1, k, j).* = g.C1(i + 1, k, j).* + g.C1(i, k, jc).*;
                g.CH(i + 1, k, jc).* = g.C1(i + 1, k, j).* - g.C1(i, k, jc).*;
            }
        };
    }
    // all in CH
    for (1..ip) |j| {
        const is = (j - 1) * (ido - 1);
        for (0..l1) |k| {
            var idij = is;
            var i: usize = 1;
            while (i <= ido - 2) : (i += 2) {
                const t1 = g.CH(i, k, j).*;
                const t2 = g.CH(i + 1, k, j).*;
                g.CH(i, k, j).* = wa[idij] * t1 - wa[idij + 1] * t2;
                g.CH(i + 1, k, j).* = wa[idij] * t2 + wa[idij + 1] * t1;
                idij += 2;
            }
        }
    }
}
