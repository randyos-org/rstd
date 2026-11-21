//! `cfftp`'s fixed-radix passes from PocketFFT (`pocketfft_hdronly.h`, cpp
//! branch, commit 4d865cba, BSD-3, LICENSE): `pass2`, `pass4`, `pass8`, and
//! `pass3`, `pass5`, `pass7`, `pass11` as one comptime-generic odd pass
//! (their macros' sums in the same order). `cc` in, `ch` out, `wa` the
//! twiddles; column `i = 0` takes none.

const std = @import("std");
const cmplx = @import("../cmplx.zig");

const C = cmplx.C;

/// `CH(a, b, c)`, `CC(a, b, c)` and `WA(x, i)`'s index arithmetic.
pub const Ix = struct {
    ido: usize,
    l1: usize,
    cdim: usize,

    pub inline fn ch(s: Ix, a: usize, b: usize, c: usize) usize {
        return a + s.ido * (b + s.l1 * c);
    }
    pub inline fn cc(s: Ix, a: usize, b: usize, c: usize) usize {
        return a + s.ido * (b + s.cdim * c);
    }
    pub inline fn wa(s: Ix, x: usize, i: usize) usize {
        return i - 1 + x * (s.ido - 1);
    }
};

pub fn Passes(comptime T: type) type {
    return struct {
        const Cx = C(T);
        const add = cmplx.add;
        const sub = cmplx.sub;

        pub fn pass2(comptime fwd: bool, ido: usize, l1: usize, cc: [*]const Cx, ch: [*]Cx, wa: [*]const Cx) void {
            const x: Ix = .{ .ido = ido, .l1 = l1, .cdim = 2 };
            for (0..l1) |k| for (0..ido) |i| {
                const a = cc[x.cc(i, 0, k)];
                const b = cc[x.cc(i, 1, k)];
                ch[x.ch(i, k, 0)] = add(T, a, b);
                ch[x.ch(i, k, 1)] = if (i == 0) sub(T, a, b) else cmplx.twiddle(T, fwd, sub(T, a, b), wa[x.wa(0, i)]);
            };
        }

        pub fn pass4(comptime fwd: bool, ido: usize, l1: usize, cc: [*]const Cx, ch: [*]Cx, wa: [*]const Cx) void {
            const x: Ix = .{ .ido = ido, .l1 = l1, .cdim = 4 };
            for (0..l1) |k| for (0..ido) |i| {
                const c0 = cc[x.cc(i, 0, k)];
                const c1 = cc[x.cc(i, 1, k)];
                const c2 = cc[x.cc(i, 2, k)];
                const c3 = cc[x.cc(i, 3, k)];
                const t2 = add(T, c0, c2);
                const t1 = sub(T, c0, c2);
                const t3 = add(T, c1, c3);
                const t4 = cmplx.rot90(T, fwd, sub(T, c1, c3));
                if (i == 0) {
                    ch[x.ch(0, k, 0)] = add(T, t2, t3);
                    ch[x.ch(0, k, 2)] = sub(T, t2, t3);
                    ch[x.ch(0, k, 1)] = add(T, t1, t4);
                    ch[x.ch(0, k, 3)] = sub(T, t1, t4);
                } else {
                    ch[x.ch(i, k, 0)] = add(T, t2, t3);
                    ch[x.ch(i, k, 1)] = cmplx.twiddle(T, fwd, add(T, t1, t4), wa[x.wa(0, i)]);
                    ch[x.ch(i, k, 2)] = cmplx.twiddle(T, fwd, sub(T, t2, t3), wa[x.wa(1, i)]);
                    ch[x.ch(i, k, 3)] = cmplx.twiddle(T, fwd, sub(T, t1, t4), wa[x.wa(2, i)]);
                }
            };
        }

        pub fn pass8(comptime fwd: bool, ido: usize, l1: usize, cc: [*]const Cx, ch: [*]Cx, wa: [*]const Cx) void {
            const x: Ix = .{ .ido = ido, .l1 = l1, .cdim = 8 };
            for (0..l1) |k| for (0..ido) |i| {
                var a1 = add(T, cc[x.cc(i, 1, k)], cc[x.cc(i, 5, k)]);
                var a5 = sub(T, cc[x.cc(i, 1, k)], cc[x.cc(i, 5, k)]);
                var a3 = add(T, cc[x.cc(i, 3, k)], cc[x.cc(i, 7, k)]);
                var a7 = cmplx.rot90(T, fwd, sub(T, cc[x.cc(i, 3, k)], cc[x.cc(i, 7, k)]));
                const t13 = a1;
                a1 = add(T, a1, a3);
                a3 = cmplx.rot90(T, fwd, sub(T, t13, a3));
                const t57 = a5;
                a5 = cmplx.rot45(T, fwd, add(T, a5, a7));
                a7 = cmplx.rot135(T, fwd, sub(T, t57, a7));
                var a0 = add(T, cc[x.cc(i, 0, k)], cc[x.cc(i, 4, k)]);
                var a4 = sub(T, cc[x.cc(i, 0, k)], cc[x.cc(i, 4, k)]);
                var a2 = add(T, cc[x.cc(i, 2, k)], cc[x.cc(i, 6, k)]);
                var a6 = sub(T, cc[x.cc(i, 2, k)], cc[x.cc(i, 6, k)]);
                if (i == 0) {
                    const p02 = add(T, a0, a2);
                    const m02 = sub(T, a0, a2);
                    ch[x.ch(0, k, 0)] = add(T, p02, a1);
                    ch[x.ch(0, k, 4)] = sub(T, p02, a1);
                    ch[x.ch(0, k, 2)] = add(T, m02, a3);
                    ch[x.ch(0, k, 6)] = sub(T, m02, a3);
                    a6 = cmplx.rot90(T, fwd, a6);
                    const p46 = add(T, a4, a6);
                    const m46 = sub(T, a4, a6);
                    ch[x.ch(0, k, 1)] = add(T, p46, a5);
                    ch[x.ch(0, k, 5)] = sub(T, p46, a5);
                    ch[x.ch(0, k, 3)] = add(T, m46, a7);
                    ch[x.ch(0, k, 7)] = sub(T, m46, a7);
                } else {
                    const t02 = a0;
                    a0 = add(T, a0, a2);
                    a2 = sub(T, t02, a2);
                    ch[x.ch(i, k, 0)] = add(T, a0, a1);
                    ch[x.ch(i, k, 4)] = cmplx.twiddle(T, fwd, sub(T, a0, a1), wa[x.wa(3, i)]);
                    ch[x.ch(i, k, 2)] = cmplx.twiddle(T, fwd, add(T, a2, a3), wa[x.wa(1, i)]);
                    ch[x.ch(i, k, 6)] = cmplx.twiddle(T, fwd, sub(T, a2, a3), wa[x.wa(5, i)]);
                    a6 = cmplx.rot90(T, fwd, a6);
                    const t46 = a4;
                    a4 = add(T, a4, a6);
                    a6 = sub(T, t46, a6);
                    ch[x.ch(i, k, 1)] = cmplx.twiddle(T, fwd, add(T, a4, a5), wa[x.wa(0, i)]);
                    ch[x.ch(i, k, 5)] = cmplx.twiddle(T, fwd, sub(T, a4, a5), wa[x.wa(4, i)]);
                    ch[x.ch(i, k, 3)] = cmplx.twiddle(T, fwd, add(T, a6, a7), wa[x.wa(2, i)]);
                    ch[x.ch(i, k, 7)] = cmplx.twiddle(T, fwd, sub(T, a6, a7), wa[x.wa(6, i)]);
                }
            };
        }

        /// `pass3`, `pass5`, `pass7` and `pass11`: with `s_m = cc_m + cc_(p-m)`
        /// and `d_m = cc_m - cc_(p-m)`, output pair `u, p - u` is `ca ± cb`,
        /// `ca = cc_0 + Σ cos(2π um/p) s_m`, `cb = i Σ ±sin(2π um/p) d_m`.
        pub fn passOdd(comptime p: usize, comptime fwd: bool, ido: usize, l1: usize, cc: [*]const Cx, ch: [*]Cx, wa: [*]const Cx) void {
            const h = (p - 1) / 2;
            const co = comptime coefficients(p, fwd);
            const x: Ix = .{ .ido = ido, .l1 = l1, .cdim = p };
            for (0..l1) |k| for (0..ido) |i| {
                const t0 = cc[x.cc(i, 0, k)];
                var s: [h]Cx = undefined;
                var d: [h]Cx = undefined;
                inline for (1..h + 1) |m| {
                    s[m - 1] = add(T, cc[x.cc(i, m, k)], cc[x.cc(i, p - m, k)]);
                    d[m - 1] = sub(T, cc[x.cc(i, m, k)], cc[x.cc(i, p - m, k)]);
                }
                var c0 = t0;
                inline for (0..h) |m| c0 = add(T, c0, s[m]);
                ch[x.ch(i, k, 0)] = c0;
                inline for (1..h + 1) |u| {
                    var ca = t0;
                    var cbi: T = 0;
                    var cbr: T = 0;
                    inline for (0..h) |m| {
                        const xr: T = co.x[u - 1][m];
                        const yi: T = co.y[u - 1][m];
                        ca.re = ca.re + xr * s[m].re;
                        ca.im = ca.im + xr * s[m].im;
                        cbi = if (m == 0) yi * d[m].re else cbi + yi * d[m].re;
                        cbr = if (m == 0) yi * d[m].im else cbr + yi * d[m].im;
                    }
                    const cb: Cx = .{ .re = -cbr, .im = cbi };
                    if (i == 0) {
                        ch[x.ch(0, k, u)] = add(T, ca, cb);
                        ch[x.ch(0, k, p - u)] = sub(T, ca, cb);
                    } else {
                        ch[x.ch(i, k, u)] = cmplx.twiddle(T, fwd, add(T, ca, cb), wa[x.wa(u - 1, i)]);
                        ch[x.ch(i, k, p - u)] = cmplx.twiddle(T, fwd, sub(T, ca, cb), wa[x.wa(p - u - 1, i)]);
                    }
                }
            };
        }
    };
}

/// PocketFFT's `twNr`, `twNi` for `k = 1 .. (p - 1) / 2`: `cos(2πk/p)`,
/// `sin(2πk/p)`.
fn table(comptime p: usize) struct { []const comptime_float, []const comptime_float } {
    return switch (p) {
        3 => .{ &.{-0.5}, &.{0.8660254037844386467637231707529362} },
        5 => .{
            &.{ 0.3090169943749474241022934171828191, -0.8090169943749474241022934171828191 },
            &.{ 0.9510565162951535721164393333793821, 0.5877852522924731291687059546390728 },
        },
        7 => .{
            &.{ 0.6234898018587335305250048840042398, -0.2225209339563144042889025644967948, -0.9009688679024191262361023195074451 },
            &.{ 0.7818314824680298087084445266740578, 0.9749279121818236070181316829939312, 0.433883739117558120475768332848359 },
        },
        11 => .{
            &.{ 0.8412535328311811688618116489193677, 0.4154150130018864255292741492296232, -0.1423148382732851404437926686163697, -0.6548607339452850640569250724662936, -0.9594929736144973898903680570663277 },
            &.{ 0.5406408174555975821076359543186917, 0.9096319953545183714117153830790285, 0.9898214418809327323760920377767188, 0.7557495743542582837740358439723444, 0.2817325568414296977114179153466169 },
        },
        else => @compileError("no fixed pass for this radix"),
    };
}

/// Pair `u`'s factors on `s_m` and `d_m`: the angle `um mod p`, folded into
/// the table's half with the sine's sign, as the macros' arguments list them.
fn coefficients(comptime p: usize, comptime fwd: bool) struct { x: [(p - 1) / 2][(p - 1) / 2]comptime_float, y: [(p - 1) / 2][(p - 1) / 2]comptime_float } {
    const h = (p - 1) / 2;
    const tw = table(p);
    const sign: comptime_float = if (fwd) -1 else 1;
    var x: [h][h]comptime_float = undefined;
    var y: [h][h]comptime_float = undefined;
    for (1..h + 1) |u| for (1..h + 1) |m| {
        const k = (u * m) % p;
        if (k <= h) {
            x[u - 1][m - 1] = tw[0][k - 1];
            y[u - 1][m - 1] = sign * tw[1][k - 1];
        } else {
            x[u - 1][m - 1] = tw[0][p - k - 1];
            y[u - 1][m - 1] = -sign * tw[1][p - k - 1];
        }
    };
    return .{ .x = x, .y = y };
}
