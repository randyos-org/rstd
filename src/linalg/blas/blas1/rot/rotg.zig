//! Setting up a rotation: `rotg` and `zrotg`, the plane rotation that zeroes
//! `b` against `a`, and `rotmg`, the modified one for `rotm`.

const std = @import("std");
const types = @import("../../types.zig");
const util = @import("../../util.zig");

const XDouble = types.XDouble;
const Float = types.Float;
const isComplex = types.isComplex;

/// `?ROTG` (`interface/rotg.c`): the rotation `[c s; -s c]` that zeroes `b`
/// against `a`. `a` becomes `r`, `b` the reconstruction value `z`.
/// Computed in `XDouble`, as the reference does in `long double` on x86 (x87
/// where `XDouble` is), or in `T` where that is wider; a narrower input is
/// rounded once on the way out.
pub fn rotg(comptime T: type, a: *T, b: *T, c: *T, s: *T) void {
    if (comptime isComplex(T)) @compileError("rotg: use zrotg for " ++ @typeName(T));
    const X = if (@bitSizeOf(T) > @bitSizeOf(XDouble)) T else XDouble;
    const da: X = a.*;
    const db: X = b.*;
    const safmin: X = std.math.floatMin(T);
    const ada = @abs(da);
    const adb = @abs(db);
    if (adb == 0) {
        c.* = 1;
        s.* = 0;
        b.* = 0;
        return;
    }
    if (ada == 0) {
        c.* = 0;
        s.* = 1;
        a.* = b.*;
        b.* = 1;
        return;
    }
    const sc = @min(@max(safmin, @max(ada, adb)), 1 / safmin);
    const sigma = std.math.copysign(@as(X, 1), if (ada > adb) da else db);
    const aa = da / sc;
    const bb = db / sc;
    const r = sigma * sc * @sqrt(aa * aa + bb * bb);
    const cx = da / r;
    const sx = db / r;
    var z: X = 1;
    if (ada > adb) z = sx;
    if (ada <= adb and cx != 0) z = 1 / cx;
    c.* = @floatCast(cx);
    s.* = @floatCast(sx);
    a.* = @floatCast(r);
    b.* = @floatCast(z);
}

/// `CROTG`/`ZROTG`: the rotation with real `c` and complex `s` that zeroes
/// `b` against `a`, `[c s; -conj(s) c] [a; b] = [r; 0]`; `a` becomes `r`.
/// LAPACK's `zlartg`, not `interface/zrotg.c`, whose scaled branches are wrong.
pub fn zrotg(comptime T: type, a: *T, b: T, c: *Float(T), s: *T) void {
    comptime util.requireComplex(T, "zrotg");
    const F = Float(T);
    const safmin = std.math.floatMin(F);
    const safmax = 1 / safmin;
    const rtmin = @sqrt(safmin);
    const f = a.*;
    const g = b;
    if (g.re == 0 and g.im == 0) {
        c.* = 1;
        s.* = .init(0, 0);
        return;
    }
    if (f.re == 0 and f.im == 0) {
        c.* = 0;
        if (g.re == 0 or g.im == 0) {
            const d = @abs(if (g.re == 0) g.im else g.re);
            s.* = scale(T, g.conjugate(), 1 / d);
            a.* = .init(d, 0);
            return;
        }
        const g1 = @max(@abs(g.re), @abs(g.im));
        if (g1 > rtmin and g1 < @sqrt(safmax / 2)) {
            const d = @sqrt(abssq(T, g));
            s.* = scale(T, g.conjugate(), 1 / d);
            a.* = .init(d, 0);
        } else {
            const u = @min(safmax, @max(safmin, g1));
            const gs = scale(T, g, 1 / u);
            const d = @sqrt(abssq(T, gs));
            s.* = scale(T, gs.conjugate(), 1 / d);
            a.* = .init(d * u, 0);
        }
        return;
    }
    const f1 = @max(@abs(f.re), @abs(f.im));
    const g1 = @max(@abs(g.re), @abs(g.im));
    var rtmax = @sqrt(safmax / 4);
    // unscaled, or scaled by `u` (and `f` by `v`, as `w = v/u`)
    var u: F = 1;
    var w: F = 1;
    var fs = f;
    var gs = g;
    if (!(f1 > rtmin and f1 < rtmax and g1 > rtmin and g1 < rtmax)) {
        u = @min(safmax, @max(safmin, @max(f1, g1)));
        gs = scale(T, g, 1 / u);
        if (f1 / u < rtmin) {
            const v = @min(safmax, @max(safmin, f1));
            w = v / u;
            fs = scale(T, f, 1 / v);
        } else {
            fs = scale(T, f, 1 / u);
        }
    }
    const f2 = abssq(T, fs);
    const h2 = f2 * w * w + abssq(T, gs);
    var r: T = undefined;
    if (f2 >= h2 * safmin) {
        c.* = @sqrt(f2 / h2);
        r = scale(T, fs, 1 / c.*);
        rtmax *= 2;
        s.* = if (f2 > rtmin and h2 < rtmax)
            gs.conjugate().mul(scale(T, fs, 1 / @sqrt(f2 * h2)))
        else
            gs.conjugate().mul(scale(T, r, 1 / h2));
    } else {
        // intermediates might overflow
        const d = @sqrt(f2 * h2);
        c.* = f2 / d;
        r = if (c.* >= safmin) scale(T, fs, 1 / c.*) else scale(T, fs, h2 / d);
        s.* = gs.conjugate().mul(scale(T, fs, 1 / d));
    }
    c.* *= w;
    a.* = scale(T, r, u);
}

fn scale(comptime T: type, z: T, f: Float(T)) T {
    return .init(z.re * f, z.im * f);
}

fn abssq(comptime T: type, z: T) Float(T) {
    return z.re * z.re + z.im * z.im;
}

/// `?ROTMG` (`interface/rotmg.c`): the modified rotation that zeroes the
/// second component of `(sqrt(d1)*x1, sqrt(d2)*y1)`, as `param` for `rotm`.
/// `d1`, `d2` and `x1` are updated.
pub fn rotmg(comptime T: type, d1: *T, d2: *T, x1: *T, y1: T, param: *[5]T) void {
    comptime util.requireReal(T, "rotmg");
    // the reference's bounds, 4096 and its square; a narrow type (`f16`)
    // keeps `d1`, `d2` within its range by 16 and 256 instead
    const narrow = @bitSizeOf(T) < 32;
    const gam: T = if (narrow) 16 else 4096;
    const gamsq: T = if (narrow) 256 else 16777216;
    const rgamsq: T = if (narrow) 1.0 / 256.0 else 5.9604645e-8;
    var h11: T = 0;
    var h21: T = 0;
    var h12: T = 0;
    var h22: T = 0;
    var flag: T = -1;
    if (d2.* == 0 or y1 == 0) {
        param[0] = -2;
        return;
    }
    if (d1.* < 0) {
        d1.* = 0;
        d2.* = 0;
        x1.* = 0;
    } else if ((d1.* == 0 or x1.* == 0) and d2.* > 0) {
        flag = 1;
        h12 = 1;
        h21 = -1;
        x1.* = y1;
        std.mem.swap(T, d1, d2);
    } else {
        const p2 = d2.* * y1;
        const p1 = d1.* * x1.*;
        const q2 = p2 * y1;
        const q1 = p1 * x1.*;
        if (@abs(q1) > @abs(q2)) {
            h11 = 1;
            h22 = 1;
            h21 = -y1 / x1.*;
            h12 = p2 / p1;
            const u = 1 - h12 * h21;
            if (u <= 0) {
                // only through rounding, as in the reference
                h11 = 0;
                h12 = 0;
                h21 = 0;
                h22 = 0;
                d1.* = 0;
                d2.* = 0;
                x1.* = 0;
            } else {
                flag = 0;
                d1.* /= u;
                d2.* /= u;
                x1.* *= u;
            }
        } else if (q2 < 0) {
            d1.* = 0;
            d2.* = 0;
            x1.* = 0;
        } else {
            flag = 1;
            h21 = -1;
            h12 = 1;
            h11 = p1 / p2;
            h22 = x1.* / y1;
            const u = 1 + h11 * h22;
            const t = d2.* / u;
            d2.* = d1.* / u;
            d1.* = t;
            x1.* = y1 * u;
        }
        while (d1.* <= rgamsq and d1.* != 0) {
            flag = -1;
            d1.* *= gam * gam;
            x1.* /= gam;
            h11 /= gam;
            h12 /= gam;
        }
        while (@abs(d1.*) > gamsq) {
            flag = -1;
            d1.* /= gam * gam;
            x1.* *= gam;
            h11 *= gam;
            h12 *= gam;
        }
        while (@abs(d2.*) <= rgamsq and d2.* != 0) {
            flag = -1;
            d2.* *= gam * gam;
            h21 /= gam;
            h22 /= gam;
        }
        while (@abs(d2.*) > gamsq) {
            flag = -1;
            d2.* /= gam * gam;
            h21 *= gam;
            h22 *= gam;
        }
    }
    if (flag < 0) {
        param[1] = h11;
        param[2] = h21;
        param[3] = h12;
        param[4] = h22;
    } else if (flag == 0) {
        param[2] = h21;
        param[3] = h12;
    } else {
        param[1] = h11;
        param[4] = h22;
    }
    param[0] = flag;
}

// ---- tests -----------------------------------------------------------------

const check = @import("../../_test/check.zig");
const reference = @import("../../_test/reference.zig");
const interface = @import("../interface.zig");
const Complex32 = types.Complex32;
const Complex64 = types.Complex64;
const reals = .{ f32, f64 } ++ check.x_real;

test "rotg zeroes b, and rotmg's param zeroes the second component" {
    inline for (reals) |T| {
        var a: T = 3;
        var b: T = 4;
        var c: T = undefined;
        var s: T = undefined;
        rotg(T, &a, &b, &c, &s);
        try check.expectClose(T, 5, a, check.tolerance(T, 4));
        try check.expectClose(T, 0.6, c, check.tolerance(T, 4));
        try check.expectClose(T, 0.8, s, check.tolerance(T, 4));
        // |a| <= |b| and c != 0: z = 1/c
        try check.expectClose(T, 1 / @as(T, 0.6), b, check.tolerance(T, 8));

        var d1: T = 2;
        var d2: T = 3;
        var x1: T = 1.5;
        const y1: T = 0.5;
        const d1_0 = d1;
        const d2_0 = d2;
        const x1_0 = x1;
        var param: [5]T = undefined;
        rotmg(T, &d1, &d2, &x1, y1, &param);
        // apply H to (x1, y1): the second component vanishes
        var hx = [1]T{x1_0};
        var hy = [1]T{y1};
        try interface.rotm(check.io, T, 1, &hx, 1, &hy, 1, &param);
        try check.expectClose(T, 0, hy[0], check.tolerance(T, 8));
        try check.expectClose(T, x1, hx[0], check.tolerance(T, 8));
        // and the weighted length is kept: d1*x1^2 + d2*y1^2
        try check.expectClose(T, d1_0 * x1_0 * x1_0 + d2_0 * y1 * y1, d1 * x1 * x1, check.tolerance(T, 16));
    }
}

test "zrotg zeroes b against a over every branch" {
    inline for (.{ Complex32, Complex64 } ++ check.x_complex) |C| {
        const F = Float(C);
        const eps: f64 = if (F == f32) 1e-5 else 1e-12;
        // magnitudes that take the unscaled and both scaled paths, with a or b zero,
        // b real or imaginary only, and a tiny against a huge b (and the reverse)
        const big: f64 = if (F == f32) 1e30 else 1e300;
        const tiny: f64 = if (F == f32) 1e-30 else 1e-300;
        const cases = [_][4]f64{
            .{ 1, 2, 3, -1 },
            .{ 1 * tiny, 2 * tiny, 3 * tiny, -1 * tiny },
            .{ 1 * big, 2 * big, 3 * big, -1 * big },
            .{ 0, 0, 3, -4 },
            .{ 0, 0, 3 * big, -4 * big },
            .{ 0, 0, 3 * tiny, -4 * tiny },
            .{ 0, 0, 0, 5 },
            .{ 0, 0, -5, 0 },
            .{ 1, 2, 0, 0 },
            .{ tiny, tiny, big, big },
            .{ big, big, tiny, tiny },
            .{ 1, -2, 0, 3 },
            .{ -2, 1, 7, 0 },
        };
        for (cases) |k| {
            const a0: C = .init(@floatCast(k[0]), @floatCast(k[1]));
            const b: C = .init(@floatCast(k[2]), @floatCast(k[3]));
            var a = a0;
            var c: F = undefined;
            var s: C = undefined;
            zrotg(C, &a, b, &c, &s);
            // [c s; -conj(s) c] [a0; b] = [r; 0], |r| = |(a0, b)|, c real, c^2 + |s|^2 = 1
            const norm = @max(@max(@abs(k[0]), @abs(k[1])), @max(@abs(k[2]), @abs(k[3])));
            const top = scale(C, a0, c).add(s.mul(b));
            const bottom = C.init(-s.re, s.im).mul(a0).add(scale(C, b, c));
            try check.expectClose(C, a, top, @floatCast(eps * norm));
            try check.expectClose(C, .init(0, 0), bottom, @floatCast(eps * norm));
            const W = reference.Wide(F);
            try check.expectClose(W, 1, @as(W, c) * c + abssq(C, s), eps);
            // |r| against the norm of (a0, b), formed in f64 on scaled parts
            const inv = 1 / norm;
            const want = @sqrt(k[0] * inv * k[0] * inv + k[1] * inv * k[1] * inv + k[2] * inv * k[2] * inv + k[3] * inv * k[3] * inv);
            const got = @sqrt(@as(W, a.re) * inv * a.re * inv + @as(W, a.im) * inv * a.im * inv);
            try check.expectClose(W, want, got, eps);
        }
    }
}

test "rotmg in f16 rescales within its range" {
    const H = f16;
    // d1 small and large enough to take each rescaling loop
    for ([_][4]H{ .{ 1.0 / 1024.0, 1, 1, 0.5 }, .{ 1000, 1, 1, 0.5 }, .{ 2, 3, 1.5, 0.5 }, .{ 1, 1.0 / 2048.0, 0.25, 4 } }) |c| {
        var d1 = c[0];
        var d2 = c[1];
        var x1 = c[2];
        const y1 = c[3];
        var param: [5]H = @splat(0);
        rotmg(H, &d1, &d2, &x1, y1, &param);
        for ([_]H{ d1, d2, x1 }) |v| try std.testing.expect(std.math.isFinite(v) and v != 0);
        var hx = [1]H{c[2]};
        var hy = [1]H{y1};
        try interface.rotm(check.io, H, 1, &hx, 1, &hy, 1, &param);
        try check.expectClose(f64, 0, hy[0], 1e-2 * @abs(@as(f64, c[2])));
        try check.expectClose(f64, x1, hx[0], 4e-3 * @abs(@as(f64, x1)));
        // the weighted length d1*x1^2 + d2*y1^2 is kept
        const was = @as(f64, c[0]) * c[2] * c[2] + @as(f64, c[1]) * y1 * y1;
        try check.expectClose(f64, was, @as(f64, d1) * x1 * x1, 1e-2 * was);
    }
}

test "rotg handles scaled and special inputs" {
    inline for (reals) |T| {
        const big: T = std.math.floatMax(T) / 4;
        const tiny: T = std.math.floatMin(T) * 4;
        for ([_][2]T{ .{ big, big }, .{ tiny, tiny }, .{ big, tiny }, .{ tiny, big }, .{ -3, 4 }, .{ 3, -4 } }) |ab| {
            var a = ab[0];
            var b = ab[1];
            var c: T = undefined;
            var s: T = undefined;
            rotg(T, &a, &b, &c, &s);
            try std.testing.expect(std.math.isFinite(a));
            const W = reference.Wide(T);
            try check.expectClose(W, 1, @as(W, c) * c + @as(W, s) * s, 1e-6);
            // the rotation zeroes b: -s*a0 + c*b0 = 0, relative to the norm
            const norm = @max(@abs(ab[0]), @abs(ab[1]));
            const resid = -@as(W, s) * ab[0] + @as(W, c) * ab[1];
            try std.testing.expect(@abs(resid) <= 1e-6 * norm);
        }
        var a: T = 0;
        var b: T = 0;
        var c: T = undefined;
        var s: T = undefined;
        rotg(T, &a, &b, &c, &s);
        try std.testing.expectEqual(@as(T, 1), c);
        try std.testing.expectEqual(@as(T, 0), s);
        a = std.math.nan(T);
        b = 1;
        rotg(T, &a, &b, &c, &s);
        try std.testing.expect(std.math.isNan(a));
    }
}
