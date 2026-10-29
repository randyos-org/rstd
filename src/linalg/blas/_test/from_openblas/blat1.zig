//! The level-1 test programs: `test/sblat1.f`, `dblat1.f`, `cblat1.f`,
//! `zblat1.f`, and the parts of their CBLAS versions `ctest/c_[sdcz]blat1.f`
//! that differ (the real `CHECK3`, `rot` and `rotm` against the programs' own
//! Fortran, and the complex `CHECK1` data), OpenBLAS v0.3.34-317. The
//! `CHECK0` to `CHECK3` logic, `STEST`, `CTEST` and `ITEST1` are ported; the
//! tables at the end are the programs' `DATA` (single and double precision
//! share them; `st7b` and `ssize3` are `sblat1.f`'s alone).

const std = @import("std");
const ob = @import("__root__.zig");
const blas = @import("../../__root__.zig");

const Float = ob.Float;
const isComplex = ob.isComplex;
const C32 = std.math.Complex(f32);
const C64 = std.math.Complex(f64);

/// `SFAC`: the factor the differences are scaled by before comparing.
const sfac = 9.765625e-4;

/// The case a failure is reported against, as the Fortran's `COMBLA` does.
const Case = struct { name: []const u8, n: isize = 0, incx: isize = 0, incy: isize = 0 };

/// `STEST`: `comp` against `want` within `|size| * eps / fac` each, in `F`.
/// `old` takes ctest's criterion, that `|size| + |fac * diff|` rounds to
/// `|size|`.
fn stest(comptime F: type, c: Case, comp: []const F, want: []const f64, size: []const f64, fac: f64, old: bool) !void {
    for (comp, 0..) |got, i| {
        const sd: F = got - @as(F, @floatCast(want[i]));
        const s: F = @abs(@as(F, @floatCast(size[i])));
        const d: F = @abs(@as(F, @floatCast(fac)) * sd);
        const ok = if (old) sdiff(F, s + d, s) == 0 else d <= s * std.math.floatEps(F);
        if (ok) continue;
        ob.say("{s} n={d} incx={d} incy={d} i={d}: got {e}, true {e}, size {e}\n", .{ c.name, c.n, c.incx, c.incy, i + 1, got, want[i], size[i] });
        return error.TestUnexpectedResult;
    }
}

/// `SDIFF`: `a - b` with `a` rounded to `F` first, through memory, as x87
/// would otherwise keep it wider.
fn sdiff(comptime F: type, a: F, b: F) F {
    var t: F = undefined;
    const v: *volatile F = &t;
    v.* = a;
    return v.* - b;
}

fn stest1(comptime F: type, c: Case, comp: F, want: f64, size: f64, fac: f64) !void {
    try stest(F, c, &.{comp}, &.{want}, &.{size}, fac, false);
}

/// `CTEST`: `STEST` over the real and imaginary parts.
fn ctest(comptime C: type, c: Case, comp: []const C, want: []const [2]f64, size: []const [2]f64, fac: f64, old: bool) !void {
    const F = Float(C);
    var sc: [40]F = undefined;
    var st: [40]f64 = undefined;
    var ss: [40]f64 = undefined;
    for (comp, 0..) |z, i| {
        sc[2 * i] = z.re;
        sc[2 * i + 1] = z.im;
        st[2 * i] = want[i][0];
        st[2 * i + 1] = want[i][1];
        ss[2 * i] = size[i][0];
        ss[2 * i + 1] = size[i][1];
    }
    const len = 2 * comp.len;
    try stest(F, c, sc[0..len], st[0..len], ss[0..len], fac, old);
}

fn itest1(c: Case, comp: isize, want: isize) !void {
    if (comp == want) return;
    ob.say("{s} n={d} incx={d}: got {d}, true {d}\n", .{ c.name, c.n, c.incx, comp, want });
    return error.TestUnexpectedResult;
}

fn load(comptime F: type, dst: []F, src: []const f64) void {
    for (dst, src[0..dst.len]) |*d, s| d.* = @floatCast(s);
}

fn loadC(comptime C: type, dst: []C, src: []const [2]f64) void {
    for (dst, src[0..dst.len]) |*d, s| d.* = .init(@floatCast(s[0]), @floatCast(s[1]));
}

// ---- real --------------------------------------------------------------------

/// `CHECK0`: `rotg` and `rotmg`.
fn check0(comptime F: type) !void {
    const p = if (F == f32) "S" else "D";
    var dbt = dbtrue;
    dbt[0] = 1.0 / 0.6;
    dbt[2] = -1.0 / 0.6;
    dbt[4] = 1.0 / 0.6;
    for (0..8) |k| {
        const c: Case = .{ .name = p ++ "ROTG", .n = @intCast(k + 1) };
        var sa: F = @floatCast(da1[k]);
        var sb: F = @floatCast(db1[k]);
        var sc: F = undefined;
        var ss: F = undefined;
        blas.rotg(F, &sa, &sb, &sc, &ss);
        try stest1(F, c, sa, datrue[k], datrue[k], sfac);
        try stest1(F, c, sb, dbt[k], dbt[k], sfac);
        try stest1(F, c, sc, dc1[k], dc1[k], sfac);
        try stest1(F, c, ss, ds1[k], ds1[k], sfac);
    }
    // DTRUE's assignments after its DATA, in the working precision
    var dt: [81]F = undefined;
    load(F, &dt, &dtrue);
    const d12: F = 4096;
    const T = struct {
        fn at(a: *[81]F, i: usize, j: usize) *F {
            return &a[(i - 1) + 9 * (j - 1)];
        }
    }.at;
    T(&dt, 1, 1).* = 12.0 / 130.0;
    T(&dt, 2, 1).* = 36.0 / 130.0;
    T(&dt, 7, 1).* = -1.0 / 6.0;
    T(&dt, 1, 2).* = 14.0 / 75.0;
    T(&dt, 2, 2).* = 49.0 / 75.0;
    T(&dt, 9, 2).* = 1.0 / 7.0;
    T(&dt, 1, 5).* = 45e-11 * (d12 * d12);
    T(&dt, 3, 5).* = 4e5 / (3 * d12);
    T(&dt, 6, 5).* = 1 / d12;
    T(&dt, 8, 5).* = 1e4 / (3 * d12);
    T(&dt, 1, 6).* = 4e10 / (1.5 * d12 * d12);
    T(&dt, 2, 6).* = 2e-2 / 1.5;
    T(&dt, 8, 6).* = 5e-7 * d12;
    T(&dt, 1, 7).* = 4.0 / 150.0;
    T(&dt, 2, 7).* = (2e-10 / 1.5) * (d12 * d12);
    T(&dt, 7, 7).* = -T(&dt, 6, 5).*;
    T(&dt, 9, 7).* = 1e4 / d12;
    T(&dt, 1, 8).* = T(&dt, 1, 7).*;
    T(&dt, 2, 8).* = 2e10 / (1.5 * d12 * d12);
    T(&dt, 1, 9).* = 32.0 / 7.0;
    T(&dt, 2, 9).* = -16.0 / 7.0;
    for (0..8) |k| {
        const c: Case = .{ .name = p ++ "ROTMG", .n = @intCast(k + 1) };
        var t: [9]F = @splat(0);
        load(F, t[0..4], dab[4 * k ..][0..4]);
        blas.rotmg(F, &t[0], &t[1], &t[2], t[3], t[4..9]);
        var want: [9]f64 = undefined;
        for (&want, dt[9 * k ..][0..9]) |*w, v| w.* = v;
        try stest(F, c, &t, &want, &want, sfac, false);
    }
}

/// `CHECK1`: `nrm2`, `asum`, `scal`, `iamax`.
fn check1(comptime F: type) !void {
    const p = if (F == f32) "S" else "D";
    for (1..3) |incx| for (0..5) |np1| {
        const n = np1;
        const len = 2 * @max(n, 1);
        const col = 8 * (np1 + 5 * (incx - 1));
        var sx: [8]F = undefined;
        const inc: isize = @intCast(incx);
        const in: isize = @intCast(n);
        load(F, sx[0..len], dv[col..]);
        try stest1(F, .{ .name = p ++ "NRM2", .n = in, .incx = inc }, try blas.nrm2(ob.io, F, in, &sx, inc), dtrue1[np1], dtrue1[np1], sfac);
        try stest1(F, .{ .name = p ++ "ASUM", .n = in, .incx = inc }, try blas.asum(ob.io, F, in, &sx, inc), dtrue3[np1], dtrue3[np1], sfac);
        try blas.scal(ob.io, F, in, @floatCast(sa1[(incx - 1) * 5 + np1]), &sx, inc);
        try stest(F, .{ .name = p ++ "SCAL", .n = in, .incx = inc }, sx[0..len], dtrue5[col..][0..len], dtrue5[col..][0..len], sfac, false);
        load(F, sx[0..len], dv[col..]);
        try itest1(.{ .name = "I" ++ p ++ "AMAX", .n = in, .incx = inc }, try blas.index(ob.io, F, .amax, in, &sx, inc), itrue2[np1]);
    };
}

const incxs = [_]isize{ 1, 2, -2, -1 };
const incys = [_]isize{ 1, -2, 1, -2 };
const ns = [_]isize{ 0, 1, 2, 4 };
/// `LENS(KN, |inc|)`: how many elements a vector of `NS(KN)` spans.
const lens = [_]usize{ 1, 1, 2, 4, 1, 1, 3, 7 };

/// `CHECK2`: `dot`, `axpy`, `copy`, `swap`, `rotm`, and `sdsdot` (single)
/// or `dsdot` (double).
fn check2(comptime F: type) !void {
    const p = if (F == f32) "S" else "D";
    const dt19x = dt19xa ++ dt19xb ++ dt19xc ++ dt19xd;
    const dt19y = dt19ya ++ dt19yb ++ dt19yc ++ dt19yd;
    const zeros: [7]f64 = @splat(0);
    for (0..4) |ki| {
        const incx = incxs[ki];
        const incy = incys[ki];
        const mx: usize = @abs(incx);
        const my: usize = @abs(incy);
        for (0..4) |kn| {
            const n = ns[kn];
            const ksize: usize = @min(2, kn + 1);
            const lenx = lens[kn + 4 * (mx - 1)];
            const leny = lens[kn + 4 * (my - 1)];
            const at = 7 * (kn + 4 * ki);
            const c: Case = .{ .name = p, .n = n, .incx = incx, .incy = incy };
            var sx: [7]F = undefined;
            var sy: [7]F = undefined;
            load(F, &sx, &dx1);
            load(F, &sy, &dy1);
            try stest1(F, c, try blas.dot(ob.io, F, n, &sx, incx, &sy, incy), dt7[kn + 4 * ki], ssize1[kn], sfac);
            try blas.axpy(ob.io, F, n, 0.3, &sx, incx, &sy, incy);
            try stest(F, c, sy[0..leny], dt8[at..], ssize2[14 * (ksize - 1) ..], sfac, false);
            load(F, &sy, &dy1);
            try blas.copy(ob.io, F, n, &sx, incx, &sy, incy);
            try stest(F, c, sy[0..leny], dt10y[at..], &zeros, 1, false);
            load(F, &sy, &dy1);
            try blas.swap(ob.io, F, n, &sx, incx, &sy, incy);
            try stest(F, c, sx[0..lenx], dt10x[at..], &zeros, 1, false);
            try stest(F, c, sy[0..leny], dt10y[at..], &zeros, 1, false);
            const kni = kn + 4 * ki;
            for (0..4) |kpar| {
                load(F, &sx, &dx1);
                load(F, &sy, &dy1);
                const stx = dt19x[7 * kpar + 28 * kni ..][0..7];
                const sty = dt19y[7 * kpar + 28 * kni ..][0..7];
                var param: [5]F = undefined;
                load(F, &param, dpar[5 * kpar ..]);
                var ssize: [7]f64 = stx.*;
                if (kpar == 1 and kni == 6) ssize[0] = 2.4;
                if (kpar == 2 and kni == 7) ssize[4] = 1.8;
                try blas.rotm(ob.io, F, n, &sx, incx, &sy, incy, &param);
                try stest(F, .{ .name = p ++ "ROTM", .n = n, .incx = incx, .incy = incy }, sx[0..lenx], stx, &ssize, sfac, false);
                try stest(F, .{ .name = p ++ "ROTM", .n = n, .incx = incx, .incy = incy }, sy[0..leny], sty, sty, sfac, false);
            }
            var fx: [7]f32 = undefined;
            var fy: [7]f32 = undefined;
            load(f32, &fx, &dx1);
            load(f32, &fy, &dy1);
            if (F == f32) {
                try stest1(f32, .{ .name = "SDSDOT", .n = n, .incx = incx, .incy = incy }, try blas.sdsdot(ob.io, n, 0.1, &fx, incx, &fy, incy), st7b[kn + 4 * ki], ssize3[kn], sfac);
            } else {
                // TESTDSDOT: compared in single precision, with its own factor
                const d: f32 = @floatCast(try blas.dsdot(ob.io, n, &fx, incx, &fy, incy));
                try stest1(f32, .{ .name = "DSDOT", .n = n, .incx = incx, .incy = incy }, d, dt7[kn + 4 * ki], ssize1[kn], 0.3125e-1);
            }
        }
    }
}

/// `CHECK3`: `rot`.
fn check3(comptime F: type) !void {
    const p = if (F == f32) "SROT" else "DROT";
    for (0..4) |ki| {
        const incx = incxs[ki];
        const incy = incys[ki];
        for (0..4) |kn| {
            const n = ns[kn];
            const ksize: usize = @min(2, kn + 1);
            const lenx = lens[kn + 4 * (@abs(incx) - 1)];
            const leny = lens[kn + 4 * (@abs(incy) - 1)];
            const at = 7 * (kn + 4 * ki);
            var sx: [7]F = undefined;
            var sy: [7]F = undefined;
            load(F, &sx, &dx1);
            load(F, &sy, &dy1);
            try blas.rot(ob.io, F, n, &sx, incx, &sy, incy, 0.8, 0.6);
            const c: Case = .{ .name = p, .n = n, .incx = incx, .incy = incy };
            try stest(F, c, sx[0..lenx], dt9x[at..], ssize2[14 * (ksize - 1) ..], sfac, false);
            try stest(F, c, sy[0..leny], dt9y[at..], ssize2[14 * (ksize - 1) ..], sfac, false);
        }
    }
    // the MWP cases: rotations by 0 and +-90 degrees, strides of every sign
    const mwpc = [11]F{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    const mwps = [11]F{ 0, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1 };
    const mwpinx = [11]isize{ 1, 1, 1, -1, 1, -1, 1, 1, -1, 1, -1 };
    const mwpiny = [11]isize{ 1, 1, -1, -1, 2, 1, 1, -1, -1, 2, 1 };
    const mwpn = [11]isize{ 5, 5, 5, 5, 3, 5, 5, 5, 5, 3, 5 };
    var tx: [11][5]f64 = undefined;
    var ty: [11][5]f64 = undefined;
    for (0..5) |ii| {
        const i: f64 = @floatFromInt(ii + 1);
        tx[0][ii] = i;
        ty[0][ii] = i;
        tx[1][ii] = i;
        ty[1][ii] = -i;
        tx[2][ii] = 6 - i;
        ty[2][ii] = i - 6;
        tx[3][ii] = i;
        ty[3][ii] = -i;
        tx[5][ii] = 6 - i;
        ty[5][ii] = i - 6;
        tx[6][ii] = -i;
        ty[6][ii] = i;
        tx[7][ii] = i - 6;
        ty[7][ii] = 6 - i;
        tx[8][ii] = -i;
        ty[8][ii] = i;
        tx[10][ii] = i - 6;
        ty[10][ii] = 6 - i;
    }
    tx[4] = .{ 1, 3, 5, 4, 5 };
    ty[4] = .{ -1, 2, -2, 4, -3 };
    tx[9] = .{ -1, -3, -5, 4, 5 };
    ty[9] = .{ 1, 2, 2, 4, 3 };
    for (0..11) |i| {
        var x = [5]F{ 1, 2, 3, 4, 5 };
        var y = [5]F{ 1, 2, 3, 4, 5 };
        try blas.rot(ob.io, F, mwpn[i], &x, mwpinx[i], &y, mwpiny[i], mwpc[i], mwps[i]);
        const c: Case = .{ .name = p, .n = mwpn[i], .incx = mwpinx[i], .incy = mwpiny[i] };
        try stest(F, c, &x, &tx[i], &tx[i], sfac, false);
        try stest(F, c, &y, &ty[i], &ty[i], sfac, false);
    }
}

/// ctest's `SROT`/`DROT`: the reference BLAS loop.
fn refRot(comptime F: type, n: isize, x: []F, incx: isize, y: []F, incy: isize, c: F, s: F) void {
    if (n <= 0) return;
    var ix: isize = if (incx < 0) (1 - n) * incx else 0;
    var iy: isize = if (incy < 0) (1 - n) * incy else 0;
    for (0..@intCast(n)) |_| {
        const xi: usize = @intCast(ix);
        const yi: usize = @intCast(iy);
        const t = c * x[xi] + s * y[yi];
        y[yi] = c * y[yi] - s * x[xi];
        x[xi] = t;
        ix += incx;
        iy += incy;
    }
}

/// ctest's `SROTM`/`DROTM`: the reference BLAS loop.
fn refRotm(comptime F: type, n: isize, x: []F, incx: isize, y: []F, incy: isize, p: *const [5]F) void {
    const flag = p[0];
    if (n <= 0 or flag + 2 == 0) return;
    var kx: isize = if (incx < 0) (1 - n) * incx else 0;
    var ky: isize = if (incy < 0) (1 - n) * incy else 0;
    for (0..@intCast(n)) |_| {
        const xi: usize = @intCast(kx);
        const yi: usize = @intCast(ky);
        const w = x[xi];
        const z = y[yi];
        if (flag < 0) {
            x[xi] = w * p[1] + z * p[3];
            y[yi] = w * p[2] + z * p[4];
        } else if (flag == 0) {
            x[xi] = w + z * p[3];
            y[yi] = w * p[2] + z;
        } else {
            x[xi] = w * p[1] + z;
            y[yi] = -w + p[4] * z;
        }
        kx += incx;
        ky += incy;
    }
}

/// ctest's `CHECK3`: `rot` and `rotm` against the reference loops, over the
/// whole of `x` and `y`.
fn ctestCheck3(comptime F: type) !void {
    const dx = if (F == f32) &c_sdx else &c_ddx;
    const dy = if (F == f32) &c_sdy else &c_ddy;
    const nss = if (F == f32) &c_sns else &c_dns;
    const size = if (F == f32) &c_sssize2 else &c_dssize2;
    const len = dx.len;
    const p = if (F == f32) "S" else "D";
    for (c_incxs, c_incys) |incx, incy| for (nss, 0..) |n, kn| {
        const ksize: usize = @min(2, kn + 1);
        const ssz = size[len * (ksize - 1) ..][0..len];
        var sx: [len]F = undefined;
        var sy: [len]F = undefined;
        load(F, &sx, dx);
        load(F, &sy, dy);
        var tx = sx;
        var ty = sy;
        const c: Case = .{ .name = p ++ "ROT", .n = n, .incx = incx, .incy = incy };
        try blas.rot(ob.io, F, n, &sx, incx, &sy, incy, 0.8, 0.6);
        refRot(F, n, &tx, incx, &ty, incy, 0.8, 0.6);
        try stest(F, c, &sx, &wideOf(F, len, tx), ssz, sfac, true);
        try stest(F, c, &sy, &wideOf(F, len, ty), ssz, sfac, true);
        load(F, &sx, dx);
        load(F, &sy, dy);
        tx = sx;
        ty = sy;
        for (0..4) |i| {
            var param: [5]F = undefined;
            load(F, &param, c_param[5 * i ..]);
            const cm: Case = .{ .name = p ++ "ROTM", .n = n, .incx = incx, .incy = incy };
            try blas.rotm(ob.io, F, n, &sx, incx, &sy, incy, &param);
            refRotm(F, n, &tx, incx, &ty, incy, &param);
            try stest(F, cm, &sx, &wideOf(F, len, tx), ssz, sfac, true);
            try stest(F, cm, &sy, &wideOf(F, len, ty), ssz, sfac, true);
        }
    };
}

fn wideOf(comptime F: type, comptime n: usize, v: [n]F) [n]f64 {
    var out: [n]f64 = undefined;
    for (&out, v) |*o, x| o.* = x;
    return out;
}

// ---- complex -----------------------------------------------------------------

/// `CHECK1`: `nrm2`, `asum`, `scal`, `scalReal`, `iamax`, from the data
/// `d` (the Fortran's or ctest's), compared by `old` as `STEST` says.
fn ccheck1(comptime C: type, d: anytype, old: bool) !void {
    const F = Float(C);
    const p = if (F == f32) "C" else "Z";
    const ca: C = .init(0.4, -0.7);
    const sa: F = 0.3;
    var cx: [8]C = undefined;
    inline for (.{ "scal", "scalReal" }) |which| {
        for (1..3) |incx| for (0..5) |np1| {
            const n = np1;
            const len = 2 * @max(n, 1);
            const col = 8 * (np1 + 5 * (incx - 1));
            const inc: isize = @intCast(incx);
            const in: isize = @intCast(n);
            const c: Case = .{ .name = p, .n = in, .incx = inc };
            loadC(C, cx[0..len], d.cv[col..]);
            try stest(F, c, &.{try blas.nrm2(ob.io, C, in, &cx, inc)}, d.strue2[np1..], d.strue2[np1..], sfac, old);
            try stest(F, c, &.{try blas.asum(ob.io, C, in, &cx, inc)}, d.strue4[np1..], d.strue4[np1..], sfac, old);
            try itest1(c, try blas.index(ob.io, C, .amax, in, &cx, inc), itrue3[np1]);
            if (comptime std.mem.eql(u8, which, "scal")) {
                try blas.scal(ob.io, C, in, ca, &cx, inc);
                try ctest(C, c, cx[0..len], d.ctrue5[col..], d.ctrue5[col..], sfac, old);
            } else {
                try blas.scalReal(ob.io, C, in, sa, &cx, inc);
                try ctest(C, c, cx[0..len], d.ctrue6[col..], d.ctrue6[col..], sfac, old);
            }
        };
        // then on the last `x`, with `incx` 1: `alpha` 0, and for the real
        // `alpha` 1 and -1
        const c: Case = .{ .name = p ++ "SCAL", .n = 5, .incx = 1 };
        const zeros: [5][2]f64 = @splat(.{ 0, 0 });
        const ones: [5][2]f64 = @splat(.{ 1, 1 });
        if (comptime std.mem.eql(u8, which, "scal")) {
            try blas.scal(ob.io, C, 5, .init(0, 0), &cx, 1);
            try ctest(C, c, cx[0..5], &zeros, &ones, sfac, old);
        } else {
            try blas.scalReal(ob.io, C, 5, 0, &cx, 1);
            try ctest(C, c, cx[0..5], &zeros, &ones, sfac, old);
            for ([_]F{ 1, -1 }) |s| {
                var want: [5][2]f64 = undefined;
                for (&want, cx[0..5]) |*w, z| w.* = .{ s * z.re, s * z.im };
                try blas.scalReal(ob.io, C, 5, s, &cx, 1);
                try ctest(C, c, cx[0..5], &want, &want, sfac, old);
            }
        }
    }
}

/// `CHECK2`: `dotc`, `dotu`, `axpy`, `copy`, `swap`.
fn ccheck2(comptime C: type) !void {
    const p = if (Float(C) == f32) "C" else "Z";
    const ca: C = .init(0.4, -0.7);
    for (0..4) |ki| {
        const incx = incxs[ki];
        const incy = incys[ki];
        for (0..4) |kn| {
            const n = ns[kn];
            const ksize: usize = @min(2, kn + 1);
            const lenx = lens[kn + 4 * (@abs(incx) - 1)];
            const leny = lens[kn + 4 * (@abs(incy) - 1)];
            const at = 7 * (kn + 4 * ki);
            const c: Case = .{ .name = p, .n = n, .incx = incx, .incy = incy };
            var cx: [7]C = undefined;
            var cy: [7]C = undefined;
            loadC(C, &cx, &cx1);
            loadC(C, &cy, &cy1);
            try ctest(C, c, &.{try blas.dotc(ob.io, C, n, &cx, incx, &cy, incy)}, ct6[kn + 4 * ki ..], csize1[kn..], sfac, false);
            try ctest(C, c, &.{try blas.dot(ob.io, C, n, &cx, incx, &cy, incy)}, ct7[kn + 4 * ki ..], csize1[kn..], sfac, false);
            try blas.axpy(ob.io, C, n, ca, &cx, incx, &cy, incy);
            try ctest(C, c, cy[0..leny], ct8[at..], csize2[7 * (ksize - 1) ..], sfac, false);
            loadC(C, &cy, &cy1);
            try blas.copy(ob.io, C, n, &cx, incx, &cy, incy);
            try ctest(C, c, cy[0..leny], ct10y[at..], &csize3, 1, false);
            loadC(C, &cy, &cy1);
            try blas.swap(ob.io, C, n, &cx, incx, &cy, incy);
            try ctest(C, c, cx[0..lenx], ct10x[at..], &csize3, 1, false);
            try ctest(C, c, cy[0..leny], ct10y[at..], &csize3, 1, false);
        }
    }
}

const fortran_c1 = .{ .cv = &cv, .strue2 = &strue2, .strue4 = &strue4, .ctrue5 = &ctrue5, .ctrue6 = &ctrue6 };
const ctest_c1 = .{ .cv = &c_cv, .strue2 = &c_strue2, .strue4 = &c_strue4, .ctrue5 = &c_ctrue5, .ctrue6 = &c_ctrue6 };

/// Every case here, for `check.eachCore`.
pub fn all() !void {
    inline for (.{ f32, f64 }) |F| {
        try check0(F);
        try check1(F);
        try check2(F);
        try check3(F);
        try ctestCheck3(F);
    }
    inline for (.{ C32, C64 }) |C| {
        try ccheck1(C, fortran_c1, false);
        try ccheck1(C, ctest_c1, true);
        try ccheck2(C);
    }
}

test "blat1: xROTG, xROTMG" {
    try check0(f32);
    try check0(f64);
}

test "blat1: xNRM2, xASUM, xSCAL, IxAMAX" {
    try check1(f32);
    try check1(f64);
}

test "blat1: xDOT, xAXPY, xCOPY, xSWAP, xROTM, SDSDOT, DSDOT" {
    try check2(f32);
    try check2(f64);
}

test "blat1: xROT" {
    try check3(f32);
    try check3(f64);
}

test "blat1: ctest xROT, xROTM against the reference loops" {
    try ctestCheck3(f32);
    try ctestCheck3(f64);
}

test "blat1: xNRM2, xASUM, xSCAL, xSSCAL, IxAMAX (complex, and ctest's data)" {
    inline for (.{ C32, C64 }) |C| {
        try ccheck1(C, fortran_c1, false);
        try ccheck1(C, ctest_c1, true);
    }
}

test "blat1: xDOTC, xDOTU, xAXPY, xCOPY, xSWAP (complex)" {
    try ccheck2(C32);
    try ccheck2(C64);
}

// ---- data --------------------------------------------------------------------

const da1 = [_]f64{
    0.3, 0.4, -0.3, -0.4, -0.3, 0.0, 0.0, 1.0,
};
const db1 = [_]f64{
    0.4, 0.3, 0.4, 0.3, -0.4, 0.0, 1.0, 0.0,
};
const dc1 = [_]f64{
    0.6, 0.8, -0.6, 0.8, 0.6, 1.0, 0.0, 1.0,
};
const ds1 = [_]f64{
    0.8, 0.6, 0.8, -0.6, 0.8, 0.0, 1.0, 0.0,
};
const datrue = [_]f64{
    0.5, 0.5, 0.5, -0.5, -0.5, 0.0, 1.0, 1.0,
};
const dbtrue = [_]f64{
    0.0, 0.6, 0.0, -0.6, 0.0, 0.0, 1.0, 0.0,
};
const dab = [_]f64{
    0.1,           0.3,  1.2,      0.2,
    0.7,           0.2,  0.6,      4.2,
    0.0,           0.0,  0.0,      0.0,
    4.0,           -1.0, 2.0,      4.0,
    6e-10,         0.02, 100000.0, 10.0,
    40000000000.0, 0.02, 1e-05,    10.0,
    2e-10,         0.04, 100000.0, 10.0,
    20000000000.0, 0.04, 1e-05,    10.0,
    4.0,           -2.0, 8.0,      4.0,
};
const dtrue = [_]f64{
    0.0, 0.0,   1.3,     0.2,  0.0,  0.0,      0.0,        0.5,   0.0,
    0.0, 0.0,   4.5,     4.2,  1.0,  0.5,      0.0,        0.0,   0.0,
    0.0, 0.0,   0.0,     0.0,  -2.0, 0.0,      0.0,        0.0,   0.0,
    0.0, 0.0,   0.0,     4.0,  -1.0, 0.0,      0.0,        0.0,   0.0,
    0.0, 0.015, 0.0,     10.0, -1.0, 0.0,      -0.0001,    0.0,   1.0,
    0.0, 0.0,   0.06144, 10.0, -1.0, 4096.0,   -1000000.0, 0.0,   1.0,
    0.0, 0.0,   15.0,    10.0, -1.0, 5e-05,    0.0,        1.0,   0.0,
    0.0, 0.0,   15.0,    10.0, -1.0, 500000.0, -4096.0,    1.0,   0.004096,
    0.0, 0.0,   7.0,     4.0,  0.0,  0.0,      -0.5,       -0.25, 0.0,
};
const sa1 = [_]f64{
    0.3, -1.0, 0.0, 1.0, 0.3,
    0.3, 0.3,  0.3, 0.3, 0.3,
};
const dv = [_]f64{
    0.1, 2.0,  2.0,  2.0,  2.0,  2.0, 2.0,  2.0,
    0.3, 3.0,  3.0,  3.0,  3.0,  3.0, 3.0,  3.0,
    0.3, -0.4, 4.0,  4.0,  4.0,  4.0, 4.0,  4.0,
    0.2, -0.6, 0.3,  5.0,  5.0,  5.0, 5.0,  5.0,
    0.1, -0.3, 0.5,  -0.1, 6.0,  6.0, 6.0,  6.0,
    0.1, 8.0,  8.0,  8.0,  8.0,  8.0, 8.0,  8.0,
    0.3, 9.0,  9.0,  9.0,  9.0,  9.0, 9.0,  9.0,
    0.3, 2.0,  -0.4, 2.0,  2.0,  2.0, 2.0,  2.0,
    0.2, 3.0,  -0.6, 5.0,  0.3,  2.0, 2.0,  2.0,
    0.1, 4.0,  -0.3, 6.0,  -0.5, 7.0, -0.1, 3.0,
};
const dtrue1 = [_]f64{
    0.0, 0.3, 0.5, 0.7, 0.6,
};
const dtrue3 = [_]f64{
    0.0, 0.3, 0.7, 1.1, 1.0,
};
const dtrue5 = [_]f64{
    0.1,  2.0,   2.0,   2.0,   2.0,   2.0, 2.0,   2.0,
    -0.3, 3.0,   3.0,   3.0,   3.0,   3.0, 3.0,   3.0,
    0.0,  0.0,   4.0,   4.0,   4.0,   4.0, 4.0,   4.0,
    0.2,  -0.6,  0.3,   5.0,   5.0,   5.0, 5.0,   5.0,
    0.03, -0.09, 0.15,  -0.03, 6.0,   6.0, 6.0,   6.0,
    0.1,  8.0,   8.0,   8.0,   8.0,   8.0, 8.0,   8.0,
    0.09, 9.0,   9.0,   9.0,   9.0,   9.0, 9.0,   9.0,
    0.09, 2.0,   -0.12, 2.0,   2.0,   2.0, 2.0,   2.0,
    0.06, 3.0,   -0.18, 5.0,   0.09,  2.0, 2.0,   2.0,
    0.03, 4.0,   -0.09, 6.0,   -0.15, 7.0, -0.03, 3.0,
};
const itrue2 = [_]isize{
    0, 1, 2, 2, 3,
};
const dx1 = [_]f64{
    0.6, 0.1, -0.5, 0.8, 0.9, -0.3, -0.4,
};
const dy1 = [_]f64{
    0.5, -0.9, 0.3, 0.7, -0.6, 0.2, 0.8,
};
const dt7 = [_]f64{
    0.0, 0.3, 0.21,  0.62,
    0.0, 0.3, -0.07, 0.85,
    0.0, 0.3, -0.79, -0.74,
    0.0, 0.3, 0.33,  1.27,
};
const st7b = [_]f64{
    0.1, 0.4, 0.31,  0.72,
    0.1, 0.4, 0.03,  0.95,
    0.1, 0.4, -0.69, -0.64,
    0.1, 0.4, 0.43,  1.37,
};
const dt8 = [_]f64{
    0.5,  0.0,   0.0,  0.0,  0.0,   0.0, 0.0,
    0.68, 0.0,   0.0,  0.0,  0.0,   0.0, 0.0,
    0.68, -0.87, 0.0,  0.0,  0.0,   0.0, 0.0,
    0.68, -0.87, 0.15, 0.94, 0.0,   0.0, 0.0,
    0.5,  0.0,   0.0,  0.0,  0.0,   0.0, 0.0,
    0.68, 0.0,   0.0,  0.0,  0.0,   0.0, 0.0,
    0.35, -0.9,  0.48, 0.0,  0.0,   0.0, 0.0,
    0.38, -0.9,  0.57, 0.7,  -0.75, 0.2, 0.98,
    0.5,  0.0,   0.0,  0.0,  0.0,   0.0, 0.0,
    0.68, 0.0,   0.0,  0.0,  0.0,   0.0, 0.0,
    0.35, -0.72, 0.0,  0.0,  0.0,   0.0, 0.0,
    0.38, -0.63, 0.15, 0.88, 0.0,   0.0, 0.0,
    0.5,  0.0,   0.0,  0.0,  0.0,   0.0, 0.0,
    0.68, 0.0,   0.0,  0.0,  0.0,   0.0, 0.0,
    0.68, -0.9,  0.33, 0.0,  0.0,   0.0, 0.0,
    0.68, -0.9,  0.33, 0.7,  -0.75, 0.2, 1.04,
};
const dt10x = [_]f64{
    0.6,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,
    0.5,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,
    0.5,  -0.9, 0.0,  0.0, 0.0,  0.0,  0.0,
    0.5,  -0.9, 0.3,  0.7, 0.0,  0.0,  0.0,
    0.6,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,
    0.5,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,
    0.3,  0.1,  0.5,  0.0, 0.0,  0.0,  0.0,
    0.8,  0.1,  -0.6, 0.8, 0.3,  -0.3, 0.5,
    0.6,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,
    0.5,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,
    -0.9, 0.1,  0.5,  0.0, 0.0,  0.0,  0.0,
    0.7,  0.1,  0.3,  0.8, -0.9, -0.3, 0.5,
    0.6,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,
    0.5,  0.0,  0.0,  0.0, 0.0,  0.0,  0.0,
    0.5,  0.3,  0.0,  0.0, 0.0,  0.0,  0.0,
    0.5,  0.3,  -0.6, 0.8, 0.0,  0.0,  0.0,
};
const dt10y = [_]f64{
    0.5,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,
    0.6,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,
    0.6,  0.1,  0.0,  0.0, 0.0,  0.0, 0.0,
    0.6,  0.1,  -0.5, 0.8, 0.0,  0.0, 0.0,
    0.5,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,
    0.6,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,
    -0.5, -0.9, 0.6,  0.0, 0.0,  0.0, 0.0,
    -0.4, -0.9, 0.9,  0.7, -0.5, 0.2, 0.6,
    0.5,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,
    0.6,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,
    -0.5, 0.6,  0.0,  0.0, 0.0,  0.0, 0.0,
    -0.4, 0.9,  -0.5, 0.6, 0.0,  0.0, 0.0,
    0.5,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,
    0.6,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,
    0.6,  -0.9, 0.1,  0.0, 0.0,  0.0, 0.0,
    0.6,  -0.9, 0.1,  0.7, -0.5, 0.2, 0.8,
};
const ssize1 = [_]f64{
    0.0, 0.3, 1.6, 3.2,
};
const ssize2 = [_]f64{
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17,
};
const ssize3 = [_]f64{
    0.1, 0.4, 1.7, 3.3,
};
const dpar = [_]f64{
    -2.0, 0.0, 0.0,  0.0,  0.0,
    -1.0, 2.0, -3.0, -4.0, 5.0,
    0.0,  0.0, 2.0,  -3.0, 0.0,
    1.0,  5.0, 2.0,  0.0,  -4.0,
};
const dt19xa = [_]f64{
    0.6,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    0.6,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    0.6,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    0.6,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    0.6,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    -0.8, 0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    -0.9, 0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    3.5,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    0.6,  0.1,  0.0,  0.0,  0.0, 0.0, 0.0,
    -0.8, 3.8,  0.0,  0.0,  0.0, 0.0, 0.0,
    -0.9, 2.8,  0.0,  0.0,  0.0, 0.0, 0.0,
    3.5,  -0.4, 0.0,  0.0,  0.0, 0.0, 0.0,
    0.6,  0.1,  -0.5, 0.8,  0.0, 0.0, 0.0,
    -0.8, 3.8,  -2.2, -1.2, 0.0, 0.0, 0.0,
    -0.9, 2.8,  -1.4, -1.3, 0.0, 0.0, 0.0,
    3.5,  -0.4, -2.2, 4.7,  0.0, 0.0, 0.0,
};
const dt19xb = [_]f64{
    0.6,  0.0, 0.0,  0.0, 0.0, 0.0,  0.0,
    0.6,  0.0, 0.0,  0.0, 0.0, 0.0,  0.0,
    0.6,  0.0, 0.0,  0.0, 0.0, 0.0,  0.0,
    0.6,  0.0, 0.0,  0.0, 0.0, 0.0,  0.0,
    0.6,  0.0, 0.0,  0.0, 0.0, 0.0,  0.0,
    -0.8, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0,
    -0.9, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0,
    3.5,  0.0, 0.0,  0.0, 0.0, 0.0,  0.0,
    0.6,  0.1, -0.5, 0.0, 0.0, 0.0,  0.0,
    0.0,  0.1, -3.0, 0.0, 0.0, 0.0,  0.0,
    -0.3, 0.1, -2.0, 0.0, 0.0, 0.0,  0.0,
    3.3,  0.1, -2.0, 0.0, 0.0, 0.0,  0.0,
    0.6,  0.1, -0.5, 0.8, 0.9, -0.3, -0.4,
    -2.0, 0.1, 1.4,  0.8, 0.6, -0.3, -2.8,
    -1.8, 0.1, 1.3,  0.8, 0.0, -0.3, -1.9,
    3.8,  0.1, -3.1, 0.8, 4.8, -0.3, -1.5,
};
const dt19xc = [_]f64{
    0.6,  0.0, 0.0,  0.0, 0.0, 0.0,  0.0,
    0.6,  0.0, 0.0,  0.0, 0.0, 0.0,  0.0,
    0.6,  0.0, 0.0,  0.0, 0.0, 0.0,  0.0,
    0.6,  0.0, 0.0,  0.0, 0.0, 0.0,  0.0,
    0.6,  0.0, 0.0,  0.0, 0.0, 0.0,  0.0,
    -0.8, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0,
    -0.9, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0,
    3.5,  0.0, 0.0,  0.0, 0.0, 0.0,  0.0,
    0.6,  0.1, -0.5, 0.0, 0.0, 0.0,  0.0,
    4.8,  0.1, -3.0, 0.0, 0.0, 0.0,  0.0,
    3.3,  0.1, -2.0, 0.0, 0.0, 0.0,  0.0,
    2.1,  0.1, -2.0, 0.0, 0.0, 0.0,  0.0,
    0.6,  0.1, -0.5, 0.8, 0.9, -0.3, -0.4,
    -1.6, 0.1, -2.2, 0.8, 5.4, -0.3, -2.8,
    -1.5, 0.1, -1.4, 0.8, 3.6, -0.3, -1.9,
    3.7,  0.1, -2.2, 0.8, 3.6, -0.3, -1.5,
};
const dt19xd = [_]f64{
    0.6,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    0.6,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    0.6,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    0.6,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    0.6,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    -0.8, 0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    -0.9, 0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    3.5,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    0.6,  0.1,  0.0,  0.0,  0.0, 0.0, 0.0,
    -0.8, -1.0, 0.0,  0.0,  0.0, 0.0, 0.0,
    -0.9, -0.8, 0.0,  0.0,  0.0, 0.0, 0.0,
    3.5,  0.8,  0.0,  0.0,  0.0, 0.0, 0.0,
    0.6,  0.1,  -0.5, 0.8,  0.0, 0.0, 0.0,
    -0.8, -1.0, 1.4,  -1.6, 0.0, 0.0, 0.0,
    -0.9, -0.8, 1.3,  -1.6, 0.0, 0.0, 0.0,
    3.5,  0.8,  -3.1, 4.8,  0.0, 0.0, 0.0,
};
const dt19ya = [_]f64{
    0.5,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    0.5,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    0.5,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    0.5,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    0.5,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    0.7,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    1.7,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    -2.6, 0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    0.5,  -0.9, 0.0,  0.0,  0.0, 0.0, 0.0,
    0.7,  -4.8, 0.0,  0.0,  0.0, 0.0, 0.0,
    1.7,  -0.7, 0.0,  0.0,  0.0, 0.0, 0.0,
    -2.6, 3.5,  0.0,  0.0,  0.0, 0.0, 0.0,
    0.5,  -0.9, 0.3,  0.7,  0.0, 0.0, 0.0,
    0.7,  -4.8, 3.0,  1.1,  0.0, 0.0, 0.0,
    1.7,  -0.7, -0.7, 2.3,  0.0, 0.0, 0.0,
    -2.6, 3.5,  -0.7, -3.6, 0.0, 0.0, 0.0,
};
const dt19yb = [_]f64{
    0.5,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,
    0.5,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,
    0.5,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,
    0.5,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,
    0.5,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,
    0.7,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,
    1.7,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,
    -2.6, 0.0,  0.0,  0.0, 0.0,  0.0, 0.0,
    0.5,  -0.9, 0.3,  0.0, 0.0,  0.0, 0.0,
    4.0,  -0.9, -0.3, 0.0, 0.0,  0.0, 0.0,
    -0.5, -0.9, 1.5,  0.0, 0.0,  0.0, 0.0,
    -1.5, -0.9, -1.8, 0.0, 0.0,  0.0, 0.0,
    0.5,  -0.9, 0.3,  0.7, -0.6, 0.2, 0.8,
    3.7,  -0.9, -1.2, 0.7, -1.5, 0.2, 2.2,
    -0.3, -0.9, 2.1,  0.7, -1.6, 0.2, 2.0,
    -1.6, -0.9, -2.1, 0.7, 2.9,  0.2, -3.8,
};
const dt19yc = [_]f64{
    0.5,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    0.5,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    0.5,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    0.5,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    0.5,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    0.7,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    1.7,  0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    -2.6, 0.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    0.5,  -0.9, 0.0,  0.0,  0.0, 0.0, 0.0,
    4.0,  -6.3, 0.0,  0.0,  0.0, 0.0, 0.0,
    -0.5, 0.3,  0.0,  0.0,  0.0, 0.0, 0.0,
    -1.5, 3.0,  0.0,  0.0,  0.0, 0.0, 0.0,
    0.5,  -0.9, 0.3,  0.7,  0.0, 0.0, 0.0,
    3.7,  -7.2, 3.0,  1.7,  0.0, 0.0, 0.0,
    -0.3, 0.9,  -0.7, 1.9,  0.0, 0.0, 0.0,
    -1.6, 2.7,  -0.7, -3.4, 0.0, 0.0, 0.0,
};
const dt19yd = [_]f64{
    0.5,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,
    0.5,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,
    0.5,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,
    0.5,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,
    0.5,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,
    0.7,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,
    1.7,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0,
    -2.6, 0.0,  0.0,  0.0, 0.0,  0.0, 0.0,
    0.5,  -0.9, 0.3,  0.0, 0.0,  0.0, 0.0,
    0.7,  -0.9, 1.2,  0.0, 0.0,  0.0, 0.0,
    1.7,  -0.9, 0.5,  0.0, 0.0,  0.0, 0.0,
    -2.6, -0.9, -1.3, 0.0, 0.0,  0.0, 0.0,
    0.5,  -0.9, 0.3,  0.7, -0.6, 0.2, 0.8,
    0.7,  -0.9, 1.2,  0.7, -1.5, 0.2, 1.6,
    1.7,  -0.9, 0.5,  0.7, -1.6, 0.2, 2.4,
    -2.6, -0.9, -1.3, 0.7, 2.9,  0.2, -4.0,
};
const dt9x = [_]f64{
    0.6,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,
    0.78,  0.0,   0.0,   0.0,  0.0,  0.0,  0.0,
    0.78,  -0.46, 0.0,   0.0,  0.0,  0.0,  0.0,
    0.78,  -0.46, -0.22, 1.06, 0.0,  0.0,  0.0,
    0.6,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,
    0.78,  0.0,   0.0,   0.0,  0.0,  0.0,  0.0,
    0.66,  0.1,   -0.1,  0.0,  0.0,  0.0,  0.0,
    0.96,  0.1,   -0.76, 0.8,  0.9,  -0.3, -0.02,
    0.6,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,
    0.78,  0.0,   0.0,   0.0,  0.0,  0.0,  0.0,
    -0.06, 0.1,   -0.1,  0.0,  0.0,  0.0,  0.0,
    0.9,   0.1,   -0.22, 0.8,  0.18, -0.3, -0.02,
    0.6,   0.0,   0.0,   0.0,  0.0,  0.0,  0.0,
    0.78,  0.0,   0.0,   0.0,  0.0,  0.0,  0.0,
    0.78,  0.26,  0.0,   0.0,  0.0,  0.0,  0.0,
    0.78,  0.26,  -0.76, 1.12, 0.0,  0.0,  0.0,
};
const dt9y = [_]f64{
    0.5,  0.0,   0.0,   0.0,  0.0,   0.0, 0.0,
    0.04, 0.0,   0.0,   0.0,  0.0,   0.0, 0.0,
    0.04, -0.78, 0.0,   0.0,  0.0,   0.0, 0.0,
    0.04, -0.78, 0.54,  0.08, 0.0,   0.0, 0.0,
    0.5,  0.0,   0.0,   0.0,  0.0,   0.0, 0.0,
    0.04, 0.0,   0.0,   0.0,  0.0,   0.0, 0.0,
    0.7,  -0.9,  -0.12, 0.0,  0.0,   0.0, 0.0,
    0.64, -0.9,  -0.3,  0.7,  -0.18, 0.2, 0.28,
    0.5,  0.0,   0.0,   0.0,  0.0,   0.0, 0.0,
    0.04, 0.0,   0.0,   0.0,  0.0,   0.0, 0.0,
    0.7,  -1.08, 0.0,   0.0,  0.0,   0.0, 0.0,
    0.64, -1.26, 0.54,  0.2,  0.0,   0.0, 0.0,
    0.5,  0.0,   0.0,   0.0,  0.0,   0.0, 0.0,
    0.04, 0.0,   0.0,   0.0,  0.0,   0.0, 0.0,
    0.04, -0.9,  0.18,  0.0,  0.0,   0.0, 0.0,
    0.04, -0.9,  0.18,  0.7,  -0.18, 0.2, 0.16,
};
const cv = [_][2]f64{
    .{ 0.1, 0.1 },  .{ 1.0, 2.0 },  .{ 1.0, 2.0 },  .{ 1.0, 2.0 }, .{ 1.0, 2.0 },  .{ 1.0, 2.0 }, .{ 1.0, 2.0 }, .{ 1.0, 2.0 },
    .{ 0.3, -0.4 }, .{ 3.0, 4.0 },  .{ 3.0, 4.0 },  .{ 3.0, 4.0 }, .{ 3.0, 4.0 },  .{ 3.0, 4.0 }, .{ 3.0, 4.0 }, .{ 3.0, 4.0 },
    .{ 0.1, -0.3 }, .{ 0.5, -0.1 }, .{ 5.0, 6.0 },  .{ 5.0, 6.0 }, .{ 5.0, 6.0 },  .{ 5.0, 6.0 }, .{ 5.0, 6.0 }, .{ 5.0, 6.0 },
    .{ 0.1, 0.1 },  .{ -0.6, 0.1 }, .{ 0.1, -0.3 }, .{ 7.0, 8.0 }, .{ 7.0, 8.0 },  .{ 7.0, 8.0 }, .{ 7.0, 8.0 }, .{ 7.0, 8.0 },
    .{ 0.3, 0.1 },  .{ 0.5, 0.0 },  .{ 0.0, 0.5 },  .{ 0.0, 0.2 }, .{ 2.0, 3.0 },  .{ 2.0, 3.0 }, .{ 2.0, 3.0 }, .{ 2.0, 3.0 },
    .{ 0.1, 0.1 },  .{ 4.0, 5.0 },  .{ 4.0, 5.0 },  .{ 4.0, 5.0 }, .{ 4.0, 5.0 },  .{ 4.0, 5.0 }, .{ 4.0, 5.0 }, .{ 4.0, 5.0 },
    .{ 0.3, -0.4 }, .{ 6.0, 7.0 },  .{ 6.0, 7.0 },  .{ 6.0, 7.0 }, .{ 6.0, 7.0 },  .{ 6.0, 7.0 }, .{ 6.0, 7.0 }, .{ 6.0, 7.0 },
    .{ 0.1, -0.3 }, .{ 8.0, 9.0 },  .{ 0.5, -0.1 }, .{ 2.0, 5.0 }, .{ 2.0, 5.0 },  .{ 2.0, 5.0 }, .{ 2.0, 5.0 }, .{ 2.0, 5.0 },
    .{ 0.1, 0.1 },  .{ 3.0, 6.0 },  .{ -0.6, 0.1 }, .{ 4.0, 7.0 }, .{ 0.1, -0.3 }, .{ 7.0, 2.0 }, .{ 7.0, 2.0 }, .{ 7.0, 2.0 },
    .{ 0.3, 0.1 },  .{ 5.0, 8.0 },  .{ 0.5, 0.0 },  .{ 6.0, 9.0 }, .{ 0.0, 0.5 },  .{ 8.0, 3.0 }, .{ 0.0, 0.2 }, .{ 9.0, 4.0 },
};
const strue2 = [_]f64{
    0.0, 0.5, 0.6, 0.7, 0.8,
};
const strue4 = [_]f64{
    0.0, 0.7, 1.0, 1.3, 1.6,
};
const ctrue5 = [_][2]f64{
    .{ 0.1, 0.1 },     .{ 1.0, 2.0 },    .{ 1.0, 2.0 },     .{ 1.0, 2.0 },   .{ 1.0, 2.0 },     .{ 1.0, 2.0 }, .{ 1.0, 2.0 },   .{ 1.0, 2.0 },
    .{ -0.16, -0.37 }, .{ 3.0, 4.0 },    .{ 3.0, 4.0 },     .{ 3.0, 4.0 },   .{ 3.0, 4.0 },     .{ 3.0, 4.0 }, .{ 3.0, 4.0 },   .{ 3.0, 4.0 },
    .{ -0.17, -0.19 }, .{ 0.13, -0.39 }, .{ 5.0, 6.0 },     .{ 5.0, 6.0 },   .{ 5.0, 6.0 },     .{ 5.0, 6.0 }, .{ 5.0, 6.0 },   .{ 5.0, 6.0 },
    .{ 0.11, -0.03 },  .{ -0.17, 0.46 }, .{ -0.17, -0.19 }, .{ 7.0, 8.0 },   .{ 7.0, 8.0 },     .{ 7.0, 8.0 }, .{ 7.0, 8.0 },   .{ 7.0, 8.0 },
    .{ 0.19, -0.17 },  .{ 0.2, -0.35 },  .{ 0.35, 0.2 },    .{ 0.14, 0.08 }, .{ 2.0, 3.0 },     .{ 2.0, 3.0 }, .{ 2.0, 3.0 },   .{ 2.0, 3.0 },
    .{ 0.1, 0.1 },     .{ 4.0, 5.0 },    .{ 4.0, 5.0 },     .{ 4.0, 5.0 },   .{ 4.0, 5.0 },     .{ 4.0, 5.0 }, .{ 4.0, 5.0 },   .{ 4.0, 5.0 },
    .{ -0.16, -0.37 }, .{ 6.0, 7.0 },    .{ 6.0, 7.0 },     .{ 6.0, 7.0 },   .{ 6.0, 7.0 },     .{ 6.0, 7.0 }, .{ 6.0, 7.0 },   .{ 6.0, 7.0 },
    .{ -0.17, -0.19 }, .{ 8.0, 9.0 },    .{ 0.13, -0.39 },  .{ 2.0, 5.0 },   .{ 2.0, 5.0 },     .{ 2.0, 5.0 }, .{ 2.0, 5.0 },   .{ 2.0, 5.0 },
    .{ 0.11, -0.03 },  .{ 3.0, 6.0 },    .{ -0.17, 0.46 },  .{ 4.0, 7.0 },   .{ -0.17, -0.19 }, .{ 7.0, 2.0 }, .{ 7.0, 2.0 },   .{ 7.0, 2.0 },
    .{ 0.19, -0.17 },  .{ 5.0, 8.0 },    .{ 0.2, -0.35 },   .{ 6.0, 9.0 },   .{ 0.35, 0.2 },    .{ 8.0, 3.0 }, .{ 0.14, 0.08 }, .{ 9.0, 4.0 },
};
const ctrue6 = [_][2]f64{
    .{ 0.1, 0.1 },    .{ 1.0, 2.0 },    .{ 1.0, 2.0 },    .{ 1.0, 2.0 },  .{ 1.0, 2.0 },    .{ 1.0, 2.0 }, .{ 1.0, 2.0 },  .{ 1.0, 2.0 },
    .{ 0.09, -0.12 }, .{ 3.0, 4.0 },    .{ 3.0, 4.0 },    .{ 3.0, 4.0 },  .{ 3.0, 4.0 },    .{ 3.0, 4.0 }, .{ 3.0, 4.0 },  .{ 3.0, 4.0 },
    .{ 0.03, -0.09 }, .{ 0.15, -0.03 }, .{ 5.0, 6.0 },    .{ 5.0, 6.0 },  .{ 5.0, 6.0 },    .{ 5.0, 6.0 }, .{ 5.0, 6.0 },  .{ 5.0, 6.0 },
    .{ 0.03, 0.03 },  .{ -0.18, 0.03 }, .{ 0.03, -0.09 }, .{ 7.0, 8.0 },  .{ 7.0, 8.0 },    .{ 7.0, 8.0 }, .{ 7.0, 8.0 },  .{ 7.0, 8.0 },
    .{ 0.09, 0.03 },  .{ 0.15, 0.0 },   .{ 0.0, 0.15 },   .{ 0.0, 0.06 }, .{ 2.0, 3.0 },    .{ 2.0, 3.0 }, .{ 2.0, 3.0 },  .{ 2.0, 3.0 },
    .{ 0.1, 0.1 },    .{ 4.0, 5.0 },    .{ 4.0, 5.0 },    .{ 4.0, 5.0 },  .{ 4.0, 5.0 },    .{ 4.0, 5.0 }, .{ 4.0, 5.0 },  .{ 4.0, 5.0 },
    .{ 0.09, -0.12 }, .{ 6.0, 7.0 },    .{ 6.0, 7.0 },    .{ 6.0, 7.0 },  .{ 6.0, 7.0 },    .{ 6.0, 7.0 }, .{ 6.0, 7.0 },  .{ 6.0, 7.0 },
    .{ 0.03, -0.09 }, .{ 8.0, 9.0 },    .{ 0.15, -0.03 }, .{ 2.0, 5.0 },  .{ 2.0, 5.0 },    .{ 2.0, 5.0 }, .{ 2.0, 5.0 },  .{ 2.0, 5.0 },
    .{ 0.03, 0.03 },  .{ 3.0, 6.0 },    .{ -0.18, 0.03 }, .{ 4.0, 7.0 },  .{ 0.03, -0.09 }, .{ 7.0, 2.0 }, .{ 7.0, 2.0 },  .{ 7.0, 2.0 },
    .{ 0.09, 0.03 },  .{ 5.0, 8.0 },    .{ 0.15, 0.0 },   .{ 6.0, 9.0 },  .{ 0.0, 0.15 },   .{ 8.0, 3.0 }, .{ 0.0, 0.06 }, .{ 9.0, 4.0 },
};
const itrue3 = [_]isize{
    0, 1, 2, 2, 2,
};
const cx1 = [_][2]f64{
    .{ 0.7, -0.8 }, .{ -0.4, -0.7 }, .{ -0.1, -0.9 }, .{ 0.2, -0.8 }, .{ -0.9, -0.4 }, .{ 0.1, 0.4 }, .{ -0.6, 0.6 },
};
const cy1 = [_][2]f64{
    .{ 0.6, -0.6 }, .{ -0.9, 0.5 }, .{ 0.7, -0.6 }, .{ 0.1, -0.5 }, .{ -0.1, -0.2 }, .{ -0.5, -0.3 }, .{ 0.8, -0.7 },
};
const ct8 = [_][2]f64{
    .{ 0.6, -0.6 },    .{ 0.0, 0.0 },     .{ 0.0, 0.0 },    .{ 0.0, 0.0 },     .{ 0.0, 0.0 },     .{ 0.0, 0.0 },   .{ 0.0, 0.0 },
    .{ 0.32, -1.41 },  .{ 0.0, 0.0 },     .{ 0.0, 0.0 },    .{ 0.0, 0.0 },     .{ 0.0, 0.0 },     .{ 0.0, 0.0 },   .{ 0.0, 0.0 },
    .{ 0.32, -1.41 },  .{ -1.55, 0.5 },   .{ 0.0, 0.0 },    .{ 0.0, 0.0 },     .{ 0.0, 0.0 },     .{ 0.0, 0.0 },   .{ 0.0, 0.0 },
    .{ 0.32, -1.41 },  .{ -1.55, 0.5 },   .{ 0.03, -0.89 }, .{ -0.38, -0.96 }, .{ 0.0, 0.0 },     .{ 0.0, 0.0 },   .{ 0.0, 0.0 },
    .{ 0.6, -0.6 },    .{ 0.0, 0.0 },     .{ 0.0, 0.0 },    .{ 0.0, 0.0 },     .{ 0.0, 0.0 },     .{ 0.0, 0.0 },   .{ 0.0, 0.0 },
    .{ 0.32, -1.41 },  .{ 0.0, 0.0 },     .{ 0.0, 0.0 },    .{ 0.0, 0.0 },     .{ 0.0, 0.0 },     .{ 0.0, 0.0 },   .{ 0.0, 0.0 },
    .{ -0.07, -0.89 }, .{ -0.9, 0.5 },    .{ 0.42, -1.41 }, .{ 0.0, 0.0 },     .{ 0.0, 0.0 },     .{ 0.0, 0.0 },   .{ 0.0, 0.0 },
    .{ 0.78, 0.06 },   .{ -0.9, 0.5 },    .{ 0.06, -0.13 }, .{ 0.1, -0.5 },    .{ -0.77, -0.49 }, .{ -0.5, -0.3 }, .{ 0.52, -1.51 },
    .{ 0.6, -0.6 },    .{ 0.0, 0.0 },     .{ 0.0, 0.0 },    .{ 0.0, 0.0 },     .{ 0.0, 0.0 },     .{ 0.0, 0.0 },   .{ 0.0, 0.0 },
    .{ 0.32, -1.41 },  .{ 0.0, 0.0 },     .{ 0.0, 0.0 },    .{ 0.0, 0.0 },     .{ 0.0, 0.0 },     .{ 0.0, 0.0 },   .{ 0.0, 0.0 },
    .{ -0.07, -0.89 }, .{ -1.18, -0.31 }, .{ 0.0, 0.0 },    .{ 0.0, 0.0 },     .{ 0.0, 0.0 },     .{ 0.0, 0.0 },   .{ 0.0, 0.0 },
    .{ 0.78, 0.06 },   .{ -1.54, 0.97 },  .{ 0.03, -0.89 }, .{ -0.18, -1.31 }, .{ 0.0, 0.0 },     .{ 0.0, 0.0 },   .{ 0.0, 0.0 },
    .{ 0.6, -0.6 },    .{ 0.0, 0.0 },     .{ 0.0, 0.0 },    .{ 0.0, 0.0 },     .{ 0.0, 0.0 },     .{ 0.0, 0.0 },   .{ 0.0, 0.0 },
    .{ 0.32, -1.41 },  .{ 0.0, 0.0 },     .{ 0.0, 0.0 },    .{ 0.0, 0.0 },     .{ 0.0, 0.0 },     .{ 0.0, 0.0 },   .{ 0.0, 0.0 },
    .{ 0.32, -1.41 },  .{ -0.9, 0.5 },    .{ 0.05, -0.6 },  .{ 0.0, 0.0 },     .{ 0.0, 0.0 },     .{ 0.0, 0.0 },   .{ 0.0, 0.0 },
    .{ 0.32, -1.41 },  .{ -0.9, 0.5 },    .{ 0.05, -0.6 },  .{ 0.1, -0.5 },    .{ -0.77, -0.49 }, .{ -0.5, -0.3 }, .{ 0.32, -1.16 },
};
const ct7 = [_][2]f64{
    .{ 0.0, 0.0 }, .{ -0.06, -0.9 }, .{ 0.65, -0.47 },  .{ -0.34, -1.22 },
    .{ 0.0, 0.0 }, .{ -0.06, -0.9 }, .{ -0.59, -1.46 }, .{ -1.04, -0.04 },
    .{ 0.0, 0.0 }, .{ -0.06, -0.9 }, .{ -0.83, 0.59 },  .{ 0.07, -0.37 },
    .{ 0.0, 0.0 }, .{ -0.06, -0.9 }, .{ -0.76, -1.15 }, .{ -1.33, -1.82 },
};
const ct6 = [_][2]f64{
    .{ 0.0, 0.0 }, .{ 0.9, 0.06 }, .{ 0.91, -0.77 }, .{ 1.8, -0.1 },
    .{ 0.0, 0.0 }, .{ 0.9, 0.06 }, .{ 1.45, 0.74 },  .{ 0.2, 0.9 },
    .{ 0.0, 0.0 }, .{ 0.9, 0.06 }, .{ -0.55, 0.23 }, .{ 0.83, -0.39 },
    .{ 0.0, 0.0 }, .{ 0.9, 0.06 }, .{ 1.04, 0.79 },  .{ 1.95, 1.22 },
};
const ct10x = [_][2]f64{
    .{ 0.7, -0.8 }, .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },  .{ 0.0, 0.0 },  .{ 0.0, 0.0 }, .{ 0.0, 0.0 },
    .{ 0.6, -0.6 }, .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },  .{ 0.0, 0.0 },  .{ 0.0, 0.0 }, .{ 0.0, 0.0 },
    .{ 0.6, -0.6 }, .{ -0.9, 0.5 },  .{ 0.0, 0.0 },   .{ 0.0, 0.0 },  .{ 0.0, 0.0 },  .{ 0.0, 0.0 }, .{ 0.0, 0.0 },
    .{ 0.6, -0.6 }, .{ -0.9, 0.5 },  .{ 0.7, -0.6 },  .{ 0.1, -0.5 }, .{ 0.0, 0.0 },  .{ 0.0, 0.0 }, .{ 0.0, 0.0 },
    .{ 0.7, -0.8 }, .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },  .{ 0.0, 0.0 },  .{ 0.0, 0.0 }, .{ 0.0, 0.0 },
    .{ 0.6, -0.6 }, .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },  .{ 0.0, 0.0 },  .{ 0.0, 0.0 }, .{ 0.0, 0.0 },
    .{ 0.7, -0.6 }, .{ -0.4, -0.7 }, .{ 0.6, -0.6 },  .{ 0.0, 0.0 },  .{ 0.0, 0.0 },  .{ 0.0, 0.0 }, .{ 0.0, 0.0 },
    .{ 0.8, -0.7 }, .{ -0.4, -0.7 }, .{ -0.1, -0.2 }, .{ 0.2, -0.8 }, .{ 0.7, -0.6 }, .{ 0.1, 0.4 }, .{ 0.6, -0.6 },
    .{ 0.7, -0.8 }, .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },  .{ 0.0, 0.0 },  .{ 0.0, 0.0 }, .{ 0.0, 0.0 },
    .{ 0.6, -0.6 }, .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },  .{ 0.0, 0.0 },  .{ 0.0, 0.0 }, .{ 0.0, 0.0 },
    .{ -0.9, 0.5 }, .{ -0.4, -0.7 }, .{ 0.6, -0.6 },  .{ 0.0, 0.0 },  .{ 0.0, 0.0 },  .{ 0.0, 0.0 }, .{ 0.0, 0.0 },
    .{ 0.1, -0.5 }, .{ -0.4, -0.7 }, .{ 0.7, -0.6 },  .{ 0.2, -0.8 }, .{ -0.9, 0.5 }, .{ 0.1, 0.4 }, .{ 0.6, -0.6 },
    .{ 0.7, -0.8 }, .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },  .{ 0.0, 0.0 },  .{ 0.0, 0.0 }, .{ 0.0, 0.0 },
    .{ 0.6, -0.6 }, .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },  .{ 0.0, 0.0 },  .{ 0.0, 0.0 }, .{ 0.0, 0.0 },
    .{ 0.6, -0.6 }, .{ 0.7, -0.6 },  .{ 0.0, 0.0 },   .{ 0.0, 0.0 },  .{ 0.0, 0.0 },  .{ 0.0, 0.0 }, .{ 0.0, 0.0 },
    .{ 0.6, -0.6 }, .{ 0.7, -0.6 },  .{ -0.1, -0.2 }, .{ 0.8, -0.7 }, .{ 0.0, 0.0 },  .{ 0.0, 0.0 }, .{ 0.0, 0.0 },
};
const ct10y = [_][2]f64{
    .{ 0.6, -0.6 },  .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },  .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },
    .{ 0.7, -0.8 },  .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },  .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },
    .{ 0.7, -0.8 },  .{ -0.4, -0.7 }, .{ 0.0, 0.0 },   .{ 0.0, 0.0 },  .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },
    .{ 0.7, -0.8 },  .{ -0.4, -0.7 }, .{ -0.1, -0.9 }, .{ 0.2, -0.8 }, .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },
    .{ 0.6, -0.6 },  .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },  .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },
    .{ 0.7, -0.8 },  .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },  .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },
    .{ -0.1, -0.9 }, .{ -0.9, 0.5 },  .{ 0.7, -0.8 },  .{ 0.0, 0.0 },  .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },
    .{ -0.6, 0.6 },  .{ -0.9, 0.5 },  .{ -0.9, -0.4 }, .{ 0.1, -0.5 }, .{ -0.1, -0.9 }, .{ -0.5, -0.3 }, .{ 0.7, -0.8 },
    .{ 0.6, -0.6 },  .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },  .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },
    .{ 0.7, -0.8 },  .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },  .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },
    .{ -0.1, -0.9 }, .{ 0.7, -0.8 },  .{ 0.0, 0.0 },   .{ 0.0, 0.0 },  .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },
    .{ -0.6, 0.6 },  .{ -0.9, -0.4 }, .{ -0.1, -0.9 }, .{ 0.7, -0.8 }, .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },
    .{ 0.6, -0.6 },  .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },  .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },
    .{ 0.7, -0.8 },  .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },  .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },
    .{ 0.7, -0.8 },  .{ -0.9, 0.5 },  .{ -0.4, -0.7 }, .{ 0.0, 0.0 },  .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },
    .{ 0.7, -0.8 },  .{ -0.9, 0.5 },  .{ -0.4, -0.7 }, .{ 0.1, -0.5 }, .{ -0.1, -0.9 }, .{ -0.5, -0.3 }, .{ 0.2, -0.8 },
};
const csize1 = [_][2]f64{
    .{ 0.0, 0.0 }, .{ 0.9, 0.9 }, .{ 1.63, 1.73 }, .{ 2.9, 2.78 },
};
const csize3 = [_][2]f64{
    .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },
    .{ 1.17, 1.17 }, .{ 1.17, 1.17 }, .{ 1.17, 1.17 }, .{ 1.17, 1.17 }, .{ 1.17, 1.17 }, .{ 1.17, 1.17 }, .{ 1.17, 1.17 },
};
const csize2 = [_][2]f64{
    .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },   .{ 0.0, 0.0 },
    .{ 1.54, 1.54 }, .{ 1.54, 1.54 }, .{ 1.54, 1.54 }, .{ 1.54, 1.54 }, .{ 1.54, 1.54 }, .{ 1.54, 1.54 }, .{ 1.54, 1.54 },
};
const c_cv = [_][2]f64{
    .{ 0.1, 0.1 },  .{ 1.0, 2.0 },  .{ 1.0, 2.0 },  .{ 1.0, 2.0 }, .{ 1.0, 2.0 },  .{ 1.0, 2.0 }, .{ 1.0, 2.0 }, .{ 1.0, 2.0 },
    .{ 0.3, -0.4 }, .{ 3.0, 4.0 },  .{ 3.0, 4.0 },  .{ 3.0, 4.0 }, .{ 3.0, 4.0 },  .{ 3.0, 4.0 }, .{ 3.0, 4.0 }, .{ 3.0, 4.0 },
    .{ 0.1, -0.3 }, .{ 0.5, -0.1 }, .{ 5.0, 6.0 },  .{ 5.0, 6.0 }, .{ 5.0, 6.0 },  .{ 5.0, 6.0 }, .{ 5.0, 6.0 }, .{ 5.0, 6.0 },
    .{ 0.1, 0.1 },  .{ -0.6, 0.1 }, .{ 0.1, -0.3 }, .{ 7.0, 8.0 }, .{ 7.0, 8.0 },  .{ 7.0, 8.0 }, .{ 7.0, 8.0 }, .{ 7.0, 8.0 },
    .{ 0.3, 0.1 },  .{ 0.1, 0.4 },  .{ 0.4, 0.1 },  .{ 0.1, 0.2 }, .{ 2.0, 3.0 },  .{ 2.0, 3.0 }, .{ 2.0, 3.0 }, .{ 2.0, 3.0 },
    .{ 0.1, 0.1 },  .{ 4.0, 5.0 },  .{ 4.0, 5.0 },  .{ 4.0, 5.0 }, .{ 4.0, 5.0 },  .{ 4.0, 5.0 }, .{ 4.0, 5.0 }, .{ 4.0, 5.0 },
    .{ 0.3, -0.4 }, .{ 6.0, 7.0 },  .{ 6.0, 7.0 },  .{ 6.0, 7.0 }, .{ 6.0, 7.0 },  .{ 6.0, 7.0 }, .{ 6.0, 7.0 }, .{ 6.0, 7.0 },
    .{ 0.1, -0.3 }, .{ 8.0, 9.0 },  .{ 0.5, -0.1 }, .{ 2.0, 5.0 }, .{ 2.0, 5.0 },  .{ 2.0, 5.0 }, .{ 2.0, 5.0 }, .{ 2.0, 5.0 },
    .{ 0.1, 0.1 },  .{ 3.0, 6.0 },  .{ -0.6, 0.1 }, .{ 4.0, 7.0 }, .{ 0.1, -0.3 }, .{ 7.0, 2.0 }, .{ 7.0, 2.0 }, .{ 7.0, 2.0 },
    .{ 0.3, 0.1 },  .{ 5.0, 8.0 },  .{ 0.1, 0.4 },  .{ 6.0, 9.0 }, .{ 0.4, 0.1 },  .{ 8.0, 3.0 }, .{ 0.1, 0.2 }, .{ 9.0, 4.0 },
};
const c_strue2 = [_]f64{
    0.0, 0.5, 0.6, 0.7, 0.7,
};
const c_strue4 = [_]f64{
    0.0, 0.7, 1.0, 1.3, 1.7,
};
const c_ctrue5 = [_][2]f64{
    .{ 0.1, 0.1 },     .{ 1.0, 2.0 },    .{ 1.0, 2.0 },     .{ 1.0, 2.0 },   .{ 1.0, 2.0 },     .{ 1.0, 2.0 }, .{ 1.0, 2.0 },   .{ 1.0, 2.0 },
    .{ -0.16, -0.37 }, .{ 3.0, 4.0 },    .{ 3.0, 4.0 },     .{ 3.0, 4.0 },   .{ 3.0, 4.0 },     .{ 3.0, 4.0 }, .{ 3.0, 4.0 },   .{ 3.0, 4.0 },
    .{ -0.17, -0.19 }, .{ 0.13, -0.39 }, .{ 5.0, 6.0 },     .{ 5.0, 6.0 },   .{ 5.0, 6.0 },     .{ 5.0, 6.0 }, .{ 5.0, 6.0 },   .{ 5.0, 6.0 },
    .{ 0.11, -0.03 },  .{ -0.17, 0.46 }, .{ -0.17, -0.19 }, .{ 7.0, 8.0 },   .{ 7.0, 8.0 },     .{ 7.0, 8.0 }, .{ 7.0, 8.0 },   .{ 7.0, 8.0 },
    .{ 0.19, -0.17 },  .{ 0.32, 0.09 },  .{ 0.23, -0.24 },  .{ 0.18, 0.01 }, .{ 2.0, 3.0 },     .{ 2.0, 3.0 }, .{ 2.0, 3.0 },   .{ 2.0, 3.0 },
    .{ 0.1, 0.1 },     .{ 4.0, 5.0 },    .{ 4.0, 5.0 },     .{ 4.0, 5.0 },   .{ 4.0, 5.0 },     .{ 4.0, 5.0 }, .{ 4.0, 5.0 },   .{ 4.0, 5.0 },
    .{ -0.16, -0.37 }, .{ 6.0, 7.0 },    .{ 6.0, 7.0 },     .{ 6.0, 7.0 },   .{ 6.0, 7.0 },     .{ 6.0, 7.0 }, .{ 6.0, 7.0 },   .{ 6.0, 7.0 },
    .{ -0.17, -0.19 }, .{ 8.0, 9.0 },    .{ 0.13, -0.39 },  .{ 2.0, 5.0 },   .{ 2.0, 5.0 },     .{ 2.0, 5.0 }, .{ 2.0, 5.0 },   .{ 2.0, 5.0 },
    .{ 0.11, -0.03 },  .{ 3.0, 6.0 },    .{ -0.17, 0.46 },  .{ 4.0, 7.0 },   .{ -0.17, -0.19 }, .{ 7.0, 2.0 }, .{ 7.0, 2.0 },   .{ 7.0, 2.0 },
    .{ 0.19, -0.17 },  .{ 5.0, 8.0 },    .{ 0.32, 0.09 },   .{ 6.0, 9.0 },   .{ 0.23, -0.24 },  .{ 8.0, 3.0 }, .{ 0.18, 0.01 }, .{ 9.0, 4.0 },
};
const c_ctrue6 = [_][2]f64{
    .{ 0.1, 0.1 },    .{ 1.0, 2.0 },    .{ 1.0, 2.0 },    .{ 1.0, 2.0 },   .{ 1.0, 2.0 },    .{ 1.0, 2.0 }, .{ 1.0, 2.0 },   .{ 1.0, 2.0 },
    .{ 0.09, -0.12 }, .{ 3.0, 4.0 },    .{ 3.0, 4.0 },    .{ 3.0, 4.0 },   .{ 3.0, 4.0 },    .{ 3.0, 4.0 }, .{ 3.0, 4.0 },   .{ 3.0, 4.0 },
    .{ 0.03, -0.09 }, .{ 0.15, -0.03 }, .{ 5.0, 6.0 },    .{ 5.0, 6.0 },   .{ 5.0, 6.0 },    .{ 5.0, 6.0 }, .{ 5.0, 6.0 },   .{ 5.0, 6.0 },
    .{ 0.03, 0.03 },  .{ -0.18, 0.03 }, .{ 0.03, -0.09 }, .{ 7.0, 8.0 },   .{ 7.0, 8.0 },    .{ 7.0, 8.0 }, .{ 7.0, 8.0 },   .{ 7.0, 8.0 },
    .{ 0.09, 0.03 },  .{ 0.03, 0.12 },  .{ 0.12, 0.03 },  .{ 0.03, 0.06 }, .{ 2.0, 3.0 },    .{ 2.0, 3.0 }, .{ 2.0, 3.0 },   .{ 2.0, 3.0 },
    .{ 0.1, 0.1 },    .{ 4.0, 5.0 },    .{ 4.0, 5.0 },    .{ 4.0, 5.0 },   .{ 4.0, 5.0 },    .{ 4.0, 5.0 }, .{ 4.0, 5.0 },   .{ 4.0, 5.0 },
    .{ 0.09, -0.12 }, .{ 6.0, 7.0 },    .{ 6.0, 7.0 },    .{ 6.0, 7.0 },   .{ 6.0, 7.0 },    .{ 6.0, 7.0 }, .{ 6.0, 7.0 },   .{ 6.0, 7.0 },
    .{ 0.03, -0.09 }, .{ 8.0, 9.0 },    .{ 0.15, -0.03 }, .{ 2.0, 5.0 },   .{ 2.0, 5.0 },    .{ 2.0, 5.0 }, .{ 2.0, 5.0 },   .{ 2.0, 5.0 },
    .{ 0.03, 0.03 },  .{ 3.0, 6.0 },    .{ -0.18, 0.03 }, .{ 4.0, 7.0 },   .{ 0.03, -0.09 }, .{ 7.0, 2.0 }, .{ 7.0, 2.0 },   .{ 7.0, 2.0 },
    .{ 0.09, 0.03 },  .{ 5.0, 8.0 },    .{ 0.03, 0.12 },  .{ 6.0, 9.0 },   .{ 0.12, 0.03 },  .{ 8.0, 3.0 }, .{ 0.03, 0.06 }, .{ 9.0, 4.0 },
};
const c_sdx = [_]f64{
    0.6, 0.1, -0.5, 0.8, 0.9, -0.3, -0.4, 0.5, -0.9, 0.3, 0.7, -0.6, 0.2, 0.8, -0.46, 0.78, -0.46, -0.22, 1.06,
};
const c_sdy = [_]f64{
    0.5, -0.9, 0.3, 0.7, -0.6, 0.2, 0.6, 0.1, -0.5, 0.8, 0.9, -0.3, 0.96, 0.1, -0.76, 0.8, 0.9, 0.66, 0.8,
};
const c_sns = [_]isize{
    0, 1, 2, 4, 5, 8, 9,
};
const c_sssize2 = [_]f64{
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17,
};
const c_param = [_]f64{
    -2.0, 1.0, 0.0,  0.0, 1.0,
    -1.0, 0.2, 0.3,  0.4, 0.5,
    0.0,  1.0, 0.3,  0.4, 1.0,
    1.0,  0.2, -1.0, 1.0, 0.5,
};
const c_incxs = [_]isize{
    1, 1, 2, 2, -2, -1, -2,
};
const c_incys = [_]isize{
    1, 2, 2, -2, 1, -2, -2,
};
const c_ddx = [_]f64{
    0.6, 0.1, -0.5, 0.8, 0.9, -0.3, -0.4, 0.7, 0.5, 0.2,
};
const c_ddy = [_]f64{
    0.5, -0.9, 0.3, 0.7, -0.6, 0.2, 0.8, -0.5, 0.1, -0.3,
};
const c_dns = [_]isize{
    0, 1, 2, 4, 5,
};
const c_dssize2 = [_]f64{
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17,
};
