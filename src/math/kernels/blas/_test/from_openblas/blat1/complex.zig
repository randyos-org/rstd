//! `cblat1.f` and `zblat1.f`'s checks: `CHECK1` and `CHECK2` on Fortran's data and `ctest`'s.

const std = @import("std");
const ob = @import("../__root__.zig");
const l1 = @import("../../../blas1/interface.zig");
const cplx = @import("complex_data.zig");
const b1 = @import("__root__.zig");

const Float = ob.Float;
const Case = b1.Case;
const ctest = b1.ctest;
const incxs = b1.incxs;
const incys = b1.incys;
const itest1 = b1.itest1;
const lens = b1.lens;
const loadC = b1.loadC;
const ns = b1.ns;
const sfac = b1.sfac;
const stest = b1.stest;

// ---- complex -----------------------------------------------------------------

/// `CHECK1`: `nrm2`, `asum`, `scal`, `scalReal`, `iamax`, from the data
/// `d` (the Fortran's or ctest's), compared by `old` as `STEST` says.
pub fn ccheck1(comptime C: type, d: anytype, old: bool) !void {
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
            try stest(F, c, &.{try l1.nrm2(ob.io, C, in, &cx, inc)}, d.strue2[np1..], d.strue2[np1..], sfac, old);
            try stest(F, c, &.{try l1.asum(ob.io, C, in, &cx, inc)}, d.strue4[np1..], d.strue4[np1..], sfac, old);
            try itest1(c, try l1.index(ob.io, C, .largest_magnitude, in, &cx, inc), cplx.itrue3[np1]);
            if (comptime std.mem.eql(u8, which, "scal")) {
                try l1.scal(ob.io, C, in, ca, &cx, inc);
                try ctest(C, c, cx[0..len], d.ctrue5[col..], d.ctrue5[col..], sfac, old);
            } else {
                try l1.scalReal(ob.io, C, in, sa, &cx, inc);
                try ctest(C, c, cx[0..len], d.ctrue6[col..], d.ctrue6[col..], sfac, old);
            }
        };
        // then on the last `x`, with `incx` 1: `alpha` 0, and for the real
        // `alpha` 1 and -1
        const c: Case = .{ .name = p ++ "SCAL", .n = 5, .incx = 1 };
        const zeros: [5][2]f64 = @splat(.{ 0, 0 });
        const ones: [5][2]f64 = @splat(.{ 1, 1 });
        if (comptime std.mem.eql(u8, which, "scal")) {
            try l1.scal(ob.io, C, 5, .init(0, 0), &cx, 1);
            try ctest(C, c, cx[0..5], &zeros, &ones, sfac, old);
        } else {
            try l1.scalReal(ob.io, C, 5, 0, &cx, 1);
            try ctest(C, c, cx[0..5], &zeros, &ones, sfac, old);
            for ([_]F{ 1, -1 }) |s| {
                var want: [5][2]f64 = undefined;
                for (&want, cx[0..5]) |*w, z| w.* = .{ s * z.re, s * z.im };
                try l1.scalReal(ob.io, C, 5, s, &cx, 1);
                try ctest(C, c, cx[0..5], &want, &want, sfac, old);
            }
        }
    }
}

/// `CHECK2`: `dotc`, `dotu`, `axpy`, `copy`, `swap`.
pub fn ccheck2(comptime C: type) !void {
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
            loadC(C, &cx, &cplx.cx1);
            loadC(C, &cy, &cplx.cy1);
            try ctest(C, c, &.{try l1.dotc(ob.io, C, n, &cx, incx, &cy, incy)}, cplx.ct6[kn + 4 * ki ..], cplx.csize1[kn..], sfac, false);
            try ctest(C, c, &.{try l1.dot(ob.io, C, n, &cx, incx, &cy, incy)}, cplx.ct7[kn + 4 * ki ..], cplx.csize1[kn..], sfac, false);
            try l1.axpy(ob.io, C, n, ca, &cx, incx, &cy, incy);
            try ctest(C, c, cy[0..leny], cplx.ct8[at..], cplx.csize2[7 * (ksize - 1) ..], sfac, false);
            loadC(C, &cy, &cplx.cy1);
            try l1.copy(ob.io, C, n, &cx, incx, &cy, incy);
            try ctest(C, c, cy[0..leny], cplx.ct10y[at..], &cplx.csize3, 1, false);
            loadC(C, &cy, &cplx.cy1);
            try l1.swap(ob.io, C, n, &cx, incx, &cy, incy);
            try ctest(C, c, cx[0..lenx], cplx.ct10x[at..], &cplx.csize3, 1, false);
            try ctest(C, c, cy[0..leny], cplx.ct10y[at..], &cplx.csize3, 1, false);
        }
    }
}

pub const fortran_c1 = .{ .cv = &cplx.cv, .strue2 = &cplx.strue2, .strue4 = &cplx.strue4, .ctrue5 = &cplx.ctrue5, .ctrue6 = &cplx.ctrue6 };
pub const ctest_c1 = .{ .cv = &cplx.c_cv, .strue2 = &cplx.c_strue2, .strue4 = &cplx.c_strue4, .ctrue5 = &cplx.c_ctrue5, .ctrue6 = &cplx.c_ctrue6 };
