//! Matrix norms and condition estimates, from Reference-LAPACK 3.12.0 as
//! OpenBLAS v0.3.34-317 ships it (`lapack-netlib/SRC/`, modified BSD):
//! `xLANGE` with `xLASSQ` here, `xLACN2` (`lacn2.zig`), `xLATRS` and `xRSCL`
//! (`latrs.zig`) and `xGECON` (`gecon.zig`), real and complex in one loop
//! (the complex branches say where `zlatrs.f` departs).

const std = @import("std");
const blas = @import("../blas/__root__.zig");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const strided = @import("../blas/strided.zig");
const interface = @import("interface.zig");
const ops = @import("ops.zig");
const lacn2_estimate = @import("lacn2.zig");
const latrs_solve = @import("latrs.zig");
const gecon_estimate = @import("gecon.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const isComplex = types.isComplex;
const ConstStrided = strided.ConstStrided;

pub const lacn2 = lacn2_estimate.lacn2;
pub const rscl = latrs_solve.rscl;
pub const latrs = latrs_solve.latrs;
pub const LatrsError = latrs_solve.LatrsError;
pub const scaledSolve = latrs_solve.scaledSolve;
pub const GeconError = gecon_estimate.GeconError;
pub const gecon = gecon_estimate.gecon;

/// `NORM`: the largest `|a(i, j)|` (`m`), the one norm (`o`, LAPACK's `1`
/// or `O`), the infinity norm (`i`) or the Frobenius norm (`f`, `F` or
/// `E`). `|.|` is the modulus for complex.
pub const Norm = enum { m, o, i, f };

// ---- xLANGE, xLASSQ ----------------------------------------------------------------

/// What `lange` returns: a slice too short or a leading dimension too
/// small, and the allocator's.
pub const LangeError = interface.ShapeError || Allocator.Error;

/// `?LANGE` (`dlange.f`, `zlange.f`): `norm` of the `m x n` `A`; 0 when
/// empty. A NaN entry gives NaN. The reference checks no argument; here a
/// slice too short or a short `lda` is an error. `i` takes `m` reals of
/// workspace.
pub fn lange(allocator: Allocator, io: Io, comptime T: type, norm: Norm, m: usize, n: usize, a: []const T, lda: usize) LangeError!Float(T) {
    comptime ops.requireFloat(T, "lange");
    _ = io;
    try util.holds(a.len, m, n, lda);
    if (m == 0 or n == 0) return 0;
    return matNorm(allocator, T, norm, m, n, a.ptr, lda);
}

/// `lange` after its checks.
pub fn matNorm(al: Allocator, comptime T: type, norm: Norm, m: usize, n: usize, a: [*]const T, lda: usize) Allocator.Error!Float(T) {
    const F = Float(T);
    if (m == 0 or n == 0) return 0;
    var value: F = 0;
    switch (norm) {
        .m => for (0..n) |j| for (0..m) |i| {
            const t = ops.abs(T, a[i + j * lda]);
            if (value < t or std.math.isNan(t)) value = t;
        },
        .o => for (0..n) |j| {
            var sum: F = 0;
            for (0..m) |i| sum += ops.abs(T, a[i + j * lda]);
            if (value < sum or std.math.isNan(sum)) value = sum;
        },
        .i => {
            const work = try al.alloc(F, m);
            defer al.free(work);
            @memset(work, 0);
            for (0..n) |j| for (0..m) |i| {
                work[i] += ops.abs(T, a[i + j * lda]);
            };
            for (work) |t| if (value < t or std.math.isNan(t)) {
                value = t;
            };
        },
        .f => {
            var scale: F = 0;
            var sum: F = 1;
            for (0..n) |j| lassq(T, .of((a + j * lda)[0..m]), &scale, &sum);
            value = scale * @sqrt(sum);
        },
    }
    return value;
}

/// Blue's thresholds and scalings (`la_constants.f90`'s `tsml`, `tbig`,
/// `ssml`, `sbig`) from Fortran's `digits`, `minexponent`, `maxexponent`.
fn Blue(comptime F: type) type {
    return struct {
        const digits = std.math.floatFractionalBits(F) + 1;
        const minexp = std.math.floatExponentMin(F) + 1;
        const maxexp = std.math.floatExponentMax(F) + 1;
        const tsml = std.math.ldexp(@as(F, 1), -@divFloor(1 - minexp, 2));
        const tbig = std.math.ldexp(@as(F, 1), @divFloor(maxexp - digits + 1, 2));
        const ssml = std.math.ldexp(@as(F, 1), -@divFloor(minexp - digits, 2));
        const sbig = std.math.ldexp(@as(F, 1), -(-@divFloor(-(maxexp + digits - 1), 2)));
    };
}

/// `?LASSQ` (`dlassq.f90`, `zlassq.f90`, Blue's algorithm): updates
/// `scale` and `sumsq` so `scale² sumsq` gains `Σ |x(i)|²` (complex: both
/// parts' squares), in `x`'s element order. A NaN in either leaves both as
/// they are.
pub fn lassq(comptime T: type, x: ConstStrided(T), scale: *Float(T), sumsq: *Float(T)) void {
    const n = x.n;
    const F = Float(T);
    const c = Blue(F);
    if (std.math.isNan(scale.*) or std.math.isNan(sumsq.*)) return;
    if (sumsq.* == 0) scale.* = 1;
    if (scale.* == 0) {
        scale.* = 1;
        sumsq.* = 0;
    }
    if (n == 0) return;
    var notbig = true;
    var asml: F = 0;
    var amed: F = 0;
    var abig: F = 0;
    for (0..n) |i| {
        const v = x.get(i);
        const parts: []const F = if (comptime isComplex(T)) &.{ v.re, v.im } else &.{v};
        for (parts) |p| {
            const ax = @abs(p);
            if (ax > c.tbig) {
                abig += (ax * c.sbig) * (ax * c.sbig);
                notbig = false;
            } else if (ax < c.tsml) {
                if (notbig) asml += (ax * c.ssml) * (ax * c.ssml);
            } else {
                amed += ax * ax;
            }
        }
    }
    // the sum so far joins the accumulator of its size
    if (sumsq.* > 0) {
        const ax = scale.* * @sqrt(sumsq.*);
        if (ax > c.tbig) {
            if (scale.* > 1) {
                scale.* *= c.sbig;
                abig += scale.* * (scale.* * sumsq.*);
            } else {
                abig += scale.* * (scale.* * (c.sbig * (c.sbig * sumsq.*)));
            }
        } else if (ax < c.tsml) {
            if (notbig) {
                if (scale.* < 1) {
                    scale.* *= c.ssml;
                    asml += scale.* * (scale.* * sumsq.*);
                } else {
                    asml += scale.* * (scale.* * (c.ssml * (c.ssml * sumsq.*)));
                }
            }
        } else {
            amed += scale.* * (scale.* * sumsq.*);
        }
    }
    if (abig > 0) {
        if (amed > 0 or std.math.isNan(amed)) abig += (amed * c.sbig) * c.sbig;
        scale.* = 1 / c.sbig;
        sumsq.* = abig;
    } else if (asml > 0) {
        if (amed > 0 or std.math.isNan(amed)) {
            const m = @sqrt(amed);
            const s = @sqrt(asml) / c.ssml;
            const ymin = if (s > m) m else s;
            const ymax = if (s > m) s else m;
            scale.* = 1;
            sumsq.* = ymax * ymax * (1 + (ymin / ymax) * (ymin / ymax));
        } else {
            scale.* = 1 / c.ssml;
            sumsq.* = asml;
        }
    } else {
        scale.* = 1;
        sumsq.* = amed;
    }
}
