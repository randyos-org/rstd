//! The least-squares drivers' scaling of `A` and `B` into `[smlnum,
//! bignum]` and back (`dgelsy.f`, `dgelss.f`, `dgelsd.f`, Reference-LAPACK
//! 3.12.0 as OpenBLAS v0.3.34-317 ships it; modified BSD).

const std = @import("std");
const types = @import("../blas/types.zig");
const ops = @import("ops.zig");
const auxiliary = @import("auxiliary.zig");
const norm = @import("norm.zig");

const Allocator = std.mem.Allocator;
const Mat = ops.Mat;
const Float = types.Float;

/// What was scaled: `iascl` and `ibscl` are the bounds `A` and `B` went to,
/// 0 where they were left alone.
pub fn Scale(comptime F: type) type {
    return struct { anrm: F, iascl: F, bnrm: F, ibscl: F };
}

/// Scales `a` and `b` (`B`'s first `m` rows) where their largest entries
/// fall outside `[smlnum, bignum]`, `smlnum = DLAMCH('S') / DLAMCH('P')`.
/// Null for a zero `A`, nothing scaled.
pub fn scaleIn(al: Allocator, comptime T: type, a: Mat(T), b: Mat(T)) Allocator.Error!?Scale(Float(T)) {
    const F = Float(T);
    const smlnum = std.math.floatMin(F) / std.math.floatEps(F);
    const bignum = 1 / smlnum;
    const anrm = try norm.matNorm(al, T, .max_abs, a.rows, a.cols, a.ptr, a.ld);
    if (anrm == 0) return null;
    var sc: Scale(F) = .{ .anrm = anrm, .iascl = 0, .bnrm = 0, .ibscl = 0 };
    sc.iascl = to(F, anrm, smlnum, bignum);
    if (sc.iascl != 0) auxiliary.lascl(T, false, anrm, sc.iascl, a);
    sc.bnrm = try norm.matNorm(al, T, .max_abs, b.rows, b.cols, b.ptr, b.ld);
    sc.ibscl = to(F, sc.bnrm, smlnum, bignum);
    if (sc.ibscl != 0) auxiliary.lascl(T, false, sc.bnrm, sc.ibscl, b);
    return sc;
}

fn to(comptime F: type, nrm: F, smlnum: F, bignum: F) F {
    if (nrm > 0 and nrm < smlnum) return smlnum;
    if (nrm > bignum) return bignum;
    return 0;
}

/// `X` (`n x nrhs`) back from the scaled problem's solution.
pub fn unscale(comptime T: type, sc: Scale(Float(T)), x: Mat(T)) void {
    if (sc.iascl != 0 and sc.iascl == sc.ibscl) {
        // the constants cancel: one step, as two may flush or overflow an
        // entry (`sgelsy.f`, OpenBLAS `6b01b70dc`; its `dgelsy.f` lags)
        auxiliary.lascl(T, false, sc.anrm, sc.bnrm, x);
        return;
    }
    if (sc.iascl != 0) auxiliary.lascl(T, false, sc.anrm, sc.iascl, x);
    if (sc.ibscl != 0) auxiliary.lascl(T, false, sc.ibscl, sc.bnrm, x);
}
