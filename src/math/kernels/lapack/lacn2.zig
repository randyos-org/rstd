//! `xLACN2`: the 1-norm estimate of a matrix known by its products, Reference-LAPACK 3.12.0
//! (`lapack-netlib/SRC/dlacn2.f`, `zlacn2.f`, OpenBLAS v0.3.34-317; modified BSD).

const std = @import("std");
const types = @import("../blas/types.zig");
const strided = @import("../blas/strided.zig");
const ops = @import("ops.zig");

const Allocator = std.mem.Allocator;
const Float = types.Float;
const isComplex = types.isComplex;

/// `ITMAX` (`dlacn2.f`): the most power-method steps.
const itmax = 5;

/// `?LACN2` (`dlacn2.f`, `zlacn2.f`; N. J. Higham, "FORTRAN codes for
/// estimating the one-norm of a real or complex matrix", ACM TOMS 14(4),
/// 1988): a lower bound on `‖B‖₁`, with `v` (`n`) holding `B w` for the
/// `w` it found, `‖v‖₁ = est`. The reference's reverse communication is a
/// callback: `apply(ctx, adjoint, x)` sets `x := B x` (its `KASE = 1`) or,
/// with `adjoint`, `x := Bᴴ x` (`KASE = 2`) on `x`'s first `n`; any error
/// it returns ends the estimate. `x` (`n`) is workspace; `v` and `x` too
/// short are `error.SliceTooShort`.
pub fn lacn2(al: Allocator, comptime T: type, comptime E: type, n: usize, v: []T, x: []T, ctx: anytype, comptime apply: fn (@TypeOf(ctx), bool, [*]T) E!void) (E || Allocator.Error || strided.SliceError)!Float(T) {
    if (v.len < n or x.len < n) return error.SliceTooShort;
    if (n == 0) return 0;
    return estimate(al, T, E, n, v.ptr, x.ptr, ctx, apply);
}

/// `lacn2` after its checks.
fn estimate(al: Allocator, comptime T: type, comptime E: type, n: usize, v: [*]T, x: [*]T, ctx: anytype, comptime apply: fn (@TypeOf(ctx), bool, [*]T) E!void) (E || Allocator.Error)!Float(T) {
    const F = Float(T);
    const nf: F = @floatFromInt(n);
    const isgn = try al.alloc(bool, if (comptime isComplex(T)) 0 else n);
    defer al.free(isgn);
    for (0..n) |i| x[i] = ops.scalar(T, 1 / nf);
    try apply(ctx, false, x);
    if (n == 1) {
        v[0] = x[0];
        return ops.abs(T, v[0]);
    }
    var est = sum1(T, n, x);
    signs(T, n, x, isgn);
    try apply(ctx, true, x);
    var j = imax1(T, n, x);
    var iter: usize = 2;
    while (true) {
        @memset(x[0..n], ops.zero(T));
        x[j] = ops.one(T);
        try apply(ctx, false, x);
        @memcpy(v[0..n], x[0..n]);
        const estold = est;
        est = sum1(T, n, v);
        if (comptime !isComplex(T)) {
            // repeated signs: converged
            for (0..n) |i| {
                if ((x[i] >= 0) != isgn[i]) break;
            } else break;
        }
        if (est <= estold) break;
        signs(T, n, x, isgn);
        try apply(ctx, true, x);
        const jlast = j;
        j = imax1(T, n, x);
        const last = if (comptime isComplex(T)) ops.abs(T, x[jlast]) else x[jlast];
        if (last != ops.abs(T, x[j]) and iter < itmax) {
            iter += 1;
            continue;
        }
        break;
    }
    // the alternating-sign test vector
    var altsgn: F = 1;
    for (0..n) |i| {
        x[i] = ops.scalar(T, altsgn * (1 + @as(F, @floatFromInt(i)) / (nf - 1)));
        altsgn = -altsgn;
    }
    try apply(ctx, false, x);
    const temp = 2 * (sum1(T, n, x) / (3 * nf));
    if (temp > est) {
        @memcpy(v[0..n], x[0..n]);
        est = temp;
    }
    return est;
}

/// `DASUM`, or `DZSUM1` (moduli) for complex.
fn sum1(comptime T: type, n: usize, x: [*]const T) Float(T) {
    var s: Float(T) = 0;
    for (0..n) |i| s += ops.abs(T, x[i]);
    return s;
}

/// `IDAMAX`, or `IZMAX1` (by modulus) for complex.
fn imax1(comptime T: type, n: usize, x: [*]const T) usize {
    var best: usize = 0;
    var m: Float(T) = -1;
    for (0..n) |i| {
        const a = ops.abs(T, x[i]);
        if (a > m) {
            m = a;
            best = i;
        }
    }
    return best;
}

/// `x := sign(x)`: real `±1` (NaN: -1), kept in `isgn`; complex `x / |x|`,
/// or 1 at or below the safe minimum.
fn signs(comptime T: type, n: usize, x: [*]T, isgn: []bool) void {
    const F = Float(T);
    for (0..n) |i| {
        if (comptime isComplex(T)) {
            const a = ops.abs(T, x[i]);
            x[i] = if (a > std.math.floatMin(F)) .init(x[i].re / a, x[i].im / a) else ops.one(T);
        } else {
            isgn[i] = x[i] >= 0;
            x[i] = if (isgn[i]) 1 else -1;
        }
    }
}
