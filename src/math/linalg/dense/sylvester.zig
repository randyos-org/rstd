//! Sylvester and Lyapunov equations (MATLAB's `sylvester` and `lyap`) by
//! Bartels and Stewart: Schur forms (`schur`), then `trsyl`.

const std = @import("std");
const lapack = @import("../../kernels/lapack/__root__.zig");
const matrix = @import("../matrix.zig");
const ops = @import("../../kernels/lapack/ops.zig");
const util = @import("util.zig");
const products = @import("products.zig");
const schur_form = @import("schur_form.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = matrix.Mat;
const ConstView = matrix.ConstView;
const ElementOf = util.ElementOf;
const constOf = util.constOf;
const mul = products.mul;
const ctranspose = products.ctranspose;

/// `Singular` when `A` and `-B` have (nearly) common eigenvalues, so the
/// solution is not unique (`trsyl`'s perturbed case); `NoConvergence` from
/// a Schur form (a NaN too).
pub const SylvesterError = error{ NotSquare, DimensionMismatch, Singular, NoConvergence } || Allocator.Error || Io.Cancelable;

/// `X` with `A X + X B = C` (MATLAB's `sylvester(A, B, C)`), `A` `m` by `m`,
/// `B` `n` by `n`, `C` `m` by `n`.
pub fn sylvester(allocator: Allocator, io: Io, a: anytype, b: anytype, c: anytype) SylvesterError!Mat(ElementOf(@TypeOf(a))) {
    return twoSided(allocator, io, a, b, c, false);
}

/// `X` with `A X + X Aᴴ + Q = 0` (MATLAB's `lyap(A, Q)`), from one Schur
/// form of `A`.
pub fn lyap(allocator: Allocator, io: Io, a: anytype, q: anytype) SylvesterError!Mat(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const ca = constOf(a);
    const cq = constOf(q);
    if (ca.rows != ca.cols) return error.NotSquare;
    if (cq.rows != ca.rows or cq.cols != ca.rows) return error.DimensionMismatch;
    const sa = try schur_form.schur(allocator, io, ca);
    defer sa.deinit(allocator);
    return bartelsStewart(T, allocator, io, sa, sa, .c, cq, true);
}

/// `X` with `A X + X B + C = 0` (MATLAB's `lyap(A, B, C)`): `sylvester(A,
/// B, -C)`.
pub fn lyapSylvester(allocator: Allocator, io: Io, a: anytype, b: anytype, c: anytype) SylvesterError!Mat(ElementOf(@TypeOf(a))) {
    return twoSided(allocator, io, a, b, c, true);
}

/// `sylvester` and `lyapSylvester`: two Schur forms, `-C` with `negate`.
fn twoSided(allocator: Allocator, io: Io, a: anytype, b: anytype, c: anytype, negate: bool) SylvesterError!Mat(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const ca = constOf(a);
    const cb = constOf(b);
    const cc = constOf(c);
    if (ca.rows != ca.cols or cb.rows != cb.cols) return error.NotSquare;
    if (cc.rows != ca.rows or cc.cols != cb.rows) return error.DimensionMismatch;
    const sa = try schur_form.schur(allocator, io, ca);
    defer sa.deinit(allocator);
    const sb = try schur_form.schur(allocator, io, cb);
    defer sb.deinit(allocator);
    return bartelsStewart(T, allocator, io, sa, sb, .n, cc, negate);
}

/// `X = Ua Y Ubᴴ / scale` for `Ta Y + Y op(Tb) = scale Uaᴴ C Ub` (`-C`
/// with `negate`).
fn bartelsStewart(comptime T: type, allocator: Allocator, io: Io, sa: schur_form.Schur(T), sb: schur_form.Schur(T), tranb: lapack.Trans, c: ConstView(T), negate: bool) SylvesterError!Mat(T) {
    const m = sa.t.rows;
    const n = sb.t.rows;
    const uah = try ctranspose(allocator, io, sa.u);
    defer uah.deinit(allocator);
    const uahc = try mul(allocator, io, uah, c);
    defer uahc.deinit(allocator);
    const f = try mul(allocator, io, uahc, sb.u);
    defer f.deinit(allocator);
    if (negate) for (f.data) |*x| {
        x.* = ops.neg(T, x.*);
    };
    const sol = lapack.trsyl(allocator, io, T, .n, tranb, 1, m, n, sa.t.data, @max(m, 1), sb.t.data, @max(n, 1), f.data, @max(m, 1)) catch |e| return ops.narrow(SylvesterError, e);
    if (sol.perturbed) return error.Singular;
    const uay = try mul(allocator, io, sa.u, f);
    defer uay.deinit(allocator);
    const ubh = try ctranspose(allocator, io, sb.u);
    defer ubh.deinit(allocator);
    const x = try mul(allocator, io, uay, ubh);
    if (sol.scale != 1) for (x.data) |*v| {
        v.* = ops.divReal(T, v.*, sol.scale);
    };
    return x;
}
