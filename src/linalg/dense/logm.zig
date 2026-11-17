//! The principal matrix logarithm by inverse scaling and squaring on the
//! Schur form: `A = U T Uᴴ`; square roots of `T` (`sqrtm`'s, so a real
//! `A` stays real) until `‖T - I‖₁ ≤ 1/4`; `log(I + X)` there by the
//! [8/8] Padé approximant in partial fractions, Gauss-Legendre's nodes and
//! weights (Higham, "Evaluating Padé approximants of the matrix logarithm",
//! SIAM J. Matrix Anal. Appl. 22(4), 2001); times `2^k`. Al-Mohy and
//! Higham's 2012 refinements (choosing the degree and the roots by norm
//! estimates) are not ported yet. Replaces `matrixfun.logm` (issue #334).

const std = @import("std");
const matrix = @import("../matrix.zig");
const vector = @import("../vector.zig");
const util = @import("util.zig");
const inverse = @import("inverse.zig");
const schur_form = @import("schur_form.zig");
const sqrtm = @import("sqrtm.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = matrix.Mat;
const ConstView = matrix.ConstView;
const isComplex = vector.isComplex;
const Real = vector.Real;
const ElementOf = util.ElementOf;
const constOf = util.constOf;
const copyOf = util.copyOf;

pub const LogmError = error{
    NotSquare,
    NoConvergence,
    /// A real `A` with a negative real eigenvalue: its principal logarithm
    /// is complex; pass `A` as complex.
    NoRealLogarithm,
    /// An eigenvalue is 0: no logarithm.
    Singular,
} || Allocator.Error || Io.Cancelable;

/// Gauss-Legendre's 8 nodes and weights, mapped to [0, 1].
const nodes = blk: {
    const t = [_]f64{ -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498, 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363 };
    var x: [8]f64 = undefined;
    for (&x, t) |*d, s| d.* = (s + 1) / 2;
    break :blk x;
};
const weights = blk: {
    const w = [_]f64{ 0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620, 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763 };
    var x: [8]f64 = undefined;
    for (&x, w) |*d, s| d.* = s / 2;
    break :blk x;
};

/// The principal logarithm (MATLAB's `logm`): `e^X = A`, the eigenvalues
/// of `X` with imaginary parts in (-π, π).
pub fn logm(allocator: Allocator, io: Io, a: anytype) LogmError!Mat(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const ca = constOf(a);
    if (ca.rows != ca.cols) return error.NotSquare;
    const n = ca.rows;
    var arena: std.heap.ArenaAllocator = .init(allocator);
    defer arena.deinit();
    const ar = arena.allocator();
    const s = try schur_form.schur(ar, io, ca);
    const t = s.t;
    // the eigenvalues on 1 x 1 blocks: none 0, none negative real
    var i: usize = 0;
    while (i < n) {
        if (comptime !isComplex(T)) {
            if (i + 1 < n and t.data[i + 1 + i * n] != 0) {
                i += 2;
                continue;
            }
        }
        const d = t.data[i * (n + 1)];
        if (comptime isComplex(T)) {
            if (d.re == 0 and d.im == 0) return error.Singular;
        } else {
            if (d == 0) return error.Singular;
            if (d < 0) return error.NoRealLogarithm;
        }
        i += 1;
    }
    var k: u32 = 0;
    while (distanceToIdentity(T, t) > 0.25) : (k += 1) {
        if (k == 64) return error.NoConvergence;
        sqrtm.rootTriangular(T, t) catch return error.NoRealLogarithm;
    }
    // X = T - I, then Σ w_j X (I + x_j X)⁻¹
    for (0..n) |d| t.data[d * (n + 1)] = sub(T, t.data[d * (n + 1)], one(T));
    const l = try Mat(T).zeros(ar, n, n);
    const m = try Mat(T).init(ar, n, n);
    for (nodes, weights) |x, w| {
        for (m.data, t.data) |*md, td| md.* = scal(T, @floatCast(x), td);
        for (0..n) |d| m.data[d * (n + 1)] = add(T, m.data[d * (n + 1)], one(T));
        const y = inverse.solve(ar, io, m, t) catch |e| switch (e) {
            error.Singular => return error.Singular,
            error.DimensionMismatch => unreachable,
            else => |x2| return x2,
        };
        for (l.data, y.data) |*ld, yd| ld.* = add(T, ld.*, scal(T, @floatCast(w), yd));
    }
    const f = std.math.ldexp(@as(Real(T), 1), @intCast(k));
    for (l.data) |*e| e.* = scal(T, f, e.*);
    const x = try sqrtm.similar(T, ar, io, s.u, l);
    const out = try Mat(T).init(allocator, n, n);
    @memcpy(out.data, x.data);
    return out;
}

fn distanceToIdentity(comptime T: type, t: Mat(T)) f64 {
    const n = t.rows;
    var worst: f64 = 0;
    for (0..n) |j| {
        var s: f64 = 0;
        for (0..n) |i| {
            const v = t.data[i + j * n];
            s += absOf(T, if (i == j) sub(T, v, one(T)) else v);
        }
        worst = @max(worst, s);
    }
    return worst;
}

fn absOf(comptime T: type, x: T) f64 {
    return @floatCast(if (comptime isComplex(T)) std.math.hypot(x.re, x.im) else @abs(x));
}

fn one(comptime T: type) T {
    return if (comptime isComplex(T)) T.init(1, 0) else 1;
}

fn add(comptime T: type, x: T, y: T) T {
    return if (comptime isComplex(T)) x.add(y) else x + y;
}

fn sub(comptime T: type, x: T, y: T) T {
    return if (comptime isComplex(T)) x.sub(y) else x - y;
}

fn scal(comptime T: type, c: Real(T), x: T) T {
    return if (comptime isComplex(T)) T.init(c * x.re, c * x.im) else c * x;
}

// ---- tests ----

const testing = std.testing;

test "logm: diagonal, rotation, complex, errors" {
    const al = testing.allocator;
    const io = testing.io;
    const V = ConstView(f64);
    const d = try logm(al, io, try V.init(&[_]f64{ std.math.e, 0, 0, 1 }, 2, 2, 2));
    defer d.deinit(al);
    try testing.expectApproxEqAbs(1, d.data[0], 1e-14);
    try testing.expectApproxEqAbs(0, d.data[3], 1e-14);
    // a rotation by 1 radian: log is [0 -1; 1 0], real
    const c = @cos(1.0);
    const s = @sin(1.0);
    const r = try logm(al, io, try V.init(&[_]f64{ c, s, -s, c }, 2, 2, 2));
    defer r.deinit(al);
    try testing.expectApproxEqAbs(0, r.data[0], 1e-13);
    try testing.expectApproxEqAbs(1, r.data[1], 1e-13);
    try testing.expectApproxEqAbs(-1, r.data[2], 1e-13);
    // [1 1; 0 1]: log is [0 1; 0 0]
    const j = try logm(al, io, try V.init(&[_]f64{ 1, 0, 1, 1 }, 2, 2, 2));
    defer j.deinit(al);
    try testing.expectApproxEqAbs(1, j.data[2], 1e-13);
    try testing.expectError(error.NoRealLogarithm, logm(al, io, try V.init(&[_]f64{ -1, 0, 0, 1 }, 2, 2, 2)));
    try testing.expectError(error.Singular, logm(al, io, try V.init(&[_]f64{ 0, 0, 0, 1 }, 2, 2, 2)));
    const C = std.math.Complex(f64);
    const z = try logm(al, io, try ConstView(C).init(&[_]C{.init(-1, 0)}, 1, 1, 1));
    defer z.deinit(al);
    try testing.expectApproxEqAbs(std.math.pi, z.data[0].im, 1e-14);
}
