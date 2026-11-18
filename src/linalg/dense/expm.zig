//! The matrix exponential by scaling and squaring with Padé approximants:
//! Al-Mohy and Higham, "A new scaling and squaring algorithm for the matrix
//! exponential", SIAM J. Matrix Anal. Appl. 31(3), 2009, Algorithm 3.1,
//! with exact 1-norms of the powers (the paper's `d_k`) in place of
//! estimates, and its `ell` from a row vector times `|A|`. The θ table and
//! unit roundoff are f64's (the paper's): f32 is then more accurate than
//! it needs, wider types stop at f64's accuracy.
//! A real `A` stays real. Replaces `matrixfun.expm` (issue #332).

const std = @import("std");
const matrix = @import("../matrix.zig");
const vector = @import("../vector.zig");
const ops = @import("../lapack/ops.zig");
const util = @import("util.zig");
const products = @import("products.zig");
const inverse = @import("inverse.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = matrix.Mat;
const ConstView = matrix.ConstView;
const isComplex = vector.isComplex;
const Real = vector.Real;
const ElementOf = util.ElementOf;
const constOf = util.constOf;
const copyOf = util.copyOf;

pub const ExpmError = error{NotSquare} || Allocator.Error || Io.Cancelable;

/// θ_m for m = 3, 5, 7, 9, 13 (the paper's Table 3.1, f64).
const theta = [_]f64{ 1.495585217958292e-2, 2.539398330063230e-1, 9.504178996162932e-1, 2.097847961257068e0, 5.371920351148152e0 };

/// The Padé numerators' coefficients `b_0 .. b_m`.
const b3 = [_]f64{ 120, 60, 12, 1 };
const b5 = [_]f64{ 30240, 15120, 3360, 420, 30, 1 };
const b7 = [_]f64{ 17297280, 8648640, 1995840, 277200, 25200, 1512, 56, 1 };
const b9 = [_]f64{ 17643225600, 8821612800, 2075673600, 302702400, 30270240, 2162160, 110880, 3960, 90, 1 };
const b13 = [_]f64{ 64764752532480000, 32382376266240000, 7771770303897600, 1187353796428800, 129060195264000, 10559470521600, 670442572800, 33522128640, 1323241920, 40840800, 960960, 16380, 182, 1 };

/// `e^A` (MATLAB's `expm`), owned by the caller. NaNs throughout when the
/// Padé denominator is singular (a NaN or infinity in `A`).
pub fn expm(allocator: Allocator, io: Io, a: anytype) ExpmError!Mat(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const ca = constOf(a);
    if (ca.rows != ca.cols) return error.NotSquare;
    const n = ca.rows;
    var arena: std.heap.ArenaAllocator = .init(allocator);
    defer arena.deinit();
    const ar = arena.allocator();
    const x = try padeScaled(T, ar, io, try copyOf(T, ar, ca));
    const out = try Mat(T).init(allocator, n, n);
    @memcpy(out.data, x.data);
    return out;
}

fn padeScaled(comptime T: type, ar: Allocator, io: Io, a: Mat(T)) ExpmError!Mat(T) {
    const n = a.rows;
    if (n == 0) return a;
    const a2 = try mul(T, ar, io, a, a);
    const a4 = try mul(T, ar, io, a2, a2);
    const a6 = try mul(T, ar, io, a2, a4);
    const d4 = std.math.pow(f64, norm1(T, a4), 1.0 / 4.0);
    const d6 = std.math.pow(f64, norm1(T, a6), 1.0 / 6.0);
    const eta1 = @max(d4, d6);
    if (eta1 <= theta[0] and try ell(T, ar, a, 3) == 0) return pade(T, ar, io, &b3, a, &.{ a2, a4, a6 });
    if (eta1 <= theta[1] and try ell(T, ar, a, 5) == 0) return pade(T, ar, io, &b5, a, &.{ a2, a4, a6 });
    const a8 = try mul(T, ar, io, a4, a4);
    const d8 = std.math.pow(f64, norm1(T, a8), 1.0 / 8.0);
    const eta3 = @max(d6, d8);
    if (eta3 <= theta[2] and try ell(T, ar, a, 7) == 0) return pade(T, ar, io, &b7, a, &.{ a2, a4, a6 });
    if (eta3 <= theta[3] and try ell(T, ar, a, 9) == 0) return pade(T, ar, io, &b9, a, &.{ a2, a4, a6, a8 });
    const a10 = try mul(T, ar, io, a4, a6);
    const d10 = std.math.pow(f64, norm1(T, a10), 1.0 / 10.0);
    const eta5 = @min(eta3, @max(d8, d10));
    var s: u32 = if (eta5 > theta[4]) @intFromFloat(@ceil(std.math.log2(eta5 / theta[4]))) else 0;
    if (s > 0) scaleAll(T, &.{ a, a2, a4, a6 }, s);
    const extra = try ell(T, ar, a, 13);
    if (extra > 0) {
        scaleAll(T, &.{ a, a2, a4, a6 }, extra);
        s += extra;
    }
    var x = try pade13(T, ar, io, a, a2, a4, a6);
    for (0..s) |_| x = try mul(T, ar, io, x, x);
    return x;
}

/// `a`, `a²`, `a⁴`, `a⁶` scaled by `2^-k`, `2^-2k`, `2^-4k`, `2^-6k`.
fn scaleAll(comptime T: type, ms: []const Mat(T), k: u32) void {
    for (ms, [_]u32{ 1, 2, 4, 6 }) |m, p| {
        const f = std.math.ldexp(@as(Real(T), 1), -@as(i32, @intCast(k * p)));
        for (m.data) |*e| e.* = scal(T, f, e.*);
    }
}

/// r_m(A) for m 3 to 9: `U = A Σ b_{2k+1} A^{2k}`, `V = Σ b_{2k} A^{2k}`,
/// then `(V - U) X = V + U`; `pows` holds `A²`, `A⁴`, ...
fn pade(comptime T: type, ar: Allocator, io: Io, b: []const f64, a: Mat(T), pows: []const Mat(T)) ExpmError!Mat(T) {
    const n = a.rows;
    const ui = try scaledIdentity(T, ar, n, b[1]);
    const v = try scaledIdentity(T, ar, n, b[0]);
    for (1..(b.len) / 2) |k| {
        addScaled(T, ui, b[2 * k + 1], pows[k - 1]);
        addScaled(T, v, b[2 * k], pows[k - 1]);
    }
    return finish(T, ar, io, try mul(T, ar, io, a, ui), v);
}

/// r_13(A), in the paper's (2.11) order.
fn pade13(comptime T: type, ar: Allocator, io: Io, a: Mat(T), a2: Mat(T), a4: Mat(T), a6: Mat(T)) ExpmError!Mat(T) {
    const n = a.rows;
    const b = b13;
    const u_in = try Mat(T).zeros(ar, n, n);
    addScaled(T, u_in, b[13], a6);
    addScaled(T, u_in, b[11], a4);
    addScaled(T, u_in, b[9], a2);
    const ui = try mul(T, ar, io, a6, u_in);
    addScaled(T, ui, b[7], a6);
    addScaled(T, ui, b[5], a4);
    addScaled(T, ui, b[3], a2);
    addIdentity(T, ui, b[1]);
    const v_in = try Mat(T).zeros(ar, n, n);
    addScaled(T, v_in, b[12], a6);
    addScaled(T, v_in, b[10], a4);
    addScaled(T, v_in, b[8], a2);
    const v = try mul(T, ar, io, a6, v_in);
    addScaled(T, v, b[6], a6);
    addScaled(T, v, b[4], a4);
    addScaled(T, v, b[2], a2);
    addIdentity(T, v, b[0]);
    return finish(T, ar, io, try mul(T, ar, io, a, ui), v);
}

/// `(V - U)⁻¹ (V + U)`.
fn finish(comptime T: type, ar: Allocator, io: Io, u: Mat(T), v: Mat(T)) ExpmError!Mat(T) {
    const p = try copyOf(T, ar, v.constView());
    const q = try copyOf(T, ar, v.constView());
    for (p.data, q.data, u.data) |*pe, *qe, ue| {
        pe.* = add(T, pe.*, ue);
        qe.* = add(T, qe.*, neg(T, ue));
    }
    return inverse.solve(ar, io, q, p) catch |e| switch (e) {
        error.Singular => {
            for (p.data) |*x| x.* = scal(T, std.math.nan(Real(T)), ops.one(T));
            return p;
        },
        error.DimensionMismatch => unreachable,
        else => |x| return x,
    };
}

/// The paper's `ℓ(A, m)`: extra squarings so the Padé error bound holds
/// for this `A`, from `‖ |A|^(2m+1) ‖₁`.
fn ell(comptime T: type, ar: Allocator, a: Mat(T), m: u32) Allocator.Error!u32 {
    const n = a.rows;
    const an = norm1(T, a);
    if (an == 0) return 0;
    // 1ᵀ |A|^p, a row at a time
    const r = try ar.alloc(f64, n);
    const t = try ar.alloc(f64, n);
    @memset(r, 1);
    for (0..2 * m + 1) |_| {
        for (t, 0..) |*tj, j| {
            var s: f64 = 0;
            for (r, 0..) |ri, i| s += ri * absOf(T, a.data[i + j * n]);
            tj.* = s;
        }
        @memcpy(r, t);
    }
    const p_norm = std.mem.max(f64, r);
    if (p_norm == 0) return 0;
    // 1 / |c_{2m+1}| = binom(2m, m) (2m + 1)!
    var recip: f64 = 1;
    for (0..m) |i| recip *= @as(f64, @floatFromInt(2 * m - i)) / @as(f64, @floatFromInt(i + 1));
    for (1..2 * m + 2) |i| recip *= @floatFromInt(i);
    const alpha = p_norm / (an * recip);
    const u = std.math.ldexp(@as(f64, 1), -53);
    const v = @ceil(std.math.log2(alpha / u) / @as(f64, @floatFromInt(2 * m)));
    return if (v > 0) @intFromFloat(v) else 0;
}

fn mul(comptime T: type, ar: Allocator, io: Io, a: Mat(T), b: Mat(T)) ExpmError!Mat(T) {
    return products.mul(ar, io, a, b) catch |e| switch (e) {
        error.DimensionMismatch => unreachable,
        else => |x| return x,
    };
}

fn norm1(comptime T: type, a: Mat(T)) f64 {
    var worst: f64 = 0;
    for (0..a.cols) |j| {
        var s: f64 = 0;
        for (a.data[j * a.rows ..][0..a.rows]) |x| s += absOf(T, x);
        worst = @max(worst, s);
    }
    return worst;
}

fn absOf(comptime T: type, x: T) f64 {
    return @floatCast(if (comptime isComplex(T)) std.math.hypot(x.re, x.im) else @abs(x));
}

fn scaledIdentity(comptime T: type, ar: Allocator, n: usize, c: f64) Allocator.Error!Mat(T) {
    const m = try Mat(T).zeros(ar, n, n);
    addIdentity(T, m, c);
    return m;
}

fn addIdentity(comptime T: type, m: Mat(T), c: f64) void {
    for (0..m.rows) |i| m.data[i * (m.rows + 1)] = add(T, m.data[i * (m.rows + 1)], scal(T, @floatCast(c), ops.one(T)));
}

fn addScaled(comptime T: type, dst: Mat(T), c: f64, src: Mat(T)) void {
    for (dst.data, src.data) |*d, s| d.* = add(T, d.*, scal(T, @floatCast(c), s));
}

fn scal(comptime T: type, c: Real(T), x: T) T {
    return if (comptime isComplex(T)) T.init(c * x.re, c * x.im) else c * x;
}

fn add(comptime T: type, x: T, y: T) T {
    return if (comptime isComplex(T)) x.add(y) else x + y;
}

fn neg(comptime T: type, x: T) T {
    return if (comptime isComplex(T)) x.neg() else -x;
}

// ---- tests ----

const testing = std.testing;

test "expm: diagonal, nilpotent, rotation, large norm, complex" {
    const al = testing.allocator;
    const io = testing.io;
    const V = ConstView(f64);
    // diag(1, 2)
    const d = try expm(al, io, try V.init(&[_]f64{ 1, 0, 0, 2 }, 2, 2, 2));
    defer d.deinit(al);
    try testing.expectApproxEqRel(std.math.e, d.data[0], 1e-15);
    try testing.expectApproxEqRel(std.math.e * std.math.e, d.data[3], 1e-15);
    // [0 1; 0 0]: I + A exactly
    const nil = try expm(al, io, try V.init(&[_]f64{ 0, 0, 1, 0 }, 2, 2, 2));
    defer nil.deinit(al);
    try testing.expectEqualSlices(f64, &.{ 1, 0, 1, 1 }, nil.data);
    // t [0 -1; 1 0]: a rotation by t, through the scaling for t = 30
    // `cos t`, `sin t` as literals: `@cos(30)` is wrong on x87 (issue #481)
    for ([_][3]f64{ .{ 0.5, 0.8775825618903728, 0.479425538604203 }, .{ 30, 0.15425144988758405, -0.9880316240928618 } }) |tcs| {
        const t = tcs[0];
        const r = try expm(al, io, try V.init(&[_]f64{ 0, t, -t, 0 }, 2, 2, 2));
        defer r.deinit(al);
        try testing.expectApproxEqAbs(tcs[1], r.data[0], 1e-12);
        try testing.expectApproxEqAbs(tcs[2], r.data[1], 1e-12);
    }
    const C = std.math.Complex(f64);
    const z = try expm(al, io, try ConstView(C).init(&[_]C{.init(0, std.math.pi)}, 1, 1, 1));
    defer z.deinit(al);
    try testing.expectApproxEqAbs(-1, z.data[0].re, 1e-15);
    try testing.expectApproxEqAbs(0, z.data[0].im, 1e-15);
}
