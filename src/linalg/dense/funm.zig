//! Functions of a matrix by the Schur-Parlett algorithm: Davies and Higham,
//! "A Schur-Parlett algorithm for computing matrix functions", SIAM J.
//! Matrix Anal. Appl. 25(2), 2003. The complex Schur form is reordered
//! (`trexc`) so eigenvalues within δ = 0.1 of each other (transitively)
//! share a diagonal block; each block's `f` by Taylor series about its
//! mean, the blocks above by the block Parlett recurrence (triangular
//! Sylvester equations). `sinm` and `cosm` go through `expm`, `atanm`
//! through `logm`; for a real `A` all three stay real by exact identities
//! (`cos A = Re e^{iA}`, `atan A = Im log(I + iA)`). Replaces
//! `matrixfun.zig`'s `applyMatrixFunction` (issue #335).

const std = @import("std");
const lapack = @import("../../math/kernels/lapack/__root__.zig");
const matrix = @import("../matrix.zig");
const vector = @import("../vector.zig");
const util = @import("util.zig");
const schur_form = @import("schur_form.zig");
const sqrtm = @import("sqrtm.zig");
const expm = @import("expm.zig");
const logm = @import("logm.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = matrix.Mat;
const ConstView = matrix.ConstView;
const isComplex = vector.isComplex;
const Real = vector.Real;
const ElementOf = util.ElementOf;
const constOf = util.constOf;

pub fn Complex(comptime T: type) type {
    return std.math.Complex(Real(T));
}

/// `f`'s `k`th derivative at `z`, for `funm`.
pub fn Derivatives(comptime C: type) type {
    return *const fn (z: C, k: usize) C;
}

pub const FunmError = error{
    NotSquare,
    NoConvergence,
    /// A block's Taylor series did not settle in 250 terms: `f` is not
    /// analytic enough there, or its eigenvalues spread too far.
    TaylorDiverged,
} || Allocator.Error || Io.Cancelable;

/// δ, the paper's blocking parameter.
const delta = 0.1;
const max_terms = 250;

/// `f(A)` for an `f` analytic on `A`'s spectrum, given every derivative
/// (MATLAB's `funm`); complex, as the Schur form it works on.
pub fn funm(allocator: Allocator, io: Io, a: anytype, f: Derivatives(Complex(ElementOf(@TypeOf(a))))) FunmError!Mat(Complex(ElementOf(@TypeOf(a)))) {
    const C = Complex(ElementOf(@TypeOf(a)));
    const F = Real(C);
    const ca = constOf(a);
    if (ca.rows != ca.cols) return error.NotSquare;
    const n = ca.rows;
    var arena: std.heap.ArenaAllocator = .init(allocator);
    defer arena.deinit();
    const ar = arena.allocator();
    const s = try schur_form.schurComplex(ar, io, ca);
    const t = s.t;
    const u = s.u;

    // clusters: eigenvalues within δ, transitively (the paper's Algorithm 4.1)
    const cluster = try ar.alloc(usize, n);
    for (cluster, 0..) |*c, i| c.* = i;
    for (0..n) |i| for (i + 1..n) |j| {
        if (dist(C, t.data[i * (n + 1)], t.data[j * (n + 1)]) <= delta) merge(cluster, i, j);
    };
    for (0..n) |i| cluster[i] = find(cluster, i);
    // reorder: each cluster contiguous, in order of first appearance
    var pos: usize = 0;
    while (pos < n) {
        const c = cluster[pos];
        var next = pos + 1;
        while (next < n) : (next += 1) {
            if (cluster[next] != c) continue;
            if (next != pos + 1) {
                var ifst = next;
                var ilst = pos + 1;
                lapack.trexc(io, C, .v, n, t.data, n, u.data, n, &ifst, &ilst) catch return error.NoConvergence;
                const moved = cluster[next];
                @memmove(cluster[pos + 2 .. next + 1], cluster[pos + 1 .. next]);
                cluster[pos + 1] = moved;
            }
            pos += 1;
        }
        pos += 1;
    }
    // block starts
    var starts: std.ArrayList(usize) = .empty;
    for (0..n) |i| if (i == 0 or cluster[i] != cluster[i - 1]) try starts.append(ar, i);
    try starts.append(ar, n);
    const b = starts.items;

    const fm = try Mat(C).zeros(ar, n, n);
    for (0..b.len - 1) |k| try taylor(C, F, ar, t, fm, b[k], b[k + 1], f);
    // the blocks above, column block by column block, bottom up
    for (1..b.len - 1) |jb| {
        var ib = jb;
        while (ib > 0) {
            ib -= 1;
            parlett(C, t, fm, b[ib], b[ib + 1], b[jb], b[jb + 1]);
        }
    }
    const x = try sqrtm.similar(C, ar, io, u, fm);
    const out = try Mat(C).init(allocator, n, n);
    @memcpy(out.data, x.data);
    return out;
}

fn dist(comptime C: type, x: C, y: C) Real(C) {
    return std.math.hypot(x.re - y.re, x.im - y.im);
}

fn find(parent: []usize, i: usize) usize {
    var r = i;
    while (parent[r] != r) r = parent[r];
    return r;
}

fn merge(parent: []usize, i: usize, j: usize) void {
    const a = find(parent, i);
    const b = find(parent, j);
    parent[@max(a, b)] = @min(a, b);
}

/// The diagonal block `lo..hi`'s `f`: `Σ f⁽ᵏ⁾(σ) / k! (T - σI)ᵏ` about the
/// mean `σ`, until two terms in a row are negligible after the block's
/// order (the nilpotent part's), `T`'s triangle kept.
fn taylor(comptime C: type, comptime F: type, ar: Allocator, t: Mat(C), fm: Mat(C), lo: usize, hi: usize, f: Derivatives(C)) FunmError!void {
    const n = t.rows;
    const m = hi - lo;
    var sigma = C.init(0, 0);
    for (lo..hi) |i| sigma = sigma.add(t.data[i * (n + 1)]);
    sigma = C.init(sigma.re / @as(F, @floatFromInt(m)), sigma.im / @as(F, @floatFromInt(m)));
    // M = T_block - σ I, P = Mᵏ / k!, F_block = Σ f⁽ᵏ⁾(σ) P
    const mm = try ar.alloc(C, m * m);
    const p = try ar.alloc(C, m * m);
    const q = try ar.alloc(C, m * m);
    for (0..m) |j| for (0..m) |i| {
        mm[i + j * m] = if (i <= j) t.data[lo + i + (lo + j) * n] else C.init(0, 0);
        if (i == j) mm[i + j * m] = mm[i + j * m].sub(sigma);
        p[i + j * m] = C.init(if (i == j) 1 else 0, 0);
    };
    const eps = std.math.floatEps(F);
    var quiet: u32 = 0;
    var k: usize = 0;
    while (k < max_terms) : (k += 1) {
        const d = f(sigma, k);
        var term: F = 0;
        var total: F = 0;
        for (0..m) |j| for (0..m) |i| {
            const v = d.mul(p[i + j * m]);
            const at = &fm.data[lo + i + (lo + j) * n];
            at.* = at.add(v);
            term = @max(term, std.math.hypot(v.re, v.im));
            total = @max(total, std.math.hypot(at.re, at.im));
        };
        if (k >= m and term <= eps * total) {
            quiet += 1;
            if (quiet == 2) return;
        } else quiet = 0;
        // P = P M / (k + 1)
        const kf: F = @floatFromInt(k + 1);
        for (0..m) |j| for (0..m) |i| {
            var s = C.init(0, 0);
            for (0..m) |l| s = s.add(p[i + l * m].mul(mm[l + j * m]));
            q[i + j * m] = C.init(s.re / kf, s.im / kf);
        };
        @memcpy(p, q);
    }
    return error.TaylorDiverged;
}

/// Block `(i, j)` of `F`, `i` above `j`: `T_ii X - X T_jj = F_ii T_ij -
/// T_ij F_jj + Σ (F_ik T_kj - T_ik F_kj)` (k between), column by column
/// (both blocks triangular).
fn parlett(comptime C: type, t: Mat(C), fm: Mat(C), ilo: usize, ihi: usize, jlo: usize, jhi: usize) void {
    const n = t.rows;
    const at = struct {
        fn f(x: Mat(C), i: usize, j: usize) C {
            return x.data[i + j * x.rows];
        }
    }.f;
    for (jlo..jhi) |c| {
        // the right side, column c, then X's earlier columns times T_jj
        for (ilo..ihi) |r| {
            var s = C.init(0, 0);
            for (ilo..ihi) |l| s = s.add(at(fm, r, l).mul(at(t, l, c)));
            for (jlo..jhi) |l| s = s.sub(at(t, r, l).mul(at(fm, l, c)));
            for (ihi..jlo) |l| s = s.add(at(fm, r, l).mul(at(t, l, c)).sub(at(t, r, l).mul(at(fm, l, c))));
            for (jlo..c) |l| s = s.add(at(fm, r, l).mul(at(t, l, c)));
            fm.data[r + c * n] = s;
        }
        // (T_ii - t_cc I) x = rhs, back substitution
        const tcc = at(t, c, c);
        var r = ihi;
        while (r > ilo) {
            r -= 1;
            var s = fm.data[r + c * n];
            for (r + 1..ihi) |l| s = s.sub(at(t, r, l).mul(fm.data[l + c * n]));
            fm.data[r + c * n] = s.div(at(t, r, r).sub(tcc));
        }
    }
}

/// The matrix cosine: `Re e^{iA}` for a real `A`, `(e^{iA} + e^{-iA}) / 2`
/// for a complex one.
pub fn cosm(allocator: Allocator, io: Io, a: anytype) FunmError!Mat(ElementOf(@TypeOf(a))) {
    return trig(.cos, allocator, io, a);
}

/// The matrix sine: `Im e^{iA}` for a real `A`, `(e^{iA} - e^{-iA}) / 2i`
/// for a complex one.
pub fn sinm(allocator: Allocator, io: Io, a: anytype) FunmError!Mat(ElementOf(@TypeOf(a))) {
    return trig(.sin, allocator, io, a);
}

fn trig(comptime which: enum { cos, sin }, allocator: Allocator, io: Io, a: anytype) FunmError!Mat(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const C = Complex(T);
    const ca = constOf(a);
    if (ca.rows != ca.cols) return error.NotSquare;
    const n = ca.rows;
    var arena: std.heap.ArenaAllocator = .init(allocator);
    defer arena.deinit();
    const ar = arena.allocator();
    // ±iA
    const ip = try Mat(C).init(ar, n, n);
    const im = try Mat(C).init(ar, n, n);
    for (0..n) |j| for (ip.data[j * n ..][0..n], im.data[j * n ..][0..n], ca.col(j)) |*p, *q, x| {
        const z: C = if (comptime isComplex(T)) x else C.init(x, 0);
        p.* = C.init(-z.im, z.re);
        q.* = C.init(z.im, -z.re);
    };
    const ep = expm.expm(ar, io, ip) catch |e| return narrowExp(e);
    const out = try Mat(T).init(allocator, n, n);
    if (comptime !isComplex(T)) {
        for (out.data, ep.data) |*o, e| o.* = if (which == .cos) e.re else e.im;
        return out;
    }
    errdefer out.deinit(allocator);
    const em = expm.expm(ar, io, im) catch |e| return narrowExp(e);
    for (out.data, ep.data, em.data) |*o, p, q| o.* = switch (which) {
        .cos => C.init((p.re + q.re) / 2, (p.im + q.im) / 2),
        // (p - q) / 2i = -i (p - q) / 2
        .sin => C.init((p.im - q.im) / 2, -(p.re - q.re) / 2),
    };
    return out;
}

fn narrowExp(e: expm.ExpmError) FunmError {
    return switch (e) {
        error.NotSquare => unreachable,
        else => |x| x,
    };
}

pub const AtanmError = FunmError || error{
    /// `±i` is an eigenvalue of `A`: `atan` has a pole there.
    Singular,
};

/// The matrix arctangent: `Im log(I + iA)` for a real `A`, `(i/2) (log(I -
/// iA) - log(I + iA))` for a complex one.
pub fn atanm(allocator: Allocator, io: Io, a: anytype) AtanmError!Mat(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const C = Complex(T);
    const ca = constOf(a);
    if (ca.rows != ca.cols) return error.NotSquare;
    const n = ca.rows;
    var arena: std.heap.ArenaAllocator = .init(allocator);
    defer arena.deinit();
    const ar = arena.allocator();
    // I + iA and I - iA
    const p = try Mat(C).init(ar, n, n);
    const q = try Mat(C).init(ar, n, n);
    for (0..n) |j| for (p.data[j * n ..][0..n], q.data[j * n ..][0..n], ca.col(j), 0..) |*x, *y, v, i| {
        const z: C = if (comptime isComplex(T)) v else C.init(v, 0);
        const d: Real(T) = if (i == j) 1 else 0;
        x.* = C.init(d - z.im, z.re);
        y.* = C.init(d + z.im, -z.re);
    };
    const lp = logm.logm(ar, io, p) catch |e| return narrowLog(e);
    const out = try Mat(T).init(allocator, n, n);
    if (comptime !isComplex(T)) {
        for (out.data, lp.data) |*o, l| o.* = l.im;
        return out;
    }
    errdefer out.deinit(allocator);
    const lq = logm.logm(ar, io, q) catch |e| return narrowLog(e);
    // (i/2) (lq - lp)
    for (out.data, lp.data, lq.data) |*o, x, y| o.* = C.init(-(y.im - x.im) / 2, (y.re - x.re) / 2);
    return out;
}

fn narrowLog(e: logm.LogmError) AtanmError {
    return switch (e) {
        error.NotSquare, error.NoRealLogarithm => unreachable,
        else => |x| x,
    };
}

// ---- tests ----

const testing = std.testing;

const Cx = std.math.Complex(f64);

fn expDerivatives(z: Cx, k: usize) Cx {
    _ = k;
    return std.math.complex.exp(z);
}

test "funm: exp by Schur-Parlett matches expm, close and repeated eigenvalues" {
    const al = testing.allocator;
    const io = testing.io;
    // eigenvalues 1, 1.05 (one block) and 3: [1 1 0; 0 1.05 1; 0 0 3]
    const a = [_]f64{ 1, 0, 0, 1, 1.05, 0, 0, 1, 3 };
    const v: ConstView(f64) = try .init(&a, 3, 3, 3);
    const got = try funm(al, io, v, expDerivatives);
    defer got.deinit(al);
    const want = try expm.expm(al, io, v);
    defer want.deinit(al);
    for (got.data, want.data) |g, w| {
        try testing.expectApproxEqAbs(w, g.re, 1e-12);
        try testing.expectApproxEqAbs(0, g.im, 1e-12);
    }
}

test "cosm, sinm and atanm: real stays real, identities hold" {
    const al = testing.allocator;
    const io = testing.io;
    const a = [_]f64{ 0.5, 0.2, -0.1, 0.3 };
    const v: ConstView(f64) = try .init(&a, 2, 2, 2);
    const c = try cosm(al, io, v);
    defer c.deinit(al);
    const s = try sinm(al, io, v);
    defer s.deinit(al);
    // cos² + sin² = I
    for (0..2) |i| for (0..2) |j| {
        var sum: f64 = 0;
        for (0..2) |k| sum += c.at(i, k) * c.at(k, j) + s.at(i, k) * s.at(k, j);
        try testing.expectApproxEqAbs(@as(f64, if (i == j) 1 else 0), sum, 1e-13);
    };
    // atan of a diagonal matrix
    const t = try atanm(al, io, try ConstView(f64).init(&[_]f64{ 1, 0, 0, 0.5 }, 2, 2, 2));
    defer t.deinit(al);
    try testing.expectApproxEqAbs(std.math.pi / 4.0, t.data[0], 1e-14);
    try testing.expectApproxEqAbs(std.math.atan(@as(f64, 0.5)), t.data[3], 1e-14);
}
