//! The generalized Lyapunov equation `A X Eᴴ + E X Aᴴ + Q = 0` (MATLAB's
//! `lyap(A, Q, [], E)`) by Penzl's generalized Bartels-Stewart method (T.
//! Penzl, "Numerical solution of generalized Lyapunov equations", Adv.
//! Comput. Math. 8, 1998; SLICOT's `SG03AD`), in the triangular form: the
//! complex QZ form `(S, T)` of `(A, E)`, then `S Y Tᴴ + T Y Sᴴ = -Q̃` a
//! column at a time.

const std = @import("std");
const matrix = @import("../matrix.zig");
const vector = @import("../vector.zig");
const util = @import("util.zig");
const products = @import("products.zig");
const schur_form = @import("schur_form.zig");
const sylvester = @import("sylvester.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = matrix.Mat;
const isComplex = vector.isComplex;
const Real = vector.Real;
const ElementOf = util.ElementOf;
const constOf = util.constOf;
const complexOf = util.complexOf;
const mul = products.mul;
const ctranspose = products.ctranspose;
const SylvesterError = sylvester.SylvesterError;

/// `X` with `A X Eᴴ + E X Aᴴ + Q = 0` (MATLAB's `lyap(A, Q, [], E)`), all
/// `n` by `n`. `Singular` when two eigenvalues of the pencil `(A, E)` have
/// `λ_i + conj(λ_j)` (nearly) zero, an infinite one included; real data is
/// solved in complex arithmetic and the real part returned.
pub fn lyapGeneralized(allocator: Allocator, io: Io, a: anytype, q: anytype, e: anytype) SylvesterError!Mat(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const ca = constOf(a);
    const cq = constOf(q);
    const ce = constOf(e);
    if (ca.rows != ca.cols or ce.rows != ce.cols) return error.NotSquare;
    const n = ca.rows;
    if (ce.rows != n or cq.rows != n or cq.cols != n) return error.DimensionMismatch;
    // `Q A Z = S`, `Q E Z = T`, so `Y = Zᴴ X Z` and `Q̃ = Q Q0 Qᴴ`
    const f = try schur_form.qzComplex(allocator, io, ca, ce);
    defer f.deinit(allocator);
    const q0 = try complexOf(T, allocator, cq);
    defer q0.deinit(allocator);
    const qq = try mul(allocator, io, f.q, q0);
    defer qq.deinit(allocator);
    const qh = try ctranspose(allocator, io, f.q);
    defer qh.deinit(allocator);
    const y = try mul(allocator, io, qq, qh);
    defer y.deinit(allocator);
    try solveTriangular(@TypeOf(y.data[0]), allocator, f.aa, f.bb, y);
    const zy = try mul(allocator, io, f.z, y);
    defer zy.deinit(allocator);
    const zh = try ctranspose(allocator, io, f.z);
    defer zh.deinit(allocator);
    const x = try mul(allocator, io, zy, zh);
    if (comptime isComplex(T)) return x;
    defer x.deinit(allocator);
    const r = try Mat(T).init(allocator, n, n);
    for (r.data, x.data) |*d, v| d.* = v.re;
    return r;
}

/// `Y` (over `R` in `y`) with `S Y Tᴴ + T Y Sᴴ + R = 0`, `S` and `T` upper
/// triangular: from the last column, `(conj(T(j, j)) S + conj(S(j, j)) T)
/// y_j = -r_j - S Σ y_l conj(T(j, l)) - T Σ y_l conj(S(j, l))`, `l > j`. A
/// pivot under `trsyl`'s `smin` (`dtrsyl.f`: `eps max(|S|, |T|)`, at least
/// the smallest normal) is `Singular`.
fn solveTriangular(comptime C: type, allocator: Allocator, s: Mat(C), t: Mat(C), y: Mat(C)) SylvesterError!void {
    const F = @FieldType(C, "re");
    const zero = C.init(0, 0);
    const n = s.rows;
    var big: F = 0;
    for (s.data, t.data) |u, v| big = @max(big, @max(u.magnitude(), v.magnitude()));
    const smin = @max(std.math.floatEps(F) * big, std.math.floatMin(F));
    const w1 = try allocator.alloc(C, n);
    defer allocator.free(w1);
    const w2 = try allocator.alloc(C, n);
    defer allocator.free(w2);
    var j = n;
    while (j > 0) {
        j -= 1;
        const yj = y.data[j * n ..][0..n];
        @memset(w1, zero);
        @memset(w2, zero);
        for (j + 1..n) |l| {
            const ct = t.data[j + l * n].conjugate();
            const cs = s.data[j + l * n].conjugate();
            for (w1, w2, y.data[l * n ..][0..n]) |*p, *r, v| {
                p.* = p.add(v.mul(ct));
                r.* = r.add(v.mul(cs));
            }
        }
        for (yj, 0..) |*v, i| {
            var acc = v.neg();
            for (i..n) |k| acc = acc.sub(s.data[i + k * n].mul(w1[k]).add(t.data[i + k * n].mul(w2[k])));
            v.* = acc;
        }
        const tj = t.data[j + j * n].conjugate();
        const sj = s.data[j + j * n].conjugate();
        var i = n;
        while (i > 0) {
            i -= 1;
            var acc = yj[i];
            for (i + 1..n) |k| acc = acc.sub(tj.mul(s.data[i + k * n]).add(sj.mul(t.data[i + k * n])).mul(yj[k]));
            const d = tj.mul(s.data[i + i * n]).add(sj.mul(t.data[i + i * n]));
            if (d.magnitude() < smin) return error.Singular;
            yj[i] = acc.div(d);
        }
    }
}
