//! The principal matrix square root by the Schur method: Björck and
//! Hammarling, "A Schur method for the square root of a matrix", Linear
//! Algebra Appl. 52/53, 1983, and for real `A` Higham's real Schur form
//! version, "Computing real square roots of a real matrix", Linear Algebra
//! Appl. 88/89, 1987: `A = U T Uᴴ`, then `R² = T` block by block (a 2 x 2
//! block per complex pair, its root in closed form; off the diagonal a
//! Sylvester equation of order at most 4), `X = U R Uᴴ`. Point by point,
//! not Deadman, Higham and Ralha's blocked form. A real `A` stays real.
//! Replaces `matrixfun.sqrtm` (issue #333).

const std = @import("std");
const mem = @import("../../../mem.zig");
const matrix = @import("../matrix.zig");
const vector = @import("../vector.zig");
const util = @import("util.zig");
const ops = @import("../../kernels/lapack/ops.zig");
const products = @import("products.zig");
const schur_form = @import("schur_form.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = matrix.Mat;
const ConstView = matrix.ConstView;
const isComplex = vector.isComplex;
const Real = vector.Real;
const ElementOf = util.ElementOf;
const constOf = util.constOf;

pub const SqrtmError = error{
    NotSquare,
    NoConvergence,
    /// A real `A` with a negative real eigenvalue: its principal root is
    /// complex; pass `A` as complex.
    NoRealSquareRoot,
} || Allocator.Error || Io.Cancelable;

/// The principal square root (MATLAB's `sqrtm`): `X² = A`, the eigenvalues
/// of `X` in the right half plane. A singular `A` may have none; then
/// infinities or NaNs come out, as MATLAB warns.
pub fn sqrtm(allocator: Allocator, io: Io, a: anytype) SqrtmError!Mat(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const ca = constOf(a);
    if (ca.rows != ca.cols) return error.NotSquare;
    var arena: std.heap.ArenaAllocator = .init(allocator);
    defer arena.deinit();
    const ar = arena.allocator();
    const s = try schur_form.schur(ar, io, ca);
    try rootTriangular(T, s.t);
    const x = try similar(T, ar, io, s.u, s.t);
    const out = try Mat(T).init(allocator, ca.rows, ca.cols);
    @memcpy(out.data, x.data);
    return out;
}

/// `U R Uᴴ`.
pub fn similar(comptime T: type, ar: Allocator, io: Io, u: Mat(T), r: Mat(T)) (Allocator.Error || Io.Cancelable)!Mat(T) {
    const ur = products.mul(ar, io, u, r) catch |e| return narrow(e);
    const uh = products.ctranspose(ar, io, u) catch |e| return narrow(e);
    return products.mul(ar, io, ur, uh) catch |e| narrow(e);
}

fn narrow(e: anytype) (Allocator.Error || Io.Cancelable) {
    return ops.narrow(Allocator.Error || Io.Cancelable, e);
}

/// `t` (upper triangular, or quasi-triangular for real `T`) becomes its
/// principal square root in place, column block by column block.
pub fn rootTriangular(comptime T: type, t: Mat(T)) error{NoRealSquareRoot}!void {
    const n = t.rows;
    // block starts: a 2 x 2 block where the subdiagonal is nonzero
    var j: usize = 0;
    while (j < n) {
        const q = blockSize(T, t, j);
        try rootDiagonal(T, t, j, q);
        // the blocks above, bottom up
        var i = j;
        while (i > 0) {
            const p = if (i >= 2 and blockSize(T, t, i - 2) == 2) @as(usize, 2) else 1;
            i -= p;
            offDiagonal(T, t, i, p, j, q);
        }
        j += q;
    }
}

/// 2 where a real `t` has a 2 x 2 block at `j`, else 1.
fn blockSize(comptime T: type, t: Mat(T), j: usize) usize {
    if (comptime isComplex(T)) return 1;
    return if (j + 1 < t.rows and t.data[j + 1 + j * t.rows] != 0) 2 else 1;
}

fn at(comptime T: type, t: Mat(T), i: usize, j: usize) *T {
    return &t.data[i + j * t.rows];
}

/// The diagonal block's root: `√t` for 1 x 1; for a standardized 2 x 2
/// block `B` with eigenvalues `θ ± iμ`, `α I + (B - θ I) / (2α)`, where
/// `α + iβ = √(θ + iμ)`.
fn rootDiagonal(comptime T: type, t: Mat(T), j: usize, q: usize) error{NoRealSquareRoot}!void {
    if (q == 1) {
        const x = at(T, t, j, j);
        if (comptime isComplex(T)) {
            x.* = std.math.complex.sqrt(x.*);
        } else {
            if (x.* < 0) return error.NoRealSquareRoot;
            x.* = @sqrt(x.*);
        }
        return;
    }
    // complex Schur forms are triangular
    if (comptime isComplex(T)) unreachable;
    const F = Real(T);
    const a = at(T, t, j, j).*;
    const b = at(T, t, j, j + 1).*;
    const c = at(T, t, j + 1, j).*;
    const d = at(T, t, j + 1, j + 1).*;
    const theta = (a + d) / 2;
    // (B - θ I)² = -μ² I for any 2 x 2 block with eigenvalues θ ± iμ
    const mu = @sqrt(-((a - d) * (a - d) / 4 + b * c));
    const root = std.math.complex.sqrt(std.math.Complex(F).init(theta, mu));
    const alpha = root.re;
    at(T, t, j, j).* = alpha + (a - theta) / (2 * alpha);
    at(T, t, j, j + 1).* = b / (2 * alpha);
    at(T, t, j + 1, j).* = c / (2 * alpha);
    at(T, t, j + 1, j + 1).* = alpha + (d - theta) / (2 * alpha);
}

/// Block `(i, j)` (`p` x `q`), all of `R`'s blocks right of it in row `i`
/// and below it in column `j` known: solves `R_ii X + X R_jj = T_ij - Σ
/// R_ik R_kj` (k strictly between) and stores `X` over `T_ij`.
fn offDiagonal(comptime T: type, t: Mat(T), i: usize, p: usize, j: usize, q: usize) void {
    var rhs: [2][2]T = undefined;
    for (0..p) |r| for (0..q) |c| {
        var s = at(T, t, i + r, j + c).*;
        for (i + p..j) |k| s = sub(T, s, mul(T, at(T, t, i + r, k).*, at(T, t, k, j + c).*));
        rhs[r][c] = s;
    };
    // (I_q ⊗ R_ii + R_jjᵀ ⊗ I_p) vec(X) = vec(rhs), vec column-major
    const m = p * q;
    var sys: [4][5]T = undefined;
    for (0..m) |row| {
        const r = row % p;
        const c = row / p;
        for (0..m) |col| {
            const r2 = col % p;
            const c2 = col / p;
            var v = zero(T);
            if (c == c2) v = add(T, v, at(T, t, i + r, i + r2).*);
            if (r == r2) v = add(T, v, at(T, t, j + c2, j + c).*);
            sys[row][col] = v;
        }
        sys[row][m] = rhs[r][c];
    }
    const x = gauss(T, m, &sys);
    for (0..m) |row| at(T, t, i + row % p, j + row / p).* = x[row];
}

/// Gaussian elimination with partial pivoting on an `m` x `m + 1`
/// augmented system, `m` at most 4.
fn gauss(comptime T: type, m: usize, sys: *[4][5]T) [4]T {
    for (0..m) |k| {
        var piv = k;
        for (k + 1..m) |r| {
            if (absOf(T, sys[r][k]) > absOf(T, sys[piv][k])) piv = r;
        }
        mem.swap([5]T, &sys[k], &sys[piv]);
        for (k + 1..m) |r| {
            const f = div(T, sys[r][k], sys[k][k]);
            for (k..m + 1) |c| sys[r][c] = sub(T, sys[r][c], mul(T, f, sys[k][c]));
        }
    }
    var x: [4]T = undefined;
    var k = m;
    while (k > 0) {
        k -= 1;
        var s = sys[k][m];
        for (k + 1..m) |c| s = sub(T, s, mul(T, sys[k][c], x[c]));
        x[k] = div(T, s, sys[k][k]);
    }
    return x;
}

fn absOf(comptime T: type, x: T) Real(T) {
    return if (comptime isComplex(T)) std.math.hypot(x.re, x.im) else @abs(x);
}

fn zero(comptime T: type) T {
    return if (comptime isComplex(T)) T.init(0, 0) else 0;
}

fn add(comptime T: type, x: T, y: T) T {
    return if (comptime isComplex(T)) x.add(y) else x + y;
}

fn sub(comptime T: type, x: T, y: T) T {
    return if (comptime isComplex(T)) x.sub(y) else x - y;
}

fn mul(comptime T: type, x: T, y: T) T {
    return if (comptime isComplex(T)) x.mul(y) else x * y;
}

fn div(comptime T: type, x: T, y: T) T {
    return if (comptime isComplex(T)) x.div(y) else x / y;
}

// ---- tests ----

const testing = std.testing;

fn expectSquare(comptime T: type, x: Mat(T), a: []const T, tol: f64) !void {
    const n = x.rows;
    for (0..n) |i| for (0..n) |j| {
        var s = zero(T);
        for (0..n) |k| s = add(T, s, mul(T, x.at(i, k), x.at(k, j)));
        const d = sub(T, s, a[i + j * n]);
        try testing.expect(absOf(T, d) < tol);
    };
}

test "sqrtm: X² = A, real with complex pairs, complex, negative eigenvalue" {
    const al = testing.allocator;
    const io = testing.io;
    // [4 1; 0 9]: triangular, root [2 0.2; 0 3]
    const t = [_]f64{ 4, 0, 1, 9 };
    const x = try sqrtm(al, io, try ConstView(f64).init(&t, 2, 2, 2));
    defer x.deinit(al);
    try expectSquare(f64, x, &t, 1e-13);
    try testing.expectApproxEqAbs(0.2, x.at(0, 1), 1e-14);
    // a rotation-scaling plus a real eigenvalue: complex pair 1 ± 2i and 3
    const r = [_]f64{ 1, -2, 0, 2, 1, 0, 1, 1, 3 };
    const y = try sqrtm(al, io, try ConstView(f64).init(&r, 3, 3, 3));
    defer y.deinit(al);
    try expectSquare(f64, y, &r, 1e-12);
    // magic(3) has a negative eigenvalue: no real root
    const m3 = [_]f64{ 8, 3, 4, 1, 5, 9, 6, 7, 2 };
    try testing.expectError(error.NoRealSquareRoot, sqrtm(al, io, try ConstView(f64).init(&m3, 3, 3, 3)));
    const C = std.math.Complex(f64);
    var mc: [9]C = undefined;
    for (&mc, m3) |*d, s| d.* = .init(s, 0);
    const z = try sqrtm(al, io, try ConstView(C).init(&mc, 3, 3, 3));
    defer z.deinit(al);
    try expectSquare(C, z, &mc, 1e-12);
}
