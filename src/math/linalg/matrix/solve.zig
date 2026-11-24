//! `det`, `inv` and `solve`: cofactors up to 3 by 3 for real types, LU past them and for complex.

const mem = @import("../../../mem.zig");
const vector = @import("../vector.zig");
const ops = @import("../../kernels/lapack/ops.zig");
const mops = @import("ops.zig");
const mplumb = @import("plumbing.zig");

const Elem = vector.Elem;
const count = vector.count;
const isComplex = vector.isComplex;
const at = mops.at;
const closedDet = mops.closedDet;
const cofactor3 = mops.cofactor3;
const eDiv = mplumb.eDiv;
const eMul = mplumb.eMul;
const eMulAdd = mplumb.eMulAdd;
const expectInline = mplumb.expectInline;
const expectSquare = mplumb.expectSquare;
const isMat = mplumb.isMat;
const spirvMat = mplumb.spirvMat;
const spvDeterminant = mplumb.spvDeterminant;
const spvInverse = mplumb.spvInverse;

/// Whether `det`, `inv` and `solve` take `T` through LU at size `n`: all
/// complex `T`, real past 3 by 3; real `det` and `inv` up to 3 by 3 take
/// cofactors.
fn byLu(comptime n: comptime_int, comptime T: type) bool {
    return isComplex(T) or n > 3;
}

/// LU with partial pivoting over `a`'s columns, in place (`L`'s unit diagonal
/// implied), as LAPACK's `getf2` but dividing by the pivot, the largest
/// `cabs1` in its column, the first of ties. Row `i` of the factors is row
/// `perm[i]` of `a`; `odd` is the swaps' parity. `error.Singular` at the first
/// zero pivot, `a` then partly factored.
pub fn luInPlace(comptime n: comptime_int, comptime T: type, a: *[n][n]T, perm: *[n]usize, odd: *bool) error{Singular}!void {
    for (perm, 0..) |*p, i| p.* = i;
    odd.* = false;
    for (0..n) |k| {
        var p = k;
        var big = ops.abs1(T, a[k][k]);
        for (k + 1..n) |i| {
            const x = ops.abs1(T, a[k][i]);
            if (x > big) {
                p = i;
                big = x;
            }
        }
        if (big == 0) return error.Singular;
        if (p != k) {
            for (a) |*c| mem.swap(T, &c[k], &c[p]);
            mem.swap(usize, &perm[k], &perm[p]);
            odd.* = !odd.*;
        }
        for (k + 1..n) |i| a[k][i] = eDiv(a[k][i], a[k][k]);
        for (k + 1..n) |j| for (k + 1..n) |i| {
            a[j][i] = eMulAdd(ops.neg(T, a[k][i]), a[j][k], a[j][i]);
        };
    }
}

/// `x` in `A x = b` from `luInPlace`'s factors, `b` in row order of `A`.
pub fn luSolve(comptime n: comptime_int, comptime T: type, lu: *const [n][n]T, perm: *const [n]usize, b: [n]T) [n]T {
    var x: [n]T = undefined;
    for (&x, perm) |*xi, p| xi.* = b[p];
    for (1..n) |i| for (0..i) |k| {
        x[i] = eMulAdd(ops.neg(T, lu[k][i]), x[k], x[i]);
    };
    var i: usize = n;
    while (i > 0) {
        i -= 1;
        for (i + 1..n) |k| x[i] = eMulAdd(ops.neg(T, lu[k][i]), x[k], x[i]);
        x[i] = eDiv(x[i], lu[i][i]);
    }
    return x;
}

/// The determinant (MATLAB's `det`, GLSL's `determinant`): real `m` up to 3
/// by 3 by cofactors, else the product of LU's pivots in order; 0 when
/// singular. Past the inline limit, a compile error naming `rstd.math.linalg.det`.
pub fn det(m: anytype) @TypeOf(m).Element {
    comptime expectInline(@TypeOf(m), @TypeOf(m), "rstd.math.linalg.det(allocator, io, &m)");
    return detAnySize(m);
}

pub fn detAnySize(m: anytype) @TypeOf(m).Element {
    const M = @TypeOf(m);
    comptime expectSquare(M);
    const n = M.rows;
    const T = M.Element;
    if (vector.on_spirv and comptime spirvMat(n, n, T)) return spvDeterminant(m);
    if (comptime !byLu(n, T)) return closedDet(m);
    var lu = m.data;
    var perm: [n]usize = undefined;
    var odd: bool = undefined;
    luInPlace(n, T, &lu, &perm, &odd) catch return ops.zero(T);
    var acc = lu[0][0];
    for (1..n) |k| acc = eMul(acc, lu[k][k]);
    return if (odd) ops.neg(T, acc) else acc;
}

/// `m⁻¹` (MATLAB's `inv`, GLSL's `inverse`): real `m` up to 3 by 3 as its
/// adjugate over `det`, `error.Singular` when that is 0; else by LU,
/// `error.Singular` at a zero pivot. Past the inline limit, a compile error
/// naming `rstd.math.linalg.invInto`.
pub fn inv(m: anytype) error{Singular}!@TypeOf(m) {
    comptime expectInline(@TypeOf(m), @TypeOf(m), "rstd.math.linalg.invInto(allocator, io, &out, &m)");
    return invAnySize(m);
}

pub fn invAnySize(m: anytype) error{Singular}!@TypeOf(m) {
    const M = @TypeOf(m);
    comptime expectSquare(M);
    const n = M.rows;
    const T = M.Element;
    if (vector.on_spirv and comptime spirvMat(n, n, T)) {
        if (spvDeterminant(m) == 0) return error.Singular;
        // a local: a literal in the payload is a SPIR-V backend TODO (issue #200)
        const r: M = .{ .data = spvInverse(m) };
        return r;
    }
    var r: M = undefined;
    if (comptime !byLu(n, T)) {
        const d = closedDet(m);
        if (d == 0) return error.Singular;
        switch (n) {
            1 => r.data[0][0] = 1 / d,
            2 => r = .{ .data = .{
                .{ at(m, 1, 1) / d, -at(m, 1, 0) / d },
                .{ -at(m, 0, 1) / d, at(m, 0, 0) / d },
            } },
            3 => inline for (0..3) |j| inline for (0..3) |i| {
                r.data[j][i] = cofactor3(m, j, i) / d;
            },
            else => unreachable,
        }
        return r;
    }
    var lu = m.data;
    var perm: [n]usize = undefined;
    var odd: bool = undefined;
    try luInPlace(n, T, &lu, &perm, &odd);
    for (0..n) |j| {
        var e: [n]T = @splat(ops.zero(T));
        e[j] = ops.one(T);
        r.data[j] = luSolve(n, T, &lu, &perm, e);
    }
    return r;
}

/// `a \ b` (MATLAB's `mldivide`) for square `a` and `b` a vector or a matrix
/// of right-hand sides, by LU at every size; `error.Singular` at a zero pivot.
/// Past the inline limit, a compile error naming `rstd.math.linalg.solveInto`.
pub fn solve(a: anytype, b: anytype) error{Singular}!@TypeOf(b) {
    comptime expectInline(@TypeOf(a), @TypeOf(b), "rstd.math.linalg.solveInto(allocator, io, &out, &a, &b)");
    return solveAnySize(a, b);
}

pub fn solveAnySize(a: anytype, b: anytype) error{Singular}!@TypeOf(b) {
    const M = @TypeOf(a);
    const B = @TypeOf(b);
    comptime expectSquare(M);
    const n = M.rows;
    const T = M.Element;
    comptime {
        const rhs_rows = if (isMat(B)) B.rows else count(B);
        const RhsT = if (isMat(B)) B.Element else Elem(B);
        if (RhsT != T or rhs_rows != n) @compileError("cannot solve " ++ @typeName(M) ++ " against " ++ @typeName(B));
    }
    var lu = a.data;
    var perm: [n]usize = undefined;
    var odd: bool = undefined;
    try luInPlace(n, T, &lu, &perm, &odd);
    if (comptime !isMat(B)) return luSolve(n, T, &lu, &perm, b);
    var r: B = undefined;
    for (0..B.cols) |j| r.data[j] = luSolve(n, T, &lu, &perm, b.data[j]);
    return r;
}
