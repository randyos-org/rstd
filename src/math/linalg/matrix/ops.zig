//! The matrix arithmetic: construction, access, transposes, products, sums, the trace and the outer product.

const rstd_matrix = @import("rstd_matrix");
const vector = @import("../vector.zig");
const ops = @import("../../kernels/lapack/ops.zig");
const mroot = @import("../matrix.zig");
const mplumb = @import("plumbing.zig");

const Vec = vector.Vec;
const Elem = vector.Elem;
const count = vector.count;
const isComplex = vector.isComplex;
const splat = vector.splat;
const mulAdd = vector.mulAdd;
const StaticMat = mroot.StaticMat;
const Column = mplumb.Column;
const eMul = mplumb.eMul;
const eMulAdd = mplumb.eMulAdd;
const expectInline = mplumb.expectInline;
const expectSquare = mplumb.expectSquare;
const isMat = mplumb.isMat;
const native = mplumb.native;
const spirvMat = mplumb.spirvMat;
const spvOuter = mplumb.spvOuter;
const toNative = mplumb.toNative;

/// `I`, or ones on the diagonal of a rectangular `M` (MATLAB's `eye`).
pub fn identity(comptime M: type) M {
    var r: M = .{ .data = @splat(@splat(ops.zero(M.Element))) };
    for (0..@min(M.rows, M.cols)) |i| r.data[i][i] = ops.one(M.Element);
    return r;
}

/// The matrix with columns `c`; vector literals coerce.
pub fn fromCols(comptime M: type, c: [M.cols][M.rows]M.Element) M {
    var r: M = undefined;
    inline for (0..M.cols) |j| r.data[j] = c[j];
    return r;
}

/// The matrix with rows `r`, arrays as `fromCols` takes.
pub fn fromRows(comptime M: type, r: [M.rows][M.cols]M.Element) M {
    var m: M = undefined;
    for (0..M.rows) |i| for (0..M.cols) |j| {
        m.data[j][i] = r[i][j];
    };
    return m;
}

/// Column `j`.
pub fn col(m: anytype, j: usize) Column(@TypeOf(m).rows, @TypeOf(m).Element) {
    return m.data[j];
}

/// Row `i`.
pub fn row(m: anytype, i: usize) Column(@TypeOf(m).cols, @TypeOf(m).Element) {
    const M = @TypeOf(m);
    var r: [M.cols]M.Element = undefined;
    for (0..M.cols) |j| r[j] = m.data[j][i];
    return r;
}

fn Transposed(comptime M: type) type {
    if (!isMat(M)) @compileError("expects a StaticMat, found " ++ @typeName(M));
    return StaticMat(M.cols, M.rows, M.Element);
}

/// `mᵀ` (MATLAB's `m.'`).
pub fn transpose(m: anytype) Transposed(@TypeOf(m)) {
    const M = @TypeOf(m);
    if (comptime native(M.Element)) return .{ .data = rstd_matrix.transpose(toNative(m)) };
    var r: Transposed(M) = undefined;
    for (0..M.cols) |j| for (0..M.rows) |i| {
        r.data[i][j] = m.data[j][i];
    };
    return r;
}

/// The `n` by `n` matrix whose columns are `v`, transposed, as columns:
/// pffft's `VTRANSPOSE4` (`src/simd/pf_sse1_float.h`, marton78/pffft
/// efc1dd2d) for any `n`, through `transpose` (`@Matrix`'s where it has one).
pub fn transposeVecs(comptime n: comptime_int, comptime T: type, v: [n]Vec(n, T)) [n]Vec(n, T) {
    var c: [n][n]T = undefined;
    inline for (0..n) |j| c[j] = v[j];
    const t = transpose(fromCols(StaticMat(n, n, T), c));
    var r: [n]Vec(n, T) = undefined;
    inline for (0..n) |j| r[j] = t.data[j];
    return r;
}

/// `mᴴ` (MATLAB's `m'`); `mᵀ` for real `T`.
pub fn ctranspose(m: anytype) Transposed(@TypeOf(m)) {
    const M = @TypeOf(m);
    if (comptime !isComplex(M.Element)) return transpose(m);
    var r: Transposed(M) = undefined;
    for (0..M.cols) |j| for (0..M.rows) |i| {
        r.data[i][j] = m.data[j][i].conjugate();
    };
    return r;
}

fn Product(comptime A: type, comptime B: type) type {
    if (isMat(A) and isMat(B)) {
        if (A.Element != B.Element or A.cols != B.rows) @compileError("cannot multiply " ++ @typeName(A) ++ " by " ++ @typeName(B));
        return StaticMat(A.rows, B.cols, A.Element);
    }
    if (isMat(A)) {
        if (A.Element != Elem(B) or A.cols != count(B)) @compileError("cannot multiply " ++ @typeName(A) ++ " by " ++ @typeName(B));
        return if (@typeInfo(B) == .vector) Vec(A.rows, A.Element) else [A.rows]A.Element;
    }
    if (isMat(B)) {
        if (B.Element != Elem(A) or B.rows != count(A)) @compileError("cannot multiply " ++ @typeName(A) ++ " by " ++ @typeName(B));
        return if (@typeInfo(A) == .vector) Vec(B.cols, B.Element) else [B.cols]B.Element;
    }
    @compileError("expects a StaticMat operand, found " ++ @typeName(A) ++ " and " ++ @typeName(B));
}

/// `a b`: a matrix by a matrix, a matrix by a vector, or a vector, as a row,
/// by a matrix (MATLAB's and GLSL's `*`). `M v` sums the columns scaled by
/// `v`'s elements in column order, `v M` is a dot per column in element order
/// (no conjugate), as `@Matrix` (decision 16).
/// Past the inline limit, a compile error naming `rstd.math.linalg.mulInto`.
pub fn mul(a: anytype, b: anytype) Product(@TypeOf(a), @TypeOf(b)) {
    comptime expectInline(@TypeOf(a), @TypeOf(b), "rstd.math.linalg.mulInto(allocator, io, &out, &a, &b)");
    return mulAnySize(a, b);
}

pub fn mulAnySize(a: anytype, b: anytype) Product(@TypeOf(a), @TypeOf(b)) {
    const A = @TypeOf(a);
    const B = @TypeOf(b);
    if (comptime isMat(A) and isMat(B)) {
        if (comptime native(A.Element)) return .{ .data = toNative(a) * toNative(b) };
        return mulCols(a, b);
    }
    if (comptime isMat(A)) {
        if (comptime native(A.Element)) return toNative(a) * @as(Vec(A.cols, A.Element), b);
        return matVec(a, b);
    }
    if (comptime native(B.Element)) return @as(Vec(B.rows, B.Element), a) * toNative(b);
    return vecMat(a, b);
}

/// `mul` of two matrices without `@Matrix`.
pub fn mulCols(a: anytype, b: anytype) Product(@TypeOf(a), @TypeOf(b)) {
    var r: Product(@TypeOf(a), @TypeOf(b)) = undefined;
    for (0..@TypeOf(b).cols) |j| r.data[j] = matVec(a, b.data[j]);
    return r;
}

/// `M v` without `@Matrix`, as a `Column`.
pub fn matVec(m: anytype, v: anytype) Column(@TypeOf(m).rows, @TypeOf(m).Element) {
    const M = @TypeOf(m);
    const T = M.Element;
    if (comptime !isComplex(T)) {
        const V = Vec(M.rows, T);
        var acc: V = @as(V, m.data[0]) * splat(V, v[0]);
        inline for (1..M.cols) |k| acc = mulAdd(V, m.data[k], splat(V, v[k]), acc);
        return acc;
    }
    var acc: [M.rows]T = undefined;
    for (0..M.rows) |i| acc[i] = eMul(m.data[0][i], v[0]);
    for (1..M.cols) |k| for (0..M.rows) |i| {
        acc[i] = eMulAdd(m.data[k][i], v[k], acc[i]);
    };
    return acc;
}

/// `v M` without `@Matrix`.
pub fn vecMat(v: anytype, m: anytype) Product(@TypeOf(v), @TypeOf(m)) {
    const M = @TypeOf(m);
    const va: [M.rows]M.Element = v;
    var r: [M.cols]M.Element = undefined;
    for (0..M.cols) |j| {
        var acc = eMul(va[0], m.data[j][0]);
        for (1..M.rows) |i| acc = eMulAdd(va[i], m.data[j][i], acc);
        r[j] = acc;
    }
    return r;
}

/// `s m`.
pub fn scale(m: anytype, s: @TypeOf(m).Element) @TypeOf(m) {
    const M = @TypeOf(m);
    if (comptime native(M.Element)) return .{ .data = toNative(m) * s };
    var r: M = undefined;
    for (0..M.cols) |j| for (0..M.rows) |i| {
        r.data[j][i] = eMul(m.data[j][i], s);
    };
    return r;
}

/// `a + b`.
pub fn add(a: anytype, b: @TypeOf(a)) @TypeOf(a) {
    const M = @TypeOf(a);
    if (comptime native(M.Element)) return .{ .data = toNative(a) + toNative(b) };
    var r: M = undefined;
    for (0..M.cols) |j| for (0..M.rows) |i| {
        r.data[j][i] = ops.add(M.Element, a.data[j][i], b.data[j][i]);
    };
    return r;
}

/// `a - b`.
pub fn sub(a: anytype, b: @TypeOf(a)) @TypeOf(a) {
    const M = @TypeOf(a);
    if (comptime native(M.Element)) return .{ .data = toNative(a) - toNative(b) };
    var r: M = undefined;
    for (0..M.cols) |j| for (0..M.rows) |i| {
        r.data[j][i] = ops.sub(M.Element, a.data[j][i], b.data[j][i]);
    };
    return r;
}

/// `Σ m[i][i]`, in order.
pub fn trace(m: anytype) @TypeOf(m).Element {
    comptime expectSquare(@TypeOf(m));
    var acc = m.data[0][0];
    for (1..@TypeOf(m).rows) |i| acc = ops.add(@TypeOf(m).Element, acc, m.data[i][i]);
    return acc;
}

fn Outer(comptime A: type, comptime B: type) type {
    if (Elem(A) != Elem(B)) @compileError("outer of " ++ @typeName(A) ++ " and " ++ @typeName(B));
    return StaticMat(count(A), count(B), Elem(A));
}

/// `a bᴴ` (GLSL's `outerProduct`; MATLAB's `a * b'`).
pub fn outer(a: anytype, b: anytype) Outer(@TypeOf(a), @TypeOf(b)) {
    const M = Outer(@TypeOf(a), @TypeOf(b));
    if (vector.on_spirv and comptime spirvMat(M.rows, M.cols, M.Element)) {
        return .{ .data = spvOuter(M.rows, M.cols, M.Element, a, b) };
    }
    var r: M = undefined;
    if (comptime !isComplex(M.Element)) {
        const V = Vec(M.rows, M.Element);
        inline for (0..M.cols) |j| r.data[j] = @as(V, a) * splat(V, b[j]);
        return r;
    }
    for (0..M.cols) |j| for (0..M.rows) |i| {
        r.data[j][i] = eMul(a[i], b[j].conjugate());
    };
    return r;
}

/// `m`'s element in row `i`, column `j`, for the closed forms.
pub inline fn at(m: anytype, comptime i: comptime_int, comptime j: comptime_int) @TypeOf(m).Element {
    return m.data[j][i];
}

/// `p s - q r`.
inline fn det2(p: anytype, q: @TypeOf(p), r: @TypeOf(p), s: @TypeOf(p)) @TypeOf(p) {
    return mulAdd(@TypeOf(p), -q, r, p * s);
}

/// The 3 by 3 cofactor of row `i`, column `j`, sign included.
pub inline fn cofactor3(m: anytype, comptime i: comptime_int, comptime j: comptime_int) @TypeOf(m).Element {
    const ia = (i + 1) % 3;
    const ib = (i + 2) % 3;
    const ja = (j + 1) % 3;
    const jb = (j + 2) % 3;
    return det2(at(m, ia, ja), at(m, ia, jb), at(m, ib, ja), at(m, ib, jb));
}

/// Real `m` up to 3 by 3: cofactors along the first row.
pub inline fn closedDet(m: anytype) @TypeOf(m).Element {
    const T = @TypeOf(m).Element;
    return switch (@TypeOf(m).rows) {
        1 => at(m, 0, 0),
        2 => det2(at(m, 0, 0), at(m, 0, 1), at(m, 1, 0), at(m, 1, 1)),
        3 => mulAdd(T, at(m, 0, 2), cofactor3(m, 0, 2), mulAdd(T, at(m, 0, 1), cofactor3(m, 0, 1), at(m, 0, 0) * cofactor3(m, 0, 0))),
        else => unreachable,
    };
}
