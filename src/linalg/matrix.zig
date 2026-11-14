//! `rstd.linalg`'s matrices (design.md, "Planned: the public `rstd.linalg`
//! layer"): `StaticMat(R, C, T)` and free functions over it, `T` a float or
//! its `std.math.Complex` (interleaved), inferring the sizes and `T`; vectors
//! are `linalg.vector`'s. Shaderlib re-exports them under GLSL names.

const std = @import("std");
const rstd_matrix = @import("rstd_matrix");
const vector = @import("vector.zig");
const ops = @import("lapack/ops.zig");

const Vec = vector.Vec;
const Elem = vector.Elem;
const count = vector.count;
const isComplex = vector.isComplex;
const Real = vector.Real;
const splat = vector.splat;
const mulAdd = vector.mulAdd;
const root = @import("root");

/// `rstd.linalg`'s build settings, from a root's `rstd_linalg_options`
/// (design.md, Dispatch, decision 3).
pub const Options = struct {
    /// `L`, per dimension: at or below it the inline kernels, past it the
    /// BLAS and LAPACK. At least 4, for shaderlib's `Mat2` to `Mat4`.
    inline_limit: usize = 8,
};

/// The root's options, or the defaults.
pub const options: Options = if (@hasDecl(root, "rstd_linalg_options")) root.rstd_linalg_options else .{};

pub const inline_limit = options.inline_limit;

comptime {
    if (inline_limit < 4) @compileError(std.fmt.comptimePrint("rstd_linalg_options.inline_limit is {d}, below 4 (shaderlib's Mat4)", .{inline_limit}));
}

/// `X`'s dimensions: a matrix's rows and columns, a vector's length.
fn dimsOf(comptime X: type) [2]comptime_int {
    return if (isMat(X)) .{ X.rows, X.cols } else .{ count(X), 1 };
}

/// Past the inline limit a `StaticMat` takes `rstd.linalg`'s `instead`,
/// with an allocator and `io` (design.md, Dispatch); SPIR-V has no limit.
fn expectInline(comptime X: type, comptime Y: type, comptime instead: []const u8) void {
    if (vector.on_spirv) return;
    for (dimsOf(X) ++ dimsOf(Y)) |d| if (d > inline_limit)
        @compileError(std.fmt.comptimePrint("{s} and {s} pass the inline limit ({d}): use {s}", .{ @typeName(X), @typeName(Y), inline_limit, instead }));
}

/// `mul`, `det`, `inv` and `solve` at every size, past the inline limit
/// too: `rstd.linalg`'s inline path (its tests run it to 16) and timing.
pub const any_size = struct {
    pub const mul = mulAnySize;
    pub const det = detAnySize;
    pub const inv = invAnySize;
    pub const solve = solveAnySize;
};

/// `R` rows by `C` columns of `T`, packed column-major: `data[j]` is column
/// `j` (design.md, decisions 1 and 14). `mul`, `transpose`, `scale`, `add`
/// and `sub` of real `T` go through `@Matrix` where the compiler has it
/// (`rstd_matrix`), else through columns in its order, to the same bits
/// (decision 16).
pub fn StaticMat(comptime R: comptime_int, comptime C: comptime_int, comptime T: type) type {
    return struct {
        data: [C][R]T,

        pub const rows = R;
        pub const cols = C;
        pub const Element = T;

        /// A view of `m`'s storage, valid while `m` is.
        pub fn view(m: *@This()) View(T) {
            const flat: *[C * R]T = @ptrCast(&m.data);
            return .{ .data = flat, .rows = R, .cols = C, .ld = R };
        }

        pub fn constView(m: *const @This()) ConstView(T) {
            const flat: *const [C * R]T = @ptrCast(&m.data);
            return .{ .data = flat, .rows = R, .cols = C, .ld = R };
        }
    };
}

/// A borrowed `rows` by `cols` matrix in `data`, column-major with column
/// stride `ld` (`data[i + j * ld]`): all or part of a `Mat`, a `StaticMat`
/// or the caller's storage. `ConstView` is its read-only form.
pub fn View(comptime T: type) type {
    return ViewOf(T, false);
}

pub fn ConstView(comptime T: type) type {
    return ViewOf(T, true);
}

fn ViewOf(comptime T: type, comptime read_only: bool) type {
    return struct {
        data: if (read_only) []const T else []T,
        rows: usize,
        cols: usize,
        ld: usize,

        const Self = @This();
        pub const Element = T;

        /// `error.DimensionMismatch` for `ld < rows`, `error.SliceTooShort`
        /// when `data` ends before the last element.
        pub fn init(data: @FieldType(Self, "data"), rows: usize, cols: usize, ld: usize) error{ DimensionMismatch, SliceTooShort }!Self {
            if (ld < @max(rows, 1)) return error.DimensionMismatch;
            if (rows > 0 and cols > 0 and data.len < (cols - 1) * ld + rows) return error.SliceTooShort;
            return .{ .data = data, .rows = rows, .cols = cols, .ld = ld };
        }

        pub fn at(v: Self, i: usize, j: usize) T {
            return v.data[i + j * v.ld];
        }

        /// Rows `i..i + rows`, columns `j..j + cols`, sharing `v`'s storage.
        pub fn sub(v: Self, i: usize, j: usize, rows: usize, cols: usize) error{DimensionMismatch}!Self {
            if (i + rows > v.rows or j + cols > v.cols) return error.DimensionMismatch;
            const start = if (rows > 0 and cols > 0) i + j * v.ld else 0;
            const end = if (rows > 0 and cols > 0) start + (cols - 1) * v.ld + rows else 0;
            return .{ .data = v.data[start..end], .rows = rows, .cols = cols, .ld = v.ld };
        }

        /// Column `j`, contiguous; empty, and in bounds, with no rows.
        pub fn col(v: Self, j: usize) @FieldType(Self, "data") {
            if (v.rows == 0) return v.data[0..0];
            return v.data[j * v.ld ..][0..v.rows];
        }

        pub fn constView(v: Self) ConstView(T) {
            return .{ .data = v.data, .rows = v.rows, .cols = v.cols, .ld = v.ld };
        }
    };
}

/// A `rows` by `cols` matrix on the heap, packed column-major (`ld = rows`),
/// unmanaged as `std.ArrayList`: freed with `deinit` on the allocator
/// that made it.
pub fn Mat(comptime T: type) type {
    return struct {
        data: []T,
        rows: usize,
        cols: usize,

        const Self = @This();
        pub const Element = T;

        /// Elements undefined.
        pub fn init(allocator: std.mem.Allocator, rows: usize, cols: usize) std.mem.Allocator.Error!Self {
            return .{ .data = try allocator.alloc(T, rows * cols), .rows = rows, .cols = cols };
        }

        pub fn zeros(allocator: std.mem.Allocator, rows: usize, cols: usize) std.mem.Allocator.Error!Self {
            const m = try init(allocator, rows, cols);
            @memset(m.data, ops.zero(T));
            return m;
        }

        pub fn deinit(m: Self, allocator: std.mem.Allocator) void {
            allocator.free(m.data);
        }

        pub fn at(m: Self, i: usize, j: usize) T {
            return m.data[i + j * m.rows];
        }

        pub fn view(m: Self) View(T) {
            return .{ .data = m.data, .rows = m.rows, .cols = m.cols, .ld = @max(m.rows, 1) };
        }

        pub fn constView(m: Self) ConstView(T) {
            return .{ .data = m.data, .rows = m.rows, .cols = m.cols, .ld = @max(m.rows, 1) };
        }
    };
}

pub const M2f32 = StaticMat(2, 2, f32);
pub const M3f32 = StaticMat(3, 3, f32);
pub const M4f32 = StaticMat(4, 4, f32);
pub const M2f64 = StaticMat(2, 2, f64);
pub const M3f64 = StaticMat(3, 3, f64);
pub const M4f64 = StaticMat(4, 4, f64);

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
/// Past the inline limit, a compile error naming `rstd.linalg.mulInto`.
pub fn mul(a: anytype, b: anytype) Product(@TypeOf(a), @TypeOf(b)) {
    comptime expectInline(@TypeOf(a), @TypeOf(b), "rstd.linalg.mulInto(allocator, io, &out, &a, &b)");
    return mulAnySize(a, b);
}

fn mulAnySize(a: anytype, b: anytype) Product(@TypeOf(a), @TypeOf(b)) {
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
fn mulCols(a: anytype, b: anytype) Product(@TypeOf(a), @TypeOf(b)) {
    var r: Product(@TypeOf(a), @TypeOf(b)) = undefined;
    for (0..@TypeOf(b).cols) |j| r.data[j] = matVec(a, b.data[j]);
    return r;
}

/// `M v` without `@Matrix`, as a `Column`.
fn matVec(m: anytype, v: anytype) Column(@TypeOf(m).rows, @TypeOf(m).Element) {
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
fn vecMat(v: anytype, m: anytype) Product(@TypeOf(v), @TypeOf(m)) {
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
inline fn at(m: anytype, comptime i: comptime_int, comptime j: comptime_int) @TypeOf(m).Element {
    return m.data[j][i];
}

/// `p s - q r`.
inline fn det2(p: anytype, q: @TypeOf(p), r: @TypeOf(p), s: @TypeOf(p)) @TypeOf(p) {
    return mulAdd(@TypeOf(p), -q, r, p * s);
}

/// The 3 by 3 cofactor of row `i`, column `j`, sign included.
inline fn cofactor3(m: anytype, comptime i: comptime_int, comptime j: comptime_int) @TypeOf(m).Element {
    const ia = (i + 1) % 3;
    const ib = (i + 2) % 3;
    const ja = (j + 1) % 3;
    const jb = (j + 2) % 3;
    return det2(at(m, ia, ja), at(m, ia, jb), at(m, ib, ja), at(m, ib, jb));
}

/// Real `m` up to 3 by 3: cofactors along the first row.
inline fn closedDet(m: anytype) @TypeOf(m).Element {
    const T = @TypeOf(m).Element;
    return switch (@TypeOf(m).rows) {
        1 => at(m, 0, 0),
        2 => det2(at(m, 0, 0), at(m, 0, 1), at(m, 1, 0), at(m, 1, 1)),
        3 => mulAdd(T, at(m, 0, 2), cofactor3(m, 0, 2), mulAdd(T, at(m, 0, 1), cofactor3(m, 0, 1), at(m, 0, 0) * cofactor3(m, 0, 0))),
        else => unreachable,
    };
}

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
            for (a) |*c| std.mem.swap(T, &c[k], &c[p]);
            std.mem.swap(usize, &perm[k], &perm[p]);
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
/// singular. Past the inline limit, a compile error naming `rstd.linalg.det`.
pub fn det(m: anytype) @TypeOf(m).Element {
    comptime expectInline(@TypeOf(m), @TypeOf(m), "rstd.linalg.det(allocator, io, &m)");
    return detAnySize(m);
}

fn detAnySize(m: anytype) @TypeOf(m).Element {
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
/// naming `rstd.linalg.invInto`.
pub fn inv(m: anytype) error{Singular}!@TypeOf(m) {
    comptime expectInline(@TypeOf(m), @TypeOf(m), "rstd.linalg.invInto(allocator, io, &out, &m)");
    return invAnySize(m);
}

fn invAnySize(m: anytype) error{Singular}!@TypeOf(m) {
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
/// Past the inline limit, a compile error naming `rstd.linalg.solveInto`.
pub fn solve(a: anytype, b: anytype) error{Singular}!@TypeOf(b) {
    comptime expectInline(@TypeOf(a), @TypeOf(b), "rstd.linalg.solveInto(allocator, io, &out, &a, &b)");
    return solveAnySize(a, b);
}

fn solveAnySize(a: anytype, b: anytype) error{Singular}!@TypeOf(b) {
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

/// `mul` into `out`.
pub fn mulInto(out: anytype, a: anytype, b: anytype) void {
    out.* = mul(a, b);
}

/// `transpose` into `out`.
pub fn transposeInto(out: anytype, m: anytype) void {
    out.* = transpose(m);
}

/// `ctranspose` into `out`.
pub fn ctransposeInto(out: anytype, m: anytype) void {
    out.* = ctranspose(m);
}

/// `outer` into `out`.
pub fn outerInto(out: anytype, a: anytype, b: anytype) void {
    out.* = outer(a, b);
}

/// `scale` into `out`.
pub fn scaleInto(out: anytype, m: anytype, s: @TypeOf(m).Element) void {
    out.* = scale(m, s);
}

/// `add` into `out`.
pub fn addInto(out: anytype, a: anytype, b: @TypeOf(a)) void {
    out.* = add(a, b);
}

/// `sub` into `out`.
pub fn subInto(out: anytype, a: anytype, b: @TypeOf(a)) void {
    out.* = sub(a, b);
}

/// `inv` into `out`, untouched on `error.Singular`.
pub fn invInto(out: anytype, m: anytype) error{Singular}!void {
    out.* = try inv(m);
}

/// `solve` into `out`, untouched on `error.Singular`.
pub fn solveInto(out: anytype, a: anytype, b: anytype) error{Singular}!void {
    out.* = try solve(a, b);
}

// ---- plumbing ---------------------------------------------------------------

fn isMat(comptime M: type) bool {
    return @typeInfo(M) == .@"struct" and @hasDecl(M, "Element") and @hasDecl(M, "rows") and
        @hasDecl(M, "cols") and M == StaticMat(M.rows, M.cols, M.Element);
}

/// A length-`n` vector of `T`: a `Vec` for real `T`, an array for complex.
fn Column(comptime n: comptime_int, comptime T: type) type {
    return if (isComplex(T)) [n]T else Vec(n, T);
}

fn expectSquare(comptime M: type) void {
    if (!isMat(M)) @compileError("expects a StaticMat, found " ++ @typeName(M));
    if (M.rows != M.cols) @compileError("expects a square matrix, found " ++ @typeName(M));
}

/// Whether `T` computes through `@Matrix`, which holds real floats only.
fn native(comptime T: type) bool {
    return rstd_matrix.has_matrix and !isComplex(T);
}

inline fn toNative(m: anytype) rstd_matrix.Matrix(@TypeOf(m).rows, @TypeOf(m).cols, @TypeOf(m).Element) {
    return m.data;
}

/// Whether an `R` by `C` matrix of `T` is a SPIR-V `OpTypeMatrix` (through
/// `@Matrix`, 2 to 4 by 2 to 4, real), which takes SPIR-V's own matrix
/// instructions (design.md, "Dispatch under SPIR-V").
fn spirvMat(comptime R: comptime_int, comptime C: comptime_int, comptime T: type) bool {
    return vector.on_spirv and native(T) and R >= 2 and R <= 4 and C >= 2 and C <= 4;
}

// GLSL.std.450 and core matrix instructions by inline asm, on `@Matrix`
// operands; literal templates, as `vector.zig`'s.

inline fn spvDeterminant(m: anytype) @TypeOf(m).Element {
    const T = @TypeOf(m).Element;
    return asm (
        \\%set = OpExtInstImport "GLSL.std.450"
        \\%ret = OpExtInst %T %set 33 %m
        : [ret] "" (-> T),
        : [T] "t" (T),
          [m] "" (toNative(m)),
    );
}

inline fn spvInverse(m: anytype) rstd_matrix.Matrix(@TypeOf(m).rows, @TypeOf(m).cols, @TypeOf(m).Element) {
    const N = rstd_matrix.Matrix(@TypeOf(m).rows, @TypeOf(m).cols, @TypeOf(m).Element);
    return asm (
        \\%set = OpExtInstImport "GLSL.std.450"
        \\%ret = OpExtInst %T %set 34 %m
        : [ret] "" (-> N),
        : [T] "t" (N),
          [m] "" (toNative(m)),
    );
}

inline fn spvOuter(comptime R: comptime_int, comptime C: comptime_int, comptime T: type, a: Vec(R, T), b: Vec(C, T)) rstd_matrix.Matrix(R, C, T) {
    const N = rstd_matrix.Matrix(R, C, T);
    return asm (
        \\%ret = OpOuterProduct %T %a %b
        : [ret] "" (-> N),
        : [T] "t" (N),
          [a] "" (a),
          [b] "" (b),
    );
}

// Element arithmetic over real and complex `T`. A sum of products rounds its
// first product and fuses the rest as `mulAdd` says; a complex product adds
// `re re` then `-im im`, and `re im` then `im re`.

pub inline fn eMul(x: anytype, y: @TypeOf(x)) @TypeOf(x) {
    const T = @TypeOf(x);
    if (comptime !isComplex(T)) return x * y;
    const R = Real(T);
    return .{
        .re = mulAdd(R, -x.im, y.im, x.re * y.re),
        .im = mulAdd(R, x.im, y.re, x.re * y.im),
    };
}

/// `acc + x y`.
pub inline fn eMulAdd(x: anytype, y: @TypeOf(x), acc: @TypeOf(x)) @TypeOf(x) {
    const T = @TypeOf(x);
    if (comptime !isComplex(T)) return mulAdd(T, x, y, acc);
    const R = Real(T);
    return .{
        .re = mulAdd(R, -x.im, y.im, mulAdd(R, x.re, y.re, acc.re)),
        .im = mulAdd(R, x.im, y.re, mulAdd(R, x.re, y.im, acc.im)),
    };
}

inline fn eDiv(x: anytype, y: @TypeOf(x)) @TypeOf(x) {
    return if (comptime isComplex(@TypeOf(x))) x.div(y) else x / y;
}

// ---- tests ------------------------------------------------------------------

const zigstd = @import("../zigstd.zig");
const V2f32 = vector.V2f32;
const V3f32 = vector.V3f32;
const V4f64 = vector.V4f64;

const bitwise = @import("blas/_test/bitwise.zig");
const expectSameBits = bitwise.expectSameBits;
const math = @import("../math.zig");

test "transposeVecs is pffft's VTRANSPOSE4 (_MM_TRANSPOSE4_PS)" {
    const V = vector.V4f32;
    const t = transposeVecs(4, f32, .{ V{ 0, 1, 2, 3 }, V{ 4, 5, 6, 7 }, V{ 8, 9, 10, 11 }, V{ 12, 13, 14, 15 } });
    try zigstd.testExpectEqual(V{ 0, 4, 8, 12 }, t[0]);
    try zigstd.testExpectEqual(V{ 3, 7, 11, 15 }, t[3]);
    const d = transposeVecs(2, f64, .{ .{ 1, 2 }, .{ 3, 4 } });
    try zigstd.testExpectEqual(vector.V2f64{ 1, 3 }, d[0]);
    try zigstd.testExpectEqual(vector.V2f64{ 2, 4 }, d[1]);
}

test "StaticMat construction, access and transposes" {
    const M = StaticMat(2, 3, f32);
    const m = fromRows(M, .{ .{ 1, 3, 5 }, .{ 2, 4, 6 } });
    try zigstd.testExpectEqual(fromCols(M, .{ .{ 1, 2 }, .{ 3, 4 }, .{ 5, 6 } }), m);
    try zigstd.testExpectEqual(V3f32{ 2, 4, 6 }, row(m, 1));
    try zigstd.testExpectEqual(V2f32{ 5, 6 }, col(m, 2));
    try zigstd.testExpectEqual(fromRows(StaticMat(3, 2, f32), .{ .{ 1, 2 }, .{ 3, 4 }, .{ 5, 6 } }), transpose(m));
    try zigstd.testExpectEqual(m, transpose(transpose(m)));
    try zigstd.testExpectEqual([3][2]f32{ .{ 1, 0 }, .{ 0, 1 }, .{ 0, 0 } }, identity(M).data);
    try zigstd.testExpectEqual(@as(f32, 4), trace(identity(M4f32)));

    const C = std.math.Complex(f64);
    const z: StaticMat(1, 2, C) = .{ .data = .{ .{.init(1, 2)}, .{.init(3, -4)} } };
    const zh = ctranspose(z);
    try zigstd.testExpectEqual(C.init(1, -2), zh.data[0][0]);
    try zigstd.testExpectEqual(C.init(3, 4), zh.data[0][1]);
    try zigstd.testExpectEqual(C.init(3, -4), transpose(z).data[0][1]);
}

test "views: init's checks, sub, col; StaticMat and Mat views share storage" {
    var buf = [_]f64{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    try std.testing.expectError(error.DimensionMismatch, View(f64).init(&buf, 3, 2, 2));
    try std.testing.expectError(error.SliceTooShort, View(f64).init(&buf, 3, 4, 3));
    const v = try View(f64).init(&buf, 3, 3, 3);
    try std.testing.expectEqual(@as(f64, 8), v.at(1, 2));
    const s = try v.sub(1, 1, 2, 2);
    try std.testing.expectEqual(@as(f64, 5), s.at(0, 0));
    try std.testing.expectEqual(@as(f64, 9), s.at(1, 1));
    try std.testing.expectEqualSlices(f64, &.{ 5, 6 }, s.col(0));
    try std.testing.expectError(error.DimensionMismatch, v.sub(2, 0, 2, 1));

    var sm = fromCols(M2f64, .{ .{ 1, 2 }, .{ 3, 4 } });
    sm.view().data[3] = 40;
    try std.testing.expectEqual(@as(f64, 40), sm.data[1][1]);
    try std.testing.expectEqual(@as(f64, 3), sm.constView().at(0, 1));
    const m = try Mat(f64).zeros(std.testing.allocator, 2, 3);
    defer m.deinit(std.testing.allocator);
    m.view().data[5] = 7;
    try std.testing.expectEqual(@as(f64, 7), m.at(1, 2));
}

test "mul, scale, add, sub and outer on exact values" {
    // [1 3 5; 2 4 6] and [7 10; 8 11; 9 12], as the fork's behavior test
    const a = fromCols(StaticMat(2, 3, f32), .{ .{ 1, 2 }, .{ 3, 4 }, .{ 5, 6 } });
    const b = fromCols(StaticMat(3, 2, f32), .{ .{ 7, 8, 9 }, .{ 10, 11, 12 } });
    try zigstd.testExpectEqual([2][2]f32{ .{ 76, 100 }, .{ 103, 136 } }, mul(a, b).data);
    try zigstd.testExpectEqual(V2f32{ 9, 12 }, mul(a, V3f32{ 1, 1, 1 }));
    try zigstd.testExpectEqual([2]f32{ 9, 12 }, mul(a, [3]f32{ 1, 1, 1 }));
    try zigstd.testExpectEqual(V3f32{ 5, 11, 17 }, mul(V2f32{ 1, 2 }, a));
    try zigstd.testExpectEqual(add(a, a), scale(a, 2));
    try zigstd.testExpectEqual(StaticMat(2, 3, f32){ .data = @splat(@splat(0)) }, sub(a, a));
    try zigstd.testExpectEqual([3][2]f32{ .{ 1, 2 }, .{ 2, 4 }, .{ 3, 6 } }, outer(V2f32{ 1, 2 }, [3]f32{ 1, 2, 3 }).data);
    var out: M2f32 = undefined;
    mulInto(&out, a, b);
    try zigstd.testExpectEqual(mul(a, b), out);

    const C = std.math.Complex(f64);
    const za: StaticMat(2, 2, C) = .{ .data = .{ .{ .init(1, 1), .init(0, 2) }, .{ .init(3, 0), .init(1, -1) } } };
    const zv = [2]C{ .init(2, 0), .init(0, 1) };
    // [1+i 3; 2i 1-i] [2; i] = [2+5i; 1+5i]
    const zr = mul(za, zv);
    try zigstd.testExpectEqual(C.init(2, 5), zr[0]);
    try zigstd.testExpectEqual(C.init(1, 5), zr[1]);
    // [2 i] [1+i 3; 2i 1-i] = [0+2i, 7+i]
    const zl = mul(zv, za);
    try zigstd.testExpectEqual(C.init(0, 2), zl[0]);
    try zigstd.testExpectEqual(C.init(7, 1), zl[1]);
    try zigstd.testExpectEqual(mul(za, za).data[1], mul(za, za.data[1]));
    // a bᴴ conjugates b
    try zigstd.testExpectEqual(C.init(0, -2), outer([1]C{.init(2, 0)}, [1]C{.init(0, 1)}).data[0][0]);
}

fn randomMat(comptime M: type, random: std.Random, s: M.Element) M {
    var m: M = undefined;
    for (&m.data) |*c| for (c) |*x| {
        x.* = (random.float(M.Element) * 2 - 1) * s;
    };
    return m;
}

fn expectSameMatBits(want: anytype, got: @TypeOf(want)) !void {
    for (want.data, got.data) |wc, gc| try expectSameOrX87(@as(Vec(@TypeOf(want).rows, @TypeOf(want).Element), wc), gc);
}

/// `expectSameBits` of two `@Vector`s, but for `f64` on x87 within
/// `bitwise.sameBitsOrX87`'s few ulps: there each code path rounds its 80-bit
/// intermediates where it stores them (issue #264).
fn expectSameOrX87(want: anytype, got: @TypeOf(want)) !void {
    const V = @typeInfo(@TypeOf(want)).vector;
    if (comptime !(math.x87_f64 and V.child == f64)) return expectSameBits(want, got);
    const wa: [V.len]f64 = want;
    const ga: [V.len]f64 = got;
    try zigstd.testExpectTrue(bitwise.sameBitsOrX87(f64, &wa, &ga));
}

/// `mul`'s `@Matrix` forms against its column code, for `R` by `K` times `K`
/// by `C`.
fn expectNativeBits(comptime T: type, comptime r: comptime_int, comptime k: comptime_int, comptime c: comptime_int, random: std.Random) !void {
    const scales = [_]f64{ 1, 1e-3, 1e3, 1e30, 1e-30 };
    for (0..50) |trial| {
        const s: T = @floatCast(scales[trial % scales.len]);
        const a = randomMat(StaticMat(r, k, T), random, s);
        const b = randomMat(StaticMat(k, c, T), random, s);
        const v: Vec(k, T) = b.data[0];
        const w: Vec(r, T) = a.data[0];
        try expectSameMatBits(mulCols(a, b), mul(a, b));
        try expectSameOrX87(@as(Vec(r, T), matVec(a, v)), mul(a, v));
        try expectSameOrX87(@as(Vec(k, T), vecMat(w, a)), mul(w, a));
        try expectSameMatBits(transpose(transpose(a)), a);
    }
}

test "mul through @Matrix has its column code's bits" {
    if (!rstd_matrix.has_matrix) return error.SkipZigTest;
    var prng: std.Random.DefaultPrng = .init(0x3a7);
    const shapes = .{ .{ 1, 1, 1 }, .{ 2, 2, 2 }, .{ 3, 3, 3 }, .{ 4, 4, 4 }, .{ 2, 3, 4 }, .{ 5, 7, 3 }, .{ 6, 6, 1 }, .{ 8, 8, 8 } };
    inline for (.{ f32, f64 }) |T| inline for (shapes) |s| try expectNativeBits(T, s[0], s[1], s[2], prng.random());
}

test "det, inv and solve on known matrices" {
    const m2 = fromRows(M2f64, .{ .{ 4, 7 }, .{ 2, 6 } });
    try zigstd.testExpectEqual(@as(f64, 10), det(m2));
    try zigstd.testExpectEqual(fromRows(M2f64, .{ .{ 0.6, -0.7 }, .{ -0.2, 0.4 } }), try inv(m2));
    const m3 = fromRows(M3f64, .{ .{ 1, 2, 3 }, .{ 0, 1, 4 }, .{ 5, 6, 0 } });
    try zigstd.testExpectEqual(@as(f64, 1), det(m3));
    try zigstd.testExpectEqual(fromRows(M3f64, .{ .{ -24, 18, 5 }, .{ 20, -15, -4 }, .{ -5, 4, 1 } }), try inv(m3)); // a permutation needs a row swap; its determinant is the swaps' sign
    const p4 = fromCols(M4f64, .{ .{ 0, 1, 0, 0 }, .{ 1, 0, 0, 0 }, .{ 0, 0, 0, 1 }, .{ 0, 0, 1, 0 } });
    try zigstd.testExpectEqual(@as(f64, 1), det(p4));
    try zigstd.testExpectEqual(transpose(p4), try inv(p4));
    try zigstd.testExpectEqual(V4f64{ 2, 1, 4, 3 }, try solve(p4, V4f64{ 1, 2, 3, 4 }));
    try zigstd.testExpectEqual(@as(f64, -2), det(fromRows(StaticMat(5, 5, f64), .{
        .{ 0, 2, 0, 0, 0 }, .{ 1, 0, 0, 0, 0 }, .{ 0, 0, 1, 0, 0 }, .{ 0, 0, 0, 1, 0 }, .{ 0, 0, 0, 0, 1 },
    })));
    try zigstd.testExpectEqual(@as(f32, -3), det(StaticMat(1, 1, f32){ .data = .{.{-3}} }));
    try zigstd.testExpectEqual(@as(f32, -0.5), (try inv(StaticMat(1, 1, f32){ .data = .{.{-2}} })).data[0][0]);

    const C = std.math.Complex(f64);
    const z: StaticMat(2, 2, C) = .{ .data = .{ .{ .init(0, 1), .init(0, 0) }, .{ .init(0, 0), .init(2, 0) } } };
    try zigstd.testExpectEqual(C.init(0, 2), det(z));
    try zigstd.testExpectEqual(C.init(0, -1), (try inv(z)).data[0][0]);
}

test "det, inv and solve on singular matrices" {
    try zigstd.testExpectEqual(@as(f64, 0), det(fromRows(M2f64, .{ .{ 1, 2 }, .{ 2, 4 } })));
    try zigstd.testExpectError(error.Singular, inv(fromRows(M2f64, .{ .{ 1, 2 }, .{ 2, 4 } })));
    try zigstd.testExpectError(error.Singular, inv(fromRows(M3f32, .{ .{ 1, 2, 3 }, .{ 4, 5, 6 }, .{ 7, 8, 9 } })));
    try zigstd.testExpectError(error.Singular, solve(fromRows(M3f32, .{ .{ 1, 2, 3 }, .{ 0, 0, 0 }, .{ 7, 8, 9 } }), V3f32{ 1, 2, 3 }));
    const s5 = StaticMat(5, 5, f64){ .data = @splat(@splat(1)) };
    try zigstd.testExpectEqual(@as(f64, 0), det(s5));
    try zigstd.testExpectError(error.Singular, inv(s5));
    var out = identity(M4f64);
    try zigstd.testExpectError(error.Singular, invInto(&out, M4f64{ .data = @splat(@splat(0)) }));
    try zigstd.testExpectEqual(identity(M4f64), out);
    const C = std.math.Complex(f32);
    try zigstd.testExpectError(error.Singular, inv(StaticMat(2, 2, C){ .data = .{ .{ .init(1, 1), .init(2, 2) }, .{ .init(1, 1), .init(2, 2) } } }));
}

fn isClose(x: anytype, y: @TypeOf(x), tol: Real(@TypeOf(x))) bool {
    return ops.abs1(@TypeOf(x), ops.sub(@TypeOf(x), x, y)) <= tol;
}

/// `inv(m) m` and `m inv(m)` are `I` and `m solve(m, b)` is `b` within `tol`,
/// for random diagonally dominant `n` by `n` matrices.
fn expectInverts(comptime T: type, comptime n: comptime_int, random: std.Random) !void {
    const M = StaticMat(n, n, T);
    const R = Real(T);
    const tol: R = 1e-12;
    for (0..20) |_| {
        var m: M = undefined;
        var b: [n]T = undefined;
        for (0..n) |j| {
            for (0..n) |i| m.data[j][i] = randomElem(T, random);
            m.data[j][j] = ops.add(T, m.data[j][j], if (comptime isComplex(T)) T.init(2 * n, 0) else 2 * n);
            b[j] = randomElem(T, random);
        }
        const mi = try inv(m);
        const eye = identity(M);
        for (mul(mi, m).data, mul(m, mi).data, eye.data) |p, q, e| for (p, q, e) |x, y, z| {
            try zigstd.testExpectTrue(isClose(x, z, tol) and isClose(y, z, tol));
        };
        const x = try solve(m, b);
        for (mul(m, x), b) |got, want| try zigstd.testExpectTrue(isClose(got, want, tol));
        const xs = try solve(m, fromCols(StaticMat(n, 1, T), .{b}));
        try zigstd.testExpectTrue(bitwise.sameBitsOrX87(T, &x, &xs.data[0]));
        // the determinant is the product of the eigenvalues: of `mᵀ` the same
        const d = det(m);
        const dt = det(transpose(m));
        try zigstd.testExpectTrue(isClose(d, dt, tol * ops.abs1(T, d)));
    }
}

fn randomElem(comptime T: type, random: std.Random) T {
    if (comptime isComplex(T)) return .init(random.float(Real(T)) * 2 - 1, random.float(Real(T)) * 2 - 1);
    return random.float(T) * 2 - 1;
}

test "inv and solve invert, real and complex, n 1 to 8" {
    var prng: std.Random.DefaultPrng = .init(0x1a);
    inline for (.{ f64, std.math.Complex(f64) }) |T| inline for (1..9) |n| try expectInverts(T, n, prng.random());
}

test {
    zigstd.testRefAllDecls(@This());
}
