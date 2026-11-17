//! Plumbing: dispatch by size, operand views, copies.

const sysinfo = @import("builtin");
const std = @import("std");
const Allocator = std.mem.Allocator;
const Io = std.Io;
const blas = @import("../blas/__root__.zig");
const ext = @import("../blas/ext/interface.zig");
const l3 = @import("../blas/blas3/core.zig");
const l2 = @import("../blas/blas2/core.zig");
const l1 = @import("../blas/blas1/interface.zig");
const lapack = @import("../lapack/__root__.zig");
const vector = @import("../vector.zig");
const matrix = @import("../matrix.zig");
const ops = @import("../lapack/ops.zig");

const View = matrix.View;
const ConstView = matrix.ConstView;
const Mat = matrix.Mat;
const isComplex = vector.isComplex;
const Real = vector.Real;
const eMul = matrix.eMul;
const eMulAdd = matrix.eMulAdd;

const inverse = @import("inverse.zig");
const det = inverse.det;
const inv = inverse.inv;
const solve = inverse.solve;
const reductions = @import("reductions.zig");
const dot = reductions.dot;
const extreme = reductions.extreme;
const sum = reductions.sum;

/// The inline limit `L` per dimension, `rstd_linalg_options.inline_limit`
/// (design.md, Dispatch).
pub const inline_limit = matrix.inline_limit;

/// The largest inline size built: `L`, and 16 in tests, which run each path
/// at every size 1 to 16 (design.md, Dispatch).
pub const fixed_max = if (sysinfo.is_test) @max(inline_limit, 16) else inline_limit;

/// The inline kernels or the BLAS and LAPACK; by size alone outside tests.
pub const Path = enum { fixed, blas };

/// `.fixed` when every dimension is at or below the inline limit.
pub fn pathOf(dims: []const usize) Path {
    for (dims) |d| if (d > inline_limit) return .blas;
    return .fixed;
}

/// The element type of an operand.
pub fn ElementOf(comptime X: type) type {
    switch (@typeInfo(X)) {
        .pointer => |p| switch (p.size) {
            .slice => return p.child,
            .one => return if (@typeInfo(p.child) == .array) @typeInfo(p.child).array.child else p.child.Element,
            else => {},
        },
        .@"struct" => if (@hasField(X, "data") and @typeInfo(@FieldType(X, "data")) == .pointer) return X.Element,
        else => {},
    }
    @compileError("expects a View, ConstView, Mat, pointer to a StaticMat or slice, found " ++ @typeName(X));
}

/// `x` read-only.
pub fn constOf(x: anytype) ConstView(ElementOf(@TypeOf(x))) {
    const T = ElementOf(@TypeOf(x));
    return switch (@typeInfo(@TypeOf(x))) {
        .pointer => |p| if (p.size == .slice)
            column(T, true, x)
        else if (@typeInfo(p.child) == .array) column(T, true, x[0..]) else x.constView(),
        else => x.constView(),
    };
}

/// `x` writable.
pub fn viewOf(x: anytype) View(ElementOf(@TypeOf(x))) {
    const X = @TypeOf(x);
    const T = ElementOf(X);
    if (X == View(T)) return x;
    if (X == ConstView(T)) @compileError("expects a writable operand, found a ConstView");
    return switch (@typeInfo(X)) {
        .pointer => |p| if (p.size == .slice)
            column(T, false, x)
        else if (@typeInfo(p.child) == .array) column(T, false, x[0..]) else x.view(),
        else => x.view(),
    };
}

pub fn column(comptime T: type, comptime read_only: bool, s: if (read_only) []const T else []T) if (read_only) ConstView(T) else View(T) {
    return .{ .data = s, .rows = s.len, .cols = 1, .ld = @max(s.len, 1) };
}

/// A vector's elements: `n` of them `inc` apart from `ptr`.
pub fn Strided(comptime T: type) type {
    return struct {
        ptr: [*]T,
        n: usize,
        inc: usize,

        fn at(s: @This(), i: usize) T {
            return s.ptr[i * s.inc];
        }

        /// The first `c` elements as an array, for the fixed-size forms.
        fn gather(s: @This(), comptime c: usize) [c]T {
            var r: [c]T = undefined;
            for (&r, 0..) |*x, i| x.* = s.at(i);
            return r;
        }
    };
}

/// `linalg.vector`'s forms at length `c` and `linalg.matrix`'s at `c` by `c`,
/// one non-inline function per size so that each gets its own comptime
/// branch quota.
pub fn Fixed(comptime c: usize) type {
    return struct {
        fn load(a: anytype) matrix.StaticMat(c, c, @TypeOf(a).Element) {
            var m: matrix.StaticMat(c, c, @TypeOf(a).Element) = undefined;
            for (&m.data, 0..) |*mc, j| @memcpy(mc, a.col(j));
            return m;
        }

        fn det(a: anytype) @TypeOf(a).Element {
            return matrix.any_size.det(load(a));
        }

        fn inv(o: anytype, a: ConstView(@TypeOf(o).Element)) error{Singular}!void {
            const r = try matrix.any_size.inv(load(a));
            for (&r.data, 0..) |*rc, j| @memcpy(o.col(j), rc);
        }

        /// `matrix.solve`'s LU, factored once for all of `b`'s columns.
        fn solve(o: anytype, a: ConstView(@TypeOf(o).Element), b: ConstView(@TypeOf(o).Element)) error{Singular}!void {
            const T = @TypeOf(o).Element;
            var f = load(a).data;
            var perm: [c]usize = undefined;
            var odd: bool = undefined;
            try matrix.luInPlace(c, T, &f, &perm, &odd);
            for (0..b.cols) |j| {
                const x = matrix.luSolve(c, T, &f, &perm, b.col(j)[0..c].*);
                @memcpy(o.col(j), &x);
            }
        }

        fn dot(x: anytype, y: @TypeOf(x)) @TypeOf(x.ptr[0]) {
            return vector.dot(x.gather(c), y.gather(c));
        }

        fn length(x: anytype) Real(@TypeOf(x.ptr[0])) {
            return vector.length(x.gather(c));
        }

        fn sum(x: anytype) @TypeOf(x.ptr[0]) {
            return vector.sum(x.gather(c));
        }

        fn extreme(comptime which: blas.Extreme, x: anytype) @TypeOf(x.ptr[0]) {
            return if (which == .max) vector.max(x.gather(c)) else vector.min(x.gather(c));
        }
    };
}

/// `Fixed(n).name(args)`, `n` from 1 to `fixed_max`.
pub inline fn fixed(comptime R: type, comptime name: []const u8, n: usize, args: anytype) R {
    switch (n) {
        inline 1...fixed_max => |c| return @call(.auto, @field(Fixed(c), name), args),
        else => unreachable,
    }
}

/// `v` as a vector: one column, or one row (stride `ld`).
pub fn strided(v: anytype) error{DimensionMismatch}!Strided(@TypeOf(v).Element) {
    const p: [*]@TypeOf(v).Element = @constCast(v.data.ptr);
    if (v.cols == 1 or v.rows == 0 or v.cols == 0) return .{ .ptr = p, .n = if (v.cols == 0) 0 else v.rows, .inc = 1 };
    if (v.rows == 1) return .{ .ptr = p, .n = v.cols, .inc = v.ld };
    return error.DimensionMismatch;
}

pub fn toIsize(n: usize) isize {
    return @intCast(n);
}

/// What the BLAS and LAPACK calls here pass on (`ops.narrow`): the arguments
/// come from views, so their argument errors are bugs.
pub const Error = Allocator.Error || Io.Cancelable;
pub const narrow = ops.narrow;

/// `a` into `dst`, column stride `ld`.
pub fn copyInto(comptime T: type, dst: []T, ld: usize, a: ConstView(T)) void {
    for (0..a.cols) |j| @memcpy(dst[j * ld ..][0..a.rows], a.col(j));
}

/// A packed copy of `a`, for LAPACK to overwrite.
pub fn copyOf(comptime T: type, allocator: Allocator, a: ConstView(T)) Allocator.Error!Mat(T) {
    const r = try Mat(T).init(allocator, a.rows, a.cols);
    copyInto(T, r.data, @max(a.rows, 1), a);
    return r;
}
