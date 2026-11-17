//! Building matrices and taking their structure apart, by MATLAB's rules:
//! `identity` (`eye`), `ones`, `tril`, `triu`, `diagMatrix` and `diagonal`
//! (the two ways of `diag`), `kron`. Diagonal `k` is above the main one for
//! `k > 0`, below for `k < 0`. Replaces `matlab/generators.zig`'s (issue
//! #325).

const std = @import("std");
const matrix = @import("matrix.zig");
const vector = @import("vector.zig");

const Allocator = std.mem.Allocator;
const Mat = matrix.Mat;
const ConstView = matrix.ConstView;
const isComplex = vector.isComplex;
const Real = vector.Real;

/// Ones on the main diagonal, rectangular as MATLAB's `eye(m, n)`.
pub fn identity(allocator: Allocator, comptime T: type, rows: usize, cols: usize) Allocator.Error!Mat(T) {
    const m = try Mat(T).zeros(allocator, rows, cols);
    for (0..@min(rows, cols)) |i| m.data[i * (rows + 1)] = scalar(T, 1);
    return m;
}

pub fn ones(allocator: Allocator, comptime T: type, rows: usize, cols: usize) Allocator.Error!Mat(T) {
    const m = try Mat(T).init(allocator, rows, cols);
    @memset(m.data, scalar(T, 1));
    return m;
}

/// `a` on and below diagonal `k`, zero above.
pub fn tril(allocator: Allocator, comptime T: type, a: ConstView(T), k: isize) Allocator.Error!Mat(T) {
    return keep(allocator, T, a, k, .lower);
}

/// `a` on and above diagonal `k`, zero below.
pub fn triu(allocator: Allocator, comptime T: type, a: ConstView(T), k: isize) Allocator.Error!Mat(T) {
    return keep(allocator, T, a, k, .upper);
}

fn keep(allocator: Allocator, comptime T: type, a: ConstView(T), k: isize, comptime side: enum { lower, upper }) Allocator.Error!Mat(T) {
    const m = try Mat(T).init(allocator, a.rows, a.cols);
    for (0..a.cols) |j| for (0..a.rows) |i| {
        const d = @as(isize, @intCast(j)) - @as(isize, @intCast(i));
        const in = switch (side) {
            .lower => d <= k,
            .upper => d >= k,
        };
        m.data[i + j * a.rows] = if (in) a.at(i, j) else scalar(T, 0);
    };
    return m;
}

/// `v` on diagonal `k` of a square matrix of order `v.len + |k|`, zero
/// elsewhere (`diag(v, k)`).
pub fn diagMatrix(allocator: Allocator, comptime T: type, v: []const T, k: isize) Allocator.Error!Mat(T) {
    const off = @abs(k);
    const n = v.len + off;
    const m = try Mat(T).zeros(allocator, n, n);
    for (v, 0..) |x, t| {
        const i = if (k < 0) t + off else t;
        const j = if (k > 0) t + off else t;
        m.data[i + j * n] = x;
    }
    return m;
}

/// Diagonal `k` of `a` as a new slice (`diag(A, k)`); empty past the corner.
pub fn diagonal(allocator: Allocator, comptime T: type, a: ConstView(T), k: isize) Allocator.Error![]T {
    const off = @abs(k);
    const row0: usize = if (k < 0) off else 0;
    const col0: usize = if (k > 0) off else 0;
    const len = if (row0 >= a.rows or col0 >= a.cols) 0 else @min(a.rows - row0, a.cols - col0);
    const d = try allocator.alloc(T, len);
    for (d, 0..) |*x, t| x.* = a.at(row0 + t, col0 + t);
    return d;
}

/// The Kronecker product: block `(i, j)` is `a(i, j) b`.
pub fn kron(allocator: Allocator, comptime T: type, a: ConstView(T), b: ConstView(T)) Allocator.Error!Mat(T) {
    const rows = a.rows * b.rows;
    const m = try Mat(T).init(allocator, rows, a.cols * b.cols);
    for (0..a.cols) |ja| for (0..b.cols) |jb| {
        const col = m.data[(ja * b.cols + jb) * rows ..][0..rows];
        for (0..a.rows) |ia| {
            const x = a.at(ia, ja);
            for (col[ia * b.rows ..][0..b.rows], b.col(jb)) |*o, y| o.* = mul(T, x, y);
        }
    };
    return m;
}

fn scalar(comptime T: type, x: Real(T)) T {
    return if (comptime isComplex(T)) T.init(x, 0) else x;
}

fn mul(comptime T: type, x: T, y: T) T {
    return if (comptime isComplex(T)) x.mul(y) else x * y;
}

// ---- tests ----

const testing = std.testing;

test "identity, ones, tril, triu" {
    const al = testing.allocator;
    const e = try identity(al, f64, 2, 3);
    defer e.deinit(al);
    try testing.expectEqualSlices(f64, &.{ 1, 0, 0, 1, 0, 0 }, e.data);
    const o = try ones(al, std.math.Complex(f32), 1, 2);
    defer o.deinit(al);
    try testing.expectEqual(@as(f32, 1), o.data[1].re);
    // [1 2 3; 4 5 6; 7 8 9]
    const a = [_]f64{ 1, 4, 7, 2, 5, 8, 3, 6, 9 };
    const v: ConstView(f64) = try .init(&a, 3, 3, 3);
    const l = try tril(al, f64, v, -1);
    defer l.deinit(al);
    try testing.expectEqualSlices(f64, &.{ 0, 4, 7, 0, 0, 8, 0, 0, 0 }, l.data);
    const u = try triu(al, f64, v, 1);
    defer u.deinit(al);
    try testing.expectEqualSlices(f64, &.{ 0, 0, 0, 2, 0, 0, 3, 6, 0 }, u.data);
}

test "diag both ways, kron" {
    const al = testing.allocator;
    const m = try diagMatrix(al, f64, &.{ 1, 2 }, 1);
    defer m.deinit(al);
    try testing.expectEqualSlices(f64, &.{ 0, 0, 0, 1, 0, 0, 0, 2, 0 }, m.data);
    const a = [_]f64{ 1, 4, 7, 2, 5, 8, 3, 6, 9 };
    const v: ConstView(f64) = try .init(&a, 3, 3, 3);
    const d = try diagonal(al, f64, v, -1);
    defer al.free(d);
    try testing.expectEqualSlices(f64, &.{ 4, 8 }, d);
    const far = try diagonal(al, f64, v, 5);
    defer al.free(far);
    try testing.expectEqual(0, far.len);
    // kron([1 2], [1; 1]) = [1 2; 1 2]
    const k = try kron(al, f64, try .init(&[_]f64{ 1, 2 }, 1, 2, 1), try .init(&[_]f64{ 1, 1 }, 2, 1, 2));
    defer k.deinit(al);
    try testing.expectEqualSlices(f64, &.{ 1, 1, 2, 2 }, k.data);
}
