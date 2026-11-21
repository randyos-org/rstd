//! Reduced row echelon form (MATLAB's `rref.m`): Gauss-Jordan elimination
//! with partial pivoting and a tolerance. Replaces matlab88's `RREF`
//! (`matlab/generators.zig`, issue #329).

const std = @import("std");
const mem = @import("../../mem.zig");
const matrix = @import("matrix.zig");
const vector = @import("vector.zig");

const Allocator = std.mem.Allocator;
const Mat = matrix.Mat;
const ConstView = matrix.ConstView;
const isComplex = vector.isComplex;
const Real = vector.Real;

pub fn Rref(comptime T: type) type {
    return struct {
        r: Mat(T),
        /// The pivot columns, increasing (MATLAB's `jb`); their count is
        /// `A`'s rank to the tolerance.
        pivots: []usize,

        pub fn deinit(self: @This(), allocator: Allocator) void {
            self.r.deinit(allocator);
            allocator.free(self.pivots);
        }
    };
}

/// `a`'s reduced row echelon form. Entries at most `tol` count as zero;
/// null is MATLAB's `max(m, n) eps ‖A‖∞`.
pub fn rref(allocator: Allocator, comptime T: type, a: ConstView(T), tol: ?Real(T)) Allocator.Error!Rref(T) {
    const F = Real(T);
    const m = a.rows;
    const n = a.cols;
    const r = try Mat(T).init(allocator, m, n);
    errdefer r.deinit(allocator);
    for (0..n) |j| @memcpy(r.data[j * m ..][0..m], a.col(j));
    var pivots: std.ArrayList(usize) = .empty;
    errdefer pivots.deinit(allocator);
    try pivots.ensureTotalCapacity(allocator, @min(m, n));

    const t = tol orelse blk: {
        var inf_norm: F = 0;
        for (0..m) |i| {
            var s: F = 0;
            for (0..n) |j| s += abs(T, r.data[i + j * m]);
            inf_norm = @max(inf_norm, s);
        }
        break :blk @as(F, @floatFromInt(@max(m, n))) * std.math.floatEps(F) * inf_norm;
    };

    const at = struct {
        fn f(x: Mat(T), i: usize, j: usize) *T {
            return &x.data[i + j * x.rows];
        }
    }.f;
    var i: usize = 0;
    var j: usize = 0;
    while (i < m and j < n) {
        // the largest in column j from row i, the first on ties
        var k = i;
        var p: F = abs(T, at(r, i, j).*);
        for (i + 1..m) |q| {
            const v = abs(T, at(r, q, j).*);
            if (v > p) {
                p = v;
                k = q;
            }
        }
        if (p <= t) {
            for (i..m) |q| at(r, q, j).* = scalar(T, 0);
            j += 1;
            continue;
        }
        pivots.appendAssumeCapacity(j);
        for (j..n) |c| mem.swap(T, at(r, i, c), at(r, k, c));
        const d = at(r, i, j).*;
        for (j..n) |c| at(r, i, c).* = div(T, at(r, i, c).*, d);
        for (0..m) |q| {
            if (q == i) continue;
            const f = at(r, q, j).*;
            for (j..n) |c| at(r, q, c).* = sub(T, at(r, q, c).*, mul(T, f, at(r, i, c).*));
        }
        i += 1;
        j += 1;
    }
    return .{ .r = r, .pivots = try pivots.toOwnedSlice(allocator) };
}

fn abs(comptime T: type, x: T) Real(T) {
    return if (comptime isComplex(T)) std.math.hypot(x.re, x.im) else @abs(x);
}

fn scalar(comptime T: type, x: Real(T)) T {
    return if (comptime isComplex(T)) T.init(x, 0) else x;
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

test "rref: MATLAB's examples, rank-deficient and wide" {
    const al = testing.allocator;
    // magic(4): rank 3; MATLAB's rref gives [1 0 0 1; 0 1 0 3; 0 0 1 -3; 0 0 0 0]
    const m4 = [_]f64{ 16, 5, 9, 4, 2, 11, 7, 14, 3, 10, 6, 15, 13, 8, 12, 1 };
    const e = try rref(al, f64, try .init(&m4, 4, 4, 4), null);
    defer e.deinit(al);
    try testing.expectEqualSlices(usize, &.{ 0, 1, 2 }, e.pivots);
    const want = [_]f64{ 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 3, -3, 0 };
    for (want, e.r.data) |w, g| try testing.expectApproxEqAbs(w, g, 1e-12);

    // [1 2 3; 2 4 7]: the second column is no pivot
    const w2 = [_]f64{ 1, 2, 2, 4, 3, 7 };
    const f = try rref(al, f64, try .init(&w2, 2, 3, 2), null);
    defer f.deinit(al);
    try testing.expectEqualSlices(usize, &.{ 0, 2 }, f.pivots);
    for ([_]f64{ 1, 0, 2, 0, 0, 1 }, f.r.data) |w, g| try testing.expectApproxEqAbs(w, g, 1e-12);

    const C = std.math.Complex(f64);
    const c = [_]C{ .init(0, 2), .init(1, 0) };
    const g = try rref(al, C, try .init(&c, 2, 1, 2), null);
    defer g.deinit(al);
    try testing.expectEqual(C.init(1, 0), g.r.data[0]);
    try testing.expectEqual(C.init(0, 0), g.r.data[1]);
}
