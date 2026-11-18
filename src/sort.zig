//! Sorting and set operations (rstd's design.md, "Sorting and sets"), over
//! what `std.sort` already has, and the set operations (MATLAB's `unique`,
//! `intersect`, `union`, `setdiff`, `ismember`). Each NaN is distinct: it
//! matches nothing, itself included.
//!
//! Floats order NaNs last ascending and first descending, as MATLAB's `sort`;
//! complex values by magnitude, then angle, a NaN part counting as a NaN;
//! -0 and +0 compare equal and keep their order. `std.sort.asc` is no strict
//! weak order once a NaN appears.

const std = @import("std");
const blas_util = @import("linalg/blas/util.zig");

const Allocator = std.mem.Allocator;

pub const Direction = enum { ascend, descend };

/// Which end of a run of equal values `searchSorted` returns.
pub const Side = enum { left, right };

/// The set operations' output order: sorted, or by first appearance (MATLAB's
/// `'stable'`).
pub const Keep = enum { sorted, stable };

/// A `lessThan` for `std.sort`, NaN-aware for floats.
pub fn order(comptime T: type, comptime direction: Direction) fn (void, T, T) bool {
    return struct {
        fn lessThan(_: void, a: T, b: T) bool {
            return precedes(T, direction, a, b);
        }
    }.lessThan;
}

/// Stable.
pub fn sort(comptime T: type, xs: []T, direction: Direction) void {
    switch (direction) {
        inline else => |d| std.sort.block(T, xs, {}, order(T, d)),
    }
}

/// The permutation that sorts `xs`, stably: MATLAB's `[~, i] = sort(xs)`.
pub fn argsort(allocator: Allocator, comptime T: type, xs: []const T, direction: Direction) error{OutOfMemory}![]usize {
    const index = try allocator.alloc(usize, xs.len);
    for (index, 0..) |*k, i| k.* = i;
    const By = struct {
        xs: []const T,
        direction: Direction,

        fn lessThan(by: @This(), a: usize, b: usize) bool {
            return precedes(T, by.direction, by.xs[a], by.xs[b]);
        }
    };
    std.sort.block(usize, index, By{ .xs = xs, .direction = direction }, By.lessThan);
    return index;
}

/// Where `x` goes in `xs`, sorted by `direction`: before its equals (`.left`)
/// or after them (`.right`). NumPy's `searchsorted`.
pub fn searchSorted(comptime T: type, xs: []const T, x: T, direction: Direction, side: Side) usize {
    const Key = struct {
        x: T,
        direction: Direction,

        fn left(key: @This(), item: T) bool {
            return precedes(T, key.direction, item, key.x);
        }

        fn right(key: @This(), item: T) bool {
            return !precedes(T, key.direction, key.x, item);
        }
    };
    const key: Key = .{ .x = x, .direction = direction };
    return switch (side) {
        .left => std.sort.partitionPoint(T, xs, key, Key.left),
        .right => std.sort.partitionPoint(T, xs, key, Key.right),
    };
}

fn precedes(comptime T: type, direction: Direction, a: T, b: T) bool {
    return switch (direction) {
        .ascend => less(T, a, b),
        .descend => less(T, b, a),
    };
}

/// NaNs sort past every number.
fn less(comptime T: type, a: T, b: T) bool {
    if (comptime isComplexType(T)) {
        // MATLAB's order: by magnitude, then angle in (-π, π]; a NaN part last
        if (isNan(T, a)) return false;
        if (isNan(T, b)) return true;
        // rounded, or the x87 compares one stored and one wider (issue #478)
        const R = @TypeOf(a.re);
        const ma = blas_util.rounded(R, std.math.hypot(a.re, a.im));
        const mb = blas_util.rounded(R, std.math.hypot(b.re, b.im));
        if (ma != mb) return ma < mb;
        return blas_util.rounded(R, std.math.atan2(a.im, a.re)) < blas_util.rounded(R, std.math.atan2(b.im, b.re));
    }
    if (comptime @typeInfo(T) == .float) {
        if (std.math.isNan(a)) return false;
        if (std.math.isNan(b)) return true;
    }
    return a < b;
}

pub const SetError = error{OutOfMemory};

/// MATLAB's `[values, first, inverse] = unique(xs)`; each NaN is distinct.
pub fn Unique(comptime T: type) type {
    return struct {
        values: []T,
        /// `values[k] == xs[first[k]]`, the first occurrence.
        first: []usize,
        /// `xs[i] == values[inverse[i]]`.
        inverse: []usize,

        pub fn deinit(self: *@This(), allocator: Allocator) void {
            allocator.free(self.values);
            allocator.free(self.first);
            allocator.free(self.inverse);
            self.* = undefined;
        }
    };
}

pub fn unique(allocator: Allocator, comptime T: type, xs: []const T, keep: Keep) SetError!Unique(T) {
    const idx = try argsort(allocator, T, xs, .ascend);
    defer allocator.free(idx);
    const inverse = try allocator.alloc(usize, xs.len);
    errdefer allocator.free(inverse);
    // the runs of equal values; the stable sort puts each run's first
    // occurrence first
    var firsts: std.ArrayList(usize) = .empty;
    defer firsts.deinit(allocator);
    for (idx, 0..) |i, k| {
        if (k == 0 or !equal(T, xs[idx[k - 1]], xs[i])) try firsts.append(allocator, i);
        inverse[i] = firsts.items.len - 1;
    }
    const n = firsts.items.len;
    const first = try allocator.dupe(usize, firsts.items);
    errdefer allocator.free(first);
    if (keep == .stable) {
        // groups by first appearance; `rank[g]` is group `g`'s new place
        const by = try argsort(allocator, usize, first, .ascend);
        defer allocator.free(by);
        const rank = try allocator.alloc(usize, n);
        defer allocator.free(rank);
        for (by, 0..) |g, r| rank[g] = r;
        for (inverse) |*v| v.* = rank[v.*];
        std.mem.sort(usize, first, {}, std.sort.asc(usize));
    }
    const values = try allocator.alloc(T, n);
    for (values, first) |*v, i| v.* = xs[i];
    return .{ .values = values, .first = first, .inverse = inverse };
}

/// Equal for the set operations: neither precedes the other, and no NaN.
fn equal(comptime T: type, a: T, b: T) bool {
    if (isNan(T, a) or isNan(T, b)) return false;
    return !precedes(T, .ascend, a, b) and !precedes(T, .ascend, b, a);
}

/// `xs`'s distinct values sorted, with their first indices.
fn distinct(allocator: Allocator, comptime T: type, xs: []const T) SetError!Unique(T) {
    return unique(allocator, T, xs, .sorted);
}

/// Distinct values with their sources in `a` and `b`, as MATLAB's
/// `[values, ia, ib]`; what the indices mean is per function.
pub fn Matched(comptime T: type) type {
    return struct {
        values: []T,
        a_index: []usize,
        b_index: []usize,

        pub fn deinit(self: *@This(), allocator: Allocator) void {
            allocator.free(self.values);
            allocator.free(self.a_index);
            allocator.free(self.b_index);
            self.* = undefined;
        }
    };
}

/// `values[k] == a[a_index[k]] == b[b_index[k]]`, first occurrences.
pub fn intersect(allocator: Allocator, comptime T: type, a: []const T, b: []const T, keep: Keep) SetError!Matched(T) {
    var ua = try distinct(allocator, T, a);
    defer ua.deinit(allocator);
    var ub = try distinct(allocator, T, b);
    defer ub.deinit(allocator);
    var hits: std.ArrayList([2]usize) = .empty;
    defer hits.deinit(allocator);
    var i: usize = 0;
    var j: usize = 0;
    while (i < ua.values.len and j < ub.values.len) {
        const x = ua.values[i];
        const y = ub.values[j];
        if (equal(T, x, y)) {
            try hits.append(allocator, .{ ua.first[i], ub.first[j] });
            i += 1;
            j += 1;
        } else if (precedes(T, .ascend, x, y) or isNan(T, y)) i += 1 else j += 1;
    }
    // `.stable`: by appearance in `a`
    if (keep == .stable) std.mem.sort([2]usize, hits.items, {}, struct {
        fn lt(_: void, p: [2]usize, q: [2]usize) bool {
            return p[0] < q[0];
        }
    }.lt);
    const out: Matched(T) = .{
        .values = try allocator.alloc(T, hits.items.len),
        .a_index = try allocator.alloc(usize, hits.items.len),
        .b_index = try allocator.alloc(usize, hits.items.len),
    };
    for (hits.items, out.values, out.a_index, out.b_index) |h, *v, *x, *y| {
        v.* = a[h[0]];
        x.* = h[0];
        y.* = h[1];
    }
    return out;
}

/// A NaN, or for complex `T` a NaN part.
fn isNan(comptime T: type, x: T) bool {
    if (comptime isComplexType(T)) return std.math.isNan(x.re) or std.math.isNan(x.im);
    return if (comptime @typeInfo(T) == .float) std.math.isNan(x) else false;
}

/// `std.math.Complex` and other `re`, `im` structs.
fn isComplexType(comptime T: type) bool {
    return @typeInfo(T) == .@"struct" and @hasField(T, "re") and @hasField(T, "im");
}

/// `values` holds `a[a_index]` and `b[b_index]`, each distinct value once,
/// from `a` where both have it. Not `union`, a keyword.
pub fn setUnion(allocator: Allocator, comptime T: type, a: []const T, b: []const T, keep: Keep) SetError!Matched(T) {
    var ua = try unique(allocator, T, a, keep);
    defer ua.deinit(allocator);
    var ub = try unique(allocator, T, b, keep);
    defer ub.deinit(allocator);
    var sa = try distinct(allocator, T, a);
    defer sa.deinit(allocator);
    // `b`'s values `a` lacks, found by search in `a`'s sorted values
    var only_b: std.ArrayList(usize) = .empty;
    defer only_b.deinit(allocator);
    for (ub.values, ub.first) |v, f| if (!contains(T, sa.values, v)) try only_b.append(allocator, f);
    const n = ua.values.len + only_b.items.len;
    var out: Matched(T) = .{
        .values = try allocator.alloc(T, n),
        .a_index = try allocator.dupe(usize, ua.first),
        .b_index = try only_b.toOwnedSlice(allocator),
    };
    @memcpy(out.values[0..ua.values.len], ua.values);
    for (out.values[ua.values.len..], out.b_index) |*v, f| v.* = b[f];
    // sorted: merge the two sorted runs (both are sorted already)
    if (keep == .sorted) std.sort.block(T, out.values, {}, order(T, .ascend));
    return out;
}

/// Whether the sorted `vs` hold a value equal to `x`.
fn contains(comptime T: type, vs: []const T, x: T) bool {
    if (isNan(T, x)) return false;
    const k = searchSorted(T, vs, x, .ascend, .left);
    return k < vs.len and equal(T, vs[k], x);
}

/// Values picked from one slice: `values[k] == xs[index[k]]`.
pub fn Selected(comptime T: type) type {
    return struct {
        values: []T,
        index: []usize,

        pub fn deinit(self: *@This(), allocator: Allocator) void {
            allocator.free(self.values);
            allocator.free(self.index);
            self.* = undefined;
        }
    };
}

/// The distinct values of `a` not in `b`, indexed into `a`.
pub fn setDiff(allocator: Allocator, comptime T: type, a: []const T, b: []const T, keep: Keep) SetError!Selected(T) {
    var ua = try unique(allocator, T, a, keep);
    defer ua.deinit(allocator);
    var sb = try distinct(allocator, T, b);
    defer sb.deinit(allocator);
    var keep_k: std.ArrayList(usize) = .empty;
    defer keep_k.deinit(allocator);
    for (ua.values, 0..) |v, k| if (!contains(T, sb.values, v)) try keep_k.append(allocator, k);
    const out: Selected(T) = .{
        .values = try allocator.alloc(T, keep_k.items.len),
        .index = try allocator.alloc(usize, keep_k.items.len),
    };
    for (keep_k.items, out.values, out.index) |k, *v, *i| {
        v.* = ua.values[k];
        i.* = ua.first[k];
    }
    return out;
}

/// For each `a[i]`, the first index of it in `set`, or null: MATLAB's
/// `[found, at] = ismember(a, set)`, `found` being `at != null`.
pub fn isMember(allocator: Allocator, comptime T: type, a: []const T, set: []const T) SetError![]?usize {
    var us = try distinct(allocator, T, set);
    defer us.deinit(allocator);
    const at = try allocator.alloc(?usize, a.len);
    for (a, at) |x, *r| {
        r.* = null;
        if (isNan(T, x)) continue;
        const k = searchSorted(T, us.values, x, .ascend, .left);
        if (k < us.values.len and equal(T, us.values[k], x)) r.* = us.first[k];
    }
    return at;
}

// ---- tests ----

const testing = std.testing;

const nan = std.math.nan(f64);

test "sort: NaNs and signed zeros" {
    var up = [_]f64{ 3, nan, 1, -0.0, 0.0, 2 };
    sort(f64, &up, .ascend);
    try testing.expectEqualSlices(f64, &.{ -0.0, 0.0, 1, 2, 3 }, up[0..5]);
    try testing.expect(std.math.signbit(up[0]) and !std.math.signbit(up[1]));
    try testing.expect(std.math.isNan(up[5]));

    var down = [_]f64{ 3, nan, 1, -0.0, 0.0, 2 };
    sort(f64, &down, .descend);
    try testing.expect(std.math.isNan(down[0]));
    try testing.expectEqualSlices(f64, &.{ 3, 2, 1, -0.0, 0.0 }, down[1..]);
    try testing.expect(std.math.signbit(down[4]) and !std.math.signbit(down[5]));

    var ints = [_]i32{ 5, -1, 3 };
    sort(i32, &ints, .descend);
    try testing.expectEqualSlices(i32, &.{ 5, 3, -1 }, &ints);
}

test "argsort is stable" {
    const index = try argsort(testing.allocator, f32, &.{ 30, 10, 20, 10 }, .ascend);
    defer testing.allocator.free(index);
    try testing.expectEqualSlices(usize, &.{ 1, 3, 2, 0 }, index);

    const down = try argsort(testing.allocator, f32, &.{ 30, 10, 20, 10 }, .descend);
    defer testing.allocator.free(down);
    try testing.expectEqualSlices(usize, &.{ 0, 2, 1, 3 }, down);
}

test searchSorted {
    const up = [_]u8{ 1, 2, 2, 3 };
    try testing.expectEqual(1, searchSorted(u8, &up, 2, .ascend, .left));
    try testing.expectEqual(3, searchSorted(u8, &up, 2, .ascend, .right));
    try testing.expectEqual(0, searchSorted(u8, &up, 0, .ascend, .left));
    try testing.expectEqual(4, searchSorted(u8, &up, 9, .ascend, .right));

    const down = [_]f64{ nan, 3, 2, 2, 1 };
    try testing.expectEqual(2, searchSorted(f64, &down, 2, .descend, .left));
    try testing.expectEqual(4, searchSorted(f64, &down, 2, .descend, .right));
    try testing.expectEqual(0, searchSorted(f64, &down, nan, .descend, .left));
    try testing.expectEqual(1, searchSorted(f64, &down, nan, .descend, .right));
}

test unique {
    const xs = [_]i32{ 3, 1, 3, 2, 1 };
    var u = try unique(testing.allocator, i32, &xs, .sorted);
    defer u.deinit(testing.allocator);
    try testing.expectEqualSlices(i32, &.{ 1, 2, 3 }, u.values);
    try testing.expectEqualSlices(usize, &.{ 1, 3, 0 }, u.first);
    try testing.expectEqualSlices(usize, &.{ 2, 0, 2, 1, 0 }, u.inverse);

    var s = try unique(testing.allocator, i32, &xs, .stable);
    defer s.deinit(testing.allocator);
    try testing.expectEqualSlices(i32, &.{ 3, 1, 2 }, s.values);

    var n = try unique(testing.allocator, f64, &.{ nan, 1, nan }, .sorted);
    defer n.deinit(testing.allocator);
    try testing.expectEqual(3, n.values.len);
}

test "intersect, setUnion, setDiff" {
    const a = [_]i32{ 5, 1, 3, 1 };
    const b = [_]i32{ 3, 4, 5 };
    var i = try intersect(testing.allocator, i32, &a, &b, .sorted);
    defer i.deinit(testing.allocator);
    try testing.expectEqualSlices(i32, &.{ 3, 5 }, i.values);
    try testing.expectEqualSlices(usize, &.{ 2, 0 }, i.a_index);
    try testing.expectEqualSlices(usize, &.{ 0, 2 }, i.b_index);

    var u = try setUnion(testing.allocator, i32, &a, &b, .sorted);
    defer u.deinit(testing.allocator);
    try testing.expectEqualSlices(i32, &.{ 1, 3, 4, 5 }, u.values);
    try testing.expectEqualSlices(usize, &.{1}, u.b_index);

    var d = try setDiff(testing.allocator, i32, &a, &b, .sorted);
    defer d.deinit(testing.allocator);
    try testing.expectEqualSlices(i32, &.{1}, d.values);
    try testing.expectEqualSlices(usize, &.{1}, d.index);
}

test isMember {
    const at = try isMember(testing.allocator, i32, &.{ 4, 7, 4 }, &.{ 1, 4, 4 });
    defer testing.allocator.free(at);
    try testing.expectEqualSlices(?usize, &.{ 1, null, 1 }, at);
}

test "set operations: stable order, NaNs distinct, -0 equal to +0" {
    const al = testing.allocator;
    const a = [_]f64{ 5, nan, 1, 3, 1, -0.0 };
    const b = [_]f64{ 3, nan, 0, 4, 5 };
    var i = try intersect(al, f64, &a, &b, .stable);
    defer i.deinit(al);
    try testing.expectEqualSlices(f64, &.{ 5, 3, -0.0 }, i.values);
    try testing.expectEqualSlices(usize, &.{ 0, 3, 5 }, i.a_index);
    var u = try setUnion(al, f64, &a, &b, .stable);
    defer u.deinit(al);
    // a's five distinct values by appearance, then b's NaN and 4 (its 0 is a's -0)
    try testing.expectEqual(7, u.values.len);
    try testing.expectEqualSlices(usize, &.{ 1, 3 }, u.b_index);
    var d = try setDiff(al, f64, &a, &b, .stable);
    defer d.deinit(al);
    try testing.expectEqual(1, d.values[1]);
    try testing.expectEqualSlices(usize, &.{ 1, 2 }, d.index);
    const at = try isMember(al, f64, &.{ nan, 0, 4 }, &b);
    defer al.free(at);
    try testing.expectEqualSlices(?usize, &.{ null, 2, 3 }, at);
}

test "Complex: MATLAB's order, by magnitude then angle, NaN last" {
    const C = std.math.Complex(f64);
    var xs = [_]C{ .init(0, 2), .init(-1, 0), .init(nan, 0), .init(1, 0), .init(0, -1), .init(2, 0) };
    sort(C, &xs, .ascend);
    // |1| ties: angles -π/2 (−i), 0 (1), π (−1); then |2|: 0 (2), π/2 (2i)
    const want = [_]C{ .init(0, -1), .init(1, 0), .init(-1, 0), .init(2, 0), .init(0, 2) };
    for (want, xs[0..5]) |w, x| try testing.expect(w.re == x.re and w.im == x.im);
    try testing.expect(std.math.isNan(xs[5].re));
    var u = try unique(testing.allocator, C, &.{ .init(1, 1), .init(1, 1), .init(nan, 1) }, .sorted);
    defer u.deinit(testing.allocator);
    try testing.expectEqual(2, u.values.len);
}
