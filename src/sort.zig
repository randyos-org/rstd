//! Sorting and set operations (rstd's design.md, "Sorting and sets"), over
//! what `std.sort` already has. The set operations are stubs: bodies return `error.NotImplemented`,
//! which leaves `SetError` when they land, and their tests skip until then.
//!
//! Floats order NaNs last ascending and first descending, as MATLAB's `sort`;
//! -0 and +0 compare equal and keep their order. `std.sort.asc` is no strict
//! weak order once a NaN appears.

const std = @import("std");

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
    if (comptime @typeInfo(T) == .float) {
        if (std.math.isNan(a)) return false;
        if (std.math.isNan(b)) return true;
    }
    return a < b;
}

pub const SetError = error{ OutOfMemory, NotImplemented };

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
    _ = allocator;
    _ = xs;
    _ = keep;
    return error.NotImplemented;
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
    _ = allocator;
    _ = a;
    _ = b;
    _ = keep;
    return error.NotImplemented;
}

/// `values` holds `a[a_index]` and `b[b_index]`, each distinct value once,
/// from `a` where both have it. Not `union`, a keyword.
pub fn setUnion(allocator: Allocator, comptime T: type, a: []const T, b: []const T, keep: Keep) SetError!Matched(T) {
    _ = allocator;
    _ = a;
    _ = b;
    _ = keep;
    return error.NotImplemented;
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
    _ = allocator;
    _ = a;
    _ = b;
    _ = keep;
    return error.NotImplemented;
}

/// For each `a[i]`, the first index of it in `set`, or null: MATLAB's
/// `[found, at] = ismember(a, set)`, `found` being `at != null`.
pub fn isMember(allocator: Allocator, comptime T: type, a: []const T, set: []const T) SetError![]?usize {
    _ = allocator;
    _ = a;
    _ = set;
    return error.NotImplemented;
}

// ---- tests ----

const testing = std.testing;

/// A stub's tests skip until its body lands.
fn skip(err: anyerror) anyerror {
    return if (err == error.NotImplemented) error.SkipZigTest else err;
}

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
    var u = unique(testing.allocator, i32, &xs, .sorted) catch |e| return skip(e);
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
    var i = intersect(testing.allocator, i32, &a, &b, .sorted) catch |e| return skip(e);
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
    const at = isMember(testing.allocator, i32, &.{ 4, 7, 4 }, &.{ 1, 4, 4 }) catch |e| return skip(e);
    defer testing.allocator.free(at);
    try testing.expectEqualSlices(?usize, &.{ 1, null, 1 }, at);
}
