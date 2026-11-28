//! Tests: `stebz` against the tridiagonal QR algorithm's eigenvalues, as
//! Reference-LAPACK 3.12.0's `TESTING/EIG/dchkst.f` checks `DSTEBZ`: every
//! range and order, `abstol` 0 and `2 safmin`, random, split, graded and
//! clustered matrices.

const std = @import("std");
const types = @import("../../blas/types.zig");
const check = @import("../../blas/_test/check.zig");
const lcheck = @import("check.zig");
const interface = @import("../interface.zig");
const tevd = @import("../tridiag_evd.zig");
const bisect = @import("../stebz.zig");

const testing = std.testing;
const SpectrumRange = interface.SpectrumRange;
const SpectrumOrder = interface.SpectrumOrder;

/// `T`'s kinds: random; random with every seventh `e` zero; graded, `d(i)`
/// and `e(i)` halving down the diagonal; Wilkinson's `W⁺`, `d(i) = |m - i|`
/// and `e = 1`, whose eigenvalues come in close pairs.
const Kind = enum { random, split, graded, wilkinson };

fn tridiagonal(comptime F: type, kind: Kind, d: []F, e: []F) void {
    const n = d.len;
    check.fill(F, d, n * 13 + 1);
    check.fill(F, e, n * 17 + 2);
    switch (kind) {
        .random => {},
        .split => for (e, 0..) |*v, i| {
            if (i % 7 == 6) v.* = 0;
        },
        .graded => {
            for (d, 0..) |*v, i| v.* = std.math.ldexp(v.*, -@as(i32, @intCast(i)));
            for (e, 0..) |*v, i| v.* = std.math.ldexp(v.*, -@as(i32, @intCast(i)) - 1);
        },
        .wilkinson => {
            for (d, 0..) |*v, i| v.* = @abs(@as(F, @floatFromInt(n / 2)) - @as(F, @floatFromInt(i)));
            @memset(e, 1);
        },
    }
}

/// `T`'s eigenvalues, ascending, by the QR algorithm.
fn eigenvalues(comptime F: type, d: []const F, e: []const F, out: []F) !void {
    const al = testing.allocator;
    @memcpy(out, d);
    const off = try al.dupe(F, e);
    defer al.free(off);
    try tevd.qrAlgorithm(F, out, off, null);
    std.mem.sort(F, out, {}, std.sort.asc(F));
}

fn fail(comptime F: type, kind: Kind, n: usize, what: []const u8, value: anytype) error{TestUnexpectedResult} {
    check.say("stebz {s} {t} n {d}: {s} {any}\n", .{ @typeName(F), kind, n, what, value });
    return error.TestUnexpectedResult;
}

fn case(comptime F: type, kind: Kind, n: usize) !void {
    const al = testing.allocator;
    const d = try al.alloc(F, n);
    defer al.free(d);
    const e = try al.alloc(F, n -| 1);
    defer al.free(e);
    tridiagonal(F, kind, d, e);
    const want = try al.alloc(F, n);
    defer al.free(want);
    try eigenvalues(F, d, e, want);
    // `dstech.f`'s tolerance, `sep.in`'s `THRESH` (60) ulps of the largest
    // eigenvalue, past underflow
    const ulp = std.math.floatEps(F);
    var mx: F = 0;
    for (want) |v| mx = @max(mx, @abs(v));
    const tol = 60 * ulp * mx + std.math.floatMin(F) / ulp;
    const w = try al.alloc(F, n);
    defer al.free(w);
    const iblock = try al.alloc(usize, n);
    defer al.free(iblock);
    const isplit = try al.alloc(usize, n);
    defer al.free(isplit);

    for ([_]F{ 0, 2 * std.math.floatMin(F) }) |abstol| {
        // every eigenvalue, both orders
        for ([_]SpectrumOrder{ .ascending, .by_block }) |order| {
            const r = try bisect.stebz(al, F, .all, order, d, e, abstol, w, iblock, isplit);
            if (r.count != n or r.unconverged or r.incomplete) return fail(F, kind, n, "all: found", r);
            try expectBlocks(F, kind, d, e, r.blocks, isplit);
            if (order == .by_block) {
                // ascending within each block, the blocks in order
                for (1..n) |i| if (iblock[i] < iblock[i - 1] or (iblock[i] == iblock[i - 1] and w[i] < w[i - 1])) return fail(F, kind, n, "by_block at", i);
                std.mem.sort(F, w, {}, std.sort.asc(F));
            }
            for (w, want) |x, y| if (!(@abs(x - y) <= tol)) return fail(F, kind, n, "all: eigenvalue", .{ x, y });
        }
        // `W⁺`'s pairs may be closer than bisection can split
        if (n < 3 or kind == .wilkinson) continue;
        // a middle run by index, then by value, its ends halfway to the next
        const first = n / 4;
        const last = n - 1 - n / 3;
        const r = try bisect.stebz(al, F, .{ .indices = .{ .first = first, .last = last } }, .ascending, d, e, abstol, w, iblock, isplit);
        if (r.count != last - first + 1 or r.incomplete) return fail(F, kind, n, "indices: found", r);
        for (w[0..r.count], want[first .. last + 1]) |x, y| if (!(@abs(x - y) <= tol)) return fail(F, kind, n, "indices: eigenvalue", .{ x, y });
        const lo = 0.5 * (want[first - 1] + want[first]);
        const hi = 0.5 * (want[last] + want[last + 1]);
        if (want[first] - want[first - 1] <= 4 * tol or want[last + 1] - want[last] <= 4 * tol) continue;
        const v = try bisect.stebz(al, F, .{ .values = .{ .lo = lo, .hi = hi } }, .ascending, d, e, abstol, w, iblock, isplit);
        if (v.count != last - first + 1) return fail(F, kind, n, "values: found", v);
        for (w[0..v.count], want[first .. last + 1]) |x, y| if (!(@abs(x - y) <= tol)) return fail(F, kind, n, "values: eigenvalue", .{ x, y });
    }
}

/// `isplit` ends a block wherever `e` is zero, and at `n`.
fn expectBlocks(comptime F: type, kind: Kind, d: []const F, e: []const F, blocks: usize, isplit: []const usize) !void {
    const n = d.len;
    // grading may split `T` where `e` is not zero
    if (kind == .graded) return;
    var b: usize = 0;
    for (1..n) |j| if (e[j - 1] == 0) {
        if (b >= blocks or isplit[b] != j) return fail(F, kind, n, "split at", j);
        b += 1;
    };
    if (b + 1 != blocks) return fail(F, kind, n, "blocks", blocks);
    if (isplit[blocks - 1] != n) return fail(F, kind, n, "last block ends at", isplit[blocks - 1]);
}

test "stebz: the QR algorithm's eigenvalues, every range and order, both abstols, every real type" {
    inline for (lcheck.test_types) |F| if (comptime !types.isComplex(F)) {
        for (std.meta.tags(Kind)) |kind| for ([_]usize{ 1, 2, 5, 23, 60 }) |n| try case(F, kind, n);
    };
}

test "stebz: one element, a range missing it, empty sizes, bad arguments" {
    const al = testing.allocator;
    var w: [4]f64 = undefined;
    var iblock: [4]usize = undefined;
    var isplit: [4]usize = undefined;
    const d = [_]f64{ 3, 1, 2, 5 };
    const e = [_]f64{ 1, 0, 1 };
    const one = try bisect.stebz(al, f64, .all, .ascending, d[0..1], e[0..0], 0, &w, &iblock, &isplit);
    try testing.expectEqual(@as(usize, 1), one.count);
    try testing.expectEqual(@as(f64, 3), w[0]);
    const none = try bisect.stebz(al, f64, .{ .values = .{ .lo = 3, .hi = 4 } }, .ascending, d[0..1], e[0..0], 0, &w, &iblock, &isplit);
    try testing.expectEqual(@as(usize, 0), none.count);
    const empty = try bisect.stebz(al, f64, .all, .ascending, d[0..0], e[0..0], 0, w[0..0], iblock[0..0], isplit[0..0]);
    try testing.expectEqual(@as(usize, 0), empty.count);
    // two blocks, by block: `{3, 1}`'s pair, then `{2, 5}`'s
    const two = try bisect.stebz(al, f64, .all, .by_block, &d, &e, 0, &w, &iblock, &isplit);
    try testing.expectEqual(@as(usize, 2), two.blocks);
    try testing.expectEqualSlices(usize, &.{ 2, 4 }, isplit[0..2]);
    try testing.expectEqualSlices(usize, &.{ 0, 0, 1, 1 }, &iblock);
    const ee = testing.expectError;
    try ee(error.EmptyInterval, bisect.stebz(al, f64, .{ .values = .{ .lo = 1, .hi = 1 } }, .ascending, &d, &e, 0, &w, &iblock, &isplit));
    try ee(error.IndexOutOfRange, bisect.stebz(al, f64, .{ .indices = .{ .first = 2, .last = 1 } }, .ascending, &d, &e, 0, &w, &iblock, &isplit));
    try ee(error.IndexOutOfRange, bisect.stebz(al, f64, .{ .indices = .{ .first = 0, .last = 4 } }, .ascending, &d, &e, 0, &w, &iblock, &isplit));
    try ee(error.SliceTooShort, bisect.stebz(al, f64, .all, .ascending, &d, e[0..2], 0, &w, &iblock, &isplit));
    try ee(error.SliceTooShort, bisect.stebz(al, f64, .all, .ascending, &d, &e, 0, w[0..3], &iblock, &isplit));
}
