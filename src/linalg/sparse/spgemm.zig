//! `C := A B`, all CSC, by Gustavson's column form (src/linalg/design.md,
//! "Sparse matrices", `spgemm`): `C(:, j) = Σ_k B(k, j) A(:, k)`, reading
//! columns only (Gustavson, "Two fast algorithms for sparse matrices", ACM
//! TOMS 4(3), 1978). The symbolic pass marks each column's rows and sorts
//! them; the numeric pass sums into `C`'s own slots, `k` ascending.

const std = @import("std");
const ops = @import("../lapack/ops.zig");
const util = @import("util.zig");
const csc = @import("csc.zig");

const Allocator = std.mem.Allocator;

/// No row: unmarked, in the symbolic and numeric scratch.
const none = std.math.maxInt(u32);

pub fn symbolic(comptime T: type, a: csc.Csc(T), allocator: Allocator, b: csc.Csc(T)) csc.SpgemmError!csc.Csc(T) {
    if (a.cols != b.rows) return error.DimensionMismatch;
    const mark = try allocator.alloc(u32, a.rows);
    defer allocator.free(mark);
    @memset(mark, none);
    const col_start = try allocator.alloc(u32, @as(usize, b.cols) + 1);
    errdefer allocator.free(col_start);
    var rows: std.ArrayList(u32) = .empty;
    defer rows.deinit(allocator);
    col_start[0] = 0;
    for (0..b.cols) |j| {
        const first = rows.items.len;
        for (b.col_start[j]..b.col_start[j + 1]) |q| {
            const k = b.row[q];
            for (a.row[a.col_start[k]..a.col_start[k + 1]]) |i| {
                if (mark[i] == j) continue;
                mark[i] = @intCast(j);
                try rows.append(allocator, i);
            }
        }
        std.mem.sort(u32, rows.items[first..], {}, std.sort.asc(u32));
        col_start[j + 1] = try util.narrow(rows.items.len);
    }
    const row = try allocator.dupe(u32, rows.items);
    errdefer allocator.free(row);
    const value = try allocator.alloc(T, row.len);
    @memset(value, ops.zero(T));
    return .{ .rows = a.rows, .cols = b.cols, .col_start = col_start, .row = row, .value = value };
}

pub fn numeric(comptime T: type, a: csc.Csc(T), allocator: Allocator, b: csc.Csc(T), c: csc.Csc(T)) csc.NumericError!void {
    if (a.cols != b.rows or c.rows != a.rows or c.cols != b.cols) return error.DimensionMismatch;
    // `slot[i]`: row `i`'s place in `c.value` for the column at hand
    const slot = try allocator.alloc(u32, a.rows);
    defer allocator.free(slot);
    @memset(slot, none);
    for (0..b.cols) |j| {
        const s = c.col_start[j];
        const e = c.col_start[j + 1];
        for (s..e) |p| slot[c.row[p]] = @intCast(p);
        defer {
            for (c.row[s..e]) |i| slot[i] = none;
        }
        @memset(c.value[s..e], ops.zero(T));
        for (b.col_start[j]..b.col_start[j + 1]) |q| {
            const k = b.row[q];
            const bkj = b.value[q];
            for (a.col_start[k]..a.col_start[k + 1]) |p| {
                const at = slot[a.row[p]];
                if (at == none) return error.PatternMismatch;
                c.value[at] = ops.add(T, c.value[at], ops.mul(T, a.value[p], bkj));
            }
        }
    }
}

pub fn spgemm(comptime T: type, a: csc.Csc(T), allocator: Allocator, b: csc.Csc(T)) csc.SpgemmError!csc.Csc(T) {
    var c = try symbolic(T, a, allocator, b);
    errdefer c.deinit(allocator);
    numeric(T, a, allocator, b, c) catch |e| return switch (e) {
        error.OutOfMemory => error.OutOfMemory,
        // `c` is `a` times `b`'s pattern
        error.DimensionMismatch, error.PatternMismatch => unreachable,
    };
    return c;
}
