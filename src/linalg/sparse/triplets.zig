//! `Triplets(T)`, the builder: `append` entries in any order, then `toCsc`.

const std = @import("std");
const ops = @import("../lapack/ops.zig");
const util = @import("util.zig");
const csc = @import("csc.zig");

const Allocator = std.mem.Allocator;

/// Entries `x[k]` at row `i[k]`, column `j[k]`, duplicates allowed.
pub fn Triplets(comptime T: type) type {
    return struct {
        rows: u32,
        cols: u32,
        i: std.ArrayList(u32) = .empty,
        j: std.ArrayList(u32) = .empty,
        x: std.ArrayList(T) = .empty,

        const Self = @This();

        pub fn init(rows: usize, cols: usize) error{TooLarge}!Self {
            return .{ .rows = try util.narrow(rows), .cols = try util.narrow(cols) };
        }

        pub fn deinit(self: *Self, allocator: Allocator) void {
            self.i.deinit(allocator);
            self.j.deinit(allocator);
            self.x.deinit(allocator);
            self.* = undefined;
        }

        pub fn append(self: *Self, allocator: Allocator, row: usize, col: usize, value: T) (Allocator.Error || error{IndexOutOfRange})!void {
            if (row >= self.rows or col >= self.cols) return error.IndexOutOfRange;
            try self.i.ensureUnusedCapacity(allocator, 1);
            try self.j.ensureUnusedCapacity(allocator, 1);
            try self.x.ensureUnusedCapacity(allocator, 1);
            self.i.appendAssumeCapacity(@intCast(row));
            self.j.appendAssumeCapacity(@intCast(col));
            self.x.appendAssumeCapacity(value);
        }

        /// The CSC, duplicates summed in insertion order, a sum of zero
        /// kept: stable counting sorts by row, then by column, then the
        /// adjacent duplicates summed.
        pub fn toCsc(self: Self, allocator: Allocator) csc.Error!csc.Csc(T) {
            const n = self.x.items.len;
            _ = try util.narrow(n);
            const ri = self.i.items;
            const cj = self.j.items;
            const next = try allocator.alloc(u32, @as(usize, @max(self.rows, self.cols)) + 1);
            defer allocator.free(next);
            const by_row = try allocator.alloc(u32, n);
            defer allocator.free(by_row);
            const order = try allocator.alloc(u32, n);
            defer allocator.free(order);
            starts(next[0 .. @as(usize, self.rows) + 1], ri);
            for (0..n) |k| {
                by_row[next[ri[k]]] = @intCast(k);
                next[ri[k]] += 1;
            }
            starts(next[0 .. @as(usize, self.cols) + 1], cj);
            for (by_row) |k| {
                order[next[cj[k]]] = k;
                next[cj[k]] += 1;
            }
            // `order` is column-major, rows ascending, duplicates adjacent
            var nnz: usize = 0;
            for (order, 0..) |k, q| nnz += @intFromBool(q == 0 or !same(self, k, order[q - 1]));
            const a: csc.Csc(T) = try .alloc(allocator, self.rows, self.cols, nnz);
            @memset(a.col_start, 0);
            var at: usize = 0;
            for (order, 0..) |k, q| {
                if (q != 0 and same(self, k, order[q - 1])) {
                    a.value[at - 1] = ops.add(T, a.value[at - 1], self.x.items[k]);
                    continue;
                }
                a.row[at] = ri[k];
                a.value[at] = self.x.items[k];
                a.col_start[cj[k] + 1] += 1;
                at += 1;
            }
            for (1..a.col_start.len) |c| a.col_start[c] += a.col_start[c - 1];
            return a;
        }

        fn same(self: Self, k: u32, prev: u32) bool {
            return self.i.items[k] == self.i.items[prev] and self.j.items[k] == self.j.items[prev];
        }
    };
}

/// `next[v]` := the first slot of key `v` in a counting sort of `keys`.
fn starts(next: []u32, keys: []const u32) void {
    @memset(next, 0);
    for (keys) |v| next[v + 1] += 1;
    for (1..next.len) |v| next[v] += next[v - 1];
}
