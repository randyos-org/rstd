//! Choosing the eigenvalues a reordering moves (MATLAB's `ordschur` and
//! `ordqz`): a logical select, a region keyword, or clusters, over the
//! diagonal blocks of a Schur or generalized Schur form.

const std = @import("std");

const Allocator = std.mem.Allocator;
const Io = std.Io;

/// `ordschur`'s and `ordqz`'s: `DimensionMismatch` unless the matrices and a
/// `.mask` or `.clusters` agree in size; `EigenvaluesTooClose` when a swap
/// would lose accuracy (`trsen`'s or `tgsen`'s `SwapRejected`, MATLAB's
/// error there).
pub const OrderError = error{ NotSquare, DimensionMismatch, EigenvaluesTooClose } || Allocator.Error || Io.Cancelable;

/// MATLAB's keywords: `Re λ < 0`, `Re λ ≥ 0`, `|λ| < 1`, `|λ| ≥ 1` (an
/// infinite `λ` is outside the unit disk; a NaN in no region).
pub const Region = enum { lhp, rhp, udi, udo };

/// Which eigenvalues move, one entry per diagonal row for `.mask` and
/// `.clusters`. In a real form a 2 x 2 block's pair goes together: either
/// row's `true` picks it, and it takes the higher of its rows' indices.
pub const Select = union(enum) {
    /// The selected eigenvalues lead, in their order.
    mask: []const bool,
    /// The eigenvalues in the region lead.
    region: Region,
    /// Clusters by index, the highest first (MATLAB's `clusters`).
    clusters: []const usize,
};

/// A block's eigenvalue as the regions read it: its real part and squared
/// modulus.
pub fn Eigen(comptime F: type) type {
    return struct { re: F, mod2: F };
}

/// `DimensionMismatch` unless a `.mask` or `.clusters` has `n` entries.
pub fn checkLength(select: Select, n: usize) error{DimensionMismatch}!void {
    switch (select) {
        .region => {},
        inline else => |s| if (s.len != n) return error.DimensionMismatch,
    }
}

/// Reorders the form `ctx` holds as `select` asks. `ctx` has `n`,
/// `opensBlock(k)` (row `k` opens a 2 x 2 block), `eigen(k, size)` (that
/// block's `Eigen`) and `reorder(mask)` (one `trsen` or `tgsen`, stable).
pub fn apply(allocator: Allocator, ctx: anytype, select: Select) OrderError!void {
    const n = ctx.n;
    if (n == 0) return;
    switch (select) {
        .mask => |m| try ctx.reorder(m),
        .region => |r| {
            const m = try allocator.alloc(bool, n);
            defer allocator.free(m);
            var k: usize = 0;
            while (k < n) {
                const size: usize = if (ctx.opensBlock(k)) 2 else 1;
                @memset(m[k..][0..size], inRegion(ctx.eigen(k, size), r));
                k += size;
            }
            try ctx.reorder(m);
        },
        .clusters => |c| try clusterOrder(allocator, ctx, c),
    }
}

fn inRegion(e: anytype, r: Region) bool {
    return switch (r) {
        .lhp => e.re < 0,
        .rhp => e.re >= 0,
        .udi => e.mod2 < 1,
        .udo => e.mod2 >= 1,
    };
}

/// Clusters by index, the highest leading: for each index `v` above the
/// lowest, one reorder brings every row of index `v` or more to the front,
/// which keeps those already there in place. `labels` follows the rows
/// through each stable reordering.
fn clusterOrder(allocator: Allocator, ctx: anytype, clusters: []const usize) OrderError!void {
    const n = ctx.n;
    const labels = try allocator.dupe(usize, clusters);
    defer allocator.free(labels);
    var k: usize = 0;
    while (k < n) : (k += 1) {
        if (!ctx.opensBlock(k)) continue;
        const v = @max(labels[k], labels[k + 1]);
        labels[k] = v;
        labels[k + 1] = v;
        k += 1;
    }
    const sorted = try allocator.dupe(usize, labels);
    defer allocator.free(sorted);
    std.mem.sortUnstable(usize, sorted, {}, std.sort.desc(usize));
    const mask = try allocator.alloc(bool, n);
    defer allocator.free(mask);
    const moved = try allocator.alloc(usize, n);
    defer allocator.free(moved);
    for (sorted[0 .. n - 1], sorted[1..]) |v, next| {
        if (v == next) continue;
        for (mask, labels) |*s, l| s.* = l >= v;
        try ctx.reorder(mask);
        // the selected rows lead, each part in its old order
        var front: usize = 0;
        for (labels) |l| {
            if (l < v) continue;
            moved[front] = l;
            front += 1;
        }
        for (labels) |l| {
            if (l >= v) continue;
            moved[front] = l;
            front += 1;
        }
        @memcpy(labels, moved);
    }
}
