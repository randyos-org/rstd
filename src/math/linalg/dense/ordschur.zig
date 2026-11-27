//! Reordering a Schur form (MATLAB's `ordschur`), by `trsen`: a logical
//! select, a region keyword, or clusters.

const std = @import("std");
const lapack = @import("../../kernels/lapack/__root__.zig");
const matrix = @import("../matrix.zig");
const vector = @import("../vector.zig");
const ops = @import("../../kernels/lapack/ops.zig");
const util = @import("util.zig");
const schur_form = @import("schur_form.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = matrix.Mat;
const isComplex = vector.isComplex;
const Real = vector.Real;
const ElementOf = util.ElementOf;
const constOf = util.constOf;
const copyOf = util.copyOf;
const Schur = schur_form.Schur;

/// `ordschur`'s: `DimensionMismatch` unless `u`, `t` and a `.mask` or
/// `.clusters` agree in size; `EigenvaluesTooClose` when a swap would lose
/// accuracy (real `T`: `trsen`'s `SwapRejected`, MATLAB's error there).
pub const OrdschurError = error{ NotSquare, DimensionMismatch, EigenvaluesTooClose } || Allocator.Error || Io.Cancelable;

/// MATLAB's keywords: `Re λ < 0`, `Re λ ≥ 0`, `|λ| < 1`, `|λ| ≥ 1`.
pub const Region = enum { lhp, rhp, udi, udo };

/// Which eigenvalues `ordschur` moves, one entry per row of `T` for
/// `.mask` and `.clusters`. For a real `T`, a 2 x 2 block's pair goes
/// together: either row's `true` picks it, and it takes the higher of its
/// rows' cluster indices.
pub const Select = union(enum) {
    /// The selected eigenvalues lead `T`, in their order.
    mask: []const bool,
    /// The eigenvalues in the region lead `T`.
    region: Region,
    /// Clusters by index, the highest first (MATLAB's `clusters`).
    clusters: []const usize,
};

/// The Schur form `A = U T Uᴴ` reordered as `select` asks (MATLAB's `[US,
/// TS] = ordschur(U, T, select)`, `keyword` or `clusters`), by `trsen`;
/// clusters by one `trsen` per index above the lowest.
pub fn ordschur(allocator: Allocator, io: Io, u: anytype, t: anytype, select: Select) OrdschurError!Schur(ElementOf(@TypeOf(t))) {
    const T = ElementOf(@TypeOf(t));
    const cu = constOf(u);
    const ct = constOf(t);
    if (ct.rows != ct.cols) return error.NotSquare;
    const n = ct.rows;
    if (cu.rows != n or cu.cols != n) return error.DimensionMismatch;
    switch (select) {
        .region => {},
        inline else => |s| if (s.len != n) return error.DimensionMismatch,
    }
    const rt = try copyOf(T, allocator, ct);
    errdefer rt.deinit(allocator);
    const ru = try copyOf(T, allocator, cu);
    errdefer ru.deinit(allocator);
    if (n == 0) return .{ .u = ru, .t = rt };
    switch (select) {
        .mask => |m| try reorder(T, allocator, io, rt, ru, m),
        .region => |r| {
            const m = try allocator.alloc(bool, n);
            defer allocator.free(m);
            regionMask(T, rt, r, m);
            try reorder(T, allocator, io, rt, ru, m);
        },
        .clusters => |c| try clusterOrder(T, allocator, io, rt, ru, c),
    }
    return .{ .u = ru, .t = rt };
}

/// One `trsen` on `rt` and `ru` in place.
fn reorder(comptime T: type, allocator: Allocator, io: Io, rt: Mat(T), ru: Mat(T), mask: []const bool) OrdschurError!void {
    const n = rt.rows;
    const w = try allocator.alloc(T, if (comptime isComplex(T)) n else 2 * n);
    defer allocator.free(w);
    const ev: lapack.Eigenvalues(T) = if (comptime isComplex(T)) w else .{ .re = w[0..n], .im = w[n..] };
    _ = lapack.trsen(allocator, io, T, .none, .update, mask, n, rt.data, n, ru.data, n, ev) catch |e| return switch (e) {
        error.SwapRejected => error.EigenvaluesTooClose,
        else => ops.narrow(OrdschurError, e),
    };
}

/// A real `T`'s row `k` opens a 2 x 2 block.
fn opensBlock(comptime T: type, t: Mat(T), k: usize) bool {
    if (comptime isComplex(T)) return false;
    return k + 1 < t.rows and t.data[k + 1 + k * t.rows] != 0;
}

/// `m[k]`: row `k`'s eigenvalue lies in `r`; a 2 x 2 block's `[a b; c d]`
/// pair has real part `(a + d) / 2` and modulus `sqrt(a d - b c)`.
fn regionMask(comptime T: type, t: Mat(T), r: Region, m: []bool) void {
    const F = Real(T);
    const n = t.rows;
    var k: usize = 0;
    while (k < n) {
        var re: F = undefined;
        var mod2: F = undefined;
        var size: usize = 1;
        if (comptime isComplex(T)) {
            const x = t.data[k + k * n];
            re = x.re;
            mod2 = x.re * x.re + x.im * x.im;
        } else if (opensBlock(T, t, k)) {
            const a = t.data[k + k * n];
            const b = t.data[k + (k + 1) * n];
            const c = t.data[k + 1 + k * n];
            const d = t.data[k + 1 + (k + 1) * n];
            re = (a + d) / 2;
            mod2 = a * d - b * c;
            size = 2;
        } else {
            re = t.data[k + k * n];
            mod2 = re * re;
        }
        const in = switch (r) {
            .lhp => re < 0,
            .rhp => re >= 0,
            .udi => mod2 < 1,
            .udo => mod2 >= 1,
        };
        @memset(m[k..][0..size], in);
        k += size;
    }
}

/// Clusters by index, the highest leading: for each index `v` above the
/// lowest, `trsen` brings every row of index `v` or more to the front, which
/// keeps those already there in place. `labels` follows the rows through
/// each stable reordering.
fn clusterOrder(comptime T: type, allocator: Allocator, io: Io, rt: Mat(T), ru: Mat(T), clusters: []const usize) OrdschurError!void {
    const n = rt.rows;
    const labels = try allocator.dupe(usize, clusters);
    defer allocator.free(labels);
    var k: usize = 0;
    while (k < n) : (k += 1) {
        if (!opensBlock(T, rt, k)) continue;
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
        try reorder(T, allocator, io, rt, ru, mask);
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
