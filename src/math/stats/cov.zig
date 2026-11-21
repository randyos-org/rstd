//! `cov` and `corrcoef` over the public layer's views.

const std = @import("std");
const blas_types = @import("../kernels/blas/types.zig");
const matrix = @import("../linalg/matrix.zig");
const l3 = @import("../kernels/blas/blas3/interface.zig");
const util = @import("util.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Elem = util.Elem;
const Wide = util.Wide;
const testing = std.testing;

pub const CovError = error{ EmptyInput, OutOfMemory, Canceled };

/// The covariance of `x`'s columns (variables) over its rows (observations),
/// as MATLAB's `cov`: divided by `rows - 1` (by 1 for one row). `syrk` over
/// the centered columns, both triangles filled. `x` is a `linalg.View` or
/// `ConstView` of `f32`, `f64` or the BLAS's extended type.
pub fn cov(allocator: Allocator, io: Io, x: anytype) CovError!matrix.Mat(Elem(@TypeOf(x.data))) {
    const F = Elem(@TypeOf(x.data));
    comptime if (!(F == f32 or F == f64 or F == blas_types.XDouble)) @compileError("rstd.stats.cov: " ++ @typeName(F) ++ " is not a BLAS float type");
    const m = x.rows;
    const n = x.cols;
    if (m == 0) return error.EmptyInput;
    // the columns, centered
    const a = try allocator.alloc(F, m * n);
    defer allocator.free(a);
    for (0..n) |j| {
        var s: Wide(F) = 0;
        for (0..m) |i| s += x.data[i + j * x.ld];
        const mu: F = @floatCast(s / @as(Wide(F), @floatFromInt(m)));
        for (0..m) |i| a[i + j * m] = x.data[i + j * x.ld] - mu;
    }
    var c = try matrix.Mat(F).init(allocator, n, n);
    errdefer c.deinit(allocator);
    if (n == 0) return c;
    const d: F = @floatFromInt(if (m > 1) m - 1 else 1);
    l3.syrk(allocator, io, F, .l, .t, @intCast(n), @intCast(m), 1 / d, a.ptr, @intCast(m), 0, c.data.ptr, @intCast(n)) catch |e| switch (e) {
        error.OutOfMemory => |z| return z,
        error.Canceled => |z| return z,
        // the sizes above are valid
        else => unreachable,
    };
    for (0..n) |j| for (j + 1..n) |i| {
        c.data[j + i * n] = c.data[i + j * n];
    };
    return c;
}

/// The correlation coefficients of `x`'s columns, as MATLAB's `corrcoef`:
/// `cov` scaled by the standard deviations, the diagonal 1 and the rest
/// clamped to [-1, 1]; NaN where a column is constant.
pub fn corrcoef(allocator: Allocator, io: Io, x: anytype) CovError!matrix.Mat(Elem(@TypeOf(x.data))) {
    const c = try cov(allocator, io, x);
    const n = c.rows;
    for (0..n) |j| for (0..n) |i| {
        if (i == j) continue;
        const r = c.data[i + j * n] / @sqrt(c.data[i + i * n] * c.data[j + j * n]);
        c.data[i + j * n] = if (std.math.isNan(r)) r else std.math.clamp(r, -1, 1);
    };
    for (0..n) |i| {
        const v = c.data[i + i * n];
        c.data[i + i * n] = if (v > 0 and std.math.isFinite(v)) 1 else std.math.nan(@TypeOf(v));
    }
    return c;
}

test "cov and corrcoef: MATLAB's on columns" {
    const al = testing.allocator;
    // columns (1, 2, 3), (2, 4, 6), (3, 1, 2)
    const data = [_]f64{ 1, 2, 3, 2, 4, 6, 3, 1, 2 };
    const x = try matrix.ConstView(f64).init(&data, 3, 3, 3);
    var c = try cov(al, testing.io, x);
    defer c.deinit(al);
    const want = [_]f64{ 1, 2, -0.5, 2, 4, -1, -0.5, -1, 1 };
    for (want, c.data) |w, g| try testing.expectApproxEqAbs(w, g, 1e-14);
    var r = try corrcoef(al, testing.io, x);
    defer r.deinit(al);
    try testing.expectEqual(1, r.data[0]);
    try testing.expectApproxEqAbs(1, r.data[1], 1e-14);
    try testing.expectApproxEqAbs(-0.5, r.data[2], 1e-14);
    try testing.expectError(error.EmptyInput, cov(al, testing.io, try matrix.ConstView(f64).init(&data, 0, 3, 1)));
}
