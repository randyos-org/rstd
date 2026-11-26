//! The iterative solver stubs: shape checks, then `NotImplemented` (N2).

const std = @import("std");
const check = @import("../../blas/_test/check.zig");
const sparse = @import("../__root__.zig");

const testing = std.testing;
const Io = std.Io;

fn identity(ctx: *const anyopaque, io: Io, x: []const f64, y: []f64) Io.Cancelable!void {
    _ = ctx;
    _ = io;
    @memcpy(y, x);
}

test "cg and gmres: DimensionMismatch, else NotImplemented" {
    const al = testing.allocator;
    const a: sparse.LinearOperator(f64) = .{ .ctx = undefined, .apply = identity };
    const b = [_]f64{ 1, 2, 3 };
    var x: [3]f64 = @splat(0);
    try testing.expectError(error.DimensionMismatch, sparse.cg(al, check.io, f64, a, null, 1e-10, 100, &b, x[0..2]));
    try testing.expectError(error.NotImplemented, sparse.cg(al, check.io, f64, a, a, 1e-10, 100, &b, &x));
    try testing.expectError(error.DimensionMismatch, sparse.gmres(al, check.io, f64, a, null, 1e-10, 100, 20, &b, x[0..2]));
    try testing.expectError(error.NotImplemented, sparse.gmres(al, check.io, f64, a, null, 1e-10, 100, 20, &b, &x));
}
