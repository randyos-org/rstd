//! CHOL over LAPACK's `potrf(.u)` (M2), in place of WPOFA (LINPACK's ZPOFA);
//! a NaN pivot fails, a diagonal's imaginary part is ignored (M1's E9).

const std = @import("std");
const scalar = @import("scalar.zig");
const CxMat = @import("CxMat.zig");
const cholesky = @import("../lapack/cholesky.zig");
const ops = @import("../lapack/ops.zig");
const Cx = scalar.Cx;
const Allocator = std.mem.Allocator;
const Io = std.Io;

pub const Error = error{NotPositiveDefinite} || Allocator.Error || Io.Cancelable;

/// CHOL: factor the Hermitian positive-definite `self` (its upper triangle
/// read) in place into upper-triangular `R`, `ctrans(R)*R == A`, with the
/// strict lower triangle zeroed. `error.NotPositiveDefinite` puts the
/// failed leading minor's 1-based order in `minor`, `self` partly factored.
/// An empty `self` is its own factor.
pub fn chol(self: CxMat, allocator: Allocator, io: Io, minor: ?*usize) Error!void {
    // `potrf` rejects an empty matrix's `lda` of 0
    if (self.n == 0) return;
    const n = self.n;
    // MATLAB's `p`: the failed pivot's index, which `potrf` does not return
    const failed = if (self.isReal())
        try factored(cholesky.potrfChecked(allocator, io, f64, .u, n, self.r, self.lda))
    else blk: {
        const z = try self.toPacked(Cx, allocator);
        defer allocator.free(z);
        defer self.fromPacked(Cx, z);
        break :blk try factored(cholesky.potrfChecked(allocator, io, Cx, .u, n, z, n));
    };
    if (failed) |k| {
        // MATLAB's order, 1-based
        if (minor) |m| m.* = k + 1;
        return error.NotPositiveDefinite;
    }
    self.triu(0);
}

fn factored(r: cholesky.CheckedError!?usize) (Allocator.Error || Io.Cancelable)!?usize {
    // the arguments come from `self`
    return r catch |e| ops.narrow(Allocator.Error || Io.Cancelable, e);
}

const zigstd = @import("../../zigstd.zig");
const testing = std.testing;

test "chol factors a known real SPD matrix and zeros the lower triangle" {
    // A = [[4, 2], [2, 3]], the lower-left entry never read
    var ar = [_]f64{ 4.0, 7.0, 2.0, 3.0 };
    var ai: [4]f64 = @splat(0.0);
    try chol(CxMat.initPacked(&ar, &ai, 2, 2), testing.allocator, testing.io, null);
    // R = [[2, 1], [0, sqrt(2)]]
    try zigstd.testExpectApproxEqAbs(2.0, ar[0], 1e-15);
    try zigstd.testExpectEqual(@as(f64, 0), ar[1]);
    try zigstd.testExpectApproxEqAbs(1.0, ar[2], 1e-15);
    try zigstd.testExpectApproxEqAbs(@sqrt(2.0), ar[3], 1e-15);
}

test "chol of a complex Hermitian matrix; a non positive-definite one fails" {
    // A = [[2, i], [-i, 2]]: R = [[sqrt2, i/sqrt2], [0, sqrt(3/2)]]
    var ar = [_]f64{ 2, 0, 0, 2 };
    var ai = [_]f64{ 0, -1, 1, 0 };
    try chol(CxMat.initPacked(&ar, &ai, 2, 2), testing.allocator, testing.io, null);
    try zigstd.testExpectApproxEqAbs(@sqrt(2.0), ar[0], 1e-15);
    try zigstd.testExpectApproxEqAbs(1.0 / @sqrt(2.0), ai[2], 1e-15);
    try zigstd.testExpectApproxEqAbs(@sqrt(1.5), ar[3], 1e-15);
    try zigstd.testExpectEqual(@as(f64, 0), ai[1]);

    var br = [_]f64{ 1.0, 0.0, 5.0, 1.0 };
    var bi: [4]f64 = @splat(0.0);
    var minor: usize = 0;
    try testing.expectError(error.NotPositiveDefinite, chol(CxMat.initPacked(&br, &bi, 2, 2), testing.allocator, testing.io, &minor));
    try zigstd.testExpectEqual(@as(usize, 2), minor);
}

test "chol of an empty matrix returns the empty factor" {
    var e: [0]f64 = .{};
    try chol(CxMat.initPacked(&e, &e, 0, 0), testing.allocator, testing.io, null);
}
