//! Orthonormal bases from the SVD (MATLAB's `orth` and `null`; `null` is a
//! Zig keyword, hence `nullSpace`), with `rank`'s default tolerance. They
//! replace `CxMat.orth`, whose unpivoted QR did not reveal rank (issue #323).

const std = @import("std");
const matrix = @import("../matrix.zig");
const vector = @import("../vector.zig");
const singular = @import("singular.zig");
const util = @import("util.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = matrix.Mat;
const Real = vector.Real;
const ElementOf = util.ElementOf;
const constOf = util.constOf;
const SvdError = singular.SvdError;

/// An orthonormal basis of `a`'s range: the first `rank` left singular
/// vectors, `m` by `rank`, owned by the caller. `tol` as `rank`'s.
pub fn orth(allocator: Allocator, io: Io, a: anytype, tol: ?Real(ElementOf(@TypeOf(a)))) SvdError!Mat(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const ca = constOf(a);
    const f = try singular.svd(allocator, io, ca, .economy);
    defer f.deinit(allocator);
    const r = singular.rankOf(Real(T), f.s, ca.rows, ca.cols, tol);
    const q = try Mat(T).init(allocator, ca.rows, r);
    @memcpy(q.data, f.u.?.data[0 .. ca.rows * r]);
    return q;
}

/// An orthonormal basis of `a`'s null space: the last `n - rank` right
/// singular vectors of the full SVD, `n` by `n - rank`, owned by the
/// caller. `tol` as `rank`'s.
pub fn nullSpace(allocator: Allocator, io: Io, a: anytype, tol: ?Real(ElementOf(@TypeOf(a)))) SvdError!Mat(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const ca = constOf(a);
    const n = ca.cols;
    const f = try singular.svd(allocator, io, ca, .full);
    defer f.deinit(allocator);
    const r = singular.rankOf(Real(T), f.s, ca.rows, n, tol);
    const z = try Mat(T).init(allocator, n, n - r);
    @memcpy(z.data, f.v.?.data[r * n .. n * n]);
    return z;
}

// ---- tests ----

const testing = std.testing;

test "orth and nullSpace: rank-deficient, orthonormal, A Z = 0" {
    const al = testing.allocator;
    const io = testing.io;
    // rank 1: [1 2; 2 4; 3 6]
    const a = [_]f64{ 1, 2, 3, 2, 4, 6 };
    const v: matrix.ConstView(f64) = try .init(&a, 3, 2, 3);
    const q = try orth(al, io, v, null);
    defer q.deinit(al);
    try testing.expectEqual(1, q.cols);
    var nq: f64 = 0;
    for (q.data) |x| nq += x * x;
    try testing.expectApproxEqAbs(1, nq, 1e-12);
    // q is ± [1 2 3] / √14
    try testing.expectApproxEqAbs(@as(f64, 2) / @sqrt(14.0), @abs(q.data[1]), 1e-12);

    const z = try nullSpace(al, io, v, null);
    defer z.deinit(al);
    try testing.expectEqual(1, z.cols);
    for (0..3) |i| try testing.expectApproxEqAbs(0, a[i] * z.data[0] + a[i + 3] * z.data[1], 1e-12);

    // full rank: an empty null space
    const id = [_]f64{ 1, 0, 0, 1 };
    const e = try nullSpace(al, io, try matrix.ConstView(f64).init(&id, 2, 2, 2), null);
    defer e.deinit(al);
    try testing.expectEqual(0, e.cols);
}
