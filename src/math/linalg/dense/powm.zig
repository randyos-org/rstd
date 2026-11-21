//! Matrix powers (MATLAB's `A^p`, `mpower`): an integer `p` by binary
//! powering (`gemm`), a negative one through `inv`; any other `p` as
//! `expm(p logm(A))`, the principal power, real for a real `A` without
//! negative real eigenvalues. Higham and Lin's Schur-Padé algorithm (SIAM
//! J. Matrix Anal. Appl. 32(3), 2011) is not ported yet. Replaces
//! `matrixfun.powm`, whose eigendecomposition failed on defective
//! matrices (issue #336).

const std = @import("std");
const matrix = @import("../matrix.zig");
const vector = @import("../vector.zig");
const util = @import("util.zig");
const products = @import("products.zig");
const inverse = @import("inverse.zig");
const expm = @import("expm.zig");
const logm = @import("logm.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = matrix.Mat;
const ConstView = matrix.ConstView;
const isComplex = vector.isComplex;
const Real = vector.Real;
const ElementOf = util.ElementOf;
const constOf = util.constOf;
const copyOf = util.copyOf;

pub const PowmError = error{
    NotSquare,
    NoConvergence,
    /// A negative `p` of a singular `A`, or a fractional one of an `A`
    /// with an eigenvalue 0.
    Singular,
    /// A fractional `p` of a real `A` with a negative real eigenvalue:
    /// pass `A` as complex.
    NoRealLogarithm,
} || Allocator.Error || Io.Cancelable;

/// `A^p`, owned by the caller; `p = 0` gives `I`.
pub fn powm(allocator: Allocator, io: Io, a: anytype, p: f64) PowmError!Mat(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const ca = constOf(a);
    if (ca.rows != ca.cols) return error.NotSquare;
    const n = ca.rows;
    if (@trunc(p) != p or @abs(p) > 1 << 53) {
        var arena: std.heap.ArenaAllocator = .init(allocator);
        defer arena.deinit();
        const ar = arena.allocator();
        const l = try logm.logm(ar, io, ca);
        for (l.data) |*e| e.* = scal(T, @floatCast(p), e.*);
        return expm.expm(allocator, io, l) catch |e| switch (e) {
            error.NotSquare => unreachable,
            else => |x| return x,
        };
    }
    var arena: std.heap.ArenaAllocator = .init(allocator);
    defer arena.deinit();
    const ar = arena.allocator();
    var base = if (p < 0) inverse.inv(ar, io, ca) catch |e| switch (e) {
        error.NotSquare => unreachable,
        else => |x| return x,
    } else try copyOf(T, ar, ca);
    var k: u64 = @intFromFloat(@abs(p));
    var acc = try Mat(T).zeros(ar, n, n);
    for (0..n) |i| acc.data[i * (n + 1)] = scal(T, 1, one(T));
    // right to left: acc *= base for each set bit, base squared between
    while (k > 0) : (k >>= 1) {
        if (k & 1 == 1) acc = try mul(T, ar, io, acc, base);
        if (k > 1) base = try mul(T, ar, io, base, base);
    }
    const out = try Mat(T).init(allocator, n, n);
    @memcpy(out.data, acc.data);
    return out;
}

fn mul(comptime T: type, ar: Allocator, io: Io, a: Mat(T), b: Mat(T)) PowmError!Mat(T) {
    return products.mul(ar, io, a, b) catch |e| switch (e) {
        error.DimensionMismatch => unreachable,
        else => |x| return x,
    };
}

fn one(comptime T: type) T {
    return if (comptime isComplex(T)) T.init(1, 0) else 1;
}

fn scal(comptime T: type, c: Real(T), x: T) T {
    return if (comptime isComplex(T)) T.init(c * x.re, c * x.im) else c * x;
}

// ---- tests ----

const testing = std.testing;

test "powm: integer, negative, zero, fractional" {
    const al = testing.allocator;
    const io = testing.io;
    const V = ConstView(f64);
    // [1 1; 0 1]^5 = [1 5; 0 1], exactly
    const j = [_]f64{ 1, 0, 1, 1 };
    const p5 = try powm(al, io, try V.init(&j, 2, 2, 2), 5);
    defer p5.deinit(al);
    try testing.expectEqualSlices(f64, &.{ 1, 0, 5, 1 }, p5.data);
    const m2 = try powm(al, io, try V.init(&j, 2, 2, 2), -2);
    defer m2.deinit(al);
    try testing.expectApproxEqAbs(-2, m2.data[2], 1e-14);
    const id = try powm(al, io, try V.init(&j, 2, 2, 2), 0);
    defer id.deinit(al);
    try testing.expectEqualSlices(f64, &.{ 1, 0, 0, 1 }, id.data);
    // diag(4, 9)^0.5 = diag(2, 3)
    const half = try powm(al, io, try V.init(&[_]f64{ 4, 0, 0, 9 }, 2, 2, 2), 0.5);
    defer half.deinit(al);
    try testing.expectApproxEqAbs(2, half.data[0], 1e-13);
    try testing.expectApproxEqAbs(3, half.data[3], 1e-13);
    try testing.expectError(error.Singular, powm(al, io, try V.init(&[_]f64{ 0, 0, 0, 1 }, 2, 2, 2), -1));
}
