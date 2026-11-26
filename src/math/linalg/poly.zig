//! Polynomials as coefficient slices, highest power first (MATLAB's
//! order): `roots` by the companion matrix's eigenvalues, `poly` from the
//! roots, `charPoly` of a square matrix. Replaces matlab88's `POLY` and
//! `ROOTS` (`matlab/matrixfun.zig`, issue #328).

const std = @import("std");
const matrix = @import("matrix.zig");
const vector = @import("vector.zig");
const linalg = @import("__root__.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = matrix.Mat;
const ConstView = matrix.ConstView;
const isComplex = vector.isComplex;
const Real = vector.Real;

pub fn Complex(comptime T: type) type {
    return std.math.Complex(Real(T));
}

pub const RootsError = error{NoConvergence} || Allocator.Error || Io.Cancelable;

/// MATLAB's `roots`: the roots of `c[0] xⁿ + ... + c[n]`, owned by the
/// caller. Leading zeros are dropped, trailing ones are roots at 0; the
/// rest are the eigenvalues of the companion matrix (`eig`, so `geev` with
/// balancing). An all-zero `c` has no roots.
pub fn roots(allocator: Allocator, io: Io, comptime T: type, c: []const T) RootsError![]Complex(T) {
    const zero = scalar(T, 0);
    const first = for (c, 0..) |x, i| {
        if (!eq(T, x, zero)) break i;
    } else return allocator.alloc(Complex(T), 0);
    var last = c.len;
    while (eq(T, c[last - 1], zero)) last -= 1;
    const nz = c.len - last;
    const p = c[first..last];
    const d = p.len - 1;
    const out = try allocator.alloc(Complex(T), d + nz);
    errdefer allocator.free(out);
    for (out[d..]) |*r| r.* = .init(0, 0);
    if (d == 0) return out;

    // first row -p[1..] / p[0], ones below the diagonal
    const a = try Mat(T).zeros(allocator, d, d);
    defer a.deinit(allocator);
    for (0..d) |j| a.data[j * d] = neg(T, div(T, p[j + 1], p[0]));
    for (1..d) |i| a.data[i + (i - 1) * d] = scalar(T, 1);
    const e = linalg.eig(allocator, io, a, false) catch |err| return switch (err) {
        error.NotSquare => unreachable,
        else => |x| x,
    };
    defer e.deinit(allocator);
    switch (e) {
        .hermitian => |h| for (out[0..d], h.values) |*r, v| {
            r.* = .init(v, 0);
        },
        .general => |g| @memcpy(out[0..d], g.values),
    }
    return out;
}

/// MATLAB's `poly` of a vector: the monic polynomial with these roots,
/// `r.len + 1` coefficients, owned by the caller.
pub fn poly(allocator: Allocator, comptime T: type, r: []const T) Allocator.Error![]T {
    const c = try allocator.alloc(T, r.len + 1);
    @memset(c, scalar(T, 0));
    c[0] = scalar(T, 1);
    for (r, 1..) |x, j| {
        var k = j;
        while (k >= 1) : (k -= 1) c[k] = sub(T, c[k], mul(T, x, c[k - 1]));
    }
    return c;
}

pub const CharPolyError = error{ NotSquare, NoConvergence } || Allocator.Error || Io.Cancelable;

/// MATLAB's `poly` of a square matrix: its characteristic polynomial,
/// `det(x I - A)`, from its eigenvalues; a real `A`'s is real, as its
/// complex roots come in conjugate pairs.
pub fn charPoly(allocator: Allocator, io: Io, comptime T: type, a: ConstView(T)) CharPolyError![]T {
    const e = try linalg.eig(allocator, io, a, false);
    defer e.deinit(allocator);
    const C = Complex(T);
    const vals: []const C = switch (e) {
        .general => |g| g.values,
        .hermitian => |h| blk: {
            const v = try allocator.alloc(C, h.values.len);
            for (v, h.values) |*d, x| d.* = .init(x, 0);
            break :blk v;
        },
    };
    defer if (e == .hermitian) allocator.free(vals);
    const cc = try poly(allocator, C, vals);
    if (comptime isComplex(T)) return cc;
    defer allocator.free(cc);
    const out = try allocator.alloc(T, cc.len);
    for (out, cc) |*o, x| o.* = x.re;
    return out;
}

fn scalar(comptime T: type, x: Real(T)) T {
    return if (comptime isComplex(T)) T.init(x, 0) else x;
}

fn eq(comptime T: type, x: T, y: T) bool {
    return if (comptime isComplex(T)) x.re == y.re and x.im == y.im else x == y;
}

fn neg(comptime T: type, x: T) T {
    return if (comptime isComplex(T)) x.neg() else -x;
}

fn sub(comptime T: type, x: T, y: T) T {
    return if (comptime isComplex(T)) x.sub(y) else x - y;
}

fn mul(comptime T: type, x: T, y: T) T {
    return if (comptime isComplex(T)) x.mul(y) else x * y;
}

fn div(comptime T: type, x: T, y: T) T {
    return if (comptime isComplex(T)) x.div(y) else x / y;
}

// ---- tests ----

const testing = std.testing;

/// `got` holds `want`'s values, in any order, within `tol`.
fn expectSameRoots(want: []const Complex(f64), got: []const Complex(f64), tol: f64) !void {
    try testing.expectEqual(want.len, got.len);
    var used: [16]bool = @splat(false);
    for (want) |w| {
        const k = for (got, 0..) |g, i| {
            if (!used[i] and @abs(g.re - w.re) < tol and @abs(g.im - w.im) < tol) break i;
        } else return error.TestExpectedEqual;
        used[k] = true;
    }
}

test "roots: real, complex pairs, leading and trailing zeros" {
    const al = testing.allocator;
    const io = testing.io;
    // (x - 1)(x - 2)(x - 3), with a leading zero and a root at 0
    const r = try roots(al, io, f64, &.{ 0, 1, -6, 11, -6, 0 });
    defer al.free(r);
    try expectSameRoots(&.{ .init(1, 0), .init(2, 0), .init(3, 0), .init(0, 0) }, r, 1e-10);
    // x² + 1
    const z = try roots(al, io, f64, &.{ 1, 0, 1 });
    defer al.free(z);
    try expectSameRoots(&.{ .init(0, 1), .init(0, -1) }, z, 1e-12);
    // x - 5: a 1 by 1 companion matrix
    const one = try roots(al, io, f64, &.{ 2, -10 });
    defer al.free(one);
    try expectSameRoots(&.{.init(5, 0)}, one, 1e-12);
    const none = try roots(al, io, f64, &.{ 0, 0 });
    defer al.free(none);
    try testing.expectEqual(0, none.len);
}

test "poly and charPoly invert roots and eig" {
    const al = testing.allocator;
    const io = testing.io;
    const c = try poly(al, f64, &.{ 1, 2, 3 });
    defer al.free(c);
    try testing.expectEqualSlices(f64, &.{ 1, -6, 11, -6 }, c);
    // MATLAB's poly([1 2; 3 4]) is [1 -5 -2]
    const a = [_]f64{ 1, 3, 2, 4 };
    const p = try charPoly(al, io, f64, try .init(&a, 2, 2, 2));
    defer al.free(p);
    for ([_]f64{ 1, -5, -2 }, p) |w, g| try testing.expectApproxEqAbs(w, g, 1e-12);
    // a rotation's: x² + 1, real from complex eigenvalues
    const rot = [_]f64{ 0, 1, -1, 0 };
    const q = try charPoly(al, io, f64, try .init(&rot, 2, 2, 2));
    defer al.free(q);
    for ([_]f64{ 1, 0, 1 }, q) |w, g| try testing.expectApproxEqAbs(w, g, 1e-12);
}
