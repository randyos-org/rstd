//! Element-wise operations over views: `add`, `sub`, `times` and `rdivide`
//! (MATLAB's names for `.*` and `./`, as `mul` is the matrix product),
//! `prod`, a `map` of any scalar function (`std.math.complex.sin`, `@round`
//! and so on), `realPart`, `imagPart` and `conj`; each with an `...Into`
//! form. Exact per element, so the inline and BLAS paths agree without a
//! dispatch (issue #326).

const std = @import("std");
const matrix = @import("../matrix.zig");
const vector = @import("../vector.zig");
const util = @import("util.zig");

const Allocator = std.mem.Allocator;
const Mat = matrix.Mat;
const ConstView = matrix.ConstView;
const View = matrix.View;
const isComplex = vector.isComplex;
const Real = vector.Real;
const ElementOf = util.ElementOf;
const constOf = util.constOf;
const viewOf = util.viewOf;

pub const ElementwiseError = error{DimensionMismatch} || Allocator.Error;
pub const ElementwiseIntoError = error{DimensionMismatch};

const Op = enum { add, sub, times, rdivide };

pub fn add(allocator: Allocator, a: anytype, b: anytype) ElementwiseError!Mat(ElementOf(@TypeOf(a))) {
    return binary(.add, allocator, a, b);
}

pub fn sub(allocator: Allocator, a: anytype, b: anytype) ElementwiseError!Mat(ElementOf(@TypeOf(a))) {
    return binary(.sub, allocator, a, b);
}

/// `a .* b`.
pub fn times(allocator: Allocator, a: anytype, b: anytype) ElementwiseError!Mat(ElementOf(@TypeOf(a))) {
    return binary(.times, allocator, a, b);
}

/// `a ./ b`.
pub fn rdivide(allocator: Allocator, a: anytype, b: anytype) ElementwiseError!Mat(ElementOf(@TypeOf(a))) {
    return binary(.rdivide, allocator, a, b);
}

pub fn addInto(out: anytype, a: anytype, b: anytype) ElementwiseIntoError!void {
    return binaryInto(.add, viewOf(out), constOf(a), constOf(b));
}

pub fn subInto(out: anytype, a: anytype, b: anytype) ElementwiseIntoError!void {
    return binaryInto(.sub, viewOf(out), constOf(a), constOf(b));
}

pub fn timesInto(out: anytype, a: anytype, b: anytype) ElementwiseIntoError!void {
    return binaryInto(.times, viewOf(out), constOf(a), constOf(b));
}

pub fn rdivideInto(out: anytype, a: anytype, b: anytype) ElementwiseIntoError!void {
    return binaryInto(.rdivide, viewOf(out), constOf(a), constOf(b));
}

fn binary(comptime op: Op, allocator: Allocator, a: anytype, b: anytype) ElementwiseError!Mat(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const ca = constOf(a);
    const cb: ConstView(T) = constOf(b);
    if (ca.rows != cb.rows or ca.cols != cb.cols) return error.DimensionMismatch;
    const m = try Mat(T).init(allocator, ca.rows, ca.cols);
    binaryFill(op, m.view(), ca, cb);
    return m;
}

fn binaryInto(comptime op: Op, o: anytype, a: anytype, b: anytype) ElementwiseIntoError!void {
    if (a.rows != o.rows or a.cols != o.cols or b.rows != o.rows or b.cols != o.cols) return error.DimensionMismatch;
    binaryFill(op, o, a, b);
}

/// `binaryInto` of views of one shape, which the caller has checked.
fn binaryFill(comptime op: Op, o: anytype, a: anytype, b: anytype) void {
    const T = @TypeOf(o).Element;
    for (0..o.cols) |j| for (o.col(j), a.col(j), b.col(j)) |*r, x, y| {
        r.* = if (comptime isComplex(T)) switch (op) {
            .add => x.add(y),
            .sub => x.sub(y),
            .times => x.mul(y),
            .rdivide => x.div(y),
        } else switch (op) {
            .add => x + y,
            .sub => x - y,
            .times => x * y,
            .rdivide => x / y,
        };
    };
}

/// The product of every element, column by column (MATLAB's `prod` of a
/// vector, or of `A(:)`); 1 when empty.
pub fn prod(a: anytype) ElementOf(@TypeOf(a)) {
    const T = ElementOf(@TypeOf(a));
    const ca = constOf(a);
    var p: T = if (comptime isComplex(T)) T.init(1, 0) else 1;
    for (0..ca.cols) |j| for (ca.col(j)) |x| {
        p = if (comptime isComplex(T)) p.mul(x) else p * x;
    };
    return p;
}

/// `f` of each element: any `fn (T) T`.
pub fn map(allocator: Allocator, a: anytype, f: anytype) Allocator.Error!Mat(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const ca = constOf(a);
    const m = try Mat(T).init(allocator, ca.rows, ca.cols);
    mapFill(m.view(), ca, f);
    return m;
}

pub fn mapInto(out: anytype, a: anytype, f: anytype) ElementwiseIntoError!void {
    const o = viewOf(out);
    const ca = constOf(a);
    if (ca.rows != o.rows or ca.cols != o.cols) return error.DimensionMismatch;
    mapFill(o, ca, f);
}

/// `mapInto` of views of one shape, which the caller has checked.
fn mapFill(o: anytype, a: anytype, f: anytype) void {
    for (0..o.cols) |j| for (o.col(j), a.col(j)) |*r, x| {
        r.* = f(x);
    };
}

/// The real parts, as a real matrix.
pub fn realPart(allocator: Allocator, a: anytype) Allocator.Error!Mat(Real(ElementOf(@TypeOf(a)))) {
    return part(.re, allocator, a);
}

/// The imaginary parts, as a real matrix; zeros for a real `a`.
pub fn imagPart(allocator: Allocator, a: anytype) Allocator.Error!Mat(Real(ElementOf(@TypeOf(a)))) {
    return part(.im, allocator, a);
}

fn part(comptime which: enum { re, im }, allocator: Allocator, a: anytype) Allocator.Error!Mat(Real(ElementOf(@TypeOf(a)))) {
    const T = ElementOf(@TypeOf(a));
    const ca = constOf(a);
    const m = try Mat(Real(T)).init(allocator, ca.rows, ca.cols);
    for (0..ca.cols) |j| for (m.data[j * ca.rows ..][0..ca.rows], ca.col(j)) |*r, x| {
        r.* = if (comptime isComplex(T)) @field(x, @tagName(which)) else if (which == .re) x else 0;
    };
    return m;
}

/// The complex conjugate; a copy for a real `a`.
pub fn conj(allocator: Allocator, a: anytype) Allocator.Error!Mat(ElementOf(@TypeOf(a))) {
    const T = ElementOf(@TypeOf(a));
    const ca = constOf(a);
    const m = try Mat(T).init(allocator, ca.rows, ca.cols);
    for (0..ca.cols) |j| for (m.data[j * ca.rows ..][0..ca.rows], ca.col(j)) |*r, x| {
        r.* = if (comptime isComplex(T)) x.conjugate() else x;
    };
    return m;
}

// ---- tests ----

const testing = std.testing;

test "add, sub, times, rdivide, prod, map, parts and conj" {
    const al = testing.allocator;
    const a = [_]f64{ 1, 2, 3, 4 };
    const b = [_]f64{ 4, 3, 2, 1 };
    const va: ConstView(f64) = try .init(&a, 2, 2, 2);
    const vb: ConstView(f64) = try .init(&b, 2, 2, 2);
    const s = try add(al, va, vb);
    defer s.deinit(al);
    try testing.expectEqualSlices(f64, &.{ 5, 5, 5, 5 }, s.data);
    const d = try sub(al, va, vb);
    defer d.deinit(al);
    try testing.expectEqualSlices(f64, &.{ -3, -1, 1, 3 }, d.data);
    const t = try times(al, va, vb);
    defer t.deinit(al);
    try testing.expectEqualSlices(f64, &.{ 4, 6, 6, 4 }, t.data);
    const q = try rdivide(al, va, vb);
    defer q.deinit(al);
    try testing.expectEqualSlices(f64, &.{ 0.25, 2.0 / 3.0, 1.5, 4 }, q.data);
    try testing.expectEqual(24, prod(va));
    const r = try map(al, va, struct {
        fn f(x: f64) f64 {
            return x * x;
        }
    }.f);
    defer r.deinit(al);
    try testing.expectEqualSlices(f64, &.{ 1, 4, 9, 16 }, r.data);
    try testing.expectError(error.DimensionMismatch, add(al, va, try ConstView(f64).init(&a, 1, 2, 1)));

    const C = std.math.Complex(f64);
    const z = [_]C{ .init(1, 2), .init(3, -4) };
    const vz: ConstView(C) = try .init(&z, 2, 1, 2);
    const im = try imagPart(al, vz);
    defer im.deinit(al);
    try testing.expectEqualSlices(f64, &.{ 2, -4 }, im.data);
    const cz = try conj(al, vz);
    defer cz.deinit(al);
    try testing.expectEqual(C.init(3, 4), cz.data[1]);
    const e = try map(al, vz, std.math.complex.exp);
    defer e.deinit(al);
    try testing.expectApproxEqAbs(std.math.exp(1.0) * @cos(2.0), e.data[0].re, 1e-12);
}
