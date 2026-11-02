//! `rstd.linalg`'s run-time-shape functions (design.md, "Operations over
//! views"), re-exported flat at `rstd.linalg`. Each operand is a `View` or
//! `ConstView`, a `Mat`, a pointer to a `StaticMat` or a slice (a column),
//! in any mix; a vector is a view of one column or one row. Shapes at or
//! below the inline limit take the fixed-size forms' operation order, so a
//! `Mat` and a `StaticMat` holding the same numbers give the same bits; past
//! it, the BLAS. Outputs must not overlap inputs.

const std = @import("std");
const Allocator = std.mem.Allocator;
const Io = std.Io;
const blas = @import("blas/__root__.zig");
const vector = @import("vector.zig");
const matrix = @import("matrix.zig");

const View = matrix.View;
const ConstView = matrix.ConstView;
const Mat = matrix.Mat;
const isComplex = vector.isComplex;
const Real = vector.Real;
const eZero = matrix.eZero;
const eMul = matrix.eMul;
const eMulAdd = matrix.eMulAdd;

/// The inline limit `L` per dimension (design.md, Dispatch); M11 makes it
/// `rstd_linalg_options.inline_limit`.
const inline_limit = 8;

/// The element type of an operand.
fn ElementOf(comptime X: type) type {
    switch (@typeInfo(X)) {
        .pointer => |p| switch (p.size) {
            .slice => return p.child,
            .one => return if (@typeInfo(p.child) == .array) @typeInfo(p.child).array.child else p.child.Element,
            else => {},
        },
        .@"struct" => if (@hasField(X, "data") and @typeInfo(@FieldType(X, "data")) == .pointer) return X.Element,
        else => {},
    }
    @compileError("expects a View, ConstView, Mat, pointer to a StaticMat or slice, found " ++ @typeName(X));
}

/// `x` read-only.
fn constOf(x: anytype) ConstView(ElementOf(@TypeOf(x))) {
    const T = ElementOf(@TypeOf(x));
    return switch (@typeInfo(@TypeOf(x))) {
        .pointer => |p| if (p.size == .slice)
            column(T, true, x)
        else if (@typeInfo(p.child) == .array) column(T, true, x[0..]) else x.constView(),
        else => x.constView(),
    };
}

/// `x` writable.
fn viewOf(x: anytype) View(ElementOf(@TypeOf(x))) {
    const X = @TypeOf(x);
    const T = ElementOf(X);
    if (X == View(T)) return x;
    if (X == ConstView(T)) @compileError("expects a writable operand, found a ConstView");
    return switch (@typeInfo(X)) {
        .pointer => |p| if (p.size == .slice)
            column(T, false, x)
        else if (@typeInfo(p.child) == .array) column(T, false, x[0..]) else x.view(),
        else => x.view(),
    };
}

fn column(comptime T: type, comptime read_only: bool, s: if (read_only) []const T else []T) if (read_only) ConstView(T) else View(T) {
    return .{ .data = s, .rows = s.len, .cols = 1, .ld = @max(s.len, 1) };
}

/// A vector's elements: `n` of them `inc` apart from `ptr`.
fn Strided(comptime T: type) type {
    return struct {
        ptr: [*]T,
        n: usize,
        inc: usize,

        fn at(s: @This(), i: usize) T {
            return s.ptr[i * s.inc];
        }

        /// The first `c` elements as an array, for the fixed-size forms.
        fn gather(s: @This(), comptime c: usize) [c]T {
            var r: [c]T = undefined;
            for (&r, 0..) |*x, i| x.* = s.at(i);
            return r;
        }
    };
}

/// `linalg.vector`'s forms at length `c`, one non-inline function per size so
/// that each gets its own comptime branch quota.
fn Fixed(comptime c: usize) type {
    return struct {
        fn dot(x: anytype, y: @TypeOf(x)) @TypeOf(x.ptr[0]) {
            return vector.dot(x.gather(c), y.gather(c));
        }

        fn length(x: anytype) Real(@TypeOf(x.ptr[0])) {
            return vector.length(x.gather(c));
        }

        fn sum(x: anytype) @TypeOf(x.ptr[0]) {
            return vector.sum(x.gather(c));
        }

        fn extreme(comptime which: blas.Extreme, x: anytype) @TypeOf(x.ptr[0]) {
            return if (which == .max) vector.max(x.gather(c)) else vector.min(x.gather(c));
        }
    };
}

/// `v` as a vector: one column, or one row (stride `ld`).
fn strided(v: anytype) error{DimensionMismatch}!Strided(@TypeOf(v).Element) {
    const p: [*]@TypeOf(v).Element = @constCast(v.data.ptr);
    if (v.cols == 1 or v.rows == 0 or v.cols == 0) return .{ .ptr = p, .n = if (v.cols == 0) 0 else v.rows, .inc = 1 };
    if (v.rows == 1) return .{ .ptr = p, .n = v.cols, .inc = v.ld };
    return error.DimensionMismatch;
}

fn len(n: usize) isize {
    return @intCast(n);
}

/// A BLAS call's error, of which only `Canceled` can arise here: the
/// arguments come from views, so the BLAS's argument errors are bugs.
fn canceled(e: anyerror) Io.Cancelable {
    return switch (e) {
        error.Canceled => error.Canceled,
        else => unreachable,
    };
}

/// As `canceled`, for the calls that may also allocate.
fn allocOrCanceled(e: anyerror) (Allocator.Error || Io.Cancelable) {
    return switch (e) {
        error.Canceled => error.Canceled,
        error.OutOfMemory => error.OutOfMemory,
        else => unreachable,
    };
}

pub const VectorError = error{DimensionMismatch} || Io.Cancelable;

/// `xᴴ y` (`Σ conj(x[i]) y[i]`, MATLAB's `dot`); `linalg.vector.dot`'s order
/// up to the inline limit, else the BLAS's `dot` or `dotc`.
pub fn dot(io: Io, x: anytype, y: anytype) VectorError!ElementOf(@TypeOf(x)) {
    const T = ElementOf(@TypeOf(x));
    const sx = try strided(constOf(x));
    const sy = try strided(constOf(y));
    if (sx.n != sy.n) return error.DimensionMismatch;
    switch (sx.n) {
        0 => return eZero(T),
        inline 1...inline_limit => |c| return Fixed(c).dot(sx, sy),
        else => {},
    }
    const f = if (comptime isComplex(T)) blas.dotc else blas.dot;
    return f(io, T, len(sx.n), sx.ptr, len(sx.inc), sy.ptr, len(sy.inc)) catch |e| return canceled(e);
}

/// The 2-norm (MATLAB's `norm` of a vector; `linalg.vector`'s `norm`
/// normalizes): `linalg.vector.length` up to the inline limit, else `nrm2`.
pub fn norm(io: Io, x: anytype) VectorError!Real(ElementOf(@TypeOf(x))) {
    const T = ElementOf(@TypeOf(x));
    const s = try strided(constOf(x));
    switch (s.n) {
        0 => return 0,
        inline 1...inline_limit => |c| return Fixed(c).length(s),
        else => {},
    }
    return blas.nrm2(io, T, len(s.n), s.ptr, len(s.inc)) catch |e| return canceled(e);
}

/// `Σ x[i]`: in order up to the inline limit, else the BLAS's `sum`, the
/// real and imaginary parts apart for complex.
pub fn sum(io: Io, x: anytype) VectorError!ElementOf(@TypeOf(x)) {
    const T = ElementOf(@TypeOf(x));
    const s = try strided(constOf(x));
    switch (s.n) {
        0 => return eZero(T),
        inline 1...inline_limit => |c| return Fixed(c).sum(s),
        else => {},
    }
    if (comptime !isComplex(T)) return blas.sum(io, T, len(s.n), s.ptr, len(s.inc)) catch |e| return canceled(e);
    const F = Real(T);
    const parts: [*]F = @ptrCast(s.ptr);
    const inc = len(2 * s.inc);
    return .{
        .re = blas.sum(io, F, len(s.n), parts, inc) catch |e| return canceled(e),
        .im = blas.sum(io, F, len(s.n), parts + 1, inc) catch |e| return canceled(e),
    };
}

/// The largest element, or the first NaN; ties keep the first; 0 for an
/// empty vector, as the BLAS's `max`. Real only.
pub fn max(io: Io, x: anytype) VectorError!ElementOf(@TypeOf(x)) {
    return extreme(.max, io, x);
}

/// The smallest element, or the first NaN, as `max`.
pub fn min(io: Io, x: anytype) VectorError!ElementOf(@TypeOf(x)) {
    return extreme(.min, io, x);
}

fn extreme(comptime which: blas.Extreme, io: Io, x: anytype) VectorError!ElementOf(@TypeOf(x)) {
    const T = ElementOf(@TypeOf(x));
    if (comptime isComplex(T)) @compileError("max and min order real elements only");
    const s = try strided(constOf(x));
    switch (s.n) {
        0 => return 0,
        inline 1...inline_limit => |c| return Fixed(c).extreme(which, s),
        else => {},
    }
    return blas.extreme(io, T, which, len(s.n), s.ptr, len(s.inc)) catch |e| return canceled(e);
}

/// `y += a x`, by the BLAS's `axpy`.
pub fn axpy(io: Io, a: anytype, x: anytype, y: anytype) VectorError!void {
    const T = ElementOf(@TypeOf(y));
    const sx = try strided(constOf(x));
    const sy = try strided(viewOf(y).constView());
    if (sx.n != sy.n) return error.DimensionMismatch;
    const alpha: T = a;
    blas.axpy(io, T, len(sx.n), alpha, sx.ptr, len(sx.inc), sy.ptr, len(sy.inc)) catch |e| return canceled(e);
}

/// `x *= a`, element by element as `linalg.matrix.scale` (design.md,
/// Dispatch: element-wise forms are inline at every size).
pub fn scale(a: anytype, x: anytype) void {
    const v = viewOf(x);
    const alpha: @TypeOf(v).Element = a;
    for (0..v.cols) |j| for (v.col(j)) |*e| {
        e.* = eMul(e.*, alpha);
    };
}

pub const MulError = error{DimensionMismatch} || Allocator.Error || Io.Cancelable;

/// `a b` (MATLAB's `*`), owned by the caller.
pub fn mul(allocator: Allocator, io: Io, a: anytype, b: anytype) MulError!Mat(ElementOf(@TypeOf(a))) {
    const ca = constOf(a);
    const cb = constOf(b);
    if (ca.cols != cb.rows) return error.DimensionMismatch;
    const r = try Mat(ElementOf(@TypeOf(a))).init(allocator, ca.rows, cb.cols);
    errdefer r.deinit(allocator);
    try mulInto(allocator, io, r.view(), ca, cb);
    return r;
}

/// `out = a b`: up to the inline limit in every dimension, `linalg.matrix.mul`'s
/// order (each column of `out` the columns of `a` scaled by `b`'s and
/// summed in order); else `gemv` for one column, `gemm` otherwise.
pub fn mulInto(allocator: Allocator, io: Io, out: anytype, a: anytype, b: anytype) MulError!void {
    const T = ElementOf(@TypeOf(a));
    const o = viewOf(out);
    const ca = constOf(a);
    const cb = constOf(b);
    if (ca.cols != cb.rows or o.rows != ca.rows or o.cols != cb.cols) return error.DimensionMismatch;
    const m = ca.rows;
    const n = cb.cols;
    const k = ca.cols;
    if (m <= inline_limit and n <= inline_limit and k <= inline_limit) {
        for (0..n) |j| for (0..m) |i| {
            var acc = eZero(T);
            if (k > 0) {
                acc = eMul(ca.at(i, 0), cb.at(0, j));
                for (1..k) |p| acc = eMulAdd(ca.at(i, p), cb.at(p, j), acc);
            }
            o.data[i + j * o.ld] = acc;
        };
        return;
    }
    const pa: [*]T = @constCast(ca.data.ptr);
    const pb: [*]T = @constCast(cb.data.ptr);
    const one: T = if (comptime isComplex(T)) T.init(1, 0) else 1;
    if (n == 1) {
        blas.gemv(allocator, io, T, .n, len(m), len(k), one, pa, len(ca.ld), pb, 1, eZero(T), o.data.ptr, 1) catch |e| return allocOrCanceled(e);
        return;
    }
    blas.gemm(allocator, io, T, .n, .n, len(m), len(n), len(k), one, pa, len(ca.ld), pb, len(cb.ld), eZero(T), o.data.ptr, len(o.ld)) catch |e| return allocOrCanceled(e);
}

pub const TransposeError = error{DimensionMismatch} || Io.Cancelable;

/// `aᵀ` (MATLAB's `a.'`), owned by the caller.
pub fn transpose(allocator: Allocator, io: Io, a: anytype) (Allocator.Error || Io.Cancelable)!Mat(ElementOf(@TypeOf(a))) {
    return transposed(.t, allocator, io, a);
}

/// `aᴴ` (MATLAB's `a'`); `aᵀ` for real `T`.
pub fn ctranspose(allocator: Allocator, io: Io, a: anytype) (Allocator.Error || Io.Cancelable)!Mat(ElementOf(@TypeOf(a))) {
    return transposed(.c, allocator, io, a);
}

fn transposed(comptime trans: blas.Trans, allocator: Allocator, io: Io, a: anytype) (Allocator.Error || Io.Cancelable)!Mat(ElementOf(@TypeOf(a))) {
    const ca = constOf(a);
    const r = try Mat(ElementOf(@TypeOf(a))).init(allocator, ca.cols, ca.rows);
    errdefer r.deinit(allocator);
    transposeWith(trans, io, r.view(), ca) catch |e| return switch (e) {
        error.DimensionMismatch => unreachable,
        error.Canceled => error.Canceled,
    };
    return r;
}

/// `out = aᵀ`, by the BLAS's `omatcopy`.
pub fn transposeInto(io: Io, out: anytype, a: anytype) TransposeError!void {
    return transposeWith(.t, io, viewOf(out), constOf(a));
}

/// `out = aᴴ`, by the BLAS's `omatcopy`.
pub fn ctransposeInto(io: Io, out: anytype, a: anytype) TransposeError!void {
    return transposeWith(.c, io, viewOf(out), constOf(a));
}

fn transposeWith(comptime trans: blas.Trans, io: Io, o: anytype, a: anytype) TransposeError!void {
    const T = @TypeOf(a).Element;
    if (o.rows != a.cols or o.cols != a.rows) return error.DimensionMismatch;
    if (a.rows == 0 or a.cols == 0) return;
    const one: T = if (comptime isComplex(T)) T.init(1, 0) else 1;
    const t: blas.Trans = if (comptime isComplex(T)) trans else .t;
    blas.omatcopy(io, T, t, len(a.rows), len(a.cols), one, a.data.ptr, len(a.ld), o.data.ptr, len(o.ld)) catch |e| return canceled(e);
}

const zigstd = @import("../zigstd.zig");
const testing = std.testing;
const C64 = std.math.Complex(f64);
const test_io = Io.failing;

test "views: init's checks, sub, col; StaticMat and Mat views share storage" {
    var buf = [_]f64{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    try testing.expectError(error.DimensionMismatch, View(f64).init(&buf, 3, 2, 2));
    try testing.expectError(error.SliceTooShort, View(f64).init(&buf, 3, 4, 3));
    const v = try View(f64).init(&buf, 3, 3, 3);
    try testing.expectEqual(@as(f64, 8), v.at(1, 2));
    const s = try v.sub(1, 1, 2, 2);
    try testing.expectEqual(@as(f64, 5), s.at(0, 0));
    try testing.expectEqual(@as(f64, 9), s.at(1, 1));
    try testing.expectEqualSlices(f64, &.{ 5, 6 }, s.col(0));
    try testing.expectError(error.DimensionMismatch, v.sub(2, 0, 2, 1));

    var sm = matrix.fromCols(matrix.M2f64, .{ .{ 1, 2 }, .{ 3, 4 } });
    sm.view().data[3] = 40;
    try testing.expectEqual(@as(f64, 40), sm.data[1][1]);
    try testing.expectEqual(@as(f64, 3), sm.constView().at(0, 1));
    const m = try Mat(f64).zeros(testing.allocator, 2, 3);
    defer m.deinit(testing.allocator);
    m.view().data[5] = 7;
    try testing.expectEqual(@as(f64, 7), m.at(1, 2));
}

fn randomElem(comptime T: type, rnd: std.Random) T {
    if (comptime isComplex(T)) return .init(rnd.float(f64) * 2 - 1, rnd.float(f64) * 2 - 1);
    return std.math.ldexp(rnd.float(T) * 2 - 1, rnd.intRangeAtMost(i32, -4, 4));
}

fn fill(comptime T: type, rnd: std.Random, s: []T) void {
    for (s) |*x| x.* = randomElem(T, rnd);
}

fn sameBits(comptime T: type, a: []const T, b: []const T) bool {
    return std.mem.eql(u8, std.mem.sliceAsBytes(a), std.mem.sliceAsBytes(b));
}

/// `mulInto` on a `Mat` against `linalg.matrix.mul` on `StaticMat`s of the
/// same numbers, `R` by `K` times `K` by `C`.
fn expectStaticBits(comptime T: type, comptime r: usize, comptime k: usize, comptime c: usize, rnd: std.Random) !void {
    const al = testing.allocator;
    for (0..20) |_| {
        var a: matrix.StaticMat(r, k, T) = undefined;
        var b: matrix.StaticMat(k, c, T) = undefined;
        fill(T, rnd, a.view().data);
        fill(T, rnd, b.view().data);
        const want = matrix.mul(a, b);
        const got = try mul(al, test_io, &a, &b);
        defer got.deinit(al);
        var w = want;
        try testing.expect(sameBits(T, w.view().data, got.data));
        // a column: `M v`
        const v: [k]T = b.data[0];
        var col: [r]T = undefined;
        try mulInto(al, test_io, &col, &a, &v);
        try testing.expect(sameBits(T, &col, &@as([r]T, w.data[0])));
    }
}

test "up to the inline limit, a Mat multiplies to a StaticMat's bits" {
    var prng: std.Random.DefaultPrng = .init(0x7e57);
    const shapes = .{ .{ 1, 1, 1 }, .{ 2, 3, 4 }, .{ 4, 4, 4 }, .{ 3, 8, 5 }, .{ 8, 8, 8 } };
    inline for (.{ f32, f64, C64 }) |T| inline for (shapes) |s| try expectStaticBits(T, s[0], s[1], s[2], prng.random());
}

/// `Σ_p a(i, p) b(p, j)` in `f64` (or its complex), against `got` within
/// `k` roundings of `T`.
fn expectNearProduct(comptime T: type, a: ConstView(T), b: ConstView(T), got: ConstView(T)) !void {
    const eps: f64 = std.math.floatEps(Real(T));
    for (0..got.cols) |j| for (0..got.rows) |i| {
        var re: f64 = 0;
        var im: f64 = 0;
        var mag: f64 = 0;
        for (0..a.cols) |p| {
            const x = a.at(i, p);
            const y = b.at(p, j);
            if (comptime isComplex(T)) {
                re += x.re * y.re - x.im * y.im;
                im += x.re * y.im + x.im * y.re;
                mag += (@abs(x.re) + @abs(x.im)) * (@abs(y.re) + @abs(y.im));
            } else {
                re += @as(f64, x) * y;
                mag += @abs(@as(f64, x) * y);
            }
        }
        const g = got.at(i, j);
        const tol = 2 * @as(f64, @floatFromInt(a.cols + 1)) * eps * mag;
        if (comptime isComplex(T)) {
            try testing.expect(@abs(g.re - re) <= tol and @abs(g.im - im) <= tol);
        } else try testing.expect(@abs(@as(f64, g) - re) <= tol);
    };
}

fn expectBlasProduct(comptime T: type, m: usize, n: usize, k: usize, rnd: std.Random) !void {
    const al = testing.allocator;
    const a = try Mat(T).init(al, m, k);
    defer a.deinit(al);
    const b = try Mat(T).init(al, k, n);
    defer b.deinit(al);
    fill(T, rnd, a.data);
    fill(T, rnd, b.data);
    // undefined storage may hold NaNs: `beta = 0` must not read them
    const out = try Mat(T).init(al, m, n);
    defer out.deinit(al);
    const nan = std.math.nan(Real(T));
    @memset(out.data, if (comptime isComplex(T)) T.init(nan, nan) else nan);
    try mulInto(al, test_io, out, a, b);
    try expectNearProduct(T, a.constView(), b.constView(), out.constView());
}

test "past the inline limit, mul is the BLAS's gemm or gemv" {
    var prng: std.Random.DefaultPrng = .init(0xb1a5);
    inline for (.{ f32, f64, C64 }) |T| {
        try expectBlasProduct(T, 9, 7, 5, prng.random());
        try expectBlasProduct(T, 3, 4, 20, prng.random());
        try expectBlasProduct(T, 40, 33, 17, prng.random());
        try expectBlasProduct(T, 30, 1, 12, prng.random());
    }
}

fn expectVectorOps(comptime T: type, n: usize, rnd: std.Random) !void {
    const al = testing.allocator;
    // `x` is a row of a 3-row matrix, stride 3; `y` a column
    const xs = try al.alloc(T, 3 * n);
    defer al.free(xs);
    fill(T, rnd, xs);
    const xv = try ConstView(T).init(xs, 1, n, 3);
    const y = try al.alloc(T, n);
    defer al.free(y);
    fill(T, rnd, y);
    var d_re: f64 = 0;
    var d_im: f64 = 0;
    var s_re: f64 = 0;
    var s_im: f64 = 0;
    var sq: f64 = 0;
    var mag: f64 = 0;
    for (0..n) |i| {
        const x = xv.at(0, i);
        if (comptime isComplex(T)) {
            d_re += x.re * y[i].re + x.im * y[i].im;
            d_im += x.re * y[i].im - x.im * y[i].re;
            s_re += x.re;
            s_im += x.im;
            sq += x.re * x.re + x.im * x.im;
            mag += (@abs(x.re) + @abs(x.im)) * (@abs(y[i].re) + @abs(y[i].im));
        } else {
            d_re += @as(f64, x) * y[i];
            s_re += x;
            sq += @as(f64, x) * x;
            mag += @abs(@as(f64, x) * y[i]);
        }
    }
    const eps: f64 = std.math.floatEps(Real(T));
    const tol = 2 * @as(f64, @floatFromInt(n + 1)) * eps;
    const d = try dot(test_io, xv, y);
    const s = try sum(test_io, xv);
    const nrm: f64 = try norm(test_io, xv);
    try testing.expect(@abs(nrm - @sqrt(sq)) <= tol * @sqrt(sq));
    if (comptime isComplex(T)) {
        try testing.expect(@abs(d.re - d_re) <= tol * mag and @abs(d.im - d_im) <= tol * mag);
        try testing.expect(@abs(s.re - s_re) <= tol * @as(f64, @floatFromInt(n)) and @abs(s.im - s_im) <= tol * @as(f64, @floatFromInt(n)));
    } else {
        try testing.expect(@abs(@as(f64, d) - d_re) <= tol * mag);
        try testing.expect(@abs(@as(f64, s) - s_re) <= tol * @as(f64, @floatFromInt(n)));
        var hi = xv.at(0, 0);
        var lo = hi;
        for (1..n) |i| {
            hi = @max(hi, xv.at(0, i));
            lo = @min(lo, xv.at(0, i));
        }
        try testing.expectEqual(hi, try max(test_io, xv));
        try testing.expectEqual(lo, try min(test_io, xv));
    }
}

test "dot, norm, sum, max and min, inline and on the BLAS, on rows and columns" {
    var prng: std.Random.DefaultPrng = .init(0xd07);
    inline for (.{ f32, f64, C64 }) |T| for ([_]usize{ 1, 3, 8, 9, 100 }) |n| try expectVectorOps(T, n, prng.random());
}

test "up to the inline limit, vector reductions keep linalg.vector's bits" {
    var prng: std.Random.DefaultPrng = .init(0x5a5e);
    var a: [6]f64 = undefined;
    var b: [6]f64 = undefined;
    fill(f64, prng.random(), &a);
    fill(f64, prng.random(), &b);
    try testing.expectEqual(vector.dot(a, b), try dot(test_io, &a, &b));
    try testing.expectEqual(vector.length(a), try norm(test_io, &a));
    try testing.expectEqual(vector.sum(a), try sum(test_io, &a));
    var z: [5]C64 = undefined;
    fill(C64, prng.random(), &z);
    try testing.expectEqual(vector.dot(z, z), try dot(test_io, &z, &z));
}

test "transposes, axpy and scale" {
    const al = testing.allocator;
    var a = matrix.fromRows(matrix.StaticMat(2, 3, f64), .{ .{ 1, 2, 3 }, .{ 4, 5, 6 } });
    const t = try transpose(al, test_io, &a);
    defer t.deinit(al);
    try testing.expectEqualSlices(f64, &.{ 1, 2, 3, 4, 5, 6 }, t.data);
    var too_small: [4]f64 = undefined;
    try testing.expectError(error.DimensionMismatch, transposeInto(test_io, &too_small, &a));

    var z: matrix.StaticMat(1, 2, C64) = .{ .data = .{ .{.init(1, 2)}, .{.init(3, -4)} } };
    const zh = try ctranspose(al, test_io, &z);
    defer zh.deinit(al);
    try testing.expectEqual(C64.init(1, -2), zh.data[0]);
    try testing.expectEqual(C64.init(3, 4), zh.data[1]);

    var y = [_]f64{ 1, 1, 1 };
    try axpy(test_io, 2, &[_]f64{ 1, 2, 3 }, &y);
    try testing.expectEqualSlices(f64, &.{ 3, 5, 7 }, &y);
    try testing.expectError(error.DimensionMismatch, axpy(test_io, 2, &[_]f64{ 1, 2 }, &y));
    scale(-1, &a);
    try testing.expectEqual(@as(f64, -6), a.data[2][1]);
    try testing.expectEqual(matrix.scale(matrix.fromCols(matrix.M2f32, .{ .{ 0.1, 0.2 }, .{ 0.3, 0.4 } }), 0.7).data[1][0], blk: {
        var s = matrix.fromCols(matrix.M2f32, .{ .{ 0.1, 0.2 }, .{ 0.3, 0.4 } });
        scale(0.7, &s);
        break :blk s.data[1][0];
    });
}

test "shape errors" {
    const al = testing.allocator;
    var a: matrix.StaticMat(2, 3, f64) = .{ .data = @splat(@splat(1)) };
    try testing.expectError(error.DimensionMismatch, dot(test_io, &a, &a));
    try testing.expectError(error.DimensionMismatch, dot(test_io, &[_]f64{1}, &[_]f64{ 1, 2 }));
    try testing.expectError(error.DimensionMismatch, mul(al, test_io, &a, &a));
    var out: matrix.M3f64 = undefined;
    try testing.expectError(error.DimensionMismatch, mulInto(al, test_io, &out, &a, &a));
    try testing.expectEqual(@as(f64, 0), try sum(test_io, &[_]f64{}));
}

test {
    zigstd.testRefAllDecls(@This());
}
