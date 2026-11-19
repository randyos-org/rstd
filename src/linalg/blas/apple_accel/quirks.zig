//! Quirks: where Accelerate's results differ from the reference's, rstd runs
//! the call (design.md, Apple Accelerate, lists them).

const std = @import("std");
const types = @import("../types.zig");
const portable = @import("../portable.zig");
const util = @import("../util.zig");
const binding = @import("binding.zig");

const isComplex = types.isComplex;
const Float = types.Float;
const isZero = util.isZero;
const x86_64 = binding.x86_64;

/// The level-1 routines `sends` covers; `dot` over `f32` is `sdot`.
pub const L1 = enum { copy, swap, scal, dot, dotc, dsdot, nrm2, asum, rot, rotm, axpy, axpby };

/// Which calls go to Accelerate: all, those whose every stride is 1, none.
const Strides = enum { all, unit, none };

/// `Strides` by letter.
const Letters = struct { s: Strides = .all, d: Strides = .all, c: Strides = .all, z: Strides = .all };

/// Where rstd's kernels are the faster, by speed alone: the arm64 flyoff of
/// 2026-10-29 (design.md, Apple Accelerate, has the rule and its rows);
/// strided there is 2. x86_64 has no data and keeps every call.
const sends: std.enums.EnumFieldStruct(L1, Letters, .{}) = if (x86_64) .{} else .{
    .nrm2 = .{ .s = .none, .d = .none, .c = .none, .z = .none },
    .dot = .{ .s = .none, .c = .unit, .z = .unit },
    .dotc = .{ .c = .unit, .z = .unit },
    .copy = .{ .s = .unit, .d = .unit },
    .swap = .{ .s = .unit, .d = .unit, .c = .unit, .z = .unit },
    .axpy = .{ .c = .unit, .z = .unit },
    .rot = .{ .c = .unit },
};

/// Whether `sends` lets `r` over `T` go to Accelerate; `unit` if every
/// stride is 1.
pub inline fn sent(comptime r: L1, comptime T: type, unit: bool) bool {
    return switch (comptime @field(@field(sends, @tagName(r)), portable.letter(T).?)) {
        .all => true,
        .unit => unit,
        .none => false,
    };
}

/// Zeroes the imaginary parts of the diagonal of the order `n` Hermitian
/// matrix at `a`, full (`ld`) or packed (`ld` 0) on the `lower` triangle.
pub fn realDiagonal(comptime T: type, n: isize, a: [*]T, ld: isize, lower: bool) void {
    if (comptime !isComplex(T)) return;
    const un: usize = @intCast(n);
    for (0..un) |j| {
        const at = if (ld != 0) j * @as(usize, @intCast(ld + 1)) else if (lower) j * (2 * un - j + 1) / 2 else j * (j + 1) / 2 + j;
        a[at].im = 0;
    }
}

/// `nonFinite` of the `lower` or upper triangle of the order `n` matrix at
/// `c`, leading dimension `ld`.
pub fn nonFiniteTri(comptime T: type, n: isize, c: [*]const T, ld: isize, lower: bool) bool {
    const un: usize = @intCast(n);
    for (0..un) |j| {
        const col_ = c + j * @as(usize, @intCast(ld));
        const bad = if (lower) nonFinite(T, n - @as(isize, @intCast(j)), col_ + j, 1) else nonFinite(T, @intCast(j + 1), col_, 1);
        if (bad) return true;
    }
    return false;
}

/// Whether complex `v` has a zero part.
pub fn zeroPart(comptime T: type, v: T) bool {
    if (comptime !isComplex(T)) return false;
    return v.re == 0 or v.im == 0;
}

/// Whether `n` elements of `x` at stride `inc`, from its lowest, hold NaN or
/// an infinity.
pub fn nonFinite(comptime T: type, n: isize, x: [*]const T, inc: isize) bool {
    const F = Float(T);
    const w: usize = if (comptime isComplex(T)) 2 else 1;
    const f: [*]const F = @ptrCast(x);
    const step = w * @abs(inc);
    // `v * 0` is NaN only for NaN and the infinities
    var acc: F = 0;
    for (0..@intCast(n)) |i| {
        for (0..w) |p| acc += f[i * step + p] * 0;
    }
    return acc != acc;
}

/// Whether some `j` has `x(j)`, and `y(j)` if given, zero: a column the
/// reference's rank updates pass over.
pub fn zeroColumn(comptime T: type, n: isize, x: [*]const T, incx: isize, y: ?[*]const T, incy: isize) bool {
    const un: usize = @intCast(n);
    for (0..un) |j| {
        if (!isZero(T, x[elem(j, un, incx)])) continue;
        const yv = y orelse return true;
        if (isZero(T, yv[elem(j, un, incy)])) return true;
    }
    return false;
}

/// Element `j`'s offset from the lowest of `n` at stride `inc`.
fn elem(j: usize, n: usize, inc: isize) usize {
    return (if (inc > 0) j else n - 1 - j) * @abs(inc);
}

/// Whether an entry on the diagonal of the order `n` matrix at `a` is too
/// small for its reciprocal: below the smallest normal (for complex, its
/// square root, as `1 / d` divides by `|d|^2`).
pub fn tinyDiagonal(comptime T: type, n: isize, a: [*]const T, lda: isize) bool {
    const F = Float(T);
    const lim = if (comptime isComplex(T)) @sqrt(std.math.floatMin(F)) else std.math.floatMin(F);
    const un: usize = @intCast(n);
    const ld: usize = @intCast(lda);
    for (0..un) |i| {
        const v = a[i * (ld + 1)];
        const m = if (comptime isComplex(T)) @max(@abs(v.re), @abs(v.im)) else @abs(v);
        if (m < lim) return true;
    }
    return false;
}

/// `nonFinite` of a `rows` by `cols` matrix at leading dimension `ld`.
pub fn nonFiniteMat(comptime T: type, rows: isize, cols: isize, b: [*]const T, ld: isize) bool {
    for (0..@intCast(cols)) |j| {
        if (nonFinite(T, rows, b + j * @as(usize, @intCast(ld)), 1)) return true;
    }
    return false;
}
