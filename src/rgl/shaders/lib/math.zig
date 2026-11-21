const std = @import("std");
const linalg = @import("../../../math/linalg/__root__.zig");
const vec = linalg.vector;
const V3f32 = vec.V3f32;

pub const radians = std.math.degreesToRadians;
pub const degrees = std.math.radiansToDegrees;
pub const sin = std.math.sin;
pub const cos = std.math.cos;
pub const tan = std.math.tan;
pub const asin = std.math.asin;
pub const acos = std.math.acos;
pub const atan = std.math.atan;
pub const pow = std.math.pow;
pub const exp = std.math.exp;
pub const log = std.math.log;
pub const exp2 = std.math.exp2;
pub const log2 = std.math.log2;
pub const sqrt = std.math.sqrt;
pub const sign = std.math.sign;
pub const floor = std.math.floor;
pub const ceil = std.math.ceil;
pub const mod = std.math.mod;

/// Scalar `s` as `T`: splatted when `T` is a vector.
inline fn like(comptime T: type, s: anytype) T {
    return if (@typeInfo(T) == .vector) @splat(s) else s;
}

/// Whether `T` is a float or a vector of floats: what GLSL.std.450 takes.
fn isFloat(comptime T: type) bool {
    return switch (@typeInfo(T)) {
        .float => true,
        .vector => |v| @typeInfo(v.child) == .float,
        else => false,
    };
}

// Under SPIR-V the operations below are GLSL.std.450 instructions by inline
// `asm` (lib/readme.md), each template a literal: a computed one crashes the
// backend (issue #199). Elsewhere they are the plain formulas.
const spirv = vec.on_spirv;

pub inline fn inversesqrt(x: anytype) @TypeOf(x) {
    if (comptime spirv) return asm (
        \\%set = OpExtInstImport "GLSL.std.450"
        \\%ret = OpExtInst %T %set 32 %x
        : [ret] "" (-> @TypeOf(x)),
        : [T] "t" (@TypeOf(x)),
          [x] "" (x),
    );
    return like(@TypeOf(x), 1) / @sqrt(x);
}

pub fn stepScalar(edge: f32, x: f32) f32 {
    if (comptime spirv) return asm (
        \\%set = OpExtInstImport "GLSL.std.450"
        \\%ret = OpExtInst %T %set 48 %edge %x
        : [ret] "" (-> f32),
        : [T] "t" (f32),
          [edge] "" (edge),
          [x] "" (x),
    );
    return if (x < edge) 0 else 1;
}

pub inline fn fract(x: anytype) @TypeOf(x) {
    if (comptime spirv) return asm (
        \\%set = OpExtInstImport "GLSL.std.450"
        \\%ret = OpExtInst %T %set 10 %x
        : [ret] "" (-> @TypeOf(x)),
        : [T] "t" (@TypeOf(x)),
          [x] "" (x),
    );
    return x - @floor(x);
}

/// `std.math.clamp`; under SPIR-V, FClamp for floats.
pub inline fn clamp(x: anytype, lo: anytype, hi: anytype) @TypeOf(x, lo, hi) {
    const T = @TypeOf(x, lo, hi);
    if (comptime spirv and isFloat(T)) return asm (
        \\%set = OpExtInstImport "GLSL.std.450"
        \\%ret = OpExtInst %T %set 43 %x %lo %hi
        : [ret] "" (-> T),
        : [T] "t" (T),
          [x] "" (@as(T, x)),
          [lo] "" (@as(T, lo)),
          [hi] "" (@as(T, hi)),
    );
    return std.math.clamp(x, lo, hi);
}

/// Median of three scalars; the standard way to recover a signed distance
/// from a 3-channel MSDF texture sample
pub inline fn median3(a: anytype, b: @TypeOf(a), c: @TypeOf(a)) @TypeOf(a) {
    return @max(@min(a, b), @min(@max(a, b), c));
}

pub const mix = vec.mix;
pub const dot = vec.dot;
pub const dot2 = vec.dot2;
pub const length = vec.length;
pub const distance = vec.distance;
pub const cross = vec.cross;
pub const sum = vec.sum;

pub inline fn normalize(v: V3f32) V3f32 {
    return vec.norm(v) catch @splat(0);
}

pub fn smoothstepScalar(edge0: f32, edge1: f32, x: f32) f32 {
    if (comptime spirv) return asm (
        \\%set = OpExtInstImport "GLSL.std.450"
        \\%ret = OpExtInst %T %set 49 %e0 %e1 %x
        : [ret] "" (-> f32),
        : [T] "t" (f32),
          [e0] "" (edge0),
          [e1] "" (edge1),
          [x] "" (x),
    );
    const t = std.math.clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
    return t * t * (3.0 - 2.0 * t);
}

pub inline fn reflect(i: anytype, n: anytype) @TypeOf(i) {
    const T = @TypeOf(i);
    if (comptime spirv) return asm (
        \\%set = OpExtInstImport "GLSL.std.450"
        \\%ret = OpExtInst %T %set 71 %i %n
        : [ret] "" (-> T),
        : [T] "t" (T),
          [i] "" (i),
          [n] "" (@as(T, n)),
    );
    return i - like(T, 2 * dot(i, n)) * n;
}

pub fn refract(i: anytype, n: anytype, eta: f32) @TypeOf(i) {
    const T = @TypeOf(i);
    if (comptime spirv) return asm (
        \\%set = OpExtInstImport "GLSL.std.450"
        \\%ret = OpExtInst %T %set 72 %i %n %eta
        : [ret] "" (-> T),
        : [T] "t" (T),
          [i] "" (i),
          [n] "" (@as(T, n)),
          [eta] "" (eta),
    );
    const dotni = dot(n, i);
    const k = 1.0 - eta * eta * (1.0 - dotni * dotni);
    return if (k < 0.0) like(T, 0) else like(T, eta) * i - like(T, eta * dotni + @sqrt(k)) * n;
}

pub inline fn faceforward(n: anytype, i: anytype, nref: anytype) @TypeOf(n) {
    const T = @TypeOf(n);
    if (comptime spirv) return asm (
        \\%set = OpExtInstImport "GLSL.std.450"
        \\%ret = OpExtInst %T %set 70 %n %i %nref
        : [ret] "" (-> T),
        : [T] "t" (T),
          [n] "" (n),
          [i] "" (@as(T, i)),
          [nref] "" (@as(T, nref)),
    );
    return if (dot(nref, i) < 0.0) n else -n;
}

// ---- GLSL 4.60's chapter 8 math, by GLSL's names ----

// angle and trigonometry; `atan(y, x)` is `atan2`, as Zig has no overloads
pub const atan2 = std.math.atan2;
pub const sinh = std.math.sinh;
pub const cosh = std.math.cosh;
pub const tanh = std.math.tanh;
pub const asinh = std.math.asinh;
pub const acosh = std.math.acosh;
pub const atanh = std.math.atanh;

// common
pub inline fn abs(x: anytype) @TypeOf(x) {
    return @abs(x);
}
pub inline fn min(x: anytype, y: anytype) @TypeOf(x, y) {
    return @min(x, y);
}
pub inline fn max(x: anytype, y: anytype) @TypeOf(x, y) {
    return @max(x, y);
}
pub inline fn trunc(x: anytype) @TypeOf(x) {
    return @trunc(x);
}
/// Halves away from zero; GLSL leaves the direction to the implementation.
pub inline fn round(x: anytype) @TypeOf(x) {
    return @round(x);
}
/// Halves to the even neighbor: there, `2 round(x / 2)`.
pub inline fn roundEven(x: anytype) @TypeOf(x) {
    const T = @TypeOf(x);
    const half = @abs(x - @trunc(x)) == like(T, 0.5);
    const even = like(T, 2) * @round(x / like(T, 2));
    if (@typeInfo(T) == .vector) return @select(Scalar(T), half, even, @round(x));
    return if (half) even else @round(x);
}
/// GLSL's `modf`: the fraction, `whole` the integer part, both with `x`'s sign.
pub inline fn modf(x: anytype, whole: *@TypeOf(x)) @TypeOf(x) {
    whole.* = @trunc(x);
    return x - whole.*;
}
pub inline fn isnan(x: anytype) Mask(@TypeOf(x)) {
    return x != x;
}
pub inline fn isinf(x: anytype) Mask(@TypeOf(x)) {
    return @abs(x) == like(@TypeOf(x), std.math.inf(Scalar(@TypeOf(x))));
}
pub inline fn fma(a: anytype, b: anytype, c: anytype) @TypeOf(a, b, c) {
    return @mulAdd(@TypeOf(a, b, c), a, b, c);
}
pub inline fn floatBitsToInt(x: anytype) Bits(@TypeOf(x), .signed) {
    return @bitCast(x);
}
pub inline fn floatBitsToUint(x: anytype) Bits(@TypeOf(x), .unsigned) {
    return @bitCast(x);
}
pub inline fn intBitsToFloat(x: anytype) Floats(@TypeOf(x)) {
    return @bitCast(x);
}
pub inline fn uintBitsToFloat(x: anytype) Floats(@TypeOf(x)) {
    return @bitCast(x);
}
/// GLSL's `frexp`: the significand in [0.5, 1) (0 for 0), `exponent` the
/// power of 2.
pub fn frexp(x: f32, exponent: *i32) f32 {
    const r = std.math.frexp(x);
    exponent.* = r.exponent;
    return r.significand;
}
pub fn ldexp(x: f32, exponent: i32) f32 {
    return std.math.ldexp(x, exponent);
}
/// Element by element: 0 below `edge`, else 1; scalars or vectors.
pub inline fn step(edge: anytype, x: anytype) @TypeOf(edge, x) {
    const T = @TypeOf(edge, x);
    const below = @as(T, x) < @as(T, edge);
    if (@typeInfo(T) == .vector) return @select(Scalar(T), below, like(T, 0), like(T, 1));
    return if (below) 0 else 1;
}
/// Element by element, scalars or vectors.
pub inline fn smoothstep(edge0: anytype, edge1: anytype, x: anytype) @TypeOf(edge0, edge1, x) {
    const T = @TypeOf(edge0, edge1, x);
    const t = std.math.clamp((@as(T, x) - @as(T, edge0)) / (@as(T, edge1) - @as(T, edge0)), like(T, 0), like(T, 1));
    return t * t * (like(T, 3) - like(T, 2) * t);
}

// matrix
pub const outerProduct = linalg.matrix.outer;
pub const transpose = linalg.matrix.transpose;
pub const determinant = linalg.matrix.det;
pub const inverse = linalg.matrix.inv;
/// Element by element.
pub inline fn matrixCompMult(a: anytype, b: @TypeOf(a)) @TypeOf(a) {
    var r = a;
    for (&r.data, b.data) |*x, y| x.* *= y;
    return r;
}

// vector relational
pub inline fn lessThan(x: anytype, y: @TypeOf(x)) Mask(@TypeOf(x)) {
    return x < y;
}
pub inline fn lessThanEqual(x: anytype, y: @TypeOf(x)) Mask(@TypeOf(x)) {
    return x <= y;
}
pub inline fn greaterThan(x: anytype, y: @TypeOf(x)) Mask(@TypeOf(x)) {
    return x > y;
}
pub inline fn greaterThanEqual(x: anytype, y: @TypeOf(x)) Mask(@TypeOf(x)) {
    return x >= y;
}
pub inline fn equal(x: anytype, y: @TypeOf(x)) Mask(@TypeOf(x)) {
    return x == y;
}
pub inline fn notEqual(x: anytype, y: @TypeOf(x)) Mask(@TypeOf(x)) {
    return x != y;
}
pub inline fn any(x: anytype) bool {
    return @reduce(.Or, x);
}
pub inline fn all(x: anytype) bool {
    return @reduce(.And, x);
}
pub inline fn not(x: anytype) @TypeOf(x) {
    return !x;
}

// packing, GLSL 4.60 section 8.4
pub fn packUnorm2x16(v: @Vector(2, f32)) u32 {
    return packN(2, 16, false, v);
}
pub fn packSnorm2x16(v: @Vector(2, f32)) u32 {
    return packN(2, 16, true, v);
}
pub fn packUnorm4x8(v: @Vector(4, f32)) u32 {
    return packN(4, 8, false, v);
}
pub fn packSnorm4x8(v: @Vector(4, f32)) u32 {
    return packN(4, 8, true, v);
}
pub fn unpackUnorm2x16(p: u32) @Vector(2, f32) {
    return unpackN(2, 16, false, p);
}
pub fn unpackSnorm2x16(p: u32) @Vector(2, f32) {
    return unpackN(2, 16, true, p);
}
pub fn unpackUnorm4x8(p: u32) @Vector(4, f32) {
    return unpackN(4, 8, false, p);
}
pub fn unpackSnorm4x8(p: u32) @Vector(4, f32) {
    return unpackN(4, 8, true, p);
}
pub fn packHalf2x16(v: @Vector(2, f32)) u32 {
    const lo: u16 = @bitCast(@as(f16, @floatCast(v[0])));
    const hi: u16 = @bitCast(@as(f16, @floatCast(v[1])));
    return @as(u32, lo) | (@as(u32, hi) << 16);
}
pub fn unpackHalf2x16(p: u32) @Vector(2, f32) {
    const lo: f16 = @bitCast(@as(u16, @truncate(p)));
    const hi: f16 = @bitCast(@as(u16, @truncate(p >> 16)));
    return .{ lo, hi };
}

/// `n` fields of `bits` each, the first lowest: `round(clamp(c) * max)`.
fn packN(comptime n: usize, comptime bits: u5, comptime signed: bool, v: @Vector(n, f32)) u32 {
    const top: f32 = @floatFromInt((@as(u32, 1) << (bits - @intFromBool(signed))) - 1);
    const mask: u32 = (@as(u32, 1) << bits) - 1;
    var r: u32 = 0;
    inline for (0..n) |k| {
        const c = std.math.clamp(v[k], if (signed) -1 else 0, 1);
        const q: i32 = @intFromFloat(@round(c * top));
        r |= (@as(u32, @bitCast(q)) & mask) << (bits * k);
    }
    return r;
}

fn unpackN(comptime n: usize, comptime bits: u5, comptime signed: bool, p: u32) @Vector(n, f32) {
    const top: f32 = @floatFromInt((@as(u32, 1) << (bits - @intFromBool(signed))) - 1);
    const mask: u32 = (@as(u32, 1) << bits) - 1;
    var r: @Vector(n, f32) = undefined;
    inline for (0..n) |k| {
        const field = (p >> (bits * k)) & mask;
        if (signed) {
            // sign-extend the field
            const s: i32 = @as(i32, @bitCast(field << (32 - @as(u6, bits)))) >> (32 - @as(u6, bits));
            r[k] = std.math.clamp(@as(f32, @floatFromInt(s)) / top, -1, 1);
        } else r[k] = @as(f32, @floatFromInt(field)) / top;
    }
    return r;
}

fn Scalar(comptime T: type) type {
    return if (@typeInfo(T) == .vector) @typeInfo(T).vector.child else T;
}

/// `bool`, or a vector of them as long as `T`.
fn Mask(comptime T: type) type {
    return if (@typeInfo(T) == .vector) @Vector(@typeInfo(T).vector.len, bool) else bool;
}

/// GLSL's bit casts are 32-bit: `i32` or `u32` lanes for `f32` ones.
fn Bits(comptime T: type, comptime signedness: std.lang.Signedness) type {
    if (Scalar(T) != f32) @compileError("GLSL's float bit casts take f32");
    const I = if (signedness == .signed) i32 else u32;
    return if (@typeInfo(T) == .vector) @Vector(@typeInfo(T).vector.len, I) else I;
}

/// `f32` lanes for `i32` or `u32` ones.
fn Floats(comptime T: type) type {
    if (@bitSizeOf(Scalar(T)) != 32) @compileError("GLSL's int bit casts take 32-bit ints");
    const F = f32;
    return if (@typeInfo(T) == .vector) @Vector(@typeInfo(T).vector.len, F) else F;
}

test "GLSL names: relational, rounding, step, packing" {
    const t = std.testing;
    const V = @Vector(3, f32);
    try t.expect(all(lessThan(V{ 1, 2, 3 }, V{ 2, 3, 4 })));
    try t.expect(!any(equal(V{ 1, 2, 3 }, V{ 0, 0, 0 })));
    try t.expectEqual(V{ 2, 4, 3 }, roundEven(V{ 2.5, 3.5, 2.7 }));
    try t.expectEqual(V{ 0, 1, 1 }, step(V{ 1, 1, 1 }, V{ 0, 1, 2 }));
    var whole: f32 = undefined;
    try t.expectEqual(@as(f32, -0.5), modf(@as(f32, -2.5), &whole));
    try t.expectEqual(@as(f32, -2), whole);
    try t.expectEqual(@as(u32, 0xffff_0000), packUnorm2x16(.{ 0, 1 }));
    try t.expectEqual(@Vector(2, f32){ 0, 1 }, unpackUnorm2x16(0xffff_0000));
    try t.expectEqual(@Vector(4, f32){ -1, 0, 1, 1 }, unpackSnorm4x8(packSnorm4x8(.{ -1, 0, 1, 2 })));
    try t.expectEqual(@Vector(2, f32){ 0.5, -2 }, unpackHalf2x16(packHalf2x16(.{ 0.5, -2 })));
    try t.expectEqual(@as(i32, 0x3f800000), floatBitsToInt(@as(f32, 1)));
}
