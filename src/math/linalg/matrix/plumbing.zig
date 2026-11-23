//! What the matrix files share: type checks, column types, the inline limit's check, and where `@Matrix` and SPIR-V's matrices apply.

const std = @import("std");
const rstd_matrix = @import("rstd_matrix");
const vector = @import("../vector.zig");
const mroot = @import("../matrix.zig");

const Vec = vector.Vec;
const count = vector.count;
const isComplex = vector.isComplex;
const Real = vector.Real;
const mulAdd = vector.mulAdd;
const StaticMat = mroot.StaticMat;
const inline_limit = mroot.inline_limit;

/// `X`'s dimensions: a matrix's rows and columns, a vector's length.
fn dimsOf(comptime X: type) [2]comptime_int {
    return if (isMat(X)) .{ X.rows, X.cols } else .{ count(X), 1 };
}

/// Past the inline limit a `StaticMat` takes `rstd.math.linalg`'s `instead`,
/// with an allocator and `io` (design.md, Dispatch); SPIR-V has no limit.
pub fn expectInline(comptime X: type, comptime Y: type, comptime instead: []const u8) void {
    if (vector.on_spirv) return;
    for (dimsOf(X) ++ dimsOf(Y)) |d| if (d > inline_limit)
        @compileError(std.fmt.comptimePrint("{s} and {s} pass the inline limit ({d}): use {s}", .{ @typeName(X), @typeName(Y), inline_limit, instead }));
}

// ---- plumbing ---------------------------------------------------------------

pub fn isMat(comptime M: type) bool {
    return @typeInfo(M) == .@"struct" and @hasDecl(M, "Element") and @hasDecl(M, "rows") and
        @hasDecl(M, "cols") and M == StaticMat(M.rows, M.cols, M.Element);
}

/// A length-`n` vector of `T`: a `Vec` for real `T`, an array for complex.
pub fn Column(comptime n: comptime_int, comptime T: type) type {
    return if (isComplex(T)) [n]T else Vec(n, T);
}

pub fn expectSquare(comptime M: type) void {
    if (!isMat(M)) @compileError("expects a StaticMat, found " ++ @typeName(M));
    if (M.rows != M.cols) @compileError("expects a square matrix, found " ++ @typeName(M));
}

/// Whether `T` computes through `@Matrix`, which holds real floats only.
pub fn native(comptime T: type) bool {
    return rstd_matrix.has_matrix and !isComplex(T);
}

pub inline fn toNative(m: anytype) rstd_matrix.Matrix(@TypeOf(m).rows, @TypeOf(m).cols, @TypeOf(m).Element) {
    return m.data;
}

/// Whether an `R` by `C` matrix of `T` is a SPIR-V `OpTypeMatrix` (through
/// `@Matrix`, 2 to 4 by 2 to 4, real), which takes SPIR-V's own matrix
/// instructions (design.md, "Dispatch under SPIR-V").
pub fn spirvMat(comptime R: comptime_int, comptime C: comptime_int, comptime T: type) bool {
    return vector.on_spirv and native(T) and R >= 2 and R <= 4 and C >= 2 and C <= 4;
}

// GLSL.std.450 and core matrix instructions by inline asm, on `@Matrix`
// operands; literal templates, as `vector.zig`'s.

pub inline fn spvDeterminant(m: anytype) @TypeOf(m).Element {
    const T = @TypeOf(m).Element;
    return asm (
        \\%set = OpExtInstImport "GLSL.std.450"
        \\%ret = OpExtInst %T %set 33 %m
        : [ret] "" (-> T),
        : [T] "t" (T),
          [m] "" (toNative(m)),
    );
}

pub inline fn spvInverse(m: anytype) rstd_matrix.Matrix(@TypeOf(m).rows, @TypeOf(m).cols, @TypeOf(m).Element) {
    const N = rstd_matrix.Matrix(@TypeOf(m).rows, @TypeOf(m).cols, @TypeOf(m).Element);
    return asm (
        \\%set = OpExtInstImport "GLSL.std.450"
        \\%ret = OpExtInst %T %set 34 %m
        : [ret] "" (-> N),
        : [T] "t" (N),
          [m] "" (toNative(m)),
    );
}

pub inline fn spvOuter(comptime R: comptime_int, comptime C: comptime_int, comptime T: type, a: Vec(R, T), b: Vec(C, T)) rstd_matrix.Matrix(R, C, T) {
    const N = rstd_matrix.Matrix(R, C, T);
    return asm (
        \\%ret = OpOuterProduct %T %a %b
        : [ret] "" (-> N),
        : [T] "t" (N),
          [a] "" (a),
          [b] "" (b),
    );
}

// Element arithmetic over real and complex `T`. A sum of products rounds its
// first product and fuses the rest as `mulAdd` says; a complex product adds
// `re re` then `-im im`, and `re im` then `im re`.

pub inline fn eMul(x: anytype, y: @TypeOf(x)) @TypeOf(x) {
    const T = @TypeOf(x);
    if (comptime !isComplex(T)) return x * y;
    const R = Real(T);
    return .{
        .re = mulAdd(R, -x.im, y.im, x.re * y.re),
        .im = mulAdd(R, x.im, y.re, x.re * y.im),
    };
}

/// `acc + x y`.
pub inline fn eMulAdd(x: anytype, y: @TypeOf(x), acc: @TypeOf(x)) @TypeOf(x) {
    const T = @TypeOf(x);
    if (comptime !isComplex(T)) return mulAdd(T, x, y, acc);
    const R = Real(T);
    return .{
        .re = mulAdd(R, -x.im, y.im, mulAdd(R, x.re, y.re, acc.re)),
        .im = mulAdd(R, x.im, y.re, mulAdd(R, x.re, y.im, acc.im)),
    };
}

pub inline fn eDiv(x: anytype, y: @TypeOf(x)) @TypeOf(x) {
    return if (comptime isComplex(@TypeOf(x))) x.div(y) else x / y;
}
