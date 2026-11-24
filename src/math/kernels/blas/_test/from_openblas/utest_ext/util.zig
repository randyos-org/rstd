//! What utest_ext's files share: the sources' tolerances, their `Case`, values and norms.

const std = @import("std");
const ob = @import("../__root__.zig");
const blas = @import("../../../__root__.zig");

pub const Float = ob.Float;
pub const C32 = std.math.Complex(f32);
pub const C64 = std.math.Complex(f64);

pub const single_eps = 1e-4;
pub const double_eps = 1e-13;
pub const single_tol = 1e-2;
pub const double_tol = 1e-10;
pub const single_eps_zgemv = 5e-3;
pub const double_eps_zgemv = 1e-11;

pub fn tol(comptime T: type, single: f64, double: f64) f64 {
    return if (Float(T) == f32) single else double;
}

// ---- the parameter cases ---------------------------------------------------------

/// One case: the source's variables by their own names.
pub const Case = struct {
    name: []const u8,
    kind: enum { run, nan, bad } = .run,
    n: isize = 0,
    N: isize = 0,
    m: isize = 0,
    k: isize = 0,
    kl: isize = 0,
    ku: isize = 0,
    lda: isize = 0,
    ldb: isize = 0,
    ldc: isize = 0,
    ldm: isize = 0,
    incx: isize = 1,
    incy: isize = 1,
    inc: isize = 1,
    inc_x: isize = 1,
    inc_y: isize = 1,
    inc_b: isize = 1,
    inc_c: isize = 1,
    trans: u8 = 'N',
    transa: u8 = 'N',
    transb: u8 = 'N',
    uplo: u8 = 'U',
    diag: u8 = 'N',
    order: u8 = 'C',
    alpha: [2]f64 = .{ 0, 0 },
    beta: [2]f64 = .{ 0, 0 },
    c: [2]f64 = .{ 0, 0 },
    s: [2]f64 = .{ 0, 0 },
    info: u8 = 0,
};

pub fn val(comptime T: type, v: [2]f64) T {
    return ob.val(T, v[0], v[1]);
}

/// `drand_generate`: values in [0, 1).
pub fn fill(comptime T: type, p: []T, r: std.Random) void {
    for (p) |*v| v.* = if (comptime isComplex(T)) .init(@floatCast(r.float(f64)), @floatCast(r.float(f64))) else @floatCast(r.float(f64));
}

pub fn span(n: isize, inc: isize) usize {
    return @intCast(@max(n, 1) * @max(@as(isize, @intCast(@abs(inc))), 1));
}

/// The norm of `got - want` over `n` elements at stride `|inc|`; elements
/// that agree exactly count as zero.
pub fn diffNorm(comptime T: type, n: isize, got: []const T, want: []const T, inc: isize) f64 {
    if (n <= 0) return 0;
    var s: f64 = 0;
    const st: usize = @abs(inc);
    for (0..@intCast(n)) |i| {
        const g = got[i * st];
        const w = want[i * st];
        if (ob.eql(g, w)) continue;
        const d = ob.sub(g, w);
        s += if (comptime isComplex(T)) @as(f64, d.re) * d.re + @as(f64, d.im) * d.im else @as(f64, d) * d;
    }
    return @sqrt(s);
}

/// The norm of `v`, `n` elements at stride `|inc|`.
pub fn norm2(comptime T: type, n: isize, v: []const T, inc: isize) f64 {
    var s: f64 = 0;
    for (0..@intCast(@max(n, 0))) |i| s += @as(f64, ob.absv(v[i * @abs(inc)])) * ob.absv(v[i * @abs(inc)]);
    return @sqrt(s);
}

pub fn expectNorm(name: []const u8, comptime T: type, norm: f64, t: f64) !void {
    if (norm <= t) return;
    ob.say("{s} {s}: norm of differences {e} over {e}\n", .{ @typeName(T), name, norm, t });
    return error.TestUnexpectedResult;
}

/// `T`'s elements of `p`, conjugated: the source's `?conjugate_vector`.
pub fn conjugate(comptime T: type, p: []T) void {
    if (comptime isComplex(T)) for (p) |*v| {
        v.* = v.conjugate();
    };
}

pub fn transOf(c: u8) blas.Trans {
    return ob.letter(blas.Trans, c) orelse .d;
}

pub fn transposed(c: u8) bool {
    return c == 'T' or c == 'C' or c == 'U' or c == 'D';
}

// ---- the matrix extensions -----------------------------------------------------

/// The sources' `?matrix_difference`: the mean over `cols` columns of `ld`
/// of the norm of `got - want` down `rows` of each; elements that agree
/// exactly count as zero.
pub fn matDiff(comptime T: type, rows: usize, cols: usize, got: []const T, want: []const T, ld: usize) f64 {
    if (rows == 0 or cols == 0) return 0;
    var sum: f64 = 0;
    for (0..cols) |j| sum += diffNorm(T, @intCast(rows), got[j * ld ..], want[j * ld ..], 1);
    return sum / @as(f64, @floatFromInt(cols));
}

pub fn transposedLetter(c: u8) bool {
    return c == 'T' or c == 'C';
}

pub const isComplex = ob.isComplex;
