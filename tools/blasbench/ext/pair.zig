//! The table races' scalars: complex pairs in f64, the Fortran ABI's types, and each entry's sizes.

const std = @import("std");
const bb = @import("../blas_bench.zig");
const kernels = @import("../kernels.zig");
const libs = @import("../libs.zig");
const xr = @import("../ext.zig");

const Prec = bb.Prec;
const gemm_sizes = kernels.gemm_sizes;
const l2_bytes = kernels.l2_bytes;
const Num = libs.Num;
const Cols = xr.Cols;
const Level = xr.Level;
const band_k = xr.band_k;
const lapack_sizes = xr.lapack_sizes;

/// A complex value, or a real one with 0 imaginary part, in f64.
pub const Pair = [2]f64;

pub fn cadd(a: Pair, b: Pair) Pair {
    return .{ a[0] + b[0], a[1] + b[1] };
}
pub fn csub(a: Pair, b: Pair) Pair {
    return .{ a[0] - b[0], a[1] - b[1] };
}
pub fn cmul(a: Pair, b: Pair) Pair {
    return .{ a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0] };
}
pub fn cdiv(a: Pair, b: Pair) Pair {
    const d = b[0] * b[0] + b[1] * b[1];
    return .{ (a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d };
}
pub fn cconj(a: Pair) Pair {
    return .{ a[0], -a[1] };
}
/// `|re| + |im|`, the size `checkGemv` gives a complex.
pub fn cmag(a: Pair) f64 {
    return @abs(a[0]) + @abs(a[1]);
}
pub fn cdist(a: Pair, b: Pair) f64 {
    return @max(@abs(a[0] - b[0]), @abs(a[1] - b[1]));
}
/// `v` rounded to `p`'s type, as the routines get it.
pub fn cnum(comptime p: Prec, v: f64) Pair {
    return .{ @as(p.F, @floatCast(v)), 0 };
}

/// Fortran's `CHARACTER*1` arguments, by reference; each call also passes
/// their hidden lengths (1) last, which a library built without them
/// ignores.
pub const fc = struct {
    pub const U: u8 = 'U';
    pub const L: u8 = 'L';
    pub const N: u8 = 'N';
    pub const S: u8 = 'S';
    pub const V: u8 = 'V';
};
pub const Ch = *const u8;
pub const In = *const c_int;
pub const one_c: c_int = 1;
pub const band_kc: c_int = band_k;

/// A Fortran scalar of `p` by reference: real, or (re, im).
pub fn Sc(comptime p: Prec) type {
    return if (p.complex) [2]p.F else p.F;
}
pub fn sc(comptime p: Prec, v: f64) Sc(p) {
    return if (p.complex) .{ @floatCast(v), 0 } else @floatCast(v);
}
/// The same scalar as rstd takes it.
pub fn num(comptime p: Prec, v: f64) Num(p) {
    return if (p.complex) .init(@floatCast(v), 0) else @floatCast(v);
}

pub fn decl(comptime E: type, comptime name: []const u8, comptime default: anytype) @TypeOf(default) {
    return if (@hasDecl(E, name)) @field(E, name) else default;
}
pub fn hasRstd(comptime E: type) bool {
    return decl(E, "has_rstd", @hasDecl(E, "rstd"));
}
pub fn bandRows(comptime E: type) usize {
    return decl(E, "band_rows", @as(usize, 2 * band_k + 1));
}
pub fn symOf(comptime E: type, comptime p: Prec) []const u8 {
    return if (@hasDecl(E, "sym")) E.sym(p) else p.name ++ E.name;
}

/// An entry's orders at level `lv`.
pub fn orders(comptime E: type, comptime lv: Level, comptime p: Prec) []const usize {
    return comptime switch (lv) {
        .l2 => blk: {
            var ns: [l2_bytes.len]usize = undefined;
            for (l2_bytes, &ns) |bytes, *n| n.* = switch (E.shape) {
                .band => bytes / (p.es() * bandRows(E)),
                .mv, .mm => std.math.sqrt(bytes / p.es()),
            };
            const out = ns;
            break :blk &out;
        },
        .l3 => &gemm_sizes,
        .lapack => &lapack_sizes,
    };
}

/// Elements of `A`, and of `B` and `C`, at order n.
pub fn lenA(comptime E: type, n: usize) usize {
    return switch (E.shape) {
        .mv, .mm => n * n,
        .band => bandRows(E) * n,
    };
}
pub fn lenV(comptime E: type, n: usize) usize {
    return switch (E.shape) {
        .mv, .band => n,
        .mm => n * n,
    };
}

/// Rows and columns of `E`'s output at order n.
pub fn outDims(comptime E: type, n: usize) [2]usize {
    return switch (E.out) {
        .a => .{ if (E.shape == .band) bandRows(E) else n, n },
        .b, .c => if (E.shape == .mm) .{ n, if (decl(E, "cols", Cols.n) == .one) 1 else n } else .{ n, 1 },
        .s => .{ n, 1 },
        else => @compileError("an output is a, b, c or s"),
    };
}

/// Floats of `E`'s output at order n; `s` is real.
pub fn outLen(comptime E: type, comptime p: Prec, n: usize) usize {
    const d = outDims(E, n);
    return d[0] * d[1] * (if (E.out == .s) 1 else p.cs());
}

/// Terms of one output entry, for the bounds.
pub fn terms(comptime E: type, n: usize) usize {
    return if (E.shape == .band) bandRows(E) else n;
}

/// A residual past this many eps is wrong, not rounding.
pub fn wrongLimit(comptime E: type, n: usize) f64 {
    return 64 * @as(f64, @floatFromInt(terms(E, n) + 8));
}
