//! The interfaces of OpenBLAS's extensions: the checks and quick returns of
//! `interface/geadd.c` and `zgeadd.c` before their kernels. Each takes the
//! element type `T` explicitly, any float type or its `Complex` (the vtable's
//! kernels for those it has slots for, the portable ones built for `T`
//! otherwise), and sizes and leading dimensions as `isize`; matrices are
//! column-major, and the row-major (`ORDER` `R`) forms are left out.
//!
//! Bad arguments return an `ArgError` in the order the reference numbers them
//! for `xerbla`, so the first bad one wins; a matrix that would not fit the
//! address space is its leading dimension's error. `io` is for the threading
//! to come (`threading.md`); every routine is serial.

const std = @import("std");
const blas = @import("../__root__.zig");
const types = @import("../types.zig");
const util = @import("../util.zig");
const portable = @import("../portable.zig");
const check = @import("../_test/blas1/check.zig");

const Io = std.Io;
const Trans = types.Trans;
const isComplex = types.isComplex;
const Float = types.Float;

/// The reference's `xerbla` argument, by name.
pub const ArgError = error{
    InvalidTrans,
    InvalidM,
    InvalidN,
    InvalidLda,
    InvalidLdc,
};
/// What an extension routine returns.
pub const Error = ArgError || Io.Cancelable;

fn floatType(comptime T: type, comptime name: []const u8) void {
    if (@typeInfo(Float(T)) != .float) @compileError(name ++ ": " ++ @typeName(T) ++ " is not a float type or its Complex");
}

/// `T`'s elements as the kernels take them: floats, two per complex.
fn raw(comptime T: type, p: [*]const T) [*]Float(T) {
    return @ptrCast(@constCast(p));
}

/// `?GEADD` (`interface/geadd.c`, `zgeadd.c`): `C := alpha*op(A) + beta*C`,
/// `op(A)` and `C` being `m x n`; with `transc` `t`, the same on `C^T`, which
/// `c` stores `n x m`. `op` transposes for `t`, and for real `c`; it never
/// conjugates, so complex `r` and `c` are `error.InvalidTrans`, as are `o` to
/// `d`. The letters come last, as in the reference (which takes any letter
/// but `T` and `C` as `N`). A zero `alpha` never reads `A`, and a zero `beta`
/// never reads `C`.
pub fn geadd(io: Io, comptime T: type, m: isize, n: isize, alpha: T, a: [*]const T, lda: isize, beta: T, c: [*]T, ldc: isize, transa: Trans, transc: Trans) Error!void {
    _ = io;
    comptime floatType(T, "geadd");
    const ta = transa.transposed();
    const tc = transc.transposed();
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, if (ta) n else m)) return error.InvalidLda;
    if (ldc < @max(1, if (tc) n else m)) return error.InvalidLdc;
    if (!plainTrans(T, transa) or !plainTrans(T, transc)) return error.InvalidTrans;
    if (m == 0 or n == 0) return;
    if (!util.fitsMatrix(if (ta) n else m, if (ta) m else n, lda, @sizeOf(T))) return error.InvalidLda;
    if (!util.fitsMatrix(if (tc) n else m, if (tc) m else n, ldc, @sizeOf(T))) return error.InvalidLdc;
    const k = portable.fetch(T, .geadd_k, blas.vtable);
    if (comptime isComplex(T)) {
        _ = k(m, n, alpha.re, alpha.im, raw(T, a), lda, beta.re, beta.im, raw(T, c), ldc, @intFromBool(ta), @intFromBool(tc));
    } else {
        _ = k(m, n, alpha, raw(T, a), lda, beta, raw(T, c), ldc, @intFromBool(ta), @intFromBool(tc));
    }
}

/// Whether `t` is a letter a routine that transposes but never conjugates
/// takes: `n` and `t`, and for a real type `r` and `c`.
fn plainTrans(comptime T: type, t: Trans) bool {
    return switch (t) {
        .n, .t => true,
        .r, .c => !isComplex(T),
        .o, .u, .s, .d => false,
    };
}

// ---- tests -----------------------------------------------------------------

const testing = std.testing;

test "geadd reports the first bad argument and returns for empty problems" {
    const nan = std.math.nan(f64);
    var a: [16]f64 = @splat(nan);
    var c: [16]f64 = @splat(nan);
    try testing.expectError(error.InvalidM, geadd(check.io, f64, -1, -1, 1, &a, 0, 1, &c, 0, .o, .n));
    try testing.expectError(error.InvalidN, geadd(check.io, f64, 0, -1, 1, &a, 0, 1, &c, 0, .n, .n));
    try testing.expectError(error.InvalidLda, geadd(check.io, f64, 3, 2, 1, &a, 2, 1, &c, 0, .n, .n));
    try testing.expectError(error.InvalidLda, geadd(check.io, f64, 3, 2, 1, &a, 1, 1, &c, 3, .t, .n));
    try testing.expectError(error.InvalidLdc, geadd(check.io, f64, 3, 2, 1, &a, 3, 1, &c, 2, .n, .n));
    try testing.expectError(error.InvalidLdc, geadd(check.io, f64, 3, 2, 1, &a, 3, 1, &c, 1, .n, .c));
    try testing.expectError(error.InvalidLdc, geadd(check.io, f64, 0, 2, 1, &a, 1, 1, &c, 0, .n, .n));
    try testing.expectError(error.InvalidTrans, geadd(check.io, f64, 1, 1, 1, &a, 1, 1, &c, 1, .s, .n));
    const Z = types.Complex64;
    var z: [4]Z = @splat(.init(0, 0));
    try testing.expectError(error.InvalidTrans, geadd(check.io, Z, 1, 1, .init(1, 0), &z, 1, .init(1, 0), &z, 1, .n, .c));
    try testing.expectError(error.InvalidTrans, geadd(check.io, Z, 1, 1, .init(1, 0), &z, 1, .init(1, 0), &z, 1, .r, .n));
    // empty: nothing read or written
    try geadd(check.io, f64, 0, 3, nan, &a, 1, nan, &c, 1, .n, .n);
    try geadd(check.io, f64, 3, 0, nan, &a, 3, nan, &c, 3, .t, .t);
    try testing.expect(std.math.isNan(c[0]));
    // a matrix past the address space
    const big = std.math.maxInt(isize) / 2;
    try testing.expectError(error.InvalidLda, geadd(check.io, f64, 2, 3, 1, &a, big, 1, &c, 2, .n, .n));
    try testing.expectError(error.InvalidLdc, geadd(check.io, f64, 2, 3, 1, &a, 2, 1, &c, big, .n, .n));
}

test "geadd matches a plain loop, real letters r and c standing in for n and t" {
    inline for (.{ f32, f64, f16, f128, types.Complex32, types.Complex64, std.math.Complex(f16) }) |T| {
        const complex = comptime isComplex(T);
        const letters: []const [2]Trans = if (complex) &.{ .{ .n, .n }, .{ .t, .n }, .{ .n, .t }, .{ .t, .t } } else &.{ .{ .n, .n }, .{ .c, .r }, .{ .r, .c }, .{ .t, .t } };
        for (letters) |l| {
            const m: usize = 5;
            const n: usize = 3;
            const ta = l[0].transposed();
            const tc = l[1].transposed();
            const lda: usize = (if (ta) n else m) + 1;
            const ldc: usize = (if (tc) n else m) + 2;
            var a: [6 * 6]T = undefined;
            var c: [7 * 7]T = undefined;
            check.fill(T, &a, 3);
            check.fill(T, &c, 4);
            var want = c;
            const alpha: T = if (complex) T.init(0.5, -1) else 0.5;
            const beta: T = if (complex) T.init(-2, 0.25) else -2;
            for (0..n) |j| for (0..m) |i| {
                const k = if (tc) j + i * ldc else i + j * ldc;
                const x = if (ta) a[j + i * lda] else a[i + j * lda];
                want[k] = check.addW(check.mulW(alpha, x), check.mulW(beta, c[k]));
            };
            try geadd(check.io, T, @intCast(m), @intCast(n), alpha, &a, @intCast(lda), beta, &c, @intCast(ldc), l[0], l[1]);
            try check.expectCloseSlices(T, &want, &c, check.tolerance(Float(T), 4));
        }
    }
}
