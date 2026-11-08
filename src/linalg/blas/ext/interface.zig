//! The interfaces of OpenBLAS's extensions: the checks and quick returns of
//! `interface/geadd.c`, `omatcopy.c`, `imatcopy.c`, `lapack/laswp.c` and
//! their `z` forms before their kernels. Each takes the
//! element type `T` explicitly, any float type or its `Complex` (the vtable's
//! kernels for those it has slots for, the portable ones built for `T`
//! otherwise), and sizes and leading dimensions as `isize`; matrices are
//! column-major, and the row-major (`ORDER` `R`) forms are left out.
//!
//! Bad arguments return an `ArgError` in the order the reference numbers them
//! for `xerbla`, so the first bad one wins; a matrix that would not fit the
//! address space is its leading dimension's error. `io` is unused; every
//! routine is serial.

const std = @import("std");
const blas = @import("../__root__.zig");
const types = @import("../types.zig");
const util = @import("../util.zig");
const portable = @import("../portable.zig");
const matcopy = @import("matcopy.zig");
const laswp_kernel = @import("laswp/kernel.zig");
const ref = @import("laswp/__root__.zig");
const check = @import("../_test/check.zig");
const bitwise = @import("../_test/bitwise.zig");
const reference = @import("../_test/reference.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Trans = types.Trans;
const isComplex = types.isComplex;
const Float = types.Float;

/// The reference's `xerbla` argument, by name.
pub const ArgError = error{
    InvalidTrans,
    InvalidM,
    InvalidN,
    InvalidRows,
    InvalidCols,
    InvalidLda,
    InvalidLdb,
    InvalidLdc,
    InvalidK1,
    InvalidIncx,
    InvalidIpiv,
};
/// What an extension routine returns; only `imatcopy` allocates.
pub const Error = ArgError || Allocator.Error || Io.Cancelable;

const raw = util.raw;

/// `?GEADD` (`interface/geadd.c`, `zgeadd.c`): `C := alpha*op(A) + beta*C`,
/// `op(A)` and `C` being `m x n`; with `transc` `t`, the same on `C^T`, which
/// `c` stores `n x m`. `op` transposes for `t`, and for real `c`; it never
/// conjugates, so complex `r` and `c` are `error.InvalidTrans`, as are `o` to
/// `d`. The letters come last, as in the reference (which takes any letter
/// but `T` and `C` as `N`). A zero `alpha` never reads `A`, and a zero `beta`
/// never reads `C`.
pub fn geadd(io: Io, comptime T: type, m: isize, n: isize, alpha: T, a: [*]const T, lda: isize, beta: T, c: [*]T, ldc: isize, transa: Trans, transc: Trans) Error!void {
    _ = io;
    comptime util.requireFloat(T, "geadd");
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

/// `?OMATCOPY` (`interface/omatcopy.c`, `zomatcopy.c`), column-major:
/// `B := alpha*op(A)`, `A` being `rows x cols`, so `B` is `cols x rows`
/// where `op` transposes. `trans` is `n`, `t`, or `r` and `c` conjugated (as
/// `n` and `t` for a real type); `o` to `d` are `error.InvalidTrans`. A zero
/// `alpha` stores zeros without reading `A`, and `alpha` 1 copies it exactly.
/// `A` and `B` must not overlap; `imatcopy` is the in-place form.
pub fn omatcopy(io: Io, comptime T: type, trans: Trans, rows: isize, cols: isize, alpha: T, a: [*]const T, lda: isize, b: [*]T, ldb: isize) Error!void {
    _ = io;
    comptime util.requireFloat(T, "omatcopy");
    const t = trans.transposed();
    try matcopyArgs(trans, rows, cols, lda, ldb);
    if (rows == 0 or cols == 0) return;
    if (!util.fitsMatrix(rows, cols, lda, @sizeOf(T))) return error.InvalidLda;
    if (!util.fitsMatrix(if (t) cols else rows, if (t) rows else cols, ldb, @sizeOf(T))) return error.InvalidLdb;
    outOf(T, trans, rows, cols, alpha, a, lda, b, ldb);
}

/// `?IMATCOPY` (`interface/imatcopy.c`, `zimatcopy.c`), column-major: `A :=
/// alpha*op(A)` in place, `A` read `rows x cols` with leading dimension
/// `lda` and written as `op(A)` with `ldb`. Where the shape and the leading
/// dimension stay (no transpose, or a square one, and `lda == ldb`) the
/// in-place kernels run, as the reference's `NEW_IMATCOPY` build does; else
/// `op(A)` goes through a buffer of `rows*cols` from `allocator`, as its
/// default build does. Both give the same bits. Arguments as `omatcopy`'s.
pub fn imatcopy(allocator: Allocator, io: Io, comptime T: type, trans: Trans, rows: isize, cols: isize, alpha: T, a: [*]T, lda: isize, ldb: isize) Error!void {
    _ = io;
    comptime util.requireFloat(T, "imatcopy");
    const t = trans.transposed();
    try matcopyArgs(trans, rows, cols, lda, ldb);
    if (rows == 0 or cols == 0) return;
    const out_rows = if (t) cols else rows;
    const out_cols = if (t) rows else cols;
    if (!util.fitsMatrix(rows, cols, lda, @sizeOf(T))) return error.InvalidLda;
    if (!util.fitsMatrix(out_rows, out_cols, ldb, @sizeOf(T))) return error.InvalidLdb;
    if (lda == ldb and (!t or rows == cols)) {
        switch (pick(T, trans)) {
            inline else => |p| if (comptime realConj(T, p)) unreachable else {
                const k = portable.fetch(T, p.inOp(), blas.vtable);
                if (comptime isComplex(T)) {
                    _ = k(rows, cols, alpha.re, alpha.im, raw(T, a), lda);
                } else {
                    _ = k(rows, cols, alpha, raw(T, a), lda);
                }
            },
        }
        return;
    }
    const buf = try allocator.alloc(T, @intCast(rows * cols));
    defer allocator.free(buf);
    outOf(T, trans, rows, cols, alpha, a, lda, buf.ptr, out_rows);
    const one: T = if (comptime isComplex(T)) .init(1, 0) else 1;
    outOf(T, .n, out_rows, out_cols, one, buf.ptr, out_rows, a, ldb);
}

/// `?LASWP` (`interface/lapack/laswp.c`, `zlaswp.c`,
/// `lapack/laswp/generic/laswp_k_1.c`), as `dlaswp.f`: rows `k1` to `k2`
/// (1-based; from `k2` for `incx < 0`) swapped with rows `ipiv(k)`; every
/// pivot checked before the first swap.
pub fn laswp(io: Io, comptime T: type, n: isize, a: [*]T, lda: isize, k1: isize, k2: isize, ipiv: [*]const isize, incx: isize) Error!void {
    _ = io;
    comptime util.requireFloat(T, "laswp");
    if (n <= 0 or incx == 0 or k2 < k1) return;
    if (lda < @max(1, k2)) return error.InvalidLda;
    if (k1 < 1) return error.InvalidK1;
    if (!util.fitsMatrix(lda, n, lda, @sizeOf(T))) return error.InvalidLda;
    if (!util.fits(k2 - k1 + 1, incx, @sizeOf(isize))) return error.InvalidIncx;
    const step: usize = @abs(incx);
    const first: usize = @intCast(k1 - 1);
    for (0..@intCast(k2 - k1 + 1)) |r| {
        const p = ipiv[first + r * step];
        if (p < 1 or p > lda) return error.InvalidIpiv;
    }
    laswp_kernel.swapRows(T, @intCast(n), a, @intCast(lda), @intCast(k1), @intCast(k2), ipiv + first, incx);
}

/// The checks of `omatcopy` and `imatcopy`, in `xerbla`'s order (`ldb` is
/// the 9th argument of one and the 8th of the other, last either way).
fn matcopyArgs(trans: Trans, rows: isize, cols: isize, lda: isize, ldb: isize) ArgError!void {
    if (trans.conjX()) return error.InvalidTrans;
    if (rows < 0) return error.InvalidRows;
    if (cols < 0) return error.InvalidCols;
    if (lda < @max(1, rows)) return error.InvalidLda;
    if (ldb < @max(1, if (trans.transposed()) cols else rows)) return error.InvalidLdb;
}

/// Which `omatcopy` or `imatcopy` kernel `trans` takes.
const Pick = matcopy.Copy;

fn pick(comptime T: type, trans: Trans) Pick {
    const c = isComplex(T) and trans.conjA();
    return if (trans.transposed()) (if (c) .ctc else .ct) else (if (c) .cnc else .cn);
}

/// Whether `p` conjugates a real `T`, which has no such kernel.
fn realConj(comptime T: type, comptime p: Pick) bool {
    return !isComplex(T) and (p == .cnc or p == .ctc);
}

/// `B := alpha*op(A)` by the `omatcopy` kernel `trans` takes.
fn outOf(comptime T: type, trans: Trans, rows: isize, cols: isize, alpha: T, a: [*]const T, lda: isize, b: [*]T, ldb: isize) void {
    switch (pick(T, trans)) {
        inline else => |p| if (comptime realConj(T, p)) unreachable else {
            const k = portable.fetch(T, p.outOp(), blas.vtable);
            if (comptime isComplex(T)) {
                _ = k(rows, cols, alpha.re, alpha.im, raw(T, a), lda, raw(T, b), ldb);
            } else {
                _ = k(rows, cols, alpha, raw(T, a), lda, raw(T, b), ldb);
            }
        },
    }
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
                want[k] = reference.addW(reference.mulW(alpha, x), reference.mulW(beta, c[k]));
            };
            try geadd(check.io, T, @intCast(m), @intCast(n), alpha, &a, @intCast(lda), beta, &c, @intCast(ldc), l[0], l[1]);
            try check.expectCloseSlices(T, &want, &c, check.tolerance(Float(T), 4));
        }
    }
}

test "omatcopy and imatcopy report the first bad argument and return for empty problems" {
    const fa = testing.failing_allocator;
    const nan = std.math.nan(f64);
    var a: [16]f64 = @splat(nan);
    var b: [16]f64 = @splat(nan);
    try testing.expectError(error.InvalidTrans, omatcopy(check.io, f64, .u, -1, -1, 1, &a, 0, &b, 0));
    try testing.expectError(error.InvalidRows, omatcopy(check.io, f64, .n, -1, -1, 1, &a, 0, &b, 0));
    try testing.expectError(error.InvalidCols, omatcopy(check.io, f64, .t, 0, -1, 1, &a, 0, &b, 0));
    try testing.expectError(error.InvalidLda, omatcopy(check.io, f64, .t, 3, 2, 1, &a, 2, &b, 0));
    try testing.expectError(error.InvalidLdb, omatcopy(check.io, f64, .n, 3, 2, 1, &a, 3, &b, 2));
    try testing.expectError(error.InvalidLdb, omatcopy(check.io, f64, .c, 3, 2, 1, &a, 3, &b, 1));
    try testing.expectError(error.InvalidLdb, imatcopy(fa, check.io, f64, .t, 3, 2, 1, &a, 3, 1));
    try testing.expectError(error.InvalidLda, imatcopy(fa, check.io, f64, .n, 0, 2, 1, &a, 0, 1));
    try omatcopy(check.io, f64, .n, 0, 3, nan, &a, 1, &b, 1);
    try imatcopy(fa, check.io, f64, .t, 3, 0, nan, &a, 3, 1);
    try testing.expect(std.math.isNan(b[0]) and std.math.isNan(a[0]));
    // in place without a buffer where the shape and leading dimension stay
    a[0] = 2;
    a[1] = 3;
    a[2] = 4;
    a[3] = 5;
    try imatcopy(fa, check.io, f64, .t, 2, 2, 2, &a, 2, 2);
    try testing.expectEqualSlices(f64, &.{ 4, 8, 6, 10 }, a[0..4]);
    try testing.expectError(error.OutOfMemory, imatcopy(fa, check.io, f64, .t, 2, 1, 2, &a, 2, 1));
    const big = std.math.maxInt(isize) / 2;
    try testing.expectError(error.InvalidLda, omatcopy(check.io, f64, .n, 2, 3, 1, &a, big, &b, 2));
    try testing.expectError(error.InvalidLdb, omatcopy(check.io, f64, .n, 2, 3, 1, &a, 2, &b, big));
}

test "omatcopy and imatcopy match the definition for every letter, square and not, in place through a buffer too" {
    const al = testing.allocator;
    inline for (.{ f32, f64, f16, f128, types.Complex32, types.Complex64, std.math.Complex(f16) }) |T| {
        const complex = comptime isComplex(T);
        for ([_]Trans{ .n, .t, .r, .c }) |tr| for ([_][2]usize{ .{ 4, 4 }, .{ 5, 3 }, .{ 2, 6 } }) |s| for ([_][2]usize{ .{ 0, 0 }, .{ 1, 0 }, .{ 0, 2 } }) |pad| {
            const rows = s[0];
            const cols = s[1];
            const t = tr.transposed();
            const out_rows = if (t) cols else rows;
            const out_cols = if (t) rows else cols;
            // equal leading dimensions run in place where the shape allows
            const lda = rows + pad[0];
            const ldb = out_rows + pad[1];
            const len = @max(lda * cols, ldb * out_cols);
            const a = try al.alloc(T, len);
            defer al.free(a);
            const b = try al.alloc(T, ldb * out_cols);
            defer al.free(b);
            check.fill(T, a, rows + 3 * cols);
            const alpha: T = if (complex) T.init(0.5, -2) else -0.5;
            const want = try al.alloc(T, ldb * out_cols);
            defer al.free(want);
            for (0..out_cols) |j| for (0..out_rows) |i| {
                var x = if (t) a[j + i * lda] else a[i + j * lda];
                if (complex and tr.conjA()) x = x.conjugate();
                want[i + j * ldb] = reference.mulW(alpha, x);
            };
            try omatcopy(check.io, T, tr, @intCast(rows), @intCast(cols), alpha, a.ptr, @intCast(lda), b.ptr, @intCast(ldb));
            for (0..out_cols) |j| for (0..out_rows) |i| {
                try check.expectClose(T, want[i + j * ldb], b[i + j * ldb], check.tolerance(Float(T), 2));
            };
            try imatcopy(al, check.io, T, tr, @intCast(rows), @intCast(cols), alpha, a.ptr, @intCast(lda), @intCast(ldb));
            for (0..out_cols) |j| for (0..out_rows) |i| {
                try testing.expect(bitwise.same(T, b[i + j * ldb], a[i + j * ldb], true));
            };
        };
    }
}

test "laswp follows dlaswp.f, and reports pivots and rows outside A" {
    const al = testing.allocator;
    inline for (.{ f32, types.Complex64, f128 }) |T| {
        for ([_]isize{ 1, -1, 2, -2 }) |incx| {
            const m = 9;
            const n = 4;
            const lda = m + 2;
            var ipiv: [2 * m]isize = undefined;
            ref.pivots(&ipiv, m, @abs(incx));
            const a = try al.alloc(T, lda * n);
            defer al.free(a);
            check.fill(T, a, 4);
            const want = try al.dupe(T, a);
            defer al.free(want);
            ref.laswpRef(T, n, want, lda, 2, 8, &ipiv, incx);
            try laswp(check.io, T, n, a.ptr, lda, 2, 8, &ipiv, incx);
            for (want, a) |w, g| try testing.expect(bitwise.same(T, w, g, true));
        }
    }
    var a: [12]f64 = @splat(1);
    var ipiv = [_]isize{ 1, 2, 3, 4 };
    try laswp(check.io, f64, 0, &a, 0, 0, 9, &ipiv, 1);
    try laswp(check.io, f64, 3, &a, 0, 0, 9, &ipiv, 0);
    try laswp(check.io, f64, 3, &a, 0, 3, 2, &ipiv, 1);
    try testing.expectError(error.InvalidLda, laswp(check.io, f64, 3, &a, 3, 1, 4, &ipiv, 1));
    try testing.expectError(error.InvalidK1, laswp(check.io, f64, 3, &a, 4, 0, 4, &ipiv, 1));
    try testing.expectError(error.InvalidIncx, laswp(check.io, f64, 3, &a, 4, 1, 4, &ipiv, std.math.maxInt(isize) / 2));
    ipiv[2] = 5;
    try testing.expectError(error.InvalidIpiv, laswp(check.io, f64, 3, &a, 4, 1, 4, &ipiv, 1));
    ipiv[2] = 0;
    try testing.expectError(error.InvalidIpiv, laswp(check.io, f64, 3, &a, 4, 1, 4, &ipiv, -1));
    for (a) |v| try testing.expectEqual(1, v);
}
