//! `gghrd`: the Hessenberg-triangular reduction of a pencil by Givens
//! rotations, faer v0.24.4's `generalized_hessenberg_unblocked` and its
//! `make_givens` and `rot` (`faer/src/linalg/gevd/gen_hessenberg/mod.rs`,
//! MIT) over Reference-LAPACK's `ILO`/`IHI` and `COMPQ`/`COMPZ`
//! (`lapack-netlib/SRC/dgghrd.f`, OpenBLAS v0.3.34-317). faer's blocked form
//! for large pencils is issue #761.

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const ops = @import("ops.zig");
const interface = @import("interface.zig");
const schur = @import("schur.zig");

const Float = types.Float;
const isComplex = types.isComplex;
const scalar = ops.scalar;
const add = ops.add;
const sub = ops.sub;
const mul = ops.mul;
const conj = ops.conj;
const mulReal = ops.mulReal;

/// LAPACK's `COMPQ`, `COMPZ`: no vectors, from the identity, or accumulated
/// into the caller's (`hseqr`'s `COMPZ`).
pub const SchurVectors = schur.SchurVectors;

/// What `gghrd` returns: `BlockOutOfRange` unless `lo <= hi <= n`, a slice
/// too short (`a`, `b`, or `q`, `z` when used) or a leading dimension too small.
pub const GghrdError = interface.BlockError || interface.ShapeError;

/// `?GGHRD`: reduces the pencil `(A, B)`, `B` upper triangular (its part
/// below the diagonal is zeroed first), to `(H, T) = (Qᴴ A Z, Qᴴ B Z)`
/// with `H` upper Hessenberg and `T` upper triangular, rotating within rows
/// and columns `lo..hi` (0-based, `hi` past the end), outside which `A` must
/// already be upper triangular (as `ggbal` leaves it). `q` and `z` (`n x n`)
/// take `Q` and `Z` (`.from_identity`) or are multiplied by them (`.accumulate`).
pub fn gghrd(comptime T: type, compq: SchurVectors, compz: SchurVectors, n: usize, lo: usize, hi: usize, a: []T, lda: usize, b: []T, ldb: usize, q: []T, ldq: usize, z: []T, ldz: usize) GghrdError!void {
    comptime ops.requireFloat(T, "gghrd");
    try interface.blockIn(lo, hi, n);
    try util.holds(a.len, n, n, lda);
    try util.holds(b.len, n, n, ldb);
    if (compq != .none) try util.holds(q.len, n, n, ldq);
    if (compz != .none) try util.holds(z.len, n, n, ldz);
    if (compq == .from_identity) identity(T, n, q, ldq);
    if (compz == .from_identity) identity(T, n, z, ldz);
    if (n <= 1) return;
    for (0..n - 1) |j| for (j + 1..n) |i| {
        b[i + j * ldb] = scalar(T, 0);
    };
    if (hi < lo + 3) return;
    for (lo..hi - 2) |jcol| {
        var jrow = hi - 1;
        while (jrow >= jcol + 2) : (jrow -= 1) {
            // rows `jrow - 1` and `jrow`: zero `A(jrow, jcol)`
            const g1 = givens(T, a[jrow - 1 + jcol * lda], a[jrow + jcol * lda]);
            a[jrow - 1 + jcol * lda] = g1.r;
            a[jrow + jcol * lda] = scalar(T, 0);
            rotRows(T, g1.c, g1.s, a, lda, jrow - 1, jrow, jcol + 1, n);
            rotRows(T, g1.c, g1.s, b, ldb, jrow - 1, jrow, jrow - 1, n);
            if (compq != .none) rotCols(T, g1.c, conj(T, g1.s), q, ldq, jrow - 1, jrow, 0, n);
            // columns `jrow` and `jrow - 1`: zero the fill-in `B(jrow, jrow - 1)`
            const g2 = givens(T, b[jrow + jrow * ldb], b[jrow + (jrow - 1) * ldb]);
            b[jrow + jrow * ldb] = g2.r;
            b[jrow + (jrow - 1) * ldb] = scalar(T, 0);
            rotCols(T, g2.c, g2.s, a, lda, jrow, jrow - 1, 0, hi);
            rotCols(T, g2.c, g2.s, b, ldb, jrow, jrow - 1, 0, jrow);
            if (compz != .none) rotCols(T, g2.c, g2.s, z, ldz, jrow, jrow - 1, 0, n);
        }
    }
}

fn identity(comptime T: type, n: usize, m: []T, ld: usize) void {
    for (0..n) |j| for (0..n) |i| {
        m[i + j * ld] = scalar(T, if (i == j) 1 else 0);
    };
}

/// A plane rotation `[c s; -conj(s) c]` taking `(f, g)` to `(r, 0)`: faer's
/// `make_givens`.
pub fn Givens(comptime T: type) type {
    return struct { c: Float(T), s: T, r: T };
}

pub fn givens(comptime T: type, f: T, g: T) Givens(T) {
    const F = Float(T);
    if (ops.isZero(T, g)) return .{ .c = 1, .s = scalar(T, 0), .r = f };
    if (ops.isZero(T, f)) {
        const d = ops.abs(T, g);
        return .{ .c = 0, .s = mulReal(T, conj(T, g), 1 / d), .r = scalar(T, d) };
    }
    const f1 = ops.abs(T, f);
    const g1 = ops.abs(T, g);
    const h1: F = std.math.hypot(f1, g1);
    const h1_inv = 1 / h1;
    const c = f1 * h1_inv;
    const r = mulReal(T, f, 1 / c);
    const s = mul(T, conj(T, g), mulReal(T, mulReal(T, r, h1_inv), h1_inv));
    return .{ .c = c, .s = s, .r = r };
}

/// faer's `rot` on rows `x` and `y` of `m`, columns `from..to`: `x := c x +
/// s y`, `y := c y - conj(s) x`.
pub fn rotRows(comptime T: type, c: Float(T), s: T, m: []T, ld: usize, x: usize, y: usize, from: usize, to: usize) void {
    for (from..to) |j| rotPair(T, c, s, &m[x + j * ld], &m[y + j * ld]);
}

/// `rot` on columns `x` and `y` of `m`, rows `from..to`.
pub fn rotCols(comptime T: type, c: Float(T), s: T, m: []T, ld: usize, x: usize, y: usize, from: usize, to: usize) void {
    for (from..to) |i| rotPair(T, c, s, &m[i + x * ld], &m[i + y * ld]);
}

/// `rot` on one pair of entries.
pub fn rotPair(comptime T: type, c: Float(T), s: T, x: *T, y: *T) void {
    const xv = x.*;
    const yv = y.*;
    x.* = add(T, mulReal(T, xv, c), mul(T, yv, s));
    y.* = sub(T, mulReal(T, yv, c), mul(T, xv, conj(T, s)));
}

// ---- tests ---------------------------------------------------------------------

const testing = std.testing;
const check = @import("../blas/_test/check.zig");
const pencil = @import("_test/pencil.zig");
const residual = pencil.residual;
const unitary = pencil.unitary;

test "gghrd: Q H Zᴴ = A and Q T Zᴴ = B, H Hessenberg, T triangular, Q and Z unitary" {
    const al = testing.allocator;
    inline for (.{ f32, f64, types.Complex64 }) |T| for ([_]usize{ 1, 2, 3, 8, 20 }) |n| {
        const a0 = try al.alloc(T, n * n);
        defer al.free(a0);
        const b0 = try al.alloc(T, n * n);
        defer al.free(b0);
        check.fill(T, a0, n + 1);
        check.fill(T, b0, n + 2);
        for (0..n) |j| for (j + 1..n) |i| {
            b0[i + j * n] = scalar(T, 0);
        };
        const a = try al.dupe(T, a0);
        defer al.free(a);
        const b = try al.dupe(T, b0);
        defer al.free(b);
        const q = try al.alloc(T, n * n);
        defer al.free(q);
        const z = try al.alloc(T, n * n);
        defer al.free(z);
        try gghrd(T, .from_identity, .from_identity, n, 0, n, a, n, b, n, q, n, z, n);
        for (0..n) |j| for (j + 1..n) |i| {
            try testing.expect(ops.isZero(T, b[i + j * n]));
            if (i > j + 1) try testing.expect(ops.isZero(T, a[i + j * n]));
        };
        const eps = std.math.floatEps(Float(T));
        const tol = 40 * @as(Float(T), @floatFromInt(n * n)) * eps;
        try testing.expect(unitary(T, n, q) <= tol);
        try testing.expect(unitary(T, n, z) <= tol);
        try testing.expect(residual(T, n, q, a, z, a0) <= tol * 4);
        try testing.expect(residual(T, n, q, b, z, b0) <= tol * 4);
    };
}

test "gghrd: an empty block is the identity's work; a bad block is an error" {
    var a = [_]f64{ 1, 2, 3, 4 };
    var b = [_]f64{ 1, 5, 0, 1 };
    var q: [4]f64 = undefined;
    var z: [4]f64 = undefined;
    try gghrd(f64, .from_identity, .from_identity, 2, 1, 1, &a, 2, &b, 2, &q, 2, &z, 2);
    try testing.expectEqualSlices(f64, &.{ 1, 0, 0, 1 }, &q);
    // `B` below its diagonal is zeroed
    try testing.expectEqual(0, b[1]);
    try testing.expectError(error.BlockOutOfRange, gghrd(f64, .none, .none, 2, 2, 1, &a, 2, &b, 2, &q, 2, &z, 2));
}

test "givens: (f, g) goes to (r, 0)" {
    inline for (.{ f64, types.Complex64 }) |T| for ([_][2]f64{ .{ 3, 4 }, .{ 0, 2 }, .{ 5, 0 }, .{ -1e-300, 1e-300 } }) |p| {
        const f = scalar(T, p[0]);
        const g = if (comptime isComplex(T)) T.init(p[1], p[1] / 3) else p[1];
        const rot = givens(T, f, g);
        var x = f;
        var y = g;
        rotPair(T, rot.c, rot.s, &x, &y);
        try testing.expect(ops.abs1(T, sub(T, x, rot.r)) <= 4 * std.math.floatEps(f64) * ops.abs1(T, rot.r));
        try testing.expect(ops.abs1(T, y) <= 4 * std.math.floatEps(f64) * ops.abs1(T, rot.r));
    };
}
