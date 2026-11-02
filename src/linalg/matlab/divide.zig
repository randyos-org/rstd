//! Matrix left/right division (REPL `\`/`/`) over `lu.zig`'s `getrf` and
//! `getrs` (M2), in place of WGECO and WGESL.

const std = @import("std");
const scalar = @import("scalar.zig");
const CxMat = @import("CxMat.zig");
const lu_mod = @import("lu.zig");
const Trans = @import("../lapack/__root__.zig").Trans;

const Allocator = std.mem.Allocator;
const Io = std.Io;
const ix = scalar.ix;

/// Right matrix division for a square divisor: solve `X*B = A` for `X` -
/// REPL semantics for `A/B` - as `Bᵀ Xᵀ = Aᵀ`. `a` (`m`x`n`) is
/// overwritten with `X`; `b` (`n`x`n`) is destroyed by factorization.
/// Returns `B`'s `rcond`; an exact zero pivot is `error.Singular`.
pub fn mrdivide(a: CxMat, b: CxMat, allocator: Allocator, io: Io) lu_mod.Error!f64 {
    const tr = try allocator.alloc(f64, a.m * a.n);
    defer allocator.free(tr);
    const ti = try allocator.alloc(f64, a.m * a.n);
    defer allocator.free(ti);
    const t = CxMat.initPacked(tr, ti, a.n, a.m);
    a.transpose(t);
    const rc = try mldivideWith(b, t, .t, allocator, io);
    t.transpose(a);
    return rc;
}

/// Left matrix division for a square `self`: solve `self*X = b` for `X` -
/// REPL semantics for `A\B` (square `A` only; a non-square least-squares
/// solve is `gels`'s, M2's QR group). `self` is destroyed by factorization;
/// `b` (`self.n`x`p`) is overwritten with `X`. Returns `self`'s `rcond`; an
/// exact zero pivot is `error.Singular`.
pub fn mldivide(self: CxMat, b: CxMat, allocator: Allocator, io: Io) lu_mod.Error!f64 {
    return mldivideWith(self, b, .n, allocator, io);
}

fn mldivideWith(a: CxMat, b: CxMat, trans: Trans, allocator: Allocator, io: Io) lu_mod.Error!f64 {
    const anorm = try lu_mod.norm1(a, allocator, io);
    const f = try lu_mod.Lu.init(a, !b.isReal(), allocator, io);
    defer f.deinit(allocator);
    if (f.singular) return error.Singular;
    const rc = try lu_mod.rcond(f, allocator, io, anorm);
    try f.solve(trans, b, allocator, io);
    return rc;
}

const zigstd = @import("../../zigstd.zig");
const testing = std.testing;

test "mldivide solves A*X=B" {
    // A = [[2,0],[1,2]], B = [[1,2],[3,4]] (both square, real).
    var ar = [_]f64{ 2.0, 1.0, 0.0, 2.0 }; // column-major
    var ai: [4]f64 = @splat(0.0);
    var br = [_]f64{ 1.0, 3.0, 2.0, 4.0 };
    var bi: [4]f64 = @splat(0.0);

    const b_orig = br;
    const a_mat = CxMat.initPacked(&ar, &ai, 2, 2);
    const b_mat = CxMat.initPacked(&br, &bi, 2, 2);
    const a_for_check = [_]f64{ 2.0, 1.0, 0.0, 2.0 };
    const rc = try mldivide(a_mat, b_mat, testing.allocator, testing.io);
    try zigstd.testExpectTrue(rc > 0.0);

    // Check A*X == original B.
    for (0..2) |row| {
        for (0..2) |col| {
            var s: f64 = 0.0;
            for (0..2) |k| s += a_for_check[ix(2, row + 1, k + 1)] * br[ix(2, k + 1, col + 1)];
            try zigstd.testExpectApproxEqAbs(b_orig[ix(2, row + 1, col + 1)], s, 1e-9);
        }
    }
}

test "mrdivide solves X*B=A" {
    // A = [[1,2],[3,4]], B = [[2,0],[1,2]] (both square, real).
    var ar = [_]f64{ 1.0, 3.0, 2.0, 4.0 }; // column-major
    var ai: [4]f64 = @splat(0.0);
    var br = [_]f64{ 2.0, 1.0, 0.0, 2.0 };
    var bi: [4]f64 = @splat(0.0);

    const b_for_check = br;
    const a_mat = CxMat.initPacked(&ar, &ai, 2, 2);
    const b_mat = CxMat.initPacked(&br, &bi, 2, 2);
    const rc = try mrdivide(a_mat, b_mat, testing.allocator, testing.io);
    try zigstd.testExpectTrue(rc > 0.0);

    // Check X*B == original A.
    const a_orig = [_]f64{ 1.0, 3.0, 2.0, 4.0 };
    for (0..2) |row| {
        for (0..2) |col| {
            var s: f64 = 0.0;
            for (0..2) |k| s += ar[ix(2, row + 1, k + 1)] * b_for_check[ix(2, k + 1, col + 1)];
            try zigstd.testExpectApproxEqAbs(a_orig[ix(2, row + 1, col + 1)], s, 1e-9);
        }
    }
}

test "divisions mix real and complex operands; a singular divisor fails" {
    // A = [[1,i],[0,2]] (2x2 complex), B = [[2,1],[0,1]] real: A\B and B/A, A*X = B
    const a0r = [_]f64{ 1, 0, 0, 2 };
    const a0i = [_]f64{ 0, 0, 1, 0 };
    const b0r = [_]f64{ 2, 0, 1, 1 };
    const b0i: [4]f64 = @splat(0.0);
    for ([_]bool{ false, true }) |right| {
        var ar = a0r;
        var ai = a0i;
        var br = b0r;
        var bi = b0i;
        const a = CxMat.initPacked(&ar, &ai, 2, 2);
        const b = CxMat.initPacked(&br, &bi, 2, 2);
        // left: X = A \ B; right: X = B / A, so X A = B
        _ = try if (right) mrdivide(b, a, testing.allocator, testing.io) else mldivide(a, b, testing.allocator, testing.io);
        const x = b;
        var cr = a0r;
        var ci = a0i;
        const a_orig = CxMat.initPacked(&cr, &ci, 2, 2);
        for (1..3) |i| {
            for (1..3) |j| {
                var s = scalar.Cx.init(0, 0);
                for (1..3) |k| s = s.add(if (right) x.get(i, k).mul(a_orig.get(k, j)) else a_orig.get(i, k).mul(x.get(k, j)));
                try zigstd.testExpectApproxEqAbs(b0r[ix(2, i, j)], s.re, 1e-12);
                try zigstd.testExpectApproxEqAbs(0.0, s.im, 1e-12);
            }
        }
    }
    var sr = [_]f64{ 1, 2, 2, 4 };
    var si: [4]f64 = @splat(0.0);
    var tr = [_]f64{ 1, 0, 0, 1 };
    var ti: [4]f64 = @splat(0.0);
    try testing.expectError(error.Singular, mldivide(CxMat.initPacked(&sr, &si, 2, 2), CxMat.initPacked(&tr, &ti, 2, 2), testing.allocator, testing.io));
}
