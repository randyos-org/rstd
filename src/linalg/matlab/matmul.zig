//! Matrix multiplication (`X*Y`), ported from `matlab88.f90`'s `STACK2`
//! (the `*` operator's `MULTIPLICATION` case, label 05) - REPL operator
//! dispatch, not a named Fortran subroutine, so there's no `WMATMUL` to
//! point at. This file is the extracted math: the row-column double loop
//! (via `wdotur`/`wdotui`) and the dimension-compatibility rule. The
//! surrounding expression-stack bookkeeping (`STACK2`'s own
//! `L`/`L2`/`LL`/`TOP`/`BOT`/`NSTK`/`ERROR`-code plumbing) is left out as
//! REPL-only, and the original's scalar-times-anything special case
//! (`M2*N2==1`/`M*N==1`) isn't reproduced here either - that's exactly
//! `blas1.wscal`/`wrscal` applied to the whole flattened matrix
//! (`wscal(m*n, s, ar, ai, 1)`), already covered, not a matmul concern.

const blas1 = @import("blas1.zig");
const CxMat = @import("CxMat.zig");

/// `C = A*B`: `A` is `m`x`n`, `B` is `n`x`p`, `C` is `m`x`p`. The shared
/// `a.n == b.m` requirement *is* the dimension-compatibility rule (checked
/// by the caller constructing matching `CxMat`s) rather than a runtime
/// check on four independent dimensions the way the original Fortran
/// needed.
///
/// `out` must not alias `a`/`b`. `scratch_r`/`scratch_i` (length >= `a.n`)
/// are caller-owned workspace, not returned data.
///
/// A row of `A` is strided (not contiguous) in this column-major layout,
/// unlike a column of `B` - so each row of `A` is gathered into
/// `scratch_r`/`scratch_i` once and reused across all `p` columns of `B`,
/// rather than re-gathering it per output element (`m` gathers total, not
/// `m*p`).
pub fn matmul(a: CxMat, b: CxMat, out: CxMat, scratch_r: []f64, scratch_i: []f64) void {
    const m = a.m;
    const n = a.n;
    const p = b.n;

    var i: usize = 1;
    while (i <= m) : (i += 1) {
        var k: usize = 1;
        while (k <= n) : (k += 1) {
            scratch_r[k - 1] = a.r[a.idx(i, k)];
            scratch_i[k - 1] = a.i[a.idx(i, k)];
        }

        var j: usize = 1;
        while (j <= p) : (j += 1) {
            const col_r = b.r[b.idx(1, j)..];
            const col_i = b.i[b.idx(1, j)..];
            out.r[out.idx(i, j)] = blas1.wdotur(n, scratch_r[0..n], scratch_i[0..n], 1, col_r, col_i, 1);
            out.i[out.idx(i, j)] = blas1.wdotui(n, scratch_r[0..n], scratch_i[0..n], 1, col_r, col_i, 1);
        }
    }
}

const std = @import("std");
const zigstd = @import("../../zigstd.zig");
const ix = blas1.ix;

test "matmul matches hand computation for a small real case" {
    // A = [[1,2],[3,4]] (2x2), B = [[5,6],[7,8]] (2x2).
    // A*B = [[1*5+2*7, 1*6+2*8], [3*5+4*7, 3*6+4*8]] = [[19,22],[43,50]].
    var ar = [_]f64{ 1.0, 3.0, 2.0, 4.0 };
    var ai: [4]f64 = @splat(0.0);
    var br = [_]f64{ 5.0, 7.0, 6.0, 8.0 };
    var bi: [4]f64 = @splat(0.0);
    var cr: [4]f64 = undefined;
    var ci: [4]f64 = undefined;
    var scratch_r: [2]f64 = undefined;
    var scratch_i: [2]f64 = undefined;

    matmul(CxMat.initPacked(&ar, &ai, 2, 2), CxMat.initPacked(&br, &bi, 2, 2), CxMat.initPacked(&cr, &ci, 2, 2), &scratch_r, &scratch_i);

    try zigstd.testExpectApproxEqAbs(19.0, cr[ix(2, 1, 1)], 1e-12);
    try zigstd.testExpectApproxEqAbs(43.0, cr[ix(2, 2, 1)], 1e-12);
    try zigstd.testExpectApproxEqAbs(22.0, cr[ix(2, 1, 2)], 1e-12);
    try zigstd.testExpectApproxEqAbs(50.0, cr[ix(2, 2, 2)], 1e-12);
    for (ci) |v| try zigstd.testExpectApproxEqAbs(0.0, v, 1e-12);
}

test "matmul handles non-square, complex operands and matches a naive reference loop" {
    const m = 3;
    const n = 4;
    const p = 2;

    var prng = std.Random.DefaultPrng.init(42);
    const rand = prng.random();
    const alloc = zigstd.testAllocator;

    const ar = try alloc.alloc(f64, m * n);
    defer alloc.free(ar);
    const ai = try alloc.alloc(f64, m * n);
    defer alloc.free(ai);
    const br = try alloc.alloc(f64, n * p);
    defer alloc.free(br);
    const bi = try alloc.alloc(f64, n * p);
    defer alloc.free(bi);
    for (ar, ai) |*r, *im| {
        r.* = rand.float(f64) * 2 - 1;
        im.* = rand.float(f64) * 2 - 1;
    }
    for (br, bi) |*r, *im| {
        r.* = rand.float(f64) * 2 - 1;
        im.* = rand.float(f64) * 2 - 1;
    }

    const cr = try alloc.alloc(f64, m * p);
    defer alloc.free(cr);
    const ci = try alloc.alloc(f64, m * p);
    defer alloc.free(ci);
    const scratch_r = try alloc.alloc(f64, n);
    defer alloc.free(scratch_r);
    const scratch_i = try alloc.alloc(f64, n);
    defer alloc.free(scratch_i);

    matmul(CxMat.initPacked(ar, ai, m, n), CxMat.initPacked(br, bi, n, p), CxMat.initPacked(cr, ci, m, p), scratch_r, scratch_i);

    var i: usize = 1;
    while (i <= m) : (i += 1) {
        var j: usize = 1;
        while (j <= p) : (j += 1) {
            var sr: f64 = 0.0;
            var si: f64 = 0.0;
            var k: usize = 1;
            while (k <= n) : (k += 1) {
                const a_re = ar[ix(m, i, k)];
                const a_im = ai[ix(m, i, k)];
                const b_re = br[ix(n, k, j)];
                const b_im = bi[ix(n, k, j)];
                sr += a_re * b_re - a_im * b_im;
                si += a_re * b_im + a_im * b_re;
            }
            try zigstd.testExpectApproxEqAbs(sr, cr[ix(m, i, j)], 1e-9);
            try zigstd.testExpectApproxEqAbs(si, ci[ix(m, i, j)], 1e-9);
        }
    }
}
