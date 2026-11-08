//! `X*Y` (`STACK2`, label 05) by `gemm` (M2); a scalar times anything is
//! `CxMat.scale`.

const std = @import("std");
const scalar = @import("scalar.zig");
const CxMat = @import("CxMat.zig");
const blas = @import("../blas/__root__.zig");
const ops = @import("../lapack/ops.zig");
const Cx = scalar.Cx;
const Allocator = std.mem.Allocator;
const Io = std.Io;
const ix = scalar.ix;

pub const Error = Allocator.Error || Io.Cancelable;

/// `out = A*B` by `gemm`, real in place for real operands, else over
/// interleaved copies; the caller matches the shapes, `out` not aliasing.
pub fn matmul(a: CxMat, b: CxMat, out: CxMat, allocator: Allocator, io: Io) Error!void {
    const m: isize = @intCast(a.m);
    const n: isize = @intCast(a.n);
    const p: isize = @intCast(b.n);
    if (a.isReal() and b.isReal()) {
        try multiplied(blas.gemm(allocator, io, f64, .n, .n, m, p, n, 1, a.r.ptr, @intCast(@max(1, a.lda)), b.r.ptr, @intCast(@max(1, b.lda)), 0, out.r.ptr, @intCast(@max(1, out.lda))));
        for (1..out.n + 1) |j| @memset(out.col(j).i, 0);
        return;
    }
    const x = try a.toPacked(Cx, allocator);
    defer allocator.free(x);
    const y = try b.toPacked(Cx, allocator);
    defer allocator.free(y);
    const z = try allocator.alloc(Cx, out.m * out.n);
    defer allocator.free(z);
    try multiplied(blas.gemm(allocator, io, Cx, .n, .n, m, p, n, Cx.init(1, 0), x.ptr, @max(1, m), y.ptr, @max(1, n), Cx.init(0, 0), z.ptr, @max(1, m)));
    out.fromPacked(Cx, z);
}

fn multiplied(r: blas.Level3Error!void) Error!void {
    // the shapes come from the `CxMat`s
    r catch |e| return ops.narrow(Error, e);
}

// ---- tests ------------------------------------------------------------------

const zigstd = @import("../../zigstd.zig");

test "matmul matches hand computation for a small real case" {
    // A = [[1,2],[3,4]] (2x2), B = [[5,6],[7,8]] (2x2).
    // A*B = [[1*5+2*7, 1*6+2*8], [3*5+4*7, 3*6+4*8]] = [[19,22],[43,50]].
    var ar = [_]f64{ 1.0, 3.0, 2.0, 4.0 };
    var ai: [4]f64 = @splat(0.0);
    var br = [_]f64{ 5.0, 7.0, 6.0, 8.0 };
    var bi: [4]f64 = @splat(0.0);
    var cr: [4]f64 = undefined;
    var ci: [4]f64 = undefined;

    try matmul(CxMat.initPacked(&ar, &ai, 2, 2), CxMat.initPacked(&br, &bi, 2, 2), CxMat.initPacked(&cr, &ci, 2, 2), zigstd.testAllocator, std.testing.io);

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

    try matmul(CxMat.initPacked(ar, ai, m, n), CxMat.initPacked(br, bi, n, p), CxMat.initPacked(cr, ci, m, p), alloc, std.testing.io);

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
            try zigstd.testExpectApproxEqAbs(sr, cr[ix(m, i, j)], 1e-14);
            try zigstd.testExpectApproxEqAbs(si, ci[ix(m, i, j)], 1e-14);
        }
    }
}
