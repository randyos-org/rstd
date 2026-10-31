//! Bidiagonalization, from faer v0.24.4
//! (<https://codeberg.org/sarah-quinones/faer>, tag `faer-v0.24.4`, MIT):
//! `faer/src/linalg/svd/bidiag.rs`. `A = U B Vᴴ` in place: `B` upper
//! bidiagonal on `A`'s diagonal and superdiagonal, `U`'s reflections below
//! the diagonal (factors in `hl`, `householder.zig`), `V`'s in the rows right
//! of the superdiagonal (factors in `hr`), each row read as a column.
//!
//! One pass per column pair, its rank-2 update deferred into the next
//! column's products (the "BLAS 2.5" scheme of Howell, Demmel, Fulton,
//! Hammarling and Marmol, "Cache efficient bidiagonalization using BLAS 2.5
//! operators", ACM TOMS 34(3), 2008): `A22 -= up y + z vp` fused with
//! `y := uᴴ A22`, faer's `bidiag_fused_op`. That pass is ported as is and
//! split by column chunks on `io` (columns are independent, so any worker
//! count gives the same bits); its other products are `gemv`.

const std = @import("std");
const ops = @import("ops.zig");
const hh = @import("householder.zig");
const threading = @import("../blas/threading.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = ops.Mat;
const Float = ops.Float;
const isComplex = ops.isComplex;
const Error = ops.Error;

/// `bidiag_in_place`: `hl` is `bl x min(m, n)`, `hr` `br x (min(m, n) - 1)`
/// (`bl`, `br > 0`).
pub fn bidiagonalize(al: Allocator, io: Io, comptime T: type, a: Mat(T), hl: Mat(T), hr: Mat(T)) Error!void {
    const m = a.rows;
    const n = a.cols;
    const size = @min(m, n);
    const y = try al.alloc(T, n);
    defer al.free(y);
    const z = try al.alloc(T, m);
    defer al.free(z);
    for (0..size) |k| {
        const a11 = a.at(k, k);
        const r = m - k - 1;
        const c = n - k - 1;
        const a21 = a.col(k) + k + 1;
        const a12 = a.col(k + 1) + k; // stride lda
        const a22 = a.sub(k + 1, k + 1, r, c);
        const y2 = y.ptr + k + 1;
        const z2 = z.ptr + k + 1;
        if (k > 0) {
            const y1 = y[k];
            const z1 = z[k];
            const up0 = a.get(k, k - 1);
            const up = a.col(k - 1) + k + 1;
            const vp = a.col(k + 1) + k - 1; // stride lda
            a11.* = ops.sub(T, a11.*, ops.add(T, ops.mul(T, up0, y1), z1));
            for (0..r) |i| a21[i] = ops.sub(T, a21[i], ops.add(T, ops.mul(T, up[i], y1), z2[i]));
            for (0..c) |j| a12[j * a.ld] = ops.sub(T, a12[j * a.ld], ops.add(T, ops.mul(T, up0, y2[j]), ops.mul(T, z1, vp[j * a.ld])));
        }
        const left = hh.make(io, T, a11, r, a21, 1, null, 0);
        const tl_inv = 1 / left.tau;
        hl.at(0, k).* = ops.scalar(T, left.tau);
        if (k > 0) {
            try fused(io, T, a22, a21, a.col(k - 1) + k + 1, z2, y2, a.col(k + 1) + k - 1, a.ld);
        } else {
            // y2 := A21ᴴ A22
            try ops.gemv(al, io, T, .u, ops.one(T), a22, a21, 1, ops.zero(T), y2, 1);
        }
        for (0..c) |j| {
            y2[j] = ops.mulReal(T, ops.add(T, y2[j], a12[j * a.ld]), tl_inv);
            a12[j * a.ld] = ops.sub(T, a12[j * a.ld], y2[j]);
        }
        const norm = ops.nrm2(io, T, c, a12, a.ld);
        if (norm >= std.math.floatMin(Float(T))) {
            const inv = 1 / norm;
            for (0..c) |j| a12[j * a.ld] = ops.mulReal(T, a12[j * a.ld], inv);
        } else if (norm != 0) {
            // faer multiplies by `1 / norm`, which overflows here
            for (0..c) |j| a12[j * a.ld] = divReal(T, a12[j * a.ld], norm);
        }
        // z2 := A22 conj(A12)
        try ops.gemv(al, io, T, .o, ops.one(T), a22, a12, a.ld, ops.zero(T), z2, 1);
        if (k + 1 == size) break;
        const right = hh.make(io, T, &a12[0], c - 1, a12 + a.ld, a.ld, null, 0);
        const tr_inv = 1 / right.tau;
        hr.at(0, k).* = ops.scalar(T, right.tau);
        const beta = a12[0];
        a12[0] = ops.mulReal(T, a12[0], norm);
        var b = y2[0];
        for (1..c) |j| b = ops.add(T, b, ops.mul(T, y2[j], ops.conj(T, a12[j * a.ld])));
        const a22a = a22.col(0);
        if (right.head_with_beta_inv) |hw| {
            for (0..r) |i| {
                var w = ops.sub(T, z2[i], ops.mul(T, a22a[i], ops.conj(T, beta)));
                w = ops.mul(T, w, ops.conj(T, hw));
                w = ops.sub(T, w, ops.mul(T, a21[i], b));
                z2[i] = ops.mulReal(T, w, tr_inv);
            }
        } else {
            for (0..r) |i| z2[i] = ops.mulReal(T, ops.sub(T, a22a[i], ops.mul(T, a21[i], b)), tr_inv);
        }
    }
    // the block factors: U's from its columns, V's from its rows, copied
    var j: usize = 0;
    while (j < size) : (j += hl.rows) {
        const b = @min(hl.rows, size - j);
        for (0..b) |k| hl.at(k, j + k).* = hl.get(0, j + k);
        try hh.blockFactor(al, io, T, a.sub(j, j, m - j, b), hl.sub(0, j, b, b));
    }
    if (size > 1) {
        const rsize = size - 1;
        const w = try ops.alloc(T, al, n - 1, @min(hr.rows, rsize));
        defer ops.free(T, al, w);
        j = 0;
        while (j < rsize) : (j += hr.rows) {
            const b = @min(hr.rows, rsize - j);
            for (0..b) |k| hr.at(k, j + k).* = hr.get(0, j + k);
            // V's columns j..j+b: rows j..j+b of A(:, 1:), from column j
            const v = w.sub(0, 0, n - 1 - j, b);
            for (0..b) |cc| for (0..n - 1 - j) |rr| {
                v.at(rr, cc).* = a.get(j + cc, 1 + j + rr);
            };
            try hh.blockFactor(al, io, T, v, hr.sub(0, j, b, b));
        }
    }
}

fn divReal(comptime T: type, x: T, d: Float(T)) T {
    return if (comptime isComplex(T)) .init(x.re / d, x.im / d) else x / d;
}

/// `bidiag_fused_op`: per column `j` of `a22`, `a := a - up y(j) - z vp(j)`,
/// then `y(j) := uᴴ a`; `vp` at stride `ldv`.
fn fused(io: Io, comptime T: type, a22: Mat(T), u: [*]const T, up: [*]const T, z: [*]const T, y: [*]T, vp: [*]const T, ldv: usize) Io.Cancelable!void {
    if (a22.cols == 0) return;
    const Ctx = struct {
        a: Mat(T),
        u: [*]const T,
        up: [*]const T,
        z: [*]const T,
        y: [*]T,
        vp: [*]const T,
        ldv: usize,

        fn run(ctx: *const @This(), ch: threading.Chunk) void {
            for (ch.start..ch.end) |j| ctx.y[j] = column(T, ctx.a.rows, ctx.a.col(j), ctx.u, ctx.up, ctx.z, ops.neg(T, ctx.y[j]), ops.neg(T, ctx.vp[j * ctx.ldv]));
        }
    };
    const ctx: Ctx = .{ .a = a22, .u = u, .up = up, .z = z, .y = y, .vp = vp, .ldv = ldv };
    const p = threading.plan(.gemv, T, threading.product(&.{ a22.rows, a22.cols }), a22.cols, .even, 1);
    try threading.forEach(io, &p, &ctx, Ctx.run);
}

/// One column of `fused`: `a += up yj + z vj`, then `uᴴ a` in four partial
/// sums, as faer's `simd_iter!` over `[indices; 4]` keeps them.
fn column(comptime T: type, m: usize, a: [*]T, u: [*]const T, up: [*]const T, z: [*]const T, yj: T, vj: T) T {
    var acc: [4]T = @splat(ops.zero(T));
    for (0..m) |i| {
        var x = ops.add(T, a[i], ops.mul(T, up[i], yj));
        x = ops.add(T, x, ops.mul(T, z[i], vj));
        a[i] = x;
        acc[i % 4] = ops.add(T, acc[i % 4], ops.mul(T, ops.conj(T, u[i]), x));
    }
    return ops.add(T, ops.add(T, acc[0], acc[1]), ops.add(T, acc[2], acc[3]));
}

// ---- tests -----------------------------------------------------------------------

const testing = std.testing;
const check = @import("../blas/_test/blas1/check.zig");
const types = @import("../blas/types.zig");

test "faer bidiag: Uᴴ A V = B, real and complex, square and tall" {
    // faer's `test_bidiag_real`, `test_bidiag_cplx`: (8, 4) and (8, 8),
    // `bl = 4`, `br = 3`; and larger, past the block sizes
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64, f32, types.Complex32 }) |T| {
        const F = Float(T);
        for ([_][2]usize{ .{ 8, 4 }, .{ 8, 8 }, .{ 37, 21 }, .{ 30, 30 } }) |mn| {
            const m = mn[0];
            const n = mn[1];
            const size = @min(m, n);
            const a0 = try ops.alloc(T, al, m, n);
            defer ops.free(T, al, a0);
            check.fill(T, a0.ptr[0 .. m * n], m + 3 * n);
            const uv = try ops.alloc(T, al, m, n);
            defer ops.free(T, al, uv);
            uv.copyFrom(a0);
            const hl = try ops.alloc(T, al, 4, size);
            defer ops.free(T, al, hl);
            const hr = try ops.alloc(T, al, 3, @max(size, 1) - 1);
            defer ops.free(T, al, hr);
            try bidiagonalize(al, check.io, T, uv, hl, hr);
            // A := Uᴴ A, then A(:, 1:size) := A(:, 1:size) V
            const x = try ops.alloc(T, al, m, n);
            defer ops.free(T, al, x);
            x.copyFrom(a0);
            try hh.applySequenceTransposed(al, check.io, T, uv.sub(0, 0, m, size), hl, true, x);
            if (size > 1) {
                // V's reflections as columns, then X Vᴴ... as (V Xᴴ)ᴴ
                const vb = try ops.alloc(T, al, n - 1, size - 1);
                defer ops.free(T, al, vb);
                for (0..size - 1) |cc| for (0..n - 1) |rr| {
                    vb.at(rr, cc).* = uv.get(cc, 1 + rr);
                };
                const xt = try ops.alloc(T, al, n - 1, m);
                defer ops.free(T, al, xt);
                for (0..m) |i| for (0..n - 1) |jj| {
                    xt.at(jj, i).* = ops.conj(T, x.get(i, 1 + jj));
                };
                // V is the rows' reflections conjugated, so (X V)ᴴ = Vᴴ Xᴴ is
                // their sequence transposed
                try hh.applySequenceTransposed(al, check.io, T, vb, hr, false, xt);
                for (0..m) |i| for (0..n - 1) |jj| {
                    x.at(i, 1 + jj).* = ops.conj(T, xt.get(jj, i));
                };
            }
            const tol = 64 * std.math.floatEps(F) * @as(F, @floatFromInt(m + n));
            for (0..n) |j| for (0..m) |i| {
                const want = if (i == j or i + 1 == j) uv.get(i, j) else ops.zero(T);
                const err = ops.abs(T, ops.sub(T, x.get(i, j), want));
                if (!(err <= tol)) {
                    check.say("bidiag {s} {d}x{d}: ({d},{d}) off by {e}\n", .{ @typeName(T), m, n, i, j, err });
                    return error.TestUnexpectedResult;
                }
            };
        }
    }
}
