//! Hessenberg reduction, from faer v0.24.4
//! (<https://codeberg.org/sarah-quinones/faer>, tag `faer-v0.24.4`, MIT):
//! `faer/src/linalg/evd/hessenberg.rs`. `A = Q H Qᴴ` in place: `H` on and
//! above the subdiagonal, `Q`'s reflections below it (factors in `h`,
//! `householder.zig`), `Q = H(0) ...` acting on rows `1..n`.
//!
//! `hessenberg_rearranged_unblocked`: one pass per column, its two-sided
//! update deferred into the next column's products, `A22 -= u yᴴ + z uᴴ`
//! fused with `x A22` and `A22 x` (faer's `hessenberg_fused_op`). Faer sums
//! that pass's `A22 x` over threads, so its bits follow the thread count;
//! here it is serial, faer's `Par::Seq` path. Faer's blocked form
//! (`hessenberg_gqvdg_blocked`, from 256 up) is not ported yet (todo.md);
//! the unblocked one serves every size.

const std = @import("std");
const ops = @import("ops.zig");
const hh = @import("householder.zig");
const tridiag = @import("tridiag.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = ops.Mat;
const Float = ops.Float;
const Error = ops.Error;

/// `hessenberg_in_place`: `a` is `n x n`, `h` `bs x (n - 1)` (`bs > 0`).
pub fn hessenberg(al: Allocator, io: Io, comptime T: type, a: Mat(T), h: Mat(T)) Error!void {
    const n = a.rows;
    if (n == 0) return;
    const buf = try al.alloc(T, 4 * n);
    defer al.free(buf);
    const y = buf[0..n];
    const z = buf[n .. 2 * n];
    const v = buf[2 * n .. 3 * n];
    const w = buf[3 * n ..];
    for (0..n) |k| {
        const rem = n - k - 1;
        const a11 = a.at(k, k);
        const a12 = a.col(k + 1) + k; // stride ld
        const a21 = a.col(k) + k + 1;
        const ld = a.ld;
        const y2 = y.ptr + k + 1;
        const z2 = z.ptr + k + 1;
        if (k > 0) {
            const y1 = y[k];
            const z1 = z[k];
            const up = a.col(k - 1) + k + 1;
            a11.* = ops.sub(T, a11.*, ops.add(T, y1, z1));
            for (0..rem) |j| a12[j * ld] = ops.sub(T, a12[j * ld], ops.add(T, y2[j], ops.mul(T, z1, ops.conj(T, up[j]))));
            for (0..rem) |i| a21[i] = ops.sub(T, a21[i], ops.add(T, ops.mul(T, up[i], y1), z2[i]));
        }
        if (k + 1 == n) break;
        const info = hh.make(io, T, &a21[0], rem - 1, a21 + 1, 1, null, 0);
        const tau_inv = 1 / info.tau;
        const beta = a21[0];
        a21[0] = ops.one(T);
        h.at(0, k).* = ops.scalar(T, info.tau);
        const x2 = a21;
        const a22 = a.sub(k + 1, k + 1, rem, rem);
        if (k > 0) {
            const up = a.col(k - 1) + k + 1;
            fused(T, a22, v.ptr + k + 1, w.ptr + k + 1, x2, up, z2, y2);
            @memcpy(y2[0..rem], v[k + 1 ..][0..rem]);
            @memcpy(z2[0..rem], w[k + 1 ..][0..rem]);
        } else {
            // z2 := A22 x2; y2 := x2ᴴ A22
            try ops.gemv(al, io, T, .n, ops.one(T), a22, x2, 1, ops.zero(T), z2, 1);
            try ops.gemv(al, io, T, .c, ops.one(T), a22, x2, 1, ops.zero(T), y2, 1);
            if (comptime ops.isComplex(T)) for (y2[0..rem]) |*e| {
                e.* = ops.conj(T, e.*);
            };
        }
        // b := x2ᴴ z2 / (2 tau); y2 := (y2 - b x2ᴴ) / tau; z2 := (z2 - b x2) / tau
        var d = ops.zero(T);
        for (0..rem) |i| d = ops.add(T, d, ops.mul(T, ops.conj(T, x2[i]), z2[i]));
        const b = ops.mulReal(T, d, 0.5 * tau_inv);
        for (0..rem) |j| y2[j] = ops.mulReal(T, ops.sub(T, y2[j], ops.mul(T, b, ops.conj(T, x2[j]))), tau_inv);
        for (0..rem) |i| z2[i] = ops.mulReal(T, ops.sub(T, z2[i], ops.mul(T, b, x2[i])), tau_inv);
        // A12 := A12 (I - x2 x2ᴴ / tau)
        var e = ops.zero(T);
        for (0..rem) |j| e = ops.add(T, e, ops.mul(T, a12[j * ld], x2[j]));
        e = ops.mulReal(T, e, tau_inv);
        for (0..rem) |j| a12[j * ld] = ops.sub(T, a12[j * ld], ops.mul(T, e, ops.conj(T, x2[j])));
        // A02 := A02 (I - x2 x2ᴴ / tau)
        if (k > 0) {
            const a02 = a.sub(0, k + 1, k, rem);
            try ops.gemv(al, io, T, .n, ops.one(T), a02, x2, 1, ops.zero(T), w.ptr, 1);
            for (0..rem) |j| {
                const c = ops.mulReal(T, ops.conj(T, x2[j]), -tau_inv);
                for (0..k) |i| a02.at(i, j).* = ops.add(T, a02.get(i, j), ops.mul(T, w[i], c));
            }
        }
        a21[0] = beta;
    }
    if (n > 1) try tridiag.sequenceFactors(al, io, T, a.sub(1, 0, n - 1, n - 1), h);
}

/// `hessenberg_fused_op`: per column `j` of `a`, `a -= l0 r0(j) + l1
/// conj(r1(j))`, then `lout(j) := l_inᴴ a` and `rout += a r_in(j)`; here
/// `l_in = r_in = x`, `l0 = r1 = u`, `l1 = z`, `r0 = y`. Partial sums in
/// fours, as faer's `simd_iter!` keeps them.
fn fused(comptime T: type, a: Mat(T), lout: [*]T, rout: [*]T, x: [*]const T, u: [*]const T, z: [*]const T, y: [*]const T) void {
    const m = a.rows;
    @memset(rout[0..m], ops.zero(T));
    for (0..a.cols) |j| {
        const col = a.col(j);
        const yj = y[j];
        const uj = ops.conj(T, u[j]);
        const xj = x[j];
        var acc: [4]T = @splat(ops.zero(T));
        for (0..m) |i| {
            const v = ops.sub(T, ops.sub(T, col[i], ops.mul(T, u[i], yj)), ops.mul(T, uj, z[i]));
            col[i] = v;
            acc[i % 4] = ops.add(T, acc[i % 4], ops.mul(T, ops.conj(T, x[i]), v));
            rout[i] = ops.add(T, rout[i], ops.mul(T, v, xj));
        }
        lout[j] = ops.add(T, ops.add(T, acc[0], acc[1]), ops.add(T, acc[2], acc[3]));
    }
}

// ---- tests -----------------------------------------------------------------------

const testing = std.testing;
const check = @import("../blas/_test/blas1/check.zig");
const types = @import("../blas/types.zig");

test "faer hessenberg: Qᴴ A Q = H, real and complex" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64, f32, types.Complex32 }) |T| {
        const F = Float(T);
        for ([_]usize{ 1, 2, 3, 5, 8, 21, 60 }) |n| {
            const a0 = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, a0);
            check.fill(T, a0.ptr[0 .. n * n], 5 * n);
            const hm = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, hm);
            hm.copyFrom(a0);
            const b = 4;
            const h = try ops.alloc(T, al, b, @max(n, 1) - 1);
            defer ops.free(T, al, h);
            try hessenberg(al, check.io, T, hm, h);
            const x = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, x);
            x.copyFrom(a0);
            if (n > 1) {
                const basis = hm.sub(1, 0, n - 1, n - 1);
                try hh.applySequenceTransposed(al, check.io, T, basis, h, true, x.sub(1, 0, n - 1, n));
                try hh.applySequenceRight(al, check.io, T, basis, h, false, false, x.sub(0, 1, n, n - 1));
            }
            const tol = 64 * std.math.floatEps(F) * @as(F, @floatFromInt(n + 1)) * @max(ops.normMax(T, a0), 1);
            for (0..n) |j| for (0..n) |i| {
                const want = if (i <= j + 1) hm.get(i, j) else ops.zero(T);
                const err = ops.abs(T, ops.sub(T, x.get(i, j), want));
                if (!(err <= tol)) {
                    check.say("hessenberg {s} n={d}: ({d},{d}) off by {e}\n", .{ @typeName(T), n, i, j, err });
                    return error.TestUnexpectedResult;
                }
            };
        }
    }
}
