//! Tridiagonalization of a self-adjoint matrix, from faer v0.24.4
//! (<https://codeberg.org/sarah-quinones/faer>, tag `faer-v0.24.4`, MIT):
//! `faer/src/linalg/evd/tridiag.rs`, its fused pass serial (faer's `Par::Seq`
//! path), so the bits never depend on the worker count.

const std = @import("std");
const ops = @import("ops.zig");
const auxiliary = @import("auxiliary.zig");
const hh = @import("householder.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = ops.Mat;
const Float = ops.Float;
const Error = ops.Error;

/// `tridiag_in_place`: `a` is `n x n` (lower triangle read), `h`
/// `bs x (n - 1)` (`bs > 0`).
pub fn tridiagonalize(al: Allocator, io: Io, comptime T: type, a: Mat(T), h: Mat(T)) Error!void {
    const n = a.rows;
    if (n == 0) return;
    const buf = try al.alloc(T, 3 * n);
    defer al.free(buf);
    const y = buf[0..n];
    const w = buf[n .. 2 * n];
    const z = buf[2 * n ..];
    for (0..n) |k| {
        const a11 = a.at(k, k);
        const a21 = a.col(k) + k + 1;
        const rem = n - k - 1;
        if (k > 0) {
            // A11 -= y1 + conj(y1); A21 -= conj(y1) up + y2
            const y1 = y[k];
            const up = a.col(k - 1) + k + 1;
            a11.* = ops.sub(T, a11.*, ops.add(T, y1, ops.conj(T, y1)));
            for (0..rem) |i| a21[i] = ops.sub(T, a21[i], ops.add(T, ops.mul(T, ops.conj(T, y1), up[i]), y[k + 1 + i]));
        }
        if (k + 1 == n) break;
        // the reflection taking A21 to beta e0; x2 its essential part
        const info = hh.make(io, T, &a21[0], rem - 1, a21 + 1, 1, null, 0);
        const tau_inv = 1 / info.tau;
        h.at(0, k).* = ops.scalar(T, info.tau);
        const x2 = a21 + 1;
        const m = rem - 1;
        const y1 = &y[k + 1];
        const y2 = y.ptr + k + 2;
        const z2 = z.ptr + k + 2;
        const w2 = w.ptr + k + 2;
        const ap11 = a.at(k + 1, k + 1);
        const ap21 = a.col(k + 1) + k + 2;
        const a22 = a.sub(k + 2, k + 2, m, m);
        if (k > 0) {
            const up1 = a.get(k + 1, k - 1);
            const up = a.col(k - 1) + k + 2;
            ap11.* = ops.sub(T, ap11.*, ops.add(T, ops.mul(T, up1, ops.conj(T, y1.*)), ops.mul(T, y1.*, ops.conj(T, up1))));
            for (0..m) |i| ap21[i] = ops.sub(T, ap21[i], ops.add(T, ops.mul(T, up[i], ops.conj(T, y1.*)), ops.mul(T, y2[i], ops.conj(T, up1))));
            @memcpy(w2[0..m], y2[0..m]);
            fused(T, a22, y2, z2, w2, up, x2, tau_inv);
        } else {
            fused(T, a22, y2, z2, null, null, x2, tau_inv);
        }
        for (0..m) |i| y2[i] = ops.add(T, ops.add(T, y2[i], z2[i]), ops.mulReal(T, ap21[i], tau_inv));
        // y1 := (a11 + A21ᴴ x2) / tau; b := (y1 + x2ᴴ y2) / (2 tau)
        var d = ap11.*;
        for (0..m) |i| d = ops.add(T, d, ops.mul(T, ops.conj(T, ap21[i]), x2[i]));
        y1.* = ops.mulReal(T, d, tau_inv);
        var e = y1.*;
        for (0..m) |i| e = ops.add(T, e, ops.mul(T, ops.conj(T, x2[i]), y2[i]));
        const b = ops.mulReal(T, e, 0.5 * tau_inv);
        y1.* = ops.sub(T, y1.*, b);
        for (0..m) |i| y2[i] = ops.sub(T, y2[i], ops.mul(T, b, x2[i]));
    }
    // the block factors, from the reflections below the subdiagonal
    if (n > 1) try hh.sequenceFactors(al, io, T, a.sub(1, 0, n - 1, n - 1), h);
}

/// `tridiag_fused_op` over the lower triangle of `a` (`m x m`): with `u`
/// and `w`, first `A -= u wᴴ + w uᴴ`; then `z := f L x` and
/// `y := f Lₛᴴ x`, `L` the lower triangle and `Lₛ` the strict one, so
/// `y + z = f A x`. Partial sums in fours, as faer's `simd_iter!` keeps
/// them.
fn fused(comptime T: type, a: Mat(T), y: [*]T, z: [*]T, w: ?[*]const T, u: ?[*]const T, x: [*]const T, f: Float(T)) void {
    const m = a.rows;
    @memset(z[0..m], ops.zero(T));
    for (0..m) |j| {
        const col = a.col(j);
        const uc = ops.mulReal(T, x[j], f);
        var acc: [4]T = @splat(ops.zero(T));
        const wj = if (w) |wv| ops.conj(T, wv[j]) else ops.zero(T);
        const uj = if (u) |uv| ops.conj(T, uv[j]) else ops.zero(T);
        for (j..m) |i| {
            var v = col[i];
            if (u) |uv| v = ops.sub(T, ops.sub(T, v, ops.mul(T, wj, uv[i])), ops.mul(T, uj, w.?[i]));
            col[i] = v;
            z[i] = ops.add(T, z[i], ops.mul(T, v, uc));
            acc[(i - j) % 4] = ops.add(T, acc[(i - j) % 4], ops.mul(T, ops.conj(T, v), x[i]));
        }
        const s = ops.add(T, ops.add(T, acc[0], acc[1]), ops.add(T, acc[2], acc[3]));
        y[j] = ops.mulReal(T, ops.sub(T, s, ops.mul(T, col[j], x[j])), f);
    }
}

// ---- tests -----------------------------------------------------------------------

const testing = std.testing;
const check = @import("../blas/_test/check.zig");
const lcheck = @import("_test/check.zig");
const types = @import("../blas/types.zig");

test "faer tridiag: Qᴴ A Q = T, real and complex (test_tridiag_real, test_tridiag_cplx)" {
    // faer's sizes 2, 3, 4, 8, 16 with block size 3, and larger
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64, f32, types.Complex32 }) |T| {
        const F = Float(T);
        for ([_]usize{ 1, 2, 3, 4, 8, 16, 45 }) |n| {
            const a0 = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, a0);
            lcheck.hermitian(T, a0, 3 * n);
            const v = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, v);
            v.copyFrom(a0);
            const b = 3;
            const h = try ops.alloc(T, al, b, @max(n, 1) - 1);
            defer ops.free(T, al, h);
            try tridiagonalize(al, check.io, T, v, h);
            // Qᴴ A Q on rows and columns 1..n
            const x = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, x);
            x.copyFrom(a0);
            if (n > 1) {
                const basis = v.sub(1, 0, n - 1, n - 1);
                try hh.applySequenceTransposed(al, check.io, T, basis, h, true, x.sub(1, 0, n - 1, n));
                try hh.applySequenceRight(al, check.io, T, basis, h, false, false, x.sub(0, 1, n, n - 1));
            }
            const tol = 64 * std.math.floatEps(F) * @as(F, @floatFromInt(n + 1)) * auxiliary.normMax(T, a0);
            for (0..n) |j| for (0..n) |i| {
                const want = if (i == j or i == j + 1) v.get(i, j) else if (j == i + 1) ops.conj(T, v.get(j, i)) else ops.zero(T);
                const err = ops.abs(T, ops.sub(T, x.get(i, j), want));
                if (!(err <= tol)) {
                    check.say("tridiag {s} n={d}: ({d},{d}) off by {e}\n", .{ @typeName(T), n, i, j, err });
                    return error.TestUnexpectedResult;
                }
            };
        }
    }
}
