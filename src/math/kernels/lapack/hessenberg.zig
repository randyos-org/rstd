//! Hessenberg reduction, from faer v0.24.4
//! (<https://codeberg.org/sarah-quinones/faer>, tag `faer-v0.24.4`, MIT):
//! `faer/src/linalg/evd/hessenberg.rs`'s `hessenberg_rearranged_unblocked`,
//! its fused pass serial (faer's `Par::Seq` path), and from order 256 the
//! blocked form (`hessenberg_blocked.zig`).

const std = @import("std");
const ops = @import("ops.zig");
const blas_calls = @import("blas_calls.zig");
const auxiliary = @import("auxiliary.zig");
const hh = @import("householder.zig");
const qr = @import("qr.zig");
const util = @import("../blas/util.zig");
const interface = @import("interface.zig");
const blocked_form = @import("hessenberg_blocked.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = ops.Mat;
const Float = ops.Float;
const Error = ops.Error;

/// `HessenbergParams`' `blocking_threshold`: the `n²` from which faer
/// reduces by blocks.
pub const blocking_threshold = 256 * 256;

/// `hessenberg_in_place`: `a` is `n x n`, `h` `bs x (n - 1)` (`bs > 0`).
pub fn hessenberg(al: Allocator, io: Io, comptime T: type, a: Mat(T), h: Mat(T)) Error!void {
    return hessenbergWith(al, io, T, a, h, blocking_threshold);
}

/// `hessenberg`, blocked from `n² >= threshold` (`h`'s rows the block
/// size), so tests reach the blocked form at small orders.
pub fn hessenbergWith(al: Allocator, io: Io, comptime T: type, a: Mat(T), h: Mat(T), threshold: usize) Error!void {
    const n = a.rows;
    if (n == 0) return;
    if (n * n >= threshold) return blocked_form.blocked(al, io, T, a, h);
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
            try blas_calls.gemv(al, io, T, .n, ops.one(T), a22, x2, 1, ops.zero(T), z2, 1);
            try blas_calls.gemv(al, io, T, .c, ops.one(T), a22, x2, 1, ops.zero(T), y2, 1);
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
            try blas_calls.gemv(al, io, T, .n, ops.one(T), a02, x2, 1, ops.zero(T), w.ptr, 1);
            for (0..rem) |j| {
                const c = ops.mulReal(T, ops.conj(T, x2[j]), -tau_inv);
                for (0..k) |i| a02.at(i, j).* = ops.add(T, a02.get(i, j), ops.mul(T, w[i], c));
            }
        }
        a21[0] = beta;
    }
    if (n > 1) try hh.sequenceFactors(al, io, T, a.sub(1, 0, n - 1, n - 1), h);
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

// ---- gehrd, orghr ------------------------------------------------------------------

/// What `gehrd` and `orghr` return: `BlockOutOfRange` unless
/// `lo <= hi <= n`, a slice too short (`a`, or `tau` under `n - 1`) or a
/// leading dimension too small, and the allocator's and cancellation.
pub const GehrdError = interface.BlockError || interface.ShapeError || Allocator.Error || Io.Cancelable;
pub const OrghrError = GehrdError;

/// The checks `gehrd` and `orghr` share.
fn block(n: usize, lo: usize, hi: usize, a_len: usize, lda: usize, tau_len: usize) GehrdError!void {
    try interface.blockIn(lo, hi, n);
    try util.holds(a_len, n, n, lda);
    if (tau_len < n -| 1) return error.SliceTooShort;
}

/// `?GEHRD` (`dgehrd.f`, `zgehrd.f`): `A := Qᴴ A Q`, upper Hessenberg, for
/// the `n x n` `A`, by faer's reduction (above) of rows and columns
/// `lo..hi` (0-based, `hi` past the end; outside them `A` is upper
/// triangular already, as `gebal` leaves it), `Q` then applied to the rest
/// of those rows and columns. `Q = H(lo) ... H(hi - 2)`, each
/// `H = I - tau v vᴴ` (the reference's `tau`, faer's reciprocal) with `v`
/// below the subdiagonal and its `tau` in `tau` (`n - 1`, 0 outside
/// `lo..hi - 1`).
pub fn gehrd(allocator: Allocator, io: Io, comptime T: type, n: usize, lo: usize, hi: usize, a: []T, lda: usize, tau: []T) GehrdError!void {
    comptime ops.requireFloat(T, "gehrd");
    try block(n, lo, hi, a.len, lda, tau.len);
    if (n == 0) return;
    const un = n;
    @memset(tau[0 .. un - 1], ops.zero(T));
    const nh = hi - lo;
    if (nh <= 1) return;
    const am: Mat(T) = .{ .ptr = a.ptr, .rows = un, .cols = un, .ld = lda };
    const bs = qr.blockSize(nh - 1, nh - 1);
    const h = try ops.alloc(T, allocator, bs, nh - 1);
    defer ops.free(T, allocator, h);
    const blk = am.sub(lo, lo, nh, nh);
    try hessenberg(allocator, io, T, blk, h);
    for (0..nh - 1) |j| tau[lo + j] = ops.scalar(T, 1 / ops.re(T, h.get(j % bs, j)));
    const basis = blk.sub(1, 0, nh - 1, nh - 1);
    // the rows above the block take `Q` from the right, the columns past it
    // `Qᴴ` from the left
    if (lo > 0) try hh.applySequenceRight(allocator, io, T, basis, h, false, false, am.sub(0, lo + 1, lo, nh - 1));
    if (hi < un) try hh.applySequenceTransposed(allocator, io, T, basis, h, true, am.sub(lo + 1, hi, nh - 1, un - hi));
}

/// `?ORGHR`, `?UNGHR` (`dorghr.f`, `zunghr.f`): `gehrd`'s `Q` (`n x n`) in
/// place of its reflections, from them and their `tau`s: shifted one column
/// right, the identity outside `lo + 1..hi`, then `orgqr`'s generation.
pub fn orghr(allocator: Allocator, io: Io, comptime T: type, n: usize, lo: usize, hi: usize, a: []T, lda: usize, tau: []const T) OrghrError!void {
    comptime ops.requireFloat(T, "orghr");
    try block(n, lo, hi, a.len, lda, tau.len);
    if (n == 0) return;
    const un = n;
    const am: Mat(T) = .{ .ptr = a.ptr, .rows = un, .cols = un, .ld = lda };
    var j = hi;
    while (j > lo + 1) {
        j -= 1;
        const c = am.col(j);
        @memset(c[0..j], ops.zero(T));
        for (j + 1..hi) |i| c[i] = am.get(i, j - 1);
        @memset(c[hi..un], ops.zero(T));
    }
    for (0..un) |k| if (k <= lo or k >= hi) {
        const c = am.col(k);
        @memset(c[0..un], ops.zero(T));
        c[k] = ops.one(T);
    };
    if (hi > lo + 1) {
        const nq = hi - lo - 1;
        try qr.generate(allocator, io, T, am.sub(lo + 1, lo + 1, nq, nq), nq, tau.ptr + lo);
    }
}

// ---- tests -----------------------------------------------------------------------

const testing = std.testing;
const check = @import("../blas/_test/check.zig");
const types = @import("../blas/types.zig");

test "faer hessenberg: Qᴴ A Q = H, real and complex, unblocked and blocked" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64, f32, types.Complex32 }) |T| {
        const F = Float(T);
        for ([_]usize{ blocking_threshold, 0 }) |threshold| for ([_]usize{ 1, 2, 3, 5, 8, 21, 60 }) |n| {
            const a0 = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, a0);
            check.fill(T, a0.ptr[0 .. n * n], 5 * n);
            const hm = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, hm);
            hm.copyFrom(a0);
            const b = 4;
            const h = try ops.alloc(T, al, b, @max(n, 1) - 1);
            defer ops.free(T, al, h);
            try hessenbergWith(al, check.io, T, hm, h, threshold);
            const x = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, x);
            x.copyFrom(a0);
            if (n > 1) {
                const basis = hm.sub(1, 0, n - 1, n - 1);
                try hh.applySequenceTransposed(al, check.io, T, basis, h, true, x.sub(1, 0, n - 1, n));
                try hh.applySequenceRight(al, check.io, T, basis, h, false, false, x.sub(0, 1, n, n - 1));
            }
            const tol = 64 * std.math.floatEps(F) * @as(F, @floatFromInt(n + 1)) * @max(auxiliary.normMax(T, a0), 1);
            for (0..n) |j| for (0..n) |i| {
                const want = if (i <= j + 1) hm.get(i, j) else ops.zero(T);
                const err = ops.abs(T, ops.sub(T, x.get(i, j), want));
                if (!(err <= tol)) {
                    check.say("hessenberg {s} n={d} threshold {d}: ({d},{d}) off by {e}\n", .{ @typeName(T), n, threshold, i, j, err });
                    return error.TestUnexpectedResult;
                }
            };
        };
    }
}

test "gehrd, orghr: Qᴴ A Q = H on the lo..hi block, Q unitary, real and complex" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64, f32, types.Complex32 }) |T| {
        const F = Float(T);
        // n, lo, hi; an empty block leaves Q = I
        for ([_][3]usize{ .{ 1, 0, 1 }, .{ 2, 0, 2 }, .{ 5, 0, 5 }, .{ 9, 2, 7 }, .{ 9, 0, 4 }, .{ 9, 5, 9 }, .{ 40, 4, 33 }, .{ 5, 2, 2 } }) |c| {
            const n = c[0];
            const lo = c[1];
            const hi = c[2];
            // upper triangular outside rows and columns lo..hi, as gebal leaves it
            const a0 = try al.alloc(T, n * n);
            defer al.free(a0);
            check.fill(T, a0, 3 * n + lo);
            for (0..n) |j| for (j + 1..n) |i| {
                if (j < lo or i >= hi) a0[i + j * n] = ops.zero(T);
            };
            const hm = try al.dupe(T, a0);
            defer al.free(hm);
            const tau = try al.alloc(T, @max(n, 2) - 1);
            defer al.free(tau);
            try gehrd(al, check.io, T, n, lo, hi, hm, n, tau);
            const q = try al.dupe(T, hm);
            defer al.free(q);
            try orghr(al, check.io, T, n, lo, hi, q, n, tau);
            const nrm = @max(auxiliary.normMax(T, .{ .ptr = a0.ptr, .rows = n, .cols = n, .ld = n }), 1);
            const tol = 64 * std.math.floatEps(F) * @as(F, @floatFromInt(n + 1)) * nrm;
            for (0..n) |j| for (0..n) |i| {
                // (Qᴴ A Q)(i, j) against H, zero below the subdiagonal
                var s = ops.zero(T);
                for (0..n) |k| for (0..n) |l| {
                    s = ops.add(T, s, ops.mul(T, ops.mul(T, ops.conj(T, q[k + i * n]), a0[k + l * n]), q[l + j * n]));
                };
                const want = if (i <= j + 1) hm[i + j * n] else ops.zero(T);
                const err = ops.abs(T, ops.sub(T, s, want));
                // Qᴴ Q = I
                var g = ops.zero(T);
                for (0..n) |k| g = ops.add(T, g, ops.mul(T, ops.conj(T, q[k + i * n]), q[k + j * n]));
                const gerr = ops.abs(T, ops.sub(T, g, if (i == j) ops.one(T) else ops.zero(T)));
                if (!(err <= tol) or !(gerr <= tol / nrm)) {
                    check.say("gehrd {s} n={d} lo={d} hi={d}: ({d},{d}) H off by {e}, QᴴQ by {e}\n", .{ @typeName(T), n, lo, hi, i, j, err, gerr });
                    return error.TestUnexpectedResult;
                }
            };
        }
    }
    var a: [4]f64 = @splat(0);
    var tau: [1]f64 = undefined;
    try testing.expectError(error.BlockOutOfRange, gehrd(al, check.io, f64, 2, 0, 3, &a, 2, &tau));
    try testing.expectError(error.BlockOutOfRange, gehrd(al, check.io, f64, 2, 2, 1, &a, 2, &tau));
    try testing.expectError(error.LeadingDimensionTooSmall, orghr(al, check.io, f64, 2, 0, 2, &a, 1, &tau));
    try testing.expectError(error.SliceTooShort, orghr(al, check.io, f64, 2, 0, 2, a[0..3], 2, &tau));
    try testing.expectError(error.SliceTooShort, gehrd(al, check.io, f64, 2, 0, 2, &a, 2, tau[0..0]));
    try gehrd(al, check.io, f64, 0, 0, 0, &a, 1, &tau);
}
