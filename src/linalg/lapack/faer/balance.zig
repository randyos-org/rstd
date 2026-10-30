//! Balancing for `geev`, from Reference-LAPACK (BSD-3, as OpenBLAS
//! v0.3.34-317 ships it in `lapack-netlib/SRC/`): `dgebal.f` and `zgebal.f`
//! (job `B`: permute to isolate eigenvalues, then scale rows and columns by
//! powers of 2 toward equal norms) and `dgebak.f`, `zgebak.f` (the vectors
//! back). Faer has no balancing; `dgeev.f` balances before the Hessenberg
//! reduction, and so does `geev` here.

const std = @import("std");
const ops = @import("ops.zig");

const Io = std.Io;
const Mat = ops.Mat;
const Float = ops.Float;
const isComplex = ops.isComplex;

/// What `balance` found: `a`'s rows and columns outside `ilo..ihi` are
/// already triangular; `perm[i]` is the row swapped with `i` there, `scale`
/// the factors inside.
pub const Info = struct { ilo: usize, ihi: usize };

fn swapRows(comptime T: type, a: Mat(T), i: usize, j: usize, from: usize, to: usize) void {
    for (from..to) |c| std.mem.swap(T, a.at(i, c), a.at(j, c));
}

fn swapCols(comptime T: type, a: Mat(T), i: usize, j: usize, rows: usize) void {
    for (0..rows) |r| std.mem.swap(T, a.at(r, i), a.at(r, j));
}

/// `IDAMAX`, `IZAMAX` (by `|re| + |im|`), then the modulus of that entry.
fn maxAbs(comptime T: type, x: [*]const T, n: usize, inc: usize) Float(T) {
    var best: usize = 0;
    var m: Float(T) = -1;
    for (0..n) |i| {
        const v = ops.abs1(T, x[i * inc]);
        if (v > m) {
            m = v;
            best = i;
        }
    }
    return ops.abs(T, x[best * inc]);
}

/// `xGEBAL` with job `B` on the `n x n` `a`, which holds no NaN.
pub fn balance(io: Io, comptime T: type, a: Mat(T), scale: []Float(T), perm: []usize) Info {
    const F = Float(T);
    const n = a.rows;
    @memset(scale, 1);
    for (perm, 0..) |*p, i| p.* = i;
    var k: usize = 0;
    var l: usize = n;
    // rows with no off-diagonal entry in `0..l` go to the bottom
    var noconv = true;
    while (noconv) {
        noconv = false;
        var i = l;
        while (i > 0) {
            i -= 1;
            var canswap = true;
            for (0..l) |j| if (i != j and !ops.isZero(T, a.get(i, j))) {
                canswap = false;
                break;
            };
            if (!canswap) continue;
            perm[l - 1] = i;
            if (i != l - 1) {
                swapCols(T, a, i, l - 1, l);
                swapRows(T, a, i, l - 1, k, n);
            }
            noconv = true;
            if (l == 1) return .{ .ilo = 0, .ihi = 1 };
            l -= 1;
        }
    }
    // columns with no off-diagonal entry in `k..l` go to the top
    noconv = true;
    while (noconv) {
        noconv = false;
        const j1 = l;
        var j = k;
        while (j < j1) : (j += 1) {
            var canswap = true;
            for (k..l) |i| if (i != j and !ops.isZero(T, a.get(i, j))) {
                canswap = false;
                break;
            };
            if (!canswap) continue;
            perm[k] = j;
            if (j != k) {
                swapCols(T, a, j, k, l);
                swapRows(T, a, j, k, k, n);
            }
            noconv = true;
            k += 1;
        }
    }
    // `SCLFAC` 2, `FACTOR` 0.95
    const sfmin1 = std.math.floatMin(F) / std.math.floatEps(F);
    const sfmax1 = 1 / sfmin1;
    const sfmin2 = sfmin1 * 2;
    const sfmax2 = 1 / sfmin2;
    noconv = true;
    while (noconv) {
        noconv = false;
        for (k..l) |i| {
            var c = ops.nrm2(io, T, l - k, a.col(i) + k, 1);
            var r = ops.nrm2(io, T, l - k, a.ptr + i + k * a.ld, a.ld);
            var ca = maxAbs(T, a.col(i), l, 1);
            var ra = maxAbs(T, a.ptr + i + k * a.ld, n - k, a.ld);
            if (c == 0 or r == 0) continue;
            if (std.math.isNan(c + ca + r + ra)) return .{ .ilo = k, .ihi = l };
            var g = r / 2;
            var f: F = 1;
            const s = c + r;
            while (c < g and @max(f, @max(c, ca)) < sfmax2 and @min(r, @min(g, ra)) > sfmin2) {
                f *= 2;
                c *= 2;
                ca *= 2;
                r /= 2;
                g /= 2;
                ra /= 2;
            }
            g = c / 2;
            while (g >= r and @max(r, ra) < sfmax2 and @min(@min(f, c), @min(g, ca)) > sfmin2) {
                f /= 2;
                c /= 2;
                g /= 2;
                ca /= 2;
                r *= 2;
                ra *= 2;
            }
            if (c + r >= 0.95 * s) continue;
            if (f < 1 and scale[i] < 1 and f * scale[i] <= sfmin1) continue;
            if (f > 1 and scale[i] > 1 and scale[i] >= sfmax1 / f) continue;
            scale[i] *= f;
            noconv = true;
            for (k..n) |cc| a.at(i, cc).* = ops.mulReal(T, a.get(i, cc), 1 / f);
            for (0..l) |rr| a.at(rr, i).* = ops.mulReal(T, a.get(rr, i), f);
        }
    }
    return .{ .ilo = k, .ihi = l };
}

/// `xGEBAK` with job `B`: `v`'s rows (the balanced matrix's right vectors,
/// or with `left` its left ones) back to `a`'s.
pub fn back(comptime T: type, info: Info, scale: []const Float(T), perm: []const usize, left: bool, v: Mat(T)) void {
    const n = v.rows;
    for (info.ilo..info.ihi) |i| {
        const s = if (left) 1 / scale[i] else scale[i];
        for (0..v.cols) |c| v.at(i, c).* = ops.mulReal(T, v.get(i, c), s);
    }
    for (0..n) |ii| {
        if (ii >= info.ilo and ii < info.ihi) continue;
        const i = if (ii < info.ilo) info.ilo - 1 - ii else ii;
        const k = perm[i];
        if (k != i) swapRows(T, v, i, k, 0, v.cols);
    }
}

// ---- tests -----------------------------------------------------------------------

const testing = std.testing;
const check = @import("../../blas/_test/blas1/check.zig");
const types = @import("../../blas/types.zig");

test "faer balance: D⁻¹ Pᵀ A P D is the balanced matrix, isolated eigenvalues outside ilo..ihi" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex32 }) |T| {
        const F = Float(T);
        const n = 9;
        const a0 = try ops.alloc(T, al, n, n);
        defer ops.free(T, al, a0);
        check.fill(T, a0.ptr[0 .. n * n], 4);
        // graded rows and columns, and two isolated eigenvalues: row 2 and
        // column 6 zero off the diagonal
        for (0..n) |j| for (0..n) |i| {
            const d = std.math.pow(F, 2, @as(F, @floatFromInt(i)) * 3 - @as(F, @floatFromInt(j)) * 3);
            a0.at(i, j).* = ops.mulReal(T, a0.get(i, j), d);
            if ((i == 2 and j != 2) or (j == 6 and i != 6)) a0.at(i, j).* = ops.zero(T);
        };
        const b = try ops.alloc(T, al, n, n);
        defer ops.free(T, al, b);
        b.copyFrom(a0);
        var scale: [n]F = undefined;
        var perm: [n]usize = undefined;
        const info = balance(check.io, T, b, &scale, &perm);
        try testing.expect(info.ilo >= 1 and info.ihi <= n - 1);
        // X = P D, the identity backed: then A X = X B
        const x = try ops.alloc(T, al, n, n);
        defer ops.free(T, al, x);
        x.identity();
        back(T, info, &scale, &perm, false, x);
        const ax = try ops.alloc(T, al, n, n);
        defer ops.free(T, al, ax);
        const xb = try ops.alloc(T, al, n, n);
        defer ops.free(T, al, xb);
        try ops.gemm(al, check.io, T, .n, .n, ops.one(T), a0, x, ops.zero(T), ax);
        try ops.gemm(al, check.io, T, .n, .n, ops.one(T), x, b, ops.zero(T), xb);
        for (0..n) |j| for (0..n) |i| {
            const err = ops.abs(T, ops.sub(T, ax.get(i, j), xb.get(i, j)));
            try testing.expect(err <= 16 * std.math.floatEps(F) * @max(ops.abs(T, ax.get(i, j)), std.math.floatMin(F)));
        };
        // outside the window: triangular
        for (0..n) |j| for (0..n) |i| {
            if (i > j and (j < info.ilo or i >= info.ihi)) try testing.expect(ops.isZero(T, b.get(i, j)));
        };
    }
}
