//! Row reduction and matrix generators: `matlab88.f90`'s RREF, HILBER and
//! MAGIC, and `MATFN6`'s EYE, ONES, DIAG, TRIL, TRIU and KRON bodies.

const std = @import("std");
const scalar = @import("scalar.zig");
const CxMat = @import("CxMat.zig");
const blas = @import("../blas/__root__.zig");
const ops = @import("../lapack/ops.zig");
const Cx = scalar.Cx;
const ix = scalar.ix;
const level1 = scalar.level1;
const Allocator = std.mem.Allocator;
const Io = std.Io;

/// RREF: reduce `self` to reduced row-echelon form in place, treating any
/// entry no larger than `eps * 2*max(M,N) * max_column_norm(A)` as a
/// numerical zero. By BLAS 1 (`i?amax`, `?swap`, `?scal`, `?axpy`; M2)
/// over a packed copy, the real routines for a real `self`.
pub fn rref(self: CxMat, allocator: Allocator, io: Io, eps: f64) (Allocator.Error || Io.Cancelable)!void {
    switch (self.isReal()) {
        inline else => |real| {
            const T = if (real) f64 else Cx;
            const x = try self.toPacked(T, allocator);
            defer allocator.free(x);
            try rrefPacked(T, io, self.m, self.n, x, eps);
            self.fromPacked(T, x);
        },
    }
}

fn rrefPacked(comptime T: type, io: Io, m: usize, n: usize, x: []T, eps: f64) Io.Cancelable!void {
    const mi: isize = @intCast(m);
    var tol: f64 = 0.0;
    for (0..n) |j| tol = @max(tol, try level1(blas.asum(io, T, mi, x[j * m ..].ptr, 1)));
    tol *= eps * @as(f64, @floatFromInt(2 * @max(m, n)));

    // 0-based; a row from column `l` on is `w` entries `m` apart
    var k: usize = 0;
    var l: usize = 0;
    while (k < m and l < n) {
        const col = x[l * m ..][0..m];
        const amax = try level1(blas.index(io, T, .amax, @intCast(m - k), col[k..].ptr, 1));
        const piv = k + @as(usize, @intCast(amax)) - 1;
        if (ops.abs1(T, col[piv]) <= tol) {
            @memset(col[k..], ops.zero(T));
            l += 1;
            continue;
        }
        const w: isize = @intCast(n - l);
        try level1(blas.swap(io, T, w, x[l * m + piv ..].ptr, mi, x[l * m + k ..].ptr, mi));
        try level1(blas.scal(io, T, w, recip(col[k]), x[l * m + k ..].ptr, mi));
        col[k] = ops.one(T);
        for (0..m) |row| {
            if (row == k) continue;
            try level1(blas.axpy(io, T, w, ops.neg(T, col[row]), x[l * m + k ..].ptr, mi, x[l * m + row ..].ptr, mi));
        }
        k += 1;
        l += 1;
    }
}

fn recip(x: anytype) @TypeOf(x) {
    return if (@TypeOf(x) == f64) 1 / x else scalar.wdiv(Cx.init(1, 0), x);
}

/// HILBER: the exact inverse of the `n` by `n` Hilbert matrix into the real
/// `a` (leading dimension `lda`), by its closed form; inverting is hopeless.
pub fn hilber(a: []f64, lda: usize, n: usize) void {
    var p: f64 = @floatFromInt(n);
    var i: usize = 1;
    while (i <= n) : (i += 1) {
        if (i != 1) {
            const im1: f64 = @floatFromInt(i - 1);
            p = @as(f64, @floatFromInt(n - i + 1)) * p * @as(f64, @floatFromInt(n + i - 1)) / (im1 * im1);
        }
        var r = p * p;
        a[ix(lda, i, i)] = r / @as(f64, @floatFromInt(2 * i - 1));
        if (i == n) continue;

        var j: usize = i + 1;
        while (j <= n) : (j += 1) {
            const jm1: f64 = @floatFromInt(j - 1);
            r = -(@as(f64, @floatFromInt(n - j + 1)) * r * @as(f64, @floatFromInt(n + j - 1))) / (jm1 * jm1);
            const v = r / @as(f64, @floatFromInt(i + j - 1));
            a[ix(lda, i, j)] = v;
            a[ix(lda, j, i)] = v;
        }
    }
}

/// EYE: fill `self` with the (rectangular) identity; 1 where `i == j`, 0
/// elsewhere.
pub fn eye(self: CxMat) void {
    var j: usize = 1;
    while (j <= self.n) : (j += 1) {
        var i: usize = 1;
        while (i <= self.m) : (i += 1) {
            self.r[self.idx(i, j)] = if (i == j) 1.0 else 0.0;
            self.i[self.idx(i, j)] = 0.0;
        }
    }
}

/// ONES: fill `self` with all ones.
pub fn ones(self: CxMat) void {
    var j: usize = 1;
    while (j <= self.n) : (j += 1) {
        const c = self.col(j);
        @memset(c.r[0..self.m], 1.0);
        @memset(c.i[0..self.m], 0.0);
    }
}

/// MAGIC: an `n` by `n` magic square into the real `a` (leading dimension
/// `lda`): Siamese for odd `n`, LUX for singly even, a pattern for doubly
/// even. None is 2 by 2; the REPL gave an empty result.
pub fn magic(io: Io, a: []f64, lda: usize, n: usize) Io.Cancelable!void {
    if (n % 4 == 0) {
        var k: usize = 1;
        var i: usize = 1;
        while (i <= n) : (i += 1) {
            var j: usize = 1;
            while (j <= n) : (j += 1) {
                var v: f64 = @floatFromInt(k);
                if ((i % 4) / 2 == (j % 4) / 2) v = @as(f64, @floatFromInt(n * n + 1)) - v;
                a[ix(lda, i, j)] = v;
                k += 1;
            }
        }
        return;
    }

    // Odd order, or the upper-left quadrant of a singly-even order.
    const m: usize = if (n % 2 == 0) n / 2 else n;
    var jz: usize = 1;
    while (jz <= m) : (jz += 1) {
        var iz: usize = 1;
        while (iz <= m) : (iz += 1) a[ix(lda, iz, jz)] = 0.0;
    }
    var i: usize = 1;
    var j: usize = (m + 1) / 2;
    const mm = m * m;
    var k: usize = 1;
    while (k <= mm) : (k += 1) {
        a[ix(lda, i, j)] = @floatFromInt(k);
        var ni = if (i == 1) m else i - 1;
        var nj = if (j == m) 1 else j + 1;
        if (a[ix(lda, ni, nj)] != 0.0) {
            ni = if (i == m) 1 else i + 1;
            nj = j;
        }
        i = ni;
        j = nj;
    }
    if (n % 2 != 0) return;

    // Rest of a singly-even order (the LUX method's quadrant fill-in and
    // column swaps).
    const t = @as(f64, @floatFromInt(m * m));
    var ie: usize = 1;
    while (ie <= m) : (ie += 1) {
        var je: usize = 1;
        while (je <= m) : (je += 1) {
            const im = ie + m;
            const jm = je + m;
            const base = a[ix(lda, ie, je)];
            a[ix(lda, ie, jm)] = base + 2.0 * t;
            a[ix(lda, im, je)] = base + 3.0 * t;
            a[ix(lda, im, jm)] = base + t;
        }
    }
    const m1 = (m - 1) / 2;
    if (m1 == 0) return;
    var jc: usize = 1;
    while (jc <= m1) : (jc += 1) {
        try swap(io, m, a[ix(lda, 1, jc)..], a[ix(lda, m + 1, jc)..]);
    }
    const m1b = (m + 1) / 2;
    const m2b = m1b + m;
    try swap(io, 1, a[ix(lda, m1b, 1)..], a[ix(lda, m2b, 1)..]);
    try swap(io, 1, a[ix(lda, m1b, m1b)..], a[ix(lda, m2b, m1b)..]);
    const m1c = n + 1 - (m - 3) / 2;
    if (m1c > n) return;
    var jd: usize = m1c;
    while (jd <= n) : (jd += 1) {
        try swap(io, m, a[ix(lda, 1, jd)..], a[ix(lda, m + 1, jd)..]);
    }
}

/// `dswap` of `n` contiguous entries.
fn swap(io: Io, n: usize, x: []f64, y: []f64) Io.Cancelable!void {
    try level1(blas.swap(io, f64, @intCast(n), x.ptr, 1, y.ptr, 1));
}

/// DIAG: `v` (length `n`) on the `k`-th diagonal (above for `k > 0`) of
/// `self`, square of order `n + |k|`, zero elsewhere.
pub fn diagFromVector(self: CxMat, vr: []const f64, vi: []const f64, n: usize, k: isize) void {
    const order = n + @as(usize, @intCast(@abs(k)));
    @memset(self.r[0 .. order * order], 0.0);
    @memset(self.i[0 .. order * order], 0.0);
    var idx: usize = 1;
    while (idx <= n) : (idx += 1) {
        const row: usize = if (k < 0) idx + @as(usize, @intCast(-k)) else idx;
        const col: usize = if (k < 0) idx else idx + @as(usize, @intCast(k));
        self.r[self.idx(row, col)] = vr[idx - 1];
        self.i[self.idx(row, col)] = vi[idx - 1];
    }
}

/// DIAG: extract the `k`-th diagonal of `self` into the length-
/// `min(m, n-k)` (`k>=0`) or `min(m+k, n)` (`k<0`) vector `dr`/`di`.
/// Returns how many entries were written.
pub fn diagToVector(self: CxMat, k: isize, dr: []f64, di: []f64) usize {
    const m = self.m;
    const n = self.n;
    const len: usize = if (k >= 0)
        @min(m, n -| @as(usize, @intCast(k)))
    else
        @min(m -| @as(usize, @intCast(-k)), n);
    var idx: usize = 1;
    while (idx <= len) : (idx += 1) {
        const row: usize = if (k < 0) idx + @as(usize, @intCast(-k)) else idx;
        const col: usize = if (k < 0) idx else idx + @as(usize, @intCast(k));
        dr[idx - 1] = self.r[self.idx(row, col)];
        di[idx - 1] = self.i[self.idx(row, col)];
    }
    return len;
}

/// TRIL: zero the strict upper triangle above the `k`-th diagonal of
/// `self` in place, leaving its lower triangular part (`k == 0`: including
/// the main diagonal; `k < 0`: below that).
pub fn tril(self: CxMat, k: isize) void {
    const m = self.m;
    const n = self.n;
    var j: usize = 1;
    while (j <= n) : (j += 1) {
        const jk = @as(isize, @intCast(j)) - k;
        if (jk <= 0) continue;
        const count: usize = @intCast(jk - 1);
        if (count > 0) {
            const c = self.col(j);
            @memset(c.r[0..@min(count, m)], 0.0);
            @memset(c.i[0..@min(count, m)], 0.0);
        }
    }
}

/// TRIU: zero the strict lower triangle below the `k`-th diagonal of
/// `self` in place, leaving its upper triangular part.
pub fn triu(self: CxMat, k: isize) void {
    const m = self.m;
    const n = self.n;
    var j: usize = 1;
    while (j <= n) : (j += 1) {
        const jk = @as(isize, @intCast(j)) - k;
        const first_zero: usize = if (jk < 1) 1 else @as(usize, @intCast(jk)) + 1;
        if (first_zero > m) continue;
        @memset(self.r[self.idx(first_zero, j)..][0 .. m - first_zero + 1], 0.0);
        @memset(self.i[self.idx(first_zero, j)..][0 .. m - first_zero + 1], 0.0);
    }
}

pub const KronOp = enum { mul, div, rdiv };

/// KRON: the Kronecker product of `a` and `b` into `out`; `.div` and `.rdiv`
/// divide each block instead (`A(ia,ja) / B`, `B / A(ia,ja)`).
pub fn kron(a: CxMat, b: CxMat, out: CxMat, op: KronOp) void {
    const ma = a.m;
    const na = a.n;
    const m = b.m;
    const n = b.n;
    const lda_result = out.lda;

    var ja: usize = 1;
    while (ja <= na) : (ja += 1) {
        var j: usize = 1;
        while (j <= n) : (j += 1) {
            var ia: usize = 1;
            while (ia <= ma) : (ia += 1) {
                const a_val = a.get(ia, ja);
                var i: usize = 1;
                while (i <= m) : (i += 1) {
                    const b_val = b.get(i, j);
                    const v = switch (op) {
                        .mul => a_val.mul(b_val),
                        .div => scalar.wdiv(a_val, b_val),
                        .rdiv => scalar.wdiv(b_val, a_val),
                    };
                    const row = (ia - 1) * m + i;
                    const col = (ja - 1) * n + j;
                    out.r[ix(lda_result, row, col)] = v.re;
                    out.i[ix(lda_result, row, col)] = v.im;
                }
            }
        }
    }
}

const zigstd = @import("../../zigstd.zig");
const testing = std.testing;

test "hilber matches the known inverse of the 3x3 Hilbert matrix" {
    // inv(hilb(3)) = [[9,-36,30],[-36,192,-180],[30,-180,180]].
    var a: [9]f64 = @splat(0.0);
    hilber(&a, 3, 3);
    try zigstd.testExpectApproxEqAbs(9.0, a[ix(3, 1, 1)], 1e-6);
    try zigstd.testExpectApproxEqAbs(-36.0, a[ix(3, 1, 2)], 1e-6);
    try zigstd.testExpectApproxEqAbs(30.0, a[ix(3, 1, 3)], 1e-6);
    try zigstd.testExpectApproxEqAbs(192.0, a[ix(3, 2, 2)], 1e-6);
    try zigstd.testExpectApproxEqAbs(-180.0, a[ix(3, 2, 3)], 1e-6);
    try zigstd.testExpectApproxEqAbs(180.0, a[ix(3, 3, 3)], 1e-6);
}

test "eye and ones produce the expected 2x3 matrices" {
    var ar: [6]f64 = @splat(0.0);
    var ai: [6]f64 = @splat(0.0);
    const mat = CxMat.init(&ar, &ai, 2, 2, 3);
    eye(mat);
    try zigstd.testExpectApproxEqAbs(1.0, ar[ix(2, 1, 1)], 1e-12);
    try zigstd.testExpectApproxEqAbs(0.0, ar[ix(2, 2, 1)], 1e-12);
    try zigstd.testExpectApproxEqAbs(1.0, ar[ix(2, 2, 2)], 1e-12);
    try zigstd.testExpectApproxEqAbs(0.0, ar[ix(2, 1, 3)], 1e-12);

    ones(mat);
    for (ar) |v| try zigstd.testExpectApproxEqAbs(1.0, v, 1e-12);
}

fn checkMagicSquare(a: []const f64, n: usize) !void {
    const target = @as(f64, @floatFromInt(n)) * (@as(f64, @floatFromInt(n * n)) + 1.0) / 2.0;
    var i: usize = 1;
    while (i <= n) : (i += 1) {
        var row_sum: f64 = 0.0;
        var col_sum: f64 = 0.0;
        var j: usize = 1;
        while (j <= n) : (j += 1) {
            row_sum += a[ix(n, i, j)];
            col_sum += a[ix(n, j, i)];
        }
        try zigstd.testExpectApproxEqAbs(target, row_sum, 1e-9);
        try zigstd.testExpectApproxEqAbs(target, col_sum, 1e-9);
    }
    var diag1: f64 = 0.0;
    var diag2: f64 = 0.0;
    var k: usize = 1;
    while (k <= n) : (k += 1) {
        diag1 += a[ix(n, k, k)];
        diag2 += a[ix(n, k, n + 1 - k)];
    }
    try zigstd.testExpectApproxEqAbs(target, diag1, 1e-9);
    try zigstd.testExpectApproxEqAbs(target, diag2, 1e-9);
}

test "magic produces valid squares for odd, singly-even and doubly-even orders" {
    var a5: [25]f64 = @splat(0.0);
    try magic(std.testing.io, &a5, 5, 5);
    try checkMagicSquare(&a5, 5);

    var a6: [36]f64 = @splat(0.0);
    try magic(std.testing.io, &a6, 6, 6);
    try checkMagicSquare(&a6, 6);

    var a8: [64]f64 = @splat(0.0);
    try magic(std.testing.io, &a8, 8, 8);
    try checkMagicSquare(&a8, 8);
}

test "diagFromVector and diagToVector round-trip" {
    const vr = [_]f64{ 1.0, 2.0, 3.0 };
    const vi = [_]f64{ 0.0, 0.0, 0.0 };
    var ar: [16]f64 = @splat(0.0); // room for order 4 (n=3, k=1)
    var ai: [16]f64 = @splat(0.0);
    diagFromVector(CxMat.initPacked(&ar, &ai, 4, 4), &vr, &vi, 3, 1);
    // order = 3+1 = 4; diagonal k=1 means A(i, i+1) = v(i).
    try zigstd.testExpectApproxEqAbs(1.0, ar[ix(4, 1, 2)], 1e-12);
    try zigstd.testExpectApproxEqAbs(2.0, ar[ix(4, 2, 3)], 1e-12);
    try zigstd.testExpectApproxEqAbs(3.0, ar[ix(4, 3, 4)], 1e-12);

    var dr: [3]f64 = @splat(0.0);
    var di: [3]f64 = @splat(0.0);
    const len = diagToVector(CxMat.initPacked(&ar, &ai, 4, 4), 1, &dr, &di);
    try zigstd.testExpectEqual(@as(usize, 3), len);
    try zigstd.testExpectApproxEqAbs(1.0, dr[0], 1e-12);
    try zigstd.testExpectApproxEqAbs(2.0, dr[1], 1e-12);
    try zigstd.testExpectApproxEqAbs(3.0, dr[2], 1e-12);
}

test "tril and triu zero the expected triangle" {
    var ar = [_]f64{ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 }; // 3x3, col-major
    var ai: [9]f64 = @splat(0.0);
    var ar2 = ar;
    var ai2 = ai;

    tril(CxMat.initPacked(&ar, &ai, 3, 3), 0);
    try zigstd.testExpectApproxEqAbs(0.0, ar[ix(3, 1, 2)], 1e-12); // above diagonal zeroed
    try zigstd.testExpectApproxEqAbs(2.0, ar[ix(3, 2, 1)], 1e-12); // below diagonal kept

    triu(CxMat.initPacked(&ar2, &ai2, 3, 3), 0);
    try zigstd.testExpectApproxEqAbs(0.0, ar2[ix(3, 2, 1)], 1e-12); // below diagonal zeroed
    try zigstd.testExpectApproxEqAbs(4.0, ar2[ix(3, 1, 2)], 1e-12); // above diagonal kept
}

test "kron matches the textbook 2x2 example" {
    // A = [[1,2],[3,4]], B = [[0,5],[6,7]].
    var ar_a = [_]f64{ 1.0, 3.0, 2.0, 4.0 };
    var ai_a: [4]f64 = @splat(0.0);
    var ar_b = [_]f64{ 0.0, 6.0, 5.0, 7.0 };
    var ai_b: [4]f64 = @splat(0.0);
    var ar: [16]f64 = @splat(0.0);
    var ai: [16]f64 = @splat(0.0);
    kron(
        CxMat.initPacked(&ar_a, &ai_a, 2, 2),
        CxMat.initPacked(&ar_b, &ai_b, 2, 2),
        CxMat.initPacked(&ar, &ai, 4, 4),
        .mul,
    );

    // Expected kron(A,B) (row-major for readability):
    // [ 0  5 |  0 10]
    // [ 6  7 | 12 14]
    // [ 0 15 |  0 20]
    // [18 21 | 24 28]
    try zigstd.testExpectApproxEqAbs(0.0, ar[ix(4, 1, 1)], 1e-12);
    try zigstd.testExpectApproxEqAbs(5.0, ar[ix(4, 1, 2)], 1e-12);
    try zigstd.testExpectApproxEqAbs(10.0, ar[ix(4, 1, 4)], 1e-12);
    try zigstd.testExpectApproxEqAbs(18.0, ar[ix(4, 4, 1)], 1e-12);
    try zigstd.testExpectApproxEqAbs(28.0, ar[ix(4, 4, 4)], 1e-12);
}

test "rref reduces a simple rank-deficient system" {
    // [[1, 2, 3], [2, 4, 6]] -> rank 1, row-reduces to [[1,2,3],[0,0,0]].
    var ar = [_]f64{ 1.0, 2.0, 2.0, 4.0, 3.0, 6.0 }; // column-major, 2x3
    var ai: [6]f64 = @splat(0.0);
    try rref(CxMat.initPacked(&ar, &ai, 2, 3), testing.allocator, testing.io, 1e-12);
    try zigstd.testExpectApproxEqAbs(1.0, ar[ix(2, 1, 1)], 1e-9);
    try zigstd.testExpectApproxEqAbs(2.0, ar[ix(2, 1, 2)], 1e-9);
    try zigstd.testExpectApproxEqAbs(3.0, ar[ix(2, 1, 3)], 1e-9);
    try zigstd.testExpectApproxEqAbs(0.0, ar[ix(2, 2, 1)], 1e-9);
    try zigstd.testExpectApproxEqAbs(0.0, ar[ix(2, 2, 2)], 1e-9);
    try zigstd.testExpectApproxEqAbs(0.0, ar[ix(2, 2, 3)], 1e-9);
}

test "rref of a complex rank-deficient matrix" {
    // [[1, i], [i, -1]]: the second row is i times the first -> [[1, i], [0, 0]]
    var ar = [_]f64{ 1, 0, 0, -1 };
    var ai = [_]f64{ 0, 1, 1, 0 };
    try rref(CxMat.initPacked(&ar, &ai, 2, 2), testing.allocator, testing.io, 1e-12);
    try zigstd.testExpectEqualSlices(f64, &.{ 1, 0, 0, 0 }, &ar);
    try zigstd.testExpectEqualSlices(f64, &.{ 0, 0, 1, 0 }, &ai);
}
