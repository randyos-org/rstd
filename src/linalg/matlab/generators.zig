//! Row reduction and matrix generation, ported from `matlab88.f90`'s RREF,
//! HILBER, MAGIC, and the EYE/ONES/DIAG/TRIL/TRIU/KRON command bodies in
//! `MATFN6`. See `blas1.zig` for the shared matrix/vector conventions this
//! port uses.

const blas1 = @import("blas1.zig");
const Cx = blas1.Cx;
const ix = blas1.ix;

/// RREF: reduce the M×N complex matrix `ar`/`ai` to reduced row-echelon
/// form in place, treating any entry no larger than
/// `eps * 2*max(M,N) * max_column_norm(A)` as a numerical zero.
pub fn rref(ar: []f64, ai: []f64, lda: usize, m: usize, n: usize, eps: f64) void {
    var tol: f64 = 0.0;
    var j: usize = 1;
    while (j <= n) : (j += 1) {
        tol = @max(tol, blas1.wasum(m, ar[ix(lda, 1, j)..], ai[ix(lda, 1, j)..], 1));
    }
    tol *= eps * @as(f64, @floatFromInt(2 * @max(m, n)));

    const lda_stride: isize = @intCast(lda);
    var k: usize = 1;
    var l: usize = 1;
    while (k <= m and l <= n) {
        const piv = blas1.iwamax(m - k + 1, ar[ix(lda, k, l)..], ai[ix(lda, k, l)..], 1) + k - 1;
        if (blas1.cabs1(ar[ix(lda, piv, l)], ai[ix(lda, piv, l)]) <= tol) {
            blas1.wset(m - k + 1, Cx.init(0.0, 0.0), ar[ix(lda, k, l)..], ai[ix(lda, k, l)..], 1);
            l += 1;
            continue;
        }
        blas1.wswap(n - l + 1, ar[ix(lda, piv, l)..], ai[ix(lda, piv, l)..], lda_stride, ar[ix(lda, k, l)..], ai[ix(lda, k, l)..], lda_stride);
        const inv = blas1.wdiv(Cx.init(1.0, 0.0), Cx.init(ar[ix(lda, k, l)], ai[ix(lda, k, l)]));
        blas1.wscal(n - l + 1, inv, ar[ix(lda, k, l)..], ai[ix(lda, k, l)..], lda_stride);
        ar[ix(lda, k, l)] = 1.0;
        ai[ix(lda, k, l)] = 0.0;

        var row: usize = 1;
        while (row <= m) : (row += 1) {
            if (row != k) {
                const t = Cx.init(-ar[ix(lda, row, l)], -ai[ix(lda, row, l)]);
                blas1.waxpy(n - l + 1, t, ar[ix(lda, k, l)..], ai[ix(lda, k, l)..], lda_stride, ar[ix(lda, row, l)..], ai[ix(lda, row, l)..], lda_stride);
            }
        }
        k += 1;
        l += 1;
    }
}

/// HILBER: generate the exact inverse of the N×N Hilbert matrix directly
/// into the real array `a` (column-major, leading dimension `lda`), using
/// the closed-form recurrence rather than computing it via LU inversion -
/// the Hilbert matrix is notoriously ill-conditioned, so an explicit
/// formula is far more accurate than inverting it numerically. `a` is
/// real-only; zero the imaginary part yourself if a complex result is
/// needed.
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

/// EYE: fill the `m`x`n` complex matrix `ar`/`ai` with the (rectangular)
/// identity - 1 where `i == j`, 0 elsewhere.
pub fn eye(ar: []f64, ai: []f64, lda: usize, m: usize, n: usize) void {
    var j: usize = 1;
    while (j <= n) : (j += 1) {
        var i: usize = 1;
        while (i <= m) : (i += 1) {
            ar[ix(lda, i, j)] = if (i == j) 1.0 else 0.0;
            ai[ix(lda, i, j)] = 0.0;
        }
    }
}

/// ONES: fill the `m`x`n` complex matrix `ar`/`ai` with all ones.
pub fn ones(ar: []f64, ai: []f64, lda: usize, m: usize, n: usize) void {
    var j: usize = 1;
    while (j <= n) : (j += 1) {
        blas1.wset(m, Cx.init(1.0, 0.0), ar[ix(lda, 1, j)..], ai[ix(lda, 1, j)..], 1);
    }
}

/// MAGIC: fill the `n`x`n` real array `a` (column-major, leading dimension
/// `lda`) with a magic square (every row, column, and both diagonals sum
/// to the same value), using the Siamese method for odd `n`, the LUX
/// method for singly-even `n`, and a direct pattern for doubly-even `n`.
/// There is no 2x2 magic square; callers should treat `n == 2` specially
/// (the original REPL substituted an empty result).
pub fn magic(a: []f64, lda: usize, n: usize) void {
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
        blas1.rswap(m, a[ix(lda, 1, jc)..], 1, a[ix(lda, m + 1, jc)..], 1);
    }
    const m1b = (m + 1) / 2;
    const m2b = m1b + m;
    blas1.rswap(1, a[ix(lda, m1b, 1)..], 1, a[ix(lda, m2b, 1)..], 1);
    blas1.rswap(1, a[ix(lda, m1b, m1b)..], 1, a[ix(lda, m2b, m1b)..], 1);
    const m1c = n + 1 - (m - 3) / 2;
    if (m1c > n) return;
    var jd: usize = m1c;
    while (jd <= n) : (jd += 1) {
        blas1.rswap(m, a[ix(lda, 1, jd)..], 1, a[ix(lda, m + 1, jd)..], 1);
    }
}

/// DIAG: if `v` (length `n`) is a vector, `diagFromVector` builds the
/// `n+|k|`-order square matrix with `v` on the `k`-th diagonal (`k == 0`:
/// main diagonal; `k > 0`: above it; `k < 0`: below it) into `ar`/`ai`
/// (must be sized for that larger order). Everything off that diagonal is
/// zeroed first.
pub fn diagFromVector(vr: []const f64, vi: []const f64, n: usize, k: isize, ar: []f64, ai: []f64) void {
    const order = n + @as(usize, @intCast(@abs(k)));
    blas1.wset(order * order, Cx.init(0.0, 0.0), ar, ai, 1);
    var idx: usize = 1;
    while (idx <= n) : (idx += 1) {
        const row: usize = if (k < 0) idx + @as(usize, @intCast(-k)) else idx;
        const col: usize = if (k < 0) idx else idx + @as(usize, @intCast(k));
        ar[ix(order, row, col)] = vr[idx - 1];
        ai[ix(order, row, col)] = vi[idx - 1];
    }
}

/// DIAG: extract the `k`-th diagonal of the `m`x`n` matrix `ar`/`ai` into
/// the length-`min(m, n-k)` (`k>=0`) or `min(m+k, n)` (`k<0`) vector
/// `dr`/`di`. Returns how many entries were written.
pub fn diagToVector(ar: []const f64, ai: []const f64, lda: usize, m: usize, n: usize, k: isize, dr: []f64, di: []f64) usize {
    const len: usize = if (k >= 0)
        @min(m, n -| @as(usize, @intCast(k)))
    else
        @min(m -| @as(usize, @intCast(-k)), n);
    var idx: usize = 1;
    while (idx <= len) : (idx += 1) {
        const row: usize = if (k < 0) idx + @as(usize, @intCast(-k)) else idx;
        const col: usize = if (k < 0) idx else idx + @as(usize, @intCast(k));
        dr[idx - 1] = ar[ix(lda, row, col)];
        di[idx - 1] = ai[ix(lda, row, col)];
    }
    return len;
}

/// TRIL: zero the strict upper triangle above the `k`-th diagonal of the
/// `m`x`n` matrix `ar`/`ai` in place, leaving its lower triangular part
/// (`k == 0`: including the main diagonal; `k < 0`: below that).
pub fn tril(ar: []f64, ai: []f64, lda: usize, m: usize, n: usize, k: isize) void {
    var j: usize = 1;
    while (j <= n) : (j += 1) {
        const jk = @as(isize, @intCast(j)) - k;
        if (jk <= 0) continue;
        const count: usize = @intCast(jk - 1);
        if (count > 0) blas1.wset(@min(count, m), Cx.init(0.0, 0.0), ar[ix(lda, 1, j)..], ai[ix(lda, 1, j)..], 1);
    }
}

/// TRIU: zero the strict lower triangle below the `k`-th diagonal of the
/// `m`x`n` matrix `ar`/`ai` in place, leaving its upper triangular part.
pub fn triu(ar: []f64, ai: []f64, lda: usize, m: usize, n: usize, k: isize) void {
    var j: usize = 1;
    while (j <= n) : (j += 1) {
        const jk = @as(isize, @intCast(j)) - k;
        const first_zero: usize = if (jk < 1) 1 else @as(usize, @intCast(jk)) + 1;
        if (first_zero > m) continue;
        blas1.wset(m - first_zero + 1, Cx.init(0.0, 0.0), ar[ix(lda, first_zero, j)..], ai[ix(lda, first_zero, j)..], 1);
    }
}

pub const KronOp = enum { mul, div, rdiv };

/// KRON: the Kronecker (tensor) product of the `ma`x`na` matrix `a` (must
/// be tightly packed - leading dimension `ma`, no padding) and the
/// `m`x`n` matrix `b` (leading dimension `ldb`), written into `ar`/`ai`
/// (`(ma*m)`x`(na*n)`, leading dimension `ma*m`). Block `(ia, ja)` of the
/// result is `A(ia,ja)` combined element-by-element with all of `B`, via
/// `op`: `.mul` is the standard Kronecker product; `.div`/`.rdiv` apply
/// the same block structure with elementwise division instead
/// (`A(ia,ja) / B(i,j)` or `B(i,j) / A(ia,ja)`), matching the original's
/// shared implementation of all three.
pub fn kron(
    ar_a: []const f64,
    ai_a: []const f64,
    ma: usize,
    na: usize,
    ar_b: []const f64,
    ai_b: []const f64,
    ldb: usize,
    m: usize,
    n: usize,
    ar: []f64,
    ai: []f64,
    op: KronOp,
) void {
    const lda_result = ma * m;
    var ja: usize = 1;
    while (ja <= na) : (ja += 1) {
        var j: usize = 1;
        while (j <= n) : (j += 1) {
            var ia: usize = 1;
            while (ia <= ma) : (ia += 1) {
                const a_val = Cx.init(ar_a[ix(ma, ia, ja)], ai_a[ix(ma, ia, ja)]);
                var i: usize = 1;
                while (i <= m) : (i += 1) {
                    const b_val = Cx.init(ar_b[ix(ldb, i, j)], ai_b[ix(ldb, i, j)]);
                    const v = switch (op) {
                        .mul => blas1.wmul(a_val, b_val),
                        .div => blas1.wdiv(a_val, b_val),
                        .rdiv => blas1.wdiv(b_val, a_val),
                    };
                    const row = (ia - 1) * m + i;
                    const col = (ja - 1) * n + j;
                    ar[ix(lda_result, row, col)] = v.re;
                    ai[ix(lda_result, row, col)] = v.im;
                }
            }
        }
    }
}

const std = @import("std");
const testing = std.testing;

test "hilber matches the known inverse of the 3x3 Hilbert matrix" {
    // inv(hilb(3)) = [[9,-36,30],[-36,192,-180],[30,-180,180]].
    var a: [9]f64 = @splat(0.0);
    hilber(&a, 3, 3);
    try testing.expectApproxEqAbs(9.0, a[ix(3, 1, 1)], 1e-6);
    try testing.expectApproxEqAbs(-36.0, a[ix(3, 1, 2)], 1e-6);
    try testing.expectApproxEqAbs(30.0, a[ix(3, 1, 3)], 1e-6);
    try testing.expectApproxEqAbs(192.0, a[ix(3, 2, 2)], 1e-6);
    try testing.expectApproxEqAbs(-180.0, a[ix(3, 2, 3)], 1e-6);
    try testing.expectApproxEqAbs(180.0, a[ix(3, 3, 3)], 1e-6);
}

test "eye and ones produce the expected 2x3 matrices" {
    var ar: [6]f64 = @splat(0.0);
    var ai: [6]f64 = @splat(0.0);
    eye(&ar, &ai, 2, 2, 3);
    try testing.expectApproxEqAbs(1.0, ar[ix(2, 1, 1)], 1e-12);
    try testing.expectApproxEqAbs(0.0, ar[ix(2, 2, 1)], 1e-12);
    try testing.expectApproxEqAbs(1.0, ar[ix(2, 2, 2)], 1e-12);
    try testing.expectApproxEqAbs(0.0, ar[ix(2, 1, 3)], 1e-12);

    ones(&ar, &ai, 2, 2, 3);
    for (ar) |v| try testing.expectApproxEqAbs(1.0, v, 1e-12);
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
        try testing.expectApproxEqAbs(target, row_sum, 1e-9);
        try testing.expectApproxEqAbs(target, col_sum, 1e-9);
    }
    var diag1: f64 = 0.0;
    var diag2: f64 = 0.0;
    var k: usize = 1;
    while (k <= n) : (k += 1) {
        diag1 += a[ix(n, k, k)];
        diag2 += a[ix(n, k, n + 1 - k)];
    }
    try testing.expectApproxEqAbs(target, diag1, 1e-9);
    try testing.expectApproxEqAbs(target, diag2, 1e-9);
}

test "magic produces valid squares for odd, singly-even and doubly-even orders" {
    var a5: [25]f64 = @splat(0.0);
    magic(&a5, 5, 5);
    try checkMagicSquare(&a5, 5);

    var a6: [36]f64 = @splat(0.0);
    magic(&a6, 6, 6);
    try checkMagicSquare(&a6, 6);

    var a8: [64]f64 = @splat(0.0);
    magic(&a8, 8, 8);
    try checkMagicSquare(&a8, 8);
}

test "diagFromVector and diagToVector round-trip" {
    const vr = [_]f64{ 1.0, 2.0, 3.0 };
    const vi = [_]f64{ 0.0, 0.0, 0.0 };
    var ar: [16]f64 = @splat(0.0); // room for order 4 (n=3, k=1)
    var ai: [16]f64 = @splat(0.0);
    diagFromVector(&vr, &vi, 3, 1, &ar, &ai);
    // order = 3+1 = 4; diagonal k=1 means A(i, i+1) = v(i).
    try testing.expectApproxEqAbs(1.0, ar[ix(4, 1, 2)], 1e-12);
    try testing.expectApproxEqAbs(2.0, ar[ix(4, 2, 3)], 1e-12);
    try testing.expectApproxEqAbs(3.0, ar[ix(4, 3, 4)], 1e-12);

    var dr: [3]f64 = @splat(0.0);
    var di: [3]f64 = @splat(0.0);
    const len = diagToVector(&ar, &ai, 4, 4, 4, 1, &dr, &di);
    try testing.expectEqual(@as(usize, 3), len);
    try testing.expectApproxEqAbs(1.0, dr[0], 1e-12);
    try testing.expectApproxEqAbs(2.0, dr[1], 1e-12);
    try testing.expectApproxEqAbs(3.0, dr[2], 1e-12);
}

test "tril and triu zero the expected triangle" {
    var ar = [_]f64{ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 }; // 3x3, col-major
    var ai: [9]f64 = @splat(0.0);
    var ar2 = ar;
    var ai2 = ai;

    tril(&ar, &ai, 3, 3, 3, 0);
    try testing.expectApproxEqAbs(0.0, ar[ix(3, 1, 2)], 1e-12); // above diagonal zeroed
    try testing.expectApproxEqAbs(2.0, ar[ix(3, 2, 1)], 1e-12); // below diagonal kept

    triu(&ar2, &ai2, 3, 3, 3, 0);
    try testing.expectApproxEqAbs(0.0, ar2[ix(3, 2, 1)], 1e-12); // below diagonal zeroed
    try testing.expectApproxEqAbs(4.0, ar2[ix(3, 1, 2)], 1e-12); // above diagonal kept
}

test "kron matches the textbook 2x2 example" {
    // A = [[1,2],[3,4]], B = [[0,5],[6,7]].
    const ar_a = [_]f64{ 1.0, 3.0, 2.0, 4.0 };
    const ai_a: [4]f64 = @splat(0.0);
    const ar_b = [_]f64{ 0.0, 6.0, 5.0, 7.0 };
    const ai_b: [4]f64 = @splat(0.0);
    var ar: [16]f64 = @splat(0.0);
    var ai: [16]f64 = @splat(0.0);
    kron(&ar_a, &ai_a, 2, 2, &ar_b, &ai_b, 2, 2, 2, &ar, &ai, .mul);

    // Expected kron(A,B) (row-major for readability):
    // [ 0  5 |  0 10]
    // [ 6  7 | 12 14]
    // [ 0 15 |  0 20]
    // [18 21 | 24 28]
    try testing.expectApproxEqAbs(0.0, ar[ix(4, 1, 1)], 1e-12);
    try testing.expectApproxEqAbs(5.0, ar[ix(4, 1, 2)], 1e-12);
    try testing.expectApproxEqAbs(10.0, ar[ix(4, 1, 4)], 1e-12);
    try testing.expectApproxEqAbs(18.0, ar[ix(4, 4, 1)], 1e-12);
    try testing.expectApproxEqAbs(28.0, ar[ix(4, 4, 4)], 1e-12);
}

test "rref reduces a simple rank-deficient system" {
    // [[1, 2, 3], [2, 4, 6]] -> rank 1, row-reduces to [[1,2,3],[0,0,0]].
    var ar = [_]f64{ 1.0, 2.0, 2.0, 4.0, 3.0, 6.0 }; // column-major, 2x3
    var ai: [6]f64 = @splat(0.0);
    rref(&ar, &ai, 2, 2, 3, 1e-12);
    try testing.expectApproxEqAbs(1.0, ar[ix(2, 1, 1)], 1e-9);
    try testing.expectApproxEqAbs(2.0, ar[ix(2, 1, 2)], 1e-9);
    try testing.expectApproxEqAbs(3.0, ar[ix(2, 1, 3)], 1e-9);
    try testing.expectApproxEqAbs(0.0, ar[ix(2, 2, 1)], 1e-9);
    try testing.expectApproxEqAbs(0.0, ar[ix(2, 2, 2)], 1e-9);
    try testing.expectApproxEqAbs(0.0, ar[ix(2, 2, 3)], 1e-9);
}
