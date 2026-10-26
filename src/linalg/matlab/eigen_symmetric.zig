//! Eigenvalues/eigenvectors of a real symmetric tridiagonal matrix by the
//! implicit QL method, ported from `matlab88.f90`'s IMTQL2 (EISPACK,
//! Martin/Wilkinson, as modified by Dubrulle). See `blas1.zig` for the
//! shared matrix/vector conventions this port uses.

const std = @import("std");
const blas1 = @import("blas1.zig");
const ix = blas1.ix;

/// IMTQL2: find the eigenvalues (always) and, if `job != 0`, eigenvectors
/// of the real symmetric tridiagonal matrix with diagonal `d` and
/// subdiagonal `e` (`e[0]` is arbitrary on input). Also finds the
/// eigenvectors of a full symmetric matrix if `htridi`+real `tred2`-style
/// reduction has been used to produce this tridiagonal form and `z` holds
/// the accumulated transformation (or the identity, for the tridiagonal
/// matrix's own eigenvectors).
///
/// On return, `d` holds the eigenvalues in ascending order and (if
/// `job != 0`) `z` (column-major, leading dimension `nm`, `n` columns)
/// holds the corresponding orthonormal eigenvectors as columns. `e` is
/// destroyed.
///
/// Returns 0 on success, or the 1-based index of an eigenvalue that failed
/// to converge after 30 iterations (eigenvalues/vectors for smaller
/// indices are still correct in that case).
pub fn imtql2(n: usize, d: []f64, e: []f64, z: []f64, nm: usize, job: usize) usize {
    if (n == 1) return 0;

    var i: usize = 2;
    while (i <= n) : (i += 1) e[i - 2] = e[i - 1];
    e[n - 1] = 0.0;

    var l: usize = 1;
    outer: while (l <= n) : (l += 1) {
        var iterations: usize = 0;
        while (true) {
            // Look for a small subdiagonal element.
            var m: usize = l;
            while (m <= n) : (m += 1) {
                if (m == n) break;
                const p_test = @abs(d[m - 1]) + @abs(d[m]);
                const s_test = p_test + @abs(e[m - 1]);
                if (p_test == s_test) break;
            }
            var p = d[l - 1];
            if (m == l) continue :outer;
            if (iterations == 30) return l;
            iterations += 1;

            // Form the shift.
            var g = (d[l] - p) / (2.0 * e[l - 1]);
            var r = @sqrt(g * g + 1.0);
            g = d[m - 1] - p + e[l - 1] / (g + std.math.copysign(r, g));
            var s: f64 = 1.0;
            var c: f64 = 1.0;
            p = 0.0;
            const mml = m - l;

            var count: usize = 1;
            while (count <= mml) : (count += 1) {
                const at = m - count;
                const f = s * e[at - 1];
                const b = c * e[at - 1];
                if (@abs(f) >= @abs(g)) {
                    c = g / f;
                    r = @sqrt(c * c + 1.0);
                    e[at] = f * r;
                    s = 1.0 / r;
                    c *= s;
                } else {
                    s = f / g;
                    r = @sqrt(s * s + 1.0);
                    e[at] = g * r;
                    c = 1.0 / r;
                    s *= c;
                }
                g = d[at] - p;
                r = (d[at - 1] - g) * s + 2.0 * c * b;
                p = s * r;
                d[at] = g + p;
                g = c * r - b;
                if (job != 0) {
                    var k: usize = 1;
                    while (k <= n) : (k += 1) {
                        const zf = z[ix(nm, k, at + 1)];
                        z[ix(nm, k, at + 1)] = s * z[ix(nm, k, at)] + c * zf;
                        z[ix(nm, k, at)] = c * z[ix(nm, k, at)] - s * zf;
                    }
                }
            }
            d[l - 1] -= p;
            e[l - 1] = g;
            e[m - 1] = 0.0;
        }
    }

    // Order eigenvalues (and eigenvectors) ascending by selection sort.
    var ii: usize = 2;
    while (ii <= n) : (ii += 1) {
        const lo = ii - 1;
        var k = lo;
        var p = d[lo - 1];
        var j: usize = ii;
        while (j <= n) : (j += 1) {
            if (d[j - 1] < p) {
                k = j;
                p = d[j - 1];
            }
        }
        if (k == lo) continue;
        d[k - 1] = d[lo - 1];
        d[lo - 1] = p;
        if (job == 0) continue;
        var row: usize = 1;
        while (row <= n) : (row += 1) {
            std.mem.swap(f64, &z[ix(nm, row, lo)], &z[ix(nm, row, k)]);
        }
    }
    return 0;
}

const zigstd = @import("../../zigstd.zig");

test "imtql2 finds eigenvalues/eigenvectors of a small tridiagonal matrix" {
    // Symmetric tridiagonal: diag [2,2,2], off-diag [1,1] -> eigenvalues
    // 2 - sqrt(2), 2, 2 + sqrt(2) (standard path-graph Laplacian-like case).
    var d = [_]f64{ 2.0, 2.0, 2.0 };
    var e = [_]f64{ 0.0, 1.0, 1.0 };
    var z: [9]f64 = @splat(0.0);
    z[ix(3, 1, 1)] = 1.0;
    z[ix(3, 2, 2)] = 1.0;
    z[ix(3, 3, 3)] = 1.0;

    const ierr = imtql2(3, &d, &e, &z, 3, 1);
    try zigstd.testExpectEqual(@as(usize, 0), ierr);

    const sqrt2 = @sqrt(2.0);
    try zigstd.testExpectApproxEqAbs(2.0 - sqrt2, d[0], 1e-10);
    try zigstd.testExpectApproxEqAbs(2.0, d[1], 1e-10);
    try zigstd.testExpectApproxEqAbs(2.0 + sqrt2, d[2], 1e-10);

    // Columns of z should be unit-norm (orthonormal eigenvectors).
    for (0..3) |col| {
        var norm2: f64 = 0.0;
        for (0..3) |row| {
            const v = z[ix(3, row + 1, col + 1)];
            norm2 += v * v;
        }
        try zigstd.testExpectApproxEqAbs(1.0, norm2, 1e-10);
    }
}
