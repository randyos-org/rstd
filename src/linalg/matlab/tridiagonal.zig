//! Householder reduction of a complex Hermitian matrix to real symmetric
//! tridiagonal form (and back-transformation of its eigenvectors), ported
//! from `matlab88.f90`'s HTRIDI/HTRIBK (EISPACK's complex analogues of
//! TRED1/TRBAK1). See `blas1.zig` for the shared matrix/vector conventions
//! this port uses.
//!
//! The original's `TAU(2,N)` work array (two rows of transformation data)
//! is split here into two separate length-N arrays `taur`/`taui`
//! (`TAU(1,*)`/`TAU(2,*)`), matching how every other multi-row real/complex
//! pair in this port is represented.

const blas1 = @import("blas1.zig");
const ix = blas1.ix;

/// HTRIDI: reduce the Hermitian matrix in the lower triangle of `ar`/`ai`
/// (column-major, leading dimension `nm`, order `n`) to a real symmetric
/// tridiagonal matrix by unitary similarity transformations.
///
/// On return: `d` holds the tridiagonal's diagonal; `e` holds its
/// subdiagonal in positions `2..n` (`e[0]` is set to 0); `e2` holds the
/// squares of `e`'s entries (may be the same array as `e` if not needed
/// separately); `ar`/`ai`'s lower triangle and `taur`/`taui` hold
/// information needed by `htribk` to recover eigenvectors.
pub fn htridi(nm: usize, n: usize, ar: []f64, ai: []f64, d: []f64, e: []f64, e2: []f64, taur: []f64, taui: []f64) void {
    taur[n - 1] = 1.0;
    taui[n - 1] = 0.0;
    var idx: usize = 1;
    while (idx <= n) : (idx += 1) d[idx - 1] = ar[ix(nm, idx, idx)];

    var ii: usize = 1;
    while (ii <= n) : (ii += 1) {
        const i = n + 1 - ii;
        const l = i - 1;
        var h: f64 = 0.0;
        var scale: f64 = 0.0;
        if (l >= 1) {
            var ksc: usize = 1;
            while (ksc <= l) : (ksc += 1) scale += @abs(ar[ix(nm, i, ksc)]) + @abs(ai[ix(nm, i, ksc)]);
        }

        if (l < 1 or scale == 0.0) {
            if (l >= 1) {
                taur[l - 1] = 1.0;
                taui[l - 1] = 0.0;
            }
            e[i - 1] = 0.0;
            e2[i - 1] = 0.0;
        } else {
            var ksc2: usize = 1;
            while (ksc2 <= l) : (ksc2 += 1) {
                ar[ix(nm, i, ksc2)] /= scale;
                ai[ix(nm, i, ksc2)] /= scale;
                h += ar[ix(nm, i, ksc2)] * ar[ix(nm, i, ksc2)] + ai[ix(nm, i, ksc2)] * ai[ix(nm, i, ksc2)];
            }
            e2[i - 1] = scale * scale * h;
            var g = @sqrt(h);
            e[i - 1] = scale * g;
            const f0 = blas1.pythag(ar[ix(nm, i, l)], ai[ix(nm, i, l)]);

            var si: f64 = undefined;
            var skip_reduction = false;
            if (f0 == 0.0) {
                taur[l - 1] = -taur[i - 1];
                si = taui[i - 1];
                ar[ix(nm, i, l)] = g;
            } else {
                taur[l - 1] = (ai[ix(nm, i, l)] * taui[i - 1] - ar[ix(nm, i, l)] * taur[i - 1]) / f0;
                si = (ar[ix(nm, i, l)] * taui[i - 1] + ai[ix(nm, i, l)] * taur[i - 1]) / f0;
                h += f0 * g;
                g = 1.0 + g / f0;
                ar[ix(nm, i, l)] *= g;
                ai[ix(nm, i, l)] *= g;
                if (l == 1) skip_reduction = true;
            }

            if (!skip_reduction) {
                var f: f64 = 0.0;
                var j: usize = 1;
                while (j <= l) : (j += 1) {
                    var gg: f64 = 0.0;
                    var gi: f64 = 0.0;
                    var k: usize = 1;
                    while (k <= j) : (k += 1) {
                        gg += ar[ix(nm, j, k)] * ar[ix(nm, i, k)] + ai[ix(nm, j, k)] * ai[ix(nm, i, k)];
                        gi += -ar[ix(nm, j, k)] * ai[ix(nm, i, k)] + ai[ix(nm, j, k)] * ar[ix(nm, i, k)];
                    }
                    if (l >= j + 1) {
                        var k2: usize = j + 1;
                        while (k2 <= l) : (k2 += 1) {
                            gg += ar[ix(nm, k2, j)] * ar[ix(nm, i, k2)] - ai[ix(nm, k2, j)] * ai[ix(nm, i, k2)];
                            gi += -ar[ix(nm, k2, j)] * ai[ix(nm, i, k2)] - ai[ix(nm, k2, j)] * ar[ix(nm, i, k2)];
                        }
                    }
                    e[j - 1] = gg / h;
                    taui[j - 1] = gi / h;
                    f += e[j - 1] * ar[ix(nm, i, j)] - taui[j - 1] * ai[ix(nm, i, j)];
                }
                const hh = f / (h + h);
                var j2: usize = 1;
                while (j2 <= l) : (j2 += 1) {
                    const fj = ar[ix(nm, i, j2)];
                    const gj = e[j2 - 1] - hh * fj;
                    e[j2 - 1] = gj;
                    const fij = -ai[ix(nm, i, j2)];
                    const gij = taui[j2 - 1] - hh * fij;
                    taui[j2 - 1] = -gij;
                    var k3: usize = 1;
                    while (k3 <= j2) : (k3 += 1) {
                        ar[ix(nm, j2, k3)] -= fj * e[k3 - 1] + gj * ar[ix(nm, i, k3)] - fij * taui[k3 - 1] - gij * ai[ix(nm, i, k3)];
                        ai[ix(nm, j2, k3)] -= fj * taui[k3 - 1] + gj * ai[ix(nm, i, k3)] + fij * e[k3 - 1] + gij * ar[ix(nm, i, k3)];
                    }
                }
            }

            var kres: usize = 1;
            while (kres <= l) : (kres += 1) {
                ar[ix(nm, i, kres)] *= scale;
                ai[ix(nm, i, kres)] *= scale;
            }
            taui[l - 1] = -si;
        }

        const hh_end = d[i - 1];
        d[i - 1] = ar[ix(nm, i, i)];
        ar[ix(nm, i, i)] = hh_end;
        ai[ix(nm, i, i)] = scale * @sqrt(h);
    }
}

/// HTRIBK: back-transform the first `m` columns of `zr` (eigenvectors of
/// the real symmetric tridiagonal matrix `htridi` produced) into the
/// eigenvectors of the original complex Hermitian matrix, using the
/// transformation data `htridi` left in `ar`/`ai`/`taur`/`taui`.
///
/// The last component of each returned eigenvector is real, and vector
/// Euclidean norms are preserved.
pub fn htribk(nm: usize, n: usize, ar: []const f64, ai: []const f64, taur: []const f64, taui: []const f64, m: usize, zr: []f64, zi: []f64) void {
    if (m == 0) return;

    var k: usize = 1;
    while (k <= n) : (k += 1) {
        var j: usize = 1;
        while (j <= m) : (j += 1) {
            zi[ix(nm, k, j)] = -(zr[ix(nm, k, j)] * taui[k - 1]);
            zr[ix(nm, k, j)] = zr[ix(nm, k, j)] * taur[k - 1];
        }
    }
    if (n == 1) return;

    var i: usize = 2;
    while (i <= n) : (i += 1) {
        const l = i - 1;
        const h = ai[ix(nm, i, i)];
        if (h == 0.0) continue;

        var j: usize = 1;
        while (j <= m) : (j += 1) {
            var s: f64 = 0.0;
            var si: f64 = 0.0;
            var k2: usize = 1;
            while (k2 <= l) : (k2 += 1) {
                s += ar[ix(nm, i, k2)] * zr[ix(nm, k2, j)] - ai[ix(nm, i, k2)] * zi[ix(nm, k2, j)];
                si += ar[ix(nm, i, k2)] * zi[ix(nm, k2, j)] + ai[ix(nm, i, k2)] * zr[ix(nm, k2, j)];
            }
            // Double divisions avoid possible underflow.
            s = (s / h) / h;
            si = (si / h) / h;
            var k3: usize = 1;
            while (k3 <= l) : (k3 += 1) {
                zr[ix(nm, k3, j)] -= s * ar[ix(nm, i, k3)] + si * ai[ix(nm, i, k3)];
                zi[ix(nm, k3, j)] -= si * ar[ix(nm, i, k3)] - s * ai[ix(nm, i, k3)];
            }
        }
    }
}

const std = @import("std");
const testing = std.testing;

test "htridi preserves trace and Frobenius norm of a real symmetric matrix" {
    // A = [[4,1,2],[1,3,0],[2,0,5]] (Hermitian similarity preserves both).
    var ar: [9]f64 = @splat(0.0);
    var ai: [9]f64 = @splat(0.0);
    ar[ix(3, 1, 1)] = 4.0;
    ar[ix(3, 2, 1)] = 1.0;
    ar[ix(3, 2, 2)] = 3.0;
    ar[ix(3, 3, 1)] = 2.0;
    ar[ix(3, 3, 2)] = 0.0;
    ar[ix(3, 3, 3)] = 5.0;

    var d: [3]f64 = undefined;
    var e: [3]f64 = undefined;
    var e2: [3]f64 = undefined;
    var taur: [3]f64 = undefined;
    var taui: [3]f64 = undefined;
    htridi(3, 3, &ar, &ai, &d, &e, &e2, &taur, &taui);

    const trace = d[0] + d[1] + d[2];
    try testing.expectApproxEqAbs(12.0, trace, 1e-10);

    const frob2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2] + 2.0 * (e[1] * e[1] + e[2] * e[2]);
    try testing.expectApproxEqAbs(60.0, frob2, 1e-10);
}
