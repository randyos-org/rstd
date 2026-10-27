//! Orthogonal (unitary) reduction of a general complex matrix to upper
//! Hessenberg form, ported from `matlab88.f90`'s CORTH (EISPACK's complex
//! analogue of ORTHES). See `blas1.zig` for the shared matrix/vector
//! conventions this port uses, and `CxMat.zig` for the `self`-shaped
//! calling convention matrix-level functions like this one now use -
//! `low`/`igh` stay loose parameters since they describe a sub-range of
//! `self`, not `self`'s own shape.

const blas1 = @import("blas1.zig");
const CxMat = @import("CxMat.zig");
const ix = blas1.ix;
const rsimd = @import("../../simd/__root__.zig");

/// CORTH: reduce the submatrix of `self` (order `self.n`, leading
/// dimension `self.lda`) in rows/columns `low..igh` to upper Hessenberg
/// form by unitary similarity transformations. If a prior balancing step
/// was not used, pass `low = 1, igh = self.n`.
///
/// On return, `self` holds the Hessenberg matrix, with information about
/// the transformations in the elements below it; `ortr`/`orti` (indexed
/// `low..igh`) hold further information needed to recover the
/// transformations (e.g. for `comqr3`'s eigenvector back-transformation).
pub fn corth(self: CxMat, low: usize, igh: usize, ortr: []f64, orti: []f64) void {
    const nm = self.lda;
    const n = self.n;
    const ar = self.r;
    const ai = self.i;

    if (igh < low + 2) return; // la = igh-1 < kp1 = low+1: nothing to do

    var m: usize = low + 1;
    while (m <= igh - 1) : (m += 1) {
        var h: f64 = 0.0;
        ortr[m - 1] = 0.0;
        orti[m - 1] = 0.0;
        // scale = sum(|ar|+|ai|) over rows m..igh of column m-1: contiguous
        // (unit-stride) since row varies with column fixed in this
        // column-major layout, so it's safe to route through mapReduce.
        const col_start = ix(nm, m, m - 1);
        const col_len = igh - m + 1;
        const scale = rsimd.reduce.mapReduce(.Add, f64, blas1.cabs1, .{ ar[col_start..][0..col_len], ai[col_start..][0..col_len] });
        if (scale == 0.0) continue;

        const mp = m + igh;
        var ii: usize = m;
        while (ii <= igh) : (ii += 1) {
            const i = mp - ii;
            ortr[i - 1] = ar[ix(nm, i, m - 1)] / scale;
            orti[i - 1] = ai[ix(nm, i, m - 1)] / scale;
            h += ortr[i - 1] * ortr[i - 1] + orti[i - 1] * orti[i - 1];
        }

        var g = @sqrt(h);
        const f0 = blas1.pythag(ortr[m - 1], orti[m - 1]);
        if (f0 != 0.0) {
            h += f0 * g;
            g /= f0;
            ortr[m - 1] = (1.0 + g) * ortr[m - 1];
            orti[m - 1] = (1.0 + g) * orti[m - 1];
        } else {
            ortr[m - 1] = g;
            ar[ix(nm, m, m - 1)] = scale;
        }

        // Form (I - (u*ut)/h) * A.
        var j: usize = m;
        while (j <= n) : (j += 1) {
            var fr: f64 = 0.0;
            var fi: f64 = 0.0;
            var ii2: usize = m;
            while (ii2 <= igh) : (ii2 += 1) {
                const i = mp - ii2;
                fr += ortr[i - 1] * ar[ix(nm, i, j)] + orti[i - 1] * ai[ix(nm, i, j)];
                fi += ortr[i - 1] * ai[ix(nm, i, j)] - orti[i - 1] * ar[ix(nm, i, j)];
            }
            fr /= h;
            fi /= h;
            var ri: usize = m;
            while (ri <= igh) : (ri += 1) {
                ar[ix(nm, ri, j)] -= fr * ortr[ri - 1] - fi * orti[ri - 1];
                ai[ix(nm, ri, j)] -= fr * orti[ri - 1] + fi * ortr[ri - 1];
            }
        }

        // Form (I - (u*ut)/h) * A * (I - (u*ut)/h).
        var row2: usize = 1;
        while (row2 <= igh) : (row2 += 1) {
            var fr: f64 = 0.0;
            var fi: f64 = 0.0;
            var jj: usize = m;
            while (jj <= igh) : (jj += 1) {
                const col = mp - jj;
                fr += ortr[col - 1] * ar[ix(nm, row2, col)] - orti[col - 1] * ai[ix(nm, row2, col)];
                fi += ortr[col - 1] * ai[ix(nm, row2, col)] + orti[col - 1] * ar[ix(nm, row2, col)];
            }
            fr /= h;
            fi /= h;
            var col2: usize = m;
            while (col2 <= igh) : (col2 += 1) {
                ar[ix(nm, row2, col2)] -= fr * ortr[col2 - 1] + fi * orti[col2 - 1];
                ai[ix(nm, row2, col2)] += fr * orti[col2 - 1] - fi * ortr[col2 - 1];
            }
        }

        ortr[m - 1] *= scale;
        orti[m - 1] *= scale;
        ar[ix(nm, m, m - 1)] = -(g * ar[ix(nm, m, m - 1)]);
        ai[ix(nm, m, m - 1)] = -(g * ai[ix(nm, m, m - 1)]);
    }
}

const std = @import("std");

const zigstd = @import("../../zigstd.zig");

test "corth preserves trace (similarity transform)" {
    // A = [[1,2,3],[4,5,6],[7,8,10]] (real, column-major).
    //
    // Note: unlike `htridi`, `corth` does not zero the strict lower
    // triangle below the subdiagonal - only `A(m, m-1)` itself becomes
    // the new Hessenberg subdiagonal entry; entries strictly below that
    // keep leftover reduction data rather than the true (zero) similarity-
    // transformed value, per the original's "information ... is stored in
    // the remaining triangles" docs. So Frobenius norm over the raw
    // storage is *not* an invariant here (only trace is, cheaply); full
    // validation of `corth` happens end-to-end once `comqr3` (which
    // consumes exactly this storage format) is in place.
    var ar = [_]f64{ 1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 10.0 };
    var ai: [9]f64 = @splat(0.0);
    var ortr: [3]f64 = @splat(0.0);
    var orti: [3]f64 = @splat(0.0);

    corth(CxMat.initPacked(&ar, &ai, 3, 3), 1, 3, &ortr, &orti);

    const trace = ar[ix(3, 1, 1)] + ar[ix(3, 2, 2)] + ar[ix(3, 3, 3)];
    try zigstd.testExpectApproxEqAbs(16.0, trace, 1e-9);
}
