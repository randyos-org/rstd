//! The unblocked kernels of `sytrf`, `hetrf` and their `_rook` forms
//! (`relapack.zig`): Bunch and Kaufman's diagonal pivoting, or the bounded
//! ("rook") form, `A = L D Lᵀ` (Hermitian: `L D Lᴴ`) with `D` of 1 x 1 and
//! 2 x 2 blocks, in LAPACK's storage and `ipiv`.
//!
//! - `tf2` is Reference-LAPACK's `xSYTF2`, `xHETF2`, `xSYTF2_ROOK` and
//!   `xHETF2_ROOK` (`lapack-netlib/SRC/dsytf2.f`, `zsytf2.f`, `zhetf2.f`,
//!   `dsytf2_rook.f`, `zsytf2_rook.f`, `zhetf2_rook.f`, as OpenBLAS
//!   v0.3.34-317 ships them).
//! - `panel` is ReLAPACK's `RELAPACK_?sytrf_rec2`, `?hetrf_rec2` and their
//!   `_rook_rec2` forms (`relapack/src/dsytrf_rec2.c`, `zsytrf_rec2.c`,
//!   `zhetrf_rec2.c`, `dsytrf_rook_rec2.c`, `zsytrf_rook_rec2.c`,
//!   `zhetrf_rook_rec2.c` in that checkout), Reference-LAPACK's `xLASYF`,
//!   `xLAHEF` and `_ROOK` forms less their trailing update, with the closing
//!   interchanges split off as `undo`: `xLASYF` updates the rest before it
//!   puts the panel in LAPACK's form, which `rec2` does first (`issues.md`).
//!
//! One body serves every form: `herm` (complex `T` only) conjugates and
//! keeps the diagonal real where the reference's `xHE` forms do, and `rook`
//! takes the bounded search and its small-pivot division. The upper
//! triangle runs the lower algorithm on the matrix read backwards (rows and
//! columns reversed), which the reference's upper loops mirror; ties in the
//! pivot search and the order of `gemv`'s sums are mirrored with it.
//! `i?amax` here takes the first NaN, else the first largest `|re| + |im|`,
//! as rstd's `iamax` kernels. The level-2 calls are loops in the reference
//! BLAS's order, `syr` and `her` skipping a zero `x(j)` as it does.

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const chol = @import("cholesky.zig");

const Float = types.Float;
const isComplex = types.isComplex;
const scalar = chol.scalar;
const add = chol.add;
const sub = chol.sub;
const mul = chol.mul;

/// `n x m` of a column-major matrix, read backwards (`flip`): `(i, j)` is
/// `(rows - 1 - i, cols - 1 - j)` of the matrix at `p`.
pub fn View(comptime T: type) type {
    return struct {
        p: [*]T,
        ld: usize,
        rows: usize,
        cols: usize,
        flip: bool,

        const Self = @This();

        pub inline fn at(v: Self, i: usize, j: usize) *T {
            const r = if (v.flip) v.rows - 1 - i else i;
            const c = if (v.flip) v.cols - 1 - j else j;
            return &v.p[r + c * v.ld];
        }

        pub inline fn get(v: Self, i: usize, j: usize) T {
            return v.at(i, j).*;
        }

        pub inline fn set(v: Self, i: usize, j: usize, x: T) void {
            v.at(i, j).* = x;
        }

        /// Rows `r` and `s` over columns `[0, cols)`.
        fn swapRows(v: Self, r: usize, s: usize, cols: usize) void {
            for (0..cols) |c| std.mem.swap(T, v.at(r, c), v.at(s, c));
        }
    };
}

/// `ipiv` of `n` pivots, as `View` reads the matrix: 1-based signed values of
/// the read order, stored as the matrix's.
pub const Piv = struct {
    p: [*]isize,
    n: usize,
    flip: bool,

    fn map(pv: Piv, v: isize) isize {
        if (!pv.flip) return v;
        const m: isize = @intCast(pv.n + 1);
        return if (v < 0) -(m + v) else m - v;
    }

    pub fn get(pv: Piv, k: usize) isize {
        return pv.map(pv.p[if (pv.flip) pv.n - 1 - k else k]);
    }

    pub fn set(pv: Piv, k: usize, v: isize) void {
        pv.p[if (pv.flip) pv.n - 1 - k else k] = pv.map(v);
    }

    /// `k` (0-based, read order) as the reference's `INFO`.
    pub fn info(pv: Piv, k: usize) usize {
        return if (pv.flip) pv.n - k else k + 1;
    }
};

/// The element arithmetic of the `herm` or symmetric forms.
fn El(comptime T: type, comptime herm_: bool) type {
    return struct {
        const herm = herm_ and isComplex(T);
        const F = Float(T);

        /// `|re| + |im|` (`CABS1`), `|x|` for real.
        fn abs1(v: T) F {
            return if (comptime isComplex(T)) @abs(v.re) + @abs(v.im) else @abs(v);
        }

        /// A diagonal element's size: `|re|` where it is real (`herm`).
        fn absd(v: T) F {
            return if (comptime herm) @abs(v.re) else abs1(v);
        }

        fn cj(v: T) T {
            return if (comptime herm) v.conjugate() else v;
        }

        /// `DBLE` where the diagonal is real.
        fn real(v: T) T {
            return if (comptime herm) .init(v.re, 0) else v;
        }

        fn div(x: T, y: T) T {
            if (comptime !isComplex(T)) return x / y;
            const z = util.cdiv(F, x.re, x.im, y.re, y.im);
            return .init(z[0], z[1]);
        }

        fn conj(v: T) T {
            return if (comptime isComplex(T)) v.conjugate() else v;
        }

        /// The first NaN, else the first largest `abs1` of column `j`'s rows
        /// `[lo, hi)`, `hi > lo`.
        fn iamaxCol(a: View(T), j: usize, lo: usize, hi: usize) usize {
            var best = lo;
            var m: F = -1;
            for (lo..hi) |i| {
                const x = abs1(a.get(i, j));
                if (std.math.isNan(x)) return i;
                if (x > m) {
                    m = x;
                    best = i;
                }
            }
            return best;
        }

        /// `iamaxCol` along row `i`'s columns `[lo, hi)`.
        fn iamaxRow(a: View(T), i: usize, lo: usize, hi: usize) usize {
            var best = lo;
            var m: F = -1;
            for (lo..hi) |j| {
                const x = abs1(a.get(i, j));
                if (std.math.isNan(x)) return j;
                if (x > m) {
                    m = x;
                    best = j;
                }
            }
            return best;
        }

        /// `xSYR` or `xHER` (`herm`) on the lower triangle of `a`'s trailing
        /// `[k, n)`: `A += alpha x xᵀ` (`x xᴴ`), `x` column `c`'s rows `[k, n)`.
        fn rank1(a: View(T), c: usize, k: usize, alpha: T) void {
            const n = a.rows;
            for (k..n) |j| {
                const xj = a.get(j, c);
                if (chol.isZero(T, xj)) {
                    if (comptime herm) a.set(j, j, real(a.get(j, j)));
                    continue;
                }
                const t = if (comptime herm) chol.mulReal(T, cj(xj), alpha.re) else mul(T, alpha, xj);
                if (comptime herm) {
                    const d = mul(T, xj, t);
                    a.set(j, j, .init(a.get(j, j).re + d.re, 0));
                } else a.set(j, j, add(T, a.get(j, j), mul(T, xj, t)));
                for (j + 1..n) |i| a.set(i, j, add(T, a.get(i, j), mul(T, a.get(i, c), t)));
            }
        }

        /// Rows and columns `kk < kp` of the lower triangle from `kk` on
        /// interchanged (`xSYTF2`'s or `xHETF2`'s swaps).
        fn swapSym(a: View(T), kk: usize, kp: usize) void {
            const n = a.rows;
            for (kp + 1..n) |i| std.mem.swap(T, a.at(i, kk), a.at(i, kp));
            for (kk + 1..kp) |j| {
                const t = cj(a.get(j, kk));
                a.set(j, kk, cj(a.get(kp, j)));
                a.set(kp, j, t);
            }
            if (comptime herm) {
                a.set(kp, kk, cj(a.get(kp, kk)));
                const r = real(a.get(kk, kk));
                a.set(kk, kk, real(a.get(kp, kp)));
                a.set(kp, kp, r);
            } else std.mem.swap(T, a.at(kk, kk), a.at(kp, kp));
        }
    };
}

/// `alpha` of the pivot rule, `(1 + sqrt 17) / 8`.
fn alphaOf(comptime F: type) F {
    return (1 + @sqrt(@as(F, 17))) / 8;
}

// ---- tf2 --------------------------------------------------------------------------

/// `xSYTF2`, `xHETF2` (`herm`) or their `_ROOK` forms on `uplo`'s triangle
/// of `n x n` `a`, `ipiv` its `n` pivots: the reference's `INFO`, the first
/// exactly zero (or NaN, outside `rook`) `D(k, k)`, or 0.
pub fn tf2(comptime T: type, comptime herm: bool, comptime rook: bool, upper: bool, n: usize, a_: [*]T, lda: usize, ipiv: [*]isize) usize {
    const E = El(T, herm);
    const F = Float(T);
    const a: View(T) = .{ .p = a_, .ld = lda, .rows = n, .cols = n, .flip = upper };
    const pv: Piv = .{ .p = ipiv, .n = n, .flip = upper };
    const alpha = alphaOf(F);
    const sfmin = std.math.floatMin(F);
    var info: usize = 0;
    var k: usize = 0;
    while (k < n) {
        var kstep: usize = 1;
        var p = k;
        var kp = k;
        const absakk = E.absd(a.get(k, k));
        var imax = k;
        var colmax: F = 0;
        if (k + 1 < n) {
            imax = E.iamaxCol(a, k, k + 1, n);
            colmax = E.abs1(a.get(imax, k));
        }
        if (@max(absakk, colmax) == 0 or (!rook and std.math.isNan(absakk))) {
            if (info == 0) info = pv.info(k);
            if (comptime E.herm) a.set(k, k, E.real(a.get(k, k)));
        } else {
            if (!rook) {
                if (!(absakk >= alpha * colmax)) {
                    var jmax = E.iamaxRow(a, imax, k, imax);
                    var rowmax = E.abs1(a.get(imax, jmax));
                    if (imax + 1 < n) {
                        jmax = E.iamaxCol(a, imax, imax + 1, n);
                        rowmax = @max(rowmax, E.abs1(a.get(jmax, imax)));
                    }
                    if (absakk >= alpha * colmax * (colmax / rowmax)) {
                        kp = k;
                    } else if (E.absd(a.get(imax, imax)) >= alpha * rowmax) {
                        kp = imax;
                    } else {
                        kp = imax;
                        kstep = 2;
                    }
                }
            } else if (absakk < alpha * colmax) while (true) {
                var jmax = k;
                var rowmax: F = 0;
                if (imax != k) {
                    jmax = E.iamaxRow(a, imax, k, imax);
                    rowmax = E.abs1(a.get(imax, jmax));
                }
                if (imax + 1 < n) {
                    const itemp = E.iamaxCol(a, imax, imax + 1, n);
                    const dtemp = E.abs1(a.get(itemp, imax));
                    if (dtemp > rowmax) {
                        rowmax = dtemp;
                        jmax = itemp;
                    }
                }
                if (!(E.absd(a.get(imax, imax)) < alpha * rowmax)) {
                    kp = imax;
                    break;
                } else if (p == jmax or rowmax <= colmax) {
                    kp = imax;
                    kstep = 2;
                    break;
                } else {
                    p = imax;
                    colmax = rowmax;
                    imax = jmax;
                }
            };
            const kk = k + kstep - 1;
            if (rook and kstep == 2 and p != k) E.swapSym(a, k, p);
            if (kp != kk) {
                E.swapSym(a, kk, kp);
                if (kstep == 2) {
                    if (comptime E.herm) a.set(k, k, E.real(a.get(k, k)));
                    std.mem.swap(T, a.at(k + 1, k), a.at(kp, k));
                }
            } else if (comptime E.herm) {
                a.set(k, k, E.real(a.get(k, k)));
                if (kstep == 2) a.set(k + 1, k + 1, E.real(a.get(k + 1, k + 1)));
            }
            if (kstep == 1) {
                if (k + 1 < n) oneByOne(T, herm, rook, a, k, sfmin);
            } else if (k + 2 < n) {
                twoByTwo(T, herm, rook, a, k);
            }
        }
        if (kstep == 1) {
            pv.set(k, @intCast(kp + 1));
        } else {
            pv.set(k, -@as(isize, @intCast((if (rook) p else kp) + 1)));
            pv.set(k + 1, -@as(isize, @intCast(kp + 1)));
        }
        k += kstep;
    }
    return info;
}

/// `tf2`'s 1 x 1 step: `L(k+1:, k)` and the rank-1 update of the rest.
fn oneByOne(comptime T: type, comptime herm: bool, comptime rook: bool, a: View(T), k: usize, sfmin: Float(T)) void {
    const E = El(T, herm);
    const n = a.rows;
    const akk = a.get(k, k);
    if (rook and !(E.absd(akk) >= sfmin)) {
        // divide first, then update by the pivot itself
        const d = E.real(akk);
        for (k + 1..n) |i| a.set(i, k, if (comptime E.herm) chol.divReal(T, a.get(i, k), d.re) else E.div(a.get(i, k), d));
        E.rank1(a, k, k + 1, chol.neg(T, d));
        return;
    }
    const r = E.div(scalar(T, 1), E.real(akk));
    E.rank1(a, k, k + 1, chol.neg(T, r));
    for (k + 1..n) |i| a.set(i, k, if (comptime E.herm) chol.mulReal(T, a.get(i, k), r.re) else mul(T, r, a.get(i, k)));
}

/// `tf2`'s 2 x 2 step: `L(k+2:, k:k+1)` and the rank-2 update of the rest.
fn twoByTwo(comptime T: type, comptime herm: bool, comptime rook: bool, a: View(T), k: usize) void {
    const E = El(T, herm);
    const F = Float(T);
    const n = a.rows;
    const one = scalar(T, 1);
    if (comptime E.herm) {
        const sub21 = a.get(k + 1, k);
        const dd = std.math.hypot(sub21.re, sub21.im);
        const d11 = a.get(k + 1, k + 1).re / dd;
        const d22 = a.get(k, k).re / dd;
        const tt = 1 / (d11 * d22 - 1);
        const d21 = chol.divReal(T, sub21, dd);
        // `xHETF2` folds `TT / D` into `D`; the rook form divides by `D` last
        const s: F = if (rook) tt else tt / dd;
        for (k + 2..n) |j| {
            const wk = chol.mulReal(T, sub(T, chol.mulReal(T, a.get(j, k), d11), mul(T, d21, a.get(j, k + 1))), s);
            const wkp1 = chol.mulReal(T, sub(T, chol.mulReal(T, a.get(j, k + 1), d22), mul(T, E.conj(d21), a.get(j, k))), s);
            for (j..n) |i| {
                const x = if (rook) chol.divReal(T, a.get(i, k), dd) else a.get(i, k);
                const y = if (rook) chol.divReal(T, a.get(i, k + 1), dd) else a.get(i, k + 1);
                a.set(i, j, sub(T, sub(T, a.get(i, j), mul(T, x, E.conj(wk))), mul(T, y, E.conj(wkp1))));
            }
            a.set(j, k, if (rook) chol.divReal(T, wk, dd) else wk);
            a.set(j, k + 1, if (rook) chol.divReal(T, wkp1, dd) else wkp1);
            a.set(j, j, E.real(a.get(j, j)));
        }
        return;
    }
    var d21 = a.get(k + 1, k);
    const d11 = E.div(a.get(k + 1, k + 1), d21);
    const d22 = E.div(a.get(k, k), d21);
    const t = E.div(one, sub(T, mul(T, d11, d22), one));
    if (!rook) d21 = E.div(t, d21);
    for (k + 2..n) |j| {
        const xk = sub(T, mul(T, d11, a.get(j, k)), a.get(j, k + 1));
        const xk1 = sub(T, mul(T, d22, a.get(j, k + 1)), a.get(j, k));
        const wk = if (rook) mul(T, t, xk) else mul(T, d21, xk);
        const wkp1 = if (rook) mul(T, t, xk1) else mul(T, d21, xk1);
        for (j..n) |i| {
            const x = if (rook) E.div(a.get(i, k), d21) else a.get(i, k);
            const y = if (rook) E.div(a.get(i, k + 1), d21) else a.get(i, k + 1);
            a.set(i, j, sub(T, sub(T, a.get(i, j), mul(T, x, wk)), mul(T, y, wkp1)));
        }
        a.set(j, k, if (rook) E.div(wk, d21) else wk);
        a.set(j, k + 1, if (rook) E.div(wkp1, d21) else wkp1);
    }
}

// ---- panel ------------------------------------------------------------------------

/// What `panel` returns: the columns it factored and the reference's `INFO`.
pub const Panel = struct { kb: usize, info: usize };

/// `RELAPACK_?sytrf_rec2` and the rest on `upper`'s side of the `n x n` `a`:
/// up to `nb` columns (lower: the first; upper: the last), one fewer where a
/// 2 x 2 pivot would cross the edge, unless `nb >= n`. `w` (`n x nb`, `ldw`)
/// takes `W = L D` (Hermitian: its conjugate) for the caller's update of the
/// rest, `A -= L Wᵀ`; the factored columns hold every interchange of the
/// panel until the caller's `undo`, which follows that update.
pub fn panel(comptime T: type, comptime herm: bool, comptime rook: bool, upper: bool, n: usize, nb: usize, a_: [*]T, lda: usize, ipiv: [*]isize, w_: [*]T, ldw: usize) Panel {
    const E = El(T, herm);
    const F = Float(T);
    const a: View(T) = .{ .p = a_, .ld = lda, .rows = n, .cols = n, .flip = upper };
    const w: View(T) = .{ .p = w_, .ld = ldw, .rows = n, .cols = nb, .flip = upper };
    const pv: Piv = .{ .p = ipiv, .n = n, .flip = upper };
    const alpha = alphaOf(F);
    const sfmin = std.math.floatMin(F);
    var info: usize = 0;
    var k: usize = 0;
    while (!((k + 1 >= nb and nb < n) or k >= n)) {
        var kstep: usize = 1;
        var p = k;
        var kp = k;
        // W(k:, k) = A(k:, k) - A(k:, :k) W(k, :k)ᵀ
        for (k..n) |i| w.set(i, k, a.get(i, k));
        if (comptime E.herm) w.set(k, k, E.real(w.get(k, k)));
        update(T, a, w, k, k, k);
        if (comptime E.herm) w.set(k, k, E.real(w.get(k, k)));
        const absakk = E.absd(w.get(k, k));
        var imax = k;
        var colmax: F = 0;
        if (k + 1 < n) {
            imax = E.iamaxCol(w, k, k + 1, n);
            colmax = E.abs1(w.get(imax, k));
        }
        if (@max(absakk, colmax) == 0) {
            if (info == 0) info = pv.info(k);
            if (rook) {
                a.set(k, k, E.real(w.get(k, k)));
                for (k + 1..n) |i| a.set(i, k, w.get(i, k));
            } else if (comptime E.herm) a.set(k, k, E.real(a.get(k, k)));
        } else {
            if (if (rook) absakk < alpha * colmax else !(absakk >= alpha * colmax)) while (true) {
                // W(k:, k+1) = row and column `imax` of A, updated
                for (k..imax) |j| w.set(j, k + 1, E.cj(a.get(imax, j)));
                w.set(imax, k + 1, E.real(a.get(imax, imax)));
                for (imax + 1..n) |i| w.set(i, k + 1, a.get(i, imax));
                update(T, a, w, k, k + 1, imax);
                if (comptime E.herm) w.set(imax, k + 1, E.real(w.get(imax, k + 1)));
                var jmax = k;
                var rowmax: F = 0;
                if (imax != k) {
                    jmax = E.iamaxCol(w, k + 1, k, imax);
                    rowmax = E.abs1(w.get(jmax, k + 1));
                }
                if (imax + 1 < n) {
                    const itemp = E.iamaxCol(w, k + 1, imax + 1, n);
                    const dtemp = E.abs1(w.get(itemp, k + 1));
                    if (!rook) {
                        rowmax = @max(rowmax, dtemp);
                    } else if (dtemp > rowmax) {
                        rowmax = dtemp;
                        jmax = itemp;
                    }
                }
                if (!rook) {
                    if (absakk >= alpha * colmax * (colmax / rowmax)) {
                        kp = k;
                    } else if (E.absd(w.get(imax, k + 1)) >= alpha * rowmax) {
                        kp = imax;
                        for (k..n) |i| w.set(i, k, w.get(i, k + 1));
                    } else {
                        kp = imax;
                        kstep = 2;
                    }
                    break;
                }
                if (!(E.absd(w.get(imax, k + 1)) < alpha * rowmax)) {
                    kp = imax;
                    for (k..n) |i| w.set(i, k, w.get(i, k + 1));
                    break;
                } else if (p == jmax or rowmax <= colmax) {
                    kp = imax;
                    kstep = 2;
                    break;
                } else {
                    p = imax;
                    colmax = rowmax;
                    imax = jmax;
                    for (k..n) |i| w.set(i, k, w.get(i, k + 1));
                }
            };
            const kk = k + kstep - 1;
            if (rook and kstep == 2 and p != k) interchange(T, herm, a, w, k, p, k, kk);
            if (kp != kk) interchange(T, herm, a, w, kk, kp, k, kk);
            if (kstep == 1) {
                for (k..n) |i| a.set(i, k, w.get(i, k));
                if (k + 1 < n) {
                    const akk = a.get(k, k);
                    if (rook and !(E.absd(akk) >= sfmin)) {
                        // `dsytrf_rook_rec2.c` divides only by a nonzero pivot;
                        // `zhetrf_rook_rec2.c`'s real one is never zero here
                        if (comptime E.herm) {
                            for (k + 1..n) |i| a.set(i, k, chol.divReal(T, a.get(i, k), akk.re));
                        } else if (!chol.isZero(T, akk)) {
                            for (k + 1..n) |i| a.set(i, k, E.div(a.get(i, k), akk));
                        }
                    } else {
                        const r = E.div(scalar(T, 1), E.real(akk));
                        for (k + 1..n) |i| a.set(i, k, if (comptime E.herm) chol.mulReal(T, a.get(i, k), r.re) else mul(T, r, a.get(i, k)));
                    }
                    if (comptime E.herm) for (k + 1..n) |i| w.set(i, k, E.conj(w.get(i, k)));
                }
            } else {
                if (k + 2 < n) panel2(T, herm, rook, a, w, k);
                a.set(k, k, w.get(k, k));
                a.set(k + 1, k, w.get(k + 1, k));
                a.set(k + 1, k + 1, w.get(k + 1, k + 1));
                if (comptime E.herm) {
                    for (k + 1..n) |i| w.set(i, k, E.conj(w.get(i, k)));
                    for (k + 2..n) |i| w.set(i, k + 1, E.conj(w.get(i, k + 1)));
                }
            }
        }
        if (kstep == 1) {
            pv.set(k, @intCast(kp + 1));
        } else {
            pv.set(k, -@as(isize, @intCast((if (rook) p else kp) + 1)));
            pv.set(k + 1, -@as(isize, @intCast(kp + 1)));
        }
        k += kstep;
    }
    return .{ .kb = k, .info = info };
}

/// The close of `xLASYF` and `rec2`: each step's interchanges undone in the
/// columns left of it (upper: right), last step first, the panel's `kb`
/// columns of the `n x n` `a` having been factored, so each column of `L`
/// holds only its own and earlier interchanges (LAPACK's form).
pub fn undoPanel(comptime T: type, comptime rook: bool, upper: bool, n: usize, a: [*]T, lda: usize, ipiv: [*]isize, kb: usize) void {
    const v: View(T) = .{ .p = a, .ld = lda, .rows = n, .cols = n, .flip = upper };
    undo(T, rook, v, .{ .p = ipiv, .n = n, .flip = upper }, kb);
}

/// `W(k:, c) -= A(k:, :k) W(r, :k)ᵀ`, `xGEMV`'s columns in order.
fn update(comptime T: type, a: View(T), w: View(T), k: usize, c: usize, r: usize) void {
    for (0..k) |j| {
        const t = w.get(r, j);
        for (k..a.rows) |i| w.set(i, c, sub(T, w.get(i, c), mul(T, a.get(i, j), t)));
    }
}

/// `panel`'s interchange of rows and columns `kk < kp` of `A` from `kk`
/// on (the step at `k` not yet written, so its column is left), of rows
/// `kk` and `kp` left of `k`, and of `W`'s rows up to column `top`.
fn interchange(comptime T: type, comptime herm: bool, a: View(T), w: View(T), kk: usize, kp: usize, k: usize, top: usize) void {
    const E = El(T, herm);
    const n = a.rows;
    a.set(kp, kp, E.real(a.get(kk, kk)));
    for (kk + 1..kp) |j| a.set(kp, j, E.cj(a.get(j, kk)));
    for (kp + 1..n) |i| a.set(i, kp, a.get(i, kk));
    a.swapRows(kk, kp, k);
    w.swapRows(kk, kp, top + 1);
}

/// `panel`'s 2 x 2 step: `L(k+2:, k:k+1)` from `W`.
fn panel2(comptime T: type, comptime herm: bool, comptime rook: bool, a: View(T), w: View(T), k: usize) void {
    const E = El(T, herm);
    const n = a.rows;
    const one = scalar(T, 1);
    var d21 = w.get(k + 1, k);
    const d11 = E.div(w.get(k + 1, k + 1), d21);
    const d22 = E.div(w.get(k, k), E.cj(d21));
    // real for the Hermitian forms (`1 / (re(d11 d22) - 1)`)
    const t: T = if (comptime E.herm) scalar(T, 1 / (mul(T, d11, d22).re - 1)) else E.div(one, sub(T, mul(T, d11, d22), one));
    if (!rook) d21 = E.div(t, d21);
    for (k + 2..n) |j| {
        const xk = sub(T, mul(T, d11, w.get(j, k)), w.get(j, k + 1));
        const xk1 = sub(T, mul(T, d22, w.get(j, k + 1)), w.get(j, k));
        if (rook) {
            const u = E.div(xk, E.cj(d21));
            const v = E.div(xk1, d21);
            a.set(j, k, if (comptime E.herm) chol.mulReal(T, u, t.re) else mul(T, t, u));
            a.set(j, k + 1, if (comptime E.herm) chol.mulReal(T, v, t.re) else mul(T, t, v));
        } else {
            a.set(j, k, mul(T, E.cj(d21), xk));
            a.set(j, k + 1, mul(T, d21, xk1));
        }
    }
}

/// `undoPanel` on a view.
fn undo(comptime T: type, comptime rook: bool, a: View(T), pv: Piv, k: usize) void {
    if (k == 0) return;
    // 1-based, as the reference's loop
    var j: isize = @intCast(k);
    while (true) {
        var jj = j;
        var jp1: isize = 1;
        var two = false;
        var jp2 = pv.get(@intCast(j - 1));
        if (jp2 < 0) {
            jp2 = -jp2;
            j -= 1;
            jp1 = -pv.get(@intCast(j - 1));
            two = true;
        }
        j -= 1;
        if (jp2 != jj and j >= 1) a.swapRows(@intCast(jp2 - 1), @intCast(jj - 1), @intCast(j));
        jj -= 1;
        if (rook and two and jp1 != jj and j >= 1) a.swapRows(@intCast(jp1 - 1), @intCast(jj - 1), @intCast(j));
        if (j <= 1) return;
    }
}
