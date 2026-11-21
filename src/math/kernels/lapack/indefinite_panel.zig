//! The blocked panel of the symmetric indefinite factorizations (`xLASYF` and its forms); see `indefinite.zig`.

const std = @import("std");
const types = @import("../blas/types.zig");
const ops = @import("ops.zig");
const indefinite = @import("indefinite.zig");

const Float = types.Float;
const scalar = ops.scalar;
const sub = ops.sub;
const mul = ops.mul;
const El = indefinite.El;
const Piv = indefinite.Piv;
const Pivot = indefinite.Pivot;
const View = indefinite.View;
const alphaOf = indefinite.alphaOf;

/// What `panel` returns: the columns it factored and the reference's `INFO`.
pub const Panel = struct { kb: usize, info: ?usize };

/// `RELAPACK_?sytrf_rec2` and the rest on `upper`'s side of the `n x n` `a`:
/// up to `nb` columns (lower: the first; upper: the last), one fewer where a
/// 2 x 2 pivot would cross the edge, unless `nb >= n`. `w` (`n x nb`, `ldw`)
/// takes `W = L D` (Hermitian: its conjugate) for the caller's update of the
/// rest, `A -= L Wᵀ`; the factored columns hold every interchange of the
/// panel until the caller's `undo`, which follows that update. A column the
/// update made zero is copied back from `W`, as `xLASYF` and `xLAHEF` do
/// (`lapack-netlib/SRC/dlasyf.f`, `zlahef.f`; OpenBLAS `ce4ece13e`), where
/// `rec2` outside `rook` leaves the stale one.
pub fn panel(comptime T: type, comptime herm: bool, comptime rook: bool, upper: bool, n: usize, nb: usize, a_: [*]T, lda: usize, ipiv: [*]Pivot, w_: [*]T, ldw: usize) Panel {
    const E = El(T, herm);
    const F = Float(T);
    const a: View(T) = .{ .p = a_, .ld = lda, .rows = n, .cols = n, .flip = upper };
    const w: View(T) = .{ .p = w_, .ld = ldw, .rows = n, .cols = nb, .flip = upper };
    const pv: Piv = .{ .p = ipiv, .n = n, .flip = upper };
    const alpha = alphaOf(F);
    const sfmin = std.math.floatMin(F);
    var info: ?usize = null;
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
            colmax = ops.abs1(T, w.get(imax, k));
        }
        if (@max(absakk, colmax) == 0) {
            if (info == null) info = pv.index(k);
            a.set(k, k, E.real(w.get(k, k)));
            for (k + 1..n) |i| a.set(i, k, w.get(i, k));
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
                    rowmax = ops.abs1(T, w.get(jmax, k + 1));
                }
                if (imax + 1 < n) {
                    const itemp = E.iamaxCol(w, k + 1, imax + 1, n);
                    const dtemp = ops.abs1(T, w.get(itemp, k + 1));
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
                            for (k + 1..n) |i| a.set(i, k, ops.divReal(T, a.get(i, k), akk.re));
                        } else if (!ops.isZero(T, akk)) {
                            for (k + 1..n) |i| a.set(i, k, ops.div(T, a.get(i, k), akk));
                        }
                    } else {
                        const r = ops.div(T, scalar(T, 1), E.real(akk));
                        for (k + 1..n) |i| a.set(i, k, if (comptime E.herm) ops.mulReal(T, a.get(i, k), r.re) else mul(T, r, a.get(i, k)));
                    }
                    if (comptime E.herm) for (k + 1..n) |i| w.set(i, k, ops.conj(T, w.get(i, k)));
                }
            } else {
                if (k + 2 < n) panel2(T, herm, a, w, k);
                a.set(k, k, w.get(k, k));
                a.set(k + 1, k, w.get(k + 1, k));
                a.set(k + 1, k + 1, w.get(k + 1, k + 1));
                if (comptime E.herm) {
                    for (k + 1..n) |i| w.set(i, k, ops.conj(T, w.get(i, k)));
                    for (k + 2..n) |i| w.set(i, k + 1, ops.conj(T, w.get(i, k + 1)));
                }
            }
        }
        if (kstep == 1) {
            pv.set(k, .of(kp, .one));
        } else {
            pv.set(k, .of(if (rook) p else kp, .two));
            pv.set(k + 1, .of(kp, .two));
        }
        k += kstep;
    }
    return .{ .kb = k, .info = info };
}

/// The close of `xLASYF` and `rec2`: each step's interchanges undone in the
/// columns left of it (upper: right), last step first, the panel's `kb`
/// columns of the `n x n` `a` having been factored, so each column of `L`
/// holds only its own and earlier interchanges (LAPACK's form).
pub fn undoPanel(comptime T: type, comptime rook: bool, upper: bool, n: usize, a: [*]T, lda: usize, ipiv: [*]Pivot, kb: usize) void {
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

/// `panel`'s 2 x 2 step: `L(k+2:, k:k+1)` from `W`, each entry divided by
/// `d21`, then scaled by `t`, as `xLASYF`, `xLAHEF` (OpenBLAS `ce4ece13e`)
/// and the `_rook` forms; `rec2` outside `rook` folds `t / d21`, which
/// overflows for a tiny `d21`.
fn panel2(comptime T: type, comptime herm: bool, a: View(T), w: View(T), k: usize) void {
    const E = El(T, herm);
    const n = a.rows;
    const one = scalar(T, 1);
    const d21 = w.get(k + 1, k);
    const d11 = ops.div(T, w.get(k + 1, k + 1), d21);
    const d22 = ops.div(T, w.get(k, k), E.cj(d21));
    // real for the Hermitian forms (`1 / (re(d11 d22) - 1)`)
    const t: T = if (comptime E.herm) scalar(T, 1 / (mul(T, d11, d22).re - 1)) else ops.div(T, one, sub(T, mul(T, d11, d22), one));
    for (k + 2..n) |j| {
        const u = ops.div(T, sub(T, mul(T, d11, w.get(j, k)), w.get(j, k + 1)), E.cj(d21));
        const v = ops.div(T, sub(T, mul(T, d22, w.get(j, k + 1)), w.get(j, k)), d21);
        a.set(j, k, if (comptime E.herm) ops.mulReal(T, u, t.re) else mul(T, t, u));
        a.set(j, k + 1, if (comptime E.herm) ops.mulReal(T, v, t.re) else mul(T, t, v));
    }
}

/// `undoPanel` on a view.
fn undo(comptime T: type, comptime rook: bool, a: View(T), pv: Piv, k: usize) void {
    if (k == 0) return;
    // `j` 1-based, as the reference's loop
    var j = k;
    while (true) {
        var jj = j;
        const p2 = pv.get(j - 1);
        const two = p2.block == .two;
        var r1: usize = 0;
        if (two) {
            j -= 1;
            r1 = pv.get(j - 1).row;
        }
        j -= 1;
        if (p2.row != jj - 1 and j >= 1) a.swapRows(p2.row, jj - 1, j);
        jj -= 1;
        if (rook and two and r1 != jj - 1 and j >= 1) a.swapRows(r1, jj - 1, j);
        if (j <= 1) return;
    }
}
