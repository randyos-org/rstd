//! The kernels of `sytrf`, `hetrf` and their `_rook` forms (`sytrf.zig`):
//! Reference-LAPACK's `xSYTF2`, `xHETF2` and `_ROOK` forms (`tf2`;
//! `lapack-netlib/SRC/dsytf2.f` and the rest, OpenBLAS v0.3.34-317) and
//! ReLAPACK's `rec2` panels (`panel`, in `indefinite_panel.zig`;
//! `relapack/src/dsytrf_rec2.c` and the rest, MIT). One body serves every
//! form; the upper triangle runs the lower algorithm on the matrix read
//! backwards (design.md's LAPACK section).

const std = @import("std");
const mem = @import("../../../mem.zig");
const types = @import("../blas/types.zig");
const ops = @import("ops.zig");
const indefinite_panel = @import("indefinite_panel.zig");

const Float = types.Float;
const isComplex = types.isComplex;
const scalar = ops.scalar;
const add = ops.add;
const sub = ops.sub;
const mul = ops.mul;

pub const Panel = indefinite_panel.Panel;
pub const panel = indefinite_panel.panel;
pub const undoPanel = indefinite_panel.undoPanel;

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
        pub fn swapRows(v: Self, r: usize, s: usize, cols: usize) void {
            for (0..cols) |c| mem.swap(T, v.at(r, c), v.at(s, c));
        }
    };
}

/// One entry of `sytrf`'s pivots (issue #314): the row (0-based) swapped
/// with this one, and whether this row belongs to a 1 x 1 or a 2 x 2 block
/// of `D`. As small as a `usize`.
pub const Pivot = packed struct(usize) {
    row: Row,
    block: Block,

    pub const Block = enum(u1) { one, two };

    /// A row: one bit short of a `usize`.
    pub const Row = switch (@bitSizeOf(usize)) {
        64 => u63,
        32 => u31,
        16 => u15,
        else => @compileError("Pivot: an unexpected usize"),
    };

    pub fn of(row: usize, block: Block) Pivot {
        return .{ .row = @intCast(row), .block = block };
    }
};

/// `ipiv` of `n` pivots, as `View` reads the matrix: rows of the read
/// order, stored as the matrix's.
pub const Piv = struct {
    p: [*]Pivot,
    n: usize,
    flip: bool,

    fn map(pv: Piv, v: Pivot) Pivot {
        return if (pv.flip) .of(pv.n - 1 - v.row, v.block) else v;
    }

    pub fn get(pv: Piv, k: usize) Pivot {
        return pv.map(pv.p[pv.index(k)]);
    }

    pub fn set(pv: Piv, k: usize, v: Pivot) void {
        pv.p[pv.index(k)] = pv.map(v);
    }

    /// `k` (read order) as stored.
    pub fn index(pv: Piv, k: usize) usize {
        return if (pv.flip) pv.n - 1 - k else k;
    }
};

/// The element arithmetic of the `herm` or symmetric forms.
pub fn El(comptime T: type, comptime herm_: bool) type {
    return struct {
        pub const herm = herm_ and isComplex(T);
        pub const F = Float(T);

        /// A diagonal element's size: `|re|` where it is real (`herm`).
        pub fn absd(v: T) F {
            return if (comptime herm) @abs(v.re) else ops.abs1(T, v);
        }

        pub fn cj(v: T) T {
            return if (comptime herm) v.conjugate() else v;
        }

        /// `DBLE` where the diagonal is real.
        pub fn real(v: T) T {
            return if (comptime herm) .init(v.re, 0) else v;
        }

        /// The first NaN, else the first largest `abs1` of column `j`'s rows
        /// `[lo, hi)`, `hi > lo`.
        pub fn iamaxCol(a: View(T), j: usize, lo: usize, hi: usize) usize {
            var best = lo;
            var m: F = -1;
            for (lo..hi) |i| {
                const x = ops.abs1(T, a.get(i, j));
                if (std.math.isNan(x)) return i;
                if (x > m) {
                    m = x;
                    best = i;
                }
            }
            return best;
        }

        /// `iamaxCol` along row `i`'s columns `[lo, hi)`.
        pub fn iamaxRow(a: View(T), i: usize, lo: usize, hi: usize) usize {
            var best = lo;
            var m: F = -1;
            for (lo..hi) |j| {
                const x = ops.abs1(T, a.get(i, j));
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
        pub fn rank1(a: View(T), c: usize, k: usize, alpha: T) void {
            const n = a.rows;
            for (k..n) |j| {
                const xj = a.get(j, c);
                if (ops.isZero(T, xj)) {
                    if (comptime herm) a.set(j, j, real(a.get(j, j)));
                    continue;
                }
                const t = if (comptime herm) ops.mulReal(T, cj(xj), alpha.re) else mul(T, alpha, xj);
                if (comptime herm) {
                    const d = mul(T, xj, t);
                    a.set(j, j, .init(a.get(j, j).re + d.re, 0));
                } else a.set(j, j, add(T, a.get(j, j), mul(T, xj, t)));
                for (j + 1..n) |i| a.set(i, j, add(T, a.get(i, j), mul(T, a.get(i, c), t)));
            }
        }

        /// Rows and columns `kk < kp` of the lower triangle from `kk` on
        /// interchanged (`xSYTF2`'s or `xHETF2`'s swaps).
        pub fn swapSym(a: View(T), kk: usize, kp: usize) void {
            const n = a.rows;
            for (kp + 1..n) |i| mem.swap(T, a.at(i, kk), a.at(i, kp));
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
            } else mem.swap(T, a.at(kk, kk), a.at(kp, kp));
        }
    };
}

/// `alpha` of the pivot rule, `(1 + sqrt 17) / 8`.
pub fn alphaOf(comptime F: type) F {
    return (1 + @sqrt(@as(F, 17))) / 8;
}

// ---- tf2 --------------------------------------------------------------------------

/// `xSYTF2`, `xHETF2` (`herm`) or their `_ROOK` forms on `uplo`'s triangle
/// of `n x n` `a`, `ipiv` its `n` pivots: the first exactly zero (or NaN,
/// outside `rook`) `D(k, k)` (0-based, as stored), or null.
pub fn tf2(comptime T: type, comptime herm: bool, comptime rook: bool, upper: bool, n: usize, a_: [*]T, lda: usize, ipiv: [*]Pivot) ?usize {
    const E = El(T, herm);
    const F = Float(T);
    const a: View(T) = .{ .p = a_, .ld = lda, .rows = n, .cols = n, .flip = upper };
    const pv: Piv = .{ .p = ipiv, .n = n, .flip = upper };
    const alpha = alphaOf(F);
    const sfmin = std.math.floatMin(F);
    var info: ?usize = null;
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
            colmax = ops.abs1(T, a.get(imax, k));
        }
        if (@max(absakk, colmax) == 0 or (!rook and std.math.isNan(absakk))) {
            if (info == null) info = pv.index(k);
            if (comptime E.herm) a.set(k, k, E.real(a.get(k, k)));
        } else {
            if (!rook) {
                if (!(absakk >= alpha * colmax)) {
                    var jmax = E.iamaxRow(a, imax, k, imax);
                    var rowmax = ops.abs1(T, a.get(imax, jmax));
                    if (imax + 1 < n) {
                        jmax = E.iamaxCol(a, imax, imax + 1, n);
                        rowmax = @max(rowmax, ops.abs1(T, a.get(jmax, imax)));
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
                    rowmax = ops.abs1(T, a.get(imax, jmax));
                }
                if (imax + 1 < n) {
                    const itemp = E.iamaxCol(a, imax, imax + 1, n);
                    const dtemp = ops.abs1(T, a.get(itemp, imax));
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
                    mem.swap(T, a.at(k + 1, k), a.at(kp, k));
                }
            } else if (comptime E.herm) {
                a.set(k, k, E.real(a.get(k, k)));
                if (kstep == 2) a.set(k + 1, k + 1, E.real(a.get(k + 1, k + 1)));
            }
            if (kstep == 1) {
                if (k + 1 < n) oneByOne(T, herm, rook, a, k, sfmin);
            } else if (k + 2 < n) {
                twoByTwo(T, herm, a, k);
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
        for (k + 1..n) |i| a.set(i, k, if (comptime E.herm) ops.divReal(T, a.get(i, k), d.re) else ops.div(T, a.get(i, k), d));
        E.rank1(a, k, k + 1, ops.neg(T, d));
        return;
    }
    const r = ops.div(T, scalar(T, 1), E.real(akk));
    E.rank1(a, k, k + 1, ops.neg(T, r));
    for (k + 1..n) |i| a.set(i, k, if (comptime E.herm) ops.mulReal(T, a.get(i, k), r.re) else mul(T, r, a.get(i, k)));
}

/// `tf2`'s 2 x 2 step: `L(k+2:, k:k+1)` and the rank-2 update of the rest.
/// Each entry is divided by `d21`, then scaled by `t`, in every form
/// (`dsytf2.f`, `zhetf2.f` and the `_rook` forms; OpenBLAS `ce4ece13e`):
/// `t / d21` overflows for a tiny `d21`.
fn twoByTwo(comptime T: type, comptime herm: bool, a: View(T), k: usize) void {
    const E = El(T, herm);
    const n = a.rows;
    const one = scalar(T, 1);
    if (comptime E.herm) {
        const sub21 = a.get(k + 1, k);
        const dd = std.math.hypot(sub21.re, sub21.im);
        const d11 = a.get(k + 1, k + 1).re / dd;
        const d22 = a.get(k, k).re / dd;
        const tt = 1 / (d11 * d22 - 1);
        const d21 = ops.divReal(T, sub21, dd);
        for (k + 2..n) |j| {
            const wk = ops.mulReal(T, ops.divReal(T, sub(T, ops.mulReal(T, a.get(j, k), d11), mul(T, d21, a.get(j, k + 1))), dd), tt);
            const wkp1 = ops.mulReal(T, ops.divReal(T, sub(T, ops.mulReal(T, a.get(j, k + 1), d22), mul(T, ops.conj(T, d21), a.get(j, k))), dd), tt);
            for (j..n) |i| {
                a.set(i, j, sub(T, sub(T, a.get(i, j), mul(T, a.get(i, k), ops.conj(T, wk))), mul(T, a.get(i, k + 1), ops.conj(T, wkp1))));
            }
            a.set(j, k, wk);
            a.set(j, k + 1, wkp1);
            a.set(j, j, E.real(a.get(j, j)));
        }
        return;
    }
    const d21 = a.get(k + 1, k);
    const d11 = ops.div(T, a.get(k + 1, k + 1), d21);
    const d22 = ops.div(T, a.get(k, k), d21);
    const t = ops.div(T, one, sub(T, mul(T, d11, d22), one));
    for (k + 2..n) |j| {
        const wk = mul(T, t, ops.div(T, sub(T, mul(T, d11, a.get(j, k)), a.get(j, k + 1)), d21));
        const wkp1 = mul(T, t, ops.div(T, sub(T, mul(T, d22, a.get(j, k + 1)), a.get(j, k)), d21));
        for (j..n) |i| {
            a.set(i, j, sub(T, sub(T, a.get(i, j), mul(T, a.get(i, k), wk)), mul(T, a.get(i, k + 1), wkp1)));
        }
        a.set(j, k, wk);
        a.set(j, k + 1, wkp1);
    }
}
