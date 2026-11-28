//! `?SYTRI` and `?HETRI`: the inverse from `sytrf`'s or `hetrf`'s
//! factor, Reference-LAPACK 3.12.0 (`lapack-netlib/SRC/dsytri.f`,
//! `zsytri.f`, `zhetri.f`, OpenBLAS v0.3.34-794; modified BSD). One body
//! serves all three.

const std = @import("std");
const types = @import("../blas/types.zig");
const interface = @import("interface.zig");
const ops = @import("ops.zig");
const blas_calls = @import("blas_calls.zig");
const indefinite = @import("indefinite.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Uplo = types.Uplo;
const Mat = ops.Mat;
const Pivot = indefinite.Pivot;

/// What `sytri` and `hetri` return: `Singular` for a 1 x 1 block of `D`
/// exactly zero (the reference's `INFO > 0`), `A` untouched; a reversed
/// `uplo`, a slice too short (`a`, or `ipiv` under `n`), a leading
/// dimension too small or a pivot that is not a row of `A`; and the
/// allocator's and cancellation.
pub const SytriError = interface.TriangleError || interface.PivotError || interface.SingularError || Allocator.Error || Io.Cancelable;

/// `?SYTRI` (`dsytri.f`, `zsytri.f`): `A⁻¹` of the symmetric (complex:
/// complex symmetric) `A`, over `sytrf`'s factor in `a`'s `uplo` triangle,
/// given its `ipiv`.
pub fn sytri(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, a: []T, lda: usize, ipiv: []const Pivot) SytriError!void {
    comptime ops.requireFloat(T, "sytri");
    return invert(allocator, io, T, false, uplo, n, a, lda, ipiv);
}

/// `?HETRI` (`zhetri.f`): `sytri` for the Hermitian `A` from `hetrf`'s
/// factor; the diagonal comes out real. Complex only.
pub fn hetri(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, a: []T, lda: usize, ipiv: []const Pivot) SytriError!void {
    comptime ops.requireComplex(T, "hetri");
    return invert(allocator, io, T, true, uplo, n, a, lda, ipiv);
}

fn invert(allocator: Allocator, io: Io, comptime T: type, comptime herm_: bool, uplo: Uplo, n: usize, a: []T, lda: usize, ipiv: []const Pivot) SytriError!void {
    const E = indefinite.El(T, herm_);
    try interface.holdsTriangle(uplo, n, a.len, lda);
    if (ipiv.len < n) return error.SliceTooShort;
    for (ipiv[0..n]) |p| if (p.row >= n) return error.PivotOutOfRange;
    if (n == 0) return;
    const m: Mat(T) = .{ .ptr = a.ptr, .rows = n, .cols = n, .ld = lda };
    for (0..n) |i| if (ipiv[i].block == .one and ops.isZero(T, m.get(i, i))) return error.Singular;
    const work = try allocator.alloc(T, n);
    defer allocator.free(work);
    const s: Step(T, E) = .{ .al = allocator, .io = io, .m = m, .work = work, .upper = uplo == .u };
    if (s.upper) {
        var k: usize = 0;
        while (k < n) {
            const two = ipiv[k].block == .two;
            // the block's rest is above it: `A(0..k, 0..k)`
            if (two) try s.block2(k, k + 1, k, 0, k) else try s.block1(k, 0, k);
            const kp = ipiv[k].row;
            if (kp != k) {
                swapCols(T, m, 0, kp, k, kp);
                for (kp + 1..k) |j| {
                    const t = E.cj(m.get(j, k));
                    m.at(j, k).* = E.cj(m.get(kp, j));
                    m.at(kp, j).* = t;
                }
                s.finishSwap(k, kp, if (two) k + 1 else null);
            }
            k += if (two) 2 else 1;
        }
        return;
    }
    var k = n;
    while (k > 0) {
        k -= 1;
        const two = ipiv[k].block == .two;
        // the block's rest is below it: `A(k + 1..n, k + 1..n)`
        if (two) try s.block2(k - 1, k, k, k + 1, n) else try s.block1(k, k + 1, n);
        const kp = ipiv[k].row;
        if (kp != k) {
            swapCols(T, m, kp + 1, n, k, kp);
            for (k + 1..kp) |j| {
                const t = E.cj(m.get(j, k));
                m.at(j, k).* = E.cj(m.get(kp, j));
                m.at(kp, j).* = t;
            }
            s.finishSwap(k, kp, if (two) k - 1 else null);
        }
        if (two) k -= 1;
    }
}

/// Rows `lo..hi` of columns `p` and `q` swapped.
fn swapCols(comptime T: type, m: Mat(T), lo: usize, hi: usize, p: usize, q: usize) void {
    for (lo..hi) |i| std.mem.swap(T, m.at(i, p), m.at(i, q));
}

/// One block's work, its rest `A(lo..hi, lo..hi)` in the factor's triangle.
fn Step(comptime T: type, comptime E: type) type {
    return struct {
        const Self = @This();

        al: Allocator,
        io: Io,
        m: Mat(T),
        work: []T,
        upper: bool,

        /// `x := -R x` for column `j`'s rows `lo..hi`, `R` the rest; returns
        /// `x`'s old value dotted with its new one.
        fn apply(s: Self, j: usize, lo: usize, hi: usize) ops.Error!T {
            const x = s.m.col(j) + lo;
            const r = hi - lo;
            @memcpy(s.work[0..r], x[0..r]);
            try blas_calls.symvAs(s.al, s.io, T, E.herm, s.upper, ops.scalar(T, -1), s.m.sub(lo, lo, r, r), s.work.ptr, x);
            return blas_calls.dot(s.io, T, E.herm, r, s.work.ptr, 1, x, 1);
        }

        /// `D(k, k) -= d`, its real part only where `D` is Hermitian.
        fn subDiag(s: Self, k: usize, d: T) void {
            s.m.at(k, k).* = ops.sub(T, s.m.get(k, k), E.real(d));
        }

        fn block1(s: Self, k: usize, lo: usize, hi: usize) ops.Error!void {
            const d = s.m.get(k, k);
            s.m.at(k, k).* = if (comptime E.herm) ops.scalar(T, 1 / d.re) else ops.recip(T, d);
            if (hi > lo) s.subDiag(k, try s.apply(k, lo, hi));
        }

        /// The 2 x 2 block on `p < q`, its off-diagonal at `(e, ...)`:
        /// `A(p, q)` upper, `A(q, p)` lower, `e` its row.
        fn block2(s: Self, p: usize, q: usize, e: usize, lo: usize, hi: usize) ops.Error!void {
            const ec = if (e == p) q else p;
            const off = s.m.at(e, ec);
            if (comptime E.herm or !types.isComplex(T)) {
                const F = types.Float(T);
                const t: F = ops.abs(T, off.*);
                const ak = ops.re(T, s.m.get(p, p)) / t;
                const akp1 = ops.re(T, s.m.get(q, q)) / t;
                const akkp1 = ops.divReal(T, off.*, t);
                const d = t * (ak * akp1 - 1);
                s.m.at(p, p).* = ops.scalar(T, akp1 / d);
                s.m.at(q, q).* = ops.scalar(T, ak / d);
                off.* = ops.neg(T, ops.divReal(T, akkp1, d));
            } else {
                const t = off.*;
                const ak = ops.div(T, s.m.get(p, p), t);
                const akp1 = ops.div(T, s.m.get(q, q), t);
                const akkp1 = ops.div(T, off.*, t);
                const d = ops.mul(T, t, ops.sub(T, ops.mul(T, ak, akp1), ops.one(T)));
                s.m.at(p, p).* = ops.div(T, akp1, d);
                s.m.at(q, q).* = ops.div(T, ak, d);
                off.* = ops.neg(T, ops.div(T, akkp1, d));
            }
            if (hi == lo) return;
            // the reference's order: the column nearer the rest first
            const near = if (s.upper) p else q;
            const far = if (s.upper) q else p;
            s.subDiag(near, try s.apply(near, lo, hi));
            const r = hi - lo;
            off.* = ops.sub(T, off.*, try blas_calls.dot(s.io, T, E.herm, r, s.m.col(near) + lo, 1, s.m.col(far) + lo, 1));
            s.subDiag(far, try s.apply(far, lo, hi));
        }

        /// The swap's diagonal and, for a 2 x 2 block, rows `k` and `kp` of
        /// its other column.
        fn finishSwap(s: Self, k: usize, kp: usize, other: ?usize) void {
            if (comptime E.herm) s.m.at(kp, k).* = s.m.get(kp, k).conjugate();
            std.mem.swap(T, s.m.at(k, k), s.m.at(kp, kp));
            if (other) |o| std.mem.swap(T, s.m.at(k, o), s.m.at(kp, o));
        }
    };
}
