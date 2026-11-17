//! `sytrf`, `hetrf` and their `_rook` forms: Reference-LAPACK's blocked
//! `xSYTRF` loop (`lapack-netlib/SRC/dsytrf.f` and the rest, OpenBLAS
//! v0.3.34-317) over `indefinite.zig`'s kernels, the trailing update by
//! ReLAPACK's `gemmt` (`relapack/src/dgemmt.c`, MIT).

const std = @import("std");
const types = @import("../blas/types.zig");
const driver = @import("../blas/blas3/driver.zig");
const ops = @import("ops.zig");
const interface = @import("interface.zig");
const relapack = @import("relapack.zig");
const indefinite = @import("indefinite.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Uplo = types.Uplo;
const isComplex = types.isComplex;
const scalar = ops.scalar;
const tf2 = indefinite.tf2;
const panel = indefinite.panel;
const undoPanel = indefinite.undoPanel;

const Error = driver.Error;

/// One entry of the pivots (`indefinite.zig`).
pub const Pivot = indefinite.Pivot;

/// What `sytrf`, `hetrf` and their `_rook` forms return: a reversed `uplo`,
/// a slice too short (`a`, or `ipiv` under `n`) or a leading dimension too
/// small, and the allocator's and cancellation.
pub const SytrfError = interface.TriangleError || Allocator.Error || Io.Cancelable;

/// `?SYTRF`: `A = U D Uᵀ` (`u`) or `L D Lᵀ` (`l`) of the symmetric (complex:
/// complex symmetric) `n x n` `A` in its `uplo` triangle, by Bunch and
/// Kaufman's diagonal pivoting. `D` and the multipliers are LAPACK's; the
/// first `n` entries of `ipiv` are 0-based (issue #314): `ipiv[k]` is `.one`
/// for a 1 x 1 block whose rows `k` and `ipiv[k].row` were swapped, `.two`
/// on both rows of a 2 x 2 block, each naming the row swapped there.
///
/// Returns the first row (0-based; the lowest for `l`, the highest for `u`)
/// of the first block of `D` met exactly zero (or a NaN pivot in `tf2`), or
/// null. A singular `D` is a result, not an error, as `getrf`'s `U`: the
/// factorization completes.
/// Reference-LAPACK's blocked `xSYTRF` (`dsytrf.f`, `zsytrf.f`, `zhetrf.f`,
/// the `_rook` forms) over ReLAPACK's `rec2` panels (issues #242, #243).
pub fn sytrf(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, a: []T, lda: usize, ipiv: []Pivot) SytrfError!?usize {
    comptime ops.requireFloat(T, "sytrf");
    return indefiniteImpl(allocator, io, T, false, false, uplo, n, a, lda, ipiv);
}

/// `?HETRF` (`zhetrf.f`, `chetrf.f`; ReLAPACK's `zhetrf_rec2.c` panels):
/// `sytrf` for the Hermitian `A`, `A = U D Uᴴ` or `L D Lᴴ`; the diagonal's
/// imaginary parts are never read and `D`'s diagonal comes out real.
/// Complex only.
pub fn hetrf(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, a: []T, lda: usize, ipiv: []Pivot) SytrfError!?usize {
    comptime ops.requireComplex(T, "hetrf");
    return indefiniteImpl(allocator, io, T, true, false, uplo, n, a, lda, ipiv);
}

/// `?SYTRF_ROOK` (`dsytrf_rook.f` and the rest; ReLAPACK's
/// `dsytrf_rook_rec2.c` panels): `sytrf` by the bounded Bunch-Kaufman
/// ("rook") pivoting, which swaps twice for a 2 x 2 block: its two entries
/// name different rows. No NaN test: only an exactly zero block is
/// singular, as the reference's.
pub fn sytrfRook(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, a: []T, lda: usize, ipiv: []Pivot) SytrfError!?usize {
    comptime ops.requireFloat(T, "sytrfRook");
    return indefiniteImpl(allocator, io, T, false, true, uplo, n, a, lda, ipiv);
}

/// `?HETRF_ROOK` (`zhetrf_rook.f`, `chetrf_rook.f`; ReLAPACK's
/// `zhetrf_rook_rec2.c` panels): `sytrfRook` for the Hermitian `A`.
/// Complex only.
pub fn hetrfRook(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, a: []T, lda: usize, ipiv: []Pivot) SytrfError!?usize {
    comptime ops.requireComplex(T, "hetrfRook");
    return indefiniteImpl(allocator, io, T, true, true, uplo, n, a, lda, ipiv);
}

/// `ILAENV`'s block for `xSYTRF` and `xHETRF` and their `_ROOK` forms
/// (`lapack-netlib/SRC/ilaenv.f`); at or past `n`, `tf2` takes all of `A`.
pub const sytrf_nb = 64;

fn indefiniteImpl(al: Allocator, io: Io, comptime T: type, comptime herm: bool, comptime rook: bool, uplo: Uplo, n: usize, a: []T, lda: usize, ipiv: []Pivot) SytrfError!?usize {
    try interface.holdsTriangle(uplo, n, a.len, lda);
    if (ipiv.len < n) return error.SliceTooShort;
    if (n == 0) return null;
    const work = try al.alloc(T, if (sytrf_nb < n) n * sytrf_nb else 0);
    defer al.free(work);
    return sytrfBlocked(al, io, T, herm, rook, uplo == .u, n, a.ptr, lda, ipiv.ptr, work.ptr, n);
}

/// `xSYTRF`'s loop: panels from the first column (upper: the last), each
/// `panel`, the update of the rest, `undoPanel`; the first zero block of `D`.
fn sytrfBlocked(al: Allocator, io: Io, comptime T: type, comptime herm: bool, comptime rook: bool, upper: bool, n: usize, a: [*]T, lda: usize, ipiv: [*]Pivot, w: [*]T, ldw: usize) Error!?usize {
    const one = scalar(T, 1);
    const minus = scalar(T, -1);
    const nb = sytrf_nb;
    var info: ?usize = null;
    if (!upper) {
        var k: usize = 0;
        while (k < n) {
            const m = n - k;
            const d = a + k * (lda + 1);
            const p = ipiv + k;
            var kb = m;
            var iinfo: ?usize = null;
            if (k + nb < n) {
                const r = panel(T, herm, rook, false, m, nb, d, lda, p, w, ldw);
                kb = r.kb;
                iinfo = r.info;
                // A22 -= L21 W21ᵀ (`xLASYF`), then L21 in LAPACK's form
                try gemmtRec(al, io, T, .l, m - kb, kb, minus, d + kb, lda, w + kb, ldw, one, d + kb * (lda + 1), lda);
                undoPanel(T, rook, false, m, d, lda, p, kb);
            } else iinfo = tf2(T, herm, rook, false, m, d, lda, p);
            info = interface.firstOf(info, iinfo, k);
            for (p[0..kb]) |*v| v.row += @intCast(k);
            k += kb;
        }
        return info;
    }
    var k = n;
    while (k > 0) {
        var kb = k;
        var iinfo: ?usize = null;
        if (k > nb) {
            const r = panel(T, herm, rook, true, k, nb, a, lda, ipiv, w, ldw);
            kb = r.kb;
            iinfo = r.info;
            // A11 -= U12 W12ᵀ (`xLASYF`), W's last `kb` columns
            try gemmtRec(al, io, T, .u, k - kb, kb, minus, a + (k - kb) * lda, lda, w + (nb - kb) * ldw, ldw, one, a, lda);
            undoPanel(T, rook, true, k, a, lda, ipiv, kb);
        } else iinfo = tf2(T, herm, rook, true, k, a, lda, ipiv);
        info = info orelse iinfo;
        k -= kb;
    }
    return info;
}

/// `RELAPACK_?gemmt` as its `sytrf` calls it (`relapack/src/dgemmt.c` and
/// the rest, `HAVE_XGEMMT` 0): `C := alpha A Bᵀ + beta C` on the `uplo`
/// triangle of `n x n` `C`, `A` and `B` `n x k`; halves by `gemm` down to the
/// relapack.crossover (`CROSSOVER_XGEMMT`, `sytrf`'s), then `gemv` per column.
fn gemmtRec(al: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, k: usize, alpha: T, a: [*]const T, lda: usize, b: [*]const T, ldb: usize, beta: T, c: [*]T, ldc: usize) Error!void {
    if (n <= relapack.crossover) {
        // `RELAPACK_?gemmt_rec2`: `xGEMV` `N` per column, `beta` first
        for (0..n) |i| {
            const lo = if (uplo == .l) i else 0;
            const hi = if (uplo == .l) n else i + 1;
            const y = c + i * ldc;
            if (!ops.isZero(T, ops.sub(T, beta, scalar(T, 1)))) for (lo..hi) |r| {
                y[r] = if (ops.isZero(T, beta)) scalar(T, 0) else ops.mul(T, beta, y[r]);
            };
            for (0..k) |j| {
                const t = ops.mul(T, alpha, b[i + j * ldb]);
                for (lo..hi) |r| y[r] = ops.add(T, y[r], ops.mul(T, t, a[r + j * lda]));
            }
        }
        return;
    }
    const n1 = relapack.split(T, n);
    const n2 = n - n1;
    try gemmtRec(al, io, T, uplo, n1, k, alpha, a, lda, b, ldb, beta, c, ldc);
    if (uplo == .l) {
        try relapack.gemm(al, io, T, .n, .t, n2, n1, k, alpha, a + n1, lda, b, ldb, beta, c + n1, ldc);
    } else {
        try relapack.gemm(al, io, T, .n, .t, n1, n2, k, alpha, a, lda, b + n1, ldb, beta, c + n1 * ldc, ldc);
    }
    try gemmtRec(al, io, T, uplo, n2, k, alpha, a + n1, lda, b + n1, ldb, beta, c + n1 * ldc + n1, ldc);
}

// ---- tests ---------------------------------------------------------------------

const testing = std.testing;
const check = @import("../blas/_test/check.zig");
const reference = @import("../blas/_test/reference.zig");
const lcheck = @import("_test/check.zig");
const symmetric = @import("_test/symmetric.zig");

/// `sytrf`, `hetrf` (`herm`) or a `_rook` form, by its flags.
pub fn indefiniteRoutine(comptime T: type, comptime herm: bool, comptime rook: bool) fn (Allocator, Io, comptime type, Uplo, usize, []T, usize, []Pivot) SytrfError!?usize {
    return if (herm) (if (rook) hetrfRook else hetrf) else (if (rook) sytrfRook else sytrf);
}

/// Around the panel width (`sytrf_nb`), past it the blocked loop.
const sytrf_sizes = [_]usize{ 1, 2, 3, 24, 57, 65, 66, 100 };

fn sytrfCases(comptime T: type, comptime herm: bool, comptime rook: bool) !void {
    const al = testing.allocator;
    const f = comptime indefiniteRoutine(T, herm, rook);
    const what = if (herm) (if (rook) "hetrfRook" else "hetrf") else (if (rook) "sytrfRook" else "sytrf");
    for ([_]bool{ true, false }) |upper| for ([_]bool{ false, true }) |zd| for (sytrf_sizes) |n| {
        // a zero 1 x 1 is singular
        if ((check.rig() and n > 66) or (zd and n == 1)) continue;
        var p = try symmetric.Indefinite(T, herm).init(al, n, upper, zd, n * 7 + @intFromBool(upper) + 2 * @as(u64, @intFromBool(zd)));
        defer p.deinit(al);
        const zero = f(al, check.io, T, lcheck.uploOf(upper), n, p.a, p.lda, p.ipiv) catch |e| {
            check.say("{s} {s} n={d}\n", .{ what, @typeName(T), n });
            return e;
        };
        try testing.expectEqual(@as(?usize, null), zero);
        try p.expectFactor(rook, al, what);
    };
}

test "relapack sytrf, sytrfRook: A = L D Lᵀ and U D Uᵀ, random and zero-diagonal, every type" {
    inline for (lcheck.test_types) |T| inline for (.{ false, true }) |rook| try sytrfCases(T, false, rook);
}

test "relapack hetrf, hetrfRook: A = L D Lᴴ and U D Uᴴ, random and zero-diagonal, complex types" {
    inline for (lcheck.test_types) |T| if (comptime isComplex(T)) inline for (.{ false, true }) |rook| try sytrfCases(T, true, rook);
}

test "relapack sytrf and the rest: three panels, f64 and c64" {
    if (check.rig()) return error.SkipZigTest;
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64 }) |T| inline for (.{ false, true }) |herm| inline for (.{ false, true }) |rook| {
        if (comptime !herm or isComplex(T)) {
            const f = comptime indefiniteRoutine(T, herm, rook);
            for ([_]bool{ true, false }) |upper| {
                var p = try symmetric.Indefinite(T, herm).init(al, 150, upper, false, 99);
                defer p.deinit(al);
                _ = try f(al, check.io, T, lcheck.uploOf(upper), 150, p.a, p.lda, p.ipiv);
                try p.expectFactor(rook, al, "sytrf 150");
            }
        }
    };
}

test "relapack sytrf and the rest: a zero row and column is singular there, the factorization completed" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64 }) |T| inline for (.{ false, true }) |herm| inline for (.{ false, true }) |rook| {
        if (comptime !herm or isComplex(T)) {
            const f = comptime indefiniteRoutine(T, herm, rook);
            for ([_]bool{ true, false }) |upper| for ([_]usize{ 0, 5, 70, 99 }) |k| {
                const n = 100;
                var p = try symmetric.Indefinite(T, herm).init(al, n, upper, false, 3 + k);
                defer p.deinit(al);
                for (0..n) |i| {
                    p.w[i + k * n] = reference.wide(T, scalar(T, 0));
                    p.w[k + i * n] = reference.wide(T, scalar(T, 0));
                }
                p.store();
                const piv = (try f(al, check.io, T, lcheck.uploOf(upper), n, p.a, p.lda, p.ipiv)) orelse return error.TestExpectedSingular;
                // the zero row may move with the interchanges (the reference's rook
                // search moves it too): `D` is zero where the result says
                try testing.expect(ops.isZero(T, p.a[piv * (p.lda + 1)]));
                try p.expectFactor(rook, al, "singular");
            };
            // all zero: the first column met, the last for `u`
            for ([_]bool{ true, false }) |upper| for ([_]usize{ 3, 70 }) |n| {
                var p = try symmetric.Indefinite(T, herm).init(al, n, upper, false, 1);
                defer p.deinit(al);
                @memset(p.w, reference.wide(T, scalar(T, 0)));
                p.store();
                try testing.expectEqual(@as(?usize, if (upper) n - 1 else 0), try f(al, check.io, T, lcheck.uploOf(upper), n, p.a, p.lda, p.ipiv));
            };
        }
    };
}

test "relapack sytrf, hetrf: a reversed uplo, short slices, a small leading dimension; nothing read for n = 0" {
    const al = testing.allocator;
    const C = types.Complex64;
    var a: [4]C = @splat(lcheck.poisoned(C));
    const unset: Pivot = .of(7, .two);
    var ipiv: [2]Pivot = @splat(unset);
    try testing.expectError(error.UnsupportedUplo, sytrf(al, check.io, C, .v, 2, &a, 2, &ipiv));
    try testing.expectError(error.LeadingDimensionTooSmall, hetrfRook(al, check.io, C, .l, 2, &a, 1, &ipiv));
    try testing.expectError(error.SliceTooShort, hetrf(al, check.io, C, .u, 3, &a, 3, &ipiv));
    try testing.expectError(error.SliceTooShort, sytrf(al, check.io, C, .u, 2, &a, 2, ipiv[0..1]));
    try testing.expectEqual(@as(?usize, null), try sytrfRook(al, check.io, C, .l, 0, &a, 1, &ipiv));
    for (a) |x| try testing.expect(lcheck.isNanT(C, x));
    try testing.expectEqual(unset, ipiv[0]);
}

const BigSym = struct {
    a: []const f64,

    fn run(self: BigSym, io: Io, out: []f64) anyerror!void {
        @memcpy(out, self.a);
        var ipiv: [lcheck.big_n]Pivot = undefined;
        _ = try sytrf(testing.allocator, io, f64, .l, lcheck.big_n, out, lcheck.big_n, &ipiv);
        @memcpy(out, self.a);
        _ = try sytrfRook(testing.allocator, io, f64, .u, lcheck.big_n, out, lcheck.big_n, &ipiv);
    }
};

test "relapack sytrf: the same bits inline and threaded" {
    const al = testing.allocator;
    const a = try al.alloc(f64, lcheck.big_n * lcheck.big_n);
    defer al.free(a);
    check.fill(f64, a, 4);
    try lcheck.expectIoBits(f64, BigSym{ .a = a }, BigSym.run, lcheck.big_n * lcheck.big_n, "sytrf");
}
