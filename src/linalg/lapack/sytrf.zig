//! `sytrf`, `hetrf` and their `_rook` forms: Reference-LAPACK's blocked
//! `xSYTRF` loop (`lapack-netlib/SRC/dsytrf.f` and the rest, OpenBLAS
//! v0.3.34-317) over `indefinite.zig`'s kernels, the trailing update by
//! ReLAPACK's `gemmt` (`relapack/src/dgemmt.c`, MIT).

const std = @import("std");
const types = @import("../blas/types.zig");
const driver = @import("../blas/blas3/driver.zig");
const ops = @import("ops.zig");
const interface = @import("interface.zig");
const lu = @import("lu.zig");
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

/// What `sytrf`, `hetrf` and their `_rook` forms return: `Singular` when a
/// diagonal block of `D` is exactly zero.
pub const SytrfError = interface.SquareError || error{Singular} || Allocator.Error || Io.Cancelable;

/// `?SYTRF`: `A = U D Uᵀ` (`u`) or `L D Lᵀ` (`l`) of the symmetric (complex:
/// complex symmetric) `n x n` `A` in its `uplo` triangle, by Bunch and
/// Kaufman's diagonal pivoting. `D`, the multipliers and `ipiv` are LAPACK's:
/// `ipiv[k]` positive for a 1 x 1 block whose rows `k + 1` and `ipiv[k]`
/// (1-based) were swapped, equal and negative on both rows of a 2 x 2 block.
/// A zero block of `D` (or a NaN pivot in `tf2`) still completes the
/// factorization, then returns `error.Singular` with the first met block's
/// first index (1-based; the lowest for `l`, the highest for `u`) in `pivot`.
/// Reference-LAPACK's blocked `xSYTRF` (`dsytrf.f`, `zsytrf.f`, `zhetrf.f`,
/// the `_rook` forms) over ReLAPACK's `rec2` panels (issues #242, #243).
pub fn sytrf(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, a: [*]T, lda: isize, ipiv: [*]isize, pivot: ?*usize) SytrfError!void {
    comptime ops.requireFloat(T, "sytrf");
    return indefiniteImpl(allocator, io, T, false, false, uplo, n, a, lda, ipiv, pivot);
}

/// `?HETRF` (`zhetrf.f`, `chetrf.f`; ReLAPACK's `zhetrf_rec2.c` panels):
/// `sytrf` for the Hermitian `A`, `A = U D Uᴴ` or `L D Lᴴ`; the diagonal's
/// imaginary parts are never read and `D`'s diagonal comes out real.
/// Complex only.
pub fn hetrf(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, a: [*]T, lda: isize, ipiv: [*]isize, pivot: ?*usize) SytrfError!void {
    comptime ops.requireComplex(T, "hetrf");
    return indefiniteImpl(allocator, io, T, true, false, uplo, n, a, lda, ipiv, pivot);
}

/// `?SYTRF_ROOK` (`dsytrf_rook.f` and the rest; ReLAPACK's
/// `dsytrf_rook_rec2.c` panels): `sytrf` by the bounded Bunch-Kaufman
/// ("rook") pivoting, whose `ipiv` names both interchanges of a 2 x 2 block
/// (`-ipiv[k]` and `-ipiv[k + 1]` for `l`). No NaN test: only an exactly
/// zero block is singular, as the reference's.
pub fn sytrfRook(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, a: [*]T, lda: isize, ipiv: [*]isize, pivot: ?*usize) SytrfError!void {
    comptime ops.requireFloat(T, "sytrfRook");
    return indefiniteImpl(allocator, io, T, false, true, uplo, n, a, lda, ipiv, pivot);
}

/// `?HETRF_ROOK` (`zhetrf_rook.f`, `chetrf_rook.f`; ReLAPACK's
/// `zhetrf_rook_rec2.c` panels): `sytrfRook` for the Hermitian `A`.
/// Complex only.
pub fn hetrfRook(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, n: isize, a: [*]T, lda: isize, ipiv: [*]isize, pivot: ?*usize) SytrfError!void {
    comptime ops.requireComplex(T, "hetrfRook");
    return indefiniteImpl(allocator, io, T, true, true, uplo, n, a, lda, ipiv, pivot);
}

/// `ILAENV`'s block for `xSYTRF` and `xHETRF` and their `_ROOK` forms
/// (`lapack-netlib/SRC/ilaenv.f`); at or past `n`, `tf2` takes all of `A`.
pub const sytrf_nb = 64;

fn indefiniteImpl(al: Allocator, io: Io, comptime T: type, comptime herm: bool, comptime rook: bool, uplo: Uplo, n: isize, a: [*]T, lda: isize, ipiv: [*]isize, pivot: ?*usize) SytrfError!void {
    if (!try interface.checkSquare(T, uplo, n, lda)) return;
    const un: usize = @intCast(n);
    const blocked = sytrf_nb < un;
    const work = try al.alloc(T, if (blocked) un * sytrf_nb else 0);
    defer al.free(work);
    const info = try sytrfBlocked(al, io, T, herm, rook, uplo == .u, un, a, @intCast(lda), ipiv, work.ptr, un);
    return lu.singular(info, pivot);
}

/// `xSYTRF`'s loop: panels from the first column (upper: the last), each
/// `panel`, the update of the rest, `undoPanel`; the reference's `INFO`.
fn sytrfBlocked(al: Allocator, io: Io, comptime T: type, comptime herm: bool, comptime rook: bool, upper: bool, n: usize, a: [*]T, lda: usize, ipiv: [*]isize, w: [*]T, ldw: usize) Error!usize {
    const one = scalar(T, 1);
    const minus = scalar(T, -1);
    const nb = sytrf_nb;
    var info: usize = 0;
    if (!upper) {
        var k: usize = 0;
        while (k < n) {
            const m = n - k;
            const d = a + k * (lda + 1);
            const p = ipiv + k;
            var kb = m;
            var iinfo: usize = 0;
            if (k + nb < n) {
                const r = panel(T, herm, rook, false, m, nb, d, lda, p, w, ldw);
                kb = r.kb;
                iinfo = r.info;
                // A22 -= L21 W21ᵀ (`xLASYF`), then L21 in LAPACK's form
                try gemmtRec(al, io, T, .l, m - kb, kb, minus, d + kb, lda, w + kb, ldw, one, d + kb * (lda + 1), lda);
                undoPanel(T, rook, false, m, d, lda, p, kb);
            } else iinfo = tf2(T, herm, rook, false, m, d, lda, p);
            if (info == 0 and iinfo != 0) info = iinfo + k;
            const off: isize = @intCast(k);
            for (p[0..kb]) |*v| v.* += if (v.* > 0) off else -off;
            k += kb;
        }
        return info;
    }
    var k = n;
    while (k > 0) {
        var kb = k;
        var iinfo: usize = 0;
        if (k > nb) {
            const r = panel(T, herm, rook, true, k, nb, a, lda, ipiv, w, ldw);
            kb = r.kb;
            iinfo = r.info;
            // A11 -= U12 W12ᵀ (`xLASYF`), W's last `kb` columns
            try gemmtRec(al, io, T, .u, k - kb, kb, minus, a + (k - kb) * lda, lda, w + (nb - kb) * ldw, ldw, one, a, lda);
            undoPanel(T, rook, true, k, a, lda, ipiv, kb);
        } else iinfo = tf2(T, herm, rook, true, k, a, lda, ipiv);
        if (info == 0 and iinfo != 0) info = iinfo;
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
pub fn indefiniteRoutine(comptime T: type, comptime herm: bool, comptime rook: bool) fn (Allocator, Io, comptime type, Uplo, isize, [*]T, isize, [*]isize, ?*usize) SytrfError!void {
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
        f(al, check.io, T, lcheck.uploOf(upper), @intCast(n), p.a.ptr, @intCast(p.lda), p.ipiv.ptr, null) catch |e| {
            check.say("{s} {s} n={d}\n", .{ what, @typeName(T), n });
            return e;
        };
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
                try f(al, check.io, T, lcheck.uploOf(upper), 150, p.a.ptr, @intCast(p.lda), p.ipiv.ptr, null);
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
                var piv: usize = 0;
                try testing.expectError(error.Singular, f(al, check.io, T, lcheck.uploOf(upper), n, p.a.ptr, @intCast(p.lda), p.ipiv.ptr, &piv));
                // the zero row may move with the interchanges (the reference's rook
                // search moves it too): `D` is zero where `pivot` says
                try testing.expect(ops.isZero(T, p.a[(piv - 1) * (p.lda + 1)]));
                try p.expectFactor(rook, al, "singular");
            };
            // all zero: the first column met, the last for `u`
            for ([_]bool{ true, false }) |upper| for ([_]usize{ 3, 70 }) |n| {
                var p = try symmetric.Indefinite(T, herm).init(al, n, upper, false, 1);
                defer p.deinit(al);
                @memset(p.w, reference.wide(T, scalar(T, 0)));
                p.store();
                var piv: usize = 0;
                try testing.expectError(error.Singular, f(al, check.io, T, lcheck.uploOf(upper), @intCast(n), p.a.ptr, @intCast(p.lda), p.ipiv.ptr, &piv));
                try testing.expectEqual(if (upper) n else 1, piv);
            };
        }
    };
}

test "relapack sytrf, hetrf: the first bad argument; nothing read for n = 0" {
    const al = testing.allocator;
    const C = types.Complex64;
    var a: [4]C = @splat(lcheck.poisoned(C));
    var ipiv: [2]isize = @splat(-7);
    try testing.expectError(error.InvalidUplo, sytrf(al, check.io, C, .v, 2, &a, 2, &ipiv, null));
    try testing.expectError(error.InvalidN, hetrf(al, check.io, C, .u, -1, &a, 2, &ipiv, null));
    try testing.expectError(error.InvalidLda, hetrfRook(al, check.io, C, .l, 2, &a, 1, &ipiv, null));
    try sytrfRook(al, check.io, C, .l, 0, &a, 1, &ipiv, null);
    for (a) |x| try testing.expect(lcheck.isNanT(C, x));
    try testing.expectEqual(-7, ipiv[0]);
}

const BigSym = struct {
    a: []const f64,

    fn run(self: BigSym, io: Io, out: []f64) anyerror!void {
        @memcpy(out, self.a);
        var ipiv: [lcheck.big_n]isize = undefined;
        try sytrf(testing.allocator, io, f64, .l, lcheck.big_n, out.ptr, lcheck.big_n, &ipiv, null);
        @memcpy(out, self.a);
        try sytrfRook(testing.allocator, io, f64, .u, lcheck.big_n, out.ptr, lcheck.big_n, &ipiv, null);
    }
};

test "relapack sytrf: the same bits inline and threaded" {
    const al = testing.allocator;
    const a = try al.alloc(f64, lcheck.big_n * lcheck.big_n);
    defer al.free(a);
    check.fill(f64, a, 4);
    try lcheck.expectIoBits(f64, BigSym{ .a = a }, BigSym.run, lcheck.big_n * lcheck.big_n, "sytrf");
}
