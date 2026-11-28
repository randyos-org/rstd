//! `sytrs`, `hetrs` and their `_rook` forms: Reference-LAPACK's `xSYTRS`,
//! `xHETRS`, `xSYTRS_ROOK` and `xHETRS_ROOK` (`lapack-netlib/SRC/dsytrs.f`,
//! `zsytrs.f`, `zhetrs.f`, `dsytrs_rook.f` and the rest, OpenBLAS
//! v0.3.34-317), their `xGER` and `xGEMV` as loops. One body serves every
//! form; the upper triangle runs the lower algorithm on `A` read backwards
//! and `B`'s rows reversed, as `sytrf` factors it (`indefinite.zig`).

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const ops = @import("ops.zig");
const interface = @import("interface.zig");
const indefinite = @import("indefinite.zig");

const Uplo = types.Uplo;
const isComplex = types.isComplex;
const scalar = ops.scalar;
const sub = ops.sub;
const mul = ops.mul;
const div = ops.div;

/// One entry of the pivots, as `sytrf` gives them.
pub const Pivot = indefinite.Pivot;

/// What the solves return: a reversed `uplo`, a slice too short (`a`, `b`,
/// or `ipiv` under `n`), a leading dimension too small, or a pivot that is
/// not a row of `A`, checked before `B` is touched.
pub const SytrsError = interface.TriangleError || interface.PivotError;

/// `?SYTRS`: solves `A X = B` for the `n x nrhs` `X`, in place of `B`, with
/// `sytrf`'s `U D Uᵀ` or `L D Lᵀ` of `A` (complex: complex symmetric) and
/// its `ipiv`. No singularity check, as in the reference.
pub fn sytrs(comptime T: type, uplo: Uplo, n: usize, nrhs: usize, a: []const T, lda: usize, ipiv: []const Pivot, b: []T, ldb: usize) SytrsError!void {
    comptime ops.requireFloat(T, "sytrs");
    return checked(T, false, false, uplo, n, nrhs, a, lda, ipiv, b, ldb);
}

/// `?HETRS`: `sytrs` with `hetrf`'s `U D Uᴴ` or `L D Lᴴ`. Complex only.
pub fn hetrs(comptime T: type, uplo: Uplo, n: usize, nrhs: usize, a: []const T, lda: usize, ipiv: []const Pivot, b: []T, ldb: usize) SytrsError!void {
    comptime ops.requireComplex(T, "hetrs");
    return checked(T, true, false, uplo, n, nrhs, a, lda, ipiv, b, ldb);
}

/// `?SYTRS_ROOK`: `sytrs` with `sytrfRook`'s factor, whose 2 x 2 blocks
/// swap both their rows.
pub fn sytrsRook(comptime T: type, uplo: Uplo, n: usize, nrhs: usize, a: []const T, lda: usize, ipiv: []const Pivot, b: []T, ldb: usize) SytrsError!void {
    comptime ops.requireFloat(T, "sytrsRook");
    return checked(T, false, true, uplo, n, nrhs, a, lda, ipiv, b, ldb);
}

/// `?HETRS_ROOK`: `hetrs` with `hetrfRook`'s factor. Complex only.
pub fn hetrsRook(comptime T: type, uplo: Uplo, n: usize, nrhs: usize, a: []const T, lda: usize, ipiv: []const Pivot, b: []T, ldb: usize) SytrsError!void {
    comptime ops.requireComplex(T, "hetrsRook");
    return checked(T, true, true, uplo, n, nrhs, a, lda, ipiv, b, ldb);
}

fn checked(comptime T: type, comptime herm: bool, comptime rook: bool, uplo: Uplo, n: usize, nrhs: usize, a: []const T, lda: usize, ipiv: []const Pivot, b: []T, ldb: usize) SytrsError!void {
    try interface.holdsTriangle(uplo, n, a.len, lda);
    try util.holds(b.len, n, nrhs, ldb);
    if (ipiv.len < n) return error.SliceTooShort;
    for (ipiv[0..n]) |p| if (p.row >= n) return error.PivotOutOfRange;
    if (n == 0 or nrhs == 0) return;
    solve(T, herm, rook, uplo == .u, n, nrhs, a, lda, ipiv, b, ldb);
}

/// `A` (the factor), `B` and the pivots in the lower algorithm's read order.
fn Sys(comptime T: type) type {
    return struct {
        a: []const T,
        lda: usize,
        b: []T,
        ldb: usize,
        ipiv: []const Pivot,
        n: usize,
        flip: bool,

        const Self = @This();

        fn r(s: Self, i: usize) usize {
            return if (s.flip) s.n - 1 - i else i;
        }

        fn fa(s: Self, i: usize, j: usize) T {
            return s.a[s.r(i) + s.r(j) * s.lda];
        }

        fn fb(s: Self, i: usize, j: usize) *T {
            return &s.b[s.r(i) + j * s.ldb];
        }

        /// The row swapped with read row `k`.
        fn piv(s: Self, k: usize) Pivot {
            const p = s.ipiv[s.r(k)];
            return .of(s.r(p.row), p.block);
        }

        fn swap(s: Self, i: usize, k: usize, nrhs: usize) void {
            if (i == k) return;
            for (0..nrhs) |j| std.mem.swap(T, s.fb(i, j), s.fb(k, j));
        }
    };
}

fn solve(comptime T: type, comptime herm_: bool, comptime rook: bool, upper: bool, n: usize, nrhs: usize, a: []const T, lda: usize, ipiv: []const Pivot, b: []T, ldb: usize) void {
    const E = indefinite.El(T, herm_);
    const herm = E.herm;
    const s: Sys(T) = .{ .a = a, .lda = lda, .b = b, .ldb = ldb, .ipiv = ipiv, .n = n, .flip = upper };

    // `L D Y = B`
    var k: usize = 0;
    while (k < n) {
        if (s.piv(k).block == .one) {
            s.swap(k, s.piv(k).row, nrhs);
            subOuter(T, s, k, k + 1, nrhs);
            if (comptime herm) {
                const d = 1 / s.fa(k, k).re;
                for (0..nrhs) |j| s.fb(k, j).* = ops.mulReal(T, s.fb(k, j).*, d);
            } else {
                const d = div(T, scalar(T, 1), s.fa(k, k));
                for (0..nrhs) |j| s.fb(k, j).* = mul(T, d, s.fb(k, j).*);
            }
            k += 1;
        } else {
            if (rook) s.swap(k, s.piv(k).row, nrhs);
            s.swap(k + 1, s.piv(k + 1).row, nrhs);
            subOuter(T, s, k, k + 2, nrhs);
            subOuter(T, s, k + 1, k + 2, nrhs);
            const akm1k = s.fa(k + 1, k);
            const akm1 = div(T, s.fa(k, k), E.cj(akm1k));
            const ak = div(T, s.fa(k + 1, k + 1), akm1k);
            const denom = sub(T, mul(T, akm1, ak), scalar(T, 1));
            for (0..nrhs) |j| {
                const bkm1 = div(T, s.fb(k, j).*, E.cj(akm1k));
                const bk = div(T, s.fb(k + 1, j).*, akm1k);
                s.fb(k, j).* = div(T, sub(T, mul(T, ak, bkm1), bk), denom);
                s.fb(k + 1, j).* = div(T, sub(T, mul(T, akm1, bk), bkm1), denom);
            }
            k += 2;
        }
    }

    // `Lᵀ X = Y` (`Lᴴ`)
    k = n;
    while (k > 0) {
        const kk = k - 1;
        if (s.piv(kk).block == .one) {
            subDot(T, herm, s, kk, kk, nrhs);
            s.swap(kk, s.piv(kk).row, nrhs);
            k -= 1;
        } else {
            // `kk` is the block's second row
            subDot(T, herm, s, kk, kk, nrhs);
            subDot(T, herm, s, kk - 1, kk, nrhs);
            s.swap(kk, s.piv(kk).row, nrhs);
            if (rook) s.swap(kk - 1, s.piv(kk - 1).row, nrhs);
            k -= 2;
        }
    }
}

/// `xGER`: rows `[from, n)` of `B` less column `c` of `A` times row `c` of `B`.
fn subOuter(comptime T: type, s: Sys(T), c: usize, from: usize, nrhs: usize) void {
    for (0..nrhs) |j| {
        const x = s.fb(c, j).*;
        for (from..s.n) |i| s.fb(i, j).* = sub(T, s.fb(i, j).*, mul(T, s.fa(i, c), x));
    }
}

/// `xGEMV` (`ZLACGV` around a conjugate transpose, `herm`): row `k` of `B`
/// less column `k` of `A` below row `c`, transposed (conjugated), times
/// those rows of `B`.
fn subDot(comptime T: type, comptime herm: bool, s: Sys(T), k: usize, c: usize, nrhs: usize) void {
    for (0..nrhs) |j| {
        var t = scalar(T, 0);
        for (c + 1..s.n) |i| {
            const p = if (comptime herm) ops.mulConj(T, s.fa(i, k), s.fb(i, j).*) else mul(T, s.fa(i, k), s.fb(i, j).*);
            t = ops.add(T, t, p);
        }
        s.fb(k, j).* = sub(T, s.fb(k, j).*, t);
    }
}

// ---- tests ---------------------------------------------------------------------

const testing = std.testing;
const check = @import("../blas/_test/check.zig");
const reference = @import("../blas/_test/reference.zig");
const lcheck = @import("_test/check.zig");
const symmetric = @import("_test/symmetric.zig");
const sytrf_mod = @import("sytrf.zig");

fn sytrsCases(comptime T: type, comptime herm: bool, comptime rook: bool) !void {
    const al = testing.allocator;
    const factor = comptime sytrf_mod.indefiniteRoutine(T, herm, rook);
    const solveWith = if (herm) (if (rook) hetrsRook else hetrs) else (if (rook) sytrsRook else sytrs);
    const what = if (herm) (if (rook) "hetrsRook" else "hetrs") else (if (rook) "sytrsRook" else "sytrs");
    const nrhs = 3;
    const eps: f128 = std.math.floatEps(types.Float(T));
    for ([_]bool{ true, false }) |upper| for ([_]bool{ false, true }) |zd| for ([_]usize{ 1, 2, 3, 24, 57 }) |n| {
        if (zd and n == 1) continue;
        const seed = n * 11 + @intFromBool(upper) + 2 * @as(u64, @intFromBool(zd));
        var p = try symmetric.Indefinite(T, herm).init(al, n, upper, zd, seed);
        defer p.deinit(al);
        const ldb = n + 1;
        const b = try al.alloc(T, ldb * nrhs);
        defer al.free(b);
        check.fill(T, b, seed + 5);
        const b0 = try al.dupe(T, b);
        defer al.free(b0);
        const uplo = lcheck.uploOf(upper);
        try testing.expectEqual(@as(?usize, null), try factor(al, check.io, T, uplo, n, p.a, p.lda, p.ipiv));
        try solveWith(T, uplo, n, nrhs, p.a, p.lda, p.ipiv, b, ldb);
        // `A X - B` within `64 (n + 2)` ulps of `|A| |X| + |B|`
        const ulps: f128 = @floatFromInt(64 * (n + 2));
        for (0..nrhs) |j| for (0..n) |i| {
            var sum = reference.wide(T, scalar(T, 0));
            var mag: f128 = 0;
            for (0..n) |k| {
                const x = reference.wide(T, b[k + j * ldb]);
                sum = reference.addW(sum, reference.mulW(p.w[i + k * n], x));
                mag += reference.abs1(p.w[i + k * n]) * reference.abs1(x);
            }
            const want = reference.wide(T, b0[i + j * ldb]);
            const d = reference.abs1(reference.subW(sum, want));
            const allowed = ulps * eps * (mag + reference.abs1(want)) + std.math.floatMin(types.Float(T));
            if (!(d <= allowed)) {
                check.say("{s} {s} {s} zd={} n={d}: (A X - B)({d},{d}) is {e}, allowed {e}\n", .{ what, @typeName(T), if (upper) "u" else "l", zd, n, i, j, @as(f64, @floatCast(d)), @as(f64, @floatCast(allowed)) });
                return error.TestUnexpectedResult;
            }
        };
    };
}

test "sytrs, sytrsRook: A X = B with sytrf's factors, random and zero-diagonal, every type" {
    inline for (lcheck.test_types) |T| inline for (.{ false, true }) |rook| try sytrsCases(T, false, rook);
}

test "hetrs, hetrsRook: A X = B with hetrf's factors, complex types" {
    inline for (lcheck.test_types) |T| if (comptime isComplex(T)) inline for (.{ false, true }) |rook| try sytrsCases(T, true, rook);
}

test "sytrs: a pivot past n, short slices and a reversed uplo are errors, B untouched" {
    const a = [_]f64{ 2, 0, 0, 3 };
    var b = [_]f64{ 1, 1 };
    const bad = [_]Pivot{ .of(0, .one), .of(2, .one) };
    try testing.expectError(error.PivotOutOfRange, sytrs(f64, .l, 2, 1, &a, 2, &bad, &b, 2));
    try testing.expectError(error.SliceTooShort, sytrs(f64, .l, 2, 1, &a, 2, bad[0..1], &b, 2));
    try testing.expectEqualSlices(f64, &.{ 1, 1 }, &b);
    const good = [_]Pivot{ .of(0, .one), .of(1, .one) };
    try sytrs(f64, .l, 2, 1, &a, 2, &good, &b, 2);
    try testing.expectEqualSlices(f64, &.{ 0.5, 1.0 / 3.0 }, &b);
}
