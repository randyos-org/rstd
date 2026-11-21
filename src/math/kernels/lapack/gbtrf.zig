//! `gbtrf`, the band LU, by Reference-LAPACK's `xGBTF2`
//! (`lapack-netlib/SRC/dgbtf2.f`, `zgbtf2.f`; OpenBLAS v0.3.34-317) in place
//! of ReLAPACK's recursion (`gbtrf`'s doc).

const std = @import("std");
const mem = @import("../../../mem.zig");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const ops = @import("ops.zig");
const interface = @import("interface.zig");

const Allocator = std.mem.Allocator;
const Float = types.Float;
const isComplex = types.isComplex;
const scalar = ops.scalar;

/// What `gbtrf` returns: a slice too short (`ab`, or `ipiv` under
/// `min(m, n)`) or a leading dimension under `2 kl + ku + 1`.
pub const GbtrfError = interface.ShapeError;

/// `?GBTRF`: `A = P L U` of the `m x n` band matrix with `kl` sub- and `ku`
/// superdiagonals in LAPACK's band storage `ab` (`ldab >= 2 kl + ku + 1`,
/// `A(i, j)` at row `kl + ku + i - j` of column `j`; the first `kl` rows
/// take `U`'s fill-in), by partial pivoting; `ipiv` (0-based) and the
/// multipliers as LAPACK stores them, each interchange applied to the
/// columns right of its step only. Returns the first `i` (0-based) whose
/// `U(i, i)` is exactly zero, or null, the factorization completed either
/// way, as `getrf`.
///
/// ReLAPACK's (`relapack/src/dgbtrf.c` and the rest; OpenBLAS v0.3.34-317)
/// runs the reference's `xGBTF2` (`gbtf2`) unless `n <= kl`, where its
/// recursion gives a wrong factor, as upstream's full recursion does on most
/// shapes (issue #244); so `gbtrf` is `gbtf2` throughout, which allocates
/// nothing and runs on one thread.
pub fn gbtrf(comptime T: type, m: usize, n: usize, kl: usize, ku: usize, ab: []T, ldab: usize, ipiv: []usize) GbtrfError!?usize {
    comptime ops.requireFloat(T, "gbtrf");
    try util.holds(ab.len, 2 * kl + ku + 1, n, ldab);
    if (ipiv.len < @min(m, n)) return error.SliceTooShort;
    if (m == 0 or n == 0) return null;
    return gbtf2(T, m, n, kl, ku, ab.ptr, ldab, ipiv.ptr);
}

/// Reference-LAPACK's `xGBTF2` (`lapack-netlib/SRC/dgbtf2.f`, `zgbtf2.f`;
/// OpenBLAS `48c7b639c`): column by column, the pivot the first NaN or else
/// the first largest `|re| + |im|` (rstd's `iamax`), the interchange over the
/// columns the step reaches (`JU`), the multipliers times the pivot's
/// reciprocal or, under the safe minimum, divided by it (`ops.scaleBelow`),
/// and `xGER` (`xGERU`) in the reference's order, skipping a zero `y(j)`.
/// The first zero pivot (0-based), or null.
fn gbtf2(comptime T: type, m: usize, n: usize, kl: usize, ku: usize, ab: [*]T, ldab: usize, ipiv: [*]usize) ?usize {
    const F = Float(T);
    const kv = ku + kl;
    const zero = scalar(T, 0);
    // `AB(r, c)`, a place of the storage
    const at = struct {
        inline fn f(p: [*]T, ld: usize, r: usize, c: usize) *T {
            return &p[r + c * ld];
        }
    }.f;
    // `A(r, c)` of the matrix, `c - kv <= r`
    const el = struct {
        inline fn f(p: [*]T, ld: usize, v: usize, r: usize, c: usize) *T {
            return &p[v + r - c + c * ld];
        }
    }.f;
    // fill-in of columns `ku + 1` to `kv - 1` (0-based), zeroed
    for (ku + 1..@max(ku + 1, @min(kv, n))) |j| for (kv - j..kl) |i| {
        at(ab, ldab, i, j).* = zero;
    };
    var info: ?usize = null;
    var ju: usize = 0;
    for (0..@min(m, n)) |j| {
        if (j + kv < n) for (0..kl) |i| {
            at(ab, ldab, i, j + kv).* = zero;
        };
        const km = @min(kl, m - j - 1);
        // the pivot among the `km + 1` rows from the diagonal
        var jp: usize = 0;
        var big: F = -1;
        for (0..km + 1) |i| {
            const v = at(ab, ldab, kv + i, j).*;
            const x = if (comptime isComplex(T)) @abs(v.re) + @abs(v.im) else @abs(v);
            if (std.math.isNan(x)) {
                jp = i;
                break;
            }
            if (x > big) {
                big = x;
                jp = i;
            }
        }
        ipiv[j] = jp + j;
        const piv = at(ab, ldab, kv + jp, j).*;
        if (ops.isZero(T, piv)) {
            if (info == null) info = j;
            continue;
        }
        ju = @max(ju, @min(j + ku + jp, n - 1));
        // rows `j` and `j + jp` over columns `j..ju`, along the band's diagonals
        if (jp != 0) for (j..ju + 1) |c| {
            mem.swap(T, el(ab, ldab, kv, j + jp, c), el(ab, ldab, kv, j, c));
        };
        if (km == 0) continue;
        ops.scaleBelow(T, (ab + kv + 1 + j * ldab)[0..km], at(ab, ldab, kv, j).*);
        // A(j+1:, j+1:ju) -= l u, `xGER` with -1
        for (j + 1..ju + 1) |c| {
            const y = el(ab, ldab, kv, j, c).*;
            if (ops.isZero(T, y)) continue;
            const t = ops.neg(T, y);
            for (1..km + 1) |i| {
                const dst = el(ab, ldab, kv, j + i, c);
                dst.* = ops.add(T, dst.*, ops.mul(T, at(ab, ldab, kv + i, j).*, t));
            }
        }
    }
    return info;
}

// ---- tests ---------------------------------------------------------------------

const testing = std.testing;
const check = @import("../blas/_test/check.zig");
const reference = @import("../blas/_test/reference.zig");
const lcheck = @import("_test/check.zig");

const W = lcheck.W;
const test_types = lcheck.test_types;
const isNanT = lcheck.isNanT;

/// `gbtrf` on a random `m x n` band matrix (`kl`, `ku`) at `ldab = 2 kl + ku
/// + 2`, the unused places NaN: `P1 L1 P2 L2 ... U` rebuilt in `W(T)` within
/// `2 (min(m, n) + 2)` ulps of the same product in magnitude, the multipliers
/// at most 1 in `|re| + |im|`. `zero_col` zeroes that column first (singular
/// there).
fn gbCase(comptime T: type, al: Allocator, m: usize, n: usize, kl: usize, ku: usize, seed: u64, zero_col: ?usize) !void {
    const kv = kl + ku;
    const ldab = 2 * kl + ku + 2;
    const ab = try al.alloc(T, ldab * n);
    defer al.free(ab);
    @memset(ab, lcheck.poisoned(T));
    const aw = try al.alloc(W(T), m * n);
    defer al.free(aw);
    const v = try al.alloc(T, m * n);
    defer al.free(v);
    check.fill(T, v, seed);
    for (0..n) |j| for (0..m) |i| {
        const inside = i <= j + kl and j <= i + ku and zero_col != j;
        const x = if (inside) v[i + j * m] else scalar(T, 0);
        aw[i + j * m] = reference.wide(T, x);
        if (i + ku >= j and i <= j + kl) ab[kv + i - j + j * ldab] = x;
    };
    // the fill-in rows, which `gbtrf` zeroes where it uses them
    for (0..n) |j| for (0..kl) |r| {
        ab[r + j * ldab] = lcheck.poisoned(T);
    };
    const ipiv = try al.alloc(usize, @max(1, @min(m, n)));
    defer al.free(ipiv);
    const piv = try gbtrf(T, m, n, kl, ku, ab, ldab, ipiv);
    const mn = @min(m, n);
    try testing.expectEqual(zero_col, piv);
    // X := U, then X := P_j L_j X from the last step back
    const x = try al.alloc(W(T), m * n);
    defer al.free(x);
    const g = try al.alloc(f128, m * n);
    defer al.free(g);
    for (0..n) |j| for (0..m) |i| {
        const u = if (i <= j and j <= i + kv and i < mn) reference.wide(T, ab[kv + i - j + j * ldab]) else reference.wide(T, scalar(T, 0));
        x[i + j * m] = u;
        g[i + j * m] = reference.abs1(u);
    };
    var s = mn;
    while (s > 0) {
        s -= 1;
        const km = @min(kl, m - s - 1);
        for (1..km + 1) |r| {
            const l = reference.wide(T, ab[kv + r + s * ldab]);
            // the pivot largest in `|re| + |im|`: `|l|` at most `sqrt 2` complex
            const l2 = if (comptime isComplex(T)) l.re * l.re + l.im * l.im else l * l;
            if (!(l2 <= (if (comptime isComplex(T)) 2 else 1) * (1 + 8 * std.math.floatEps(Float(T))))) {
                check.say("gbtrf {s}: |L({d},{d})|^2 = {e}\n", .{ @typeName(T), s + r, s, @as(f64, @floatCast(l2)) });
                return error.TestUnexpectedResult;
            }
            for (0..n) |j| {
                x[s + r + j * m] = reference.addW(x[s + r + j * m], reference.mulW(l, x[s + j * m]));
                g[s + r + j * m] += reference.abs1(l) * g[s + j * m];
            }
        }
        const p = ipiv[s];
        if (p != s) for (0..n) |j| {
            mem.swap(W(T), &x[s + j * m], &x[p + j * m]);
            mem.swap(f128, &g[s + j * m], &g[p + j * m]);
        };
    }
    const eps: f128 = std.math.floatEps(Float(T));
    const ulps: f128 = @floatFromInt(2 * (mn + 2));
    for (0..n) |j| for (0..m) |i| {
        const d = reference.abs1(reference.subW(x[i + j * m], aw[i + j * m]));
        if (!(d <= ulps * eps * g[i + j * m] + std.math.floatMin(Float(T)))) {
            check.say("gbtrf {s} {d}x{d} kl={d} ku={d}: (P L U)({d},{d}) off by {e}\n", .{ @typeName(T), m, n, kl, ku, i, j, @as(f64, @floatCast(d)) });
            return error.TestUnexpectedResult;
        }
    };
    // the padding row, never used
    for (0..n) |j| if (!isNanT(T, ab[ldab - 1 + j * ldab])) {
        check.say("gbtrf {s}: wrote the padding of column {d}\n", .{ @typeName(T), j });
        return error.TestUnexpectedResult;
    };
}

const gb_cases = [_][4]usize{ .{ 1, 1, 0, 0 }, .{ 5, 5, 0, 0 }, .{ 30, 30, 2, 3 }, .{ 60, 60, 10, 5 }, .{ 50, 40, 45, 3 }, .{ 40, 50, 3, 45 }, .{ 100, 100, 30, 20 }, .{ 80, 100, 90, 7 }, .{ 3, 7, 1, 1 } };

test "gbtrf: P L U = A in band storage, tall, wide and square, every type" {
    const al = testing.allocator;
    inline for (test_types) |T| for (gb_cases) |c| {
        if (check.rig() and c[0] * c[1] > 4000) continue;
        gbCase(T, al, c[0], c[1], c[2], c[3], c[0] + c[2] * 7, null) catch |e| {
            check.say("gbtrf {s} {d}x{d} kl={d} ku={d}\n", .{ @typeName(T), c[0], c[1], c[2], c[3] });
            return e;
        };
    };
}

test "gbtrf: a zero column is singular there, the factorization completed" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64 }) |T| for ([_]usize{ 0, 7, 39 }) |z| {
        try gbCase(T, al, 40, 40, 6, 4, 3, z);
    };
}

test "gbtrf: a denormal pivot is swapped in and divided by" {
    inline for (test_types) |T| {
        const tm = std.math.floatTrueMin(Float(T));
        const nan = lcheck.poisoned(T);
        var ipiv: [2]usize = @splat(std.math.maxInt(usize));
        // [[4tm, 1], [8tm, 3]], kl = ku = 1: rows swapped, L(2, 1) = 1/2,
        // U(2, 2) = -1/2
        var ab: [8]T = .{ nan, nan, scalar(T, 4 * tm), scalar(T, 8 * tm), nan, scalar(T, 1), scalar(T, 3), nan };
        try testing.expectEqual(@as(?usize, null), try gbtrf(T, 2, 2, 1, 1, &ab, 4, &ipiv));
        try testing.expectEqual([2]usize{ 1, 1 }, ipiv);
        try testing.expectEqual([4]T{ scalar(T, 8 * tm), scalar(T, 0.5), scalar(T, 3), scalar(T, -0.5) }, .{ ab[2], ab[3], ab[5], ab[6] });
        if (comptime isComplex(T)) {
            // 4tm over 8tm i: -i/2
            ab = .{ nan, nan, .init(4 * tm, 0), .init(0, 8 * tm), nan, scalar(T, 1), scalar(T, 3), nan };
            try testing.expectEqual(@as(?usize, null), try gbtrf(T, 2, 2, 1, 1, &ab, 4, &ipiv));
            try testing.expectEqual(T.init(0, -0.5), ab[3]);
        }
    }
}

test "gbtrf: a slice too short or leading dimension too small, nothing read when empty" {
    var ab: [8]f64 = @splat(lcheck.poisoned(f64));
    var ipiv: [2]usize = @splat(std.math.maxInt(usize));
    try testing.expectError(error.LeadingDimensionTooSmall, gbtrf(f64, 2, 2, 1, 1, &ab, 3, &ipiv));
    try testing.expectError(error.SliceTooShort, gbtrf(f64, 2, 3, 1, 1, &ab, 4, &ipiv));
    try testing.expectError(error.SliceTooShort, gbtrf(f64, 2, 2, 1, 1, &ab, 4, ipiv[0..1]));
    try testing.expectEqual(@as(?usize, null), try gbtrf(f64, 0, 2, 1, 1, &ab, 4, &ipiv));
    for (ab) |x| try testing.expect(std.math.isNan(x));
    try testing.expectEqual(std.math.maxInt(usize), ipiv[0]);
}
