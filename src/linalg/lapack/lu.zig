//! LU factorization with partial pivoting. `getrf` is OpenBLAS's
//! `interface/lapack/getrf.c` and `zgetrf.c` over
//! `lapack/getrf/getrf_single.c`, whose base case is `lapack/getf2/getf2_k.c`
//! and `zgetf2_k.c` (OpenBLAS v0.3.34-317). The semantics are
//! Reference-LAPACK's `xGETRF2` (`lapack-netlib/SRC/zgetrf2.f` in that
//! checkout).
//!
//! Where they differ, the reference wins:
//! - a nonzero pivot below the safe minimum (denormal, or NaN) is swapped in
//!   and divided by, where upstream's `getf2` neither swaps nor scales yet
//!   records the swap in `ipiv`, so `P L U` is not `A`;
//! - complex: that test is on `|p|`, as `zgetrf2.f`'s `ABS`, where upstream
//!   tests each part.
//!
//! The updates are the level-3 drivers' `trsm` and `gemm`, called directly as
//! upstream calls its kernels, so Accelerate never runs them; upstream fuses
//! the interchanges, `trsm` and `gemm` per `GEMM_R` columns, the same
//! arithmetic. Its interchanges are `laswp_plus` then `gemm_oncopy` (the
//! `#if 1`), not `laswp_ncopy`, and nothing calls `neg_tcopy`; nor here.
//! Serial: `getrf_parallel.c` is not ported. `getf2`'s `dot`, `gemv`, `swap`
//! and `scal` are loops, its `iamax` the core's kernel: the pivot is the
//! first NaN, else the first largest `|re| + |im|`, as `izamax`.

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const portable = @import("../blas/portable.zig");
const driver = @import("../blas/blas3/driver.zig");
const laswp = @import("../blas/ext/laswp/kernel.zig");
const Vtable = @import("../blas/Vtable.zig");
const chol = @import("cholesky.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const isComplex = types.isComplex;
const scalar = chol.scalar;

/// What `getrf` returns: `Singular` when `U` has a zero on its diagonal.
pub const GetrfError = error{ InvalidM, InvalidN, InvalidLda, Singular } || Allocator.Error || Io.Cancelable;

/// `?GETRF`: `A = P L U` for the `m x n` `A`, in place: `L` unit lower
/// trapezoidal below the diagonal, `U` upper trapezoidal on and above it.
/// Row `k` (1-based) was interchanged with row `ipiv[k - 1]`, for the
/// `min(m, n)` entries of `ipiv`.
///
/// When `U(i, i)` is exactly zero, the factorization still completes, as the
/// reference's does, then returns `error.Singular` with `pivot.*` set to the
/// first such `i` (the reference's `INFO`); `A` and `ipiv` hold the whole
/// factorization. `pivot` is written only then. A NaN pivot is not zero:
/// its column below fills with NaN, and nothing fails.
pub fn getrf(allocator: Allocator, io: Io, comptime T: type, m: isize, n: isize, a: [*]T, lda: isize, ipiv: [*]isize, pivot: ?*usize) GetrfError!void {
    comptime chol.floatType(T, "getrf");
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, m)) return error.InvalidLda;
    if (m == 0 or n == 0) return;
    if (!util.fitsMatrix(m, n, lda, @sizeOf(T))) return error.InvalidLda;
    const vt = chol.table(T);
    const blk = driver.Blocking.of(T, vt);
    const k = try factor(allocator, io, T, vt, blk.un, blk.q, @intCast(m), @intCast(n), a, @intCast(lda), ipiv);
    if (k == 0) return;
    if (pivot) |p| p.* = k;
    return error.Singular;
}

/// `getrf_single.c`: the first zero pivot's index, or 0. Pivots are 1-based
/// from `a`'s first row. Blocks of `min(mn / 2, Q)` columns, rounded up to
/// `unroll_n`, each factored recursively, until they would be at most
/// `2 * unroll_n` (upstream's POWER8 build stops at `unroll_n`).
fn factor(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, un: usize, q: usize, m: usize, n: usize, a: [*]T, lda: usize, ipiv: [*]isize) driver.Error!usize {
    const mn = @min(m, n);
    const blocking = @min((mn / 2 + un - 1) / un * un, q);
    if (blocking <= 2 * un) return getf2(T, vt, m, n, a, lda, ipiv);
    const one = scalar(T, 1);
    const minus = scalar(T, -1);
    var info: usize = 0;
    var j: usize = 0;
    while (j < mn) : (j += blocking) {
        const jb = @min(mn - j, blocking);
        const d = a + j * (lda + 1);
        const k = try factor(allocator, io, T, vt, un, q, m - j, jb, d, lda, ipiv + j);
        if (k != 0 and info == 0) info = k + j;
        for (ipiv[j..][0..jb]) |*p| p.* += @intCast(j);
        const rest = n - j - jb;
        if (rest == 0) continue;
        // A12 := L11^-1 P A12; A22 -= A21 A12
        const c = a + (j + jb) * lda;
        laswp.swapRows(T, rest, c, lda, j + 1, j + jb, ipiv + j, 1);
        try driver.trsm(allocator, io, T, vt, true, false, false, false, true, jb, rest, one, d, lda, c + j, lda, 1);
        if (m > j + jb) try driver.gemm(allocator, io, T, vt, false, false, .n, m - j - jb, rest, jb, minus, d + jb, lda, c + j, lda, one, c + j + jb, lda, 1);
    }
    // the later interchanges, into each block's columns
    j = 0;
    while (j < mn) : (j += blocking) {
        const jb = @min(mn - j, blocking);
        laswp.swapRows(T, jb, a + j * lda, lda, j + jb + 1, mn, ipiv + j + jb, 1);
    }
    return info;
}

/// `getf2_k.c`, `zgetf2_k.c`: left-looking, a column at a time. Column `j`
/// takes the earlier interchanges, then `U(:j, j)` by forward substitution,
/// then `A(j:, j) -= L(j:, :j) U(:j, j)`, its pivot, and `L(j + 1:, j)`.
fn getf2(comptime T: type, vt: *const Vtable, m: usize, n: usize, a: [*]T, lda: usize, ipiv: [*]isize) usize {
    const iamax = portable.kernel(T, vt, .iamax_k);
    var info: usize = 0;
    for (0..n) |j| {
        const b = a + j * lda;
        const len = @min(j, m);
        for (0..len) |i| {
            const ip: usize = @intCast(ipiv[i] - 1);
            if (ip != i) std.mem.swap(T, &b[i], &b[ip]);
        }
        for (1..@max(1, len)) |i| {
            var s = scalar(T, 0);
            for (0..i) |k| s = add(T, s, mul(T, a[i + k * lda], b[k]));
            b[i] = sub(T, b[i], s);
        }
        if (j >= m) continue;
        for (0..j) |k| {
            const x = b[k];
            const l = a + k * lda;
            for (j..m) |i| b[i] = sub(T, b[i], mul(T, l[i], x));
        }
        const p = j + @as(usize, @intCast(iamax(@intCast(m - j), @ptrCast(b + j), 1))) - 1;
        ipiv[j] = @intCast(p + 1);
        const piv = b[p];
        if (isZero(T, piv)) {
            if (info == 0) info = j + 1;
            continue;
        }
        if (p != j) for (0..j + 1) |k| std.mem.swap(T, &a[j + k * lda], &a[p + k * lda]);
        scaleBelow(T, b[j + 1 .. m], piv);
    }
    return info;
}

/// `x /= p`, as `xGETRF2`'s base case: times `1 / p` where `|p|` is at least
/// the safe minimum, else divided by `p`, whose reciprocal may overflow.
/// Complex quotients by Smith's algorithm, as `util.cdiv` and upstream.
fn scaleBelow(comptime T: type, x: []T, p: T) void {
    const F = Float(T);
    const sfmin = std.math.floatMin(F);
    if (comptime isComplex(T)) {
        if (std.math.hypot(p.re, p.im) >= sfmin) {
            const r = util.cdiv(F, 1, 0, p.re, p.im);
            for (x) |*v| {
                const z = util.cmul(F, r[0], r[1], v.re, v.im);
                v.* = .init(z[0], z[1]);
            }
        } else for (x) |*v| {
            const z = util.cdiv(F, v.re, v.im, p.re, p.im);
            v.* = .init(z[0], z[1]);
        }
    } else if (@abs(p) >= sfmin) {
        const r = 1 / p;
        for (x) |*v| v.* *= r;
    } else for (x) |*v| {
        v.* /= p;
    }
}

// ---- element arithmetic ---------------------------------------------------------

fn isZero(comptime T: type, v: T) bool {
    return if (comptime isComplex(T)) v.re == 0 and v.im == 0 else v == 0;
}

fn add(comptime T: type, x: T, y: T) T {
    return if (comptime isComplex(T)) .init(x.re + y.re, x.im + y.im) else x + y;
}

fn sub(comptime T: type, x: T, y: T) T {
    return if (comptime isComplex(T)) .init(x.re - y.re, x.im - y.im) else x - y;
}

fn mul(comptime T: type, x: T, y: T) T {
    if (comptime !isComplex(T)) return x * y;
    return .init(x.re * y.re - x.im * y.im, x.re * y.im + x.im * y.re);
}

// ---- tests ---------------------------------------------------------------------

const testing = std.testing;
const blas = @import("../blas/__root__.zig");
const check = @import("../blas/_test/blas1/check.zig");
const guard = @import("../blas/_test/blas2/guard.zig");

const W = chol.W;
const test_types = chol.test_types;
const poisoned = chol.poisoned;
const isNanT = chol.isNanT;
const fromF64 = chol.fromF64;

const nan = std.math.nan(f64);
const inf = std.math.inf(f64);

/// A random `m x n` `A` in [-2, 2), stored with `lda = m + 2`, the padding
/// rows NaN; `aw` keeps it in `W(T)`.
fn Problem(comptime T: type) type {
    return struct {
        m: usize,
        n: usize,
        lda: usize,
        a: []T,
        aw: []W(T),
        ipiv: []isize,

        const Self = @This();

        fn init(al: Allocator, m: usize, n: usize, seed: u64) !Self {
            const lda = m + 2;
            const a = try al.alloc(T, @max(1, lda * n));
            errdefer al.free(a);
            const aw = try al.alloc(W(T), m * n);
            errdefer al.free(aw);
            const ipiv = try al.alloc(isize, @max(1, @min(m, n)));
            const v = try al.alloc(T, m * n);
            defer al.free(v);
            check.fill(T, v, seed);
            for (aw, v) |*w, x| w.* = check.wide(T, x);
            var self: Self = .{ .m = m, .n = n, .lda = lda, .a = a, .aw = aw, .ipiv = ipiv };
            self.store();
            return self;
        }

        /// `aw` into `a`, the padding poisoned; `ipiv` poisoned too.
        fn store(self: *Self) void {
            @memset(self.a, poisoned(T));
            for (0..self.n) |j| for (0..self.m) |i| {
                self.a[i + j * self.lda] = chol.toT(T, self.aw[i + j * self.m]);
            };
            @memset(self.ipiv, -1);
        }

        fn deinit(self: *Self, al: Allocator) void {
            al.free(self.a);
            al.free(self.aw);
            al.free(self.ipiv);
        }

        fn factorIt(self: *Self, al: Allocator, pivot: ?*usize) GetrfError!void {
            return getrf(al, check.io, T, @intCast(self.m), @intCast(self.n), self.a.ptr, @intCast(self.lda), self.ipiv.ptr, pivot);
        }

        fn l(self: *const Self, i: usize, k: usize) W(T) {
            return if (i == k) check.wide(T, scalar(T, 1)) else check.wide(T, self.a[i + k * self.lda]);
        }

        fn u(self: *const Self, k: usize, j: usize) W(T) {
            return check.wide(T, self.a[k + j * self.lda]);
        }

        /// `A` with `ipiv`'s interchanges applied, in `W(T)`.
        fn permuted(self: *const Self, al: Allocator) ![]W(T) {
            const pa = try al.dupe(W(T), self.aw);
            for (0..@min(self.m, self.n)) |k| {
                const p: usize = @intCast(self.ipiv[k] - 1);
                for (0..self.n) |j| std.mem.swap(W(T), &pa[k + j * self.m], &pa[p + j * self.m]);
            }
            return pa;
        }

        /// Pivots inside their rows, `|L|` at most 1 (complex: `sqrt 2`, the
        /// pivot being largest in `|re| + |im|`), `L U` against `P A` within
        /// `2 (mn + 2)` ulps of `|L| |U|`, and the padding untouched.
        fn expectFactor(self: *const Self, al: Allocator) !void {
            const m = self.m;
            const n = self.n;
            const mn = @min(m, n);
            const F = Float(T);
            const eps: f128 = std.math.floatEps(F);
            for (self.ipiv[0..mn], 0..) |p, k| if (p < k + 1 or p > m) {
                check.say("getrf {d}x{d}: ipiv[{d}] = {d}\n", .{ m, n, k, p });
                return error.TestUnexpectedResult;
            };
            const bound: f128 = (if (comptime isComplex(T)) 2 else 1) * (1 + 8 * eps);
            for (0..mn) |k| for (k + 1..m) |i| {
                const x = self.l(i, k);
                const s2 = if (comptime isComplex(T)) x.re * x.re + x.im * x.im else x * x;
                if (!(s2 <= bound)) {
                    check.say("getrf {d}x{d}: |L({d},{d})|^2 = {e}\n", .{ m, n, i, k, @as(f64, @floatCast(s2)) });
                    return error.TestUnexpectedResult;
                }
            };
            const pa = try self.permuted(al);
            defer al.free(pa);
            const ulps: f128 = @floatFromInt(2 * (mn + 2));
            for (0..n) |j| for (0..m) |i| {
                var s: W(T) = check.wide(T, scalar(T, 0));
                var mag: f128 = 0;
                for (0..@min(i, j, mn - 1) + 1) |k| {
                    const x = self.l(i, k);
                    const y = self.u(k, j);
                    s = check.addW(s, check.mulW(x, y));
                    mag += check.abs1(x) * check.abs1(y);
                }
                const d = check.abs1(check.subW(s, pa[i + j * m]));
                if (!(d <= ulps * eps * mag + std.math.floatMin(F))) {
                    check.say("getrf {d}x{d}: (L U)({d},{d}) off by {e}, allowed {e}\n", .{ m, n, i, j, @as(f64, @floatCast(d)), @as(f64, @floatCast(ulps * eps * mag)) });
                    return error.TestUnexpectedResult;
                }
            };
            for (0..n) |j| for (m..self.lda) |i| if (!isNanT(T, self.a[i + j * self.lda])) {
                check.say("getrf {d}x{d}: wrote padding ({d},{d})\n", .{ m, n, i, j });
                return error.TestUnexpectedResult;
            };
        }
    };
}

/// `vt` with every precision's `q` at `2 * unroll_n + 1`, so the blocks are
/// capped at an odd width just past the base case's; static, as a vtable is
/// large for a rig's stack.
var shrunk: Vtable = undefined;

fn shrink(vt: *const Vtable) *const Vtable {
    shrunk = vt.*;
    inline for (.{ "s", "d", "q", "c", "z", "x" }) |p| {
        const un = @max(1, @field(shrunk, p ++ "gemm_unroll_n"));
        @field(shrunk, p ++ "gemm_q") = 2 * un + 1;
    }
    return &shrunk;
}

/// `m x n` shapes: square, tall, wide, a row and a column, across both
/// blockings' block edges.
const shapes = [_][2]usize{
    .{ 1, 1 },   .{ 2, 2 },   .{ 3, 3 },   .{ 7, 7 },   .{ 1, 9 },   .{ 9, 1 },   .{ 17, 17 },   .{ 33, 33 },
    .{ 40, 25 }, .{ 25, 40 }, .{ 50, 50 }, .{ 70, 70 }, .{ 90, 61 }, .{ 61, 90 }, .{ 130, 130 },
};

fn factorCases(vt: *const Vtable) anyerror!void {
    const al = testing.allocator;
    inline for (test_types) |T| for ([_]*const Vtable{ shrink(vt), vt }, 0..) |v, vi| {
        const saved = blas.testUse(v);
        defer blas.testRestore(saved);
        for (shapes) |s| {
            if (check.rig() and s[0] * s[1] > 70 * 70) continue;
            var p = try Problem(T).init(al, s[0], s[1], s[0] * 31 + s[1]);
            defer p.deinit(al);
            var poison: guard.Poison = .{ .child = al };
            var pivot: usize = 99;
            p.factorIt(poison.allocator(), &pivot) catch |e| {
                check.say("getrf {s} {d}x{d}\n", .{ @typeName(T), s[0], s[1] });
                return e;
            };
            try testing.expectEqual(@as(usize, 99), pivot);
            p.expectFactor(al) catch |e| {
                check.say("{s}, blocking {s}\n", .{ @typeName(T), if (vi == 0) "shrunk" else "own" });
                return e;
            };
        }
    };
}

test "getrf: P L U = A, square, tall and wide, on every linked core" {
    try check.eachCore(factorCases);
}

/// `getrf` on `p` must fail at pivot `want`, its factorization complete.
fn expectSingular(comptime T: type, al: Allocator, p: *Problem(T), want: usize, what: []const u8) !void {
    var pivot: usize = 0;
    var poison: guard.Poison = .{ .child = al };
    if (p.factorIt(poison.allocator(), &pivot)) |_| {
        check.say("{s} {s} {d}x{d}: no error, want pivot {d}\n", .{ what, @typeName(T), p.m, p.n, want });
        return error.TestExpectedError;
    } else |e| if (e != error.Singular or pivot != want) {
        check.say("{s} {s} {d}x{d}: {s} at pivot {d}, want {d}\n", .{ what, @typeName(T), p.m, p.n, @errorName(e), pivot, want });
        return error.TestUnexpectedResult;
    }
    p.expectFactor(al) catch |e| {
        check.say("{s} {s} {d}x{d}, singular at {d}\n", .{ what, @typeName(T), p.m, p.n, want });
        return e;
    };
}

/// Columns made zero: the first and last, and either side of block edges on
/// both blockings.
const zero_at = [_]usize{ 0, 1, 2, 5, 12, 13, 17, 25, 32, 33, 69 };

fn singularCases(vt: *const Vtable) anyerror!void {
    const al = testing.allocator;
    inline for (test_types) |T| for ([_]*const Vtable{ shrink(vt), vt }, 0..) |v, vi| {
        const saved = blas.testUse(v);
        defer blas.testRestore(saved);
        const n: usize = if (vi == 0) 40 else 70;
        for ([_]usize{ n, n + 9 }) |m| {
            var p = try Problem(T).init(al, m, n, 5 + m);
            defer p.deinit(al);
            const aw = try al.dupe(W(T), p.aw);
            defer al.free(aw);
            for (zero_at) |k| {
                if (k >= n) continue;
                // a zero column `k`, then a second one after it
                @memcpy(p.aw, aw);
                for (0..m) |i| p.aw[i + k * m] = check.wide(T, scalar(T, 0));
                if (k + 3 < n) for (0..m) |i| {
                    p.aw[i + (k + 3) * m] = check.wide(T, scalar(T, 0));
                };
                p.store();
                try expectSingular(T, al, &p, k + 1, "zero column");
            }
            // rank one, `x yᵀ` with `x` of powers of two and `y(0) = 1`: every
            // step is exact, so `U(2, 2)` is zero
            for (0..m) |i| for (0..n) |j| {
                const xi = std.math.ldexp(@as(f64, if (i % 3 == 1) -1 else 1), -@as(i32, @intCast((i * 7) % 5)));
                const yj: f64 = if (j == 0) 1 else @floatFromInt(@as(i32, @intCast(j % 7)) - 3);
                p.aw[i + j * m] = check.wide(T, fromF64(T, xi * yj));
            };
            p.store();
            try expectSingular(T, al, &p, 2, "rank one");
        }
    };
}

test "getrf: the first zero pivot, the factorization completed, on every linked core" {
    try check.eachCore(singularCases);
}

/// A NaN anywhere is a pivot's or reaches one, and is never zero: no error,
/// even where it reaches a zero column, as the reference's `INFO`.
fn nanCases(vt: *const Vtable) anyerror!void {
    const al = testing.allocator;
    inline for (test_types) |T| for ([_]*const Vtable{ shrink(vt), vt }) |v| {
        const saved = blas.testUse(v);
        defer blas.testRestore(saved);
        const n = 40;
        var p = try Problem(T).init(al, n, n, 11);
        defer p.deinit(al);
        for ([_][2]usize{ .{ 0, 0 }, .{ 7, 3 }, .{ 39, 20 }, .{ 0, 39 } }) |at| {
            p.store();
            p.a[at[0] + at[1] * p.lda] = fromF64(T, nan);
            // a zero column past the NaN's
            const z = @min(at[1] + 2, n - 1);
            if (z != at[1]) for (0..n) |i| {
                p.a[i + z * p.lda] = scalar(T, 0);
            };
            var pivot: usize = 99;
            var poison: guard.Poison = .{ .child = al };
            p.factorIt(poison.allocator(), &pivot) catch |e| {
                check.say("getrf {s}: NaN at ({d},{d}): {s} at {d}\n", .{ @typeName(T), at[0], at[1], @errorName(e), pivot });
                return e;
            };
            // the last pivot is NaN
            try testing.expect(isNanT(T, p.a[(n - 1) * (p.lda + 1)]));
        }
    };
}

test "getrf: NaN is never a zero pivot, on every linked core" {
    try check.eachCore(nanCases);
}

/// The `m x 1` column `col` factored: its pivot (1-based) and `L`.
fn column(comptime T: type, comptime m: usize, col: [m]T) !struct { usize, [m]T } {
    var a = col;
    var ipiv: [1]isize = .{-1};
    try getrf(testing.failing_allocator, check.io, T, m, 1, &a, m, &ipiv, null);
    return .{ @intCast(ipiv[0]), a };
}

test "getrf: the first largest |re| + |im| pivots, the first NaN before it" {
    inline for (test_types) |T| {
        const c = struct {
            fn f(re: f64, im: f64) T {
                return if (comptime isComplex(T)) .init(@floatCast(re), @floatCast(im)) else @floatCast(re);
            }
        }.f;
        // ties keep the first; no swap for the diagonal's
        try testing.expectEqual(@as(usize, 1), (try column(T, 3, .{ c(1, 0), c(-1, 0), c(1, 0) }))[0]);
        try testing.expectEqual(@as(usize, 2), (try column(T, 4, .{ c(-1, 0), c(2, 0), c(-2, 0), c(2, 0) }))[0]);
        if (comptime isComplex(T)) {
            // 7, 5, 7, 7: the first; and 3 - 4i over 5, the same modulus
            try testing.expectEqual(@as(usize, 1), (try column(T, 4, .{ c(3, 4), c(5, 0), c(-4, 3), c(0, -7) }))[0]);
            try testing.expectEqual(@as(usize, 3), (try column(T, 3, .{ c(1, 0), c(5, 0), c(3, -4) }))[0]);
        }
        // the first NaN wins, and fills L
        const r = try column(T, 4, .{ c(1, 0), c(nan, 0), c(5, 0), c(nan, 0) });
        try testing.expectEqual(@as(usize, 2), r[0]);
        for (r[1][1..]) |x| try testing.expect(if (comptime isComplex(T)) std.math.isNan(x.re) else std.math.isNan(x));
        // an infinite pivot: finite elements over it are zero, infinite ones NaN
        const s = try column(T, 3, .{ c(1, 0), c(inf, 0), c(-inf, 0) });
        try testing.expectEqual(@as(usize, 2), s[0]);
        try testing.expectEqual(c(0, 0), s[1][1]);
        try testing.expect(if (comptime isComplex(T)) std.math.isNan(s[1][2].re) else std.math.isNan(s[1][2]));
    }
}

test "getrf: a denormal pivot is swapped in and divided by" {
    inline for (test_types) |T| {
        const F = Float(T);
        const tm = std.math.floatTrueMin(F);
        // [[4tm, 1], [8tm, 3]]: rows swapped, L(2, 1) = 1/2, U(2, 2) = -1/2
        var a: [4]T = .{ scalar(T, 4 * tm), scalar(T, 8 * tm), scalar(T, 1), scalar(T, 3) };
        var ipiv: [2]isize = @splat(-1);
        try getrf(testing.failing_allocator, check.io, T, 2, 2, &a, 2, &ipiv, null);
        try testing.expectEqual([2]isize{ 2, 2 }, ipiv);
        try testing.expectEqual([4]T{ scalar(T, 8 * tm), scalar(T, 0.5), scalar(T, 3), scalar(T, -0.5) }, a);
        if (comptime isComplex(T)) {
            // 4tm over 8tm i: -i/2
            const r = try column(T, 2, .{ .init(4 * tm, 0), .init(0, 8 * tm) });
            try testing.expectEqual(@as(usize, 2), r[0]);
            try testing.expectEqual(T.init(0, -0.5), r[1][1]);
        }
    }
}

test "getrf: an exactly singular U; pivot is written only on failure" {
    inline for (.{ f64, types.Complex64 }) |T| {
        const al = testing.failing_allocator;
        var pivot: usize = 99;
        var ipiv: [3]isize = @splat(-1);
        // a zero matrix fails at 1, a zero row's `U(1, 1)` too
        var z: [9]T = @splat(scalar(T, 0));
        try testing.expectError(error.Singular, getrf(al, check.io, T, 3, 3, &z, 3, &ipiv, &pivot));
        try testing.expectEqual(@as(usize, 1), pivot);
        try testing.expectEqual([3]isize{ 1, 2, 3 }, ipiv);
        var row: [3]T = .{ scalar(T, 0), scalar(T, 2), scalar(T, 3) };
        try testing.expectError(error.Singular, getrf(al, check.io, T, 1, 3, &row, 1, &ipiv, &pivot));
        try testing.expectEqual(@as(usize, 1), pivot);
        // [[1, 2], [2, 4]]: rows swapped, then U(2, 2) = 2 - 1/2 4 = 0
        var a: [4]T = .{ scalar(T, 1), scalar(T, 2), scalar(T, 2), scalar(T, 4) };
        try testing.expectError(error.Singular, getrf(al, check.io, T, 2, 2, &a, 2, &ipiv, &pivot));
        try testing.expectEqual(@as(usize, 2), pivot);
        try testing.expectEqual([4]T{ scalar(T, 2), scalar(T, 0.5), scalar(T, 4), scalar(T, 0) }, a);
        a = .{ scalar(T, 1), scalar(T, 2), scalar(T, 2), scalar(T, 5) };
        pivot = 99;
        try getrf(al, check.io, T, 2, 2, &a, 2, &ipiv, &pivot);
        try testing.expectEqual(@as(usize, 99), pivot);
        try testing.expectEqual([4]T{ scalar(T, 2), scalar(T, 0.5), scalar(T, 5), scalar(T, -0.5) }, a);
    }
}

test "getrf: reports the first bad argument and touches nothing when empty" {
    const fa = testing.failing_allocator;
    var a: [4]f64 = @splat(nan);
    var ipiv: [2]isize = @splat(-1);
    try testing.expectError(error.InvalidM, getrf(fa, check.io, f64, -1, -1, &a, 0, &ipiv, null));
    try testing.expectError(error.InvalidN, getrf(fa, check.io, f64, 2, -1, &a, 0, &ipiv, null));
    try testing.expectError(error.InvalidLda, getrf(fa, check.io, f64, 2, 2, &a, 1, &ipiv, null));
    try testing.expectError(error.InvalidLda, getrf(fa, check.io, f64, 0, 2, &a, 0, &ipiv, null));
    try testing.expectError(error.InvalidLda, getrf(fa, check.io, f64, 2, 2, &a, std.math.maxInt(isize) / 2, &ipiv, null));
    try getrf(fa, check.io, f64, 0, 2, &a, 1, &ipiv, null);
    try getrf(fa, check.io, f64, 2, 0, &a, 2, &ipiv, null);
    for (a) |v| try testing.expect(std.math.isNan(v));
    try testing.expectEqual([2]isize{ -1, -1 }, ipiv);
}
