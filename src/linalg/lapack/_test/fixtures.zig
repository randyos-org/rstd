//! Test problems shared by the LAPACK tests: `potrf`'s, `getrf`'s and
//! `trtri`'s with known answers, ReLAPACK's comparisons, and the Schur forms
//! `trsyl` and `tgsyl` solve on.

const std = @import("std");
const types = @import("../../blas/types.zig");
const check = @import("../../blas/_test/check.zig");
const reference = @import("../../blas/_test/reference.zig");
const ops = @import("../ops.zig");
const chol = @import("../cholesky.zig");
const lu = @import("../lu.zig");
const lcheck = @import("check.zig");

const Allocator = std.mem.Allocator;
const Float = types.Float;
const isComplex = types.isComplex;
const scalar = ops.scalar;
const re = ops.re;
const W = lcheck.W;
const toT = lcheck.toT;
const poisoned = lcheck.poisoned;
const isNanT = lcheck.isNanT;

/// A test matrix: `L`, lower with a diagonal in [2, 4) and the rest within
/// `2 / n` (so `L` and `A = L Lᴴ` are well conditioned), and `A` built from
/// it in `W(T)`, stored in `uplo`'s triangle of `n x n` with `lda = n + 2`.
/// The other triangle, the padding rows and a complex diagonal's imaginary
/// parts hold NaN.
pub fn CholProblem(comptime T: type) type {
    return struct {
        n: usize,
        lda: usize,
        upper: bool,
        l: []T,
        a: []T,
        aw: []W(T),

        const Self = @This();

        pub fn init(al: Allocator, n: usize, upper: bool, seed: u64) !Self {
            const lda = n + 2;
            const l = try al.alloc(T, n * n);
            errdefer al.free(l);
            const a = try al.alloc(T, @max(1, lda * n));
            errdefer al.free(a);
            const aw = try al.alloc(W(T), n * n);
            check.fill(T, l, seed);
            const scale: Float(T) = @floatCast(1 / @as(f64, @floatFromInt(@max(n, 1))));
            for (0..n) |j| for (0..n) |i| {
                const v = &l[i + j * n];
                if (i < j) {
                    v.* = scalar(T, 0);
                } else if (i == j) {
                    v.* = scalar(T, 3 + re(T, v.*) / 2);
                } else if (comptime isComplex(T)) {
                    v.* = .init(v.re * scale, v.im * scale);
                } else {
                    v.* = v.* * scale;
                }
            };
            for (0..n) |j| for (0..n) |i| {
                var s: W(T) = reference.wide(T, scalar(T, 0));
                for (0..@min(i, j) + 1) |k| s = reference.addW(s, reference.mulW(reference.wide(T, l[i + k * n]), reference.conjW(reference.wide(T, l[j + k * n]))));
                aw[i + j * n] = s;
            };
            var self: Self = .{ .n = n, .lda = lda, .upper = upper, .l = l, .a = a, .aw = aw };
            self.store();
            return self;
        }

        /// `A` into `a`, poisoned where `potrf` must not read.
        pub fn store(self: *Self) void {
            @memset(self.a, poisoned(T));
            for (0..self.n) |j| for (0..self.n) |i| {
                if (if (self.upper) i > j else i < j) continue;
                const v = toT(T, self.aw[i + j * self.n]);
                self.a[i + j * self.lda] = if (comptime isComplex(T)) (if (i == j) .init(v.re, std.math.nan(Float(T))) else v) else v;
            };
        }

        pub fn deinit(self: *Self, al: Allocator) void {
            al.free(self.l);
            al.free(self.a);
            al.free(self.aw);
        }

        /// The computed factor's `(i, j)` as lower, `i >= j`: `L(i, j)`, or
        /// `conj(U(j, i))`.
        pub fn got(self: *const Self, i: usize, j: usize) T {
            if (!self.upper) return self.a[i + j * self.lda];
            const v = self.a[j + i * self.lda];
            return if (comptime isComplex(T)) v.conjugate() else v;
        }

        pub fn factorIt(self: *Self, al: Allocator, minor: ?*usize) chol.PotrfError!void {
            return chol.potrf(al, check.io, T, if (self.upper) .u else .l, @intCast(self.n), self.a.ptr, @intCast(self.lda), minor);
        }

        /// The factor against `A` (backward, within `2 (n + 2)` ulps of
        /// `|L| |L|ᴴ`) and against the known `L` (forward, loosely), with
        /// the poisoned elements untouched.
        pub fn expectFactor(self: *const Self) !void {
            const n = self.n;
            const F = Float(T);
            const eps: f128 = std.math.floatEps(F);
            const ulps: f128 = @floatFromInt(2 * (n + 2));
            for (0..n) |j| for (j..n) |i| {
                var s: W(T) = reference.wide(T, scalar(T, 0));
                var mag: f128 = 0;
                for (0..j + 1) |k| {
                    const x = reference.wide(T, self.got(i, k));
                    const y = reference.wide(T, self.got(j, k));
                    s = reference.addW(s, reference.mulW(x, reference.conjW(y)));
                    mag += reference.abs1(x) * reference.abs1(y);
                }
                const d = reference.abs1(reference.subW(s, self.aw[i + j * n]));
                if (!(d <= ulps * eps * mag)) {
                    check.say("potrf {s} n={d}: (L Lᴴ)({d},{d}) off by {e}, allowed {e}\n", .{ if (self.upper) "u" else "l", n, i, j, @as(f64, @floatCast(d)), @as(f64, @floatCast(ulps * eps * mag)) });
                    return error.TestUnexpectedResult;
                }
                const f = reference.abs1(reference.subW(reference.wide(T, self.got(i, j)), reference.wide(T, self.l[i + j * n])));
                if (!(f <= 64 * @as(f128, @floatFromInt(n)) * eps * 4)) {
                    check.say("potrf {s} n={d}: L({d},{d}) off by {e}\n", .{ if (self.upper) "u" else "l", n, i, j, @as(f64, @floatCast(f)) });
                    return error.TestUnexpectedResult;
                }
                if (comptime isComplex(T)) if (i == j and self.got(i, i).im != 0) {
                    check.say("potrf n={d}: L({d},{d}) has imaginary part {e}\n", .{ n, i, i, self.got(i, i).im });
                    return error.TestUnexpectedResult;
                };
            };
            for (0..n) |j| for (0..self.lda) |i| {
                const stored = i < n and (if (self.upper) i <= j else i >= j);
                if (!stored and !isNanT(T, self.a[i + j * self.lda])) {
                    check.say("potrf {s} n={d}: wrote ({d},{d}), outside its triangle\n", .{ if (self.upper) "u" else "l", n, i, j });
                    return error.TestUnexpectedResult;
                }
            };
        }
    };
}

/// A random `m x n` `A` in [-2, 2), stored with `lda = m + 2`, the padding
/// rows NaN; `aw` keeps it in `W(T)`.
pub fn LuProblem(comptime T: type) type {
    return struct {
        m: usize,
        n: usize,
        lda: usize,
        a: []T,
        aw: []W(T),
        ipiv: []isize,

        const Self = @This();

        pub fn init(al: Allocator, m: usize, n: usize, seed: u64) !Self {
            const lda = m + 2;
            const a = try al.alloc(T, @max(1, lda * n));
            errdefer al.free(a);
            const aw = try al.alloc(W(T), m * n);
            errdefer al.free(aw);
            const ipiv = try al.alloc(isize, @max(1, @min(m, n)));
            const v = try al.alloc(T, m * n);
            defer al.free(v);
            check.fill(T, v, seed);
            for (aw, v) |*w, x| w.* = reference.wide(T, x);
            var self: Self = .{ .m = m, .n = n, .lda = lda, .a = a, .aw = aw, .ipiv = ipiv };
            self.store();
            return self;
        }

        /// `aw` into `a`, the padding poisoned; `ipiv` poisoned too.
        pub fn store(self: *Self) void {
            @memset(self.a, poisoned(T));
            for (0..self.n) |j| for (0..self.m) |i| {
                self.a[i + j * self.lda] = toT(T, self.aw[i + j * self.m]);
            };
            @memset(self.ipiv, -1);
        }

        pub fn deinit(self: *Self, al: Allocator) void {
            al.free(self.a);
            al.free(self.aw);
            al.free(self.ipiv);
        }

        pub fn factorIt(self: *Self, al: Allocator, pivot: ?*usize) lu.GetrfError!void {
            return lu.getrf(al, check.io, T, @intCast(self.m), @intCast(self.n), self.a.ptr, @intCast(self.lda), self.ipiv.ptr, pivot);
        }

        pub fn l(self: *const Self, i: usize, k: usize) W(T) {
            return if (i == k) reference.wide(T, scalar(T, 1)) else reference.wide(T, self.a[i + k * self.lda]);
        }

        pub fn u(self: *const Self, k: usize, j: usize) W(T) {
            return reference.wide(T, self.a[k + j * self.lda]);
        }

        /// `A` with `ipiv`'s interchanges applied, in `W(T)`.
        pub fn permuted(self: *const Self, al: Allocator) ![]W(T) {
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
        pub fn expectFactor(self: *const Self, al: Allocator) !void {
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
                var s: W(T) = reference.wide(T, scalar(T, 0));
                var mag: f128 = 0;
                for (0..@min(i, j, mn - 1) + 1) |k| {
                    const x = self.l(i, k);
                    const y = self.u(k, j);
                    s = reference.addW(s, reference.mulW(x, y));
                    mag += reference.abs1(x) * reference.abs1(y);
                }
                const d = reference.abs1(reference.subW(s, pa[i + j * m]));
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

/// A test triangle in `upper`'s side of `n x n`, stored at `lda = n + 2`:
/// the diagonal in [2.5, 3.5), the rest within `1 / n`, so it and its
/// inverse are well conditioned; `w` keeps it whole in `W(T)`, a unit
/// diagonal as ones. The other triangle, the padding and a unit diagonal
/// hold NaN.
pub fn Tri(comptime T: type) type {
    return struct {
        n: usize,
        lda: usize,
        upper: bool,
        unit: bool,
        a: []T,
        w: []W(T),

        const Self = @This();

        pub fn init(al: Allocator, n: usize, upper: bool, unit: bool, seed: u64) !Self {
            const lda = n + 2;
            const a = try al.alloc(T, @max(1, lda * n));
            errdefer al.free(a);
            const w = try al.alloc(W(T), n * n);
            errdefer al.free(w);
            const v = try al.alloc(T, n * n);
            defer al.free(v);
            check.fill(T, v, seed);
            const scale: Float(T) = @floatCast(1 / @as(f64, @floatFromInt(@max(n, 1))));
            for (0..n) |j| for (0..n) |i| {
                const x = v[i + j * n];
                const e = if (if (upper) i > j else i < j)
                    scalar(T, 0)
                else if (i == j)
                    scalar(T, if (unit) 1 else 3 + ops.re(T, x) / 2)
                else
                    ops.mulReal(T, x, scale);
                w[i + j * n] = reference.wide(T, e);
            };
            var self: Self = .{ .n = n, .lda = lda, .upper = upper, .unit = unit, .a = a, .w = w };
            self.store();
            return self;
        }

        pub fn deinit(self: *Self, al: Allocator) void {
            al.free(self.a);
            al.free(self.w);
        }

        /// Whether the routines read and write `(i, j)` of `a`.
        pub fn stored(self: *const Self, i: usize, j: usize) bool {
            if (i >= self.n) return false;
            if (i == j) return !self.unit;
            return if (self.upper) i < j else i > j;
        }

        pub fn store(self: *Self) void {
            @memset(self.a, poisoned(T));
            for (0..self.n) |j| for (0..self.n) |i| if (self.stored(i, j)) {
                self.a[i + j * self.lda] = toT(T, self.w[i + j * self.n]);
            };
        }

        /// The computed matrix's `(i, j)`: ones on a unit diagonal, zeros
        /// outside the triangle.
        pub fn got(self: *const Self, i: usize, j: usize) W(T) {
            if (i == j and self.unit) return reference.wide(T, scalar(T, 1));
            if (!self.stored(i, j)) return reference.wide(T, scalar(T, 0));
            return reference.wide(T, self.a[i + j * self.lda]);
        }

        pub fn expectUntouched(self: *const Self, what: []const u8) !void {
            for (0..self.n) |j| for (0..self.lda) |i| if (!self.stored(i, j) and !isNanT(T, self.a[i + j * self.lda])) {
                check.say("{s} {s} {s} n={d}: wrote ({d},{d}), outside its triangle\n", .{ what, @typeName(T), if (self.upper) "u" else "l", self.n, i, j });
                return error.TestUnexpectedResult;
            };
        }

        /// `A X = I` within `8 (n + 2)` ulps of `|A| |X|`, `X` the computed
        /// inverse, nothing outside the triangle written.
        pub fn expectInverse(self: *const Self) !void {
            const n = self.n;
            const eps: f128 = std.math.floatEps(Float(T));
            const ulps: f128 = @floatFromInt(8 * (n + 2));
            for (0..n) |j| for (0..n) |i| {
                if (if (self.upper) i > j else i < j) continue;
                var s: W(T) = reference.wide(T, scalar(T, if (i == j) -1 else 0));
                var mag: f128 = 0;
                for (@min(i, j)..@max(i, j) + 1) |k| {
                    const x = self.w[i + k * n];
                    const y = self.got(k, j);
                    s = reference.addW(s, reference.mulW(x, y));
                    mag += reference.abs1(x) * reference.abs1(y);
                }
                const d = reference.abs1(s);
                if (!(d <= ulps * eps * mag)) {
                    check.say("trtri {s} {s}{s} n={d}: (A X - I)({d},{d}) = {e}, allowed {e}\n", .{ @typeName(T), if (self.upper) "u" else "l", if (self.unit) "u" else "n", n, i, j, @as(f64, @floatCast(d)), @as(f64, @floatCast(ulps * eps * mag)) });
                    return error.TestUnexpectedResult;
                }
            };
            try self.expectUntouched("trtri");
        }

        /// The stored triangle against `U Uᴴ` or `Lᴴ L` of `w` within
        /// `2 (n + 2)` ulps of the product in magnitude, a complex diagonal
        /// exactly real, nothing outside the triangle written.
        pub fn expectProduct(self: *const Self) !void {
            const n = self.n;
            const eps: f128 = std.math.floatEps(Float(T));
            const ulps: f128 = @floatFromInt(2 * (n + 2));
            for (0..n) |j| for (0..n) |i| {
                if (!self.stored(i, j)) continue;
                var s: W(T) = reference.wide(T, scalar(T, 0));
                var mag: f128 = 0;
                for (@max(i, j)..n) |k| {
                    // `U(i, k) conj(U(j, k))`, or `conj(L(k, i)) L(k, j)`
                    const x = if (self.upper) self.w[i + k * n] else reference.conjW(self.w[k + i * n]);
                    const y = if (self.upper) reference.conjW(self.w[j + k * n]) else self.w[k + j * n];
                    s = reference.addW(s, reference.mulW(x, y));
                    mag += reference.abs1(x) * reference.abs1(y);
                }
                const g = self.a[i + j * self.lda];
                const d = reference.abs1(reference.subW(reference.wide(T, g), s));
                if (!(d <= ulps * eps * mag)) {
                    check.say("lauum {s} {s} n={d}: ({d},{d}) off by {e}, allowed {e}\n", .{ @typeName(T), if (self.upper) "u" else "l", n, i, j, @as(f64, @floatCast(d)), @as(f64, @floatCast(ulps * eps * mag)) });
                    return error.TestUnexpectedResult;
                }
                if (comptime isComplex(T)) if (i == j and g.im != 0) {
                    check.say("lauum {s} n={d}: ({d},{d}) has imaginary part {e}\n", .{ @typeName(T), n, i, i, g.im });
                    return error.TestUnexpectedResult;
                };
            };
            try self.expectUntouched("lauum");
        }
    };
}

/// `x` and `y` within `64 n` ulps of `y`'s largest element, NaN where both
/// are.
pub fn expectAgree(comptime T: type, what: []const u8, n: usize, x: []const T, y: []const T) !void {
    var big: f128 = 0;
    for (y) |v| if (!isNanT(T, v)) {
        big = @max(big, reference.abs1(reference.wide(T, v)));
    };
    const tol = 64 * @as(f128, @floatFromInt(n)) * std.math.floatEps(Float(T)) * big;
    for (x, y, 0..) |u, v, i| {
        if (isNanT(T, u) and isNanT(T, v)) continue;
        const d = reference.abs1(reference.subW(reference.wide(T, u), reference.wide(T, v)));
        if (!(d <= tol)) {
            check.say("{s} {s} n={d}: element {d} differs by {e}, allowed {e}\n", .{ what, @typeName(T), n, i, @as(f64, @floatCast(d)), @as(f64, @floatCast(tol)) });
            return error.TestUnexpectedResult;
        }
    }
}

/// A diagonally dominant symmetric `n x n`, positive definite.
pub fn dominant(al: Allocator, n: usize, seed: u64) ![]f64 {
    const a = try al.alloc(f64, n * n);
    check.fill(f64, a, seed);
    for (0..n) |j| for (0..n) |i| {
        if (i > j) a[i + j * n] = a[j + i * n];
        if (i == j) a[i + j * n] += @floatFromInt(2 * n);
    };
    return a;
}

/// `trsyl`'s and `tgsyl`'s: an `n x n` upper triangular (real:
/// quasi-triangular, a 2 x 2 block every `every` rows) matrix whose
/// eigenvalues' real parts lie in `[lo, lo + 2)`, the rest within `1 / n`;
/// zero below.
pub fn schur(comptime T: type, al: Allocator, n: usize, lo: f64, every: usize, seed: u64) ![]T {
    const a = try al.alloc(T, n * n);
    check.fill(T, a, seed);
    const s: Float(T) = @floatCast(1 / @as(f64, @floatFromInt(@max(n, 1))));
    for (0..n) |j| for (0..n) |i| {
        const v = &a[i + j * n];
        if (i > j) {
            v.* = scalar(T, 0);
        } else if (i == j) {
            v.* = ops.add(T, scalar(T, @floatCast(lo + 1)), ops.mulReal(T, v.*, 0.5));
        } else v.* = ops.mulReal(T, v.*, s);
    };
    if (comptime !isComplex(T)) {
        var i: usize = 0;
        while (i + 1 < n) : (i += every) {
            // `[d, 0.7; -x, d]`, eigenvalues `d +- i sqrt(0.7 x)`
            a[i + 1 + i * n] = -0.5;
            a[i + (i + 1) * n] = 0.7;
            a[i + 1 + (i + 1) * n] = a[i + i * n];
        }
    }
    return a;
}
