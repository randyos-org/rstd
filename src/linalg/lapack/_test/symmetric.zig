//! Symmetric and Hermitian test problems: `sytrf`'s `Indefinite` and its
//! `ldlProduct`, and `sygst`'s `sandwich`.

const std = @import("std");
const types = @import("../../blas/types.zig");
const check = @import("../../blas/_test/check.zig");
const reference = @import("../../blas/_test/reference.zig");
const ops = @import("../ops.zig");
const lcheck = @import("check.zig");

const Allocator = std.mem.Allocator;
const Float = types.Float;
const scalar = ops.scalar;
const W = lcheck.W;
const toT = lcheck.toT;
const poisoned = lcheck.poisoned;
const isNanT = lcheck.isNanT;

/// `sytrf`'s factor in `a` (`n x n`, `lda`) and `ipiv` multiplied out in
/// `W(T)`: `L D Lᵀ` or `U D Uᵀ` (`ᴴ` with `herm`), as LAPACK's `xLAVSY` and
/// `xLAVHE` read the format, applied to the identity; `mag` gets the same
/// products in magnitude (`|L| |D| |L|ᵀ`, permuted).
pub fn ldlProduct(comptime T: type, comptime herm: bool, comptime rook: bool, al: Allocator, upper: bool, n: usize, a: []const T, lda: usize, ipiv: []const isize, out: []W(T), mag: []f128) !void {
    // `P = S(p1) S(p2)` per step, in the order the factorization took them
    const Step = struct { b: usize, s: usize, p1: [2]usize, p2: ?[2]usize };
    var steps: std.ArrayList(Step) = .empty;
    defer steps.deinit(al);
    const at = struct {
        fn neg(v: isize) usize {
            return @intCast(-v - 1);
        }
    };
    if (!upper) {
        var k: usize = 0;
        while (k < n) {
            if (ipiv[k] > 0) {
                try steps.append(al, .{ .b = k, .s = 1, .p1 = .{ k, @intCast(ipiv[k] - 1) }, .p2 = null });
                k += 1;
            } else {
                try steps.append(al, if (rook)
                    .{ .b = k, .s = 2, .p1 = .{ k, at.neg(ipiv[k]) }, .p2 = .{ k + 1, at.neg(ipiv[k + 1]) } }
                else
                    .{ .b = k, .s = 2, .p1 = .{ k + 1, at.neg(ipiv[k]) }, .p2 = null });
                k += 2;
            }
        }
    } else {
        var k: usize = n;
        while (k > 0) {
            const j = k - 1;
            if (ipiv[j] > 0) {
                try steps.append(al, .{ .b = j, .s = 1, .p1 = .{ j, @intCast(ipiv[j] - 1) }, .p2 = null });
                k -= 1;
            } else {
                try steps.append(al, if (rook)
                    .{ .b = j - 1, .s = 2, .p1 = .{ j, at.neg(ipiv[j]) }, .p2 = .{ j - 1, at.neg(ipiv[j - 1]) } }
                else
                    .{ .b = j - 1, .s = 2, .p1 = .{ j - 1, at.neg(ipiv[j]) }, .p2 = null });
                k -= 2;
            }
        }
    }
    const F = reference.Wide(T);
    const zero = reference.wide(T, scalar(T, 0));
    const Ctx = struct {
        n: usize,
        x: []F,
        m: []f128,
        fn swap(c: @This(), r: usize, s: usize) void {
            if (r == s) return;
            for (0..c.n) |j| {
                std.mem.swap(F, &c.x[r + j * c.n], &c.x[s + j * c.n]);
                std.mem.swap(f128, &c.m[r + j * c.n], &c.m[s + j * c.n]);
            }
        }
    };
    for (out, mag, 0..) |*o, *g, i| {
        const one = i % (n + 1) == 0;
        o.* = if (one) reference.wide(T, scalar(T, 1)) else zero;
        g.* = if (one) 1 else 0;
    }
    const c: Ctx = .{ .n = n, .x = out, .m = mag };
    const el = struct {
        fn of(v: T) F {
            return reference.wide(T, v);
        }
        fn cjw(v: F) F {
            return if (comptime herm) reference.conjW(v) else v;
        }
    };
    // the multiplier rows of step `st`
    const Rows = struct {
        fn lo(st: Step, up: bool) usize {
            return if (up) 0 else st.b + st.s;
        }
        fn hi(st: Step, up: bool, nn: usize) usize {
            return if (up) st.b else nn;
        }
    };
    // X := Lᵀ X (`ᴴ`): each step's `Pᵀ` (p1, then p2), then its multipliers
    for (steps.items) |st| {
        c.swap(st.p1[0], st.p1[1]);
        if (st.p2) |q| c.swap(q[0], q[1]);
        for (st.b..st.b + st.s) |col| for (Rows.lo(st, upper)..Rows.hi(st, upper, n)) |r| {
            const l = el.cjw(el.of(a[r + col * lda]));
            for (0..n) |j| {
                c.x[col + j * n] = reference.addW(c.x[col + j * n], reference.mulW(l, c.x[r + j * n]));
                c.m[col + j * n] += reference.abs1(l) * c.m[r + j * n];
            }
        };
    }
    // X := D X
    for (steps.items) |st| {
        const b = st.b;
        if (st.s == 1) {
            const d = el.of(a[b + b * lda]);
            const dd = if (comptime herm) reference.wide(T, scalar(T, a[b + b * lda].re)) else d;
            for (0..n) |j| {
                c.x[b + j * n] = reference.mulW(dd, c.x[b + j * n]);
                c.m[b + j * n] *= reference.abs1(dd);
            }
            continue;
        }
        const d11 = if (comptime herm) reference.wide(T, scalar(T, a[b + b * lda].re)) else el.of(a[b + b * lda]);
        const d22 = if (comptime herm) reference.wide(T, scalar(T, a[b + 1 + (b + 1) * lda].re)) else el.of(a[b + 1 + (b + 1) * lda]);
        const d21 = if (upper) el.cjw(el.of(a[b + (b + 1) * lda])) else el.of(a[b + 1 + b * lda]);
        const d12 = el.cjw(d21);
        for (0..n) |j| {
            const x0 = c.x[b + j * n];
            const x1 = c.x[b + 1 + j * n];
            const m0 = c.m[b + j * n];
            const m1 = c.m[b + 1 + j * n];
            c.x[b + j * n] = reference.addW(reference.mulW(d11, x0), reference.mulW(d12, x1));
            c.x[b + 1 + j * n] = reference.addW(reference.mulW(d21, x0), reference.mulW(d22, x1));
            c.m[b + j * n] = reference.abs1(d11) * m0 + reference.abs1(d12) * m1;
            c.m[b + 1 + j * n] = reference.abs1(d21) * m0 + reference.abs1(d22) * m1;
        }
    }
    // X := L X: steps in reverse, multipliers then `P` (p2, then p1)
    var i = steps.items.len;
    while (i > 0) {
        i -= 1;
        const st = steps.items[i];
        for (st.b..st.b + st.s) |col| for (Rows.lo(st, upper)..Rows.hi(st, upper, n)) |r| {
            const l = el.of(a[r + col * lda]);
            for (0..n) |j| {
                c.x[r + j * n] = reference.addW(c.x[r + j * n], reference.mulW(l, c.x[col + j * n]));
                c.m[r + j * n] += reference.abs1(l) * c.m[col + j * n];
            }
        };
        if (st.p2) |q| c.swap(q[0], q[1]);
        c.swap(st.p1[0], st.p1[1]);
    }
}

/// A random symmetric (`herm`: Hermitian) `n x n` in its `upper` or lower
/// triangle at `lda = n + 2`, `zero_diag` with a zero diagonal (2 x 2
/// pivots); the other triangle, the padding and a Hermitian diagonal's
/// imaginary parts NaN. `w` keeps the whole matrix.
pub fn Indefinite(comptime T: type, comptime herm: bool) type {
    return struct {
        n: usize,
        lda: usize,
        upper: bool,
        a: []T,
        w: []W(T),
        ipiv: []isize,

        const Self = @This();

        pub fn init(al: Allocator, n: usize, upper: bool, zero_diag: bool, seed: u64) !Self {
            const lda = n + 2;
            const v = try al.alloc(T, n * n);
            defer al.free(v);
            check.fill(T, v, seed);
            const w = try al.alloc(W(T), n * n);
            errdefer al.free(w);
            for (0..n) |j| for (0..n) |i| {
                var x = reference.wide(T, v[@max(i, j) + @min(i, j) * n]);
                if (i < j and herm) x = reference.conjW(x);
                if (i == j) x = if (zero_diag) reference.wide(T, scalar(T, 0)) else if (comptime herm) reference.wide(T, scalar(T, v[i + i * n].re)) else x;
                w[i + j * n] = x;
            };
            const a = try al.alloc(T, @max(1, lda * n));
            errdefer al.free(a);
            const ipiv = try al.alloc(isize, @max(1, n));
            var self: Self = .{ .n = n, .lda = lda, .upper = upper, .a = a, .w = w, .ipiv = ipiv };
            self.store();
            return self;
        }

        pub fn store(self: *Self) void {
            @memset(self.a, poisoned(T));
            for (0..self.n) |j| for (0..self.n) |i| {
                if (if (self.upper) i > j else i < j) continue;
                var x = toT(T, self.w[i + j * self.n]);
                if (comptime herm) if (i == j) {
                    x.im = std.math.nan(Float(T));
                };
                self.a[i + j * self.lda] = x;
            };
            @memset(self.ipiv, 0);
        }

        pub fn deinit(self: *Self, al: Allocator) void {
            al.free(self.a);
            al.free(self.w);
            al.free(self.ipiv);
        }

        /// The factor against `A` within `32 (n + 2)` ulps of `|A| + |L| |D|
        /// |L|ᵀ`, a Hermitian `D`'s diagonal real, nothing outside the
        /// triangle written.
        pub fn expectFactor(self: *const Self, comptime rook: bool, al: Allocator, what: []const u8) !void {
            const n = self.n;
            const got = try al.alloc(W(T), n * n);
            defer al.free(got);
            const mag = try al.alloc(f128, n * n);
            defer al.free(mag);
            try ldlProduct(T, herm, rook, al, self.upper, n, self.a, self.lda, self.ipiv[0..n], got, mag);
            const eps: f128 = std.math.floatEps(Float(T));
            const ulps: f128 = @floatFromInt(32 * (n + 2));
            for (0..n) |j| for (0..n) |i| {
                const want = self.w[i + j * n];
                const d = reference.abs1(reference.subW(got[i + j * n], want));
                const allowed = ulps * eps * (reference.abs1(want) + mag[i + j * n]) + std.math.floatMin(Float(T));
                if (!(d <= allowed)) {
                    check.say("{s} {s} {s} n={d}: (L D Lᵀ)({d},{d}) off by {e}, allowed {e}\n", .{ what, @typeName(T), if (self.upper) "u" else "l", n, i, j, @as(f64, @floatCast(d)), @as(f64, @floatCast(allowed)) });
                    return error.TestUnexpectedResult;
                }
            };
            for (0..n) |j| for (0..self.lda) |i| {
                const x = self.a[i + j * self.lda];
                const inside = i < n and (if (self.upper) i <= j else i >= j);
                if (!inside and !isNanT(T, x)) {
                    check.say("{s} {s} n={d}: wrote ({d},{d}), outside its triangle\n", .{ what, @typeName(T), n, i, j });
                    return error.TestUnexpectedResult;
                }
                if (comptime herm) if (i == j and i < n and x.im != 0) {
                    check.say("{s} {s} n={d}: D({d},{d}) has imaginary part {e}\n", .{ what, @typeName(T), n, i, i, x.im });
                    return error.TestUnexpectedResult;
                };
            };
        }
    };
}

/// `M X Mᴴ` of `n x n` wide matrices, and its magnitude `|M| |X| |M|ᴴ`.
pub fn sandwich(comptime T: type, al: Allocator, n: usize, m: []const W(T), x: []const W(T), out: []W(T), mag: []f128) !void {
    const t = try al.alloc(W(T), n * n);
    defer al.free(t);
    const tm = try al.alloc(f128, n * n);
    defer al.free(tm);
    const zero = reference.wide(T, scalar(T, 0));
    for (0..n) |j| for (0..n) |i| {
        var s = zero;
        var g: f128 = 0;
        for (0..n) |k| {
            s = reference.addW(s, reference.mulW(m[i + k * n], x[k + j * n]));
            g += reference.abs1(m[i + k * n]) * reference.abs1(x[k + j * n]);
        }
        t[i + j * n] = s;
        tm[i + j * n] = g;
    };
    for (0..n) |j| for (0..n) |i| {
        var s = zero;
        var g: f128 = 0;
        for (0..n) |k| {
            const c = reference.conjW(m[j + k * n]);
            s = reference.addW(s, reference.mulW(t[i + k * n], c));
            g += tm[i + k * n] * reference.abs1(c);
        }
        out[i + j * n] = s;
        mag[i + j * n] = g;
    };
}
