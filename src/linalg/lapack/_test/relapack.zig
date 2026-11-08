//! Tests of `relapack.zig`'s routines: `potrf`, `getrf`, `trtri` and
//! `lauum` (`-Drelapack`), `sygst`, `hegst` and `pbtrf`.

const std = @import("std");
const types = @import("../../blas/types.zig");
const check = @import("../../blas/_test/check.zig");
const guard = @import("../../blas/_test/guard.zig");
const reference = @import("../../blas/_test/reference.zig");
const Vtable = @import("../../blas/Vtable.zig");
const options = @import("lapack_options");
const lcheck = @import("check.zig");
const fixtures = @import("fixtures.zig");
const symmetric = @import("symmetric.zig");
const ops = @import("../ops.zig");
const chol = @import("../cholesky.zig");
const lu = @import("../lu.zig");
const tri = @import("../triangular.zig");
const relapack = @import("../relapack.zig");

const testing = std.testing;
const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const isComplex = types.isComplex;
const scalar = ops.scalar;
const test_types = lcheck.test_types;
const isNanT = lcheck.isNanT;
const W = lcheck.W;

/// Around the crossover and `split`'s rounding for every element size.
pub const test_sizes = [_]usize{ 1, 2, 24, 25, 33, 49, 70, 100 };

fn potrfCases(vt: *const Vtable) anyerror!void {
    const al = testing.allocator;
    guard.use(vt);
    defer guard.restore();
    inline for (test_types) |T| for ([_]bool{ true, false }) |upper| for (test_sizes) |n| {
        if (check.rig() and n > 70) continue;
        var p = try fixtures.CholProblem(T).init(al, n, upper, n * 37 + @intFromBool(upper));
        defer p.deinit(al);
        relapack.potrf(al, check.io, T, lcheck.uploOf(upper), @intCast(n), p.a.ptr, @intCast(p.lda), null) catch |e| {
            check.say("relapack potrf {s} n={d} upper={}\n", .{ @typeName(T), n, upper });
            return e;
        };
        p.expectFactor() catch |e| {
            check.say("relapack potrf {s} n={d}\n", .{ @typeName(T), n });
            return e;
        };
    };
}

test "relapack potrf: A = L Lᴴ from a known L, both triangles, on every linked core" {
    try check.eachCore(potrfCases);
}

test "relapack potrf: the first minor that is not positive definite, either half" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex32 }) |T| for ([_]bool{ true, false }) |upper| {
        const n = 49;
        var p = try fixtures.CholProblem(T).init(al, n, upper, 5);
        defer p.deinit(al);
        for ([_]usize{ 0, 23, 24, 25, 40, 48 }) |k| {
            p.store();
            p.a[k * (p.lda + 1)] = scalar(T, -100);
            var minor: usize = 0;
            try testing.expectError(error.NotPositiveDefinite, relapack.potrf(al, check.io, T, lcheck.uploOf(upper), n, p.a.ptr, @intCast(p.lda), &minor));
            try testing.expectEqual(k + 1, minor);
        }
    };
}

const getrf_shapes = [_][2]usize{ .{ 1, 1 }, .{ 24, 24 }, .{ 25, 25 }, .{ 50, 50 }, .{ 70, 33 }, .{ 33, 70 }, .{ 100, 60 }, .{ 26, 100 } };

fn getrfCases(vt: *const Vtable) anyerror!void {
    const al = testing.allocator;
    guard.use(vt);
    defer guard.restore();
    inline for (test_types) |T| for (getrf_shapes) |s| {
        if (check.rig() and s[0] * s[1] > 3000) continue;
        var p = try fixtures.LuProblem(T).init(al, s[0], s[1], s[0] * 101 + s[1]);
        defer p.deinit(al);
        try relapack.getrf(al, check.io, T, @intCast(s[0]), @intCast(s[1]), p.a.ptr, @intCast(p.lda), p.ipiv.ptr, null);
        p.expectFactor(al) catch |e| {
            check.say("relapack getrf {s}\n", .{@typeName(T)});
            return e;
        };
    };
}

test "relapack getrf: P L U = A, square, tall and wide, on every linked core" {
    try check.eachCore(getrfCases);
}

test "relapack getrf: a zero column gives that pivot, the factorization completed, either half" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64 }) |T| for ([_]usize{ 0, 10, 30, 45 }) |k| {
        var p = try fixtures.LuProblem(T).init(al, 50, 50, 9);
        defer p.deinit(al);
        for (0..50) |i| p.aw[i + k * 50] = reference.wide(T, scalar(T, 0));
        p.store();
        var piv: usize = 0;
        try testing.expectError(error.Singular, relapack.getrf(al, check.io, T, 50, 50, p.a.ptr, @intCast(p.lda), p.ipiv.ptr, &piv));
        try testing.expectEqual(k + 1, piv);
        try p.expectFactor(al);
    };
}

fn triCases(vt: *const Vtable) anyerror!void {
    const al = testing.allocator;
    guard.use(vt);
    defer guard.restore();
    inline for (test_types) |T| for ([_]bool{ true, false }) |upper| for ([_]bool{ false, true }) |unit| for (test_sizes) |n| {
        if (check.rig() and n > 70) continue;
        var t = try fixtures.Tri(T).init(al, n, upper, unit, n * 13 + @intFromBool(upper));
        defer t.deinit(al);
        try relapack.trtri(al, check.io, T, lcheck.uploOf(upper), if (unit) .u else .n, @intCast(n), t.a.ptr, @intCast(t.lda), null);
        try t.expectInverse();
        if (unit) continue;
        t.store();
        try relapack.lauum(al, check.io, T, lcheck.uploOf(upper), @intCast(n), t.a.ptr, @intCast(t.lda));
        try t.expectProduct();
    };
}

test "relapack trtri and lauum: A X = I and U Uᴴ, Lᴴ L, both triangles, on every linked core" {
    try check.eachCore(triCases);
}

test "relapack trtri: a zero first diagonal element is singular at 1, nothing written" {
    const al = testing.allocator;
    var t = try fixtures.Tri(f64).init(al, 30, true, false, 3);
    defer t.deinit(al);
    t.a[0] = 0;
    const before = try al.dupe(f64, t.a);
    defer al.free(before);
    var piv: usize = 0;
    try testing.expectError(error.Singular, relapack.trtri(al, check.io, f64, .u, .n, 30, t.a.ptr, @intCast(t.lda), &piv));
    try testing.expectEqual(1, piv);
    for (before, t.a) |x, y| try testing.expect(lcheck.sameBits(f64, x, y));
}

test "relapack potrf, getrf, trtri, lauum: the first bad argument, nothing read when empty" {
    const al = testing.allocator;
    var a: [4]f64 = @splat(lcheck.poisoned(f64));
    var ipiv: [2]isize = @splat(-1);
    try testing.expectError(error.InvalidUplo, relapack.potrf(al, check.io, f64, .v, 2, &a, 2, null));
    try testing.expectError(error.InvalidN, relapack.potrf(al, check.io, f64, .u, -1, &a, 2, null));
    try testing.expectError(error.InvalidLda, relapack.potrf(al, check.io, f64, .u, 2, &a, 1, null));
    try testing.expectError(error.InvalidM, relapack.getrf(al, check.io, f64, -1, 2, &a, 2, &ipiv, null));
    try testing.expectError(error.InvalidN, relapack.getrf(al, check.io, f64, 2, -1, &a, 2, &ipiv, null));
    try testing.expectError(error.InvalidLda, relapack.getrf(al, check.io, f64, 2, 2, &a, 1, &ipiv, null));
    try testing.expectError(error.InvalidUplo, relapack.trtri(al, check.io, f64, .m, .n, 2, &a, 2, null));
    try testing.expectError(error.InvalidLda, relapack.lauum(al, check.io, f64, .l, 2, &a, 1));
    try relapack.potrf(al, check.io, f64, .u, 0, &a, 1, null);
    try relapack.getrf(al, check.io, f64, 0, 2, &a, 1, &ipiv, null);
    try relapack.trtri(al, check.io, f64, .l, .n, 0, &a, 1, null);
    try relapack.lauum(al, check.io, f64, .l, 0, &a, 1);
    for (a) |x| try testing.expect(std.math.isNan(x));
    try testing.expectEqual(-1, ipiv[0]);
}

test "with -Drelapack, ReLAPACK's potrf, getrf, trtri and lauum agree with OpenBLAS's" {
    if (!options.relapack) return error.SkipZigTest;
    const root = @import("../__root__.zig");
    const al = testing.allocator;
    const io = check.io;
    inline for (.{ f32, f64, types.Complex32, types.Complex64 }) |T| for ([_]bool{ true, false }) |upper| {
        const n = 100;
        const uplo = lcheck.uploOf(upper);
        var p = try fixtures.CholProblem(T).init(al, n, upper, 77);
        defer p.deinit(al);
        const x = try al.dupe(T, p.a);
        defer al.free(x);
        try root.potrf(al, io, T, uplo, n, x.ptr, @intCast(p.lda), null);
        try chol.potrf(al, io, T, uplo, n, p.a.ptr, @intCast(p.lda), null);
        try fixtures.expectAgree(T, "potrf", n, x, p.a);

        var t = try fixtures.Tri(T).init(al, n, upper, false, 78);
        defer t.deinit(al);
        const z = try al.dupe(T, t.a);
        defer al.free(z);
        try root.trtri(al, io, T, uplo, .n, n, z.ptr, @intCast(t.lda), null);
        try tri.trtri(al, io, T, uplo, .n, n, t.a.ptr, @intCast(t.lda), null);
        try fixtures.expectAgree(T, "trtri", n, z, t.a);
        t.store();
        @memcpy(z, t.a);
        try root.lauum(al, io, T, uplo, n, z.ptr, @intCast(t.lda));
        try tri.lauum(al, io, T, uplo, n, t.a.ptr, @intCast(t.lda));
        try fixtures.expectAgree(T, "lauum", n, z, t.a);

        var g = try fixtures.LuProblem(T).init(al, n, n, 79 + @as(u64, @intFromBool(upper)));
        defer g.deinit(al);
        const y = try al.dupe(T, g.a);
        defer al.free(y);
        const ip = try al.alloc(isize, n);
        defer al.free(ip);
        try root.getrf(al, io, T, n, n, y.ptr, @intCast(g.lda), ip.ptr, null);
        try lu.getrf(al, io, T, n, n, g.a.ptr, @intCast(g.lda), g.ipiv.ptr, null);
        try testing.expectEqualSlices(isize, g.ipiv, ip);
        try fixtures.expectAgree(T, "getrf", n, y, g.a);
    };
}

const Big = struct {
    a: []const f64,

    fn potrfRun(self: Big, io: Io, out: []f64) anyerror!void {
        @memcpy(out, self.a);
        try relapack.potrf(testing.allocator, io, f64, .l, lcheck.big_n, out.ptr, lcheck.big_n, null);
    }

    fn getrfRun(self: Big, io: Io, out: []f64) anyerror!void {
        @memcpy(out, self.a);
        var ipiv: [lcheck.big_n]isize = undefined;
        try relapack.getrf(testing.allocator, io, f64, lcheck.big_n, lcheck.big_n, out.ptr, lcheck.big_n, &ipiv, null);
    }
};

test "relapack potrf and getrf: the same bits inline and threaded" {
    const al = testing.allocator;
    const a = try fixtures.dominant(al, lcheck.big_n, 1);
    defer al.free(a);
    try lcheck.expectIoBits(f64, Big{ .a = a }, Big.potrfRun, lcheck.big_n * lcheck.big_n, "potrf");
    try lcheck.expectIoBits(f64, Big{ .a = a }, Big.getrfRun, lcheck.big_n * lcheck.big_n, "getrf");
}

// ---- sygst tests ------------------------------------------------------------------

fn gstCase(comptime T: type, al: Allocator, itype: isize, upper: bool, n: usize, seed: u64) !void {
    // `B` from a known Cholesky factor, `A` Hermitian, both in `uplo`'s triangle
    var p = try fixtures.CholProblem(T).init(al, n, upper, seed);
    defer p.deinit(al);
    var h = try symmetric.Indefinite(T, isComplex(T)).init(al, n, upper, false, seed + 1);
    defer h.deinit(al);
    const ldb = n + 1;
    const b = try al.alloc(T, @max(1, ldb * n));
    defer al.free(b);
    @memset(b, lcheck.poisoned(T));
    // `L`, or `U = Lᴴ`, whole in `W(T)`
    const f = try al.alloc(W(T), n * n);
    defer al.free(f);
    for (0..n) |j| for (0..n) |i| {
        const l = reference.wide(T, p.l[i + j * n]);
        f[i + j * n] = if (upper) reference.conjW(reference.wide(T, p.l[j + i * n])) else l;
        if (if (upper) i <= j else i >= j) b[i + j * ldb] = lcheck.toT(T, f[i + j * n]);
    };
    const routine = if (comptime isComplex(T)) relapack.hegst else relapack.sygst;
    try routine(al, check.io, T, itype, lcheck.uploOf(upper), @intCast(n), h.a.ptr, @intCast(h.lda), b.ptr, @intCast(ldb));
    // the result, whole
    const c = try al.alloc(W(T), n * n);
    defer al.free(c);
    for (0..n) |j| for (0..n) |i| {
        const s = if (upper) i <= j else i >= j;
        const v = reference.wide(T, if (s) h.a[i + j * h.lda] else h.a[j + i * h.lda]);
        c[i + j * n] = if (s) v else reference.conjW(v);
    };
    // `M X Mᴴ` against the other side: `L C Lᴴ` or `Uᴴ C U` against `A`
    // (itype 1); `C` against `Lᴴ A L` or `U A Uᴴ`
    const m = try al.alloc(W(T), n * n);
    defer al.free(m);
    for (0..n) |j| for (0..n) |i| {
        const flip = (itype == 1) == upper;
        m[i + j * n] = if (flip) reference.conjW(f[j + i * n]) else f[i + j * n];
    };
    const r = try al.alloc(W(T), n * n);
    defer al.free(r);
    const mag = try al.alloc(f128, n * n);
    defer al.free(mag);
    try symmetric.sandwich(T, al, n, m, if (itype == 1) c else h.w, r, mag);
    const other = if (itype == 1) h.w else c;
    const eps: f128 = std.math.floatEps(Float(T));
    const ulps: f128 = @floatFromInt(32 * (n + 2));
    for (0..n) |j| for (0..n) |i| {
        const d = reference.abs1(reference.subW(r[i + j * n], other[i + j * n]));
        const allowed = ulps * eps * (mag[i + j * n] + reference.abs1(other[i + j * n])) + std.math.floatMin(Float(T));
        if (!(d <= allowed)) {
            check.say("sygst {s} itype {d} {s} n={d}: ({d},{d}) off by {e}, allowed {e}\n", .{ @typeName(T), itype, if (upper) "u" else "l", n, i, j, @as(f64, @floatCast(d)), @as(f64, @floatCast(allowed)) });
            return error.TestUnexpectedResult;
        }
    };
    for (0..n) |j| for (0..h.lda) |i| {
        const x = h.a[i + j * h.lda];
        const inside = i < n and (if (upper) i <= j else i >= j);
        if (!inside and !isNanT(T, x)) {
            check.say("sygst {s} n={d}: wrote ({d},{d}), outside its triangle\n", .{ @typeName(T), n, i, j });
            return error.TestUnexpectedResult;
        }
        if (comptime isComplex(T)) if (i == j and i < n and x.im != 0) {
            check.say("hegst {s} n={d}: ({d},{d}) has imaginary part {e}\n", .{ @typeName(T), n, i, i, x.im });
            return error.TestUnexpectedResult;
        };
    };
}

test "relapack sygst, hegst: every itype, both triangles, against the products, every type" {
    const al = testing.allocator;
    inline for (test_types) |T| for ([_]isize{ 1, 2, 3 }) |itype| for ([_]bool{ true, false }) |upper| {
        for ([_]usize{ 1, 2, 24, 25, 40, 57 }) |n| {
            if (check.rig() and n > 40) continue;
            gstCase(T, al, itype, upper, n, n * 5 + @as(u64, @intCast(itype))) catch |e| {
                check.say("sygst {s} itype {d} n={d}\n", .{ @typeName(T), itype, n });
                return e;
            };
        }
    };
}

test "relapack sygst, hegst: past the crossover twice over, f64 and c64" {
    if (check.rig()) return error.SkipZigTest;
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64 }) |T| for ([_]isize{ 1, 2 }) |itype| for ([_]bool{ true, false }) |upper| {
        try gstCase(T, al, itype, upper, 100, 7);
    };
}

test "relapack sygst, hegst: the first bad argument, nothing read for n = 0" {
    const al = testing.allocator;
    var a: [4]f64 = @splat(lcheck.poisoned(f64));
    const b: [4]f64 = @splat(lcheck.poisoned(f64));
    try testing.expectError(error.InvalidItype, relapack.sygst(al, check.io, f64, 0, .v, 2, &a, 2, &b, 2));
    try testing.expectError(error.InvalidItype, relapack.sygst(al, check.io, f64, 4, .u, 2, &a, 2, &b, 2));
    try testing.expectError(error.InvalidUplo, relapack.sygst(al, check.io, f64, 1, .m, 2, &a, 2, &b, 2));
    try testing.expectError(error.InvalidN, relapack.sygst(al, check.io, f64, 1, .u, -1, &a, 2, &b, 2));
    try testing.expectError(error.InvalidLda, relapack.sygst(al, check.io, f64, 2, .u, 2, &a, 1, &b, 2));
    try testing.expectError(error.InvalidLdb, relapack.sygst(al, check.io, f64, 3, .l, 2, &a, 2, &b, 1));
    try relapack.sygst(al, check.io, f64, 1, .l, 0, &a, 1, &b, 1);
    for (a) |x| try testing.expect(std.math.isNan(x));
    var z: [1]types.Complex32 = @splat(lcheck.poisoned(types.Complex32));
    try testing.expectError(error.InvalidLdb, relapack.hegst(al, check.io, types.Complex32, 1, .u, 1, &z, 1, &z, 0));
}

const BigGst = struct {
    a: []const f64,
    b: []const f64,

    fn run(self: BigGst, io: Io, out: []f64) anyerror!void {
        @memcpy(out, self.a);
        try relapack.sygst(testing.allocator, io, f64, 1, .l, lcheck.big_n, out.ptr, lcheck.big_n, self.b.ptr, lcheck.big_n);
    }
};

test "relapack sygst: the same bits inline and threaded" {
    const al = testing.allocator;
    const a = try fixtures.dominant(al, lcheck.big_n, 5);
    defer al.free(a);
    const b = try fixtures.dominant(al, lcheck.big_n, 6);
    defer al.free(b);
    try relapack.potrf(al, check.io, f64, .l, lcheck.big_n, b.ptr, lcheck.big_n, null);
    try lcheck.expectIoBits(f64, BigGst{ .a = a, .b = b }, BigGst.run, lcheck.big_n * lcheck.big_n, "sygst");
}

// ---- pbtrf tests ------------------------------------------------------------------

/// A Hermitian (real: symmetric) positive-definite band matrix `A = L Lᴴ`,
/// `L` lower with `kd` subdiagonals, its diagonal in [2, 4) and the rest
/// within `1 / (kd + 1)`; stored in `upper`'s band form at `ldab = kd + 2`,
/// every unused place NaN.
fn Band(comptime T: type) type {
    return struct {
        n: usize,
        kd: usize,
        ldab: usize,
        upper: bool,
        ab: []T,
        aw: []W(T),

        const Self = @This();

        fn init(al: Allocator, n: usize, kd: usize, upper: bool, seed: u64) !Self {
            const v = try al.alloc(T, n * n);
            defer al.free(v);
            check.fill(T, v, seed);
            const scale: Float(T) = @floatCast(1 / @as(f64, @floatFromInt(kd + 1)));
            const l = try al.alloc(W(T), n * n);
            defer al.free(l);
            for (0..n) |j| for (0..n) |i| {
                const x = v[i + j * n];
                l[i + j * n] = reference.wide(T, if (i < j or i > j + kd) scalar(T, 0) else if (i == j) scalar(T, 3 + ops.re(T, x) / 2) else ops.mulReal(T, x, scale));
            };
            const aw = try al.alloc(W(T), n * n);
            errdefer al.free(aw);
            for (0..n) |j| for (0..n) |i| {
                var s = reference.wide(T, scalar(T, 0));
                for (0..@min(i, j) + 1) |k| s = reference.addW(s, reference.mulW(l[i + k * n], reference.conjW(l[j + k * n])));
                aw[i + j * n] = s;
            };
            const ldab = kd + 2;
            const ab = try al.alloc(T, @max(1, ldab * n));
            var self: Self = .{ .n = n, .kd = kd, .ldab = ldab, .upper = upper, .ab = ab, .aw = aw };
            self.store();
            return self;
        }

        fn deinit(self: *Self, al: Allocator) void {
            al.free(self.ab);
            al.free(self.aw);
        }

        /// The storage place of `(i, j)`, `i >= j` (`l`) or `i <= j` (`u`),
        /// inside the band, else null.
        fn place(self: *const Self, i: usize, j: usize) ?usize {
            if (self.upper) {
                if (i > j or j - i > self.kd) return null;
                return self.kd + i - j + j * self.ldab;
            }
            if (i < j or i - j > self.kd) return null;
            return i - j + j * self.ldab;
        }

        fn store(self: *Self) void {
            @memset(self.ab, lcheck.poisoned(T));
            for (0..self.n) |j| for (0..self.n) |i| if (self.place(i, j)) |at| {
                self.ab[at] = lcheck.toT(T, self.aw[i + j * self.n]);
            };
        }

        /// The factor's `(i, j)` as lower, `i >= j`, zero outside the band.
        fn got(self: *const Self, i: usize, j: usize) W(T) {
            const at = (if (self.upper) self.place(j, i) else self.place(i, j)) orelse return reference.wide(T, scalar(T, 0));
            const v = reference.wide(T, self.ab[at]);
            return if (self.upper) reference.conjW(v) else v;
        }

        /// `F Fᴴ` against `A` within `2 (n + 2)` ulps of `|F| |F|ᴴ`, the
        /// unused places untouched.
        fn expectFactor(self: *const Self) !void {
            const n = self.n;
            const eps: f128 = std.math.floatEps(Float(T));
            const ulps: f128 = @floatFromInt(2 * (n + 2));
            for (0..n) |j| for (j..n) |i| {
                var s = reference.wide(T, scalar(T, 0));
                var mag: f128 = 0;
                const lo = @min(if (i > self.kd) i - self.kd else 0, j + 1);
                for (lo..j + 1) |k| {
                    const x = self.got(i, k);
                    const y = reference.conjW(self.got(j, k));
                    s = reference.addW(s, reference.mulW(x, y));
                    mag += reference.abs1(x) * reference.abs1(y);
                }
                const d = reference.abs1(reference.subW(s, self.aw[i + j * n]));
                if (!(d <= ulps * eps * mag + std.math.floatMin(Float(T)))) {
                    check.say("pbtrf {s} {s} n={d} kd={d}: (F Fᴴ)({d},{d}) off by {e}\n", .{ @typeName(T), if (self.upper) "u" else "l", n, self.kd, i, j, @as(f64, @floatCast(d)) });
                    return error.TestUnexpectedResult;
                }
            };
            for (0..n) |j| for (0..self.ldab) |r| {
                const used = r <= self.kd and (if (self.upper) r + j >= self.kd else r + j < n);
                if (!used and !isNanT(T, self.ab[r + j * self.ldab])) {
                    check.say("pbtrf {s} n={d} kd={d}: wrote unused ({d},{d})\n", .{ @typeName(T), n, self.kd, r, j });
                    return error.TestUnexpectedResult;
                }
            };
        }
    };
}

const band_cases = [_][2]usize{ .{ 1, 0 }, .{ 5, 0 }, .{ 30, 1 }, .{ 30, 5 }, .{ 26, 25 }, .{ 50, 10 }, .{ 60, 30 }, .{ 40, 60 } };

test "relapack pbtrf: A = F Fᴴ in band storage, both triangles, every type" {
    const al = testing.allocator;
    inline for (test_types) |T| for ([_]bool{ true, false }) |upper| for (band_cases) |c| {
        var p = try Band(T).init(al, c[0], c[1], upper, c[0] * 3 + c[1]);
        defer p.deinit(al);
        relapack.pbtrf(al, check.io, T, lcheck.uploOf(upper), @intCast(c[0]), @intCast(c[1]), p.ab.ptr, @intCast(p.ldab), null) catch |e| {
            check.say("pbtrf {s} n={d} kd={d}\n", .{ @typeName(T), c[0], c[1] });
            return e;
        };
        try p.expectFactor();
    };
}

test "relapack pbtrf: a rest wider than the band past the block (ReLAPACK's dense case), f64 and c64" {
    if (check.rig()) return error.SkipZigTest;
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64 }) |T| for ([_]bool{ true, false }) |upper| for ([_][2]usize{ .{ 100, 50 }, .{ 100, 60 }, .{ 120, 30 } }) |c| {
        var p = try Band(T).init(al, c[0], c[1], upper, 11);
        defer p.deinit(al);
        try relapack.pbtrf(al, check.io, T, lcheck.uploOf(upper), @intCast(c[0]), @intCast(c[1]), p.ab.ptr, @intCast(p.ldab), null);
        try p.expectFactor();
    };
}

test "relapack pbtrf: the first minor that is not positive definite, NaN too" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex32 }) |T| for ([_]bool{ true, false }) |upper| {
        var p = try Band(T).init(al, 60, 10, upper, 2);
        defer p.deinit(al);
        for ([_]usize{ 0, 23, 30, 59 }) |k| for ([_]f64{ -100, std.math.nan(f64) }) |x| {
            p.store();
            p.ab[p.place(k, k).?] = lcheck.fromF64(T, x);
            var minor: usize = 0;
            try testing.expectError(error.NotPositiveDefinite, relapack.pbtrf(al, check.io, T, lcheck.uploOf(upper), 60, 10, p.ab.ptr, @intCast(p.ldab), &minor));
            try testing.expectEqual(k + 1, minor);
        };
    };
}

test "relapack pbtrf: the first bad argument, nothing read for n = 0" {
    const al = testing.allocator;
    var ab: [6]f64 = @splat(lcheck.poisoned(f64));
    try testing.expectError(error.InvalidUplo, relapack.pbtrf(al, check.io, f64, .v, 2, 1, &ab, 2, null));
    try testing.expectError(error.InvalidN, relapack.pbtrf(al, check.io, f64, .u, -1, 1, &ab, 2, null));
    try testing.expectError(error.InvalidKd, relapack.pbtrf(al, check.io, f64, .u, 2, -1, &ab, 2, null));
    try testing.expectError(error.InvalidLdab, relapack.pbtrf(al, check.io, f64, .l, 2, 1, &ab, 1, null));
    try relapack.pbtrf(al, check.io, f64, .l, 0, 1, &ab, 2, null);
    for (ab) |x| try testing.expect(std.math.isNan(x));
}
