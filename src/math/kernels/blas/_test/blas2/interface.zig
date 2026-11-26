//! `blas2/interface.zig`'s tests: `gemv`, `ger`, `symv`, `hemv` and the half-precision gemvs.

const std = @import("std");
const blas = @import("../../__root__.zig");
const types = @import("../../types.zig");
const util = @import("../../util.zig");
const sym = @import("../../blas2/symmetric.zig");
const portable = @import("../../portable.zig");
const bf = @import("../../ext/bf16/__root__.zig");
const interface = @import("../../blas2/interface.zig");

const Trans = types.Trans;
const Uplo = types.Uplo;
const isComplex = types.isComplex;
const Float = types.Float;
const isZero = util.isZero;
const gemv = interface.gemv;
const HalfPrecision = interface.HalfPrecision;
const gemvHalf = interface.gemvHalf;
const ger = interface.ger;
const gerc = interface.gerc;
const symv = interface.symv;
const hemv = interface.hemv;

const Vtable = @import("../../Vtable.zig");
const check = @import("../check.zig");
const guard = @import("../guard.zig");
const reference = @import("../reference.zig");
const Complex32 = types.Complex32;
const Complex64 = types.Complex64;

const tri = @import("../../blas2/triangular.zig");
const band = @import("../../blas2/banded.zig");
const testing = std.testing;
/// Scratch from the heap holds NaN, so a missing clear shows.
var poison: guard.Poison = .{ .child = testing.allocator };
const all = check.all_types;

/// Points the global vtable at the portable kernels until
/// `defer blas.testRestore(saved)`.
fn usePortable() blas.TestSaved {
    var vt: Vtable = .{};
    portable.install(&vt);
    return blas.testUse(&vt);
}

fn mulT(comptime T: type, p: T, q: T) T {
    return if (comptime isComplex(T)) p.mul(q) else p * q;
}

pub fn addT(comptime T: type, p: T, q: T) T {
    return if (comptime isComplex(T)) p.add(q) else p + q;
}

fn conjT(comptime T: type, p: T) T {
    return if (comptime isComplex(T)) p.conjugate() else p;
}

pub fn zeroT(comptime T: type) T {
    return if (comptime isComplex(T)) .init(0, 0) else 0;
}

pub fn val(comptime T: type, re: Float(T), im: Float(T)) T {
    if (comptime isComplex(T)) return .init(re, im);
    // an unused f80 argument leaks an x87 register on the x86_64 backend (issue #272)
    std.mem.doNotOptimizeAway(im);
    return re;
}

/// Element `i` of a vector of `n` at stride `inc`, as an index into its
/// storage: stride -2 puts element 0 last.
fn vi(n: usize, inc: isize, i: usize) usize {
    const k: usize = @abs(inc);
    return if (inc < 0) (n - 1 - i) * k else i * k;
}

fn randomVec(comptime T: type, n: usize, inc: isize, seed: u64) ![]T {
    const v = try testing.allocator.alloc(T, guard.span(n, @abs(inc)));
    check.fill(T, v, seed);
    return v;
}

fn expectSame(comptime T: type, want: []const T, got: []const T, terms: usize) !void {
    try check.expectCloseSlices(T, want, got, check.tolerance(Float(T), 8 * (terms + 2)));
}

const Case = struct { m: usize, n: usize };
const shapes = [_]Case{ .{ .m = 1, .n = 1 }, .{ .m = 5, .n = 3 }, .{ .m = 3, .n = 5 }, .{ .m = 9, .n = 7 }, .{ .m = 17, .n = 4 } };
const stride_pairs = [_][2]isize{ .{ 1, 1 }, .{ 2, -3 }, .{ -1, 2 }, .{ -2, -1 } };

/// `alpha`, `beta` pairs: general, and each special case.
fn scalars(comptime T: type) [6][2]T {
    return .{
        .{ val(T, 1.25, -0.5), val(T, 0.5, 0.25) },
        .{ val(T, 1, 0), val(T, 1, 0) },
        .{ val(T, 0, 0), val(T, 0.5, 0) },
        .{ val(T, 1.5, 0.5), val(T, 0, 0) },
        .{ val(T, 0, 0), val(T, 1, 0) },
        .{ val(T, -1, 0), val(T, 2, 0) },
    };
}

pub fn testGemv(comptime T: type, trans: Trans, m: usize, n: usize, incx: isize, incy: isize, alpha: T, beta: T) !void {
    const tr = trans.transposed();
    const lenx = if (tr) m else n;
    const leny = if (tr) n else m;
    const lda = m + 2;
    const a = try testing.allocator.alloc(T, lda * n);
    defer testing.allocator.free(a);
    check.fill(T, a, m * 7 + n);
    const x = try randomVec(T, lenx, incx, 3);
    defer testing.allocator.free(x);
    const y = try randomVec(T, leny, incy, 4);
    defer testing.allocator.free(y);
    const want = try testing.allocator.dupe(T, y);
    defer testing.allocator.free(want);
    const conj_a = trans.conjA();
    const conj_x = trans.conjX();
    for (0..leny) |r| {
        var acc = zeroT(T);
        for (0..lenx) |q| {
            var e = if (tr) a[q + r * lda] else a[r + q * lda];
            var xv = x[vi(lenx, incx, q)];
            if (conj_a) e = conjT(T, e);
            if (conj_x) xv = conjT(T, xv);
            acc = addT(T, acc, mulT(T, e, xv));
        }
        const k = vi(leny, incy, r);
        const scaled = if (isZero(T, beta)) zeroT(T) else mulT(T, beta, y[k]);
        want[k] = if (isZero(T, alpha)) scaled else addT(T, scaled, mulT(T, alpha, acc));
    }
    try gemv(poison.allocator(), check.io, T, trans, @intCast(m), @intCast(n), alpha, a.ptr, @intCast(lda), x.ptr, incx, beta, y.ptr, incy);
    try expectSame(T, want, y, lenx);
}

test "interface: gemv matches plain arithmetic, for every trans, stride and scalar case" {
    const saved = usePortable();
    defer blas.testRestore(saved);
    inline for (all) |T| {
        const transes = if (comptime isComplex(T)) std.enums.values(Trans) else &[_]Trans{ .n, .t, .r, .c };
        for (transes) |trans| for (shapes) |s| for (stride_pairs) |inc| for (scalars(T)) |ab| {
            testGemv(T, trans, s.m, s.n, inc[0], inc[1], ab[0], ab[1]) catch |e| {
                check.say("gemv {s} {s} m={d} n={d} inc={d},{d}\n", .{ @typeName(T), @tagName(trans), s.m, s.n, inc[0], inc[1] });
                return e;
            };
        };
    }
}

test "interface: gemv beta 0 overwrites NaN in y, and alpha 0 never reads x or a" {
    const saved = usePortable();
    defer blas.testRestore(saved);
    inline for (all) |T| {
        const nan = val(T, std.math.nan(Float(T)), std.math.nan(Float(T)));
        var a = [_]T{ nan, nan, nan, nan };
        var x = [_]T{ nan, nan };
        var y = [_]T{ nan, nan };
        try gemv(testing.allocator, check.io, T, .n, 2, 2, val(T, 0, 0), &a, 2, &x, 1, val(T, 0, 0), &y, 1);
        try expectSame(T, &[_]T{ zeroT(T), zeroT(T) }, &y, 1);
        y = .{ val(T, 1, 0), val(T, 2, 0) };
        try gemv(testing.allocator, check.io, T, .t, 2, 2, val(T, 0, 0), &a, 2, &x, 1, val(T, 2, 0), &y, 1);
        try expectSame(T, &[_]T{ val(T, 2, 0), val(T, 4, 0) }, &y, 1);
    }
}

test "interface: gemv, ger, symv and hemv report the first bad argument and do nothing for empty problems" {
    const saved = usePortable();
    defer blas.testRestore(saved);
    const al = testing.allocator;
    var a = [_]f64{ 1, 2, 3, 4 };
    var x = [_]f64{ 1, 1 };
    var y = [_]f64{ 5, 6 };
    try testing.expectError(error.InvalidTrans, gemv(al, check.io, f64, .o, 2, 2, 1, &a, 2, &x, 1, 1, &y, 1));
    try testing.expectError(error.InvalidM, gemv(al, check.io, f64, .n, -1, 2, 1, &a, 2, &x, 1, 1, &y, 1));
    try testing.expectError(error.InvalidN, gemv(al, check.io, f64, .n, 2, -1, 1, &a, 2, &x, 1, 1, &y, 1));
    try testing.expectError(error.InvalidLda, gemv(al, check.io, f64, .n, 2, 2, 1, &a, 1, &x, 1, 1, &y, 1));
    try testing.expectError(error.InvalidIncx, gemv(al, check.io, f64, .n, 2, 2, 1, &a, 2, &x, 0, 1, &y, 1));
    try testing.expectError(error.InvalidIncy, gemv(al, check.io, f64, .n, 2, 2, 1, &a, 2, &x, 1, 1, &y, 0));
    // m wins over a bad stride
    try testing.expectError(error.InvalidM, gemv(al, check.io, f64, .n, -1, 2, 1, &a, 2, &x, 0, 1, &y, 0));
    try testing.expectError(error.InvalidM, ger(al, check.io, f64, -1, 2, 1, &x, 1, &y, 1, &a, 2));
    try testing.expectError(error.InvalidN, ger(al, check.io, f64, 2, -1, 1, &x, 1, &y, 1, &a, 2));
    try testing.expectError(error.InvalidIncx, ger(al, check.io, f64, 2, 2, 1, &x, 0, &y, 1, &a, 2));
    try testing.expectError(error.InvalidIncy, ger(al, check.io, f64, 2, 2, 1, &x, 1, &y, 0, &a, 2));
    try testing.expectError(error.InvalidLda, ger(al, check.io, f64, 2, 2, 1, &x, 1, &y, 1, &a, 1));
    try testing.expectError(error.InvalidUplo, symv(al, check.io, f64, .v, 2, 1, &a, 2, &x, 1, 1, &y, 1));
    try testing.expectError(error.InvalidN, symv(al, check.io, f64, .u, -1, 1, &a, 2, &x, 1, 1, &y, 1));
    try testing.expectError(error.InvalidLda, symv(al, check.io, f64, .u, 2, 1, &a, 1, &x, 1, 1, &y, 1));
    try testing.expectError(error.InvalidIncx, symv(al, check.io, f64, .l, 2, 1, &a, 2, &x, 0, 1, &y, 1));
    try testing.expectError(error.InvalidIncy, symv(al, check.io, f64, .l, 2, 1, &a, 2, &x, 1, 1, &y, 0));
    var za = [_]Complex64{ .init(1, 0), .init(0, 0), .init(0, 0), .init(1, 0) };
    var zx = [_]Complex64{ .init(1, 0), .init(1, 0) };
    var zy = [_]Complex64{ .init(1, 0), .init(1, 0) };
    try testing.expectError(error.InvalidLda, hemv(al, check.io, Complex64, .m, 2, .init(1, 0), &za, 1, &zx, 1, .init(1, 0), &zy, 1));
    // empty problems leave y alone, even with a scaling beta
    try gemv(al, check.io, f64, .n, 0, 2, 1, &a, 1, &x, 1, 2, &y, 1);
    try gemv(al, check.io, f64, .t, 2, 0, 1, &a, 2, &x, 1, 2, &y, 1);
    try ger(al, check.io, f64, 0, 2, 1, &x, 1, &y, 1, &a, 1);
    try symv(al, check.io, f64, .u, 0, 1, &a, 1, &x, 1, 2, &y, 1);
    try testing.expectEqualSlices(f64, &.{ 5, 6 }, &y);
    try testing.expectEqualSlices(f64, &.{ 1, 2, 3, 4 }, &a);
}

/// `ger` or `gerc` against plain arithmetic; with `holes`, every other `y` is
/// zero, so the interface makes a run of each column between.
pub fn testGer(comptime T: type, conj: bool, m: usize, n: usize, incx: isize, incy: isize, alpha: T, holes: bool) !void {
    const lda = m + 2;
    const a = try testing.allocator.alloc(T, lda * n);
    defer testing.allocator.free(a);
    check.fill(T, a, m * 5 + n);
    const x = try randomVec(T, m, incx, 6);
    defer testing.allocator.free(x);
    const y = try randomVec(T, n, incy, 7);
    defer testing.allocator.free(y);
    if (holes) for (0..n) |j| if (j % 2 == 1) {
        y[vi(n, incy, j)] = zeroT(T);
    };
    const want = try testing.allocator.dupe(T, a);
    defer testing.allocator.free(want);
    for (0..n) |j| for (0..m) |i| {
        const yv = y[vi(n, incy, j)];
        const t = mulT(T, x[vi(m, incx, i)], if (conj) conjT(T, yv) else yv);
        want[i + j * lda] = addT(T, want[i + j * lda], mulT(T, alpha, t));
    };
    if (conj) {
        if (comptime !isComplex(T)) unreachable;
        try gerc(poison.allocator(), check.io, T, @intCast(m), @intCast(n), alpha, x.ptr, incx, y.ptr, incy, a.ptr, @intCast(lda));
    } else {
        try ger(poison.allocator(), check.io, T, @intCast(m), @intCast(n), alpha, x.ptr, incx, y.ptr, incy, a.ptr, @intCast(lda));
    }
    try expectSame(T, want, a, 1);
}

test "interface: ger and gerc match plain arithmetic, for every stride and an alpha of 0" {
    const saved = usePortable();
    defer blas.testRestore(saved);
    inline for (all) |T| {
        for (shapes) |s| for (stride_pairs) |inc| for ([_]T{ val(T, 1.25, -0.5), val(T, 0, 0), val(T, -1, 0) }) |alpha| {
            for ([_]bool{ false, true }) |holes| {
                testGer(T, false, s.m, s.n, inc[0], inc[1], alpha, holes) catch |e| {
                    check.say("ger {s} m={d} n={d} inc={d},{d} holes={}\n", .{ @typeName(T), s.m, s.n, inc[0], inc[1], holes });
                    return e;
                };
                if (comptime isComplex(T)) try testGer(T, true, s.m, s.n, inc[0], inc[1], alpha, holes);
            }
        };
    }
}

/// `S[i][j]` of the symmetric or Hermitian matrix `uplo` stores in `a`.
pub fn symAt(comptime T: type, kind: sym.Kind, upper: bool, a: []const T, lda: usize, i: usize, j: usize) T {
    const stored = if (upper) i <= j else i >= j;
    const e = if (stored) a[i + j * lda] else a[j + i * lda];
    if (i == j) return if (comptime !isComplex(T)) e else if (kind == .symmetric) e else val(T, e.re, 0);
    return switch (kind) {
        .symmetric => e,
        .hermitian => if (stored) e else conjT(T, e),
        .hermitian_rev => if (stored) conjT(T, e) else e,
    };
}

pub fn testMv(comptime T: type, herm: bool, uplo: Uplo, n: usize, incx: isize, incy: isize, alpha: T, beta: T) !void {
    const lda = n + 1;
    const a = try testing.allocator.alloc(T, lda * n);
    defer testing.allocator.free(a);
    check.fill(T, a, n * 3 + 1);
    const x = try randomVec(T, n, incx, 8);
    defer testing.allocator.free(x);
    const y = try randomVec(T, n, incy, 9);
    defer testing.allocator.free(y);
    const want = try testing.allocator.dupe(T, y);
    defer testing.allocator.free(want);
    const upper = uplo == .u or uplo == .v;
    const kind: sym.Kind = if (uplo == .v or uplo == .m) .hermitian_rev else if (herm) .hermitian else .symmetric;
    for (0..n) |i| {
        var acc = zeroT(T);
        for (0..n) |j| acc = addT(T, acc, mulT(T, symAt(T, kind, upper, a, lda, i, j), x[vi(n, incx, j)]));
        const k = vi(n, incy, i);
        const scaled = if (isZero(T, beta)) zeroT(T) else mulT(T, beta, y[k]);
        want[k] = if (isZero(T, alpha)) scaled else addT(T, scaled, mulT(T, alpha, acc));
    }
    if (kind == .symmetric) {
        try symv(poison.allocator(), check.io, T, uplo, @intCast(n), alpha, a.ptr, @intCast(lda), x.ptr, incx, beta, y.ptr, incy);
    } else if (comptime isComplex(T)) {
        try hemv(poison.allocator(), check.io, T, uplo, @intCast(n), alpha, a.ptr, @intCast(lda), x.ptr, incx, beta, y.ptr, incy);
    }
    try expectSame(T, want, y, n);
}

const mv_sizes = [_]usize{ 1, 2, 5, 8, 9, 17 };

test "interface: symv and hemv match plain arithmetic, for every triangle, stride and scalar case" {
    const saved = usePortable();
    defer blas.testRestore(saved);
    inline for (all) |T| {
        for ([_]bool{ false, true }) |herm| {
            if (herm and !comptime isComplex(T)) continue;
            const uplos = if (herm) std.enums.values(Uplo) else &[_]Uplo{ .u, .l };
            for (uplos) |uplo| for (mv_sizes) |n| for (stride_pairs) |inc| for (scalars(T)) |ab| {
                testMv(T, herm, uplo, n, inc[0], inc[1], ab[0], ab[1]) catch |e| {
                    check.say("{s} {s} {s} n={d} inc={d},{d}\n", .{ if (herm) "hemv" else "symv", @typeName(T), @tagName(uplo), n, inc[0], inc[1] });
                    return e;
                };
            };
        }
    }
}

test "interface: problems over the stack buffer allocate, and small gemv and ger never do" {
    const saved = usePortable();
    defer blas.testRestore(saved);
    inline for (all) |T| {
        try testGemv(T, .n, 300, 300, 1, 1, val(T, 1.5, 0.5), val(T, 0.5, 0));
        try testGemv(T, if (isComplex(T)) .c else .t, 300, 300, -2, 3, val(T, 1.5, 0.5), val(T, 1, 0));
        try testGer(T, false, 600, 3, 2, -1, val(T, 0.75, 0.25), true);
        try testMv(T, false, .l, 300, -1, 2, val(T, 1.5, 0.5), val(T, 0.5, 0));
        if (comptime isComplex(T)) {
            try testMv(T, true, .u, 300, 2, 1, val(T, 1.5, 0.5), val(T, 0.5, 0));
        }
        // a failing allocator is never asked for tiny problems
        var a = [_]T{ val(T, 1, 0), val(T, 2, 0), val(T, 3, 0), val(T, 4, 0) };
        var x = [_]T{ val(T, 1, 0), val(T, 1, 0) };
        var y = [_]T{ val(T, 1, 0), val(T, 1, 0) };
        const fa = testing.failing_allocator;
        try gemv(fa, check.io, T, .n, 2, 2, val(T, 1, 0), &a, 2, &x, 1, val(T, 1, 0), &y, 1);
        try ger(fa, check.io, T, 2, 2, val(T, 1, 0), &x, 2, &y, 1, &a, 2);
        // nor for symv on the portable gemv (plain loops), or a unit-stride ger
        try symv(fa, check.io, T, .u, 2, val(T, 1, 0), &a, 2, &x, 1, val(T, 1, 0), &y, 1);
        const big = 600;
        const bx = try testing.allocator.alloc(T, 2 * big);
        defer testing.allocator.free(bx);
        check.fill(T, bx, 3);
        try ger(fa, check.io, T, big, 1, val(T, 1, 0), bx.ptr, 1, &y, 2, bx.ptr, big);
        // and a big one that needs a buffer reports the failure
        try testing.expectError(error.OutOfMemory, ger(fa, check.io, T, big, 1, val(T, 1, 0), bx.ptr, 2, &y, 2, bx.ptr, big));
    }
}

// ---- half-precision gemv tests ----------------------------------------------

fn halfIn(comptime E: type, v: f32) E {
    return if (E == types.BFloat16) bf.toBf16(v) else @floatCast(v);
}

fn halfVal(comptime E: type, v: E) f64 {
    return if (E == types.BFloat16) bf.fromBf16(v) else @floatCast(v);
}

pub fn testGemvHalf(comptime kind: HalfPrecision, trans: Trans, m: usize, n: usize, incx: isize, incy: isize, alpha: f32, beta: f32) !void {
    const In = kind.In();
    const Out = kind.Out();
    const tr = trans.transposed();
    const lenx = if (tr) m else n;
    const leny = if (tr) n else m;
    const lda = m + 2;
    const a = try testing.allocator.alloc(In, lda * n);
    defer testing.allocator.free(a);
    const x = try testing.allocator.alloc(In, guard.span(lenx, @abs(incx)));
    defer testing.allocator.free(x);
    const y = try testing.allocator.alloc(Out, guard.span(leny, @abs(incy)));
    defer testing.allocator.free(y);
    var f: [1024]f32 = undefined;
    check.fill(f32, f[0..a.len], m * 7 + n);
    for (a, f[0..a.len]) |*d, v| d.* = halfIn(In, v);
    check.fill(f32, f[0..x.len], 3);
    for (x, f[0..x.len]) |*d, v| d.* = halfIn(In, v);
    check.fill(f32, f[0..y.len], 4);
    for (y, f[0..y.len]) |*d, v| d.* = halfIn(Out, v);
    const want = try testing.allocator.dupe(Out, y);
    defer testing.allocator.free(want);
    for (0..leny) |r| {
        var acc: f64 = 0;
        for (0..lenx) |q| {
            const e = if (tr) a[q + r * lda] else a[r + q * lda];
            acc += halfVal(In, e) * halfVal(In, x[vi(lenx, incx, q)]);
        }
        const k = vi(leny, incy, r);
        const af: f64 = halfVal(Out, halfIn(Out, alpha));
        const bt: f64 = halfVal(Out, halfIn(Out, beta));
        const prior: f64 = if (bt == 0) 0 else bt * halfVal(Out, y[k]);
        want[k] = halfIn(Out, @floatCast(if (af == 0) prior else af * acc + prior));
    }
    try gemvHalf(check.io, kind, trans, @intCast(m), @intCast(n), halfIn(Out, alpha), a.ptr, @intCast(lda), x.ptr, incx, halfIn(Out, beta), y.ptr, incy);
    for (want, y, 0..) |w, g, i| {
        const wv = halfVal(Out, w);
        const tol: f64 = if (Out == types.BFloat16) 0.012 * @abs(wv) + 1e-3 else 1e-5 * @as(f64, @floatFromInt(lenx + 2)) + 1e-6 * @abs(wv);
        if (!(@abs(wv - halfVal(Out, g)) <= tol)) {
            check.say("{s} {s} m={d} n={d} inc={d},{d} alpha={d} beta={d} y[{d}]: got {d}, want {d}\n", .{ @tagName(kind), @tagName(trans), m, n, incx, incy, alpha, beta, i, halfVal(Out, g), wv });
            return error.TestExpectedApproxEqAbs;
        }
    }
}

const half_scalars = [_][2]f32{ .{ 1.25, 0.5 }, .{ 1, 1 }, .{ 0, 0.5 }, .{ 1.5, 0 }, .{ 0, 1 }, .{ -1, 2 }, .{ 0, 0 } };

test "interface: the half-precision gemvs match plain arithmetic, for every trans, stride and scalar case" {
    const saved = usePortable();
    defer blas.testRestore(saved);
    inline for (std.enums.values(HalfPrecision)) |kind| {
        for ([_]Trans{ .n, .t, .r, .c }) |trans| for (shapes) |s| for (stride_pairs) |inc| for (half_scalars) |ab| {
            testGemvHalf(kind, trans, s.m, s.n, inc[0], inc[1], ab[0], ab[1]) catch |e| return e;
        };
    }
}

test "interface: half-precision gemv beta 0 overwrites NaN, alpha 0 never reads x or a" {
    const saved = usePortable();
    defer blas.testRestore(saved);
    inline for (std.enums.values(HalfPrecision)) |kind| {
        const In = kind.In();
        const Out = kind.Out();
        const nan_in = halfIn(In, std.math.nan(f32));
        var a = [_]In{ nan_in, nan_in, nan_in, nan_in };
        var x = [_]In{ nan_in, nan_in };
        var y = [_]Out{ halfIn(Out, std.math.nan(f32)), halfIn(Out, std.math.nan(f32)) };
        try gemvHalf(check.io, kind, .n, 2, 2, halfIn(Out, 0), &a, 2, &x, 1, halfIn(Out, 0), &y, 1);
        try testing.expectEqual(@as(f64, 0), halfVal(Out, y[0]));
        try testing.expectEqual(@as(f64, 0), halfVal(Out, y[1]));
        y = .{ halfIn(Out, 1), halfIn(Out, 2) };
        try gemvHalf(check.io, kind, .t, 2, 2, halfIn(Out, 0), &a, 2, &x, -1, halfIn(Out, 2), &y, 1);
        try testing.expectEqual(@as(f64, 2), halfVal(Out, y[0]));
        try testing.expectEqual(@as(f64, 4), halfVal(Out, y[1]));
        // beta 1 with alpha 0 touches nothing, even through a stride
        y = .{ halfIn(Out, 1), halfIn(Out, 2) };
        try gemvHalf(check.io, kind, .n, 2, 2, halfIn(Out, 0), &a, 2, &x, 1, halfIn(Out, 1), &y, -1);
        try testing.expectEqual(@as(f64, 1), halfVal(Out, y[0]));
    }
}

test "interface: half-precision gemv reports the first bad argument and ignores empty problems" {
    const saved = usePortable();
    defer blas.testRestore(saved);
    inline for (std.enums.values(HalfPrecision)) |kind| {
        const In = kind.In();
        const Out = kind.Out();
        var a = [_]In{ halfIn(In, 1), halfIn(In, 2), halfIn(In, 3), halfIn(In, 4) };
        var x = [_]In{ halfIn(In, 1), halfIn(In, 1) };
        var y = [_]Out{ halfIn(Out, 5), halfIn(Out, 6) };
        const one = halfIn(Out, 1);
        try testing.expectError(error.InvalidTrans, gemvHalf(check.io, kind, .o, 2, 2, one, &a, 2, &x, 1, one, &y, 1));
        try testing.expectError(error.InvalidTrans, gemvHalf(check.io, kind, .d, -1, -1, one, &a, 0, &x, 0, one, &y, 0));
        try testing.expectError(error.InvalidM, gemvHalf(check.io, kind, .n, -1, -1, one, &a, 0, &x, 0, one, &y, 0));
        try testing.expectError(error.InvalidN, gemvHalf(check.io, kind, .n, 2, -1, one, &a, 0, &x, 0, one, &y, 0));
        try testing.expectError(error.InvalidLda, gemvHalf(check.io, kind, .n, 2, 2, one, &a, 1, &x, 0, one, &y, 0));
        try testing.expectError(error.InvalidIncx, gemvHalf(check.io, kind, .t, 2, 2, one, &a, 2, &x, 0, one, &y, 0));
        try testing.expectError(error.InvalidIncy, gemvHalf(check.io, kind, .t, 2, 2, one, &a, 2, &x, 1, one, &y, 0));
        // an empty problem does not even apply beta
        try gemvHalf(check.io, kind, .n, 0, 2, one, &a, 1, &x, 1, halfIn(Out, 0), &y, 1);
        try gemvHalf(check.io, kind, .n, 2, 0, one, &a, 2, &x, 1, halfIn(Out, 0), &y, 1);
        try testing.expectEqual(@as(f64, 5), halfVal(Out, y[0]));
    }
}
