//! Tests of `../matrix.zig`.

const std = @import("std");
const rstd_matrix = @import("rstd_matrix");
const vector = @import("../vector.zig");
const ops = @import("../../kernels/lapack/ops.zig");
const mroot = @import("../matrix.zig");
const mops = @import("../matrix/ops.zig");
const msol = @import("../matrix/solve.zig");

const Vec = vector.Vec;
const isComplex = vector.isComplex;
const Real = vector.Real;
const M2f32 = mroot.M2f32;
const M2f64 = mroot.M2f64;
const M3f32 = mroot.M3f32;
const M3f64 = mroot.M3f64;
const M4f32 = mroot.M4f32;
const M4f64 = mroot.M4f64;
const Mat = mroot.Mat;
const StaticMat = mroot.StaticMat;
const View = mroot.View;
const invInto = mroot.invInto;
const mulInto = mroot.mulInto;
const add = mops.add;
const col = mops.col;
const ctranspose = mops.ctranspose;
const fromCols = mops.fromCols;
const fromRows = mops.fromRows;
const identity = mops.identity;
const matVec = mops.matVec;
const mul = mops.mul;
const mulCols = mops.mulCols;
const outer = mops.outer;
const row = mops.row;
const scale = mops.scale;
const sub = mops.sub;
const trace = mops.trace;
const transpose = mops.transpose;
const transposeVecs = mops.transposeVecs;
const vecMat = mops.vecMat;
const det = msol.det;
const inv = msol.inv;
const solve = msol.solve;

// ---- tests ------------------------------------------------------------------

const zigstd = @import("../../../zigstd.zig");
const V2f32 = vector.V2f32;
const V3f32 = vector.V3f32;
const V4f64 = vector.V4f64;

const bitwise = @import("../../kernels/blas/_test/bitwise.zig");
const expectSameBits = bitwise.expectSameBits;
const math = @import("../../__root__.zig");

test "transposeVecs is pffft's VTRANSPOSE4 (_MM_TRANSPOSE4_PS)" {
    const V = vector.V4f32;
    const t = transposeVecs(4, f32, .{ V{ 0, 1, 2, 3 }, V{ 4, 5, 6, 7 }, V{ 8, 9, 10, 11 }, V{ 12, 13, 14, 15 } });
    try zigstd.testExpectEqual(V{ 0, 4, 8, 12 }, t[0]);
    try zigstd.testExpectEqual(V{ 3, 7, 11, 15 }, t[3]);
    const d = transposeVecs(2, f64, .{ .{ 1, 2 }, .{ 3, 4 } });
    try zigstd.testExpectEqual(vector.V2f64{ 1, 3 }, d[0]);
    try zigstd.testExpectEqual(vector.V2f64{ 2, 4 }, d[1]);
}

test "StaticMat construction, access and transposes" {
    const M = StaticMat(2, 3, f32);
    const m = fromRows(M, .{ .{ 1, 3, 5 }, .{ 2, 4, 6 } });
    try zigstd.testExpectEqual(fromCols(M, .{ .{ 1, 2 }, .{ 3, 4 }, .{ 5, 6 } }), m);
    try zigstd.testExpectEqual(V3f32{ 2, 4, 6 }, row(m, 1));
    try zigstd.testExpectEqual(V2f32{ 5, 6 }, col(m, 2));
    try zigstd.testExpectEqual(fromRows(StaticMat(3, 2, f32), .{ .{ 1, 2 }, .{ 3, 4 }, .{ 5, 6 } }), transpose(m));
    try zigstd.testExpectEqual(m, transpose(transpose(m)));
    try zigstd.testExpectEqual([3][2]f32{ .{ 1, 0 }, .{ 0, 1 }, .{ 0, 0 } }, identity(M).data);
    try zigstd.testExpectEqual(@as(f32, 4), trace(identity(M4f32)));

    const C = std.math.Complex(f64);
    const z: StaticMat(1, 2, C) = .{ .data = .{ .{.init(1, 2)}, .{.init(3, -4)} } };
    const zh = ctranspose(z);
    try zigstd.testExpectEqual(C.init(1, -2), zh.data[0][0]);
    try zigstd.testExpectEqual(C.init(3, 4), zh.data[0][1]);
    try zigstd.testExpectEqual(C.init(3, -4), transpose(z).data[0][1]);
}

test "views: init's checks, sub, col; StaticMat and Mat views share storage" {
    var buf = [_]f64{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    try std.testing.expectError(error.DimensionMismatch, View(f64).init(&buf, 3, 2, 2));
    try std.testing.expectError(error.SliceTooShort, View(f64).init(&buf, 3, 4, 3));
    const v = try View(f64).init(&buf, 3, 3, 3);
    try std.testing.expectEqual(@as(f64, 8), v.at(1, 2));
    const s = try v.sub(1, 1, 2, 2);
    try std.testing.expectEqual(@as(f64, 5), s.at(0, 0));
    try std.testing.expectEqual(@as(f64, 9), s.at(1, 1));
    try std.testing.expectEqualSlices(f64, &.{ 5, 6 }, s.col(0));
    try std.testing.expectError(error.DimensionMismatch, v.sub(2, 0, 2, 1));

    var sm = fromCols(M2f64, .{ .{ 1, 2 }, .{ 3, 4 } });
    sm.view().data[3] = 40;
    try std.testing.expectEqual(@as(f64, 40), sm.data[1][1]);
    try std.testing.expectEqual(@as(f64, 3), sm.constView().at(0, 1));
    const m = try Mat(f64).zeros(std.testing.allocator, 2, 3);
    defer m.deinit(std.testing.allocator);
    m.view().data[5] = 7;
    try std.testing.expectEqual(@as(f64, 7), m.at(1, 2));
}

test "mul, scale, add, sub and outer on exact values" {
    // [1 3 5; 2 4 6] and [7 10; 8 11; 9 12], as the fork's behavior test
    const a = fromCols(StaticMat(2, 3, f32), .{ .{ 1, 2 }, .{ 3, 4 }, .{ 5, 6 } });
    const b = fromCols(StaticMat(3, 2, f32), .{ .{ 7, 8, 9 }, .{ 10, 11, 12 } });
    try zigstd.testExpectEqual([2][2]f32{ .{ 76, 100 }, .{ 103, 136 } }, mul(a, b).data);
    try zigstd.testExpectEqual(V2f32{ 9, 12 }, mul(a, V3f32{ 1, 1, 1 }));
    try zigstd.testExpectEqual([2]f32{ 9, 12 }, mul(a, [3]f32{ 1, 1, 1 }));
    try zigstd.testExpectEqual(V3f32{ 5, 11, 17 }, mul(V2f32{ 1, 2 }, a));
    try zigstd.testExpectEqual(add(a, a), scale(a, 2));
    try zigstd.testExpectEqual(StaticMat(2, 3, f32){ .data = @splat(@splat(0)) }, sub(a, a));
    try zigstd.testExpectEqual([3][2]f32{ .{ 1, 2 }, .{ 2, 4 }, .{ 3, 6 } }, outer(V2f32{ 1, 2 }, [3]f32{ 1, 2, 3 }).data);
    var out: M2f32 = undefined;
    mulInto(&out, a, b);
    try zigstd.testExpectEqual(mul(a, b), out);

    const C = std.math.Complex(f64);
    const za: StaticMat(2, 2, C) = .{ .data = .{ .{ .init(1, 1), .init(0, 2) }, .{ .init(3, 0), .init(1, -1) } } };
    const zv = [2]C{ .init(2, 0), .init(0, 1) };
    // [1+i 3; 2i 1-i] [2; i] = [2+5i; 1+5i]
    const zr = mul(za, zv);
    try zigstd.testExpectEqual(C.init(2, 5), zr[0]);
    try zigstd.testExpectEqual(C.init(1, 5), zr[1]);
    // [2 i] [1+i 3; 2i 1-i] = [0+2i, 7+i]
    const zl = mul(zv, za);
    try zigstd.testExpectEqual(C.init(0, 2), zl[0]);
    try zigstd.testExpectEqual(C.init(7, 1), zl[1]);
    try zigstd.testExpectEqual(mul(za, za).data[1], mul(za, za.data[1]));
    // a bᴴ conjugates b
    try zigstd.testExpectEqual(C.init(0, -2), outer([1]C{.init(2, 0)}, [1]C{.init(0, 1)}).data[0][0]);
}

fn randomMat(comptime M: type, random: std.Random, s: M.Element) M {
    var m: M = undefined;
    for (&m.data) |*c| for (c) |*x| {
        x.* = (random.float(M.Element) * 2 - 1) * s;
    };
    return m;
}

fn expectSameMatBits(want: anytype, got: @TypeOf(want)) !void {
    for (want.data, got.data) |wc, gc| try expectSameOrX87(@as(Vec(@TypeOf(want).rows, @TypeOf(want).Element), wc), gc);
}

/// `expectSameBits` of two `@Vector`s, but for `f64` on x87 within
/// `bitwise.sameBitsOrX87`'s few ulps: there each code path rounds its 80-bit
/// intermediates where it stores them (issue #264).
fn expectSameOrX87(want: anytype, got: @TypeOf(want)) !void {
    const V = @typeInfo(@TypeOf(want)).vector;
    if (comptime !(math.x87_f64 and V.child == f64)) return expectSameBits(want, got);
    const wa: [V.len]f64 = want;
    const ga: [V.len]f64 = got;
    try zigstd.testExpectTrue(bitwise.sameBitsOrX87(f64, &wa, &ga));
}

/// `mul`'s `@Matrix` forms against its column code, for `R` by `K` times `K`
/// by `C`.
fn expectNativeBits(comptime T: type, comptime r: comptime_int, comptime k: comptime_int, comptime c: comptime_int, random: std.Random) !void {
    const scales = [_]f64{ 1, 1e-3, 1e3, 1e30, 1e-30 };
    for (0..50) |trial| {
        const s: T = @floatCast(scales[trial % scales.len]);
        const a = randomMat(StaticMat(r, k, T), random, s);
        const b = randomMat(StaticMat(k, c, T), random, s);
        const v: Vec(k, T) = b.data[0];
        const w: Vec(r, T) = a.data[0];
        try expectSameMatBits(mulCols(a, b), mul(a, b));
        try expectSameOrX87(@as(Vec(r, T), matVec(a, v)), mul(a, v));
        try expectSameOrX87(@as(Vec(k, T), vecMat(w, a)), mul(w, a));
        try expectSameMatBits(transpose(transpose(a)), a);
    }
}

test "mul through @Matrix has its column code's bits" {
    if (!rstd_matrix.has_matrix) return error.SkipZigTest;
    var prng: std.Random.DefaultPrng = .init(0x3a7);
    const shapes = .{ .{ 1, 1, 1 }, .{ 2, 2, 2 }, .{ 3, 3, 3 }, .{ 4, 4, 4 }, .{ 2, 3, 4 }, .{ 5, 7, 3 }, .{ 6, 6, 1 }, .{ 8, 8, 8 } };
    inline for (.{ f32, f64 }) |T| inline for (shapes) |s| try expectNativeBits(T, s[0], s[1], s[2], prng.random());
}

test "det, inv and solve on known matrices" {
    const m2 = fromRows(M2f64, .{ .{ 4, 7 }, .{ 2, 6 } });
    try zigstd.testExpectEqual(@as(f64, 10), det(m2));
    try zigstd.testExpectEqual(fromRows(M2f64, .{ .{ 0.6, -0.7 }, .{ -0.2, 0.4 } }), try inv(m2));
    const m3 = fromRows(M3f64, .{ .{ 1, 2, 3 }, .{ 0, 1, 4 }, .{ 5, 6, 0 } });
    try zigstd.testExpectEqual(@as(f64, 1), det(m3));
    try zigstd.testExpectEqual(fromRows(M3f64, .{ .{ -24, 18, 5 }, .{ 20, -15, -4 }, .{ -5, 4, 1 } }), try inv(m3)); // a permutation needs a row swap; its determinant is the swaps' sign
    const p4 = fromCols(M4f64, .{ .{ 0, 1, 0, 0 }, .{ 1, 0, 0, 0 }, .{ 0, 0, 0, 1 }, .{ 0, 0, 1, 0 } });
    try zigstd.testExpectEqual(@as(f64, 1), det(p4));
    try zigstd.testExpectEqual(transpose(p4), try inv(p4));
    try zigstd.testExpectEqual(V4f64{ 2, 1, 4, 3 }, try solve(p4, V4f64{ 1, 2, 3, 4 }));
    try zigstd.testExpectEqual(@as(f64, -2), det(fromRows(StaticMat(5, 5, f64), .{
        .{ 0, 2, 0, 0, 0 }, .{ 1, 0, 0, 0, 0 }, .{ 0, 0, 1, 0, 0 }, .{ 0, 0, 0, 1, 0 }, .{ 0, 0, 0, 0, 1 },
    })));
    try zigstd.testExpectEqual(@as(f32, -3), det(StaticMat(1, 1, f32){ .data = .{.{-3}} }));
    try zigstd.testExpectEqual(@as(f32, -0.5), (try inv(StaticMat(1, 1, f32){ .data = .{.{-2}} })).data[0][0]);

    const C = std.math.Complex(f64);
    const z: StaticMat(2, 2, C) = .{ .data = .{ .{ .init(0, 1), .init(0, 0) }, .{ .init(0, 0), .init(2, 0) } } };
    try zigstd.testExpectEqual(C.init(0, 2), det(z));
    try zigstd.testExpectEqual(C.init(0, -1), (try inv(z)).data[0][0]);
}

test "det, inv and solve on singular matrices" {
    try zigstd.testExpectEqual(@as(f64, 0), det(fromRows(M2f64, .{ .{ 1, 2 }, .{ 2, 4 } })));
    try zigstd.testExpectError(error.Singular, inv(fromRows(M2f64, .{ .{ 1, 2 }, .{ 2, 4 } })));
    try zigstd.testExpectError(error.Singular, inv(fromRows(M3f32, .{ .{ 1, 2, 3 }, .{ 4, 5, 6 }, .{ 7, 8, 9 } })));
    try zigstd.testExpectError(error.Singular, solve(fromRows(M3f32, .{ .{ 1, 2, 3 }, .{ 0, 0, 0 }, .{ 7, 8, 9 } }), V3f32{ 1, 2, 3 }));
    const s5 = StaticMat(5, 5, f64){ .data = @splat(@splat(1)) };
    try zigstd.testExpectEqual(@as(f64, 0), det(s5));
    try zigstd.testExpectError(error.Singular, inv(s5));
    var out = identity(M4f64);
    try zigstd.testExpectError(error.Singular, invInto(&out, M4f64{ .data = @splat(@splat(0)) }));
    try zigstd.testExpectEqual(identity(M4f64), out);
    const C = std.math.Complex(f32);
    try zigstd.testExpectError(error.Singular, inv(StaticMat(2, 2, C){ .data = .{ .{ .init(1, 1), .init(2, 2) }, .{ .init(1, 1), .init(2, 2) } } }));
}

fn isClose(x: anytype, y: @TypeOf(x), tol: Real(@TypeOf(x))) bool {
    return ops.abs1(@TypeOf(x), ops.sub(@TypeOf(x), x, y)) <= tol;
}

/// `inv(m) m` and `m inv(m)` are `I` and `m solve(m, b)` is `b` within `tol`,
/// for random diagonally dominant `n` by `n` matrices.
fn expectInverts(comptime T: type, comptime n: comptime_int, random: std.Random) !void {
    const M = StaticMat(n, n, T);
    const R = Real(T);
    const tol: R = 1e-12;
    for (0..20) |_| {
        var m: M = undefined;
        var b: [n]T = undefined;
        for (0..n) |j| {
            for (0..n) |i| m.data[j][i] = randomElem(T, random);
            m.data[j][j] = ops.add(T, m.data[j][j], if (comptime isComplex(T)) T.init(2 * n, 0) else 2 * n);
            b[j] = randomElem(T, random);
        }
        const mi = try inv(m);
        const eye = identity(M);
        for (mul(mi, m).data, mul(m, mi).data, eye.data) |p, q, e| for (p, q, e) |x, y, z| {
            try zigstd.testExpectTrue(isClose(x, z, tol) and isClose(y, z, tol));
        };
        const x = try solve(m, b);
        for (mul(m, x), b) |got, want| try zigstd.testExpectTrue(isClose(got, want, tol));
        const xs = try solve(m, fromCols(StaticMat(n, 1, T), .{b}));
        try zigstd.testExpectTrue(bitwise.sameBitsOrX87(T, &x, &xs.data[0]));
        // the determinant is the product of the eigenvalues: of `mᵀ` the same
        const d = det(m);
        const dt = det(transpose(m));
        try zigstd.testExpectTrue(isClose(d, dt, tol * ops.abs1(T, d)));
    }
}

fn randomElem(comptime T: type, random: std.Random) T {
    if (comptime isComplex(T)) return .init(random.float(Real(T)) * 2 - 1, random.float(Real(T)) * 2 - 1);
    return random.float(T) * 2 - 1;
}

test "inv and solve invert, real and complex, n 1 to 8" {
    var prng: std.Random.DefaultPrng = .init(0x1a);
    inline for (.{ f64, std.math.Complex(f64) }) |T| inline for (1..9) |n| try expectInverts(T, n, prng.random());
}
