//! `trsyl`: ReLAPACK's recursion (`relapack/src/dtrsyl.c` and the rest,
//! OpenBLAS v0.3.34-317, MIT) over `sylvester.zig`'s `trsy2`.

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const driver = @import("../blas/blas3/driver.zig");
const ops = @import("ops.zig");
const relapack = @import("relapack.zig");
const sylvester = @import("sylvester.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const isComplex = types.isComplex;
const Trans = types.Trans;
const scalar = ops.scalar;
const Solved = sylvester.Solved;
const scaleBlock = sylvester.scaleBlock;
const trsy2 = sylvester.trsy2;

const Error = driver.Error;

/// What `trsyl` returns: bad arguments in the reference's order;
/// `Perturbed` (the reference's `INFO = 1`: `A` and `-isgn B` have close or
/// common eigenvalues, so perturbed values were used), with `C` and `scale`
/// written; `error.OutOfMemory` or `error.Canceled`.
pub const TrsylError = error{ InvalidTrana, InvalidTranb, InvalidIsgn, InvalidM, InvalidN, InvalidLda, InvalidLdb, InvalidLdc, Perturbed } || Allocator.Error || Io.Cancelable;

/// `?TRSYL` by ReLAPACK's recursion (`relapack/src/dtrsyl.c`, `ztrsyl.c`
/// and the rest; OpenBLAS v0.3.34-317): `op(A) X + isgn X op(B) = scale C`
/// for `X`, into `C` (`m x n`), `A` (`m x m`) and `B` (`n x n`) upper
/// triangular, real ones in Schur form (2 x 2 diagonal blocks); `op` is `n`
/// or the (conjugate) transpose, `t` or `c` for real `T` and `c` only for
/// complex, as the reference. `scale <= 1` keeps `X` from overflowing.
///
/// The larger side splits in two (never inside a 2 x 2 block), the halves
/// joined by `gemm`; the base case is `trsy2`, the reference's
/// `xTRSYL`, which reads the whole of `A` and `B` (`xLANGE`).
///
/// The reference wins where ReLAPACK, solving the top half first (`op(A)`
/// transposed), rescales the bottom half by the second solve's `scale`; the
/// top half, solved first, takes it here.
pub fn trsyl(allocator: Allocator, io: Io, comptime T: type, trana: Trans, tranb: Trans, isgn: isize, m: isize, n: isize, a: [*]const T, lda: isize, b: [*]const T, ldb: isize, c: [*]T, ldc: isize, scale: *Float(T)) TrsylError!void {
    comptime ops.requireFloat(T, "trsyl");
    const ta = sylTrans(T, trana) orelse return error.InvalidTrana;
    const tb = sylTrans(T, tranb) orelse return error.InvalidTranb;
    if (isgn != 1 and isgn != -1) return error.InvalidIsgn;
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, m)) return error.InvalidLda;
    if (ldb < @max(1, n)) return error.InvalidLdb;
    if (ldc < @max(1, m)) return error.InvalidLdc;
    scale.* = 1;
    if (m == 0 or n == 0) return;
    if (!util.fitsMatrix(m, m, lda, @sizeOf(T))) return error.InvalidLda;
    if (!util.fitsMatrix(n, n, ldb, @sizeOf(T))) return error.InvalidLdb;
    if (!util.fitsMatrix(m, n, ldc, @sizeOf(T))) return error.InvalidLdc;
    const r = try trsylRec(allocator, io, T, ta, tb, isgn, @intCast(m), @intCast(n), a, @intCast(lda), b, @intCast(ldb), c, @intCast(ldc));
    scale.* = r.scale;
    if (r.perturbed) return error.Perturbed;
}

/// Whether `t` transposes `A` or `B` for `trsyl`, or null where the
/// reference rejects it (complex `t`).
fn sylTrans(comptime T: type, t: Trans) ?bool {
    return switch (t) {
        .n => false,
        .c => true,
        .t => if (comptime isComplex(T)) null else true,
        .r, .o, .u, .s, .d => null,
    };
}

/// `RELAPACK_?trsyl_rec`.
fn trsylRec(al: Allocator, io: Io, comptime T: type, ta: bool, tb: bool, isgn: isize, m: usize, n: usize, a: [*]const T, lda: usize, b: [*]const T, ldb: usize, c: [*]T, ldc: usize) Error!Solved(Float(T)) {
    if (m <= relapack.crossover and n <= relapack.crossover) return trsy2(T, ta, tb, isgn, m, n, a, lda, b, ldb, c, ldc);
    const minus = scalar(T, -1);
    const msgn = scalar(T, @floatFromInt(-isgn));
    const opa: Trans = if (ta) relapack.ctrans(T) else .n;
    const opb: Trans = if (tb) relapack.ctrans(T) else .n;
    var r1: Solved(Float(T)) = undefined;
    var r2: Solved(Float(T)) = undefined;
    if (m > n) {
        var m1 = relapack.split(T, m);
        // never inside a 2 x 2 block
        if (comptime !isComplex(T)) if (a[m1 + lda * (m1 - 1)] != 0) {
            m1 += 1;
        };
        const m2 = m - m1;
        const a_tr = a + lda * m1;
        const a_br = a_tr + m1;
        const c_b = c + m1;
        if (!ta) {
            // C_T = scale1 C_T - A_TR C_B
            r1 = try trsylRec(al, io, T, ta, tb, isgn, m2, n, a_br, lda, b, ldb, c_b, ldc);
            try relapack.gemm(al, io, T, .n, .n, m1, n, m2, minus, a_tr, lda, c_b, ldc, scalar(T, r1.scale), c, ldc);
            r2 = try trsylRec(al, io, T, ta, tb, isgn, m1, n, a, lda, b, ldb, c, ldc);
            if (r2.scale != 1) scaleBlock(T, m2, n, r2.scale, c_b, ldc);
        } else {
            // C_B = scale1 C_B - A_TR' C_T
            r1 = try trsylRec(al, io, T, ta, tb, isgn, m1, n, a, lda, b, ldb, c, ldc);
            try relapack.gemm(al, io, T, opa, .n, m2, n, m1, minus, a_tr, lda, c, ldc, scalar(T, r1.scale), c_b, ldc);
            r2 = try trsylRec(al, io, T, ta, tb, isgn, m2, n, a_br, lda, b, ldb, c_b, ldc);
            if (r2.scale != 1) scaleBlock(T, m1, n, r2.scale, c, ldc);
        }
    } else {
        var n1 = relapack.split(T, n);
        if (comptime !isComplex(T)) if (b[n1 + ldb * (n1 - 1)] != 0) {
            n1 += 1;
        };
        const n2 = n - n1;
        const b_tr = b + ldb * n1;
        const b_br = b_tr + n1;
        const c_r = c + ldc * n1;
        if (!tb) {
            // C_R = scale1 C_R -/+ C_L B_TR
            r1 = try trsylRec(al, io, T, ta, tb, isgn, m, n1, a, lda, b, ldb, c, ldc);
            try relapack.gemm(al, io, T, .n, .n, m, n2, n1, msgn, c, ldc, b_tr, ldb, scalar(T, r1.scale), c_r, ldc);
            r2 = try trsylRec(al, io, T, ta, tb, isgn, m, n2, a, lda, b_br, ldb, c_r, ldc);
            if (r2.scale != 1) scaleBlock(T, m, n1, r2.scale, c, ldc);
        } else {
            // C_L = scale1 C_L -/+ C_R B_TR'
            r1 = try trsylRec(al, io, T, ta, tb, isgn, m, n2, a, lda, b_br, ldb, c_r, ldc);
            try relapack.gemm(al, io, T, .n, opb, m, n1, n2, msgn, c_r, ldc, b_tr, ldb, scalar(T, r1.scale), c, ldc);
            r2 = try trsylRec(al, io, T, ta, tb, isgn, m, n1, a, lda, b, ldb, c, ldc);
            if (r2.scale != 1) scaleBlock(T, m, n2, r2.scale, c_r, ldc);
        }
    }
    return .{ .scale = r1.scale * r2.scale, .perturbed = r1.perturbed or r2.perturbed };
}

// ---- tests ---------------------------------------------------------------------

const testing = std.testing;
const check = @import("../blas/_test/check.zig");
const reference = @import("../blas/_test/reference.zig");
const lcheck = @import("_test/check.zig");
const fixtures = @import("_test/fixtures.zig");

const schur = fixtures.schur;

fn trsylCase(comptime T: type, al: Allocator, ta: bool, tb: bool, isgn: isize, m: usize, n: usize, seed: u64) !void {
    const a = try schur(T, al, m, 2, 3, seed);
    defer al.free(a);
    // `-isgn B`'s spectrum away from `A`'s
    const b = try schur(T, al, n, if (isgn == 1) 2 else -4, 4, seed + 1);
    defer al.free(b);
    const c0 = try al.alloc(T, m * n);
    defer al.free(c0);
    check.fill(T, c0, seed + 2);
    const x = try al.dupe(T, c0);
    defer al.free(x);
    const op = struct {
        fn of(comptime U: type, t: bool) Trans {
            return if (!t) .n else if (comptime isComplex(U)) .c else .t;
        }
    };
    var scale: Float(T) = 0;
    try trsyl(al, check.io, T, op.of(T, ta), op.of(T, tb), isgn, @intCast(m), @intCast(n), a.ptr, @intCast(m), b.ptr, @intCast(n), x.ptr, @intCast(m), &scale);
    try testing.expect(scale > 0 and scale <= 1);
    const eps: f128 = std.math.floatEps(Float(T));
    const ulps: f128 = @floatFromInt(32 * (m + n + 2));
    const sg: f128 = @floatFromInt(isgn);
    for (0..n) |j| for (0..m) |i| {
        var r = reference.scaleW(-@as(f128, scale), reference.wide(T, c0[i + j * m]));
        var mag: f128 = scale * reference.abs1(reference.wide(T, c0[i + j * m]));
        for (0..m) |k| {
            const av = reference.wide(T, if (ta) a[k + i * m] else a[i + k * m]);
            const aw = if (ta) reference.conjW(av) else av;
            const xv = reference.wide(T, x[k + j * m]);
            r = reference.addW(r, reference.mulW(aw, xv));
            mag += reference.abs1(aw) * reference.abs1(xv);
        }
        for (0..n) |k| {
            const bv = reference.wide(T, if (tb) b[j + k * n] else b[k + j * n]);
            const bw = if (tb) reference.conjW(bv) else bv;
            const xv = reference.wide(T, x[i + k * m]);
            r = reference.addW(r, reference.scaleW(sg, reference.mulW(xv, bw)));
            mag += reference.abs1(bw) * reference.abs1(xv);
        }
        const d = reference.abs1(r);
        if (!(d <= ulps * eps * mag + std.math.floatMin(Float(T)))) {
            check.say("trsyl {s} {d}x{d} {any}{any} isgn {d}: residual ({d},{d}) {e}, allowed {e}\n", .{ @typeName(T), m, n, ta, tb, isgn, i, j, @as(f64, @floatCast(d)), @as(f64, @floatCast(ulps * eps * mag)) });
            return error.TestUnexpectedResult;
        }
    };
}

const trsyl_shapes = [_][2]usize{ .{ 1, 1 }, .{ 5, 3 }, .{ 24, 24 }, .{ 30, 10 }, .{ 10, 40 }, .{ 60, 50 }, .{ 100, 70 } };

test "relapack trsyl: op(A) X + isgn X op(B) = scale C, every op and sign, quasi-triangular real, every type" {
    const al = testing.allocator;
    inline for (lcheck.test_types) |T| for (trsyl_shapes) |s| for ([_]bool{ false, true }) |ta| for ([_]bool{ false, true }) |tb| for ([_]isize{ 1, -1 }) |sg| {
        if (check.rig() and s[0] * s[1] > 1000) continue;
        trsylCase(T, al, ta, tb, sg, s[0], s[1], s[0] * 13 + s[1]) catch |e| {
            check.say("trsyl {s} {d}x{d}\n", .{ @typeName(T), s[0], s[1] });
            return e;
        };
    };
}

test "relapack trsyl: common eigenvalues are perturbed, C and scale written; the first bad argument" {
    const al = testing.allocator;
    // A X - X A = C, A = 2 I: singular
    var a = [_]f64{ 2, 0, 0, 2 };
    var c = [_]f64{ 1, 2, 3, 4 };
    var scale: f64 = 0;
    try testing.expectError(error.Perturbed, trsyl(al, check.io, f64, .n, .n, -1, 2, 2, &a, 2, &a, 2, &c, 2, &scale));
    try testing.expect(scale > 0 and scale <= 1 and std.math.isFinite(c[0]));
    const z: [4]types.Complex64 = @splat(lcheck.poisoned(types.Complex64));
    var zc = z;
    try testing.expectError(error.InvalidTrana, trsyl(al, check.io, types.Complex64, .t, .n, 1, 2, 2, &z, 2, &z, 2, &zc, 2, &scale));
    try testing.expectError(error.InvalidTranb, trsyl(al, check.io, f64, .n, .r, 1, 2, 2, &a, 2, &a, 2, &c, 2, &scale));
    try testing.expectError(error.InvalidIsgn, trsyl(al, check.io, f64, .n, .n, 0, 2, 2, &a, 2, &a, 2, &c, 2, &scale));
    try testing.expectError(error.InvalidM, trsyl(al, check.io, f64, .n, .n, 1, -1, 2, &a, 2, &a, 2, &c, 2, &scale));
    try testing.expectError(error.InvalidN, trsyl(al, check.io, f64, .n, .n, 1, 2, -1, &a, 2, &a, 2, &c, 2, &scale));
    try testing.expectError(error.InvalidLda, trsyl(al, check.io, f64, .n, .n, 1, 2, 2, &a, 1, &a, 2, &c, 2, &scale));
    try testing.expectError(error.InvalidLdb, trsyl(al, check.io, f64, .n, .n, 1, 2, 2, &a, 2, &a, 1, &c, 2, &scale));
    try testing.expectError(error.InvalidLdc, trsyl(al, check.io, f64, .n, .n, 1, 2, 2, &a, 2, &a, 2, &c, 1, &scale));
    scale = 0;
    try trsyl(al, check.io, types.Complex64, .c, .c, 1, 0, 2, &z, 1, &z, 2, &zc, 1, &scale);
    try testing.expectEqual(1, scale);
    for (zc) |v| try testing.expect(lcheck.isNanT(types.Complex64, v));
}

const BigSyl = struct {
    a: []const f64,
    b: []const f64,
    c: []const f64,

    fn run(self: BigSyl, io: Io, out: []f64) anyerror!void {
        @memcpy(out, self.c);
        var scale: f64 = 0;
        try trsyl(testing.allocator, io, f64, .n, .t, 1, lcheck.big_n, lcheck.big_n, self.a.ptr, lcheck.big_n, self.b.ptr, lcheck.big_n, out.ptr, lcheck.big_n, &scale);
    }
};

test "relapack trsyl: the same bits inline and threaded" {
    const al = testing.allocator;
    const a = try schur(f64, al, lcheck.big_n, 2, 3, 1);
    defer al.free(a);
    const b = try schur(f64, al, lcheck.big_n, 2, 5, 2);
    defer al.free(b);
    const c = try al.alloc(f64, lcheck.big_n * lcheck.big_n);
    defer al.free(c);
    check.fill(f64, c, 3);
    try lcheck.expectIoBits(f64, BigSyl{ .a = a, .b = b, .c = c }, BigSyl.run, lcheck.big_n * lcheck.big_n, "trsyl");
}
