//! Batched `gemm`, column-major: OpenBLAS's `interface/gemm_batch.c` (groups
//! of problems that share their sizes and scalars) and
//! `gemm_batch_strided.c` (one problem repeated at fixed strides), for any
//! type `gemm` takes and for the half gemms (upstream has `sbgemm`'s; `bgemm`
//! and `shgemm` come along).
//!
//! As upstream, every group's arguments are checked before any problem
//! runs, so a bad one leaves every `C` as it was. Upstream hands whole
//! problems to its threads (`gemm_batch_thread.c`), each run serially; here
//! the problems run in order, each threading as `gemm` does. The bits are
//! `gemm`'s either way. Upstream also takes its small-matrix kernels for
//! problems up to `100^3` multiply-adds; rstd has none yet (T30).

const std = @import("std");
const interface = @import("interface.zig");
const half = @import("half.zig");
const types = @import("../types.zig");
const util = @import("../util.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Trans = types.Trans;
const HalfGemm = half.HalfGemm;

/// What the batched routines return: `gemm`'s errors, and `InvalidCount`
/// for a group whose `a`, `b` and `c` differ in length.
pub const Error = interface.Error || error{InvalidCount};

/// One group of `gemmBatch`: `C := alpha*op(A)*op(B) + beta*C` for each
/// `a[i]`, `b[i]`, `c[i]`, all with these sizes, scalars and leading
/// dimensions. `In` is the type of `A` and `B`, `Out` that of `C` and the
/// scalars.
pub fn Group(comptime In: type, comptime Out: type) type {
    return struct {
        transa: Trans = .n,
        transb: Trans = .n,
        m: isize,
        n: isize,
        k: isize,
        alpha: Out,
        a: []const [*]const In,
        lda: isize,
        b: []const [*]const In,
        ldb: isize,
        beta: Out,
        c: []const [*]Out,
        ldc: isize,
    };
}

/// `gemm`'s checks for one problem shape, in `xerbla`'s order.
fn checkOne(transa: Trans, transb: Trans, m: isize, n: isize, k: isize, lda: isize, ldb: isize, ldc: isize, size_in: usize, size_out: usize) interface.ArgError!void {
    if (transa.conjX() or transb.conjX()) return error.InvalidTrans;
    const ta = transa.transposed();
    const tb = transb.transposed();
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    if (k < 0) return error.InvalidK;
    if (lda < @max(1, if (ta) k else m)) return error.InvalidLda;
    if (ldb < @max(1, if (tb) n else k)) return error.InvalidLdb;
    if (ldc < @max(1, m)) return error.InvalidLdc;
    if (!util.fitsMatrix(if (ta) k else m, if (ta) m else k, lda, size_in)) return error.InvalidLda;
    if (!util.fitsMatrix(if (tb) n else k, if (tb) k else n, ldb, size_in)) return error.InvalidLdb;
    if (!util.fitsMatrix(m, n, ldc, size_out)) return error.InvalidLdc;
}

/// The routine one problem runs: `gemm` over `T`, or a half gemm.
fn Runner(comptime In: type, comptime Out: type, comptime kind: ?HalfGemm) type {
    return struct {
        fn one(allocator: Allocator, io: Io, transa: Trans, transb: Trans, m: isize, n: isize, k: isize, alpha: Out, a: [*]const In, lda: isize, b: [*]const In, ldb: isize, beta: Out, c: [*]Out, ldc: isize) interface.Error!void {
            if (comptime kind) |h| return half.gemmHalf(allocator, io, h, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
            return interface.gemm(allocator, io, In, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
        }

        fn groups(allocator: Allocator, io: Io, gs: []const Group(In, Out)) Error!void {
            for (gs) |g| {
                try checkOne(g.transa, g.transb, g.m, g.n, g.k, g.lda, g.ldb, g.ldc, @sizeOf(In), @sizeOf(Out));
                if (g.b.len != g.a.len or g.c.len != g.a.len) return error.InvalidCount;
            }
            for (gs) |g| for (g.a, g.b, g.c) |a, b, c| {
                try one(allocator, io, g.transa, g.transb, g.m, g.n, g.k, g.alpha, a, g.lda, b, g.ldb, g.beta, c, g.ldc);
            };
        }

        fn strided(allocator: Allocator, io: Io, transa: Trans, transb: Trans, m: isize, n: isize, k: isize, alpha: Out, a: [*]const In, lda: isize, stride_a: isize, b: [*]const In, ldb: isize, stride_b: isize, beta: Out, c: [*]Out, ldc: isize, stride_c: isize, count: isize) Error!void {
            try checkOne(transa, transb, m, n, k, lda, ldb, ldc, @sizeOf(In), @sizeOf(Out));
            var i: isize = 0;
            while (i < count) : (i += 1) {
                const pa = util.at(In, @constCast(a), i * stride_a);
                const pb = util.at(In, @constCast(b), i * stride_b);
                const pc = util.at(Out, c, i * stride_c);
                try one(allocator, io, transa, transb, m, n, k, alpha, @ptrCast(pa), lda, @ptrCast(pb), ldb, beta, @ptrCast(pc), ldc);
            }
        }
    };
}

/// `?GEMM_BATCH` (`interface/gemm_batch.c`): each group's problems, after
/// every group is checked.
pub fn gemmBatch(allocator: Allocator, io: Io, comptime T: type, groups: []const Group(T, T)) Error!void {
    return Runner(T, T, null).groups(allocator, io, groups);
}

/// `?GEMM_BATCH_STRIDED` (`interface/gemm_batch_strided.c`): `count`
/// problems, the `i`th at `a + i*stride_a`, `b + i*stride_b` and `c +
/// i*stride_c`; nothing for a `count` below 1. The strides are not checked,
/// as upstream's are not.
pub fn gemmBatchStrided(allocator: Allocator, io: Io, comptime T: type, transa: Trans, transb: Trans, m: isize, n: isize, k: isize, alpha: T, a: [*]const T, lda: isize, stride_a: isize, b: [*]const T, ldb: isize, stride_b: isize, beta: T, c: [*]T, ldc: isize, stride_c: isize, count: isize) Error!void {
    return Runner(T, T, null).strided(allocator, io, transa, transb, m, n, k, alpha, a, lda, stride_a, b, ldb, stride_b, beta, c, ldc, stride_c, count);
}

/// `SBGEMM_BATCH`, and the same for `bgemm` and `shgemm`: `gemmBatch` over
/// `gemmHalf`.
pub fn gemmHalfBatch(allocator: Allocator, io: Io, comptime kind: HalfGemm, groups: []const Group(kind.In(), kind.Out())) Error!void {
    return Runner(kind.In(), kind.Out(), kind).groups(allocator, io, groups);
}

/// `SBGEMM_BATCH_STRIDED`, and the same for `bgemm` and `shgemm`.
pub fn gemmHalfBatchStrided(allocator: Allocator, io: Io, comptime kind: HalfGemm, transa: Trans, transb: Trans, m: isize, n: isize, k: isize, alpha: kind.Out(), a: [*]const kind.In(), lda: isize, stride_a: isize, b: [*]const kind.In(), ldb: isize, stride_b: isize, beta: kind.Out(), c: [*]kind.Out(), ldc: isize, stride_c: isize, count: isize) Error!void {
    return Runner(kind.In(), kind.Out(), kind).strided(allocator, io, transa, transb, m, n, k, alpha, a, lda, stride_a, b, ldb, stride_b, beta, c, ldc, stride_c, count);
}

// ---- tests -----------------------------------------------------------------------

const testing = std.testing;
const check = @import("../_test/blas1/check.zig");
const bf16_kernel = @import("../ext/bf16/kernel.zig");

fn val(comptime T: type, v: f32) T {
    if (comptime types.isComplex(T)) return .init(@floatCast(v), @floatCast(v / 2));
    if (T == types.BF16 or T == types.HF16) return bf16_kernel.narrow(T, v);
    return @floatCast(v);
}

fn sameBits(comptime T: type, x: T, y: T) bool {
    return std.mem.eql(u8, std.mem.asBytes(&x), std.mem.asBytes(&y));
}

/// Two groups, of 3 and 2 problems, against each problem's own call.
fn groupsCase(comptime In: type, comptime Out: type, comptime kind: ?HalfGemm) !void {
    const R = Runner(In, Out, kind);
    const al = testing.allocator;
    const m = 5;
    const n = 4;
    const k = 6;
    var a: [5][m * k]In = undefined;
    var b: [5][k * n]In = undefined;
    var c: [5][m * n]Out = undefined;
    for (0..5) |p| {
        for (&a[p], 0..) |*v, i| v.* = val(In, @floatFromInt((i * 7 + p * 3) % 11));
        for (&b[p], 0..) |*v, i| v.* = val(In, @floatFromInt((i * 5 + p) % 9));
        for (&c[p], 0..) |*v, i| v.* = val(Out, @floatFromInt((i + p) % 5));
    }
    const orig = c;
    var want = c;
    // group 0: `A` `m x k`, `B` `k x n`; group 1: `A^T`, `B^T` on the same storage
    try R.one(al, check.io, .n, .n, m, n, k, val(Out, 1.5), &a[0], m, &b[0], k, val(Out, 0.5), &want[0], m);
    try R.one(al, check.io, .n, .n, m, n, k, val(Out, 1.5), &a[1], m, &b[1], k, val(Out, 0.5), &want[1], m);
    try R.one(al, check.io, .n, .n, m, n, k, val(Out, 1.5), &a[2], m, &b[2], k, val(Out, 0.5), &want[2], m);
    try R.one(al, check.io, .t, .t, 4, 4, 5, val(Out, -1), &a[3], 5, &b[3], 4, val(Out, 0), &want[3], 5);
    try R.one(al, check.io, .t, .t, 4, 4, 5, val(Out, -1), &a[4], 5, &b[4], 4, val(Out, 0), &want[4], 5);
    const gs = [_]Group(In, Out){
        .{ .m = m, .n = n, .k = k, .alpha = val(Out, 1.5), .a = &.{ &a[0], &a[1], &a[2] }, .lda = m, .b = &.{ &b[0], &b[1], &b[2] }, .ldb = k, .beta = val(Out, 0.5), .c = &.{ &c[0], &c[1], &c[2] }, .ldc = m },
        .{ .transa = .t, .transb = .t, .m = 4, .n = 4, .k = 5, .alpha = val(Out, -1), .a = &.{ &a[3], &a[4] }, .lda = 5, .b = &.{ &b[3], &b[4] }, .ldb = 4, .beta = val(Out, 0), .c = &.{ &c[3], &c[4] }, .ldc = 5 },
    };
    try R.groups(al, check.io, &gs);
    for (want, c) |w, g| for (w, g) |x, y| try testing.expect(sameBits(Out, x, y));
    // a bad later group runs nothing, nor does a group of unequal lengths
    const before = c;
    var bad = gs;
    bad[1].ldc = 3;
    try testing.expectError(error.InvalidLdc, R.groups(al, check.io, &bad));
    bad = gs;
    bad[1].c = bad[1].c[0..1];
    try testing.expectError(error.InvalidCount, R.groups(al, check.io, &bad));
    for (before, c) |w, g| for (w, g) |x, y| try testing.expect(sameBits(Out, x, y));
    // strided: the three problems of group 0 laid end to end, `C` by a
    // negative stride
    var sa: [3 * m * k]In = undefined;
    var sb: [3 * k * n]In = undefined;
    var sc: [3 * m * n]Out = undefined;
    for (0..3) |p| {
        @memcpy(sa[p * m * k ..][0 .. m * k], &a[p]);
        @memcpy(sb[p * k * n ..][0 .. k * n], &b[p]);
        @memcpy(sc[(2 - p) * m * n ..][0 .. m * n], &orig[p]);
    }
    try R.strided(al, check.io, .n, .n, m, n, k, val(Out, 1.5), &sa, m, m * k, &sb, k, k * n, val(Out, 0.5), sc[2 * m * n ..].ptr, m, -(m * n), 3);
    for (0..3) |p| for (sc[(2 - p) * m * n ..][0 .. m * n], c[p]) |x, y| try testing.expect(sameBits(Out, x, y));
    try testing.expectError(error.InvalidK, R.strided(al, check.io, .n, .n, m, n, -1, val(Out, 1), &sa, m, 0, &sb, k, 0, val(Out, 0), &sc, m, 0, 3));
    try R.strided(al, check.io, .n, .n, m, n, k, val(Out, 1), &sa, m, 0, &sb, k, 0, val(Out, 0), &sc, m, 0, 0);
}

test "gemm batches: groups and strides give each problem's own gemm bits, and check every group first" {
    inline for (.{ f32, f64, types.Complex32, types.Complex64 }) |T| try groupsCase(T, T, null);
    inline for (std.enums.values(HalfGemm)) |kind| try groupsCase(kind.In(), kind.Out(), kind);
}
