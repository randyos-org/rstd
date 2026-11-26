//! The public level-3 routines (issue #311): matrices as a slice and its
//! leading dimension, column-major, sizes `usize`. Each checks that every
//! slice holds its matrix (`util.holds`), by the shape its `trans` and `side`
//! letters give it, then delegates to the reference-shaped core
//! (`interface.zig`, `half.zig`, `batch.zig`), which the LAPACK routines and
//! the ported tests call.

const std = @import("std");
const types = @import("../types.zig");
const strided = @import("../strided.zig");
const util = @import("../util.zig");
const iface = @import("interface.zig");
const half = @import("half.zig");
const batch = @import("batch.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const Trans = types.Trans;
const Uplo = types.Uplo;
const Diag = types.Diag;
const Side = types.Side;
const holds = util.holds;

pub const HalfPrecision = half.HalfPrecision;

/// A level-3 routine's bad input: a slice too short for its matrix, a leading
/// dimension below its rows; and `error.OutOfMemory` and `error.Canceled`.
pub const Error = strided.SliceError || strided.LeadingDimensionError || Allocator.Error || Io.Cancelable;
/// `gemmtHalf`'s, which allocates nothing.
pub const HalfError = strided.SliceError || strided.LeadingDimensionError || Io.Cancelable;
/// The batched routines': `Error`, and `error.LengthMismatch` for a group whose
/// lists of `A`, `B` and `C` differ in length.
pub const BatchError = Error || strided.LengthError;

fn len(n: usize) isize {
    return @intCast(n);
}

/// The stored rows and columns of a matrix `op` makes `rows` by `cols`.
fn stored(trans: Trans, rows: usize, cols: usize) [2]usize {
    return if (trans.transposed()) .{ cols, rows } else .{ rows, cols };
}

/// `a` holds `op(A)`, `rows` by `cols`.
fn holdsOp(a_len: usize, trans: Trans, rows: usize, cols: usize, lda: usize) (strided.SliceError || strided.LeadingDimensionError)!void {
    const s = stored(trans, rows, cols);
    return holds(a_len, s[0], s[1], lda);
}

/// The core's errors past these checks: the allocator's, cancellation, and
/// the leading dimensions should the core's rule be stricter than `holds`.
fn mapped(e: iface.Error) Error {
    return switch (e) {
        error.OutOfMemory => error.OutOfMemory,
        error.Canceled => error.Canceled,
        error.InvalidLda, error.InvalidLdb, error.InvalidLdc => error.LeadingDimensionTooSmall,
        error.InvalidUplo, error.InvalidTrans, error.InvalidM, error.InvalidN, error.InvalidK => unreachable,
    };
}

fn mappedHalf(e: half.HalfError) HalfError {
    return switch (e) {
        error.Canceled => error.Canceled,
        error.InvalidLda, error.InvalidLdb, error.InvalidLdc => error.LeadingDimensionTooSmall,
        error.InvalidUplo, error.InvalidTrans, error.InvalidM, error.InvalidN, error.InvalidK => unreachable,
    };
}

fn mappedBatch(e: batch.Error) BatchError {
    return switch (e) {
        error.InvalidCount => error.LengthMismatch,
        else => |f| mapped(@errorCast(f)),
    };
}

// ---- gemm --------------------------------------------------------------------

/// `C := alpha op(A) op(B) + beta C`, `op(A)` `m` by `k`, `op(B)` `k` by `n`.
/// A zero `beta` never reads `C`; a zero `alpha` or `k` reads neither `A` nor `B`.
pub fn gemm(allocator: Allocator, io: Io, comptime T: type, transa: Trans, transb: Trans, m: usize, n: usize, k: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: T, c: []T, ldc: usize) Error!void {
    try holdsOp(a.len, transa, m, k, lda);
    try holdsOp(b.len, transb, k, n, ldb);
    try holds(c.len, m, n, ldc);
    iface.gemm(allocator, io, T, transa, transb, len(m), len(n), len(k), alpha, a.ptr, len(lda), b.ptr, len(ldb), beta, c.ptr, len(ldc)) catch |e| return mapped(e);
}

/// `gemm` by three real products (complex only).
pub fn gemm3m(allocator: Allocator, io: Io, comptime T: type, transa: Trans, transb: Trans, m: usize, n: usize, k: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: T, c: []T, ldc: usize) Error!void {
    try holdsOp(a.len, transa, m, k, lda);
    try holdsOp(b.len, transb, k, n, ldb);
    try holds(c.len, m, n, ldc);
    iface.gemm3m(allocator, io, T, transa, transb, len(m), len(n), len(k), alpha, a.ptr, len(lda), b.ptr, len(ldb), beta, c.ptr, len(ldc)) catch |e| return mapped(e);
}

/// `gemm` into `uplo`'s triangle of the `m` by `m` `C` only.
pub fn gemmt(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, transa: Trans, transb: Trans, m: usize, k: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: T, c: []T, ldc: usize) Error!void {
    try holdsOp(a.len, transa, m, k, lda);
    try holdsOp(b.len, transb, k, m, ldb);
    try holds(c.len, m, m, ldc);
    iface.gemmt(allocator, io, T, uplo, transa, transb, len(m), len(k), alpha, a.ptr, len(lda), b.ptr, len(ldb), beta, c.ptr, len(ldc)) catch |e| return mapped(e);
}

/// `gemm` in half precision: `kind` picks bfloat16 or `f16` in, `f32` or
/// bfloat16 out.
pub fn gemmHalf(allocator: Allocator, io: Io, comptime kind: HalfPrecision, transa: Trans, transb: Trans, m: usize, n: usize, k: usize, alpha: kind.Out(), a: []const kind.In(), lda: usize, b: []const kind.In(), ldb: usize, beta: kind.Out(), c: []kind.Out(), ldc: usize) Error!void {
    try holdsOp(a.len, transa, m, k, lda);
    try holdsOp(b.len, transb, k, n, ldb);
    try holds(c.len, m, n, ldc);
    half.gemmHalf(allocator, io, kind, transa, transb, len(m), len(n), len(k), alpha, a.ptr, len(lda), b.ptr, len(ldb), beta, c.ptr, len(ldc)) catch |e| return mapped(e);
}

/// `gemmt` in half precision.
pub fn gemmtHalf(io: Io, comptime kind: HalfPrecision, uplo: Uplo, transa: Trans, transb: Trans, m: usize, k: usize, alpha: kind.Out(), a: []const kind.In(), lda: usize, b: []const kind.In(), ldb: usize, beta: kind.Out(), c: []kind.Out(), ldc: usize) HalfError!void {
    try holdsOp(a.len, transa, m, k, lda);
    try holdsOp(b.len, transb, k, m, ldb);
    try holds(c.len, m, m, ldc);
    half.gemmtHalf(io, kind, uplo, transa, transb, len(m), len(k), alpha, a.ptr, len(lda), b.ptr, len(ldb), beta, c.ptr, len(ldc)) catch |e| return mappedHalf(e);
}

// ---- symmetric and Hermitian -------------------------------------------------

/// `A`'s order: `m` on the left, `n` on the right.
fn order(side: Side, m: usize, n: usize) usize {
    return if (side == .l) m else n;
}

/// `C := alpha A B + beta C` (`side` `l`) or `alpha B A + beta C` (`r`), `A`
/// symmetric (`uplo`'s triangle), `B` and `C` `m` by `n`.
pub fn symm(allocator: Allocator, io: Io, comptime T: type, side: Side, uplo: Uplo, m: usize, n: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: T, c: []T, ldc: usize) Error!void {
    try holds(a.len, order(side, m, n), order(side, m, n), lda);
    try holds(b.len, m, n, ldb);
    try holds(c.len, m, n, ldc);
    iface.symm(allocator, io, T, side, uplo, len(m), len(n), alpha, a.ptr, len(lda), b.ptr, len(ldb), beta, c.ptr, len(ldc)) catch |e| return mapped(e);
}

/// `symm` with `A` Hermitian. Complex only.
pub fn hemm(allocator: Allocator, io: Io, comptime T: type, side: Side, uplo: Uplo, m: usize, n: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: T, c: []T, ldc: usize) Error!void {
    try holds(a.len, order(side, m, n), order(side, m, n), lda);
    try holds(b.len, m, n, ldb);
    try holds(c.len, m, n, ldc);
    iface.hemm(allocator, io, T, side, uplo, len(m), len(n), alpha, a.ptr, len(lda), b.ptr, len(ldb), beta, c.ptr, len(ldc)) catch |e| return mapped(e);
}

/// `symm` by three real products (complex only).
pub fn symm3m(allocator: Allocator, io: Io, comptime T: type, side: Side, uplo: Uplo, m: usize, n: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: T, c: []T, ldc: usize) Error!void {
    try holds(a.len, order(side, m, n), order(side, m, n), lda);
    try holds(b.len, m, n, ldb);
    try holds(c.len, m, n, ldc);
    iface.symm3m(allocator, io, T, side, uplo, len(m), len(n), alpha, a.ptr, len(lda), b.ptr, len(ldb), beta, c.ptr, len(ldc)) catch |e| return mapped(e);
}

/// `hemm` by three real products. Complex only.
pub fn hemm3m(allocator: Allocator, io: Io, comptime T: type, side: Side, uplo: Uplo, m: usize, n: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: T, c: []T, ldc: usize) Error!void {
    try holds(a.len, order(side, m, n), order(side, m, n), lda);
    try holds(b.len, m, n, ldb);
    try holds(c.len, m, n, ldc);
    iface.hemm3m(allocator, io, T, side, uplo, len(m), len(n), alpha, a.ptr, len(lda), b.ptr, len(ldb), beta, c.ptr, len(ldc)) catch |e| return mapped(e);
}

/// `C := alpha op(A) op(A)ᵀ + beta C`, `op(A)` `n` by `k`, `C` symmetric.
pub fn syrk(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, n: usize, k: usize, alpha: T, a: []const T, lda: usize, beta: T, c: []T, ldc: usize) Error!void {
    try holdsOp(a.len, trans, n, k, lda);
    try holds(c.len, n, n, ldc);
    iface.syrk(allocator, io, T, uplo, trans, len(n), len(k), alpha, a.ptr, len(lda), beta, c.ptr, len(ldc)) catch |e| return mapped(e);
}

/// `C := alpha op(A) op(A)ᴴ + beta C`, real `alpha` and `beta`, `C`
/// Hermitian. Complex only.
pub fn herk(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, n: usize, k: usize, alpha: Float(T), a: []const T, lda: usize, beta: Float(T), c: []T, ldc: usize) Error!void {
    try holdsOp(a.len, trans, n, k, lda);
    try holds(c.len, n, n, ldc);
    iface.herk(allocator, io, T, uplo, trans, len(n), len(k), alpha, a.ptr, len(lda), beta, c.ptr, len(ldc)) catch |e| return mapped(e);
}

/// `C := alpha op(A) op(B)ᵀ + alpha op(B) op(A)ᵀ + beta C`, `C` symmetric.
pub fn syr2k(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, n: usize, k: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: T, c: []T, ldc: usize) Error!void {
    try holdsOp(a.len, trans, n, k, lda);
    try holdsOp(b.len, trans, n, k, ldb);
    try holds(c.len, n, n, ldc);
    iface.syr2k(allocator, io, T, uplo, trans, len(n), len(k), alpha, a.ptr, len(lda), b.ptr, len(ldb), beta, c.ptr, len(ldc)) catch |e| return mapped(e);
}

/// `C := alpha op(A) op(B)ᴴ + conj(alpha) op(B) op(A)ᴴ + beta C`, real
/// `beta`, `C` Hermitian. Complex only.
pub fn her2k(allocator: Allocator, io: Io, comptime T: type, uplo: Uplo, trans: Trans, n: usize, k: usize, alpha: T, a: []const T, lda: usize, b: []const T, ldb: usize, beta: Float(T), c: []T, ldc: usize) Error!void {
    try holdsOp(a.len, trans, n, k, lda);
    try holdsOp(b.len, trans, n, k, ldb);
    try holds(c.len, n, n, ldc);
    iface.her2k(allocator, io, T, uplo, trans, len(n), len(k), alpha, a.ptr, len(lda), b.ptr, len(ldb), beta, c.ptr, len(ldc)) catch |e| return mapped(e);
}

// ---- triangular ----------------------------------------------------------------

/// `B := alpha op(A) B` (`side` `l`) or `alpha B op(A)` (`r`), `A`
/// triangular, `B` `m` by `n`.
pub fn trmm(allocator: Allocator, io: Io, comptime T: type, side: Side, uplo: Uplo, transa: Trans, diag: Diag, m: usize, n: usize, alpha: T, a: []const T, lda: usize, b: []T, ldb: usize) Error!void {
    try holds(a.len, order(side, m, n), order(side, m, n), lda);
    try holds(b.len, m, n, ldb);
    iface.trmm(allocator, io, T, side, uplo, transa, diag, len(m), len(n), alpha, a.ptr, len(lda), b.ptr, len(ldb)) catch |e| return mapped(e);
}

/// `B := alpha op(A)⁻¹ B` (`side` `l`) or `alpha B op(A)⁻¹` (`r`).
pub fn trsm(allocator: Allocator, io: Io, comptime T: type, side: Side, uplo: Uplo, transa: Trans, diag: Diag, m: usize, n: usize, alpha: T, a: []const T, lda: usize, b: []T, ldb: usize) Error!void {
    try holds(a.len, order(side, m, n), order(side, m, n), lda);
    try holds(b.len, m, n, ldb);
    iface.trsm(allocator, io, T, side, uplo, transa, diag, len(m), len(n), alpha, a.ptr, len(lda), b.ptr, len(ldb)) catch |e| return mapped(e);
}

// ---- batched -----------------------------------------------------------------

/// One group of `gemmBatch` or `gemmHalfBatch`: one shape and scalars, and
/// one `A`, `B` and `C` per problem, the three lists as long.
pub fn GemmGroup(comptime In: type, comptime Out: type) type {
    return struct {
        transa: Trans = .n,
        transb: Trans = .n,
        m: usize,
        n: usize,
        k: usize,
        alpha: Out,
        a: []const []const In,
        lda: usize,
        b: []const []const In,
        ldb: usize,
        beta: Out,
        c: []const []Out,
        ldc: usize,
    };
}

/// `groups` as the core's, every matrix checked; the pointer lists are one
/// allocation, which the caller frees.
fn coreGroups(comptime In: type, comptime Out: type, allocator: Allocator, groups: []const GemmGroup(In, Out)) BatchError!struct { []batch.GemmGroup(In, Out), [][*]const In, [][*]Out } {
    var problems: usize = 0;
    for (groups) |g| {
        if (g.b.len != g.a.len or g.c.len != g.a.len) return error.LengthMismatch;
        for (g.a, g.b, g.c) |a, b, c| {
            try holdsOp(a.len, g.transa, g.m, g.k, g.lda);
            try holdsOp(b.len, g.transb, g.k, g.n, g.ldb);
            try holds(c.len, g.m, g.n, g.ldc);
        }
        problems += g.a.len;
    }
    const out = try allocator.alloc(batch.GemmGroup(In, Out), groups.len);
    errdefer allocator.free(out);
    const ins = try allocator.alloc([*]const In, 2 * problems);
    errdefer allocator.free(ins);
    const outs = try allocator.alloc([*]Out, problems);
    var at: usize = 0;
    for (groups, out) |g, *o| {
        const p = g.a.len;
        const ga = ins[2 * at ..][0..p];
        const gb = ins[2 * at + p ..][0..p];
        const gc = outs[at..][0..p];
        for (g.a, g.b, g.c, ga, gb, gc) |a, b, c, *pa, *pb, *pc| {
            pa.* = a.ptr;
            pb.* = b.ptr;
            pc.* = c.ptr;
        }
        o.* = .{ .transa = g.transa, .transb = g.transb, .m = len(g.m), .n = len(g.n), .k = len(g.k), .alpha = g.alpha, .a = ga, .lda = len(g.lda), .b = gb, .ldb = len(g.ldb), .beta = g.beta, .c = gc, .ldc = len(g.ldc) };
        at += p;
    }
    return .{ out, ins, outs };
}

/// Many `gemm`s, by group: each group one shape and scalars over its lists of
/// matrices, the problems spread over `io`.
pub fn gemmBatch(allocator: Allocator, io: Io, comptime T: type, groups: []const GemmGroup(T, T)) BatchError!void {
    const g, const ins, const outs = try coreGroups(T, T, allocator, groups);
    defer allocator.free(g);
    defer allocator.free(ins);
    defer allocator.free(outs);
    batch.gemmBatch(allocator, io, T, g) catch |e| return mappedBatch(e);
}

/// `gemmBatch` in half precision.
pub fn gemmHalfBatch(allocator: Allocator, io: Io, comptime kind: HalfPrecision, groups: []const GemmGroup(kind.In(), kind.Out())) BatchError!void {
    const g, const ins, const outs = try coreGroups(kind.In(), kind.Out(), allocator, groups);
    defer allocator.free(g);
    defer allocator.free(ins);
    defer allocator.free(outs);
    batch.gemmHalfBatch(allocator, io, kind, g) catch |e| return mappedBatch(e);
}

/// `x` holds `count` matrices `stride` apart, each `rows` by `cols` at `ld`.
fn holdsStrided(x_len: usize, rows: usize, cols: usize, ld: usize, stride: usize, count: usize) (strided.SliceError || strided.LeadingDimensionError)!void {
    try holds(x_len, rows, cols, ld);
    if (count <= 1 or rows == 0 or cols == 0) return;
    const one = (cols - 1) * ld + rows;
    const last = std.math.mul(usize, count - 1, stride) catch return error.SliceTooShort;
    const need = std.math.add(usize, last, one) catch return error.SliceTooShort;
    if (need > x_len) return error.SliceTooShort;
}

/// `count` `gemm`s of one shape, the matrices `stride_a`, `stride_b` and
/// `stride_c` apart in their slices.
pub fn gemmBatchStrided(allocator: Allocator, io: Io, comptime T: type, transa: Trans, transb: Trans, m: usize, n: usize, k: usize, alpha: T, a: []const T, lda: usize, stride_a: usize, b: []const T, ldb: usize, stride_b: usize, beta: T, c: []T, ldc: usize, stride_c: usize, count: usize) Error!void {
    const sa = stored(transa, m, k);
    const sb = stored(transb, k, n);
    try holdsStrided(a.len, sa[0], sa[1], lda, stride_a, count);
    try holdsStrided(b.len, sb[0], sb[1], ldb, stride_b, count);
    try holdsStrided(c.len, m, n, ldc, stride_c, count);
    batch.gemmBatchStrided(allocator, io, T, transa, transb, len(m), len(n), len(k), alpha, a.ptr, len(lda), len(stride_a), b.ptr, len(ldb), len(stride_b), beta, c.ptr, len(ldc), len(stride_c), len(count)) catch |e| return mapped(@errorCast(e));
}

/// `gemmBatchStrided` in half precision.
pub fn gemmHalfBatchStrided(allocator: Allocator, io: Io, comptime kind: HalfPrecision, transa: Trans, transb: Trans, m: usize, n: usize, k: usize, alpha: kind.Out(), a: []const kind.In(), lda: usize, stride_a: usize, b: []const kind.In(), ldb: usize, stride_b: usize, beta: kind.Out(), c: []kind.Out(), ldc: usize, stride_c: usize, count: usize) Error!void {
    const sa = stored(transa, m, k);
    const sb = stored(transb, k, n);
    try holdsStrided(a.len, sa[0], sa[1], lda, stride_a, count);
    try holdsStrided(b.len, sb[0], sb[1], ldb, stride_b, count);
    try holdsStrided(c.len, m, n, ldc, stride_c, count);
    batch.gemmHalfBatchStrided(allocator, io, kind, transa, transb, len(m), len(n), len(k), alpha, a.ptr, len(lda), len(stride_a), b.ptr, len(ldb), len(stride_b), beta, c.ptr, len(ldc), len(stride_c), len(count)) catch |e| return mapped(@errorCast(e));
}

// ---- tests -------------------------------------------------------------------

const check = @import("../_test/check.zig");
const testing = std.testing;

test "api: gemm on slices gives the core's bits, transposed too" {
    const al = testing.allocator;
    var a: [12]f64 = undefined;
    check.fill(f64, &a, 1);
    var b: [15]f64 = undefined;
    check.fill(f64, &b, 2);
    var c1: [20]f64 = undefined;
    check.fill(f64, &c1, 3);
    var c2 = c1;
    // op(A) 4 x 3, A stored 3 x 4 (transposed); B 3 x 5
    try gemm(al, check.io, f64, .t, .n, 4, 5, 3, 1.5, &a, 3, &b, 3, 0.5, &c1, 4);
    try iface.gemm(al, check.io, f64, .t, .n, 4, 5, 3, 1.5, &a, 3, &b, 3, 0.5, &c2, 4);
    try testing.expectEqualSlices(f64, &c2, &c1);
}

test "api: a short slice or a small leading dimension is an error, by side and trans" {
    const al = testing.allocator;
    var a: [16]f64 = @splat(1);
    var b: [12]f64 = @splat(1);
    var c: [12]f64 = @splat(0);
    try testing.expectError(error.SliceTooShort, gemm(al, check.io, f64, .n, .n, 4, 3, 4, 1, a[0..15], 4, &b, 4, 0, &c, 4));
    try testing.expectError(error.LeadingDimensionTooSmall, gemm(al, check.io, f64, .n, .n, 4, 3, 4, 1, &a, 3, &b, 4, 0, &c, 4));
    // B 4 x 3 on the right needs a 3 x 3 A, on the left a 4 x 4 one
    try trsm(al, check.io, f64, .r, .u, .n, .n, 4, 3, 1, a[0..9], 3, &b, 4);
    try testing.expectError(error.SliceTooShort, trsm(al, check.io, f64, .l, .u, .n, .n, 4, 3, 1, a[0..9], 4, &b, 4));
    try testing.expectError(error.SliceTooShort, syrk(al, check.io, f64, .u, .t, 3, 4, 1, a[0..11], 4, 0, c[0..9], 3));
}

test "api: gemmBatch is each problem's gemm, and mismatched lists are an error" {
    const al = testing.allocator;
    var a: [2][4]f64 = undefined;
    var b: [2][4]f64 = undefined;
    var c: [2][4]f64 = @splat(@splat(0));
    for (&a, &b, 0..) |*x, *y, i| {
        check.fill(f64, x, @intCast(2 * i + 1));
        check.fill(f64, y, @intCast(2 * i + 2));
    }
    const as = [_][]const f64{ &a[0], &a[1] };
    const bs = [_][]const f64{ &b[0], &b[1] };
    const cs = [_][]f64{ &c[0], &c[1] };
    const g = [_]GemmGroup(f64, f64){.{ .m = 2, .n = 2, .k = 2, .alpha = 1, .a = &as, .lda = 2, .b = &bs, .ldb = 2, .beta = 0, .c = &cs, .ldc = 2 }};
    try gemmBatch(al, check.io, f64, &g);
    for (0..2) |i| {
        var want: [4]f64 = @splat(0);
        try iface.gemm(al, check.io, f64, .n, .n, 2, 2, 2, 1, &a[i], 2, &b[i], 2, 0, &want, 2);
        try testing.expectEqualSlices(f64, &want, &c[i]);
    }
    const short = [_]GemmGroup(f64, f64){.{ .m = 2, .n = 2, .k = 2, .alpha = 1, .a = &as, .lda = 2, .b = bs[0..1], .ldb = 2, .beta = 0, .c = &cs, .ldc = 2 }};
    try testing.expectError(error.LengthMismatch, gemmBatch(al, check.io, f64, &short));
    // two A 4 apart, one B (stride 0); C one element short of its second matrix
    const af: [8]f64 = @bitCast(a);
    var flat: [8]f64 = @splat(0);
    try testing.expectError(error.SliceTooShort, gemmBatchStrided(al, check.io, f64, .n, .n, 2, 2, 2, 1, &af, 2, 4, &b[0], 2, 0, 0, flat[0..7], 2, 4, 2));
    try gemmBatchStrided(al, check.io, f64, .n, .n, 2, 2, 2, 1, &af, 2, 4, &b[0], 2, 0, 0, &flat, 2, 4, 2);
}
