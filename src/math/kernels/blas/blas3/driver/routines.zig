//! The level-3 routines on the drivers: `gemm`, `gemmt`, `symm` and `syrk` (and their Hermitian forms), checked by the interface.

const std = @import("std");
const types = @import("../../types.zig");
const Vtable = @import("../../Vtable.zig");
const portable = @import("../../portable.zig");
const gemm_kernels = @import("../gemm/kernel.zig");
const small = @import("../gemm/small.zig");
const sme2 = @import("../gemm/sme2/driver.zig");
const droot = @import("../driver.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const isComplex = types.isComplex;
pub const Conj = gemm_kernels.Conj;
const Desc = droot.Desc;
const Error = droot.Error;
const Matrix = droot.Matrix;
const Pack = droot.Pack;
const Region = droot.Region;
const isOne = droot.isOne;
const isZero = droot.isZero;
const raw = droot.raw;
const run = droot.run;
const scale = droot.scale;

// ---- the routines ----------------------------------------------------------------

/// `C := alpha*op(A)*op(B) + beta*C` with the kernels of `vt`; `ta`, `tb`
/// say whether `op` transposes, `conj` which it conjugates (complex).
/// Arguments are checked by the caller (`interface.zig`), and `m`, `n` are
/// positive. Up to `want` workers run it (`threaded.gemmWorkers`).
pub fn gemm(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, ta: bool, tb: bool, conj: Conj, m: usize, n: usize, k: usize, alpha: T, a: [*]const T, lda: usize, b: [*]const T, ldb: usize, beta: T, c: [*]T, ldc: usize, want: u32) Error!void {
    if (k == 0 or isZero(T, alpha)) {
        if (!isOne(T, beta)) scale(T, vt, m, n, beta, c, ldc);
        return;
    }
    const cj = comptime isComplex(T);
    // upstream takes the SME2 kernels before the small-matrix ones
    if (comptime !cj) if (try sme2.gemm(allocator, io, T, vt, ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, want)) return;
    if (small.gemm(T, vt, ta, tb, cj and conj.ofA(), cj and conj.ofB(), m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)) return;
    const x: Matrix(T) = .{ .p = @constCast(a), .ld = lda };
    const y: Matrix(T) = .{ .p = @constCast(b), .ld = ldb };
    const d: Desc = .{ .rows = if (ta) .t else .n, .cols = if (tb) .t else .n, .conj = if (comptime isComplex(T)) conj else .n };
    try run(allocator, io, T, vt, d, m, n, k, alpha, x, y, null, unlessOne(T, beta), c, ldc, want);
}

/// `C += alpha*A*B` on the blocked loops alone, one worker, never the
/// small-matrix path: upstream's LAPACK drivers call the kernels
/// (`getrf_single.c`), so the bits match `drive` on any part of `C`.
pub fn gemmBlocked(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, m: usize, n: usize, k: usize, alpha: T, a: [*]const T, lda: usize, b: [*]const T, ldb: usize, c: [*]T, ldc: usize) Error!void {
    const x: Matrix(T) = .{ .p = @constCast(a), .ld = lda };
    const y: Matrix(T) = .{ .p = @constCast(b), .ld = ldb };
    try run(allocator, io, T, vt, .{ .rows = .n, .cols = .n }, m, n, k, alpha, x, y, null, null, c, ldc, 1);
}

/// `C := alpha*op(A)*op(B) + beta*C` on the `upper` or lower triangle of the
/// `m x m` `C`, `op(A)` being `m x k`: `gemm`'s operands over `syrk`'s
/// triangle. Upstream's `interface/gemmt.c` took a `gemv` per column
/// instead at commit a69062524, conjugating `B` in place for it; this
/// reads `B` only. `m` positive. Up to `want` workers run it, over the rank-k
/// updates' split, as `interface/gemmt.c` takes `syrk_thread` since
/// OpenBLAS PR 6060.
pub fn gemmt(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, upper: bool, ta: bool, tb: bool, conj: Conj, m: usize, k: usize, alpha: T, a: [*]const T, lda: usize, b: [*]const T, ldb: usize, beta: T, c: [*]T, ldc: usize, want: u32) Error!void {
    const region: Region = if (upper) .upper else .lower;
    if (k == 0 or isZero(T, alpha)) {
        if (!isOne(T, beta)) scaleBlock(T, vt, region, false, 0, m, 0, m, beta, c, ldc);
        return;
    }
    const x: Matrix(T) = .{ .p = @constCast(a), .ld = lda };
    const y: Matrix(T) = .{ .p = @constCast(b), .ld = ldb };
    const d: Desc = .{ .rows = if (ta) .t else .n, .cols = if (tb) .t else .n, .conj = if (comptime isComplex(T)) conj else .n, .region = region };
    try run(allocator, io, T, vt, d, m, m, k, alpha, x, y, null, unlessOne(T, beta), c, ldc, want);
}

/// `C := alpha*S*B + beta*C` (`left`) or `alpha*B*S + beta*C`, `S` the
/// symmetric (`herm`: Hermitian) matrix whose `upper` or lower triangle `a`
/// holds; `symm_k.c` and `zhemm_k.c` over `level3.c`. `m`, `n` positive.
/// Up to `want` workers run it.
pub fn symm(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, comptime herm: bool, left: bool, upper: bool, m: usize, n: usize, alpha: T, a: [*]const T, lda: usize, b: [*]const T, ldb: usize, beta: T, c: [*]T, ldc: usize, want: u32) Error!void {
    if (isZero(T, alpha)) {
        if (!isOne(T, beta)) scale(T, vt, m, n, beta, c, ldc);
        return;
    }
    const s: Matrix(T) = .{ .p = @constCast(a), .ld = lda };
    const g: Matrix(T) = .{ .p = @constCast(b), .ld = ldb };
    const t: Pack = if (upper) .upper else .lower;
    // `hemm`'s right side packs `S^T`: the kernel conjugates it (`NC`)
    const d: Desc = if (left) .{ .rows = t, .cols = .n, .herm = herm } else .{ .rows = .n, .cols = t, .herm = herm, .conj = if (herm) .r else .n };
    try run(allocator, io, T, vt, d, m, n, if (left) m else n, alpha, if (left) s else g, if (left) g else s, null, unlessOne(T, beta), c, ldc, want);
}

/// `C := beta*C` on rows `[ib, ie)` and columns `[jb, je)` of `C`, the part
/// of them in `region`; with `real_diag` `beta` is real and the diagonal
/// comes out real (`zherk_k.c`'s `syrk_beta`). A zero `beta` never reads `C`.
pub fn scaleBlock(comptime T: type, vt: *const Vtable, region: Region, real_diag: bool, ib: usize, ie: usize, jb: usize, je: usize, beta: T, c: [*]T, ldc: usize) void {
    if (region == .full) return scale(T, vt, ie - ib, je - jb, beta, c + ib + jb * ldc, ldc);
    for (jb..je) |j| {
        const lo = if (region == .lower) @max(ib, j) else ib;
        const hi = if (region == .upper) @min(ie, j + 1) else ie;
        if (lo >= hi) continue;
        const col = c + lo + j * ldc;
        if (comptime isComplex(T)) {
            if (real_diag) {
                // a real `beta` scales each part, as the reference's real times complex
                const rows = hi - lo;
                _ = portable.kernel(Float(T), vt, .gemm_beta)(@intCast(2 * rows), 1, 0, beta.re, raw(T, col), 0, raw(T, col), 0, raw(T, col), @intCast(2 * rows));
                if (lo <= j and j < hi) c[j + j * ldc].im = 0;
                continue;
            }
        }
        scale(T, vt, hi - lo, 1, beta, col, ldc);
    }
}

/// `beta` as `run` takes it: null for 1.
pub fn unlessOne(comptime T: type, beta: T) ?T {
    return if (isOne(T, beta)) null else beta;
}

/// `C := alpha*op(A)*op(A)' + beta*C` on the `upper` or lower triangle of
/// the `n x n` `C`, `op(A)` being `n x k`; `'` is the transpose, or with
/// `herm` the conjugate transpose (`alpha`, `beta` then real, given as
/// complex with zero imaginary parts). `trans` takes `op(A) = A'`. With `b`,
/// `syr2k`: `alpha*op(A)*op(B)' + alpha2*op(B)*op(A)' + beta*C`, `alpha2`
/// being `alpha` (or its conjugate for `herm`). `level3_syrk.c` and
/// `level3_syr2k.c`, over `syrk_k.c`'s triangle. `n` positive. Up to `want`
/// workers run it.
pub fn syrk(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, comptime herm: bool, upper: bool, trans: bool, n: usize, k: usize, alpha: T, a: [*]const T, lda: usize, b: ?[*]const T, ldb: usize, beta: T, c: [*]T, ldc: usize, want: u32) Error!void {
    const region: Region = if (upper) .upper else .lower;
    if (k == 0 or isZero(T, alpha)) {
        if (!isOne(T, beta)) scaleBlock(T, vt, region, herm, 0, n, 0, n, beta, c, ldc);
        return;
    }
    const x: Matrix(T) = .{ .p = @constCast(a), .ld = lda };
    const y: Matrix(T) = if (b) |p| .{ .p = @constCast(p), .ld = ldb } else x;
    const alpha2: ?T = if (b == null) null else if (comptime herm) alpha.conjugate() else alpha;
    // `herk`'s `op(A)'` conjugates the columns (N) or `op(A)` the rows (C)
    try run(allocator, io, T, vt, .{
        .rows = if (trans) .t else .n,
        .cols = if (trans) .n else .t,
        .conj = if (herm) (if (trans) .l else .r) else .n,
        .region = region,
        .real_diag = herm,
    }, n, n, k, alpha, x, y, alpha2, unlessOne(T, beta), c, ldc, want);
}
