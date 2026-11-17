//! The public extension routines (issue #312): matrices as a slice and its
//! leading dimension, column-major, sizes `usize`, each slice checked to hold
//! its matrix (`util.holds`), over the reference-shaped core
//! (`interface.zig`). Pivots are 0-based row indices (user, 2026-11-15); the
//! core keeps the reference's 1-based ones, and `laswp` converts.

const std = @import("std");
const types = @import("../types.zig");
const strided = @import("../strided.zig");
const util = @import("../util.zig");
const core = @import("interface.zig");
const laswp_kernel = @import("laswp/kernel.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Trans = types.Trans;
const holds = util.holds;
const ConstStrided = strided.ConstStrided;

/// A letter the routine does not take: a conjugating one for `geadd`, `o` to
/// `d` for the matrix copies.
pub const TransError = error{UnsupportedTranspose};
/// `omatcopy`'s `A` and `B` sharing memory, which it cannot copy between.
pub const OverlapError = error{Overlap};
/// What `geadd` returns.
pub const GeaddError = TransError || strided.SliceError || strided.LeadingDimensionError;
/// What `omatcopy` returns.
pub const OmatcopyError = GeaddError || OverlapError;
/// What `imatcopy` returns: `geadd`'s, and the allocator's.
pub const ImatcopyError = GeaddError || Allocator.Error;
/// What `laswp` returns: a row or pivot outside the matrix, a slice too
/// short or a leading dimension too small.
pub const LaswpError = error{PivotOutOfRange} || strided.SliceError || strided.LeadingDimensionError;

fn len(n: usize) isize {
    return @intCast(n);
}

/// The core's errors past these checks.
fn mapped(e: core.Error) ImatcopyError {
    return switch (e) {
        error.OutOfMemory => error.OutOfMemory,
        error.InvalidTrans => error.UnsupportedTranspose,
        error.InvalidLda, error.InvalidLdb, error.InvalidLdc => error.LeadingDimensionTooSmall,
        error.Canceled, error.InvalidM, error.InvalidN, error.InvalidRows, error.InvalidCols, error.InvalidK1, error.InvalidIncx, error.InvalidIpiv => unreachable,
    };
}

/// The stored shape of a matrix `op` makes `rows` by `cols`.
fn stored(trans: Trans, rows: usize, cols: usize) [2]usize {
    return if (trans.transposed()) .{ cols, rows } else .{ rows, cols };
}

/// `C := alpha op(A) + beta C`, `op(A)` and `C` `m` by `n`; with `transc`
/// `t`, the same on `Cᵀ`, which `c` stores `n` by `m`. `op` transposes, never
/// conjugates. A zero `alpha` never reads `A`, a zero `beta` never reads `C`.
pub fn geadd(io: Io, comptime T: type, m: usize, n: usize, alpha: T, a: []const T, lda: usize, beta: T, c: []T, ldc: usize, transa: Trans, transc: Trans) GeaddError!void {
    const sa = stored(transa, m, n);
    const sc = stored(transc, m, n);
    try holds(a.len, sa[0], sa[1], lda);
    try holds(c.len, sc[0], sc[1], ldc);
    core.geadd(io, T, len(m), len(n), alpha, a.ptr, len(lda), beta, c.ptr, len(ldc), transa, transc) catch |e| return @errorCast(mapped(e));
}

/// Whether two slices share any memory.
fn overlaps(comptime T: type, a: []const T, b: []const T) bool {
    if (a.len == 0 or b.len == 0) return false;
    const a0 = @intFromPtr(a.ptr);
    const b0 = @intFromPtr(b.ptr);
    return a0 < b0 + b.len * @sizeOf(T) and b0 < a0 + a.len * @sizeOf(T);
}

/// `B := alpha op(A)`, `A` `rows` by `cols`, so `B` is `cols` by `rows`
/// where `op` transposes; `r` and `c` conjugate. A zero `alpha` stores zeros
/// without reading `A`, and `alpha` 1 copies it exactly. `A` and `B` must not
/// share memory; `imatcopy` works in place.
pub fn omatcopy(io: Io, comptime T: type, trans: Trans, rows: usize, cols: usize, alpha: T, a: []const T, lda: usize, b: []T, ldb: usize) OmatcopyError!void {
    const sb = stored(trans, rows, cols);
    try holds(a.len, rows, cols, lda);
    try holds(b.len, sb[0], sb[1], ldb);
    if (overlaps(T, a, b)) return error.Overlap;
    core.omatcopy(io, T, trans, len(rows), len(cols), alpha, a.ptr, len(lda), b.ptr, len(ldb)) catch |e| return @errorCast(mapped(e));
}

/// `A := alpha op(A)` in place: `A` read `rows` by `cols` at `lda`, written as
/// `op(A)` at `ldb`, so `a` must hold both. Through a buffer from `allocator`
/// where the shape or the leading dimension changes.
pub fn imatcopy(allocator: Allocator, io: Io, comptime T: type, trans: Trans, rows: usize, cols: usize, alpha: T, a: []T, lda: usize, ldb: usize) ImatcopyError!void {
    const sb = stored(trans, rows, cols);
    try holds(a.len, rows, cols, lda);
    try holds(a.len, sb[0], sb[1], ldb);
    core.imatcopy(allocator, io, T, trans, len(rows), len(cols), alpha, a.ptr, len(lda), len(ldb)) catch |e| return mapped(e);
}

/// Row interchanges, as `getrf` records them: for each pivot in `ipiv`'s
/// order, row `first + p` of the `m` by `n` matrix swaps with row `ipiv(p)`,
/// `p` being the pivot's place in the slice. Every index is 0-based. A
/// backward view (a negative stride) applies the same swaps last to first,
/// which undoes them. Every pivot is checked before the first swap.
pub fn laswp(io: Io, comptime T: type, m: usize, n: usize, a: []T, lda: usize, first: usize, ipiv: ConstStrided(usize)) LaswpError!void {
    _ = io;
    try holds(a.len, m, n, lda);
    const count = ipiv.n;
    if (count == 0 or n == 0) return;
    if (first > m or count > m - first) return error.PivotOutOfRange;
    for (0..count) |i| if (ipiv.get(i) >= m) return error.PivotOutOfRange;
    // the reference's 1-based pivots, a chunk at a time, through its kernel
    var buf: [256]isize = undefined;
    const backward = ipiv.stride < 0 and count > 1;
    const step: usize = @abs(ipiv.stride);
    var done: usize = 0;
    while (done < count) {
        const c = @min(buf.len, count - done);
        // the chunk's places in the slice, from the end when backward
        const p0 = if (backward) count - done - c else done;
        for (buf[0..c], 0..) |*b, j| b.* = @intCast(ipiv.data[(p0 + j) * step] + 1);
        const k1 = first + p0 + 1;
        laswp_kernel.swapRows(T, n, a.ptr, lda, k1, k1 + c - 1, &buf, if (backward) -1 else 1);
        done += c;
    }
}

// ---- tests -------------------------------------------------------------------

const check = @import("../_test/check.zig");
const testing = std.testing;

test "api: omatcopy transposes into a slice, and refuses overlap and short slices" {
    const a = [_]f64{ 1, 2, 3, 4, 5, 6 };
    var b: [6]f64 = @splat(0);
    try omatcopy(check.io, f64, .t, 2, 3, 1, &a, 2, &b, 3);
    try testing.expectEqualSlices(f64, &.{ 1, 3, 5, 2, 4, 6 }, &b);
    var both = [_]f64{ 1, 2, 3, 4, 5, 6, 0, 0 };
    try testing.expectError(error.Overlap, omatcopy(check.io, f64, .n, 2, 2, 1, both[0..4], 2, both[2..6], 2));
    try testing.expectError(error.SliceTooShort, omatcopy(check.io, f64, .n, 2, 3, 1, a[0..5], 2, &b, 2));
    try testing.expectError(error.LeadingDimensionTooSmall, omatcopy(check.io, f64, .t, 2, 3, 1, &a, 2, &b, 2));
}

test "api: geadd's conjugating letters are an error; imatcopy in place" {
    const a = [_]types.Complex64{ .init(1, 1), .init(2, 2) };
    var c = [_]types.Complex64{ .init(0, 0), .init(0, 0) };
    try testing.expectError(error.UnsupportedTranspose, geadd(check.io, types.Complex64, 2, 1, .init(1, 0), &a, 2, .init(0, 0), &c, 2, .r, .n));
    var x = [_]f64{ 1, 2, 3, 4 };
    try imatcopy(testing.allocator, check.io, f64, .t, 2, 2, 2, &x, 2, 2);
    try testing.expectEqualSlices(f64, &.{ 2, 6, 4, 8 }, &x);
}

test "api: laswp's 0-based pivots match the core's 1-based ones, and a backward view undoes them" {
    // 4 x 2, column-major: rows 0..3 are 10, 20, 30, 40 (and 11, 21, 31, 41)
    var a = [_]f64{ 10, 20, 30, 40, 11, 21, 31, 41 };
    var b = a;
    const piv = [_]usize{ 2, 3, 3 };
    try laswp(check.io, f64, 4, 2, &a, 4, 0, .of(&piv));
    const piv1 = [_]isize{ 3, 4, 4 };
    try core.laswp(check.io, f64, 2, &b, 4, 1, 3, &piv1, 1);
    try testing.expectEqualSlices(f64, &b, &a);
    try laswp(check.io, f64, 4, 2, &a, 4, 0, try .init(&piv, 3, -1));
    try testing.expectEqualSlices(f64, &.{ 10, 20, 30, 40, 11, 21, 31, 41 }, &a);
    try testing.expectError(error.PivotOutOfRange, laswp(check.io, f64, 4, 2, &a, 4, 0, .of(&[_]usize{4})));
    try testing.expectError(error.PivotOutOfRange, laswp(check.io, f64, 4, 2, &a, 4, 2, .of(&piv)));
}
