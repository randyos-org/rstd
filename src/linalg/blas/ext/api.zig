//! The public extension routines (issue #312): matrices as a slice and its
//! leading dimension, column-major, sizes `usize`, each slice checked to hold
//! its matrix (`util.holds`), over the reference-shaped core
//! (`interface.zig`). Pivots are 0-based row indices (user, 2026-11-15, issue
//! #398), inside too: `laswp` hands its view to the row-swap kernel as it is.

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
        error.Canceled, error.InvalidM, error.InvalidN, error.InvalidRows, error.InvalidCols => unreachable,
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
    const backward = ipiv.stride < 0 and count > 1;
    const step: usize = if (count > 1) @abs(ipiv.stride) else 1;
    if (step == 1) return laswp_kernel.swapRows(T, n, a.ptr, lda, first, count, ipiv.data.ptr, backward);
    // pivots `step` apart: one swap at a time, in the view's order
    for (0..count) |s| {
        const p = if (backward) count - 1 - s else s;
        laswp_kernel.swapRows(T, n, a.ptr, lda, first + p, 1, ipiv.data.ptr + p * step, false);
    }
}

// ---- tests -------------------------------------------------------------------

const check = @import("../_test/check.zig");
const testing = std.testing;
const ref = @import("laswp/__root__.zig");
const bitwise = @import("../_test/bitwise.zig");

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

test "api: laswp follows dlaswp.f on 0-based pivots, any stride, and checks every pivot first" {
    const al = testing.allocator;
    inline for (.{ f32, types.Complex64, f128 }) |T| {
        for ([_]isize{ 1, -1, 2, -2 }) |incx| {
            const m = 9;
            const n = 4;
            const lda = m + 2;
            var ipiv: [2 * m]isize = undefined;
            ref.pivots(&ipiv, m, @abs(incx));
            const a = try al.alloc(T, lda * n);
            defer al.free(a);
            check.fill(T, a, 4);
            const want = try al.dupe(T, a);
            defer al.free(want);
            // rows 2 to 8, 1-based, as `dlaswp.f` takes them
            ref.laswpRef(T, n, want, lda, 2, 8, &ipiv, incx);
            var zero: [2 * m]usize = undefined;
            for (&zero, ipiv) |*z, p| z.* = @intCast(p - 1);
            try laswp(check.io, T, m, n, a, lda, 1, try .init(zero[1..], 7, incx));
            for (want, a) |w, g| try testing.expect(bitwise.same(T, w, g, true));
        }
    }
    var a: [12]f64 = @splat(1);
    var piv = [_]usize{ 0, 1, 2, 3 };
    try laswp(check.io, f64, 4, 3, &a, 4, 0, try .init(&piv, 0, 1));
    try testing.expectError(error.LeadingDimensionTooSmall, laswp(check.io, f64, 4, 3, &a, 3, 0, .of(&piv)));
    try testing.expectError(error.PivotOutOfRange, laswp(check.io, f64, 4, 3, &a, 4, 1, .of(&piv)));
    piv[2] = 4;
    try testing.expectError(error.PivotOutOfRange, laswp(check.io, f64, 4, 3, &a, 4, 0, .of(&piv)));
    try testing.expectError(error.PivotOutOfRange, laswp(check.io, f64, 4, 3, &a, 4, 0, try .init(&piv, 4, -1)));
    for (a) |v| try testing.expectEqual(1, v);
}
