//! QR without pivoting, from faer v0.24.4
//! (<https://codeberg.org/sarah-quinones/faer>, tag `faer-v0.24.4`, MIT):
//! `faer/src/linalg/qr/no_pivoting/factor.rs`'s staircase QR, behind
//! Reference-LAPACK's `?GEQRF`, `?ORGQR` and `?ORMQR` (`dgeqrf.f`,
//! `dorgqr.f`, `dormqr.f` and their `z` forms). `tau` is LAPACK's, `1 /
//! faer's`; departures: design.md's LAPACK section.

const std = @import("std");
const ops = @import("ops.zig");
const blas_calls = @import("blas_calls.zig");
const hh = @import("householder.zig");
const util = @import("../blas/util.zig");
const interface = @import("interface.zig");
const types = @import("../blas/types.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = ops.Mat;
const Float = ops.Float;
const Error = ops.Error;
const Trans = types.Trans;
const Side = types.Side;

/// `QrParams::auto`'s `blocking_threshold`: below `48 * 48` elements a
/// block's halves are single reflections.
const blocking_threshold = 48 * 48;

/// `recommended_block_size`: the block size for an `m x n` factorization.
pub fn blockSize(m: usize, n: usize) usize {
    const prod = m *| n;
    const bs: usize = if (prod > 8192 * 8192) 256 else if (prod > 2048 * 2048) 128 else if (prod > 1024 * 1024) 64 else if (prod > 512 * 512) 48 else if (prod > 128 * 128) 32 else if (prod > 32 * 32) 8 else if (prod > 16 * 16) 4 else 1;
    return @max(@min(bs, @min(m, n)), 1);
}

/// `qr_in_place`: factors `a` in place, `h` (`bs x min(m, n)`, `bs > 0`)
/// taking the block factors; returns the rank faer estimates. Asserts `h`'s
/// shape (internal callers size it).
pub fn factor(al: Allocator, io: Io, comptime T: type, a: Mat(T), h: Mat(T)) Error!usize {
    std.debug.assert(h.rows > 0 and h.cols == @min(a.rows, a.cols));
    const rank = try blocked(al, io, T, a, h, 0, 0);
    const bs = h.rows;
    const n = h.cols;
    for (rank..n) |j| @memset(h.col(j)[0..bs], ops.zero(T));
    var c = rank / bs * bs;
    while (c < n) : (c += bs) {
        const b = @min(bs, n - c);
        for (@max(rank, c)..c + b) |j| h.at(j - c, j).* = ops.scalar(T, std.math.inf(Float(T)));
    }
    return rank;
}

/// `qr_in_place_unblocked`: reflections from column `col` and row `row` on,
/// `tau`s into row `h` (stride `ldh`), each applied only to `a`'s columns.
fn unblocked(al: Allocator, io: Io, comptime T: type, a: Mat(T), h: [*]T, ldh: usize, size: usize, row_start: usize, col_start: usize) Error!usize {
    const F = Float(T);
    const m = a.rows;
    const n = a.cols;
    var row = row_start;
    var c = col_start;
    while (row < @min(size, m) and c < n) : (c += 1) {
        const above = blas_calls.nrm2(io, T, row, a.col(c), 1);
        const tail = m - row - 1;
        const info = if (row == c)
            hh.make(io, T, a.at(row, c), tail, a.col(c) + row + 1, 1, null, 0)
        else blk: {
            const info = hh.make(io, T, a.at(row, c), tail, a.col(row) + row + 1, 1, a.col(c) + row + 1, 1);
            @memset(a.col(c)[row + 1 ..][0..@min(tail, c - row)], ops.zero(T));
            break :blk info;
        };
        const v = a.col(row) + row + 1;
        const norm = std.math.hypot(info.norm, above);
        const threshold = std.math.floatEps(F) * (@as(F, @floatFromInt(m - row)) * 16) * norm;
        const tau_inv = 1 / info.tau;
        h[row * ldh] = ops.scalar(T, info.tau);
        if (tau_inv < std.math.floatMin(F)) {
            if (info.norm > 0) row += 1;
        } else if (!(info.norm <= threshold)) {
            // faer tests `norm > threshold`, false for NaN, so a NaN column
            // took no reflection and read as rank deficient (issue #256)
            // w := A01 + vᴴ A11, k := -w / tau; A01 += k, A11 += v kᵀ
            const rest = n - c - 1;
            if (rest > 0) {
                const w = try al.alloc(T, rest);
                defer al.free(w);
                for (w, 0..) |*x, j| x.* = a.get(row, c + 1 + j);
                try blas_calls.gemv(al, io, T, .u, ops.one(T), a.sub(row + 1, c + 1, tail, rest), v, 1, ops.one(T), w.ptr, 1);
                for (w, 0..) |*x, j| {
                    x.* = ops.neg(T, ops.mulReal(T, x.*, tau_inv));
                    const p = a.at(row, c + 1 + j);
                    p.* = ops.add(T, p.*, x.*);
                }
                try blas_calls.ger(al, io, T, ops.one(T), v, 1, w.ptr, 1, a.sub(row + 1, c + 1, tail, rest));
            }
            row += 1;
        }
    }
    return row;
}

/// `qr_in_place_blocked`: blocks of `h.rows` reflections (not columns: a
/// block takes columns until it has them all), each from halves recursively;
/// returns the next row.
fn blocked(al: Allocator, io: Io, comptime T: type, a: Mat(T), h: Mat(T), row_start: usize, col_start: usize) Error!usize {
    const m = a.rows;
    const n = a.cols;
    const size = @min(m, n);
    const bs = h.rows;
    if (bs == 1) return unblocked(al, io, T, a, h.ptr, h.ld, h.cols, row_start, col_start);
    const sub = if (m * n < blocking_threshold) 1 else bs / 2;
    var row = row_start;
    var c = col_start;
    while (row < size and c < n) {
        const block = @min(bs, @min(size - row, n - c));
        const sb = @min(block, sub);
        const start = row;
        var offset: usize = 0;
        while (offset < block and c < n) {
            const width = @min(n - c, block - offset);
            const sw = @min(width, sb);
            const new_row = try blocked(al, io, T, a.sub(0, 0, m, c + width), h.sub(offset, 0, sw, h.cols), row, c);
            const local = new_row - row;
            if (local > 0) {
                // the halves' factors onto the block's diagonal
                var k: usize = sw;
                while (k < local) : (k += sw) {
                    const s = @min(sw, local - k);
                    for (0..s) |j| for (0..j + 1) |i| {
                        h.at(offset + k + i, row + k + j).* = h.get(offset + i, row + k + j);
                    };
                }
                if (local > sw) try hh.blockFactor(al, io, T, a.sub(row, row, m - row, local), h.sub(offset, row, local, local));
                if (offset > 0) try hh.blockFactor(al, io, T, a.sub(start, start, m - start, offset + local), h.sub(0, start, offset + local, offset + local));
                if (c + width < n) try hh.applyBlock(al, io, T, a.sub(row, row, m - row, local), h.sub(offset, row, local, local), false, true, a.sub(row, c + width, m - row, n - c - width));
            }
            offset += local;
            row += local;
            c += width;
        }
    }
    return row;
}

// ---- LAPACK ----------------------------------------------------------------------

/// What `geqrf` returns: a slice too short (`a`, or `tau` under
/// `min(m, n)`) or a leading dimension too small, and the allocator's and
/// cancellation.
pub const GeqrfError = interface.ShapeError || Allocator.Error || Io.Cancelable;
/// What `orgqr` returns: `TooManyColumns` for `n > m`, `TooManyReflections`
/// for `k > n`, and `GeqrfError`'s.
pub const OrgqrError = error{TooManyColumns} || interface.ReflectionsError || GeqrfError;
/// What `ormqr` returns: a `trans` other than `n`, `t` or `c` (complex: `n`
/// or `c`), `TooManyReflections` for `k` past the side's order, and
/// `GeqrfError`'s.
pub const OrmqrError = interface.TransError || interface.ReflectionsError || GeqrfError;

/// `?GEQRF`: `A = Q R` for the `m x n` `A`: `R` on and above the diagonal,
/// `Q = H(0) ... H(k-1)` (`k = min(m, n)`) as reflections below it and
/// their `tau`s in `tau` (`k`).
pub fn geqrf(allocator: Allocator, io: Io, comptime T: type, m: usize, n: usize, a: []T, lda: usize, tau: []T) GeqrfError!void {
    comptime ops.requireFloat(T, "geqrf");
    try util.holds(a.len, m, n, lda);
    if (tau.len < @min(m, n)) return error.SliceTooShort;
    if (m == 0 or n == 0) return;
    const am: Mat(T) = .{ .ptr = a.ptr, .rows = m, .cols = n, .ld = lda };
    try factorLapack(allocator, io, T, am, tau.ptr);
}

/// `geqrf` after its checks.
pub fn factorLapack(al: Allocator, io: Io, comptime T: type, a: Mat(T), tau: [*]T) Error!void {
    const k = @min(a.rows, a.cols);
    if (k == 0) return;
    const bs = blockSize(a.rows, a.cols);
    const h = try ops.alloc(T, al, bs, k);
    defer ops.free(T, al, h);
    try factorZeroed(al, io, T, a, h);
    for (0..k) |j| tau[j] = ops.scalar(T, 1 / ops.re(T, h.get(j % bs, j)));
}

/// `factor`, with the identity reflections' `v` zeroed, as `xLARFG` leaves
/// a zero column.
pub fn factorZeroed(al: Allocator, io: Io, comptime T: type, a: Mat(T), h: Mat(T)) Error!void {
    const rank = try factor(al, io, T, a, h);
    for (rank..h.cols) |j| @memset(a.col(j)[j + 1 .. a.rows], ops.zero(T));
}

/// faer's `T` diagonal entry for LAPACK's `tau`: `1 / tau`, infinite for 0.
fn factorOf(comptime T: type, tau: T) T {
    return if (ops.isZero(T, tau)) ops.scalar(T, std.math.inf(Float(T))) else ops.recip(T, tau);
}

/// The block factors (`bs x k`, allocated) of the `k` reflections in `v`'s
/// columns with LAPACK `tau`s.
pub fn factorsOf(al: Allocator, io: Io, comptime T: type, v: Mat(T), tau: [*]const T) Error!Mat(T) {
    const k = v.cols;
    const bs = blockSize(v.rows, k);
    const h = try ops.alloc(T, al, bs, k);
    errdefer ops.free(T, al, h);
    for (0..k) |j| h.at(j % bs, j).* = factorOf(T, tau[j]);
    try hh.sequenceFactor(al, io, T, v, h);
    return h;
}

/// `?ORGQR`, `?UNGQR`: the `m x n` `Q` with orthonormal columns, the first
/// `n` of `H(0) ... H(k-1)`, over the reflections `geqrf` left in `A`'s
/// first `k` columns.
pub fn orgqr(allocator: Allocator, io: Io, comptime T: type, m: usize, n: usize, k: usize, a: []T, lda: usize, tau: []const T) OrgqrError!void {
    comptime ops.requireFloat(T, "orgqr");
    if (n > m) return error.TooManyColumns;
    if (k > n) return error.TooManyReflections;
    try util.holds(a.len, m, n, lda);
    if (tau.len < k) return error.SliceTooShort;
    if (n == 0) return;
    const am: Mat(T) = .{ .ptr = a.ptr, .rows = m, .cols = n, .ld = lda };
    try generate(allocator, io, T, am, k, tau.ptr);
}

/// `orgqr` after its checks: `a`'s first `k` columns hold the reflections.
pub fn generate(al: Allocator, io: Io, comptime T: type, a: Mat(T), k: usize, tau: [*]const T) Error!void {
    if (k == 0) return a.identity();
    const v = try ops.alloc(T, al, a.rows, k);
    defer ops.free(T, al, v);
    v.copyFrom(a.sub(0, 0, a.rows, k));
    const h = try factorsOf(al, io, T, v, tau);
    defer ops.free(T, al, h);
    a.identity();
    try hh.applySequence(al, io, T, v, h, false, a);
}

/// `?ORMQR`, `?UNMQR`: `C := op(Q) C` (side `l`) or `C op(Q)` (side `r`),
/// `op` by `trans` (`t` and `c` both `Qᴴ`), for the `m x n` `C` and the
/// `Q` of `k` reflections in `A`'s first `k` columns (`m` rows for `l`, `n`
/// for `r`). Only `A`'s strictly lower part is read.
pub fn ormqr(allocator: Allocator, io: Io, comptime T: type, side: Side, trans: Trans, m: usize, n: usize, k: usize, a: []const T, lda: usize, tau: []const T, c: []T, ldc: usize) OrmqrError!void {
    comptime ops.requireFloat(T, "ormqr");
    try interface.takesAdjoint(T, trans);
    const nq = if (side == .l) m else n;
    if (k > nq) return error.TooManyReflections;
    try util.holds(a.len, nq, k, lda);
    try util.holds(c.len, m, n, ldc);
    if (tau.len < k) return error.SliceTooShort;
    if (m == 0 or n == 0 or k == 0) return;
    const v: Mat(T) = .{ .ptr = @constCast(a.ptr), .rows = nq, .cols = k, .ld = lda };
    const cm: Mat(T) = .{ .ptr = c.ptr, .rows = m, .cols = n, .ld = ldc };
    try multiply(allocator, io, T, side, trans != .n, v, tau.ptr, cm);
}

/// `ormqr` after its checks: `adjoint` takes `Qᴴ`.
pub fn multiply(al: Allocator, io: Io, comptime T: type, side: Side, adjoint: bool, v: Mat(T), tau: [*]const T, c: Mat(T)) Error!void {
    const h = try factorsOf(al, io, T, v, tau);
    defer ops.free(T, al, h);
    switch (side) {
        .l => if (adjoint) try hh.applySequenceTransposed(al, io, T, v, h, true, c) else try hh.applySequence(al, io, T, v, h, false, c),
        .r => try hh.applySequenceRight(al, io, T, v, h, false, adjoint, c),
    }
}
