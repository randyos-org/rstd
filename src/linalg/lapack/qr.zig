//! QR without pivoting, from faer v0.24.4
//! (<https://codeberg.org/sarah-quinones/faer>, tag `faer-v0.24.4`, MIT):
//! `faer/src/linalg/qr/no_pivoting/factor.rs`'s staircase QR, behind
//! Reference-LAPACK's `?GEQRF`, `?ORGQR` and `?ORMQR` (`dgeqrf.f`,
//! `dorgqr.f`, `dormqr.f` and their `z` forms). `tau` is LAPACK's, `1 /
//! faer's`; departures: design.md's LAPACK section.

const std = @import("std");
const ops = @import("ops.zig");
const blas_calls = @import("blas_calls.zig");
const auxiliary = @import("auxiliary.zig");
const hh = @import("householder.zig");
const util = @import("../blas/util.zig");
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
/// taking the block factors; returns the rank faer estimates.
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

/// What `geqrf` returns.
pub const GeqrfError = error{ InvalidM, InvalidN, InvalidLda } || Allocator.Error || Io.Cancelable;
/// What `orgqr` returns: `InvalidN` unless `0 <= n <= m`, `InvalidK` unless
/// `0 <= k <= n`.
pub const OrgqrError = error{ InvalidM, InvalidN, InvalidK, InvalidLda } || Allocator.Error || Io.Cancelable;
/// What `ormqr` returns: `InvalidTrans` for other than `n`, `t` or `c`
/// (complex: `n` or `c`), `InvalidK` unless `0 <= k <=` the side's order.
pub const OrmqrError = error{ InvalidTrans, InvalidM, InvalidN, InvalidK, InvalidLda, InvalidLdc } || Allocator.Error || Io.Cancelable;

/// `?GEQRF`: `A = Q R` for the `m x n` `A`: `R` on and above the diagonal,
/// `Q = H(0) ... H(k-1)` (`k = min(m, n)`) as reflections below it and
/// their `tau`s in `tau` (`k`).
pub fn geqrf(allocator: Allocator, io: Io, comptime T: type, m: isize, n: isize, a: [*]T, lda: isize, tau: [*]T) GeqrfError!void {
    comptime ops.requireFloat(T, "geqrf");
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, m)) return error.InvalidLda;
    if (m == 0 or n == 0) return;
    if (!util.fitsMatrix(m, n, lda, @sizeOf(T))) return error.InvalidLda;
    const am: Mat(T) = .{ .ptr = a, .rows = @intCast(m), .cols = @intCast(n), .ld = @intCast(lda) };
    try factorLapack(allocator, io, T, am, tau);
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
pub fn orgqr(allocator: Allocator, io: Io, comptime T: type, m: isize, n: isize, k: isize, a: [*]T, lda: isize, tau: [*]const T) OrgqrError!void {
    comptime ops.requireFloat(T, "orgqr");
    if (m < 0) return error.InvalidM;
    if (n < 0 or n > m) return error.InvalidN;
    if (k < 0 or k > n) return error.InvalidK;
    if (lda < @max(1, m)) return error.InvalidLda;
    if (n == 0) return;
    if (!util.fitsMatrix(m, n, lda, @sizeOf(T))) return error.InvalidLda;
    const am: Mat(T) = .{ .ptr = a, .rows = @intCast(m), .cols = @intCast(n), .ld = @intCast(lda) };
    try generate(allocator, io, T, am, @intCast(k), tau);
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
pub fn ormqr(allocator: Allocator, io: Io, comptime T: type, side: Side, trans: Trans, m: isize, n: isize, k: isize, a: [*]const T, lda: isize, tau: [*]const T, c: [*]T, ldc: isize) OrmqrError!void {
    comptime ops.requireFloat(T, "ormqr");
    switch (trans) {
        .n, .c => {},
        .t => if (comptime ops.isComplex(T)) return error.InvalidTrans,
        .r, .o, .u, .s, .d => return error.InvalidTrans,
    }
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    const nq = if (side == .l) m else n;
    if (k < 0 or k > nq) return error.InvalidK;
    if (lda < @max(1, nq)) return error.InvalidLda;
    if (ldc < @max(1, m)) return error.InvalidLdc;
    if (m == 0 or n == 0 or k == 0) return;
    if (!util.fitsMatrix(nq, k, lda, @sizeOf(T))) return error.InvalidLda;
    if (!util.fitsMatrix(m, n, ldc, @sizeOf(T))) return error.InvalidLdc;
    const v: Mat(T) = .{ .ptr = @constCast(a), .rows = @intCast(nq), .cols = @intCast(k), .ld = @intCast(lda) };
    const cm: Mat(T) = .{ .ptr = c, .rows = @intCast(m), .cols = @intCast(n), .ld = @intCast(ldc) };
    try multiply(allocator, io, T, side, trans != .n, v, tau, cm);
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

// ---- tests -----------------------------------------------------------------------

const testing = std.testing;
const check = @import("../blas/_test/check.zig");
const lcheck = @import("_test/check.zig");
const blas = @import("../blas/__root__.zig");
const test_io = @import("../blas/_test/io.zig");

/// `Q R` against `A`, within `n eps ‖A‖` per element (faer's `test_qr` takes
/// `1e-10` at `c64`).
fn expectQr(comptime T: type, al: Allocator, a0: Mat(T), qr: Mat(T), h: Mat(T)) !void {
    const F = Float(T);
    const m = a0.rows;
    const n = a0.cols;
    const q = try ops.alloc(T, al, m, m);
    defer ops.free(T, al, q);
    q.identity();
    try hh.applySequence(al, check.io, T, qr, h, false, q);
    var amax: F = 0;
    for (0..n) |j| for (0..m) |i| {
        amax = @max(amax, ops.abs(T, a0.get(i, j)));
    };
    const tol = @as(F, @floatFromInt(4 * (m + n))) * std.math.floatEps(F) * @max(amax, 1);
    for (0..n) |j| for (0..m) |i| {
        var s = ops.zero(T);
        for (0..@min(j + 1, m)) |k| s = ops.add(T, s, ops.mul(T, q.get(i, k), qr.get(k, j)));
        const err = ops.abs(T, ops.sub(T, s, a0.get(i, j)));
        if (!(err <= tol)) {
            check.say("qr {s} {d}x{d}: ({d},{d}) off by {e}, tol {e}\n", .{ @typeName(T), m, n, i, j, err, tol });
            return error.TestUnexpectedResult;
        }
    };
}

test "faer qr: Q R = A, square, tall, wide, blocked and rank deficient" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex64, f32, types.Complex32 }) |T| {
        for ([_][3]usize{ .{ 1, 1, 0 }, .{ 2, 2, 0 }, .{ 8, 8, 0 }, .{ 40, 17, 0 }, .{ 17, 40, 0 }, .{ 70, 70, 0 }, .{ 130, 60, 0 }, .{ 64, 64, 3 }, .{ 80, 50, 1 } }) |c| {
            const m = c[0];
            const n = c[1];
            const a0 = try ops.alloc(T, al, m, n);
            defer ops.free(T, al, a0);
            if (c[2] == 0) {
                check.fill(T, a0.ptr[0 .. m * n], m * 31 + n);
            } else try lcheck.rankProduct(T, al, a0, c[2], 3);
            const qr = try ops.alloc(T, al, m, n);
            defer ops.free(T, al, qr);
            qr.copyFrom(a0);
            const h = try ops.alloc(T, al, blockSize(m, n), @min(m, n));
            defer ops.free(T, al, h);
            const rank = try factor(al, check.io, T, qr, h);
            try testing.expect(rank >= @min(c[2], @min(m, n)));
            try expectQr(T, al, a0, qr, h);
        }
    }
}

const test_types = .{ f32, f64, types.Complex32, types.Complex64 };

/// `expectQr`'s tolerance, `4 (m + n) eps max(1, ‖A‖max)` per element.
fn qrTol(comptime F: type, m: usize, n: usize, amax: F) F {
    return @as(F, @floatFromInt(4 * (m + n))) * std.math.floatEps(F) * @max(amax, 1);
}

test "faer qr: geqrf and orgqr, Q R = A and Qᴴ Q = I in LAPACK's layout" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        // m, n, rank (0 full, else a product; n + 1 a zero column)
        for ([_][3]usize{ .{ 1, 1, 0 }, .{ 5, 3, 0 }, .{ 3, 5, 0 }, .{ 40, 17, 0 }, .{ 17, 40, 0 }, .{ 70, 70, 0 }, .{ 130, 60, 0 }, .{ 64, 64, 3 }, .{ 30, 20, 21 } }) |c| {
            const m = c[0];
            const n = c[1];
            const k = @min(m, n);
            const a0 = try ops.alloc(T, al, m, n);
            defer ops.free(T, al, a0);
            if (c[2] == n + 1) {
                check.fill(T, a0.ptr[0 .. m * n], 9);
                @memset(a0.col(4)[0..m], ops.zero(T));
            } else if (c[2] == 0) {
                check.fill(T, a0.ptr[0 .. m * n], m + 7 * n);
            } else try lcheck.rankProduct(T, al, a0, c[2], m + 7 * n);
            const lda = m + 1;
            const a = try ops.alloc(T, al, lda, n);
            defer ops.free(T, al, a);
            a.fill(lcheck.poisoned(T));
            a.sub(0, 0, m, n).copyFrom(a0);
            const tau = try al.alloc(T, k);
            defer al.free(tau);
            try geqrf(al, check.io, T, @intCast(m), @intCast(n), a.ptr, @intCast(lda), tau.ptr);
            for (0..n) |j| try testing.expect(lcheck.isNanT(T, a.get(m, j)));
            // R, then Q's first k columns
            const r = try ops.alloc(T, al, k, n);
            defer ops.free(T, al, r);
            for (0..n) |j| for (0..k) |i| {
                r.at(i, j).* = if (i <= j) a.get(i, j) else ops.zero(T);
            };
            try orgqr(al, check.io, T, @intCast(m), @intCast(k), @intCast(k), a.ptr, @intCast(lda), tau.ptr);
            const q = a.sub(0, 0, m, k);
            const qr_ = try ops.alloc(T, al, m, n);
            defer ops.free(T, al, qr_);
            try blas_calls.gemm(al, check.io, T, .n, .n, ops.one(T), q, r, ops.zero(T), qr_);
            try lcheck.expectNear(T, "geqrf Q R", qr_, a0, qrTol(F, m, n, auxiliary.normMax(T, a0)));
            const g = try ops.alloc(T, al, k, k);
            defer ops.free(T, al, g);
            try blas_calls.gemm(al, check.io, T, .c, .n, ops.one(T), q, q, ops.zero(T), g);
            const id = try ops.alloc(T, al, k, k);
            defer ops.free(T, al, id);
            id.identity();
            try lcheck.expectNear(T, "orgqr Qh Q", g, id, qrTol(F, m, n, 1));
            if (c[2] == n + 1) try testing.expect(ops.isZero(T, r.get(4, 4)));
        }
    }
}

/// `C := op(H(0) ... H(k-1)) C` or `C op(...)`, `H = I - tau v vᴴ`, one
/// reflection at a time, `v` unit lower in `v`'s columns.
fn refApply(comptime T: type, side: Side, adjoint: bool, v: Mat(T), tau: []const T, c: Mat(T)) void {
    const k = tau.len;
    const left = side == .l;
    // left: Q C applies the last first; Qᴴ C the first first; right mirrors
    for (0..k) |s| {
        const j = if (left != adjoint) k - 1 - s else s;
        const t = if (adjoint) ops.conj(T, tau[j]) else tau[j];
        const nq = v.rows;
        const others = if (left) c.cols else c.rows;
        for (0..others) |o| {
            // d := vᴴ x (left, x a column) or x v (right, x a row)
            var d = ops.zero(T);
            for (j..nq) |i| {
                const vi = if (i == j) ops.one(T) else v.get(i, j);
                d = if (left) ops.add(T, d, ops.mul(T, ops.conj(T, vi), c.get(i, o))) else ops.add(T, d, ops.mul(T, c.get(o, i), vi));
            }
            d = ops.mul(T, t, d);
            for (j..nq) |i| {
                const vi = if (i == j) ops.one(T) else v.get(i, j);
                if (left) {
                    c.at(i, o).* = ops.sub(T, c.get(i, o), ops.mul(T, vi, d));
                } else {
                    c.at(o, i).* = ops.sub(T, c.get(o, i), ops.mul(T, d, ops.conj(T, vi)));
                }
            }
        }
    }
}

test "faer qr: ormqr and orgqr against one reflection at a time, faer's tau and the reference's" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        const nq = 37;
        const k = 20;
        const other = 9;
        for (0..2) |kind| {
            const v = try ops.alloc(T, al, nq, nq);
            defer ops.free(T, al, v);
            const tau = try al.alloc(T, k);
            defer al.free(tau);
            if (kind == 0) {
                // geqrf's
                check.fill(T, v.ptr[0 .. nq * nq], 21);
                try geqrf(al, check.io, T, nq, k, v.ptr, nq, tau.ptr);
            } else {
                // the reference's form: `tau = (1 - e^(i phi)) / ‖v‖²`
                // (`xLARFG`'s `(beta - alpha) / beta`), real 2 / ‖v‖²; one 0
                check.fill(T, v.ptr[0 .. nq * nq], 22);
                var prng = std.Random.DefaultPrng.init(23);
                for (tau, 0..) |*t, j| {
                    var nv: F = 1;
                    for (j + 1..nq) |i| nv += ops.abs2(T, v.get(i, j));
                    const phi = prng.random().float(F) * 2 * std.math.pi;
                    t.* = if (comptime ops.isComplex(T)) .init((1 - @cos(phi)) / nv, -@sin(phi) / nv) else 2 / nv;
                }
                tau[5] = ops.zero(T);
            }
            for ([_]Side{ .l, .r }) |side| for ([_]bool{ false, true }) |adj| {
                const rows: usize = if (side == .l) nq else other;
                const cols: usize = if (side == .l) other else nq;
                const c = try ops.alloc(T, al, rows, cols);
                defer ops.free(T, al, c);
                check.fill(T, c.ptr[0 .. rows * cols], 24);
                const want = try ops.alloc(T, al, rows, cols);
                defer ops.free(T, al, want);
                want.copyFrom(c);
                refApply(T, side, adj, v.sub(0, 0, nq, k), tau, want);
                const tr: Trans = if (adj) .c else .n;
                try ormqr(al, check.io, T, side, tr, @intCast(rows), @intCast(cols), k, v.ptr, nq, tau.ptr, c.ptr, @intCast(rows));
                try lcheck.expectNear(T, "ormqr", c, want, qrTol(F, nq, nq, auxiliary.normMax(T, want)));
            };
            // orgqr: Q's first n columns
            const n = 30;
            const q = try ops.alloc(T, al, nq, nq);
            defer ops.free(T, al, q);
            q.copyFrom(v);
            try orgqr(al, check.io, T, nq, n, k, q.ptr, nq, tau.ptr);
            const want = try ops.alloc(T, al, nq, n);
            defer ops.free(T, al, want);
            want.identity();
            refApply(T, .l, false, v.sub(0, 0, nq, k), tau, want);
            try lcheck.expectNear(T, "orgqr", q.sub(0, 0, nq, n), want, qrTol(F, nq, nq, 1));
        }
    }
}

test "faer qr: geqrf, orgqr and ormqr report the first bad argument and read nothing when empty" {
    const fa = testing.failing_allocator;
    const nan = std.math.nan(f64);
    var a: [4]f64 = @splat(nan);
    var tau: [2]f64 = @splat(nan);
    var c: [4]f64 = @splat(nan);
    try testing.expectError(error.InvalidM, geqrf(fa, check.io, f64, -1, -1, &a, 0, &tau));
    try testing.expectError(error.InvalidN, geqrf(fa, check.io, f64, 2, -1, &a, 0, &tau));
    try testing.expectError(error.InvalidLda, geqrf(fa, check.io, f64, 2, 2, &a, 1, &tau));
    try testing.expectError(error.InvalidM, orgqr(fa, check.io, f64, -1, -1, -1, &a, 0, &tau));
    try testing.expectError(error.InvalidN, orgqr(fa, check.io, f64, 1, 2, -1, &a, 0, &tau));
    try testing.expectError(error.InvalidK, orgqr(fa, check.io, f64, 2, 1, 2, &a, 0, &tau));
    try testing.expectError(error.InvalidLda, orgqr(fa, check.io, f64, 2, 2, 2, &a, 1, &tau));
    try testing.expectError(error.InvalidTrans, ormqr(fa, check.io, f64, .l, .r, -1, -1, -1, &a, 0, &tau, &c, 0));
    try testing.expectError(error.InvalidTrans, ormqr(fa, check.io, types.Complex64, .l, .t, 1, 1, 1, @ptrCast(&a), 1, @ptrCast(&tau), @ptrCast(&c), 1));
    try testing.expectError(error.InvalidM, ormqr(fa, check.io, f64, .l, .t, -1, -1, -1, &a, 0, &tau, &c, 0));
    try testing.expectError(error.InvalidN, ormqr(fa, check.io, f64, .l, .t, 2, -1, -1, &a, 0, &tau, &c, 0));
    try testing.expectError(error.InvalidK, ormqr(fa, check.io, f64, .r, .n, 3, 2, 3, &a, 0, &tau, &c, 0));
    try testing.expectError(error.InvalidLda, ormqr(fa, check.io, f64, .l, .n, 2, 1, 1, &a, 1, &tau, &c, 0));
    try testing.expectError(error.InvalidLdc, ormqr(fa, check.io, f64, .r, .n, 2, 1, 1, &a, 1, &tau, &c, 1));
    // empty: nothing read or written
    try geqrf(fa, check.io, f64, 0, 2, &a, 1, &tau);
    try orgqr(fa, check.io, f64, 2, 0, 0, &a, 2, &tau);
    try ormqr(fa, check.io, f64, .l, .n, 2, 2, 0, &a, 2, &tau, &c, 2);
    for ([_][]const f64{ &a, &tau, &c }) |x| for (x) |e| try testing.expect(std.math.isNan(e));
}

test "faer qr: geqrf and ormqr give one worker's bits at any worker count" {
    if (comptime !test_io.hosted) return error.SkipZigTest;
    _ = blas.coreName();
    const al = testing.allocator;
    var runs: test_io.WorkerRuns = .init(al);
    defer runs.deinit();
    inline for (.{ f64, types.Complex32 }) |T| {
        const m = 300;
        const n = 140;
        const a0 = try ops.alloc(T, al, m, n);
        defer ops.free(T, al, a0);
        check.fill(T, a0.ptr[0 .. m * n], 77);
        var outs: [3][]T = undefined;
        for (&outs) |*o| o.* = try al.alloc(T, m * n + n + m * 40);
        defer for (outs) |o| al.free(o);
        for (0..test_io.worker_runs.len) |ci| {
            const io = runs.io(ci);
            const o = outs[ci];
            const a: Mat(T) = .{ .ptr = o.ptr, .rows = m, .cols = n, .ld = m };
            a.copyFrom(a0);
            try geqrf(al, io, T, m, n, a.ptr, m, o.ptr + m * n);
            const c: Mat(T) = .{ .ptr = o.ptr + m * n + n, .rows = m, .cols = 40, .ld = m };
            check.fill(T, c.ptr[0 .. m * 40], 78);
            try ormqr(al, io, T, .l, .c, m, 40, n, a.ptr, m, o.ptr + m * n, c.ptr, m);
        }
        try test_io.expectSameBytes(outs);
    }
}
