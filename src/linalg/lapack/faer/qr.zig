//! QR without pivoting, from faer v0.24.4
//! (<https://codeberg.org/sarah-quinones/faer>, tag `faer-v0.24.4`, MIT):
//! `faer/src/linalg/qr/no_pivoting/factor.rs`. `A = Q R` in place: `R` on and
//! above the diagonal, `Q`'s reflections below it (`householder.zig`), their
//! block factors in a `bs x min(m, n)` matrix.
//!
//! Faer's QR is a staircase: a column whose remaining norm is within
//! `16 (m - row) eps` of its whole norm takes no reflection of its own, so
//! the rank it returns can fall short of `min(m, n)`; the reflections past it
//! are the identity (`tau` infinite). Blocks are recursive halves down to
//! single reflections, each block's trailing columns updated by `gemm`
//! (`householder.applyBlock`). The SVD factors tall matrices with it; T111's
//! `geqrf` will too.

const std = @import("std");
const ops = @import("ops.zig");
const hh = @import("householder.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = ops.Mat;
const Float = ops.Float;
const Error = ops.Error;

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
        const above = ops.nrm2(io, T, row, a.col(c), 1);
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
        } else if (info.norm > threshold) {
            // w := A01 + vᴴ A11, k := -w / tau; A01 += k, A11 += v kᵀ
            const rest = n - c - 1;
            if (rest > 0) {
                const w = try al.alloc(T, rest);
                defer al.free(w);
                for (w, 0..) |*x, j| x.* = a.get(row, c + 1 + j);
                try ops.gemv(al, io, T, .u, ops.one(T), a.sub(row + 1, c + 1, tail, rest), v, 1, ops.one(T), w.ptr, 1);
                for (w, 0..) |*x, j| {
                    x.* = ops.neg(T, ops.mulReal(T, x.*, tau_inv));
                    const p = a.at(row, c + 1 + j);
                    p.* = ops.add(T, p.*, x.*);
                }
                try ops.ger(al, io, T, ops.one(T), v, 1, w.ptr, 1, a.sub(row + 1, c + 1, tail, rest));
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

// ---- tests -----------------------------------------------------------------------

const testing = std.testing;
const check = @import("../../blas/_test/blas1/check.zig");
const types = @import("../../blas/types.zig");

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
            } else {
                // rank c[2]: a product of random factors
                const x = try ops.alloc(T, al, m, c[2]);
                defer ops.free(T, al, x);
                const y = try ops.alloc(T, al, c[2], n);
                defer ops.free(T, al, y);
                check.fill(T, x.ptr[0 .. m * c[2]], 3);
                check.fill(T, y.ptr[0 .. c[2] * n], 5);
                try ops.gemm(al, check.io, T, .n, .n, ops.one(T), x, y, ops.zero(T), a0);
            }
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
