//! LQ: `?GELQF`, `?ORGLQ` and `?ORMLQ` behind Reference-LAPACK's arguments
//! (`dgelqf.f`, `dorglq.f`, `dormlq.f` and their `z` forms, 3.12.0 as
//! OpenBLAS v0.3.34-317 ships it), by libflame 5.2.0's `FLA_LQ_UT` method
//! (`src/map/lapack2flamec/FLA_gelqf.c`; no code ported) on faer's
//! reflections (`householder.zig`, MIT) run on `Aᴴ`.

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const ops = @import("ops.zig");
const blas_calls = @import("blas_calls.zig");
const hh = @import("householder.zig");
const qr = @import("qr.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = ops.Mat;
const Float = types.Float;
const isComplex = types.isComplex;
const Trans = types.Trans;
const Side = types.Side;
const Error = ops.Error;

/// What `gelqf` returns.
pub const GelqfError = error{ InvalidM, InvalidN, InvalidLda } || Allocator.Error || Io.Cancelable;
/// What `orglq` returns: `InvalidN` unless `n >= m`, `InvalidK` unless
/// `0 <= k <= m`.
pub const OrglqError = error{ InvalidM, InvalidN, InvalidK, InvalidLda } || Allocator.Error || Io.Cancelable;
/// What `ormlq` returns: `InvalidTrans` for other than `n`, `t` or `c`
/// (complex: `n` or `c`), `InvalidK` unless `0 <= k <=` the side's order.
pub const OrmlqError = error{ InvalidTrans, InvalidM, InvalidN, InvalidK, InvalidLda, InvalidLdc } || Allocator.Error || Io.Cancelable;

/// `?GELQF`: `A = L Q` for the `m x n` `A`: `L` on and below the diagonal,
/// `Q = H(k-1)ᴴ ... H(0)ᴴ` (`k = min(m, n)`) as reflections right of it
/// (`conj(v)` in row `i` past column `i`) and their `tau`s in `tau` (`k`).
pub fn gelqf(allocator: Allocator, io: Io, comptime T: type, m: isize, n: isize, a: [*]T, lda: isize, tau: [*]T) GelqfError!void {
    comptime ops.requireFloat(T, "gelqf");
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, m)) return error.InvalidLda;
    if (m == 0 or n == 0) return;
    if (!util.fitsMatrix(m, n, lda, @sizeOf(T))) return error.InvalidLda;
    const am: Mat(T) = .{ .ptr = a, .rows = @intCast(m), .cols = @intCast(n), .ld = @intCast(lda) };
    const b = try ops.alloc(T, allocator, am.cols, am.rows);
    defer ops.free(T, allocator, b);
    adjoint(T, am, b);
    try factor(allocator, io, T, b, tau);
    adjoint(T, b, am);
}

/// `dst := srcᴴ`.
fn adjoint(comptime T: type, src: Mat(T), dst: Mat(T)) void {
    for (0..src.cols) |j| for (0..src.rows) |i| {
        dst.at(j, i).* = ops.conj(T, src.get(i, j));
    };
}

/// QR of `a` in place by blocks with the UT transform, LAPACK's `tau`s
/// (`min(rows, cols)`) out.
pub fn factor(al: Allocator, io: Io, comptime T: type, a: Mat(T), tau: [*]T) Error!void {
    const F = Float(T);
    const m = a.rows;
    const n = a.cols;
    const k = @min(m, n);
    if (k == 0) return;
    const bs = qr.blockSize(m, n);
    const h = try ops.alloc(T, al, bs, bs);
    defer ops.free(T, al, h);
    const w = try al.alloc(T, n);
    defer al.free(w);
    var j: usize = 0;
    while (j < k) : (j += bs) {
        const jb = @min(bs, k - j);
        // the panel, one reflection at a time
        for (j..j + jb) |c| {
            const tail = m - c - 1;
            const v = a.col(c) + c + 1;
            const info = hh.make(io, T, a.at(c, c), tail, v, 1, null, 0);
            h.at(c - j, c - j).* = ops.scalar(T, info.tau);
            const tau_inv = 1 / info.tau;
            tau[c] = ops.scalar(T, tau_inv);
            const rest = j + jb - c - 1;
            if (tau_inv < std.math.floatMin(F) or rest == 0) continue;
            // w := A01 + vᴴ A11, w := -w / tau; A01 += w, A11 += v wᵀ
            for (w[0..rest], 0..) |*x, i| x.* = a.get(c, c + 1 + i);
            try blas_calls.gemv(al, io, T, .u, ops.one(T), a.sub(c + 1, c + 1, tail, rest), v, 1, ops.one(T), w.ptr, 1);
            for (w[0..rest], 0..) |*x, i| {
                x.* = ops.neg(T, ops.mulReal(T, x.*, tau_inv));
                const p = a.at(c, c + 1 + i);
                p.* = ops.add(T, p.*, x.*);
            }
            try blas_calls.ger(al, io, T, ops.one(T), v, 1, w.ptr, 1, a.sub(c + 1, c + 1, tail, rest));
        }
        if (j + jb >= n) continue;
        // the block's factor, then the rest: M := Hᴴ M
        const vb = a.sub(j, j, m - j, jb);
        const t = h.sub(0, 0, jb, jb);
        try hh.blockFactor(al, io, T, vb, t);
        try hh.applyBlock(al, io, T, vb, t, false, true, a.sub(j, j + jb, m - j, n - j - jb));
    }
}

/// `V = (A's first k rows)ᴴ`, `nq x k`, allocated.
fn basisOf(al: Allocator, comptime T: type, a: Mat(T), k: usize) Allocator.Error!Mat(T) {
    const v = try ops.alloc(T, al, a.cols, k);
    adjoint(T, a.sub(0, 0, k, a.cols), v);
    return v;
}

/// `?ORGLQ`, `?UNGLQ`: the `m x n` `Q` with orthonormal rows, the first `m`
/// of `H(k-1)ᴴ ... H(0)ᴴ`, over the reflections `gelqf` left in `A`'s first
/// `k` rows.
pub fn orglq(allocator: Allocator, io: Io, comptime T: type, m: isize, n: isize, k: isize, a: [*]T, lda: isize, tau: [*]const T) OrglqError!void {
    comptime ops.requireFloat(T, "orglq");
    if (m < 0) return error.InvalidM;
    if (n < m) return error.InvalidN;
    if (k < 0 or k > m) return error.InvalidK;
    if (lda < @max(1, m)) return error.InvalidLda;
    if (m == 0) return;
    if (!util.fitsMatrix(m, n, lda, @sizeOf(T))) return error.InvalidLda;
    const am: Mat(T) = .{ .ptr = a, .rows = @intCast(m), .cols = @intCast(n), .ld = @intCast(lda) };
    const uk: usize = @intCast(k);
    if (uk == 0) return am.identity();
    const v = try basisOf(allocator, T, am, uk);
    defer ops.free(T, allocator, v);
    const h = try qr.factorsOf(allocator, io, T, v, tau);
    defer ops.free(T, allocator, h);
    // Qᴴ's first m columns: Q_b I, Q_b = H(0) ... H(k-1)
    const g = try ops.alloc(T, allocator, am.cols, am.rows);
    defer ops.free(T, allocator, g);
    g.identity();
    try hh.applySequence(allocator, io, T, v, h, false, g);
    adjoint(T, g, am);
}

/// `?ORMLQ`, `?UNMLQ`: `C := op(Q) C` (side `l`) or `C op(Q)` (side `r`),
/// `op` by `trans` (`t` and `c` both `Qᴴ`), for the `m x n` `C` and the
/// `Q` of `k` reflections in `A`'s first `k` rows (`m` columns for `l`, `n`
/// for `r`). Only `A`'s strictly upper part is read. A right-side product
/// works on `Cᴴ`, so takes `m n` elements of workspace.
pub fn ormlq(allocator: Allocator, io: Io, comptime T: type, side: Side, trans: Trans, m: isize, n: isize, k: isize, a: [*]const T, lda: isize, tau: [*]const T, c: [*]T, ldc: isize) OrmlqError!void {
    comptime ops.requireFloat(T, "ormlq");
    switch (trans) {
        .n, .c => {},
        .t => if (comptime isComplex(T)) return error.InvalidTrans,
        .r, .o, .u, .s, .d => return error.InvalidTrans,
    }
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    const nq = if (side == .l) m else n;
    if (k < 0 or k > nq) return error.InvalidK;
    if (lda < @max(1, k)) return error.InvalidLda;
    if (ldc < @max(1, m)) return error.InvalidLdc;
    if (m == 0 or n == 0 or k == 0) return;
    if (!util.fitsMatrix(k, nq, lda, @sizeOf(T))) return error.InvalidLda;
    if (!util.fitsMatrix(m, n, ldc, @sizeOf(T))) return error.InvalidLdc;
    const am: Mat(T) = .{ .ptr = @constCast(a), .rows = @intCast(k), .cols = @intCast(nq), .ld = @intCast(lda) };
    const cm: Mat(T) = .{ .ptr = c, .rows = @intCast(m), .cols = @intCast(n), .ld = @intCast(ldc) };
    const v = try basisOf(allocator, T, am, @intCast(k));
    defer ops.free(T, allocator, v);
    const h = try qr.factorsOf(allocator, io, T, v, tau);
    defer ops.free(T, allocator, h);
    // Q = Q_bᴴ; Q C = Q_bᴴ C; C Q = (Q_b Cᴴ)ᴴ
    const adj = trans != .n;
    switch (side) {
        .l => if (adj) try hh.applySequence(allocator, io, T, v, h, false, cm) else try hh.applySequenceTransposed(allocator, io, T, v, h, true, cm),
        .r => {
            const wk = try ops.alloc(T, allocator, cm.cols, cm.rows);
            defer ops.free(T, allocator, wk);
            adjoint(T, cm, wk);
            if (adj) try hh.applySequenceTransposed(allocator, io, T, v, h, true, wk) else try hh.applySequence(allocator, io, T, v, h, false, wk);
            adjoint(T, wk, cm);
        },
    }
}

// ---- tests -----------------------------------------------------------------------

const testing = std.testing;
const check = @import("../blas/_test/check.zig");
const lcheck = @import("_test/check.zig");
const test_io = @import("../blas/_test/io.zig");
const blas = @import("../blas/__root__.zig");

const test_types = lcheck.test_types;

const thresh = lcheck.thresh;

/// `C := op(A) op(B)` by loops, in `T`.
fn mul(comptime T: type, ta: bool, a: Mat(T), tb: bool, b: Mat(T), c: Mat(T)) void {
    const k = if (ta) a.rows else a.cols;
    for (0..c.cols) |j| for (0..c.rows) |i| {
        var s = ops.zero(T);
        for (0..k) |p| {
            const x = if (ta) ops.conj(T, a.get(p, i)) else a.get(i, p);
            const y = if (tb) ops.conj(T, b.get(j, p)) else b.get(p, j);
            s = ops.add(T, s, ops.mul(T, x, y));
        }
        c.at(i, j).* = s;
    };
}

/// `dlqt01.f`'s ratios: `‖L - A Qᴴ‖ / (n ‖A‖ eps)` and `‖I - Q Qᴴ‖ /
/// (n eps)` below the threshold, `Q` (`n x n`) from `orglq`; then
/// `dlqt03.f`'s: `ormlq` against `Q` for every side and `op`.
fn lqCase(comptime T: type, al: Allocator, m: usize, n: usize, seed: u64) !void {
    const F = Float(T);
    const eps = std.math.floatEps(F);
    const k = @min(m, n);
    const a0 = try ops.alloc(T, al, m, n);
    defer ops.free(T, al, a0);
    check.fill(T, a0.ptr[0 .. m * n], seed);
    const a = try ops.alloc(T, al, m, n);
    defer ops.free(T, al, a);
    a.copyFrom(a0);
    const tau = try al.alloc(T, @max(k, 1));
    defer al.free(tau);
    const mi: isize = @intCast(m);
    const ni: isize = @intCast(n);
    const ki: isize = @intCast(k);
    try gelqf(al, check.io, T, mi, ni, a.ptr, @intCast(a.ld), tau.ptr);
    var anorm: F = 0;
    for (a0.ptr[0 .. m * n]) |v| anorm = @max(anorm, ops.abs(T, v));
    const q = try ops.alloc(T, al, n, n);
    defer ops.free(T, al, q);
    q.fill(ops.zero(T));
    for (0..n) |j| for (0..k) |i| {
        q.at(i, j).* = a.get(i, j);
    };
    try orglq(al, check.io, T, ni, ni, ki, q.ptr, ni, tau.ptr);
    const nf: F = @floatFromInt(@max(n, 1));
    // L = A Qᴴ, zero right of its diagonal
    const l = try ops.alloc(T, al, m, n);
    defer ops.free(T, al, l);
    mul(T, false, a0, true, q, l);
    const want = try ops.alloc(T, al, m, n);
    defer ops.free(T, al, want);
    for (0..n) |j| for (0..m) |i| {
        want.at(i, j).* = if (j <= i) a.get(i, j) else ops.zero(T);
    };
    try lcheck.expectNear(T, "gelqf L", l, want, thresh * nf * eps * @max(anorm, 1));
    const g = try ops.alloc(T, al, n, n);
    defer ops.free(T, al, g);
    mul(T, false, q, true, q, g);
    const id = try ops.alloc(T, al, n, n);
    defer ops.free(T, al, id);
    id.identity();
    try lcheck.expectNear(T, "orglq Q Qᴴ", g, id, thresh * nf * eps);
    // ormlq on random C, both sides, against Q's products
    if (k == 0) return;
    const ops_list: []const types.Trans = if (comptime isComplex(T)) &.{ .n, .c } else &.{ .n, .t };
    for ([_]Side{ .l, .r }) |side| for (ops_list) |tr| {
        const cr = if (side == .l) n else 5;
        const cc = if (side == .l) 5 else n;
        const c = try ops.alloc(T, al, cr, cc);
        defer ops.free(T, al, c);
        check.fill(T, c.ptr[0 .. cr * cc], seed + 9);
        const r = try ops.alloc(T, al, cr, cc);
        defer ops.free(T, al, r);
        // Q here is orglq's full n x n, op(Q) = Q or Qᴴ
        if (side == .l) mul(T, tr != .n, q, false, c, r) else mul(T, false, c, tr != .n, q, r);
        try ormlq(al, check.io, T, side, tr, @intCast(cr), @intCast(cc), ki, a.ptr, @intCast(a.ld), tau.ptr, c.ptr, @intCast(c.ld));
        try lcheck.expectNear(T, "ormlq", c, r, thresh * nf * eps * 4);
    };
}

test "lapack gelqf: L = A Qᴴ and Q Qᴴ = I (orglq), ormlq on both sides, square, wide, tall, blocked" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        for ([_][2]usize{ .{ 1, 1 }, .{ 1, 4 }, .{ 3, 5 }, .{ 5, 3 }, .{ 10, 10 }, .{ 40, 70 }, .{ 70, 40 } }) |mn| {
            try lqCase(T, al, mn[0], mn[1], mn[0] * 97 + mn[1]);
        }
    }
}

test "lapack gelqf: a zero row and a zero matrix give tau = 0; empty; bad arguments" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        var a = [_]T{ ops.zero(T), ops.zero(T), ops.zero(T), ops.zero(T) };
        var tau = [_]T{ ops.one(T), ops.one(T) };
        try gelqf(al, check.io, T, 2, 2, &a, 2, &tau);
        for (tau) |t| try testing.expect(ops.isZero(T, t));
        for (a) |v| try testing.expect(ops.isZero(T, v));
        try gelqf(al, check.io, T, 0, 2, &a, 1, &tau);
    }
    var a: [4]f64 = undefined;
    var tau: [2]f64 = undefined;
    try testing.expectError(error.InvalidM, gelqf(al, check.io, f64, -1, 1, &a, 1, &tau));
    try testing.expectError(error.InvalidN, gelqf(al, check.io, f64, 1, -1, &a, 1, &tau));
    try testing.expectError(error.InvalidLda, gelqf(al, check.io, f64, 2, 2, &a, 1, &tau));
    try testing.expectError(error.InvalidN, orglq(al, check.io, f64, 2, 1, 1, &a, 2, &tau));
    try testing.expectError(error.InvalidK, orglq(al, check.io, f64, 2, 2, 3, &a, 2, &tau));
    try testing.expectError(error.InvalidLda, orglq(al, check.io, f64, 2, 2, 1, &a, 1, &tau));
    try testing.expectError(error.InvalidTrans, ormlq(al, check.io, f64, .l, .r, 2, 2, 1, &a, 1, &tau, &a, 2));
    try testing.expectError(error.InvalidTrans, ormlq(al, check.io, types.Complex64, .l, .t, 0, 0, 0, @ptrCast(&a), 1, @ptrCast(&tau), @ptrCast(&a), 1));
    try testing.expectError(error.InvalidK, ormlq(al, check.io, f64, .r, .n, 2, 1, 2, &a, 2, &tau, &a, 2));
    try testing.expectError(error.InvalidLda, ormlq(al, check.io, f64, .l, .n, 2, 2, 2, &a, 1, &tau, &a, 2));
    try testing.expectError(error.InvalidLdc, ormlq(al, check.io, f64, .l, .n, 2, 2, 1, &a, 1, &tau, &a, 1));
}

test "lapack gelqf: one worker's bits at any worker count" {
    if (comptime !test_io.hosted) return error.SkipZigTest;
    _ = blas.coreName();
    const al = testing.allocator;
    var runs: test_io.WorkerRuns = .init(al);
    defer runs.deinit();
    inline for (.{ f64, types.Complex32 }) |T| {
        const m = 90;
        const n = 300;
        const a0 = try al.alloc(T, m * n);
        defer al.free(a0);
        check.fill(T, a0, 5);
        var outs: [3][]T = undefined;
        for (&outs) |*o| o.* = try al.alloc(T, m * n + m);
        defer for (outs) |o| al.free(o);
        for (0..test_io.worker_runs.len) |c| {
            const io = runs.io(c);
            @memcpy(outs[c][0 .. m * n], a0);
            try gelqf(al, io, T, m, n, outs[c].ptr, m, outs[c].ptr + m * n);
        }
        try test_io.expectSameBytes(outs);
    }
}
