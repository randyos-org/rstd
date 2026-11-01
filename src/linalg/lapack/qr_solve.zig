//! Least squares by QR, from faer v0.24.4
//! (<https://codeberg.org/sarah-quinones/faer>, tag `faer-v0.24.4`, MIT):
//! `faer/src/linalg/qr/no_pivoting/solve.rs` (`solve_lstsq_in_place`:
//! `Qᴴ B`, then `R`'s solve), behind LAPACK's `?GELS` (Reference-LAPACK's
//! `dgels.f` and its `s`, `c`, `z` forms: arguments, `xerbla` order, the
//! scaling of `A` and `B`, the four cases and where the answer lands).
//!
//! Faer solves only the tall, untransposed case. The reference factors a
//! wide `A` as `L Q` (`xGELQF`); faer has no LQ, so a wide `A` is factored
//! as `Aᴴ = Q R` (`qr.zig`) and written back conjugate transposed, which is
//! `xGELQF`'s layout: `L = Rᴴ` on and below the diagonal, the reflections
//! conjugated in the rows right of it. The triangular solves are rstd's
//! `trtrs`, so a zero on `R`'s diagonal is `error.Singular`, the reference's
//! `INFO > 0`.

const std = @import("std");
const ops = @import("ops.zig");
const hh = @import("householder.zig");
const qr = @import("qr.zig");
const util = @import("../blas/util.zig");
const chol = @import("cholesky.zig");
const tri = @import("triangular.zig");
const types = @import("../blas/types.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = ops.Mat;
const Float = ops.Float;
const isComplex = ops.isComplex;
const Trans = types.Trans;

/// What `gels` returns: `InvalidTrans` for other than `n`, `t` or `c`
/// (complex: `n` or `c`); `Singular` when `A`'s triangular factor has a
/// zero on its diagonal (`A` is rank deficient).
pub const GelsError = error{ InvalidTrans, InvalidM, InvalidN, InvalidNrhs, InvalidLda, InvalidLdb, Singular } || Allocator.Error || Io.Cancelable;

/// `?GELS`: for the `m x n` `A` of full rank, the least-squares solution of
/// `op(A) X = B` where it is overdetermined, the minimum-norm one where
/// underdetermined (`trans` `t` and `c` both `Aᴴ`). `B` is `max(m, n) x
/// nrhs` and takes `X` in its first `n` rows (`m` for `Aᴴ`); for an
/// overdetermined system the rows past them hold each column's residual,
/// whose norm is the residual's. `A` takes `geqrf`'s factor (`m >= n`) or
/// `gelqf`'s (`m < n`). `Singular` (with the zero's 1-based index in
/// `pivot`) leaves `B` part solved, as the reference's `INFO > 0`.
pub fn gels(allocator: Allocator, io: Io, comptime T: type, trans: Trans, m: isize, n: isize, nrhs: isize, a: [*]T, lda: isize, b: [*]T, ldb: isize, pivot: ?*usize) GelsError!void {
    comptime chol.floatType(T, "gels");
    switch (trans) {
        .n, .c => {},
        .t => if (comptime isComplex(T)) return error.InvalidTrans,
        .r, .o, .u, .s, .d => return error.InvalidTrans,
    }
    if (m < 0) return error.InvalidM;
    if (n < 0) return error.InvalidN;
    if (nrhs < 0) return error.InvalidNrhs;
    if (lda < @max(1, m)) return error.InvalidLda;
    if (ldb < @max(1, @max(m, n))) return error.InvalidLdb;
    if (!util.fitsMatrix(m, n, lda, @sizeOf(T))) return error.InvalidLda;
    if (!util.fitsMatrix(@max(m, n), nrhs, ldb, @sizeOf(T))) return error.InvalidLdb;
    const um: usize = @intCast(m);
    const un: usize = @intCast(n);
    const bm: Mat(T) = .{ .ptr = b, .rows = @max(um, un), .cols = @intCast(nrhs), .ld = @intCast(ldb) };
    if (um == 0 or un == 0 or bm.cols == 0) return bm.fill(ops.zero(T));
    const am: Mat(T) = .{ .ptr = a, .rows = um, .cols = un, .ld = @intCast(lda) };
    solve(allocator, io, T, trans != .n, am, bm, pivot) catch |e| return switch (e) {
        error.Singular => error.Singular,
        error.OutOfMemory => error.OutOfMemory,
        error.Canceled => error.Canceled,
    };
}

const SolveError = error{Singular} || ops.Error;

/// `gels` after its checks; `adjoint` solves with `Aᴴ`.
fn solve(al: Allocator, io: Io, comptime T: type, adjoint: bool, a: Mat(T), b: Mat(T), pivot: ?*usize) SolveError!void {
    const F = Float(T);
    const m = a.rows;
    const n = a.cols;
    const nrhs = b.cols;
    const smlnum = std.math.floatMin(F) / std.math.floatEps(F);
    const anrm = ops.normMax(T, a);
    if (anrm == 0) return b.fill(ops.zero(T));
    const to_a = ops.scaleTo(F, anrm, smlnum);
    if (to_a) |t| ops.lascl(T, anrm, t, a);
    const brow = if (adjoint) n else m;
    const bb = b.sub(0, 0, brow, nrhs);
    const bnrm = ops.normMax(T, bb);
    const to_b = ops.scaleTo(F, bnrm, smlnum);
    if (to_b) |t| ops.lascl(T, bnrm, t, bb);
    // `Aᴴ` (wide) or `A` (tall), factored `Q R`
    const tall = m >= n;
    const fm = if (tall) a else try ops.alloc(T, al, n, m);
    defer if (!tall) ops.free(T, al, fm);
    if (!tall) transposeInto(T, a, fm);
    const k = fm.cols;
    const h = try ops.alloc(T, al, qr.blockSize(fm.rows, k), k);
    defer ops.free(T, al, h);
    try qr.factorZeroed(al, io, T, fm, h);
    if (!tall) transposeInto(T, fm, a);
    const ct: Trans = if (comptime isComplex(T)) .c else .t;
    // `A X = B` tall, or `Aᴴ X = B` wide: overdetermined, `R X = Qᴴ B`;
    // the others underdetermined, `X = Q [R⁻ᴴ B; 0]`
    const over = tall != adjoint;
    const big = fm.rows;
    if (over) {
        try hh.applySequenceTransposed(al, io, T, fm, h, true, b.sub(0, 0, big, nrhs));
        try trtrs(al, io, T, .n, fm, b, pivot);
    } else {
        try trtrs(al, io, T, ct, fm, b, pivot);
        b.sub(k, 0, big - k, nrhs).fill(ops.zero(T));
        try hh.applySequence(al, io, T, fm, h, false, b.sub(0, 0, big, nrhs));
    }
    const scllen = if (over) k else big;
    const x = b.sub(0, 0, scllen, nrhs);
    if (to_a) |t| ops.lascl(T, anrm, t, x);
    if (to_b) |t| ops.lascl(T, t, bnrm, x);
}

/// `dst := srcᴴ`.
fn transposeInto(comptime T: type, src: Mat(T), dst: Mat(T)) void {
    for (0..src.cols) |j| for (0..src.rows) |i| {
        dst.at(j, i).* = ops.conj(T, src.get(i, j));
    };
}

/// `trtrs` on `R` (`r`'s upper `k x k`) with `B`'s first `k` rows.
fn trtrs(al: Allocator, io: Io, comptime T: type, trans: Trans, r: Mat(T), b: Mat(T), pivot: ?*usize) SolveError!void {
    const k = r.cols;
    tri.trtrs(al, io, T, .u, trans, .n, @intCast(k), @intCast(b.cols), r.ptr, @intCast(r.ld), b.ptr, @intCast(b.ld), pivot) catch |e| return switch (e) {
        error.Singular => error.Singular,
        error.OutOfMemory => error.OutOfMemory,
        error.Canceled => error.Canceled,
        else => unreachable,
    };
}

// ---- tests -----------------------------------------------------------------------

const testing = std.testing;
const check = @import("../blas/_test/blas1/check.zig");

const test_types = .{ f32, f64, types.Complex32, types.Complex64 };

/// Reference-LAPACK's test threshold (`TESTING/dtest.in`, `THRESH` 30) on
/// its residual ratios (`dqrt16.f`, `dqrt17.f`), `‖·‖` the max norm here,
/// times 10 for `κ(A)` of a random matrix at aspect ratio up to 1 / 2 where
/// an answer is compared directly.
const thresh = 30;

/// `x := op(A) y` (`n` or `c`).
fn apply(comptime T: type, al: Allocator, adjoint: bool, a: Mat(T), y: Mat(T), x: Mat(T)) !void {
    try ops.gemm(al, check.io, T, if (adjoint) .c else .n, .n, ops.one(T), a, y, ops.zero(T), x);
}

test "faer gels: tall and wide, A and Aᴴ, least squares and minimum norm" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        const eps = std.math.floatEps(F);
        for ([_][2]usize{ .{ 1, 1 }, .{ 7, 3 }, .{ 3, 7 }, .{ 60, 25 }, .{ 25, 60 }, .{ 40, 40 }, .{ 130, 60 } }) |mn| for ([_]bool{ false, true }) |adj| {
            const m = mn[0];
            const n = mn[1];
            const nrhs = 3;
            // op(A) is rows x cols
            const rows = if (adj) n else m;
            const cols = if (adj) m else n;
            const a0 = try ops.alloc(T, al, m, n);
            defer ops.free(T, al, a0);
            check.fill(T, a0.ptr[0 .. m * n], 31 * m + n);
            const big = @max(m, n);
            const x0 = try ops.alloc(T, al, cols, nrhs);
            defer ops.free(T, al, x0);
            const b0 = try ops.alloc(T, al, big, nrhs);
            defer ops.free(T, al, b0);
            b0.fill(ops.zero(T));
            const rhs = b0.sub(0, 0, rows, nrhs);
            var rnorm: [nrhs]F = @splat(0);
            if (rows >= cols) {
                // B = op(A) x0 + r0, r0 orthogonal to op(A)'s range: x0 is the answer
                check.fill(T, x0.ptr[0 .. cols * nrhs], 5);
                const opa = try ops.alloc(T, al, rows, cols);
                defer ops.free(T, al, opa);
                for (0..cols) |j| for (0..rows) |i| {
                    opa.at(i, j).* = if (adj) ops.conj(T, a0.get(j, i)) else a0.get(i, j);
                };
                const tau = try al.alloc(T, cols);
                defer al.free(tau);
                const q = try ops.alloc(T, al, rows, rows);
                defer ops.free(T, al, q);
                q.sub(0, 0, rows, cols).copyFrom(opa);
                try qr.factorLapack(al, check.io, T, q.sub(0, 0, rows, cols), tau.ptr);
                try qr.generate(al, check.io, T, q, cols, tau.ptr);
                // r0: random in Q's trailing columns
                const w = try ops.alloc(T, al, rows - cols, nrhs);
                defer ops.free(T, al, w);
                check.fill(T, w.ptr[0 .. (rows - cols) * nrhs], 6);
                try ops.gemm(al, check.io, T, .n, .n, ops.one(T), q.sub(0, cols, rows, rows - cols), w, ops.zero(T), rhs);
                for (0..nrhs) |j| rnorm[j] = ops.nrm2(check.io, T, rows - cols, w.col(j), 1);
                try ops.gemm(al, check.io, T, .n, .n, ops.one(T), opa, x0, ops.one(T), rhs);
            } else {
                // x0 in op(A)ᴴ's range: the minimum-norm answer
                const z = try ops.alloc(T, al, rows, nrhs);
                defer ops.free(T, al, z);
                check.fill(T, z.ptr[0 .. rows * nrhs], 7);
                try apply(T, al, !adj, a0, z, x0);
                try apply(T, al, adj, a0, x0, rhs);
            }
            const a = try ops.alloc(T, al, m, n);
            defer ops.free(T, al, a);
            a.copyFrom(a0);
            const b = try ops.alloc(T, al, big, nrhs);
            defer ops.free(T, al, b);
            b.copyFrom(b0);
            try gels(al, check.io, T, if (adj) .c else .n, @intCast(m), @intCast(n), nrhs, a.ptr, @intCast(m), b.ptr, @intCast(big), null);
            const x = b.sub(0, 0, cols, nrhs);
            const tol = 10 * thresh * @as(F, @floatFromInt(big)) * eps * @max(ops.normMax(T, x0), 1);
            for (0..nrhs) |j| for (0..cols) |i| {
                const err = ops.abs(T, ops.sub(T, x.get(i, j), x0.get(i, j)));
                if (!(err <= tol)) {
                    check.say("gels {s} {d}x{d} adj {}: x({d},{d}) off by {e}, tol {e}\n", .{ @typeName(T), m, n, adj, i, j, err, tol });
                    return error.TestUnexpectedResult;
                }
            };
            // the residual's norm in the rows past the answer
            if (rows > cols) for (0..nrhs) |j| {
                const got = ops.nrm2(check.io, T, rows - cols, b.col(j) + cols, 1);
                try testing.expect(@abs(got - rnorm[j]) <= tol);
            };
        };
    }
}

test "faer gels: an inconsistent tall system's residual is orthogonal to A (dqrt17)" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        const m = 50;
        const n = 20;
        const nrhs = 2;
        const a0 = try ops.alloc(T, al, m, n);
        defer ops.free(T, al, a0);
        check.fill(T, a0.ptr[0 .. m * n], 41);
        const b0 = try ops.alloc(T, al, m, nrhs);
        defer ops.free(T, al, b0);
        check.fill(T, b0.ptr[0 .. m * nrhs], 42);
        const a = try ops.alloc(T, al, m, n);
        defer ops.free(T, al, a);
        a.copyFrom(a0);
        const b = try ops.alloc(T, al, m, nrhs);
        defer ops.free(T, al, b);
        b.copyFrom(b0);
        try gels(al, check.io, T, .n, m, n, nrhs, a.ptr, m, b.ptr, m, null);
        // r = B - A X; ‖Aᴴ r‖ / (max(m, n, nrhs) ‖A‖ ‖B‖ eps) < thresh
        const r = b0;
        try ops.gemm(al, check.io, T, .n, .n, ops.neg(T, ops.one(T)), a0, b.sub(0, 0, n, nrhs), ops.one(T), r);
        const g = try ops.alloc(T, al, n, nrhs);
        defer ops.free(T, al, g);
        try ops.gemm(al, check.io, T, .c, .n, ops.one(T), a0, r, ops.zero(T), g);
        const bscale = ops.normMax(T, b.sub(0, 0, m, nrhs));
        const ratio = ops.normMax(T, g) / (@as(F, m) * ops.normMax(T, a0) * @max(bscale, 1) * std.math.floatEps(F));
        try testing.expect(ratio < thresh);
    }
}

test "faer gels: rank deficient is Singular, a zero A zeroes B, empty zeroes B, tiny and huge scale" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex32 }) |T| {
        const F = Float(T);
        const m = 9;
        const n = 5;
        var a: [m * n]T = undefined;
        var b: [m * 2]T = undefined;
        check.fill(T, &a, 1);
        check.fill(T, &b, 2);
        @memset(a[2 * m .. 3 * m], ops.zero(T));
        var piv: usize = 0;
        try testing.expectError(error.Singular, gels(al, check.io, T, .n, m, n, 2, &a, m, &b, m, &piv));
        try testing.expectEqual(@as(usize, 3), piv);
        @memset(&a, ops.zero(T));
        try gels(al, check.io, T, .n, m, n, 2, &a, m, &b, m, null);
        for (b) |x| try testing.expect(ops.isZero(T, x));
        // m = 0: B's max(m, n) rows zeroed
        check.fill(T, &b, 3);
        try gels(al, check.io, T, .n, 0, 3, 2, &a, 1, &b, 3, null);
        for (b[0..6]) |x| try testing.expect(ops.isZero(T, x));
        // scaled past SMLNUM and BIGNUM: the unscaled answer, scaled back
        for ([_]F{ std.math.floatMin(F) * 16, std.math.floatMax(F) / 16 }) |sc| {
            var a1: [m * n]T = undefined;
            var b1: [m]T = undefined;
            check.fill(T, &a1, 4);
            check.fill(T, &b1, 5);
            var a2 = a1;
            var b2 = b1;
            for (&a2) |*x| x.* = ops.mulReal(T, x.*, sc);
            try gels(al, check.io, T, .n, m, n, 1, &a1, m, &b1, m, null);
            try gels(al, check.io, T, .n, m, n, 1, &a2, m, &b2, m, null);
            // relative to the answer's largest entry: 10 thresh eps
            var big: F = 0;
            for (b1[0..n]) |x| big = @max(big, ops.abs(T, x));
            for (b1[0..n], b2[0..n]) |x, y| {
                const want = ops.mulReal(T, x, 1 / sc);
                try testing.expect(ops.abs(T, ops.sub(T, y, want)) <= 10 * thresh * std.math.floatEps(F) * (big / sc));
            }
        }
        // NaN in A: it reaches X, nothing fails
        check.fill(T, &a, 6);
        check.fill(T, &b, 7);
        a[4] = ops.scalar(T, std.math.nan(F));
        try gels(al, check.io, T, .n, m, n, 2, &a, m, &b, m, null);
        try testing.expect(std.math.isNan(ops.re(T, b[0])));
    }
}

test "faer gels: reports the first bad argument" {
    const fa = testing.failing_allocator;
    var a: [4]f64 = @splat(0);
    var b: [4]f64 = @splat(0);
    try testing.expectError(error.InvalidTrans, gels(fa, check.io, f64, .r, -1, -1, -1, &a, 0, &b, 0, null));
    try testing.expectError(error.InvalidTrans, gels(fa, check.io, types.Complex64, .t, 1, 1, 1, @ptrCast(&a), 1, @ptrCast(&b), 1, null));
    try testing.expectError(error.InvalidM, gels(fa, check.io, f64, .n, -1, -1, -1, &a, 0, &b, 0, null));
    try testing.expectError(error.InvalidN, gels(fa, check.io, f64, .n, 1, -1, -1, &a, 0, &b, 0, null));
    try testing.expectError(error.InvalidNrhs, gels(fa, check.io, f64, .n, 1, 1, -1, &a, 0, &b, 0, null));
    try testing.expectError(error.InvalidLda, gels(fa, check.io, f64, .n, 2, 1, 1, &a, 1, &b, 0, null));
    try testing.expectError(error.InvalidLdb, gels(fa, check.io, f64, .n, 1, 2, 1, &a, 1, &b, 1, null));
}
