//! Eigendecompositions: LAPACK's `?SYEVD`, `?HEEVD`, `?SYEV` and `?HEEV`
//! (Reference-LAPACK's `dsyevd.f`, `dsyev.f` and their `s`, `c`, `z` forms,
//! as OpenBLAS v0.3.34-317 ships them in `lapack-netlib/`: arguments,
//! `xerbla` order, the scaling of `A`) over faer v0.24.4's
//! (<https://codeberg.org/sarah-quinones/faer>, tag `faer-v0.24.4`, MIT)
//! `faer/src/linalg/evd/mod.rs`.
//!
//! Self-adjoint: tridiagonalize (`tridiag.zig`); take the real tridiagonal
//! (a complex one's subdiagonal by its moduli, the phases moved into the
//! vectors); decompose it (`tridiag_evd.zig`) by divide and conquer from
//! 128 up for `syevd` (the reference's `xSTEDC`), by QR iteration for
//! `syev` (its `xSTEQR`) and for values alone; apply the reflections.
//!
//! Departures, all in design.md's LAPACK section: no `WORK`, `RWORK` or
//! `IWORK` (the allocator gives workspace); `NoConvergence` for the
//! reference's `INFO > 0`, and for a NaN or infinity in `A`, which faer's
//! solvers cannot take; with `jobz` `n`, `A` is left as it was.

const std = @import("std");
const ops = @import("ops.zig");
const hh = @import("householder.zig");
const qr = @import("qr.zig");
const tridiag = @import("tridiag.zig");
const tevd = @import("tridiag_evd.zig");
const util = @import("../../blas/util.zig");
const chol = @import("../cholesky.zig");
const types = @import("../../blas/types.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = ops.Mat;
const Float = ops.Float;
const isComplex = ops.isComplex;
const Uplo = types.Uplo;

/// `JOBZ`, `JOBVL`, `JOBVR`: vectors (`v`) or not (`n`).
pub const Jobz = enum { n, v };

/// What `syevd`, `heevd`, `syev` and `heev` return.
pub const SyevError = error{ InvalidUplo, InvalidN, InvalidLda, NoConvergence } || Allocator.Error || Io.Cancelable;

/// `SelfAdjointEvdParams::auto`'s `recursion_threshold`: below it, divide
/// and conquer is QR iteration.
pub const dc_threshold = 128;

/// `?SYEVD`, `?HEEVD`: the eigenvalues of the `n x n` symmetric (Hermitian)
/// `A`, from its `uplo` triangle, into `w` (`n`, nondecreasing); with
/// `jobz` `v` the orthonormal eigenvectors over `A`. Divide and conquer.
pub fn syevd(allocator: Allocator, io: Io, comptime T: type, jobz: Jobz, uplo: Uplo, n: isize, a: [*]T, lda: isize, w: [*]Float(T)) SyevError!void {
    return driver(allocator, io, T, jobz, uplo, n, a, lda, w, dc_threshold);
}

/// `?SYEV`, `?HEEV`: as `syevd`, the tridiagonal by QR iteration.
pub fn syev(allocator: Allocator, io: Io, comptime T: type, jobz: Jobz, uplo: Uplo, n: isize, a: [*]T, lda: isize, w: [*]Float(T)) SyevError!void {
    return driver(allocator, io, T, jobz, uplo, n, a, lda, w, std.math.maxInt(usize));
}

fn driver(al: Allocator, io: Io, comptime T: type, jobz: Jobz, uplo: Uplo, n: isize, a: [*]T, lda: isize, w: [*]Float(T), threshold: usize) SyevError!void {
    comptime chol.floatType(T, "syevd");
    if (uplo.isRev()) return error.InvalidUplo;
    if (n < 0) return error.InvalidN;
    if (lda < @max(1, n)) return error.InvalidLda;
    if (n == 0) return;
    if (!util.fitsMatrix(n, n, lda, @sizeOf(T))) return error.InvalidLda;
    const un: usize = @intCast(n);
    const am: Mat(T) = .{ .ptr = a, .rows = un, .cols = un, .ld = @intCast(lda) };
    return symmetric(al, io, T, jobz == .v, uplo.isLower(), am, w[0..un], threshold);
}

/// The drivers after their checks: `dsyevd.f:292-307`'s scaling into
/// `[sqrt(SMLNUM), sqrt(BIGNUM)]`, faer's `self_adjoint_evd`, the values
/// scaled back.
pub fn symmetric(al: Allocator, io: Io, comptime T: type, vectors: bool, lower: bool, a: Mat(T), w: []Float(T), threshold: usize) (ops.Error || error{NoConvergence})!void {
    const F = Float(T);
    const n = a.rows;
    // the triangle, as the lower one (`dlansy('M')` on the way)
    const trid = try ops.alloc(T, al, n, n);
    defer ops.free(T, al, trid);
    var anrm: F = 0;
    for (0..n) |j| for (j..n) |i| {
        const x = if (lower) a.get(i, j) else ops.conj(T, a.get(j, i));
        trid.at(i, j).* = x;
        const v = ops.abs(T, x);
        if (!std.math.isFinite(v)) return error.NoConvergence;
        anrm = @max(anrm, v);
    };
    const smlnum = std.math.floatMin(F) / std.math.floatEps(F);
    const rmin = @sqrt(smlnum);
    const rmax = @sqrt(1 / smlnum);
    const sigma: ?F = if (anrm > 0 and anrm < rmin) rmin / anrm else if (anrm > rmax) rmax / anrm else null;
    if (sigma) |s| ops.lascl(T, 1, s, trid);
    try selfAdjointEvd(al, io, T, trid, w, if (vectors) a else null, threshold);
    if (sigma) |s| for (w) |*x| {
        x.* *= 1 / s;
    };
}

/// `self_adjoint_evd`: `trid`'s lower triangle (overwritten) to the
/// eigenvalues in `s`, nondecreasing, and the vectors in `u` if given.
pub fn selfAdjointEvd(al: Allocator, io: Io, comptime T: type, trid: Mat(T), s: []Float(T), u: ?Mat(T), threshold: usize) (ops.Error || error{NoConvergence})!void {
    const F = Float(T);
    const n = trid.rows;
    if (n == 0) return;
    const bs = qr.blockSize(n, n);
    const h = try ops.alloc(T, al, bs, @max(n, 2) - 1);
    defer ops.free(T, al, h);
    try tridiag.tridiagonalize(al, io, T, trid, h.sub(0, 0, bs, n - 1));
    const real = try al.alloc(F, 2 * n);
    defer al.free(real);
    const diag = real[0..n];
    const offdiag = real[n..];
    for (0..n) |i| {
        diag[i] = ops.re(T, trid.get(i, i));
        offdiag[i] = if (i + 1 < n) (if (comptime isComplex(T)) ops.abs(T, trid.get(i + 1, i)) else trid.get(i + 1, i)) else 0;
    }
    const um = u orelse {
        try tevd.qrAlgorithm(F, diag, offdiag, null);
        @memcpy(s, diag);
        return;
    };
    const ur: Mat(F) = if (comptime isComplex(T)) try ops.alloc(F, al, n, n) else um;
    defer if (comptime isComplex(T)) ops.free(F, al, ur);
    if (n < threshold) {
        try tevd.qrAlgorithm(F, diag, offdiag, ur);
    } else {
        try tevd.divideAndConquer(al, io, F, diag, offdiag, ur, threshold);
    }
    if (comptime isComplex(T)) {
        // U = diag(phases) Ur, the phases the subdiagonal's
        var x = ops.one(T);
        for (0..n) |i| {
            if (i > 0) x = normalized(T, ops.mul(T, trid.get(i, i - 1), x));
            for (0..n) |j| um.at(i, j).* = ops.mulReal(T, x, ur.get(i, j));
        }
    }
    if (n > 1) try hh.applySequence(al, io, T, trid.sub(1, 0, n - 1, n - 1), h.sub(0, 0, bs, n - 1), false, um.sub(1, 0, n - 1, n));
    @memcpy(s, diag);
}

/// `x / |x|` (1 for 0). Faer multiplies by `1 / |x|`, infinite for a
/// denormal `x`; that divides (as `svd.zig`'s).
fn normalized(comptime T: type, x: T) T {
    if (ops.isZero(T, x)) return ops.one(T);
    const norm1 = @max(@abs(x.re), @abs(x.im));
    const y: T = if (norm1 >= std.math.floatMin(Float(T))) ops.mulReal(T, x, 1 / norm1) else .init(x.re / norm1, x.im / norm1);
    return ops.mulReal(T, y, 1 / ops.abs(T, y));
}

// ---- tests -----------------------------------------------------------------------

const testing = std.testing;
const check = @import("../../blas/_test/blas1/check.zig");
const blas = @import("../../blas/__root__.zig");
const threading = @import("../../blas/threading.zig");

const test_types = .{ f32, f64, types.Complex32, types.Complex64 };

/// Faer's `test_self_adjoint_evd` tolerance: `ApproxEq::eps()` (`128 eps`)
/// times `sqrt(8 n)`, against `max(1, ‖A‖max)` per element.
fn tolerance(comptime F: type, n: usize, amax: F) F {
    return 128 * std.math.floatEps(F) * @sqrt(8 * @as(F, @floatFromInt(n))) * @max(amax, 1);
}

/// `‖A V - V Λ‖`, `‖Vᴴ V - I‖` per element within `tolerance`, `w`
/// nondecreasing.
fn expectSym(comptime T: type, al: Allocator, a0: Mat(T), w: []const Float(T), v: Mat(T), what: []const u8) !void {
    const F = Float(T);
    const n = a0.rows;
    for (1..n) |i| if (!(w[i - 1] <= w[i])) {
        check.say("{s} {s} n={d}: w out of order at {d}\n", .{ what, @typeName(T), n, i });
        return error.TestUnexpectedResult;
    };
    const r = try ops.alloc(T, al, n, n);
    defer ops.free(T, al, r);
    for (0..n) |j| for (0..n) |i| {
        r.at(i, j).* = ops.mulReal(T, v.get(i, j), w[j]);
    };
    try ops.gemm(al, check.io, T, .n, .n, ops.neg(T, ops.one(T)), a0, v, ops.one(T), r);
    const tol = tolerance(F, n, ops.normMax(T, a0));
    if (!(ops.normMax(T, r) <= tol)) {
        check.say("{s} {s} n={d}: A V - V L = {e}, tol {e}\n", .{ what, @typeName(T), n, ops.normMax(T, r), tol });
        return error.TestUnexpectedResult;
    }
    try ops.gemm(al, check.io, T, .c, .n, ops.one(T), v, v, ops.zero(T), r);
    for (0..n) |i| r.at(i, i).* = ops.sub(T, r.get(i, i), ops.one(T));
    if (!(ops.normMax(T, r) <= tolerance(F, n, 1))) {
        check.say("{s} {s} n={d}: Vh V - I = {e}\n", .{ what, @typeName(T), n, ops.normMax(T, r) });
        return error.TestUnexpectedResult;
    }
}

/// A random Hermitian `A + Aᴴ`.
fn hermitian(comptime T: type, a: Mat(T), seed: u64) void {
    const n = a.rows;
    check.fill(T, a.ptr[0 .. a.ld * n], seed);
    for (0..n) |j| for (0..j + 1) |i| {
        const s = ops.add(T, a.get(i, j), ops.conj(T, a.get(j, i)));
        a.at(i, j).* = s;
        a.at(j, i).* = ops.conj(T, s);
    };
}

const Routine = enum { syevd, syev };

/// The routine on a copy of `a0` (padded `lda`), both triangles hidden
/// behind NaN but the one `uplo` names; the vectors into `v`.
fn decompose(comptime T: type, al: Allocator, io: Io, r: Routine, jobz: Jobz, uplo: Uplo, a0: Mat(T), w: []Float(T), v: ?Mat(T)) !void {
    const n = a0.rows;
    const lda = n + 1;
    const a = try ops.alloc(T, al, lda, n);
    defer ops.free(T, al, a);
    a.fill(chol.poisoned(T));
    for (0..n) |j| for (0..n) |i| {
        if (if (uplo == .l) i >= j else i <= j) a.at(i, j).* = a0.get(i, j);
    };
    switch (r) {
        .syevd => try syevd(al, io, T, jobz, uplo, @intCast(n), a.ptr, @intCast(lda), w.ptr),
        .syev => try syev(al, io, T, jobz, uplo, @intCast(n), a.ptr, @intCast(lda), w.ptr),
    }
    for (0..n) |j| try testing.expect(chol.isNanT(T, a.get(n, j)));
    if (v) |x| x.copyFrom(a.sub(0, 0, n, n));
}

test "faer evd: syevd and syev, A V = V L with V orthonormal, both triangles, faer's sizes and past 128" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        for ([_]usize{ 1, 2, 4, 10, 15, 20, 50, 130, 150 }) |n| {
            const a0 = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, a0);
            hermitian(T, a0, 7 * n + 1);
            const w = try al.alloc(F, n);
            defer al.free(w);
            const w2 = try al.alloc(F, n);
            defer al.free(w2);
            const v = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, v);
            for ([_]Routine{ .syevd, .syev }) |r| for ([_]Uplo{ .l, .u }) |uplo| {
                try decompose(T, al, check.io, r, .v, uplo, a0, w, v);
                try expectSym(T, al, a0, w, v, @tagName(r));
                try decompose(T, al, check.io, r, .n, uplo, a0, w2, null);
                for (w, w2) |x, y| try testing.expect(@abs(x - y) <= tolerance(F, n, ops.normMax(T, a0)));
            };
        }
    }
}

test "faer evd: special matrices, faer's test_special: zero, ones, identity; rank one, clustered, graded" {
    const al = testing.allocator;
    inline for (test_types) |T| {
        const F = Float(T);
        for ([_]usize{ 1, 2, 4, 15, 50, 150 }) |n| for (0..6) |kind| {
            const a0 = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, a0);
            switch (kind) {
                0 => a0.fill(ops.zero(T)),
                1 => a0.fill(ops.one(T)),
                2 => a0.identity(),
                // clustered: I + eps-sized noise
                3 => {
                    hermitian(T, a0, n);
                    for (0..n) |j| for (a0.col(j)[0..n]) |*x| {
                        x.* = ops.mulReal(T, x.*, std.math.floatEps(F));
                    };
                    for (0..n) |i| a0.at(i, i).* = ops.add(T, a0.get(i, i), ops.one(T));
                },
                // graded: D H D, D from 1 down by 10 per two rows
                4 => {
                    hermitian(T, a0, n + 1);
                    for (0..n) |j| for (0..n) |i| {
                        const s = std.math.pow(F, 10, -@as(F, @floatFromInt(((i + j) / 4) % 8)));
                        a0.at(i, j).* = ops.mulReal(T, a0.get(i, j), s);
                    };
                },
                // tridiagonal Wilkinson W+ (a close pair at the top)
                else => {
                    a0.fill(ops.zero(T));
                    for (0..n) |i| {
                        const c = @as(F, @floatFromInt(n / 2)) - @as(F, @floatFromInt(i));
                        a0.at(i, i).* = ops.scalar(T, @abs(c));
                        if (i + 1 < n) {
                            a0.at(i + 1, i).* = ops.one(T);
                            a0.at(i, i + 1).* = ops.one(T);
                        }
                    }
                },
            }
            const w = try al.alloc(F, n);
            defer al.free(w);
            const v = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, v);
            for ([_]Routine{ .syevd, .syev }) |r| {
                decompose(T, al, check.io, r, .v, .l, a0, w, v) catch |e| {
                    check.say("{s} {s} n={d} kind {d}: {t}\n", .{ @tagName(r), @typeName(T), n, kind, e });
                    return e;
                };
                try expectSym(T, al, a0, w, v, @tagName(r));
            }
        };
    }
}

test "faer evd: OpenBLAS's test_laed3_nan.c: dsyevd on I + 10 v vᵀ, nearly degenerate, no NaN" {
    // OpenBLAS `test/test_laed3_nan.c` (71f22601): n = 8, v = (1, 1 + 1e-8,
    // 1, ...), jobz 'V', uplo 'U'; here through `syevd` (n below 128: QR
    // iteration), through divide and conquer at faer's test threshold 8 and
    // at 4, and at n = 200, where `syevd` divides
    const al = testing.allocator;
    for ([_]usize{ 8, 200 }) |n| {
        const a0 = try ops.alloc(f64, al, n, n);
        defer ops.free(f64, al, a0);
        const v = try al.alloc(f64, n);
        defer al.free(v);
        @memset(v, 1);
        v[1] = 1 + 1e-8;
        for (0..n) |j| for (0..n) |i| {
            a0.at(i, j).* = (if (i == j) @as(f64, 1) else 0) + 10 * v[i] * v[j];
        };
        const w = try al.alloc(f64, n);
        defer al.free(w);
        const a = try ops.alloc(f64, al, n, n);
        defer ops.free(f64, al, a);
        a.copyFrom(a0);
        try syevd(al, check.io, f64, .v, .u, @intCast(n), a.ptr, @intCast(n), w.ptr);
        for (w) |x| try testing.expect(!std.math.isNan(x));
        for (a.ptr[0 .. n * n]) |x| try testing.expect(!std.math.isNan(x));
        try expectSym(f64, al, a0, w, a, "laed3");
        for ([_]usize{ 4, 8 }) |thr| {
            a.copyFrom(a0);
            try symmetric(al, check.io, f64, true, false, a, w, thr);
            for (a.ptr[0 .. n * n]) |x| try testing.expect(!std.math.isNan(x));
            try expectSym(f64, al, a0, w, a, "laed3 dc");
        }
    }
}

test "faer evd: NaN and infinity are NoConvergence; scaled past the bounds; 1 x 1 and empty" {
    const al = testing.allocator;
    inline for (.{ f64, types.Complex32 }) |T| {
        const F = Float(T);
        for ([_]F{ std.math.nan(F), std.math.inf(F) }) |bad| {
            var a: [9]T = @splat(ops.one(T));
            a[1] = ops.scalar(T, bad);
            var w: [3]F = undefined;
            try testing.expectError(error.NoConvergence, syevd(al, check.io, T, .v, .l, 3, &a, 3, &w));
            // the other triangle is not read
            a[1] = ops.one(T);
            a[3] = ops.scalar(T, bad);
            try syevd(al, check.io, T, .n, .l, 3, &a, 3, &w);
        }
        const n = 20;
        for ([_]F{ std.math.floatMin(F) * 4, std.math.floatMax(F) / 4 }) |sc| {
            const a0 = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, a0);
            hermitian(T, a0, 3);
            var w0: [n]F = undefined;
            var w1: [n]F = undefined;
            const a = try ops.alloc(T, al, n, n);
            defer ops.free(T, al, a);
            a.copyFrom(a0);
            try syevd(al, check.io, T, .n, .l, n, a.ptr, n, &w0);
            for (0..n) |j| for (a0.col(j)[0..n]) |*x| {
                x.* = ops.mulReal(T, x.*, sc / 8);
            };
            try syevd(al, check.io, T, .n, .l, n, a0.ptr, n, &w1);
            const tol = tolerance(F, n, 1) * (sc / 8) * 8;
            for (w0, w1) |x, y| try testing.expect(@abs(x * (sc / 8) - y) <= tol);
        }
        var one = [_]T{ops.scalar(T, -3)};
        var w: [1]F = undefined;
        try syevd(al, check.io, T, .v, .u, 1, &one, 1, &w);
        try testing.expectEqual(@as(F, -3), w[0]);
        try testing.expectEqual(@as(F, 1), ops.re(T, one[0]));
    }
    const fa = testing.failing_allocator;
    var a: [4]f64 = @splat(std.math.nan(f64));
    var w: [2]f64 = @splat(std.math.nan(f64));
    try testing.expectError(error.InvalidUplo, syevd(fa, check.io, f64, .v, .v, -1, &a, 0, &w));
    try testing.expectError(error.InvalidN, syevd(fa, check.io, f64, .v, .l, -1, &a, 0, &w));
    try testing.expectError(error.InvalidLda, syev(fa, check.io, f64, .v, .u, 2, &a, 1, &w));
    try syevd(fa, check.io, f64, .v, .l, 0, &a, 1, &w);
    for (a ++ w) |x| try testing.expect(std.math.isNan(x));
}

test "faer evd: syevd gives one worker's bits at any worker count" {
    if (comptime !threading.hosted) return error.SkipZigTest;
    _ = blas.coreName();
    const al = testing.allocator;
    const cpus = blas.config.cpu_count;
    defer blas.config.cpu_count = cpus;
    defer threading.test_serial_max = null;
    var threaded: Io.Threaded = .init(al, .{ .async_limit = .limited(8) });
    defer threaded.deinit();
    inline for (.{ f64, types.Complex32 }) |T| {
        const n = 260;
        const a0 = try ops.alloc(T, al, n, n);
        defer ops.free(T, al, a0);
        hermitian(T, a0, 5);
        var outs: [3][]T = undefined;
        var ws: [3][]Float(T) = undefined;
        for (0..3) |c| {
            outs[c] = try al.alloc(T, n * n);
            ws[c] = try al.alloc(Float(T), n);
        }
        defer for (0..3) |c| {
            al.free(outs[c]);
            al.free(ws[c]);
        };
        for ([_]u32{ 1, 3, 8 }, 0..) |wk, c| {
            blas.config.cpu_count = wk;
            threading.test_serial_max = if (wk == 1) null else 4096;
            @memcpy(outs[c], a0.ptr[0 .. n * n]);
            try syevd(al, if (wk == 1) Io.failing else threaded.io(), T, .v, .l, n, outs[c].ptr, n, ws[c].ptr);
        }
        for (1..3) |c| {
            try testing.expectEqualSlices(u8, std.mem.sliceAsBytes(outs[0]), std.mem.sliceAsBytes(outs[c]));
            try testing.expectEqualSlices(u8, std.mem.sliceAsBytes(ws[0]), std.mem.sliceAsBytes(ws[c]));
        }
    }
}
