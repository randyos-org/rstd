//! Least squares by the SVD, `gelss` and `gelsd`: Reference-LAPACK 3.12.0
//! (`lapack-netlib/SRC/dgelss.f`, `zgelss.f`, `dgelsd.f`, `zgelsd.f` with
//! `dlalsd.f`, OpenBLAS v0.3.34-794; modified BSD), over rstd's SVD (faer's,
//! `svd.zig`): `U` and `V` are formed and applied, where the reference
//! applies them to `B` as it goes; the solution and rank are the same.

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const interface = @import("interface.zig");
const ops = @import("ops.zig");
const blas_calls = @import("blas_calls.zig");
const auxiliary = @import("auxiliary.zig");
const svd = @import("svd.zig");
const lsq_scale = @import("lsq_scale.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Mat = ops.Mat;
const Float = types.Float;

/// What `gelss` and `gelsd` return: a slice too short (`a`, `b`, or `s`
/// under `min(m, n)`) or a leading dimension too small, `NoConvergence`
/// from the SVD (the reference's `INFO > 0`), and the allocator's and
/// cancellation.
pub const GelssError = interface.ShapeError || error{NoConvergence} || Allocator.Error || Io.Cancelable;
pub const GelsdError = GelssError;

/// `?GELSS` (`dgelss.f`, `zgelss.f`): the minimum-norm solution of
/// `min ‖A X - B‖` for the `m x n` `A` of any rank, by the SVD, its
/// bidiagonal by QR iteration. `s` (`min(m, n)`) gets `A`'s singular
/// values, nonincreasing; those at or below `max(rcond s[0], sfmin)` count
/// as zero (`rcond < 0`: `eps`), and the rest give the rank returned. `B`
/// is `max(m, n) x nrhs`; `X` replaces its first `n` rows. `A` is
/// destroyed.
pub fn gelss(allocator: Allocator, io: Io, comptime T: type, m: usize, n: usize, nrhs: usize, a: []T, lda: usize, b: []T, ldb: usize, s: []Float(T), rcond: Float(T)) GelssError!usize {
    comptime ops.requireFloat(T, "gelss");
    return solve(allocator, io, T, .gelss, m, n, nrhs, a, lda, b, ldb, s, rcond);
}

/// `?GELSD` (`dgelsd.f`, `zgelsd.f`): `gelss` with the bidiagonal by divide
/// and conquer. Values at or below `rcond s[0]` count as zero (`rcond`
/// outside `(0, 1)`: `eps`), as `dlalsd.f` has it.
pub fn gelsd(allocator: Allocator, io: Io, comptime T: type, m: usize, n: usize, nrhs: usize, a: []T, lda: usize, b: []T, ldb: usize, s: []Float(T), rcond: Float(T)) GelsdError!usize {
    comptime ops.requireFloat(T, "gelsd");
    return solve(allocator, io, T, .gelsd, m, n, nrhs, a, lda, b, ldb, s, rcond);
}

const Driver = enum { gelss, gelsd };

fn solve(al: Allocator, io: Io, comptime T: type, comptime driver: Driver, m: usize, n: usize, nrhs: usize, a: []T, lda: usize, b: []T, ldb: usize, s: []Float(T), rcond: Float(T)) GelssError!usize {
    const F = Float(T);
    const mn = @min(m, n);
    try util.holds(a.len, m, n, lda);
    try util.holds(b.len, @max(m, n), nrhs, ldb);
    if (s.len < mn) return error.SliceTooShort;
    if (mn == 0) return 0;
    const am: Mat(T) = .{ .ptr = a.ptr, .rows = m, .cols = n, .ld = lda };
    const bm: Mat(T) = .{ .ptr = b.ptr, .rows = @max(m, n), .cols = nrhs, .ld = ldb };
    const sv = s[0..mn];
    const sc = try lsq_scale.scaleIn(al, T, am, bm.sub(0, 0, m, nrhs)) orelse {
        bm.fill(ops.zero(T));
        @memset(sv, 0);
        return 0;
    };
    const u = try ops.alloc(T, al, m, mn);
    defer ops.free(T, al, u);
    const v = try ops.alloc(T, al, n, mn);
    defer ops.free(T, al, v);
    try svd.svd(al, io, T, am, sv, u, v, if (driver == .gelss) std.math.maxInt(usize) else svd.dc_threshold, null);
    const thr: F = switch (driver) {
        // `dgelss.f`: `EPS = DLAMCH('P')`, `SFMIN = DLAMCH('S')`
        .gelss => @max((if (rcond < 0) std.math.floatEps(F) else rcond) * sv[0], std.math.floatMin(F)),
        // `dlalsd.f`: `EPS = DLAMCH('Epsilon')`
        .gelsd => (if (rcond <= 0 or rcond >= 1) std.math.floatEps(F) / 2 else rcond) * sv[0],
    };
    var rank: usize = 0;
    while (rank < mn and sv[rank] > thr) rank += 1;
    // X = V Σ⁺ Uᴴ B over the rank's leading values
    if (nrhs > 0) {
        const c = try ops.alloc(T, al, rank, nrhs);
        defer ops.free(T, al, c);
        try blas_calls.gemm(al, io, T, .c, .n, ops.one(T), u.sub(0, 0, m, rank), bm.sub(0, 0, m, nrhs), ops.zero(T), c);
        for (0..nrhs) |j| for (0..rank) |i| {
            c.at(i, j).* = ops.mulReal(T, c.get(i, j), 1 / sv[i]);
        };
        const x = bm.sub(0, 0, n, nrhs);
        if (rank == 0) x.fill(ops.zero(T)) else try blas_calls.gemm(al, io, T, .n, .n, ops.one(T), v.sub(0, 0, n, rank), c, ops.zero(T), x);
        lsq_scale.unscale(T, sc, x);
    }
    if (sc.iascl != 0) auxiliary.lascl(F, false, sc.iascl, sc.anrm, .{ .ptr = sv.ptr, .rows = mn, .cols = 1, .ld = mn });
    return rank;
}
