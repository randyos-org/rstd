//! `?SYEVX`, `?HEEVX`: selected eigenvalues and eigenvectors of a symmetric
//! (Hermitian) matrix, Reference-LAPACK 3.12.0 (`lapack-netlib/SRC/dsyevx.f`,
//! `zheevx.f`, OpenBLAS v0.3.34-794-gdbf2c24c9; modified BSD), over `sytrd`,
//! `stebz`, `stein` and `ormtr`. The reference's `DSTERF` and `DSTEQR` fast
//! path is `tridiag_evd.zig`'s QR algorithm, as in `stebr.zig`, with `Q`
//! applied by `ormtr` where the reference forms it by `DORGTR`.

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const interface = @import("interface.zig");
const ops = @import("ops.zig");
const norm_structured = @import("norm_structured.zig");
const tevd = @import("tridiag_evd.zig");
const tridiagonalize = @import("sytrd.zig");
const bisect = @import("stebz.zig");
const inverse_iteration = @import("stein.zig");
const tridiagonal_q = @import("ormtr.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const isComplex = types.isComplex;
const Uplo = types.Uplo;
const VectorJob = interface.VectorJob;
const SpectrumRange = interface.SpectrumRange;

/// What `syevx` found: `count` eigenvalues in `w`, ascending, and with
/// `jobz = .compute` as many eigenvectors in `z`. `unconverged` vectors did
/// not converge (the reference's `INFO` from 1 to `n`), their columns in
/// `ifail`; `values_unconverged` and `incomplete` are `stebz`'s flags (`INFO
/// = n + 1` to `n + 3`), and with `incomplete` no vectors are computed, as
/// the reference.
pub const SelectedEigenpairs = struct { count: usize, unconverged: usize = 0, values_unconverged: bool = false, incomplete: bool = false };

/// What `syevx` returns: a reversed `uplo`; `stebz`'s range errors and
/// `IntervalTooSmall`; a slice too short (`w` and, for vectors, `ifail`
/// under `n`, `z` under its columns: `n`, or `last - first + 1` for
/// `.indices`) or a leading dimension too small; the allocator's and
/// cancellation.
pub const SyevxError = error{ EmptyInterval, IndexOutOfRange, IntervalTooSmall } || interface.TriangleError || Allocator.Error || Io.Cancelable;

/// `?SYEVX`, `?HEEVX` (`dsyevx.f`, `zheevx.f`): the eigenvalues of the
/// `n x n` symmetric (Hermitian) `A`, its `uplo` triangle read, that
/// `range` selects, each to within `abstol` (as `stebz`), and with `jobz`
/// their orthonormal eigenvectors. `A` is overwritten by `sytrd`'s
/// reduction. A matrix whose largest entry is outside `[sqrt(smlnum),
/// min(sqrt(bignum), safmin^(-1/4))]` is scaled into it first (`ISCALE`).
/// Unlike the reference, which swaps `IFAIL`'s entries with the columns
/// they do not index, `ifail` lists the unconverged columns where they end
/// up after sorting.
pub fn syevx(allocator: Allocator, io: Io, comptime T: type, jobz: VectorJob, range: SpectrumRange(Float(T)), uplo: Uplo, n: usize, a: []T, lda: usize, abstol: Float(T), w: []Float(T), z: []T, ldz: usize, ifail: []usize) SyevxError!SelectedEigenpairs {
    comptime ops.requireFloat(T, "syevx");
    const F = Float(T);
    const wantz = jobz == .compute;
    try interface.holdsTriangle(uplo, n, a.len, lda);
    switch (range) {
        .all => {},
        .values => |v| if (n > 0 and !(v.lo < v.hi)) return error.EmptyInterval,
        .indices => |ix| if (n > 0 and !(ix.first <= ix.last and ix.last < n)) return error.IndexOutOfRange,
    }
    if (w.len < n) return error.SliceTooShort;
    const most = switch (range) {
        .indices => |ix| if (n > 0) ix.last - ix.first + 1 else 0,
        else => n,
    };
    if (wantz) {
        try util.holds(z.len, n, most, ldz);
        if (ifail.len < n) return error.SliceTooShort;
    }
    if (n == 0) return .{ .count = 0 };
    if (n == 1) {
        const a00 = ops.re(T, a[0]);
        const found = switch (range) {
            .values => |v| v.lo < a00 and v.hi >= a00,
            else => true,
        };
        if (found) w[0] = a00;
        if (wantz) z[0] = ops.one(T);
        return .{ .count = @intFromBool(found) };
    }

    // `ISCALE`: into `[rmin, rmax]`
    const safmin = std.math.floatMin(F);
    // `DLAMCH('P')`
    const eps = std.math.floatEps(F);
    const smlnum = safmin / eps;
    const bignum = 1 / smlnum;
    const rmin = @sqrt(smlnum);
    const rmax = @min(@sqrt(bignum), 1 / @sqrt(@sqrt(safmin)));
    const anrm = (if (comptime isComplex(T)) norm_structured.lanhe else norm_structured.lansy)(allocator, io, T, .max_abs, uplo, n, a, lda) catch |e| switch (e) {
        error.UnsupportedUplo, error.SliceTooShort, error.LeadingDimensionTooSmall => unreachable,
        error.OutOfMemory => |x| return x,
    };
    var sigma: F = 1;
    if (anrm > 0 and anrm < rmin) {
        sigma = rmin / anrm;
    } else if (anrm > rmax) {
        sigma = rmax / anrm;
    }
    const iscale = sigma != 1 and sigma > 0;
    var abstll = abstol;
    var sel = range;
    if (iscale) {
        for (0..n) |j| {
            const col = if (uplo == .l) a[j + j * lda .. n + j * lda] else a[j * lda .. j + 1 + j * lda];
            for (col) |*x| x.* = ops.mulReal(T, x.*, sigma);
        }
        if (abstol > 0) abstll = abstol * sigma;
        if (sel == .values) sel = .{ .values = .{ .lo = sel.values.lo * sigma, .hi = sel.values.hi * sigma } };
    }

    // `d`, `e`, the reflections' `tau`; a copy of `e` for the QR algorithm
    const d = try allocator.alloc(F, 3 * n);
    defer allocator.free(d);
    const e = d[n..][0 .. n - 1];
    const ee = d[2 * n ..][0 .. n - 1];
    const tau = try allocator.alloc(T, n - 1);
    defer allocator.free(tau);
    try tridiagonalize.sytrd(allocator, io, T, uplo, n, a, lda, d[0..n], e, tau);
    const iblock = try allocator.alloc(usize, 2 * n);
    defer allocator.free(iblock);
    const isplit = iblock[n..];

    var out: SelectedEigenpairs = .{ .count = 0 };
    const every = sel == .all or (sel == .indices and sel.indices.first == 0 and sel.indices.last == n - 1);
    var fast = false;
    if (every and abstol <= 0) {
        @memcpy(w[0..n], d[0..n]);
        @memcpy(ee, e);
        if (!wantz) {
            fast = if (tevd.qrAlgorithm(F, w[0..n], ee, null)) true else |_| false;
        } else {
            const v = try ops.alloc(F, allocator, n, n);
            defer ops.free(F, allocator, v);
            fast = if (tevd.qrAlgorithm(F, w[0..n], ee, v)) true else |_| false;
            if (fast) {
                for (0..n) |j| for (0..n) |i| {
                    z[i + j * ldz] = ops.scalar(T, v.get(i, j));
                };
                try applyQ(allocator, io, T, uplo, n, n, a, lda, tau, z, ldz);
            }
        }
        if (fast) out.count = n;
    }
    if (!fast) {
        const found = bisect.stebz(allocator, F, sel, if (wantz) .by_block else .ascending, d[0..n], e, abstll, w[0..n], iblock[0..n], isplit) catch |err| switch (err) {
            // checked above, or sized above
            error.EmptyInterval, error.IndexOutOfRange, error.SliceTooShort => unreachable,
            error.IntervalTooSmall, error.OutOfMemory => |x| return x,
        };
        out = .{ .count = found.count, .values_unconverged = found.unconverged, .incomplete = found.incomplete };
        if (wantz and !found.incomplete) {
            out.unconverged = inverse_iteration.stein(allocator, io, T, d[0..n], e, w[0..found.count], iblock[0..found.count], isplit, z, ldz, ifail) catch |err| switch (err) {
                // `stebz`'s order, sized above
                error.TooManyEigenvalues, error.BlocksOutOfOrder, error.EigenvaluesOutOfOrder, error.SplitOutOfOrder, error.SliceTooShort, error.LeadingDimensionTooSmall => unreachable,
                error.OutOfMemory => |x| return x,
            };
            try applyQ(allocator, io, T, uplo, n, found.count, a, lda, tau, z, ldz);
        }
    }
    const m = out.count;
    if (iscale) for (w[0..m]) |*x| {
        x.* /= sigma;
    };
    // ascending, the vectors with them (`ORDER = 'B'` left them by block)
    if (wantz and !fast and !out.incomplete) {
        const failed = try allocator.alloc(bool, m);
        defer allocator.free(failed);
        @memset(failed, false);
        for (ifail[0..out.unconverged]) |j| failed[j] = true;
        for (0..m -| 1) |j| {
            var best: ?usize = null;
            var t = w[j];
            for (j + 1..m) |jj| if (w[jj] < t) {
                best = jj;
                t = w[jj];
            };
            if (best) |i| {
                std.mem.swap(F, &w[i], &w[j]);
                std.mem.swap(bool, &failed[i], &failed[j]);
                for (0..n) |r| std.mem.swap(T, &z[r + i * ldz], &z[r + j * ldz]);
            }
        }
        var k: usize = 0;
        for (failed, 0..) |f, j| if (f) {
            ifail[k] = j;
            k += 1;
        };
    }
    return out;
}

/// `Z := Q Z` for `sytrd`'s `Q` and the `n x m` `Z`.
fn applyQ(al: Allocator, io: Io, comptime T: type, uplo: Uplo, n: usize, m: usize, a: []const T, lda: usize, tau: []const T, z: []T, ldz: usize) (Allocator.Error || Io.Cancelable)!void {
    if (m == 0) return;
    tridiagonal_q.ormtr(al, io, T, .l, uplo, .n, n, m, a, lda, tau, z, ldz) catch |err| switch (err) {
        // `.n` and the shapes checked above
        error.UnsupportedUplo, error.UnsupportedTranspose, error.TooManyReflections, error.SliceTooShort, error.LeadingDimensionTooSmall => unreachable,
        else => |x| return x,
    };
}
