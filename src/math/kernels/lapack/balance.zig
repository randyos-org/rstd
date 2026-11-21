//! Balancing, `gebal` and `gebak`, from Reference-LAPACK 3.12.0 as OpenBLAS
//! v0.3.34-317 ships it (`lapack-netlib/SRC/`, modified BSD): `xGEBAL` and
//! `xGEBAK`, their block `lo..hi` and `scale`'s rows 0-based.

const std = @import("std");
const blas = @import("../blas/__root__.zig");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const interface = @import("interface.zig");
const ops = @import("ops.zig");
const blas_calls = @import("blas_calls.zig");

const Io = std.Io;
const Float = types.Float;
const Side = types.Side;
const Error = ops.Error;
const Mat = ops.Mat;

/// LAPACK's `JOB`: nothing (`N`), permute only (`P`), scale only (`S`), or
/// both (`B`).
pub const Balancing = enum { none, permute, scale, both };

/// `SCLFAC`: scalings are powers of 2.
const sclfac = 2;
/// `FACTOR`: a scaling must cut `c + r` below this fraction of it.
const factor = 0.95;

/// The rows and columns `lo..hi` (0-based, `hi` past the end) that
/// `gebal` leaves to reduce; outside them `A` is upper triangular.
pub const BalancedRange = struct { lo: usize, hi: usize };

/// What `gebal` returns: a slice too short (`a`, or `scale` under `n`) or a
/// leading dimension too small, `NanInA` for a NaN met while scaling (the
/// reference's `INFO = -3`, `A` then part balanced).
pub const GebalError = interface.ShapeError || interface.NanError || Io.Cancelable;

/// `?GEBAL` (`dgebal.f`, `zgebal.f`): balances the `n x n` `A` in place,
/// `D⁻¹ Pᵀ A P D`, and returns the block `lo..hi` outside which `A` is then
/// upper triangular (`0..0` for `n = 0`). `scale` (`n`) gets, outside the
/// block, the row (0-based) each row was swapped with (in the reference's
/// order: from `n - 1` down to `hi`, then up from 0), and inside it `D`'s
/// entries.
pub fn gebal(io: Io, comptime T: type, job: Balancing, n: usize, a: []T, lda: usize, scale: []Float(T)) GebalError!BalancedRange {
    comptime ops.requireFloat(T, "gebal");
    try util.holds(a.len, n, n, lda);
    if (scale.len < n) return error.SliceTooShort;
    if (n == 0) return .{ .lo = 0, .hi = 0 };
    return balanceJob(io, T, job, n, a.ptr, lda, scale.ptr) catch |e| switch (e) {
        error.NanInA, error.Canceled => |x| x,
        error.OutOfMemory => unreachable,
    };
}

/// `gebal` after its checks.
pub fn balanceJob(io: Io, comptime T: type, job: Balancing, n: usize, a: [*]T, lda: usize, scale: [*]Float(T)) (Error || interface.NanError)!BalancedRange {
    const F = Float(T);
    if (job == .none) {
        for (0..n) |i| scale[i] = 1;
        return .{ .lo = 0, .hi = n };
    }
    // `k`, `l`: the 1-based active window
    var k: usize = 1;
    var l: usize = n;
    if (job != .scale) {
        // rows isolating an eigenvalue go to the bottom
        var noconv = true;
        while (noconv) {
            noconv = false;
            var i = l;
            while (i >= 1) : (i -= 1) {
                const can = for (1..l + 1) |j| {
                    if (i != j and !ops.isZero(T, at(a, lda, i, j))) break false;
                } else true;
                if (can) {
                    scale[l - 1] = @floatFromInt(i - 1);
                    if (i != l) {
                        try blas_calls.swap(io, T, l, a + (i - 1) * lda, 1, a + (l - 1) * lda, 1);
                        try blas_calls.swap(io, T, n - k + 1, a + (i - 1) + (k - 1) * lda, lda, a + (l - 1) + (k - 1) * lda, lda);
                    }
                    noconv = true;
                    // the reference's `SCALE(1) = 1` is both the row it
                    // swapped with and `D(1)`; 0-based it is `D(0)`
                    if (l == 1) {
                        scale[0] = 1;
                        return .{ .lo = 0, .hi = 1 };
                    }
                    l -= 1;
                }
                if (i == 1) break;
            }
        }
        // columns isolating one go to the left
        noconv = true;
        while (noconv) {
            noconv = false;
            const last = l;
            var j = k;
            while (j <= last) : (j += 1) {
                const can = for (k..l + 1) |i| {
                    if (i != j and !ops.isZero(T, at(a, lda, i, j))) break false;
                } else true;
                if (can) {
                    scale[k - 1] = @floatFromInt(j - 1);
                    if (j != k) {
                        try blas_calls.swap(io, T, l, a + (j - 1) * lda, 1, a + (k - 1) * lda, 1);
                        try blas_calls.swap(io, T, n - k + 1, a + (j - 1) + (k - 1) * lda, lda, a + (k - 1) + (k - 1) * lda, lda);
                    }
                    noconv = true;
                    k += 1;
                }
            }
        }
    }
    for (k - 1..l) |i| scale[i] = 1;
    if (job == .permute) return .{ .lo = k - 1, .hi = l };
    // scale toward equal row and column norms
    const sfmin1 = std.math.floatMin(F) / std.math.floatEps(F);
    const sfmax1 = 1 / sfmin1;
    const sfmin2 = sfmin1 * sclfac;
    const sfmax2 = 1 / sfmin2;
    var noconv = true;
    while (noconv) {
        noconv = false;
        for (k..l + 1) |i| {
            const w = l - k + 1;
            var c = blas_calls.nrm2(io, T, w, a + (k - 1) + (i - 1) * lda, 1);
            var r = blas_calls.nrm2(io, T, w, a + (i - 1) + (k - 1) * lda, lda);
            const ica = blas_calls.iamax(T, l, a + (i - 1) * lda, 1);
            var ca = ops.abs(T, a[ica + (i - 1) * lda]);
            const ira = blas_calls.iamax(T, n - k + 1, a + (i - 1) + (k - 1) * lda, lda);
            var ra = ops.abs(T, a[(i - 1) + (ira + k - 1) * lda]);
            if (c == 0 or r == 0) continue;
            if (std.math.isNan(c + ca + r + ra)) return error.NanInA;
            var g = r / sclfac;
            var f: F = 1;
            const s = c + r;
            while (c < g and @max(f, c, ca) < sfmax2 and @min(r, g, ra) > sfmin2) {
                f *= sclfac;
                c *= sclfac;
                ca *= sclfac;
                r /= sclfac;
                g /= sclfac;
                ra /= sclfac;
            }
            g = c / sclfac;
            while (g >= r and @max(r, ra) < sfmax2 and @min(f, c, g, ca) > sfmin2) {
                f /= sclfac;
                c /= sclfac;
                g /= sclfac;
                ca /= sclfac;
                r *= sclfac;
                ra *= sclfac;
            }
            if (c + r >= factor * s) continue;
            if (f < 1 and scale[i - 1] < 1 and f * scale[i - 1] <= sfmin1) continue;
            if (f > 1 and scale[i - 1] > 1 and scale[i - 1] >= sfmax1 / f) continue;
            g = 1 / f;
            scale[i - 1] *= f;
            noconv = true;
            try blas_calls.scal(io, T, n - k + 1, g, a + (i - 1) + (k - 1) * lda, lda);
            try blas_calls.scal(io, T, l, f, a + (i - 1) * lda, 1);
        }
    }
    return .{ .lo = k - 1, .hi = l };
}

/// `A(i, j)`, 1-based.
fn at(a: anytype, lda: usize, i: usize, j: usize) @typeInfo(@TypeOf(a)).pointer.child {
    return a[(i - 1) + (j - 1) * lda];
}

/// What `gebak` returns: `BlockOutOfRange` unless `lo <= hi <= n`, a slice
/// too short (`v`, or `scale` under `n`) or a leading dimension too small,
/// `PivotOutOfRange` for a `scale` entry outside the block (with `job`
/// `.permute` or `.both`) that is not a row of `V`.
pub const GebakError = interface.BlockError || interface.ShapeError || interface.PivotError || Io.Cancelable;

/// `?GEBAK` (`dgebak.f`, `zgebak.f`): `V := P D V` for right eigenvectors
/// (`side` `r`) or `V := P D⁻¹ V` for left ones (`l`), the `n x m` `V`'s
/// rows, from `gebal`'s `job`, block `lo..hi` and `scale`.
pub fn gebak(io: Io, comptime T: type, job: Balancing, side: Side, n: usize, lo: usize, hi: usize, scale: []const Float(T), m: usize, v: []T, ldv: usize) GebakError!void {
    comptime ops.requireFloat(T, "gebak");
    try interface.blockIn(lo, hi, n);
    try util.holds(v.len, n, m, ldv);
    if (scale.len < n) return error.SliceTooShort;
    if (job == .permute or job == .both) for (0..n) |i| if (i < lo or i >= hi) {
        const k = scale[i];
        if (!(k >= 0 and k < @as(Float(T), @floatFromInt(n)) and k == @trunc(k))) return error.PivotOutOfRange;
    };
    if (n == 0 or m == 0 or job == .none) return;
    backJob(io, T, job, side, n, lo, hi, scale.ptr, m, v.ptr, ldv) catch |e| return switch (e) {
        error.OutOfMemory => unreachable,
        error.Canceled => error.Canceled,
    };
}

/// `gebak` after its checks.
pub fn backJob(io: Io, comptime T: type, job: Balancing, side: Side, n: usize, lo: usize, hi: usize, scale: [*]const Float(T), m: usize, v: [*]T, ldv: usize) Error!void {
    // `dgebak.f`: a block of one is not scaled
    if (hi > lo + 1 and (job == .scale or job == .both)) {
        for (lo..hi) |i| {
            const s = if (side == .r) scale[i] else 1 / scale[i];
            try blas_calls.scal(io, T, m, s, v + i, ldv);
        }
    }
    if (job == .permute or job == .both) {
        // from `lo - 1` down to 0, then from `hi` up
        for (0..n) |ii| {
            var i = ii;
            if (i >= lo and i < hi) continue;
            if (i < lo) i = lo - 1 - ii;
            const k: usize = @intFromFloat(scale[i]);
            if (k == i) continue;
            try blas_calls.swap(io, T, m, v + i, ldv, v + k, ldv);
        }
    }
}

// ---- geev's job B -------------------------------------------------------------------

/// `xGEBAL` with job `B` on the `n x n` `a`, which holds no NaN; `scale`
/// (`n`) as `gebal`'s.
pub fn balance(io: Io, comptime T: type, a: Mat(T), scale: []Float(T)) Error!BalancedRange {
    return balanceJob(io, T, .both, a.rows, a.ptr, a.ld, scale.ptr) catch |e| switch (e) {
        error.NanInA => unreachable,
        error.OutOfMemory, error.Canceled => |x| return x,
    };
}

/// `xGEBAK` with job `B`: `v`'s rows (the balanced matrix's right vectors,
/// or with `left` its left ones) back to `a`'s.
pub fn back(io: Io, comptime T: type, block: BalancedRange, scale: []const Float(T), left: bool, v: Mat(T)) Error!void {
    if (v.rows == 0 or v.cols == 0) return;
    try backJob(io, T, .both, if (left) .l else .r, v.rows, block.lo, block.hi, scale.ptr, v.cols, v.ptr, v.ld);
}
