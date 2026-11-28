//! `lange` for structured matrices: `?LANSY`, `?LANHE` (a symmetric or
//! Hermitian matrix from one triangle) and `?LANTR` (a trapezoid),
//! Reference-LAPACK 3.12.0 (`lapack-netlib/SRC/dlansy.f`, `zlansy.f`,
//! `zlanhe.f`, `dlantr.f`, `zlantr.f`, OpenBLAS v0.3.34-317; modified BSD).

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const interface = @import("interface.zig");
const ops = @import("ops.zig");
const norm = @import("norm.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;
const Uplo = types.Uplo;
const Diag = types.Diag;
const NormType = norm.NormType;
const lassq = norm.lassq;

/// What `lansy`, `lanhe` and `lantr` return: a `uplo` that is neither
/// triangle, a slice too short or a leading dimension too small, and the
/// allocator's for `.one` and `.infinity`'s workspace.
pub const LansyError = interface.TriangleError || Allocator.Error;
pub const LantrError = LansyError;

/// `VALUE := SUM` where the reference's `VALUE .LT. SUM .OR. DISNAN(SUM)`.
fn keep(comptime F: type, value: *F, sum: F) void {
    if (value.* < sum or std.math.isNan(sum)) value.* = sum;
}

/// `?LANSY` (`dlansy.f`, `zlansy.f`): `norm` of the symmetric (complex:
/// complex symmetric) `n x n` `A` from its `uplo` triangle; 0 when empty. A
/// NaN entry gives NaN; `.one` and `.infinity` agree.
pub fn lansy(allocator: Allocator, io: Io, comptime T: type, norm_type: NormType, uplo: Uplo, n: usize, a: []const T, lda: usize) LansyError!Float(T) {
    comptime ops.requireFloat(T, "lansy");
    _ = io;
    try interface.holdsTriangle(uplo, n, a.len, lda);
    return symmetric(allocator, T, false, norm_type, uplo == .u, n, a, lda);
}

/// `?LANHE` (`zlanhe.f`): `lansy` for the Hermitian `A`, its diagonal's
/// real parts alone read. Complex only.
pub fn lanhe(allocator: Allocator, io: Io, comptime T: type, norm_type: NormType, uplo: Uplo, n: usize, a: []const T, lda: usize) LansyError!Float(T) {
    comptime ops.requireComplex(T, "lanhe");
    _ = io;
    try interface.holdsTriangle(uplo, n, a.len, lda);
    return symmetric(allocator, T, true, norm_type, uplo == .u, n, a, lda);
}

fn symmetric(al: Allocator, comptime T: type, comptime hermitian: bool, norm_type: NormType, upper: bool, n: usize, a: []const T, lda: usize) Allocator.Error!Float(T) {
    const F = Float(T);
    if (n == 0) return 0;
    const at = struct {
        fn f(x: []const T, ld: usize, i: usize, j: usize) F {
            return ops.abs(T, x[i + j * ld]);
        }
    }.f;
    // the diagonal's modulus, its real part's for a Hermitian `A`
    const diag = struct {
        fn f(x: []const T, ld: usize, j: usize) F {
            return if (hermitian) @abs(ops.re(T, x[j + j * ld])) else ops.abs(T, x[j + j * ld]);
        }
    }.f;
    var value: F = 0;
    switch (norm_type) {
        .max_abs => for (0..n) |j| {
            if (upper) {
                for (0..j) |i| keep(F, &value, at(a, lda, i, j));
                keep(F, &value, diag(a, lda, j));
            } else {
                keep(F, &value, diag(a, lda, j));
                for (j + 1..n) |i| keep(F, &value, at(a, lda, i, j));
            }
        },
        .one, .infinity => {
            const work = try al.alloc(F, n);
            defer al.free(work);
            if (upper) {
                for (0..n) |j| {
                    var sum: F = 0;
                    for (0..j) |i| {
                        const absa = at(a, lda, i, j);
                        sum += absa;
                        work[i] += absa;
                    }
                    work[j] = sum + diag(a, lda, j);
                }
                for (work) |w| keep(F, &value, w);
            } else {
                @memset(work, 0);
                for (0..n) |j| {
                    var sum = work[j] + diag(a, lda, j);
                    for (j + 1..n) |i| {
                        const absa = at(a, lda, i, j);
                        sum += absa;
                        work[i] += absa;
                    }
                    keep(F, &value, sum);
                }
            }
        },
        .frobenius => {
            var scale: F = 0;
            var sum: F = 1;
            if (upper) {
                for (1..n) |j| lassq(T, .of(a[j * lda ..][0..j]), &scale, &sum);
            } else {
                for (0..n - 1) |j| lassq(T, .of(a[j + 1 + j * lda ..][0 .. n - 1 - j]), &scale, &sum);
            }
            sum *= 2;
            if (hermitian) {
                // `zlanhe.f`'s own update on the diagonal's real parts
                for (0..n) |i| {
                    const re = ops.re(T, a[i + i * lda]);
                    if (re == 0) continue;
                    const absa = @abs(re);
                    if (scale < absa) {
                        sum = 1 + sum * (scale / absa) * (scale / absa);
                        scale = absa;
                    } else sum += (absa / scale) * (absa / scale);
                }
            } else lassq(T, .{ .data = a[0 .. (n - 1) * (lda + 1) + 1], .n = n, .stride = @intCast(lda + 1) }, &scale, &sum);
            value = scale * @sqrt(sum);
        },
    }
    return value;
}

/// `?LANTR` (`dlantr.f`, `zlantr.f`): `norm` of the `m x n` upper or lower
/// trapezoid in `a`'s `uplo` part, its diagonal 1 for `diag` `.u` (not
/// read); 0 when empty. A NaN entry gives NaN.
pub fn lantr(allocator: Allocator, io: Io, comptime T: type, norm_type: NormType, uplo: Uplo, diag: Diag, m: usize, n: usize, a: []const T, lda: usize) LantrError!Float(T) {
    comptime ops.requireFloat(T, "lantr");
    _ = io;
    if (uplo.isRev()) return error.UnsupportedUplo;
    try util.holds(a.len, m, n, lda);
    const F = Float(T);
    if (@min(m, n) == 0) return 0;
    const upper = uplo == .u;
    const unit = diag == .u;
    const at = struct {
        fn f(x: []const T, ld: usize, i: usize, j: usize) F {
            return ops.abs(T, x[i + j * ld]);
        }
    }.f;
    // rows `lo..hi` of column `j` inside the trapezoid, the diagonal left out when unit
    const Rows = struct { lo: usize, hi: usize };
    const rows = struct {
        fn f(up: bool, u: bool, mm: usize, j: usize) Rows {
            if (up) return .{ .lo = 0, .hi = @min(mm, if (u) j else j + 1) };
            return .{ .lo = @min(mm, if (u) j + 1 else j), .hi = mm };
        }
    }.f;
    const cols = if (upper) n else @min(m, n);
    var value: F = 0;
    switch (norm_type) {
        .max_abs => {
            value = if (unit) 1 else 0;
            for (0..cols) |j| {
                const r = rows(upper, unit, m, j);
                for (r.lo..r.hi) |i| keep(F, &value, at(a, lda, i, j));
            }
        },
        .one => for (0..cols) |j| {
            const r = rows(upper, unit, m, j);
            // a unit diagonal counts where the column holds it
            var sum: F = if (unit and (if (upper) j < m else true)) 1 else 0;
            for (r.lo..r.hi) |i| sum += at(a, lda, i, j);
            keep(F, &value, sum);
        },
        .infinity => {
            const work = try allocator.alloc(F, m);
            defer allocator.free(work);
            @memset(work, 0);
            if (unit) @memset(work[0..@min(m, n)], 1);
            for (0..cols) |j| {
                const r = rows(upper, unit, m, j);
                for (r.lo..r.hi) |i| work[i] += at(a, lda, i, j);
            }
            for (work) |w| keep(F, &value, w);
        },
        .frobenius => {
            var scale: F = if (unit) 1 else 0;
            var sum: F = if (unit) @floatFromInt(@min(m, n)) else 1;
            for (0..cols) |j| {
                const r = rows(upper, unit, m, j);
                if (r.hi > r.lo) lassq(T, .of(a[r.lo + j * lda ..][0 .. r.hi - r.lo]), &scale, &sum);
            }
            value = scale * @sqrt(sum);
        },
    }
    return value;
}
