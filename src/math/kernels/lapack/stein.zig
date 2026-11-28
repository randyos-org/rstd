//! `xSTEIN`: eigenvectors of a symmetric tridiagonal matrix by inverse
//! iteration, Reference-LAPACK 3.12.0 (`lapack-netlib/SRC/dstein.f`,
//! `zstein.f`, OpenBLAS v0.3.34-794-gdbf2c24c9; modified BSD), over
//! `lagtf.zig` and `laruv.zig`. Indices are 0-based. Serial.

const std = @import("std");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const interface = @import("interface.zig");
const ops = @import("ops.zig");
const blas_calls = @import("blas_calls.zig");
const tridiagonal = @import("lagtf.zig");
const random = @import("laruv.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Float = types.Float;

/// `MAXITS`: the most inverse iterations per eigenvector.
const maxits = 5;
/// `EXTRA`: iterations past the first whose growth passes the test.
const extra = 2;

/// What `stein` returns: more eigenvalues than rows (`TooManyEigenvalues`,
/// the reference's `INFO = -4`), `iblock` decreasing (`BlocksOutOfOrder`,
/// -6) or a block's eigenvalues not ascending (`EigenvaluesOutOfOrder`,
/// -5), `isplit` not increasing to at most `n` (`SplitOutOfOrder`), a slice
/// too short or a leading dimension too small, and the allocator's.
pub const SteinError = error{ TooManyEigenvalues, BlocksOutOfOrder, EigenvaluesOutOfOrder, SplitOutOfOrder } || interface.ShapeError || Allocator.Error;

/// `?STEIN` (`dstein.f`, `zstein.f`): for each eigenvalue `w[j]` of the
/// symmetric tridiagonal `T` (diagonal `d`, `n`; off-diagonal `e`, `n - 1`),
/// in block `iblock[j]` of `stebz`'s `isplit`, the eigenvector of its block
/// by inverse iteration from the reference's random start, into column `j`
/// of `z` (`n x m`, `m = w.len`; real in a complex `z`), zero outside the
/// block and unit in the 2-norm, its largest entry positive. Vectors of
/// eigenvalues closer than `1e-3 ‖T_b‖₁` are reorthogonalized against each
/// other. Returns how many did not converge in `maxits` iterations (the
/// reference's `INFO`), their columns in `ifail` (`m` long) in order; their
/// `z` columns hold the last iterate.
pub fn stein(allocator: Allocator, io: Io, comptime T: type, d: []const Float(T), e: []const Float(T), w: []const Float(T), iblock: []const usize, isplit: []const usize, z: []T, ldz: usize, ifail: []usize) SteinError!usize {
    comptime ops.requireFloat(T, "stein");
    const F = Float(T);
    const n = d.len;
    const m = w.len;
    if (m > n) return error.TooManyEigenvalues;
    if ((n > 0 and e.len < n - 1) or iblock.len < m or ifail.len < m) return error.SliceTooShort;
    try util.holds(z.len, n, m, ldz);
    for (1..@max(m, 1)) |j| {
        if (iblock[j] < iblock[j - 1]) return error.BlocksOutOfOrder;
        if (iblock[j] == iblock[j - 1] and w[j] < w[j - 1]) return error.EigenvaluesOutOfOrder;
    }
    if (n == 0 or m == 0) return 0;
    const nblocks = iblock[m - 1] + 1;
    if (isplit.len < nblocks) return error.SliceTooShort;
    for (0..nblocks) |b| {
        if (isplit[b] > n or isplit[b] <= (if (b > 0) isplit[b - 1] else 0)) return error.SplitOutOfOrder;
    }
    if (n == 1) {
        z[0] = ops.one(T);
        return 0;
    }
    // `DLAMCH('P')`
    const eps = std.math.floatEps(F);
    var seed: [4]u32 = @splat(1);
    // the iterate, then `lagtf`'s factor
    const work = try allocator.alloc(F, 5 * n);
    defer allocator.free(work);
    const in = try allocator.alloc(usize, n);
    defer allocator.free(in);

    var failed: usize = 0;
    var j1: usize = 0;
    for (0..nblocks) |nblk| {
        const b1 = if (nblk == 0) 0 else isplit[nblk - 1];
        const bn = isplit[nblk];
        const blksiz = bn - b1;
        var gpind = j1;
        var ortol: F = 0;
        var dtpcrt: F = 0;
        var onenrm: F = 0;
        if (blksiz > 1) {
            onenrm = @max(@abs(d[b1]) + @abs(e[b1]), @abs(d[bn - 1]) + @abs(e[bn - 2]));
            for (b1 + 1..bn - 1) |i| onenrm = @max(onenrm, @abs(d[i]) + @abs(e[i - 1]) + @abs(e[i]));
            // `ODM3`, `ODM1`
            ortol = 1e-3 * onenrm;
            dtpcrt = @sqrt(0.1 / @as(F, @floatFromInt(blksiz)));
        }
        const v = work[0..blksiz];
        const f: tridiagonal.Factor(F) = .{
            .a = work[n..][0..blksiz],
            .b = work[2 * n ..][0..blksiz -| 1],
            .c = work[3 * n ..][0..blksiz -| 1],
            .d = work[4 * n ..][0..blksiz -| 2],
            .in = in[0..blksiz],
        };
        var jblk: usize = 0;
        var xjm: F = 0;
        var j = j1;
        while (j < m and iblock[j] == nblk) : (j += 1) {
            jblk += 1;
            var xj = w[j];
            if (blksiz == 1) {
                v[0] = 1;
            } else {
                if (jblk > 1) {
                    // `PERTOL`: eigenvalues ten ulps apart at least
                    const pertol = 10 * @abs(eps * xj);
                    if (xj - xjm < pertol) xj = xjm + pertol;
                }
                random.larnvUniform(F, &seed, v);
                @memcpy(f.a, d[b1..bn]);
                @memcpy(f.b, e[b1 .. bn - 1]);
                @memcpy(f.c, e[b1 .. bn - 1]);
                var tol: F = 0;
                tridiagonal.lagtf(F, xj, tol, f);
                var its: usize = 0;
                var nrmchk: usize = 0;
                const converged = while (true) {
                    its += 1;
                    if (its > maxits) break false;
                    const scl = @as(F, @floatFromInt(blksiz)) * onenrm * @max(eps, @abs(f.a[blksiz - 1])) / @abs(v[iamax(F, v)]);
                    for (v) |*x| x.* *= scl;
                    tridiagonal.lagts(F, f, v, &tol);
                    // against the cluster's earlier vectors
                    if (jblk > 1) {
                        if (@abs(xj - xjm) > ortol) gpind = j;
                        for (gpind..j) |i| {
                            const col = z[b1 + i * ldz ..][0..blksiz];
                            var ztr: F = 0;
                            for (v, col) |x, y| ztr += x * ops.re(T, y);
                            for (v, col) |*x, y| x.* -= ztr * ops.re(T, y);
                        }
                    }
                    if (@abs(v[iamax(F, v)]) < dtpcrt) continue;
                    nrmchk += 1;
                    if (nrmchk < extra + 1) continue;
                    break true;
                };
                if (!converged) {
                    ifail[failed] = j;
                    failed += 1;
                }
                var scl = 1 / blas_calls.nrm2(io, F, blksiz, v.ptr, 1);
                if (v[iamax(F, v)] < 0) scl = -scl;
                for (v) |*x| x.* *= scl;
            }
            const col = z[j * ldz ..][0..n];
            @memset(col, ops.zero(T));
            for (col[b1..bn], v) |*o, x| o.* = ops.scalar(T, x);
            xjm = xj;
        }
        j1 = j;
    }
    return failed;
}

/// `IDAMAX`: the first entry of largest magnitude.
fn iamax(comptime F: type, v: []const F) usize {
    var k: usize = 0;
    for (v, 0..) |x, i| if (@abs(x) > @abs(v[k])) {
        k = i;
    };
    return k;
}
