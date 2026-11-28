//! `xSTEBZ`: selected eigenvalues of a symmetric tridiagonal matrix by
//! bisection, Reference-LAPACK 3.12.0 (`lapack-netlib/SRC/dstebz.f`,
//! OpenBLAS v0.3.34-794-gdbf2c24c9; modified BSD), over `laebz.zig`. Indices are
//! 0-based. Serial.

const std = @import("std");
const strided = @import("../blas/strided.zig");
const interface = @import("interface.zig");
const ops = @import("ops.zig");
const bisection = @import("laebz.zig");

const Allocator = std.mem.Allocator;
pub const SpectrumRange = interface.SpectrumRange;
pub const SpectrumOrder = interface.SpectrumOrder;
const laebz = bisection.laebz;

/// `dstebz.f`'s `FUDGE`: how far the Gershgorin bounds are widened, in
/// units of their rounding.
const fudge = 2.1;
/// `RELFAC`: the relative tolerance in units of `ulp`.
const relfac = 2;

/// What `stebz` found: `count` eigenvalues in `w`, `blocks` diagonal
/// blocks. `unconverged` is the reference's `INFO = 1` or 3: some
/// eigenvalue's interval did not narrow to the tolerance in the iterations
/// allowed, its `w` the midpoint. `incomplete` is `INFO = 2` or 3: for
/// `.indices`, eigenvalues so close together that bisection could not tell
/// them apart, so `count` is not the number asked for.
pub const BisectedEigenvalues = struct { count: usize, blocks: usize, unconverged: bool = false, incomplete: bool = false };

/// What `stebz` returns: `EmptyInterval` for `.values` with `lo >= hi` (or
/// NaN), `IndexOutOfRange` for `.indices` not `first <= last < n`, a slice
/// too short, `IntervalTooSmall` for the reference's `INFO = 4` (the
/// widened Gershgorin interval missed eigenvalues: `T` is too close to
/// overflow), and the allocator's.
pub const StebzError = error{ EmptyInterval, IndexOutOfRange, IntervalTooSmall } || strided.SliceError || Allocator.Error;

/// `?STEBZ` (`dstebz.f`): the eigenvalues of the symmetric tridiagonal
/// `T` (diagonal `d`, `n` long; off-diagonal `e`, `n - 1`) that `range`
/// selects, each to within `abstol` (`abstol <= 0`: `ulp ‖T‖₁`; `2 safmin`
/// for the most accurate). `T` splits into blocks where `e(j)² <
/// ulp² |d(j)| |d(j + 1)|`: `isplit[b]` is one past block `b`'s last row
/// and `iblock[i]` the block of `w[i]`. `order` sorts `w` ascending or
/// leaves it by block, each block ascending. `w`, `iblock` and `isplit`
/// hold `n`. `T`'s entries must stay below `sqrt(overflow) underflow^(1/4)`
/// in magnitude, as the reference asks.
pub fn stebz(allocator: Allocator, comptime F: type, range: SpectrumRange(F), order: SpectrumOrder, d: []const F, e: []const F, abstol: F, w: []F, iblock: []usize, isplit: []usize) StebzError!BisectedEigenvalues {
    comptime ops.requireFloat(F, "stebz");
    comptime ops.requireReal(F, "stebz");
    const n = d.len;
    switch (range) {
        .all => {},
        .values => |v| if (!(v.lo < v.hi)) return error.EmptyInterval,
        .indices => |ix| if (n > 0 and !(ix.first <= ix.last and ix.last < n)) return error.IndexOutOfRange,
    }
    if (n == 0) return .{ .count = 0, .blocks = 0 };
    if (e.len < n - 1 or w.len < n or iblock.len < n or isplit.len < n) return error.SliceTooShort;
    var sel = range;
    if (sel == .indices and sel.indices.first == 0 and sel.indices.last == n - 1) sel = .all;
    const safemn = std.math.floatMin(F);
    // `DLAMCH('P')`
    const ulp = std.math.floatEps(F);
    const rtoli = ulp * relfac;

    if (n == 1) {
        isplit[0] = 1;
        if (sel == .values and (sel.values.lo >= d[0] or sel.values.hi < d[0])) return .{ .count = 0, .blocks = 1 };
        w[0] = d[0];
        iblock[0] = 0;
        return .{ .count = 1, .blocks = 1 };
    }

    // the squared off-diagonal, zero where `T` splits; `pivmin`
    const e2 = try allocator.alloc(F, n);
    defer allocator.free(e2);
    e2[n - 1] = 0;
    var nsplit: usize = 0;
    var pivmin: F = 1;
    for (1..n) |j| {
        const t1 = e[j - 1] * e[j - 1];
        const t2 = (@abs(d[j]) * ulp) * (@abs(d[j - 1]) * ulp);
        if (t1 < t2 or t1 == 0) {
            isplit[nsplit] = j;
            nsplit += 1;
            e2[j - 1] = 0;
        } else {
            e2[j - 1] = t1;
            pivmin = @max(pivmin, t1);
        }
    }
    isplit[nsplit] = n;
    nsplit += 1;
    pivmin *= safemn;

    // `MMAX` intervals for one block, as the reference's workspace holds
    const ab = try allocator.alloc([2]F, n);
    defer allocator.free(ab);
    const nab = try allocator.alloc([2]isize, n);
    defer allocator.free(nab);
    const c = try allocator.alloc(F, n);
    defer allocator.free(c);

    // `(wl, wu]`, and for `.indices` the inner ends `wlu`, `wul` and the
    // counts at `wl` and `wu`
    var wl: F = 0;
    var wu: F = 0;
    var wlu: F = 0;
    var wul: F = 0;
    switch (sel) {
        .all => {},
        .values => |v| {
            wl = v.lo;
            wu = v.hi;
        },
        .indices => |ix| {
            var gl = d[0];
            var gu = d[0];
            var t1: F = 0;
            for (0..n - 1) |j| {
                const t2 = @sqrt(e2[j]);
                gu = @max(gu, d[j] + t1 + t2);
                gl = @min(gl, d[j] - t1 - t2);
                t1 = t2;
            }
            gu = @max(gu, d[n - 1] + t1);
            gl = @min(gl, d[n - 1] - t1);
            const tnorm = @max(@abs(gl), @abs(gu));
            const nf: F = @floatFromInt(n);
            gl = gl - fudge * tnorm * ulp * nf - fudge * 2 * pivmin;
            gu = gu + fudge * tnorm * ulp * nf + fudge * pivmin;
            const atoli = if (abstol <= 0) ulp * tnorm else abstol;
            var fab = [2][2]F{ .{ gl, gu }, .{ gl, gu } };
            const top: isize = @intCast(n + 1);
            var fnab = [2][2]isize{ .{ -1, top }, .{ -1, top } };
            var fc = [2]F{ gl, gu };
            const last: isize = @intCast(ix.last + 1);
            var nval = [2]isize{ @intCast(ix.first), last };
            _ = laebz(F, .find, iterations(F, tnorm, pivmin), 2, atoli, rtoli, pivmin, d, e2, .{ .ab = &fab, .nab = &fnab, .c = &fc, .nval = &nval });
            // `laebz` may have swapped the two
            const lo: usize, const hi: usize = if (nval[1] == last) .{ 0, 1 } else .{ 1, 0 };
            wl = fab[lo][0];
            wlu = fab[lo][1];
            wu = fab[hi][1];
            wul = fab[hi][0];
            const nwl = fnab[lo][0];
            const nwu = fnab[hi][1];
            const ni: isize = @intCast(n);
            if (nwl < 0 or nwl >= ni or nwu < 1 or nwu > ni) return error.IntervalTooSmall;
        },
    }

    var m: usize = 0;
    var nwl: isize = 0;
    var nwu: isize = 0;
    var unconverged = false;
    var iend: usize = 0;
    for (0..nsplit) |jb| {
        const ibegin = iend;
        iend = isplit[jb];
        const in = iend - ibegin;
        if (in == 1) {
            const v = d[ibegin] - pivmin;
            if (sel == .all or wl >= v) nwl += 1;
            if (sel == .all or wu >= v) nwu += 1;
            if (sel == .all or (wl < v and wu >= v)) {
                w[m] = d[ibegin];
                iblock[m] = jb;
                m += 1;
            }
            continue;
        }
        var gl = d[ibegin];
        var gu = d[ibegin];
        var t1: F = 0;
        for (ibegin..iend - 1) |j| {
            const t2 = @abs(e[j]);
            gu = @max(gu, d[j] + t1 + t2);
            gl = @min(gl, d[j] - t1 - t2);
            t1 = t2;
        }
        gu = @max(gu, d[iend - 1] + t1);
        gl = @min(gl, d[iend - 1] - t1);
        const bnorm = @max(@abs(gl), @abs(gu));
        const widen = fudge * bnorm * ulp * @as(F, @floatFromInt(in)) + fudge * pivmin;
        gl -= widen;
        gu += widen;
        const atoli = if (abstol <= 0) ulp * @max(@abs(gl), @abs(gu)) else abstol;
        if (sel != .all) {
            if (gu <= wl) {
                nwl += @intCast(in);
                nwu += @intCast(in);
                continue;
            }
            gl = @max(gl, wl);
            gu = @min(gu, wu);
            if (gl >= gu) continue;
        }
        const db = d[ibegin..iend];
        const eb = e2[ibegin .. iend - 1];
        const iv: bisection.Intervals(F) = .{ .ab = ab[0..in], .nab = nab[0..in], .c = c[0..in] };
        ab[0] = .{ gl, gu };
        const im = laebz(F, .count, 0, 1, atoli, rtoli, pivmin, db, eb, iv).mout;
        nwl += nab[0][0];
        nwu += nab[0][1];
        const iwoff: isize = @as(isize, @intCast(m)) - nab[0][0];
        const out = laebz(F, .bisect, iterations(F, gu - gl, pivmin), 1, atoli, rtoli, pivmin, db, eb, iv);
        for (0..out.mout) |j| {
            const mid = 0.5 * (ab[j][0] + ab[j][1]);
            if (j >= out.mout - out.unconverged) unconverged = true;
            const from: usize = @intCast(nab[j][0] + iwoff);
            const to: usize = @intCast(nab[j][1] + iwoff);
            @memset(w[from..to], mid);
            @memset(iblock[from..to], jb);
        }
        m += im;
    }

    var incomplete = false;
    if (sel == .indices) {
        const ix = sel.indices;
        var idiscl: isize = @as(isize, @intCast(ix.first)) - nwl;
        var idiscu: isize = nwu - @as(isize, @intCast(ix.last + 1));
        if (idiscl > 0 or idiscu > 0) {
            var im: usize = 0;
            for (0..m) |je| {
                if (w[je] <= wlu and idiscl > 0) {
                    idiscl -= 1;
                } else if (w[je] >= wul and idiscu > 0) {
                    idiscu -= 1;
                } else {
                    w[im] = w[je];
                    iblock[im] = iblock[je];
                    im += 1;
                }
            }
            m = im;
        }
        if (idiscl > 0 or idiscu > 0) {
            // the reference's `IBLOCK = 0`: dropped
            const killed = std.math.maxInt(usize);
            if (idiscl > 0) {
                var wkill = wu;
                for (0..@as(usize, @intCast(idiscl))) |_| {
                    var iw: ?usize = null;
                    for (0..m) |je| if (iblock[je] != killed and (w[je] < wkill or iw == null)) {
                        iw = je;
                        wkill = w[je];
                    };
                    if (iw) |k| iblock[k] = killed;
                }
            }
            if (idiscu > 0) {
                var wkill = wl;
                for (0..@as(usize, @intCast(idiscu))) |_| {
                    var iw: ?usize = null;
                    for (0..m) |je| if (iblock[je] != killed and (w[je] > wkill or iw == null)) {
                        iw = je;
                        wkill = w[je];
                    };
                    if (iw) |k| iblock[k] = killed;
                }
            }
            var im: usize = 0;
            for (0..m) |je| if (iblock[je] != killed) {
                w[im] = w[je];
                iblock[im] = iblock[je];
                im += 1;
            };
            m = im;
        }
        incomplete = idiscl < 0 or idiscu < 0;
    }

    // `ORDER = 'E'`: a selection sort, as the reference's
    if (order == .ascending and nsplit > 1) for (0..m -| 1) |je| {
        var ie: ?usize = null;
        var t = w[je];
        for (je + 1..m) |j| if (w[j] < t) {
            ie = j;
            t = w[j];
        };
        if (ie) |k| {
            std.mem.swap(F, &w[k], &w[je]);
            std.mem.swap(usize, &iblock[k], &iblock[je]);
        }
    };
    return .{ .count = m, .blocks = nsplit, .unconverged = unconverged, .incomplete = incomplete };
}

/// `ITMAX`: `INT(LOG2((width + pivmin) / pivmin)) + 2`, the halvings that
/// take `width` down to `pivmin`.
fn iterations(comptime F: type, width: F, pivmin: F) usize {
    return @as(usize, @intFromFloat((@log(width + pivmin) - @log(pivmin)) / @log(@as(F, 2)))) + 2;
}
