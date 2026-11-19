//! The symmetric tridiagonal's eigenvalues by divide and conquer carrying
//! only each block's first and last eigenvector rows: Reference-LAPACK's
//! `xSTEBR` (`lapack-netlib/SRC/dstebr.f`, `dlaed0_br.f`,
//! `dlaed0_merge_pair.f`, `dlaed7_br.f`, `dsteqr_boundary_br.f`;
//! Reference-LAPACK PR 1424, by Ruiyi Zhan and Shaoshuai Zhang, OpenBLAS
//! v0.3.34-794-gdbf2c24c9; modified BSD). Deflation in `laed8_br.zig`, the
//! secular solve in `laed9_br.zig` and `laed4.zig`. The `s` and `d`
//! routines are the same; indices are 0-based. Serial.
//!
//! Reused rather than ported: `xSTERF` (small matrices and blocks) is
//! `tridiag_evd.zig`'s QR iteration without vectors, and the leaves' `xSTEQR`
//! (`COMPZ = 'I'`) the same with them; `xLANST('M')`, `xLASCL` and `xLAMRG`
//! are `auxiliary.zig`'s, `xLASRT('I')` a sort.

const std = @import("std");
const ops = @import("ops.zig");
const auxiliary = @import("auxiliary.zig");
const tevd = @import("tridiag_evd.zig");
const deflate = @import("laed8_br.zig");
const secular = @import("laed9_br.zig");

const Allocator = std.mem.Allocator;
const Mat = ops.Mat;

/// What `stebr` returns: `NoConvergence` when a leaf's QR iteration or a
/// secular root fails (the reference's `INFO > 0`).
pub const Error = Allocator.Error || error{NoConvergence};

/// `SMLSIZ`, `ILAENV(9, ...)` (`ilaenv.f:744`): the leaves' largest order,
/// and the order up to which `xSTERF` takes the whole matrix.
pub const smlsiz = 25;

/// `DSTEBR`: the eigenvalues of the symmetric tridiagonal with diagonal `d`
/// and off-diagonal `e` (`n - 1`) into `d`, increasing; `e` is destroyed.
/// The input is finite (the drivers check). Blocks split where `|e(i)| <=
/// eps sqrt(|d(i)|) sqrt(|d(i + 1)|)`; those past `smlsiz` scaled to unit
/// norm and divided, the rest by `xSTERF`.
pub fn stebr(al: Allocator, comptime F: type, d: []F, e: []F) Error!void {
    const n = d.len;
    if (n <= 1) return;
    if (n <= smlsiz) return tevd.qrAlgorithm(F, d, e[0 .. n - 1], null);
    if (lanst(F, d, e[0 .. n - 1]) <= 0) return;
    const eps = std.math.floatEps(F) / 2; // `DLAMCH('Epsilon')`
    const ws = try Work(F).init(al, n);
    defer ws.deinit(al);
    var start: usize = 0;
    while (start < n) {
        var finish = start;
        while (finish < n - 1) {
            const tiny = eps * @sqrt(@abs(d[finish])) * @sqrt(@abs(d[finish + 1]));
            if (!(@abs(e[finish]) > tiny)) break;
            finish += 1;
        }
        const m = finish - start + 1;
        const db = d[start .. finish + 1];
        const eb = e[start..finish];
        if (m > smlsiz) {
            const orgnrm = lanst(F, db, eb);
            auxiliary.lascl(F, false, orgnrm, 1, column(F, db));
            auxiliary.lascl(F, false, orgnrm, 1, column(F, eb));
            try laed0Br(F, db, eb, ws);
            auxiliary.lascl(F, false, 1, orgnrm, column(F, db));
        } else if (m > 1) {
            try tevd.qrAlgorithm(F, db, eb, null);
        }
        start = finish + 1;
    }
    std.mem.sort(F, d, {}, std.sort.asc(F));
}

/// `xLANST('M')`: the largest `|d(i)|`, `|e(i)|`, NaN if any is NaN.
fn lanst(comptime F: type, d: []const F, e: []const F) F {
    const a = auxiliary.normMax(F, column(F, @constCast(d)));
    if (std.math.isNan(a) or e.len == 0) return a;
    const b = auxiliary.normMax(F, column(F, @constCast(e)));
    return if (std.math.isNan(b)) b else @max(a, b);
}

fn column(comptime F: type, x: []F) Mat(F) {
    return .{ .ptr = x.ptr, .rows = x.len, .cols = 1, .ld = @max(x.len, 1) };
}

/// The merge tree's state and every merge's scratch, sized for the whole.
fn Work(comptime F: type) type {
    return struct {
        idx: []usize,
        vals: []F,
        /// the blocks' ends (`dlaed0_br.f`'s `IWORK(1:N)`)
        ends: []usize,
        /// each block's sorting permutation (`INDXQ`)
        indxq: []usize,
        /// the first and last eigenvector rows per eigenvalue (`BLO`, `BHI`)
        blo: []F,
        bhi: []F,
        /// the rank-one update (`dlaed7_br.f`'s `Z`)
        z: []F,
        /// the parent's boundary rows in the children's bases (`QBR`)
        q: []F,
        /// a leaf's eigenvectors
        leaf: []F,
        deflate: deflate.Scratch(F),
        secular: secular.Scratch(F),

        const Self = @This();

        fn init(al: Allocator, n: usize) Allocator.Error!Self {
            const idx = try al.alloc(usize, 5 * n);
            errdefer al.free(idx);
            const vals = try al.alloc(F, 16 * n + smlsiz * smlsiz);
            var w: Self = undefined;
            w.idx = idx;
            w.vals = vals;
            var ii: usize = 0;
            var vi: usize = 0;
            inline for (.{ "ends", "indxq" }) |f| {
                @field(w, f) = idx[ii..][0..n];
                ii += n;
            }
            inline for (.{ .{ "blo", 1 }, .{ "bhi", 1 }, .{ "z", 1 }, .{ "q", 2 } }) |f| {
                @field(w, f[0]) = vals[vi..][0 .. f[1] * n];
                vi += f[1] * n;
            }
            inline for (.{ "indxp", "indx", "perm" }) |f| {
                @field(w.deflate, f) = idx[ii..][0..n];
                ii += n;
            }
            inline for (.{ .{ "dlambda", 1 }, .{ "w", 1 }, .{ "q2", 2 } }) |f| {
                @field(w.deflate, f[0]) = vals[vi..][0 .. f[1] * n];
                vi += f[1] * n;
            }
            inline for (.{ .{ "delta", 1 }, .{ "wprod", 1 }, .{ "wsgn", 1 }, .{ "tau", 1 }, .{ "near", 2 } }) |f| {
                @field(w.secular, f[0]) = vals[vi..][0 .. f[1] * n];
                vi += f[1] * n;
            }
            w.leaf = vals[vi..][0 .. smlsiz * smlsiz];
            return w;
        }

        fn deinit(w: Self, al: Allocator) void {
            al.free(w.idx);
            al.free(w.vals);
        }
    };
}

/// `DLAED0_BR` (with `DLAED0_MERGE_PAIR` and `DSTEQR_BOUNDARY_BR`): `d`, `e`
/// of one unsplit block past `smlsiz`, scaled to unit norm, to its
/// eigenvalues, increasing. Halved into leaves of at most `smlsiz` by
/// rank-one cuts, the leaves solved with vectors, keeping their first and
/// last rows, then merged pairwise up the tree.
fn laed0Br(comptime F: type, d: []F, e: []F, ws: Work(F)) error{NoConvergence}!void {
    const n = d.len;
    const ends = ws.ends;
    ends[0] = n;
    var subpbs: usize = 1;
    while (ends[subpbs - 1] > smlsiz) : (subpbs *= 2) {
        var j = subpbs;
        while (j > 0) : (j -= 1) {
            ends[2 * j - 1] = (ends[j - 1] + 1) / 2;
            ends[2 * j - 2] = ends[j - 1] / 2;
        }
    }
    for (1..subpbs) |j| ends[j] += ends[j - 1];
    for (ends[0 .. subpbs - 1]) |c| {
        d[c - 1] -= @abs(e[c - 1]);
        d[c] -= @abs(e[c - 1]);
    }
    // the leaves, and each one's rows of its eigenvectors
    for (0..subpbs) |i| {
        const lo = if (i == 0) 0 else ends[i - 1];
        const m = ends[i] - lo;
        const q: Mat(F) = .{ .ptr = ws.leaf.ptr, .rows = m, .cols = m, .ld = m };
        try tevd.qrAlgorithm(F, d[lo..ends[i]], e[lo .. ends[i] - 1], q);
        for (0..m) |j| {
            ws.blo[lo + j] = q.get(0, j);
            ws.bhi[lo + j] = q.get(m - 1, j);
            ws.indxq[lo + j] = j;
        }
    }
    // pairs merged, the rows updated in place in the merged order
    while (subpbs > 1) : (subpbs /= 2) {
        const wantq = subpbs > 2;
        for (0..subpbs / 2) |mi| {
            const lo = if (mi == 0) 0 else ends[2 * mi - 1];
            const hi = ends[2 * mi + 1];
            const cut = ends[2 * mi] - lo;
            try laed7Br(F, wantq, d[lo..hi], ws.indxq[lo..hi], &e[lo + cut - 1], ws.blo[lo..hi], ws.bhi[lo..hi], cut, ws);
        }
        for (0..subpbs / 2) |mi| ends[mi] = ends[2 * mi + 1];
    }
    // from the last merge's order
    const tmp = ws.z[0..n];
    for (tmp, ws.indxq[0..n]) |*x, j| x.* = d[j];
    @memcpy(d, tmp);
}

/// `DLAED7_BR`: merges the eigenvalues of two halves (`d[0..cut]`,
/// `d[cut..]`, each in its `indxq` order) across the rank-one cut `rho`,
/// into `d` with `indxq` the merged increasing order. `wantq` (below the
/// root) carries the parent's boundary rows `blo`, `bhi`, in `d`'s order.
fn laed7Br(comptime F: type, wantq: bool, d: []F, indxq: []usize, rho: *F, blo: []F, bhi: []F, cut: usize, ws: Work(F)) error{NoConvergence}!void {
    const n = d.len;
    const z = ws.z[0..n];
    // z: the left's last rows, the right's first; q the parent's first and
    // last rows of blockdiag(Q_left, Q_right)
    const q: ?[]F = if (wantq) ws.q[0 .. 2 * n] else null;
    for (0..n) |j| {
        const left = j < cut;
        z[j] = if (left) bhi[j] else blo[j];
        if (q) |qq| {
            qq[2 * j] = if (left) blo[j] else 0;
            qq[2 * j + 1] = if (left) 0 else bhi[j];
        }
    }
    const k = deflate.laed8Br(F, d, q, indxq, rho, cut, z, ws.deflate);
    if (k != 0) {
        try secular.laed9Br(F, d[0..k], ws.deflate.q2, q, rho.*, ws.deflate.dlambda, ws.deflate.w, ws.secular);
        auxiliary.lamrg(F, k, n - k, d, true, false, indxq);
    } else for (indxq, 0..) |*x, i| {
        x.* = i;
    }
    if (q) |qq| for (0..n) |j| {
        blo[j] = qq[2 * j];
        bhi[j] = qq[2 * j + 1];
    };
}
