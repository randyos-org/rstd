//! LU with partial pivoting and its solves: OpenBLAS v0.3.34-317's `getrf`
//! (`interface/lapack/getrf.c`, `zgetrf.c`, `lapack/getrf/getrf_single.c`,
//! `getrf_parallel.c`, `lapack/getf2/getf2_k.c`, `zgetf2_k.c`), `getrs`
//! (`lapack/getrs/getrs_single.c`, `zgetrs_single.c`) and `gesv`
//! (`interface/lapack/gesv.c`), Reference-LAPACK's `xGETRF2` semantics
//! (`lapack-netlib/SRC/zgetrf2.f`). Departures and threading: design.md's
//! LAPACK section.

const std = @import("std");
const mem = @import("../../../mem.zig");
const types = @import("../blas/types.zig");
const util = @import("../blas/util.zig");
const portable = @import("../blas/portable.zig");
const threading = @import("../blas/threading.zig");
const driver = @import("../blas/blas3/driver.zig");
const triangular = @import("../blas/blas2/triangular.zig");
const laswp = @import("../blas/ext/laswp/kernel.zig");
const Vtable = @import("../blas/Vtable.zig");
const ops = @import("ops.zig");
const interface = @import("interface.zig");
const lapack_tri = @import("triangular.zig");
const accel = @import("../blas/apple_accel.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Trans = types.Trans;
const isComplex = types.isComplex;
const scalar = ops.scalar;
const isZero = ops.isZero;
const add = ops.add;
const sub = ops.sub;
const mul = ops.mul;
const Work = driver.Work;

/// What `getrf` returns: a slice too short (`a`, or `ipiv` under
/// `min(m, n)`), a leading dimension too small, and the allocator's and
/// cancellation.
pub const GetrfError = interface.ShapeError || Allocator.Error || Io.Cancelable;

/// `?GETRF`: `A = P L U` for the `m x n` `A`, in place: `L` unit lower
/// trapezoidal below the diagonal, `U` upper trapezoidal on and above it.
/// Row `k` was interchanged with row `ipiv[k]`, both 0-based (issue #398),
/// for the first `min(m, n)` entries of `ipiv`.
///
/// Returns the first `i` (0-based) whose `U(i, i)` is exactly zero, or null.
/// A singular `U` is a result, not an error (issue #313): the factorization
/// completes, as the reference's does, and the caller decides. A NaN pivot
/// is not zero: its column below fills with NaN.
pub fn getrf(allocator: Allocator, io: Io, comptime T: type, m: usize, n: usize, a: []T, lda: usize, ipiv: []usize) GetrfError!?usize {
    comptime ops.requireFloat(T, "getrf");
    try util.holds(a.len, m, n, lda);
    if (ipiv.len < @min(m, n)) return error.SliceTooShort;
    if (m == 0 or n == 0) return null;
    if (getrfAccel(T, m, n, a, lda, ipiv)) |k| return k;
    const want = interface.threads(.getrf, T, @as(u64, m) * n) orelse 1;
    return factorWith(allocator, io, T, m, n, a.ptr, lda, ipiv.ptr, want);
}

/// `accel.getrf` (macOS on arm64; issue #840), `INFO` as `getrf`'s result,
/// or null where rstd runs it.
fn getrfAccel(comptime T: type, m: usize, n: usize, a: []T, lda: usize, ipiv: []usize) ??usize {
    const info = accel.getrf(T, @intCast(m), @intCast(n), a.ptr, @intCast(lda), ipiv[0..@min(m, n)]) orelse return null;
    // 1-based at its ABI
    return if (info == 0) @as(?usize, null) else @as(usize, @intCast(info - 1));
}

/// The factorization on up to `want` workers: the first zero pivot, or null.
fn factorWith(allocator: Allocator, io: Io, comptime T: type, m: usize, n: usize, a: [*]T, lda: usize, ipiv: [*]usize, want: u32) driver.Error!?usize {
    const vt = interface.table(T);
    const blk = driver.Blocking.of(T, vt);
    if (want > 1) return factorAhead(allocator, io, T, vt, blk, m, n, a, lda, ipiv, want);
    return factor(allocator, io, T, vt, blk.un, blk.q, m, n, a, lda, ipiv);
}

/// `getrf_single.c`'s block width for `min(m, n) = mn`: half of it rounded up
/// to `unroll_n`, at most `Q`; at most `2 * unroll_n` is `getf2`'s
/// (upstream's POWER8 build stops at `unroll_n`).
fn blockWidth(mn: usize, un: usize, q: usize) ?usize {
    const blocking = @min((mn / 2 + un - 1) / un * un, q);
    return if (blocking <= 2 * un) null else blocking;
}

/// `getrf_single.c`: the first zero pivot's index, or null. Pivots are 0-based
/// from `a`'s first row. Blocks of `blockWidth`, each factored recursively.
pub fn factor(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, un: usize, q: usize, m: usize, n: usize, a: [*]T, lda: usize, ipiv: [*]usize) driver.Error!?usize {
    const mn = @min(m, n);
    const blocking = blockWidth(mn, un, q) orelse return getf2(T, vt, m, n, a, lda, ipiv);
    const one = scalar(T, 1);
    const minus = scalar(T, -1);
    var info: ?usize = null;
    var j: usize = 0;
    while (j < mn) : (j += blocking) {
        const jb = @min(mn - j, blocking);
        const d = a + j * (lda + 1);
        const k = try factor(allocator, io, T, vt, un, q, m - j, jb, d, lda, ipiv + j);
        info = interface.firstOf(info, k, j);
        for (ipiv[j..][0..jb]) |*p| p.* += j;
        const rest = n - j - jb;
        if (rest == 0) continue;
        // A12 := L11^-1 P A12; A22 -= A21 A12, never on the small-matrix
        // path, which `factorAhead`'s parts do not take
        const c = a + (j + jb) * lda;
        laswp.swapRows(T, rest, c, lda, j, jb, ipiv + j, false);
        try driver.trsm(allocator, io, T, vt, true, false, false, false, true, jb, rest, one, d, lda, c + j, lda, 1);
        if (m > j + jb) try driver.gemmBlocked(allocator, io, T, vt, m - j - jb, rest, jb, minus, d + jb, lda, c + j, lda, c + j + jb, lda);
    }
    swapBack(T, mn, blocking, a, lda, ipiv);
    return info;
}

/// The later interchanges, into each block's columns.
fn swapBack(comptime T: type, mn: usize, blocking: usize, a: [*]T, lda: usize, ipiv: [*]const usize) void {
    var j: usize = 0;
    while (j < mn) : (j += blocking) {
        const jb = @min(mn - j, blocking);
        laswp.swapRows(T, jb, a + j * lda, lda, j + jb, mn - (j + jb), ipiv + j + jb, false);
    }
}

/// `getrf_parallel.c` over `factor`'s blocks: the first block recursively,
/// then for each block its interchanges, `trsm` and `gemm` on the columns
/// right of it, which the caller does for the next block's columns and then
/// factors that block (`GETRF_SINGLE`), while the other workers do the rest
/// (`inner_basic_thread`). Every column gets `factor`'s kernel calls, cut at
/// whole `unroll_n` panels, so the bits are `factor`'s.
pub fn factorAhead(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, blk: driver.Blocking, m: usize, n: usize, a: [*]T, lda: usize, ipiv: [*]usize, want: u32) driver.Error!?usize {
    const mn = @min(m, n);
    const blocking = blockWidth(mn, blk.un, blk.q) orelse return getf2(T, vt, m, n, a, lda, ipiv);
    var jb = @min(mn, blocking);
    var info = try factorAhead(allocator, io, T, vt, blk, m, jb, a, lda, ipiv, want);
    var j: usize = 0;
    while (j < mn) {
        for (ipiv[j..][0..jb]) |*p| p.* += j;
        const next = @min(mn - j - jb, blocking);
        if (n > j + jb) {
            var step: Step(T) = .{ .vt = vt, .blk = blk, .m = m, .a = a, .lda = lda, .ipiv = ipiv, .j = j, .jb = jb, .next = next, .allocator = allocator };
            const k = try step.run(io, n - j - jb, want);
            info = interface.firstOf(info, k, j + jb);
        }
        j += jb;
        jb = next;
    }
    swapBack(T, mn, blocking, a, lda, ipiv);
    return info;
}

/// One step of `factorAhead`: block `j` (`jb` wide) applied to the columns
/// right of it, and the next block (`next` wide) factored. Part 0 is the
/// caller's: the next block's columns, rounded up to whole `unroll_n`
/// panels; the rest is cut into chunks at panel edges, which the other
/// workers claim, and the caller too once its block is done. Each part has
/// packing space of its own, allocated before any task starts.
fn Step(comptime T: type) type {
    return struct {
        const Self = @This();

        vt: *const Vtable,
        blk: driver.Blocking,
        m: usize,
        a: [*]T,
        lda: usize,
        ipiv: [*]usize,
        j: usize,
        jb: usize,
        next: usize,
        allocator: Allocator,
        /// Part 0's width, the chunks of the rest, and where each part's
        /// packing space starts in `buf` (part `i + 1` is chunk `i`).
        w0: usize = 0,
        chunks: threading.Chunks = undefined,
        at: [threading.chunks_max + 2]usize = undefined,
        buf: [*]T = undefined,
        claimed: std.atomic.Value(usize) = .init(0),
        /// The next block's first zero pivot (relative), and the caller's
        /// allocation failure there.
        info: ?usize = null,
        failed: bool = false,

        pub fn run(s: *Self, io: Io, rest: usize, want: u32) driver.Error!?usize {
            const un = s.blk.un;
            s.w0 = @min(rest, (s.next + un - 1) / un * un);
            s.chunks = .init(rest - s.w0, want - 1, .even, un);
            const parts = 1 + if (rest > s.w0) s.chunks.count else 0;
            var total: usize = 0;
            for (0..parts) |i| {
                s.at[i] = total;
                total += Work(T).up(s.space(s.width(i)));
            }
            const buf = try s.allocator.alignedAlloc(T, .fromByteUnits(Work(T).alignment), @max(total, 1));
            defer s.allocator.free(buf);
            s.buf = buf.ptr;
            if (parts == 1) {
                try s.work(io, 0, null);
            } else {
                try threading.team(io, s, want);
            }
            if (s.failed) return error.OutOfMemory;
            return s.info;
        }

        /// Part `i`'s first column and width, from the first column right of
        /// block `j`.
        pub fn first(s: *const Self, i: usize) usize {
            return if (i == 0) 0 else s.w0 + s.chunks.get(i - 1).start;
        }

        pub fn width(s: *const Self, i: usize) usize {
            if (i == 0) return s.w0;
            const c = s.chunks.get(i - 1);
            return c.end - c.start;
        }

        /// Elements of packing space for `cols` columns: the `trsm`'s, then
        /// the `gemm`'s, one after the other.
        pub fn space(s: *const Self, cols: usize) usize {
            if (cols == 0) return 0;
            const t = Work(T).total(driver.triLens(s.blk, s.jb, cols, s.jb));
            const g = if (s.m > s.j + s.jb) Work(T).total(Work(T).lens(s.blk, s.m - s.j - s.jb, cols, s.jb, false)) else 0;
            return @max(t, g);
        }

        /// `team`'s plan: every worker that started takes part.
        pub fn plan(_: *Self, started: u32) u32 {
            return started;
        }

        /// Worker 0, the caller, does part 0 and the next block, then claims
        /// chunks as the others do.
        pub fn work(s: *Self, io: Io, i: u32, lent: ?*threading.Lent) Io.Cancelable!void {
            _ = lent;
            if (i == 0) {
                s.update(0);
                if (s.next > 0) {
                    const off = s.j + s.jb;
                    const d = s.a + off * (s.lda + 1);
                    s.info = factor(s.allocator, io, T, s.vt, s.blk.un, s.blk.q, s.m - off, s.next, d, s.lda, s.ipiv + off) catch |e| switch (e) {
                        error.Canceled => return error.Canceled,
                        error.OutOfMemory => {
                            s.failed = true;
                            return;
                        },
                    };
                }
            }
            const count = if (s.chunks.bounds[s.chunks.count] > 0) s.chunks.count else 0;
            while (true) {
                try threading.checkCancel(io);
                const c = s.claimed.fetchAdd(1, .monotonic);
                if (c >= count) return;
                s.update(c + 1);
            }
        }

        /// `inner_basic_thread` on part `i`: block `j`'s interchanges, then
        /// `A12 := L11^-1 A12` and `A22 -= A21 A12` on its columns, with
        /// `factor`'s driver loops on the part's own packing space.
        pub fn update(s: *const Self, i: usize) void {
            const cols = s.width(i);
            if (cols == 0) return;
            const j = s.j;
            const jb = s.jb;
            const lda = s.lda;
            const d = s.a + j * (lda + 1);
            const b = s.a + (j + jb + s.first(i)) * lda;
            const base = s.buf + s.at[i];
            laswp.swapRows(T, cols, b, lda, j, jb, s.ipiv + j, false);
            const o: driver.Tri = .{ .left = true, .upper = false, .trans = false, .conj = false, .unit = true };
            driver.triDrive(T, s.vt, true, o, s.blk, jb, cols, d, lda, null, b + j, lda, Work(T).at(base, driver.triLens(s.blk, jb, cols, jb)));
            if (s.m <= j + jb) return;
            const rows = s.m - j - jb;
            const work_g = Work(T).at(base, Work(T).lens(s.blk, rows, cols, jb, false));
            driver.drive(T, s.vt, .{ .rows = .n, .cols = .n }, s.blk, rows, cols, jb, scalar(T, -1), .{ .p = d + jb, .ld = lda }, .{ .p = b + j, .ld = lda }, null, b + j + jb, lda, work_g);
        }
    };
}

/// What `getrs` returns: a `trans` it does not take, a pivot that is not a
/// row of `A`, a slice too short or a leading dimension too small, and the
/// allocator's and cancellation.
pub const GetrsError = interface.TransError || interface.PivotError || interface.ShapeError || Allocator.Error || Io.Cancelable;

/// `?GETRS`: solves `op(A) X = B` for the `n x nrhs` `X`, in place of `B`,
/// given `getrf`'s `P L U` of `A` and its `ipiv`. `trans` is `n`, `t`, `r`
/// (`conj(A)`) or `c` (`Aᴴ`); real `r` and `c` are `n` and `t`, as in the
/// level-3 routines, where the reference takes no `r`. No singularity check,
/// as in the reference. Every pivot must lie in `0..n - 1`, else
/// `error.PivotOutOfRange` before `B` is touched; the reference checks none.
pub fn getrs(allocator: Allocator, io: Io, comptime T: type, trans: Trans, n: usize, nrhs: usize, a: []const T, lda: usize, ipiv: []const usize, b: []T, ldb: usize) GetrsError!void {
    comptime ops.requireFloat(T, "getrs");
    try interface.takes(trans);
    try util.holds(a.len, n, n, lda);
    try util.holds(b.len, n, nrhs, ldb);
    try interface.pivotsIn(ipiv, n);
    if (n == 0 or nrhs == 0) return;
    if (try accel.getrs(allocator, T, trans, @intCast(n), @intCast(nrhs), a.ptr, @intCast(lda), ipiv[0..n], b.ptr, @intCast(ldb))) return;
    const want = interface.threads(.getrs, T, n) orelse 1;
    try solve(allocator, io, T, trans, n, nrhs, a.ptr, lda, ipiv.ptr, b.ptr, ldb, want);
}

/// `getrs_single.c`, `zgetrs_single.c`, the `trsm`s on up to `want` workers
/// (`getrs_parallel.c`).
pub fn solve(allocator: Allocator, io: Io, comptime T: type, trans: Trans, m: usize, r: usize, a: [*]const T, la: usize, ipiv: [*]const usize, b: [*]T, lb: usize, want: u32) driver.Error!void {
    const vt = interface.table(T);
    const t = trans.transposed();
    // `L Y = P B`, then `U X = Y`; transposed `op(U) Y = B`, then
    // `op(L) Z = Y` and `X = Pᵀ Z`: the first solve is upper when transposed,
    // and `L` is the unit one
    if (!t) laswp.swapRows(T, r, b, lb, 0, m, ipiv, false);
    for ([_]bool{ t, !t }) |upper| {
        if (!isComplex(T) and r == 1) {
            try triangular.trsvTable(allocator, io, T, vt, if (upper) .u else .l, trans, if (upper) .n else .u, m, a, la, b);
        } else {
            try driver.trsm(allocator, io, T, vt, true, upper, t, trans.conjA(), !upper, m, r, scalar(T, 1), a, la, b, lb, want);
        }
    }
    if (t) laswp.swapRows(T, r, b, lb, 0, m, ipiv, true);
}

/// What `gesv` returns: `Singular` when `U` has a zero on its diagonal, a
/// slice too short or a leading dimension too small, and the allocator's
/// and cancellation.
pub const GesvError = interface.SingularError || interface.ShapeError || Allocator.Error || Io.Cancelable;

/// `?GESV`: solves `A X = B` for the `n x nrhs` `X`, in place of `B`, by
/// `getrf` (into `A` and `ipiv`) and `getrs`. When some `U(i, i)` is
/// exactly zero, returns `error.Singular`, the factorization complete and
/// `B` untouched; `getrf` on its own reports which `i`.
pub fn gesv(allocator: Allocator, io: Io, comptime T: type, n: usize, nrhs: usize, a: []T, lda: usize, ipiv: []usize, b: []T, ldb: usize) GesvError!void {
    comptime ops.requireFloat(T, "gesv");
    try util.holds(a.len, n, n, lda);
    try util.holds(b.len, n, nrhs, ldb);
    if (ipiv.len < n) return error.SliceTooShort;
    if (n == 0) return;
    const want = interface.threads(.gesv, T, @as(u64, n) * n) orelse 1;
    if (getrfAccel(T, n, n, a, lda, ipiv)) |k| {
        if (k != null) return error.Singular;
        if (nrhs == 0) return;
        if (try accel.getrs(allocator, T, .n, @intCast(n), @intCast(nrhs), a.ptr, @intCast(lda), ipiv[0..n], b.ptr, @intCast(ldb))) return;
        return solve(allocator, io, T, .n, n, nrhs, a.ptr, lda, ipiv.ptr, b.ptr, ldb, want);
    }
    if (try factorWith(allocator, io, T, n, n, a.ptr, lda, ipiv.ptr, want) != null) return error.Singular;
    if (nrhs == 0) return;
    try solve(allocator, io, T, .n, n, nrhs, a.ptr, lda, ipiv.ptr, b.ptr, ldb, want);
}

/// What `getri` returns: `Singular` when `U` has a zero on its diagonal, a
/// pivot that is not a row of `A`, a slice too short or a leading dimension
/// too small, and the allocator's and cancellation.
pub const GetriError = interface.SingularError || interface.PivotError || interface.ShapeError || Allocator.Error || Io.Cancelable;

/// `?GETRI`: `A⁻¹` in place of `getrf`'s `P L U` of `A` and its `ipiv`, by
/// faer's `lu/partial_pivoting/inverse.rs` (faer v0.24.4,
/// <https://codeberg.org/sarah-quinones/faer>, MIT): `U⁻¹` and `L⁻¹` in
/// place (`trtri`), `U⁻¹ L⁻¹` by `trmm` over a copy of `L⁻¹`, then the
/// columns swapped back, last pivot first (`dgetri.f`'s order). A zero on
/// `U`'s diagonal returns `error.Singular` before `A` is touched. Every
/// pivot must lie in `0..n - 1`, as `getrs` checks.
pub fn getri(allocator: Allocator, io: Io, comptime T: type, n: usize, a: []T, lda: usize, ipiv: []const usize) GetriError!void {
    comptime ops.requireFloat(T, "getri");
    try util.holds(a.len, n, n, lda);
    try interface.pivotsIn(ipiv, n);
    if (n == 0) return;
    const m = n;
    const la = lda;
    // `U` first: its singularity check leaves `A` as `getrf` left it
    for ([_]types.Uplo{ .u, .l }) |uplo| {
        lapack_tri.trtri(allocator, io, T, uplo, if (uplo == .u) .n else .u, n, a, lda) catch |e| switch (e) {
            error.Singular, error.OutOfMemory, error.Canceled => |x| return x,
            // the arguments are checked above
            error.UnsupportedUplo, error.SliceTooShort, error.LeadingDimensionTooSmall => unreachable,
        };
    }
    const w = try allocator.alloc(T, m * m);
    defer allocator.free(w);
    for (0..m) |j| for (0..m) |i| {
        w[i + j * m] = if (i > j) a[i + j * la] else scalar(T, if (i == j) 1 else 0);
    };
    const vt = interface.table(T);
    const want = interface.threads(.trtri, T, m) orelse 1;
    try driver.trmm(allocator, io, T, vt, true, true, false, false, false, m, m, scalar(T, 1), a.ptr, la, w.ptr, m, want);
    for (0..m) |j| @memcpy(a[j * la ..][0..m], w[j * m ..][0..m]);
    var j = m - 1;
    while (j > 0) {
        j -= 1;
        const p = ipiv[j];
        if (p != j) for (0..m) |i| mem.swap(T, &a[i + j * la], &a[i + p * la]);
    }
}

/// `getf2_k.c`, `zgetf2_k.c`: left-looking, a column at a time. Column `j`
/// takes the earlier interchanges, then `U(:j, j)` by forward substitution,
/// then `A(j:, j) -= L(j:, :j) U(:j, j)`, its pivot, and `L(j + 1:, j)`.
pub fn getf2(comptime T: type, vt: *const Vtable, m: usize, n: usize, a: [*]T, lda: usize, ipiv: [*]usize) ?usize {
    const iamax = portable.kernel(T, vt, .iamax_k);
    var info: ?usize = null;
    for (0..n) |j| {
        const b = a + j * lda;
        const len = @min(j, m);
        for (0..len) |i| {
            const ip = ipiv[i];
            if (ip != i) mem.swap(T, &b[i], &b[ip]);
        }
        for (1..@max(1, len)) |i| {
            var s = scalar(T, 0);
            for (0..i) |k| s = add(T, s, mul(T, a[i + k * lda], b[k]));
            b[i] = sub(T, b[i], s);
        }
        if (j >= m) continue;
        for (0..j) |k| {
            const x = b[k];
            const l = a + k * lda;
            for (j..m) |i| b[i] = sub(T, b[i], mul(T, l[i], x));
        }
        const p = j + @as(usize, @intCast(iamax(@intCast(m - j), @ptrCast(b + j), 1))) - 1;
        ipiv[j] = p;
        const piv = b[p];
        if (isZero(T, piv)) {
            if (info == null) info = j;
            continue;
        }
        if (p != j) for (0..j + 1) |k| mem.swap(T, &a[j + k * lda], &a[p + k * lda]);
        ops.scaleBelow(T, b[j + 1 .. m], piv);
    }
    return info;
}
