//! Threaded level 3 in BLIS's model (Van Zee and van de Geijn; BLIS 2.1's
//! `docs/Multithreading.md`), written from its design. Workers split `C`
//! `ic x jc` ways. The `ic` workers of a `jc` group pack each depth slice of
//! their columns of `op(B)` together into one buffer, between a barrier that
//! frees it and one that fills it (`bli_packm_int.c`); each packs its own
//! blocks of `op(A)`. Loops, packing, kernels and depth slices are the
//! serial driver's and row ranges start on row panels, so every element of
//! `C` gets the serial kernel calls in the serial order: the same bits for
//! any worker count, `Io` or scheduling.

const std = @import("std");
const sysinfo = @import("builtin");
const blas = @import("../__root__.zig");
const Vtable = @import("../Vtable.zig");
const driver = @import("driver.zig");
const threading = @import("../threading.zig");
const types = @import("../types.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Blocking = driver.Blocking;
const Chunks = threading.Chunks;
const Desc = driver.Desc;
const Matrix = driver.Matrix;

const Error = driver.Error;

/// The workers `family`'s rule (`threading.md`) gives `work` on the
/// selected core, for the drivers' `want`.
pub fn workers(comptime family: threading.Family, comptime T: type, work: u64) u32 {
    // selects the core first, which fills `blas.config`
    const r = threading.rule(family, T, blas.coreName());
    return @min(threading.workers(&blas.config, r.sizeMax(work)), threading.chunks_max);
}

/// `driver.run`'s loops on up to `want` workers.
pub fn run(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, d: Desc, m: usize, n: usize, k: usize, alpha: T, x: Matrix(T), y: Matrix(T), twice: ?T, beta: ?T, c: [*]T, ldc: usize, want: u32) Error!void {
    return Team(T, true).run(allocator, io, vt, d, m, n, k, alpha, x, y, twice, beta, c, ldc, @min(want, threading.chunks_max));
}

fn ceilDiv(a: usize, b: usize) usize {
    return (a + b - 1) / b;
}

/// A chunk of `Chunks`, or an empty range at `len` past its count.
fn part(c: *const Chunks, i: usize, len: usize) threading.Chunk {
    return if (i < c.count) c.get(i) else .{ .index = i, .start = len, .end = len };
}

/// `ic x jc` ways for `w` workers: the most workers that leave each part a
/// row panel and a column panel, then the least packing (`n * ic + m * jc`,
/// as BLIS's `bli_thread_partition_2x2` balances `m / ic` against `n / jc`),
/// more `ic` on a tie. Returns `ic`, `jc`.
fn ways(m: usize, n: usize, w: usize, um: usize, un: usize) [2]usize {
    const max_ic = ceilDiv(m, um);
    const max_jc = ceilDiv(n, un);
    var best: [2]usize = .{ 1, 1 };
    var best_cost: usize = std.math.maxInt(usize);
    var ic: usize = 1;
    while (ic <= @min(w, max_ic)) : (ic += 1) {
        const jc = @min(w / ic, max_jc);
        const cost = n * ic + m * jc;
        const more = ic * jc > best[0] * best[1];
        if (more or (ic * jc == best[0] * best[1] and cost <= best_cost)) {
            best = .{ ic, jc };
            best_cost = cost;
        }
    }
    return best;
}

/// One call's shared state; `fused` multiplies a worker's own chunks of `B`
/// with its first row block as it packs them.
fn Team(comptime T: type, comptime fused: bool) type {
    return struct {
        const Self = @This();
        const W = driver.Work(T);

        vt: *const Vtable,
        d: Desc,
        blk: Blocking,
        m: usize,
        n: usize,
        k: usize,
        alpha: T,
        x: Matrix(T),
        y: Matrix(T),
        /// `syr2k`'s second pass: its `alpha`, the operands swapped.
        twice: ?T,
        beta: ?T,
        c: [*]T,
        ldc: usize,
        /// Per worker a row block of `op(A)` and a tile; per group a slice of `B`.
        sa: [*]T,
        sa_len: usize,
        sb: [*]T,
        sb_len: usize = 0,
        barriers: []threading.Barrier,
        ic: usize = 1,
        jc: usize = 1,
        cols: Chunks = undefined,
        /// Columns of a group's slice of `B`.
        nc: usize = 0,

        /// The groups' columns for `jc` groups; a triangle's balanced by cost
        /// and cut where the serial driver's chunks start (`snap`).
        fn groups(d: Desc, blk: Blocking, n: usize, jc: usize) Chunks {
            const cost: threading.Shape = switch (d.region) {
                .full => .even,
                .upper => .rising,
                .lower => .falling,
            };
            const c: Chunks = .init(n, jc, cost, blk.un);
            if (d.region == .full) return c;
            var out: Chunks = .{ .count = 0, .bounds = undefined };
            out.bounds[0] = 0;
            for (c.bounds[1 .. c.count + 1]) |b| {
                const s = snap(blk, n, b);
                if (s > out.bounds[out.count]) {
                    out.count += 1;
                    out.bounds[out.count] = s;
                }
            }
            return out;
        }

        /// The start of the serial driver's chunk at or before column `x`:
        /// `drive` steps `chunk` from each `R` block's start. A tile across a
        /// triangle's diagonal goes through `tmp`, so a triangle's tiles must
        /// be the serial ones; a full region's columns may split anywhere.
        fn snap(blk: Blocking, n: usize, x: usize) usize {
            const js = x / blk.r * blk.r;
            if (x == @min(js + blk.r, n)) return x;
            return js + (x - js) / (3 * blk.un) * (3 * blk.un);
        }

        /// The end of the slice from `j0` in a group ending at `end`.
        fn sliceEnd(j: *const Self, j0: usize, end: usize) usize {
            if (j.d.region == .full) return @min(j0 + j.nc, end);
            const be = blockEnd(j.blk, j.n, j0);
            const e = snap(j.blk, j.n, @min(j0 + j.nc, end, be));
            return if (e > j0) e else j0 + j.blk.chunk(be - j0);
        }

        /// The end of the serial driver's column block holding `jj`.
        fn blockEnd(blk: Blocking, n: usize, jj: usize) usize {
            return @min(jj / blk.r * blk.r + blk.r, n);
        }

        /// The chunk at `jj` of the slice ending at `j1`: the serial one for a
        /// triangle.
        fn step(j: *const Self, jj: usize, j1: usize) usize {
            return j.blk.chunk(if (j.d.region == .full) j1 - jj else blockEnd(j.blk, j.n, jj) - jj);
        }

        /// Columns of a group's slice of `B` for `jc` groups.
        fn sliceCols(d: Desc, blk: Blocking, n: usize, jc: usize) usize {
            const g = groups(d, blk, n, jc);
            var widest: usize = 0;
            for (0..g.count) |i| widest = @max(widest, g.get(i).end - g.get(i).start);
            return @min(blk.r, widest);
        }

        fn depth(blk: Blocking, k: usize) usize {
            return @min(k, blk.q + blk.um);
        }

        fn run(allocator: Allocator, io: Io, vt: *const Vtable, d: Desc, m: usize, n: usize, k: usize, alpha: T, x: Matrix(T), y: Matrix(T), twice: ?T, beta: ?T, c: [*]T, ldc: usize, want: u32) Error!void {
            const blk: Blocking = .of(T, vt);
            const sa_len = W.up(@min(m, blk.p) * depth(blk, k)) + W.up(if (d.region != .full) blk.um * blk.un else 0);
            // the most that any number of groups `plan` may pick takes
            var sb_len: usize = 0;
            for (1..want + 1) |jc| sb_len = @max(sb_len, jc * W.up(sliceCols(d, blk, n, jc) * depth(blk, k)));
            const buf = try allocator.alignedAlloc(T, .fromByteUnits(W.alignment), want * sa_len + sb_len);
            defer allocator.free(buf);
            const barriers = try allocator.alloc(threading.Barrier, want);
            defer allocator.free(barriers);
            var job: Self = .{
                .vt = vt,
                .d = d,
                .blk = blk,
                .m = m,
                .n = n,
                .k = k,
                .alpha = alpha,
                .x = x,
                .y = y,
                .twice = twice,
                .beta = beta,
                .c = c,
                .ldc = ldc,
                .sa = buf.ptr,
                .sa_len = sa_len,
                .sb = buf.ptr + want * sa_len,
                .barriers = barriers,
            };
            try threading.team(io, &job, want);
        }

        pub fn plan(j: *Self, started: u32) u32 {
            const blk = j.blk;
            const w = ways(j.m, j.n, started, blk.um, blk.un);
            j.ic = w[0];
            j.jc = w[1];
            j.cols = groups(j.d, blk, j.n, j.jc);
            j.nc = sliceCols(j.d, blk, j.n, j.jc);
            j.sb_len = W.up(j.nc * depth(blk, j.k));
            for (j.barriers[0..j.jc]) |*b| b.* = .{ .parties = @intCast(j.ic) };
            return @intCast(j.ic * j.jc);
        }

        /// Worker `r`'s rows for columns `[j0, j1)` of the serial driver's
        /// column block at `js`: the rows `drive` gives that block, less the
        /// row panels the columns' triangle misses, split on those panels.
        fn rowsOf(j: *const Self, r: usize, js: usize, j0: usize, j1: usize) threading.Chunk {
            const um = j.blk.um;
            var lo: usize = 0;
            var hi = j.m;
            switch (j.d.region) {
                .full => {},
                .lower => lo = js + (j0 - js) / um * um,
                .upper => hi = @min(j.m, js + j.blk.r, (j1 + um - 1) / um * um),
            }
            const parts: Chunks = .init(hi - lo, j.ic, .even, um);
            const p = part(&parts, r, hi - lo);
            return .{ .index = r, .start = lo + p.start, .end = lo + p.end };
        }

        /// The kernel on rows `[is, is + mi)` packed at `sa` and chunk `[jj,
        /// jj + w)` of the slice at `j0` packed at `sb`; `s` scales `C` first.
        fn tile(j: *const Self, alpha: T, min_l: usize, s: ?T, is: usize, mi: usize, jj: usize, w: usize, j0: usize, sa: [*]T, sb: [*]T, tmp: [*]T) void {
            driver.block(T, j.vt, j.d, j.blk, is, mi, jj, w, min_l, alpha, s, sa, sb + min_l * (jj - j0), j.c, j.ldc, tmp);
        }

        pub fn work(j: *Self, io: Io, index: u32) Io.Cancelable!void {
            const blk = j.blk;
            const r = index % j.ic;
            const g = index / j.ic;
            const cols = part(&j.cols, g, j.n);
            const sa = j.sa + index * j.sa_len;
            const tmp = sa + W.up(@min(j.m, blk.p) * depth(blk, j.k));
            const sb = j.sb + g * j.sb_len;
            const bar = &j.barriers[g];
            var first = true;
            for (0..@as(usize, if (j.twice == null) 1 else 2)) |pass| {
                const x = if (pass == 0) j.x else j.y;
                const y = if (pass == 0) j.y else j.x;
                const alpha = if (pass == 0) j.alpha else j.twice.?;
                const beta = if (pass == 0) j.beta else null;
                var j0 = cols.start;
                while (j0 < cols.end) {
                    // a triangle's rows follow the serial column block
                    const js = j0 / blk.r * blk.r;
                    const j1 = j.sliceEnd(j0, cols.end);
                    defer j0 = j1;
                    const rows = j.rowsOf(r, js, j0, j1);
                    var chunks: usize = 0;
                    var jj = j0;
                    while (jj < j1) : (chunks += 1) jj += j.step(jj, j1);
                    // this worker's chunks of the slice
                    const c0 = chunks * r / j.ic;
                    const c1 = chunks * (r + 1) / j.ic;
                    var ls: usize = 0;
                    var min_l: usize = 0;
                    while (ls < j.k) : (ls += min_l) {
                        min_l = Blocking.split(j.k - ls, blk.q, blk.um);
                        const s = if (ls == 0) beta else null;
                        try threading.checkCancel(io);
                        // the buffer is free once every peer is done with the last slice
                        if (!first) try bar.wait(io);
                        first = false;
                        var min_i: usize = if (fused and rows.end > rows.start) Blocking.split(rows.end - rows.start, blk.p, blk.um) else 0;
                        if (min_i > 0) driver.packRows(T, j.vt, j.d, x, min_l, min_i, ls, rows.start, sa);
                        jj = j0;
                        var i: usize = 0;
                        while (jj < j1) : (i += 1) {
                            const w = j.step(jj, j1);
                            if (i >= c0 and i < c1) {
                                driver.packCols(T, j.vt, j.d, y, min_l, w, ls, jj, sb + min_l * (jj - j0));
                                if (min_i > 0) j.tile(alpha, min_l, s, rows.start, min_i, jj, w, j0, sa, sb, tmp);
                            }
                            jj += w;
                        }
                        // the slice is whole
                        try bar.wait(io);
                        var is = rows.start;
                        if (min_i > 0) {
                            // the first row block's other chunks
                            jj = j0;
                            i = 0;
                            while (jj < j1) : (i += 1) {
                                const w = j.step(jj, j1);
                                if (i < c0 or i >= c1) j.tile(alpha, min_l, s, rows.start, min_i, jj, w, j0, sa, sb, tmp);
                                jj += w;
                            }
                            is += min_i;
                        }
                        while (is < rows.end) : (is += min_i) {
                            min_i = Blocking.split(rows.end - is, blk.p, blk.um);
                            driver.packRows(T, j.vt, j.d, x, min_l, min_i, ls, is, sa);
                            jj = j0;
                            while (jj < j1) {
                                const w = j.step(jj, j1);
                                j.tile(alpha, min_l, s, is, min_i, jj, w, j0, sa, sb, tmp);
                                jj += w;
                            }
                        }
                    }
                }
            }
        }
    };
}

// ---- tests -------------------------------------------------------------------------

const zigstd = @import("../../../zigstd.zig");
const portable = @import("../portable.zig");
const gemm_kernels = @import("gemm/kernel.zig");
const guard = @import("../_test/blas2/guard.zig");
const interface = @import("interface.zig");

/// Mixed magnitudes, so a sum's bits follow its order.
fn fill(comptime T: type, x: []T, seed: u64) void {
    var prng: std.Random.DefaultPrng = .init(seed);
    const rand = prng.random();
    const F = types.Float(T);
    for (x) |*v| {
        const re = std.math.ldexp(rand.float(F) * 2 - 1, rand.intRangeAtMost(i32, -12, 12));
        v.* = if (comptime types.isComplex(T)) .init(re, std.math.ldexp(rand.float(F) * 2 - 1, rand.intRangeAtMost(i32, -12, 12))) else re;
    }
}

fn scalar(comptime T: type, re: f64, im: f64) T {
    return if (comptime types.isComplex(T)) .init(@floatCast(re), @floatCast(im)) else @floatCast(re);
}

fn sameBits(comptime T: type, a: []const T, b: []const T) bool {
    return std.mem.eql(u8, std.mem.sliceAsBytes(a), std.mem.sliceAsBytes(b));
}

/// The portable kernels, blocking shrunk as `check.zig` shrinks it, so every
/// block and panel edge comes at small sizes; `r` leaves column tails inside.
fn shrunk(comptime T: type) Vtable {
    var vt: Vtable = .{};
    gemm_kernels.installGeneric(&vt);
    const blk = Blocking.of(T, &vt);
    const pre = comptime portable.letter(T).? ++ "gemm_";
    @field(vt, pre ++ "p") = @intCast(2 * blk.um);
    @field(vt, pre ++ "q") = 5;
    @field(vt, pre ++ "r") = @intCast(6 * blk.un + 1);
    return vt;
}

/// The `op` pairs a type has: transposes, and for complex conjugates too.
fn transPairs(comptime T: type) []const [2]types.Trans {
    return comptime blk: {
        const all: []const types.Trans = if (types.isComplex(T)) &.{ .n, .t, .r, .c } else &.{ .n, .t };
        var out: [all.len * all.len][2]types.Trans = undefined;
        for (all, 0..) |x, i| for (all, 0..) |y, j| {
            out[i * all.len + j] = .{ x, y };
        };
        const final = out;
        break :blk &final;
    };
}

/// One shape against `driver.gemm`, bit for bit, at 1 to 9 workers on `io`.
/// Operand padding holds NaN; `C`'s, a sentinel that must survive; scratch
/// comes poisoned.
fn shape(comptime T: type, vt: *const Vtable, io: Io, tr: [2]types.Trans, m: usize, n: usize, k: usize, alpha: T, beta: T) !void {
    const al = std.testing.allocator;
    const ta = tr[0].transposed();
    const tb = tr[1].transposed();
    const conj: driver.Conj = .of(types.isComplex(T) and tr[0].conjA(), types.isComplex(T) and tr[1].conjA());
    const pad = 3;
    const ra = if (ta) k else m;
    const ca = if (ta) m else k;
    const rb = if (tb) n else k;
    const cb = if (tb) k else n;
    const a = try al.alloc(T, (ra + pad) * ca);
    defer al.free(a);
    const b = try al.alloc(T, (rb + pad) * cb);
    defer al.free(b);
    const c0 = try al.alloc(T, (m + pad) * n + pad);
    defer al.free(c0);
    const want = try al.alloc(T, c0.len);
    defer al.free(want);
    const got = try al.alloc(T, c0.len);
    defer al.free(got);
    fill(T, a, m * 31 + k);
    fill(T, b, n * 17 + k);
    fill(T, c0, m + n);
    const nan = scalar(T, std.math.nan(f64), std.math.nan(f64));
    for (0..ca) |j| @memset(a[j * (ra + pad) + ra ..][0..pad], nan);
    for (0..cb) |j| @memset(b[j * (rb + pad) + rb ..][0..pad], nan);
    const sentinel = scalar(T, 1234.5, -6.5);
    for (0..n) |j| @memset(c0[j * (m + pad) + m ..][0..pad], sentinel);
    @memset(c0[n * (m + pad) ..], sentinel);

    var poison: guard.Poison = .{ .child = al };
    const pa = poison.allocator();
    @memcpy(want, c0);
    try driver.gemm(pa, Io.failing, T, vt, ta, tb, conj, m, n, k, alpha, a.ptr, ra + pad, b.ptr, rb + pad, beta, want.ptr, m + pad, 1);
    for (1..10) |w| {
        @memcpy(got, c0);
        try driver.gemm(pa, io, T, vt, ta, tb, conj, m, n, k, alpha, a.ptr, ra + pad, b.ptr, rb + pad, beta, got.ptr, m + pad, @intCast(w));
        if (!sameBits(T, want, got)) {
            std.debug.print("{s} w={d} {t}{t} m={d} n={d} k={d}\n", .{ @typeName(T), w, tr[0], tr[1], m, n, k });
            return error.TestUnexpectedResult;
        }
    }
}

/// Every `op` pair over shapes around every block and panel edge (complex,
/// fewer shapes over its sixteen pairs).
fn shapes(comptime T: type, io: Io) !void {
    const vt = shrunk(T);
    const blk = Blocking.of(T, &vt);
    const um = blk.um;
    const un = blk.un;
    const rr = blk.r;
    const ms = [_]usize{ 1, um + 1, 2 * blk.p + um + 1, 4 * blk.p + 3, 9 * um + 5 };
    const ns = [_]usize{ 1, un + 1, 3 * un + 1, rr + 1, 2 * rr + 3, 5 * rr };
    const ks = [_]usize{ 1, blk.q + 1, 2 * blk.q + 3 };
    const step: usize = if (comptime types.isComplex(T)) 5 else 2;
    const alpha = scalar(T, -0.75, 0.25);
    const beta = scalar(T, 0.5, -0.125);
    for (transPairs(T)) |tr| {
        var x: usize = 0;
        for (ms) |m| for (ns) |n| for (ks) |k| {
            x += 1;
            if (x % step != 0) continue;
            try shape(T, &vt, io, tr, m, n, k, alpha, beta);
        };
        // beta 0 and 1, which the first slice folds differently
        try shape(T, &vt, io, tr, 9 * um + 5, 2 * rr + 3, 2 * blk.q + 3, scalar(T, 1.5, 0), scalar(T, 0, 0));
        try shape(T, &vt, io, tr, 9 * um + 5, 2 * rr + 3, 2 * blk.q + 3, scalar(T, 1, 0), scalar(T, 1, 0));
    }
}

/// `async_limit` and `concurrent_limit` of the threaded tests: none, fewer
/// than asked, oversubscribed.
const limits = [_]Io.Limit{ .nothing, .limited(3), .limited(8) };

const precisions = .{ f32, f64, types.Complex32, types.Complex64 };

test "threaded gemm matches the serial driver bit for bit at 1 to 9 workers, on every Io" {
    inline for (precisions) |T| {
        try shapes(T, Io.failing);
        if (comptime threading.hosted) for (limits) |limit| {
            var threaded: Io.Threaded = .init(std.testing.allocator, .{ .async_limit = limit, .concurrent_limit = limit });
            defer threaded.deinit();
            try shapes(T, threaded.io());
        };
    }
}

/// The routines over the threaded loops besides `gemm`.
const Routine = enum { symm, hemm, syrk, herk, syr2k, her2k };

/// `r` on `vt` at up to `want` workers: `C` is `m x n` (`n x n` for the rank
/// updates), the triangle or symmetric operand `lo`, the depth `k`; `opt`
/// is the side or `trans`, then `upper`.
fn routine(comptime T: type, vt: *const Vtable, al: Allocator, io: Io, r: Routine, opt: [2]bool, m: usize, n: usize, k: usize, alpha: T, beta: T, a: [*]const T, b: [*]const T, ld: usize, c: [*]T, want: u32) !void {
    const re = if (comptime types.isComplex(T)) scalar(T, alpha.re, 0) else alpha;
    switch (r) {
        .symm => try driver.symm(al, io, T, vt, false, opt[0], opt[1], m, n, alpha, a, ld, b, ld, beta, c, ld, want),
        .syrk => try driver.syrk(al, io, T, vt, false, opt[1], opt[0], n, k, alpha, a, ld, null, 0, beta, c, ld, want),
        .syr2k => try driver.syrk(al, io, T, vt, false, opt[1], opt[0], n, k, alpha, a, ld, b, ld, beta, c, ld, want),
        .hemm, .herk, .her2k => if (comptime types.isComplex(T)) switch (r) {
            .hemm => try driver.symm(al, io, T, vt, true, opt[0], opt[1], m, n, alpha, a, ld, b, ld, beta, c, ld, want),
            .herk => try driver.syrk(al, io, T, vt, true, opt[1], opt[0], n, k, re, a, ld, null, 0, beta, c, ld, want),
            else => try driver.syrk(al, io, T, vt, true, opt[1], opt[0], n, k, alpha, a, ld, b, ld, beta, c, ld, want),
        } else unreachable,
    }
}

/// `r`'s shapes around every block, panel and triangle edge, all options,
/// against the serial loops bit for bit at 1 to 9 workers on `io`.
fn routineShapes(comptime T: type, vt: *const Vtable, io: Io, r: Routine) !void {
    const al = std.testing.allocator;
    const blk = Blocking.of(T, vt);
    const herm = r == .herk or r == .her2k;
    const rank = r != .symm and r != .hemm;
    const ns = [_]usize{ 1, blk.un + 1, 2 * blk.p + blk.um + 1, blk.r + 1, 2 * blk.r + 3 };
    const ks = [_]usize{ 1, blk.q + 1, 2 * blk.q + 3 };
    const ld = 2 * blk.r + 2 * blk.q + 8;
    const a = try al.alloc(T, ld * ld);
    defer al.free(a);
    const b = try al.alloc(T, ld * ld);
    defer al.free(b);
    const c0 = try al.alloc(T, ld * ld);
    defer al.free(c0);
    const want = try al.alloc(T, c0.len);
    defer al.free(want);
    const got = try al.alloc(T, c0.len);
    defer al.free(got);
    fill(T, a, 11);
    fill(T, b, 12);
    fill(T, c0, 13);
    var poison: guard.Poison = .{ .child = al };
    const pa = poison.allocator();
    const alpha = scalar(T, -0.75, 0.25);
    const beta = if (herm) scalar(T, 0.5, 0) else scalar(T, 0.5, -0.125);
    var x: usize = 0;
    for ([_][2]bool{ .{ false, false }, .{ true, false }, .{ false, true }, .{ true, true } }) |opt| {
        for (ns) |m| for (ns) |n| for (ks) |k| {
            // the rank updates' `C` is square; `symm` has no depth
            if ((rank and m != n) or (!rank and k != 1)) continue;
            x += 1;
            if (x % 2 != 0) continue;
            @memcpy(want, c0);
            try routine(T, vt, pa, Io.failing, r, opt, m, n, k, alpha, beta, a.ptr, b.ptr, ld, want.ptr, 1);
            for (1..10) |w| {
                @memcpy(got, c0);
                try routine(T, vt, pa, io, r, opt, m, n, k, alpha, beta, a.ptr, b.ptr, ld, got.ptr, @intCast(w));
                if (!sameBits(T, want, got)) {
                    std.debug.print("{s} {t} w={d} opt={any} m={d} n={d} k={d}\n", .{ @typeName(T), r, w, opt, m, n, k });
                    return error.TestUnexpectedResult;
                }
            }
        };
    }
}

/// The global vtable with `check.zig`'s shrunk blocking: the core's own
/// copies and kernels, every edge at small sizes.
fn shrunkCore(comptime T: type) Vtable {
    var vt = blas.vtable().*;
    const blk = Blocking.of(T, &vt);
    const pre = comptime portable.letter(T).? ++ "gemm_";
    @field(vt, pre ++ "p") = @intCast(2 * blk.um);
    @field(vt, pre ++ "q") = 5;
    @field(vt, pre ++ "r") = @intCast(6 * blk.un + 1);
    return vt;
}

test "threaded symm, hemm and rank-k updates match the serial loops bit for bit at 1 to 9 workers" {
    inline for (precisions) |T| {
        const vt = shrunkCore(T);
        for (std.enums.values(Routine)) |r| {
            if (!types.isComplex(T) and (r == .hemm or r == .herk or r == .her2k)) continue;
            try routineShapes(T, &vt, Io.failing, r);
            if (comptime threading.hosted) for (limits) |limit| {
                var threaded: Io.Threaded = .init(std.testing.allocator, .{ .async_limit = limit, .concurrent_limit = limit });
                defer threaded.deinit();
                try routineShapes(T, &vt, threaded.io(), r);
            };
        }
    }
}

test "threaded gemm keeps the serial bits on a fiber Io, peers on one thread included" {
    if (comptime !threading.hosted or Io.Evented == void) return error.SkipZigTest;
    const Options = @typeInfo(@TypeOf(Io.Evented.init)).@"fn".params[2].type.?;
    const has_limit = @hasField(Options, "thread_limit");
    // no thread beside the caller's: a waiter that never slept would hang
    const thread_limits: []const ?usize = if (has_limit) &.{ 0, 3 } else &.{null};
    for (thread_limits) |tl| {
        var opts: Options = .{};
        if (has_limit) opts.thread_limit = tl;
        var ev: Io.Evented = undefined;
        ev.init(std.testing.allocator, opts) catch return error.SkipZigTest;
        defer ev.deinit();
        try shapes(f64, ev.io());
    }
}

test "threaded gemm matches the serial driver on the selected core's own blocking" {
    if (comptime !threading.hosted) return error.SkipZigTest;
    var threaded: Io.Threaded = .init(std.testing.allocator, .{});
    defer threaded.deinit();
    const vt = blas.vtable();
    inline for (precisions) |T| {
        const blk = Blocking.of(T, vt);
        try shape(T, vt, threaded.io(), .{ .n, .t }, 2 * blk.p + blk.um + 3, 3 * blk.un * 8 + 5, blk.q + 7, scalar(T, -0.75, 0.5), scalar(T, 0.5, 0));
        try shape(T, vt, threaded.io(), .{ .t, .n }, 3 * blk.p + 1, 97, 2 * blk.q + 1, scalar(T, 1.25, -1), scalar(T, 1, 0));
    }
}

test "the public gemm threads above the threshold and keeps the serial bits" {
    if (comptime !threading.hosted) return error.SkipZigTest;
    const saved = blas.config.cpu_count;
    defer blas.config.cpu_count = saved;
    var threaded: Io.Threaded = .init(std.testing.allocator, .{});
    defer threaded.deinit();
    const al = std.testing.allocator;
    inline for (precisions) |T| {
        // past `m * n * k` of 8 serial units for every type's rule
        const m = 211;
        const n = 197;
        const k = 173;
        const a = try al.alloc(T, m * k);
        defer al.free(a);
        const b = try al.alloc(T, k * n);
        defer al.free(b);
        const want = try al.alloc(T, m * n);
        defer al.free(want);
        const got = try al.alloc(T, m * n);
        defer al.free(got);
        fill(T, a, 5);
        fill(T, b, 6);
        fill(T, want, 7);
        @memcpy(got, want);
        const alpha = scalar(T, 0.5, 2);
        const beta = scalar(T, -1.5, 0.25);
        try driver.gemm(al, Io.failing, T, blas.vtable(), true, false, .n, m, n, k, alpha, a.ptr, k, b.ptr, k, beta, want.ptr, m, 1);
        for ([_]u32{ 1, 3, 8 }) |cpus| {
            blas.config.cpu_count = cpus;
            try zigstd.testExpectEqual(cpus, workers(.gemm, T, m * n * k));
            const c = try al.dupe(T, got);
            defer al.free(c);
            try interface.gemm(al, threaded.io(), T, .t, .n, m, n, k, alpha, a.ptr, k, b.ptr, k, beta, c.ptr, m);
            try zigstd.testExpectTrue(sameBits(T, want, c));
        }
    }
}

test "the public symm, herk and syr2k thread above their thresholds and keep the serial bits" {
    if (comptime !threading.hosted) return error.SkipZigTest;
    const saved = blas.config.cpu_count;
    defer blas.config.cpu_count = saved;
    blas.config.cpu_count = 5;
    var threaded: Io.Threaded = .init(std.testing.allocator, .{});
    defer threaded.deinit();
    const io = threaded.io();
    const al = std.testing.allocator;
    const T = types.Complex64;
    const n = 157;
    const k = 141;
    const a = try al.alloc(T, n * n);
    defer al.free(a);
    const b = try al.alloc(T, n * n);
    defer al.free(b);
    const want = try al.alloc(T, n * n);
    defer al.free(want);
    const got = try al.alloc(T, n * n);
    defer al.free(got);
    fill(T, a, 21);
    fill(T, b, 22);
    fill(T, want, 23);
    @memcpy(got, want);
    const vt = blas.vtable();
    const alpha = scalar(T, 0.5, -1);
    try zigstd.testExpectEqual(@as(u32, 5), workers(.symm, T, 2 * n * n * n));
    try driver.symm(al, Io.failing, T, vt, false, false, true, n, n, alpha, a.ptr, n, b.ptr, n, alpha, want.ptr, n, 1);
    try interface.symm(al, io, T, .r, .u, n, n, alpha, a.ptr, n, b.ptr, n, alpha, got.ptr, n);
    try zigstd.testExpectTrue(sameBits(T, want, got));
    try zigstd.testExpectEqual(@as(u32, 5), workers(.syrk, T, (n + 1) * n * k));
    try driver.syrk(al, Io.failing, T, vt, true, false, true, n, k, scalar(T, 0.25, 0), a.ptr, k, null, 0, scalar(T, 2, 0), want.ptr, n, 1);
    try interface.herk(al, io, T, .l, .c, n, k, 0.25, a.ptr, k, 2, got.ptr, n);
    try zigstd.testExpectTrue(sameBits(T, want, got));
    try zigstd.testExpectEqual(@as(u32, 5), workers(.syr2k, T, n * k));
    try driver.syrk(al, Io.failing, T, vt, false, true, false, n, k, alpha, a.ptr, n, b.ptr, n, alpha, want.ptr, n, 1);
    try interface.syr2k(al, io, T, .u, .n, n, k, alpha, a.ptr, n, b.ptr, n, alpha, got.ptr, n);
    try zigstd.testExpectTrue(sameBits(T, want, got));
}

test "ways use what they may and keep a panel per part" {
    try zigstd.testExpectEqual([2]usize{ 4, 4 }, ways(1000, 1000, 16, 8, 4));
    try zigstd.testExpectEqual([2]usize{ 1, 7 }, ways(8, 1000, 7, 8, 4));
    try zigstd.testExpectEqual([2]usize{ 7, 1 }, ways(1000, 4, 7, 8, 4));
    try zigstd.testExpectEqual([2]usize{ 1, 1 }, ways(1, 1, 9, 8, 4));
}

/// A kernel that holds every call until released, counting them and the
/// threads that make them.
const Hold = struct {
    const Kernel = @typeInfo(@FieldType(Vtable, "dgemm_kernel")).optional.child;
    const cc = types.kernelCc(f64);
    const BlasLong = types.BlasLong;
    const Atomic = std.atomic.Value(u32);

    var real: Kernel = undefined;
    var calls: Atomic = .init(0);
    var inside: Atomic = .init(0);
    var held: std.atomic.Value(bool) = .init(true);
    var first: std.atomic.Value(std.Thread.Id) = .init(0);
    var others: Atomic = .init(0);

    fn kernel(m: BlasLong, n: BlasLong, k: BlasLong, alpha: f64, a: [*]f64, b: [*]f64, c: [*]f64, ldc: BlasLong) callconv(cc) i32 {
        _ = calls.fetchAdd(1, .acq_rel);
        _ = inside.fetchAdd(1, .acq_rel);
        const id = std.Thread.getCurrentId();
        if (first.cmpxchgStrong(0, id, .acq_rel, .acquire)) |f| {
            if (f != id) _ = others.fetchAdd(1, .acq_rel);
        }
        defer _ = inside.fetchSub(1, .acq_rel);
        while (held.load(.acquire)) std.Thread.yield() catch {};
        return real(m, n, k, alpha, a, b, c, ldc);
    }

    fn release(io: Io) void {
        while (calls.load(.acquire) == 0) std.Thread.yield() catch {};
        io.sleep(.fromMilliseconds(50), .awake) catch {};
        held.store(false, .release);
    }
};

fn canceledGemm(io: Io, vt: *const Vtable, a: []f64, b: []f64, c: []f64, n: usize, k: usize) Error!void {
    return driver.gemm(std.testing.allocator, io, f64, vt, false, false, .n, n, n, k, 1, a.ptr, n, b.ptr, k, 1, c.ptr, n, 4);
}

test "a canceled threaded gemm returns error.Canceled with no worker left running" {
    if (comptime !threading.hosted) return error.SkipZigTest;
    var vt = shrunk(f64);
    Hold.real = vt.dgemm_kernel.?;
    vt.dgemm_kernel = &Hold.kernel;
    const n = 64;
    const k = 2000;
    const al = std.testing.allocator;
    const a = try al.alloc(f64, n * k);
    defer al.free(a);
    const b = try al.alloc(f64, k * n);
    defer al.free(b);
    const c = try al.alloc(f64, n * n);
    defer al.free(c);
    fill(f64, a, 1);
    fill(f64, b, 2);
    for ([_]Io.Limit{ .limited(3), .limited(8) }) |limit| {
        var threaded: Io.Threaded = .init(al, .{ .async_limit = limit, .concurrent_limit = limit });
        defer threaded.deinit();
        const io = threaded.io();
        Hold.calls.store(0, .release);
        Hold.held.store(true, .release);
        Hold.first.store(0, .release);
        Hold.others.store(0, .release);
        const releaser = try std.Thread.spawn(.{}, Hold.release, .{io});
        var future = io.concurrent(canceledGemm, .{ io, &vt, a, b, c, n, k }) catch return error.SkipZigTest;
        while (Hold.calls.load(.acquire) == 0) std.Thread.yield() catch {};
        const result = future.cancel(io);
        releaser.join();
        try zigstd.testExpectError(error.Canceled, result);
        // workers did run beside the caller
        try zigstd.testExpectTrue(Hold.others.load(.acquire) > 0);
        // nothing calls the kernel once the call has returned
        const after = Hold.calls.load(.acquire);
        try zigstd.testExpectEqual(@as(u32, 0), Hold.inside.load(.acquire));
        io.sleep(.fromMilliseconds(20), .awake) catch {};
        try zigstd.testExpectEqual(after, Hold.calls.load(.acquire));
    }
}

// `L3THREAD_BENCH` set, a ReleaseFast Windows build: medians of serial and
// the threaded driver (plain and fused) at 2, 4, 8 and 16 workers.
test "l3thread bench" {
    if (comptime sysinfo.target.os.tag != .windows or sysinfo.mode != @field(@TypeOf(sysinfo.mode), "ReleaseFast")) return error.SkipZigTest;
    if (std.process.Environ.getWindows(.{ .block = .global }, std.unicode.wtf8ToWtf16LeStringLiteral("L3THREAD_BENCH")) == null) return error.SkipZigTest;
    const al = std.heap.page_allocator;
    var threaded: Io.Threaded = .init(al, .{});
    defer threaded.deinit();
    const io = threaded.io();
    const vt = blas.vtable();
    std.debug.print("core {s}\n", .{blas.coreName()});
    // scratch from a buffer touched once, so no call pays for fresh pages
    const scratch = try al.alloc(u8, 256 << 20);
    defer al.free(scratch);
    @memset(scratch, 0);
    var fba: std.heap.FixedBufferAllocator = .init(scratch);
    inline for (.{ f64, f32 }) |T| for ([_]usize{ 256, 512, 1000, 2000 }) |n| {
        const a = try al.alloc(T, n * n);
        defer al.free(a);
        const b = try al.alloc(T, n * n);
        defer al.free(b);
        const c = try al.alloc(T, n * n);
        defer al.free(c);
        fill(T, a, 1);
        fill(T, b, 2);
        @memset(c, 0);
        const reps: usize = if (n >= 2000) 9 else if (n >= 1000) 15 else 25;
        const counts = [_]u32{ 2, 4, 8, 16 };
        // serial, then plain and fused at each count
        var times: [1 + 2 * counts.len][25]f64 = undefined;
        // the first round warms up and is not kept
        for (0..reps + 1) |round| for (0..times.len) |x| {
            const w: u32 = if (x == 0) 1 else counts[(x - 1) / 2];
            const fused = x != 0 and x % 2 == 0;
            fba.reset();
            const t0 = Io.Clock.awake.now(io).toNanoseconds();
            if (w == 1) {
                try driver.gemm(fba.allocator(), io, T, vt, false, false, .n, n, n, n, 1, a.ptr, n, b.ptr, n, 0, c.ptr, n, 1);
            } else if (fused) {
                try Team(T, true).run(fba.allocator(), io, vt, .{ .rows = .n, .cols = .n }, n, n, n, 1, .{ .p = a.ptr, .ld = n }, .{ .p = b.ptr, .ld = n }, null, 0, c.ptr, n, w);
            } else {
                try Team(T, false).run(fba.allocator(), io, vt, .{ .rows = .n, .cols = .n }, n, n, n, 1, .{ .p = a.ptr, .ld = n }, .{ .p = b.ptr, .ld = n }, null, 0, c.ptr, n, w);
            }
            if (round > 0) times[x][round - 1] = @floatFromInt(Io.Clock.awake.now(io).toNanoseconds() - t0);
        };
        const flops = 2.0 * @as(f64, @floatFromInt(n * n * n));
        std.debug.print("{s} n={d}: serial {d:.1}", .{ @typeName(T), n, flops / median(times[0][0..reps]) });
        for (counts, 0..) |w, i| std.debug.print(" | w={d} plain {d:.1} fused {d:.1}", .{ w, flops / median(times[1 + 2 * i][0..reps]), flops / median(times[2 + 2 * i][0..reps]) });
        std.debug.print(" GFLOPS\n", .{});
    };
}

fn median(x: []f64) f64 {
    std.mem.sort(f64, x, {}, std.sort.asc(f64));
    return x[x.len / 2];
}
