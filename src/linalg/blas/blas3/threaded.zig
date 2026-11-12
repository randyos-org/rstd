//! Threaded level 3 in BLIS's model (IPDPS 2014, <https://doi.org/10.1109/IPDPS.2014.110>;
//! BLIS 2.1's `docs/Multithreading.md`): `C` split `ic x jc` ways, the serial
//! loops, packing and kernel calls, so the same bits at any worker count.

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
    const r = threading.rule(family, T, blas.coreType());
    return @min(threading.workers(&blas.config, r.sizeMax(work)), threading.chunks_max);
}

/// `driver.run`'s loops on up to `want` workers.
pub fn run(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, d: Desc, m: usize, n: usize, k: usize, alpha: T, x: Matrix(T), y: Matrix(T), twice: ?T, beta: ?T, c: [*]T, ldc: usize, want: u32) Error!void {
    return Team(T).run(allocator, io, vt, d, m, n, k, alpha, x, y, twice, beta, c, ldc, @min(want, threading.chunks_max));
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

/// One call's shared state. A worker multiplies its own chunks of `B` with
/// its first row block as it packs them, while they are in cache, as
/// `level3_thread.c` does; measured even with packing first, or up to 4% ahead.
fn Team(comptime T: type) type {
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

        pub fn work(j: *Self, io: Io, index: u32, lent: ?*threading.Lent) Io.Cancelable!void {
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
                        if (!first) try bar.wait(io, lent);
                        first = false;
                        var min_i: usize = if (rows.end > rows.start) Blocking.split(rows.end - rows.start, blk.p, blk.um) else 0;
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
                        try bar.wait(io, lent);
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
const check = @import("../_test/check.zig");

const fill = check.fillMixed;

test "ways use what they may and keep a panel per part" {
    try zigstd.testExpectEqual([2]usize{ 4, 4 }, ways(1000, 1000, 16, 8, 4));
    try zigstd.testExpectEqual([2]usize{ 1, 7 }, ways(8, 1000, 7, 8, 4));
    try zigstd.testExpectEqual([2]usize{ 7, 1 }, ways(1000, 4, 7, 8, 4));
    try zigstd.testExpectEqual([2]usize{ 1, 1 }, ways(1, 1, 9, 8, 4));
}

// `L3THREAD_BENCH` set, a fast Windows build (`rstd-test -Doptimize=fast`):
// medians of serial and threaded `gemm` at 2, 4, 8 and 16 workers.
test "l3thread bench" {
    if (comptime sysinfo.target.os.tag != .windows or sysinfo.mode != .fast) return error.SkipZigTest;
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
        // serial, then each count
        var times: [1 + counts.len][25]f64 = undefined;
        // the first round warms up and is not kept
        for (0..reps + 1) |round| for (0..times.len) |x| {
            fba.reset();
            const t0 = Io.Clock.awake.now(io).toNanoseconds();
            try driver.gemm(fba.allocator(), io, T, vt, false, false, .n, n, n, n, 1, a.ptr, n, b.ptr, n, 0, c.ptr, n, if (x == 0) 1 else counts[x - 1]);
            if (round > 0) times[x][round - 1] = @floatFromInt(Io.Clock.awake.now(io).toNanoseconds() - t0);
        };
        const flops = 2.0 * @as(f64, @floatFromInt(n * n * n));
        std.debug.print("{s} n={d}: serial {d:.1}", .{ @typeName(T), n, flops / median(times[0][0..reps]) });
        for (counts, 0..) |w, i| std.debug.print(" | w={d} {d:.1}", .{ w, flops / median(times[1 + i][0..reps]) });
        std.debug.print(" GFLOPS\n", .{});
    };
}

fn median(x: []f64) f64 {
    std.mem.sort(f64, x, {}, std.sort.asc(f64));
    return x[x.len / 2];
}
