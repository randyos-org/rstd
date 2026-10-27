//! Two prototypes of threaded real `gemm`, compared before step 2 of
//! `threading.md`; internal, no public routine calls them. Both keep the
//! serial driver's depth slices (`Blocking.split`) and start every row range
//! on a row panel, so each element of `C` is computed by one worker in the
//! serial order and the bits are `driver.gemm`'s. Both wait on peers, so
//! they start through `threading.team` (`Group.concurrent`, all or none).
//!
//! - `.openblas`: a port of OpenBLAS's `driver/level3/level3_thread.c`
//!   (v0.3.34, BSD-3): the `CNAME` grid; each worker packs its share of its
//!   N group's `B` in two halves and publishes them through readiness flags;
//!   its peers multiply them with their own packed `A` and release them.
//! - `.blis`: BLIS's `jc x ic` ways (`docs/Multithreading.md`, tag 2.1),
//!   written from the docs: the `ic` workers of a `jc` group pack each slice
//!   of `B` together into one of two buffers and meet at one barrier per
//!   slice, where BLIS meets twice around one buffer (`bli_packm_int.c`).

const std = @import("std");
const sysinfo = @import("builtin");
const Vtable = @import("../Vtable.zig");
const driver = @import("driver.zig");
const threading = @import("../threading.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const Atomic = std.atomic.Value(u32);
const Blocking = driver.Blocking;
const Chunks = threading.Chunks;

pub const Design = enum { openblas, blis };

pub const Error = Allocator.Error || Io.Cancelable;

/// Flags of different workers sit on lines of their own.
const line = 64;

/// `C := alpha*op(A)*op(B) + beta*C` for real `T` on up to `want` workers
/// (at most `threading.chunks_max`), divided as `design` says. `m`, `n`
/// positive, as `driver.gemm` takes them.
pub fn gemm(allocator: Allocator, io: Io, comptime design: Design, comptime T: type, vt: *const Vtable, ta: bool, tb: bool, m: usize, n: usize, k: usize, alpha: T, a: [*]const T, lda: usize, b: [*]const T, ldb: usize, beta: T, c: [*]T, ldc: usize, want: u32) Error!void {
    const w = @min(want, threading.chunks_max);
    if (w <= 1 or k == 0 or driver.isZero(T, alpha)) return driver.gemm(allocator, io, T, vt, ta, tb, .n, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    const p: Problem(T) = .{
        .vt = vt,
        .d = .{ .rows = if (ta) .t else .n, .cols = if (tb) .t else .n },
        .blk = .of(T, vt),
        .m = m,
        .n = n,
        .k = k,
        .alpha = alpha,
        .beta = driver.unlessOne(T, beta),
        .x = .{ .p = @constCast(a), .ld = lda },
        .y = .{ .p = @constCast(b), .ld = ldb },
        .c = c,
        .ldc = ldc,
    };
    return switch (design) {
        .openblas => OpenBlas(T).run(allocator, io, &p, w),
        .blis => Blis(T).run(allocator, io, &p, w),
    };
}

fn ceilDiv(a: usize, b: usize) usize {
    return (a + b - 1) / b;
}

/// `x` rounded up to whole cache lines of `T`.
fn lines(comptime T: type, x: usize) usize {
    return std.mem.alignForward(usize, x, line / @sizeOf(T));
}

/// A chunk of `Chunks`, or an empty range at `len` past its count.
fn part(c: *const Chunks, i: usize, len: usize) threading.Chunk {
    return if (i < c.count) c.get(i) else .{ .index = i, .start = len, .end = len };
}

/// One call's operands and the loops both designs share.
fn Problem(comptime T: type) type {
    return struct {
        const Self = @This();

        vt: *const Vtable,
        d: driver.Desc,
        blk: Blocking,
        m: usize,
        n: usize,
        k: usize,
        alpha: T,
        beta: ?T,
        x: driver.Matrix(T),
        y: driver.Matrix(T),
        c: [*]T,
        ldc: usize,

        /// The deepest slice: `split` can pass `q` by its rounding.
        fn depth(p: *const Self) usize {
            return @min(p.k, p.blk.q + p.blk.um);
        }

        fn sliceAt(p: *const Self, ls: usize) usize {
            return Blocking.split(p.k - ls, p.blk.q, p.blk.um);
        }

        fn rowsAt(p: *const Self, is: usize, end: usize) usize {
            return Blocking.split(end - is, p.blk.p, p.blk.um);
        }

        /// Chunks of columns `[j0, j1)`, as `blk.chunk` steps them.
        fn chunkCount(p: *const Self, j0: usize, j1: usize) usize {
            var count: usize = 0;
            var jj = j0;
            while (jj < j1) : (count += 1) jj += p.blk.chunk(j1 - jj);
            return count;
        }

        /// Packs chunks `[c0, c1)` of columns `[j0, j1)` of `op(B)` over the
        /// slice at `ls`; a chunk at `jj` goes to `sb + min_l * (jj - j0)`.
        fn packB(p: *const Self, ls: usize, min_l: usize, j0: usize, j1: usize, sb: [*]T, c0: usize, c1: usize) void {
            var jj = j0;
            var i: usize = 0;
            while (jj < j1 and i < c1) : (i += 1) {
                const w = p.blk.chunk(j1 - jj);
                if (i >= c0) driver.packCols(T, p.vt, p.d, p.y, min_l, w, ls, jj, sb + min_l * (jj - j0));
                jj += w;
            }
        }

        /// The kernel on rows `[is, is + mi)` packed at `sa` and one chunk,
        /// columns `[jj, jj + w)` packed at `pb`; the first slice scales the
        /// block of `C` by `beta` first, as the serial driver does.
        fn tile(p: *const Self, ls: usize, min_l: usize, is: usize, mi: usize, jj: usize, w: usize, sa: [*]T, pb: [*]T) void {
            const cb = p.c + is + jj * p.ldc;
            if (ls == 0) if (p.beta) |s| driver.scale(T, p.vt, mi, w, s, cb, p.ldc);
            driver.kernel(T, p.vt, .n, mi, w, min_l, p.alpha, sa, pb, cb, p.ldc);
        }

        /// `tile` over every chunk of columns `[j0, j1)` packed at `sb`.
        fn multiply(p: *const Self, ls: usize, min_l: usize, is: usize, mi: usize, j0: usize, j1: usize, sa: [*]T, sb: [*]T) void {
            var jj = j0;
            while (jj < j1) {
                const w = p.blk.chunk(j1 - jj);
                p.tile(ls, min_l, is, mi, jj, w, sa, sb + min_l * (jj - j0));
                jj += w;
            }
        }
    };
}

// ---- (A) OpenBLAS: shared packed B, readiness flags ------------------------------

/// `CNAME`'s grid for `w` workers: rows split while each part keeps
/// `ratio` rows, columns to at most `ratio * tm` per part, then moved
/// toward square parts. Returns `tm`, `tn`.
fn grid(m: usize, n: usize, w: usize, ratio: usize) [2]usize {
    var tm: usize = w;
    if (m < 2 * ratio) tm = 1 else while (m < tm * ratio) tm /= 2;
    var tn: usize = 1;
    if (n >= ratio * tm) {
        tn = ceilDiv(n, ratio * tm);
        if (tm * tn > w) tn = w / tm;
        var cost: usize = 0;
        var div: usize = 0;
        var i: usize = 1;
        while (i * i <= tm) : (i += 1) {
            if (tm % i != 0) continue;
            const j = tm / i;
            const cost_i = n * j + m * tn * i;
            const cost_j = n * i + m * tn * j;
            if (cost == 0 or cost_i < cost) {
                cost = cost_i;
                div = i;
            }
            if (cost_j < cost) {
                cost = cost_j;
                div = j;
            }
        }
        if (div > 1) {
            tm /= div;
            tn *= div;
        }
    }
    return .{ tm, tn };
}

fn OpenBlas(comptime T: type) type {
    return struct {
        const Self = @This();
        const P = Problem(T);
        /// `DIVIDE_RATE`: halves of a worker's share, published one by one.
        const sides = 2;

        /// One half of one worker's share.
        const Slot = struct {
            /// The step it was last packed for.
            ready: Atomic align(line) = .init(0),
            /// Peers yet to release it.
            pending: Atomic align(line) = .init(0),
        };

        p: *const P,
        /// Per worker: a row block of `op(A)`, then the two halves of its
        /// share of `op(B)`.
        buf: [*]T,
        sa_len: usize,
        half_len: usize,
        /// Columns a half has room for.
        half_cols: usize,
        flags: [][sides]Slot,
        /// Columns per worker per round, so a share fits its halves.
        per: usize,
        tm: usize = 1,
        tn: usize = 1,
        rows: Chunks = undefined,

        fn run(allocator: Allocator, io: Io, p: *const P, want: u32) Error!void {
            const blk = p.blk;
            const depth = p.depth();
            const per = @min(blk.r, ceilDiv(p.n, want));
            // `Chunks` rounds each of two splits down to `un`
            const cap = @min(p.n, per + 2 * blk.un + 2);
            const half_cols = @min(cap, ceilDiv(cap, 2) + blk.un);
            const sa_len = lines(T, @min(p.m, blk.p) * depth);
            const half_len = lines(T, half_cols * depth);
            const buf = try allocator.alignedAlloc(T, .fromByteUnits(line), want * (sa_len + sides * half_len));
            defer allocator.free(buf);
            const flags = try allocator.alloc([sides]Slot, want);
            defer allocator.free(flags);
            var job: Self = .{ .p = p, .buf = buf.ptr, .sa_len = sa_len, .half_len = half_len, .half_cols = half_cols, .flags = flags, .per = per };
            try threading.team(io, &job, want);
        }

        pub fn plan(j: *Self, started: u32) u32 {
            const g = grid(j.p.m, j.p.n, started, @intCast(@max(1, j.p.vt.switch_ratio)));
            j.tm = g[0];
            j.tn = g[1];
            j.rows = .init(j.p.m, j.tm, .even, j.p.blk.um);
            const active = j.tm * j.tn;
            for (j.flags[0..active]) |*f| f.* = @splat(.{});
            return @intCast(active);
        }

        /// Worker `im` of group `jn`'s columns in the round `[js, js + width)`.
        fn share(j: *const Self, js: usize, width: usize, jn: usize, im: usize) [2]usize {
            const un = j.p.blk.un;
            const groups: Chunks = .init(width, j.tn, .even, un);
            const g = part(&groups, jn, width);
            const parts: Chunks = .init(g.end - g.start, j.tm, .even, un);
            const q = part(&parts, im, g.end - g.start);
            return .{ js + g.start + q.start, js + g.start + q.end };
        }

        /// Half `s` of a share, the first rounded up to `un`.
        fn half(j: *const Self, sh: [2]usize, s: usize) [2]usize {
            const un = j.p.blk.un;
            const mid = @min(sh[1], sh[0] + ceilDiv(ceilDiv(sh[1] - sh[0], 2), un) * un);
            return if (s == 0) .{ sh[0], mid } else .{ mid, sh[1] };
        }

        fn halfBuf(j: *const Self, worker: usize, s: usize) [*]T {
            return j.buf + worker * (j.sa_len + sides * j.half_len) + j.sa_len + s * j.half_len;
        }

        fn release(io: Io, f: *Slot) void {
            if (f.pending.fetchSub(1, .acq_rel) == 1) io.futexWake(u32, &f.pending.raw, 1);
        }

        pub fn work(j: *Self, io: Io, index: u32) Io.Cancelable!void {
            const p = j.p;
            const im = index % j.tm;
            const jn = index / j.tm;
            const g0 = jn * j.tm;
            const rows = part(&j.rows, im, p.m);
            const sa = j.buf + index * (j.sa_len + sides * j.half_len);
            // peers that read this worker's halves: the group's others with rows
            const readers: u32 = @intCast(j.rows.count - @intFromBool(im < j.rows.count));
            const round = j.per * j.tm * j.tn;
            var step: u32 = 0;
            var js: usize = 0;
            while (js < p.n) : (js += round) {
                const width = @min(round, p.n - js);
                const own = j.share(js, width, jn, im);
                var ls: usize = 0;
                var min_l: usize = 0;
                while (ls < p.k) : (ls += min_l) {
                    min_l = p.sliceAt(ls);
                    step += 1;
                    try threading.checkCancel(io);
                    var min_i: usize = if (rows.end > rows.start) p.rowsAt(rows.start, rows.end) else 0;
                    if (min_i > 0) driver.packRows(T, p.vt, p.d, p.x, min_l, min_i, ls, rows.start, sa);

                    // pack the own halves, each chunk multiplied as it is packed
                    for (0..sides) |s| {
                        const h = j.half(own, s);
                        if (h[0] == h[1]) continue;
                        std.debug.assert(h[1] - h[0] <= j.half_cols);
                        const f = &j.flags[index][s];
                        try threading.waitFor(io, &f.pending, 0);
                        const pb = j.halfBuf(index, s);
                        var jj = h[0];
                        while (jj < h[1]) {
                            const w = p.blk.chunk(h[1] - jj);
                            const at = pb + min_l * (jj - h[0]);
                            driver.packCols(T, p.vt, p.d, p.y, min_l, w, ls, jj, at);
                            if (min_i > 0) p.tile(ls, min_l, rows.start, min_i, jj, w, sa, at);
                            jj += w;
                        }
                        if (readers > 0) {
                            f.pending.store(readers, .monotonic);
                            f.ready.store(step, .release);
                            io.futexWake(u32, &f.ready.raw, std.math.maxInt(u32));
                        }
                    }
                    if (min_i == 0) continue;
                    const single = rows.start + min_i == rows.end;

                    // the peers' halves, from the next peer on
                    for (1..j.tm) |d| {
                        const q = (im + d) % j.tm;
                        const sh = j.share(js, width, jn, q);
                        for (0..sides) |s| {
                            const h = j.half(sh, s);
                            if (h[0] == h[1]) continue;
                            const f = &j.flags[g0 + q][s];
                            try threading.waitFor(io, &f.ready, step);
                            p.multiply(ls, min_l, rows.start, min_i, h[0], h[1], sa, j.halfBuf(g0 + q, s));
                            if (single) release(io, f);
                        }
                    }

                    // the later row blocks over every half of the group
                    var is = rows.start + min_i;
                    while (is < rows.end) : (is += min_i) {
                        try threading.checkCancel(io);
                        min_i = p.rowsAt(is, rows.end);
                        driver.packRows(T, p.vt, p.d, p.x, min_l, min_i, ls, is, sa);
                        const last = is + min_i == rows.end;
                        for (0..j.tm) |d| {
                            const q = (im + d) % j.tm;
                            const sh = j.share(js, width, jn, q);
                            for (0..sides) |s| {
                                const h = j.half(sh, s);
                                if (h[0] == h[1]) continue;
                                p.multiply(ls, min_l, is, min_i, h[0], h[1], sa, j.halfBuf(g0 + q, s));
                                if (last and q != im) release(io, &j.flags[g0 + q][s]);
                            }
                        }
                    }
                }
            }
        }
    };
}

// ---- (B) BLIS: jc x ic ways, cooperative packing, barriers -------------------------

/// `ic x jc` ways for `w` workers: the most workers that leave each part a
/// row panel and a column panel, then the squarest parts (least `n * ic +
/// m * jc`), more `ic` on a tie. Returns `ic`, `jc`.
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

fn Blis(comptime T: type) type {
    return struct {
        const Self = @This();
        const P = Problem(T);

        p: *const P,
        /// Per worker a row block of `op(A)`; per group two `B` buffers.
        sa: [*]T,
        sa_len: usize,
        sb: [*]T,
        /// Columns of `B` the call has room for, shared by its groups.
        budget: usize,
        sb_len: usize = 0,
        nc: usize = 0,
        ic: usize = 1,
        jc: usize = 1,
        rows: Chunks = undefined,
        cols: Chunks = undefined,
        barriers: []threading.Barrier,

        fn run(allocator: Allocator, io: Io, p: *const P, want: u32) Error!void {
            const blk = p.blk;
            const depth = p.depth();
            // the columns of `want` workers' shares, as `OpenBlas` has room for
            const budget = want * (@min(blk.r, ceilDiv(p.n, want)) + blk.un);
            const sa_len = lines(T, @min(p.m, blk.p) * depth);
            // two buffers per group, each on its own lines
            const sb_len = 2 * (budget * depth + want * line);
            const buf = try allocator.alignedAlloc(T, .fromByteUnits(line), want * sa_len + sb_len);
            defer allocator.free(buf);
            const barriers = try allocator.alloc(threading.Barrier, want);
            defer allocator.free(barriers);
            var job: Self = .{ .p = p, .sa = buf.ptr, .sa_len = sa_len, .sb = buf.ptr + want * sa_len, .budget = budget, .barriers = barriers };
            try threading.team(io, &job, want);
        }

        pub fn plan(j: *Self, started: u32) u32 {
            const blk = j.p.blk;
            const w = ways(j.p.m, j.p.n, started, blk.um, blk.un);
            j.ic = w[0];
            j.jc = w[1];
            j.rows = .init(j.p.m, j.ic, .even, blk.um);
            j.cols = .init(j.p.n, j.jc, .even, blk.un);
            // each group's share of the budget, whole column panels
            j.nc = @max(blk.un, @min(blk.r, j.budget / j.jc / blk.un * blk.un));
            j.sb_len = lines(T, j.nc * j.p.depth());
            for (j.barriers[0..j.jc]) |*b| b.* = .{ .parties = @intCast(j.ic) };
            return @intCast(j.ic * j.jc);
        }

        pub fn work(j: *Self, io: Io, index: u32) Io.Cancelable!void {
            const p = j.p;
            const r = index % j.ic;
            const g = index / j.ic;
            const cols = part(&j.cols, g, p.n);
            const rows = part(&j.rows, r, p.m);
            const sa = j.sa + index * j.sa_len;
            const bufs = j.sb + g * 2 * j.sb_len;
            const bar = &j.barriers[g];
            var t: usize = 0;
            var j0 = cols.start;
            while (j0 < cols.end) : (j0 += j.nc) {
                const j1 = @min(j0 + j.nc, cols.end);
                const chunks = p.chunkCount(j0, j1);
                var ls: usize = 0;
                var min_l: usize = 0;
                while (ls < p.k) : (ls += min_l) {
                    min_l = p.sliceAt(ls);
                    try threading.checkCancel(io);
                    // the barrier before last let every peer finish with this buffer
                    const sb = bufs + (t % 2) * j.sb_len;
                    t += 1;
                    p.packB(ls, min_l, j0, j1, sb, chunks * r / j.ic, chunks * (r + 1) / j.ic);
                    try bar.wait(io);
                    var is = rows.start;
                    var min_i: usize = 0;
                    while (is < rows.end) : (is += min_i) {
                        min_i = p.rowsAt(is, rows.end);
                        driver.packRows(T, p.vt, p.d, p.x, min_l, min_i, ls, is, sa);
                        p.multiply(ls, min_l, is, min_i, j0, j1, sa, sb);
                    }
                }
            }
        }
    };
}

// ---- tests -------------------------------------------------------------------------

const zigstd = @import("../../../zigstd.zig");
const types = @import("../types.zig");
const blas = @import("../__root__.zig");
const gemm_kernels = @import("gemm/kernel.zig");
const guard = @import("../_test/blas2/guard.zig");

const designs = [_]Design{ .openblas, .blis };

/// Mixed magnitudes, so a sum's bits follow its order.
fn fill(comptime T: type, x: []T, seed: u64) void {
    var prng: std.Random.DefaultPrng = .init(seed);
    const rand = prng.random();
    for (x) |*v| v.* = std.math.ldexp(rand.float(T) * 2 - 1, rand.intRangeAtMost(i32, -12, 12));
}

fn sameBits(comptime T: type, a: []const T, b: []const T) bool {
    return std.mem.eql(u8, std.mem.sliceAsBytes(a), std.mem.sliceAsBytes(b));
}

/// The portable kernels, blocking shrunk as `check.zig` shrinks it, so every
/// block and panel edge comes at small sizes.
fn shrunk(comptime T: type) Vtable {
    var vt: Vtable = .{};
    gemm_kernels.installGeneric(&vt);
    const blk = Blocking.of(T, &vt);
    const pre = if (T == f32) "sgemm_" else "dgemm_";
    @field(vt, pre ++ "p") = @intCast(2 * blk.um);
    @field(vt, pre ++ "q") = 5;
    @field(vt, pre ++ "r") = @intCast(6 * blk.un + 1);
    vt.switch_ratio = 2;
    return vt;
}

/// One shape against `driver.gemm`, bit for bit, at 1 to 9 workers of each
/// design on `io`. Operand padding holds NaN; `C`'s, a sentinel that must
/// survive; scratch comes poisoned.
fn shape(comptime T: type, vt: *const Vtable, io: Io, ta: bool, tb: bool, m: usize, n: usize, k: usize, alpha: T, beta: T) !void {
    const al = std.testing.allocator;
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
    for (0..ca) |j| @memset(a[j * (ra + pad) + ra ..][0..pad], std.math.nan(T));
    for (0..cb) |j| @memset(b[j * (rb + pad) + rb ..][0..pad], std.math.nan(T));
    for (0..n) |j| @memset(c0[j * (m + pad) + m ..][0..pad], 1234.5);
    @memset(c0[n * (m + pad) ..], 1234.5);

    var poison: guard.Poison = .{ .child = al };
    const pa = poison.allocator();
    @memcpy(want, c0);
    try driver.gemm(pa, Io.failing, T, vt, ta, tb, .n, m, n, k, alpha, a.ptr, ra + pad, b.ptr, rb + pad, beta, want.ptr, m + pad);
    inline for (designs) |design| for (1..10) |w| {
        @memcpy(got, c0);
        try gemm(pa, io, design, T, vt, ta, tb, m, n, k, alpha, a.ptr, ra + pad, b.ptr, rb + pad, beta, got.ptr, m + pad, @intCast(w));
        if (!sameBits(T, want, got)) {
            std.debug.print("{t} {s} w={d} ta={} tb={} m={d} n={d} k={d}\n", .{ design, @typeName(T), w, ta, tb, m, n, k });
            return error.TestUnexpectedResult;
        }
    };
}

/// Every transpose over shapes around every block and panel edge.
fn shapes(comptime T: type, io: Io) !void {
    const vt = shrunk(T);
    const blk = Blocking.of(T, &vt);
    const um = blk.um;
    const un = blk.un;
    const rr = blk.r;
    const ms = [_]usize{ 1, um + 1, 2 * blk.p + um + 1, 4 * blk.p + 3, 9 * um + 5 };
    const ns = [_]usize{ 1, un + 1, 3 * un + 1, rr + 1, 2 * rr + 3, 5 * rr };
    const ks = [_]usize{ 1, blk.q + 1, 2 * blk.q + 3 };
    for ([_][2]bool{ .{ false, false }, .{ true, false }, .{ false, true }, .{ true, true } }) |t| {
        for (ms, 0..) |m, mi| for (ns, 0..) |n, ni| for (ks, 0..) |k, ki| {
            if ((mi + ni + ki) % 2 != 0) continue;
            try shape(T, &vt, io, t[0], t[1], m, n, k, -0.75, 0.5);
        };
        // beta 0 and 1, which the first slice folds differently
        try shape(T, &vt, io, t[0], t[1], 9 * um + 5, 2 * rr + 3, 2 * blk.q + 3, 1.5, 0);
        try shape(T, &vt, io, t[0], t[1], 9 * um + 5, 2 * rr + 3, 2 * blk.q + 3, 1, 1);
    }
}

/// `async_limit` and `concurrent_limit` of the threaded tests: none, fewer
/// than asked, oversubscribed.
const limits = [_]Io.Limit{ .nothing, .limited(3), .limited(8) };

test "prototypes match the serial driver bit for bit at 1 to 9 workers, on every Io" {
    inline for (.{ f32, f64 }) |T| {
        try shapes(T, Io.failing);
        if (comptime threading.hosted) for (limits) |limit| {
            var threaded: Io.Threaded = .init(std.testing.allocator, .{ .async_limit = limit, .concurrent_limit = limit });
            defer threaded.deinit();
            try shapes(T, threaded.io());
        };
    }
}

test "prototypes match the serial driver on the selected core's own blocking" {
    if (comptime !threading.hosted) return error.SkipZigTest;
    var threaded: Io.Threaded = .init(std.testing.allocator, .{});
    defer threaded.deinit();
    const vt = blas.vtable();
    inline for (.{ f32, f64 }) |T| {
        const blk = Blocking.of(T, vt);
        try shape(T, vt, threaded.io(), false, true, 2 * blk.p + blk.um + 3, 3 * blk.un * 8 + 5, blk.q + 7, -0.75, 0.5);
        try shape(T, vt, threaded.io(), true, false, 3 * blk.p + 1, 97, 2 * blk.q + 1, 1.25, 1);
    }
}

test "grid and ways use what they may and keep a panel per part" {
    try zigstd.testExpectEqual([2]usize{ 1, 8 }, grid(3, 1000, 8, 2));
    try zigstd.testExpectEqual([2]usize{ 1, 1 }, grid(1, 1, 8, 2));
    const g = grid(1000, 1000, 16, 4);
    try zigstd.testExpectTrue(g[0] * g[1] <= 16 and g[0] * g[1] >= 8);
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

fn canceledGemm(io: Io, design: Design, vt: *const Vtable, a: []f64, b: []f64, c: []f64, n: usize, k: usize) Error!void {
    return switch (design) {
        inline else => |d| gemm(std.testing.allocator, io, d, f64, vt, false, false, n, n, k, 1, a.ptr, n, b.ptr, k, 1, c.ptr, n, 4),
    };
}

test "a canceled call returns error.Canceled with no worker left running" {
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
    for ([_]Io.Limit{ .limited(3), .limited(8) }) |limit| for (designs) |design| {
        var threaded: Io.Threaded = .init(al, .{ .async_limit = limit, .concurrent_limit = limit });
        defer threaded.deinit();
        const io = threaded.io();
        Hold.calls.store(0, .release);
        Hold.held.store(true, .release);
        Hold.first.store(0, .release);
        Hold.others.store(0, .release);
        const releaser = try std.Thread.spawn(.{}, Hold.release, .{io});
        var future = io.concurrent(canceledGemm, .{ io, design, &vt, a, b, c, n, k }) catch return error.SkipZigTest;
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
    };
}

// `L3PROTO_BENCH` set, a ReleaseFast Windows build: medians of serial and
// both designs at 1, 2, 4, 8 and 16 workers.
test "l3proto bench" {
    if (comptime sysinfo.target.os.tag != .windows or sysinfo.mode != @field(@TypeOf(sysinfo.mode), "ReleaseFast")) return error.SkipZigTest;
    if (std.process.Environ.getWindows(.{ .block = .global }, std.unicode.wtf8ToWtf16LeStringLiteral("L3PROTO_BENCH")) == null) return error.SkipZigTest;
    const al = std.heap.page_allocator;
    var threaded: Io.Threaded = .init(al, .{});
    defer threaded.deinit();
    const io = threaded.io();
    const vt = blas.vtable();
    std.debug.print("core {s}\n", .{blas.coreName()});
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
        const reps: usize = if (n >= 2000) 5 else if (n >= 1000) 9 else 15;
        const counts = [_]u32{ 2, 4, 8, 16 };
        // serial, then each design at each count
        var times: [1 + 2 * counts.len][15]f64 = undefined;
        for (0..reps) |r| for (0..times.len) |x| {
            const design: Design = if (x == 0 or x % 2 == 1) .openblas else .blis;
            const w: u32 = if (x == 0) 1 else counts[(x - 1) / 2];
            const t0 = Io.Clock.awake.now(io).toNanoseconds();
            switch (design) {
                inline else => |d| try gemm(al, io, d, T, vt, false, false, n, n, n, 1, a.ptr, n, b.ptr, n, 0, c.ptr, n, w),
            }
            times[x][r] = @floatFromInt(Io.Clock.awake.now(io).toNanoseconds() - t0);
        };
        const flops = 2.0 * @as(f64, @floatFromInt(n * n * n));
        std.debug.print("{s} n={d}: serial {d:.1}", .{ @typeName(T), n, flops / median(times[0][0..reps]) });
        for (counts, 0..) |w, i| std.debug.print(" | w={d} A {d:.1} B {d:.1}", .{ w, flops / median(times[1 + 2 * i][0..reps]), flops / median(times[2 + 2 * i][0..reps]) });
        std.debug.print(" GFLOPS\n", .{});
    };
}

fn median(x: []f64) f64 {
    std.mem.sort(f64, x, {}, std.sort.asc(f64));
    return x[x.len / 2];
}
