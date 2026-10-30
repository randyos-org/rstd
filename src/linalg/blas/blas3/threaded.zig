//! Threaded level 3 in BLIS's model (Smith, van de Geijn, Smelyanskiy,
//! Hammond, Van Zee, "Anatomy of High-Performance Many-Threaded Matrix
//! Multiplication", IPDPS 2014, 1049-1059, <https://doi.org/10.1109/IPDPS.2014.110>;
//! BLIS 2.1's `docs/Multithreading.md`), written from its design. Workers split `C`
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
    return Team(T).run(allocator, io, vt, d, m, n, k, alpha, x, y, twice, beta, c, ldc, @min(want, threading.chunks_max));
}

/// `driver.trmm` or, with `solve`, `trsm` on up to `want` workers: column
/// panels of `B` on the left, row panels on the right, as OpenBLAS's
/// `gemm_thread_n` and `gemm_thread_m` divide it and BLIS splits its `jc` and
/// `ic` loops. Each part is an independent serial problem on packing space
/// of its own, so they run through `threading.forEach`; parts start on
/// `unroll_n` columns or `unroll_m` rows, so tiles are the serial ones.
pub fn triangular(allocator: Allocator, io: Io, comptime T: type, vt: *const Vtable, comptime solve: bool, o: driver.Tri, blk: Blocking, m: usize, n: usize, a: [*]T, lda: usize, s: ?T, b: [*]T, ldb: usize, want: u32) Error!void {
    const W = driver.Work(T);
    const Parts = struct {
        vt: *const Vtable,
        o: driver.Tri,
        blk: Blocking,
        m: usize,
        n: usize,
        a: [*]T,
        lda: usize,
        s: ?T,
        b: [*]T,
        ldb: usize,
        buf: [*]T,
        /// Each part's packing space, from `buf`.
        at: [threading.chunks_max + 1]usize,

        /// Part `c`'s `m`, `n` and `B`.
        fn of(p: *const @This(), c: threading.Chunk) struct { usize, usize, [*]T } {
            return if (p.o.left) .{ p.m, c.end - c.start, p.b + c.start * p.ldb } else .{ c.end - c.start, p.n, p.b + c.start };
        }

        fn lens(p: *const @This(), c: threading.Chunk) [3]usize {
            const d = p.of(c);
            return driver.triLens(p.blk, d[0], d[1], if (p.o.left) p.m else p.n);
        }

        fn one(p: *const @This(), c: threading.Chunk) void {
            const d = p.of(c);
            driver.triDrive(T, p.vt, solve, p.o, p.blk, d[0], d[1], p.a, p.lda, p.s, d[2], p.ldb, W.at(p.buf + p.at[c.index], p.lens(c)));
        }
    };
    const plan: threading.Plan = .{
        .chunks = .init(if (o.left) n else m, want, .even, if (o.left) blk.un else blk.um),
        .workers = want,
    };
    var parts: Parts = .{ .vt = vt, .o = o, .blk = blk, .m = m, .n = n, .a = a, .lda = lda, .s = s, .b = b, .ldb = ldb, .buf = undefined, .at = undefined };
    parts.at[0] = 0;
    for (0..plan.chunks.count) |i| parts.at[i + 1] = parts.at[i] + W.total(parts.lens(plan.chunks.get(i)));
    const buf = try allocator.alignedAlloc(T, .fromByteUnits(W.alignment), parts.at[plan.chunks.count]);
    defer allocator.free(buf);
    parts.buf = buf.ptr;
    const ctx: *const Parts = &parts;
    try threading.forEach(io, &plan, ctx, Parts.one);
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
const util = @import("../util.zig");

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

/// The depths the level 3 kernels see, a bit each (past 255 in the last):
/// `Q`'s slices, which must not follow the worker count (threading.md,
/// Guardrails).
const Depths = struct {
    // 32-bit words: 32-bit targets may lack 64-bit atomics (powerpc-power9)
    var seen: [8]std.atomic.Value(u32) = @splat(.init(0));

    fn note(k: types.BlasLong) void {
        const d: usize = @intCast(std.math.clamp(k, 0, 255));
        _ = seen[d / 32].fetchOr(@as(u32, 1) << @intCast(d % 32), .monotonic);
    }

    /// Whether the depths seen since the last `take` are `want`; clears them.
    fn same(want: [8]u32) bool {
        return std.mem.eql(u32, &want, &take());
    }

    /// The depths seen since the last `take`.
    fn take() [8]u32 {
        var out: [8]u32 = undefined;
        for (&seen, &out) |*s, *o| o.* = s.swap(0, .monotonic);
        return out;
    }

    /// Kernel `name` noting its depth, its third argument, then calling `real`.
    fn Spy(comptime name: []const u8) type {
        const F = @typeInfo(@FieldType(Vtable, name)).optional.child;
        const Fn = @typeInfo(F).pointer.child;
        const A = std.meta.ArgsTuple(Fn);
        const R = @typeInfo(Fn).@"fn".return_type.?;
        const cc = types.kernelCc(f64);
        // Zig 0.17 names it `fields`
        const info = @typeInfo(A).@"struct";
        const n = if (@hasField(@TypeOf(info), "field_names")) info.field_names.len else info.fields.len;
        return struct {
            var real: F = undefined;
            fn t(comptime i: usize) type {
                return @FieldType(A, std.fmt.comptimePrint("{d}", .{i}));
            }
            const wrap: F = switch (n) {
                8 => &struct {
                    fn f(a0: t(0), a1: t(1), k: t(2), a3: t(3), a4: t(4), a5: t(5), a6: t(6), a7: t(7)) callconv(cc) R {
                        note(k);
                        return real(a0, a1, k, a3, a4, a5, a6, a7);
                    }
                }.f,
                9 => &struct {
                    fn f(a0: t(0), a1: t(1), k: t(2), a3: t(3), a4: t(4), a5: t(5), a6: t(6), a7: t(7), a8: t(8)) callconv(cc) R {
                        note(k);
                        return real(a0, a1, k, a3, a4, a5, a6, a7, a8);
                    }
                }.f,
                10 => &struct {
                    fn f(a0: t(0), a1: t(1), k: t(2), a3: t(3), a4: t(4), a5: t(5), a6: t(6), a7: t(7), a8: t(8), a9: t(9)) callconv(cc) R {
                        note(k);
                        return real(a0, a1, k, a3, a4, a5, a6, a7, a8, a9);
                    }
                }.f,
                else => @compileError(name),
            };
        };
    }

    /// `vt` with `T`'s `gemm`, `trsm` and `trmm` kernels noting their depths.
    fn spied(comptime T: type, vt: Vtable) Vtable {
        const l = comptime portable.letter(T).?;
        const complex = comptime types.isComplex(T);
        const sides = .{ "LN", "LT", "RN", "RT" } ++ if (complex) .{ "LR", "LC", "RR", "RC" } else .{};
        const gemms = if (complex) .{ "gemm_kernel_n", "gemm_kernel_l", "gemm_kernel_r", "gemm_kernel_b" } else .{"gemm_kernel"};
        var out = vt;
        inline for (gemms) |s| patch(&out, l ++ s);
        inline for (sides) |s| {
            patch(&out, l ++ "trsm_kernel_" ++ s);
            patch(&out, l ++ "trmm_kernel_" ++ s);
        }
        return out;
    }

    fn patch(vt: *Vtable, comptime name: []const u8) void {
        if (@field(vt, name)) |real| {
            Spy(name).real = real;
            @field(vt, name) = Spy(name).wrap;
        }
    }
};

/// One shape against `driver.gemm`, bit for bit, at 1 to 9 workers on `io`,
/// and with the kernels' depths if `vt` is `Depths.spied`. Operand padding
/// holds NaN; `C`'s, a sentinel that must survive; scratch comes poisoned.
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
    const depths = Depths.take();
    for (1..10) |w| {
        @memcpy(got, c0);
        try driver.gemm(pa, io, T, vt, ta, tb, conj, m, n, k, alpha, a.ptr, ra + pad, b.ptr, rb + pad, beta, got.ptr, m + pad, @intCast(w));
        if (!sameBits(T, want, got) or !Depths.same(depths)) {
            std.debug.print("{s} w={d} {t}{t} m={d} n={d} k={d}\n", .{ @typeName(T), w, tr[0], tr[1], m, n, k });
            return error.TestUnexpectedResult;
        }
    }
}

/// Every `op` pair over shapes around every block and panel edge (complex,
/// fewer shapes over its sixteen pairs).
fn shapes(comptime T: type, io: Io) !void {
    const vt = Depths.spied(T, shrunk(T));
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

test "threaded gemm matches the serial driver's bits and kernel depths at 1 to 9 workers, on every Io" {
    inline for (precisions) |T| {
        try shapes(T, Io.failing);
        if (comptime threading.hosted) for (limits) |limit| {
            var threaded: Io.Threaded = .init(std.testing.allocator, .{ .async_limit = limit, .concurrent_limit = limit });
            defer threaded.deinit();
            try shapes(T, threaded.io());
        };
        if (comptime threading.evented_workers) {
            var ev: threading.Evented = undefined;
            try zigstd.testExpectTrue(threading.initEvented(&ev, 3));
            defer ev.deinit();
            try shapes(T, ev.io());
        }
    }
}

/// The threaded routines besides `gemm`.
const Routine = enum { symm, hemm, syrk, herk, syr2k, her2k, trmm, trsm };

/// `r` on `vt` at up to `want` workers: `C` (`trmm`'s and `trsm`'s `B`) is
/// `m x n` (`n x n` for the rank updates), the depth `k`; `opt` is the side
/// or `trans`, then `upper`, then for `trmm` and `trsm` `trans` (and for
/// complex, conjugating).
fn routine(comptime T: type, vt: *const Vtable, al: Allocator, io: Io, r: Routine, opt: [3]bool, m: usize, n: usize, k: usize, alpha: T, beta: T, a: [*]const T, b: [*]const T, ld: usize, c: [*]T, want: u32) !void {
    const re = if (comptime types.isComplex(T)) scalar(T, alpha.re, 0) else alpha;
    switch (r) {
        .trmm => try driver.trmm(al, io, T, vt, opt[0], opt[1], opt[2], opt[2], false, m, n, alpha, a, ld, c, ld, want),
        .trsm => try driver.trsm(al, io, T, vt, opt[0], opt[1], opt[2], opt[2], false, m, n, alpha, a, ld, c, ld, want),
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

/// A tame triangle for `trsm` in `a` at `ld`: small off the diagonal, at
/// least 1 on it.
fn tame(comptime T: type, a: []T, ld: usize) void {
    for (0..ld) |j| for (0..ld) |i| {
        const v = &a[i + j * ld];
        const f = scalar(T, if (i == j) 1.5 else 0.25 / @as(f64, @floatFromInt(ld)), 0);
        v.* = if (comptime types.isComplex(T)) v.mul(f).add(if (i == j) scalar(T, 1, 0) else scalar(T, 0, 0)) else v.* * f + @as(T, if (i == j) 1 else 0);
    };
}

/// `r`'s shapes around every block, panel and triangle edge, all options,
/// against the serial loops bit for bit at 1 to 9 workers on `io`, and with
/// the kernels' depths if `vt` is `Depths.spied`.
fn routineShapes(comptime T: type, vt: *const Vtable, io: Io, r: Routine) !void {
    const al = std.testing.allocator;
    const blk = Blocking.of(T, vt);
    const herm = r == .herk or r == .her2k;
    const rank = r == .syrk or r == .herk or r == .syr2k or r == .her2k;
    const tri = r == .trmm or r == .trsm;
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
    if (tri) tame(T, a, ld);
    var poison: guard.Poison = .{ .child = al };
    const pa = poison.allocator();
    const alpha = scalar(T, -0.75, 0.25);
    const beta = if (herm) scalar(T, 0.5, 0) else scalar(T, 0.5, -0.125);
    var x: usize = 0;
    for (0..8) |bits| {
        const opt: [3]bool = .{ bits & 1 != 0, bits & 2 != 0, bits & 4 != 0 };
        if (!tri and opt[2]) continue;
        for (ns) |m| for (ns) |n| for (ks) |k| {
            // the rank updates' `C` is square; the others have no depth
            if ((rank and m != n) or (!rank and k != 1)) continue;
            x += 1;
            if (x % 2 != 0) continue;
            @memcpy(want, c0);
            try routine(T, vt, pa, Io.failing, r, opt, m, n, k, alpha, beta, a.ptr, b.ptr, ld, want.ptr, 1);
            const depths = Depths.take();
            for (1..10) |w| {
                @memcpy(got, c0);
                try routine(T, vt, pa, io, r, opt, m, n, k, alpha, beta, a.ptr, b.ptr, ld, got.ptr, @intCast(w));
                if (!sameBits(T, want, got) or !Depths.same(depths)) {
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

test "threaded symm, hemm, rank-k updates, trmm and trsm match the serial loops' bits and kernel depths at 1 to 9 workers" {
    inline for (precisions) |T| {
        const vt = Depths.spied(T, shrunkCore(T));
        for (std.enums.values(Routine)) |r| {
            if (!types.isComplex(T) and (r == .hemm or r == .herk or r == .her2k)) continue;
            try routineShapes(T, &vt, Io.failing, r);
            if (comptime threading.hosted) for (limits) |limit| {
                var threaded: Io.Threaded = .init(std.testing.allocator, .{ .async_limit = limit, .concurrent_limit = limit });
                defer threaded.deinit();
                try routineShapes(T, &vt, threaded.io(), r);
            };
            if (comptime threading.evented_workers) {
                var ev: threading.Evented = undefined;
                try zigstd.testExpectTrue(threading.initEvented(&ev, 3));
                defer ev.deinit();
                try routineShapes(T, &vt, ev.io(), r);
            }
        }
    }
}

test "threaded gemm keeps the serial bits under the caller's rounding toward zero and flush to zero" {
    if (comptime !threading.hosted or !util.FpState.governs) return error.SkipZigTest;
    const own: util.FpState = .capture();
    defer own.apply();
    var threaded: Io.Threaded = .init(std.testing.allocator, .{});
    defer threaded.deinit();
    try fpShapes(own, threaded.io());
    if (comptime threading.evented_workers) {
        var ev: threading.Evented = undefined;
        try zigstd.testExpectTrue(threading.initEvented(&ev, 3));
        defer ev.deinit();
        try fpShapes(own, ev.io());
    }
}

/// `gemm` on `io` under `own`, then rounding toward zero.
fn fpShapes(own: util.FpState, io: Io) !void {
    inline for (precisions) |T| {
        const vt = shrunk(T);
        const blk = Blocking.of(T, &vt);
        // the workers start under the default state first
        for ([_]util.FpState{ own, own.towardZero() }) |s| {
            s.apply();
            try shape(T, &vt, io, .{ .n, .t }, 9 * blk.um + 5, 2 * blk.r + 3, 2 * blk.q + 3, scalar(T, -0.75, 0.25), scalar(T, 0.5, -0.125));
            try zigstd.testExpectTrue(util.FpState.capture().eql(s));
        }
    }
}

test "threaded gemm keeps the serial bits on a fiber Io with every peer on the caller's thread" {
    const Evented = threading.Evented;
    if (comptime !threading.hosted or Evented == void) return error.SkipZigTest;
    // one thread, so a waiter that never yielded would hang; with worker
    // threads std's `Uring` and `Dispatch` crash in `Group` (issues.md)
    var ev: Evented = undefined;
    if (!threading.initEvented(&ev, 0)) return error.SkipZigTest;
    defer ev.deinit();
    const vt = shrunk(f64);
    const blk = Blocking.of(f64, &vt);
    // `ic` and `jc` both past 1, several slices
    try shape(f64, &vt, ev.io(), .{ .n, .t }, 9 * blk.um + 5, 2 * blk.r + 3, 2 * blk.q + 3, -0.75, 0.5);
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
    _ = blas.coreName();
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

test "the public symm, herk, syr2k and trmm thread above their thresholds and keep the serial bits" {
    if (comptime !threading.hosted) return error.SkipZigTest;
    // selecting the core fills `blas.config` once; override after
    _ = blas.coreName();
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
    try zigstd.testExpectEqual(@as(u32, 5), workers(.trsm, T, n * k));
    try driver.trmm(al, Io.failing, T, vt, false, false, true, true, true, k, n, alpha, a.ptr, n, want.ptr, n, 1);
    try interface.trmm(al, io, T, .r, .l, .c, .u, k, n, alpha, a.ptr, n, got.ptr, n);
    try zigstd.testExpectTrue(sameBits(T, want, got));
}

/// `r` (`gemm` if null) on `io` at 2 to 8 workers, round after round,
/// against one worker's bits on `Io.failing`. `C` (`B`) is `m x n` at
/// `ld`; every element past it is a guard band that must keep its value.
/// Scratch comes NaN-poisoned, so a read before its packing shows.
fn stress(comptime T: type, io: Io, r: ?Routine, opt: [3]bool, m: usize, n: usize, k: usize, rounds: usize) !void {
    const al = std.testing.allocator;
    const vt = shrunkCore(T);
    const ld = @max(m, n, k) + 3;
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
    fill(T, a, 41);
    fill(T, b, 42);
    fill(T, c0, 43);
    if (r == .trsm) tame(T, a, ld);
    var poison: guard.Poison = .{ .child = al };
    const pa = poison.allocator();
    const alpha = scalar(T, -0.75, 0.25);
    const beta = scalar(T, 0.5, 0);
    for (0..rounds + 1) |round| {
        const w: u32 = if (round == 0) 1 else @intCast(2 + round % 7);
        const out = if (round == 0) want else got;
        @memcpy(out, c0);
        const rio = if (round == 0) Io.failing else io;
        if (r) |rr| {
            try routine(T, &vt, pa, rio, rr, opt, m, n, k, alpha, beta, a.ptr, b.ptr, ld, out.ptr, w);
        } else {
            try driver.gemm(pa, rio, T, &vt, opt[0], opt[1], .n, m, n, k, alpha, a.ptr, ld, b.ptr, ld, beta, out.ptr, ld, w);
        }
        if (round > 0 and !sameBits(T, want, got)) {
            std.debug.print("{s} {?t} round {d} w={d}\n", .{ @typeName(T), r, round, w });
            return error.TestUnexpectedResult;
        }
    }
}

test "race stress: threaded level 3 keeps one worker's bits and its guard bands, round after round" {
    if (comptime !threading.hosted) return error.SkipZigTest;
    var threaded: Io.Threaded = .init(std.testing.allocator, .{});
    defer threaded.deinit();
    try stresses(threaded.io());
    if (comptime threading.evented_workers) {
        var ev: threading.Evented = undefined;
        try zigstd.testExpectTrue(threading.initEvented(&ev, 3));
        defer ev.deinit();
        try stresses(ev.io());
    }
}

/// The race stress test's routines on `io`.
fn stresses(io: Io) !void {
    const rounds = threading.stress_rounds;
    // shared packed `B` behind its barriers: `ic` and `jc` past 1, several slices
    {
        const blk = Blocking.of(f64, &shrunkCore(f64));
        try stress(f64, io, null, .{ false, true, false }, 9 * blk.um + 5, 2 * blk.r + 3, 2 * blk.q + 3, rounds);
    }
    // a triangle's balanced groups, both passes
    {
        const blk = Blocking.of(types.Complex32, &shrunkCore(types.Complex32));
        try stress(types.Complex32, io, .syr2k, .{ false, true, false }, 2 * blk.r + 3, 2 * blk.r + 3, 2 * blk.q + 3, rounds / 2);
    }
    // independent panels through `forEach`
    {
        const blk = Blocking.of(f32, &shrunkCore(f32));
        try stress(f32, io, .trsm, .{ true, false, false }, 2 * blk.r + 3, 2 * blk.r + 3, 1, rounds / 2);
    }
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
        try cancelGemm(threaded.io(), &vt, a, b, c, n, k);
    }
    if (comptime threading.evented_workers) {
        var ev: threading.Evented = undefined;
        try zigstd.testExpectTrue(threading.initEvented(&ev, 3));
        defer ev.deinit();
        try cancelGemm(ev.io(), &vt, a, b, c, n, k);
    }
}

/// Cancels a held `gemm` on `io` once a kernel call has started.
fn cancelGemm(io: Io, vt: *const Vtable, a: []f64, b: []f64, c: []f64, n: usize, k: usize) !void {
    Hold.calls.store(0, .release);
    Hold.held.store(true, .release);
    Hold.first.store(0, .release);
    Hold.others.store(0, .release);
    const releaser = try std.Thread.spawn(.{}, Hold.release, .{io});
    var future = io.concurrent(canceledGemm, .{ io, vt, a, b, c, n, k }) catch return error.SkipZigTest;
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
