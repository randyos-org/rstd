//! The `sparse` section (T151; only when named, any host): CSC `A x` and
//! `Aᵀ y` (`rstd.math.kernels.sparse`) against dense `gemv` on the same matrix,
//! `f64`, one worker, over fill (1% to 60%), size (100 to 3000, 10 000 with
//! `--max-n=10000`: its dense matrix is 800 MB) and pattern (random,
//! banded): the crossover for a dense-or-sparse switch; and `Aᵀ y` in four
//! sums, with scalar loads and (x86_64 with AVX2) a gather (`sparse_gather.zig`).

const std = @import("std");
const Io = std.Io;
const rstd = @import("rstd");
const bench = @import("bench.zig");
const bb = @import("blas_bench.zig");
const libs = @import("libs.zig");
const sg = @import("sparse_gather.zig");

const blas = rstd.linalg.blas;
const CompressedColumns = rstd.math.kernels.sparse.CompressedColumns(f64);
const fail = libs.fail;

const timing: bench.Timing = .{ .min_ns = 20 * std.time.ns_per_ms, .batches = 7 };
const sizes = [_]usize{ 100, 300, 1000, 3000, 10_000 };
const fills = [_]f64{ 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6 };
const Pattern = enum { random, banded };

pub fn run(io: Io, gpa: std.mem.Allocator) !void {
    try bb.clock(io, "sparse start");
    // 10 000 only when asked: its dense matrix is 800 MB
    const cap = if (libs.n_max == std.math.maxInt(usize)) 3000 else libs.n_max;
    std.debug.print(
        \\
        \\µs per call, best of {d} batches of at least {d} ms; one worker; `f64`; ratios are CSC over dense, the variants over CSC Aᵀ y (gather: x86_64 with AVX2)
        \\
        \\| pattern | n | fill | nnz | gemv N | CSC A x | ratio | gemv T | CSC Aᵀ y | ratio | 4 sums | ratio | gather | ratio |
        \\| ------- | -: | ---: | --: | -----: | ------: | ----: | -----: | -------: | ----: | -----: | ----: | -----: | ----: |
        \\
    , .{ timing.batches, @divTrunc(timing.min_ns, std.time.ns_per_ms) });
    for (std.meta.tags(Pattern)) |pat| {
        for (sizes) |n| {
            if (n > cap) continue;
            for (fills) |f| try race(io, gpa, pat, n, f);
        }
    }
}

fn race(io: Io, gpa: std.mem.Allocator, pat: Pattern, n: usize, fill: f64) !void {
    const a = try gpa.alloc(f64, n * n);
    defer gpa.free(a);
    var nnz: usize = 0;
    var rng = std.Random.DefaultPrng.init(n);
    const rand = rng.random();
    // banded: half-width `w`, so the fill is about `(2 w + 1) / n`
    const w: usize = @intFromFloat(@max(0, (fill * @as(f64, @floatFromInt(n)) - 1) / 2));
    for (0..n) |j| for (0..n) |i| {
        const on = switch (pat) {
            .random => rand.float(f64) < fill,
            .banded => (if (i > j) i - j else j - i) <= w,
        };
        a[i + j * n] = if (on) rand.float(f64) * 2 - 1 else 0;
        if (on) nnz += 1;
    };
    var s = try CompressedColumns.alloc(gpa, @intCast(n), @intCast(n), nnz);
    defer s.deinit(gpa);
    var k: u32 = 0;
    for (0..n) |j| {
        s.col_start[j] = k;
        for (0..n) |i| if (a[i + j * n] != 0) {
            s.row[k] = @intCast(i);
            s.value[k] = a[i + j * n];
            k += 1;
        };
    }
    s.col_start[n] = k;

    const xy = try gpa.alloc(f64, 2 * n);
    defer gpa.free(xy);
    for (xy, 0..) |*v, i| v.* = @floatFromInt(@as(i64, @intCast(i % 7)) - 3);
    // the gather's bits are the scalar four sums': same order, other loads
    if (sg.has_gather) {
        const ys = try gpa.alloc(f64, 2 * n);
        defer gpa.free(ys);
        sg.dots4(s, xy[n..], ys[0..n]);
        sg.dotsGather(s, xy[n..], ys[n..]);
        if (!std.mem.eql(f64, ys[0..n], ys[n..])) return error.GatherMismatch;
    }
    var ctx: [n_calls]Call = undefined;
    for (&ctx, 0..) |*c, i| c.* = .{ .gpa = gpa, .a = a, .s = s, .n = n, .x = xy[0..n], .y = xy[n..], .which = @fromBackingInt(@intCast(i)) };
    var calls: [n_calls]bench.Call = undefined;
    for (&calls, &ctx) |*c, *x| c.* = .of(x);
    var secs: [n_calls]f64 = undefined;
    try bench.race(io, timing, &calls, &secs);
    if (libs.lib_err) |err| return err;
    var us: [6]f64 = @splat(std.math.nan(f64));
    for (us[0..n_calls], secs) |*u, t| u.* = t * 1e6;
    const real_fill = @as(f64, @floatFromInt(nnz)) / @as(f64, @floatFromInt(n * n));
    std.debug.print("| {t} | {d} | {d:.1}% | {d} | {d:.2} | {d:.2} | {d:.2} | {d:.2} | {d:.2} | {d:.2} | {d:.2} | {d:.2} | {d:.2} | {d:.2} |\n", .{ pat, n, 100 * real_fill, nnz, us[0], us[1], us[1] / us[0], us[2], us[3], us[3] / us[2], us[4], us[4] / us[3], us[5], us[5] / us[3] });
}

/// The contestants, the gather last, raced only where built.
const Which = enum { dense_n, csc_n, dense_t, csc_t, csc_t4, csc_tg };
const n_calls = if (sg.has_gather) 6 else 5;

const Call = struct {
    gpa: std.mem.Allocator,
    a: []const f64,
    s: CompressedColumns,
    n: usize,
    x: []f64,
    y: []f64,
    which: Which,

    pub fn run(o: *anyopaque) void {
        const c: *Call = @ptrCast(@alignCast(o));
        const one = Io.failing;
        const r = switch (c.which) {
            .dense_n => blas.gemv(c.gpa, one, f64, .n, 1, c.a, c.n, .of(c.x), 0, .of(c.y)),
            .csc_n => c.s.mulVec(c.gpa, one, .n, 1, c.x, 0, c.y),
            .dense_t => blas.gemv(c.gpa, one, f64, .t, 1, c.a, c.n, .of(c.y), 0, .of(c.x)),
            .csc_t => c.s.mulVec(c.gpa, one, .t, 1, c.y, 0, c.x),
            .csc_t4 => sg.dots4(c.s, c.y, c.x),
            .csc_tg => sg.dotsGather(c.s, c.y, c.x),
        };
        r catch |e| fail(e);
    }
};
