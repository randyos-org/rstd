//! `lib3`'s `gemm` race.

const std = @import("std");
const rstd = @import("rstd");
const bb = @import("../blas_bench.zig");
const kernels = @import("../kernels.zig");
const ext = @import("../ext.zig");
const lb = @import("../libs.zig");
const lv = @import("via.zig");
const tb = @import("table.zig");
const l2r = @import("level2.zig");

const Io = std.Io;
const blas = rstd.math.kernels.blas;
const clock = bb.clock;
const l3_timing = bb.l3_timing;
const ob_asm = bb.ob_asm;
const Prec = bb.Prec;
const precs = bb.precs;
const asmGemm = kernels.asmGemm;
const Blocking = kernels.Blocking;
const checkGemm = kernels.checkGemm;
const gemm_sizes = kernels.gemm_sizes;
const GemmCall = kernels.GemmCall;
const ext_l3 = ext.ext_l3;
const extTable = ext.extTable;
const rstd_alloc = lb.rstd_alloc;
const GemmFn = lv.GemmFn;
const libGemm = lv.libGemm;
const Entrants = tb.Entrants;
const diff_note = tb.diff_note;
const lib3_cols = tb.lib3_cols;
const printHeader = tb.printHeader;
const rstdFirst = tb.rstdFirst;
const rstd_col = tb.rstd_col;
const sectionMark = tb.sectionMark;
const sumBound = tb.sumBound;
const diffScaled = l2r.diffScaled;
const magAt = l2r.magAt;

/// `GemmCall`, which is real only and needs the `.S`.
pub fn Loops(comptime p: Prec) type {
    return if (p.complex or !ob_asm) void else GemmCall(p);
}

/// One gemm contestant: a library's, or OpenBLAS's `.S` through the bench's
/// copy of `level3.c`'s loops (`GemmCall`).
pub fn LibGemmCall(comptime p: Prec) type {
    return struct {
        f: ?GemmFn(p),
        loops: ?Loops(p),
        n: usize,
        a: [*]p.F,
        b: [*]p.F,
        c: [*]p.F,
        pub fn call(s: *const @This()) void {
            if (s.f) |f| return f(s.n, s.a, s.b, s.c);
            if (comptime Loops(p) != void) s.loops.?.call();
        }
        pub fn run(o: *anyopaque) void {
            const s: *@This() = @ptrCast(@alignCast(o));
            s.call();
        }
    };
}

/// The size of the terms of each float of `c0 + A*B`, square `n`:
/// `|c0| + |A| |B|`, a complex's size `|re| + |im|`, by rstd's dgemm.
/// `abs_a` and `abs_b` are scratch.
pub fn gemmScale(comptime p: Prec, n: usize, a: []const p.F, b: []const p.F, c0: []const p.F, abs_a: []f64, abs_b: []f64, out: []f64) !void {
    const cs = p.cs();
    const sq = n * n;
    for (0..sq) |i| {
        abs_a[i] = magAt(p, a, i * cs);
        abs_b[i] = magAt(p, b, i * cs);
        out[i] = magAt(p, c0, i * cs);
    }
    try blas.gemm(rstd_alloc, lb.rstd_io, f64, .n, .n, n, n, n, 1, abs_a, n, abs_b, n, 1, out, n);
    // a complex's two floats share its scale, spread from the back
    if (p.complex) {
        var i = sq;
        while (i > 0) {
            i -= 1;
            out[2 * i + 1] = out[i];
            out[2 * i] = out[i];
        }
    }
}

pub fn libLevel3(io: Io, gpa: std.mem.Allocator) !void {
    sectionMark("lib3");
    try clock(io, "library level 3 start");
    std.debug.print(
        \\
        \\gemm NN, square, `C += A*B`: GFLOP/s, best of {d} batches of at least {d} ms; ratios are over OpenBLAS{s}. {s} {s}
        \\
        \\
    , .{ l3_timing.batches, @divTrunc(l3_timing.min_ns, std.time.ns_per_ms), if (ob_asm) ". `ob .S` is OpenBLAS's vortex `.S` through the bench's loops" else "", lb.io_note, diff_note });
    printHeader("p | n", "- | -:", &lib3_cols);
    inline for (precs) |p| try raceLibGemm(io, gpa, p);
    try extTable(io, gpa, .l3, ext_l3);
}

pub fn raceLibGemm(io: Io, gpa: std.mem.Allocator, comptime p: Prec) !void {
    const F = p.F;
    const cs = p.cs();
    const max_n = gemm_sizes[gemm_sizes.len - 1];
    const len = max_n * max_n * cs;
    const a = try gpa.alignedAlloc(F, .@"64", len);
    defer gpa.free(a);
    const b = try gpa.alignedAlloc(F, .@"64", len);
    defer gpa.free(b);
    const c = try gpa.alignedAlloc(F, .@"64", len);
    defer gpa.free(c);
    const c0 = try gpa.alloc(F, len);
    defer gpa.free(c0);
    // rstd's output, and the size of each float's terms
    const kc = try gpa.alloc(F, len);
    defer gpa.free(kc);
    const scale = try gpa.alloc(f64, len);
    defer gpa.free(scale);
    const abs_a = try gpa.alloc(f64, max_n * max_n);
    defer gpa.free(abs_a);
    const abs_b = try gpa.alloc(f64, max_n * max_n);
    defer gpa.free(abs_b);
    var rng = std.Random.DefaultPrng.init(4);
    const rand = rng.random();
    for (a, b, c0) |*x, *y, *z| {
        x.* = @floatCast(rand.float(f64) * 2 - 1);
        y.* = @floatCast(rand.float(f64) * 2 - 1);
        z.* = @floatCast(rand.float(f64) * 2 - 1);
    }

    // the `.S` through the bench's loops, real only
    var sa: []align(64) F = &.{};
    var sb: []align(64) F = &.{};
    var blk: Blocking = undefined;
    if (comptime Loops(p) != void) {
        blk = .of(p, asmGemm(p));
        sa = try gpa.alignedAlloc(F, .@"64", blk.p * (blk.q + blk.um));
        sb = try gpa.alignedAlloc(F, .@"64", (blk.q + blk.um) * blk.r);
    }
    defer if (comptime Loops(p) != void) {
        gpa.free(sa);
        gpa.free(sb);
    };

    const fns = [lib3_cols.len]?GemmFn(p){ libGemm(.openblas, p), libGemm(.rstd, p), libGemm(.blis, p), libGemm(.accelerate, p), null };
    const LC = LibGemmCall(p);
    for (gemm_sizes) |n| {
        if (n > lb.n_max) continue;
        const m = n * n * cs;
        // rstd's dgemm; nothing of rstd's runs in `alone`, where nothing is compared
        if (!bb.alone) try gemmScale(p, n, a[0..m], b[0..m], c0[0..m], abs_a, abs_b, scale[0..m]);
        var es: Entrants(LC, lib3_cols.len) = .{};
        for (0..lib3_cols.len) |o| {
            const i = rstdFirst(o);
            if (!lib3_cols[i].shown()) continue;
            const loops: ?Loops(p) = if (comptime Loops(p) != void) (if (i == 4) .{ .im = asmGemm(p), .blk = blk, .n = n, .a = a.ptr, .b = b.ptr, .c = c.ptr, .sa = sa.ptr, .sb = sb.ptr } else null) else null;
            if (fns[i] == null and loops == null) continue;
            const ctx = es.add(i, .{ .f = fns[i], .loops = loops, .n = n, .a = a.ptr, .b = b.ptr, .c = c.ptr });
            @memcpy(c[0..m], c0[0..m]);
            ctx.call();
            if (lb.lib_err) |err| return err;
            const e = checkGemm(p, n, n, n, a[0..m], b[0..m], c[0..m], c0[0..m], 256);
            if (e > 1) {
                std.debug.print("WRONG: {s}gemm {s} n={d}: {d:.2} of tolerance\n", .{ p.name, lib3_cols[i].label, n, e });
                return error.WrongGemm;
            }
            if (i == rstd_col) {
                @memcpy(kc[0..m], c[0..m]);
                es.mine = .{ .err = e * @as(f64, @floatFromInt(n + 8)), .bits = std.hash.Wyhash.hash(0, std.mem.sliceAsBytes(c[0..m])) };
            } else if (lib3_cols[i].lib != null) es.diffs[i] = .{ .eps = diffScaled(F, kc, c, scale[0..m]), .bound = sumBound(n) };
        }
        const nf: f64 = @floatFromInt(n);
        const rates = try es.race(io, l3_timing, nf * nf * nf * @as(f64, if (p.complex) 8 else 2));
        std.debug.print("| {s} | {d} | ", .{ p.name, n });
        es.print(&lib3_cols, &rates);
    }
}
