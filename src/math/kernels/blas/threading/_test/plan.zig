//! `threading.zig`'s tests of the rules, the chunks, the worker count and
//! the plan.

const std = @import("std");
const sysinfo = @import("builtin");
const blas = @import("../../__root__.zig");
const Config = @import("../../Config.zig");
const threading = @import("../../threading.zig");
const zigstd = @import("../../../../../zigstd.zig");
const coretypes = @import("../../cpu/coretypes.zig");

const chunks_max = threading.chunks_max;
const Family = threading.Family;
const rule = threading.rule;
const product = threading.product;
const cpuCount = threading.cpuCount;
const workers = threading.workers;
const Shape = threading.Shape;
const Chunks = threading.Chunks;
const serialFloor = threading.serialFloor;
const serial = threading.serial;
const plan = threading.plan;
const planWith = threading.planWith;

const X86 = coretypes.X86CoreType;
const Arm64 = coretypes.Arm64CoreType;

test "chunks cover the range once, in order, on quantum bounds" {
    const lens = [_]usize{ 0, 1, 2, 3, 10, 63, 64, 65, 100, 1000, 12345, std.math.maxInt(usize) };
    const wants = [_]usize{ 0, 1, 2, 3, 7, 64, 65, 1000 };
    const quanta = [_]usize{ 0, 1, 4, 16 };
    for (lens) |len| for (wants) |want| for (quanta) |q| inline for (.{ Shape.even, Shape.rising, Shape.falling }) |shape| {
        const c: Chunks = .init(len, want, shape, q);
        try zigstd.testExpectTrue(c.count >= 1 and c.count <= @max(@min(want, chunks_max), 1));
        try zigstd.testExpectEqual(@as(usize, 0), c.bounds[0]);
        try zigstd.testExpectEqual(len, c.bounds[c.count]);
        for (0..c.count) |i| {
            const ch = c.get(i);
            try zigstd.testExpectEqual(i, ch.index);
            if (len > 0) try zigstd.testExpectTrue(ch.start < ch.end);
            if (i > 0) try zigstd.testExpectEqual(@as(usize, 0), ch.start % @max(q, 1));
        }
    };
}

test "even chunks differ by at most a quantum; rising and falling ones balance a triangle" {
    const even: Chunks = .init(1000, 7, .even, 1);
    try zigstd.testExpectEqual(@as(usize, 7), even.count);
    for (0..even.count) |i| {
        const w = even.get(i).end - even.get(i).start;
        try zigstd.testExpectTrue(w == 142 or w == 143);
    }
    const len = 10000;
    inline for (.{ Shape.rising, Shape.falling }) |shape| {
        const c: Chunks = .init(len, 8, shape, 1);
        try zigstd.testExpectEqual(@as(usize, 8), c.count);
        for (0..c.count) |i| {
            var cost: u64 = 0;
            for (c.get(i).start..c.get(i).end) |u| cost += if (shape == .rising) u + 1 else len - u;
            const share = @as(u64, len) * (len + 1) / 2 / 8;
            // a bound is off by at most a row or two of `len`
            try zigstd.testExpectTrue(cost + 3 * len >= share and cost <= share + 3 * len);
        }
    }
}

test "rules: thresholds at their edges, tiers and the chunk cap" {
    const axpy = rule(.axpy, f64, X86.haswell);
    try zigstd.testExpectEqual(@as(usize, 1), axpy.chunkCount(10000));
    try zigstd.testExpectEqual(@as(u32, 1), axpy.sizeMax(10000));
    try zigstd.testExpectEqual(@as(usize, 2), axpy.chunkCount(10001));
    try zigstd.testExpectEqual(@as(u32, std.math.maxInt(u32)), axpy.sizeMax(10001));
    try zigstd.testExpectEqual(@as(usize, chunks_max), axpy.chunkCount(1 << 40));

    const gemm = rule(.gemm, f64, X86.haswell);
    try zigstd.testExpectEqual(@as(u32, 1), gemm.sizeMax(262144));
    try zigstd.testExpectEqual(@as(u32, 5), gemm.sizeMax(262144 * 5 + 1));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.gemm, std.math.Complex(f32), X86.haswell).sizeMax(32768));
    try zigstd.testExpectEqual(@as(u32, 2), rule(.gemm, std.math.Complex(f32), X86.haswell).sizeMax(65536));

    const ger = rule(.ger, f32, X86.haswell);
    try zigstd.testExpectEqual(@as(u32, 1), ger.sizeMax(8192));
    try zigstd.testExpectEqual(@as(u32, 3), ger.sizeMax(8193));
    try zigstd.testExpectEqual(@as(u32, 3), ger.sizeMax(81920));
    try zigstd.testExpectEqual(@as(u32, std.math.maxInt(u32)), ger.sizeMax(81921));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.ger, std.math.Complex(f64), X86.haswell).sizeMax(9216));
    try zigstd.testExpectEqual(@as(u32, std.math.maxInt(u32)), rule(.ger, std.math.Complex(f64), X86.haswell).sizeMax(9217));

    // `n < 8 MiB / sizeof` serial; complex counts its elements
    try zigstd.testExpectEqual(@as(u32, 1), rule(.swap, f64, X86.haswell).sizeMax(1048575));
    try zigstd.testExpectEqual(@as(usize, 2), rule(.swap, f64, X86.haswell).chunkCount(1048576));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.swap, std.math.Complex(f64), X86.haswell).sizeMax(524287));
    try zigstd.testExpectEqual(@as(usize, 2), rule(.swap, std.math.Complex(f64), X86.haswell).chunkCount(524288));

    // OpenBLAS's level 2, on a core without rstd's own (haswell's, below)
    const ztrmv = rule(.trmv, std.math.Complex(f64), X86.zen);
    try zigstd.testExpectEqual(@as(u32, 1), ztrmv.sizeMax(96));
    try zigstd.testExpectEqual(@as(u32, 2), ztrmv.sizeMax(97));
    try zigstd.testExpectEqual(@as(u32, 2), ztrmv.sizeMax(127));
    try zigstd.testExpectEqual(@as(u32, std.math.maxInt(u32)), ztrmv.sizeMax(128));
    try zigstd.testExpectEqual(@as(u32, 2), rule(.trmv, f64, X86.zen).sizeMax(499));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.symv, f64, X86.haswell).sizeMax(199));
    try zigstd.testExpectEqual(@as(usize, 2), rule(.symv, f64, X86.haswell).chunkCount(200));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.hemv, std.math.Complex(f64), X86.haswell).sizeMax(361));
    try zigstd.testExpectEqual(@as(usize, 8), rule(.hemv, std.math.Complex(f64), X86.haswell).chunkCount(362));
    try zigstd.testExpectEqual(@as(usize, 21), rule(.hemv, std.math.Complex(f64), X86.haswell).chunkCount(1024));
    try zigstd.testExpectEqual(@as(u32, 2), rule(.tpmv, std.math.Complex(f32), X86.zen).sizeMax(499));

    try zigstd.testExpectEqual(@as(u32, 1), rule(.asum, f32, X86.haswell).sizeMax(100000));
    try zigstd.testExpectEqual(@as(u32, 3), rule(.asum, f32, X86.haswell).sizeMax(300000));
    try zigstd.testExpectEqual(@as(u32, 3), rule(.asum, std.math.Complex(f32), X86.haswell).sizeMax(30000));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.asum, f64, Arm64.vortexm4).sizeMax(10000));
    try zigstd.testExpectEqual(@as(u32, std.math.maxInt(u32)), rule(.asum, f64, Arm64.vortexm4).sizeMax(10001));
    try zigstd.testExpectEqual(@as(u32, 3), rule(.asum, f32, Arm64.cortexa53).sizeMax(300000));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.trsm, f64, X86.haswell).sizeMax(1023));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.syrk, f64, X86.haswell).sizeMax(product(&.{ 101, 100, 43 })));
    try zigstd.testExpectEqual(@as(u32, std.math.maxInt(u32)), rule(.syrk, f64, X86.haswell).sizeMax(product(&.{ 101, 100, 44 })));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.syr2k, f64, X86.haswell).sizeMax(999));
    try zigstd.testExpectEqual(@as(u32, std.math.maxInt(u32)), rule(.syr2k, std.math.Complex(f32), X86.haswell).sizeMax(1000));
}

test "rules: haswell's own level-2 thresholds (issue #845), the other cores keeping OpenBLAS's" {
    const C = std.math.Complex;
    inline for (.{ .syr, .tpmv }) |f| {
        try zigstd.testExpectEqual(@as(u32, 1), rule(f, f64, X86.haswell).sizeMax(511));
        try zigstd.testExpectEqual(@as(usize, 2), rule(f, f64, X86.haswell).chunkCount(512));
    }
    try zigstd.testExpectEqual(@as(u32, 1), rule(.her, C(f32), X86.haswell).sizeMax(255));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.trmv, f32, X86.haswell).sizeMax(511));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.trmv, C(f32), X86.haswell).sizeMax(383));
    try zigstd.testExpectEqual(@as(usize, 2), rule(.trmv, C(f64), X86.haswell).chunkCount(256));
    // `tbmv` (issue #854), chunked as `tpmv`
    try zigstd.testExpectEqual(@as(u32, 1), rule(.tbmv, f32, X86.haswell).sizeMax(1535));
    try zigstd.testExpectEqual(@as(usize, 4), rule(.tbmv, f32, X86.haswell).chunkCount(1536));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.tbmv, C(f32), X86.haswell).sizeMax(1023));
    try zigstd.testExpectEqual(@as(usize, 3), rule(.tbmv, f64, X86.haswell).chunkCount(1024));
    try zigstd.testExpectEqual(@as(usize, 2), rule(.tbmv, f64, X86.zen).chunkCount(50));
    // the complex symmetric `syr` and other cores keep OpenBLAS's
    try zigstd.testExpectEqual(@as(usize, 2), rule(.syr, C(f64), X86.haswell).chunkCount(50));
    try zigstd.testExpectEqual(@as(usize, 2), rule(.syr, f64, X86.zen).chunkCount(50));
    try zigstd.testExpectEqual(@as(usize, 2), rule(.trmv, f64, X86.skylakex).chunkCount(50));
}

test "rules: vortexm4's own level-2 thresholds (issue #853), the other cores keeping OpenBLAS's" {
    const C = std.math.Complex;
    const A = Arm64.vortexm4;
    try zigstd.testExpectEqual(@as(u32, 1), rule(.syr, f32, A).sizeMax(767));
    try zigstd.testExpectEqual(@as(usize, 2), rule(.syr, f64, A).chunkCount(768));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.her, C(f32), A).sizeMax(383));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.her, C(f64), A).sizeMax(255));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.trmv, f32, A).sizeMax(1023));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.trmv, f64, A).sizeMax(511));
    try zigstd.testExpectEqual(@as(usize, 2), rule(.trmv, C(f64), A).chunkCount(256));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.tpmv, f64, A).sizeMax(767));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.tpmv, C(f32), A).sizeMax(383));
    try zigstd.testExpectEqual(@as(u32, 2), rule(.tpmv, C(f32), A).sizeMax(499));
    try zigstd.testExpectEqual(@as(u32, std.math.maxInt(u32)), rule(.tpmv, C(f64), A).sizeMax(500));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.tbmv, f32, A).sizeMax(2047));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.tbmv, C(f64), A).sizeMax(1535));
    // `symv` (the fused tridiagonalization's too, #885), the complex symmetric
    // `syr` and other arm64 cores keep OpenBLAS's
    try zigstd.testExpectEqual(@as(usize, 2), rule(.symv, f64, A).chunkCount(200));
    try zigstd.testExpectEqual(@as(u32, 2), rule(.symv, C(f32), A).sizeMax(50));
    try zigstd.testExpectEqual(@as(usize, 2), rule(.syr, C(f64), A).chunkCount(50));
    try zigstd.testExpectEqual(@as(usize, 2), rule(.syr, f64, Arm64.vortex).chunkCount(50));
}

test "rules: vortexm4's own axpy and dot split (issue #728), the other cores keeping OpenBLAS's" {
    const C = std.math.Complex;
    inline for (.{ .axpy, .dot }) |f| inline for (.{ f32, f64, C(f32), C(f64) }) |T| {
        const r = rule(f, T, Arm64.vortexm4);
        try zigstd.testExpectEqual(@as(u32, 1), r.sizeMax(10000));
        try zigstd.testExpectEqual(@as(u32, 2), r.sizeMax(20000));
        try zigstd.testExpectEqual(@as(u32, 4), r.sizeMax(32768));
        try zigstd.testExpectEqual(@as(u32, 6), r.sizeMax(65536));
        try zigstd.testExpectEqual(@as(u32, std.math.maxInt(u32)), r.sizeMax(70001));
        try zigstd.testExpectEqual(@as(usize, 7), r.chunkCount(32768));
        try zigstd.testExpectEqual(@as(usize, 14), r.chunkCount(65536));
        try zigstd.testExpectEqual(@as(usize, 2), rule(f, T, Arm64.vortex).chunkCount(10001));
        try zigstd.testExpectEqual(@as(u32, std.math.maxInt(u32)), rule(f, T, Arm64.vortex).sizeMax(65536));
    };
}

test "rules: OpenBLAS's per-core overrides apply to their core and precision only" {
    try zigstd.testExpectEqual(@as(u32, 1), rule(.dot, f64, Arm64.neoversev1).sizeMax(64500));
    try zigstd.testExpectEqual(@as(u32, 2), rule(.dot, f64, Arm64.neoversev1).sizeMax(64501));
    try zigstd.testExpectEqual(@as(u32, 56), rule(.dot, f64, Arm64.neoversev1).sizeMax(1010000));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.dot, f32, Arm64.neoversev1).sizeMax(110000));
    try zigstd.testExpectEqual(@as(u32, 8), rule(.dot, f32, Arm64.a64fx).sizeMax(1020000));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.dot, std.math.Complex(f64), Arm64.neoversev1).sizeMax(10000));
    try zigstd.testExpectEqual(@as(u32, std.math.maxInt(u32)), rule(.dot, std.math.Complex(f64), Arm64.neoversev1).sizeMax(10001));

    try zigstd.testExpectEqual(@as(u32, 1), rule(.gemm, f32, Arm64.neoversev1).sizeMax(262143));
    try zigstd.testExpectEqual(@as(u32, 6), rule(.gemm, f32, Arm64.neoversev1).sizeMax(262144));
    try zigstd.testExpectEqual(@as(u32, 64), rule(.gemm, f32, Arm64.neoversev2).sizeMax(1073741823));
    // double gemm keeps the generic rule
    try zigstd.testExpectEqual(@as(u32, 1), rule(.gemm, f64, Arm64.neoversev1).sizeMax(262144));

    try zigstd.testExpectEqual(@as(u32, 1), rule(.gemv, f64, Arm64.neoversev1).sizeMax(8099));
    try zigstd.testExpectEqual(@as(u32, 2), rule(.gemv, f64, Arm64.neoversev1).sizeMax(8100));
    try zigstd.testExpectEqual(@as(u32, 4), rule(.gemv, f32, Arm64.neoversev1).sizeMax(25600));
    try zigstd.testExpectEqual(@as(u32, 8), rule(.gemv, f64, Arm64.neoversev2).sizeMax(24964));
    try zigstd.testExpectEqual(@as(u32, 1), rule(.gemv, f64, X86.haswell).sizeMax(460799));
}

test "serial is exact: its floor is one chunk and one worker on every core" {
    inline for (comptime std.meta.tags(Family)) |f| {
        inline for (.{ f32, f64, std.math.Complex(f32), std.math.Complex(f64) }) |T| {
            const floor = comptime serialFloor(f, T);
            try zigstd.testExpectTrue(serial(f, T, floor));
            try zigstd.testExpectTrue(!serial(f, T, floor + 1));
            inline for (coretypes.all) |C| {
                for (std.meta.tags(C)) |c| {
                    const r = rule(f, T, c);
                    try zigstd.testExpectEqual(@as(usize, 1), r.chunkCount(floor));
                    try zigstd.testExpectEqual(@as(u32, 1), r.sizeMax(floor));
                }
            }
        }
    }
    // the per-core overrides lower the floor below a default
    try zigstd.testExpectEqual(@as(u64, 8099), comptime serialFloor(.gemv, f64));
}

test "the worker count is the least of the CPUs, the ceiling and the size" {
    var cfg: Config = .{ .cpu_count = 8 };
    try zigstd.testExpectTrue(cfg.threadCeiling() == null);
    const expect = struct {
        fn f(c: *const Config, size_max: u32, want: u32) !void {
            try zigstd.testExpectEqual(if (sysinfo.single_threaded) 1 else want, workers(c, size_max));
        }
    }.f;
    try expect(&cfg, std.math.maxInt(u32), 8);
    try expect(&cfg, 5, 5);
    try expect(&cfg, 0, 1);
    cfg.goto_num_threads = 6;
    try expect(&cfg, 7, 6);
    cfg.num_threads = 3;
    try expect(&cfg, 7, 3);
    try expect(&cfg, 2, 2);
    cfg.cpu_count = 1;
    try expect(&cfg, 7, 1);
    cfg = .{ .cpu_count = 2, .num_threads = 16 };
    try expect(&cfg, std.math.maxInt(u32), 2);
    try zigstd.testExpectTrue(cpuCount() >= 1);
}

test "chunks follow the size and the core, never the CPU count or the ceiling" {
    const cfgs = [_]Config{ .{}, .{ .cpu_count = 4 }, .{ .cpu_count = 64 }, .{ .cpu_count = 64, .num_threads = 2 }, .{ .cpu_count = 3, .goto_num_threads = 5 } };
    const r = rule(.dot, f64, X86.haswell);
    for ([_]usize{ 1, 10000, 10001, 123457, 1 << 30 }) |n| {
        const base = planWith(&cfgs[0], r, n, n, .even, 8);
        for (cfgs[1..]) |*cfg| {
            const p = planWith(cfg, r, n, n, .even, 8);
            try zigstd.testExpectEqual(base.chunks.count, p.chunks.count);
            try zigstd.testExpectEqualSlices(usize, base.chunks.bounds[0 .. base.chunks.count + 1], p.chunks.bounds[0 .. p.chunks.count + 1]);
            try zigstd.testExpectTrue(p.workers >= 1 and p.workers <= p.chunks.count);
        }
    }
}

test "single-threaded builds resolve to one worker over the same chunks" {
    if (!sysinfo.single_threaded) return error.SkipZigTest;
    _ = blas.vtable();
    try zigstd.testExpectEqual(@as(u32, 1), blas.config.cpu_count);
    const p = plan(.dot, f64, 1 << 20, 1 << 20, .even, 1);
    try zigstd.testExpectEqual(@as(u32, 1), p.workers);
    try zigstd.testExpectEqual(@as(usize, chunks_max), p.chunks.count);
}

test "product saturates" {
    try zigstd.testExpectEqual(@as(u64, 24), product(&.{ 2, 3, 4 }));
    try zigstd.testExpectEqual(@as(u64, std.math.maxInt(u64)), product(&.{ 1 << 40, 1 << 40 }));
}
