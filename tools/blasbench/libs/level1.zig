//! `lib1`: each library's level 1 against rstd's.

const std = @import("std");
const rstd = @import("rstd");
const bb = @import("../blas_bench.zig");
const kernels = @import("../kernels.zig");
const lb = @import("../libs.zig");
const lv = @import("via.zig");
const tb = @import("table.zig");

const Io = std.Io;
const arm = bb.arm;
const clock = bb.clock;
const Fn = bb.Fn;
const l1_timing = bb.l1_timing;
const Prec = bb.Prec;
const precs = bb.precs;
const asmFor = kernels.asmFor;
const Bufs = kernels.Bufs;
const candidate = kernels.candidate;
const checkL1 = kernels.checkL1;
const l1_cases = kernels.l1_cases;
const l1Bytes = kernels.l1Bytes;
const L1 = kernels.L1;
const L1Call = kernels.L1Call;
const l1Elem = kernels.l1Elem;
const l1Result = kernels.l1Result;
const l1Scalar = kernels.l1Scalar;
const l1Slot = kernels.l1Slot;
const reduces = kernels.reduces;
const libL1 = lv.libL1;
const Diff = tb.Diff;
const Entrants = tb.Entrants;
const diff_note = tb.diff_note;
const lib_cols = tb.lib_cols;
const printHeader = tb.printHeader;
const printLead = tb.printLead;
const rstdFirst = tb.rstdFirst;
const rstdRuns = tb.rstdRuns;
const rstd_col = tb.rstd_col;
const sectionMark = tb.sectionMark;
const sumBound = tb.sumBound;

pub fn libLevel1(io: Io, gpa: std.mem.Allocator) !void {
    sectionMark("lib1");
    try clock(io, "library level 1 start");
    std.debug.print(
        \\
        \\GB/s, best of {d} batches of at least {d} ms; ratios are over OpenBLAS. `rstd` is the public routine{s}. {s} {s}
        \\
        \\
    , .{ l1_timing.batches, @divTrunc(l1_timing.min_ns, std.time.ns_per_ms), if (arm) ", `asm` the converted kernel alone (`from`)" else "", lb.io_note, diff_note });
    if (arm)
        printHeader("routine | p | n | inc | rstd runs | from", "------- | - | -: | --: | --------- | ----", &lib_cols)
    else
        printHeader("routine | p | n | inc | rstd runs", "------- | - | -: | --: | ---------", &lib_cols);
    inline for (precs) |p| {
        const len = l1_cases[3].bytes / @sizeOf(p.F);
        const b = try Bufs(p.F).init(gpa, len, 1);
        defer b.deinit(gpa);
        // rstd's output
        const kx = try gpa.alloc(p.F, b.x.len);
        defer gpa.free(kx);
        const ky = try gpa.alloc(p.F, b.y.len);
        defer gpa.free(ky);
        inline for (comptime std.meta.tags(L1)) |r| {
            if (kernels.l1Wanted(r)) try raceLibL1(io, r, p, b, kx, ky);
        }
    }
}

/// rstd's largest difference from a library's result: `b` holds the
/// library's output, `kx` and `ky` rstd's; `rs` and `ls` are reductions'.
pub fn diffL1(comptime r: L1, comptime p: Prec, n: usize, inc: usize, b: Bufs(p.F), kx: []const p.F, ky: []const p.F, rs: [2]f64, ls: [2]f64) Diff {
    const eps: f64 = std.math.floatEps(p.F);
    switch (r) {
        .iamax => return .{ .eps = if (rs[0] == ls[0]) 0 else std.math.inf(f64), .bound = 0 },
        .dot, .nrm2, .asum, .sum => {
            const s = l1Scalar(r, p, b, n, inc).scale;
            return .{ .eps = @max(@abs(rs[0] - ls[0]), @abs(rs[1] - ls[1])) / (eps * @max(s, 1e-300)), .bound = sumBound(n) };
        },
        .axpy, .scal, .copy, .swap, .rot => {
            var worst: f64 = 0;
            for (0..n) |k| {
                for (0..p.cs()) |q| {
                    const w = l1Elem(r, p, b, k, inc, q);
                    const at = k * inc * p.cs() + q;
                    const dx = @abs(@as(f64, kx[at]) - b.x[at]) / (eps * @max(w[2], 1e-300));
                    const dy = @abs(@as(f64, ky[at]) - b.y[at]) / (eps * @max(w[3], 1e-300));
                    worst = @max(worst, dx, dy);
                }
            }
            return .{ .eps = worst, .bound = if (r == .copy or r == .swap) 0 else 8 };
        },
    }
}

pub fn raceLibL1(io: Io, comptime r: L1, comptime p: Prec, b: Bufs(p.F), kx: []p.F, ky: []p.F) !void {
    const name = comptime l1Slot(r, p);
    const C = L1Call(r, p);
    const asm_k = if (arm and !bb.alone) asmFor(name, r != .sum) else null;
    const fns = [lib_cols.len]?Fn(name){
        libL1(.openblas, r, p),
        libL1(.rstd, r, p),
        libL1(.blis, r, p),
        libL1(.accelerate, r, p),
        if (asm_k) |a| a.f else null,
        if (arm) candidate(r, p, 4) else null,
    };
    for (l1_cases) |case| {
        const n = case.bytes / p.es() / case.inc;
        if (n > lb.n_max) continue;
        const span = n * case.inc * p.cs();
        var es: Entrants(C, lib_cols.len) = .{};
        var ok = true;
        var rs: [2]f64 = .{ 0, 0 };
        for (0..lib_cols.len) |o| {
            const i = rstdFirst(o);
            if (!lib_cols[i].shown()) continue;
            const f = fns[i] orelse continue;
            const e = checkL1(r, p, f, n, case.inc, b);
            if (lb.lib_err) |err| return err;
            if (e > 1) {
                std.debug.print("WRONG: {s} {s} {s} n={d} inc={d}: {d:.2} of tolerance\n", .{ @tagName(r), p.name, lib_cols[i].label, n, case.inc, e });
                ok = false;
            }
            const res: [2]f64 = if (comptime reduces(r)) l1Result(r, p, f, n, case.inc, b) else .{ 0, 0 };
            if (i == rstd_col) {
                rs = res;
                var h: std.hash.Wyhash = .init(0);
                if (comptime reduces(r)) h.update(std.mem.asBytes(&res)) else {
                    @memcpy(kx[0..span], b.x[0..span]);
                    @memcpy(ky[0..span], b.y[0..span]);
                    h.update(std.mem.sliceAsBytes(b.x[0..span]));
                    h.update(std.mem.sliceAsBytes(b.y[0..span]));
                }
                es.mine = .{ .err = e * @as(f64, if (comptime reduces(r)) @floatFromInt(n + 8) else 8), .bits = h.final() };
            } else if (lib_cols[i].lib != null) es.diffs[i] = diffL1(r, p, n, case.inc, b, kx, ky, rs, res);
            _ = es.add(i, .{ .f = f, .n = @intCast(n), .inc = @intCast(case.inc), .x = b.x.ptr, .y = b.y.ptr });
        }
        if (!ok) continue;
        b.reset();
        const rates = try es.race(io, l1_timing, l1Bytes(r, n, p.es()));
        printLead("| {s} | {s} | {d} | {d} | {s} | ", .{ @tagName(r), p.name, n, case.inc, rstdRuns(name) });
        if (arm) std.debug.print("{s} | ", .{if (asm_k) |a| a.from else "none"});
        es.print(&lib_cols, &rates);
    }
}
