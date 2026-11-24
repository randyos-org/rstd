//! `lib2`'s `gemv` race.

const std = @import("std");
const bb = @import("../blas_bench.zig");
const kernels = @import("../kernels.zig");
const ext = @import("../ext.zig");
const lb = @import("../libs.zig");
const lv = @import("via.zig");
const tb = @import("table.zig");

const Io = std.Io;
const arm = bb.arm;
const clock = bb.clock;
const fillRandom = bb.fillRandom;
const Fn = bb.Fn;
const l1_timing = bb.l1_timing;
const Prec = bb.Prec;
const precs = bb.precs;
const asmFor = kernels.asmFor;
const checkGemv = kernels.checkGemv;
const gemv_alpha = kernels.gemv_alpha;
const GemvCall = kernels.GemvCall;
const gemvCandidate = kernels.gemvCandidate;
const l2_bytes = kernels.l2_bytes;
const rounded = kernels.rounded;
const ext_l2 = ext.ext_l2;
const extTable = ext.extTable;
const libGemv = lv.libGemv;
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

pub fn libLevel2(io: Io, gpa: std.mem.Allocator) !void {
    sectionMark("lib2");
    try clock(io, "library level 2 start");
    std.debug.print(
        \\
        \\GFLOP/s, square, unit strides, `beta` 1, best of {d} batches of at least {d} ms; ratios are over OpenBLAS. {s} {s}
        \\
        \\
    , .{ l1_timing.batches, @divTrunc(l1_timing.min_ns, std.time.ns_per_ms), lb.io_note, diff_note });
    if (arm)
        printHeader("routine | p | n | rstd runs | from", "------- | - | -: | --------- | ----", &lib_cols)
    else
        printHeader("routine | p | n | rstd runs", "------- | - | -: | ---------", &lib_cols);
    inline for (precs) |p| {
        inline for (.{ "n", "t" }) |t| try raceLibGemv(io, gpa, p, t);
    }
    try extTable(io, gpa, .l2, ext_l2);
}

/// The size of the terms of each float of `y0 + alpha*op(A)*x`, as
/// `checkGemv` has it.
pub fn gemvScale(comptime p: Prec, trans: bool, n: usize, a: []const p.F, x: []const p.F, y0: []const p.F, out: []f64) void {
    const cs = p.cs();
    const g = @abs(rounded(p, gemv_alpha));
    const al = if (p.complex) 2 * g else g;
    for (0..n) |i| {
        var mag: f64 = 0;
        for (0..n) |l| {
            const ai = if (trans) (l + i * n) * cs else (i + l * n) * cs;
            mag += magAt(p, a, ai) * magAt(p, x, l * cs);
        }
        for (0..cs) |q| out[i * cs + q] = @abs(y0[i * cs + q]) + al * mag;
    }
}

/// `|re| + |im|` of the element at float `i`.
pub fn magAt(comptime p: Prec, v: []const p.F, i: usize) f64 {
    return if (p.complex) @as(f64, @abs(v[i])) + @abs(v[i + 1]) else @abs(v[i]);
}

/// The largest `|want - got| / (eps * scale)` over `scale`'s floats.
pub fn diffScaled(comptime F: type, want: []const F, got: []const F, scale: []const f64) f64 {
    const eps: f64 = std.math.floatEps(F);
    var worst: f64 = 0;
    for (want[0..scale.len], got[0..scale.len], scale) |w, g, s| worst = @max(worst, @abs(@as(f64, w) - g) / (eps * @max(s, 1e-300)));
    return worst;
}

pub fn raceLibGemv(io: Io, gpa: std.mem.Allocator, comptime p: Prec, comptime t: []const u8) !void {
    const F = p.F;
    const name = p.name ++ "gemv_" ++ t;
    const C = GemvCall(p, name);
    const cs = p.cs();
    const max_n: usize = std.math.sqrt(l2_bytes[2] / p.es());
    const a = try gpa.alignedAlloc(F, .@"64", max_n * max_n * cs);
    defer gpa.free(a);
    const x = try gpa.alignedAlloc(F, .@"64", max_n * cs);
    defer gpa.free(x);
    const y = try gpa.alignedAlloc(F, .@"64", max_n * cs + 64);
    defer gpa.free(y);
    const y0 = try gpa.alloc(F, y.len);
    defer gpa.free(y0);
    const ky = try gpa.alloc(F, y.len);
    defer gpa.free(ky);
    const scale = try gpa.alloc(f64, max_n * cs);
    defer gpa.free(scale);
    const buf = try gpa.alignedAlloc(F, .@"64", 4096);
    defer gpa.free(buf);
    var rng = std.Random.DefaultPrng.init(2);
    const rand = rng.random();
    fillRandom(F, a, rand);
    fillRandom(F, x, rand);
    fillRandom(F, y0, rand);

    const asm_k = if (arm and !bb.alone) asmFor(name, true) else null;
    const fns = [lib_cols.len]?Fn(name){
        libGemv(.openblas, p, t),
        libGemv(.rstd, p, t),
        libGemv(.blis, p, t),
        libGemv(.accelerate, p, t),
        if (asm_k) |k| k.f else null,
        if (arm) gemvCandidate(p, t, 4) else null,
    };
    for (l2_bytes) |bytes| {
        const n: usize = std.math.sqrt(bytes / p.es());
        if (n > lb.n_max) continue;
        const out = n * cs;
        gemvScale(p, t[0] == 't', n, a, x, y0, scale[0..out]);
        var es: Entrants(C, lib_cols.len) = .{};
        var ok = true;
        for (0..lib_cols.len) |o| {
            const i = rstdFirst(o);
            if (!lib_cols[i].shown()) continue;
            const f = fns[i] orelse continue;
            const ctx = es.add(i, .{ .f = f, .n = @intCast(n), .a = a.ptr, .x = x.ptr, .y = y.ptr, .buf = buf.ptr });
            @memcpy(y, y0);
            ctx.call();
            if (lb.lib_err) |err| return err;
            const e = checkGemv(p, t[0] == 't', n, a, x, y, y0);
            if (e > 1) {
                std.debug.print("WRONG: gemv_{s} {s} {s} n={d}: {d:.2} of tolerance\n", .{ t, p.name, lib_cols[i].label, n, e });
                ok = false;
            }
            if (i == rstd_col) {
                @memcpy(ky, y);
                es.mine = .{ .err = e * @as(f64, @floatFromInt(n + 8)), .bits = std.hash.Wyhash.hash(0, std.mem.sliceAsBytes(y[0..out])) };
            } else if (lib_cols[i].lib != null) es.diffs[i] = .{ .eps = diffScaled(F, ky, y, scale[0..out]), .bound = sumBound(n) };
        }
        if (!ok) continue;
        @memcpy(y, y0);
        const nf: f64 = @floatFromInt(n);
        const rates = try es.race(io, l1_timing, nf * nf * @as(f64, if (p.complex) 8 else 2));
        printLead("| gemv_{s} | {s} | {d} | {s} | ", .{ t, p.name, n, rstdRuns(name) });
        if (arm) std.debug.print("{s} | ", .{if (asm_k) |k| k.from else "none"});
        es.print(&lib_cols, &rates);
    }
}
