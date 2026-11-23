//! One entry's race: its contestants, their residuals and agreement.

const std = @import("std");
const bb = @import("../blas_bench.zig");
const libs = @import("../libs.zig");
const xr = @import("../ext.zig");
const pr = @import("pair.zig");
const wk = @import("work.zig");

const Io = std.Io;
const l1_timing = bb.l1_timing;
const l3_timing = bb.l3_timing;
const Prec = bb.Prec;
const fail = libs.fail;
const lookupF = libs.lookupF;
const rstd_col = libs.rstd_col;
const rstdFirst = libs.rstdFirst;
const sumBound = libs.sumBound;
const Buf = xr.Buf;
const Level = xr.Level;
const ext_cols = xr.ext_cols;
const Pair = pr.Pair;
const cconj = pr.cconj;
const decl = pr.decl;
const hasRstd = pr.hasRstd;
const orders = pr.orders;
const outDims = pr.outDims;
const outLen = pr.outLen;
const symOf = pr.symOf;
const terms = pr.terms;
const wrongLimit = pr.wrongLimit;
const Work = wk.Work;

/// A contestant's largest residual over `E`'s output in eps of the type
/// relative to the size of its terms: every entry of a small output, else
/// 256 sampled; null without `residual`.
pub fn checkExt(comptime E: type, comptime p: Prec, w: *const Work(p)) ?f64 {
    if (!@hasDecl(E, "residual")) return null;
    const eps: f64 = std.math.floatEps(p.F);
    const d = outDims(E, w.n);
    var worst: f64 = 0;
    if (d[0] * d[1] <= 4096) {
        for (0..d[1]) |j| for (0..d[0]) |i| {
            const r = E.residual(p, w, i, j) orelse continue;
            worst = @max(worst, r[0] / (eps * @max(r[1], 1e-300)));
        };
    } else {
        var rng = std.Random.DefaultPrng.init(7);
        const rand = rng.random();
        var hits: usize = 0;
        var tries: usize = 0;
        while (hits < 256 and tries < 4096) : (tries += 1) {
            const r = E.residual(p, w, rand.uintLessThan(usize, d[0]), rand.uintLessThan(usize, d[1])) orelse continue;
            hits += 1;
            worst = @max(worst, r[0] / (eps * @max(r[1], 1e-300)));
        }
    }
    return if (std.math.isNan(worst)) std.math.inf(f64) else worst;
}

/// rstd's largest difference from `got` in eps relative to its own largest
/// entry.
pub fn diffNorm(comptime F: type, mine: []const F, got: []const F) f64 {
    var big: f64 = 0;
    var worst: f64 = 0;
    for (mine, got) |a, b| {
        big = @max(big, @abs(@as(f64, a)));
        worst = @max(worst, @abs(@as(f64, a) - b));
    }
    const r = worst / (std.math.floatEps(F) * @max(big, 1e-300));
    return if (std.math.isNan(r)) std.math.inf(f64) else r;
}

pub fn ExtCall(comptime E: type, comptime p: Prec) type {
    return struct {
        /// null: rstd's
        f: ?*const E.Sig(p),
        w: *Work(p),
        pub fn call(s: *const @This()) void {
            s.w.restore(E);
            if (s.f) |f| return E.lib(p, f, s.w);
            if (comptime hasRstd(E)) E.rstd(p, s.w) catch |e| fail(e);
        }
        pub fn run(o: *anyopaque) void {
            const s: *@This() = @ptrCast(@alignCast(o));
            s.call();
        }
    };
}

pub fn raceExt(io: Io, gpa: std.mem.Allocator, comptime E: type, comptime table_level: Level, comptime p: Prec) !void {
    const lv: Level = comptime decl(E, "level", table_level);
    const ns = comptime orders(E, lv, p);
    var w = try Work(p).init(gpa, E, lv, ns[ns.len - 1]);
    defer w.deinit(gpa);
    const kout = try gpa.alloc(p.F, outLen(E, p, ns[ns.len - 1]));
    defer gpa.free(kout);
    const Sig = *const E.Sig(p);
    const sym = comptime symOf(E, p);
    const fns = [ext_cols.len]?Sig{
        lookupF(.openblas, Sig, sym),
        null,
        lookupF(.blis, Sig, sym),
        lookupF(.accelerate, Sig, sym),
        lookupF(.reference, Sig, sym),
        lookupF(.eigen, Sig, sym),
        lookupF(.aoclflame, Sig, sym),
    };
    const EC = ExtCall(E, p);
    const timing = if (lv == .l2) l1_timing else l3_timing;
    for (ns) |n| {
        if (n > libs.n_max) continue;
        try w.setup(gpa, E, n);
        const len = outLen(E, p, n);
        var es: libs.Entrants(EC, ext_cols.len) = .{};
        var ok = true;
        var have_mine = false;
        for (0..ext_cols.len) |o| {
            const i = rstdFirst(o);
            const col = ext_cols[i];
            if (!col.shown() or (lv != .l2 and n > col.max_n)) continue;
            if (i == rstd_col) {
                if (!comptime hasRstd(E)) continue;
            } else if (fns[i] == null) continue;
            const ctx = es.add(i, .{ .f = fns[i], .w = &w });
            w.reset(E);
            ctx.call();
            if (libs.lib_err) |err| return err;
            const e = checkExt(E, p, &w);
            if (w.info != 0 or (e != null and e.? > wrongLimit(E, n))) {
                std.debug.print("WRONG: {s} {s} {s} n={d}: info {d}, residual {?d:.1} eps\n", .{ E.name, p.name, col.label, n, w.info, e });
                ok = false;
            }
            const got = w.out(E);
            if (i == rstd_col) {
                @memcpy(kout[0..len], got);
                es.mine = .{ .err = e, .bits = std.hash.Wyhash.hash(0, std.mem.sliceAsBytes(got)) };
                have_mine = true;
            } else if (have_mine and col.lib != null) es.diffs[i] = .{ .eps = diffNorm(p.F, kout[0..len], got), .bound = sumBound(terms(E, n)) };
        }
        if (!ok or es.count == 0) continue;
        w.reset(E);
        const rates = try es.race(io, timing, E.flops(@floatFromInt(n)) * @as(f64, if (p.complex) 4 else 1));
        std.debug.print("| {s} | {s} | {d} | ", .{ E.name, p.name, n });
        es.print(&ext_cols, &rates);
    }
}

/// `A[i, k]` of a symmetric (`herm` false) or Hermitian matrix whose upper
/// triangle `b` stores, order n.
pub fn symAt(comptime p: Prec, w: *const Work(p), comptime b: Buf, comptime herm: bool, i: usize, k: usize) Pair {
    const n = w.n;
    if (i <= k) return w.at(b, i + k * n);
    const v = w.at(b, k + i * n);
    return if (herm) cconj(v) else v;
}
