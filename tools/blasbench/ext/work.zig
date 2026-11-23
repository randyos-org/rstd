//! An entry's operands at its largest order, and their preparation.

const std = @import("std");
const rstd = @import("rstd");
const bb = @import("../blas_bench.zig");
const xr = @import("../ext.zig");
const pr = @import("pair.zig");

const lapack = rstd.math.kernels.lapack;
const Prec = bb.Prec;
const Buf = xr.Buf;
const Fill = xr.Fill;
const Level = xr.Level;
const Prep = xr.Prep;
const Restore = xr.Restore;
const Pair = pr.Pair;
const cconj = pr.cconj;
const cdiv = pr.cdiv;
const cmul = pr.cmul;
const csub = pr.csub;
const decl = pr.decl;
const lenA = pr.lenA;
const lenV = pr.lenV;
const outLen = pr.outLen;

/// An entry's operands at its largest order, shared by its contestants.
pub fn Work(comptime p: Prec) type {
    return struct {
        const W = @This();
        const F = p.F;
        const cs = p.cs();

        n: usize = 0,
        a: []align(64) F,
        a0: []align(64) F,
        orig: []align(64) F,
        b: []align(64) F,
        b0: []align(64) F,
        c: []align(64) F,
        c0: []align(64) F,
        s: []align(64) F,
        u: []align(64) F,
        vt: []align(64) F,
        work: []align(64) F,
        rwork: []align(64) F,
        ipiv: []c_int,
        ipiv0: []c_int,
        /// rstd's pivots, 0-based (issue #398)
        rpiv: []usize,
        rpiv0: []usize,
        /// `sytrf`'s, 0-based with their block of `D` (issue #314)
        spiv: []lapack.Pivot,
        iwork: []c_int,
        info: c_int = 0,

        pub fn init(gpa: std.mem.Allocator, comptime E: type, comptime lv: Level, n: usize) !W {
            const la = lenA(E, n) * cs;
            const lb = lenV(E, n) * cs;
            const big = lv == .lapack;
            var w: W = undefined;
            w.n = 0;
            w.info = 0;
            inline for (.{ .{ "a", la }, .{ "a0", la }, .{ "orig", la }, .{ "b", lb }, .{ "b0", lb }, .{ "c", lb }, .{ "c0", lb } }) |f|
                @field(w, f[0]) = try gpa.alignedAlloc(F, .@"64", f[1]);
            w.s = try gpa.alignedAlloc(F, .@"64", n * cs);
            const sq = if (big) n * n * cs else 0;
            w.u = try gpa.alignedAlloc(F, .@"64", sq);
            w.vt = try gpa.alignedAlloc(F, .@"64", sq);
            w.work = try gpa.alignedAlloc(F, .@"64", if (big) lwork(n) * cs else 0);
            w.rwork = try gpa.alignedAlloc(F, .@"64", if (big) lrwork(n) else 0);
            w.iwork = try gpa.alloc(c_int, if (big) liwork(n) else 0);
            w.ipiv = try gpa.alloc(c_int, n);
            w.ipiv0 = try gpa.alloc(c_int, n);
            w.rpiv = try gpa.alloc(usize, n);
            w.rpiv0 = try gpa.alloc(usize, n);
            w.spiv = try gpa.alloc(lapack.Pivot, n);
            return w;
        }

        pub fn deinit(w: *W, gpa: std.mem.Allocator) void {
            inline for (.{ "a", "a0", "orig", "b", "b0", "c", "c0", "s", "u", "vt", "work", "rwork", "iwork", "ipiv", "ipiv0", "rpiv", "rpiv0", "spiv" }) |f| gpa.free(@field(w, f));
        }

        /// Workspace for every routine's blocked path: past each one's
        /// documented need (`gesdd`'s 4n² + 7n, `syevd`'s 2n² + 6n + 1,
        /// `n nb` with `nb` up to 64).
        pub fn lwork(n: usize) usize {
            return 5 * n * n + 128 * n;
        }
        pub fn lrwork(n: usize) usize {
            return 5 * n * n + 64 * n;
        }
        pub fn liwork(n: usize) usize {
            return 8 * n + 64;
        }

        pub fn ni(w: *const W) c_int {
            return @intCast(w.n);
        }
        pub fn in(w: *const W) isize {
            return @intCast(w.n);
        }
        pub fn lw(w: *const W) c_int {
            return @intCast(lwork(w.n));
        }
        pub fn lrw(w: *const W) c_int {
            return @intCast(lrwork(w.n));
        }
        pub fn liw(w: *const W) c_int {
            return @intCast(liwork(w.n));
        }

        pub fn buf(w: *const W, comptime b: Buf) []const F {
            return @field(w, @tagName(b));
        }

        /// Element k of `b`, as (re, im); `s` is real.
        pub fn at(w: *const W, comptime b: Buf, k: usize) Pair {
            const v = w.buf(b);
            if (b == .s or !p.complex) return .{ v[k], 0 };
            return .{ v[2 * k], v[2 * k + 1] };
        }

        pub fn set(v: []F, k: usize, x: Pair) void {
            if (p.complex) {
                v[2 * k] = @floatCast(x[0]);
                v[2 * k + 1] = @floatCast(x[1]);
            } else v[k] = @floatCast(x[0]);
        }

        /// `E`'s output at the current order.
        pub fn out(w: *const W, comptime E: type) []const F {
            return w.buf(E.out)[0..outLen(E, p, w.n)];
        }

        /// Every operand fresh, before a check.
        pub fn reset(w: *W, comptime E: type) void {
            const la = lenA(E, w.n) * cs;
            const lb = lenV(E, w.n) * cs;
            @memcpy(w.a[0..la], w.a0[0..la]);
            @memcpy(w.b[0..lb], w.b0[0..lb]);
            @memcpy(w.c[0..lb], w.c0[0..lb]);
            @memcpy(w.ipiv[0..w.n], w.ipiv0[0..w.n]);
            @memcpy(w.rpiv[0..w.n], w.rpiv0[0..w.n]);
            w.info = 0;
        }

        /// What a call of `E` overwrites, fresh; inside its timing.
        pub fn restore(w: *W, comptime E: type) void {
            const r = decl(E, "restore", Restore{});
            if (r.a) {
                const la = lenA(E, w.n) * cs;
                @memcpy(w.a[0..la], w.a0[0..la]);
            }
            const lb = lenV(E, w.n) * cs;
            if (r.b) @memcpy(w.b[0..lb], w.b0[0..lb]);
            if (r.c) @memcpy(w.c[0..lb], w.c0[0..lb]);
        }

        /// The operands at order n: random, as `E`'s fill, and prepared.
        pub fn setup(w: *W, gpa: std.mem.Allocator, comptime E: type, n: usize) !void {
            w.n = n;
            var rng = std.Random.DefaultPrng.init(5);
            const r = rng.random();
            const la = lenA(E, n) * cs;
            const lb = lenV(E, n) * cs;
            for (w.orig[0..la]) |*v| v.* = @floatCast(r.float(f64) * 2 - 1);
            for (w.b0[0..lb]) |*v| v.* = @floatCast(r.float(f64) * 2 - 1);
            for (w.c0[0..lb]) |*v| v.* = @floatCast(r.float(f64) * 2 - 1);
            const fill = comptime decl(E, "fill", Fill.general);
            if (fill != .general) {
                if (E.shape == .band) @compileError("a band is filled general");
                if (fill == .spd) {
                    for (0..n) |j| for (0..j) |i| set(w.orig, j + i * n, cconj(w.at(.orig, i + j * n)));
                    for (0..n) |i| set(w.orig, i + i * n, .{ w.at(.orig, i + i * n)[0], 0 });
                }
                const two_n: F = @floatFromInt(2 * n);
                for (0..n) |i| w.orig[(i + i * n) * cs] += two_n;
            }
            @memcpy(w.a0[0..la], w.orig[0..la]);
            for (w.ipiv0[0..n], w.rpiv0[0..n], 0..) |*d, *e, k| {
                d.* = @intCast(k + 1);
                e.* = k;
            }
            const prep = comptime decl(E, "prep", Prep.none);
            if (prep != .none) try prepare(gpa, p, w, prep);
        }
    };
}

/// `prep`'s factorization of `orig` into `a0`, unblocked in f64: `potrf`
/// writes `U` over the upper triangle, `getrf` `L` and `U` over all of it.
pub fn prepare(gpa: std.mem.Allocator, comptime p: Prec, w: *Work(p), comptime what: Prep) !void {
    const n = w.n;
    const m = try gpa.alloc(Pair, n * n);
    defer gpa.free(m);
    for (m, 0..) |*v, k| v.* = w.at(.orig, k);
    switch (what) {
        .none => return,
        .potrf => for (0..n) |j| {
            for (0..j + 1) |i| {
                var s = m[i + j * n];
                for (0..i) |k| s = csub(s, cmul(cconj(m[k + i * n]), m[k + j * n]));
                m[i + j * n] = if (i < j) cdiv(s, m[i + i * n]) else .{ @sqrt(s[0]), 0 };
            }
        },
        .getrf => for (0..n) |k| {
            const piv = m[k + k * n];
            for (k + 1..n) |i| m[i + k * n] = cdiv(m[i + k * n], piv);
            for (k + 1..n) |j| {
                const u = m[k + j * n];
                for (k + 1..n) |i| m[i + j * n] = csub(m[i + j * n], cmul(m[i + k * n], u));
            }
        },
    }
    for (0..n) |j| for (0..n) |i| {
        if (what == .potrf and i > j) continue;
        Work(p).set(w.a0, i + j * n, m[i + j * n]);
    };
}
