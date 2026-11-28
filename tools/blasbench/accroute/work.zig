//! The `accroute` section's operands: one case's vectors and matrices at one
//! size, their fill, and its sizes.

const std = @import("std");
const rstd = @import("rstd");
const bb = @import("../blas_bench.zig");
const kernels = @import("../kernels.zig");
const libs = @import("../libs.zig");

const blas = rstd.math.kernels.blas;
const Prec = bb.Prec;
const gemm_sizes = kernels.gemm_sizes;
const l2_bytes = kernels.l2_bytes;
const Num = libs.Num;

/// `vec`: `x` and `y` alone; `full`, `band` (upper, `band_k` above the
/// diagonal, `band_k + 1` rows), `packed_` (upper): `A` and `x`, `y` n long;
/// `mm`: `A`, `B` and `C` (in `y`) n by n.
pub const Shape = enum { vec, full, band, packed_, mm };

/// `dominant` adds twice the terms of a row to the diagonal, for the solves.
/// Every fill keeps the diagonal real, as a Hermitian `A` has it.
pub const Fill = enum { general, dominant };

/// What a case's result is: `A`, `x`, `y` (`C` in `mm`), or a scalar `r`.
pub const Out = enum { a, x, y, r };

/// What each timed call starts from fresh: `x` for the triangular level 2,
/// whose `x` would grow or shrink without bound.
pub const Restore = enum { none, x };

pub const Size = struct { n: usize, inc: usize = 1 };

/// Off-diagonals of a band.
pub const band_k = 32;
pub const band_rows = band_k + 1;

/// Level 1's vector bytes, unit strides (L1, L1, L2, memory), then stride 2.
const vec_unit = [_]usize{ 4 << 10, 32 << 10, 1 << 20, 64 << 20 };
const vec_strided = [_]usize{ 32 << 10, 1 << 20 };
/// Level 3's orders: `gemm_sizes` to 1024, where one complex call alone on
/// one core takes about a second.
const mm_sizes = gemm_sizes[0 .. gemm_sizes.len - 1];

/// A shape's sizes for `p`: level 1's by vector bytes, level 2's by `A`'s
/// (`l2_bytes`, as `lib2`).
pub fn sizes(comptime shape: Shape, comptime p: Prec) []const Size {
    return comptime blk: {
        const es = p.es();
        var out: []const Size = &.{};
        switch (shape) {
            .vec => {
                for (vec_unit) |b| out = out ++ [_]Size{.{ .n = b / es }};
                for (vec_strided) |b| out = out ++ [_]Size{.{ .n = b / es / 2, .inc = 2 }};
            },
            .full => for (l2_bytes) |b| {
                out = out ++ [_]Size{.{ .n = std.math.sqrt(b / es) }};
            },
            .band => for (l2_bytes) |b| {
                out = out ++ [_]Size{.{ .n = b / (es * band_rows) }};
            },
            .packed_ => for (l2_bytes) |b| {
                out = out ++ [_]Size{.{ .n = std.math.sqrt(2 * b / es) }};
            },
            .mm => for (mm_sizes) |n| {
                out = out ++ [_]Size{.{ .n = n }};
            },
        }
        const final = out[0..out.len].*;
        break :blk &final;
    };
}

/// `p`'s scalar (re, im); real types take `re`.
pub fn scalar(comptime p: Prec, re: f64, im: f64) Num(p) {
    return if (p.complex) .init(@floatCast(re), @floatCast(im)) else @floatCast(re);
}

pub fn Work(comptime p: Prec) type {
    return struct {
        const W = @This();
        pub const T = Num(p);

        shape: Shape,
        n: usize,
        inc: usize,
        a: []T,
        a0: []T,
        b: []T,
        x: []T,
        x0: []T,
        y: []T,
        y0: []T,
        r: p.F = 0,

        pub fn init(gpa: std.mem.Allocator, shape: Shape, sz: Size) !W {
            const n = sz.n;
            const la: usize = switch (shape) {
                .vec => 0,
                .full, .mm => n * n,
                .band => band_rows * n,
                .packed_ => n * (n + 1) / 2,
            };
            const lv: usize = switch (shape) {
                .vec => (n - 1) * sz.inc + 1,
                .full, .band, .packed_ => n,
                .mm => 0,
            };
            const lc: usize = if (shape == .mm) n * n else lv;
            var w: W = .{ .shape = shape, .n = n, .inc = sz.inc, .a = &.{}, .a0 = &.{}, .b = &.{}, .x = &.{}, .x0 = &.{}, .y = &.{}, .y0 = &.{} };
            errdefer w.deinit(gpa);
            w.a = try gpa.alloc(T, la);
            w.a0 = try gpa.alloc(T, la);
            w.b = try gpa.alloc(T, if (shape == .mm) n * n else 0);
            w.x = try gpa.alloc(T, lv);
            w.x0 = try gpa.alloc(T, lv);
            w.y = try gpa.alloc(T, lc);
            w.y0 = try gpa.alloc(T, lc);
            return w;
        }

        pub fn deinit(w: *W, gpa: std.mem.Allocator) void {
            inline for (.{ "a", "a0", "b", "x", "x0", "y", "y0" }) |f| gpa.free(@field(w, f));
        }

        /// Random entries in [-1, 1), the diagonal real, as `f` says.
        pub fn fill(w: *W, f: Fill) void {
            var rng = std.Random.DefaultPrng.init(w.n);
            const r = rng.random();
            for ([_][]T{ w.a0, w.b, w.x0, w.y0 }) |v| for (v) |*e| {
                e.* = scalar(p, r.float(f64) * 2 - 1, r.float(f64) * 2 - 1);
            };
            if (w.shape == .vec) return;
            const terms: usize = if (w.shape == .band) band_rows else w.n;
            const add: p.F = if (f == .dominant) @floatFromInt(2 * terms) else 0;
            for (0..w.n) |i| {
                const d = &w.a0[w.diagonal(i)];
                if (p.complex) d.* = .init(d.re + add, 0) else d.* += add;
            }
        }

        /// Where `A(i, i)` is stored.
        fn diagonal(w: *const W, i: usize) usize {
            return switch (w.shape) {
                .vec => unreachable,
                .full, .mm => i + i * w.n,
                .band => band_k + i * band_rows,
                .packed_ => i * (i + 1) / 2 + i,
            };
        }

        /// Every operand fresh.
        pub fn reset(w: *W) void {
            @memcpy(w.a, w.a0);
            @memcpy(w.x, w.x0);
            @memcpy(w.y, w.y0);
            w.r = 0;
        }

        pub fn restore(w: *W, what: Restore) void {
            if (what == .x) @memcpy(w.x, w.x0);
        }

        /// The result as floats of `p`.
        pub fn out(w: *const W, what: Out) []const p.F {
            const v: []const T = switch (what) {
                .a => w.a,
                .x => w.x,
                .y => w.y,
                .r => return @as([*]const p.F, @ptrCast(&w.r))[0..1],
            };
            return @as([*]const p.F, @ptrCast(v.ptr))[0 .. v.len * p.cs()];
        }

        pub fn xv(w: *W) blas.Strided(T) {
            return blas.Strided(T).init(w.x, w.n, @intCast(w.inc)) catch unreachable;
        }
        pub fn xc(w: *W) blas.ConstStrided(T) {
            return blas.ConstStrided(T).init(w.x, w.n, @intCast(w.inc)) catch unreachable;
        }
        pub fn yv(w: *W) blas.Strided(T) {
            return blas.Strided(T).init(w.y, w.n, @intCast(w.inc)) catch unreachable;
        }
        pub fn yc(w: *W) blas.ConstStrided(T) {
            return blas.ConstStrided(T).init(w.y, w.n, @intCast(w.inc)) catch unreachable;
        }
    };
}
