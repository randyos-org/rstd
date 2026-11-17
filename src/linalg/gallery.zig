//! Test matrices (MATLAB's `magic`, `hilb`, `invhilb`, `rand` and `randn`),
//! for tests and examples, as `Mat(T)` of a float or complex `T`. Random
//! fills draw from the caller's `std.Random`. Replaces matlab88's
//! `generators.magic`, `hilber` and `random.Urand` (issue #330).

const std = @import("std");
const matrix = @import("matrix.zig");
const vector = @import("vector.zig");

const Allocator = std.mem.Allocator;
const Mat = matrix.Mat;
const View = matrix.View;
const isComplex = vector.isComplex;
const Real = vector.Real;

/// MATLAB's `magic(n)` (`magic.m`): rows, columns and both diagonals sum to
/// `n (n² + 1) / 2`, for `n` other than 2. Odd `n` by de la Loubère's
/// method, doubly even by the complement pattern, singly even by LUX.
pub fn magic(allocator: Allocator, comptime T: type, n: usize) Allocator.Error!Mat(T) {
    const m = try Mat(T).init(allocator, n, n);
    fillMagic(T, m, n);
    return m;
}

fn fillMagic(comptime T: type, m: Mat(T), n: usize) void {
    const set = struct {
        fn f(a: Mat(T), i: usize, j: usize, v: usize) void {
            a.data[i + j * a.rows] = scalar(T, @floatFromInt(v));
        }
    }.f;
    if (n == 0) return;
    if (n % 2 == 1) {
        // 1-based: n A + B + 1, A = (i + j - (n + 3) / 2) mod n, B = (i + 2 j - 2) mod n
        for (0..n) |j| for (0..n) |i| {
            const a = (i + j + 2 + n - (n + 3) / 2) % n;
            const b = (i + 2 * j + 1) % n;
            set(m, i, j, n * a + b + 1);
        };
    } else if (n % 4 == 0) {
        for (0..n) |j| for (0..n) |i| {
            const v = i * n + j + 1;
            const keep = ((i + 1) % 4) / 2 == ((j + 1) % 4) / 2;
            set(m, i, j, if (keep) n * n + 1 - v else v);
        };
    } else {
        // the odd square of order p in each quadrant, offset by 0, 2p², 3p², p²
        const p = n / 2;
        for (0..p) |j| for (0..p) |i| {
            const a = (i + j + 2 + p - (p + 3) / 2) % p;
            const b = (i + 2 * j + 1) % p;
            const v = p * a + b + 1;
            set(m, i, j, v);
            set(m, i, j + p, v + 2 * p * p);
            set(m, i + p, j, v + 3 * p * p);
            set(m, i + p, j + p, v + p * p);
        };
        if (n == 2) return;
        // swap the halves of the first k and last k - 1 columns, then two entries
        const k = (n - 2) / 4;
        for (0..n) |j| {
            if (j >= k and j < n - k + 1) continue;
            for (0..p) |i| std.mem.swap(T, &m.data[i + j * n], &m.data[i + p + j * n]);
        }
        for ([_]usize{ 0, k }) |j| std.mem.swap(T, &m.data[k + j * n], &m.data[k + p + j * n]);
    }
}

/// The Hilbert matrix, `1 / (i + j + 1)` 0-based: famously ill-conditioned.
pub fn hilbert(allocator: Allocator, comptime T: type, n: usize) Allocator.Error!Mat(T) {
    const m = try Mat(T).init(allocator, n, n);
    for (0..n) |j| for (0..n) |i| {
        m.data[i + j * n] = scalar(T, 1 / @as(Real(T), @floatFromInt(i + j + 1)));
    };
    return m;
}

/// The exact inverse of `hilbert(n)` by its closed form (MATLAB's `invhilb`,
/// matlab88's `HILBER`): integers, exact while they fit `T`'s significand
/// (to `n = 13` in f64).
pub fn invHilbert(allocator: Allocator, comptime T: type, n: usize) Allocator.Error!Mat(T) {
    const F = Real(T);
    const m = try Mat(T).init(allocator, n, n);
    const nf: F = @floatFromInt(n);
    var p = nf;
    for (1..n + 1) |i| {
        const fi: F = @floatFromInt(i);
        if (i > 1) p = (nf - fi + 1) * p * (nf + fi - 1) / ((fi - 1) * (fi - 1));
        var r = p * p;
        m.data[(i - 1) * (n + 1)] = scalar(T, r / (2 * fi - 1));
        for (i + 1..n + 1) |j| {
            const fj: F = @floatFromInt(j);
            r = -((nf - fj + 1) * r * (nf + fj - 1)) / ((fj - 1) * (fj - 1));
            const v = scalar(T, r / (fi + fj - 1));
            m.data[(i - 1) + (j - 1) * n] = v;
            m.data[(j - 1) + (i - 1) * n] = v;
        }
    }
    return m;
}

/// Uniform in [0, 1), each part of a complex value drawn alike (MATLAB's
/// `rand`).
pub fn uniform(allocator: Allocator, comptime T: type, rows: usize, cols: usize, rnd: std.Random) Allocator.Error!Mat(T) {
    const m = try Mat(T).init(allocator, rows, cols);
    fillUniform(T, m.view(), rnd);
    return m;
}

/// Standard normal, each part of a complex value drawn alike (MATLAB's
/// `randn` per part).
pub fn normal(allocator: Allocator, comptime T: type, rows: usize, cols: usize, rnd: std.Random) Allocator.Error!Mat(T) {
    const m = try Mat(T).init(allocator, rows, cols);
    fillNormal(T, m.view(), rnd);
    return m;
}

/// `uniform`'s values into a view, column by column.
pub fn fillUniform(comptime T: type, v: View(T), rnd: std.Random) void {
    for (0..v.cols) |j| for (v.col(j)) |*x| {
        x.* = if (comptime isComplex(T)) T.init(draw(Real(T), rnd, false), draw(Real(T), rnd, false)) else draw(T, rnd, false);
    };
}

/// `normal`'s values into a view, column by column.
pub fn fillNormal(comptime T: type, v: View(T), rnd: std.Random) void {
    for (0..v.cols) |j| for (v.col(j)) |*x| {
        x.* = if (comptime isComplex(T)) T.init(draw(Real(T), rnd, true), draw(Real(T), rnd, true)) else draw(T, rnd, true);
    };
}

/// `std.Random` serves f32 and f64; other widths take an f64 draw.
fn draw(comptime F: type, rnd: std.Random, comptime norm: bool) F {
    const G = if (F == f32) f32 else f64;
    const x = if (norm) rnd.floatNorm(G) else rnd.float(G);
    return @floatCast(x);
}

fn scalar(comptime T: type, x: Real(T)) T {
    return if (comptime isComplex(T)) T.init(x, 0) else x;
}

// ---- tests ----

const testing = std.testing;

test "magic: MATLAB's squares, every sum equal" {
    const al = testing.allocator;
    const m3 = try magic(al, f64, 3);
    defer m3.deinit(al);
    try testing.expectEqualSlices(f64, &.{ 8, 3, 4, 1, 5, 9, 6, 7, 2 }, m3.data);
    const m4 = try magic(al, f64, 4);
    defer m4.deinit(al);
    try testing.expectEqualSlices(f64, &.{ 16, 5, 9, 4, 2, 11, 7, 14, 3, 10, 6, 15, 13, 8, 12, 1 }, m4.data);
    const m6 = try magic(al, f64, 6);
    defer m6.deinit(al);
    try testing.expectEqualSlices(f64, &.{ 35, 3, 31, 8, 30, 4 }, m6.data[0..6]);

    for (1..15) |n| {
        if (n == 2) continue;
        const m = try magic(al, f64, n);
        defer m.deinit(al);
        const want: f64 = @floatFromInt(n * (n * n + 1) / 2);
        var d1: f64 = 0;
        var d2: f64 = 0;
        var seen = try al.alloc(bool, n * n);
        defer al.free(seen);
        @memset(seen, false);
        for (0..n) |i| {
            var r: f64 = 0;
            var c: f64 = 0;
            for (0..n) |j| {
                r += m.at(i, j);
                c += m.at(j, i);
                const v: usize = @intFromFloat(m.at(i, j));
                try testing.expect(v >= 1 and v <= n * n and !seen[v - 1]);
                seen[v - 1] = true;
            }
            try testing.expectEqual(want, r);
            try testing.expectEqual(want, c);
            d1 += m.at(i, i);
            d2 += m.at(i, n - 1 - i);
        }
        try testing.expectEqual(want, d1);
        try testing.expectEqual(want, d2);
    }
}

test "hilbert and invHilbert: exact inverse, H H⁻¹ = I" {
    const al = testing.allocator;
    const inv3 = try invHilbert(al, f64, 3);
    defer inv3.deinit(al);
    try testing.expectEqualSlices(f64, &.{ 9, -36, 30, -36, 192, -180, 30, -180, 180 }, inv3.data);
    const n = 5;
    const h = try hilbert(al, f64, n);
    defer h.deinit(al);
    const hi = try invHilbert(al, f64, n);
    defer hi.deinit(al);
    for (0..n) |i| for (0..n) |j| {
        var s: f64 = 0;
        for (0..n) |k| s += h.at(i, k) * hi.at(k, j);
        try testing.expectApproxEqAbs(@as(f64, if (i == j) 1 else 0), s, 1e-9);
    };
    const hc = try hilbert(al, std.math.Complex(f32), 2);
    defer hc.deinit(al);
    try testing.expectEqual(@as(f32, 0.5), hc.at(0, 1).re);
}

test "uniform and normal: the caller's generator, every element drawn" {
    const al = testing.allocator;
    var prng: std.Random.DefaultPrng = .init(27);
    const u = try uniform(al, f64, 40, 25, prng.random());
    defer u.deinit(al);
    var mean: f64 = 0;
    for (u.data) |x| {
        try testing.expect(x >= 0 and x < 1);
        mean += x;
    }
    try testing.expectApproxEqAbs(0.5, mean / 1000, 0.05);
    const z = try normal(al, std.math.Complex(f64), 30, 30, prng.random());
    defer z.deinit(al);
    var m2: f64 = 0;
    for (z.data) |x| m2 += x.re * x.re + x.im * x.im;
    try testing.expectApproxEqAbs(2, m2 / 900, 0.3);
    // the same seed, the same matrix; f16 through an f64 draw
    var a: std.Random.DefaultPrng = .init(5);
    var b: std.Random.DefaultPrng = .init(5);
    const x = try normal(al, f16, 3, 3, a.random());
    defer x.deinit(al);
    const y = try normal(al, f16, 3, 3, b.random());
    defer y.deinit(al);
    try testing.expectEqualSlices(f16, x.data, y.data);
}
