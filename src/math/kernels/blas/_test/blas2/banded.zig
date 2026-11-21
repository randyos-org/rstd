//! `blas2/banded.zig`'s tests, the rank updates' too (`rank.zig`).

const std = @import("std");
const blas = @import("../../__root__.zig");
const types = @import("../../types.zig");
const util = @import("../../util.zig");
const accel = @import("../../apple_accel.zig");
const ops = @import("../../blas2/ops.zig");
const Kind = @import("../../blas2/symmetric.zig").Kind;
const banded = @import("../../blas2/banded.zig");

const Num = ops.Num;
const Trans = types.Trans;
const Uplo = types.Uplo;
const isComplex = types.isComplex;
const Float = types.Float;
const zero = util.zero;
const gbmv = banded.gbmv;
const sbmv = banded.sbmv;
const hbmv = banded.hbmv;
const spmv = banded.spmv;
const hpmv = banded.hpmv;

const rank = @import("../../blas2/rank.zig");
const spr = rank.spr;
const syr = rank.syr;
const hpr = rank.hpr;
const her = rank.her;
const spr2 = rank.spr2;
const syr2 = rank.syr2;
const hpr2 = rank.hpr2;
const her2 = rank.her2;

const Vtable = @import("../../Vtable.zig");
const check = @import("../check.zig");
const Complex32 = types.Complex32;
const Complex64 = types.Complex64;

const all = check.all_types;
const talloc = std.testing.allocator;
const max_mem = 256;

/// Where logical element `i` of a vector of `n` with stride `inc` lives.
fn vecAt(i: usize, n: usize, inc: isize) usize {
    const s: usize = @abs(inc);
    return if (inc > 0) i * s else (n - 1 - i) * s;
}

const Store = enum { full, packed_, band };

/// Memory index of the stored entry `(i, j)`, or null when it is not in the
/// stored triangle or band.
fn storeIdx(comptime st: Store, lower: bool, n: usize, ld: usize, k: usize, i: usize, j: usize) ?usize {
    if (lower) {
        if (i < j) return null;
    } else if (i > j) return null;
    const d = if (lower) i - j else j - i;
    return switch (st) {
        .full => i + j * ld,
        .packed_ => if (lower) j * (2 * n - j + 1) / 2 + d else j * (j + 1) / 2 + i,
        .band => if (d > k) null else (if (lower) d else k - d) + j * ld,
    };
}

/// Entry `(i, j)` of the matrix a stored triangle stands for.
fn dense(comptime T: type, comptime st: Store, comptime kind: Kind, lower: bool, n: usize, ld: usize, k: usize, a: []const T, i: usize, j: usize) T {
    if (i == j) {
        const d = a[storeIdx(st, lower, n, ld, k, i, i).?];
        return if (kind == .symmetric) d else if (comptime isComplex(T)) T.init(d.re, 0) else d;
    }
    if (storeIdx(st, lower, n, ld, k, i, j)) |p| {
        return if (kind == .hermitian_rev) Num(T).conj(a[p], true) else a[p];
    }
    if (storeIdx(st, lower, n, ld, k, j, i)) |p| {
        return if (kind == .hermitian) Num(T).conj(a[p], true) else a[p];
    }
    return zero(T);
}

fn tol(comptime T: type, terms: usize) Float(T) {
    return check.tolerance(Float(T), 8 * terms);
}

const cases = [_]struct { n: usize, k: usize }{ .{ .n = 1, .k = 0 }, .{ .n = 2, .k = 1 }, .{ .n = 3, .k = 0 }, .{ .n = 5, .k = 2 }, .{ .n = 9, .k = 1 }, .{ .n = 9, .k = 4 }, .{ .n = 9, .k = 12 } };
const inc_pairs = [_][2]isize{ .{ 1, 1 }, .{ 2, 3 }, .{ -1, 1 }, .{ 1, -2 }, .{ -3, -2 } };

const scalarOf = check.val;

/// Checks `y := alpha*S*x + beta*y` against dense arithmetic.
fn checkMv(comptime T: type, comptime st: Store, comptime kind: Kind, lower: bool) !void {
    var seed: u64 = 1;
    for (cases) |c| for (inc_pairs) |inc| for ([_]f64{ 1, 0, -0.5 }) |bre| {
        seed += 1;
        const n = c.n;
        const k = if (st == .band) c.k else 0;
        const ld: usize = k + 1 + @as(usize, @intCast(seed & 1));
        var a: [max_mem]T = undefined;
        var xs: [max_mem]T = undefined;
        var ys: [max_mem]T = undefined;
        check.fill(T, &a, seed);
        check.fill(T, &xs, seed + 100);
        check.fill(T, &ys, seed + 200);
        const alpha = scalarOf(T, 1.25, -0.5);
        const beta = scalarOf(T, bre, if (bre == -0.5) 0.75 else 0);
        var want = ys;
        for (0..n) |i| {
            var acc = zero(T);
            for (0..n) |j| acc = Num(T).add(acc, Num(T).mul(dense(T, st, kind, lower, n, ld, k, &a, i, j), xs[vecAt(j, n, inc[0])], false));
            const yp = vecAt(i, n, inc[1]);
            want[yp] = Num(T).add(if (bre == 0) zero(T) else Num(T).mul(beta, ys[yp], false), Num(T).mul(alpha, acc, false));
        }
        if (bre == 0) for (0..n) |i| {
            ys[vecAt(i, n, inc[1])] = scalarOf(T, std.math.nan(f64), std.math.nan(f64));
        };
        const un: Uplo = if (lower) .l else .u;
        const hu: Uplo = switch (kind) {
            .hermitian_rev => if (lower) .m else .v,
            else => un,
        };
        const nn: isize = @intCast(n);
        const lda: isize = @intCast(ld);
        switch (st) {
            .band => if (kind == .symmetric)
                try sbmv(talloc, check.io, T, un, nn, @intCast(k), alpha, &a, lda, &xs, inc[0], beta, &ys, inc[1])
            else
                try hbmv(talloc, check.io, T, hu, nn, @intCast(k), alpha, &a, lda, &xs, inc[0], beta, &ys, inc[1]),
            .packed_ => if (kind == .symmetric)
                try spmv(talloc, check.io, T, un, nn, alpha, &a, &xs, inc[0], beta, &ys, inc[1])
            else
                try hpmv(talloc, check.io, T, un, nn, alpha, &a, &xs, inc[0], beta, &ys, inc[1]),
            .full => unreachable,
        }
        check.expectCloseSlices(T, &want, &ys, tol(T, n)) catch |e| {
            check.say("mv {s} {s} n={d} k={d} inc={d},{d} lower={}: ", .{ @tagName(st), @tagName(kind), n, k, inc[0], inc[1], lower });
            return e;
        };
    };
}

fn kindsOf(comptime T: type) []const Kind {
    return if (comptime isComplex(T)) &.{ .symmetric, .hermitian, .hermitian_rev } else &.{.symmetric};
}

/// Checks the rank-1 and rank-2 updates against dense arithmetic.
fn checkRank(comptime T: type, comptime st: Store, comptime herm: bool, lower: bool) !void {
    var seed: u64 = 1000;
    for (cases) |c| for (inc_pairs) |inc| for ([_]bool{ false, true }) |two| {
        seed += 1;
        const n = c.n;
        const ld: usize = n + @as(usize, @intCast(seed & 1));
        var a: [max_mem]T = undefined;
        var xs: [max_mem]T = undefined;
        var ys: [max_mem]T = undefined;
        check.fill(T, &a, seed);
        check.fill(T, &xs, seed + 100);
        check.fill(T, &ys, seed + 200);
        if (c.n > 2) xs[vecAt(1, n, inc[0])] = zero(T); // the zero-skip path
        const alpha = scalarOf(T, 1.25, -0.5);
        const alpha_f: Float(T) = 1.25;
        const a0 = alpha;
        var want = a;
        for (0..n) |c_| for (0..n) |r| {
            const p = storeIdx(st, lower, n, ld, 0, r, c_) orelse continue;
            const xr = xs[vecAt(r, n, inc[0])];
            const xc = xs[vecAt(c_, n, inc[0])];
            const yr = ys[vecAt(r, n, inc[1])];
            const yc = ys[vecAt(c_, n, inc[1])];
            const v = if (herm and two)
                Num(T).add(Num(T).mul(alpha, Num(T).mul(xr, Num(T).conj(yc, true), false), false), Num(T).mul(Num(T).conj(alpha, true), Num(T).mul(yr, Num(T).conj(xc, true), false), false))
            else if (herm)
                Num(T).scale(Num(T).mul(xr, Num(T).conj(xc, true), false), alpha_f)
            else if (two)
                Num(T).mul(alpha, Num(T).add(Num(T).mul(xr, yc, false), Num(T).mul(yr, xc, false)), false)
            else
                Num(T).mul(alpha, Num(T).mul(xr, xc, false), false);
            want[p] = Num(T).add(want[p], v);
            if (herm) {
                if (r == c_) want[p].im = 0;
            }
        };
        const un: Uplo = if (lower) .l else .u;
        const nn: isize = @intCast(n);
        const lda: isize = @intCast(ld);
        switch (st) {
            .packed_ => if (two) {
                if (!herm)
                    try spr2(talloc, check.io, T, un, nn, a0, &xs, inc[0], &ys, inc[1], &a)
                else
                    try hpr2(talloc, check.io, T, un, nn, a0, &xs, inc[0], &ys, inc[1], &a);
            } else {
                if (!herm)
                    try spr(talloc, check.io, T, un, nn, a0, &xs, inc[0], &a)
                else
                    try hpr(talloc, check.io, T, un, nn, alpha_f, &xs, inc[0], &a);
            },
            .full => if (two) {
                if (!herm)
                    try syr2(talloc, check.io, T, un, nn, a0, &xs, inc[0], &ys, inc[1], &a, lda)
                else
                    try her2(talloc, check.io, T, un, nn, a0, &xs, inc[0], &ys, inc[1], &a, lda);
            } else {
                if (!herm)
                    try syr(talloc, check.io, T, un, nn, a0, &xs, inc[0], &a, lda)
                else
                    try her(talloc, check.io, T, un, nn, alpha_f, &xs, inc[0], &a, lda);
            },
            .band => unreachable,
        }
        check.expectCloseSlices(T, &want, &a, tol(T, 4)) catch |e| {
            check.say("rank{s} {s} herm={} n={d} inc={d},{d} lower={}: ", .{ if (two) "2" else "1", @tagName(st), herm, n, inc[0], inc[1], lower });
            return e;
        };
    };
}

/// `gbmv` against dense arithmetic over every `trans`.
fn checkGbmv(comptime T: type) !void {
    const shapes = [_][4]usize{ .{ 1, 1, 0, 0 }, .{ 4, 4, 1, 1 }, .{ 5, 3, 2, 0 }, .{ 3, 6, 0, 2 }, .{ 7, 5, 3, 2 }, .{ 2, 9, 4, 4 } };
    var seed: u64 = 5000;
    inline for (std.enums.values(Trans)) |tr| {
        if (comptime !isComplex(T) and (tr.conjA() or tr.conjX())) continue;
        for (shapes) |s| for (inc_pairs) |inc| for ([_]f64{ 1, 0, 0.5 }) |bre| {
            seed += 1;
            const m, const n, const kl, const ku = s;
            const ld: usize = kl + ku + 1 + @as(usize, @intCast(seed & 1));
            const lenx = if (tr.transposed()) m else n;
            const leny = if (tr.transposed()) n else m;
            var a: [max_mem]T = undefined;
            var xs: [max_mem]T = undefined;
            var ys: [max_mem]T = undefined;
            check.fill(T, &a, seed);
            check.fill(T, &xs, seed + 100);
            check.fill(T, &ys, seed + 200);
            const alpha = scalarOf(T, 0.75, 1.5);
            const beta = scalarOf(T, bre, if (bre == 0.5) -1 else 0);
            var want = ys;
            for (0..leny) |o| {
                var acc = zero(T);
                for (0..lenx) |p| {
                    // A(i, j) with i the row, j the column of the m by n matrix
                    const i = if (tr.transposed()) p else o;
                    const j = if (tr.transposed()) o else p;
                    if (i + ku < j or i > j + kl) continue;
                    var e = a[ku + i - j + j * ld];
                    if (tr.conjA()) e = Num(T).conj(e, true);
                    var xv = xs[vecAt(p, lenx, inc[0])];
                    if (tr.conjX()) xv = Num(T).conj(xv, true);
                    acc = Num(T).add(acc, Num(T).mul(e, xv, false));
                }
                const yp = vecAt(o, leny, inc[1]);
                want[yp] = Num(T).add(if (bre == 0) zero(T) else Num(T).mul(beta, ys[yp], false), Num(T).mul(alpha, acc, false));
            }
            if (bre == 0) for (0..leny) |o| {
                ys[vecAt(o, leny, inc[1])] = scalarOf(T, std.math.nan(f64), std.math.nan(f64));
            };
            try gbmv(talloc, check.io, T, tr, @intCast(m), @intCast(n), @intCast(kl), @intCast(ku), alpha, &a, @intCast(ld), &xs, inc[0], beta, &ys, inc[1]);
            check.expectCloseSlices(T, &want, &ys, tol(T, lenx)) catch |e| {
                check.say("gbmv {s} {s} m={d} n={d} kl={d} ku={d} inc={d},{d}: ", .{ @typeName(T), @tagName(tr), m, n, kl, ku, inc[0], inc[1] });
                return e;
            };
        };
    }
}

fn runAll(comptime list: anytype) !void {
    inline for (list) |T| {
        try checkGbmv(T);
        inline for (.{ Store.band, Store.packed_ }) |st| {
            inline for (comptime kindsOf(T)) |kind| {
                // `hpmv` takes no `HEMVREV` triangle
                if (st == .packed_ and kind == .hermitian_rev) continue;
                for ([_]bool{ false, true }) |lower| try checkMv(T, st, kind, lower);
            }
        }
        inline for (.{ Store.full, Store.packed_ }) |st| {
            inline for (if (comptime isComplex(T)) .{ false, true } else .{false}) |herm| for ([_]bool{ false, true }) |lower| try checkRank(T, st, herm, lower);
        }
    }
}

/// Runs the cases with `vt`'s kernels in the global vtable.
fn withCore(vt: *const Vtable) anyerror!void {
    const saved = blas.testUse(vt);
    defer blas.testRestore(saved);
    try runAll(all);
}

test "banded, packed and rank updates match dense arithmetic (portable kernels)" {
    try runAll(all);
}

test "f80 without slots: banded, packed and rank updates match dense arithmetic (portable kernels)" {
    if (comptime check.f80_unslotted.len == 0) return error.SkipZigTest;
    try runAll(check.f80_unslotted);
}

test "banded, packed and rank updates match dense arithmetic on every core" {
    try check.eachCore(withCore);
}

test "strided vectors beyond the stack buffer are copied from the allocator" {
    const n = 1200; // > stack_bytes / 8
    const k = 3;
    const ld = k + 1;
    const a = try talloc.alloc(f64, ld * n);
    defer talloc.free(a);
    const xs = try talloc.alloc(f64, 2 * n);
    defer talloc.free(xs);
    const ys = try talloc.alloc(f64, 3 * n);
    defer talloc.free(ys);
    check.fill(f64, a, 1);
    check.fill(f64, xs, 2);
    check.fill(f64, ys, 3);
    const want = try talloc.dupe(f64, ys);
    defer talloc.free(want);
    for (0..n) |i| {
        var acc: f64 = 0;
        for (0..n) |j| {
            const d = if (j > i) j - i else i - j;
            if (d > k) continue;
            acc += a[k - d + @max(i, j) * ld] * xs[vecAt(j, n, -2)];
        }
        want[vecAt(i, n, 3)] += 2 * acc;
    }
    try sbmv(talloc, check.io, f64, .u, n, k, 2, a.ptr, ld, xs.ptr, -2, 1, ys.ptr, 3);
    try check.expectCloseSlices(f64, want, ys, tol(f64, k));
}

test "an allocation failure comes back as OutOfMemory" {
    // rstd's own path, which allocates
    _ = blas.vtable();
    const was = accel.routing();
    _ = accel.route(false);
    defer _ = accel.route(was);
    var failing = std.testing.FailingAllocator.init(talloc, .{ .fail_index = 0 });
    const n = 1200;
    const a = try talloc.alloc(f32, n);
    defer talloc.free(a);
    const xs = try talloc.alloc(f32, 2 * n);
    defer talloc.free(xs);
    @memset(a, 1);
    @memset(xs, 1);
    try std.testing.expectError(error.OutOfMemory, syr(failing.allocator(), check.io, f32, .l, n, 1, xs.ptr, 2, a.ptr, n));
}

test "invalid arguments and quick returns" {
    var a = [_]f64{ 1, 2, 3, 4 };
    var x = [_]f64{ 1, 2 };
    var y = [_]f64{ 5, 6 };
    try std.testing.expectError(error.InvalidLda, sbmv(talloc, check.io, f64, .u, 2, 1, 1, &a, 1, &x, 1, 1, &y, 1));
    try std.testing.expectError(error.InvalidIncx, sbmv(talloc, check.io, f64, .u, 2, 0, 1, &a, 1, &x, 0, 1, &y, 1));
    try std.testing.expectError(error.InvalidN, spmv(talloc, check.io, f64, .u, -1, 1, &a, &x, 1, 1, &y, 1));
    try std.testing.expectError(error.InvalidLda, syr(talloc, check.io, f64, .l, 2, 1, &x, 1, &a, 1));
    try std.testing.expectError(error.InvalidLda, gbmv(talloc, check.io, f64, .n, 2, 2, 1, 1, 1, &a, 2, &x, 1, 1, &y, 1));
    // n == 0 and alpha == 0 leave A alone; beta still scales y
    try syr(talloc, check.io, f64, .u, 0, 1, &x, 1, &a, 1);
    try syr(talloc, check.io, f64, .u, 2, 0, &x, 1, &a, 2);
    try std.testing.expectEqualSlices(f64, &.{ 1, 2, 3, 4 }, &a);
    try spmv(talloc, check.io, f64, .u, 2, 0, &a, &x, 1, 2, &y, 1);
    try std.testing.expectEqualSlices(f64, &.{ 10, 12 }, &y);
    // the symmetric routines take no Hermitian `uplo`, real `gbmv` no `o` to `d`
    try std.testing.expectError(error.InvalidUplo, sbmv(talloc, check.io, f64, .v, -1, -1, 1, &a, 0, &x, 0, 1, &y, 0));
    try std.testing.expectError(error.InvalidUplo, syr2(talloc, check.io, f64, .m, -1, 1, &x, 0, &y, 0, &a, 0));
    // only `hbmv` of these takes `v` and `m`; upstream's others get them from
    // row-major CBLAS alone
    var ca: [4]Complex64 = @splat(.init(1, 0));
    var cx: [2]Complex64 = @splat(.init(1, 0));
    var cy = cx;
    const c1: Complex64 = .init(1, 0);
    for ([_]Uplo{ .v, .m }) |u| {
        try std.testing.expectError(error.InvalidUplo, hpmv(talloc, check.io, Complex64, u, -1, c1, &ca, &cx, 0, c1, &cy, 0));
        try std.testing.expectError(error.InvalidUplo, her(talloc, check.io, Complex64, u, -1, 1, &cx, 0, &ca, 0));
        try std.testing.expectError(error.InvalidUplo, hpr(talloc, check.io, Complex64, u, -1, 1, &cx, 0, &ca));
        try std.testing.expectError(error.InvalidUplo, her2(talloc, check.io, Complex64, u, -1, c1, &cx, 0, &cy, 0, &ca, 0));
        try std.testing.expectError(error.InvalidUplo, hpr2(talloc, check.io, Complex64, u, -1, c1, &cx, 0, &cy, 0, &ca));
        try hbmv(talloc, check.io, Complex64, u, 2, 1, c1, &ca, 2, &cx, 1, c1, &cy, 1);
    }
    try std.testing.expectError(error.InvalidTrans, gbmv(talloc, check.io, f64, .d, -1, 2, 1, 1, 1, &a, 2, &x, 1, 1, &y, 1));
    // real `r` and `c` are `n` and `t`
    var yn = [_]f64{ 0, 0 };
    var yr = yn;
    try gbmv(talloc, check.io, f64, .t, 2, 2, 1, 0, 1, &a, 2, &x, 1, 0, &yn, 1);
    try gbmv(talloc, check.io, f64, .c, 2, 2, 1, 0, 1, &a, 2, &x, 1, 0, &yr, 1);
    try std.testing.expectEqualSlices(f64, &yn, &yr);
}
