//! `blas2/interface.zig`'s tests at any float type and on every linked core.

const std = @import("std");
const blas = @import("../../__root__.zig");
const types = @import("../../types.zig");
const sym = @import("../../blas2/symmetric.zig");
const Vtable = @import("../../Vtable.zig");
const check = @import("../check.zig");
const reference = @import("../reference.zig");
const tri = @import("../../blas2/triangular.zig");
const band = @import("../../blas2/banded.zig");
const interface = @import("../../blas2/interface.zig");
const basic = @import("interface.zig");

const Trans = types.Trans;
const Uplo = types.Uplo;
const isComplex = types.isComplex;
const testing = std.testing;
const all = check.all_types;
const gemv = interface.gemv;
const HalfGemv = interface.HalfGemv;
const ger = interface.ger;
const gerc = interface.gerc;
const symv = interface.symv;
const hemv = interface.hemv;
const addT = basic.addT;
const zeroT = basic.zeroT;
const val = basic.val;
const testGemv = basic.testGemv;
const testGer = basic.testGer;
const symAt = basic.symAt;
const testMv = basic.testMv;
const testGemvHalf = basic.testGemvHalf;

/// `op(A)` element `(r, q)` of the column-major `a`, widened: `A`'s `(q, r)`
/// for `tr`, conjugated for `cj`.
fn opAt(comptime T: type, a: []const T, lda: usize, tr: bool, cj: bool, r: usize, q: usize) reference.Wide(T) {
    const e = reference.wide(T, if (tr) a[q + r * lda] else a[r + q * lda]);
    return if (cj) reference.conjW(e) else e;
}

/// `y := alpha*M*x + beta*y` against `got`, `M` given by `ctx.at(r, q)`; `x`
/// at stride 2, `y` and `got` at stride -1.
fn expectMv(comptime T: type, rows: usize, cols: usize, ctx: anytype, alpha: T, x: []const T, beta: T, y: []const T, got: []const T) !void {
    const wa = reference.wide(T, alpha);
    const wb = reference.wide(T, beta);
    for (0..rows) |r| {
        var acc = reference.wide(T, zeroT(T));
        var mag: f128 = 0;
        for (0..cols) |q| {
            const e = ctx.at(r, q);
            const xv = reference.wide(T, x[2 * q]);
            acc = reference.addW(acc, reference.mulW(e, xv));
            mag += reference.abs1(e) * reference.abs1(xv);
        }
        const yv = reference.wide(T, y[rows - 1 - r]);
        const want = reference.addW(reference.mulW(wa, acc), reference.mulW(wb, yv));
        try reference.expectNear(T, want, got[rows - 1 - r], @floatFromInt(2 * (cols + 2)), reference.abs1(wa) * mag + reference.abs1(wb) * reference.abs1(yv));
    }
}

/// Level 2 over `T`, which has no vtable slots, against an `f128` reference:
/// `gemv`, `ger`, `symv`, `hemv`, `trmv`, `trsv` and `gbmv`.
fn anyFloat(comptime T: type) !void {
    const W = reference.Wide(T);
    const al = testing.allocator;
    const complex = comptime isComplex(T);
    const m = 7;
    const n = 5;
    const lda = m + 1;
    var a: [lda * n]T = undefined;
    var x: [2 * m]T = undefined;
    var y: [m]T = undefined;
    check.fill(T, &a, 21);
    check.fill(T, &x, 22);
    check.fill(T, &y, 23);
    const alpha = val(T, 1.25, -0.5);
    const beta = val(T, -0.75, 0.25);

    const OpCtx = struct {
        a: []const T,
        tr: bool,
        cj: bool,
        fn at(c: @This(), r: usize, q: usize) W {
            return opAt(T, c.a, lda, c.tr, c.cj, r, q);
        }
    };
    const transes = if (complex) [_]Trans{ .n, .t, .c } else [_]Trans{ .n, .t };
    for (transes) |t| {
        const tr = t.transposed();
        var got = y;
        try gemv(al, check.io, T, t, m, n, alpha, &a, lda, &x, 2, beta, &got, -1);
        const ctx: OpCtx = .{ .a = &a, .tr = tr, .cj = t == .c };
        try expectMv(T, if (tr) n else m, if (tr) m else n, ctx, alpha, &x, beta, &y, &got);
    }

    const wa = reference.wide(T, alpha);
    for ([_]bool{ false, true }) |conj| {
        if (conj and !complex) continue;
        var got = a;
        if (!conj) {
            try ger(al, check.io, T, m, n, alpha, &x, 2, &y, -1, &got, lda);
        } else if (comptime complex) {
            try gerc(al, check.io, T, m, n, alpha, &x, 2, &y, -1, &got, lda);
        }
        for (0..n) |j| for (0..m) |i| {
            const xv = reference.wide(T, x[2 * i]);
            const yv = reference.wide(T, y[n - 1 - j]);
            const e = reference.wide(T, a[i + j * lda]);
            const want = reference.addW(e, reference.mulW(wa, reference.mulW(xv, if (conj) reference.conjW(yv) else yv)));
            const mag = reference.abs1(wa) * reference.abs1(xv) * reference.abs1(yv) + reference.abs1(e);
            try reference.expectNear(T, want, got[i + j * lda], 8, mag);
        };
    }

    // `symv` and `hemv` on the leading `n` by `n` block
    const SymCtx = struct {
        a: []const T,
        kind: sym.Kind,
        upper: bool,
        fn at(c: @This(), r: usize, q: usize) W {
            return reference.wide(T, symAt(T, c.kind, c.upper, c.a, lda, r, q));
        }
    };
    for ([_]bool{ false, true }) |herm| {
        if (herm and !complex) continue;
        for ([_]Uplo{ .u, .l }) |uplo| {
            var got = y;
            if (!herm) {
                try symv(al, check.io, T, uplo, n, alpha, &a, lda, &x, 2, beta, &got, -1);
            } else if (comptime complex) {
                try hemv(al, check.io, T, uplo, n, alpha, &a, lda, &x, 2, beta, &got, -1);
            }
            const ctx: SymCtx = .{ .a = &a, .kind = if (herm) .hermitian else .symmetric, .upper = uplo == .u };
            try expectMv(T, n, n, ctx, alpha, &x, beta, &y, &got);
        }
    }

    // `trmv`, and `trsv` on a matrix kept well away from singular
    var d = a;
    for (0..n) |i| d[i + i * lda] = addT(T, d[i + i * lda], val(T, 4, 0));
    const TriCtx = struct {
        a: []const T,
        upper: bool,
        tr: bool,
        cj: bool,
        fn at(c: @This(), r: usize, q: usize) W {
            // `A(i, j)`, zero outside its triangle
            const i = if (c.tr) q else r;
            const j = if (c.tr) r else q;
            if (if (c.upper) i > j else i < j) return reference.wide(T, zeroT(T));
            return opAt(T, c.a, lda, c.tr, c.cj, r, q);
        }
    };
    for ([_]Uplo{ .u, .l }) |uplo| for ([_]Trans{ .n, .t, .c }) |t| {
        const ctx: TriCtx = .{ .a = &d, .upper = uplo == .u, .tr = t != .n, .cj = t == .c and complex };
        var got = x;
        try tri.trmv(al, check.io, T, uplo, t, .n, n, &d, lda, &got, 2);
        // `op(A) x` is `expectMv` with `alpha` 1 and `beta` 0, read back at
        // stride -1
        var rev: [n]T = undefined;
        for (0..n) |i| rev[n - 1 - i] = got[2 * i];
        try expectMv(T, n, n, ctx, val(T, 1, 0), &x, zeroT(T), &rev, &rev);
        got = x;
        try tri.trsv(al, check.io, T, uplo, t, .n, n, &d, lda, &got, 2);
        // the residual of `op(A) z = x`, relative to `|op(A)| |z|`
        for (0..n) |r| {
            var acc = reference.wide(T, zeroT(T));
            var mag: f128 = 0;
            for (0..n) |q| {
                const z = reference.wide(T, got[2 * q]);
                acc = reference.addW(acc, reference.mulW(ctx.at(r, q), z));
                mag += reference.abs1(ctx.at(r, q)) * reference.abs1(z);
            }
            try reference.expectNear(T, acc, x[2 * r], @floatFromInt(4 * (n + 2)), mag);
        }
    };

    // `gbmv`, one sub- and two superdiagonals in band storage
    const kl = 1;
    const ku = 2;
    const ldb = kl + ku + 1;
    var ab: [ldb * n]T = undefined;
    check.fill(T, &ab, 24);
    const BandCtx = struct {
        ab: []const T,
        tr: bool,
        fn at(c: @This(), r: usize, q: usize) W {
            const i = if (c.tr) q else r;
            const j = if (c.tr) r else q;
            if (i + ku < j or i > j + kl) return reference.wide(T, zeroT(T));
            return reference.wide(T, c.ab[ku + i - j + j * ldb]);
        }
    };
    for ([_]Trans{ .n, .t }) |t| {
        const tr = t == .t;
        var got = y;
        try band.gbmv(al, check.io, T, t, m, n, kl, ku, alpha, &ab, ldb, &x, 2, beta, &got, -1);
        const ctx: BandCtx = .{ .ab = &ab, .tr = tr };
        try expectMv(T, if (tr) n else m, if (tr) m else n, ctx, alpha, &x, beta, &y, &got);
    }
}

test "interface: level 2 at any float type runs the portable kernels built for it" {
    inline for (check.any_floats) |F| {
        inline for (.{ F, std.math.Complex(F) }) |T| anyFloat(T) catch |e| {
            check.say("{s}\n", .{@typeName(T)});
            return e;
        };
    }
}

fn checkAll(vt: *const Vtable) !void {
    // the interface reaches kernels through the global vtable
    const saved = blas.testUse(vt);
    defer blas.testRestore(saved);
    inline for (std.enums.values(HalfGemv)) |kind| {
        try testGemvHalf(kind, .n, 9, 7, 2, 1, 1.25, 0.5);
        try testGemvHalf(kind, .t, 9, 7, -1, 3, 1.25, 0.5);
        try testGemvHalf(kind, .n, 33, 4, 1, 1, 1.25, 0);
        try testGemvHalf(kind, .t, 4, 33, 1, 1, 1.25, 1);
        try testGemvHalf(kind, .r, 17, 3, 1, -2, 0, 0.5);
    }
    inline for (all) |T| {
        try testGemv(T, .n, 9, 7, 2, 1, val(T, 1.25, -0.5), val(T, 0.5, 0.25));
        try testGemv(T, .t, 9, 7, -1, 3, val(T, 1.25, -0.5), val(T, 0.5, 0.25));
        try testGemv(T, .n, 33, 4, 1, 1, val(T, 1.25, -0.5), val(T, 0, 0));
        try testGemv(T, .t, 4, 33, 1, 1, val(T, 1.25, -0.5), val(T, 1, 0));
        if (comptime isComplex(T)) {
            inline for (.{ Trans.r, Trans.c, Trans.o, Trans.u, Trans.s, Trans.d }) |t| try testGemv(T, t, 9, 7, 2, -1, val(T, 1.25, -0.5), val(T, 0.5, 0.25));
        }
        try testGer(T, false, 17, 9, 1, 1, val(T, 1.25, -0.5), false);
        try testGer(T, false, 9, 17, -2, 3, val(T, 1.25, -0.5), true);
        if (comptime isComplex(T)) try testGer(T, true, 9, 17, 2, 1, val(T, 1.25, -0.5), true);
        for ([_]bool{ false, true }) |herm| {
            if (herm and !comptime isComplex(T)) continue;
            for (if (herm) std.enums.values(Uplo) else &[_]Uplo{ .u, .l }) |uplo| {
                for ([_]usize{ 1, 7, 8, 17, 33 }) |n| {
                    try testMv(T, herm, uplo, n, 1, 1, val(T, 1.25, -0.5), val(T, 0.5, 0.25));
                    try testMv(T, herm, uplo, n, -2, 3, val(T, 1.25, -0.5), val(T, 1, 0));
                }
            }
        }
    }
}

test "every linked core's kernels work through the level-2 interface, where this machine runs them" {
    try check.eachCore(checkAll);
}
