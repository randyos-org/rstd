//! Tests of `../../blas3/half.zig`.

const std = @import("std");
const l2 = @import("../../blas2/core.zig");
const driver = @import("../../blas3/driver.zig");
const test_io = @import("../io.zig");
const types = @import("../../types.zig");
const portable = @import("../../portable.zig");
const util = @import("../../util.zig");
const Vtable = @import("../../Vtable.zig");
const bf = @import("../../ext/bf16/__root__.zig");
const bf16_kernel = @import("../../ext/bf16/kernel.zig");
const hroot = @import("../../blas3/half.zig");

const Allocator = std.mem.Allocator;
const Io = std.Io;
const BlasLong = types.BlasLong;
const Blocking = driver.Blocking;
const Trans = types.Trans;
const Uplo = types.Uplo;
const HalfPrecision = hroot.HalfPrecision;
const blocking = hroot.blocking;
const gemm = hroot.gemm;
const gemmHalf = hroot.gemmHalf;
const gemmtHalf = hroot.gemmtHalf;

// ---- tests -----------------------------------------------------------------------

const testing = std.testing;
const check = @import("../check.zig");
const params = @import("../../blas3/gemm/params.zig");
const zigstd = @import("../../../../../zigstd.zig");

const kinds = std.enums.values(HalfPrecision);

fn halfOf(comptime E: type, v: f32) E {
    return bf16_kernel.narrow(E, v);
}

fn valOf(comptime E: type, v: E) f32 {
    return bf16_kernel.wide(E, v);
}

/// Operands of one test problem, `A` and `B` stored with leading dimensions
/// past their rows, `C` with NaN past its rows (never read nor written).
fn Case(comptime kind: HalfPrecision) type {
    const I = kind.In();
    const O = kind.Out();
    return struct {
        const Self = @This();
        a: []I,
        b: []I,
        c: []O,
        lda: usize,
        ldb: usize,
        ldc: usize,

        fn init(al: Allocator, ta: bool, tb: bool, m: usize, n: usize, k: usize, seed: u64) !Self {
            const lda = (if (ta) k else m) + 2;
            const ldb = (if (tb) n else k) + 1;
            const ldc = m + 3;
            const a = try al.alloc(I, lda * @max(1, if (ta) m else k));
            errdefer al.free(a);
            const b = try al.alloc(I, ldb * @max(1, if (tb) k else n));
            errdefer al.free(b);
            const c = try al.alloc(O, ldc * n);
            var prng = std.Random.DefaultPrng.init(seed);
            const r = prng.random();
            for (a) |*v| v.* = halfOf(I, r.float(f32) * 2 - 1);
            for (b) |*v| v.* = halfOf(I, r.float(f32) * 2 - 1);
            for (c, 0..) |*v, i| v.* = halfOf(O, if (i % ldc < m) r.float(f32) * 2 - 1 else std.math.nan(f32));
            return .{ .a = a, .b = b, .c = c, .lda = lda, .ldb = ldb, .ldc = ldc };
        }

        fn deinit(s: Self, al: Allocator) void {
            al.free(s.a);
            al.free(s.b);
            al.free(s.c);
        }

        fn opA(s: Self, ta: bool, i: usize, l: usize) f32 {
            return valOf(I, if (ta) s.a[l + i * s.lda] else s.a[i + l * s.lda]);
        }

        fn opB(s: Self, tb: bool, l: usize, j: usize) f32 {
            return valOf(I, if (tb) s.b[j + l * s.ldb] else s.b[l + j * s.ldb]);
        }
    };
}

/// The kernels' arithmetic for one element, slice by slice: `beta` first
/// (as the beta kernel), then per depth slice the `f32` sum in order and
/// `C + alpha*sum`, rounded into `O`. The products of half values are exact
/// in `f32`, so only the sums and the update round.
fn model(comptime kind: HalfPrecision, blk: Blocking, slice: ?usize, cs: Case(kind), ta: bool, tb: bool, k: usize, alpha: f32, beta: f32, i: usize, j: usize, c0: kind.Out()) kind.Out() {
    const O = kind.Out();
    const one = struct {
        fn madd(x: f32, y: f32, z: f32) f32 {
            return util.rounded(f32, util.rounded(f32, x * y) + z);
        }
    };
    if (k == 0 or alpha == 0) return if (beta == 1) c0 else if (beta == 0) halfOf(O, 0) else halfOf(O, beta * valOf(O, c0));
    var cv: O = if (beta == 1) c0 else if (beta == 0) halfOf(O, 0) else halfOf(O, beta * valOf(O, c0));
    var ls: usize = 0;
    while (ls < k) {
        const min_l = if (slice) |q| @min(q, k - ls) else Blocking.split(k - ls, blk.q, blk.um);
        var acc: f32 = 0;
        for (ls..ls + min_l) |l| acc = one.madd(cs.opA(ta, i, l), cs.opB(tb, l, j), acc);
        cv = halfOf(O, one.madd(acc, alpha, valOf(O, cv)));
        ls += min_l;
    }
    return cv;
}

fn sameBits(comptime T: type, x: T, y: T) bool {
    const U = @Int(.unsigned, @bitSizeOf(T));
    const bx: U = if (T == types.BFloat16) @backingInt(x) else @bitCast(x);
    const by: U = if (T == types.BFloat16) @backingInt(y) else @bitCast(y);
    return bx == by;
}

/// `kind` on `vt` at `want` workers against `model`, bit for bit, and the
/// rows past `m` untouched; `slice`: the small path's depth slices, where
/// it takes the problem.
fn exact(comptime kind: HalfPrecision, vt: *const Vtable, io: Io, slice: ?usize, ta: bool, tb: bool, m: usize, n: usize, k: usize, alpha: f32, beta: f32, want: u32) !void {
    const O = kind.Out();
    const al = testing.allocator;
    const cs = try Case(kind).init(al, ta, tb, m, n, k, m * 131 + n * 17 + k);
    defer cs.deinit(al);
    const orig = try al.dupe(O, cs.c);
    defer al.free(orig);
    const af = valOf(O, halfOf(O, alpha));
    const bt = valOf(O, halfOf(O, beta));
    try gemm(al, io, kind, vt, ta, tb, m, n, k, halfOf(O, alpha), cs.a.ptr, cs.lda, cs.b.ptr, cs.ldb, halfOf(O, beta), cs.c.ptr, cs.ldc, want);
    const blk = blocking(kind, vt);
    for (0..n) |j| for (0..cs.ldc) |i| {
        const at = i + j * cs.ldc;
        const want_v = if (i < m) model(kind, blk, slice, cs, ta, tb, k, af, bt, i, j, orig[at]) else orig[at];
        if (!sameBits(O, want_v, cs.c[at])) {
            check.say("{s} ta={} tb={} m={d} n={d} k={d} alpha={d} beta={d} want={d} C[{d},{d}]: got {d}, want {d}\n", .{ @tagName(kind), ta, tb, m, n, k, alpha, beta, want, i, j, valOf(O, cs.c[at]), valOf(O, want_v) });
            return error.TestExpectedEqual;
        }
    };
}

/// A portable vtable with small blocks, so a few dozen rows and columns
/// cross every block edge.
fn smallBlocks() Vtable {
    var vt: Vtable = .{};
    portable.install(&vt);
    inline for (kinds) |kind| {
        const g = comptime kind.gemm();
        portable.paramPtr(&vt, g, .p).* = 2 * portable.paramPtr(&vt, g, .unroll_m).*;
        portable.paramPtr(&vt, g, .q).* = 5;
        portable.paramPtr(&vt, g, .r).* = 3 * portable.paramPtr(&vt, g, .unroll_n).* + 1;
    }
    return vt;
}

const shapes = [_][3]usize{ .{ 1, 1, 1 }, .{ 3, 2, 7 }, .{ 17, 9, 7 }, .{ 33, 20, 13 }, .{ 8, 40, 26 }, .{ 50, 3, 4 } };
const scalars = [_][2]f32{ .{ 1, 0 }, .{ 1.5, 2 }, .{ -0.5, 1 }, .{ 0, 0.5 }, .{ 2, -1 } };

test "half gemm: the portable kernels give the model's bits, every trans, shape and scalar case" {
    inline for (.{ false, true }) |small| {
        var vt: Vtable = if (small) smallBlocks() else .{};
        if (!small) portable.install(&vt);
        inline for (kinds) |kind| for (shapes) |s| for (scalars) |ab| for ([_]bool{ false, true }) |ta| for ([_]bool{ false, true }) |tb| {
            try exact(kind, &vt, check.io, null, ta, tb, s[0], s[1], s[2], ab[0], ab[1], 1);
        };
    }
}

fn alwaysHalf(_: i32, _: i32, _: BlasLong, _: BlasLong, _: BlasLong, _: f32, _: f32) callconv(.c) i32 {
    return 1;
}

test "half gemm: Cooper Lake's small-matrix sbgemm gives the portable bits in its depth slices" {
    const Cooperlake = @import("../../blas3/gemm/small_x86_64.zig").Cooperlake;
    var vt: Vtable = .{};
    portable.install(&vt);
    vt.sbgemm_small_matrix_permit = &alwaysHalf;
    inline for (.{ false, true }) |b0| inline for (0..2) |xa| inline for (0..2) |xb| {
        const s = comptime portable.GemmSlot.small(b0, xa, xb);
        portable.gemmSlotPtr(&vt, .sb, s).* = Cooperlake.entry(comptime portable.gemmSlotTag(.sb, s));
    };
    // every row block and column tail, and depths past one and two slices
    const deep = [_][3]usize{ .{ 70, 17, 5 }, .{ 15, 8, Cooperlake.slice + 3 }, .{ 33, 2, 2 * Cooperlake.slice + 1 } };
    for (shapes ++ deep) |s| for (scalars) |ab| for ([_]bool{ false, true }) |ta| for ([_]bool{ false, true }) |tb| {
        try exact(.bf16_to_f32, &vt, check.io, Cooperlake.slice, ta, tb, s[0], s[1], s[2], ab[0], ab[1], 1);
    };
    // Cooper Lake's own permit: 8, 16 or 32 columns up to `256^3`
    params.installCore(&vt, .x86_64, .cooperlake, .{});
    try exact(.bf16_to_f32, &vt, check.io, Cooperlake.slice, false, true, 21, 16, 9, 1.5, 2, 1);
    try exact(.bf16_to_f32, &vt, check.io, null, false, true, 21, 15, 9, 1.5, 2, 1);
}

test "half gemm: the same bits at 1 to 9 workers, on every Io" {
    var vt = smallBlocks();
    if (comptime test_io.hosted) {
        var t: Io.Threaded = .init(testing.allocator, .{ .async_limit = .limited(3), .concurrent_limit = .limited(3) });
        defer t.deinit();
        inline for (kinds) |kind| for ([_][3]usize{ .{ 41, 37, 12 }, .{ 70, 9, 11 }, .{ 5, 64, 6 } }) |s| for (1..10) |w| {
            try exact(kind, &vt, t.io(), null, true, false, s[0], s[1], s[2], 1.25, 0.5, @intCast(w));
            try exact(kind, &vt, t.io(), null, false, true, s[0], s[1], s[2], 1.25, 0.5, @intCast(w));
        };
    }
    inline for (kinds) |kind| for (1..5) |w| try exact(kind, &vt, Io.failing, null, false, false, 30, 29, 8, 1, 0, @intCast(w));
}

test "half gemm: a vector tile's lanes give the scalar bits, special values included" {
    // `widen` against `fromBf16` on every class: zeros, denormals, NaNs
    const bits = [_]u16{ 0x0000, 0x8000, 0x0001, 0x8075, 0x3f80, 0xc049, 0x7f80, 0xff80, 0x7f81, 0xffc1, 0x0080, 0x7f7f, 0x1234, 0xabcd, 0x0100, 0x8100 };
    var h: [bits.len]types.BFloat16 = undefined;
    for (&h, bits) |*d, s| d.* = @fromBackingInt(@intCast(s));
    inline for (.{ 2, 4, 8, 16 }) |n| {
        var at: usize = 0;
        while (at + n <= bits.len) : (at += n) {
            const v: [n]f32 = bf16_kernel.widen(types.BFloat16, n, h[at..].ptr);
            for (0..n) |l| try testing.expectEqual(@as(u32, @bitCast(bf.fromBf16(h[at + l]))), @as(u32, @bitCast(v[l])));
        }
    }
}

test "half gemm: alpha 0 or k 0 only scales C; beta 0 overwrites NaN; beta 1 leaves C" {
    const al = testing.allocator;
    var vt: Vtable = .{};
    portable.install(&vt);
    inline for (kinds) |kind| {
        const I = kind.In();
        const O = kind.Out();
        const nan_i = halfOf(I, std.math.nan(f32));
        var a = @as([4]I, @splat(nan_i));
        var c = [_]O{ halfOf(O, 1), halfOf(O, 2), halfOf(O, 3), halfOf(O, 4) };
        try gemm(al, Io.failing, kind, &vt, false, false, 2, 2, 2, halfOf(O, 0), &a, 2, &a, 2, halfOf(O, 2), &c, 2, 1);
        try testing.expectEqual(@as(f32, 8), valOf(O, c[3]));
        try gemm(al, Io.failing, kind, &vt, false, false, 2, 2, 0, halfOf(O, 1), &a, 2, &a, 2, halfOf(O, 1), &c, 2, 1);
        try testing.expectEqual(@as(f32, 8), valOf(O, c[3]));
        c = @splat(halfOf(O, std.math.nan(f32)));
        var one = @as([4]I, @splat(halfOf(I, 1)));
        try gemm(al, Io.failing, kind, &vt, true, true, 2, 2, 2, halfOf(O, 1), &one, 2, &one, 2, halfOf(O, 0), &c, 2, 1);
        for (c) |v| try testing.expectEqual(@as(f32, 2), valOf(O, v));
    }
}

test "interface: half gemm reports the first bad argument and does nothing for empty problems" {
    const fa = testing.failing_allocator;
    inline for (kinds) |kind| {
        const I = kind.In();
        const O = kind.Out();
        var a = @as([16]I, @splat(halfOf(I, 1)));
        var c = @as([16]O, @splat(halfOf(O, std.math.nan(f32))));
        const one = halfOf(O, 1);
        const zero = halfOf(O, 0);
        try testing.expectError(error.InvalidTrans, gemmHalf(fa, check.io, kind, .n, .o, -1, -1, -1, one, &a, 0, &a, 0, zero, &c, 0));
        try testing.expectError(error.InvalidM, gemmHalf(fa, check.io, kind, .n, .n, -1, -1, -1, one, &a, 0, &a, 0, zero, &c, 0));
        try testing.expectError(error.InvalidN, gemmHalf(fa, check.io, kind, .n, .n, 0, -1, -1, one, &a, 0, &a, 0, zero, &c, 0));
        try testing.expectError(error.InvalidK, gemmHalf(fa, check.io, kind, .n, .n, 0, 0, -1, one, &a, 0, &a, 0, zero, &c, 0));
        try testing.expectError(error.InvalidLda, gemmHalf(fa, check.io, kind, .n, .n, 3, 2, 4, one, &a, 2, &a, 1, zero, &c, 1));
        try testing.expectError(error.InvalidLda, gemmHalf(fa, check.io, kind, .t, .n, 3, 2, 4, one, &a, 3, &a, 1, zero, &c, 1));
        try testing.expectError(error.InvalidLdb, gemmHalf(fa, check.io, kind, .n, .n, 3, 2, 4, one, &a, 3, &a, 3, zero, &c, 1));
        try testing.expectError(error.InvalidLdb, gemmHalf(fa, check.io, kind, .n, .c, 3, 2, 4, one, &a, 3, &a, 1, zero, &c, 1));
        try testing.expectError(error.InvalidLdc, gemmHalf(fa, check.io, kind, .n, .n, 3, 2, 4, one, &a, 3, &a, 4, zero, &c, 2));
        const big = std.math.maxInt(isize) / 2;
        try testing.expectError(error.InvalidLda, gemmHalf(fa, check.io, kind, .n, .n, 2, 2, 3, one, &a, big, &a, 3, zero, &c, 2));
        // empty: nothing read or written, nothing allocated
        try gemmHalf(fa, check.io, kind, .n, .n, 0, 3, 2, one, &a, 1, &a, 2, zero, &c, 1);
        try gemmHalf(fa, check.io, kind, .n, .n, 3, 0, 2, one, &a, 3, &a, 2, zero, &c, 3);
        try testing.expect(std.math.isNan(valOf(O, c[0])));
        // `r` and `c` are `n` and `t`
        var x = [_]I{ halfOf(I, 1), halfOf(I, 2), halfOf(I, 3), halfOf(I, 4) };
        var y = @as([4]O, @splat(halfOf(O, 0)));
        try gemmHalf(fa, check.io, kind, .c, .r, 2, 2, 2, one, &x, 2, &x, 2, zero, &y, 2);
        // A^T B for A = B = [1 3; 2 4]
        try testing.expectEqual(@as(f32, 5), valOf(O, y[0]));
        try testing.expectEqual(@as(f32, 11), valOf(O, y[1]));
        try testing.expectEqual(@as(f32, 25), valOf(O, y[3]));
    }
}

test "interface: half gemm on a problem past the stack, at the global vtable" {
    const al = testing.allocator;
    inline for (kinds) |kind| {
        const O = kind.Out();
        const m = 70;
        const n = 45;
        const k = 300;
        const cs = try Case(kind).init(al, false, true, m, n, k, 7);
        defer cs.deinit(al);
        const orig = try al.dupe(O, cs.c);
        defer al.free(orig);
        try gemmHalf(al, check.io, kind, .n, .t, m, n, k, halfOf(O, 0.75), cs.a.ptr, @intCast(cs.lda), cs.b.ptr, @intCast(cs.ldb), halfOf(O, -1), cs.c.ptr, @intCast(cs.ldc));
        // `bgemm` rounds to bfloat16 per depth slice: 8 bits each, at
        // most `k` slices; `f32` out by the sums alone
        const rel: f32 = if (O == types.BFloat16) 0.004 * @as(f32, @floatFromInt(k)) else 2e-7 * @as(f32, @floatFromInt(k));
        for (0..n) |j| for (0..m) |i| {
            var want: f64 = -@as(f64, valOf(O, orig[i + j * cs.ldc]));
            var mag: f64 = @abs(want);
            for (0..k) |l| {
                const p = @as(f64, cs.opA(false, i, l)) * cs.opB(true, l, j) * 0.75;
                want += p;
                mag += @abs(p);
            }
            const got = valOf(O, cs.c[i + j * cs.ldc]);
            if (!(@abs(got - want) <= rel * mag)) {
                check.say("{s} C[{d},{d}]: got {d}, want {d}\n", .{ @tagName(kind), i, j, got, want });
                return error.TestExpectedApproxEqAbs;
            }
        };
    }
}

test "interface: gemmtHalf gives gemvHalf's bits per column, and leaves the other triangle" {
    const al = testing.allocator;
    inline for (kinds) |kind| {
        const I = kind.In();
        const O = kind.Out();
        for ([_]Uplo{ .u, .l }) |uplo| for ([_]Trans{ .n, .t }) |transa| for ([_]Trans{ .n, .c }) |transb| for ([_]usize{ 0, 1, 6 }) |k| for (scalars) |ab| {
            const m = 9;
            const ta = transa.transposed();
            const tb = transb.transposed();
            const cs = try Case(kind).init(al, ta, tb, m, m, k, k * 3 + m);
            defer cs.deinit(al);
            const want = try al.dupe(O, cs.c);
            defer al.free(want);
            const alpha = halfOf(O, ab[0]);
            const beta = halfOf(O, ab[1]);
            try gemmtHalf(check.io, kind, uplo, transa, transb, m, @intCast(k), alpha, cs.a.ptr, @intCast(cs.lda), cs.b.ptr, @intCast(cs.ldb), beta, cs.c.ptr, @intCast(cs.ldc));
            // the reference: one `gemvHalf` per column
            for (0..m) |i| {
                const r0 = if (uplo == .u) 0 else i;
                const len = if (uplo == .u) i + 1 else m - i;
                const bb: [*]I = if (tb) cs.b.ptr + i else cs.b.ptr + i * cs.ldb;
                const incb: isize = @intCast(if (tb) cs.ldb else 1);
                const col = want.ptr + r0 + i * cs.ldc;
                if (k == 0) {
                    if (ab[1] != 1) for (col[0..len]) |*v| {
                        v.* = if (ab[1] == 0) halfOf(O, 0) else halfOf(O, ab[1] * valOf(O, v.*));
                    };
                } else if (ta) {
                    try l2.gemvHalf(check.io, kind, .t, @intCast(k), @intCast(len), alpha, cs.a.ptr + r0 * cs.lda, @intCast(cs.lda), bb, incb, beta, col, 1);
                } else {
                    try l2.gemvHalf(check.io, kind, .n, @intCast(len), @intCast(k), alpha, cs.a.ptr + r0, @intCast(cs.lda), bb, incb, beta, col, 1);
                }
            }
            for (want, cs.c, 0..) |w, g, at| {
                if (!sameBits(O, w, g) and !(std.math.isNan(valOf(O, w)) and std.math.isNan(valOf(O, g)))) {
                    check.say("{s} {s} {s} {s} k={d} alpha={d} beta={d} C[{d}]: got {d}, want {d}\n", .{ @tagName(kind), @tagName(uplo), @tagName(transa), @tagName(transb), k, ab[0], ab[1], at, valOf(O, g), valOf(O, w) });
                    return error.TestExpectedEqual;
                }
            }
        };
    }
}

test "interface: gemmtHalf reports the first bad argument, and returns for empty problems" {
    inline for (kinds) |kind| {
        const I = kind.In();
        const O = kind.Out();
        var a = @as([16]I, @splat(halfOf(I, std.math.nan(f32))));
        var c = @as([16]O, @splat(halfOf(O, std.math.nan(f32))));
        const one = halfOf(O, 1);
        try testing.expectError(error.InvalidUplo, gemmtHalf(check.io, kind, .v, .o, .n, -1, -1, one, &a, 0, &a, 0, one, &c, 0));
        try testing.expectError(error.InvalidTrans, gemmtHalf(check.io, kind, .u, .o, .n, -1, -1, one, &a, 0, &a, 0, one, &c, 0));
        try testing.expectError(error.InvalidM, gemmtHalf(check.io, kind, .u, .n, .n, -1, -1, one, &a, 0, &a, 0, one, &c, 0));
        try testing.expectError(error.InvalidK, gemmtHalf(check.io, kind, .l, .n, .n, 0, -1, one, &a, 0, &a, 0, one, &c, 0));
        try testing.expectError(error.InvalidLda, gemmtHalf(check.io, kind, .l, .t, .n, 3, 4, one, &a, 3, &a, 4, one, &c, 3));
        try testing.expectError(error.InvalidLdb, gemmtHalf(check.io, kind, .u, .n, .c, 3, 4, one, &a, 3, &a, 2, one, &c, 3));
        try testing.expectError(error.InvalidLdc, gemmtHalf(check.io, kind, .u, .n, .n, 3, 4, one, &a, 3, &a, 4, one, &c, 2));
        try gemmtHalf(check.io, kind, .u, .n, .n, 0, 3, one, &a, 1, &a, 3, halfOf(O, 0), &c, 1);
        try testing.expect(std.math.isNan(valOf(O, c[0])));
        // a zero `alpha` scales the triangle without reading `A` or `B`
        try gemmtHalf(check.io, kind, .l, .n, .n, 2, 2, halfOf(O, 0), &a, 2, &a, 2, halfOf(O, 0), &c, 2);
        try testing.expectEqual(@as(f32, 0), valOf(O, c[0]));
        try testing.expectEqual(@as(f32, 0), valOf(O, c[1]));
        try testing.expect(std.math.isNan(valOf(O, c[2])));
        try testing.expectEqual(@as(f32, 0), valOf(O, c[3]));
    }
}

test "every core keeps the portable half gemm kernels and their unrolls" {
    const Each = struct {
        fn run(vt: *const Vtable) anyerror!void {
            inline for (bf16_kernel.families) |f| {
                const K = bf16_kernel.GemmKernels(f[1], f[2]);
                inline for (bf16_kernel.gemm_slots) |s| try zigstd.testExpectTrue(portable.gemmSlotPtr(vt, f[0], s).* == &@field(K, @tagName(s)));
                try zigstd.testExpectEqual(K.unroll_m, portable.paramPtr(vt, f[0], .unroll_m).*);
                try zigstd.testExpectEqual(K.unroll_n, portable.paramPtr(vt, f[0], .unroll_n).*);
            }
        }
    };
    try check.eachCore(Each.run);
}
