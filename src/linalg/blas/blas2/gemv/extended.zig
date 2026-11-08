//! The sums of OpenBLAS's x87 `q`/`x` gemv kernels (`kernel/x86_64/qgemv_n.S`,
//! `qgemv_t.S`, `xgemv_n.S`, `xgemv_t.S`; on x86, `kernel/x86/`'s, which differ
//! only in their blocks), which the portable gemv takes for floats wider than
//! `f64`. At `f80` they give the assembly's bits: the same blocks and
//! accumulators, filled and folded in its order, with nothing fused.

const std = @import("std");
const sysinfo = @import("builtin");
const types = @import("../../types.zig");
const portable = @import("../../portable.zig");
const util = @import("../../util.zig");
const Vtable = @import("../../Vtable.zig");
const gemm_params = @import("../../blas3/gemm/params.zig");

const BlasLong = types.BlasLong;
const XDouble = types.XDouble;

pub const X86CoreType = gemm_params.CoreOf(.x86);

const x86 = sysinfo.target.cpu.arch == .x86;

/// Rows a block of `qgemv_t.S` (`complex`: `xgemv_t.S`): 4096 on x86_64,
/// 1000 and 400 on x86.
pub fn extendedRowBlock(comptime complex: bool) BlasLong {
    return if (!x86) 4096 else if (complex) 400 else 1000;
}

/// Columns a block of `qgemv_n.S` (`complex`: `xgemv_n.S`): 32 on x86_64. On
/// x86, `core`'s `DTB_DEFAULT_ENTRIES`, halved for `xgemv_n` on `PENTIUM4` and
/// `ATHLON` (`qgemv_n`'s 32 on `ATHLON` and `OPTERON` is their DTB anyway).
pub fn extendedColBlock(comptime complex: bool, comptime core: X86CoreType) BlasLong {
    if (!x86) return 32;
    const dtb = gemm_params.dtbFor(.x86, core);
    const half = complex and (core == .northwood or core == .prescott or core == .athlon);
    return if (half) @divExact(dtb, 2) else dtb;
}

/// `qgemv_n.S` and `qgemv_t.S`: per block of `p` along the sum, each element
/// of `y` in its whole fours is one sequential sum; the one to three left take
/// four accumulators by index % 4, eight a trip, the tail into the first,
/// folded ((a0 + a1) + a2) + a3. Then `y += sum * alpha`, block by block.
pub fn extendedGemv(comptime F: type, comptime trans: bool, comptime p: BlasLong, m: BlasLong, n: BlasLong, alpha: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong) void {
    const outs = if (trans) n else m;
    const inner = if (trans) m else n;
    const seq = outs & -4;
    // `A`'s step along the sum and along `y`
    const step = if (trans) 1 else lda;
    const across = if (trans) lda else 1;
    var s: BlasLong = 0;
    while (s < inner) : (s += p) {
        const len = @min(p, inner - s);
        var o: BlasLong = 0;
        while (o < outs) : (o += 1) {
            const ao = s * step + o * across;
            var c: [4]F = @splat(0);
            var k: BlasLong = 0;
            if (o < seq) {
                while (k < len) : (k += 1) c[0] += util.at(F, a, ao + k * step).* * util.at(F, x, (s + k) * inc_x).*;
            } else {
                const main = len & -8;
                while (k < main) : (k += 8) {
                    inline for (0..8) |r| c[r % 4] += util.at(F, a, ao + (k + r) * step).* * util.at(F, x, (s + k + r) * inc_x).*;
                }
                while (k < len) : (k += 1) c[0] += util.at(F, a, ao + k * step).* * util.at(F, x, (s + k) * inc_x).*;
                c[0] = ((c[0] + c[1]) + c[2]) + c[3];
            }
            util.at(F, y, o * inc_y).* += c[0] * alpha;
        }
    }
}

/// `xgemv_n.S` and `xgemv_t.S`: per block of `p` along the sum, each element
/// of `y` sums `ar*xr`, `ai*xr` (subtracted for `conj`), `ar*xi` and `ai*xi`
/// apart, combines them as the variant's signs say, and adds `alpha` times it.
pub fn extendedZgemv(comptime F: type, comptime trans: bool, comptime conj: bool, comptime xconj: bool, comptime p: BlasLong, m: BlasLong, n: BlasLong, alpha_r: F, alpha_i: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong) void {
    const outs = if (trans) n else m;
    const inner = if (trans) m else n;
    const step = if (trans) 1 else lda;
    const across = if (trans) lda else 1;
    var s: BlasLong = 0;
    while (s < inner) : (s += p) {
        const len = @min(p, inner - s);
        var o: BlasLong = 0;
        while (o < outs) : (o += 1) {
            var c: [4]F = @splat(0);
            var k: BlasLong = 0;
            while (k < len) : (k += 1) {
                const e = 2 * (s * step + o * across + k * step);
                const ar = util.at(F, a, e).*;
                const ai = util.at(F, a, e + 1).*;
                const xr = util.at(F, x, 2 * (s + k) * inc_x).*;
                const xi = util.at(F, x, 2 * (s + k) * inc_x + 1).*;
                c[0] += ar * xr;
                if (conj) c[1] -= ai * xr else c[1] += ai * xr;
                c[2] += ar * xi;
                c[3] += ai * xi;
            }
            const re = if (conj == xconj) c[0] - c[3] else c[0] + c[3];
            const im = if (xconj) c[1] - c[2] else c[1] + c[2];
            util.at(F, y, 2 * o * inc_y).* += re * alpha_r - im * alpha_i;
            util.at(F, y, 2 * o * inc_y + 1).* += im * alpha_r + re * alpha_i;
        }
    }
}

/// The `n`-form kernels over `XDouble` with `p` columns a block.
fn ExtendedBlocked(comptime p: BlasLong) type {
    const F = XDouble;
    return struct {
        fn q(m: BlasLong, n: BlasLong, _: BlasLong, alpha: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, _: [*]F) callconv(types.kernelCc(F)) i32 {
            extendedGemv(F, false, p, m, n, alpha, a, lda, x, inc_x, y, inc_y);
            return 0;
        }
        fn Z(comptime conj: bool, comptime xconj: bool) type {
            return struct {
                fn k(m: BlasLong, n: BlasLong, _: BlasLong, alpha_r: F, alpha_i: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, _: [*]F) callconv(types.kernelCc(F)) i32 {
                    extendedZgemv(F, false, conj, xconj, p, m, n, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y);
                    return 0;
                }
            };
        }
    };
}

/// On x86, points the extended `n`-form slots at kernels with `core`'s column
/// block (`extendedColBlock`); the portable ones have the unknown core's.
pub fn installExtendedCore(vt: *Vtable, core: X86CoreType) void {
    if (comptime !x86 or !types.isExtended(XDouble)) return;
    switch (core) {
        inline else => |c| {
            vt.qgemv_n = &ExtendedBlocked(extendedColBlock(false, c)).q;
            const Z = ExtendedBlocked(extendedColBlock(true, c)).Z;
            inline for (.{ .{ .gemv_n, false, false }, .{ .gemv_r, true, false }, .{ .gemv_o, false, true }, .{ .gemv_s, true, true } }) |v| {
                portable.slotPtr(vt, .x, v[0]).* = &Z(v[1], v[2]).k;
            }
        },
    }
}

/// The portable kernels over `T`, any float or its `Complex`, by the slot

// ---- tests -----------------------------------------------------------------

const check = @import("../../_test/check.zig");
const zigstd = @import("../../../../zigstd.zig");
const x86cpu = @import("../../cpu/__root__.zig").x86;

/// The same sums written out plainly: by element of `y`, then by block.
const PlainExtended = struct {
    fn gemv(comptime F: type, trans: bool, p: usize, m: usize, n: usize, alpha: F, a: []const F, lda: usize, x: []const F, ix: usize, y: []F, iy: usize) void {
        const outs = if (trans) n else m;
        const inner = if (trans) m else n;
        for (0..outs) |o| {
            var b: usize = 0;
            while (b < inner) : (b += p) {
                const len = @min(p, inner - b);
                const seq = o < outs / 4 * 4;
                var acc: [4]F = @splat(0);
                for (0..len) |k| {
                    const e = if (trans) a[b + k + o * lda] else a[o + (b + k) * lda];
                    acc[if (seq or k >= len / 8 * 8) 0 else k % 4] += e * x[(b + k) * ix];
                }
                const sum = if (seq) acc[0] else ((acc[0] + acc[1]) + acc[2]) + acc[3];
                y[o * iy] += sum * alpha;
            }
        }
    }

    fn zgemv(comptime F: type, trans: bool, conj: bool, xconj: bool, p: usize, m: usize, n: usize, alpha: [2]F, a: []const F, lda: usize, x: []const F, ix: usize, y: []F, iy: usize) void {
        const outs = if (trans) n else m;
        const inner = if (trans) m else n;
        for (0..outs) |o| {
            var b: usize = 0;
            while (b < inner) : (b += p) {
                var c: [4]F = @splat(0);
                for (b..@min(b + p, inner)) |i| {
                    const e = 2 * (if (trans) i + o * lda else o + i * lda);
                    // accumulator 2q + r takes part r of `A` times part q of `x`
                    for (0..2) |q| for (0..2) |r| {
                        const t = a[e + r] * x[2 * i * ix + q];
                        if (conj and r == 1 and q == 0) c[2 * q + r] -= t else c[2 * q + r] += t;
                    };
                }
                const re = if (conj == xconj) c[0] - c[3] else c[0] + c[3];
                const im = if (xconj) c[1] - c[2] else c[1] + c[2];
                y[2 * o * iy] += re * alpha[0] - im * alpha[1];
                y[2 * o * iy + 1] += im * alpha[0] + re * alpha[1];
            }
        }
    }
};

fn expectExtendedBits(comptime F: type, want: []const F, got: []const F) !void {
    const U = @Int(.unsigned, @bitSizeOf(F));
    for (want, got, 0..) |w, g, i| if (@as(U, @bitCast(w)) != @as(U, @bitCast(g))) {
        check.say("y[{d}]: got {e}, want {e}\n", .{ i, g, w });
        return error.TestExpectedEqual;
    };
}

const extended_variants = .{
    .{ false, false }, .{ true, false }, .{ false, true }, .{ true, true },
};

/// One shape at one block size, real and every complex variant, against `PlainExtended`.
fn checkExtendedShape(comptime F: type, comptime trans: bool, comptime p: BlasLong, m: usize, n: usize, inc: [2]usize) !void {
    const al = std.testing.allocator;
    const lda = m + 1;
    const outs = if (trans) n else m;
    const inner = if (trans) m else n;
    const a = try al.alloc(F, 2 * lda * n);
    defer al.free(a);
    const x = try al.alloc(F, 2 * ((inner - 1) * inc[0] + 1));
    defer al.free(x);
    const y = try al.alloc(F, 2 * ((outs - 1) * inc[1] + 1));
    defer al.free(y);
    const want = try al.alloc(F, y.len);
    defer al.free(want);
    check.fill(F, a, m * 7 + n);
    check.fill(F, x, m + 3);
    const alpha: F = 1.25;
    const ml: BlasLong = @intCast(m);
    const nl: BlasLong = @intCast(n);
    const ll: BlasLong = @intCast(lda);
    const ix: BlasLong = @intCast(inc[0]);
    const iy: BlasLong = @intCast(inc[1]);
    errdefer check.say("{s} trans={} p={d} m={d} n={d} inc={d},{d}: ", .{ @typeName(F), trans, p, m, n, inc[0], inc[1] });

    check.fill(F, y, n + 5);
    @memcpy(want, y);
    PlainExtended.gemv(F, trans, @intCast(p), m, n, alpha, a, lda, x, inc[0], want, inc[1]);
    extendedGemv(F, trans, p, ml, nl, alpha, a.ptr, ll, x.ptr, ix, y.ptr, iy);
    try expectExtendedBits(F, want, y);

    inline for (extended_variants) |v| {
        errdefer check.say("conj={} xconj={} ", .{ v[0], v[1] });
        check.fill(F, y, n + 6);
        @memcpy(want, y);
        PlainExtended.zgemv(F, trans, v[0], v[1], @intCast(p), m, n, .{ alpha, -0.5 }, a, lda, x, inc[0], want, inc[1]);
        extendedZgemv(F, trans, v[0], v[1], p, ml, nl, alpha, -0.5, a.ptr, ll, x.ptr, ix, y.ptr, iy);
        try expectExtendedBits(F, want, y);
    }
}

test "extended gemv: the sums follow the assembly's order, bit for bit" {
    const shapes = [_][2]usize{
        .{ 1, 1 },   .{ 2, 3 },   .{ 3, 2 },   .{ 4, 4 },   .{ 5, 7 }, .{ 7, 5 }, .{ 8, 9 }, .{ 9, 8 },
        .{ 17, 13 }, .{ 13, 17 }, .{ 37, 70 }, .{ 70, 37 },
    };
    inline for (.{ f80, f128 }) |F| {
        for (shapes) |s| for ([_][2]usize{ .{ 1, 1 }, .{ 2, 3 } }) |inc| {
            inline for (.{ comptime extendedRowBlock(false), comptime extendedRowBlock(true), 20 }) |p| try checkExtendedShape(F, true, p, s[0], s[1], inc);
            inline for (.{ comptime extendedColBlock(false, .unknown), comptime extendedColBlock(true, .unknown), 20 }) |p| try checkExtendedShape(F, false, p, s[0], s[1], inc);
        };
        // past one block of rows, the transposed kernels' blocking
        try checkExtendedShape(F, true, extendedRowBlock(false), 4100, 6, .{ 2, 3 });
        try checkExtendedShape(F, true, extendedRowBlock(true), 4100, 6, .{ 1, 1 });
        try checkExtendedShape(F, false, extendedColBlock(false, .unknown), 6, 4100, .{ 2, 3 });
    }
}

test "extended gemv: the sums at f80 carry the full significand" {
    // 1 + 2^-63 squared is 1 + 2^-62 in f80 and 1 in f64
    const e = 1 + std.math.ldexp(@as(f80, 1), -63);
    const want = 1 + std.math.ldexp(@as(f80, 1), -62);
    var a = [_]f80{ e, 0 };
    var y = [_]f80{ 0, 0 };
    extendedGemv(f80, true, 4096, 1, 1, 1, &a, 1, &a, 1, &y, 1);
    try zigstd.testExpectEqual(want, y[0]);
    y[0] = 0;
    extendedZgemv(f80, false, false, false, 32, 1, 1, 1, 0, &a, 1, &a, 1, &y, 1);
    try zigstd.testExpectEqual(want, y[0]);
}

test "extended gemv: the sums give the values traced by hand from the .S files" {
    // 2^64 + 1 rounds to 2^64 in f80 (ties to even), so each value below
    // depends on which accumulator a term lands in, on the fold order and on
    // where blocks end. Traced once from kernel/x86_64/{q,x}gemv_{n,t}.S;
    // the plain loop of the portable gemv gives the value in parentheses.
    const F = f80;
    const t = std.math.ldexp(@as(F, 1), 64);
    var ones: [16]F = @splat(1);

    // qgemv_t: columns in fours, each one sum (`faddp` into st(1)..st(4)); the
    // n % 4 left go into four sums by row % 4, eight a trip (.L65), the tail
    // into the first (.L71), folded ((a0 + a1) + a2) + a3 (.L70). Each column
    // is t, 1, 0, 0, 0, 1, 0, 0, 1 and y starts at -t: one sum gives t, so 0;
    // column 4's a0 = t + 1 = t, a1 = 2 give t + 2, so 2 (0). A tail into a1
    // would make a1 = 3 and give t + 4.
    const col = [9]F{ t, 1, 0, 0, 0, 1, 0, 0, 1 };
    var at: [45]F = undefined;
    for (0..5) |j| @memcpy(at[9 * j ..][0..9], &col);
    var y: [5]F = @splat(-t);
    extendedGemv(F, true, extendedRowBlock(false), 9, 5, 1, &at, 9, &ones, 1, &y, 1);
    try zigstd.testExpectEqualSlices(F, &.{ 0, 0, 0, 0, 2 }, &y);

    // qgemv_n: rows in fours, each one sum over the columns; the m % 4 left as
    // qgemv_t's leftover columns. Each row is `col`: 0 for rows 0 to 3, 2 for
    // row 4. The portable `y += (alpha*x)*a` gives 3 for every row.
    var an: [45]F = undefined;
    for (0..9) |j| an[5 * j ..][0..5].* = @splat(col[j]);
    y = @splat(-t);
    extendedGemv(F, false, extendedColBlock(false, .unknown), 5, 9, 1, &an, 5, &ones, 1, &y, 1);
    try zigstd.testExpectEqualSlices(F, &.{ 0, 0, 0, 0, 2 }, &y);

    // xgemv_n and xgemv_t: four sums c1..c4 = ar*xr, ai*xr, ar*xi, ai*xi
    // (`CONJ`: c2 -= ai*xr, `fsubrp`), then re = c1 - c4, im = c2 + c3; `CONJ`
    // re = c1 + c4; `XCONJ` re = c1 + c4, im = c2 - c3; both, re = c1 - c4,
    // im = c2 - c3 (`fsubp` is fsubr in AT&T). A = (t, t), (1, 1) and
    // x = (1, 1), (1, 0): c1 = t + 1 = t, c2 = t (`CONJ`: -t - 1 = -t),
    // c3 = c4 = t. With alpha = 1, y = alpha * (re, im).
    var az = [_]F{ t, t, 1, 1 };
    var xz = [_]F{ 1, 1, 1, 0 };
    const want = [4][2]F{
        .{ 0, 2 * t }, // n, t (1, 2t)
        .{ 2 * t, 0 }, // r, c (2t, -1)
        .{ 2 * t, 0 }, // o, u (2t, 1)
        .{ 0, -2 * t }, // s, d (1, -2t)
    };
    inline for (extended_variants, want) |v, w| {
        var yz: [2]F = @splat(0);
        extendedZgemv(F, true, v[0], v[1], extendedRowBlock(true), 2, 1, 1, 0, &az, 2, &xz, 1, &yz, 1);
        try zigstd.testExpectEqualSlices(F, &w, &yz);
        yz = @splat(0);
        extendedZgemv(F, false, v[0], v[1], extendedColBlock(true, .unknown), 1, 2, 1, 0, &az, 1, &xz, 1, &yz, 1);
        try zigstd.testExpectEqualSlices(F, &w, &yz);
    }
}

test "extended gemv: the sums add each block to y as it ends" {
    // y starts at -t; t at the first index of the sum and 1 just past the
    // first block. ExtendedBlocked: -t + t, then + 1, so 1; one sum: t + 1 = t, so 0.
    const al = std.testing.allocator;
    const F = f80;
    const t = std.math.ldexp(@as(F, 1), 64);
    inline for (.{ false, true }) |complex| {
        const w = if (complex) 2 else 1;
        // the transposed kernels: rows in blocks, four columns
        {
            const p = comptime extendedRowBlock(complex);
            const m: usize = @intCast(p + 1);
            const a = try al.alloc(F, w * m * 4);
            defer al.free(a);
            const x = try al.alloc(F, w * m);
            defer al.free(x);
            @memset(a, 0);
            for (0..m) |i| {
                x[w * i] = 1;
                if (complex) x[w * i + 1] = 0;
            }
            for (0..4) |j| {
                a[w * j * m] = t;
                a[w * (j * m + m - 1)] = 1;
            }
            var y: [8]F = @splat(0);
            for (0..4) |j| y[w * j] = -t;
            if (complex) {
                extendedZgemv(F, true, false, false, p, @intCast(m), 4, 1, 0, a.ptr, @intCast(m), x.ptr, 1, &y, 1);
                try zigstd.testExpectEqualSlices(F, &.{ 1, 0, 1, 0, 1, 0, 1, 0 }, &y);
            } else {
                extendedGemv(F, true, p, @intCast(m), 4, 1, a.ptr, @intCast(m), x.ptr, 1, &y, 1);
                try zigstd.testExpectEqualSlices(F, &.{ 1, 1, 1, 1 }, y[0..4]);
            }
        }
        // the plain kernels: columns in blocks, four rows
        {
            const p = comptime extendedColBlock(complex, .unknown);
            const n: usize = @intCast(p + 1);
            const a = try al.alloc(F, w * 4 * n);
            defer al.free(a);
            const x = try al.alloc(F, w * n);
            defer al.free(x);
            @memset(a, 0);
            for (0..n) |j| {
                x[w * j] = 1;
                if (complex) x[w * j + 1] = 0;
            }
            for (0..4) |i| {
                a[w * i] = t;
                a[w * (i + 4 * (n - 1))] = 1;
            }
            var y: [8]F = @splat(0);
            for (0..4) |i| y[w * i] = -t;
            if (complex) {
                extendedZgemv(F, false, false, false, p, 4, @intCast(n), 1, 0, a.ptr, 4, x.ptr, 1, &y, 1);
                try zigstd.testExpectEqualSlices(F, &.{ 1, 0, 1, 0, 1, 0, 1, 0 }, &y);
            } else {
                extendedGemv(F, false, p, 4, @intCast(n), 1, a.ptr, 4, x.ptr, 1, &y, 1);
                try zigstd.testExpectEqualSlices(F, &.{ 1, 1, 1, 1 }, y[0..4]);
            }
        }
    }
}

test "extended gemv: the blocks by target and core" {
    if (x86) {
        try zigstd.testExpectEqual(@as(BlasLong, 1000), extendedRowBlock(false));
        try zigstd.testExpectEqual(@as(BlasLong, 400), extendedRowBlock(true));
        try zigstd.testExpectEqual(@as(BlasLong, 256), extendedColBlock(false, .core2));
        try zigstd.testExpectEqual(@as(BlasLong, 256), extendedColBlock(true, .core2));
        try zigstd.testExpectEqual(@as(BlasLong, 64), extendedColBlock(false, .northwood));
        try zigstd.testExpectEqual(@as(BlasLong, 32), extendedColBlock(true, .northwood));
        try zigstd.testExpectEqual(@as(BlasLong, 32), extendedColBlock(false, .athlon));
        try zigstd.testExpectEqual(@as(BlasLong, 16), extendedColBlock(true, .athlon));
        try zigstd.testExpectEqual(@as(BlasLong, 32), extendedColBlock(true, .opteron));
    } else {
        try zigstd.testExpectEqual(@as(BlasLong, 4096), extendedRowBlock(false));
        try zigstd.testExpectEqual(@as(BlasLong, 4096), extendedRowBlock(true));
        try zigstd.testExpectEqual(@as(BlasLong, 32), extendedColBlock(false, .core2));
        try zigstd.testExpectEqual(@as(BlasLong, 32), extendedColBlock(true, .athlon));
    }
}

/// x86's `qgemv_n.S` and `xgemv_n.S` columns a block (`P`) by core:
/// `getarch.c`'s `DTB_DEFAULT_ENTRIES`, which `xgemv_n.S` halves for
/// `PENTIUM4` (northwood, prescott) and `ATHLON`; `qgemv_n.S` takes 32 on
/// `ATHLON` and `OPTERON`, their own anyway.
const x87_col_blocks = [_]struct { core: X86CoreType, q: usize, x: usize }{
    .{ .core = .katmai, .q = 64, .x = 64 },
    .{ .core = .northwood, .q = 64, .x = 32 },
    .{ .core = .prescott, .q = 64, .x = 32 },
    .{ .core = .athlon, .q = 32, .x = 16 },
    .{ .core = .opteron, .q = 32, .x = 32 },
    .{ .core = .barcelona, .q = 48, .x = 48 },
    .{ .core = .bobcat, .q = 40, .x = 40 },
    .{ .core = .core2, .q = 256, .x = 256 },
    .{ .core = .nano, .q = 128, .x = 128 },
};

/// A core's vtable, static: a vtable is large for a rig's stack.
var x87_core_vt: Vtable = undefined;

/// `y := -t + A*x` on four rows, `t` = 2^64 in column 0 and 1 in column
/// `at`, the rest 0, `x` ones: a block ending before `at` gives 1, one
/// holding both gives `t + 1 = t`, so 0.
fn x87BlockProbe(comptime complex: bool, at: usize, n: usize) !XDouble {
    const al = std.testing.allocator;
    const F = XDouble;
    const w = if (complex) 2 else 1;
    const t = std.math.ldexp(@as(F, 1), 64);
    const a = try al.alloc(F, w * 4 * n);
    defer al.free(a);
    const x = try al.alloc(F, w * n);
    defer al.free(x);
    @memset(a, 0);
    @memset(x, 0);
    for (0..n) |j| x[w * j] = 1;
    for (0..4) |i| {
        a[w * i] = t;
        a[w * (i + 4 * at)] = 1;
    }
    var y: [8]F = @splat(0);
    for (0..4) |i| y[w * i] = -t;
    const ln: BlasLong = @intCast(n);
    if (complex) {
        _ = x87_core_vt.xgemv_n.?(4, ln, 0, 1, 0, a.ptr, 4, x.ptr, 1, &y, 1, undefined);
    } else {
        _ = x87_core_vt.qgemv_n.?(4, ln, 0, 1, a.ptr, 4, x.ptr, 1, &y, 1, undefined);
    }
    for (1..4) |i| try zigstd.testExpectEqual(y[0], y[w * i]);
    return y[0];
}

test "extended gemv: each x86 core's qgemv_n and xgemv_n slots take its column block" {
    if (comptime !x86 or !types.isExtended(XDouble)) return error.SkipZigTest;
    for (x87_col_blocks) |c| {
        x87_core_vt = x86cpu.vtableFromCoreType(c.core);
        inline for (.{ false, true }) |complex| {
            const p = if (complex) c.x else c.q;
            errdefer check.say("{s} {s}gemv_n, block {d}: ", .{ @tagName(c.core), if (complex) "x" else "q", p });
            // the 1 last in the first block, then first in the second
            try zigstd.testExpectEqual(@as(XDouble, 0), try x87BlockProbe(complex, p - 1, 2 * p + 1));
            try zigstd.testExpectEqual(@as(XDouble, 1), try x87BlockProbe(complex, p, 2 * p + 1));
        }
    }
}
