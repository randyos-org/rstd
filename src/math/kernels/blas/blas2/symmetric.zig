//! The portable symv and hemv kernels and their tests' reference, shared by
//! `symv/` and `hemv/`: `kernel/generic/symv_k.c`, `zsymv_k.c` and
//! `zhemv_k.c`, composing the core's own copy and gemv kernels in the same
//! order and blocking.

const std = @import("std");
const sysinfo = @import("builtin");
const blas = @import("../__root__.zig");
const portable = @import("../portable.zig");
const types = @import("../types.zig");
const util = @import("../util.zig");
const Vtable = @import("../Vtable.zig");
const gemv_kernels = @import("gemv/__root__.zig");
const check = @import("../_test/check.zig");
const guard = @import("../_test/guard.zig");
const copy_cases = @import("../_test/blas2/symcopy.zig").copy_cases;

/// The portable gemv slots, to tell a native one from them.
const portable_gemv: Vtable = blk: {
    var vt: Vtable = .{};
    gemv_kernels.installGeneric(&vt);
    break :blk vt;
};

const BlasLong = types.BlasLong;

/// How `S` comes from the stored triangle: symmetric (`S = S^T`), Hermitian
/// (`S = S^H`, the diagonal's imaginary part ignored), or the transpose of
/// the Hermitian one (`HEMVREV`: the `V` and `M` kernels).
pub const Kind = enum { symmetric, hermitian, hermitian_rev };

/// `SYMV_P` of the running core; the generic target's value until a core's
/// parameters are installed.
pub fn symvP() usize {
    const p = blas.table().symv_p;
    return if (p > 0) @intCast(p) else 16;
}

const page = 4096;

fn pageUp(addr: usize) usize {
    return (addr + page - 1) & ~@as(usize, page - 1);
}

/// Elements of the buffer a gemv kernel of `m` by `n` gets
/// (`interface/gemv.c`'s `buffer_size`; `cs` is 2 for complex).
pub fn gemvLen(comptime F: type, comptime cs: usize, m: usize, n: usize) usize {
    // the reference's pads are bytes; whole elements cover them at 12-byte f80
    const pad_bytes = 128 + if (sysinfo.target.os.tag == .windows) 160 else 0;
    const pad = (pad_bytes + @sizeOf(F) - 1) / @sizeOf(F);
    return (cs * (m + n) + pad + 3) & ~@as(usize, 3);
}

/// Bytes the generic symv and hemv need from a page-aligned buffer: the block
/// copy, then `x` and `y` copies of `m` on a page of their own where the
/// strides are not 1, then the gemv buffer.
pub fn composeBytes(comptime F: type, comptime cs: usize, m: usize, unit_x: bool, unit_y: bool) usize {
    const p = symvP();
    const es = cs * @sizeOf(F);
    var at = pageUp(p * p * es);
    if (!unit_y) at = pageUp(at + m * es);
    if (!unit_x) at = pageUp(at + m * es);
    return at + gemvLen(F, cs, m, p) * @sizeOf(F);
}

/// What the SSE `symv` and `hemv` assembly of x86_64 needs, whatever the
/// buffer's alignment: `x` copied, a gap of up to 512 bytes, then `y` copied
/// (`symv_L_sse2.S`, `zsymv_U_sse2.S` and the rest), and 4 elements of
/// slack past each. The 1280 bytes are ours: the 512 with room to spare,
/// not derived.
pub fn sseBytes(comptime F: type, comptime cs: usize, m: usize) usize {
    return 2 * cs * m * @sizeOf(F) + 1280;
}

/// The `symcopy.h` copies: the stored triangle of the `m` by `m` block at `a`
/// as the full block `b` (leading dimension `m`). `V` and `M` store the
/// Hermitian matrix's conjugate; its diagonal's imaginary part is 0.
fn copyBlock(comptime T: type, comptime lower: bool, comptime kind: Kind, m: usize, a: [*]const types.Float(T), lda: usize, b: [*]types.Float(T)) void {
    const complex = comptime types.isComplex(T);
    const cs = if (complex) 2 else 1;
    for (0..m) |j| for (0..m) |i| {
        const stored = if (lower) i >= j else i <= j;
        const src = a + cs * (if (stored) i + j * lda else j + i * lda);
        const dst = b + cs * (i + j * m);
        dst[0] = src[0];
        if (complex) {
            dst[1] = switch (kind) {
                .symmetric => src[1],
                .hermitian => if (i == j) 0 else if (stored) src[1] else -src[1],
                .hermitian_rev => if (i == j) 0 else if (stored) -src[1] else src[1],
            };
        }
    };
}

/// The same product as plain loops that call no kernel and need no buffer,
/// for the columns `offset` selects: the last `offset` of them for the upper triangle, the first
/// `offset` for the lower (`offset == m` is the whole product; the threaded
/// drivers split the columns).
pub fn Direct(comptime T: type, comptime lower: bool, comptime kind: Kind) type {
    const complex = comptime types.isComplex(T);
    const F = types.Float(T);
    // the row sums accumulate wide (`types.Acc`)
    const A = types.Acc(F);
    const W = if (complex) std.math.Complex(A) else A;
    return struct {
        fn wide(v: T) W {
            return if (complex) .init(v.re, v.im) else v;
        }
        fn narrow(v: W) T {
            return if (complex) .init(@floatCast(v.re), @floatCast(v.im)) else @floatCast(v);
        }
        fn get(a: [*]F, i: BlasLong) T {
            if (complex) return .init(util.at(F, a, 2 * i).*, util.at(F, a, 2 * i + 1).*);
            return util.at(F, a, i).*;
        }
        fn add(y: [*]F, i: BlasLong, v: T) void {
            if (complex) {
                util.at(F, y, 2 * i).* += v.re;
                util.at(F, y, 2 * i + 1).* += v.im;
            } else util.at(F, y, i).* += v;
        }
        fn mul(p: T, q: T) T {
            return if (complex) p.mul(q) else p * q;
        }
        /// The stored `a` as row `i`'s entry, and as its mirror's.
        fn pair(a: T) [2]T {
            return switch (kind) {
                .symmetric => .{ a, a },
                .hermitian => .{ a, a.conjugate() },
                .hermitian_rev => .{ a.conjugate(), a },
            };
        }
        fn diag(a: T) T {
            return if (kind == .symmetric) a else .init(a.re, 0);
        }
        fn run(m: BlasLong, offset: BlasLong, alpha: T, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong) void {
            var j: BlasLong = if (lower) 0 else m - offset;
            const end: BlasLong = if (lower) offset else m;
            while (j < end) : (j += 1) {
                const xj = mul(alpha, get(x, j * inc_x));
                var acc: W = if (complex) .init(0, 0) else 0;
                var i: BlasLong = if (lower) j + 1 else 0;
                const stop: BlasLong = if (lower) m else j;
                while (i < stop) : (i += 1) {
                    const e = pair(get(a, i + j * lda));
                    add(y, i * inc_y, mul(e[0], xj));
                    const p = wide(e[1]);
                    const q = wide(get(x, i * inc_x));
                    acc = if (complex) acc.add(p.mul(q)) else acc + p * q;
                }
                add(y, j * inc_y, narrow(if (complex) wide(alpha).mul(acc) else wide(alpha) * acc));
                add(y, j * inc_y, mul(diag(get(a, j + j * lda)), xj));
            }
        }
        pub fn real(m: BlasLong, offset: BlasLong, alpha: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, _: [*]F) callconv(types.kernelCc(F)) i32 {
            run(m, offset, alpha, a, lda, x, inc_x, y, inc_y);
            return 0;
        }
        pub fn cplx(m: BlasLong, offset: BlasLong, alpha_r: F, alpha_i: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, _: [*]F) callconv(types.kernelCc(F)) i32 {
            run(m, offset, .init(alpha_r, alpha_i), a, lda, x, inc_x, y, inc_y);
            return 0;
        }
    };
}

/// `symv_k.c`, `zsymv_k.c` and `zhemv_k.c` for the columns `offset` selects: the
/// last `offset` of them for the upper triangle, the first `offset` for the
/// lower (`offset == m` is the whole product; the threaded drivers split the
/// columns). `pre` is the precision letter, whose kernels it calls.
pub fn Compose(comptime pre: []const u8, comptime T: type, comptime lower: bool, comptime kind: Kind) type {
    const complex = comptime types.isComplex(T);
    const F = types.Float(T);
    const cs: usize = if (complex) 2 else 1;
    const es = cs * @sizeOf(F);
    const prec: portable.Prec = .from(pre);
    // the gemv that reads `A` transposed, and the one that adds `A x`
    const first: portable.Op = switch (kind) {
        .symmetric, .hermitian_rev => .gemv_t,
        .hermitian => .gemv_c,
    };
    const second: portable.Op = if (kind == .hermitian_rev) .gemv_r else .gemv_n;
    return struct {
        pub const first_gemv = first;
        pub const second_gemv = second;
        fn gemv(comptime op: portable.Op, m: BlasLong, n: BlasLong, ar: F, ai: F, a: [*]F, lda: BlasLong, x: [*]F, y: [*]F, buf: [*]F) void {
            const f = blas.kernel(comptime portable.slotTag(prec, op));
            if (comptime complex) {
                _ = f(m, n, 0, ar, ai, a, lda, x, 1, y, 1, buf);
            } else {
                _ = f(m, n, 0, ar, a, lda, x, 1, y, 1, buf);
            }
        }
        /// Element `i` of the vector at `p`.
        fn at(p: [*]F, i: BlasLong) [*]F {
            return p + cs * @as(usize, @intCast(i));
        }
        fn pageAt(addr: usize) [*]F {
            return @ptrFromInt(addr);
        }
        fn run(m: BlasLong, offset: BlasLong, ar: F, ai: F, a: [*]F, lda: BlasLong, x: [*]F, incx: BlasLong, y: [*]F, incy: BlasLong, buffer: [*]F) void {
            const copy = blas.kernel(portable.slotTag(.from(pre), .copy_k));
            const p: BlasLong = @intCast(symvP());
            const symbuffer = buffer;
            var gemvbuffer = pageUp(@intFromPtr(buffer) + @as(usize, @intCast(p * p)) * es);
            var buffer_x = gemvbuffer;
            var xs = x;
            var ys = y;
            if (incy != 1) {
                ys = pageAt(gemvbuffer);
                buffer_x = pageUp(gemvbuffer + @as(usize, @intCast(m)) * es);
                gemvbuffer = buffer_x;
                _ = copy(m, y, incy, ys, 1);
            }
            if (incx != 1) {
                xs = pageAt(buffer_x);
                gemvbuffer = pageUp(buffer_x + @as(usize, @intCast(m)) * es);
                _ = copy(m, x, incx, xs, 1);
            }
            const gb = pageAt(gemvbuffer);
            var is: BlasLong = if (lower) 0 else m - offset;
            const end: BlasLong = if (lower) offset else m;
            while (is < end) : (is += p) {
                const min_i = @min(end - is, p);
                const col = at(a, is * lda);
                if (!lower and is > 0) {
                    gemv(first, is, min_i, ar, ai, col, lda, xs, at(ys, is), gb);
                    gemv(second, is, min_i, ar, ai, col, lda, at(xs, is), ys, gb);
                }
                copyBlock(T, lower, kind, @intCast(min_i), at(a, is + is * lda), @intCast(lda), symbuffer);
                gemv(.gemv_n, min_i, min_i, ar, ai, symbuffer, min_i, at(xs, is), at(ys, is), gb);
                if (lower and m - is > min_i) {
                    const below = at(a, is + min_i + is * lda);
                    gemv(first, m - is - min_i, min_i, ar, ai, below, lda, at(xs, is + min_i), at(ys, is), gb);
                    gemv(second, m - is - min_i, min_i, ar, ai, below, lda, at(xs, is), at(ys, is + min_i), gb);
                }
            }
            if (incy != 1) _ = copy(m, ys, 1, y, incy);
        }
        pub fn real(m: BlasLong, offset: BlasLong, alpha: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, buffer: [*]F) callconv(.c) i32 {
            run(m, offset, alpha, 0, a, lda, x, inc_x, y, inc_y, buffer);
            return 0;
        }
        pub fn cplx(m: BlasLong, offset: BlasLong, alpha_r: F, alpha_i: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, buffer: [*]F) callconv(.c) i32 {
            run(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
            return 0;
        }
    };
}

/// Which kernels the product runs on: the composed form, the reference's,
/// with a native gemv; else the direct loops, 5 times faster over the
/// portable gemv (haswell, n to 1024, our timing for daa7a2d).
pub fn Auto(comptime pre: []const u8, comptime T: type, comptime lower: bool, comptime kind: Kind) type {
    const F = types.Float(T);
    const C = Compose(pre, T, lower, kind);
    const D = Direct(T, lower, kind);
    return struct {
        /// Whether the running core's gemv is its own, so the composed form
        /// runs; the direct loops need no buffer.
        pub fn native() bool {
            inline for ([_]portable.Op{ C.first_gemv, C.second_gemv, .gemv_n }) |op| {
                if (portable.slotPtr(blas.table(), .from(pre), op).* != portable.slotPtr(&portable_gemv, .from(pre), op).*) return true;
            }
            return false;
        }
        pub fn real(m: BlasLong, offset: BlasLong, alpha: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, buffer: [*]F) callconv(.c) i32 {
            return if (native()) C.real(m, offset, alpha, a, lda, x, inc_x, y, inc_y, buffer) else D.real(m, offset, alpha, a, lda, x, inc_x, y, inc_y, buffer);
        }
        pub fn cplx(m: BlasLong, offset: BlasLong, alpha_r: F, alpha_i: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, buffer: [*]F) callconv(.c) i32 {
            return if (native()) C.cplx(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer) else D.cplx(m, offset, alpha_r, alpha_i, a, lda, x, inc_x, y, inc_y, buffer);
        }
    };
}

// ---- the tests' reference --------------------------------------------------

/// `checkKernel` of the slot's kernel and of the composed one, which the slot
/// may not be on this core (`Auto`).
pub fn checkBoth(comptime pre: []const u8, comptime T: type, comptime lower: bool, comptime kind: Kind, kern: anytype, name: []const u8) !void {
    // x86_64's SSE kernels read past their copies, so their callers clear the
    // buffer; the composed kernel clears what it reads itself
    try checkKernel(T, lower, kind, kern, name, .junk);
    const C = Compose(pre, T, lower, kind);
    try checkKernel(T, lower, kind, if (comptime types.isComplex(T)) &C.cplx else &C.real, name, .nan);
}

/// Orders across the kernels' blocks, up to several blocks of each size.
const sizes = [_]usize{ 1, 2, 3, 4, 7, 8, 9, 16, 17, 33, 65 };
/// The `SYMV_P` of `param.h`'s cores.
const blocks = [_]i32{ 4, 8, 16 };

const strides = [_][2]usize{ .{ 1, 1 }, .{ 2, 3 } };

/// `kern` against plain arithmetic over the full matrix, for the whole
/// product and for part of it, split as the threaded driver
/// (`driver/level2/symv_thread.c`) splits it: on a multiple of 4 columns,
/// which the assembly kernels rely on. Under each block size, for the kernels
/// that compose others (`guard.use` puts the vtable they compose in place).
pub fn checkKernel(comptime T: type, comptime lower: bool, comptime kind: Kind, kern: anytype, name: []const u8, fill: guard.Fill) !void {
    const F = types.Float(T);
    const complex = comptime types.isComplex(T);
    const alpha: T = if (complex) .init(1.25, -0.5) else 1.25;
    const cs = if (complex) 2 else 1;
    for (blocks) |p| for (sizes) |m| for (strides) |inc| for ([_]usize{ m, if (m > 4) (if (lower) 4 else m - 4) else m }) |offset| {
        blas.testTable().symv_p = p;
        const lda = m + 1;
        const a = guard.operand(T, .a, lda * m);
        const x = guard.operand(T, .x, guard.span(m, inc[0]));
        const y = guard.operand(T, .y, guard.span(m, inc[1]));
        const want = guard.operand(T, .want, y.len);
        check.fill(T, a, m * 5 + offset);
        check.fill(T, x, m + 3);
        check.fill(T, y, m + 7);
        @memcpy(want, y);
        const lo = if (lower) 0 else m - offset;
        const hi = if (lower) offset else m;
        for (lo..hi) |j| for (0..m) |i| {
            const stored = if (lower) i >= j else i <= j;
            if (!stored) continue;
            const e = a[i + j * lda];
            // `S[i][j]` from the stored entry, and `S[j][i]` if not on the diagonal
            const sij: T, const sji: T = switch (kind) {
                .symmetric => .{ e, e },
                .hermitian => .{ e, if (complex) e.conjugate() else e },
                .hermitian_rev => .{ if (complex) e.conjugate() else e, e },
            };
            if (i == j) {
                const d: T = if (kind == .symmetric) e else if (complex) .init(e.re, 0) else e;
                want[j * inc[1]] = if (complex) want[j * inc[1]].add(alpha.mul(d).mul(x[j * inc[0]])) else want[j * inc[1]] + alpha * d * x[j * inc[0]];
            } else if (complex) {
                want[i * inc[1]] = want[i * inc[1]].add(alpha.mul(sij).mul(x[j * inc[0]]));
                want[j * inc[1]] = want[j * inc[1]].add(alpha.mul(sji).mul(x[i * inc[0]]));
            } else {
                want[i * inc[1]] += alpha * sij * x[j * inc[0]];
                want[j * inc[1]] += alpha * sji * x[i * inc[0]];
            }
        };
        const bytes = @max(composeBytes(F, cs, m, inc[0] == 1, inc[1] == 1), sseBytes(F, cs, m));
        const buf = guard.open(F, bytes, 0, fill);
        if (complex) {
            _ = kern(@intCast(m), @intCast(offset), alpha.re, alpha.im, @ptrCast(a.ptr), @intCast(lda), @ptrCast(x.ptr), @intCast(inc[0]), @ptrCast(y.ptr), @intCast(inc[1]), buf);
        } else {
            _ = kern(@intCast(m), @intCast(offset), alpha, a.ptr, @intCast(lda), x.ptr, @intCast(inc[0]), y.ptr, @intCast(inc[1]), buf);
        }
        check.expectCloseSlices(T, want, y, check.tolerance(F, 8 * m)) catch |e| {
            check.say("{s} p={d} m={d} offset={d}: ", .{ name, p, m, offset });
            return check.case(e, "symv", T, m, inc);
        };
        guard.intact(bytes, 0) catch |e| {
            check.say("{s} p={d} m={d} offset={d}: ", .{ name, p, m, offset });
            return check.case(e, "symv buffer", T, m, inc);
        };
    };
}

test "gemvLen pads whole elements covering the reference's bytes" {
    const win = sysinfo.target.os.tag == .windows;
    // f32 and f64 as the reference has them
    try std.testing.expectEqual((2 * 7 + 128 / 8 + (if (win) 160 / 8 else 0) + 3) & ~@as(usize, 3), gemvLen(f64, 2, 3, 4));
    try std.testing.expectEqual((5 + 128 / 4 + (if (win) 160 / 4 else 0) + 3) & ~@as(usize, 3), gemvLen(f32, 1, 2, 3));
    // a 12-byte element (f80 on x86) is not truncated
    inline for (.{ [3]u32, f80, f128 }) |F| {
        try std.testing.expect(gemvLen(F, 1, 0, 0) * @sizeOf(F) >= 128 + if (win) 160 else 0);
    }
}

test "block copies match symcopy.h" {
    const nan = std.math.nan(f64);
    for (copy_cases) |c| {
        const lda = c.m + 2;
        var a: [2 * 13 * 11]f64 = undefined;
        for (&a, 0..) |*v, i| v.* = @floatFromInt(1 + i);
        var b: [2 * 11 * 11]f64 = @splat(nan);
        switch (c.lower) {
            inline else => |lower| switch (c.kind) {
                inline else => |kind| if (c.complex) {
                    copyBlock(std.math.Complex(f64), lower, kind, c.m, &a, lda, &b);
                } else {
                    copyBlock(f64, lower, kind, c.m, &a, lda, &b);
                },
            },
        }
        try std.testing.expectEqualSlices(f64, c.want, b[0..c.want.len]);
        for (b[c.want.len..]) |v| try std.testing.expect(std.math.isNan(v));
    }
}
