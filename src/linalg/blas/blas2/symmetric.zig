//! The portable symv and hemv kernels and their tests' reference, shared by
//! `symv/` and `hemv/`: `y += alpha*S*x` for an `m` by `m` matrix `S` given by
//! one stored triangle of `A` (column-major).
//!
//! The kernels are `kernel/generic/symv_k.c`, `zsymv_k.c` and `zhemv_k.c`, which
//! compose the core's own copy and gemv kernels (`blas.kernel`, so they follow
//! the running core), in the same order and with the same blocking.

const std = @import("std");
const sysinfo = @import("builtin");
const blas = @import("../__root__.zig");
const types = @import("../types.zig");
const check = @import("../blas1/check.zig");
const guard = @import("guard.zig");
const util = @import("../util.zig");
const Vtable = @import("../Vtable.zig");
const gemv_kernels = @import("gemv/__root__.zig");

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
    const pad = if (sysinfo.os.tag == .windows) 160 / @sizeOf(F) else 0;
    return (cs * (m + n) + 128 / @sizeOf(F) + pad + 3) & ~@as(usize, 3);
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
/// slack past each.
pub fn sseBytes(comptime F: type, comptime cs: usize, m: usize) usize {
    return 2 * cs * m * @sizeOf(F) + 1280;
}

/// The `symcopy.h` copies: the stored triangle of the `m` by `m` block at `a`
/// as the full block `b` (leading dimension `m`). `V` and `M` store the
/// Hermitian matrix's conjugate; its diagonal's imaginary part is 0.
fn copyBlock(comptime T: type, comptime lower: bool, comptime kind: Kind, m: usize, a: [*]const check.FloatOf(T), lda: usize, b: [*]check.FloatOf(T)) void {
    const complex = comptime check.isComplex(T);
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
    const complex = comptime check.isComplex(T);
    const F = check.FloatOf(T);
    return struct {
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
                var acc: T = if (complex) .init(0, 0) else 0;
                var i: BlasLong = if (lower) j + 1 else 0;
                const stop: BlasLong = if (lower) m else j;
                while (i < stop) : (i += 1) {
                    const e = pair(get(a, i + j * lda));
                    add(y, i * inc_y, mul(e[0], xj));
                    acc = if (complex) acc.add(mul(e[1], get(x, i * inc_x))) else acc + e[1] * get(x, i * inc_x);
                }
                add(y, j * inc_y, mul(alpha, acc));
                add(y, j * inc_y, mul(diag(get(a, j + j * lda)), xj));
            }
        }
        pub fn real(m: BlasLong, offset: BlasLong, alpha: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, _: [*]F) callconv(.c) i32 {
            run(m, offset, alpha, a, lda, x, inc_x, y, inc_y);
            return 0;
        }
        pub fn cplx(m: BlasLong, offset: BlasLong, alpha_r: F, alpha_i: F, a: [*]F, lda: BlasLong, x: [*]F, inc_x: BlasLong, y: [*]F, inc_y: BlasLong, _: [*]F) callconv(.c) i32 {
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
    const complex = comptime check.isComplex(T);
    const F = check.FloatOf(T);
    const cs: usize = if (complex) 2 else 1;
    const es = cs * @sizeOf(F);
    // the gemv that reads `A` transposed, and the one that adds `A x`
    const first = pre ++ "gemv_" ++ switch (kind) {
        .symmetric, .hermitian_rev => "t",
        .hermitian => "c",
    };
    const second = pre ++ "gemv_" ++ if (kind == .hermitian_rev) "r" else "n";
    return struct {
        pub const first_gemv = first;
        pub const second_gemv = second;
        fn gemv(comptime name: []const u8, m: BlasLong, n: BlasLong, ar: F, ai: F, a: [*]F, lda: BlasLong, x: [*]F, y: [*]F, buf: [*]F) void {
            const f = blas.kernel(@field(blas.Slot, name));
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
            const copy = blas.kernel(@field(blas.Slot, pre ++ "copy_k"));
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
            // `zgemv_t.S` and its siblings return NaN if the buffer holds NaN
            // or infinity; cleared as `interface.gemv` does
            if (comptime complex and sysinfo.cpu.arch == .x86_64) @memset(gb[0..gemvLen(F, cs, @intCast(m), @intCast(p))], 0);
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
                gemv(pre ++ "gemv_n", min_i, min_i, ar, ai, symbuffer, min_i, at(xs, is), at(ys, is), gb);
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

/// Which kernels the product runs on. Measured on a core whose gemv is the
/// portable one (haswell, a Zig loop), composing is 5 times slower than the
/// direct loops (18 against 3.7 ns per element, n up to 1024), the block
/// copies and short gemvs costing more than they save; with a native gemv the
/// composed form is the reference's, so it is used exactly then.
pub fn Auto(comptime pre: []const u8, comptime T: type, comptime lower: bool, comptime kind: Kind) type {
    const F = check.FloatOf(T);
    const C = Compose(pre, T, lower, kind);
    const D = Direct(T, lower, kind);
    return struct {
        /// Whether the running core's gemv is its own, so the composed form
        /// runs; the direct loops need no buffer.
        pub fn native() bool {
            inline for (.{ C.first_gemv, C.second_gemv, pre ++ "gemv_n" }) |name| {
                if (@field(blas.table().*, name) != @field(portable_gemv, name)) return true;
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

/// `checkKernel` of the slot's kernel and of the composed one, which the slot
/// may not be on this core (`Auto`).
pub fn checkBoth(comptime pre: []const u8, comptime T: type, comptime lower: bool, comptime kind: Kind, kern: anytype, name: []const u8) !void {
    // x86_64's SSE kernels read past their copies, so their callers clear the
    // buffer; the composed kernel clears what it reads itself
    try checkKernel(T, lower, kind, kern, name, .junk);
    const C = Compose(pre, T, lower, kind);
    try checkKernel(T, lower, kind, if (comptime check.isComplex(T)) &C.cplx else &C.real, name, .nan);
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
    const F = check.FloatOf(T);
    const complex = comptime check.isComplex(T);
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

/// The block copies' results from `symcopy.h` itself (OpenBLAS v0.3.34-317, built
/// with `a[i] = 1 + i`, `lda = m + 2`), for odd and even `m`; re-derive them
/// from the header rather than editing.
const copy_cases = [_]struct { lower: bool, kind: Kind, complex: bool, m: usize, want: []const f64 }{
    .{ .lower = true, .kind = .symmetric, .complex = false, .m = 1, .want = &.{1} },
    .{ .lower = false, .kind = .symmetric, .complex = false, .m = 1, .want = &.{1} },
    .{ .lower = true, .kind = .symmetric, .complex = true, .m = 1, .want = &.{ 1, 2 } },
    .{ .lower = false, .kind = .symmetric, .complex = true, .m = 1, .want = &.{ 1, 2 } },
    .{ .lower = true, .kind = .hermitian, .complex = true, .m = 1, .want = &.{ 1, 0 } },
    .{ .lower = false, .kind = .hermitian, .complex = true, .m = 1, .want = &.{ 1, 0 } },
    .{ .lower = true, .kind = .hermitian_rev, .complex = true, .m = 1, .want = &.{ 1, 0 } },
    .{ .lower = false, .kind = .hermitian_rev, .complex = true, .m = 1, .want = &.{ 1, 0 } },
    .{ .lower = true, .kind = .symmetric, .complex = false, .m = 2, .want = &.{ 1, 2, 2, 6 } },
    .{ .lower = false, .kind = .symmetric, .complex = false, .m = 2, .want = &.{ 1, 5, 5, 6 } },
    .{ .lower = true, .kind = .symmetric, .complex = true, .m = 2, .want = &.{ 1, 2, 3, 4, 3, 4, 11, 12 } },
    .{ .lower = false, .kind = .symmetric, .complex = true, .m = 2, .want = &.{ 1, 2, 9, 10, 9, 10, 11, 12 } },
    .{ .lower = true, .kind = .hermitian, .complex = true, .m = 2, .want = &.{ 1, 0, 3, 4, 3, -4, 11, 0 } },
    .{ .lower = false, .kind = .hermitian, .complex = true, .m = 2, .want = &.{ 1, 0, 9, -10, 9, 10, 11, 0 } },
    .{ .lower = true, .kind = .hermitian_rev, .complex = true, .m = 2, .want = &.{ 1, 0, 3, -4, 3, 4, 11, 0 } },
    .{ .lower = false, .kind = .hermitian_rev, .complex = true, .m = 2, .want = &.{ 1, 0, 9, 10, 9, -10, 11, 0 } },
    .{ .lower = true, .kind = .symmetric, .complex = false, .m = 3, .want = &.{ 1, 2, 3, 2, 7, 8, 3, 8, 13 } },
    .{ .lower = false, .kind = .symmetric, .complex = false, .m = 3, .want = &.{ 1, 6, 11, 6, 7, 12, 11, 12, 13 } },
    .{ .lower = true, .kind = .symmetric, .complex = true, .m = 3, .want = &.{ 1, 2, 3, 4, 5, 6, 3, 4, 13, 14, 15, 16, 5, 6, 15, 16, 25, 26 } },
    .{ .lower = false, .kind = .symmetric, .complex = true, .m = 3, .want = &.{ 1, 2, 11, 12, 21, 22, 11, 12, 13, 14, 23, 24, 21, 22, 23, 24, 25, 26 } },
    .{ .lower = true, .kind = .hermitian, .complex = true, .m = 3, .want = &.{ 1, 0, 3, 4, 5, 6, 3, -4, 13, 0, 15, 16, 5, -6, 15, -16, 25, 0 } },
    .{ .lower = false, .kind = .hermitian, .complex = true, .m = 3, .want = &.{ 1, 0, 11, -12, 21, -22, 11, 12, 13, 0, 23, -24, 21, 22, 23, 24, 25, 0 } },
    .{ .lower = true, .kind = .hermitian_rev, .complex = true, .m = 3, .want = &.{ 1, 0, 3, -4, 5, -6, 3, 4, 13, 0, 15, -16, 5, 6, 15, 16, 25, 0 } },
    .{ .lower = false, .kind = .hermitian_rev, .complex = true, .m = 3, .want = &.{ 1, 0, 11, 12, 21, 22, 11, -12, 13, 0, 23, 24, 21, -22, 23, -24, 25, 0 } },
    .{ .lower = true, .kind = .symmetric, .complex = false, .m = 4, .want = &.{ 1, 2, 3, 4, 2, 8, 9, 10, 3, 9, 15, 16, 4, 10, 16, 22 } },
    .{ .lower = false, .kind = .symmetric, .complex = false, .m = 4, .want = &.{ 1, 7, 13, 19, 7, 8, 14, 20, 13, 14, 15, 21, 19, 20, 21, 22 } },
    .{ .lower = true, .kind = .symmetric, .complex = true, .m = 4, .want = &.{ 1, 2, 3, 4, 5, 6, 7, 8, 3, 4, 15, 16, 17, 18, 19, 20, 5, 6, 17, 18, 29, 30, 31, 32, 7, 8, 19, 20, 31, 32, 43, 44 } },
    .{ .lower = false, .kind = .symmetric, .complex = true, .m = 4, .want = &.{ 1, 2, 13, 14, 25, 26, 37, 38, 13, 14, 15, 16, 27, 28, 39, 40, 25, 26, 27, 28, 29, 30, 41, 42, 37, 38, 39, 40, 41, 42, 43, 44 } },
    .{ .lower = true, .kind = .hermitian, .complex = true, .m = 4, .want = &.{ 1, 0, 3, 4, 5, 6, 7, 8, 3, -4, 15, 0, 17, 18, 19, 20, 5, -6, 17, -18, 29, 0, 31, 32, 7, -8, 19, -20, 31, -32, 43, 0 } },
    .{ .lower = false, .kind = .hermitian, .complex = true, .m = 4, .want = &.{ 1, 0, 13, -14, 25, -26, 37, -38, 13, 14, 15, 0, 27, -28, 39, -40, 25, 26, 27, 28, 29, 0, 41, -42, 37, 38, 39, 40, 41, 42, 43, 0 } },
    .{ .lower = true, .kind = .hermitian_rev, .complex = true, .m = 4, .want = &.{ 1, 0, 3, -4, 5, -6, 7, -8, 3, 4, 15, 0, 17, -18, 19, -20, 5, 6, 17, 18, 29, 0, 31, -32, 7, 8, 19, 20, 31, 32, 43, 0 } },
    .{ .lower = false, .kind = .hermitian_rev, .complex = true, .m = 4, .want = &.{ 1, 0, 13, 14, 25, 26, 37, 38, 13, -14, 15, 0, 27, 28, 39, 40, 25, -26, 27, -28, 29, 0, 41, 42, 37, -38, 39, -40, 41, -42, 43, 0 } },
    .{ .lower = true, .kind = .symmetric, .complex = false, .m = 5, .want = &.{ 1, 2, 3, 4, 5, 2, 9, 10, 11, 12, 3, 10, 17, 18, 19, 4, 11, 18, 25, 26, 5, 12, 19, 26, 33 } },
    .{ .lower = false, .kind = .symmetric, .complex = false, .m = 5, .want = &.{ 1, 8, 15, 22, 29, 8, 9, 16, 23, 30, 15, 16, 17, 24, 31, 22, 23, 24, 25, 32, 29, 30, 31, 32, 33 } },
    .{ .lower = true, .kind = .symmetric, .complex = true, .m = 5, .want = &.{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 3, 4, 17, 18, 19, 20, 21, 22, 23, 24, 5, 6, 19, 20, 33, 34, 35, 36, 37, 38, 7, 8, 21, 22, 35, 36, 49, 50, 51, 52, 9, 10, 23, 24, 37, 38, 51, 52, 65, 66 } },
    .{ .lower = false, .kind = .symmetric, .complex = true, .m = 5, .want = &.{ 1, 2, 15, 16, 29, 30, 43, 44, 57, 58, 15, 16, 17, 18, 31, 32, 45, 46, 59, 60, 29, 30, 31, 32, 33, 34, 47, 48, 61, 62, 43, 44, 45, 46, 47, 48, 49, 50, 63, 64, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66 } },
    .{ .lower = true, .kind = .hermitian, .complex = true, .m = 5, .want = &.{ 1, 0, 3, 4, 5, 6, 7, 8, 9, 10, 3, -4, 17, 0, 19, 20, 21, 22, 23, 24, 5, -6, 19, -20, 33, 0, 35, 36, 37, 38, 7, -8, 21, -22, 35, -36, 49, 0, 51, 52, 9, -10, 23, -24, 37, -38, 51, -52, 65, 0 } },
    .{ .lower = false, .kind = .hermitian, .complex = true, .m = 5, .want = &.{ 1, 0, 15, -16, 29, -30, 43, -44, 57, -58, 15, 16, 17, 0, 31, -32, 45, -46, 59, -60, 29, 30, 31, 32, 33, 0, 47, -48, 61, -62, 43, 44, 45, 46, 47, 48, 49, 0, 63, -64, 57, 58, 59, 60, 61, 62, 63, 64, 65, 0 } },
    .{ .lower = true, .kind = .hermitian_rev, .complex = true, .m = 5, .want = &.{ 1, 0, 3, -4, 5, -6, 7, -8, 9, -10, 3, 4, 17, 0, 19, -20, 21, -22, 23, -24, 5, 6, 19, 20, 33, 0, 35, -36, 37, -38, 7, 8, 21, 22, 35, 36, 49, 0, 51, -52, 9, 10, 23, 24, 37, 38, 51, 52, 65, 0 } },
    .{ .lower = false, .kind = .hermitian_rev, .complex = true, .m = 5, .want = &.{ 1, 0, 15, 16, 29, 30, 43, 44, 57, 58, 15, -16, 17, 0, 31, 32, 45, 46, 59, 60, 29, -30, 31, -32, 33, 0, 47, 48, 61, 62, 43, -44, 45, -46, 47, -48, 49, 0, 63, 64, 57, -58, 59, -60, 61, -62, 63, -64, 65, 0 } },
    .{ .lower = true, .kind = .symmetric, .complex = false, .m = 6, .want = &.{ 1, 2, 3, 4, 5, 6, 2, 10, 11, 12, 13, 14, 3, 11, 19, 20, 21, 22, 4, 12, 20, 28, 29, 30, 5, 13, 21, 29, 37, 38, 6, 14, 22, 30, 38, 46 } },
    .{ .lower = false, .kind = .symmetric, .complex = false, .m = 6, .want = &.{ 1, 9, 17, 25, 33, 41, 9, 10, 18, 26, 34, 42, 17, 18, 19, 27, 35, 43, 25, 26, 27, 28, 36, 44, 33, 34, 35, 36, 37, 45, 41, 42, 43, 44, 45, 46 } },
    .{ .lower = true, .kind = .symmetric, .complex = true, .m = 6, .want = &.{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 3, 4, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 5, 6, 21, 22, 37, 38, 39, 40, 41, 42, 43, 44, 7, 8, 23, 24, 39, 40, 55, 56, 57, 58, 59, 60, 9, 10, 25, 26, 41, 42, 57, 58, 73, 74, 75, 76, 11, 12, 27, 28, 43, 44, 59, 60, 75, 76, 91, 92 } },
    .{ .lower = false, .kind = .symmetric, .complex = true, .m = 6, .want = &.{ 1, 2, 17, 18, 33, 34, 49, 50, 65, 66, 81, 82, 17, 18, 19, 20, 35, 36, 51, 52, 67, 68, 83, 84, 33, 34, 35, 36, 37, 38, 53, 54, 69, 70, 85, 86, 49, 50, 51, 52, 53, 54, 55, 56, 71, 72, 87, 88, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 89, 90, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92 } },
    .{ .lower = true, .kind = .hermitian, .complex = true, .m = 6, .want = &.{ 1, 0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 3, -4, 19, 0, 21, 22, 23, 24, 25, 26, 27, 28, 5, -6, 21, -22, 37, 0, 39, 40, 41, 42, 43, 44, 7, -8, 23, -24, 39, -40, 55, 0, 57, 58, 59, 60, 9, -10, 25, -26, 41, -42, 57, -58, 73, 0, 75, 76, 11, -12, 27, -28, 43, -44, 59, -60, 75, -76, 91, 0 } },
    .{ .lower = false, .kind = .hermitian, .complex = true, .m = 6, .want = &.{ 1, 0, 17, -18, 33, -34, 49, -50, 65, -66, 81, -82, 17, 18, 19, 0, 35, -36, 51, -52, 67, -68, 83, -84, 33, 34, 35, 36, 37, 0, 53, -54, 69, -70, 85, -86, 49, 50, 51, 52, 53, 54, 55, 0, 71, -72, 87, -88, 65, 66, 67, 68, 69, 70, 71, 72, 73, 0, 89, -90, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 0 } },
    .{ .lower = true, .kind = .hermitian_rev, .complex = true, .m = 6, .want = &.{ 1, 0, 3, -4, 5, -6, 7, -8, 9, -10, 11, -12, 3, 4, 19, 0, 21, -22, 23, -24, 25, -26, 27, -28, 5, 6, 21, 22, 37, 0, 39, -40, 41, -42, 43, -44, 7, 8, 23, 24, 39, 40, 55, 0, 57, -58, 59, -60, 9, 10, 25, 26, 41, 42, 57, 58, 73, 0, 75, -76, 11, 12, 27, 28, 43, 44, 59, 60, 75, 76, 91, 0 } },
    .{ .lower = false, .kind = .hermitian_rev, .complex = true, .m = 6, .want = &.{ 1, 0, 17, 18, 33, 34, 49, 50, 65, 66, 81, 82, 17, -18, 19, 0, 35, 36, 51, 52, 67, 68, 83, 84, 33, -34, 35, -36, 37, 0, 53, 54, 69, 70, 85, 86, 49, -50, 51, -52, 53, -54, 55, 0, 71, 72, 87, 88, 65, -66, 67, -68, 69, -70, 71, -72, 73, 0, 89, 90, 81, -82, 83, -84, 85, -86, 87, -88, 89, -90, 91, 0 } },
    .{ .lower = true, .kind = .symmetric, .complex = false, .m = 7, .want = &.{ 1, 2, 3, 4, 5, 6, 7, 2, 11, 12, 13, 14, 15, 16, 3, 12, 21, 22, 23, 24, 25, 4, 13, 22, 31, 32, 33, 34, 5, 14, 23, 32, 41, 42, 43, 6, 15, 24, 33, 42, 51, 52, 7, 16, 25, 34, 43, 52, 61 } },
    .{ .lower = false, .kind = .symmetric, .complex = false, .m = 7, .want = &.{ 1, 10, 19, 28, 37, 46, 55, 10, 11, 20, 29, 38, 47, 56, 19, 20, 21, 30, 39, 48, 57, 28, 29, 30, 31, 40, 49, 58, 37, 38, 39, 40, 41, 50, 59, 46, 47, 48, 49, 50, 51, 60, 55, 56, 57, 58, 59, 60, 61 } },
    .{ .lower = true, .kind = .symmetric, .complex = true, .m = 7, .want = &.{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 3, 4, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 5, 6, 23, 24, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 7, 8, 25, 26, 43, 44, 61, 62, 63, 64, 65, 66, 67, 68, 9, 10, 27, 28, 45, 46, 63, 64, 81, 82, 83, 84, 85, 86, 11, 12, 29, 30, 47, 48, 65, 66, 83, 84, 101, 102, 103, 104, 13, 14, 31, 32, 49, 50, 67, 68, 85, 86, 103, 104, 121, 122 } },
    .{ .lower = false, .kind = .symmetric, .complex = true, .m = 7, .want = &.{ 1, 2, 19, 20, 37, 38, 55, 56, 73, 74, 91, 92, 109, 110, 19, 20, 21, 22, 39, 40, 57, 58, 75, 76, 93, 94, 111, 112, 37, 38, 39, 40, 41, 42, 59, 60, 77, 78, 95, 96, 113, 114, 55, 56, 57, 58, 59, 60, 61, 62, 79, 80, 97, 98, 115, 116, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 99, 100, 117, 118, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 119, 120, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122 } },
    .{ .lower = true, .kind = .hermitian, .complex = true, .m = 7, .want = &.{ 1, 0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 3, -4, 21, 0, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 5, -6, 23, -24, 41, 0, 43, 44, 45, 46, 47, 48, 49, 50, 7, -8, 25, -26, 43, -44, 61, 0, 63, 64, 65, 66, 67, 68, 9, -10, 27, -28, 45, -46, 63, -64, 81, 0, 83, 84, 85, 86, 11, -12, 29, -30, 47, -48, 65, -66, 83, -84, 101, 0, 103, 104, 13, -14, 31, -32, 49, -50, 67, -68, 85, -86, 103, -104, 121, 0 } },
    .{ .lower = false, .kind = .hermitian, .complex = true, .m = 7, .want = &.{ 1, 0, 19, -20, 37, -38, 55, -56, 73, -74, 91, -92, 109, -110, 19, 20, 21, 0, 39, -40, 57, -58, 75, -76, 93, -94, 111, -112, 37, 38, 39, 40, 41, 0, 59, -60, 77, -78, 95, -96, 113, -114, 55, 56, 57, 58, 59, 60, 61, 0, 79, -80, 97, -98, 115, -116, 73, 74, 75, 76, 77, 78, 79, 80, 81, 0, 99, -100, 117, -118, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 0, 119, -120, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 0 } },
    .{ .lower = true, .kind = .hermitian_rev, .complex = true, .m = 7, .want = &.{ 1, 0, 3, -4, 5, -6, 7, -8, 9, -10, 11, -12, 13, -14, 3, 4, 21, 0, 23, -24, 25, -26, 27, -28, 29, -30, 31, -32, 5, 6, 23, 24, 41, 0, 43, -44, 45, -46, 47, -48, 49, -50, 7, 8, 25, 26, 43, 44, 61, 0, 63, -64, 65, -66, 67, -68, 9, 10, 27, 28, 45, 46, 63, 64, 81, 0, 83, -84, 85, -86, 11, 12, 29, 30, 47, 48, 65, 66, 83, 84, 101, 0, 103, -104, 13, 14, 31, 32, 49, 50, 67, 68, 85, 86, 103, 104, 121, 0 } },
    .{ .lower = false, .kind = .hermitian_rev, .complex = true, .m = 7, .want = &.{ 1, 0, 19, 20, 37, 38, 55, 56, 73, 74, 91, 92, 109, 110, 19, -20, 21, 0, 39, 40, 57, 58, 75, 76, 93, 94, 111, 112, 37, -38, 39, -40, 41, 0, 59, 60, 77, 78, 95, 96, 113, 114, 55, -56, 57, -58, 59, -60, 61, 0, 79, 80, 97, 98, 115, 116, 73, -74, 75, -76, 77, -78, 79, -80, 81, 0, 99, 100, 117, 118, 91, -92, 93, -94, 95, -96, 97, -98, 99, -100, 101, 0, 119, 120, 109, -110, 111, -112, 113, -114, 115, -116, 117, -118, 119, -120, 121, 0 } },
    .{ .lower = true, .kind = .symmetric, .complex = false, .m = 8, .want = &.{ 1, 2, 3, 4, 5, 6, 7, 8, 2, 12, 13, 14, 15, 16, 17, 18, 3, 13, 23, 24, 25, 26, 27, 28, 4, 14, 24, 34, 35, 36, 37, 38, 5, 15, 25, 35, 45, 46, 47, 48, 6, 16, 26, 36, 46, 56, 57, 58, 7, 17, 27, 37, 47, 57, 67, 68, 8, 18, 28, 38, 48, 58, 68, 78 } },
    .{ .lower = false, .kind = .symmetric, .complex = false, .m = 8, .want = &.{ 1, 11, 21, 31, 41, 51, 61, 71, 11, 12, 22, 32, 42, 52, 62, 72, 21, 22, 23, 33, 43, 53, 63, 73, 31, 32, 33, 34, 44, 54, 64, 74, 41, 42, 43, 44, 45, 55, 65, 75, 51, 52, 53, 54, 55, 56, 66, 76, 61, 62, 63, 64, 65, 66, 67, 77, 71, 72, 73, 74, 75, 76, 77, 78 } },
    .{ .lower = true, .kind = .symmetric, .complex = true, .m = 8, .want = &.{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 3, 4, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 5, 6, 25, 26, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 7, 8, 27, 28, 47, 48, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 9, 10, 29, 30, 49, 50, 69, 70, 89, 90, 91, 92, 93, 94, 95, 96, 11, 12, 31, 32, 51, 52, 71, 72, 91, 92, 111, 112, 113, 114, 115, 116, 13, 14, 33, 34, 53, 54, 73, 74, 93, 94, 113, 114, 133, 134, 135, 136, 15, 16, 35, 36, 55, 56, 75, 76, 95, 96, 115, 116, 135, 136, 155, 156 } },
    .{ .lower = false, .kind = .symmetric, .complex = true, .m = 8, .want = &.{ 1, 2, 21, 22, 41, 42, 61, 62, 81, 82, 101, 102, 121, 122, 141, 142, 21, 22, 23, 24, 43, 44, 63, 64, 83, 84, 103, 104, 123, 124, 143, 144, 41, 42, 43, 44, 45, 46, 65, 66, 85, 86, 105, 106, 125, 126, 145, 146, 61, 62, 63, 64, 65, 66, 67, 68, 87, 88, 107, 108, 127, 128, 147, 148, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 109, 110, 129, 130, 149, 150, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 131, 132, 151, 152, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 153, 154, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156 } },
    .{ .lower = true, .kind = .hermitian, .complex = true, .m = 8, .want = &.{ 1, 0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 3, -4, 23, 0, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 5, -6, 25, -26, 45, 0, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 7, -8, 27, -28, 47, -48, 67, 0, 69, 70, 71, 72, 73, 74, 75, 76, 9, -10, 29, -30, 49, -50, 69, -70, 89, 0, 91, 92, 93, 94, 95, 96, 11, -12, 31, -32, 51, -52, 71, -72, 91, -92, 111, 0, 113, 114, 115, 116, 13, -14, 33, -34, 53, -54, 73, -74, 93, -94, 113, -114, 133, 0, 135, 136, 15, -16, 35, -36, 55, -56, 75, -76, 95, -96, 115, -116, 135, -136, 155, 0 } },
    .{ .lower = false, .kind = .hermitian, .complex = true, .m = 8, .want = &.{ 1, 0, 21, -22, 41, -42, 61, -62, 81, -82, 101, -102, 121, -122, 141, -142, 21, 22, 23, 0, 43, -44, 63, -64, 83, -84, 103, -104, 123, -124, 143, -144, 41, 42, 43, 44, 45, 0, 65, -66, 85, -86, 105, -106, 125, -126, 145, -146, 61, 62, 63, 64, 65, 66, 67, 0, 87, -88, 107, -108, 127, -128, 147, -148, 81, 82, 83, 84, 85, 86, 87, 88, 89, 0, 109, -110, 129, -130, 149, -150, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 0, 131, -132, 151, -152, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 0, 153, -154, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 0 } },
    .{ .lower = true, .kind = .hermitian_rev, .complex = true, .m = 8, .want = &.{ 1, 0, 3, -4, 5, -6, 7, -8, 9, -10, 11, -12, 13, -14, 15, -16, 3, 4, 23, 0, 25, -26, 27, -28, 29, -30, 31, -32, 33, -34, 35, -36, 5, 6, 25, 26, 45, 0, 47, -48, 49, -50, 51, -52, 53, -54, 55, -56, 7, 8, 27, 28, 47, 48, 67, 0, 69, -70, 71, -72, 73, -74, 75, -76, 9, 10, 29, 30, 49, 50, 69, 70, 89, 0, 91, -92, 93, -94, 95, -96, 11, 12, 31, 32, 51, 52, 71, 72, 91, 92, 111, 0, 113, -114, 115, -116, 13, 14, 33, 34, 53, 54, 73, 74, 93, 94, 113, 114, 133, 0, 135, -136, 15, 16, 35, 36, 55, 56, 75, 76, 95, 96, 115, 116, 135, 136, 155, 0 } },
    .{ .lower = false, .kind = .hermitian_rev, .complex = true, .m = 8, .want = &.{ 1, 0, 21, 22, 41, 42, 61, 62, 81, 82, 101, 102, 121, 122, 141, 142, 21, -22, 23, 0, 43, 44, 63, 64, 83, 84, 103, 104, 123, 124, 143, 144, 41, -42, 43, -44, 45, 0, 65, 66, 85, 86, 105, 106, 125, 126, 145, 146, 61, -62, 63, -64, 65, -66, 67, 0, 87, 88, 107, 108, 127, 128, 147, 148, 81, -82, 83, -84, 85, -86, 87, -88, 89, 0, 109, 110, 129, 130, 149, 150, 101, -102, 103, -104, 105, -106, 107, -108, 109, -110, 111, 0, 131, 132, 151, 152, 121, -122, 123, -124, 125, -126, 127, -128, 129, -130, 131, -132, 133, 0, 153, 154, 141, -142, 143, -144, 145, -146, 147, -148, 149, -150, 151, -152, 153, -154, 155, 0 } },
    .{ .lower = true, .kind = .symmetric, .complex = false, .m = 9, .want = &.{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 2, 13, 14, 15, 16, 17, 18, 19, 20, 3, 14, 25, 26, 27, 28, 29, 30, 31, 4, 15, 26, 37, 38, 39, 40, 41, 42, 5, 16, 27, 38, 49, 50, 51, 52, 53, 6, 17, 28, 39, 50, 61, 62, 63, 64, 7, 18, 29, 40, 51, 62, 73, 74, 75, 8, 19, 30, 41, 52, 63, 74, 85, 86, 9, 20, 31, 42, 53, 64, 75, 86, 97 } },
    .{ .lower = false, .kind = .symmetric, .complex = false, .m = 9, .want = &.{ 1, 12, 23, 34, 45, 56, 67, 78, 89, 12, 13, 24, 35, 46, 57, 68, 79, 90, 23, 24, 25, 36, 47, 58, 69, 80, 91, 34, 35, 36, 37, 48, 59, 70, 81, 92, 45, 46, 47, 48, 49, 60, 71, 82, 93, 56, 57, 58, 59, 60, 61, 72, 83, 94, 67, 68, 69, 70, 71, 72, 73, 84, 95, 78, 79, 80, 81, 82, 83, 84, 85, 96, 89, 90, 91, 92, 93, 94, 95, 96, 97 } },
    .{ .lower = true, .kind = .symmetric, .complex = true, .m = 9, .want = &.{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 3, 4, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 5, 6, 27, 28, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 7, 8, 29, 30, 51, 52, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 9, 10, 31, 32, 53, 54, 75, 76, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 11, 12, 33, 34, 55, 56, 77, 78, 99, 100, 121, 122, 123, 124, 125, 126, 127, 128, 13, 14, 35, 36, 57, 58, 79, 80, 101, 102, 123, 124, 145, 146, 147, 148, 149, 150, 15, 16, 37, 38, 59, 60, 81, 82, 103, 104, 125, 126, 147, 148, 169, 170, 171, 172, 17, 18, 39, 40, 61, 62, 83, 84, 105, 106, 127, 128, 149, 150, 171, 172, 193, 194 } },
    .{ .lower = false, .kind = .symmetric, .complex = true, .m = 9, .want = &.{ 1, 2, 23, 24, 45, 46, 67, 68, 89, 90, 111, 112, 133, 134, 155, 156, 177, 178, 23, 24, 25, 26, 47, 48, 69, 70, 91, 92, 113, 114, 135, 136, 157, 158, 179, 180, 45, 46, 47, 48, 49, 50, 71, 72, 93, 94, 115, 116, 137, 138, 159, 160, 181, 182, 67, 68, 69, 70, 71, 72, 73, 74, 95, 96, 117, 118, 139, 140, 161, 162, 183, 184, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 119, 120, 141, 142, 163, 164, 185, 186, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 143, 144, 165, 166, 187, 188, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 167, 168, 189, 190, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 191, 192, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194 } },
    .{ .lower = true, .kind = .hermitian, .complex = true, .m = 9, .want = &.{ 1, 0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 3, -4, 25, 0, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 5, -6, 27, -28, 49, 0, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 7, -8, 29, -30, 51, -52, 73, 0, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 9, -10, 31, -32, 53, -54, 75, -76, 97, 0, 99, 100, 101, 102, 103, 104, 105, 106, 11, -12, 33, -34, 55, -56, 77, -78, 99, -100, 121, 0, 123, 124, 125, 126, 127, 128, 13, -14, 35, -36, 57, -58, 79, -80, 101, -102, 123, -124, 145, 0, 147, 148, 149, 150, 15, -16, 37, -38, 59, -60, 81, -82, 103, -104, 125, -126, 147, -148, 169, 0, 171, 172, 17, -18, 39, -40, 61, -62, 83, -84, 105, -106, 127, -128, 149, -150, 171, -172, 193, 0 } },
    .{ .lower = false, .kind = .hermitian, .complex = true, .m = 9, .want = &.{ 1, 0, 23, -24, 45, -46, 67, -68, 89, -90, 111, -112, 133, -134, 155, -156, 177, -178, 23, 24, 25, 0, 47, -48, 69, -70, 91, -92, 113, -114, 135, -136, 157, -158, 179, -180, 45, 46, 47, 48, 49, 0, 71, -72, 93, -94, 115, -116, 137, -138, 159, -160, 181, -182, 67, 68, 69, 70, 71, 72, 73, 0, 95, -96, 117, -118, 139, -140, 161, -162, 183, -184, 89, 90, 91, 92, 93, 94, 95, 96, 97, 0, 119, -120, 141, -142, 163, -164, 185, -186, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 0, 143, -144, 165, -166, 187, -188, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 0, 167, -168, 189, -190, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 0, 191, -192, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 0 } },
    .{ .lower = true, .kind = .hermitian_rev, .complex = true, .m = 9, .want = &.{ 1, 0, 3, -4, 5, -6, 7, -8, 9, -10, 11, -12, 13, -14, 15, -16, 17, -18, 3, 4, 25, 0, 27, -28, 29, -30, 31, -32, 33, -34, 35, -36, 37, -38, 39, -40, 5, 6, 27, 28, 49, 0, 51, -52, 53, -54, 55, -56, 57, -58, 59, -60, 61, -62, 7, 8, 29, 30, 51, 52, 73, 0, 75, -76, 77, -78, 79, -80, 81, -82, 83, -84, 9, 10, 31, 32, 53, 54, 75, 76, 97, 0, 99, -100, 101, -102, 103, -104, 105, -106, 11, 12, 33, 34, 55, 56, 77, 78, 99, 100, 121, 0, 123, -124, 125, -126, 127, -128, 13, 14, 35, 36, 57, 58, 79, 80, 101, 102, 123, 124, 145, 0, 147, -148, 149, -150, 15, 16, 37, 38, 59, 60, 81, 82, 103, 104, 125, 126, 147, 148, 169, 0, 171, -172, 17, 18, 39, 40, 61, 62, 83, 84, 105, 106, 127, 128, 149, 150, 171, 172, 193, 0 } },
    .{ .lower = false, .kind = .hermitian_rev, .complex = true, .m = 9, .want = &.{ 1, 0, 23, 24, 45, 46, 67, 68, 89, 90, 111, 112, 133, 134, 155, 156, 177, 178, 23, -24, 25, 0, 47, 48, 69, 70, 91, 92, 113, 114, 135, 136, 157, 158, 179, 180, 45, -46, 47, -48, 49, 0, 71, 72, 93, 94, 115, 116, 137, 138, 159, 160, 181, 182, 67, -68, 69, -70, 71, -72, 73, 0, 95, 96, 117, 118, 139, 140, 161, 162, 183, 184, 89, -90, 91, -92, 93, -94, 95, -96, 97, 0, 119, 120, 141, 142, 163, 164, 185, 186, 111, -112, 113, -114, 115, -116, 117, -118, 119, -120, 121, 0, 143, 144, 165, 166, 187, 188, 133, -134, 135, -136, 137, -138, 139, -140, 141, -142, 143, -144, 145, 0, 167, 168, 189, 190, 155, -156, 157, -158, 159, -160, 161, -162, 163, -164, 165, -166, 167, -168, 169, 0, 191, 192, 177, -178, 179, -180, 181, -182, 183, -184, 185, -186, 187, -188, 189, -190, 191, -192, 193, 0 } },
    .{ .lower = true, .kind = .symmetric, .complex = false, .m = 10, .want = &.{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 2, 14, 15, 16, 17, 18, 19, 20, 21, 22, 3, 15, 27, 28, 29, 30, 31, 32, 33, 34, 4, 16, 28, 40, 41, 42, 43, 44, 45, 46, 5, 17, 29, 41, 53, 54, 55, 56, 57, 58, 6, 18, 30, 42, 54, 66, 67, 68, 69, 70, 7, 19, 31, 43, 55, 67, 79, 80, 81, 82, 8, 20, 32, 44, 56, 68, 80, 92, 93, 94, 9, 21, 33, 45, 57, 69, 81, 93, 105, 106, 10, 22, 34, 46, 58, 70, 82, 94, 106, 118 } },
    .{ .lower = false, .kind = .symmetric, .complex = false, .m = 10, .want = &.{ 1, 13, 25, 37, 49, 61, 73, 85, 97, 109, 13, 14, 26, 38, 50, 62, 74, 86, 98, 110, 25, 26, 27, 39, 51, 63, 75, 87, 99, 111, 37, 38, 39, 40, 52, 64, 76, 88, 100, 112, 49, 50, 51, 52, 53, 65, 77, 89, 101, 113, 61, 62, 63, 64, 65, 66, 78, 90, 102, 114, 73, 74, 75, 76, 77, 78, 79, 91, 103, 115, 85, 86, 87, 88, 89, 90, 91, 92, 104, 116, 97, 98, 99, 100, 101, 102, 103, 104, 105, 117, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118 } },
    .{ .lower = true, .kind = .symmetric, .complex = true, .m = 10, .want = &.{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 3, 4, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 5, 6, 29, 30, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 7, 8, 31, 32, 55, 56, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 9, 10, 33, 34, 57, 58, 81, 82, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 11, 12, 35, 36, 59, 60, 83, 84, 107, 108, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 13, 14, 37, 38, 61, 62, 85, 86, 109, 110, 133, 134, 157, 158, 159, 160, 161, 162, 163, 164, 15, 16, 39, 40, 63, 64, 87, 88, 111, 112, 135, 136, 159, 160, 183, 184, 185, 186, 187, 188, 17, 18, 41, 42, 65, 66, 89, 90, 113, 114, 137, 138, 161, 162, 185, 186, 209, 210, 211, 212, 19, 20, 43, 44, 67, 68, 91, 92, 115, 116, 139, 140, 163, 164, 187, 188, 211, 212, 235, 236 } },
    .{ .lower = false, .kind = .symmetric, .complex = true, .m = 10, .want = &.{ 1, 2, 25, 26, 49, 50, 73, 74, 97, 98, 121, 122, 145, 146, 169, 170, 193, 194, 217, 218, 25, 26, 27, 28, 51, 52, 75, 76, 99, 100, 123, 124, 147, 148, 171, 172, 195, 196, 219, 220, 49, 50, 51, 52, 53, 54, 77, 78, 101, 102, 125, 126, 149, 150, 173, 174, 197, 198, 221, 222, 73, 74, 75, 76, 77, 78, 79, 80, 103, 104, 127, 128, 151, 152, 175, 176, 199, 200, 223, 224, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 129, 130, 153, 154, 177, 178, 201, 202, 225, 226, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 155, 156, 179, 180, 203, 204, 227, 228, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 181, 182, 205, 206, 229, 230, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 207, 208, 231, 232, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 233, 234, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236 } },
    .{ .lower = true, .kind = .hermitian, .complex = true, .m = 10, .want = &.{ 1, 0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 3, -4, 27, 0, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 5, -6, 29, -30, 53, 0, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 7, -8, 31, -32, 55, -56, 79, 0, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 9, -10, 33, -34, 57, -58, 81, -82, 105, 0, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 11, -12, 35, -36, 59, -60, 83, -84, 107, -108, 131, 0, 133, 134, 135, 136, 137, 138, 139, 140, 13, -14, 37, -38, 61, -62, 85, -86, 109, -110, 133, -134, 157, 0, 159, 160, 161, 162, 163, 164, 15, -16, 39, -40, 63, -64, 87, -88, 111, -112, 135, -136, 159, -160, 183, 0, 185, 186, 187, 188, 17, -18, 41, -42, 65, -66, 89, -90, 113, -114, 137, -138, 161, -162, 185, -186, 209, 0, 211, 212, 19, -20, 43, -44, 67, -68, 91, -92, 115, -116, 139, -140, 163, -164, 187, -188, 211, -212, 235, 0 } },
    .{ .lower = false, .kind = .hermitian, .complex = true, .m = 10, .want = &.{ 1, 0, 25, -26, 49, -50, 73, -74, 97, -98, 121, -122, 145, -146, 169, -170, 193, -194, 217, -218, 25, 26, 27, 0, 51, -52, 75, -76, 99, -100, 123, -124, 147, -148, 171, -172, 195, -196, 219, -220, 49, 50, 51, 52, 53, 0, 77, -78, 101, -102, 125, -126, 149, -150, 173, -174, 197, -198, 221, -222, 73, 74, 75, 76, 77, 78, 79, 0, 103, -104, 127, -128, 151, -152, 175, -176, 199, -200, 223, -224, 97, 98, 99, 100, 101, 102, 103, 104, 105, 0, 129, -130, 153, -154, 177, -178, 201, -202, 225, -226, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 0, 155, -156, 179, -180, 203, -204, 227, -228, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 0, 181, -182, 205, -206, 229, -230, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 0, 207, -208, 231, -232, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 0, 233, -234, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 0 } },
    .{ .lower = true, .kind = .hermitian_rev, .complex = true, .m = 10, .want = &.{ 1, 0, 3, -4, 5, -6, 7, -8, 9, -10, 11, -12, 13, -14, 15, -16, 17, -18, 19, -20, 3, 4, 27, 0, 29, -30, 31, -32, 33, -34, 35, -36, 37, -38, 39, -40, 41, -42, 43, -44, 5, 6, 29, 30, 53, 0, 55, -56, 57, -58, 59, -60, 61, -62, 63, -64, 65, -66, 67, -68, 7, 8, 31, 32, 55, 56, 79, 0, 81, -82, 83, -84, 85, -86, 87, -88, 89, -90, 91, -92, 9, 10, 33, 34, 57, 58, 81, 82, 105, 0, 107, -108, 109, -110, 111, -112, 113, -114, 115, -116, 11, 12, 35, 36, 59, 60, 83, 84, 107, 108, 131, 0, 133, -134, 135, -136, 137, -138, 139, -140, 13, 14, 37, 38, 61, 62, 85, 86, 109, 110, 133, 134, 157, 0, 159, -160, 161, -162, 163, -164, 15, 16, 39, 40, 63, 64, 87, 88, 111, 112, 135, 136, 159, 160, 183, 0, 185, -186, 187, -188, 17, 18, 41, 42, 65, 66, 89, 90, 113, 114, 137, 138, 161, 162, 185, 186, 209, 0, 211, -212, 19, 20, 43, 44, 67, 68, 91, 92, 115, 116, 139, 140, 163, 164, 187, 188, 211, 212, 235, 0 } },
    .{ .lower = false, .kind = .hermitian_rev, .complex = true, .m = 10, .want = &.{ 1, 0, 25, 26, 49, 50, 73, 74, 97, 98, 121, 122, 145, 146, 169, 170, 193, 194, 217, 218, 25, -26, 27, 0, 51, 52, 75, 76, 99, 100, 123, 124, 147, 148, 171, 172, 195, 196, 219, 220, 49, -50, 51, -52, 53, 0, 77, 78, 101, 102, 125, 126, 149, 150, 173, 174, 197, 198, 221, 222, 73, -74, 75, -76, 77, -78, 79, 0, 103, 104, 127, 128, 151, 152, 175, 176, 199, 200, 223, 224, 97, -98, 99, -100, 101, -102, 103, -104, 105, 0, 129, 130, 153, 154, 177, 178, 201, 202, 225, 226, 121, -122, 123, -124, 125, -126, 127, -128, 129, -130, 131, 0, 155, 156, 179, 180, 203, 204, 227, 228, 145, -146, 147, -148, 149, -150, 151, -152, 153, -154, 155, -156, 157, 0, 181, 182, 205, 206, 229, 230, 169, -170, 171, -172, 173, -174, 175, -176, 177, -178, 179, -180, 181, -182, 183, 0, 207, 208, 231, 232, 193, -194, 195, -196, 197, -198, 199, -200, 201, -202, 203, -204, 205, -206, 207, -208, 209, 0, 233, 234, 217, -218, 219, -220, 221, -222, 223, -224, 225, -226, 227, -228, 229, -230, 231, -232, 233, -234, 235, 0 } },
    .{ .lower = true, .kind = .symmetric, .complex = false, .m = 11, .want = &.{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 2, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 3, 16, 29, 30, 31, 32, 33, 34, 35, 36, 37, 4, 17, 30, 43, 44, 45, 46, 47, 48, 49, 50, 5, 18, 31, 44, 57, 58, 59, 60, 61, 62, 63, 6, 19, 32, 45, 58, 71, 72, 73, 74, 75, 76, 7, 20, 33, 46, 59, 72, 85, 86, 87, 88, 89, 8, 21, 34, 47, 60, 73, 86, 99, 100, 101, 102, 9, 22, 35, 48, 61, 74, 87, 100, 113, 114, 115, 10, 23, 36, 49, 62, 75, 88, 101, 114, 127, 128, 11, 24, 37, 50, 63, 76, 89, 102, 115, 128, 141 } },
    .{ .lower = false, .kind = .symmetric, .complex = false, .m = 11, .want = &.{ 1, 14, 27, 40, 53, 66, 79, 92, 105, 118, 131, 14, 15, 28, 41, 54, 67, 80, 93, 106, 119, 132, 27, 28, 29, 42, 55, 68, 81, 94, 107, 120, 133, 40, 41, 42, 43, 56, 69, 82, 95, 108, 121, 134, 53, 54, 55, 56, 57, 70, 83, 96, 109, 122, 135, 66, 67, 68, 69, 70, 71, 84, 97, 110, 123, 136, 79, 80, 81, 82, 83, 84, 85, 98, 111, 124, 137, 92, 93, 94, 95, 96, 97, 98, 99, 112, 125, 138, 105, 106, 107, 108, 109, 110, 111, 112, 113, 126, 139, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 140, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141 } },
    .{ .lower = true, .kind = .symmetric, .complex = true, .m = 11, .want = &.{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 3, 4, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 5, 6, 31, 32, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 7, 8, 33, 34, 59, 60, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 9, 10, 35, 36, 61, 62, 87, 88, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 11, 12, 37, 38, 63, 64, 89, 90, 115, 116, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 13, 14, 39, 40, 65, 66, 91, 92, 117, 118, 143, 144, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 15, 16, 41, 42, 67, 68, 93, 94, 119, 120, 145, 146, 171, 172, 197, 198, 199, 200, 201, 202, 203, 204, 17, 18, 43, 44, 69, 70, 95, 96, 121, 122, 147, 148, 173, 174, 199, 200, 225, 226, 227, 228, 229, 230, 19, 20, 45, 46, 71, 72, 97, 98, 123, 124, 149, 150, 175, 176, 201, 202, 227, 228, 253, 254, 255, 256, 21, 22, 47, 48, 73, 74, 99, 100, 125, 126, 151, 152, 177, 178, 203, 204, 229, 230, 255, 256, 281, 282 } },
    .{ .lower = false, .kind = .symmetric, .complex = true, .m = 11, .want = &.{ 1, 2, 27, 28, 53, 54, 79, 80, 105, 106, 131, 132, 157, 158, 183, 184, 209, 210, 235, 236, 261, 262, 27, 28, 29, 30, 55, 56, 81, 82, 107, 108, 133, 134, 159, 160, 185, 186, 211, 212, 237, 238, 263, 264, 53, 54, 55, 56, 57, 58, 83, 84, 109, 110, 135, 136, 161, 162, 187, 188, 213, 214, 239, 240, 265, 266, 79, 80, 81, 82, 83, 84, 85, 86, 111, 112, 137, 138, 163, 164, 189, 190, 215, 216, 241, 242, 267, 268, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 139, 140, 165, 166, 191, 192, 217, 218, 243, 244, 269, 270, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 167, 168, 193, 194, 219, 220, 245, 246, 271, 272, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 195, 196, 221, 222, 247, 248, 273, 274, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 223, 224, 249, 250, 275, 276, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 251, 252, 277, 278, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 279, 280, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282 } },
    .{ .lower = true, .kind = .hermitian, .complex = true, .m = 11, .want = &.{ 1, 0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 3, -4, 29, 0, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 5, -6, 31, -32, 57, 0, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 7, -8, 33, -34, 59, -60, 85, 0, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 9, -10, 35, -36, 61, -62, 87, -88, 113, 0, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 11, -12, 37, -38, 63, -64, 89, -90, 115, -116, 141, 0, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 13, -14, 39, -40, 65, -66, 91, -92, 117, -118, 143, -144, 169, 0, 171, 172, 173, 174, 175, 176, 177, 178, 15, -16, 41, -42, 67, -68, 93, -94, 119, -120, 145, -146, 171, -172, 197, 0, 199, 200, 201, 202, 203, 204, 17, -18, 43, -44, 69, -70, 95, -96, 121, -122, 147, -148, 173, -174, 199, -200, 225, 0, 227, 228, 229, 230, 19, -20, 45, -46, 71, -72, 97, -98, 123, -124, 149, -150, 175, -176, 201, -202, 227, -228, 253, 0, 255, 256, 21, -22, 47, -48, 73, -74, 99, -100, 125, -126, 151, -152, 177, -178, 203, -204, 229, -230, 255, -256, 281, 0 } },
    .{ .lower = false, .kind = .hermitian, .complex = true, .m = 11, .want = &.{ 1, 0, 27, -28, 53, -54, 79, -80, 105, -106, 131, -132, 157, -158, 183, -184, 209, -210, 235, -236, 261, -262, 27, 28, 29, 0, 55, -56, 81, -82, 107, -108, 133, -134, 159, -160, 185, -186, 211, -212, 237, -238, 263, -264, 53, 54, 55, 56, 57, 0, 83, -84, 109, -110, 135, -136, 161, -162, 187, -188, 213, -214, 239, -240, 265, -266, 79, 80, 81, 82, 83, 84, 85, 0, 111, -112, 137, -138, 163, -164, 189, -190, 215, -216, 241, -242, 267, -268, 105, 106, 107, 108, 109, 110, 111, 112, 113, 0, 139, -140, 165, -166, 191, -192, 217, -218, 243, -244, 269, -270, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 0, 167, -168, 193, -194, 219, -220, 245, -246, 271, -272, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 0, 195, -196, 221, -222, 247, -248, 273, -274, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 0, 223, -224, 249, -250, 275, -276, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 0, 251, -252, 277, -278, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 0, 279, -280, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 0 } },
    .{ .lower = true, .kind = .hermitian_rev, .complex = true, .m = 11, .want = &.{ 1, 0, 3, -4, 5, -6, 7, -8, 9, -10, 11, -12, 13, -14, 15, -16, 17, -18, 19, -20, 21, -22, 3, 4, 29, 0, 31, -32, 33, -34, 35, -36, 37, -38, 39, -40, 41, -42, 43, -44, 45, -46, 47, -48, 5, 6, 31, 32, 57, 0, 59, -60, 61, -62, 63, -64, 65, -66, 67, -68, 69, -70, 71, -72, 73, -74, 7, 8, 33, 34, 59, 60, 85, 0, 87, -88, 89, -90, 91, -92, 93, -94, 95, -96, 97, -98, 99, -100, 9, 10, 35, 36, 61, 62, 87, 88, 113, 0, 115, -116, 117, -118, 119, -120, 121, -122, 123, -124, 125, -126, 11, 12, 37, 38, 63, 64, 89, 90, 115, 116, 141, 0, 143, -144, 145, -146, 147, -148, 149, -150, 151, -152, 13, 14, 39, 40, 65, 66, 91, 92, 117, 118, 143, 144, 169, 0, 171, -172, 173, -174, 175, -176, 177, -178, 15, 16, 41, 42, 67, 68, 93, 94, 119, 120, 145, 146, 171, 172, 197, 0, 199, -200, 201, -202, 203, -204, 17, 18, 43, 44, 69, 70, 95, 96, 121, 122, 147, 148, 173, 174, 199, 200, 225, 0, 227, -228, 229, -230, 19, 20, 45, 46, 71, 72, 97, 98, 123, 124, 149, 150, 175, 176, 201, 202, 227, 228, 253, 0, 255, -256, 21, 22, 47, 48, 73, 74, 99, 100, 125, 126, 151, 152, 177, 178, 203, 204, 229, 230, 255, 256, 281, 0 } },
    .{ .lower = false, .kind = .hermitian_rev, .complex = true, .m = 11, .want = &.{ 1, 0, 27, 28, 53, 54, 79, 80, 105, 106, 131, 132, 157, 158, 183, 184, 209, 210, 235, 236, 261, 262, 27, -28, 29, 0, 55, 56, 81, 82, 107, 108, 133, 134, 159, 160, 185, 186, 211, 212, 237, 238, 263, 264, 53, -54, 55, -56, 57, 0, 83, 84, 109, 110, 135, 136, 161, 162, 187, 188, 213, 214, 239, 240, 265, 266, 79, -80, 81, -82, 83, -84, 85, 0, 111, 112, 137, 138, 163, 164, 189, 190, 215, 216, 241, 242, 267, 268, 105, -106, 107, -108, 109, -110, 111, -112, 113, 0, 139, 140, 165, 166, 191, 192, 217, 218, 243, 244, 269, 270, 131, -132, 133, -134, 135, -136, 137, -138, 139, -140, 141, 0, 167, 168, 193, 194, 219, 220, 245, 246, 271, 272, 157, -158, 159, -160, 161, -162, 163, -164, 165, -166, 167, -168, 169, 0, 195, 196, 221, 222, 247, 248, 273, 274, 183, -184, 185, -186, 187, -188, 189, -190, 191, -192, 193, -194, 195, -196, 197, 0, 223, 224, 249, 250, 275, 276, 209, -210, 211, -212, 213, -214, 215, -216, 217, -218, 219, -220, 221, -222, 223, -224, 225, 0, 251, 252, 277, 278, 235, -236, 237, -238, 239, -240, 241, -242, 243, -244, 245, -246, 247, -248, 249, -250, 251, -252, 253, 0, 279, 280, 261, -262, 263, -264, 265, -266, 267, -268, 269, -270, 271, -272, 273, -274, 275, -276, 277, -278, 279, -280, 281, 0 } },
};

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
