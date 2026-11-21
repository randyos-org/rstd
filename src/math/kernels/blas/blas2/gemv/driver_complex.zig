//! The complex gemv drivers (`driver.zig`): `{c,z}gemv_n_4.c` and `{c,z}gemv_t_4.c`.

const util = @import("../../util.zig");
const types = @import("../../types.zig");
const fma = @import("../../fma.zig");
const driver = @import("driver.zig");

const BlasLong = types.BlasLong;
const Slot = util.Slot;
const off = driver.off;
const choose = driver.choose;

/// What a complex slot is built as: `trans` is `TRANS`, `conj` `CONJ`
/// (`A` conjugated), `xconj` `XCONJ` (`x` conjugated, written as the product
/// conjugated, then `alpha` applied to its conjugate).
const Cform = struct { trans: bool, conj: bool, xconj: bool };

fn cform(comptime slot: Slot) Cform {
    const t: types.Trans = switch (slot) {
        .cgemv_n, .zgemv_n => .n,
        .cgemv_t, .zgemv_t => .t,
        .cgemv_r, .zgemv_r => .r,
        .cgemv_c, .zgemv_c => .c,
        .cgemv_o, .zgemv_o => .o,
        .cgemv_u, .zgemv_u => .u,
        .cgemv_s, .zgemv_s => .s,
        .cgemv_d, .zgemv_d => .d,
        else => @compileError(@tagName(slot) ++ " is not a complex gemv kernel"),
    };
    return .{ .trans = t.transposed(), .conj = t.conjA(), .xconj = t.conjX() };
}

/// The float type of complex `slot`.
fn Cfloat(comptime slot: Slot) type {
    return if (@tagName(slot)[0] == 'c') f32 else f64;
}

/// The suffix of the converted blocks' form (`uk.zig`'s `Form.suffix`).
fn suffix(comptime f: Cform) []const u8 {
    return if (f.conj and f.xconj) "_cx" else if (f.conj) "_c" else if (f.xconj) "_x" else "";
}

/// `a*b + c*d` and `a*b - c*d` as clang contracts them: the first product
/// fused, the second rounded.
inline fn sum2(comptime F: type, comptime fuse: bool, a: F, b: F, c: F, d: F) F {
    return fma.madd(fuse, F, a, b, c * d);
}
inline fn diff2(comptime F: type, comptime fuse: bool, a: F, b: F, c: F, d: F) F {
    return fma.madd(fuse, F, a, b, -(c * d));
}

/// One element's product terms as the C writes them: `(ar*xr - ai*xi,
/// ar*xi + ai*xr)` where `CONJ` and `XCONJ` agree, else `(ar*xr + ai*xi,
/// ar*xi - ai*xr)`.
inline fn cterm(comptime F: type, comptime fuse: bool, comptime same: bool, ar: F, ai: F, xr: F, xi: F) [2]F {
    return if (same)
        .{ diff2(F, fuse, ar, xr, ai, xi), sum2(F, fuse, ar, xi, ai, xr) }
    else
        .{ sum2(F, fuse, ar, xr, ai, xi), diff2(F, fuse, ar, xi, ai, xr) };
}

/// `y += alpha * t` as the C writes it: plain, or with `XCONJ`
/// `y_r += ar*tr + ai*ti`, `y_i -= ar*ti - ai*tr`.
inline fn cupdate(comptime F: type, comptime fuse: bool, comptime xconj: bool, y: [*]F, ar: F, ai: F, tr: F, ti: F) void {
    if (xconj) {
        y[0] += sum2(F, fuse, ar, tr, ai, ti);
        y[1] -= diff2(F, fuse, ar, ti, ai, tr);
    } else {
        y[0] += diff2(F, fuse, ar, tr, ai, ti);
        y[1] += sum2(F, fuse, ar, ti, ai, tr);
    }
}

/// The C blocks of `{c,z}gemv_n_4.c` and `{c,z}gemv_t_4.c` and `n`'s `add_y`.
fn CFallback(comptime F: type, comptime f: Cform, comptime fuse: bool) type {
    const same = f.conj == f.xconj;
    return struct {
        /// `kernel_4x<cols>` of `n`: `y += A x` over `cols` columns, column by
        /// column for each row.
        fn nBlock(comptime cols: comptime_int, rows: BlasLong, ap: [*]const [*]F, x: [*]F, y: [*]F) void {
            for (0..@intCast(rows)) |r| {
                const i = 2 * r;
                inline for (0..cols) |c| {
                    const t = cterm(F, fuse, same, ap[c][i], ap[c][i + 1], x[2 * c], x[2 * c + 1]);
                    y[i] += t[0];
                    y[i + 1] += t[1];
                }
            }
        }
        fn n4x4(rows: BlasLong, ap: [*]const [*]F, x: [*]F, y: [*]F) void {
            nBlock(4, rows, ap, x, y);
        }
        fn n4x2(rows: BlasLong, ap: [*]const [*]F, x: [*]F, y: [*]F) void {
            nBlock(2, rows, ap, x, y);
        }
        fn n4x1(rows: BlasLong, ap: [*]F, x: [*]F, y: [*]F) void {
            const one = [1][*]F{ap};
            nBlock(1, rows, &one, x, y);
        }

        /// `kernel_4x<cols>` of `t`: `cols` dot products, then `alpha` into `y`.
        fn tBlock(comptime cols: comptime_int, rows: BlasLong, ap: [*]const [*]F, x: [*]F, y: [*]F, alpha: [*]const F) void {
            var tr: [cols]F = @splat(0);
            var ti: [cols]F = @splat(0);
            for (0..@intCast(rows)) |r| {
                const i = 2 * r;
                inline for (0..cols) |c| {
                    const v = cterm(F, fuse, same, ap[c][i], ap[c][i + 1], x[i], x[i + 1]);
                    tr[c] += v[0];
                    ti[c] += v[1];
                }
            }
            inline for (0..cols) |c| cupdate(F, fuse, f.xconj, y + 2 * c, alpha[0], alpha[1], tr[c], ti[c]);
        }
        fn t4x4(rows: BlasLong, ap: [*]const [*]F, x: [*]F, y: [*]F, alpha: [*]const F) void {
            tBlock(4, rows, ap, x, y, alpha);
        }
        fn t4x2(rows: BlasLong, ap: [*]const [*]F, x: [*]F, y: [*]F, alpha: [*]const F) void {
            tBlock(2, rows, ap, x, y, alpha);
        }
        fn t4x1(rows: BlasLong, ap: [*]F, x: [*]F, y: [*]F, alpha: [*]const F) void {
            const one = [1][*]F{ap};
            tBlock(1, rows, &one, x, y, alpha);
        }

        /// `n`'s `add_y`: `dest += alpha * src`, or with `XCONJ` the product's
        /// conjugate, `(ar*sr + ai*si, -ar*si + ai*sr)`.
        fn addY(count: BlasLong, src: [*]F, dest: [*]F, inc_dest: BlasLong, ar: F, ai: F) void {
            var d: BlasLong = 0;
            for (0..@intCast(count)) |k| {
                const sr = src[2 * k];
                const si = src[2 * k + 1];
                const tr = if (f.xconj) sum2(F, fuse, ar, sr, ai, si) else diff2(F, fuse, ar, sr, ai, si);
                const ti = if (f.xconj) sum2(F, fuse, -ar, si, ai, sr) else sum2(F, fuse, ar, si, ai, sr);
                util.at(F, dest, d).* += tr;
                util.at(F, dest, d + 1).* += ti;
                d += inc_dest;
            }
        }
    };
}

fn CN4(comptime F: type) type {
    return *const fn (BlasLong, [*]const [*]F, [*]F, [*]F) void;
}
fn CN1(comptime F: type) type {
    return *const fn (BlasLong, [*]F, [*]F, [*]F) void;
}
fn CT4(comptime F: type) type {
    return *const fn (BlasLong, [*]const [*]F, [*]F, [*]F, [*]const F) void;
}
fn CT1(comptime F: type) type {
    return *const fn (BlasLong, [*]F, [*]F, [*]F, [*]const F) void;
}
fn AddY(comptime F: type) type {
    return *const fn (BlasLong, [*]F, [*]F, BlasLong, F, F) void;
}

/// `{c,z}gemv_n_4.c`'s `CNAME` of `slot`, around the microkernel `M`.
pub fn ZgemvN(comptime slot: Slot, comptime M: type, comptime fuse: bool) type {
    const F = Cfloat(slot);
    const f = cform(slot);
    const p = @tagName(slot)[0..1];
    const sx = suffix(f);
    const C = CFallback(F, f, fuse);
    const k4x4 = choose(CN4(F), M, void, p ++ "gemv_kernel_4x4" ++ sx, &C.n4x4);
    const k4x2 = choose(CN4(F), M, void, p ++ "gemv_kernel_4x2" ++ sx, &C.n4x2);
    const k4x1 = choose(CN1(F), M, void, p ++ "gemv_kernel_4x1" ++ sx, &C.n4x1);
    // the microkernel's `add_y` past its strided C, which is `C.addY`
    const add_asm: ?AddY(F) = if (@hasDecl(M, "add_y" ++ sx ++ "_asm")) &@field(M, "add_y" ++ sx ++ "_asm") else null;
    const nb_max: BlasLong = if (F == f32) 2048 else 1024;
    return struct {
        fn addY(count: BlasLong, src: [*]F, dest: [*]F, inc_dest: BlasLong, ar: F, ai: F) void {
            if (add_asm) |k| if (inc_dest == 2) return k(count, src, dest, inc_dest, ar, ai);
            C.addY(count, src, dest, inc_dest, ar, ai);
        }

        pub fn call(m: BlasLong, n: BlasLong, _: BlasLong, alpha_r: F, alpha_i: F, a_in: [*]F, lda_in: BlasLong, x: [*]F, inc_x_in: BlasLong, y: [*]F, inc_y_in: BlasLong, buffer: [*]F) callconv(.c) i32 {
            if (m < 1 or n < 1) return 0;
            const ybuffer = buffer;
            const inc_x = inc_x_in * 2;
            const inc_y = inc_y_in * 2;
            const lda = lda_in * 2;
            const lda4 = 4 * lda;
            var a = a_in;
            var xbuffer: [8]F = undefined;

            const n1 = @divTrunc(n, 4);
            const n2 = @rem(n, 4);
            const m3 = @rem(m, 4);
            var m1 = m - m3;
            const m2 = @rem(m, nb_max) - m3;

            var y_ptr = y;
            var nb = nb_max;
            while (nb == nb_max) {
                m1 -= nb;
                if (m1 < 0) {
                    if (m2 == 0) break;
                    nb = m2;
                }
                var a_ptr = a;
                var ap: [4][*]F = undefined;
                for (&ap, 0..) |*q, k| q.* = a_ptr + k * @as(usize, @intCast(lda));
                var x_ptr = x;
                @memset(ybuffer[0..@intCast(2 * nb)], 0);

                if (inc_x == 2) {
                    var i: BlasLong = 0;
                    while (i < n1) : (i += 1) {
                        k4x4(nb, &ap, x_ptr, ybuffer);
                        for (&ap) |*q| q.* += @intCast(lda4);
                        a_ptr += @intCast(lda4);
                        x_ptr += 8;
                    }
                    if (n2 & 2 != 0) {
                        k4x2(nb, &ap, x_ptr, ybuffer);
                        x_ptr += 4;
                        a_ptr += @intCast(2 * lda);
                    }
                    if (n2 & 1 != 0) k4x1(nb, a_ptr, x_ptr, ybuffer);
                } else {
                    var i: BlasLong = 0;
                    while (i < n1) : (i += 1) {
                        for (0..4) |k| {
                            xbuffer[2 * k] = x_ptr[0];
                            xbuffer[2 * k + 1] = x_ptr[1];
                            x_ptr = off(F, x_ptr, inc_x);
                        }
                        k4x4(nb, &ap, &xbuffer, ybuffer);
                        for (&ap) |*q| q.* += @intCast(lda4);
                        a_ptr += @intCast(lda4);
                    }
                    i = 0;
                    while (i < n2) : (i += 1) {
                        xbuffer[0] = x_ptr[0];
                        xbuffer[1] = x_ptr[1];
                        x_ptr = off(F, x_ptr, inc_x);
                        k4x1(nb, a_ptr, &xbuffer, ybuffer);
                        a_ptr += @intCast(lda);
                    }
                }
                addY(nb, ybuffer, y_ptr, inc_y, alpha_r, alpha_i);
                a += @intCast(2 * nb);
                y_ptr = off(F, y_ptr, nb * inc_y);
            }

            if (m3 == 0) return 0;
            // the last rows: each a sum over the columns from 0, in the
            // C's order whichever of its loops (packed or strided) runs
            const rows: usize = @intCast(m3);
            var tr: [3]F = @splat(0);
            var ti: [3]F = @splat(0);
            var a_ptr = a;
            var x_ptr = x;
            for (0..@intCast(n)) |_| {
                for (0..rows) |r| {
                    const v = cterm(F, fuse, f.conj == f.xconj, a_ptr[2 * r], a_ptr[2 * r + 1], x_ptr[0], x_ptr[1]);
                    tr[r] += v[0];
                    ti[r] += v[1];
                }
                a_ptr = off(F, a_ptr, lda);
                x_ptr = off(F, x_ptr, inc_x);
            }
            for (0..rows) |r| {
                cupdate(F, fuse, f.xconj, y_ptr, alpha_r, alpha_i, tr[r], ti[r]);
                if (r + 1 < rows) y_ptr = off(F, y_ptr, inc_y);
            }
            return 0;
        }
    };
}

/// `{c,z}gemv_t_4.c`'s `CNAME` of `slot`, around the microkernel `M`.
pub fn ZgemvT(comptime slot: Slot, comptime M: type, comptime fuse: bool) type {
    const F = Cfloat(slot);
    const f = cform(slot);
    const p = @tagName(slot)[0..1];
    const sx = suffix(f);
    const C = CFallback(F, f, fuse);
    const k4x4 = choose(CT4(F), M, void, p ++ "gemv_kernel_4x4" ++ sx, &C.t4x4);
    const k4x2 = choose(CT4(F), M, void, p ++ "gemv_kernel_4x2" ++ sx, &C.t4x2);
    const k4x1 = choose(CT1(F), M, void, p ++ "gemv_kernel_4x1" ++ sx, &C.t4x1);
    const nb_max: BlasLong = if (F == f32) 2048 else 1024;
    return struct {
        pub fn call(m: BlasLong, n: BlasLong, _: BlasLong, alpha_r: F, alpha_i: F, a_in: [*]F, lda_in: BlasLong, x_in: [*]F, inc_x_in: BlasLong, y: [*]F, inc_y_in: BlasLong, buffer: [*]F) callconv(.c) i32 {
            if (m < 1 or n < 1) return 0;
            const inc_x = inc_x_in * 2;
            const inc_y = inc_y_in * 2;
            const lda = lda_in * 2;
            const lda4 = lda * 4;
            var a = a_in;
            var x = x_in;
            var xbuffer = buffer;
            var ybuffer: [8]F = undefined;
            const alpha = [2]F{ alpha_r, alpha_i };

            const n1 = n >> 2;
            const n2 = n & 3;
            const m3 = m & 3;
            var m1 = m - m3;
            const m2 = (m & (nb_max - 1)) - m3;

            var nb = nb_max;
            while (nb == nb_max) {
                m1 -= nb;
                if (m1 < 0) {
                    if (m2 == 0) break;
                    nb = m2;
                }
                var y_ptr = y;
                var a_ptr = a;
                var ap: [4][*]F = undefined;
                for (&ap, 0..) |*q, k| q.* = a_ptr + k * @as(usize, @intCast(lda));
                if (inc_x != 2) {
                    var s: BlasLong = 0;
                    for (0..@intCast(nb)) |k| {
                        xbuffer[2 * k] = util.at(F, x, s).*;
                        xbuffer[2 * k + 1] = util.at(F, x, s + 1).*;
                        s += inc_x;
                    }
                } else xbuffer = x;

                if (inc_y == 2) {
                    var i: BlasLong = 0;
                    while (i < n1) : (i += 1) {
                        k4x4(nb, &ap, xbuffer, y_ptr, &alpha);
                        for (&ap) |*q| q.* += @intCast(lda4);
                        a_ptr += @intCast(lda4);
                        y_ptr += 8;
                    }
                    if (n2 & 2 != 0) {
                        k4x2(nb, &ap, xbuffer, y_ptr, &alpha);
                        a_ptr += @intCast(lda * 2);
                        y_ptr += 4;
                    }
                    if (n2 & 1 != 0) k4x1(nb, a_ptr, xbuffer, y_ptr, &alpha);
                } else {
                    var i: BlasLong = 0;
                    while (i < n1) : (i += 1) {
                        ybuffer = @splat(0);
                        k4x4(nb, &ap, xbuffer, &ybuffer, &alpha);
                        for (&ap) |*q| q.* += @intCast(lda4);
                        a_ptr += @intCast(lda4);
                        for (0..4) |k| {
                            y_ptr[0] += ybuffer[2 * k];
                            y_ptr[1] += ybuffer[2 * k + 1];
                            y_ptr = off(F, y_ptr, inc_y);
                        }
                    }
                    i = 0;
                    while (i < n2) : (i += 1) {
                        ybuffer = @splat(0);
                        k4x1(nb, a_ptr, xbuffer, &ybuffer, &alpha);
                        a_ptr += @intCast(lda);
                        y_ptr[0] += ybuffer[0];
                        y_ptr[1] += ybuffer[1];
                        y_ptr = off(F, y_ptr, inc_y);
                    }
                }
                a += @intCast(2 * nb);
                x = off(F, x, nb * inc_x);
            }

            if (m3 == 0) return 0;
            // the last rows: per column, the first row's term, the others added
            const rows: usize = @intCast(m3);
            var xs: [3][2]F = undefined;
            var x_ptr = x;
            for (0..rows) |r| {
                xs[r] = .{ x_ptr[0], x_ptr[1] };
                if (r + 1 < rows) x_ptr = off(F, x_ptr, inc_x);
            }
            var a_ptr = a;
            var y_ptr = y;
            for (0..@intCast(n)) |_| {
                var t = cterm(F, fuse, f.conj == f.xconj, a_ptr[0], a_ptr[1], xs[0][0], xs[0][1]);
                for (1..rows) |r| {
                    const v = cterm(F, fuse, f.conj == f.xconj, a_ptr[2 * r], a_ptr[2 * r + 1], xs[r][0], xs[r][1]);
                    t[0] += v[0];
                    t[1] += v[1];
                }
                cupdate(F, fuse, f.xconj, y_ptr, alpha_r, alpha_i, t[0], t[1]);
                a_ptr = off(F, a_ptr, lda);
                y_ptr = off(F, y_ptr, inc_y);
            }
            return 0;
        }
    };
}
