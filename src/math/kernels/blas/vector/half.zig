//! The half-precision `@Vector` gemvs (`vector.zig`).

const types = @import("../types.zig");
const util = @import("../util.zig");
const bf16 = @import("../ext/bf16/kernel.zig");
const shared = @import("util.zig");

const BlasLong = types.BlasLong;
const lanes = shared.lanes;
const splatTo = shared.splatTo;

/// The half-precision gemvs (`sbgemv`, `bgemv`, `shgemv`) in
/// `ext/bf16/kernel.zig`'s bits: `f32` sums in index order, unfused, one
/// rounding on the store. Only the shape of upstream's kernels is kept
/// (`kernel/riscv64/sbgemv_?_vector.c`, `kernel/x86_64/sbgemv_n_microk_cooperlake.c`,
/// `kernel/arm64/sbgemv_n_neon.c`, `sbgemv_t_bfdot.c`, `bgemv_n_sve_v3x4.c`),
/// whose fused or paired products round otherwise.
pub fn Half(comptime E: type, comptime Y: type) type {
    return struct {
        const N = lanes(f32);
        const V = @Vector(N, f32);

        inline fn wide(e: E) f32 {
            return bf16.wide(E, e);
        }

        /// `N` elements from `p`, widened as `bf16.widen` does: a bfloat16
        /// denormal flushes to a signed zero and a NaN comes out quiet.
        inline fn widen(p: [*]const E) V {
            const w = bf16.widen(E, N, p);
            return if (N == 1) @splat(w) else w;
        }

        /// `alpha*acc + beta*y`, rounded once each; `beta == 0` does not read `y`.
        inline fn result(comptime T: type, alpha: f32, beta: f32, acc: T, yv: anytype) T {
            if (beta == 0) return splatTo(T, alpha) * acc;
            return splatTo(T, alpha) * acc + splatTo(T, beta) * @as(T, yv.*);
        }

        /// `y[at] = alpha*acc + beta*y[at]` in `f32`, rounded to `Y`.
        inline fn put(yp: *Y, alpha: f32, beta: f32, acc: f32) void {
            const wy = if (beta == 0) 0 else bf16.wide(Y, yp.*);
            yp.* = bf16.narrow(Y, result(f32, alpha, beta, acc, &wy));
        }

        /// Rows `r` to `r + B` from column 0 on, `B / N` vectors.
        fn rowBlock(comptime B: usize, r: usize, cols: usize, alpha: f32, a: [*]E, ld_a: usize, x: [*]E, incx: BlasLong, beta: f32, y: [*]Y, incy: BlasLong) void {
            var acc: [B / N]V = @splat(@splat(0));
            for (0..cols) |q| {
                const xv: V = @splat(wide(util.at(E, x, @as(BlasLong, @intCast(q)) * incx).*));
                inline for (&acc, 0..) |*s, v| s.* = s.* + widen(a + r + v * N + q * ld_a) * xv;
            }
            inline for (acc, 0..) |s, v| {
                const r0 = r + v * N;
                if (Y == f32 and incy == 1) {
                    const yp: *[N]f32 = (y + r0)[0..N];
                    yp.* = result(V, alpha, beta, s, yp);
                } else {
                    const sa: [N]f32 = s;
                    inline for (sa, 0..) |e, l| {
                        put(util.at(Y, y, @as(BlasLong, @intCast(r0 + l)) * incy), alpha, beta, e);
                    }
                }
            }
        }

        /// `y = alpha*A*x + beta*y`, a lane per row of `A`.
        pub fn gemv_n(m: BlasLong, n: BlasLong, alpha_y: Y, a: [*]E, lda: BlasLong, x: [*]E, incx: BlasLong, beta_y: Y, y: [*]Y, incy: BlasLong) callconv(.c) i32 {
            const alpha = bf16.wide(Y, alpha_y);
            const beta = bf16.wide(Y, beta_y);
            const rows: usize = @intCast(@max(m, 0));
            const cols: usize = @intCast(@max(n, 0));
            const ld_a: usize = @intCast(lda);
            var r: usize = 0;
            while (r + 4 * N <= rows) : (r += 4 * N) rowBlock(4 * N, r, cols, alpha, a, ld_a, x, incx, beta, y, incy);
            while (r + N <= rows) : (r += N) rowBlock(N, r, cols, alpha, a, ld_a, x, incx, beta, y, incy);
            while (r < rows) : (r += 1) {
                var acc: f32 = 0;
                for (0..cols) |q| acc = acc + wide(a[r + q * ld_a]) * wide(util.at(E, x, @as(BlasLong, @intCast(q)) * incx).*);
                put(util.at(Y, y, @as(BlasLong, @intCast(r)) * incy), alpha, beta, acc);
            }
            return 0;
        }

        /// `y = alpha*A'*x + beta*y`, a lane per column of `A`: `N` columns'
        /// `N` rows at a time, transposed, so each column sums in row order.
        pub fn gemv_t(m: BlasLong, n: BlasLong, alpha_y: Y, a: [*]E, lda: BlasLong, x: [*]E, incx: BlasLong, beta_y: Y, y: [*]Y, incy: BlasLong) callconv(.c) i32 {
            const alpha = bf16.wide(Y, alpha_y);
            const beta = bf16.wide(Y, beta_y);
            const rows: usize = @intCast(@max(m, 0));
            const cols: usize = @intCast(@max(n, 0));
            const ld_a: usize = @intCast(lda);
            var c: usize = 0;
            while (c + N <= cols) : (c += N) {
                var acc: V = @splat(0);
                var q: usize = 0;
                while (q + N <= rows) : (q += N) {
                    var t: [N]V = undefined;
                    inline for (&t, 0..) |*v, l| v.* = widen(a + q + (c + l) * ld_a);
                    t = transpose(t);
                    inline for (t, 0..) |v, l| acc = acc + v * @as(V, @splat(wide(util.at(E, x, @as(BlasLong, @intCast(q + l)) * incx).*)));
                }
                while (q < rows) : (q += 1) {
                    var row: [N]f32 = undefined;
                    inline for (&row, 0..) |*e, l| e.* = wide(a[q + (c + l) * ld_a]);
                    acc = acc + @as(V, row) * @as(V, @splat(wide(util.at(E, x, @as(BlasLong, @intCast(q)) * incx).*)));
                }
                const sa: [N]f32 = acc;
                inline for (sa, 0..) |e, l| {
                    put(util.at(Y, y, @as(BlasLong, @intCast(c + l)) * incy), alpha, beta, e);
                }
            }
            while (c < cols) : (c += 1) {
                var acc: f32 = 0;
                for (0..rows) |q| acc = acc + wide(a[q + c * ld_a]) * wide(util.at(E, x, @as(BlasLong, @intCast(q)) * incx).*);
                put(util.at(Y, y, @as(BlasLong, @intCast(c)) * incy), alpha, beta, acc);
            }
            return 0;
        }

        /// `t` transposed: halves of rows trade their off-diagonal halves,
        /// then quarters, down to single lanes.
        inline fn transpose(t: [N]V) [N]V {
            var v = t;
            comptime var h = N / 2;
            inline while (h >= 1) : (h /= 2) {
                const lo = comptime mask(h, false);
                const hi = comptime mask(h, true);
                inline for (0..N) |r| {
                    if (r & h != 0) continue;
                    const p = v[r];
                    const q = v[r + h];
                    v[r] = @shuffle(f32, p, q, lo);
                    v[r + h] = @shuffle(f32, p, q, hi);
                }
            }
            return v;
        }

        /// The upper (`second`) or lower row of a pair after trading the
        /// lanes with bit `h` set in the first for those clear in the second.
        fn mask(comptime h: usize, comptime second: bool) [N]i32 {
            var m: [N]i32 = undefined;
            for (0..N) |l| {
                const own = l & h == 0;
                m[l] = if (!second)
                    (if (own) @intCast(l) else ~@as(i32, @intCast(l - h)))
                else
                    (if (own) @intCast(l + h) else ~@as(i32, @intCast(l)));
            }
            return m;
        }
    };
}
