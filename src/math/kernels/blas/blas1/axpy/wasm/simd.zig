//! WebAssembly `axpy` on SIMD128: OpenBLAS `kernel/wasm/axpy.c`, for the core
//! whose `KERNEL` file names it (`WASM128_GENERIC`); the complex kernels there
//! are the plain C loops, which the portable kernels are.
//!
//! The unit-stride loop keeps eight independent vectors in flight so the
//! loads, multiplies, adds and stores can overlap, then a vector at a time,
//! then scalars. A multiply and an add, not a fused one: WebAssembly has none
//! outside the relaxed extension, and the reference notes that `axpy` is
//! checked to machine epsilon. Strided access stays scalar, since SIMD128 has
//! no gather.

const types = @import("../../../types.zig");
const util = @import("../../../util.zig");
const kernel = @import("../driver.zig");

const BlasLong = types.BlasLong;

/// `y := alpha*x + y`, the reference's kernel ABI.
pub fn Real(comptime T: type) type {
    // a 128-bit vector
    const vlen = 16 / @sizeOf(T);
    const V = @Vector(vlen, T);
    const unroll = 8;
    const chunk = vlen * unroll;
    return struct {
        /// The contiguous case.
        fn unit(len: usize, x: [*]const T, y: [*]T, da: T) void {
            const va: V = @splat(da);
            var i: usize = 0;
            const n_main = len & ~@as(usize, chunk - 1);
            while (i < n_main) : (i += chunk) {
                var xs: [unroll]V = undefined;
                var ys: [unroll]V = undefined;
                inline for (0..unroll) |k| xs[k] = util.load(V, x + i + k * vlen);
                inline for (0..unroll) |k| ys[k] = util.load(V, y + i + k * vlen);
                inline for (0..unroll) |k| util.store(V, y + i + k * vlen, ys[k] + va * xs[k]);
            }
            while (i + vlen <= len) : (i += vlen) {
                util.store(V, y + i, util.load(V, y + i) + va * util.load(V, x + i));
            }
            while (i < len) : (i += 1) y[i] += da * x[i];
        }

        pub fn call(
            n: BlasLong,
            _: BlasLong,
            _: BlasLong,
            da: T,
            x: [*]T,
            inc_x: BlasLong,
            y: [*]T,
            inc_y: BlasLong,
            _: [*]T,
            _: BlasLong,
        ) callconv(.c) i32 {
            if (n <= 0) return 0;
            const len: usize = @intCast(n);

            if (inc_x == 1 and inc_y == 1) {
                unit(len, x, y, da);
                return 0;
            }

            // scalar, four at a time as the reference does: the driver's loop
            kernel.Real(T, 1, null, false).call(n, da, x, inc_x, y, inc_y);
            return 0;
        }
    };
}
