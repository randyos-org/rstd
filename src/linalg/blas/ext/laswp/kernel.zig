//! The row interchanges of `?LASWP`, and `laswp_ncopy` over the gemm tile.
//!
//! `swapRows` is `lapack/laswp/generic/laswp_k_1.c` (`laswp_plus` and
//! `laswp_minus`), as Reference-LAPACK's `dlaswp.f` defines it: rows `k1` to
//! `k2` in turn, each swapped with row `ipiv(k)`. The unrolled upstream
//! variants give the same matrix, so this is all of them.
//!
//! `laswp_ncopy` is `kernel/generic/laswp_ncopy_*.c`: the interchanges, then
//! rows `k1` to `k2` packed as `gemm_oncopy` packs them, so it is built with
//! the gemm kernels (`blas3/kernels.zig`'s `sets`). Upstream leaves rows `k1` to
//! `k2` of `A` unswapped, having packed them; here they are swapped too, as
//! `lapack/getrf/getrf_single.c` does where it calls `laswp_plus` and
//! `gemm_oncopy` in its place.

const types = @import("../../types.zig");
const gk = @import("../../blas3/gemm/kernel.zig");

const BlasLong = types.BlasLong;

/// Swaps rows `k1` to `k2` (1-based) of the `n` columns of `a` with rows
/// `ipiv(k1)` to `ipiv(k2)`, `ipiv(k)` at `ipiv[(k - k1) * |incx|]`; for a
/// negative `incx` from `k2` back, so from the far end of `ipiv`, as the
/// reference's `IX0` has it. The caller checks every index.
pub fn swapRows(comptime E: type, n: usize, a: [*]E, lda: usize, k1: usize, k2: usize, ipiv: [*]const isize, incx: isize) void {
    if (k2 < k1) return;
    const rows = k2 - k1 + 1;
    const step: usize = @abs(incx);
    for (0..n) |j| {
        const col = a + j * lda;
        for (0..rows) |r| {
            const ix = if (incx > 0) r else rows - 1 - r;
            const i = if (incx > 0) k1 + r else k2 - r;
            const ip: usize = @intCast(ipiv[ix * step]);
            if (ip == i) continue;
            const t = col[i - 1];
            col[i - 1] = col[ip - 1];
            col[ip - 1] = t;
        }
    }
}

pub fn Kernels(comptime T: type, comptime geo: gk.Geometry) type {
    const complex = types.isComplex(T);
    const F = types.Float(T);
    const cc = types.kernelCc(F);
    const E = if (complex) [2]F else F;
    return struct {
        /// Line `j` at depth `l` is `A(l, j)`, as `gemm_oncopy` packs it.
        const Source = struct {
            a: [*]const E,
            ld: usize,

            pub fn at(s: Source, j: usize, l: usize) E {
                return s.a[j * s.ld + l];
            }
        };

        /// `ipiv` holds `ipiv(k1)` to `ipiv(k2)` from `ipiv[k1 - 1]`, as the
        /// reference passes the whole vector.
        pub fn laswp_ncopy(n: BlasLong, k1: BlasLong, k2: BlasLong, a: [*]F, lda: BlasLong, ipiv: [*]isize, buffer: [*]F) callconv(cc) i32 {
            if (n <= 0 or k2 < k1) return 0;
            const ld: usize = @intCast(lda);
            const first: usize = @intCast(k1);
            const ea: [*]E = @ptrCast(a);
            swapRows(E, @intCast(n), ea, ld, first, @intCast(k2), ipiv + (first - 1), 1);
            const src: Source = .{ .a = ea + (first - 1), .ld = ld };
            gk.packWith(E, geo.nr, @intCast(k2 - k1 + 1), @intCast(n), @ptrCast(buffer), src);
            return 0;
        }
    };
}
