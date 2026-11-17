//! The row interchanges of `?LASWP`, and `laswp_ncopy` over the gemm tile.
//!
//! `swapRows` is `lapack/laswp/generic/laswp_k_1.c` (`laswp_plus` and
//! `laswp_minus`), as Reference-LAPACK's `dlaswp.f` defines it: rows `first` to
//! `first + count - 1` in turn, each swapped with row `ipiv(k)`; 0-based rows and
//! pivots (issue #398), where the reference counts from 1. The unrolled upstream
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

/// Swaps rows `first` to `first + count - 1` of the `n` columns of `a`, row
/// `first + r` with row `ipiv[r]`, 0-based; with `reverse` from the last back,
/// which undoes the forward order. The caller checks every index.
pub fn swapRows(comptime E: type, n: usize, a: [*]E, lda: usize, first: usize, count: usize, ipiv: [*]const usize, reverse: bool) void {
    for (0..n) |j| {
        const col = a + j * lda;
        for (0..count) |s| {
            const r = if (reverse) count - 1 - s else s;
            const i = first + r;
            const ip = ipiv[r];
            if (ip == i) continue;
            const t = col[i];
            col[i] = col[ip];
            col[ip] = t;
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
        /// reference passes the whole vector. The slot keeps OpenBLAS's kernel
        /// ABI (`Vtable.zig`): 1-based rows and pivots, here only.
        pub fn laswp_ncopy(n: BlasLong, k1: BlasLong, k2: BlasLong, a: [*]F, lda: BlasLong, ipiv: [*]isize, buffer: [*]F) callconv(cc) i32 {
            if (n <= 0 or k2 < k1) return 0;
            const ld: usize = @intCast(lda);
            const first: usize = @intCast(k1 - 1);
            const last: usize = @intCast(k2 - 1);
            const ea: [*]E = @ptrCast(a);
            for (0..@intCast(n)) |j| {
                const col = ea + j * ld;
                for (first..last + 1) |i| {
                    const ip: usize = @intCast(ipiv[i] - 1);
                    if (ip == i) continue;
                    const t = col[i];
                    col[i] = col[ip];
                    col[ip] = t;
                }
            }
            const src: Source = .{ .a = ea + first, .ld = ld };
            gk.packWith(E, geo.nr, @intCast(k2 - k1 + 1), @intCast(n), @ptrCast(buffer), src);
            return 0;
        }
    };
}
