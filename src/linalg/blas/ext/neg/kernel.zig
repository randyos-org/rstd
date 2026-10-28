//! `neg_tcopy` over the gemm tile: `kernel/generic/neg_tcopy_*.c` and
//! `zneg_tcopy_*.c`, which are `gemm_tcopy_*.c` with every element negated.
//! Built with the `gemm_itcopy` it mirrors (`gemm/kernel.zig`'s `sets`), so
//! its layout is the installed gemm kernels' on every core. Nothing upstream
//! calls it.

const types = @import("../../types.zig");
const gk = @import("../../blas3/gemm/kernel.zig");

const BlasLong = types.BlasLong;

pub fn Kernels(comptime T: type, comptime geo: gk.Geometry) type {
    const complex = types.isComplex(T);
    const F = types.Float(T);
    const cc = types.kernelCc(F);
    const E = if (complex) [2]F else F;
    return struct {
        /// Line `j` at depth `l` is `-A(l, j)`: `gemm_itcopy` negated.
        const Source = struct {
            a: [*]const E,
            ld: usize,

            pub fn at(s: Source, j: usize, l: usize) E {
                const v = s.a[l * s.ld + j];
                return if (comptime complex) .{ -v[0], -v[1] } else -v;
            }
        };

        pub fn neg_tcopy(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, b: [*]F) callconv(cc) i32 {
            const src: Source = .{ .a = @ptrCast(a), .ld = @intCast(lda) };
            gk.packWith(E, geo.mr(), @intCast(m), @intCast(n), @ptrCast(b), src);
            return 0;
        }
    };
}
