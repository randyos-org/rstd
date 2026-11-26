//! The real gemm kernels: `beta` and the register tile over packed panels.

const types = @import("../../../types.zig");
const util = @import("../../../util.zig");
const loops = @import("../loops.zig");
const kroot = @import("../kernel.zig");
const kpack = @import("pack.zig");

const BlasLong = types.BlasLong;
const Geometry = kroot.Geometry;
const Lanes = kroot.Lanes;
const madd = kroot.madd;
const splat = kroot.splat;
const store = kroot.store;
const pack = kpack.pack;

pub fn RealKernels(comptime T: type, comptime geo: Geometry) type {
    const vl = geo.vl;
    // C for the vtable's types, Zig's own for the rest
    const cc = types.kernelCc(T);
    // `f16` sums in `f32`; vectors are only ever `f32` or `f64`, summed as is
    const A = types.Acc(T);
    if (A != T and vl > 1) @compileError("a widened sum needs scalar lanes");
    const mr = geo.mr();
    const nr = geo.nr;
    return struct {
        pub const unroll_m: i32 = mr;
        pub const unroll_n: i32 = nr;

        /// `C := beta*C` on `m x n`; a zero `beta` writes zeros without
        /// reading `C`.
        pub fn gemm_beta(m: BlasLong, n: BlasLong, _: BlasLong, b: T, _: [*]T, _: BlasLong, _: [*]T, _: BlasLong, c: [*]T, ldc: BlasLong) callconv(cc) i32 {
            const height: usize = @intCast(m);
            const ld: usize = @intCast(ldc);
            for (0..@intCast(n)) |j| {
                const col = c[j * ld ..][0..height];
                if (b == 0) {
                    @memset(col, 0);
                } else {
                    for (col) |*v| v.* = util.scalar(T, v.* * b);
                }
            }
            return 0;
        }

        /// Rows of `op(A) = A^T` (`m` deep, `n` rows), `mr` a panel.
        pub fn gemm_incopy(m: BlasLong, n: BlasLong, a: [*]T, lda: BlasLong, b: [*]T) callconv(cc) i32 {
            pack(T, mr, false, m, n, a, lda, b);
            return 0;
        }
        /// Rows of `op(A) = A`.
        pub fn gemm_itcopy(m: BlasLong, n: BlasLong, a: [*]T, lda: BlasLong, b: [*]T) callconv(cc) i32 {
            pack(T, mr, true, m, n, a, lda, b);
            return 0;
        }
        /// Columns of `op(B) = B`, `nr` a panel.
        pub fn gemm_oncopy(m: BlasLong, n: BlasLong, a: [*]T, lda: BlasLong, b: [*]T) callconv(cc) i32 {
            pack(T, nr, false, m, n, a, lda, b);
            return 0;
        }
        /// Columns of `op(B) = B^T`.
        pub fn gemm_otcopy(m: BlasLong, n: BlasLong, a: [*]T, lda: BlasLong, b: [*]T) callconv(cc) i32 {
            pack(T, nr, true, m, n, a, lda, b);
            return 0;
        }

        /// `C += alpha * A * B` for `A` packed `bm x bk` (rows in panels of
        /// `mr`) and `B` packed `bk x bn` (columns in panels of `nr`), through
        /// `loops.zig`.
        pub fn gemm_kernel(bm: BlasLong, bn: BlasLong, bk: BlasLong, alpha: T, ba: [*]T, bb: [*]T, c: [*]T, ldc: BlasLong) callconv(cc) i32 {
            L.run(@intCast(bm), @intCast(bn), @intCast(bk), alpha, ba, bb, c, @intCast(ldc));
            return 0;
        }

        const L = loops.Loops(T, geo, Out);

        /// `C += alpha * sum`, rounded back to `T`.
        const Out = struct {
            pub const C = T;
            pub const scale = 1;
            pub const Alpha = T;

            pub inline fn prep(comptime lanes: usize, alpha: T) Lanes(A, lanes) {
                return splat(A, lanes, alpha);
            }

            pub inline fn put(p: [*]T, sum: A, alpha: T) void {
                const r = madd(A, 1, sum, alpha, p[0]);
                p[0] = if (A == T) r else @floatCast(r);
            }

            pub inline fn putLanes(comptime lanes: usize, p: [*]T, sums: Lanes(A, lanes), al: Lanes(A, lanes)) void {
                const r = madd(A, lanes, sums, al, L.wide(lanes, p));
                if (A == T) store(T, lanes, p, r) else p[0] = @floatCast(r);
            }
        };
    };
}
