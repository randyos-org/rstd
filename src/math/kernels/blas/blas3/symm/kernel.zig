//! The `symm` and `hemm` copies (`kernel/generic/symm_{u,l}copy_*.c`,
//! `zsymm_*copy_*.c`, `zhemm_{ut,lt}copy_*.c`): a block of the whole matrix
//! from its stored triangle, in the gemm layout (`gemm/kernel.zig`).

const types = @import("../../types.zig");
const gk = @import("../gemm/kernel.zig");

const BlasLong = types.BlasLong;

/// Element `(x + j, y + l)` of the matrix whose `upper` (or lower) triangle
/// `a` stores, conjugated for `herm` where it comes from the other side, with
/// a real diagonal; `E` is `F`, or `[2]F` for complex. 3M's copies read it too.
pub fn Source(comptime E: type) type {
    const complex = @typeInfo(E) == .array;
    return struct {
        a: [*]const E,
        ld: usize,
        x: usize,
        y: usize,
        upper: bool,
        herm: bool,

        pub fn at(s: @This(), j: usize, l: usize) E {
            const x = s.x + j;
            const y = s.y + l;
            const mirrored = if (s.upper) x > y else x < y;
            const v = if (mirrored) s.a[y + x * s.ld] else s.a[x + y * s.ld];
            if (comptime !complex) return v;
            if (!s.herm) return v;
            if (x == y) return .{ v[0], 0 };
            return if (mirrored) .{ v[0], -v[1] } else v;
        }
    };
}

pub fn Kernels(comptime T: type, comptime geo: gk.Geometry) type {
    const complex = types.isComplex(T);
    const F = types.Float(T);
    const cc = types.kernelCc(F);
    // one element as the copies move it
    const E = if (complex) [2]F else F;
    const mr = geo.mr();
    const nr = geo.nr;
    return struct {
        const Src = Source(E);

        /// One function for every copy of a precision: the options are values.
        fn copy(u: usize, upper: bool, herm: bool, m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]F) i32 {
            const src: Src = .{ .a = @ptrCast(a), .ld = @intCast(lda), .x = @intCast(posX), .y = @intCast(posY), .upper = upper, .herm = herm };
            return gk.copyFrom(E, u, m, n, b, src);
        }

        pub fn symm_iutcopy(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]F) callconv(cc) i32 {
            return copy(mr, true, false, m, n, a, lda, posX, posY, b);
        }
        pub fn symm_iltcopy(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]F) callconv(cc) i32 {
            return copy(mr, false, false, m, n, a, lda, posX, posY, b);
        }
        pub fn symm_outcopy(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]F) callconv(cc) i32 {
            return copy(nr, true, false, m, n, a, lda, posX, posY, b);
        }
        pub fn symm_oltcopy(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]F) callconv(cc) i32 {
            return copy(nr, false, false, m, n, a, lda, posX, posY, b);
        }

        pub const hemm_iutcopy = if (complex) hemmIut else {};
        pub const hemm_iltcopy = if (complex) hemmIlt else {};
        pub const hemm_outcopy = if (complex) hemmOut else {};
        pub const hemm_oltcopy = if (complex) hemmOlt else {};

        fn hemmIut(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]F) callconv(cc) i32 {
            return copy(mr, true, true, m, n, a, lda, posX, posY, b);
        }
        fn hemmIlt(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]F) callconv(cc) i32 {
            return copy(mr, false, true, m, n, a, lda, posX, posY, b);
        }
        fn hemmOut(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]F) callconv(cc) i32 {
            return copy(nr, true, true, m, n, a, lda, posX, posY, b);
        }
        fn hemmOlt(m: BlasLong, n: BlasLong, a: [*]F, lda: BlasLong, posX: BlasLong, posY: BlasLong, b: [*]F) callconv(cc) i32 {
            return copy(nr, false, true, m, n, a, lda, posX, posY, b);
        }
    };
}
