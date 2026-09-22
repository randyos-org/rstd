//! AArch64 `axpy` on SVE: OpenBLAS `kernel/arm64/axpy_sve.c`, for the cores
//! whose `KERNEL` file names it (`A64FX`, `NEOVERSEV1`).
//!
//! The C is a loop of vector-length-agnostic intrinsics: full vectors under an
//! all-true predicate, then one predicated tail. Zig has no SVE intrinsics, so
//! the contiguous case is the loop the intrinsics stand for, in one `asm`
//! block: a `whilelt` predicate covers the tail without a separate step, and
//! `incd`/`incw` steps by the vector length. The strided case is scalar, as in
//! the reference.

const std = @import("std");
const types = @import("../../../types.zig");
const util = @import("../../../util.zig");

const BlasLong = types.BlasLong;

/// `y := alpha*x + y`, the reference's kernel ABI.
pub fn Real(comptime T: type) type {
    // the element suffix, its size shift, and the increment-by-vector-length
    const e = if (T == f64) "d" else "s";
    const s = if (T == f64) "3" else "2";
    const incr = if (T == f64) "incd" else "incw";
    const template = std.fmt.comptimePrint(
        \\ mov z0.{0s}, %[a:{0s}]
        \\ mov %[i], #0
        \\ whilelt p0.{0s}, %[i], %[n]
        \\ b.none 2f
        \\1:
        \\ ld1{3s} {{z1.{0s}}}, p0/z, [%[x], %[i], lsl #{1s}]
        \\ ld1{3s} {{z2.{0s}}}, p0/z, [%[y], %[i], lsl #{1s}]
        \\ fmla z2.{0s}, p0/m, z1.{0s}, z0.{0s}
        \\ st1{3s} {{z2.{0s}}}, p0, [%[y], %[i], lsl #{1s}]
        \\ {2s} %[i]
        \\ whilelt p0.{0s}, %[i], %[n]
        \\ b.first 1b
        \\2:
    , .{ e, s, incr, if (T == f64) "d" else "w" });
    return struct {
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
            if (n < 0) return 0;
            if (da == 0.0) return 0;

            if (inc_x == 1 and inc_y == 1) {
                var i: i64 = 0;
                asm volatile (template
                    : [i] "=&r" (i),
                    : [a] "w" (da),
                      [n] "r" (@as(i64, n)),
                      [x] "r" (x),
                      [y] "r" (y),
                    : .{ .memory = true, .nzcv = true, .z0 = true, .z1 = true, .z2 = true, .p0 = true });
                return 0;
            }

            var i: BlasLong = 0;
            var ix: BlasLong = 0;
            var iy: BlasLong = 0;
            while (i < n) : (i += 1) {
                util.at(T, y, iy).* += da * util.at(T, x, ix).*;
                ix += inc_x;
                iy += inc_y;
            }
            return 0;
        }
    };
}
