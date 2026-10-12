//! Vtable wrappers for the `swap` kernels.
//!
//! Argument names are copied from the OpenBLAS kernel sources (see the
//! `CNAME(...)` definitions under `kernel/<arch>/*.c`). Several are named
//! `dummy*` in the reference: kernel-ABI padding slots this shape does not use.

const blas = @import("../../__root__.zig");
const types = @import("../../types.zig");

const BlasLong = types.BlasLong;
const BlasInt = types.BlasInt;
const HF16 = types.HF16;
const BF16 = types.BF16;
const XDouble = types.XDouble;
const Complex32 = types.Complex32;
const Complex64 = types.Complex64;
const ComplexXD = types.ComplexXD;

pub fn sbswap_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, dummy3: f32, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, dummy: [*]f32, dummy2: BlasLong) i32 {
    return blas.kernel(.sbswap_k)(n, dummy0, dummy1, dummy3, x, inc_x, y, inc_y, dummy, dummy2);
}
pub fn sswap_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, dummy3: f32, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, dummy: [*]f32, dummy2: BlasLong) i32 {
    return blas.kernel(.sswap_k)(n, dummy0, dummy1, dummy3, x, inc_x, y, inc_y, dummy, dummy2);
}
pub fn dswap_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, dummy3: f64, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, dummy: [*]f64, dummy2: BlasLong) i32 {
    return blas.kernel(.dswap_k)(n, dummy0, dummy1, dummy3, x, inc_x, y, inc_y, dummy, dummy2);
}
pub fn qswap_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, dummy3: XDouble, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, dummy: [*]XDouble, dummy2: BlasLong) i32 {
    return blas.kernel(.qswap_k)(n, dummy0, dummy1, dummy3, x, inc_x, y, inc_y, dummy, dummy2);
}
pub fn cswap_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, dummy3: f32, dummy4: f32, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, dummy: [*]f32, dummy2: BlasLong) i32 {
    return blas.kernel(.cswap_k)(n, dummy0, dummy1, dummy3, dummy4, x, inc_x, y, inc_y, dummy, dummy2);
}
pub fn zswap_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, dummy3: f64, dummy4: f64, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, dummy: [*]f64, dummy2: BlasLong) i32 {
    return blas.kernel(.zswap_k)(n, dummy0, dummy1, dummy3, dummy4, x, inc_x, y, inc_y, dummy, dummy2);
}
pub fn xswap_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, dummy3: XDouble, dummy4: XDouble, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, dummy: [*]XDouble, dummy2: BlasLong) i32 {
    return blas.kernel(.xswap_k)(n, dummy0, dummy1, dummy3, dummy4, x, inc_x, y, inc_y, dummy, dummy2);
}
