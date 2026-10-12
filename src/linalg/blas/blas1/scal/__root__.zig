//! Vtable wrappers for the `scal` kernels.
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

pub fn bscal_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, da: BF16, x: [*]BF16, inc_x: BlasLong, y: [*]BF16, inc_y: BlasLong, dummy: [*]BF16, dummy2: BlasLong) i32 {
    return blas.kernel(.bscal_k)(n, dummy0, dummy1, da, x, inc_x, y, inc_y, dummy, dummy2);
}
pub fn sbscal_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, da: f32, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, dummy: [*]f32, dummy2: BlasLong) i32 {
    return blas.kernel(.sbscal_k)(n, dummy0, dummy1, da, x, inc_x, y, inc_y, dummy, dummy2);
}
pub fn sscal_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, da: f32, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, dummy: [*]f32, dummy2: BlasLong) i32 {
    return blas.kernel(.sscal_k)(n, dummy0, dummy1, da, x, inc_x, y, inc_y, dummy, dummy2);
}
pub fn dscal_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, da: f64, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, dummy: [*]f64, dummy2: BlasLong) i32 {
    return blas.kernel(.dscal_k)(n, dummy0, dummy1, da, x, inc_x, y, inc_y, dummy, dummy2);
}
pub fn qscal_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, da: XDouble, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, dummy: [*]XDouble, dummy2: BlasLong) i32 {
    return blas.kernel(.qscal_k)(n, dummy0, dummy1, da, x, inc_x, y, inc_y, dummy, dummy2);
}
pub fn cscal_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, da_r: f32, da_i: f32, x: [*]f32, inc_x: BlasLong, y: [*]f32, inc_y: BlasLong, dummy: [*]f32, dummy2: BlasLong) i32 {
    return blas.kernel(.cscal_k)(n, dummy0, dummy1, da_r, da_i, x, inc_x, y, inc_y, dummy, dummy2);
}
pub fn zscal_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, da_r: f64, da_i: f64, x: [*]f64, inc_x: BlasLong, y: [*]f64, inc_y: BlasLong, dummy: [*]f64, dummy2: BlasLong) i32 {
    return blas.kernel(.zscal_k)(n, dummy0, dummy1, da_r, da_i, x, inc_x, y, inc_y, dummy, dummy2);
}
pub fn xscal_k(n: BlasLong, dummy0: BlasLong, dummy1: BlasLong, da_r: XDouble, da_i: XDouble, x: [*]XDouble, inc_x: BlasLong, y: [*]XDouble, inc_y: BlasLong, dummy: [*]XDouble, dummy2: BlasLong) i32 {
    return blas.kernel(.xscal_k)(n, dummy0, dummy1, da_r, da_i, x, inc_x, y, inc_y, dummy, dummy2);
}
